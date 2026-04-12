#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"

using namespace hcm;

// ============================================================================
// TageAhead — Ahead-pipelined TAGE with SINGLE combined RAM.
//
// Architecture: ONE RAM holding all tables' tag+pred data per entry.
// All tables share the same lineaddr-based index. Per-table history folding
// goes into tags only (not index). This ensures a single RAM read per
// predict1, keeping the floorplan compact → P2 < 1 cycle.
//
// Pipeline pattern (gshareN_ahead_best):
//   execute_if(true_block) {
//     shift data[0] → data[1]          // fast: reg→reg, old [0] value
//     data[0] = combined_ram.read(idx)  // RAM read for current block
//     tag compare on data[1] → pred     // prediction from previous read
//   }
//   pred[R] = ...  // outside execute_if (unconditional)
// ============================================================================

namespace ta {

constexpr u64 clog2(u64 x) {
  u64 r = 0, v = x - 1;
  while (v > 0) { v >>= 1; r++; }
  return r;
}

constexpr u64 next_pow2(u64 x) {
  u64 r = 1;
  while (r < x) r <<= 1;
  return r;
}

constexpr double constexpr_pow(double base, double exp) {
  if (exp == 0.0) return 1.0;
  if (base == 0.0) return 0.0;
  return std::exp(exp * std::log(base));
}

template <std::size_t N>
constexpr std::array<u64, N> geometric_hist(u64 min_h, u64 max_h) {
  static_assert(N >= 2);
  std::array<u64, N> h{};
  u64 prev = 0;
  double ratio = static_cast<double>(max_h) / static_cast<double>(min_h);
  for (std::size_t i = 0; i < N; i++) {
    double e = static_cast<double>(i) / static_cast<double>(N - 1);
    u64 hl = static_cast<u64>(static_cast<double>(min_h) * constexpr_pow(ratio, e));
    hl = (hl > prev + 1) ? hl : prev + 1;
    h[N - 1 - i] = hl;
    prev = hl;
  }
  return h;
}

// XOR Fold (incremental folded history hash)
template <u64 F> struct XORFold {
  reg<F> folded;

  val<F> get() { return folded; }
  void fanout(auto fo) { folded.fanout(fo); }

  template <u64 MAXL> void update(auto &gh, u64 ghlen, auto in) {
    constexpr u64 inbits = std::min(F, u64(decltype(in)::size));
    val<inbits> input = in;
    auto f = folded.make_array(val<1>{});
    val<1> outbit = gh[ghlen - 1];
    u64 outpos = ghlen % F;
    arr<val<1>, F> ff = [&](u64 i) {
      if (i == 0)
        return (outpos == 0) ? f[F - 1].fo1() ^ outbit.fo1() : f[F - 1].fo1();
      else
        return (outpos == i) ? f[i - 1].fo1() ^ outbit.fo1() : f[i - 1].fo1();
    };
    auto x = input.fo1().make_array(val<1>{});
    arr<val<1>, F> y = [&](u64 i) {
      return (i < x.size) ? x[i].fo1() ^ ff[i].fo1() : ff[i].fo1();
    };
    folded = y.fo1().concat();
  }
};

// Global History Register
template <u64 MAXLEN> struct TAGlobalHistory {
  arr<reg<1>, MAXLEN> hist;

  val<1> operator[](u64 pos) { return hist[pos]; }

  void update(auto in) {
    auto input = in.fo1().make_array(val<1>{});
    static_assert(input.size <= MAXLEN);
    for (u64 i = MAXLEN - 1; i >= input.size; i--)
      hist[i] = hist[i - 1];
    for (u64 i = input.size - 1; i >= 1; i--)
      hist[i] = hist[i - 1] ^ input[i].fo1();
    hist[0] = input[0].fo1();
  }

  void fanout(auto fo) {
    for (u64 i = 0; i < MAXLEN; i++) hist[i].fanout(fo);
  }
};

// History Folder
template <u64 NUM_TABLES, u64 MAXHIST, u64 TAG_FOLD_W,
          std::array<u64, NUM_TABLES> HIST_LEN = std::array<u64, NUM_TABLES>{}>
struct TAHistoryFolder {
  TAGlobalHistory<MAXHIST> gh;
  std::array<XORFold<TAG_FOLD_W>, NUM_TABLES> tag_folds;

  val<TAG_FOLD_W> get_tag_fold(u64 i) { return tag_folds[i].get(); }

  void fanout(auto fo) {
    static_loop<NUM_TABLES>([&]<u64 I>() {
      tag_folds[I].fanout(fo);
    });
  }

  void update(auto branchbits) {
    branchbits.fanout(hard<NUM_TABLES + 1>{});
    gh.fanout(hard<std::max(u64(2), NUM_TABLES + 1)>{});
    static_loop<NUM_TABLES>([&]<u64 I>() {
      static constexpr u64 HL = HIST_LEN[I];
      tag_folds[I].template update<MAXHIST>(gh, HL, branchbits);
    });
    gh.update(branchbits);
  }
};

} // namespace ta

// ============================================================================
// TageAheadImpl — Single combined RAM architecture
// ============================================================================

template <u64 NUM_TABLES_V = 8, u64 TAG_WIDTH_V = 11, u64 TABLE_SIZE_V = 4096,
          u64 LINEINST_V = 1024, u64 N_V = 7, u64 PATH_BITS_V = 6,
          u64 MINH_V = 8, u64 MAXH_V = 400>
struct TageAheadImpl : predictor {

  // ======== Constants ========
  static constexpr u64 NUM_TABLES = NUM_TABLES_V;
  static constexpr u64 TAG_WIDTH = TAG_WIDTH_V;
  static constexpr u64 TABLE_SIZE = TABLE_SIZE_V;
  static constexpr u64 N = N_V;
  static constexpr u64 LANES = ta::next_pow2(N);
  static constexpr u64 LOG_LANES = ta::clog2(LANES);
  static constexpr u64 LINEINST = LINEINST_V;
  static constexpr u64 LOG_LINEINST = ta::clog2(LINEINST);
  static constexpr u64 PATHBITS = PATH_BITS_V;
  static constexpr u64 MAXHIST = MAXH_V;
  static constexpr u64 IDX_BITS = ta::clog2(TABLE_SIZE);
  static constexpr u64 HTAG_BITS = TAG_WIDTH - LOG_LANES; // tag minus lane bits

  static constexpr auto HIST_LEN = ta::geometric_hist<NUM_TABLES>(MINH_V, MAXH_V);

  // Combined entry: per-table [tag TAG_WIDTH | pred LANES]
  static constexpr u64 ENTRY_PER_TABLE = TAG_WIDTH + LANES;
  // Total entry across all tables
  static constexpr u64 LOGBANKS = std::bit_width(NUM_TABLES - 1);
  static constexpr u64 BANKS = 1 << LOGBANKS;

  static_assert(NUM_TABLES > 0 && NUM_TABLES <= BANKS);
  static_assert(N >= 1 && LANES >= N);
  static_assert(std::has_single_bit(LINEINST));
  static_assert(TAG_WIDTH > LOG_LANES, "TAG_WIDTH must be > LOG_LANES");

  // ======== Storage ========

  // SINGLE combined RAM: each entry holds tag+pred for all tables.
  // Entry per bank: [tag TAG_WIDTH | pred LANES] = ENTRY_PER_TABLE bits.
  // BANKS banks (one per table, padded to power of 2).
  hcm::ram<arr<val<ENTRY_PER_TABLE>, BANKS>, TABLE_SIZE> combined_hi{"ta_combined"};

  // Pipeline registers: depth 2
  // data[0/1]: per-table combined tag+pred
  arr<reg<ENTRY_PER_TABLE>, BANKS> data[2];
  reg<IDX_BITS> data_idx[2]; // stored index for writes

  // History
  ta::TAHistoryFolder<NUM_TABLES, MAXHIST, HTAG_BITS, HIST_LEN> gfolds;
  bool gfolds_inited = false;

  // Computed htag per table (from folded history, stored for comparison)
  arr<reg<HTAG_BITS>, NUM_TABLES> data_htag[2];

  // Prediction state
  reg<LOGBANKS> path;
  reg<LANES> unordered_pred;
  arr<reg<1>, LANES> pred;

  // Block tracking (gshareN_ahead_best pattern)
  reg<1> true_block = 1;
  reg<1> last_condbr_dir = 1;
  reg<LOG_LINEINST> block_entry;
  reg<LANES> XL;
  reg<N + 1> rank;
  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<1>, N> branch_dir;

  bool params_printed = false;

  // ======== Helpers ========
  val<1> line_end() { return (block_entry + block_size) >= hard<LINEINST>{}; }
  val<1> last_pred() { return rank >> (N - num_branch); }

  // ======== Predictor Interface ========

  val<1> predict1(val<64> inst_pc) override {
    if (!gfolds_inited) { gfolds_inited = true; }

    if (!params_printed) {
      std::cerr << "=== TageAhead (Single-RAM) ===\n"
                << "NUM_TABLES=" << NUM_TABLES << "  TABLE_SIZE=" << TABLE_SIZE
                << "  TAG=" << TAG_WIDTH << "  LANES=" << LANES
                << "  N=" << N << "  LINEINST=" << LINEINST << "\n"
                << "Hist lens: ";
      for (u64 i = 0; i < NUM_TABLES; i++)
        std::cerr << (i ? "," : "") << HIST_LEN[i];
      std::cerr << "\n\n";
      params_printed = true;
    }

    inst_pc.fanout(hard<4>{});
    true_block.fanout(hard<8 + BANKS * 2>{});

    block_entry = select(true_block, val<LOG_LINEINST>{inst_pc >> 2},
                         val<LOG_LINEINST>{block_entry + block_size});
    block_entry.fanout(hard<LINEINST + N + 1>{});

    rank = select(true_block, val<N + 1>{1}, rank << num_branch);
    rank.fanout(hard<N + 2>{});

    XL = select(true_block, val<LOG_LANES>{inst_pc >> 6}.decode().concat(),
                XL.rotate_left(num_branch));
    XL.fanout(hard<LANES>{});

    gfolds.fanout(hard<2>{});

    execute_if(true_block, [&]() {
      // 1. Pipeline shift
      data[1] = data[0].fo1();
      data_idx[1] = data_idx[0].fo1();
      for (u64 i = 0; i < NUM_TABLES; i++)
        data_htag[1][i] = data_htag[0][i].fo1();

      // 2. Compute index & read SINGLE combined RAM
      val<IDX_BITS> lineaddr = inst_pc >> 2;
      lineaddr.fanout(hard<3>{});
      data_idx[0] = lineaddr;
      data[0] = combined_hi.read(lineaddr);

      // 3. Compute per-table htags from folded history
      static_loop<NUM_TABLES>([&]<u64 I>() {
        data_htag[0][I] =
            val<HTAG_BITS>{lineaddr}.reverse() ^ gfolds.get_tag_fold(I);
      });

      // 4. TEMPORARY: simple bank select on prediction portion of data[1]
      // Extract pred bits (low LANES bits of each combined entry)
      // For timing test: just use bank select like gshareN_ahead_best
      path = val<LOGBANKS>{num_branch} + ~last_condbr_dir;
      // Select bank → get combined entry (TAG_WIDTH + LANES bits)
      // Extract prediction portion (low LANES bits)
      unordered_pred = data[1].select(path);
      unordered_pred.fanout(hard<LANES>{});
    });

    // Lane scatter (outside execute_if)
    for (u64 i = 0; i < LANES; i++) {
      pred[i] = (unordered_pred & XL.rotate_left(i)) != hard<0>{};
    }
    pred.fanout(hard<LINEINST * 2>{});

    block_size = 1;
    num_branch = 0;
    reuse_prediction(~line_end());
    return pred[0];
  }

  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) override {
    block_size++;
    reuse_prediction(~line_end());
    return pred[num_branch];
  }

  val<1> predict2([[maybe_unused]] val<64> inst_pc) override {
    return pred[num_branch];
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) override {
    return pred[num_branch];
  }

  void update_condbr([[maybe_unused]] val<64> branch_pc, val<1> taken,
                     [[maybe_unused]] val<64> next_pc) override {
    assert(num_branch < N);
    branch_dir[num_branch] = taken.fo1();
    num_branch++;
    reuse_prediction(~(line_end() | last_pred()));
  }

  void update_cycle(instruction_info &block_end_info) override {
    val<1> &mispredict = block_end_info.is_mispredict;
    val<64> &next_pc = block_end_info.next_pc;

    if (num_branch == 0) {
      last_condbr_dir = 0;
      execute_if(true_block, [&]() {
        next_pc.fanout(hard<2>{});
        gfolds.update(val<PATHBITS>{next_pc >> 2});
      });
      true_block = 1;
      num_branch = 0;
      return;
    }

    branch_dir.fanout(hard<3>{});
    last_condbr_dir = branch_dir[num_branch - 1];

    need_extra_cycle(mispredict);

    // No writes for Phase 1 timing test — just update true_block + history
    true_block = arr<val<1>, 4>{~mispredict, branch_dir[num_branch - 1],
                                last_pred(), line_end()}
                     .fold_or();
    true_block.fanout(hard<MAXHIST + 4>{});
    execute_if(true_block, [&]() {
      next_pc.fanout(hard<2>{});
      gfolds.update(val<PATHBITS>{next_pc >> 2});
    });

    num_branch = 0;
  }
};

// User-facing alias
template <u64 NUM_TABLES = 8, u64 TAG_WIDTH = 11, u64 TABLE_SIZE = 4096,
          u64 TA_LINEINST = 1024, u64 TA_N = 7, u64 TA_PATH_BITS = 6,
          u64 TA_MINH = 8, u64 TA_MAXH = 400>
using TageAhead = TageAheadImpl<NUM_TABLES, TAG_WIDTH, TABLE_SIZE,
                                TA_LINEINST, TA_N, TA_PATH_BITS,
                                TA_MINH, TA_MAXH>;
