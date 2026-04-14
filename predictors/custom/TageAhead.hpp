#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"
#include "../common.hpp"
#include "custom_common.hpp"

using namespace hcm;

// ============================================================================
// Table Config
// ============================================================================

template <u64 N = 4, u64 SIZE = 1024, u64 TAG = 11, u64 MINH = 8,
          u64 MAXH = 100, u64 SIZE_RATIO = 1,
          ta::HistSeries HIST = ta::HistSeries::GEOMETRIC,
          typename TagFn = ta::GradedTag<TAG, TAG - 1>,
          typename SizeFn = ta::GeoSize<SIZE, SIZE_RATIO>>
struct TATableConfig {
  static constexpr u64 NUM_TABLES = N;
  static constexpr u64 MINHIST = MINH;
  static constexpr u64 MAXHIST = MAXH;
  static constexpr auto TABLE_SIZE = ta::generate_table_sizes<N>(SizeFn{});
  static constexpr auto TAG_WIDTH = ta::generate_tag_widths<N>(TagFn{});
  static constexpr auto HIST_LEN = []() {
    if constexpr (HIST == ta::HistSeries::GEOMETRIC)
      return ta::geometric_hist<N>(MINH, MAXH);
    else if constexpr (HIST == ta::HistSeries::QUADRATIC)
      return ta::quadratic_hist<N>(MINH, MAXH);
    else if constexpr (HIST == ta::HistSeries::SUPEREXP)
      return ta::superexp_hist<N>(MINH, 0.1, 1.1);
    else if constexpr (HIST == ta::HistSeries::ROS)
      return ta::ros_hist<N>(MINH, MAXH);
  }();
};

// ============================================================================
// TageAhead — Ahead-pipelined TAGE predictor
//
// P1 = P2: both return the same reg-based prediction (P2 < 1 cycle).
// predict1 reads TAGE tables for the NEXT block (ahead pipeline).
// predict2/predict1 return prediction from regs written in previous predict1.
// Fallback on secondary tag miss: bimodal (ahead-pipelined).
//
// Self-contained: no dependency on common.hpp or TageTable.hpp.
// Follows gshareN_ahead_best pipeline pattern.
// ============================================================================

template <typename TableCfg = TATableConfig<>,
          u64 N = 8,                   // max conditional branches per block
          u64 PATHBITS = 6,            // bits of next_pc injected into history
          u64 SEC_TAG_BITS = 3,        // secondary tag width (ahead ambiguity)
          u64 CTR_WIDTH = 1,           // prediction counter width per lane
          u64 HYST_WIDTH = 2,          // hysteresis width (separate from ctr)
          u64 U_WIDTH = 1,             // usefulness counter width
          u64 BIMODAL_CAPACITY = 8192, // bimodal fallback table size
          u64 META_WIDTH = 4,          // meta counter width (provider vs alt)
          u64 META_CAPACITY = 256,     // meta table entries
          u64 META_PIPE = 2,           // meta pipeline depth
          u64 LINEINST = 128,          // line size in instructions
          HistUpdate HIST_MODE =
              HistUpdate::PATH // what goes into history: PATH, DIR, or BOTH
          >
struct TageAhead : predictor {

  // ======== Derived Constants ========

  static constexpr u64 NT = TableCfg::NUM_TABLES;
  static constexpr u64 LOGLINEINST = ta::clog2(LINEINST);
  static constexpr u64 MAXHIST = TableCfg::MAXHIST;
  static constexpr u64 MAX_TAG_WIDTH = ta::array_max(TableCfg::TAG_WIDTH);
  static constexpr u64 MAX_TABLE_SIZE = ta::array_max(TableCfg::TABLE_SIZE);
  static constexpr u64 MAX_IDX_BITS = ta::clog2(MAX_TABLE_SIZE);

  // Prediction bits per entry: one CTR_WIDTH counter per branch
  static constexpr u64 PRED_BITS = N * CTR_WIDTH;

  // ======================================================================
  // Storage
  // ======================================================================
  // ---- Table tuple (per-table tag width and table size) ----
  using Tables = typename TAMakeTableTuple<TableCfg, CTR_WIDTH, HYST_WIDTH,
                                           U_WIDTH, SEC_TAG_BITS, N,
                                           std::make_index_sequence<NT>>::type;
  Tables tables;
  hcm::ram<val<N>, BIMODAL_CAPACITY> bim_ctr{"bim"};
  hcm::ram<val<META_WIDTH>, META_CAPACITY> meta_ctr{"meta"};

  // ======================================================================
  // Registers
  // ======================================================================
  // ---- Global History (shared, folds live in per-table TATable) ----
  global_history<MAXHIST> gh;

  // ---- Bimodal Fallback (ahead-pipelined) ----

  // ---- Block Tracking ----
  reg<1> true_block = 1;
  reg<1> last_condbr_dir = 1;

  // track where we enter the block
  reg<LOGLINEINST> block_entry;

  // ---- Simulation Artifacts (free in hardware) ----
  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<1>, N> branch_dir;

  // ---- Secondary tag (precomputed in update_cycle from next_pc) ----
  // update_cycle for block B computes curr_sec_tag = hash(next_pc).
  // next_pc is the PC of block B+1, which is exactly the block whose
  // prefetched data gets shifted into current_* in the same update_cycle.
  // So when predict1 for B+1 runs, curr_sec_tag already identifies
  // which successor path the current_* entries were trained for.
  reg<SEC_TAG_BITS> curr_sec_tag;

  // ---- Pipeline Regs [NT] ----
  // prefetch_*: written in predict1 (ahead reads for next block)
  reg<MAX_TAG_WIDTH> prefetch_tag[NT];
  reg<1> prefetch_tag_hit[NT]; // primary tag match (computed off crit path)
  reg<PRED_BITS> prefetch_pred[NT];
  reg<SEC_TAG_BITS> prefetch_sec[NT];
  reg<MAX_IDX_BITS> prefetch_idx[NT];
  reg<std::max(u64(1), HYST_WIDTH)> prefetch_hyst[NT];
  reg<U_WIDTH> prefetch_u[NT];

  // current_*: shifted from prefetch, used for prediction
  reg<MAX_TAG_WIDTH> current_tag[NT];
  reg<1> current_tag_hit[NT];
  reg<PRED_BITS> current_pred[NT];
  reg<SEC_TAG_BITS> current_sec[NT];
  reg<MAX_IDX_BITS> current_idx[NT];
  reg<std::max(u64(1), HYST_WIDTH)> current_hyst[NT];
  reg<U_WIDTH> current_u[NT];

  // ---- Prediction reg (shared by
  // predict1/reuse_predict1/predict2/reuse_predict2) ----
  arr<reg<1>, N> pred;

  // ---- Meta pipeline (shifted each update_cycle) ----
  reg<META_WIDTH, i64> meta_pipe[META_PIPE];

  // ======================================================================
  // Helpers
  // ======================================================================

  val<1> line_end() { return (block_entry + block_size) == hard<LINEINST>{}; }

  // ======================================================================
  // predict1/2, reuse_predict1/2, update_condbr, update_cycle — TODO
  // ======================================================================

  val<1> predict1([[maybe_unused]] val<64> inst_pc) {
    // Ahead reads for next block (off crit path, needs inst_pc)
    execute_if(true_block, [&]() {
      static_loop<NT>([&]<u64 I>() {
        auto &t = std::get<I>(tables);
        auto idx = t.fold_idx.get() ^ val<t.IDX_BITS>{inst_pc >> 2};
        auto computed_tag = t.fold_tag.get() ^ val<t.tag_width>{inst_pc >> 4};

        auto stored_tag = t.tag_ram.read(idx);
        stored_tag.fanout(hard<2>{});
        prefetch_tag[I] = stored_tag;
        prefetch_tag_hit[I] =
            val<MAX_TAG_WIDTH>{stored_tag} == val<MAX_TAG_WIDTH>{computed_tag};
        prefetch_pred[I] = t.pred_ram.read(idx);
        prefetch_sec[I] = t.sec_ram.read(idx);
        prefetch_idx[I] = idx;
        prefetch_hyst[I] = t.hyst_ram.read(idx);
        prefetch_u[I] = t.u_ram.read(idx);
      });
    });

    // Crit path: just read precomputed prediction from reg
    block_size = 1;
    num_branch = 0;
    reuse_prediction(~line_end());
    return pred[num_branch];
  }
  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) {
    block_size++;
    reuse_prediction(~line_end());
    return pred[num_branch];
  }

  val<1> predict2([[maybe_unused]] val<64> inst_pc) { return pred[num_branch]; }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) {
    return pred[num_branch];
  }

  void update_condbr([[maybe_unused]] val<64> branch_pc, val<1> taken,
                     [[maybe_unused]] val<64> next_pc) {
    assert(num_branch < N);
    branch_dir[num_branch] = taken.fo1();
    num_branch++;
    reuse_prediction(~line_end() & val<1>{num_branch < N});
  }

  void update_cycle([[maybe_unused]] instruction_info &block_end_info) {
    // ---- Prefetch part ------
    // 1. Pipeline shift: prefetch → current (unconditional — only true blocks
    // reach here)
    static_loop<NT>([&]<u64 I>() {
      current_tag[I] = prefetch_tag[I];
      current_tag_hit[I] = prefetch_tag_hit[I];
      current_pred[I] = prefetch_pred[I];
      current_sec[I] = prefetch_sec[I];
      current_idx[I] = prefetch_idx[I];
      current_hyst[I] = prefetch_hyst[I];
      current_u[I] = prefetch_u[I];
    });

    // Precompute secondary tag for next block
    curr_sec_tag = val<SEC_TAG_BITS>{block_end_info.next_pc >> 2};

    // 2. Provider / altpred resolution (shortest → longest, longest wins)
    static constexpr u64 PROV_IDX_BITS = std::max(u64(1), ta::clog2(NT));
    val<PRED_BITS> provider_pred = val<PRED_BITS>{0};
    val<PRED_BITS> alt_pred = val<PRED_BITS>{0};
    val<1> has_provider = val<1>{0};
    val<1> has_alt = val<1>{0};
    val<1> provider_weak = val<1>{0};
    val<PROV_IDX_BITS> provider_idx = val<PROV_IDX_BITS>{0};

    static_loop<NT>([&]<u64 I>() {
      constexpr u64 J = NT - 1 - I;
      val<1> sec_hit =
          val<SEC_TAG_BITS>{current_sec[J]} == val<SEC_TAG_BITS>{curr_sec_tag};
      val<1> full_hit = current_tag_hit[J] & sec_hit;

      // on hit: previous provider becomes alt, this becomes provider
      has_alt = select(full_hit, has_provider, has_alt);
      alt_pred = select(full_hit, provider_pred, alt_pred);
      provider_pred =
          select(full_hit, val<PRED_BITS>{current_pred[J]}, provider_pred);
      provider_weak =
          select(full_hit,
                 val<1>{val<std::max(u64(1), HYST_WIDTH)>{current_hyst[J]} ==
                        hard<0>{}} &
                     val<1>{val<U_WIDTH>{current_u[J]} == hard<0>{}},
                 provider_weak);
      provider_idx =
          select(full_hit, val<PROV_IDX_BITS>{hard<J>{}}, provider_idx);
      has_provider = has_provider | full_hit;
    });

    // 3. Bimodal fallback read
    auto bim_idx =
        val<ta::clog2(BIMODAL_CAPACITY)>{block_end_info.next_pc >> 2};
    val<N> bim_pred = bim_ctr.read(bim_idx);

    // Alt falls back to bimodal if no second TAGE hit
    alt_pred = select(has_alt, alt_pred, val<PRED_BITS>{bim_pred});

    // 4. Meta read + pipeline shift
    auto meta_idx = val<ta::clog2(META_CAPACITY)>{block_end_info.next_pc >> 2};
    meta_pipe[0] = meta_ctr.read(meta_idx);
    for (u64 i = META_PIPE - 1; i > 0; i--)
      meta_pipe[i] = meta_pipe[i - 1];
    val<1> meta_use_alt =
        val<META_WIDTH, i64>{meta_pipe[META_PIPE - 1]} >= hard<0>{};

    // 5. Final prediction: provider, alt, or bimodal
    // newly_alloc = weak hyst AND u==0
    // If provider is newly allocated and meta says use alt → alt_pred
    // If no provider → bim_pred (already in alt_pred via fallback)
    val<1> use_alt = provider_weak & meta_use_alt & has_alt;
    val<PRED_BITS> tage_pred = select(use_alt, alt_pred, provider_pred);
    val<PRED_BITS> final_pred =
        select(has_provider, tage_pred, val<PRED_BITS>{bim_pred});

    static_loop<N>([&]<u64 I>() { pred[I] = final_pred >> I; });
  }

  // ---- resolution part ----
};
