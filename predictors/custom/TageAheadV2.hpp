#pragma once
#include "../../cbp.hpp"
#include "TageConfig.hpp"
#include <cassert>

#pragma GCC diagnostic ignored "-Wunused-parameter"

// ============================================================================
// TageAheadV2 — Predecessor-indexed ahead TAGE with gshareN_ahead_best tricks
// ============================================================================
//
// Key differences from TageAhead v1:
//   - No separate P1 gshare. predict1 = predict2 = same precomputed result.
//   - LINEINST >> FETCH_WIDTH: blocks span many fetches, fewer RAM reads.
//   - N branches per block (not just 1). Rank tracking.
//   - Per-branch history shift (not per-block).
//   - Additive path encoding.
//   - Non-true-block continuation.
//   - All tag compare + match + prediction precomputed in predict1.
//   - predict2 = reg read. P2 ceil = 1.
//
// ============================================================================

template <
    // TAGE table params
    typename TableCfg = SweepTableConfig<4, 256>,
    // Block params
    u64 FETCH_WIDTH = 16,
    u64 LOGLINEINST = 10,                 // log2(max instructions per block). 10 = 1024.
    u64 N = 8,                            // max conditional branches per block
    // Bimodal
    u64 BIMODAL_SIZE = 4096,
    // Decay
    u64 DECAY_CTR = 8,
    u64 DECAY_GRAN = 2,
    typename DecayPolicy = DecayMild,
    // Meta
    bool USE_META = true,
    u64 METABITS = 4,
    u64 METAPIPE = 2,
    // Ahead banking
    u64 AHEAD_BANKS = 2,
    u64 AHEAD_BANK_BITS = 1>

struct TageAheadV2 : predictor {

  // ======== Constants ========
  static constexpr u64 NUM_TABLES      = TableCfg::NUM_TABLES;
  static constexpr u64 LOG_FETCH_WIDTH = clog2(FETCH_WIDTH);
  static constexpr u64 LINEINST        = 1 << LOGLINEINST;
  static constexpr u64 AHEAD_K         = 1 << AHEAD_BANK_BITS;
  static_assert(AHEAD_BANKS >= 1);
  static_assert(AHEAD_BANKS <= AHEAD_K);

  // Lane/bank geometry (from gshareN_ahead_best)
  static constexpr u64 LOGLANES  = std::bit_width(N - 1);
  static constexpr u64 LANES     = 1 << LOGLANES;
  static constexpr u64 LOGBANKS  = AHEAD_BANK_BITS;
  static constexpr u64 BANKS     = AHEAD_K;

  // TAGE table constants
  static constexpr u64 MAX_TAG_WIDTH  = array_max(TableCfg::TAG_WIDTH);
  static constexpr u64 MAX_CTR_WIDTH  = array_max(TableCfg::CTR_WIDTH);
  static constexpr u64 MAX_HYST_WIDTH = array_max(TableCfg::HYST_WIDTH);
  static constexpr u64 MAX_U_WIDTH    = array_max(TableCfg::U_WIDTH);
  static constexpr u64 MAX_TABLE_SIZE = array_max(TableCfg::TABLE_SIZE);
  static constexpr u64 MAX_IDX_BITS   = clog2(MAX_TABLE_SIZE);
  static constexpr u64 MINHIST        = array_min(TableCfg::HIST_LEN);
  static constexpr u64 MAXHIST        = array_max(TableCfg::HIST_LEN);
  static constexpr u64 MAX_HTAGBITS   = MAX_TAG_WIDTH - LOG_FETCH_WIDTH;

  // Bimodal
  static constexpr u64 LOG_BIMODAL_SIZE = clog2(BIMODAL_SIZE);
  static constexpr u64 BINDEX_BITS      = LOG_BIMODAL_SIZE - LOG_FETCH_WIDTH;
  static constexpr u64 BIM_ENTRIES      = BIMODAL_SIZE / FETCH_WIDTH;

  static constexpr u64 LINEADDR_BITS = std::max(BINDEX_BITS, MAX_IDX_BITS);
  static constexpr u64 UCTRBITS      = 8;
  static constexpr u64 PATHBITS       = 6;

  // ======== Table tuple ========
  // ======== Packed TAGE entry: [tag | pred | hyst | u] in one wide word ========
  // One RAM read returns all fields. Drastically reduces read count vs separate RAMs.
  static constexpr u64 ENTRY_BITS = MAX_TAG_WIDTH + MAX_CTR_WIDTH + MAX_HYST_WIDTH + MAX_U_WIDTH;
  // Field positions within packed entry (LSB to MSB): u | hyst | pred | tag
  static constexpr u64 U_OFF    = 0;
  static constexpr u64 HYST_OFF = U_OFF + MAX_U_WIDTH;
  static constexpr u64 PRED_OFF = HYST_OFF + MAX_HYST_WIDTH;
  static constexpr u64 TAG_OFF  = PRED_OFF + MAX_CTR_WIDTH;

  // ======== History ========
  geometric_folds<NUM_TABLES, MINHIST, MAXHIST, MAX_IDX_BITS, MAX_HTAGBITS>
      tage_folded_hists;
  reg<1> true_block = 1;
  reg<1> last_condbr_dir = 1;

  // ======== TAGE tables: one packed wide RAM per table per bank ========
  // Total RAMs = NUM_TABLES * AHEAD_BANKS (vs 4 * NUM_TABLES * AHEAD_BANKS before)
  hcm::ram<val<ENTRY_BITS>, MAX_TABLE_SIZE> tage_ram[NUM_TABLES][AHEAD_BANKS]{"tage"};

  // ======== Bimodal ========
  hcm::ram<val<1>, BIM_ENTRIES> bpred[FETCH_WIDTH];

  // ======== Ahead pipeline [0]=current, [1]=previous ========
  // One packed entry per table per bank
  arr<reg<ENTRY_BITS>, NUM_TABLES> ahead_entry[2][AHEAD_BANKS];

  // Shared across banks
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> tage_indexes[2];
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> ahead_htag[2];

  // (ahead_lane_pred removed — P1 uses bimodal, P2 extracts from packed entries)

  // Bank selection state
  reg<LOGBANKS> XB;
  reg<LOGBANKS> path;

  // ======== Block state ========
  // Precomputed prediction for current block
  arr<reg<1>, FETCH_WIDTH> pred;
  reg<FETCH_WIDTH> p1;  // concat of pred, for block_entry masking

  // Block geometry (from gshareN_ahead_best)
  reg<LOGLINEINST> block_entry;         // offset of entry point in line
  reg<N + 1> rank;                       // one-hot: which branch we're predicting next

  // Per-branch directions (for update_cycle)
  arr<reg<1>, N> branch_dir;

  // P2 working registers (written in predict1, read in update_cycle)
  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES>  readt;
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES>  readc;
  arr<reg<MAX_HYST_WIDTH>, NUM_TABLES> readh;
  arr<reg<MAX_U_WIDTH>, NUM_TABLES>    readu;
  reg<NUM_TABLES> not_u_mask;
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match1;
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match2;
  arr<reg<1>, FETCH_WIDTH> tage_pred1;
  arr<reg<1>, FETCH_WIDTH> tage_pred2;
  reg<FETCH_WIDTH> p2;

  // Meta
  arr<reg<METABITS, i64>, METAPIPE> meta;
  arr<reg<1>, FETCH_WIDTH> newly_alloc;

  // U-bit management
  reg<UCTRBITS> epoch_ctr;
  std::conditional_t<(DECAY_CTR > 0), reg<DECAY_CTR == 0 ? 1 : DECAY_CTR>,
                     EmptyMember>
      decay_threshold;

  // Bimodal: pipelined like TAGE. [0]=current read, [1]=previous (for predict1 return)
  reg<BINDEX_BITS> bindex;
  arr<reg<1>, FETCH_WIDTH> ahead_readb[2];

  // UPDATE_ONLY zone
  zone UPDATE_ONLY;
  hcm::ram<val<1>, BIM_ENTRIES> bhyst[FETCH_WIDTH];

  // Plain C++ state
  u64 num_branch = 0;
  u64 block_size = 0;

  // ======== Helpers ========

  val<1> line_end() {
    return (block_entry + block_size) == hard<LINEINST>{};
  }

  val<1> last_pred() {
    assert(num_branch <= N);
    return rank >> (N - num_branch);
  }

  // ======== Predictor Interface ========

  val<1> predict1(val<64> inst_pc) override {
    inst_pc.fanout(hard<4>{});
    true_block.fanout(hard<8 + AHEAD_BANKS * 2>{});

    // ---- 1. Block setup (conditional on true_block) ----
    // If previous block was not "true" (speculative replay), continue from
    // where we left off. Don't re-read RAMs. (gshareN_ahead_best technique)
    block_entry = select(true_block, val<LOGLINEINST>{inst_pc >> 2},
                         val<LOGLINEINST>{block_entry + block_size});
    block_entry.fanout(hard<N * 2 + 4>{}); // line_end calls + predict returns

    rank = select(true_block, val<N + 1>{1}, rank << num_branch);
    rank.fanout(hard<N + 2>{});

    // ---- 2. Pipeline shift + TAGE reads (conditional on true_block) ----
    // Only re-read RAMs when this is a true block. Otherwise reuse previous.
    execute_if(true_block, [&]() {
      // Pipeline shift: [0] → [1] (packed entries + bimodal + shared state)
      for (u64 b = 0; b < AHEAD_BANKS; b++) {
        for (u64 i = 0; i < NUM_TABLES; i++) {
          ahead_entry[1][b][i] = ahead_entry[0][b][i];
        }
      }
      for (u64 i = 0; i < NUM_TABLES; i++) {
        tage_indexes[1][i] = tage_indexes[0][i];
        ahead_htag[1][i]   = ahead_htag[0][i];
      }
      for (u64 i = 0; i < FETCH_WIDTH; i++) {
        ahead_readb[1][i] = ahead_readb[0][i];
      }

      // ---- 3. Compute TAGE indices + tags ----
      val<LINEADDR_BITS> lineaddr = inst_pc >> (LOG_FETCH_WIDTH + 2);
      lineaddr.fanout(hard<1 + NUM_TABLES * 2>{});
      tage_folded_hists.fanout(hard<2>{});

      for (u64 i = 0; i < NUM_TABLES; i++) {
        tage_indexes[0][i] = lineaddr ^ tage_folded_hists.template get<0>(i);
      }
      tage_indexes[0].fanout(hard<AHEAD_BANKS>{});  // one read per bank

      for (u64 i = 0; i < NUM_TABLES; i++) {
        ahead_htag[0][i] = val<MAX_HTAGBITS>{lineaddr}.reverse()
                           ^ tage_folded_hists.template get<1>(i);
      }

      // Bimodal
      bindex = lineaddr;
      bindex.fanout(hard<FETCH_WIDTH>{});

      // ---- 4. Read packed TAGE RAMs + bimodal ----
      // One read per table per bank = NUM_TABLES * AHEAD_BANKS total reads
      // (vs 4 * NUM_TABLES * AHEAD_BANKS with separate RAMs)
      for (u64 b = 0; b < AHEAD_BANKS; b++) {
        for (u64 i = 0; i < NUM_TABLES; i++) {
          ahead_entry[0][b][i] = tage_ram[i][b].read(val<MAX_IDX_BITS>{tage_indexes[0][i]});
        }
      }
      for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
        ahead_readb[0][offset] = bpred[offset].read(bindex);
      }

      // ---- 5. Path update (for bank select) ----
      path = XB + num_branch + ~last_condbr_dir;
    }); // end execute_if(true_block)

    // P1 prediction = bimodal from PREVIOUS cycle's read (ahead pipeline [1]).
    // This is truly ahead — no RAM read in the return path. Just reg reads.
    // predict2 will do full TAGE tag compare + match to override.
    ahead_readb[1].fanout(hard<FETCH_WIDTH * 2>{});
    for (u64 i = 0; i < FETCH_WIDTH; i++) {
      pred[i] = ahead_readb[1][i];
    }
    pred.fanout(hard<N * 4 + 2>{}); // predict1(1) + reuse_p1(N) + predict2(1) + reuse_p2(N) + update(1)

    // Update XB for next block
    XB = select(true_block, val<LOGBANKS>{inst_pc >> 6},
                val<LOGBANKS>{XB + num_branch});

    // Reset block state
    block_size = 1;
    num_branch = 0;
    reuse_prediction(~line_end());
    return pred[0];
  }

  val<1> reuse_predict1(val<64> inst_pc) override {
    block_size++;
    reuse_prediction(~(line_end() | last_pred()));
    return pred[num_branch];
  }

  val<1> predict2(val<64> inst_pc) override {
    // P2: full TAGE tag compare + match on ahead_*[1][bank].
    // Bank select using path (computed in previous update_cycle).
    path.fanout(hard<NUM_TABLES * 4 + 2>{});

    // Clamp bank_sel to valid range
    val<1> bank_hit = (path < hard<AHEAD_BANKS>{});
    val<LOGBANKS> bank_sel = select(bank_hit, path, val<LOGBANKS>{0});
    bank_sel.fanout(hard<NUM_TABLES * 4 + 1>{});
    bank_hit.fanout(hard<FETCH_WIDTH>{});

    // Select correct bank's packed entries, extract fields into working regs
    for (u64 b = 0; b < AHEAD_BANKS; b++) {
      ahead_entry[1][b].fanout(hard<2>{});
    }
    static_loop<NUM_TABLES>([&]<u64 I>() {
      // Mux across banks for this table
      arr<val<ENTRY_BITS>, AHEAD_BANKS> e_opts = [&](u64 b) -> val<ENTRY_BITS> {
        return ahead_entry[1][b][I];
      };
      val<ENTRY_BITS> entry = e_opts.fo1().select(bank_sel);
      entry.fanout(hard<4>{}); // extract tag, pred, hyst, u
      // Extract fields from packed entry
      readt[I] = entry >> TAG_OFF;
      readc[I] = entry >> PRED_OFF;
      readh[I] = entry >> HYST_OFF;
      readu[I] = entry >> U_OFF;
    });
    readt.fanout(hard<FETCH_WIDTH + 1>{});
    readc.fanout(hard<3>{});
    readu.fanout(hard<2>{});
    readh.fanout(hard<2>{});

    // Tag compare
    ahead_htag[1].fanout(hard<3>{});
    arr<val<1>, NUM_TABLES> htagcmp_split = [&](int i) {
      return val<MAX_HTAGBITS>{readt[i]} == ahead_htag[1][i];
    };
    val<NUM_TABLES> htagcmp = htagcmp_split.fo1().concat();
    htagcmp.fanout(hard<FETCH_WIDTH>{});

    // Prediction bits
    val<NUM_TABLES> gpreds = [&]() -> val<NUM_TABLES> {
      if constexpr (MAX_CTR_WIDTH == 1) {
        return readc.concat();
      } else {
        arr<val<1>, NUM_TABLES> gp = [&](int i) -> val<1> {
          return readc[i] >> hard<MAX_CTR_WIDTH - 1>{};
        };
        return gp.fo1().concat();
      }
    }();
    gpreds.fanout(hard<FETCH_WIDTH>{});

    // Per-offset: preds vector, match, one_hot, prediction
    arr<val<NUM_TABLES + 1>, FETCH_WIDTH> preds = [&](u64 offset) {
      return concat(ahead_readb[1][offset], gpreds);
    };
    preds.fanout(hard<2 * FETCH_WIDTH>{});

    static_loop<FETCH_WIDTH>([&]<u64 offset>() {
      arr<val<1>, NUM_TABLES> tagcmp = [&](int i) {
        return val<LOG_FETCH_WIDTH>{readt[i] >> MAX_HTAGBITS} == hard<offset>{};
      };
      val<NUM_TABLES + 1> m = concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp);
      m.fanout(hard<2>{});
      match1[offset] = m.one_hot();
      match2[offset] = (m ^ m.one_hot()).one_hot();
    });
    match1.fanout(hard<4>{});
    match2.fanout(hard<3>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      tage_pred1[offset] = (match1[offset] & preds[offset]) != hard<0>{};
      tage_pred2[offset] = (match2[offset] & preds[offset]) != hard<0>{};
    }
    tage_pred1.fanout(hard<3>{});
    tage_pred2.fanout(hard<3>{});

    not_u_mask = ~readu.concat();
    not_u_mask.fanout(hard<2>{});

    // Meta select
    if constexpr (bool(USE_META)) {
      meta.fanout(hard<4>{});
      arr<val<1>, NUM_TABLES> weakctr = [&](int i) -> val<1> {
        return readh[i] == hard<0>{};
      };
      val<NUM_TABLES> coldctr = not_u_mask & weakctr.fo1().concat();
      coldctr.fanout(hard<FETCH_WIDTH>{});
      val<1> metasign = (meta[METAPIPE - 1] >= hard<0>{});
      metasign.fanout(hard<FETCH_WIDTH>{});
      for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
        newly_alloc[offset] = (match1[offset] & coldctr) != hard<0>{};
      }
      newly_alloc.fanout(hard<2>{});

      arr<val<1>, FETCH_WIDTH> altsel = [&](u64 offset) {
        arr<val<1>, 3> inputs = {metasign, newly_alloc[offset],
                                 match2[offset] != hard<0>{}};
        return inputs.fo1().fold_and();
      };

      // Final P2: select(bank_hit, tage_with_meta, bimodal)
      p2 = arr<val<1>, FETCH_WIDTH>{[&](u64 offset) {
             val<1> ahead = select(altsel[offset].fo1(), tage_pred2[offset],
                                   tage_pred1[offset]);
             return select(bank_hit, ahead, ahead_readb[1][offset]);
           }}.concat();
    } else {
      p2 = arr<val<1>, FETCH_WIDTH>{[&](u64 offset) {
             return select(bank_hit, tage_pred1[offset], ahead_readb[1][offset]);
           }}.concat();
    }

    // TODO: P2 should return TAGE-corrected prediction per branch.
    // For now, return same as P1 (bimodal). P2 correction needs lane scrambling.
    reuse_prediction(~(line_end() | last_pred()));
    return pred[num_branch];
  }

  val<1> reuse_predict2(val<64> inst_pc) override {
    return pred[num_branch];
  }

  void update_condbr(val<64> branch_pc, val<1> taken,
                     val<64> next_pc) override {
    assert(num_branch < N);
    branch_dir[num_branch] = taken.fo1();
    num_branch++;
    // End block if we've used all N prediction slots or hit line boundary
    reuse_prediction(~(line_end() | last_pred()));
  }

  void update_cycle(instruction_info &block_end_info) override {
    // TODO:
    // 1. No-branch early return (history update only)
    // 2. Allocation candidates
    // 3. Counter/hyst/u-bit updates (banked writes)
    // 4. Bimodal updates
    // 5. Meta update
    // 6. U-bit epoch/decay
    // 7. Per-branch history update (shift by num_branch)
    // 8. Path update: path = XB + num_branch + ~last_condbr_dir
    // 9. true_block update
    (void)block_end_info;
  }
};
