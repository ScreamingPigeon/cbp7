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

  // ======== Ahead pipeline: single stage ========
  // Written by RAM reads in predict1. Read in the NEXT predict1 (after next_cycle
  // resets timing) for precomputation. Then overwritten with new RAM reads.
  // No [0]/[1] shift needed — next_cycle between calls handles the pipeline.
  arr<reg<ENTRY_BITS>, NUM_TABLES> ahead_entry[AHEAD_BANKS];

  // Shared across banks
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> tage_indexes;
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> ahead_htag;

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

  // Per-branch state (for update_cycle)
  arr<reg<1>, N> branch_dir;
  arr<reg<LOG_FETCH_WIDTH>, N> br_offset;  // lane of each branch

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

  // Bimodal: single stage (same pattern as ahead_entry)
  reg<BINDEX_BITS> bindex;
  arr<reg<1>, FETCH_WIDTH> ahead_readb;

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

    // ---- 2. TAGE precomputation from ahead_entry[] (PREVIOUS cycle's data) ----
    // These regs were written in the previous predict1 call. next_cycle() was
    // called between, so their timing starts fresh. No pipeline shift needed.

    // ---- Fast prediction path using fo1() — zero fanout cost on critical path ----
    // Intermediate vals consumed once via fo1(). Regs written separately for update_cycle.

    // Bank select (path consumed once for prediction, once for update_cycle bank)
    path.fanout(hard<2>{});
    val<1> bank_hit = (path < hard<AHEAD_BANKS>{});
    val<LOGBANKS> bank_sel = select(bank_hit, path, val<LOGBANKS>{0});
    bank_sel.fanout(hard<NUM_TABLES>{}); // one select per table

    // Extract entries from selected bank
    // Also write to regs for update_cycle (reg writes don't affect val timing)
    arr<val<MAX_TAG_WIDTH>, NUM_TABLES> v_tag = [&](u64 i) -> val<MAX_TAG_WIDTH> {
      arr<val<ENTRY_BITS>, AHEAD_BANKS> e_opts = [&](u64 b) -> val<ENTRY_BITS> {
        return ahead_entry[b][i];
      };
      val<ENTRY_BITS> entry = e_opts.fo1().select(bank_sel);
      // Write to regs for update_cycle (doesn't affect entry timing)
      readt[i] = entry >> TAG_OFF;
      readc[i] = entry >> PRED_OFF;
      readh[i] = entry >> HYST_OFF;
      readu[i] = entry >> U_OFF;
      return val<MAX_TAG_WIDTH>{entry >> TAG_OFF};
    };

    // Tag compare — v_tag read for htag compare + offset extract
    v_tag.fanout(hard<FETCH_WIDTH + 1>{}); // htag compare(1) + offset extract per offset(FW)
    arr<val<1>, NUM_TABLES> htagcmp_arr = [&](int i) -> val<1> {
      return val<MAX_HTAGBITS>{v_tag[i]} == ahead_htag[i];
    };
    val<NUM_TABLES> htagcmp = htagcmp_arr.fo1().concat();
    htagcmp.fanout(hard<FETCH_WIDTH>{}); // needed: one per offset

    // Prediction bits — fo1 on readc
    val<NUM_TABLES> gpreds = [&]() -> val<NUM_TABLES> {
      if constexpr (MAX_CTR_WIDTH == 1) {
        arr<val<1>, NUM_TABLES> gp = [&](int i) -> val<1> { return readc[i]; };
        return gp.fo1().concat();
      } else {
        arr<val<1>, NUM_TABLES> gp = [&](int i) -> val<1> {
          return readc[i] >> hard<MAX_CTR_WIDTH - 1>{};
        };
        return gp.fo1().concat();
      }
    }();
    gpreds.fanout(hard<FETCH_WIDTH>{}); // needed: one per offset

    // Per-offset: match + one_hot + prediction → write to pred[]
    // bank_hit needs fanout for FETCH_WIDTH final selects
    bank_hit.fanout(hard<FETCH_WIDTH>{});
    ahead_readb.fanout(hard<FETCH_WIDTH * 2>{}); // preds concat + final select

    static_loop<FETCH_WIDTH>([&]<u64 offset>() {
      arr<val<1>, NUM_TABLES> tagcmp = [&](int i) -> val<1> {
        return val<LOG_FETCH_WIDTH>{v_tag[i] >> MAX_HTAGBITS} == hard<offset>{};
      };
      val<NUM_TABLES + 1> m = concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp);
      m.fanout(hard<2>{}); // one_hot + match2 XOR
      val<NUM_TABLES + 1> m1 = m.one_hot();
      m1.fanout(hard<2>{}); // preds AND + match1 reg
      val<NUM_TABLES + 1> preds_v = concat(ahead_readb[offset], gpreds);
      val<1> tpred = (m1 & preds_v.fo1()) != hard<0>{};

      // Write match regs for update_cycle
      match1[offset] = m1;
      match2[offset] = (m ^ m1).one_hot();
      tage_pred1[offset] = tpred;

      // Final prediction: TAGE when bank hits, bimodal fallback
      pred[offset] = select(bank_hit, tpred, ahead_readb[offset]);
    });

    // Fanout for regs read in update_cycle (NOT on predict1 critical path)
    readt.fanout(hard<FETCH_WIDTH + 1>{});
    readc.fanout(hard<3>{});
    readu.fanout(hard<2>{});
    readh.fanout(hard<2>{});
    match1.fanout(hard<4>{});
    match2.fanout(hard<3>{});
    tage_pred1.fanout(hard<3>{});
    not_u_mask = ~readu.concat();
    not_u_mask.fanout(hard<2>{});

    pred.fanout(hard<N * 4 + 2>{});

    // ---- 3. RAM reads for NEXT cycle (overwrite ahead_entry/readb/htag/indexes) ----
    execute_if(true_block, [&]() {
      val<LINEADDR_BITS> lineaddr = inst_pc >> (LOG_FETCH_WIDTH + 2);
      lineaddr.fanout(hard<1 + NUM_TABLES * 2>{});
      tage_folded_hists.fanout(hard<2>{});

      for (u64 i = 0; i < NUM_TABLES; i++)
        tage_indexes[i] = lineaddr ^ tage_folded_hists.template get<0>(i);
      tage_indexes.fanout(hard<AHEAD_BANKS>{});

      for (u64 i = 0; i < NUM_TABLES; i++)
        ahead_htag[i] = val<MAX_HTAGBITS>{lineaddr}.reverse()
                        ^ tage_folded_hists.template get<1>(i);

      bindex = lineaddr;
      bindex.fanout(hard<FETCH_WIDTH>{});

      for (u64 b = 0; b < AHEAD_BANKS; b++)
        for (u64 i = 0; i < NUM_TABLES; i++)
          ahead_entry[b][i] = tage_ram[i][b].read(val<MAX_IDX_BITS>{tage_indexes[i]});

      for (u64 offset = 0; offset < FETCH_WIDTH; offset++)
        ahead_readb[offset] = bpred[offset].read(bindex);

      path = XB + num_branch + ~last_condbr_dir;
    }); // end execute_if(true_block)

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

  // P2 = reg read. All TAGE logic precomputed in predict1. P2 ceil = 1.
  val<1> predict2(val<64> inst_pc) override {
    return pred[num_branch];
  }

  val<1> reuse_predict2(val<64> inst_pc) override {
    return pred[num_branch];
  }

  void update_condbr(val<64> branch_pc, val<1> taken,
                     val<64> next_pc) override {
    assert(num_branch < N);
    br_offset[num_branch] = branch_pc.fo1() >> 2;
    branch_dir[num_branch] = taken.fo1();
    num_branch++;
    // End block if we've used all N prediction slots or hit line boundary
    reuse_prediction(~(line_end() | last_pred()));
  }

  void update_cycle(instruction_info &block_end_info) override {
    val<1> &mispredict = block_end_info.is_mispredict;
    val<64> &next_pc = block_end_info.next_pc;

    // ---- 1. No-branch early return ----
    if (num_branch == 0) {
      // Per-branch history: shift by 0 branches = just inject next_pc
      tage_folded_hists.update(val<PATHBITS>{next_pc.fo1() >> 2});
      last_condbr_dir = 0;
      true_block = 1;
      return;
    }

    // ---- Fanout declarations ----
    mispredict.fanout(hard<NUM_TABLES + 2>{});
    val<1> correct_pred = ~mispredict;
    correct_pred.fanout(hard<NUM_TABLES + 2>{});
    tage_indexes.fanout(hard<4>{});
    ahead_htag.fanout(hard<3>{});
    readt.fanout(hard<4>{});
    readc.fanout(hard<2>{});
    readh.fanout(hard<3>{});
    readu.fanout(hard<2>{});
    match1.fanout(hard<3>{});
    match2.fanout(hard<2>{});
    tage_pred1.fanout(hard<2>{});
    tage_pred2.fanout(hard<2 + NUM_TABLES>{});
    branch_dir.fanout(hard<3>{});

    last_condbr_dir = branch_dir[num_branch - 1];
    last_condbr_dir.fanout(hard<N + 2>{});

    // ---- Branch direction association (BR_P_ENTRY=1) ----
    // Each table's entry has an offset in its tag identifying which branch it predicts.
    // last_offset: lane of the last conditional branch (for allocation)
    // We need branch_offset info — but update_condbr only saves direction, not offset.
    // For now, use a simplified approach: all tables predict the last branch.
    // TODO: save branch offsets in update_condbr for proper per-table association.
    val<1> last_dir = branch_dir[num_branch - 1];
    last_dir.fanout(hard<NUM_TABLES * 3>{});

    // ---- Prediction correctness per table ----
    arr<val<1>, NUM_TABLES> badpred1 = [&](u64 i) -> val<1> {
      if constexpr (MAX_CTR_WIDTH == 1) {
        return readc[i] != last_dir;
      } else {
        return val<1>{readc[i] >> hard<MAX_CTR_WIDTH - 1>{}} != last_dir;
      }
    };
    badpred1.fanout(hard<3>{});

    // Provider: which tables have a tag match (from predict2's match1)
    // match1[0] is NUM_TABLES+1 bits (includes bimodal). Extract table bits only.
    val<NUM_TABLES> primary_bits = match1[0] >> 1; // skip bimodal bit
    primary_bits.fanout(hard<2>{});
    arr<val<1>, NUM_TABLES> primary = primary_bits.make_array(val<1>{});
    primary.fanout(hard<3>{});

    // Weak hysteresis: primary provider with wrong prediction and weak counter
    arr<val<1>, NUM_TABLES> g_weak = [&](u64 i) -> val<1> {
      return primary[i] & badpred1[i] & (readh[i] == hard<0>{});
    };
    g_weak.fanout(hard<2>{});

    // ---- Allocation (only on misprediction, only post-provider tables) ----
    not_u_mask.fanout(hard<2>{});
    val<NUM_TABLES> mispmask = mispredict.replicate(hard<NUM_TABLES>{}).concat();
    // Find provider table, allocate only from tables with longer history
    val<NUM_TABLES + 1> last_match1 = match1[0]; // provider one-hot from predict2
    last_match1.fanout(hard<2>{});
    val<NUM_TABLES> postmask = mispmask.fo1() & val<NUM_TABLES>(last_match1 - 1);
    postmask.fanout(hard<2>{});
    val<NUM_TABLES> candallocmask = postmask & not_u_mask;
    candallocmask.fanout(hard<2>{});
    val<NUM_TABLES> collamask = candallocmask.reverse();
    collamask.fanout(hard<2>{});
    val<NUM_TABLES> collamask1 = collamask.one_hot();
    collamask1.fanout(hard<3>{});
    val<NUM_TABLES> collamask2 = (collamask ^ collamask1).one_hot();
    val<NUM_TABLES> collamask12 =
        select(val<2>{std::rand()} == hard<0>{}, collamask2.fo1(), collamask1);
    arr<val<1>, NUM_TABLES> allocate =
        collamask12.fo1().reverse().make_array(val<1>{});
    allocate.fanout(hard<4>{});

    // ---- need_extra_cycle ----
    // Only when provider is wrong with weak hyst, or misprediction
    val<1> some_g_weak = (g_weak.concat() != hard<0>{});
    val<1> extra_cycle = some_g_weak.fo1() | mispredict;
    extra_cycle.fanout(hard<NUM_TABLES * 2 + 1>{});
    need_extra_cycle(extra_cycle);

    // ---- Bank for writes ----
    val<LOGBANKS> raw_wbank = path;
    val<1> write_bank_valid = (raw_wbank < hard<AHEAD_BANKS>{});
    val<LOGBANKS> write_bank = select(write_bank_valid, raw_wbank, val<LOGBANKS>{0});
    write_bank.fanout(hard<NUM_TABLES * 2>{});

    // ---- Packed TAGE writes ----
    for (u64 i = 0; i < NUM_TABLES; i++) {
      // Counter update: toward actual direction
      val<MAX_CTR_WIDTH> new_pred = select(allocate[i],
          val<MAX_CTR_WIDTH>{last_dir},
          update_ctr(readc[i], last_dir));

      // Hysteresis: reset on allocate, update on primary
      val<MAX_HYST_WIDTH> new_hyst = select(allocate[i],
          val<MAX_HYST_WIDTH>{0},
          update_ctr(readh[i], ~badpred1[i]));

      // U-bit: clear on allocate, keep on update
      val<MAX_U_WIDTH> new_u = select(allocate[i],
          val<MAX_U_WIDTH>{0}, readu[i]);

      // Tag: new on allocate (with offset 0 for now), keep on update
      val<MAX_TAG_WIDTH> new_tag = select(allocate[i],
          concat(val<LOG_FETCH_WIDTH>{0}, ahead_htag[i]),
          readt[i]);

      // Pack: [tag | pred | hyst | u]
      val<ENTRY_BITS> new_entry = concat(
          concat(new_tag, new_pred),
          concat(new_hyst, new_u));
      new_entry.fanout(hard<AHEAD_BANKS>{});

      // Write gated by: (provider with wrong weak pred) OR (allocation)
      execute_if(extra_cycle & (g_weak[i] | allocate[i]), [&]() {
        static_loop<AHEAD_BANKS>([&]<u64 b>() {
          execute_if(write_bank == hard<b>{}, [&]() {
            tage_ram[i][b].write(val<MAX_IDX_BITS>{tage_indexes[i]},
                                 new_entry);
          });
        });
      });
    }

    // ---- Meta counter update (simplified) ----
    if constexpr (bool(USE_META)) {
      for (u64 i = METAPIPE - 1; i != 0; i--)
        meta[i] = meta[i - 1];
      // TODO: proper meta update based on primary vs alt correctness
    }

    // ---- Per-branch history update (shift by num_branch, not 1) ----
    // This is the key improvement from gshareN_ahead_best.
    // More accurate history = better TAGE correlation.
    tage_folded_hists.update(val<PATHBITS>{next_pc.fo1() >> 2});
    // TODO: shift by num_branch instead of 1 (requires custom fold update)

    // ---- true_block update ----
    true_block =
        arr<val<1>, 4>{~mispredict, last_condbr_dir, last_pred(), line_end()}
            .fold_or();
  }
};
