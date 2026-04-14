#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"
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
          u64 LINEINST = 8,            // line size in instructions
          bool SHARED_HYS = true,      // shared hyst: 2 entries share 1 counter
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
                                           U_WIDTH, SEC_TAG_BITS, N, SHARED_HYS,
                                           std::make_index_sequence<NT>>::type;
  Tables tables;
  hcm::ram<val<N>, BIMODAL_CAPACITY> bim_ctr{"bim"};
  hcm::ram<val<META_WIDTH>, META_CAPACITY> meta_ctr{"meta"};

  // ======================================================================
  // Registers
  // ======================================================================
  // ---- Global History (shared, folds live in per-table TATable) ----
  ta_global_history<MAXHIST> gh;

  // ---- Bimodal Fallback (ahead-pipelined) ----
  reg<PRED_BITS> prefetch_bim;
  reg<PRED_BITS> current_bim;

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

// ---- Timing debug taps (zero cost in normal builds) ----
#ifdef TIMING_DEBUG
  reg<1> dbg_full_hits;     // after per-table hit computation
  reg<1> dbg_bim_pred;      // after bimodal read
  reg<1> dbg_match;         // after concat into match bitmask
  reg<1> dbg_match1;        // after one_hot (provider)
  reg<1> dbg_match2;        // after one_hot (alt)
  reg<1> dbg_provider_pred; // after replicate-mask-fold
  reg<1> dbg_alt_pred;      // after replicate-mask-fold
  reg<1> dbg_provider_weak; // after weakness check
  reg<1> dbg_has_alt;       // after has_alt mask
  reg<1> dbg_meta_use_alt;  // after meta pipeline read
  reg<1> dbg_use_alt;       // after final use_alt AND
  reg<1> dbg_final_pred;    // after select mux
  reg<1> dbg_mispredict;    // framework's is_mispredict
  reg<1> dbg_true_block;    // true_block after computation
  reg<1> dbg_inst_pc;       // inst_pc arrival in predict1
  reg<1> dbg_fold_idx;      // worst fold_idx.get() across tables
  reg<1> dbg_fold_tag;      // worst fold_tag.get() across tables
  reg<1> dbg_p1_return;     // predict1 return value timing
  reg<1> dbg_hist_input;    // hist_input timing in update_cycle
  reg<1> dbg_gh_fanout;     // gh after fanout in update_cycle
  reg<1> dbg_fold_compute;  // compute_update result timing
  reg<1> dbg_fold_apply;    // after apply_update (fold write)
#endif

  // ======================================================================
  // Helpers
  // ======================================================================

  val<1> line_end() { return (block_entry + block_size) == hard<LINEINST>{}; }

  // ======================================================================
  // predict1/2, reuse_predict1/2, update_condbr, update_cycle — TODO
  // ======================================================================

  val<1> predict1([[maybe_unused]] val<64> inst_pc) {
    inst_pc.fanout(
        hard<2 * NT + 1>{}); // 2 reads per table (>>2, >>4) + bimodal (>>2)

#ifdef TIMING_DEBUG
    dbg_inst_pc = val<1>{inst_pc};
#endif

    // Ahead reads for next block (off crit path, needs inst_pc)
    execute_if(true_block, [&]() {
      static_loop<NT>([&]<u64 I>() {
        auto &t = std::get<I>(tables);
        auto idx = t.fold_idx.get() ^ val<t.IDX_BITS>{inst_pc >> 2};
        idx.fanout(hard<6>{}); // 5 RAM reads + prefetch_idx write
        auto computed_tag = t.fold_tag.get() ^ val<t.tag_width>{inst_pc >> 4};

#ifdef TIMING_DEBUG
        if constexpr (I == 0) {
          dbg_fold_idx = val<1>{t.fold_idx.get()};
          dbg_fold_tag = val<1>{t.fold_tag.get()};
        }
#endif

        auto stored_tag = t.tag_ram.read(idx);
        stored_tag.fanout(hard<2>{});
        prefetch_tag[I] = stored_tag;
        prefetch_tag_hit[I] =
            val<MAX_TAG_WIDTH>{stored_tag} == val<MAX_TAG_WIDTH>{computed_tag};
        prefetch_pred[I] = t.pred_ram.read(idx);
        prefetch_sec[I] = t.sec_ram.read(idx);
        prefetch_idx[I] = idx;
        prefetch_hyst[I] = t.hyst_ram.read(val<t.HYST_IDX_BITS>{idx});
        prefetch_u[I] = t.u_ram.read(idx);
      });

      // Bimodal ahead read (direct-mapped, no tag match needed)
      auto bim_idx = val<ta::clog2(BIMODAL_CAPACITY)>{inst_pc >> 2};
      prefetch_bim = bim_ctr.read(bim_idx);
    });

    // Crit path: just read precomputed prediction from reg
    block_entry.fanout(hard<2 * LINEINST>{}); // read in line_end() across
                                              // predict + reuse + update_condbr
    pred.fanout(
        hard<2 * LINEINST + 1>{}); // predict1/2, reuse + 1 training old_pred read
    block_size = 1;
    num_branch = 0;
    reuse_prediction((block_entry + block_size) == hard<LINEINST>{});
#ifdef TIMING_DEBUG
    dbg_p1_return = val<1>{pred[num_branch]};
#endif
    return pred[num_branch];
  }
  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) {
    block_size++;
    reuse_prediction((block_entry + block_size) == hard<LINEINST>{});
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
    current_bim = prefetch_bim;

    // Precompute secondary tag for next block
    block_end_info.next_pc.fanout(hard<3>{}); // curr_sec_tag + meta_idx + hist path_bits
    curr_sec_tag = val<SEC_TAG_BITS>{block_end_info.next_pc >> 2};

    // ================================================================
    // Provider / altpred resolution via bitmask + one_hot
    //
    // We cannot reassign val<N> (private operator=), so we avoid
    // accumulation loops. Instead:
    //   1. Compute per-table hit bits → arr<val<1>, NT>
    //   2. Concat into val<NT+1> with bimodal as MSB always-hit
    //   3. one_hot() → lowest set bit = longest-history hit = provider
    //   4. one_hot() on remainder → alt
    //   5. Replicate one-hot bits to PRED_BITS width, AND with each
    //      table's prediction, fold_or → extract provider/alt pred
    // ================================================================

    // 2. Per-table hit: primary tag matched (off crit path) AND
    //    stored secondary tag matches curr_sec_tag (same successor path)
    curr_sec_tag.fanout(hard<NT>{}); // compared once per table
    arr<val<1>, NT> full_hits = [&](u64 i) {
      return val<1>{current_tag_hit[i]} & (val<SEC_TAG_BITS>{current_sec[i]} ==
                                           val<SEC_TAG_BITS>{curr_sec_tag});
    };

#ifdef TIMING_DEBUG
    dbg_full_hits = full_hits[0];
#endif

    // 3. Bimodal — ahead-pipelined, already in current_bim from pipe shift.
    val<PRED_BITS> bim_pred = val<PRED_BITS>{current_bim};

#ifdef TIMING_DEBUG
    dbg_bim_pred = val<1>{bim_pred};
#endif

    // 4. Build match bitmask.
    //    Bit layout of val<NT+1>:
    //      bit 0     = table 0 (longest history)
    //      bit NT-1  = table NT-1 (shortest history)
    //      bit NT    = bimodal (always 1 — fallback)
    //    one_hot() returns lowest set bit → longest-history hit = provider.
    //    Second one_hot() on (match ^ match1) → alt provider.
    val<NT + 1> match = concat(val<1>{1}, full_hits.concat());
    match.fanout(hard<2>{}); // one_hot + XOR
    val<NT + 1> match1 = match.one_hot();
    match1.fanout(hard<2>{}); // make_array + XOR with match
    val<NT + 1> match2 = (match ^ match1).one_hot();
    match2.fanout(hard<2>{}); // make_array + has_alt mask

#ifdef TIMING_DEBUG
    dbg_match = val<1>{match};
    dbg_match1 = val<1>{match1};
    dbg_match2 = val<1>{match2};
#endif

    // 5. Prediction array: one PRED_BITS-wide entry per table + bimodal.
    //    table_preds[0..NT-1] = TAGE tables, table_preds[NT] = bimodal.
    arr<val<PRED_BITS>, NT + 1> table_preds = [&](u64 i) -> val<PRED_BITS> {
      if (i < NT)
        return val<PRED_BITS>{current_pred[i]};
      return val<PRED_BITS>{bim_pred};
    };

    // 6. Extract provider and alt predictions.
    //    For each table: replicate its one-hot match bit to PRED_BITS width,
    //    AND with that table's prediction. Since match1 is one-hot, exactly
    //    one table contributes non-zero bits. fold_or collapses to that pred.
    arr<val<1>, NT + 1> m1_bits = match1.make_array(val<1>{});
    arr<val<1>, NT + 1> m2_bits = match2.make_array(val<1>{});
    m1_bits.fanout(hard<2>{});
    m2_bits.fanout(hard<2>{});

    arr<val<PRED_BITS>, NT + 1> provider_masked = [&](u64 i) {
      return m1_bits[i].replicate(hard<PRED_BITS>{}).concat() & table_preds[i];
    };
    arr<val<PRED_BITS>, NT + 1> alt_masked = [&](u64 i) {
      return m2_bits[i].replicate(hard<PRED_BITS>{}).concat() & table_preds[i];
    };
    val<PRED_BITS> provider_pred = provider_masked.fold_or();
    val<PRED_BITS> alt_pred = alt_masked.fold_or();

#ifdef TIMING_DEBUG
    dbg_provider_pred = val<1>{provider_pred};
    dbg_alt_pred = val<1>{alt_pred};
#endif

    // 7. Provider weakness: newly allocated entry = hyst==0 AND u==0.
    //    Only check the provider table (mask by match1 bit). Bimodal
    //    (index NT) is never considered "weak".
    arr<val<1>, NT + 1> weak_arr = [&](u64 i) -> val<1> {
      if (i < NT)
        return m1_bits[i] &
               val<1>{val<std::max(u64(1), HYST_WIDTH)>{current_hyst[i]} ==
                      hard<0>{}} &
               val<1>{val<U_WIDTH>{current_u[i]} == hard<0>{}};
      return val<1>{0};
    };
    val<1> provider_weak = weak_arr.fold_or();

    // 8. has_alt: does the alt match point to a TAGE table (not bimodal)?
    //    Mask out the bimodal bit (MSB) and check if anything remains.
    val<1> has_alt = (match2 & val<NT + 1>{(u64(1) << NT) - 1}) != hard<0>{};

#ifdef TIMING_DEBUG
    dbg_provider_weak = provider_weak;
    dbg_has_alt = has_alt;
#endif

    // 9. Meta counter: predicts whether to trust a newly-allocated provider
    //    or fall back to alt. Read from PC-indexed RAM, shifted through a
    //    META_PIPE-stage pipeline. Sign bit of oldest stage decides.
    auto meta_idx = val<ta::clog2(META_CAPACITY)>{block_end_info.next_pc >> 2};
    meta_pipe[0] = meta_ctr.read(meta_idx);
    for (u64 i = META_PIPE - 1; i > 0; i--)
      meta_pipe[i] = meta_pipe[i - 1];
    val<1> meta_use_alt =
        val<META_WIDTH, i64>{meta_pipe[META_PIPE - 1]} >= hard<0>{};

#ifdef TIMING_DEBUG
    dbg_meta_use_alt = meta_use_alt;
#endif

    // 10. Final prediction mux.
    //     If provider is newly allocated AND meta says use alt AND
    //     alt is a real TAGE hit (not just bimodal) → use alt_pred.
    //     Otherwise → provider_pred (which is bimodal if no TAGE hit,
    //     since match1 falls through to the bimodal bit).
    val<1> use_alt = provider_weak & meta_use_alt & has_alt;
    val<PRED_BITS> final_pred = select(use_alt, alt_pred, provider_pred);

#ifdef TIMING_DEBUG
    dbg_use_alt = use_alt;
    dbg_final_pred = val<1>{final_pred};
#endif

    // 11. Save old prediction (block B) before scatter overwrites with B+1.
    branch_dir.fanout(hard<3>{}); // badpred + true_block + hist_input(BOTH)
    arr<val<1>, N> old_pred = [&](u64 i) -> val<1> {
      return val<1>{pred[i]};
    };

    // 12. Scatter PRED_BITS into per-branch prediction regs.
    //     pred[0] = LSB = branch 0's prediction, pred[I] = bit I.
    static_loop<N>([&]<u64 I>() { pred[I] = final_pred >> I; });

    // ================================================================
    // Training (for block B, using old_pred[] and branch_dir[])
    //
    // old_pred[i] = what we predicted for branch i of block B
    // branch_dir[i] = actual direction of branch i (from update_condbr)
    // ================================================================

    // ---- No conditional branches: update history, skip training ----
    // (matches gshareN_ahead_best pattern: history always advances,
    //  true_block=1 since there's nothing to mispredict)
    if (num_branch == 0) {
      val<PATHBITS> path_bits = val<PATHBITS>{block_end_info.next_pc >> 2};
      path_bits.fanout(hard<NT * 2 + 1>{});
      gh.fanout(hard<NT * 2 + 1>{});
      static_loop<NT>([&]<u64 I>() {
        auto &t = std::get<I>(tables);
        t.fold_idx.update(gh, hard<TableCfg::HIST_LEN[I]>{}, path_bits);
        t.fold_tag.update(gh, hard<TableCfg::HIST_LEN[I]>{}, path_bits);
      });
      gh.update(path_bits);
      last_condbr_dir = 0;
      true_block = 1;
      return;
    }

    // ---- Step 1: Correctness ----
    // badpred[i] = 1 when prediction for branch i was wrong
    // (used by counter/alloc updates in later steps)
    val<1> &mispredict = block_end_info.is_mispredict;
    need_extra_cycle(mispredict);
    [[maybe_unused]] arr<val<1>, N> badpred = [&](u64 i) -> val<1> {
      return old_pred[i] ^ val<1>{branch_dir[i]};
    };

    // ---- Step 8: History update ----
    // true_block uses framework's mispredict signal (not our computed
    // correct_pred) to avoid timing bleed from old_pred reg reads
    true_block =
        ~mispredict | val<1>{branch_dir[num_branch - 1]} | line_end();
    true_block.fanout(hard<NT * 2 + 3>{});

#ifdef TIMING_DEBUG
    dbg_mispredict = mispredict;
    dbg_true_block = true_block;
#endif

    // Compute new fold values OUTSIDE execute_if — runs in parallel with
    // true_block gate, so timing is max(fold_computation, true_block)
    // instead of additive (true_block + fold_computation).
    auto hist_input = [&]() {
      if constexpr (HIST_MODE == HistUpdate::PATH)
        return val<PATHBITS>{block_end_info.next_pc >> 2};
      else if constexpr (HIST_MODE == HistUpdate::DIR)
        return val<1>{branch_dir[num_branch - 1]};
      else // BOTH: concat(direction, path)
        return concat(val<1>{branch_dir[num_branch - 1]},
                      val<PATHBITS>{block_end_info.next_pc >> 2});
    }();

    hist_input.fanout(hard<NT * 2 + 1>{});
    gh.fanout(hard<NT * 2 + 1>{});

#ifdef TIMING_DEBUG
    dbg_hist_input = val<1>{hist_input};
    dbg_gh_fanout = val<1>{gh[0]};
#endif

    // Per-table: compute new fold values, apply with mux (no execute_if gate).
    // select(true_block, new, old) avoids the execute_if timing bleed —
    // both paths resolve in parallel, mux adds ~10ps constant overhead.
    static_loop<NT>([&]<u64 I>() {
      auto &t = std::get<I>(tables);
      auto new_idx = t.fold_idx.compute_update(
          gh, hard<TableCfg::HIST_LEN[I]>{}, hist_input);
      auto new_tag = t.fold_tag.compute_update(
          gh, hard<TableCfg::HIST_LEN[I]>{}, hist_input);
#ifdef TIMING_DEBUG
      if constexpr (I == 0) {
        dbg_fold_compute = val<1>{new_idx};
      }
#endif
      t.fold_idx.apply_update(new_idx, true_block);
      t.fold_tag.apply_update(new_tag, true_block);
#ifdef TIMING_DEBUG
      if constexpr (I == 0) {
        dbg_fold_apply = val<1>{t.fold_idx.get()};
      }
#endif
    });
    gh.update(hist_input, true_block);
  }
};
