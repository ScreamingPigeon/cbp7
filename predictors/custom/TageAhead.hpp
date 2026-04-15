#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"
#include "custom_common.hpp"
#ifdef TAGE_MONITOR
#include "TAMonitor.hpp"
#endif

using namespace hcm;

// ============================================================================
// Table Config
// ============================================================================

template <u64 N = 8, u64 SIZE = 1024, u64 TAG = 11, u64 MINH = 4,
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
// Fallback on secondary tag miss: bimodal or gshare (ahead-pipelined).
//
// Self-contained: no dependency on common.hpp or TageTable.hpp.
// Follows gshareN_ahead_best pipeline pattern.
// ============================================================================

template <
    typename TableCfg = TATableConfig<>,
    u64 N = 8,                   // max conditional branches per block
    u64 PATHBITS = 6,            // bits of next_pc injected into history
    u64 SEC_TAG_BITS = 3,        // secondary tag width (ahead ambiguity)
    bool USE_SEC_TAG = true,     // enable secondary tag matching
    u64 CTR_WIDTH = 1,           // prediction counter width per lane
    u64 HYST_WIDTH = 2,          // hysteresis width (separate from ctr)
    u64 U_WIDTH = 1,             // usefulness counter width
    u64 FB_CAPACITY = 8192, // fallback table size (bimodal or gshare)
    bool USE_GSHARE = false,     // use gshare base (PC^history) vs bimodal (PC)
    u64 GS_HIST = 16,           // gshare history length (only when USE_GSHARE)
    u64 META_WIDTH = 4,          // meta counter width (provider vs alt)
    u64 META_CAPACITY = 256,     // meta table entries
    u64 META_PIPE = 2,           // meta pipeline depth
    u64 LINEINST = 1024,         // line size in instructions
    bool SHARED_HYS = true,      // shared hyst: 2 entries share 1 counter
    HistUpdate HIST_MODE =
        HistUpdate::PATH, // what goes into history: PATH, DIR, or BOTH
    // ---- Global pressure counters ----
    u64 ACC_WIDTH = 4,   // accuracy counter width
    u64 ALLOC_WIDTH = 4, // alloc pressure counter width
    // ---- Probabilistic u-bit decay ----
    bool DECAY_ENABLE = false, DecayMiss DECAY_MISS = DecayMiss::TAG_OR_SEC,
    DecayOp DECAY_OP = DecayOp::DECREMENT,
    auto DECAY_LFSR_WIDTHS = ta::uniform_array<u64, TableCfg::NUM_TABLES>(8),
    typename DecayThreshFn = ta::DefaultDecayThresh,
    // ---- Epoch-based u-bit reset ----
    bool EPOCH_ENABLE = true, typename EpochTriggerFn = ta::DefaultEpochTrigger,
    bool EPOCH_RESET_ACC = false, // reset acc_ctr on epoch fire
    bool EPOCH_RESET_ALLOC = true // reset alloc_ctr on epoch fire
    >
struct TageAhead : predictor {

  // ======== Derived Constants ========

  static_assert(!(DECAY_ENABLE && EPOCH_ENABLE),
                "DECAY_ENABLE and EPOCH_ENABLE are mutually exclusive");

  static constexpr u64 NT = TableCfg::NUM_TABLES;
  static constexpr u64 LOGLINEINST = ta::clog2(LINEINST);
  static constexpr u64 MAXHIST = TableCfg::MAXHIST;
  static constexpr u64 MAX_TAG_WIDTH = ta::array_max(TableCfg::TAG_WIDTH);
  static constexpr u64 MAX_TABLE_SIZE = ta::array_max(TableCfg::TABLE_SIZE);
  static constexpr u64 MAX_IDX_BITS = ta::clog2(MAX_TABLE_SIZE);

  static constexpr u64 MATCH_BITS = NT + 1; // NT tables + fallback

  // Per-bit gh fanout: only fanout bits each table actually reads
  static constexpr auto GH_FANOUT = []() {
    auto fo = ta::gh_per_bit_fanout<MAXHIST, NT, TableCfg::HIST_LEN>();
    if constexpr (USE_GSHARE)
      fo[GS_HIST - 1] += 1; // fb_fold reads gh[GS_HIST-1]
    return fo;
  }();

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
  hcm::ram<val<N>, FB_CAPACITY> fb_ctr{"fb"};
  ta_rwram<META_WIDTH, META_CAPACITY, 2> meta_ctr{"meta"};

  // ======================================================================
  // Registers
  // ======================================================================
  // ---- Global History (shared, folds live in per-table TATable) ----
  ta_global_history<MAXHIST> gh;

  // ---- Fallback Predictor (ahead-pipelined) ----
  // USE_GSHARE=false: bimodal (PC-indexed)
  // USE_GSHARE=true:  gshare (PC ^ folded_history indexed)
  reg<PRED_BITS> prefetch_fb;
  reg<PRED_BITS> current_fb;
  static constexpr u64 FB_IDX_BITS = ta::clog2(FB_CAPACITY);
  reg<FB_IDX_BITS> prefetch_fb_idx;
  reg<FB_IDX_BITS> current_fb_idx;

  // Gshare fold register — folds GS_HIST bits of global history into
  // FB_IDX_BITS for the fallback index. Zero cost when USE_GSHARE=false
  // (constexpr-gated, no reads/writes occur).
  ta_folded_gh<FB_IDX_BITS> fb_fold;

  // Piped PC for allocation tag recomputation (stores inst_pc >> 2)
  static constexpr u64 ALLOC_PC_BITS = MAX_TAG_WIDTH + 2;
  reg<ALLOC_PC_BITS> prefetch_pc;
  reg<ALLOC_PC_BITS> current_pc;

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
  reg<MAX_TAG_WIDTH> prefetch_ctag[NT]; // computed tag for allocation piping

  // current_*: shifted from prefetch, used for prediction
  reg<MAX_TAG_WIDTH> current_tag[NT];
  reg<1> current_tag_hit[NT];
  reg<PRED_BITS> current_pred[NT];
  reg<SEC_TAG_BITS> current_sec[NT];
  reg<MAX_IDX_BITS> current_idx[NT];
  reg<std::max(u64(1), HYST_WIDTH)> current_hyst[NT];
  reg<U_WIDTH> current_u[NT];
  reg<MAX_TAG_WIDTH> current_ctag[NT]; // piped computed tag for allocation

  // train_*: saved from current_* before pipeline shift, used for training
  // (training is for block A, resolution is for block B)
  reg<MAX_IDX_BITS> train_idx[NT];
  reg<std::max(u64(1), HYST_WIDTH)> train_hyst[NT];
  reg<U_WIDTH> train_u[NT];
  reg<PRED_BITS> train_fb;
  reg<FB_IDX_BITS> train_fb_idx;
  reg<ALLOC_PC_BITS> train_pc;
  reg<MAX_TAG_WIDTH> train_ctag[NT]; // piped computed tag for allocation

  // Piped resolution values from previous update_cycle (block A's resolution)
  reg<MATCH_BITS> train_match1;
  reg<PRED_BITS> train_provider_pred;
  reg<1> train_provider_weak;
  reg<1> train_altdiff;

  // Guard: skip training until piped resolution regs have been populated
  reg<1> train_valid = 0;

  // ---- U-bit decay (probabilistic) ----
  // Piped tag/sec hit for decay miss detection
  reg<1> train_tag_hit[NT];
  reg<1> train_sec_hit[NT];

  // Global pressure counters (always declared, zero-cost when unused)
  reg<ACC_WIDTH> acc_ctr;
  reg<ALLOC_WIDTH> alloc_ctr;

  // Per-table LFSRs (varying widths via tuple)
  static constexpr u64 MAX_LFSR_WIDTH =
      DECAY_ENABLE ? ta::array_max(DECAY_LFSR_WIDTHS) : 1;
  typename ta::TALfsrTuple<DECAY_LFSR_WIDTHS,
                           std::make_index_sequence<NT>>::type decay_lfsrs;

  // ---- Prediction reg (shared by
  // predict1/reuse_predict1/predict2/reuse_predict2) ----
  arr<reg<1>, N> pred;

  // ---- Meta pipeline (shifted each update_cycle) ----
  static constexpr u64 META_IDX_BITS = ta::clog2(META_CAPACITY);
  reg<META_WIDTH, i64> meta_pipe[META_PIPE];
  reg<META_IDX_BITS> meta_idx_pipe[META_PIPE];

// ---- Timing debug taps (zero cost in normal builds) ----
#ifdef TIMING_DEBUG
  reg<1> dbg_full_hits;     // after per-table hit computation
  reg<1> dbg_fb_pred;      // after fallback read
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
  // Resolution gaps
  reg<1> dbg_altdiff;        // provider_pred ^ alt_pred
  reg<1> dbg_actual_dir;     // branch_dir scatter
  reg<1> dbg_any_prov_wrong; // t_pp ^ actual_dir
  // Allocation chain
  reg<1> dbg_alloc_base;    // t_match1 - 1
  reg<1> dbg_notumask;      // u_zero concat
  reg<1> dbg_candallocmask; // misp_alloc & notumask
  reg<1> dbg_alloc_target;  // one_hot of candidates
  reg<1> dbg_noalloc;       // candallocmask == 0
  reg<1> dbg_uclearmask;    // u-clear fallback
  // Train write chain (per-table)
  reg<1> dbg_fb_write;      // fb_ctr.write gate
  reg<1> dbg_meta_write;     // meta_ctr.write gate
  reg<1> dbg_pred_write[NT]; // pred_ram.write gate per table
  reg<1> dbg_hyst_write[NT]; // hyst_ram.write gate per table
  reg<1> dbg_tag_write[NT];  // tag_ram.write gate per table
  reg<1> dbg_u_write[NT];    // u_ram.write gate per table
  // Decay (per-table)
  reg<1> dbg_decay_fire[NT];   // decay_fire per table
  reg<1> dbg_decay_merged[NT]; // merged u value per table
  reg<1> dbg_epoch_fire;       // epoch trigger
#endif

#ifdef TAGE_MONITOR
  TAMonitor<NT, N, MAX_TABLE_SIZE, USE_GSHARE, FB_CAPACITY> mon;
  ~TageAhead() { mon.print_summary(); }
#endif

  // ======================================================================
  // Helpers
  // ======================================================================

  val<1> line_end() { return (block_entry + block_size) == hard<LINEINST>{}; }

  // ======================================================================
  // predict1/2, reuse_predict1/2, update_condbr, update_cycle — TODO
  // ======================================================================

  val<1> predict1([[maybe_unused]] val<64> inst_pc) {
#ifdef TAGE_MONITOR
    mon.record_predict1();
#endif
    inst_pc.fanout(
        hard<2 * NT + 1>{}); // 2 reads per table (>>2, >>4) + fb (>>2)

#ifdef TIMING_DEBUG
    dbg_inst_pc = val<1>{inst_pc};
#endif

    // Ahead reads for next block — run unconditionally (no true_block gate).
    //
    // Why no execute_if(true_block)?
    //   true_block is computed in update_cycle at +90ps (depends on mispredict,
    //   branch_dir, and line_end). Gating the RAM reads on it would delay the
    //   entire prefetch path: RAM addresses (fold_idx at -74, inst_pc at -255)
    //   are ready long before true_block, but execute_if holds the reads until
    //   the gate resolves. This adds ~165ps to pipe_shift, pushing the
    //   downstream resolution chain (one_hot → fold_or → select → scatter)
    //   well past the 1-cycle target.
    //
    // Why is this safe?
    //   When true_block=0 (mispredicted taken branch mid-line), the framework
    //   fires extra_cycle, which re-invokes predict1 with the corrected PC.
    //   The prefetch regs from the spurious read are unconditionally
    //   overwritten in that corrected predict1 call before they ever shift into
    //   current_*. In hardware, the RAM reads happen regardless — execute_if
    //   only gates the output latch, so removing it has zero area/power impact.
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
      prefetch_ctag[I] = computed_tag;
      prefetch_pred[I] = t.pred_ram.read(idx);
      if constexpr (USE_SEC_TAG)
        prefetch_sec[I] = t.sec_ram.read(idx);
      prefetch_idx[I] = idx;
      prefetch_hyst[I] = t.hyst_ram.read(val<t.HYST_IDX_BITS>{idx});
      prefetch_u[I] = t.u_ram.read(idx);
    });

    // Fallback ahead read (direct-mapped, no tag match needed)
    // USE_GSHARE: index = PC ^ folded_history; bimodal: index = PC
    auto fb_idx = [&]() {
      if constexpr (USE_GSHARE)
        return val<FB_IDX_BITS>{inst_pc >> 2} ^ fb_fold.get();
      else
        return val<FB_IDX_BITS>{inst_pc >> 2};
    }();
    prefetch_fb_idx = fb_idx;
    prefetch_fb = fb_ctr.read(fb_idx);
    prefetch_pc = val<ALLOC_PC_BITS>{inst_pc >> 2};

    // Crit path: just read precomputed prediction from reg
    block_entry.fanout(hard<2 * LINEINST>{}); // read in line_end() across
                                              // predict + reuse + update_condbr
    pred.fanout(hard<2 * LINEINST + 1>{});    // predict1/2, reuse + 1 training
                                              // old_pred read
    block_size = 1;
    num_branch = 0;
    reuse_prediction(~line_end());
#ifdef TIMING_DEBUG
    dbg_p1_return = val<1>{pred[num_branch]};
#endif
    return pred[num_branch];
  }
  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) {
#ifdef TAGE_MONITOR
    mon.record_reuse_predict1();
#endif
    block_size++;
    reuse_prediction(~line_end());
    return pred[num_branch];
  }

  val<1> predict2([[maybe_unused]] val<64> inst_pc) {
#ifdef TAGE_MONITOR
    mon.record_predict2();
#endif
    return pred[num_branch];
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) {
#ifdef TAGE_MONITOR
    mon.record_reuse_predict2();
#endif
    return pred[num_branch];
  }

  void update_condbr([[maybe_unused]] val<64> branch_pc, val<1> taken,
                     [[maybe_unused]] val<64> next_pc) {
    assert(num_branch < N);
    branch_dir[num_branch] = taken.fo1();
#ifdef TAGE_MONITOR
    mon.record_branch_pc(static_cast<u64>(branch_pc));
#endif
    num_branch++;
    reuse_prediction(~line_end() & val<1>{num_branch < N});
  }

  void update_cycle([[maybe_unused]] instruction_info &block_end_info) {
    // std::cerr << "=== ENTER update_cycle (num_branch=" << num_branch << ")
    // ===\n";
    // ---- Prefetch part ------
    // 1. Pipeline shift: prefetch → current (unconditional — only true blocks
    // reach here)

    // Save current_* into train_* before shift (block A's data for training)
    static_loop<NT>([&]<u64 I>() {
      train_idx[I] = current_idx[I];
      train_hyst[I] = current_hyst[I];
      train_u[I] = current_u[I];
      train_ctag[I] = current_ctag[I];
      train_tag_hit[I] = current_tag_hit[I];
      if constexpr (USE_SEC_TAG)
        train_sec_hit[I] = (val<SEC_TAG_BITS>{current_sec[I]} ==
                            val<SEC_TAG_BITS>{curr_sec_tag});
    });
    train_fb = current_fb;
    train_fb_idx = current_fb_idx;
    train_pc = current_pc;

    static_loop<NT>([&]<u64 I>() {
      current_tag[I] = prefetch_tag[I];
      current_tag_hit[I] = prefetch_tag_hit[I];
      current_pred[I] = prefetch_pred[I];
      if constexpr (USE_SEC_TAG)
        current_sec[I] = prefetch_sec[I];
      current_idx[I] = prefetch_idx[I];
      current_hyst[I] = prefetch_hyst[I];
      current_u[I] = prefetch_u[I];
      current_ctag[I] = prefetch_ctag[I];
    });
    current_fb = prefetch_fb;
    current_fb_idx = prefetch_fb_idx;
    current_pc = prefetch_pc;

    // Precompute secondary tag for next block
    if constexpr (USE_SEC_TAG) {
      block_end_info.next_pc.fanout(
          hard<3>{}); // curr_sec_tag + meta_idx + hist path_bits
      curr_sec_tag = val<SEC_TAG_BITS>{block_end_info.next_pc >> 2};
    } else {
      block_end_info.next_pc.fanout(hard<2>{}); // meta_idx + hist path_bits
    }

    // ================================================================
    // Provider / altpred resolution via bitmask + one_hot
    //
    // We cannot reassign val<N> (private operator=), so we avoid
    // accumulation loops. Instead:
    //   1. Compute per-table hit bits → arr<val<1>, NT>
    //   2. Concat into val<NT+1> with fallback as MSB always-hit
    //   3. one_hot() → lowest set bit = longest-history hit = provider
    //   4. one_hot() on remainder → alt
    //   5. Replicate one-hot bits to PRED_BITS width, AND with each
    //      table's prediction, fold_or → extract provider/alt pred
    // ================================================================

    // 2. Per-table hit: primary tag matched (off crit path) AND
    //    optionally stored secondary tag matches curr_sec_tag
    if constexpr (USE_SEC_TAG)
      curr_sec_tag.fanout(hard<NT>{}); // compared once per table
    arr<val<1>, NT> full_hits = [&](u64 i) {
      if constexpr (USE_SEC_TAG) {
        return val<1>{current_tag_hit[i]} &
               (val<SEC_TAG_BITS>{current_sec[i]} ==
                val<SEC_TAG_BITS>{curr_sec_tag});
      } else {
        return val<1>{current_tag_hit[i]};
      }
    };

#ifdef TAGE_MONITOR
    for (u64 i = 0; i < NT; i++) {
      if constexpr (USE_SEC_TAG) {
        mon.record_tag_lookup(
            i, static_cast<u64>(current_tag_hit[i]),
            static_cast<u64>(val<SEC_TAG_BITS>{current_sec[i]} ==
                             val<SEC_TAG_BITS>{curr_sec_tag}));
      } else {
        mon.record_tag_lookup(i, static_cast<u64>(current_tag_hit[i]), false);
      }
    }
#endif

#ifdef TIMING_DEBUG
    dbg_full_hits = full_hits[0];
#endif

    // 3. Fallback — ahead-pipelined, already in current_fb from pipe shift.
    val<PRED_BITS> fb_pred = val<PRED_BITS>{current_fb};

#ifdef TIMING_DEBUG
    dbg_fb_pred = val<1>{fb_pred};
#endif

    // 4. Build match bitmask.
    //    Bit layout of val<NT+1>:
    //      bit 0     = table 0 (longest history)
    //      bit NT-1  = table NT-1 (shortest history)
    //      bit NT    = fallback (always 1)
    //    one_hot() returns lowest set bit → longest-history hit = provider.
    //    Second one_hot() on (match ^ match1) → alt provider.
    val<MATCH_BITS> match = concat(val<1>{1}, full_hits.concat());
    match.fanout(hard<2>{}); // one_hot + XOR
    val<MATCH_BITS> match1 = match.one_hot();
    match1.fanout(hard<3>{}); // make_array + XOR with match + alloc_base
    val<MATCH_BITS> match2 = (match ^ match1).one_hot();
    match2.fanout(hard<2>{}); // make_array + has_alt mask

#ifdef TIMING_DEBUG
    dbg_match = val<1>{match};
    dbg_match1 = val<1>{match1};
    dbg_match2 = val<1>{match2};
#endif

    // 5. Prediction array: one PRED_BITS-wide entry per table + fallback.
    //    table_preds[0..NT-1] = TAGE tables, table_preds[NT] = fallback.
    arr<val<PRED_BITS>, NT + 1> table_preds = [&](u64 i) -> val<PRED_BITS> {
      if (i < NT)
        return val<PRED_BITS>{current_pred[i]};
      return val<PRED_BITS>{fb_pred};
    };

    // 6. Extract provider and alt predictions.
    //    For each table: replicate its one-hot match bit to PRED_BITS width,
    //    AND with that table's prediction. Since match1 is one-hot, exactly
    //    one table contributes non-zero bits. fold_or collapses to that pred.
    arr<val<1>, NT + 1> m1_bits = match1.make_array(val<1>{});
    arr<val<1>, NT + 1> m2_bits = match2.make_array(val<1>{});
    m1_bits.fanout(hard<6>{});
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
    //    Only check the provider table (mask by match1 bit). Fallback
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

    // 8. has_alt: does the alt match point to a TAGE table (not fallback)?
    //    Mask out the fallback bit (MSB) and check if anything remains.
    val<1> has_alt =
        (match2 & val<MATCH_BITS>{(u64(1) << NT) - 1}) != hard<0>{};

#ifdef TIMING_DEBUG
    dbg_provider_weak = provider_weak;
    dbg_has_alt = has_alt;
#endif

    // 9. Meta counter: predicts whether to trust a newly-allocated provider
    //    or fall back to alt. Read from PC-indexed RAM, shifted through a
    //    META_PIPE-stage pipeline. Sign bit of oldest stage decides.
    auto meta_idx = val<META_IDX_BITS>{block_end_info.next_pc >> 2};
    meta_pipe[0] = meta_ctr.read(meta_idx);
    meta_idx_pipe[0] = meta_idx;
    for (u64 i = META_PIPE - 1; i > 0; i--) {
      meta_pipe[i] = meta_pipe[i - 1];
      meta_idx_pipe[i] = meta_idx_pipe[i - 1];
    }
    val<1> meta_use_alt =
        val<META_WIDTH, i64>{meta_pipe[META_PIPE - 1]} >= hard<0>{};

#ifdef TIMING_DEBUG
    dbg_meta_use_alt = meta_use_alt;
#endif

    // 10. Final prediction mux.
    //     If provider is newly allocated AND meta says use alt AND
    //     alt is a real TAGE hit (not just fallback) → use alt_pred.
    //     Otherwise → provider_pred (which is fallback if no TAGE hit,
    //     since match1 falls through to the fallback bit).
    val<1> use_alt = provider_weak & meta_use_alt & has_alt;
    val<PRED_BITS> final_pred = select(use_alt, alt_pred, provider_pred);

#ifdef TIMING_DEBUG
    dbg_use_alt = use_alt;
    dbg_final_pred = val<1>{final_pred};
#endif

    // 11. Save old prediction (block B) before scatter overwrites with B+1.
    branch_dir.fanout(
        hard<4>{}); // badpred + true_block + hist_input + actual_dir
    arr<val<1>, N> old_pred = [&](u64 i) -> val<1> { return val<1>{pred[i]}; };

    // 12. Scatter PRED_BITS into per-branch prediction regs.
    //     pred[0] = LSB = branch 0's prediction, pred[I] = bit I.
    static_loop<N>([&]<u64 I>() { pred[I] = final_pred >> I; });

#ifdef TAGE_MONITOR
    // Record resolution for each branch in this block
    for (u64 r = 0; r < num_branch; r++) {
      bool meta_overrode =
          static_cast<u64>(provider_weak) && static_cast<u64>(has_alt);
      bool meta_chose = static_cast<u64>(use_alt);
      bool pred_taken = (static_cast<u64>(final_pred) >> r) & 1;
      mon.record_prediction(r, static_cast<u64>(match1),
                            static_cast<u64>(match2), meta_overrode, meta_chose,
                            pred_taken);
    }
#endif

    // ================================================================
    // Training (for block A, using old_pred[], branch_dir[],
    // and piped resolution from previous cycle's train_* regs)
    //
    // old_pred[i] = what we predicted for branch i of block A
    // branch_dir[i] = actual direction of branch i (from update_condbr)
    // t_m1/t_pp/t_pw/t_ad = block A's resolution (piped from prev cycle)
    // train_idx/hyst/u/bim = block A's prefetched data
    // ================================================================

    // Read OLD piped resolution values BEFORE overwriting with current
    // resolution
    val<MATCH_BITS> t_match1 = train_match1;
    arr<val<1>, NT + 1> t_m1 = t_match1.make_array(val<1>{});
    val<PRED_BITS> t_pp = train_provider_pred;
    val<1> t_pw = train_provider_weak;
    val<1> t_ad = train_altdiff;

    // Save current resolution → train regs (for NEXT cycle's training)
    alt_pred.fanout(hard<2>{});
    val<1> altdiff = (provider_pred ^ alt_pred) != hard<0>{};
    train_match1 = match1;
    train_provider_pred = provider_pred;
    train_provider_weak = provider_weak;
    train_altdiff = altdiff;
#ifdef TIMING_DEBUG
    dbg_altdiff = val<1>{altdiff};
#endif

    // Read train_valid BEFORE setting it to 1 (regs may be immediate-write)
    val<1> do_train = train_valid;
    train_valid = 1;

    // ---- No conditional branches: update history, skip training ----
    if (num_branch == 0) {
      val<PATHBITS> path_bits = val<PATHBITS>{block_end_info.next_pc >> 2};
      path_bits.fanout(hard<NT * 2 + 1 + (USE_GSHARE ? 1 : 0)>{});
      gh.template fanout_per_bit<GH_FANOUT>();
      static_loop<NT>([&]<u64 I>() {
        auto &t = std::get<I>(tables);
        t.fold_idx.update(gh, hard<TableCfg::HIST_LEN[I]>{}, path_bits);
        t.fold_tag.update(gh, hard<TableCfg::HIST_LEN[I]>{}, path_bits);
      });
      if constexpr (USE_GSHARE)
        fb_fold.update(gh, hard<GS_HIST>{}, path_bits);
      gh.update(path_bits);
      last_condbr_dir = 0;
      true_block = 1;
      return;
    }

    // ---- Step 1: Correctness ----
    val<1> &mispredict = block_end_info.is_mispredict;
    mispredict.fanout(
        hard<5>{}); // extra_cycle + bim + true_block + dbg + alloc
    need_extra_cycle(mispredict);
    do_train.fanout(hard<4 * NT + 3>{}); // gates bim, meta, 4 writes per table

#ifdef TAGE_MONITOR
    mon.record_block(static_cast<u64>(val<LOGLINEINST>{block_entry}),
                     block_size, num_branch, static_cast<u64>(mispredict));
    for (u64 r = 0; r < num_branch; r++)
      mon.record_outcome(r, static_cast<u64>(branch_dir[r]),
                         r == num_branch - 1 && static_cast<u64>(mispredict));
    if (!static_cast<u64>(do_train))
      mon.record_train_skip();
#endif
    [[maybe_unused]] arr<val<1>, N> badpred = [&](u64 i) -> val<1> {
      return old_pred[i] ^ val<1>{branch_dir[i]};
    };

    // ---- Compute training signals ----
    val<PRED_BITS> actual_dir = arr<val<1>, N>{[&](u64 i) -> val<1> {
                                  return val<1>{branch_dir[i]};
                                }}.concat();

    // Provider wrong on any branch? (uses piped provider_pred)
    t_pp.fanout(hard<2>{});
    actual_dir.fanout(hard<NT + 1>{});
    val<1> any_provider_wrong = (t_pp ^ actual_dir) != hard<0>{};
    any_provider_wrong.fanout(hard<3 * NT + 1>{});
#ifdef TIMING_DEBUG
    dbg_actual_dir = val<1>{actual_dir};
    dbg_any_prov_wrong = any_provider_wrong;
#endif
    t_pw.fanout(hard<NT + 1>{});
    t_m1.fanout(hard<6>{});
    t_ad.fanout(hard<NT + 1>{});

    // ---- Allocation masks (needed before merged loop) ----
    val<NT> alloc_base = val<NT>{t_match1 - 1};
    arr<val<1>, NT> u_zero = [&](u64 i) -> val<1> {
      return val<U_WIDTH>{train_u[i]} == hard<0>{};
    };
    val<NT> notumask = u_zero.concat();
    notumask.fanout(hard<2>{});
    val<NT> misp_alloc = alloc_base & mispredict.replicate(hard<NT>{}).concat();
    misp_alloc.fanout(hard<2>{});
    val<NT> candallocmask = misp_alloc & notumask;
    candallocmask.fanout(hard<2>{});
    val<NT> alloc_target = candallocmask.reverse().one_hot().reverse();
    arr<val<1>, NT> allocate = alloc_target.make_array(val<1>{});
    allocate.fanout(hard<6>{});
    val<1> noalloc = (candallocmask == hard<0>{});
    val<NT> uclearmask = misp_alloc & noalloc.replicate(hard<NT>{}).concat();
    arr<val<1>, NT> uclear = uclearmask.make_array(val<1>{});
    uclear.fanout(hard<2>{});
#ifdef TIMING_DEBUG
    dbg_alloc_base = val<1>{alloc_base};
    dbg_notumask = val<1>{notumask};
    dbg_candallocmask = val<1>{candallocmask};
    dbg_alloc_target = val<1>{alloc_target};
    dbg_noalloc = noalloc;
    dbg_uclearmask = val<1>{uclearmask};
#endif

#ifdef TAGE_MONITOR
    {
      u64 at = static_cast<u64>(alloc_target);
      mon.record_allocation(at != 0, at);
      // Provider index for cascade tracking
      u64 prov_idx = TAMonitor<NT, N, MAX_TABLE_SIZE, USE_GSHARE, FB_CAPACITY>::decode_provider(
          static_cast<u64>(t_match1));
      if (static_cast<u64>(mispredict)) {
        if (static_cast<u64>(misp_alloc) == 0)
          mon.record_alloc_blocked(); // fallback was provider, no candidates
        mon.record_alloc_cascade(prov_idx, at);
      }
    }
#endif

    // ---- Step 4: Fallback update (mispredict + fallback is provider only) ----
    val<1> fb_changed = actual_dir != val<PRED_BITS>{train_fb};
    val<1> fb_gate = do_train & t_m1[NT] & mispredict & fb_changed;
    execute_if(fb_gate, [&]() {
      fb_ctr.write(val<FB_IDX_BITS>{train_fb_idx}, actual_dir);
    });
#ifdef TAGE_MONITOR
    if (static_cast<u64>(fb_gate))
      mon.record_fb_write(static_cast<u64>(val<FB_IDX_BITS>{train_fb_idx}));
#endif

    // ---- Step 5: Meta counter update ----
    auto old_meta = val<META_WIDTH, i64>{meta_pipe[META_PIPE - 1]};
    auto new_meta = ta_update_ctr(old_meta, any_provider_wrong);
    val<1> meta_gate = do_train & t_pw & t_ad & (new_meta != old_meta);
    execute_if(meta_gate, [&]() {
      meta_ctr.write(val<META_IDX_BITS>{meta_idx_pipe[META_PIPE - 1]}, new_meta,
                     hard<0>{});
    });
#ifdef TIMING_DEBUG
    dbg_fb_write = fb_gate;
    dbg_meta_write = meta_gate;
#endif

    // ---- Merged per-table writes (one write per RAM per table) ----
    // For each table: alloc takes priority over update. Mux selects data.
    train_pc.fanout(hard<NT>{});
    static_loop<NT>([&]<u64 I>() {
      auto &t = std::get<I>(tables);
      val<1> do_alloc = allocate[I];

      // pred_ram: alloc writes actual_dir, update writes actual_dir → same data
      val<1> do_pred_update = t_m1[I] & t_pw & any_provider_wrong;
      execute_if(do_train & (do_alloc | do_pred_update), [&]() {
        t.pred_ram.write(val<t.IDX_BITS>{train_idx[I]}, actual_dir, hard<0>{});
      });

      // hyst_ram: alloc writes 0, update writes new_hyst → mux on do_alloc
      constexpr u64 HW = std::max(u64(1), HYST_WIDTH);
      auto old_hyst = val<HW>{train_hyst[I]};
      auto new_hyst = ta_update_ctr(old_hyst, ~any_provider_wrong);
      auto hyst_data = select(do_alloc, val<HW>{0}, new_hyst);
      val<1> do_hyst_update = t_m1[I] & (new_hyst != old_hyst);
      execute_if(do_train & (do_alloc | do_hyst_update), [&]() {
        t.hyst_ram.write(val<t.HYST_IDX_BITS>{train_idx[I]}, hyst_data,
                         hard<0>{});
      });

      // tag_ram + sec_ram: alloc only (plain RAM, protected by extra_cycle)
      execute_if(do_train & do_alloc, [&]() {
        // Use piped computed tag (from predict1 time, not current folds)
        t.tag_ram.write(val<t.IDX_BITS>{train_idx[I]},
                        val<t.tag_width>{train_ctag[I]});
        if constexpr (USE_SEC_TAG)
          t.sec_ram.write(val<t.IDX_BITS>{train_idx[I]},
                          val<SEC_TAG_BITS>{curr_sec_tag});
      });
#ifdef TAGE_MONITOR
      if (static_cast<u64>(do_train & do_alloc))
        mon.record_tage_write(I,
                              static_cast<u64>(val<t.IDX_BITS>{train_idx[I]}));
#endif

      // u_ram: combined provider update + allocation + uclear + decay
      val<U_WIDTH> base_newu = val<U_WIDTH>{~any_provider_wrong} &
                               val<U_WIDTH>{~allocate[I]} &
                               val<U_WIDTH>{~uclear[I]};
      val<1> base_u_write = (t_m1[I] & t_ad) | allocate[I] | uclear[I];

      // Probabilistic decay: on tag/sec miss, LFSR < threshold → decay u
      auto [newu, u_write] = [&]() -> std::pair<val<U_WIDTH>, val<1>> {
        if constexpr (!DECAY_ENABLE) {
          return {base_newu, base_u_write};
        } else {
          constexpr u64 LW = DECAY_LFSR_WIDTHS[I];
          auto &lfsr = std::get<I>(decay_lfsrs);

          // Miss condition from piped tag/sec hit
          val<1> tag_missed = ~val<1>{train_tag_hit[I]};
          val<1> decay_miss = [&]() {
            if constexpr (!USE_SEC_TAG || DECAY_MISS == DecayMiss::TAG)
              return tag_missed;
            else {
              val<1> sec_missed = ~val<1>{train_sec_hit[I]};
              if constexpr (DECAY_MISS == DecayMiss::SEC)
                return sec_missed;
              else if constexpr (DECAY_MISS == DecayMiss::TAG_OR_SEC)
                return tag_missed | sec_missed;
              else
                return tag_missed & sec_missed;
            }
          }();

          // Threshold from global counters
          auto thresh = DecayThreshFn::template compute<I, LW>(
              val<ACC_WIDTH>{acc_ctr}, val<ALLOC_WIDTH>{alloc_ctr});

          // LFSR fires when below threshold, on miss, not allocating
          val<1> decay_fire =
              decay_miss & ~allocate[I] & (val<LW>{lfsr} < thresh);

          // Apply decay op to train_u
          val<U_WIDTH> old_u = val<U_WIDTH>{train_u[I]};
          val<U_WIDTH> decayed_u = [&]() {
            if constexpr (DECAY_OP == DecayOp::DECREMENT)
              return select(old_u == hard<0>{}, old_u, val<U_WIDTH>{old_u - 1});
            else if constexpr (DECAY_OP == DecayOp::HALVE)
              return val<U_WIDTH>{old_u >> 1};
            else // CLEAR
              return val<U_WIDTH>{0};
          }();

          // Mux: if base write active, use base_newu; else if decay, use
          // decayed_u
          val<U_WIDTH> merged = select(base_u_write, base_newu,
                                       select(decay_fire, decayed_u, old_u));
          val<1> merged_write = base_u_write | decay_fire;
          val<1> merged_changed = merged != old_u;
          return {merged, merged_write & merged_changed};
        }
      }();

      execute_if(do_train & u_write, [&]() {
        t.u_ram.write(val<t.IDX_BITS>{train_idx[I]}, newu, hard<0>{});
      });

#ifdef TAGE_MONITOR
      if (static_cast<u64>(do_train & u_write))
        mon.record_u_write(I, static_cast<u64>(newu) != 0);
      if constexpr (DECAY_ENABLE) {
        // decay_fire = u_write & ~base_u_write (only non-base u write is decay)
        if (static_cast<u64>(do_train & u_write & ~base_u_write))
          mon.record_decay_fire();
      }
#endif
#ifdef TIMING_DEBUG
      dbg_pred_write[I] = do_train & (do_alloc | do_pred_update);
      dbg_hyst_write[I] = do_train & (do_alloc | do_hyst_update);
      dbg_tag_write[I] = do_train & do_alloc;
      dbg_u_write[I] = do_train & u_write;
      if constexpr (DECAY_ENABLE) {
        dbg_decay_fire[I] = u_write & ~base_u_write;
        dbg_decay_merged[I] = val<1>{newu};
      }
#endif
    });

    // ---- Global pressure counter updates ----
    if constexpr (DECAY_ENABLE || EPOCH_ENABLE) {
      // Accuracy counter: increment on correct, decrement on mispredict
      auto new_acc = ta_update_ctr(val<ACC_WIDTH>{acc_ctr}, ~mispredict);

      // Alloc pressure counter: increment when no alloc slot found
      val<1> any_alloc = alloc_target != hard<0>{};
      auto new_alloc = ta_update_ctr(val<ALLOC_WIDTH>{alloc_ctr}, ~any_alloc);

      if constexpr (EPOCH_ENABLE) {
        // Epoch: bulk reset u_ram when trigger fires
        val<1> epoch_fire =
            EpochTriggerFn::template should_fire<ACC_WIDTH, ALLOC_WIDTH>(
                val<ACC_WIDTH>{acc_ctr}, val<ALLOC_WIDTH>{alloc_ctr});
        execute_if(epoch_fire, [&]() {
          static_loop<NT>([&]<u64 I>() { std::get<I>(tables).u_ram.reset(); });
        });
#ifdef TIMING_DEBUG
        dbg_epoch_fire = epoch_fire;
#endif
#ifdef TAGE_MONITOR
        if (static_cast<u64>(epoch_fire))
          mon.record_epoch_reset();
#endif
        // Optionally reset counters on epoch fire
        if constexpr (EPOCH_RESET_ACC)
          acc_ctr = select(epoch_fire, val<ACC_WIDTH>{0}, new_acc);
        else
          acc_ctr = new_acc;
        if constexpr (EPOCH_RESET_ALLOC)
          alloc_ctr = select(epoch_fire, val<ALLOC_WIDTH>{0}, new_alloc);
        else
          alloc_ctr = new_alloc;
      } else {
        acc_ctr = new_acc;
        alloc_ctr = new_alloc;
      }
    }

#ifdef TAGE_MONITOR
    if constexpr (DECAY_ENABLE || EPOCH_ENABLE)
      mon.record_pressure(static_cast<u64>(acc_ctr),
                          static_cast<u64>(alloc_ctr));
#endif

    // ---- Decay: LFSR tick ----
    if constexpr (DECAY_ENABLE) {
      static_loop<NT>([&]<u64 I>() {
        constexpr u64 LW = DECAY_LFSR_WIDTHS[I];
        auto &lfsr = std::get<I>(decay_lfsrs);
        val<LW> old_lfsr = val<LW>{lfsr};
        val<1> feedback = old_lfsr & hard<1>{};
        val<LW> shifted = val<LW>{old_lfsr >> 1};
        val<LW> tap_mask = val<LW>{u64(1) << (LW - 1)};
        lfsr = shifted ^ (tap_mask & feedback.replicate(hard<LW>{}).concat());
      });
    }

    // ---- Step 8: History update ----
    // true_block uses framework's mispredict signal (not our computed
    // correct_pred) to avoid timing bleed from old_pred reg reads
    true_block = ~mispredict | val<1>{branch_dir[num_branch - 1]} | line_end();
    true_block.fanout(
        hard<NT * 2 + 2 + (USE_GSHARE ? 1 : 0)>{}); // NT fold_idx + NT fold_tag apply_update muxes +
                             // gh.update + monitor + (fb_fold if gshare)

#ifdef TAGE_MONITOR
    if (static_cast<u64>(true_block))
      mon.record_true_block();
#endif

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

    hist_input.fanout(hard<NT * 2 + 1 + (USE_GSHARE ? 1 : 0)>{});
    gh.template fanout_per_bit<GH_FANOUT>();

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
    if constexpr (USE_GSHARE) {
      auto new_fb = fb_fold.compute_update(gh, hard<GS_HIST>{}, hist_input);
      fb_fold.apply_update(new_fb, true_block);
    }
    gh.update(hist_input, true_block);
    // std::cerr << "=== EXIT update_cycle ===\n";
  }
};
