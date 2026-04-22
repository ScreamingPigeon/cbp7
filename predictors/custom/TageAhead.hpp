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

template <u64 N = 7, u64 SIZE = 1024, u64 TAG = 8, u64 MINH = 10,
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
    u64 N = 8,               // max conditional branches per block
    u64 PATHBITS = 6,        // bits of next_pc injected into history
    u64 SEC_TAG_BITS = 3,    // secondary tag width (ahead ambiguity)
    bool USE_SEC_TAG = true, // enable secondary tag matching
    u64 NUM_PATHS = 1, // parallel resolution chains (paper: 1 << SEC_TAG_BITS)
    typename SecTagHashFn =
        ta::DefaultSecTagHash, // sec_tag hash: PC → val<SEC_TAG_BITS>
    u64 CTR_WIDTH = 1,         // prediction counter width per lane
    u64 HYST_WIDTH = 2,        // hysteresis width (separate from ctr)
    u64 U_WIDTH = 1,           // usefulness counter width
    u64 FB_CAPACITY = 8192,    // fallback table size (bimodal or gshare)
    bool USE_GSHARE = true,    // use gshare base (PC^history) vs bimodal (PC)
    u64 GS_HIST = 6,           // gshare history length (only when USE_GSHARE)
    u64 META_WIDTH = 6,        // meta counter width (provider vs alt)
    u64 META_CAPACITY = 1024,  // meta table entries
    u64 META_PIPE = 2,         // meta pipeline depth
    u64 LINEINST = 1024,       // line size in instructions
    bool SHARED_HYS = true,    // shared hyst: 2 entries share 1 counter
    HistUpdate HIST_MODE =
        HistUpdate::PATH, // what goes into history: PATH, DIR, or BOTH
    // ---- Allocation policy ----
    typename AllocCfg = TADefaultAllocConfig,
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
    bool EPOCH_RESET_ACC = false,  // reset acc_ctr on epoch fire
    bool EPOCH_RESET_ALLOC = true, // reset alloc_ctr on epoch fire
    // ---- Fallback reconciliation ----
    // When enabled, tracks agreement between fb and TAGE via fb_hyst RAM.
    // Overwrites fb pred when they persistently disagree (hyst weak).
    // Mirrors Tage.hpp's P1/P2 reconciliation — keeps fb aligned with TAGE.
    bool FB_RECONCILE = false,
    // ---- Far-allocation pressure ----
    // When > 0, allocation distance >= FARALLOC_DIST from provider biases
    // alloc_ctr harder toward epoch/decay (extra decrement).
    // Mirrors Tage.hpp's faralloc-based uctr adaptation.
    u64 FARALLOC_DIST = 0>
struct TageAhead : predictor {

  // ======== Derived Constants ========

  static_assert(!(DECAY_ENABLE && EPOCH_ENABLE),
                "DECAY_ENABLE and EPOCH_ENABLE are mutually exclusive");
  static_assert(NUM_PATHS == 1 ||
                    (USE_SEC_TAG && NUM_PATHS == (u64(1) << SEC_TAG_BITS)),
                "NUM_PATHS must be 1 or 2^SEC_TAG_BITS with USE_SEC_TAG");
  static_assert(NUM_PATHS <= 4, "NUM_PATHS > 4 not yet supported");

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
  // Fallback hysteresis: tracks agreement between fb and TAGE.
  // hyst=1 → agree, hyst=0 → disagree (weak → eligible for reconciliation).
  // Only accessed when FB_RECONCILE=true; zero cost otherwise.
  hcm::ram<val<1>, FB_CAPACITY> fb_hyst{"fb_hyst"};
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
  // Piped fb hyst for reconciliation (only accessed when FB_RECONCILE=true)
  reg<1> prefetch_fb_hyst;
  reg<1> current_fb_hyst;

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
  reg<PRED_BITS>
      train_pred[NT]; // per-table pred for per-table training signals
  reg<std::max(u64(1), HYST_WIDTH)> train_hyst[NT];
  reg<U_WIDTH> train_u[NT];
  reg<PRED_BITS> train_fb;
  reg<FB_IDX_BITS> train_fb_idx;
  reg<1> train_fb_hyst; // piped fb hyst for reconciliation
  reg<ALLOC_PC_BITS> train_pc;
  reg<MAX_TAG_WIDTH> train_ctag[NT]; // piped computed tag for allocation

  // Piped resolution values from previous update_cycle (block A's resolution)
  reg<MATCH_BITS> train_match1;
  reg<PRED_BITS> train_provider_pred;
  reg<1>
      train_provider_weak; // newly-allocated weakness (hyst==0 & u==0) — meta
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

  // Allocation LFSR (8-bit, hardware randomness for target policy + action
  // gating)
  reg<8> alloc_lfsr;

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
  reg<1> dbg_full_hits;            // after per-table hit computation
  reg<1> dbg_fb_pred;              // after fallback read
  reg<1> dbg_match;                // after concat into match bitmask
  reg<1> dbg_match1;               // after one_hot (provider)
  reg<1> dbg_match2;               // after one_hot (alt)
  reg<1> dbg_provider_pred;        // after replicate-mask-fold
  reg<1> dbg_alt_pred;             // after replicate-mask-fold
  reg<1> dbg_provider_weak;        // after weakness check
  reg<1> dbg_has_alt;              // after has_alt mask
  reg<1> dbg_meta_use_alt;         // after meta pipeline read
  reg<1> dbg_use_alt;              // after final use_alt AND
  reg<1> dbg_final_pred;           // after select mux
  reg<1> dbg_mispredict;           // framework's is_mispredict
  reg<1> dbg_true_block;           // true_block after computation
  reg<1> dbg_inst_pc;              // inst_pc arrival in predict1
  reg<1> dbg_fold_idx;             // worst fold_idx.get() across tables
  reg<1> dbg_fold_tag;             // worst fold_tag.get() across tables
  reg<1> dbg_p1_return;            // predict1 return value timing
  reg<1> dbg_hist_input;           // hist_input timing in update_cycle
  reg<1> dbg_gh_fanout;            // gh after fanout in update_cycle
  reg<1> dbg_fold_compute;         // compute_update result timing
  reg<1> dbg_fold_apply;           // after apply_update (fold write)
  reg<1> dbg_fold_early_write;     // fold write in num_branch==0 path
  reg<1> dbg_fold_read_in_compute; // folded reg read inside compute_update
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
  reg<1> dbg_fb_write;       // fb_ctr.write gate
  reg<1> dbg_meta_write;     // meta_ctr.write gate
  reg<1> dbg_pred_write[NT]; // pred_ram.write gate per table
  reg<1> dbg_hyst_write[NT]; // hyst_ram.write gate per table
  reg<1> dbg_tag_write[NT];  // tag_ram.write gate per table
  reg<1> dbg_u_write[NT];    // u_ram.write gate per table
  // Decay (per-table)
  reg<1> dbg_decay_fire[NT];   // decay_fire per table
  reg<1> dbg_decay_merged[NT]; // merged u value per table
  reg<1> dbg_epoch_fire;       // epoch trigger
  reg<1> dbg_next_pc;          // raw next_pc timing from block_end_info
  reg<1> dbg_curr_sec_tag;     // curr_sec_tag after hash
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
    mon.shadow_block_pc = static_cast<u64>(inst_pc);
#endif
    inst_pc.fanout(
        hard<2 * NT + 2>{}); // 2 reads per table (>>2, >>4) + fb (>>2)
                              // + prefetch_pc (>>2)

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
      t.fold_idx.fanout(hard<2>{}); // get() + compute_update
      t.fold_tag.fanout(hard<2>{}); // get() + compute_update
      auto fold_idx_val = t.fold_idx.get();
      auto idx = fold_idx_val.fo1() ^ val<t.IDX_BITS>{inst_pc >> 2};
      idx.fanout(hard<6>{}); // 5 RAM reads + prefetch_idx write
      auto fold_tag_val = t.fold_tag.get();
      auto computed_tag = fold_tag_val.fo1() ^ val<t.tag_width>{inst_pc >> 4};
      computed_tag.fanout(hard<2>{}); // tag comparison + prefetch_ctag write

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
      if constexpr (USE_GSHARE) {
        fb_fold.fanout(hard<2>{}); // get() + compute_update
        auto fb_fold_val = fb_fold.get();
        return val<FB_IDX_BITS>{inst_pc >> 2} ^ fb_fold_val.fo1();
      } else
        return val<FB_IDX_BITS>{inst_pc >> 2};
    }();
    prefetch_fb_idx = fb_idx;
    prefetch_fb = fb_ctr.read(fb_idx);
    if constexpr (FB_RECONCILE)
      prefetch_fb_hyst = fb_hyst.read(fb_idx);
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
    mon.record_branch_info(num_branch, static_cast<u64>(branch_pc),
                           static_cast<u64>(taken));
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
      train_pred[I] =
          current_pred[I]; // per-table pred for per-table wrong signals
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
    if constexpr (FB_RECONCILE)
      train_fb_hyst = current_fb_hyst;
    train_pc = current_pc;

    // Fanout on prefetch_* regs: read twice (shift + resolution chain bypass)
    static_loop<NT>([&]<u64 I>() {
      prefetch_tag_hit[I].fanout(hard<2>{}); // shift + full_hits
      prefetch_pred[I].fanout(hard<2>{});    // shift + table_preds
      prefetch_hyst[I].fanout(hard<2>{});    // shift + weak_mask
      prefetch_u[I].fanout(hard<2>{});       // shift + weak_mask
      if constexpr (USE_SEC_TAG)
        prefetch_sec[I].fanout(hard<2>{});   // shift + full_hits
    });
    prefetch_fb.fanout(hard<2>{}); // shift + table_preds

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
    if constexpr (FB_RECONCILE)
      current_fb_hyst = prefetch_fb_hyst;
    current_pc = prefetch_pc;

#ifdef TAGE_MONITOR
    mon.shift_block_pc();
#endif

    // Precompute secondary tag for next block
    // Use a local val to bypass the reg's transparent-latch penalty (119ps).
    // The reg is only needed for cross-cycle reads; same-cycle resolution
    // uses sec_tag_now directly.
    if constexpr (!USE_SEC_TAG)
      block_end_info.next_pc.fanout(hard<2>{}); // meta_idx + hist path_bits
    else
      // +1 for sec_tag_alloc duplicate hash (split to reduce crit-path fanout)
      block_end_info.next_pc.fanout(hard<4 + (NUM_PATHS > 1 ? 1 : 0)>{});
    auto sec_tag_now = [&]() {
      if constexpr (USE_SEC_TAG)
        return SecTagHashFn::template apply<SEC_TAG_BITS>(
            block_end_info.next_pc);
      else
        return hard<0>{};
    }();
    // Separate alloc hash: same computation, independent val with its own
    // fanout tree. Keeps alloc's NT writes off sec_tag_now's critical-path
    // buffer tree (hard<NT+1>=3 FO2 vs old hard<2*NT+1>=4 FO2).
    auto sec_tag_alloc = [&]() {
      if constexpr (USE_SEC_TAG)
        return SecTagHashFn::template apply<SEC_TAG_BITS>(
            block_end_info.next_pc);
      else
        return hard<0>{};
    }();
    if constexpr (USE_SEC_TAG) {
      // Critical-path fanout: reg write (1) + full_hits comparisons (NT) = NT+1
      sec_tag_now.fanout(hard<NT + 1>{});
      curr_sec_tag = sec_tag_now; // store for next-cycle use
      // Alloc-path fanout: NT sec_ram writes (off critical path)
      sec_tag_alloc.fanout(hard<NT>{});
#ifdef TIMING_DEBUG
      dbg_next_pc = val<1>{block_end_info.next_pc};
      dbg_curr_sec_tag = val<1>{sec_tag_now};
#endif
    }

    // ================================================================
    // Meta pipeline: shift FIRST, then write [0].
    // Fixes stale-value bug where meta_pipe[META_PIPE-1] reads the
    // new RAM value instead of the properly delayed old value.
    // ================================================================
    for (u64 i = META_PIPE - 1; i > 0; i--) {
      meta_pipe[i] = meta_pipe[i - 1];
      meta_idx_pipe[i] = meta_idx_pipe[i - 1];
    }
    {
      auto meta_idx = val<META_IDX_BITS>{block_end_info.next_pc >> 2};
      meta_pipe[0] = meta_ctr.read(meta_idx);
      meta_idx_pipe[0] = meta_idx;
    }
    val<1> meta_use_alt =
        val<META_WIDTH, i64>{meta_pipe[META_PIPE - 1]} >= hard<0>{};
#ifdef TIMING_DEBUG
    dbg_meta_use_alt = meta_use_alt;
#endif

    // ================================================================
    // Provider / altpred resolution via bitmask + one_hot
    //
    // NUM_PATHS > 1 (paper: Cai/Deshmukh/Patt ISCA'25 Sec 4.2):
    //   Run NUM_PATHS independent parallel chains. Each chain compares
    //   stored sec_tag against a compile-time constant (0, 1, ...).
    //   Each chain independently finds its own provider and alt.
    //   Final select with curr_sec_tag is a regular MUX OFF crit path.
    //
    // NUM_PATHS == 1: single chain with runtime curr_sec_tag comparison
    //   in full_hits (original behavior, curr_sec_tag on crit path).
    // ================================================================

    // Save old prediction (block B) before scatter overwrites with B+1.
    branch_dir.fanout(
        hard<4>{}); // badpred + true_block + hist_input + actual_dir
    arr<val<1>, N> old_pred = [&](u64 i) -> val<1> { return val<1>{pred[i]}; };

#ifdef TAGE_MONITOR
    for (u64 i = 0; i < NT; i++) {
      if constexpr (USE_SEC_TAG) {
        mon.record_tag_lookup(
            i, static_cast<u64>(current_tag_hit[i]),
            static_cast<u64>(val<SEC_TAG_BITS>{current_sec[i]} ==
                             val<SEC_TAG_BITS>{sec_tag_now}));
      } else {
        mon.record_tag_lookup(i, static_cast<u64>(current_tag_hit[i]), false);
      }
      mon.record_collision_check(
          i, static_cast<u64>(val<MAX_IDX_BITS>{current_idx[i]}),
          static_cast<u64>(current_tag_hit[i]));
    }
#endif

    // -- Resolve one chain: full_hits → match → provider/alt → final pred --
    // Critical path: full_hits → match → one_hot → pp/ap → ua → select
    // hyst_weak (phw) is NOT computed here — it's training-only, so it's
    // computed from piped regs in the training section to stay off crit path.

    // Precompute per-table weakness — read prefetch_* directly to bypass
    // the current_* transparent-latch penalty (saves one reg hop on crit path).
    val<NT> weak_mask = [&]() {
      constexpr u64 HW = std::max(u64(1), HYST_WIDTH);
      arr<val<1>, NT> w = [&](u64 i) -> val<1> {
        return val<1>{val<HW>{prefetch_hyst[i]} == hard<0>{}} &
               val<1>{val<U_WIDTH>{prefetch_u[i]} == hard<0>{}};
      };
      return w.concat();
    }();

    auto resolve_chain = [&](arr<val<1>, NT> full_hits) {
      val<MATCH_BITS> match = concat(val<1>{1}, full_hits.concat());
      match.fanout(hard<2>{}); // one_hot + ha

      // has_alt: true iff any TAGE table matched. If so, provider is a
      // table and fallback (always bit 0) is below it as alt.
      // Depends only on match — no match1 dependency.
      val<1> ha = (match >> 1) != hard<0>{};

      val<MATCH_BITS> match1 = match.one_hot();
      // make_array + XOR(remainder) + pw + train_match1 save
      match1.fanout(hard<4>{});
      val<MATCH_BITS> remainder = match ^ match1;
      val<MATCH_BITS> match2 = remainder.fo1().one_hot();
      match2.fanout(hard<2>{}); // make_array (alt pred extraction)

      arr<val<PRED_BITS>, NT + 1> table_preds = [&](u64 i) -> val<PRED_BITS> {
        if (i < NT)
          return val<PRED_BITS>{prefetch_pred[i]};
        return val<PRED_BITS>{prefetch_fb};
      };
      table_preds.fanout(hard<2>{}); // pmask + amask

      arr<val<1>, NT + 1> m1_bits = match1.make_array(val<1>{});
      arr<val<1>, NT + 1> m2_bits = match2.make_array(val<1>{});

      arr<val<PRED_BITS>, NT + 1> pmask = [&](u64 i) {
        return m1_bits[i].fo1().replicate(hard<PRED_BITS>{}).concat() &
               table_preds[i];
      };
      arr<val<PRED_BITS>, NT + 1> amask = [&](u64 i) {
        return m2_bits[i].fo1().replicate(hard<PRED_BITS>{}).concat() &
               table_preds[i];
      };
      val<PRED_BITS> pp = pmask.fold_or();
      val<PRED_BITS> ap = amask.fold_or();
      // pp: make_array + ad XOR + train_provider_pred save
      pp.fanout(hard<3>{});
      // ap: make_array + ad XOR
      ap.fanout(hard<2>{});

      // Provider weakness: AND match1's table bits with precomputed
      // weak_mask, then check nonzero. Single gate after match1 instead
      // of fold_or over NT+1 terms.
      val<1> pw = (val<NT>{match1 >> 1} & weak_mask) != hard<0>{};
      pw.fanout(hard<2>{}); // ua computation + train_provider_weak save
      val<1> ua = pw & meta_use_alt.fo1() & ha.fo1();
      // ua: N scatter reads
      ua.fanout(hard<N>{});

      // Final select is NOT done here — it's fused with the per-branch
      // scatter outside, using 1-bit selects to avoid PRED_BITS-wide
      // replication of ua (saves ~45ps ctrl-to-out delay).
      auto pp_xor_ap = pp ^ ap;
      val<1> ad = pp_xor_ap.fo1() != hard<0>{};

      return std::tuple{std::move(match1), std::move(match2), std::move(pp),
                        std::move(ap),     std::move(pw),     std::move(ad),
                        std::move(ua)};
    };

    // Run resolution chain(s) and select active result.
    // Returns {match1, match2, pp, ap, pw, ad, ua} — final select is
    // deferred to per-branch 1-bit scatter below.
    auto [match1, match2, provider_pred, alt_pred, provider_weak, altdiff,
          use_alt] = [&]() {
      if constexpr (NUM_PATHS > 1) {
        // Multi-chain reads: read prefetch_* directly to bypass current_*
        // transparent-latch penalty (predict1 wrote prefetch, shift copies
        // to current — reading prefetch skips the second reg hop).

        // Chain 0: stored sec_tag == 0 (compile-time constant)
        arr<val<1>, NT> fh0 = [&](u64 i) {
          return val<1>{prefetch_tag_hit[i]} &
                 (val<SEC_TAG_BITS>{prefetch_sec[i]} == hard<0>{});
        };
        auto [m1_0, m2_0, pp0, ap0, pw0, ad0, ua0] = resolve_chain(fh0);

        // Chain 1: stored sec_tag == 1 (compile-time constant)
        arr<val<1>, NT> fh1 = [&](u64 i) {
          return val<1>{prefetch_tag_hit[i]} &
                 (val<SEC_TAG_BITS>{prefetch_sec[i]} == hard<1>{});
        };
        auto [m1_1, m2_1, pp1, ap1, pw1, ad1, ua1] = resolve_chain(fh1);

        // Fields muxed through NP select: m1, m2, pp, ap, pw, ad, ua.
        static constexpr u64 MUX_FIELDS = 7;

        if constexpr (NUM_PATHS == 2) {
          // 2-to-1 mux: derive sel directly from next_pc (bypass reg chain)
          val<1> sec_sel = val<1>{SecTagHashFn::template apply<SEC_TAG_BITS>(
              block_end_info.next_pc)};
          sec_sel.fanout(hard<MUX_FIELDS>{});
          return std::tuple{
              select(sec_sel, m1_1, m1_0), select(sec_sel, m2_1, m2_0),
              select(sec_sel, pp1, pp0),   select(sec_sel, ap1, ap0),
              select(sec_sel, pw1, pw0),   select(sec_sel, ad1, ad0),
              select(sec_sel, ua1, ua0)};
        } else if constexpr (NUM_PATHS == 4) {
          // Chains 2-3
          arr<val<1>, NT> fh2 = [&](u64 i) {
            return val<1>{prefetch_tag_hit[i]} &
                   (val<SEC_TAG_BITS>{prefetch_sec[i]} == hard<2>{});
          };
          auto [m1_2, m2_2, pp2, ap2, pw2, ad2, ua2] = resolve_chain(fh2);

          arr<val<1>, NT> fh3 = [&](u64 i) {
            return val<1>{prefetch_tag_hit[i]} &
                   (val<SEC_TAG_BITS>{prefetch_sec[i]} == hard<3>{});
          };
          auto [m1_3, m2_3, pp3, ap3, pw3, ad3, ua3] = resolve_chain(fh3);

          // 4-to-1 mux tree: derive sel directly from next_pc (bypass reg)
          val<SEC_TAG_BITS> sec_idx =
              SecTagHashFn::template apply<SEC_TAG_BITS>(
                  block_end_info.next_pc);
          val<1> lo = sec_idx & hard<1>{};
          val<1> hi = (sec_idx >> 1) & hard<1>{};
          lo.fanout(hard<MUX_FIELDS * 2>{}); // each field needs lo for 2 pairs
          hi.fanout(hard<MUX_FIELDS>{});     // each field needs hi for final

          auto mux4 = [&](auto a, auto b, auto c, auto d) {
            auto ab = select(lo, b, a);
            auto cd = select(lo, d, c);
            return select(hi, cd, ab);
          };
          return std::tuple{
              mux4(m1_0, m1_1, m1_2, m1_3), mux4(m2_0, m2_1, m2_2, m2_3),
              mux4(pp0, pp1, pp2, pp3),     mux4(ap0, ap1, ap2, ap3),
              mux4(pw0, pw1, pw2, pw3),     mux4(ad0, ad1, ad2, ad3),
              mux4(ua0, ua1, ua2, ua3)};
        }
      } else {
        // Single chain: original behavior
        // sec_tag_now fanout already declared above (covers all reads)
        arr<val<1>, NT> full_hits = [&](u64 i) {
          if constexpr (USE_SEC_TAG) {
            return val<1>{prefetch_tag_hit[i]} &
                   (val<SEC_TAG_BITS>{prefetch_sec[i]} ==
                    val<SEC_TAG_BITS>{sec_tag_now});
          } else {
            return val<1>{prefetch_tag_hit[i]};
          }
        };
        return resolve_chain(full_hits);
      }
    }();

    // Per-branch 1-bit scatter: split pp/ap into per-bit arrays to reduce
    // fanout on the wide values (fanout(3) + N×fo1 vs fanout(N+2) on each).
    arr<val<1>, PRED_BITS> pp_bits = provider_pred.make_array(val<1>{});
    arr<val<1>, PRED_BITS> ap_bits = alt_pred.make_array(val<1>{});
    static_loop<N>([&]<u64 I>() {
      pred[I] = select(use_alt, ap_bits[I].fo1(), pp_bits[I].fo1());
    });

#ifdef TIMING_DEBUG
    dbg_full_hits = val<1>{pred[0]}; // timing proxy
    dbg_fb_pred = val<1>{pred[0]};   // timing proxy
    dbg_match = val<1>{match1};      // timing proxy
    dbg_match1 = val<1>{match1};
    dbg_match2 = val<1>{match2};
    dbg_provider_pred = val<1>{provider_pred};
    dbg_alt_pred = val<1>{alt_pred};
    dbg_provider_weak = provider_weak;
    dbg_has_alt = provider_weak; // timing proxy
    dbg_use_alt = use_alt;
    dbg_final_pred = val<1>{pred[0]};
#endif

#ifdef TAGE_MONITOR
    for (u64 r = 0; r < num_branch; r++) {
      bool meta_overrode =
          static_cast<u64>(provider_weak) &&
          ((static_cast<u64>(match2) & ((1ULL << NT) - 1)) != 0);
      bool meta_chose = static_cast<u64>(use_alt);
      bool pred_taken = static_cast<u64>(pred[r]);
      mon.record_prediction(r, static_cast<u64>(match1),
                            static_cast<u64>(match2), meta_overrode, meta_chose,
                            pred_taken);
    }
    {
      u64 m1v = static_cast<u64>(match1);
      u64 prov = decltype(mon)::decode_provider(m1v);
      if (prov < NT) {
        u64 prov_index = static_cast<u64>(val<MAX_IDX_BITS>{current_idx[prov]});
        u64 nc = 0;
        for (u64 r = 0; r < num_branch; r++)
          if (static_cast<u64>(pred[r]) == static_cast<u64>(branch_dir[r]))
            nc++;
        mon.record_provider_entry(prov, prov_index, num_branch, nc);
      }
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
    val<1> t_pw = train_provider_weak; // newly-alloc (hyst==0 & u==0) — meta
    val<1> t_ad = train_altdiff;

    // Hyst-only weakness: computed from piped regs (not resolution chain)
    // to keep it off the resolution critical path. Gates pred/hyst counter
    // updates — a useful entry (u>0) with weak hyst should still flip.
    // Tage.hpp equivalent: g_weak = primary & badpred & (hyst==0).
    constexpr u64 HW_T = std::max(u64(1), HYST_WIDTH);
    arr<val<1>, NT + 1> t_hyst_weak_arr = [&](u64 i) -> val<1> {
      if (i < NT)
        return t_m1[i] & val<1>{val<HW_T>{train_hyst[i]} == hard<0>{}};
      return val<1>{0};
    };
    val<1> t_phw = t_hyst_weak_arr.fold_or();

    // Save current resolution → train regs (for NEXT cycle's training)
    train_match1 = match1;
    train_provider_pred = provider_pred;
    train_provider_weak = provider_weak;
    train_altdiff = altdiff;
#ifdef TIMING_DEBUG
    dbg_altdiff = altdiff;
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
#ifdef TIMING_DEBUG
      dbg_fold_early_write = val<1>{std::get<0>(tables).fold_idx.get()};
#endif
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
        hard<5 + (AllocCfg::ALLOC_TRIGGER == AllocTrigger::MISPREDICT
                      ? 1
                      : 0)>{}); // extra_cycle + fb + true_block + dbg + acc_ctr
                                // + (alloc if MISPREDICT)
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
    // any_provider_wrong(1) + per-table table_wrong in loop(NT) + fb_changed(1)
    actual_dir.fanout(hard<NT + 2>{});
    val<1> any_provider_wrong = (t_pp ^ actual_dir) != hard<0>{};
    // Now only consumed by: meta update(1) + alloc trigger if TAGE_WRONG(1)
    any_provider_wrong.fanout(
        hard<1 + (AllocCfg::ALLOC_TRIGGER == AllocTrigger::TAGE_WRONG
                      ? 1
                      : 1)>{}); // patch for now to avoid
                                // diddling fo1()
#ifdef TIMING_DEBUG
    dbg_actual_dir = val<1>{actual_dir};
    dbg_any_prov_wrong = any_provider_wrong;
#endif
    t_pw.fanout(hard<2>{});       // meta gate + meta update direction
    t_phw.fanout(hard<NT + 1>{}); // per-table pred/hyst update gate
    t_m1.fanout(hard<6>{});
    t_ad.fanout(hard<NT + 1>{});

    // ---- Step 3: Allocation ----

    // 3a. Allocation trigger
    val<1> alloc_trigger = [&]() -> val<1> {
      if constexpr (AllocCfg::ALLOC_TRIGGER == AllocTrigger::MISPREDICT)
        return mispredict;
      else if constexpr (AllocCfg::ALLOC_TRIGGER == AllocTrigger::TAGE_WRONG)
        return any_provider_wrong;
      else {
        static_assert(AllocCfg::ALLOC_TRIGGER == AllocTrigger::ALWAYS);
        return val<1>{1};
      }
    }();
    val<NT> triggermask = alloc_trigger.replicate(hard<NT>{}).concat();

    // 3b. Allocation action (probabilistic gating)
    val<8> alloc_rng = val<8>{alloc_lfsr};
    val<NT> gated_triggermask = [&]() -> val<NT> {
      if constexpr (AllocCfg::ALLOC_ACTION == AllocAction::STANDARD) {
        return triggermask;
      } else if constexpr (AllocCfg::ALLOC_ACTION == AllocAction::FILTERED) {
        static_assert(ACC_WIDTH > 0, "FILTERED requires ACC_WIDTH > 0");
        val<1> allow = (val<ACC_WIDTH>{acc_ctr} >= val<ACC_WIDTH>{alloc_rng});
        return triggermask & allow.replicate(hard<NT>{}).concat();
      } else {
        static_assert(AllocCfg::ALLOC_ACTION == AllocAction::THROTTLED);
        static_assert(ALLOC_WIDTH > 0, "THROTTLED requires ALLOC_WIDTH > 0");
        val<1> allow =
            (val<ALLOC_WIDTH>{alloc_ctr} >= val<ALLOC_WIDTH>{alloc_rng});
        return triggermask & allow.replicate(hard<NT>{}).concat();
      }
    }();

    // 3c. Candidate mask: tables above provider with u=0
    val<NT> alloc_base = val<NT>{t_match1 - 1};
    arr<val<1>, NT> u_zero = [&](u64 i) -> val<1> {
      return val<U_WIDTH>{train_u[i]} == hard<0>{};
    };
    val<NT> notumask = u_zero.concat();
    notumask.fanout(hard<2>{});
    val<NT> postmask = alloc_base & gated_triggermask;
    postmask.fanout(hard<2>{});
    val<NT> candallocmask = [&]() {
      val<NT> base = postmask & notumask;
      if constexpr (USE_SEC_TAG) {
        // Allocation promotion (Cai/Deshmukh/Patt ISCA'25 Sec 4.3):
        // Skip entries with same primary tag but different sec_tag (siblings).
        // Promotes allocation to the next higher table.
        arr<val<1>, NT> not_sibling = [&](u64 i) -> val<1> {
          return ~(val<1>{train_tag_hit[i]} & ~val<1>{train_sec_hit[i]});
        };
        return base & not_sibling.concat();
      } else {
        return base;
      }
    }();
    candallocmask.fanout(hard<2>{});

    // 3d. Target policy (may skip closest candidates)
    val<NT> collamask = AllocCfg::TARGET_POLICY::template apply<NT>(
        candallocmask.reverse(), val<ALLOC_WIDTH>{alloc_ctr},
        val<ACC_WIDTH>{acc_ctr}, alloc_rng);

    // 3e. Final allocation decision (one-hot or two-hot)
    arr<val<1>, NT> allocate = [&]() -> arr<val<1>, NT> {
      if constexpr (AllocCfg::MAX_ALLOC >= 2) {
        collamask.fanout(hard<3>{});
        val<NT> pick1 = collamask.one_hot();
        pick1.fanout(hard<3>{});
        val<NT> pick2 = [&]() -> val<NT> {
          val<NT> basic2 = (collamask ^ pick1).one_hot();
          if constexpr (AllocCfg::NON_CONSECUTIVE) {
            val<NT> neighbors = (pick1 << 1) | (pick1 >> 1);
            val<NT> nc_mask = (collamask ^ pick1) & ~neighbors;
            val<NT> nc_pick = nc_mask.reverse().one_hot();
            return select(nc_mask != hard<0>{}, nc_pick, basic2);
          } else {
            return basic2;
          }
        }();
        return (pick1 | pick2).reverse().make_array(val<1>{});
      } else {
        val<NT> alloc_target_rev = collamask.one_hot();
        return alloc_target_rev.reverse().make_array(val<1>{});
      }
    }();
    allocate.fanout(hard<6>{});
    val<NT> alloc_target = [&]() {
      arr<val<1>, NT> a = allocate;
      return a.concat();
    }();
    val<1> noalloc = (candallocmask == hard<0>{});
    val<NT> uclearmask = postmask & noalloc.replicate(hard<NT>{}).concat();
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
      u64 prov_idx =
          TAMonitor<NT, N, MAX_TABLE_SIZE, USE_GSHARE,
                    FB_CAPACITY>::decode_provider(static_cast<u64>(t_match1));
      if (static_cast<u64>(alloc_trigger)) {
        if (static_cast<u64>(postmask) == 0)
          mon.record_alloc_blocked(); // fallback was provider, no candidates
        mon.record_alloc_cascade(prov_idx, at);
        if constexpr (USE_SEC_TAG) {
          u64 sibling = 0;
          for (u64 i = 0; i < NT; i++)
            if (static_cast<u64>(train_tag_hit[i]) &&
                !static_cast<u64>(train_sec_hit[i]))
              sibling |= (u64(1) << i);
          sibling &= static_cast<u64>(postmask); // only above provider
          if (sibling)
            mon.record_alloc_sibling_skip(sibling);
        }
      }
    }
#endif

    // ---- Step 4: Fallback update ----
    // 4a. Direct update: mispredict + fallback is provider → write actual_dir.
    val<1> fb_changed = actual_dir != val<PRED_BITS>{train_fb};
    val<1> fb_gate = do_train & t_m1[NT] & mispredict & fb_changed;
    execute_if(fb_gate, [&]() {
      fb_ctr.write(val<FB_IDX_BITS>{train_fb_idx}, actual_dir);
    });

    // 4b. Reconciliation (Tage.hpp P1/P2 pattern): when TAGE and fb disagree
    // persistently (fb_hyst weak), overwrite fb with TAGE's prediction.
    // Also update fb_hyst every cycle to track agreement.
    if constexpr (FB_RECONCILE) {
      val<1> fb_not_provider = ~t_m1[NT];
      val<1> fb_tage_disagree = (val<PRED_BITS>{train_fb} ^ t_pp) != hard<0>{};
      val<1> fb_hyst_weak = ~val<1>{train_fb_hyst};

      // Overwrite fb pred when hyst weak and they disagree (non-provider only;
      // provider case handled above by fb_gate)
      val<1> reconcile_gate =
          do_train & fb_not_provider & fb_hyst_weak & fb_tage_disagree;
      execute_if(reconcile_gate,
                 [&]() { fb_ctr.write(val<FB_IDX_BITS>{train_fb_idx}, t_pp); });

      // Update fb_hyst: 1 if agree, 0 if disagree — tracks recent consensus.
      execute_if(do_train, [&]() {
        fb_hyst.write(val<FB_IDX_BITS>{train_fb_idx}, ~fb_tage_disagree);
      });
    }

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
    // Uses per-table wrong signal (this table's own prediction vs actual)
    // instead of blanket any_provider_wrong. For the provider table these
    // are identical; for non-provider tables this enables future per-table
    // update policies (e.g. partial update).
    train_pc.fanout(hard<NT>{});
    static_loop<NT>([&]<u64 I>() {
      auto &t = std::get<I>(tables);
      val<1> do_alloc = allocate[I];

      // Per-table wrong: does THIS table's stored pred disagree with actual?
      val<1> table_wrong =
          (val<PRED_BITS>{train_pred[I]} ^ actual_dir) != hard<0>{};

      // pred_ram: alloc writes actual_dir, update writes actual_dir → same
      // data. Gated on hyst-only weakness (t_phw), not newly-alloc (t_pw). A
      // useful entry (u>0) with weak hyst that's wrong should still flip.
      val<1> do_pred_update = t_m1[I] & t_phw & table_wrong;
      execute_if(do_train & (do_alloc | do_pred_update), [&]() {
        t.pred_ram.write(val<t.IDX_BITS>{train_idx[I]}, actual_dir, hard<0>{});
      });

      // hyst_ram: alloc writes 0, update writes new_hyst → mux on do_alloc.
      // Direction based on this table's own prediction accuracy.
      constexpr u64 HW = std::max(u64(1), HYST_WIDTH);
      auto old_hyst = val<HW>{train_hyst[I]};
      auto new_hyst = ta_update_ctr(old_hyst, ~table_wrong);
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
                          val<SEC_TAG_BITS>{sec_tag_alloc});
      });
#ifdef TAGE_MONITOR
      if (static_cast<u64>(do_train & do_alloc)) {
        u64 tidx = static_cast<u64>(val<t.IDX_BITS>{train_idx[I]});
        mon.record_tage_write(I, tidx);
        mon.record_entry_alloc_diag(I, tidx);
      }
#endif

      // u_ram: combined provider update + allocation + uclear + decay.
      // Uses per-table wrong signal: u set when THIS table predicted correctly.
      val<U_WIDTH> base_newu = val<U_WIDTH>{~table_wrong} &
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

      // Alloc pressure counter: increment when no alloc slot found,
      // with optional extra decrement for far allocations.
      val<1> any_alloc = alloc_target != hard<0>{};
      auto new_alloc = [&]() {
        if constexpr (FARALLOC_DIST > 0) {
          // Far allocation: alloc target is >= FARALLOC_DIST tables from
          // provider. Shift provider mask down by FARALLOC_DIST; if alloc
          // target is still the one_hot winner, it's far → extra decrement.
          val<NT> allocmask1 = alloc_target;
          val<NT> provider_bits = val<NT>{t_match1 >> 1};
          val<1> faralloc =
              (((provider_bits >> FARALLOC_DIST) | allocmask1).one_hot() ^
               allocmask1) == hard<0>{};
          // Two-step update: first normal, then extra decrement if far
          auto step1 =
              ta_update_ctr(val<ALLOC_WIDTH>{alloc_ctr}, ~any_alloc);
          return ta_update_ctr(step1, ~(faralloc & any_alloc));
        } else {
          return ta_update_ctr(val<ALLOC_WIDTH>{alloc_ctr}, ~any_alloc);
        }
      }();

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

    // ---- Allocation LFSR tick (8-bit, polynomial x^8+x^6+x^5+x^4+1) ----
    {
      val<8> old = val<8>{alloc_lfsr};
      val<1> fb = old & hard<1>{};
      val<8> shifted = val<8>{old >> 1};
      val<8> tap = val<8>{u64(1) << 7};
      alloc_lfsr = shifted ^ (tap & fb.replicate(hard<8>{}).concat());
    }

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
        hard<NT * 2 + 2 + (USE_GSHARE ? 1 : 0)>{}); // NT fold_idx + NT fold_tag
                                                    // apply_update muxes +
                                                    // gh.update + monitor +
                                                    // (fb_fold if gshare)

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
#ifdef TIMING_DEBUG
      if constexpr (I == 0) {
        dbg_fold_read_in_compute = val<1>{t.fold_idx.get()};
      }
#endif
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
