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

template <u64 N = 8, u64 SIZE = 1024, u64 TAG = 10, u64 MINH = 10,
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
    u64 CTR_WIDTH = 1,       // prediction counter width per lane
    u64 HYST_WIDTH = 3,      // hysteresis width (separate from ctr)
    u64 U_WIDTH = 2,         // usefulness counter width
    u64 FB_CAPACITY = 8192 * 2, // fallback table size (bimodal or gshare)
    bool USE_GSHARE = false, // use gshare base (PC^history) vs bimodal (PC)
    u64 GS_HIST = 6,         // gshare history length (only when USE_GSHARE)
    u64 META_WIDTH = 5,      // meta counter width (provider vs alt)
    u64 META_CAPACITY = 256, // meta table entries
    u64 META_PIPE = 2,       // meta pipeline depth
    u64 LINEINST = 1024,     // line size in instructions
    bool SHARED_HYS = true,  // shared hyst: 2 entries share 1 counter
    HistUpdate HIST_MODE =
        HistUpdate::PATH, // what goes into history: PATH, DIR, or BOTH
    // ---- Allocation policy ----
    typename AllocCfg = TAAlloc2,
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
    u64 EPOCH_CTR_WIDTH =
        16, // epoch counter width (for interval-based triggers)
    bool EPOCH_RESET_ACC = false,  // reset acc_ctr on epoch fire
    bool EPOCH_RESET_ALLOC = true, // reset alloc_ctr on epoch fire
    // ---- Replacement policy (Technique 4) ----
    typename ReplacePolicyFn = ReplaceUZero,
    // ---- Alt bank promotion (Technique 3) ----
    typename AltPromoteFn = AltPromoteOff,
    // ---- DIP-like allocation (Technique 6) ----
    u64 DIP_PROB_256 =
        0, // probability of setting u>0 on alloc (0=off, 256=always)
    u64 DIP_INIT_U = 1, // initial u value when DIP fires
    // ---- Provider u-bit update policy ----
    UProvUpdate U_PROV_UPDATE = UProvUpdate::INC_DEC,
    // ---- Fallback trains toward P2 (gshare only) ----
    bool FB_TRAIN_P2 = false, // when TAGE is provider and gshare disagrees,
                              // write TAGE pred into gshare
    // ---- Far-allocation epoch pressure ----
    u64 FARALLOC_DIST =
        3, // 0=off; when >0, alloc_ctr uses faralloc logic (Tage.hpp style)
    // ---- Per-branch hyst/u banking ----
    u64 HYST_BANKS = 2, // hyst banks per table (1=shared, N=per-branch)
    u64 HYST_BANK_BIT =
        0,             // starting bit of branch rank for hyst bank selection
    u64 U_BANKS = 4,   // u banks per table (1=shared, N=per-branch)
    u64 U_BANK_BIT = 0 // starting bit of branch rank for u bank selection
    >
struct TageAhead : predictor {

  // ======== Derived Constants ========

  static_assert(!(DECAY_ENABLE && EPOCH_ENABLE),
                "DECAY_ENABLE and EPOCH_ENABLE are mutually exclusive");
  static_assert(CTR_WIDTH == 1);
  static_assert(!FB_TRAIN_P2 || USE_GSHARE,
                "FB_TRAIN_P2 requires USE_GSHARE=true");

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

  // Per-branch banking helpers
  static_assert(HYST_BANKS >= 1 && HYST_BANKS <= N &&
                std::has_single_bit(HYST_BANKS));
  static_assert(U_BANKS >= 1 && U_BANKS <= N && std::has_single_bit(U_BANKS));
  static constexpr u64 HYST_BANK_BITS = std::bit_width(HYST_BANKS) - 1;
  static constexpr u64 U_BANK_BITS = std::bit_width(U_BANKS) - 1;
  // Compile-time bank index for branch b
  static constexpr u64 hyst_bank_of(u64 b) {
    return (b >> HYST_BANK_BIT) & (HYST_BANKS - 1);
  }
  static constexpr u64 u_bank_of(u64 b) {
    return (b >> U_BANK_BIT) & (U_BANKS - 1);
  }
  // Compile-time bitmask of branches in a given hyst/u bank
  static constexpr u64 hyst_bank_mask(u64 bank) {
    u64 mask = 0;
    for (u64 b = 0; b < N; b++)
      if (hyst_bank_of(b) == bank)
        mask |= (1ULL << b);
    return mask;
  }
  static constexpr u64 u_bank_mask(u64 bank) {
    u64 mask = 0;
    for (u64 b = 0; b < N; b++)
      if (u_bank_of(b) == bank)
        mask |= (1ULL << b);
    return mask;
  }

  // ======================================================================
  // Storage
  // ======================================================================
  // ---- Table tuple (per-table tag width and table size) ----
  using Tables =
      typename TAMakeTableTuple<TableCfg, CTR_WIDTH, HYST_WIDTH, U_WIDTH,
                                SEC_TAG_BITS, N, SHARED_HYS, HYST_BANKS,
                                U_BANKS, std::make_index_sequence<NT>>::type;
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
  reg<std::max(u64(1), HYST_WIDTH)> prefetch_hyst[NT][HYST_BANKS];
  reg<U_WIDTH> prefetch_u[NT][U_BANKS];
  reg<MAX_TAG_WIDTH> prefetch_ctag[NT]; // computed tag for allocation piping

  // current_*: shifted from prefetch, used for prediction
  reg<MAX_TAG_WIDTH> current_tag[NT];
  reg<1> current_tag_hit[NT];
  reg<PRED_BITS> current_pred[NT];
  reg<SEC_TAG_BITS> current_sec[NT];
  reg<MAX_IDX_BITS> current_idx[NT];
  reg<std::max(u64(1), HYST_WIDTH)> current_hyst[NT][HYST_BANKS];
  reg<U_WIDTH> current_u[NT][U_BANKS];
  reg<MAX_TAG_WIDTH> current_ctag[NT]; // piped computed tag for allocation

  // train_*: saved from current_* before pipeline shift, used for training
  // (training is for block A, resolution is for block B)
  reg<MAX_IDX_BITS> train_idx[NT];
  reg<std::max(u64(1), HYST_WIDTH)> train_hyst[NT][HYST_BANKS];
  reg<U_WIDTH> train_u[NT][U_BANKS];
  reg<PRED_BITS> train_fb;
  reg<FB_IDX_BITS> train_fb_idx;
  reg<ALLOC_PC_BITS> train_pc;
  reg<MAX_TAG_WIDTH> train_ctag[NT]; // piped computed tag for allocation

  // Piped resolution values from previous update_cycle (block A's resolution)
  reg<MATCH_BITS> train_match1;
  reg<MATCH_BITS> train_match2; // alt match (for alt promotion)
  reg<PRED_BITS> train_provider_pred;
  reg<PRED_BITS> train_alt_pred; // alt prediction (for alt promotion)
  reg<1> train_provider_weak;
  reg<1> train_provider_weakconf; // hyst==0 only (for pred flip)
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
  reg<EPOCH_CTR_WIDTH> epoch_ctr;

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
  reg<1> dbg_full_hits;     // after per-table hit computation
  reg<1> dbg_fb_pred;       // after fallback read
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
#endif

#ifdef TAGE_MONITOR
  TAMonitor<NT, N, MAX_TABLE_SIZE, USE_GSHARE, FB_CAPACITY, HYST_BANKS, U_BANKS>
      mon;
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
      static_loop<HYST_BANKS>([&]<u64 HB>() {
        prefetch_hyst[I][HB] = t.hyst_ram[HB].read(val<t.HYST_IDX_BITS>{idx});
      });
      static_loop<U_BANKS>(
          [&]<u64 UB>() { prefetch_u[I][UB] = t.u_ram[UB].read(idx); });
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
      static_loop<HYST_BANKS>(
          [&]<u64 HB>() { train_hyst[I][HB] = current_hyst[I][HB]; });
      static_loop<U_BANKS>(
          [&]<u64 UB>() { train_u[I][UB] = current_u[I][UB]; });
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
      static_loop<HYST_BANKS>(
          [&]<u64 HB>() { current_hyst[I][HB] = prefetch_hyst[I][HB]; });
      static_loop<U_BANKS>(
          [&]<u64 UB>() { current_u[I][UB] = prefetch_u[I][UB]; });
      current_ctag[I] = prefetch_ctag[I];
    });
    current_fb = prefetch_fb;
    current_fb_idx = prefetch_fb_idx;
    current_pc = prefetch_pc;

#ifdef TAGE_MONITOR
    mon.shift_block_pc();
#endif

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
      bool th = static_cast<u64>(current_tag_hit[i]);
      bool sh = USE_SEC_TAG
                    ? static_cast<u64>(val<SEC_TAG_BITS>{current_sec[i]} ==
                                       val<SEC_TAG_BITS>{curr_sec_tag})
                    : false;
      mon.record_tag_lookup(i, th, sh);
      mon.record_table_pred(
          i, static_cast<u64>(val<PRED_BITS>{current_pred[i]}), th, sh);
      // Feature 7: tag collision check
      mon.record_collision_check(
          i, static_cast<u64>(val<MAX_IDX_BITS>{current_idx[i]}), th);
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

    // 7. Provider weakness (per-branch via banking):
    //    provider_weak = any branch has (hyst==0 AND u==0) → meta override
    //    provider_weakconf = any branch has hyst==0 → pred direction flip
    constexpr u64 HW_RES = std::max(u64(1), HYST_WIDTH);
    arr<val<1>, NT + 1> weak_arr = [&](u64 i) -> val<1> {
      if (i < NT) {
        if constexpr (HYST_BANKS == 1 && U_BANKS == 1) {
          return m1_bits[i] &
                 val<1>{val<HW_RES>{current_hyst[i][0]} == hard<0>{}} &
                 val<1>{val<U_WIDTH>{current_u[i][0]} == hard<0>{}};
        } else {
          // OR across branches: any branch weak in this table?
          arr<val<1>, N> per_br = [&](u64 b) -> val<1> {
            u64 hb = hyst_bank_of(b);
            u64 ub = u_bank_of(b);
            return val<1>{val<HW_RES>{current_hyst[i][hb]} == hard<0>{}} &
                   val<1>{val<U_WIDTH>{current_u[i][ub]} == hard<0>{}};
          };
          return m1_bits[i] & (per_br.concat() != hard<0>{});
        }
      }
      return val<1>{0};
    };
    val<1> provider_weak = weak_arr.fold_or();
    arr<val<1>, NT + 1> weakconf_arr = [&](u64 i) -> val<1> {
      if (i < NT) {
        if constexpr (HYST_BANKS == 1) {
          return m1_bits[i] &
                 val<1>{val<HW_RES>{current_hyst[i][0]} == hard<0>{}};
        } else {
          arr<val<1>, N> per_br = [&](u64 b) -> val<1> {
            u64 hb = hyst_bank_of(b);
            return val<1>{val<HW_RES>{current_hyst[i][hb]} == hard<0>{}};
          };
          return m1_bits[i] & (per_br.concat() != hard<0>{});
        }
      }
      return val<1>{0};
    };
    val<1> provider_weakconf = weakconf_arr.fold_or();

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
    // Feature 5: Record provider entry usage for lifetime tracking
    {
      u64 m1v = static_cast<u64>(match1);
      u64 prov = decltype(mon)::decode_provider(m1v);
      if (prov < NT) {
        u64 prov_index = static_cast<u64>(val<MAX_IDX_BITS>{current_idx[prov]});
        u64 fp = static_cast<u64>(final_pred);
        u64 nc = 0;
        for (u64 r = 0; r < num_branch; r++)
          if (((fp >> r) & 1) == static_cast<u64>(branch_dir[r]))
            nc++;
        mon.record_provider_entry(prov, prov_index, num_branch, nc);
      }
    }
    // Record num_branch for sec rejection analysis pipeline
    mon.record_table_pred_num_branch(num_branch);
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
    val<1> t_pwc = train_provider_weakconf;
    val<1> t_ad = train_altdiff;

    // Save current resolution → train regs (for NEXT cycle's training)
    constexpr bool ALT_PROMOTE_ACTIVE =
        !std::is_same_v<AltPromoteFn, AltPromoteOff>;
    alt_pred.fanout(hard<2 + (ALT_PROMOTE_ACTIVE ? 1 : 0)>{});
    val<1> altdiff = (provider_pred ^ alt_pred) != hard<0>{};
    train_match1 = match1;
    train_provider_pred = provider_pred;
    train_provider_weak = provider_weak;
    train_provider_weakconf = provider_weakconf;
    train_altdiff = altdiff;
    if constexpr (ALT_PROMOTE_ACTIVE) {
      train_match2 = match2;
      train_alt_pred = alt_pred;
    }
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
        hard<5 + (AllocCfg::ALLOC_TRIGGER == AllocTrigger::MISPREDICT ? 1 : 0) +
             (FB_TRAIN_P2 ? 1 : 0) +
             (FARALLOC_DIST > 0 ? 1 : 0)>{}); // extra_cycle + fb + true_block +
                                              // dbg + acc_ctr
                                              // + (alloc if MISPREDICT) + fb_p2
                                              // + faralloc
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
    // Sec tag rejection analysis: compare shadowed per-table preds vs actual
    if (static_cast<u64>(do_train)) {
      u64 actual_bits = 0;
      for (u64 r = 0; r < num_branch; r++)
        actual_bits |= (static_cast<u64>(branch_dir[r]) & 1) << r;
      mon.eval_sec_rejections(actual_bits);
    }
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
    actual_dir.fanout(hard<NT + 1 + (ALT_PROMOTE_ACTIVE ? NT : 0)>{});
    val<1> any_provider_wrong = (t_pp ^ actual_dir) != hard<0>{};
    any_provider_wrong.fanout(
        hard<3 * NT + 1 +
             (AllocCfg::ALLOC_TRIGGER == AllocTrigger::TAGE_WRONG ? 1 : 0) +
             (ALT_PROMOTE_ACTIVE ? NT : 0)>{});

#ifdef TIMING_DEBUG
    dbg_actual_dir = val<1>{actual_dir};
    dbg_any_prov_wrong = any_provider_wrong;
#endif
    t_pw.fanout(hard<2>{});   // meta gate only
    t_pwc.fanout(hard<NT>{}); // pred flip per table
    t_m1.fanout(hard<6 + (FB_TRAIN_P2 ? 1 : 0)>{});
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

    // 3b. Alt bank promotion (Technique 3) — computed here so alloc_rng is
    // available
    val<8> alloc_rng = val<8>{alloc_lfsr};

    // 3c. Allocation action (probabilistic gating)
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

    // 3c. Candidate mask: tables above provider with replaceable entries
    val<NT> alloc_base = val<NT>{t_match1 - 1};
    arr<val<1>, NT> replaceable = [&](u64 i) -> val<1> {
      // Use bank 0 for replacement check (heuristic — entry-level decision)
      return ReplacePolicyFn::template is_replaceable<
          U_WIDTH, std::max(u64(1), HYST_WIDTH)>(
          val<U_WIDTH>{train_u[i][0]},
          val<std::max(u64(1), HYST_WIDTH)>{train_hyst[i][0]},
          val<ALLOC_WIDTH>{alloc_ctr}, val<ACC_WIDTH>{acc_ctr});
    };
    val<NT> notumask = replaceable.concat();
    notumask.fanout(hard<2>{});
    val<NT> postmask = alloc_base & gated_triggermask;
    postmask.fanout(hard<2>{});
    val<NT> candallocmask = postmask & notumask;
    candallocmask.fanout(hard<2>{});

    // 3d. Target policy (may skip closest candidates)
    val<NT> collamask = AllocCfg::TARGET_POLICY::template apply<NT>(
        candallocmask.reverse(), val<ALLOC_WIDTH>{alloc_ctr},
        val<ACC_WIDTH>{acc_ctr}, alloc_rng);

    // 3e. Final allocation decision (up to MAX_ALLOC picks with prob decay)
    // Helper: recursive pick from candidate mask. Each successive pick
    // has probability 1/2^(K*DECAY_SHIFT). Returns OR of all picks.
    constexpr auto multi_pick = []<u64 K>(auto self, val<NT> cands,
                                          val<8> rng) -> val<NT> {
      constexpr u64 MA = AllocCfg::MAX_ALLOC;
      constexpr u64 DS = AllocCfg::ALLOC_DECAY_SHIFT;
      val<NT> pick = cands.one_hot();
      // Gate this pick by probability (pick 0 always fires)
      auto gated_pick = [&]() {
        if constexpr (K == 0 || DS == 0) {
          return pick;
        } else {
          constexpr u64 SHIFT = std::min(K * DS, u64(7));
          constexpr u64 THRESH = u64(1) << (8 - SHIFT);
          val<1> prob_ok = (rng < hard<THRESH>{});
          return val<NT>{pick & prob_ok.replicate(hard<NT>{}).concat()};
        }
      }();
      if constexpr (K + 1 >= MA) {
        return gated_pick;
      } else {
        pick.fanout(hard<2>{});
        val<NT> next_cands = [&]() -> val<NT> {
          if constexpr (AllocCfg::NON_CONSECUTIVE) {
            val<NT> neighbors = (pick << 1) | (pick >> 1);
            return (cands ^ pick) & ~neighbors;
          } else {
            return cands ^ pick;
          }
        }();
        return gated_pick |
               self.template operator()<K + 1>(self, next_cands, rng);
      }
    };
    arr<val<1>, NT> allocate = [&]() -> arr<val<1>, NT> {
      if constexpr (AllocCfg::MAX_ALLOC == 1) {
        val<NT> alloc_target_rev = collamask.one_hot();
        return alloc_target_rev.reverse().make_array(val<1>{});
      } else {
        collamask.fanout(hard<2>{});
        val<NT> result =
            multi_pick.template operator()<0>(multi_pick, collamask, alloc_rng);
        return result.reverse().make_array(val<1>{});
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
          TAMonitor<NT, N, MAX_TABLE_SIZE, USE_GSHARE, FB_CAPACITY, HYST_BANKS,
                    U_BANKS>::decode_provider(static_cast<u64>(t_match1));
      if (static_cast<u64>(alloc_trigger)) {
        if (static_cast<u64>(postmask) == 0)
          mon.record_alloc_blocked(); // fallback was provider, no candidates
        mon.record_alloc_cascade(prov_idx, at);
      }
    }
#endif

    // ---- Step 4: Fallback update ----
    // 4a: Standard — mispredict when fallback is provider → write actual dir
    val<1> fb_changed = actual_dir != val<PRED_BITS>{train_fb};
    val<1> fb_gate = do_train & t_m1[NT] & mispredict & fb_changed;
    // 4b: Train toward P2 — on mispredict when TAGE is provider,
    //     overwrite gshare with actual dir (like Tage.hpp P1→P2)
    if constexpr (FB_TRAIN_P2) {
      val<1> tage_is_prov = ~t_m1[NT];
      val<1> fb_p2_gate =
          do_train & tage_is_prov & mispredict & fb_changed & ~fb_gate;
      val<1> fb_any_write = fb_gate | fb_p2_gate;
      execute_if(fb_any_write, [&]() {
        fb_ctr.write(val<FB_IDX_BITS>{train_fb_idx}, actual_dir);
      });
    } else {
      execute_if(fb_gate, [&]() {
        fb_ctr.write(val<FB_IDX_BITS>{train_fb_idx}, actual_dir);
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
    train_pc.fanout(hard<NT>{});
    static_loop<NT>([&]<u64 I>() {
      auto &t = std::get<I>(tables);
      val<1> do_alloc = allocate[I];

      // Per-branch wrong bits for banked training
      val<N> per_branch_wrong = t_pp ^ actual_dir;

      // pred_ram: with banking, only flip pred bits where branch is weakconf
      // AND wrong
      if constexpr (HYST_BANKS == 1) {
        // Original: block-wide weakconf gate
        val<1> do_pred_update = t_m1[I] & t_pwc & any_provider_wrong;
        execute_if(do_train & (do_alloc | do_pred_update), [&]() {
          t.pred_ram.write(val<t.IDX_BITS>{train_idx[I]}, actual_dir,
                           hard<0>{});
        });
      } else {
        // Per-branch: build weakconf mask from banked hyst, only flip
        // weak+wrong bits
        constexpr u64 HW_P = std::max(u64(1), HYST_WIDTH);
        arr<val<1>, N> per_br_weakconf = [&](u64 b) -> val<1> {
          u64 hb = hyst_bank_of(b);
          return val<1>{val<HW_P>{train_hyst[I][hb]} == hard<0>{}};
        };
        val<N> weakconf_mask = per_br_weakconf.concat();
        val<N> flip_mask = per_branch_wrong & weakconf_mask;
        val<1> any_flip = t_m1[I] & (flip_mask != hard<0>{});
        val<PRED_BITS> new_pred = t_pp ^ flip_mask;
        execute_if(do_train & (do_alloc | any_flip), [&]() {
          t.pred_ram.write(val<t.IDX_BITS>{train_idx[I]},
                           select(do_alloc, actual_dir, new_pred), hard<0>{});
        });
      }

      // hyst_ram: per-bank training
      constexpr u64 HW = std::max(u64(1), HYST_WIDTH);
      static_loop<HYST_BANKS>([&]<u64 HB>() {
        auto old_hyst = val<HW>{train_hyst[I][HB]};
        // Bank wrong: OR of per_branch_wrong for branches in this bank
        val<1> bank_wrong =
            (per_branch_wrong & hard<hyst_bank_mask(HB)>{}) != hard<0>{};
        auto new_hyst = ta_update_ctr(old_hyst, ~bank_wrong);
        auto hyst_data = select(do_alloc, val<HW>{0}, new_hyst);
        val<1> do_hyst_update = t_m1[I] & (new_hyst != old_hyst);
        execute_if(do_train & (do_alloc | do_hyst_update), [&]() {
          t.hyst_ram[HB].write(val<t.HYST_IDX_BITS>{train_idx[I]}, hyst_data,
                               hard<0>{});
        });
#ifdef TAGE_MONITOR
        if (static_cast<u64>(do_train & do_hyst_update))
          mon.record_hyst_update(I, HB, !static_cast<u64>(bank_wrong));
#endif
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
      if (static_cast<u64>(do_train & do_alloc)) {
        u64 tidx = static_cast<u64>(val<t.IDX_BITS>{train_idx[I]});
        mon.record_tage_write(I, tidx);
        mon.record_entry_alloc_diag(I, tidx);
      }
#endif

      // u_ram: combined provider update + allocation + uclear + decay + DIP
      // Banked: each u bank gets per-bank wrong signal
      val<U_WIDTH> alloc_u = [&]() -> val<U_WIDTH> {
        if constexpr (DIP_PROB_256 == 0) {
          return val<U_WIDTH>{0};
        } else {
          val<1> dip_fire = (alloc_rng < hard<DIP_PROB_256>{});
          return select(dip_fire, val<U_WIDTH>{hard<DIP_INIT_U>{}},
                        val<U_WIDTH>{0});
        }
      }();

      static_loop<U_BANKS>([&]<u64 UB>() {
        // Per-bank provider_u: depends on U_PROV_UPDATE policy
        val<1> bank_u_wrong =
            (per_branch_wrong & hard<u_bank_mask(UB)>{}) != hard<0>{};
        val<1> bank_u_correct = ~bank_u_wrong;
        val<U_WIDTH> old_u = val<U_WIDTH>{train_u[I][UB]};
        constexpr u64 MAX_U = (1ULL << U_WIDTH) - 1;

        // provider_u: value to write; prov_u_write: should we write at all
        auto [provider_u,
              prov_u_write] = [&]() -> std::pair<val<U_WIDTH>, val<1>> {
          if constexpr (U_PROV_UPDATE == UProvUpdate::SET_OR_CLEAR) {
            // correct → max_u, wrong → 0, always write
            return {val<U_WIDTH>{
                        bank_u_correct.replicate(hard<U_WIDTH>{}).concat()},
                    val<1>{hard<1>{}}};
          } else if constexpr (U_PROV_UPDATE == UProvUpdate::SET_ON_CORRECT) {
            // correct → max_u, wrong → no write
            return {val<U_WIDTH>{hard<MAX_U>{}}, bank_u_correct};
          } else if constexpr (U_PROV_UPDATE == UProvUpdate::INC_DEC) {
            // correct → saturating u+1, wrong → saturating u-1
            val<U_WIDTH> inc =
                select(old_u == hard<MAX_U>{}, old_u, val<U_WIDTH>{old_u + 1});
            val<U_WIDTH> dec =
                select(old_u == hard<0>{}, old_u, val<U_WIDTH>{old_u - 1});
            val<U_WIDTH> new_val = select(bank_u_correct, inc, dec);
            return {new_val, new_val != old_u};
          } else { // INC_ONLY
            // correct → saturating u+1, wrong → no write
            val<U_WIDTH> inc =
                select(old_u == hard<MAX_U>{}, old_u, val<U_WIDTH>{old_u + 1});
            return {inc, bank_u_correct & (inc != old_u)};
          }
        }();

        struct UResult {
          val<U_WIDTH> newu;
          val<1> uw;
        };
        auto [base_newu, base_u_write] = [&]() -> UResult {
          if constexpr (ALT_PROMOTE_ACTIVE) {
            val<MATCH_BITS> t_match2_loc = train_match2;
            arr<val<1>, NT + 1> t_m2 = t_match2_loc.make_array(val<1>{});
            val<PRED_BITS> t_ap = train_alt_pred;
            val<1> any_alt_correct = (t_ap ^ actual_dir) == hard<0>{};
            val<1> alt_prom = AltPromoteFn::should_promote(
                any_provider_wrong, any_alt_correct,
                val<ALLOC_WIDTH>{alloc_ctr}, val<ACC_WIDTH>{acc_ctr},
                alloc_rng);
            val<1> do_alt_promote = alt_prom & t_m2[I];
            return {
                select(do_alt_promote, val<U_WIDTH>{1},
                       select(allocate[I], alloc_u,
                              select(uclear[I], val<U_WIDTH>{0}, provider_u))),
                (t_m1[I] & prov_u_write & t_ad) | allocate[I] | uclear[I] |
                    do_alt_promote};
          } else {
            return {select(allocate[I], alloc_u,
                           select(uclear[I], val<U_WIDTH>{0}, provider_u)),
                    (t_m1[I] & prov_u_write & t_ad) | allocate[I] | uclear[I]};
          }
        }();

        // Probabilistic decay
        auto [newu, u_write] = [&]() -> std::pair<val<U_WIDTH>, val<1>> {
          if constexpr (!DECAY_ENABLE) {
            return {base_newu, base_u_write};
          } else {
            constexpr u64 LW = DECAY_LFSR_WIDTHS[I];
            auto &lfsr = std::get<I>(decay_lfsrs);
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
            auto thresh = DecayThreshFn::template compute<I, LW>(
                val<ACC_WIDTH>{acc_ctr}, val<ALLOC_WIDTH>{alloc_ctr});
            val<1> decay_fire =
                decay_miss & ~allocate[I] & (val<LW>{lfsr} < thresh);
            val<U_WIDTH> old_u = val<U_WIDTH>{train_u[I][UB]};
            val<U_WIDTH> decayed_u = [&]() {
              if constexpr (DECAY_OP == DecayOp::DECREMENT)
                return select(old_u == hard<0>{}, old_u,
                              val<U_WIDTH>{old_u - 1});
              else if constexpr (DECAY_OP == DecayOp::HALVE)
                return val<U_WIDTH>{old_u >> 1};
              else
                return val<U_WIDTH>{0};
            }();
            val<U_WIDTH> merged = select(base_u_write, base_newu,
                                         select(decay_fire, decayed_u, old_u));
            val<1> merged_write = base_u_write | decay_fire;
            val<1> merged_changed = merged != old_u;
            return {merged, merged_write & merged_changed};
          }
        }();

        execute_if(do_train & u_write, [&]() {
          t.u_ram[UB].write(val<t.IDX_BITS>{train_idx[I]}, newu, hard<0>{});
        });
#ifdef TAGE_MONITOR
        if (static_cast<u64>(do_train & u_write))
          mon.record_u_write(I, UB,
                             static_cast<u64>(val<t.IDX_BITS>{train_idx[I]}),
                             static_cast<u64>(newu));
#endif
      }); // end static_loop<U_BANKS>
#ifdef TIMING_DEBUG
      dbg_tag_write[I] = do_train & do_alloc;
      // pred_write and hyst_write not easily tracked with banking; use alloc
      // gate
      dbg_pred_write[I] = do_train & do_alloc;
      dbg_hyst_write[I] = do_train & do_alloc;
      dbg_u_write[I] = do_train & do_alloc;
#endif
    });

    // ---- Global pressure counter updates ----
    if constexpr (DECAY_ENABLE || EPOCH_ENABLE) {
      // Accuracy counter: increment on correct, decrement on mispredict
      auto new_acc = ta_update_ctr(val<ACC_WIDTH>{acc_ctr}, ~mispredict);

      // Alloc pressure counter
      auto new_alloc = [&]() {
        if constexpr (FARALLOC_DIST > 0) {
          // Tage.hpp-style: only update on mispredict; inc on far alloc, dec on
          // close
          val<NT> at = alloc_target;
          at.fanout(hard<2>{});
          val<NT> shifted_prov = val<NT>{t_match1} >> FARALLOC_DIST;
          val<1> faralloc = ((shifted_prov | at).one_hot() ^ at) == hard<0>{};
          return select(mispredict,
                        ta_update_ctr(val<ALLOC_WIDTH>{alloc_ctr}, faralloc),
                        val<ALLOC_WIDTH>{alloc_ctr});
        } else {
          // Default: increment when no alloc slot found
          val<1> any_alloc = alloc_target != hard<0>{};
          return ta_update_ctr(val<ALLOC_WIDTH>{alloc_ctr}, ~any_alloc);
        }
      }();

      if constexpr (EPOCH_ENABLE) {
        // Epoch counter: free-running, incremented each update_cycle
        val<EPOCH_CTR_WIDTH> ectr = val<EPOCH_CTR_WIDTH>{epoch_ctr};

        // Epoch: bulk reset u_ram when trigger fires
        val<1> epoch_fire =
            EpochTriggerFn::template should_fire<ACC_WIDTH, ALLOC_WIDTH,
                                                 EPOCH_CTR_WIDTH>(
                val<ACC_WIDTH>{acc_ctr}, val<ALLOC_WIDTH>{alloc_ctr}, ectr);
        execute_if(epoch_fire, [&]() {
          static_loop<NT>([&]<u64 I>() {
            static_loop<U_BANKS>(
                [&]<u64 UB>() { std::get<I>(tables).u_ram[UB].reset(); });
          });
        });
        epoch_ctr = val<EPOCH_CTR_WIDTH>{ectr + 1};
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

// ---- Convenience typedefs for technique sweeps ----
// Technique 3: Alt bank promotion variants
using TageAhead_AltPromAlways =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 3, 2, 2, 8192, true, 6, 5, 256, 2,
              1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4, false,
              DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplaceUZero, AltPromoteAlways>;
using TageAhead_AltPromProb128 =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 3, 2, 2, 8192, true, 6, 5, 256, 2,
              1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4, false,
              DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplaceUZero, AltPromoteProb<128>>;
using TageAhead_AltPromPress =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 3, 2, 2, 8192, true, 6, 5, 256, 2,
              1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4, false,
              DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplaceUZero, AltPromotePressure>;

// Technique 4: Replacement policy variants
using TageAhead_ReplWeakConf1 =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 3, 2, 2, 8192, true, 6, 5, 256, 2,
              1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4, false,
              DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplaceUZeroWeakConf<1>>;
using TageAhead_ReplPressAdapt =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 3, 2, 2, 8192, true, 6, 5, 256, 2,
              1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4, false,
              DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplacePressureAdaptive<1, 4>>;

// Technique 6: DIP-like allocation variants
using TageAhead_DIP64 =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 3, 2, 2, 8192, true, 6, 5, 256, 2,
              1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4, false,
              DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplaceUZero, AltPromoteOff, 64, 1,
              UProvUpdate::SET_OR_CLEAR>;
using TageAhead_DIP128 =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 3, 2, 2, 8192, true, 6, 5, 256, 2,
              1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4, false,
              DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplaceUZero, AltPromoteOff, 128, 1,
              UProvUpdate::SET_OR_CLEAR>;

// Combined: best of each (alt promote + confidence gate + DIP)
using TageAhead_AllTechniques =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 3, 2, 2, 8192, true, 6, 5, 256, 2,
              1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4, false,
              DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplaceUZeroWeakConf<1>, AltPromoteAlways, 64, 1,
              UProvUpdate::SET_OR_CLEAR>;

// Per-branch banking variants (HYST_BANKS, HYST_BANK_BIT, U_BANKS, U_BANK_BIT)
// 2 hyst banks, bit 1 (pairs {0,1},{2,3},{4,5},{6,7} interleaved)
using TageAhead_HystBank2 =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 1, 3, 2, 8192 * 2, false, 6, 5,
              256, 2, 1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4,
              false, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplaceUZero, AltPromoteOff, 0, 1,
              UProvUpdate::SET_OR_CLEAR, false, 0, 2, 1, 1, 0>;
// 4 hyst+u banks, bit 0 (interleaved)
using TageAhead_Bank4 =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 1, 3, 2, 8192 * 2, false, 6, 5,
              256, 2, 1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4,
              false, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplaceUZero, AltPromoteOff, 0, 1,
              UProvUpdate::SET_OR_CLEAR, false, 0, 4, 0, 4, 0>;
// 8 hyst+u banks (fully per-branch)
using TageAhead_Bank8 =
    TageAhead<TATableConfig<>, 8, 6, 3, true, 1, 3, 2, 8192 * 2, false, 6, 5,
              256, 2, 1024, true, HistUpdate::PATH, TADefaultAllocConfig, 4, 4,
              false, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, TATableConfig<>::NUM_TABLES>(8),
              ta::DefaultDecayThresh, true, ta::DefaultEpochTrigger, 16, false,
              true, ReplaceUZero, AltPromoteOff, 0, 1,
              UProvUpdate::SET_OR_CLEAR, false, 0, 8, 0, 8, 0>;
