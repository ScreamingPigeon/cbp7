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

template <u64 N = 12, u64 SIZE = 1024, u64 TAG = 11, u64 MINH = 4,
          u64 MAXH = 500, u64 SIZE_RATIO = 4,
          ta::HistSeries HIST = ta::HistSeries::GEOMETRIC,
          typename TagFn = ta::UniformTag<TAG>,
          typename SizeFn = ta::GeoSize<SIZE, SIZE_RATIO>,
          auto BANK_SHIFTS = ta::uniform_array<u64, N>(1),
          auto BANK_COUNTS = ta::uniform_array<u64, N>(2)>
struct TATableConfig {
  static constexpr u64 NUM_TABLES = N;
  static constexpr u64 MINHIST = MINH;
  static constexpr u64 MAXHIST = MAXH;
  static constexpr auto TABLE_SIZE = ta::generate_table_sizes<N>(SizeFn{});
  static constexpr auto TAG_WIDTH = ta::generate_tag_widths<N>(TagFn{});
  static constexpr auto RWRAM_BANK_SHIFT = BANK_SHIFTS;
  static constexpr auto RWRAM_BANKS = BANK_COUNTS;
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
    u64 BR_P_ENTRY = N,        // branches per entry (N=all share, 1=each independent)
    bool INTERLEAVED = false,  // grouping: false=contiguous, true=interleaved
    u64 HYST_WIDTH = 2,        // hysteresis width (separate from ctr)
    u64 U_WIDTH = 2,           // usefulness counter width
    UMispPolicy U_MISP = UMispPolicy::UNTOUCHED, // u-bit on provider mispredict
    UClearPolicy U_CLEAR = UClearPolicy::DECREMENT, // u-bit on alloc failure
    u64 FB_CAPACITY = 8192,   // fallback table size (bimodal or gshare)
    bool USE_GSHARE = false,  // use gshare base (PC^history) vs bimodal (PC)
    u64 GS_HIST = 6,          // gshare history length (only when USE_GSHARE)
    u64 META_WIDTH = 2,       // meta counter width (provider vs alt)
    u64 META_CAPACITY = 1024, // meta table entries
    u64 META_PIPE = 2,        // meta pipeline depth
    u64 LINEINST = 512,       // line size in instructions
    bool SHARED_HYS = true,   // shared hyst: 2 entries share 1 counter
    HistUpdate HIST_MODE =
        HistUpdate::PATH, // what goes into history: PATH, DIR, or BOTH
    // ---- Allocation policy ----
    typename AllocCfg = TADefaultAllocConfig,
    // ---- Sibling skip policy ----
    SiblingPolicy SIBLING_POLICY =
        SiblingPolicy::ALL,      // NONE=never skip, ALL=skip siblings
    u64 SIBLING_TABLE_FLOOR = 0, // skip siblings only for tables >= this index
    // ---- Global pressure counters ----
    u64 ACC_WIDTH = 10,   // accuracy counter width
    u64 ALLOC_WIDTH = 10, // alloc pressure counter width
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
    // ---- Fallback bimodal hysteresis ----
    // Traditional 2-bit saturating counter for fb: pred only flips when hyst
    // is weak (wrong twice). Mirrors Tage.hpp's bhyst. Half-sized RAM (shared).
    bool FB_BIM_HYST = false,
    // ---- Far-allocation pressure ----
    // When > 0, allocation distance >= FARALLOC_DIST from provider biases
    // alloc_ctr harder toward epoch/decay (extra decrement).
    // Mirrors Tage.hpp's faralloc-based uctr adaptation.
    u64 FARALLOC_DIST = 0,
    // ---- Floorplan tuning ----
    // Reverse RAM declaration order so biggest tables are declared first.
    bool REVERSE_TABLE_ORDER = false,
    // ---- Per-table sec-tag policy ----
    // Controls which tables require sec-tag match. Only used when
    // USE_SEC_TAG=true and NUM_PATHS==1 (single-chain).
    typename SecTagPolicy = ta::SecTagAll,
    // ---- Fallback banking ----
    u64 FB_BANKS = 1,         // address-banked fb: 1=single RAM, 2/4/8=banked
    u64 FB_BANK_SEL_SHIFT = 0, // bank select from fb_idx bits [SHIFT+SEL_BITS-1:SHIFT]
    bool FB_SPLIT_WIDTH = false, // split fb into per-group 1-bit RAMs
    // ---- Meta banking ----
    u64 META_BANKS = 1,       // separate meta RAMs (each serves N/META_BANKS branches)
    // ---- PC index shift ----
    u64 PC_IDX_SHIFT = 2>     // inst_pc >> PC_IDX_SHIFT for table index (2=word, 5=block)
struct TageAhead : predictor {

  // ======== Derived Constants ========

  static_assert(!(DECAY_ENABLE && EPOCH_ENABLE),
                "DECAY_ENABLE and EPOCH_ENABLE are mutually exclusive");
  // static_assert(NUM_PATHS == 1 ||
  //                   (USE_SEC_TAG && NUM_PATHS == (u64(1) << SEC_TAG_BITS)),
  //               "NUM_PATHS must be 1 or 2^SEC_TAG_BITS with USE_SEC_TAG");
  static_assert(NUM_PATHS <= 4, "NUM_PATHS > 4 not yet supported");
  static_assert(BR_P_ENTRY >= 1 && BR_P_ENTRY <= N,
                "BR_P_ENTRY must be in [1, N]");

  static constexpr u64 NT = TableCfg::NUM_TABLES;
  static constexpr u64 LOGLINEINST = ta::clog2(LINEINST);
  static constexpr u64 MAXHIST = TableCfg::MAXHIST;
  static constexpr u64 MAX_TAG_WIDTH = ta::array_max(TableCfg::TAG_WIDTH);
  static constexpr u64 MAX_TABLE_SIZE = ta::array_max(TableCfg::TABLE_SIZE);
  static constexpr u64 MAX_IDX_BITS = ta::clog2(MAX_TABLE_SIZE);
  static constexpr u64 MAX_BANKS = ta::array_max(TableCfg::RWRAM_BANKS);

  static constexpr u64 MATCH_BITS = NT + 1; // NT tables + fallback

  // ---- Per-group tag encoding (BR_P_ENTRY) ----
  static constexpr u64 NUM_GROUPS =
      (BR_P_ENTRY == N) ? 1 : (N + BR_P_ENTRY - 1) / BR_P_ENTRY;
  static constexpr u64 GROUP_BITS =
      (NUM_GROUPS > 1) ? ta::clog2(NUM_GROUPS) : 0;
  static constexpr u64 MAX_HTAG_WIDTH = MAX_TAG_WIDTH - GROUP_BITS;
  // Ensure enough tag bits remain for htag after group_id
  static_assert(GROUP_BITS == 0 || MAX_TAG_WIDTH > GROUP_BITS,
                "TAG_WIDTH must be > GROUP_BITS when BR_P_ENTRY < N");
  // Multi-path × multi-group not yet supported
  static_assert(NUM_PATHS == 1 || NUM_GROUPS == 1,
                "NUM_PATHS > 1 with NUM_GROUPS > 1 not yet supported");

  // Map branch position → group index
  static constexpr auto BRANCH_TO_GROUP = []() {
    std::array<u64, N> b2g{};
    for (u64 i = 0; i < N; i++) {
      if constexpr (INTERLEAVED)
        b2g[i] = i % NUM_GROUPS;
      else
        b2g[i] = i / BR_P_ENTRY;
    }
    return b2g;
  }();

  // Map branch position → position within its group
  static constexpr auto BRANCH_IN_GROUP = []() {
    std::array<u64, N> big{};
    for (u64 i = 0; i < N; i++) {
      if constexpr (INTERLEAVED)
        big[i] = i / NUM_GROUPS;
      else
        big[i] = i % BR_P_ENTRY;
    }
    return big;
  }();

  // Number of branches in each group (last group may be smaller)
  static constexpr auto GROUP_SIZE = []() {
    std::array<u64, NUM_GROUPS> gs{};
    for (u64 i = 0; i < N; i++)
      gs[BRANCH_TO_GROUP[i]]++;
    return gs;
  }();

  // First branch index in each group (for contiguous fb extraction)
  static constexpr auto GROUP_START = []() {
    std::array<u64, NUM_GROUPS> s{};
    for (u64 g = 0; g < NUM_GROUPS; g++) {
      for (u64 i = 0; i < N; i++) {
        if (BRANCH_TO_GROUP[i] == g) { s[g] = i; break; }
      }
    }
    return s;
  }();

  // ---- Meta banking (per-branch meta counters) ----
  static_assert(META_BANKS >= 1 && META_BANKS <= N,
                "META_BANKS must be in [1, N]");
  static_assert(META_BANKS == 1 || NUM_GROUPS > 1,
                "META_BANKS > 1 requires NUM_GROUPS > 1 (BPE < N)");
  // Map branch position → meta bank index
  static constexpr auto BRANCH_TO_META_BANK = []() {
    std::array<u64, N> b2m{};
    for (u64 i = 0; i < N; i++)
      b2m[i] = i * META_BANKS / N;
    return b2m;
  }();
  // Map group → meta bank (for NUM_GROUPS>1: group's first branch determines)
  static constexpr auto GROUP_TO_META_BANK = []() {
    std::array<u64, NUM_GROUPS> g2m{};
    for (u64 g = 0; g < NUM_GROUPS; g++)
      g2m[g] = GROUP_START[g] * META_BANKS / N;
    return g2m;
  }();

  // ---- Sec-tag adaptive benefit counter ----
  static constexpr u64 BENEFIT_WIDTH = [] {
    if constexpr (requires { SecTagPolicy::BENEFIT_WIDTH; })
      return SecTagPolicy::BENEFIT_WIDTH;
    else
      return u64(0);
  }();
  static constexpr bool HAS_BENEFIT_CTR = BENEFIT_WIDTH > 0;
  static constexpr u64 BENEFIT_REG_WIDTH = std::max(u64(1), BENEFIT_WIDTH);

  // Per-bit gh fanout: only fanout bits each table actually reads
  static constexpr auto GH_FANOUT = []() {
    auto fo = ta::gh_per_bit_fanout<MAXHIST, NT, TableCfg::HIST_LEN>();
    if constexpr (USE_GSHARE)
      fo[GS_HIST - 1] += 1; // fb_fold reads gh[GS_HIST-1]
    return fo;
  }();

  // Prediction bits per entry: one CTR_WIDTH counter per branch-in-group
  static constexpr u64 PRED_BITS = BR_P_ENTRY * CTR_WIDTH;
  // Fallback stays N-wide (direct-mapped, no tags/groups)
  static constexpr u64 FB_PRED_BITS = N * CTR_WIDTH;

  // ---- Fallback banking ----
  static_assert(FB_BANKS >= 1 && (FB_BANKS & (FB_BANKS - 1)) == 0,
                "FB_BANKS must be a power of 2");
  static_assert(FB_CAPACITY % FB_BANKS == 0,
                "FB_CAPACITY must be divisible by FB_BANKS");
  static constexpr u64 FB_BANK_SEL_BITS = (FB_BANKS > 1) ? ta::clog2(FB_BANKS) : 0;
  static constexpr u64 FB_BANK_CAPACITY = FB_CAPACITY / FB_BANKS;
  static constexpr u64 FB_BANK_IDX_BITS = ta::clog2(FB_BANK_CAPACITY);
  // Upper bits above the bank-select field: fb_idx[IDX-1 : SHIFT+SEL]
  static constexpr u64 FB_BANK_UPPER_BITS =
      (FB_BANKS > 1) ? (ta::clog2(FB_CAPACITY) - FB_BANK_SEL_SHIFT - FB_BANK_SEL_BITS) : 0;
  static_assert(FB_BANKS == 1 ||
                FB_BANK_SEL_SHIFT + FB_BANK_SEL_BITS <= ta::clog2(FB_CAPACITY),
                "FB_BANK_SEL_SHIFT + SEL_BITS exceeds FB_IDX_BITS");
  static_assert(FB_BANKS == 1 ||
                FB_BANK_UPPER_BITS + FB_BANK_SEL_SHIFT == FB_BANK_IDX_BITS,
                "bank index bits must equal upper + lower bits around select field");
  // Width per fb RAM: 1-bit per group when FB_SPLIT_WIDTH, else N-wide
  static constexpr u64 FB_RAM_WIDTH = FB_SPLIT_WIDTH ? CTR_WIDTH : FB_PRED_BITS;
  static constexpr u64 FB_NUM_WIDTH_RAMS = FB_SPLIT_WIDTH ? NUM_GROUPS : 1;

  // ---- Fallback bimodal hysteresis ----
  static constexpr u64 FB_BH_CAPACITY = FB_CAPACITY / 2; // half-sized (shared hyst)
  static constexpr u64 FB_BH_IDX_BITS = ta::clog2(FB_BH_CAPACITY);
  static_assert(!FB_SPLIT_WIDTH || NUM_GROUPS > 1,
                "FB_SPLIT_WIDTH requires NUM_GROUPS > 1 (BPE < N)");

  // ======================================================================
  // Predict-path Registers (declared before tables → at first table loc)
  // ======================================================================
  // ---- Global History (shared, folds live in per-table TATable) ----
  ta_global_history<MAXHIST> gh;

  // ---- Fallback Predictor (ahead-pipelined) ----
  // USE_GSHARE=false: bimodal (PC-indexed)
  // USE_GSHARE=true:  gshare (PC ^ folded_history indexed)
  reg<FB_PRED_BITS> prefetch_fb;
  reg<FB_PRED_BITS> current_fb;
  static constexpr u64 FB_IDX_BITS = ta::clog2(FB_CAPACITY);
  reg<FB_IDX_BITS> prefetch_fb_idx;
  reg<FB_IDX_BITS> current_fb_idx;
  // Piped fb hyst for reconciliation (only accessed when FB_RECONCILE=true)
  reg<1> prefetch_fb_hyst;
  reg<1> current_fb_hyst;
  // fb bimodal hyst read value (read in update_cycle before extra_cycle,
  // written after extra_cycle — no predict1 pipeline needed)
  reg<FB_PRED_BITS> fb_bh_read;

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
  reg<SEC_TAG_BITS> curr_sec_tag;

  // ---- Pipeline Regs [NT] ----
  // prefetch_*: written in predict1 (ahead reads for next block)
  reg<MAX_TAG_WIDTH> prefetch_tag[NT];
  reg<1> prefetch_tag_hit[NT]; // full tag match (NUM_GROUPS==1) or htag match
  reg<PRED_BITS> prefetch_pred[NT];
  reg<SEC_TAG_BITS> prefetch_sec[NT];
  reg<MAX_IDX_BITS> prefetch_idx[NT];
  reg<std::max(u64(1), HYST_WIDTH)> prefetch_hyst[NT];
  reg<U_WIDTH> prefetch_u[NT];
  reg<MAX_TAG_WIDTH> prefetch_ctag[NT]; // computed tag for allocation piping
  // Per-group tag encoding: stored group_id extracted from tag upper bits
  reg<std::max(u64(1), GROUP_BITS)> prefetch_group_id[NT];

  // current_*: shifted from prefetch, used for prediction
  reg<MAX_TAG_WIDTH> current_tag[NT];
  reg<1> current_tag_hit[NT];
  reg<PRED_BITS> current_pred[NT];
  reg<SEC_TAG_BITS> current_sec[NT];
  reg<MAX_IDX_BITS> current_idx[NT];
  reg<std::max(u64(1), HYST_WIDTH)> current_hyst[NT];
  reg<U_WIDTH> current_u[NT];
  reg<MAX_TAG_WIDTH> current_ctag[NT]; // piped computed tag for allocation
  reg<std::max(u64(1), GROUP_BITS)> current_group_id[NT];

  // ---- Prediction reg (shared by
  // predict1/reuse_predict1/predict2/reuse_predict2) ----
  arr<reg<1>, N> pred;

  // ---- Meta pipeline (shifted each update_cycle) ----
  static constexpr u64 META_IDX_BITS = ta::clog2(META_CAPACITY);
  reg<META_WIDTH, i64> meta_pipe[META_BANKS][META_PIPE];
  reg<META_IDX_BITS> meta_idx_pipe[META_BANKS][META_PIPE];

  // ======================================================================
  // Storage
  // ======================================================================
  // ---- Table tuple (per-table tag width and table size) ----
  using Tables = typename TAMakeTableTuple<
      TableCfg, CTR_WIDTH, HYST_WIDTH, U_WIDTH, SEC_TAG_BITS, BR_P_ENTRY,
      SHARED_HYS, REVERSE_TABLE_ORDER, std::make_index_sequence<NT>>::type;
  Tables tables;

  // Access logical table I (remaps through reversed storage when needed)
  template <u64 I> auto &table() {
    return std::get<REVERSE_TABLE_ORDER ? (NT - 1 - I) : I>(tables);
  }
  template <u64 I> const auto &table() const {
    return std::get<REVERSE_TABLE_ORDER ? (NT - 1 - I) : I>(tables);
  }
  // Fallback RAMs: FB_NUM_WIDTH_RAMS × FB_BANKS, each FB_BANK_CAPACITY entries
  // When FB_BANKS=1 && !FB_SPLIT_WIDTH: single N-wide RAM (baseline)
  // When FB_BANKS>1: address-banked, only one bank active per access
  // When FB_SPLIT_WIDTH: per-group 1-bit RAMs
  hcm::ram<val<FB_RAM_WIDTH>, FB_BANK_CAPACITY> fb_ctr[FB_NUM_WIDTH_RAMS][FB_BANKS];

  // ======================================================================
  // Update-only RAMs and Registers (in update zone)
  // ======================================================================
  hcm::zone update_zone;
  // Fallback hysteresis: tracks agreement between fb and TAGE.
  // hyst=1 → agree, hyst=0 → disagree (weak → eligible for reconciliation).
  // Only accessed when FB_RECONCILE=true; zero cost otherwise.
  hcm::ram<val<1>, FB_CAPACITY> fb_hyst{"fb_hyst"};
  // Fallback bimodal hyst: per-branch confidence (1=strong, 0=weak).
  // Only flip fb pred when hyst is weak. Half-sized (shared).
  hcm::ram<val<FB_PRED_BITS>, FB_BH_CAPACITY> fb_bim_hyst{"fb_bim_hyst"};
  ta_rwram<META_WIDTH, META_CAPACITY, 2> meta_ctr[META_BANKS];

  // train_*: saved from current_* before pipeline shift, used for training
  // (training is for block A, resolution is for block B)
  reg<MAX_IDX_BITS> train_idx[NT];
  reg<PRED_BITS>
      train_pred[NT]; // per-table pred for per-table training signals
  reg<std::max(u64(1), HYST_WIDTH)> train_hyst[NT];
  reg<U_WIDTH> train_u[NT];
  reg<FB_PRED_BITS> train_fb;
  reg<FB_IDX_BITS> train_fb_idx;
  reg<1> train_fb_hyst; // piped fb hyst for reconciliation
  reg<ALLOC_PC_BITS> train_pc;
  reg<MAX_TAG_WIDTH> train_ctag[NT]; // piped computed tag for allocation
  reg<SEC_TAG_BITS> train_sec_tag;   // piped sec_tag for allocation writes
  reg<std::max(u64(1), GROUP_BITS)> train_group_id[NT]; // per-group tag encoding

  // Piped resolution values from previous update_cycle (block A's resolution)
  // When NUM_GROUPS > 1, these are per-group arrays.
  reg<MATCH_BITS> train_match1[NUM_GROUPS];
  reg<PRED_BITS> train_provider_pred[NUM_GROUPS];
  reg<1> train_provider_weak[NUM_GROUPS];
  reg<1> train_altdiff[NUM_GROUPS];

  // Guard: skip training until piped resolution regs have been populated
  reg<1> train_valid = 0;

  // ---- U-bit decay (probabilistic) ----
  // Piped tag/sec hit for decay miss detection
  reg<1> train_tag_hit[NT];
  reg<1> train_sec_hit[NT];

  // Global pressure counters (always declared, zero-cost when unused)
  reg<ACC_WIDTH> acc_ctr;
  reg<ALLOC_WIDTH> alloc_ctr;

  // Sec-tag adaptive benefit counter (dead 1-bit reg when policy doesn't use it)
  reg<BENEFIT_REG_WIDTH> benefit_ctr;

  // Decay LFSR widths (used to parameterize std::rand() masking)
  static constexpr u64 MAX_LFSR_WIDTH =
      DECAY_ENABLE ? ta::array_max(DECAY_LFSR_WIDTHS) : 1;

#ifdef TAGE_MONITOR
  TAMonitor<NT, N, MAX_TABLE_SIZE, USE_GSHARE, FB_CAPACITY, NUM_GROUPS, MAX_BANKS, META_BANKS> mon;
  ~TageAhead() {
    mon.print_summary();
    print_params(std::cerr);
    // rwram bank conflict stats
    std::cerr << "\n=== rwram bank conflict stats ===\n";
    static_loop<NT>([&]<u64 I>() {
      auto &t = table<I>();
      t.pred_ram.print_conflict_stats();
      t.hyst_ram.print_conflict_stats();
      t.u_ram.print_conflict_stats();
    });
    for (u64 mb = 0; mb < META_BANKS; mb++)
      meta_ctr[mb].print_conflict_stats();
    std::cerr << "=================================\n";
  }
  void print_params(std::ostream &os) const {
    os << "\n=== TageAhead Parameters ===\n";
    os << "Tables: " << NT << "  N: " << N << "  PathBits: " << PATHBITS
       << "  HistMode: "
       << (HIST_MODE == HistUpdate::PATH  ? "PATH"
           : HIST_MODE == HistUpdate::DIR ? "DIR"
                                          : "BOTH")
       << "\n";
    os << "HistLen: [";
    for (u64 i = 0; i < NT; i++)
      os << (i ? "," : "") << TableCfg::HIST_LEN[i];
    os << "]\n";
    os << "TableSize: [";
    for (u64 i = 0; i < NT; i++)
      os << (i ? "," : "") << TableCfg::TABLE_SIZE[i];
    os << "]\n";
    os << "TagWidth: [";
    for (u64 i = 0; i < NT; i++)
      os << (i ? "," : "") << TableCfg::TAG_WIDTH[i];
    os << "]\n";
    os << "CtrWidth: " << CTR_WIDTH << "  HystWidth: " << HYST_WIDTH
       << "  SharedHyst: " << (SHARED_HYS ? "yes" : "no") << "\n";
    os << "SecTag: " << (USE_SEC_TAG ? "yes" : "no")
       << "  SecTagBits: " << SEC_TAG_BITS << "  NumPaths: " << NUM_PATHS
       << "\n";
    os << "UWidth: " << U_WIDTH << "  UMisp: "
       << (U_MISP == UMispPolicy::UNTOUCHED ? "UNTOUCHED"
           : U_MISP == UMispPolicy::ZERO    ? "ZERO"
                                            : "DECREMENT")
       << "  UClear: "
       << (U_CLEAR == UClearPolicy::DISABLED ? "DISABLED"
           : U_CLEAR == UClearPolicy::ZERO   ? "ZERO"
                                             : "DECREMENT")
       << "\n";
    os << "Fallback: " << (USE_GSHARE ? "Gshare" : "Bimodal")
       << "  Capacity: " << FB_CAPACITY;
    if (USE_GSHARE)
      os << "  GsHist: " << GS_HIST;
    os << "  Reconcile: " << (FB_RECONCILE ? "yes" : "no") << "\n";
    os << "Meta: width=" << META_WIDTH << " cap=" << META_CAPACITY
       << " pipe=" << META_PIPE << " banks=" << META_BANKS << "\n";
    os << "Alloc: trigger="
       << (AllocCfg::ALLOC_TRIGGER == AllocTrigger::MISPREDICT ? "MISPREDICT"
           : AllocCfg::ALLOC_TRIGGER == AllocTrigger::TAGE_WRONG ? "TAGE_WRONG"
                                                                  : "ALWAYS")
       << "  action="
       << (AllocCfg::ALLOC_ACTION == AllocAction::STANDARD   ? "STANDARD"
           : AllocCfg::ALLOC_ACTION == AllocAction::FILTERED ? "FILTERED"
                                                             : "THROTTLED")
       << "  target=" << AllocCfg::TARGET_POLICY::name()
       << "  maxAlloc=" << AllocCfg::MAX_ALLOC
       << "  nonConsec=" << (AllocCfg::NON_CONSECUTIVE ? "yes" : "no")
       << "\n";
    os << "  SiblingPolicy="
       << (SIBLING_POLICY == SiblingPolicy::NONE ? "NONE" : "ALL")
       << "  SiblingFloor=" << SIBLING_TABLE_FLOOR
       << "  FarAllocDist=" << FARALLOC_DIST << "\n";
    os << "Pressure: AccWidth=" << ACC_WIDTH << "  AllocWidth=" << ALLOC_WIDTH
       << "\n";
    os << "Epoch: " << (EPOCH_ENABLE ? "yes" : "no");
    if (EPOCH_ENABLE)
      os << "  ResetAcc=" << (EPOCH_RESET_ACC ? "yes" : "no")
         << "  ResetAlloc=" << (EPOCH_RESET_ALLOC ? "yes" : "no");
    os << "\n";
    os << "Decay: " << (DECAY_ENABLE ? "yes" : "no");
    if (DECAY_ENABLE) {
      os << "  Miss="
         << (DECAY_MISS == DecayMiss::TAG            ? "TAG"
             : DECAY_MISS == DecayMiss::SEC          ? "SEC"
             : DECAY_MISS == DecayMiss::TAG_OR_SEC   ? "TAG_OR_SEC"
                                                     : "TAG_AND_SEC")
         << "  Op="
         << (DECAY_OP == DecayOp::DECREMENT ? "DECREMENT"
             : DECAY_OP == DecayOp::HALVE   ? "HALVE"
                                            : "CLEAR")
         << "  LfsrWidths=[";
      for (u64 i = 0; i < NT; i++)
        os << (i ? "," : "") << DECAY_LFSR_WIDTHS[i];
      os << "]";
    }
    os << "\n";
    os << "LineInst: " << LINEINST
       << "  ReverseTableOrder: " << (REVERSE_TABLE_ORDER ? "yes" : "no")
       << "\n";
    os << "=== End Parameters ===\n\n";
  }
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
    mon.table_sizes = TableCfg::TABLE_SIZE;
    mon.table_bank_shift = TableCfg::RWRAM_BANK_SHIFT;
    mon.record_predict1();
    mon.shadow_block_pc = static_cast<u64>(inst_pc);
#endif
    // NOTE: @prakhar @claude audit
    inst_pc.fanout(
        hard<2 * NT + 2>{}); // 2 reads per table (>>2, >>4) + fb (>>2)
                             // + prefetch_pc (>>2)

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
      auto &t = table<I>();
      // NOTE: @prakhar @claude audit
      t.fold_idx.fanout(hard<2>{}); // get() + compute_update
      // NOTE: @prakhar @claude audit
      t.fold_tag.fanout(hard<2>{}); // get() + compute_update
      auto fold_idx_val = t.fold_idx.get();
      auto idx = fold_idx_val.fo1() ^ val<t.IDX_BITS>{inst_pc >> PC_IDX_SHIFT};
      // NOTE: @prakhar @claude audit
      idx.fanout(hard<6>{}); // 5 RAM reads + prefetch_idx write
      auto fold_tag_val = t.fold_tag.get();
      auto computed_tag = fold_tag_val.fo1() ^ val<t.tag_width>{inst_pc >> (PC_IDX_SHIFT + 2)};
      // NOTE: @prakhar @claude audit
      computed_tag.fanout(hard<2>{}); // tag comparison + prefetch_ctag write

      auto stored_tag = t.tag_ram.read(idx);
      // NOTE: @prakhar @claude audit
      stored_tag.fanout(hard<(NUM_GROUPS > 1) ? 3 : 2>{});
      prefetch_tag[I] = stored_tag;
      if constexpr (NUM_GROUPS > 1) {
        // Split tag comparison: compare only htag portion (lower bits)
        static constexpr u64 PER_HTAG = t.tag_width - GROUP_BITS;
        prefetch_tag_hit[I] =
            val<PER_HTAG>{stored_tag} == val<PER_HTAG>{computed_tag};
        // Extract stored group_id from upper bits of tag
        prefetch_group_id[I] =
            val<GROUP_BITS>{stored_tag >> PER_HTAG};
      } else {
        prefetch_tag_hit[I] =
            val<MAX_TAG_WIDTH>{stored_tag} == val<MAX_TAG_WIDTH>{computed_tag};
      }
      prefetch_ctag[I] = computed_tag;
      prefetch_pred[I] = t.pred_ram.read(idx);
      if constexpr (USE_SEC_TAG)
        prefetch_sec[I] = t.sec_ram.read(idx);
      prefetch_idx[I] = idx;
      prefetch_hyst[I] = t.hyst_ram.read(val<t.HYST_IDX_BITS>{idx});
      prefetch_u[I] = t.u_ram.read(idx);

#ifdef DEBUG_PRINT
      if constexpr (I == 0) {
        std::cerr << "\n───────────────────────────────────────\n";
        std::cerr << "=== predict1 table[0] ===\n";
        inst_pc.print(" pc=", "\n", true, std::cerr);
        fold_idx_val.print("  fold_idx_val=", "\n", true, std::cerr);
        idx.print("  idx=", "\n", true, std::cerr);
        fold_tag_val.print("  fold_tag_val=", "\n", true, std::cerr);
        computed_tag.print("  computed_tag=", "\n", true, std::cerr);
        stored_tag.print("  stored_tag=", "\n", true, std::cerr);
        prefetch_tag_hit[I].print("  tag_hit=", "\n", true, std::cerr);
        prefetch_sec[I].print("  sec=", "\n", true, std::cerr);
        prefetch_pred[I].print("  pred=", "\n", true, std::cerr);
        prefetch_hyst[I].print("  hyst=", "\n", true, std::cerr);
        prefetch_u[I].print("  u=", "\n", true, std::cerr);
      }
#endif
    });

    // Fallback ahead read (direct-mapped, no tag match needed)
    // USE_GSHARE: index = PC ^ folded_history; bimodal: index = PC
    auto fb_idx = [&]() {
      if constexpr (USE_GSHARE) {
        // NOTE: @prakhar @claude audit
        fb_fold.fanout(hard<2>{}); // get() + compute_update
        auto fb_fold_val = fb_fold.get();
        return val<FB_IDX_BITS>{inst_pc >> PC_IDX_SHIFT} ^ fb_fold_val.fo1();
      } else
        return val<FB_IDX_BITS>{inst_pc >> PC_IDX_SHIFT};
    }();
    // NOTE: @prakhar @claude audit
    fb_idx.fanout(
        hard<2 + (FB_RECONCILE ? 1 : 0)>{}); // prefetch_fb_idx + fb_ctr.read [+
                                             // fb_hyst.read]
    prefetch_fb_idx = fb_idx;

    // FB read: dispatch based on banking config
    if constexpr (FB_BANKS == 1 && !FB_SPLIT_WIDTH) {
      // Baseline: single N-wide RAM
      prefetch_fb = fb_ctr[0][0].read(fb_idx);
    } else if constexpr (FB_BANKS == 1 && FB_SPLIT_WIDTH) {
      // Width-split only: NUM_GROUPS per-group 1-bit RAMs, no address banking
      arr<val<1>, NUM_GROUPS> fb_bits = [&](u64 g) -> val<1> {
        return fb_ctr[g][0].read(fb_idx);
      };
      prefetch_fb = fb_bits.concat();
    } else if constexpr (FB_BANKS > 1 && !FB_SPLIT_WIDTH) {
      // Address-banked only: read all banks, MUX by bank_sel
      auto fb_bank_sel = val<FB_BANK_SEL_BITS>{fb_idx >> FB_BANK_SEL_SHIFT};
      // Remove bank-select bits: idx = {upper, lower} around sel field
      auto fb_bank_idx = [&]() {
        if constexpr (FB_BANK_SEL_SHIFT == 0)
          return val<FB_BANK_IDX_BITS>{fb_idx >> FB_BANK_SEL_BITS};
        else
          return val<FB_BANK_IDX_BITS>{
            concat(val<FB_BANK_UPPER_BITS>{fb_idx >> (FB_BANK_SEL_SHIFT + FB_BANK_SEL_BITS)},
                   val<FB_BANK_SEL_SHIFT>{fb_idx})};
      }();
      arr<val<FB_PRED_BITS>, FB_BANKS> bank_reads = [&](u64 b) -> val<FB_PRED_BITS> {
        return fb_ctr[0][b].read(fb_bank_idx);
      };
      prefetch_fb = bank_reads.select(fb_bank_sel);
    } else {
      // Both: read all banks per group, MUX by bank_sel
      auto fb_bank_sel = val<FB_BANK_SEL_BITS>{fb_idx >> FB_BANK_SEL_SHIFT};
      auto fb_bank_idx = [&]() {
        if constexpr (FB_BANK_SEL_SHIFT == 0)
          return val<FB_BANK_IDX_BITS>{fb_idx >> FB_BANK_SEL_BITS};
        else
          return val<FB_BANK_IDX_BITS>{
            concat(val<FB_BANK_UPPER_BITS>{fb_idx >> (FB_BANK_SEL_SHIFT + FB_BANK_SEL_BITS)},
                   val<FB_BANK_SEL_SHIFT>{fb_idx})};
      }();
      arr<val<1>, NUM_GROUPS> fb_bits = [&](u64 g) -> val<1> {
        arr<val<1>, FB_BANKS> grp_bank_reads = [&](u64 b) -> val<1> {
          return fb_ctr[g][b].read(fb_bank_idx);
        };
        return grp_bank_reads.select(fb_bank_sel);
      };
      prefetch_fb = fb_bits.concat();
    }

    if constexpr (FB_RECONCILE)
      prefetch_fb_hyst = fb_hyst.read(fb_idx);
    prefetch_pc = val<ALLOC_PC_BITS>{inst_pc >> 2};

#ifdef DEBUG_PRINT
    prefetch_fb.print("  fb_pred=", "\n", true, std::cerr);
    fb_idx.print("  fb_idx=", "\n", true, std::cerr);
#endif

    // Crit path: just read precomputed prediction from reg
    // NOTE: @prakhar @claude audit
    block_entry.fanout(hard<2 * LINEINST>{}); // line_end() across predict +
                                              // reuse + update_condbr
    // NOTE: @prakhar @claude audit
    // predict1/2, reuse calls [+ 1 old_pred read under TAGE_MONITOR]
    pred.fanout(hard<2 * LINEINST
#ifdef TAGE_MONITOR
                     + 1
#endif
                     >{});
    block_size = 1;
    num_branch = 0;
    reuse_prediction(~line_end());
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
    // curr_sec_tag read NT times (train_sec_hit) + 1 (train_sec_tag save)
    // NOTE: @prakhar @claude audit
    if constexpr (USE_SEC_TAG)
      curr_sec_tag.fanout(
          hard<NT + 1>{}); // NT train_sec_hit comparisons + train_sec_tag save
#ifdef TAGE_MONITOR
    // Capture raw sec values and tag hits before they're consumed
    u64 mon_saved_sec[NT]{};
    bool mon_tag_hit[NT]{};
    bool mon_sec_hit[NT]{};
#endif
    static_loop<NT>([&]<u64 I>() {
      train_idx[I] = current_idx[I].fo1();
      train_pred[I] =
          current_pred[I].fo1(); // per-table pred for per-table wrong signals
      train_hyst[I] = current_hyst[I].fo1();
      train_u[I] = current_u[I].fo1();
      train_ctag[I] = current_ctag[I].fo1();
      train_tag_hit[I] = current_tag_hit[I].fo1();
      if constexpr (NUM_GROUPS > 1)
        train_group_id[I] = current_group_id[I].fo1();
      if constexpr (USE_SEC_TAG) {
        train_sec_hit[I] = (val<SEC_TAG_BITS>{current_sec[I].fo1()} ==
                            val<SEC_TAG_BITS>{curr_sec_tag});
#ifdef TAGE_MONITOR
        mon_saved_sec[I] = static_cast<u64>(current_sec[I]);
        mon_tag_hit[I] = static_cast<u64>(current_tag_hit[I]);
        mon_sec_hit[I] = static_cast<u64>(train_sec_hit[I]);
#endif
      }
    });
    train_fb = current_fb.fo1();
    train_fb_idx = current_fb_idx.fo1();
    if constexpr (FB_RECONCILE)
      train_fb_hyst = current_fb_hyst.fo1();
    train_pc = current_pc.fo1();
    if constexpr (USE_SEC_TAG)
      train_sec_tag = curr_sec_tag;

    // Fanout on prefetch_* regs: read twice (shift + resolution chain bypass)
    // Resolution chains: NUM_PATHS when NUM_GROUPS==1, NUM_GROUPS when NUM_PATHS==1
    static constexpr u64 NUM_CHAINS = std::max(NUM_PATHS, NUM_GROUPS);
    static_loop<NT>([&]<u64 I>() {
      // NOTE: @prakhar @claude audit
      prefetch_tag_hit[I].fanout(
          hard<1 + NUM_CHAINS>{}); // shift + full_hits (1 per chain)
      prefetch_pred[I].fanout(
          hard<1 + NUM_CHAINS>{}); // shift + table_preds (1 per chain)
      // NOTE: @prakhar @claude audit
      prefetch_hyst[I].fanout(hard<2>{}); // shift + weak_mask
      // NOTE: @prakhar @claude audit
      prefetch_u[I].fanout(hard<2>{}); // shift + weak_mask
      if constexpr (USE_SEC_TAG)
        prefetch_sec[I].fanout(
            hard<1 + NUM_CHAINS>{}); // shift + full_hits (1 per chain)
      if constexpr (NUM_GROUPS > 1)
        prefetch_group_id[I].fanout(
            hard<1 + NUM_GROUPS>{}); // shift + group_id cmp (1 per group)
    });
    // NOTE: @prakhar @claude audit
    prefetch_fb.fanout(
        hard<1 + NUM_CHAINS>{}); // shift + table_preds (1 per chain)

    static_loop<NT>([&]<u64 I>() {
      current_tag[I] = prefetch_tag[I].fo1();
      current_tag_hit[I] = prefetch_tag_hit[I];
      current_pred[I] = prefetch_pred[I];
      if constexpr (USE_SEC_TAG)
        current_sec[I] = prefetch_sec[I];
      current_idx[I] = prefetch_idx[I].fo1();
      current_hyst[I] = prefetch_hyst[I];
      current_u[I] = prefetch_u[I];
      current_ctag[I] = prefetch_ctag[I].fo1();
      if constexpr (NUM_GROUPS > 1)
        current_group_id[I] = prefetch_group_id[I];
    });
    current_fb = prefetch_fb;
    current_fb_idx = prefetch_fb_idx.fo1();
    if constexpr (FB_RECONCILE)
      current_fb_hyst = prefetch_fb_hyst.fo1();
    current_pc = prefetch_pc.fo1();

#ifdef TAGE_MONITOR
    mon.shift_block_pc(num_branch);
#endif

    // Precompute secondary tag for next block
    // Use a local val to bypass the reg's transparent-latch penalty (119ps).
    // The reg is only needed for cross-cycle reads; same-cycle resolution
    // uses sec_tag_now directly.
    if constexpr (!USE_SEC_TAG)
      // NOTE: @prakhar @claude audit
      block_end_info.next_pc.fanout(hard<2>{}); // meta_idx + hist path_bits
    else
      // NOTE: @prakhar @claude audit
      // sec_tag_now (1) + meta_idx (1) + hist path_bits (1) [+ multi-path sel
      // (1)]
      block_end_info.next_pc.fanout(hard<3 + (NUM_PATHS > 1 ? 1 : 0)>{});
    auto sec_tag_now = [&]() {
      if constexpr (USE_SEC_TAG)
        return SecTagHashFn::template apply<SEC_TAG_BITS>(
            block_end_info.next_pc);
      else
        return hard<0>{};
    }();
    if constexpr (USE_SEC_TAG) {
      // NOTE: @prakhar @claude audit
      // Critical-path fanout: reg write (1) + full_hits comparisons (NT) = NT+1
      sec_tag_now.fanout(hard<NT + 1>{});
      curr_sec_tag = sec_tag_now; // store for next-cycle use
      // NOTE: @prakhar @claude audit
      // Alloc-path fanout: NT sec_ram writes (off critical path)
      // train_sec_tag = curr_sec_tag saved before overwrite = hash(predicted
      // block's PC). Old sec_tag_alloc was hash(successor), one stage off.
      train_sec_tag.fanout(hard<NT>{});
    }

#ifdef TAGE_MONITOR
    // sec_tag consistency check: for tag-hit entries, does stored sec_tag
    // match predict-time (curr_sec_tag from prev cycle) vs update-time
    // (sec_tag_now from this cycle)?
    if constexpr (USE_SEC_TAG) {
      u64 now_val = static_cast<u64>(sec_tag_now);
      for (u64 i = 0; i < NT; i++) {
        if (mon_tag_hit[i]) {
          auto &c = mon.cum;
          c.sec_tag_checks++;
          bool match_curr = mon_sec_hit[i]; // predict-time match
          bool match_now = (mon_saved_sec[i] == now_val);
          if (match_curr)
            c.sec_tag_match_curr++;
          if (match_now)
            c.sec_tag_match_now++;
          if (match_curr && match_now)
            c.sec_tag_match_both++;
          if (!match_curr && !match_now)
            c.sec_tag_match_neither++;
        }
      }
    }
#endif

    // ================================================================
    // Meta pipeline: shift FIRST, then write [0].
    // Fixes stale-value bug where meta_pipe[META_PIPE-1] reads the
    // new RAM value instead of the properly delayed old value.
    // META_BANKS separate RAMs: each bank has its own pipe.
    // ================================================================
    for (u64 mb = 0; mb < META_BANKS; mb++) {
      for (u64 i = META_PIPE - 1; i > 0; i--) {
        meta_pipe[mb][i] = meta_pipe[mb][i - 1].fo1();
        meta_idx_pipe[mb][i] = meta_idx_pipe[mb][i - 1].fo1();
      }
    }
    {
      auto meta_idx = val<META_IDX_BITS>{block_end_info.next_pc >> 2};
      // NOTE: @prakhar @claude audit
      meta_idx.fanout(hard<1 + META_BANKS>{}); // meta_idx_pipe[0] + K reads
      for (u64 mb = 0; mb < META_BANKS; mb++) {
        meta_pipe[mb][0] = meta_ctr[mb].read(meta_idx);
        meta_idx_pipe[mb][0] = meta_idx;
      }
    }
    // NOTE: @prakhar @claude audit
    // meta_pipe[mb][META_PIPE-1] read twice: meta_use_alt + old_meta in training
    arr<val<1>, META_BANKS> meta_use_alt = [&](u64 mb) -> val<1> {
      meta_pipe[mb][META_PIPE - 1].fanout(hard<2>{});
      return val<META_WIDTH, i64>{meta_pipe[mb][META_PIPE - 1]} >= hard<0>{};
    };
    // For multi-path configs, each meta_use_alt[mb] needs NUM_PATHS reads.
    // For multi-group (BPE1), each is read once by its group's resolve_chain.
    if constexpr (NUM_PATHS >= 2) {
      for (u64 mb = 0; mb < META_BANKS; mb++)
        meta_use_alt[mb].fanout(hard<NUM_PATHS>{});
    }

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
    // NOTE: @prakhar @claude audit
    branch_dir.fanout(hard<3>{}); // true_block + hist_input + actual_dir
#ifdef TAGE_MONITOR
    arr<val<1>, N> old_pred = [&](u64 i) -> val<1> { return val<1>{pred[i]}; };
#endif

#ifdef TAGE_MONITOR
    for (u64 i = 0; i < NT; i++) {
      u64 idx_val = static_cast<u64>(val<MAX_IDX_BITS>{current_idx[i]});
      if constexpr (USE_SEC_TAG) {
        // Use pre-shift values (mon_tag_hit/mon_sec_hit) to match the
        // pipeline stage being trained.  Post-shift current_sec[i] holds
        // the *next* block's prefetch data — comparing it against
        // sec_tag_now was off by one stage.
        mon.record_tag_lookup(i, mon_tag_hit[i], mon_sec_hit[i]);
        // Per-bank: which bank did this hit land in?
        mon.record_bank_tag_hit(i, idx_val, mon_tag_hit[i], mon_sec_hit[i]);
      } else {
        bool th = static_cast<u64>(current_tag_hit[i]);
        mon.record_tag_lookup(i, th, false);
        mon.record_bank_tag_hit(i, idx_val, th, false);
      }
      // Per-group tag lookup: attribute hit to the group stored in this entry
      if constexpr (NUM_GROUPS > 1) {
        bool th = USE_SEC_TAG ? mon_tag_hit[i] : (bool)static_cast<u64>(current_tag_hit[i]);
        bool sh = USE_SEC_TAG ? mon_sec_hit[i] : false;
        if (th) {
          u64 gid = static_cast<u64>(val<GROUP_BITS>{current_group_id[i]});
          mon.record_group_tag_lookup(gid, i, th, sh);
        }
      }
      mon.record_collision_check(i, idx_val,
          static_cast<u64>(current_tag_hit[i]), mon.current_block_pc);
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
      return w.fo1().concat(); // rvalue → fo1 per element (each used once)
    }();
    // weak_mask read once per resolve_chain call (NUM_CHAINS calls total)
    // NOTE: @prakhar @claude audit
    if constexpr (NUM_CHAINS >= 2)
      weak_mask.fanout(hard<NUM_CHAINS>{});

    auto resolve_chain = [&](arr<val<1>, NT> full_hits,
                             val<PRED_BITS> fb_pred,
                             u64 meta_bank = 0) {
      val<MATCH_BITS> match = concat(val<1>{1}, full_hits.fo1().concat());
      // NOTE: @prakhar @claude audit
      match.fanout(hard<3>{}); // ha(>>1) + one_hot + remainder(^match1)

      // has_alt: true iff any TAGE table matched. If so, provider is a
      // table and fallback (always bit 0) is below it as alt.
      // Depends only on match — no match1 dependency.
      val<1> ha = val<NT>{match} != hard<0>{};

      val<MATCH_BITS> match1 = match.one_hot();
      // NOTE: @prakhar @claude audit
      // make_array + XOR(remainder) + pw(>>1) + train_match1 save
      match1.fanout(hard<4>{});
      val<MATCH_BITS> remainder = match ^ match1;
      val<MATCH_BITS> match2 = remainder.fo1().one_hot();

      arr<val<PRED_BITS>, NT + 1> table_preds = [&](u64 i) -> val<PRED_BITS> {
        if (i < NT)
          return val<PRED_BITS>{prefetch_pred[i]};
        return fb_pred;
      };
      // NOTE: @prakhar @claude audit: pmask + amask = 2
      table_preds.fanout(hard<2>{});

      arr<val<1>, NT + 1> m1_bits = match1.make_array(val<1>{});
      // NOTE: @prakhar @claude audit: make_array(1) + return(1) = 2
      match2.fanout(hard<2>{});
      arr<val<1>, NT + 1> m2_bits = match2.make_array(val<1>{});

      arr<val<PRED_BITS>, NT + 1> pmask = [&](u64 i) {
        return m1_bits[i].fo1().replicate(hard<PRED_BITS>{}).concat() &
               table_preds[i];
      };
      arr<val<PRED_BITS>, NT + 1> amask = [&](u64 i) {
        return m2_bits[i].fo1().replicate(hard<PRED_BITS>{}).concat() &
               table_preds[i];
      };
      val<PRED_BITS> pp = pmask.fo1().fold_or();
      val<PRED_BITS> ap = amask.fo1().fold_or();
      // NOTE: @prakhar @claude audit: pp^ap(1) + return(1) = 2
      pp.fanout(hard<2>{});
      // NOTE: @prakhar @claude audit: pp^ap(1) + return(1) = 2
      ap.fanout(hard<2>{});

      // Provider weakness: AND match1's table bits with precomputed
      // weak_mask, then check nonzero. Single gate after match1 instead
      // of fold_or over NT+1 terms.
      auto pw_mask = [&]() -> val<NT> {
        if constexpr (NUM_CHAINS == 1)
          return weak_mask.fo1();
        else
          return weak_mask;
      }();
      val<1> pw = (val<NT>{match1} & pw_mask.fo1()) != hard<0>{};
      // NOTE: @prakhar @claude audit: ua_comp(1) + return(1) = 2
      pw.fanout(hard<2>{});
      // meta_use_alt: select the appropriate bank's counter.
      // For single-bank configs (META_BANKS==1), bank is always 0.
      auto mua = [&]() -> val<1> {
        if constexpr (META_BANKS == 1 && NUM_CHAINS == 1)
          return meta_use_alt[0].fo1();
        else
          return val<1>{meta_use_alt[meta_bank]};
      }();
      val<1> ua = pw & mua.fo1() & ha.fo1();
      // ua: no internal reads, returned via fo1; outside binding gets fanout(N)

      // Final select is NOT done here — it's fused with the per-branch
      // scatter outside, using 1-bit selects to avoid PRED_BITS-wide
      // replication of ua (saves ~45ps ctrl-to-out delay).
      auto pp_xor_ap = pp ^ ap;
      val<1> ad = pp_xor_ap.fo1() != hard<0>{};

      // fanout'd vars use lvalue copy (last credit); non-fanout'd use fo1
      return std::tuple{val<MATCH_BITS>{match1},
                        val<MATCH_BITS>{match2},
                        val<PRED_BITS>{pp},
                        val<PRED_BITS>{ap},
                        val<1>{pw},
                        ad.fo1(),
                        ua.fo1()};
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
        auto [m1_0, m2_0, pp0, ap0, pw0, ad0, ua0] =
            resolve_chain(fh0.fo1(), val<PRED_BITS>{prefetch_fb});

        // Chain 1: stored sec_tag == 1 (compile-time constant)
        arr<val<1>, NT> fh1 = [&](u64 i) {
          return val<1>{prefetch_tag_hit[i]} &
                 (val<SEC_TAG_BITS>{prefetch_sec[i]} == hard<1>{});
        };
        auto [m1_1, m2_1, pp1, ap1, pw1, ad1, ua1] =
            resolve_chain(fh1.fo1(), val<PRED_BITS>{prefetch_fb});

        // Fields muxed through NP select: m1, m2, pp, ap, pw, ad, ua.
        static constexpr u64 MUX_FIELDS = 7;

        if constexpr (NUM_PATHS == 2) {
          // 2-to-1 mux: derive sel directly from next_pc (bypass reg chain)
          val<1> sec_sel = val<1>{SecTagHashFn::template apply<SEC_TAG_BITS>(
              block_end_info.next_pc)};
          // NOTE: @prakhar @claude audit
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
          auto [m1_2, m2_2, pp2, ap2, pw2, ad2, ua2] =
              resolve_chain(fh2.fo1(), val<PRED_BITS>{prefetch_fb});

          arr<val<1>, NT> fh3 = [&](u64 i) {
            return val<1>{prefetch_tag_hit[i]} &
                   (val<SEC_TAG_BITS>{prefetch_sec[i]} == hard<3>{});
          };
          auto [m1_3, m2_3, pp3, ap3, pw3, ad3, ua3] =
              resolve_chain(fh3.fo1(), val<PRED_BITS>{prefetch_fb});

          // 4-to-1 mux tree: derive sel directly from next_pc (bypass reg)
          val<SEC_TAG_BITS> sec_idx =
              SecTagHashFn::template apply<SEC_TAG_BITS>(
                  block_end_info.next_pc);
          sec_idx.fanout(hard<2>{}); // lo + hi
          val<1> lo = sec_idx & hard<1>{};
          val<1> hi = (sec_idx >> 1) & hard<1>{};
          // NOTE: @prakhar @claude audit
          lo.fanout(hard<MUX_FIELDS * 2>{}); // each field needs lo for 2 pairs
          // NOTE: @prakhar @claude audit
          hi.fanout(hard<MUX_FIELDS>{}); // each field needs hi for final

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
        // Single-path: sec-tag policy controls per-table enforcement.
        // sec_tag_now fanout already declared above (covers all reads)

        // Precompute per-table sec-tag enforcement (compile-time).
        static constexpr auto sec_enforce = []() {
          std::array<bool, NT> a{};
          // Use a helper to expand the parameter pack
          auto fill = [&]<u64... Is>(std::index_sequence<Is...>) {
            ((a[Is] = SecTagPolicy::template apply<Is, NT>()), ...);
          };
          fill(std::make_index_sequence<NT>{});
          return a;
        }();
        // Count how many tables enforce (for fanout on runtime gate)
        static constexpr u64 SEC_ENFORCE_COUNT = []() {
          u64 n = 0;
          for (u64 i = 0; i < NT; i++) if (sec_enforce[i]) n++;
          return n;
        }();

        // Runtime-gated policies: single shared gate signal.
        [[maybe_unused]] val<1> sec_rt_gate = [&]() -> val<1> {
          if constexpr (USE_SEC_TAG && SecTagPolicy::RUNTIME &&
                        SEC_ENFORCE_COUNT > 0) {
            return SecTagPolicy::template gate<0, NT>(
                val<ACC_WIDTH>{acc_ctr}, val<ALLOC_WIDTH>{alloc_ctr},
                val<BENEFIT_REG_WIDTH>{benefit_ctr});
          } else {
            return val<1>{hard<1>{}};
          }
        }();
        if constexpr (USE_SEC_TAG && SecTagPolicy::RUNTIME &&
                      SEC_ENFORCE_COUNT > 1)
          sec_rt_gate.fanout(hard<SEC_ENFORCE_COUNT * NUM_GROUPS>{});

        // Helper: build per-table sec-tag-filtered hit, shared by all groups
        auto sec_filtered_hit = [&](u64 i) -> val<1> {
          if constexpr (!USE_SEC_TAG || SEC_ENFORCE_COUNT == 0) {
            return val<1>{prefetch_tag_hit[i]};
          } else {
            val<1> sec_match =
                val<SEC_TAG_BITS>{prefetch_sec[i]} ==
                val<SEC_TAG_BITS>{sec_tag_now};
            if constexpr (!SecTagPolicy::RUNTIME) {
              if (sec_enforce[i])
                return val<1>{prefetch_tag_hit[i]} & sec_match;
              else
                return val<1>{prefetch_tag_hit[i]};
            } else {
              if (sec_enforce[i])
                return val<1>{prefetch_tag_hit[i]} &
                       (sec_match | ~sec_rt_gate);
              else
                return val<1>{prefetch_tag_hit[i]};
            }
          }
        };

        if constexpr (NUM_GROUPS == 1) {
          // Original single-chain path
          arr<val<1>, NT> full_hits = [&](u64 i) {
            return sec_filtered_hit(i);
          };
#ifdef DEBUG_PRINT
          std::cerr << "--- sec-tag path ---\n";
          sec_tag_now.print("  sec_tag_now=", "\n", true, std::cerr);
          if constexpr (USE_SEC_TAG && SecTagPolicy::RUNTIME)
            sec_rt_gate.print("  sec_rt_gate=", "\n", true, std::cerr);
          for (u64 i = 0; i < NT; i++) {
            std::cerr << "  sec_enforce[" << i << "]=" << sec_enforce[i];
            full_hits[i].print("  full_hit=", "\n", true, std::cerr);
          }
#endif
          return resolve_chain(full_hits.fo1(), val<PRED_BITS>{prefetch_fb});
        }
        // NUM_GROUPS > 1: return dummy values (per-group block below does real work)
        return std::tuple{val<MATCH_BITS>{hard<0>{}}, val<MATCH_BITS>{hard<0>{}},
                          val<PRED_BITS>{hard<0>{}}, val<PRED_BITS>{hard<0>{}},
                          val<1>{hard<0>{}}, val<1>{hard<0>{}}, val<1>{hard<0>{}}};
      }
    }();

    // NUM_GROUPS>1: pre-read ALL old train reg values BEFORE per-group resolution
    // overwrites them. Training and any_provider_wrong use these pre-read values.
    // (For NUM_GROUPS==1, reads happen after resolution — see below.)
    u64 train_group = (NUM_GROUPS > 1)
        ? BRANCH_TO_GROUP[num_branch > 0 ? num_branch - 1 : 0]
        : 0;
    // Pre-read all groups' provider predictions (needed by any_provider_wrong)
    arr<val<PRED_BITS>, NUM_GROUPS> old_provider_pred = [&](u64 g) -> val<PRED_BITS> {
      if constexpr (NUM_GROUPS > 1)
        return val<PRED_BITS>{train_provider_pred[g].fo1()};
      else
        return val<PRED_BITS>{hard<0>{}}; // placeholder
    };
    val<MATCH_BITS> pre_read_match1 = [&]() {
      if constexpr (NUM_GROUPS > 1)
        return train_match1[train_group].fo1();
      else
        return val<MATCH_BITS>{hard<0>{}}; // placeholder, overwritten below
    }();
    // Pre-read pw/ad for all groups (META_BANKS>1 needs per-group training).
    // For META_BANKS==1, only train_group is needed (backward compat).
    arr<val<1>, NUM_GROUPS> pre_read_pw = [&](u64 g) -> val<1> {
      if constexpr (NUM_GROUPS > 1) {
        if constexpr (META_BANKS > 1)
          return train_provider_weak[g].fo1();
        else
          return (g == u64(train_group)) ? train_provider_weak[g].fo1() : val<1>{hard<0>{}};
      } else
        return val<1>{hard<0>{}}; // placeholder
    };
    arr<val<1>, NUM_GROUPS> pre_read_ad = [&](u64 g) -> val<1> {
      if constexpr (NUM_GROUPS > 1) {
        if constexpr (META_BANKS > 1)
          return train_altdiff[g].fo1();
        else
          return (g == u64(train_group)) ? train_altdiff[g].fo1() : val<1>{hard<0>{}};
      } else
        return val<1>{hard<0>{}}; // placeholder
    };

    if constexpr (NUM_GROUPS > 1) {
      // ---- Per-group resolution chains ----
      // Each group gets its own full_hits (htag_hit & group_id match & sec_tag)
      // and its own resolve_chain call. Results stored in per-group arrays.

      // Precompute sec-tag enforcement (same as single-chain path above)
      static constexpr auto sec_enforce = []() {
        std::array<bool, NT> a{};
        auto fill = [&]<u64... Is>(std::index_sequence<Is...>) {
          ((a[Is] = SecTagPolicy::template apply<Is, NT>()), ...);
        };
        fill(std::make_index_sequence<NT>{});
        return a;
      }();
      static constexpr u64 SEC_ENFORCE_COUNT = []() {
        u64 n = 0;
        for (u64 i = 0; i < NT; i++) if (sec_enforce[i]) n++;
        return n;
      }();

      [[maybe_unused]] val<1> sec_rt_gate = [&]() -> val<1> {
        if constexpr (USE_SEC_TAG && SecTagPolicy::RUNTIME &&
                      SEC_ENFORCE_COUNT > 0) {
          return SecTagPolicy::template gate<0, NT>(
              val<ACC_WIDTH>{acc_ctr}, val<ALLOC_WIDTH>{alloc_ctr},
              val<BENEFIT_REG_WIDTH>{benefit_ctr});
        } else {
          return val<1>{hard<1>{}};
        }
      }();
      if constexpr (USE_SEC_TAG && SecTagPolicy::RUNTIME &&
                    SEC_ENFORCE_COUNT > 1)
        sec_rt_gate.fanout(hard<SEC_ENFORCE_COUNT * NUM_GROUPS>{});

      // Extract per-group fb bits from N-wide fallback
      auto extract_group_fb = [&]<u64 G>() -> val<PRED_BITS> {
        // Gather BR_P_ENTRY bits for group G from the N-wide fb
        arr<val<CTR_WIDTH>, BR_P_ENTRY> fb_bits = [&](u64 j) {
          u64 branch_idx = (j < GROUP_SIZE[G])
              ? [&]() {
                  // Find the j-th branch in group G
                  u64 count = 0;
                  for (u64 b = 0; b < N; b++) {
                    if (BRANCH_TO_GROUP[b] == G) {
                      if (count == j) return b;
                      count++;
                    }
                  }
                  return u64(0); // unreachable
                }()
              : u64(0);
          return val<CTR_WIDTH>{prefetch_fb >> (branch_idx * CTR_WIDTH)};
        };
        return fb_bits.concat();
      };

      static_loop<NUM_GROUPS>([&]<u64 G>() {
        // Per-group full_hits: htag_hit & (group_id == G) & sec_tag
        arr<val<1>, NT> group_hits = [&](u64 i) -> val<1> {
          val<1> htag_ok = val<1>{prefetch_tag_hit[i]};
          val<1> group_ok =
              val<GROUP_BITS>{prefetch_group_id[i]} == hard<G>{};
          val<1> base_hit = htag_ok & group_ok;
          // Apply sec-tag filtering
          if constexpr (!USE_SEC_TAG || SEC_ENFORCE_COUNT == 0) {
            return base_hit;
          } else {
            val<1> sec_match =
                val<SEC_TAG_BITS>{prefetch_sec[i]} ==
                val<SEC_TAG_BITS>{sec_tag_now};
            if constexpr (!SecTagPolicy::RUNTIME) {
              if (sec_enforce[i])
                return base_hit & sec_match;
              else
                return base_hit;
            } else {
              if (sec_enforce[i])
                return base_hit & (sec_match | ~sec_rt_gate);
              else
                return base_hit;
            }
          }
        };

        val<PRED_BITS> fb_g = extract_group_fb.template operator()<G>();
        auto [m1, m2, pp, ap, pw, ad, ua] =
            resolve_chain(group_hits.fo1(), fb_g, GROUP_TO_META_BANK[G]);

        // Scatter: assign pred[I] for branches belonging to this group
        // NOTE: @prakhar @claude audit
        pp.fanout(hard<2>{}); // make_array + train save
        if constexpr (GROUP_SIZE[G] >= 2)
          ua.fanout(hard<GROUP_SIZE[G]>{}); // one per branch in group
        arr<val<1>, PRED_BITS> pp_bits = pp.make_array(val<1>{});
        arr<val<1>, PRED_BITS> ap_bits = ap.fo1().make_array(val<1>{});
        static_loop<N>([&]<u64 I>() {
          if constexpr (BRANCH_TO_GROUP[I] == G) {
            static constexpr u64 POS = BRANCH_IN_GROUP[I];
            pred[I] = select(ua, ap_bits[POS].fo1(), pp_bits[POS].fo1());
          }
        });

        // Save per-group resolution → train regs (OK because training pre-reads old values above)
        train_match1[G] = m1.fo1();
        train_provider_pred[G] = pp;
        train_provider_weak[G] = pw.fo1();
        train_altdiff[G] = ad.fo1();

#ifdef TAGE_MONITOR
        // Record per-group provider info for monitor
        {
          u64 m1v = static_cast<u64>(m1);
          u64 m2v = static_cast<u64>(m2);
          bool meta_overrode = static_cast<u64>(pw) &&
              ((m2v & ((1ULL << NT) - 1)) != 0);
          bool meta_chose = static_cast<u64>(ua);
          u64 prov = decltype(mon)::decode_provider(m1v);
          bool has_tage = (prov < NT);

          // any_tag_hit for this group: htag matched AND group_id==G
          // (ignoring sec_tag — that's the counterfactual question)
          bool any_tag_hit_g = false;
          u64 tag_only_provider_g = NT; // highest matching table ignoring sec
          for (u64 i = 0; i < NT; i++) {
            if (mon_tag_hit[i]) {
              // Check if this entry's group_id matches G
              u64 gid = static_cast<u64>(val<GROUP_BITS>{current_group_id[i]});
              if (gid == G) {
                if (!any_tag_hit_g) tag_only_provider_g = i;
                any_tag_hit_g = true;
              }
            }
          }

          // tag_only_pred: what would tag-only provider (ignoring sec) predict?
          // Use software-only reads to avoid creating hardware
          u64 tag_only_pred_val = 0;
          if (tag_only_provider_g < NT) {
            tag_only_pred_val = static_cast<u64>(prefetch_pred[tag_only_provider_g]);
          } else {
            // No tag match at all — use fb prediction (software read)
            tag_only_pred_val = static_cast<u64>(prefetch_fb);
          }

          // Attribute to each branch in this group
          static_loop<N>([&]<u64 I>() {
            if constexpr (BRANCH_TO_GROUP[I] == G) {
              static constexpr u64 POS = BRANCH_IN_GROUP[I];
              bool pred_taken = static_cast<u64>(pred[I]);
              bool tag_only_taken = (tag_only_pred_val >> POS) & 1;
              mon.record_prediction(I, m1v, m2v, meta_overrode, meta_chose,
                                    pred_taken, any_tag_hit_g, has_tage,
                                    tag_only_taken, GROUP_TO_META_BANK[G]);
            }
          });

          // Record provider entry lifetime + bank
          if (has_tage) {
            u64 prov_index = static_cast<u64>(val<MAX_IDX_BITS>{current_idx[prov]});
            // Count correct predictions for this group's branches
            u64 nc = 0, nb = 0;
            static_loop<N>([&]<u64 I>() {
              if constexpr (BRANCH_TO_GROUP[I] == G) {
                nb++;
                if (static_cast<u64>(pred[I]) == static_cast<u64>(branch_dir[I]))
                  nc++;
              }
            });
            mon.record_provider_entry(prov, prov_index, nb, nc);
            bool all_correct = (nc == nb);
            mon.record_bank_provider(prov, prov_index, all_correct);
          }
        }
#endif
      });
    }

    if constexpr (NUM_GROUPS == 1) {

    // Per-branch 1-bit scatter: split pp/ap into per-bit arrays to reduce
    // fanout on the wide values (fanout(3) + N×fo1 vs fanout(N+2) on each).
    // NOTE: @prakhar @claude audit: make_array(1) + train_provider_pred save(1)
    // = 2
    provider_pred.fanout(hard<2>{});
    // NOTE: @prakhar @claude audit: N reads in static_loop scatter
    use_alt.fanout(hard<N>{});
    arr<val<1>, PRED_BITS> pp_bits = provider_pred.make_array(val<1>{});
    arr<val<1>, PRED_BITS> ap_bits = alt_pred.fo1().make_array(val<1>{});
    static_loop<N>([&]<u64 I>() {
      pred[I] = select(use_alt, ap_bits[I].fo1(), pp_bits[I].fo1());
    });

#ifdef DEBUG_PRINT
    std::cerr << "--- resolve_chain output ---\n";
    block_end_info.next_pc.print("  blk_end_info.next_pc= ", "\n", true,
                                 std::cerr);
    match1.print("  match1=", "\n", true, std::cerr);
    match2.print("  match2=", "\n", true, std::cerr);
    provider_pred.print("  provider_pred=", "\n", true, std::cerr);
    alt_pred.print("  alt_pred=", "\n", true, std::cerr);
    provider_weak.print("  provider_weak=", "\n", true, std::cerr);
    altdiff.print("  altdiff=", "\n", true, std::cerr);
    use_alt.print("  use_alt=", "\n", true, std::cerr);
    meta_use_alt.print("  meta_use_alt=", "\n", true, std::cerr);
    weak_mask.print("  weak_mask=", "\n", true, std::cerr);
    for (u64 i = 0; i < N; i++)
      pred[i].print(("  pred[" + std::to_string(i) + "]=").c_str(), "\n", true,
                    std::cerr);
#endif

#ifdef TAGE_MONITOR
    // Provider source: any_tag_hit = at least one table matched on primary tag
    // has_tage_provider = provider is a TAGE table (not fallback)
    bool mon_any_tag_hit = false;
    u64 mon_tag_only_provider = NT; // fallback if no tag match
    for (u64 i = 0; i < NT; i++) {
      if (mon_tag_hit[i]) {
        if (!mon_any_tag_hit) mon_tag_only_provider = i;
        mon_any_tag_hit = true;
      }
    }
    u64 m1v_check = static_cast<u64>(match1);
    bool mon_has_tage = (decltype(mon)::decode_provider(m1v_check) < NT);
    // Counterfactual: what would the tag-only provider have predicted?
    // (ignoring sec-tag filter — for branches where sec-tag rejected all)
    u64 mon_tag_only_pred = (mon_tag_only_provider < NT)
        ? static_cast<u64>(prefetch_pred[mon_tag_only_provider])
        : static_cast<u64>(prefetch_fb);
    for (u64 r = 0; r < num_branch; r++) {
      bool meta_overrode =
          static_cast<u64>(provider_weak) &&
          ((static_cast<u64>(match2) & ((1ULL << NT) - 1)) != 0);
      bool meta_chose = static_cast<u64>(use_alt);
      bool pred_taken = static_cast<u64>(pred[r]);
      bool tag_only_taken = (mon_tag_only_pred >> r) & 1;
      mon.record_prediction(r, static_cast<u64>(match1),
                            static_cast<u64>(match2), meta_overrode, meta_chose,
                            pred_taken, mon_any_tag_hit, mon_has_tage,
                            tag_only_taken);
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
        // Per-bank: which bank provided?
        bool all_correct = (nc == num_branch);
        mon.record_bank_provider(prov, prov_index, all_correct);
      }
    }
#endif

    } // end if constexpr (NUM_GROUPS == 1)

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
    // resolution.
    // Declare fanout on training regs before any reads:
    //   train_hyst[I]:     2 reads (t_hyst_weak_arr + training loop line ~1188)
    //   train_u[I]:        2 reads (u_zero lambda + decay line ~1254)
    //   train_tag_hit[I]:  reads: (sibling lambda if active for I) + (decay if
    //   enabled) train_sec_hit[I]:  reads: (sibling lambda if active for I) +
    //   (decay if enabled) train_idx[I]:      5 reads (pred/hyst/tag/sec/u RAM
    //   writes) single-read regs
    //   (match1/provider_pred/provider_weak/altdiff/pred/ctag) use fo1() at the
    //   point of use below.
    static_loop<NT>([&]<u64 I>() {
      train_hyst[I].fanout(hard<2>{}); // NOTE: @prakhar @claude audit
      // NOTE: @prakhar @claude audit
      // u_zero(1) + base old_u(1) [+ decay old_u(1)]
      train_u[I].fanout(hard<2 + (DECAY_ENABLE ? 1 : 0)>{});
      // sibling skip reads train_tag_hit/train_sec_hit for ALL tables
      // when SIBLING_POLICY==ALL (runtime lambda iterates all indices).
      // SIBLING_TABLE_FLOOR masks the result post-hoc, not the read.
      static constexpr bool SIBLING_ACTIVE =
          USE_SEC_TAG && SIBLING_POLICY == SiblingPolicy::ALL;
      static constexpr u64 TAG_HIT_READS =
          (SIBLING_ACTIVE ? 1 : 0) + (DECAY_ENABLE ? 1 : 0);
      static constexpr u64 SEC_HIT_READS =
          (SIBLING_ACTIVE ? 1 : 0) + (DECAY_ENABLE ? 1 : 0);
      // NOTE: @prakhar @claude audit
      if constexpr (TAG_HIT_READS > 1)
        train_tag_hit[I].fanout(hard<TAG_HIT_READS>{});
      // NOTE: @prakhar @claude audit
      if constexpr (SEC_HIT_READS > 1)
        train_sec_hit[I].fanout(hard<SEC_HIT_READS>{});
      // TAG/SEC_HIT_READS==1 → use fo1() at point of use; ==0 → unreferenced
      // train_group_id: sibling detection only (entry_actual_dir no longer needs it).
      // When SIBLING_ACTIVE: 1 read in candallocmask. Otherwise unused in HW.
      // (Monitor peeks via CHEATING_MODE are free.)
      train_idx[I].fanout(hard<5
#ifdef TAGE_MONITOR
                                  + 2
#endif
                                  >{}); // NOTE: @prakhar @claude audit
    });

    // Read old train reg values for training.
    // NUM_GROUPS>1: use pre-read values (read before per-group loop overwrote train regs).
    // NUM_GROUPS==1: read now (resolution hasn't overwritten train regs yet).
    val<MATCH_BITS> t_match1 = [&]() {
      if constexpr (NUM_GROUPS > 1)
        return pre_read_match1;
      else
        return train_match1[train_group].fo1();
    }();
    t_match1.fanout(
        hard<2 + (FARALLOC_DIST > 0 ? 1 : 0)>{}); // make_array(1) +
                                                  // alloc_base(1) [+ FARALLOC
                                                  // provider_bits(1)]
    arr<val<1>, NT + 1> t_m1 = t_match1.make_array(val<1>{});
    // t_m1[I<NT]: 4 reads — t_hyst_weak_arr(1) + do_pred_update(1) +
    // do_hyst_update(1) + base_u_write(1) t_m1[NT]:   1 read normally, 2 if
    // FB_RECONCILE — declared at use site
    static_loop<NT>([&]<u64 I>() {
      t_m1[I].fanout(hard<4>{}); // NOTE: @prakhar @claude audit
    });
    if constexpr (FB_RECONCILE || FB_BIM_HYST)
      t_m1[NT].fanout(hard<1 + (FB_RECONCILE ? 1 : 0) + (FB_BIM_HYST ? 1 : 0)>{}); // NOTE: @prakhar @claude audit
    val<PRED_BITS> t_pp = [&]() -> val<PRED_BITS> {
      if constexpr (NUM_GROUPS > 1)
        return val<PRED_BITS>{old_provider_pred[train_group]};
      else
        return train_provider_pred[train_group].fo1();
    }();
    val<1> t_pw = [&]() -> val<1> {
      if constexpr (NUM_GROUPS > 1)
        return pre_read_pw[train_group];
      else
        return train_provider_weak[train_group].fo1();
    }();
    val<1> t_ad = [&]() -> val<1> {
      if constexpr (NUM_GROUPS > 1)
        return pre_read_ad[train_group];
      else
        return train_altdiff[train_group].fo1();
    }();

    // Save current resolution → train regs (for NEXT cycle's training)
    // NUM_GROUPS==1: save here (after reads above).
    // NUM_GROUPS>1: already saved in per-group loop (safe because pre-read above).
    if constexpr (NUM_GROUPS == 1) {
      train_match1[0] = match1.fo1();
      train_provider_pred[0] = provider_pred;
      train_provider_weak[0] = provider_weak.fo1();
      train_altdiff[0] = altdiff.fo1();
    }

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
    val<1> t_phw = t_hyst_weak_arr.fo1().fold_or();

    // Read train_valid BEFORE setting it to 1 (regs may be immediate-write)
    val<1> do_train = train_valid.fo1();
    train_valid = 1;

    // ---- No conditional branches: update history, skip training ----
    if (num_branch == 0) {
      val<PATHBITS> path_bits =
          val<PATHBITS>{block_end_info.next_pc.fo1() >> 2};
      // NOTE: @prakhar @claude audit
      path_bits.fanout(hard<NT * 2 + 1 + (USE_GSHARE ? 1 : 0)>{});
      gh.template fanout_per_bit<GH_FANOUT>();
      static_loop<NT>([&]<u64 I>() {
        auto &t = table<I>();
        t.fold_idx.apply_update(t.fold_idx.compute_update(
            gh, hard<TableCfg::HIST_LEN[I]>{}, path_bits));
        t.fold_tag.apply_update(t.fold_tag.compute_update(
            gh, hard<TableCfg::HIST_LEN[I]>{}, path_bits));
      });
      if constexpr (USE_GSHARE)
        fb_fold.apply_update(
            fb_fold.compute_update(gh, hard<GS_HIST>{}, path_bits));
      gh.update(path_bits);
      last_condbr_dir = 0;
      true_block = 1;
      return;
    }

    // ---- Step 1: Correctness ----
    val<1> &mispredict = block_end_info.is_mispredict;
    // NOTE: @prakhar @claude audit: extra_cycle(1) + fb_gate(1) + true_block(1)
    //   + acc_ctr(1 if DECAY||EPOCH) + alloc_trigger(1 if MISPREDICT)
    //   [+ print(1) if DEBUG_PRINT]
    static constexpr u64 MISP_READS =
        3 + ((DECAY_ENABLE || EPOCH_ENABLE) ? 1 : 0) +
        (AllocCfg::ALLOC_TRIGGER == AllocTrigger::MISPREDICT ? 1 : 0);
    mispredict.fanout(hard<MISP_READS>{});
    // Read fb bimodal hyst BEFORE extra_cycle (cycle 1); write happens AFTER (cycle 2)
    // Declare train_fb_idx fanout early since bh read precedes other uses.
    // Reads: bh_read(1) + fb_write(1) + bh_write(1) [+ reconcile(2) if FB_RECONCILE]
    if constexpr (FB_BIM_HYST) {
      train_fb_idx.fanout(hard<3 + (FB_RECONCILE ? 2 : 0)>{});
      fb_bh_read = execute_if(mispredict, [&]() {
        return fb_bim_hyst.read(val<FB_BH_IDX_BITS>{train_fb_idx});
      });
    }
    need_extra_cycle(mispredict);
#ifdef DEBUG_PRINT
    std::cerr << "--- update_cycle ---\n";
    mispredict.print("  mispredict=", "\n", true, std::cerr);
#endif
    // NOTE: @prakhar @claude audit
    do_train.fanout(
        hard<4 * NT + 2 + (FB_RECONCILE ? 2 : 0) + (FB_BIM_HYST ? 1 : 0)>{}); // fb + meta + 4 per table
                                                      // [+ reconcile] [+ bim_hyst]

#ifdef TAGE_MONITOR
    mon.record_block(static_cast<u64>(val<LOGLINEINST>{block_entry}),
                     block_size, num_branch, static_cast<u64>(mispredict));
    for (u64 r = 0; r < num_branch; r++) {
      bool actual = static_cast<u64>(branch_dir[r]);
      bool is_misp = r == num_branch - 1 && static_cast<u64>(mispredict);
      mon.record_outcome(r, actual, is_misp);
      // Per-group: attribute this branch's outcome to its group's chain
      if constexpr (NUM_GROUPS > 1) {
        u64 g = BRANCH_TO_GROUP[r];
        bool pred_correct = (mon.shadow_pred[r] == actual);
        u64 prov = mon.shadow_provider[r];
        mon.record_group_outcome(g, pred_correct, is_misp, prov);
      }
    }
    if (!static_cast<u64>(do_train))
      mon.record_train_skip();
#endif
#ifdef TAGE_MONITOR
    // TODO: No idea why this even is here
    // arr<val<1>, N> badpred = [&](u64 i) -> val<1> {
    //   return old_pred[i].fo1() ^ val<1>{branch_dir[i]};
    // };
#endif

    // ---- Compute training signals ----
    // Per-group actual direction: each group's branches concatenated into
    // PRED_BITS = BR_P_ENTRY * CTR_WIDTH bits.
    // NUM_GROUPS==1: single val<PRED_BITS> from all N branches (same as before).
    // NUM_GROUPS>1: arr of NUM_GROUPS vals, each BR_P_ENTRY*CTR_WIDTH bits.
    arr<val<PRED_BITS>, NUM_GROUPS> actual_dir_g = [&](u64 g) -> val<PRED_BITS> {
      arr<val<CTR_WIDTH>, BR_P_ENTRY> bits = [&](u64 j) -> val<CTR_WIDTH> {
        // Find the j-th branch in group g
        u64 branch_idx = 0;
        u64 count = 0;
        for (u64 b = 0; b < N; b++) {
          if (BRANCH_TO_GROUP[b] == g) {
            if (count == j) { branch_idx = b; break; }
            count++;
          }
        }
        return val<CTR_WIDTH>{branch_dir[branch_idx]};
      };
      return bits.concat();
    };
    // actual_dir_g[0] used directly (no fo1 alias — multiple reads needed)

    // Provider wrong on any branch? (uses piped provider_pred)
    // NUM_GROUPS==1: single comparison. NUM_GROUPS>1: per-group, fold_or.
    // NOTE: @prakhar @claude audit
    if constexpr (NUM_GROUPS == 1) {
      if constexpr (FB_RECONCILE)
        t_pp.fanout(hard<2>{}); // any_provider_wrong + fb_tage_disagree
      // NOTE: @prakhar @claude audit
      // actual_dir_g[0] reads: any_provider_wrong(1) + actual_dir_full(1) + entry_actual_dir×NT
      actual_dir_g.fanout(hard<2 + NT>{});
    } else {
      // Per-group: each actual_dir_g[G] read for any_provider_wrong + table_wrong selects
      actual_dir_g.fanout(hard<2 + NT>{}); // any_prov_wrong(1) + fb_changed(1) + table_wrong×NT
    }
    // Per-group wrong signals (kept for META_BANKS>1 per-bank training)
    arr<val<1>, NUM_GROUPS> per_group_wrong = [&](u64 g) -> val<1> {
      if constexpr (NUM_GROUPS == 1) {
        if constexpr (FB_RECONCILE)
          return (t_pp ^ actual_dir_g[0]) != hard<0>{};
        else
          return (t_pp.fo1() ^ actual_dir_g[0]) != hard<0>{};
      } else {
        return (val<PRED_BITS>{old_provider_pred[g]} ^
                actual_dir_g[g]) != hard<0>{};
      }
    };
    val<1> any_provider_wrong = [&]() {
      if constexpr (NUM_GROUPS == 1)
        return per_group_wrong[0];
      else
        return per_group_wrong.fo1().fold_or();
    }();
    // NOTE: @prakhar @claude audit: meta_update(1) [+ alloc_trigger(1) if
    // TAGE_WRONG]
    // NOTE: @prakhar @claude audit
    // META_BANKS==1: any_provider_wrong used by meta_update(1) + alloc_trigger(1)
    // META_BANKS>1:  meta uses per_group_wrong; any_provider_wrong only for alloc
    if constexpr (AllocCfg::ALLOC_TRIGGER == AllocTrigger::TAGE_WRONG) {
      if constexpr (META_BANKS == 1)
        any_provider_wrong.fanout(hard<2>{});
      // META_BANKS>1: no extra fanout needed — fo1 suffices for alloc_trigger
    }
    // NOTE: @prakhar @claude audit: meta_gate(1) + meta_update_dir(1) = 2
    t_pw.fanout(hard<2>{});
    // NOTE: @prakhar @claude audit: do_pred_update per table = NT
    t_phw.fanout(hard<NT>{});
    // t_m1[NT] fanout declared above (fo1 or fanout per
    // FB_RECONCILE/FB_BIM_HYST)
    // NOTE: @prakhar @claude audit: u_correct×NT(NT) + meta_gate(1)
    //   [+ u_wrong×NT(NT) if U_MISP != UNTOUCHED]
    t_ad.fanout(hard<NT + (U_MISP != UMispPolicy::UNTOUCHED ? NT : 0) + 1>{});

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
    val<NT> triggermask = alloc_trigger.fo1().replicate(hard<NT>{}).concat();

    // 3b. Allocation action (probabilistic gating)
    val<8> alloc_rng = val<8>{static_cast<u64>(std::rand()) & 0xFF};
    val<NT> gated_triggermask = [&]() -> val<NT> {
      if constexpr (AllocCfg::ALLOC_ACTION == AllocAction::STANDARD) {
        return triggermask.fo1();
      } else if constexpr (AllocCfg::ALLOC_ACTION == AllocAction::FILTERED) {
        static_assert(ACC_WIDTH > 0, "FILTERED requires ACC_WIDTH > 0");
        val<1> allow = (val<ACC_WIDTH>{acc_ctr} >= val<ACC_WIDTH>{alloc_rng});
        return triggermask.fo1() & allow.fo1().replicate(hard<NT>{}).concat();
      } else {
        static_assert(AllocCfg::ALLOC_ACTION == AllocAction::THROTTLED);
        static_assert(ALLOC_WIDTH > 0, "THROTTLED requires ALLOC_WIDTH > 0");
        val<1> allow =
            (val<ALLOC_WIDTH>{alloc_ctr} >= val<ALLOC_WIDTH>{alloc_rng});
        return triggermask.fo1() & allow.fo1().replicate(hard<NT>{}).concat();
      }
    }();

    // 3c. Candidate mask: tables above provider with u=0
    val<NT> alloc_base = val<NT>{t_match1 - 1};
    arr<val<1>, NT> u_zero = [&](u64 i) -> val<1> {
      return val<U_WIDTH>{train_u[i]} == hard<0>{};
    };
    val<NT> notumask = u_zero.fo1().concat();
    // NOTE: @prakhar @claude audit
    static constexpr bool UCLEAR_ACTIVE = U_CLEAR != UClearPolicy::DISABLED;
    val<NT> postmask = alloc_base.fo1() & gated_triggermask.fo1();
    if constexpr (UCLEAR_ACTIVE) {
      // NOTE: @prakhar @claude audit
      // postmask: candallocmask(1) + uclearmask(1)
      postmask.fanout(hard<2>{});
    }
    val<NT> candallocmask = [&]() {
      val<NT> base = postmask & notumask.fo1();
      if constexpr (USE_SEC_TAG && SIBLING_POLICY == SiblingPolicy::ALL) {
        // Allocation promotion (Cai/Deshmukh/Patt ISCA'25 Sec 4.3):
        // Skip entries with same primary tag but different sec_tag (siblings).
        // NUM_GROUPS>1: also skip if same htag but different group_id.
        // Sibling = htag_hit & ~(group_match & sec_hit).
        // Promotes allocation to the next higher table.
        // SIBLING_TABLE_FLOOR: only skip for tables >= floor index.
        // Below floor: always 1 (never skip). At/above floor: check.
        static constexpr u64 FLOOR_MASK =
            SIBLING_TABLE_FLOOR >= NT ? ~u64(0)
                                      : ~((u64(1) << SIBLING_TABLE_FLOOR) - 1);
        // For NUM_GROUPS>1, compute alloc_group_id: the group of the
        // mispredicting (last) branch, used to check group_match.
        [[maybe_unused]] val<std::max(u64(1), GROUP_BITS)> alloc_gid_sib =
            val<std::max(u64(1), GROUP_BITS)>{
                NUM_GROUPS > 1
                    ? BRANCH_TO_GROUP[num_branch > 0 ? num_branch - 1 : 0]
                    : u64(0)};
        arr<val<1>, NT> not_sibling = [&](u64 i) -> val<1> {
          // Read train_tag_hit/train_sec_hit for ALL tables (runtime i).
          // For tables below FLOOR_MASK, the result is masked out below.
          // When DECAY also reads these, fanout(2) was declared and lvalue
          // reads are fine. When sibling is the only reader, use fo1.
          auto th = [&]() -> val<1> {
            if constexpr (DECAY_ENABLE)
              return val<1>{train_tag_hit[i]};
            else
              return train_tag_hit[i].fo1();
          }();
          auto sh = [&]() -> val<1> {
            if constexpr (DECAY_ENABLE)
              return val<1>{train_sec_hit[i]};
            else
              return train_sec_hit[i].fo1();
          }();
          if constexpr (NUM_GROUPS > 1) {
            // Group-aware sibling: sibling if htag matches but
            // (group differs OR sec_tag differs).
            // not_sibling = ~htag_hit | (group_match & sec_hit)
            auto gid = val<GROUP_BITS>{train_group_id[i].fo1()};
            val<1> group_match =
                gid.fo1() == val<GROUP_BITS>{alloc_gid_sib};
            return ~th.fo1() | (group_match.fo1() & sh.fo1());
          } else {
            return ~(th.fo1() & ~sh.fo1());
          }
        };
        // Apply floor: force not_sibling=1 for tables below floor.
        // FLOOR_MASK has bits set only for tables >= SIBLING_TABLE_FLOOR.
        val<NT> sibling_mask =
            not_sibling.fo1().concat() | ~val<NT>{hard<FLOOR_MASK>{}};
        return base.fo1() & sibling_mask.fo1();
      } else {
        return base.fo1();
      }
    }();
    // NOTE: @prakhar @claude audit
    candallocmask.fanout(hard<2>{}); // collamask(1) + noalloc(1)

    // acc_ctr / alloc_ctr fanout: read in collamask(1) +
    // gated_triggermask(FILTERED/THROTTLED)
    // + decay/epoch loop(NT each) + new_acc/new_alloc(1 each) + epoch_check(1
    // each) NOTE: @prakhar @claude audit
    static constexpr bool SEC_TAG_RUNTIME =
        USE_SEC_TAG && SecTagPolicy::RUNTIME && NUM_PATHS == 1;
    static constexpr u64 ACC_CTR_READS =
        1 + // collamask lambda
        (AllocCfg::ALLOC_ACTION == AllocAction::FILTERED
             ? 1
             : 0) +                                     // gated_triggermask
        ((DECAY_ENABLE || EPOCH_ENABLE) ? NT + 1 : 0) + // loop(NT) + new_acc
        (EPOCH_ENABLE ? 1 : 0) +                        // epoch_check
        (SEC_TAG_RUNTIME ? 1 : 0);                      // sec-tag gate
    static constexpr u64 ALLOC_CTR_READS =
        1 + // collamask lambda
        (AllocCfg::ALLOC_ACTION == AllocAction::THROTTLED
             ? 1
             : 0) +                                     // gated_triggermask
        ((DECAY_ENABLE || EPOCH_ENABLE) ? NT + 1 : 0) + // loop(NT) + new_alloc
        (EPOCH_ENABLE ? 1 : 0) +                        // epoch_check
        (SEC_TAG_RUNTIME ? 1 : 0);                      // sec-tag gate
    // NOTE: @prakhar @claude audit
    acc_ctr.fanout(hard<ACC_CTR_READS>{});
    // NOTE: @prakhar @claude audit
    alloc_ctr.fanout(hard<ALLOC_CTR_READS>{});
    // benefit_ctr fanout: gate read (1) + ta_update_ctr input (1) = 2
    // When HAS_BENEFIT_CTR is false, benefit_ctr is a dead 1-bit reg, skip.
    static constexpr u64 BENEFIT_CTR_READS =
        (SEC_TAG_RUNTIME && HAS_BENEFIT_CTR) ? 2 : 0;
    if constexpr (BENEFIT_CTR_READS > 0)
      benefit_ctr.fanout(hard<BENEFIT_CTR_READS>{});

    // 3d. Target policy (may skip closest candidates)
    // alloc_rng: STANDARD uses it once (apply only); FILTERED/THROTTLED use
    // it twice (lambda + apply), so declare fanout<2> in those paths.
    // NOTE: @prakhar @claude audit
    if constexpr (AllocCfg::ALLOC_ACTION != AllocAction::STANDARD)
      alloc_rng.fanout(hard<2>{}); // lambda read + apply read
    val<NT> collamask = [&]() -> val<NT> {
      if constexpr (AllocCfg::ALLOC_ACTION == AllocAction::STANDARD)
        return AllocCfg::TARGET_POLICY::template apply<NT>(
            candallocmask.reverse(), val<ALLOC_WIDTH>{alloc_ctr},
            val<ACC_WIDTH>{acc_ctr}, alloc_rng.fo1());
      else
        return AllocCfg::TARGET_POLICY::template apply<NT>(
            candallocmask.reverse(), val<ALLOC_WIDTH>{alloc_ctr},
            val<ACC_WIDTH>{acc_ctr}, alloc_rng);
    }();

    // 3e. Final allocation decision (one-hot or two-hot)
    arr<val<1>, NT> allocate = [&]() -> arr<val<1>, NT> {
      if constexpr (AllocCfg::MAX_ALLOC >= 2) {
        // NOTE: @prakhar @claude audit
        collamask.fanout(hard<AllocCfg::NON_CONSECUTIVE ? 3 : 2>{});
        val<NT> pick1 = collamask.one_hot();
        // NOTE: @prakhar @claude audit
        pick1.fanout(hard<AllocCfg::NON_CONSECUTIVE ? 5 : 2>{});
        val<NT> pick2 = [&]() -> val<NT> {
          val<NT> basic2 = (collamask ^ pick1).one_hot();
          if constexpr (AllocCfg::NON_CONSECUTIVE) {
            val<NT> neighbors = (pick1 << 1) | (pick1 >> 1);
            val<NT> nc_mask = (collamask ^ pick1) & ~neighbors.fo1();
            nc_mask.fanout(hard<2>{}); // != hard<0> + nc_pick
            val<NT> nc_pick = nc_mask.reverse().one_hot();
            return select(nc_mask != hard<0>{}, nc_pick.fo1(), basic2.fo1());
          } else {
            return basic2;
          }
        }();
        return (pick1 | pick2.fo1()).reverse().make_array(val<1>{});
      } else {
        return collamask.fo1().one_hot().reverse().make_array(val<1>{});
      }
    }();
    // NOTE: @prakhar @claude audit
    allocate.fanout(
        hard<3 + (DECAY_ENABLE ? 1 : 0)
#ifdef TAGE_MONITOR
             + (DECAY_ENABLE ? 1 : 0)
#endif
             >{}); // do_alloc(1) + alloc_target(1) +
                   // base_u_write(1) [+ decay(1)] [+ monitor(1)]
    val<NT> alloc_target = [&]() {
      arr<val<1>, NT> a = allocate;
      return a.fo1().concat();
    }();
    // NOTE: @prakhar @claude audit
    alloc_target.fanout(hard<2>{}); // line 1395(!= hard<0>) + line
                                    // 1401(allocmask1, if FARALLOC_DIST>0)
    // uclear: on alloc failure, optionally modify u-bits of tables above
    // provider. DISABLED = no-op, ZERO = clear to 0, DECREMENT = sat dec.
    arr<val<1>, NT> uclear = [&]() -> arr<val<1>, NT> {
      if constexpr (U_CLEAR == UClearPolicy::DISABLED) {
        return arr<val<1>, NT>{[](u64) -> val<1> { return hard<0>{}; }};
      } else {
        val<1> noalloc = (candallocmask == hard<0>{});
        val<NT> uclearmask =
            postmask & noalloc.fo1().replicate(hard<NT>{}).concat();
        return uclearmask.fo1().make_array(val<1>{});
      }
    }();
    // uclear: single read per table (fo1 in mask). No fanout needed.

#ifdef DEBUG_PRINT
    std::cerr << "--- training signals ---\n";
    t_match1.print("  t_match1=", "\n", true, std::cerr);
    t_pp.print("  t_pp=", "\n", true, std::cerr);
    t_pw.print("  t_pw=", "\n", true, std::cerr);
    t_ad.print("  t_ad=", "\n", true, std::cerr);
    t_phw.print("  t_phw=", "\n", true, std::cerr);
    actual_dir_g[0].print("  actual_dir=", "\n", true, std::cerr);
    any_provider_wrong.print("  any_prov_wrong=", "\n", true, std::cerr);
    do_train.print("  do_train=", "\n", true, std::cerr);
    std::cerr << "--- allocation path ---\n";
    alloc_trigger.print("  alloc_trigger=", "\n", true, std::cerr);
    gated_triggermask.print("  gated_trigmask=", "\n", true, std::cerr);
    alloc_base.print("  alloc_base=", "\n", true, std::cerr);
    notumask.print("  notumask=", "\n", true, std::cerr);
    postmask.print("  postmask=", "\n", true, std::cerr);
    candallocmask.print("  candallocmask=", "\n", true, std::cerr);
    alloc_target.print("  alloc_target=", "\n", true, std::cerr);
    std::cerr << "--- sec-tag train ---\n";
    train_sec_tag.print("  train_sec_tag=", "\n", true, std::cerr);
    for (u64 i = 0; i < NT; i++) {
      train_tag_hit[i].print(("  train_tag_hit[" + std::to_string(i) + "]=").c_str(), "\n", true, std::cerr);
      train_sec_hit[i].print(("  train_sec_hit[" + std::to_string(i) + "]=").c_str(), "\n", true, std::cerr);
    }
#endif

#ifdef TAGE_MONITOR
    {
      u64 at = static_cast<u64>(alloc_target);
      mon.record_allocation(at != 0, at);
      // Per-group allocation: attribute to the mispredicting branch's group
      if constexpr (NUM_GROUPS > 1) {
        u64 alloc_group = num_branch > 0 ? BRANCH_TO_GROUP[num_branch - 1] : 0;
        if (at != 0) {
          for (u64 i = 0; i < NT; i++)
            if (at & (u64(1) << i))
              mon.record_group_allocation(alloc_group, true, i);
        } else if (static_cast<u64>(alloc_trigger)) {
          mon.record_group_allocation(alloc_group, false, 0);
        }
      }
      // Provider index for cascade tracking
      u64 prov_idx =
          decltype(mon)::decode_provider(static_cast<u64>(t_match1));
      if (static_cast<u64>(alloc_trigger)) {
        if (static_cast<u64>(postmask) == 0)
          mon.record_alloc_blocked(); // fallback was provider, no candidates
        mon.record_alloc_cascade(prov_idx, at);
        if constexpr (USE_SEC_TAG) {
          u64 sibling = 0;
          u64 alloc_g = NUM_GROUPS > 1
              ? BRANCH_TO_GROUP[num_branch > 0 ? num_branch - 1 : 0]
              : 0;
          for (u64 i = 0; i < NT; i++) {
            bool th = static_cast<u64>(train_tag_hit[i]);
            bool sh = static_cast<u64>(train_sec_hit[i]);
            if constexpr (NUM_GROUPS > 1) {
              u64 gid = static_cast<u64>(
                  val<GROUP_BITS>{train_group_id[i]});
              bool gm = (gid == alloc_g);
              // Sibling: htag matches but (group or sec) differs
              if (th && !(gm && sh))
                sibling |= (u64(1) << i);
            } else {
              if (th && !sh)
                sibling |= (u64(1) << i);
            }
          }
          sibling &= static_cast<u64>(postmask); // only above provider
          if (sibling)
            mon.record_alloc_sibling_skip(sibling);
        }
      }
      // u-clear source tracking: per-table breakdown
      if constexpr (UCLEAR_ACTIVE) {
        // Reconstruct uclearmask from postmask/notumask (same logic as alloc).
        u64 pm = static_cast<u64>(postmask);
        u64 num = static_cast<u64>(notumask);
        u64 cam = static_cast<u64>(candallocmask);
        bool had_u0 = (pm & num) != 0;
        bool had_cand = cam != 0;
        bool noalloc = (pm & num) == 0;
        u64 ucm = noalloc ? pm : 0;
        for (u64 i = 0; i < NT; i++) {
          if (ucm & (u64(1) << i))
            mon.record_uclear_source(i, 0); // alloc_fail (genuine)
          // Would-have-been sibling uclear: u=0 existed but candallocmask
          // filtered it out (sibling skip). Under old code this triggered
          // uclear.
          if (!(ucm & (u64(1) << i)) && (pm & (u64(1) << i)) && had_u0 &&
              !had_cand && static_cast<u64>(alloc_trigger))
            mon.record_uclear_source(i, 1); // sibling (counterfactual)
        }
      }
    }
#endif

    // ---- Step 4: Fallback update ----
    // 4a. Direct update: mispredict + fallback is provider → write actual_dir.
    // NOTE: @prakhar @claude audit — fanout counts vary by FB_RECONCILE/FB_BIM_HYST
    static constexpr u64 FB_EXTRA = (FB_RECONCILE ? 1 : 0) + (FB_BIM_HYST ? 1 : 0);
    if constexpr (FB_RECONCILE) {
      // NOTE: @prakhar @claude audit
      train_fb_hyst.fanout(hard<2>{}); // fb_hyst_weak + silent update check
    }
    if constexpr (FB_EXTRA > 0) {
      train_fb.fanout(hard<1 + FB_EXTRA>{});   // NOTE: @prakhar @claude audit
      if constexpr (!FB_BIM_HYST) // BIM_HYST declares train_fb_idx fanout early (before need_extra_cycle)
        train_fb_idx.fanout(hard<1 + 2 * FB_EXTRA>{}); // NOTE: @prakhar @claude audit
    }
    if constexpr (FB_BIM_HYST)
      fb_bh_read.fanout(hard<2>{}); // fb_bh_weak + bh_changed
    // Full N-wide actual direction for fb path (fb_ctr stores FB_PRED_BITS = N*CTR_WIDTH)
    val<FB_PRED_BITS> actual_dir_full = [&]() -> val<FB_PRED_BITS> {
      if constexpr (NUM_GROUPS == 1) {
        return val<FB_PRED_BITS>{actual_dir_g[0]};
      } else {
        return arr<val<1>, N>{[&](u64 i) -> val<1> {
          return val<1>{branch_dir[i]};
        }}.concat();
      }
    }();
    if constexpr (NUM_GROUPS > 1 || FB_BIM_HYST)
      actual_dir_full.fanout(hard<1 + (NUM_GROUPS > 1 ? 1 : 0) + (FB_BIM_HYST ? 1 : 0)>{});
    val<1> fb_changed = [&]() -> val<1> {
      if constexpr (FB_EXTRA > 0)
        return actual_dir_full != val<FB_PRED_BITS>{train_fb};
      else
        return actual_dir_full != val<FB_PRED_BITS>{train_fb.fo1()};
    }();
    // fb_gate + fb write data depend on FB_BIM_HYST
    val<FB_PRED_BITS> fb_write_data = [&]() -> val<FB_PRED_BITS> {
      if constexpr (FB_BIM_HYST) {
        // Per-branch: only flip prediction where hyst is weak AND wrong
        val<FB_PRED_BITS> fb_wrong = actual_dir_full ^ val<FB_PRED_BITS>{train_fb};
        val<FB_PRED_BITS> fb_bh_weak = ~val<FB_PRED_BITS>{fb_bh_read};
        val<FB_PRED_BITS> fb_flip = fb_wrong & fb_bh_weak;
        return val<FB_PRED_BITS>{train_fb} ^ fb_flip;
      } else {
        return actual_dir_full;
      }
    }();
    val<1> fb_gate = [&]() -> val<1> {
      if constexpr (FB_EXTRA > 0)
        return do_train & t_m1[NT] & mispredict & fb_changed.fo1();
      else
        return do_train & t_m1[NT].fo1() & mispredict & fb_changed.fo1();
    }();
    // Helper: write to banked fb RAMs
    auto fb_write = [&](auto idx_val, auto data_val, val<1> gate) {
      if constexpr (FB_BANKS == 1 && !FB_SPLIT_WIDTH) {
        execute_if(gate, [&]() {
          fb_ctr[0][0].write(idx_val, data_val);
        });
      } else if constexpr (FB_BANKS == 1 && FB_SPLIT_WIDTH) {
        execute_if(gate, [&]() {
          for (u64 g = 0; g < NUM_GROUPS; g++) {
            fb_ctr[g][0].write(idx_val, val<1>{data_val >> g});
          }
        });
      } else if constexpr (FB_BANKS > 1 && !FB_SPLIT_WIDTH) {
        auto bank_sel = val<FB_BANK_SEL_BITS>{idx_val >> FB_BANK_SEL_SHIFT};
        auto bank_idx = [&]() {
          if constexpr (FB_BANK_SEL_SHIFT == 0)
            return val<FB_BANK_IDX_BITS>{idx_val >> FB_BANK_SEL_BITS};
          else
            return val<FB_BANK_IDX_BITS>{
              concat(val<FB_BANK_UPPER_BITS>{idx_val >> (FB_BANK_SEL_SHIFT + FB_BANK_SEL_BITS)},
                     val<FB_BANK_SEL_SHIFT>{idx_val})};
        }();
        static_loop<FB_BANKS>([&]<u64 B>() {
          val<1> bank_en = gate & (bank_sel == hard<B>{});
          execute_if(bank_en, [&]() {
            fb_ctr[0][B].write(bank_idx, data_val);
          });
        });
      } else {
        auto bank_sel = val<FB_BANK_SEL_BITS>{idx_val >> FB_BANK_SEL_SHIFT};
        auto bank_idx = [&]() {
          if constexpr (FB_BANK_SEL_SHIFT == 0)
            return val<FB_BANK_IDX_BITS>{idx_val >> FB_BANK_SEL_BITS};
          else
            return val<FB_BANK_IDX_BITS>{
              concat(val<FB_BANK_UPPER_BITS>{idx_val >> (FB_BANK_SEL_SHIFT + FB_BANK_SEL_BITS)},
                     val<FB_BANK_SEL_SHIFT>{idx_val})};
        }();
        static_loop<FB_BANKS>([&]<u64 B>() {
          val<1> bank_en = gate & (bank_sel == hard<B>{});
          execute_if(bank_en, [&]() {
            for (u64 g = 0; g < NUM_GROUPS; g++) {
              fb_ctr[g][B].write(bank_idx, val<1>{data_val >> g});
            }
          });
        });
      }
    };

    // NOTE: @prakhar @claude audit
    fb_gate.fanout(hard<5>{}); // fb_ctr.write: B=2 bank writes + write_bank +
                               // write_localaddr + write_data
    if constexpr (FB_EXTRA > 0)
      fb_write(val<FB_IDX_BITS>{train_fb_idx}, fb_write_data, fb_gate);
    else
      fb_write(val<FB_IDX_BITS>{train_fb_idx.fo1()}, fb_write_data, fb_gate);

    // 4b. Bimodal hysteresis update: correct→strong(1), wrong→weak(0)
    if constexpr (FB_BIM_HYST) {
      val<FB_PRED_BITS> fb_wrong = actual_dir_full ^ val<FB_PRED_BITS>{train_fb};
      val<FB_PRED_BITS> new_bh = ~fb_wrong; // 1=correct, 0=wrong
      val<1> bh_changed = new_bh != val<FB_PRED_BITS>{fb_bh_read};
      val<1> bh_gate = do_train & t_m1[NT] & mispredict & bh_changed;
      execute_if(bh_gate, [&]() {
        fb_bim_hyst.write(val<FB_BH_IDX_BITS>{train_fb_idx}, new_bh);
      });
    }

    // 4c. Reconciliation (Tage.hpp P1/P2 pattern): when TAGE and fb disagree
    // persistently (fb_hyst weak), overwrite fb with TAGE's prediction.
    // Also update fb_hyst every cycle to track agreement.
    if constexpr (FB_RECONCILE) {
      val<1> fb_not_provider = ~t_m1[NT];
      val<1> fb_tage_disagree = (val<PRED_BITS>{train_fb} ^ t_pp) != hard<0>{};
      val<1> fb_hyst_weak = ~val<1>{train_fb_hyst.fo1()};

      // Overwrite fb pred when hyst weak and they disagree (non-provider only;
      // provider case handled above by fb_gate)
      val<1> reconcile_gate =
          do_train & fb_not_provider & fb_hyst_weak & fb_tage_disagree;
      fb_write(val<FB_IDX_BITS>{train_fb_idx}, t_pp, reconcile_gate);

      // Update fb_hyst: 1 if agree, 0 if disagree — tracks recent consensus.
      // Silent update elimination: skip write when new == old.
      val<1> new_fb_hyst = ~fb_tage_disagree;
      val<1> fb_hyst_changed = new_fb_hyst != val<1>{train_fb_hyst};
      execute_if(do_train & fb_hyst_changed, [&]() {
        fb_hyst.write(val<FB_IDX_BITS>{train_fb_idx}, new_fb_hyst);
      });
    }

#ifdef TAGE_MONITOR
    if (static_cast<u64>(fb_gate))
      mon.record_fb_write(static_cast<u64>(val<FB_IDX_BITS>{train_fb_idx}));
#endif

    // ---- Step 5: Meta counter update ----
    if constexpr (META_BANKS == 1) {
      // Single meta bank: original behavior — uses any_provider_wrong
      auto old_meta = val<META_WIDTH, i64>{meta_pipe[0][META_PIPE - 1]};
      old_meta.fanout(hard<2>{});
      auto new_meta = [&]() {
        if constexpr (AllocCfg::ALLOC_TRIGGER == AllocTrigger::TAGE_WRONG)
          return ta_update_ctr(old_meta, any_provider_wrong);
        else
          return ta_update_ctr(old_meta, any_provider_wrong.fo1());
      }();
      new_meta.fanout(hard<2>{});
      val<1> meta_gate = do_train & t_pw & t_ad & (new_meta != old_meta);
      meta_gate.fanout(hard<5>{});
      execute_if(meta_gate, [&]() {
        meta_ctr[0].write(
            val<META_IDX_BITS>{meta_idx_pipe[0][META_PIPE - 1].fo1()},
            new_meta, hard<0>{});
      });
    } else {
      // Per-bank meta: each bank trained with its own group's wrong signal,
      // pw, and ad. Multiple groups may map to same bank — fold_or them.
      // Constexpr: which groups map to each bank?
      static constexpr auto META_BANK_GROUPS = []() {
        // groups_for_bank[mb] = bitmask of groups mapping to bank mb
        std::array<u64, META_BANKS> mask{};
        for (u64 g = 0; g < NUM_GROUPS; g++)
          mask[GROUP_TO_META_BANK[g]] |= (1ULL << g);
        return mask;
      }();

      // Number of groups per meta bank (constexpr)
      static constexpr auto META_BANK_GROUP_COUNT = []() {
        std::array<u64, META_BANKS> cnt{};
        for (u64 g = 0; g < NUM_GROUPS; g++)
          cnt[GROUP_TO_META_BANK[g]]++;
        return cnt;
      }();
      static constexpr u64 MAX_GROUPS_PER_BANK = []() {
        u64 mx = 0;
        for (u64 mb = 0; mb < META_BANKS; mb++)
          if (META_BANK_GROUP_COUNT[mb] > mx) mx = META_BANK_GROUP_COUNT[mb];
        return mx;
      }();

      for (u64 mb = 0; mb < META_BANKS; mb++) {
        // Per-group gating: gate fires only when pw AND ad are true for the
        // SAME group. Direction comes only from groups where gate fired.
        // This avoids the pairing bug where OR(pw) & OR(ad) could fire from
        // different groups, and ensures correct decrement (groups where gate
        // fires but prediction was correct push counter back toward provider).
        u64 gi = 0;
        // bank_gate_arr[j] = pw_g & ad_g  (group g has weak provider with differing alt)
        arr<val<1>, MAX_GROUPS_PER_BANK> bank_gate_arr = [&]([[maybe_unused]] u64 j) -> val<1> {
          for (u64 g = gi; g < NUM_GROUPS; g++) {
            if (GROUP_TO_META_BANK[g] == mb) {
              gi = g + 1;
              return pre_read_pw[g] & pre_read_ad[g];
            }
          }
          return hard<0>{}; // padding (OR identity)
        };
        gi = 0;
        // bank_dir_arr[j] = pw_g & ad_g & wrong_g  (direction: wrong only from gated groups)
        arr<val<1>, MAX_GROUPS_PER_BANK> bank_dir_arr = [&]([[maybe_unused]] u64 j) -> val<1> {
          for (u64 g = gi; g < NUM_GROUPS; g++) {
            if (GROUP_TO_META_BANK[g] == mb) {
              gi = g + 1;
              return pre_read_pw[g] & pre_read_ad[g] & per_group_wrong[g];
            }
          }
          return hard<0>{};
        };
        val<1> bank_gate = bank_gate_arr.fo1().fold_or();
        val<1> bank_wrong_dir = bank_dir_arr.fo1().fold_or();

        auto old_meta = val<META_WIDTH, i64>{meta_pipe[mb][META_PIPE - 1]};
        old_meta.fanout(hard<2>{});
        auto new_meta = ta_update_ctr(old_meta, bank_wrong_dir);
        new_meta.fanout(hard<2>{});
        val<1> meta_gate = do_train & bank_gate & (new_meta != old_meta);
        meta_gate.fanout(hard<5>{});
        execute_if(meta_gate, [&]() {
          meta_ctr[mb].write(
              val<META_IDX_BITS>{meta_idx_pipe[mb][META_PIPE - 1].fo1()},
              new_meta, hard<0>{});
        });
      }
    }

    // ---- Merged per-table writes (one write per RAM per table) ----
    // For each table: alloc takes priority over update. Mux selects data.
    // Uses per-table wrong signal (this table's own prediction vs actual)
    // instead of blanket any_provider_wrong. For the provider table these
    // are identical; for non-provider tables this enables future per-table
    // update policies (e.g. partial update).
    // NOTE: @prakhar @claude audit
    train_pc.fanout(hard<NT>{});
    static_loop<NT>([&]<u64 I>() {
      auto &t = table<I>();
      val<1> do_alloc = allocate[I];
      // NOTE: @prakhar @claude audit
      do_alloc.fanout(hard<4>{}); // 1252(pred gate) + 1261(hyst select) +
                                  // 1263(hyst gate) + 1269(alloc gate)

      // Per-table wrong: does THIS table's stored pred disagree with actual?
      // We always compare against the training group's actual directions.
      // Provider: stored group == train_group (group-aware resolution ensures this).
      // Alloc: new entry is for train_group.
      // No-write tables: value is unused.
      auto entry_actual_dir = val<PRED_BITS>{actual_dir_g[train_group]};
      entry_actual_dir.fanout(hard<2>{}); // table_wrong(1) + pred_ram.write(1)
      val<1> table_wrong =
          (val<PRED_BITS>{train_pred[I].fo1()} ^ entry_actual_dir) != hard<0>{};
      // NOTE: @prakhar @claude audit
      // do_pred_update(1) + new_hyst(1) + u_correct(1) [+ u_wrong(1) if
      // !UNTOUCHED]
      table_wrong.fanout(
          hard<3 + (U_MISP != UMispPolicy::UNTOUCHED ? 1 : 0)>{});

      // pred_ram: alloc writes entry_actual_dir, update writes entry_actual_dir.
      // Gated on hyst-only weakness (t_phw), not newly-alloc (t_pw). A
      // useful entry (u>0) with weak hyst that's wrong should still flip.
      // Silent: do_pred_update already requires table_wrong (pred != actual),
      // so non-alloc writes only fire when value changes.
      val<1> do_pred_update = t_m1[I] & t_phw & table_wrong;
      // NOTE: @prakhar @claude audit
      val<1> gate_pred = do_train & (do_alloc | do_pred_update.fo1());
      // NOTE: @prakhar @claude audit
      gate_pred.fanout(hard<5>{}); // pred_ram ta_rwram B=2: B+3=5 writes
      execute_if(gate_pred, [&]() {
        t.pred_ram.write(val<t.IDX_BITS>{train_idx[I]}, entry_actual_dir, hard<0>{});
      });

      // hyst_ram: alloc writes 0, update writes new_hyst → mux on do_alloc.
      // Direction based on this table's own prediction accuracy.
      constexpr u64 HW = std::max(u64(1), HYST_WIDTH);
      auto old_hyst = val<HW>{train_hyst[I]};
      // NOTE: @prakhar @claude audit
      old_hyst.fanout(
          hard<2>{}); // ta_update_ctr(1) + (new_hyst != old_hyst)(1)
      auto new_hyst = ta_update_ctr(old_hyst, ~table_wrong);
      // NOTE: @prakhar @claude audit
      new_hyst.fanout(hard<2>{}); // select(1) + (new_hyst != old_hyst)(1)
      auto hyst_data = select(do_alloc, val<HW>{0}, new_hyst);
      val<1> do_hyst_update = t_m1[I] & (new_hyst != old_hyst);
      // NOTE: @prakhar @claude audit
      val<1> gate_hyst = do_train & (do_alloc | do_hyst_update.fo1());
      // NOTE: @prakhar @claude audit
      gate_hyst.fanout(hard<5>{}); // hyst_ram ta_rwram B=2: B+3=5 writes
      execute_if(gate_hyst, [&]() {
        t.hyst_ram.write(val<t.HYST_IDX_BITS>{train_idx[I]}, hyst_data.fo1(),
                         hard<0>{});
      });

      // tag_ram + sec_ram: alloc only (plain RAM, protected by extra_cycle)
      // NOTE: @prakhar @claude audit
      val<1> gate_alloc = do_train & do_alloc;
      gate_alloc.fanout(
          hard<2>{}); // tag_ram hcm::ram(1) + sec_ram hcm::ram(1, USE_SEC_TAG)
      execute_if(gate_alloc, [&]() {
        // Use piped computed tag (from predict1 time, not current folds)
        // NUM_GROUPS>1: prepend alloc_group_id to htag
        if constexpr (NUM_GROUPS > 1) {
          // alloc_group_id = group of the last (mispredicting) branch
          static constexpr u64 PER_HTAG = t.tag_width - GROUP_BITS;
          val<GROUP_BITS> alloc_gid =
              val<GROUP_BITS>{BRANCH_TO_GROUP[num_branch > 0 ? num_branch - 1 : 0]};
          auto full_tag = concat(alloc_gid, val<PER_HTAG>{train_ctag[I].fo1()});
          t.tag_ram.write(val<t.IDX_BITS>{train_idx[I]},
                          val<t.tag_ram_width>{full_tag});
        } else {
          t.tag_ram.write(val<t.IDX_BITS>{train_idx[I]},
                          val<t.tag_ram_width>{train_ctag[I].fo1()});
        }
        if constexpr (USE_SEC_TAG)
          t.sec_ram.write(val<t.IDX_BITS>{train_idx[I]},
                          val<SEC_TAG_BITS>{train_sec_tag});
      });
#ifdef TAGE_MONITOR
      {
        u64 tidx = static_cast<u64>(val<t.IDX_BITS>{train_idx[I]});
        if (static_cast<u64>(do_train & do_alloc)) {
          mon.record_tage_write(I, tidx);
          mon.record_entry_alloc_diag(I, tidx);
          mon.record_bank_write(I, tidx, decltype(mon)::BankWriteType::ALLOC);
          if constexpr (NUM_GROUPS > 1) {
            u64 gid = static_cast<u64>(val<GROUP_BITS>{train_group_id[I]});
            mon.record_group_occupy(gid, I, tidx);
          }
        }
        if (static_cast<u64>(gate_pred))
          mon.record_bank_write(I, tidx, decltype(mon)::BankWriteType::PRED_UPDATE);
        if (static_cast<u64>(gate_hyst))
          mon.record_bank_write(I, tidx, decltype(mon)::BankWriteType::HYST_UPDATE);
      }
#endif

      // u_ram: combined provider update + allocation + uclear + decay.
      // Standard TAGE: increment u when correct & alt differs, optionally
      // modify on wrong. U_MISP controls wrong-prediction behavior:
      //   UNTOUCHED: leave u alone (default, preserves accumulated usefulness)
      //   ZERO:      clear u to 0 (aggressive, original buggy behavior)
      //   DECREMENT: saturating decrement (gradual degradation)
      val<1> u_correct = t_m1[I] & t_ad & ~table_wrong;
      // Single fo1() read — high-drive to mask. No fanout needed.
      val<1> u_wrong = [&]() -> val<1> {
        if constexpr (U_MISP == UMispPolicy::UNTOUCHED)
          return hard<0>{};
        else
          return t_m1[I] & t_ad & table_wrong;
      }();
      // Single fo1() read — high-drive to mask. No fanout needed.
      val<U_WIDTH> old_u = val<U_WIDTH>{train_u[I]};
      // Need u_dec when U_MISP==DECREMENT or U_CLEAR==DECREMENT
      static constexpr bool NEED_DEC = U_MISP == UMispPolicy::DECREMENT ||
                                       U_CLEAR == UClearPolicy::DECREMENT;
      // NOTE: @prakhar @claude audit:
      //   u_inc: !=maxval(1) + +(1) = 2
      //   u_dec (NEED_DEC): !=0(1) + -(1) = 2
      //   silent update (!DECAY): base_newu != old_u = 1
      //   DECAY fallback: old_u in select = 1
      old_u.fanout(hard<2 + (NEED_DEC ? 2 : 0) + 1 + (DECAY_ENABLE ? 1 : 0)>{});
      // Saturating increment: old_u + 1, clamped at maxval.
      // Shallower than select-mux: compare + add, no mux in path.
      auto u_inc =
          val<U_WIDTH>{old_u + val<U_WIDTH>{old_u != hard<old_u.maxval>{}}};
      auto u_dec = [&]() -> val<U_WIDTH> {
        if constexpr (NEED_DEC)
          // Saturating decrement: old_u - 1, clamped at 0.
          return val<U_WIDTH>{old_u - val<U_WIDTH>{old_u != hard<0>{}}};
        else
          return val<U_WIDTH>{0}; // unused placeholder
      }();

      // uclear value: DISABLED → old_u (no-op), ZERO → 0, DECREMENT → dec
      auto uclear_val = [&]() -> val<U_WIDTH> {
        if constexpr (U_CLEAR == UClearPolicy::DISABLED)
          return val<U_WIDTH>{0}; // never reached (uclear always 0)
        else if constexpr (U_CLEAR == UClearPolicy::ZERO)
          return val<U_WIDTH>{0};
        else
          return u_dec.fo1();
      }();

      // Flat OR-of-ANDs: conditions are mutually exclusive per table
      // (provider can't be allocated, uclear can't fire on allocated table).
      // allocate → 0 (vanishes in OR), no-op → don't write (data irrelevant).
      // Depth: 1 AND + 1 OR, vs 3 nested selects (3 mux levels).
      // u_correct/u_wrong use fo1() — high-drive, single read each.
      val<U_WIDTH> base_newu = [&]() -> val<U_WIDTH> {
        auto uc_mask =
            val<U_WIDTH>{u_correct.fo1().replicate(hard<U_WIDTH>{}).concat()};
        auto inc_term = val<U_WIDTH>{uc_mask.fo1() & u_inc.fo1()};
        if constexpr (U_MISP == UMispPolicy::UNTOUCHED) {
          if constexpr (U_CLEAR == UClearPolicy::DISABLED)
            return inc_term;
          else {
            auto ucl_mask = val<U_WIDTH>{
                uclear[I].fo1().replicate(hard<U_WIDTH>{}).concat()};
            return inc_term.fo1() |
                   val<U_WIDTH>{ucl_mask.fo1() & uclear_val.fo1()};
          }
        } else if constexpr (U_MISP == UMispPolicy::ZERO) {
          // u_wrong → 0 (vanishes in OR, same as allocate)
          if constexpr (U_CLEAR == UClearPolicy::DISABLED)
            return inc_term;
          else {
            auto ucl_mask = val<U_WIDTH>{
                uclear[I].fo1().replicate(hard<U_WIDTH>{}).concat()};
            return inc_term.fo1() |
                   val<U_WIDTH>{ucl_mask.fo1() & uclear_val.fo1()};
          }
        } else { // DECREMENT
          auto uw_mask =
              val<U_WIDTH>{u_wrong.fo1().replicate(hard<U_WIDTH>{}).concat()};
          auto wrong_term = val<U_WIDTH>{uw_mask.fo1() & u_dec.fo1()};
          if constexpr (U_CLEAR == UClearPolicy::DISABLED)
            return inc_term.fo1() | wrong_term;
          else {
            auto ucl_mask = val<U_WIDTH>{
                uclear[I].fo1().replicate(hard<U_WIDTH>{}).concat()};
            return inc_term.fo1() | wrong_term.fo1() |
                   val<U_WIDTH>{ucl_mask.fo1() & uclear_val.fo1()};
          }
        }
      }();
      // Probabilistic decay: on tag/sec miss, random < threshold → decay u
      // Silent update elimination: only write when new != old (both paths).
#ifdef TAGE_MONITOR
      bool mon_is_base_write = false; // set inside lambda for decay fire detection
#endif
      auto [newu, u_write] = [&]() -> std::pair<val<U_WIDTH>, val<1>> {
        if constexpr (!DECAY_ENABLE) {
          // NOTE: @prakhar @claude audit
          // base_newu reads: != 0 (base_u_write) + != old_u (changed) + return
          // = 3
          base_newu.fanout(hard<3>{});
          // base_u_write: derived from base_newu to avoid second read of
          // u_correct/u_wrong. Since u_inc >= 1, (base_newu != 0) captures
          // u_correct/u_wrong/uclear. allocate writes 0, OR it separately.
          val<1> base_u_write = (base_newu != hard<0>{}) | allocate[I];
          val<1> base_changed = base_newu != old_u;
          return {base_newu, base_u_write.fo1() & base_changed.fo1()};
        } else {
          constexpr u64 LW = DECAY_LFSR_WIDTHS[I];

          // Miss condition from piped tag/sec hit
          val<1> tag_missed = ~val<1>{train_tag_hit[I]};
          val<1> decay_miss = [&]() {
            if constexpr (!USE_SEC_TAG || DECAY_MISS == DecayMiss::TAG)
              return tag_missed.fo1();
            else {
              val<1> sec_missed = ~val<1>{train_sec_hit[I]};
              if constexpr (DECAY_MISS == DecayMiss::SEC)
                return sec_missed.fo1();
              else if constexpr (DECAY_MISS == DecayMiss::TAG_OR_SEC)
                return tag_missed.fo1() | sec_missed.fo1();
              else
                return tag_missed.fo1() & sec_missed.fo1();
            }
          }();

          // Threshold from global counters
          auto thresh = DecayThreshFn::template compute<I, LW>(
              val<ACC_WIDTH>{acc_ctr}, val<ALLOC_WIDTH>{alloc_ctr});

          // Random fires when below threshold, on miss, not allocating
          val<LW> rng = val<LW>{static_cast<u64>(std::rand())};
          val<1> decay_fire =
              decay_miss.fo1() & ~allocate[I] & (rng.fo1() < thresh.fo1());

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

          // base_u_write for DECAY path: same derivation as !DECAY.
          // base_newu reads: != 0 (base_u_write) + select (merged) + return = 3
          base_newu.fanout(hard<3>{});
          val<1> base_u_write = (base_newu != hard<0>{}) | allocate[I];
          base_u_write.fanout(hard<2
#ifdef TAGE_MONITOR
                                    + 1
#endif
                                    >{}); // select(1) + merged_write(1) [+ monitor(1)]
#ifdef TAGE_MONITOR
          mon_is_base_write = static_cast<u64>(base_u_write);
#endif

          // Mux: if base write active, use base_newu; else if decay, use
          // decayed_u
          decay_fire.fanout(hard<2>{}); // select(1) + merged_write(1)
          val<U_WIDTH> merged =
              select(base_u_write, base_newu,
                     select(decay_fire, decayed_u.fo1(), old_u));
          val<1> merged_write = base_u_write | decay_fire;
          merged.fanout(hard<2>{}); // merged_changed(1) + return(1)
          val<1> merged_changed = merged != old_u;
#ifdef DEBUG_PRINT
          if constexpr (I == 0) {
            std::cerr << "--- decay path [T0] ---\n";
            decay_miss.print("  decay_miss=", "\n", true, std::cerr);
            decay_fire.print("  decay_fire=", "\n", true, std::cerr);
            decayed_u.print("  decayed_u=", "\n", true, std::cerr);
            merged.print("  merged=", "\n", true, std::cerr);
            merged_write.print("  merged_write=", "\n", true, std::cerr);
          }
#endif
          return {merged, merged_write.fo1() & merged_changed.fo1()};
        }
      }();

      // NOTE: @prakhar @claude audit
#ifdef TAGE_MONITOR
      u_write.fanout(hard<2>{}); // gate_u(1) + monitor peek(1)
      newu.fanout(hard<2>{});    // u_ram.write(1) + monitor peek(1)
      val<1> gate_u = do_train & u_write;
#else
      val<1> gate_u = do_train & u_write.fo1();
#endif
      gate_u.fanout(hard<5>{}); // u_ram ta_rwram B=2: B+3=5 writes
      execute_if(gate_u, [&]() {
#ifdef TAGE_MONITOR
        t.u_ram.write(val<t.IDX_BITS>{train_idx[I]}, newu, hard<0>{});
#else
        t.u_ram.write(val<t.IDX_BITS>{train_idx[I]}, newu.fo1(), hard<0>{});
#endif
      });
#ifdef DEBUG_PRINT
      if constexpr (I == 0) {
        std::cerr << "--- u-bit path [T0] ---\n";
        old_u.print("  old_u=", "\n", true, std::cerr);
        u_correct.print("  u_correct=", "\n", true, std::cerr);
        u_wrong.print("  u_wrong=", "\n", true, std::cerr);
        uclear[I].print("  uclear[0]=", "\n", true, std::cerr);
        allocate[I].print("  allocate[0]=", "\n", true, std::cerr);
        base_newu.print("  base_newu=", "\n", true, std::cerr);
        u_write.print("  u_write=", "\n", true, std::cerr);
        gate_u.print("  gate_u=", "\n", true, std::cerr);
      }
#endif

#ifdef TAGE_MONITOR
      if (static_cast<u64>(gate_u)) {
        mon.record_u_write(I, static_cast<u64>(newu) != 0);
        u64 tidx_u = static_cast<u64>(val<t.IDX_BITS>{train_idx[I]});
        mon.record_bank_write(I, tidx_u, decltype(mon)::BankWriteType::U_WRITE);
      }
      if constexpr (DECAY_ENABLE) {
        // Decay fire: u_write active but not from base logic.
        if (static_cast<u64>(gate_u) && !mon_is_base_write)
          mon.record_decay_fire();
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
      // NOTE: @prakhar @claude audit
      any_alloc.fanout(hard<2>{}); // FARALLOC_DIST=0: 1 read; FARALLOC_DIST>0:
                                   // 2 reads (1407+1408)
      auto new_alloc = [&]() {
        if constexpr (FARALLOC_DIST > 0) {
          // Far allocation: alloc target is >= FARALLOC_DIST tables from
          // provider. Shift provider mask down by FARALLOC_DIST; if alloc
          // target is still the one_hot winner, it's far → extra decrement.
          val<NT> allocmask1 = alloc_target;
          // NOTE: @prakhar @claude audit
          allocmask1.fanout(hard<2>{}); // used twice in line below: |
                                        // allocmask1, ^ allocmask1
          val<NT> provider_bits = val<NT>{t_match1};
          val<1> faralloc =
              (((provider_bits >> FARALLOC_DIST) | allocmask1).one_hot() ^
               allocmask1) == hard<0>{};
          // Two-step update: first normal, then extra decrement if far
          auto step1 = ta_update_ctr(val<ALLOC_WIDTH>{alloc_ctr}, ~any_alloc);
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
        // NOTE: @prakhar @claude audit
        epoch_fire.fanout(
            hard<1 + (EPOCH_RESET_ACC ? 1 : 0) +
                 (EPOCH_RESET_ALLOC ? 1 : 0)>{}); // execute_if + selects:
                                                  // EPOCH_RESET_ACC(1) +
                                                  // EPOCH_RESET_ALLOC(1)
        execute_if(epoch_fire, [&]() {
          static_loop<NT>([&]<u64 I>() { table<I>().u_ram.reset(); });
        });
#ifdef TAGE_MONITOR
        if (static_cast<u64>(epoch_fire)) {
          mon.record_epoch_reset();
          for (u64 i = 0; i < NT; i++)
            mon.record_uclear_source(i, 2); // epoch
        }
#endif
        // Optionally reset counters on epoch fire
        if constexpr (EPOCH_RESET_ACC)
          acc_ctr = select(epoch_fire, val<ACC_WIDTH>{0}, new_acc.fo1());
        else
          acc_ctr = new_acc.fo1();
        if constexpr (EPOCH_RESET_ALLOC)
          alloc_ctr = select(epoch_fire, val<ALLOC_WIDTH>{0}, new_alloc.fo1());
        else
          alloc_ctr = new_alloc.fo1();
      } else {
        acc_ctr = new_acc.fo1();
        alloc_ctr = new_alloc.fo1();
      }
    }

    // ---- Sec-tag adaptive benefit counter update ----
    if constexpr (HAS_BENEFIT_CTR && USE_SEC_TAG) {
      // Counterfactual sec-tag benefit signal.
      // Update only when enforcement WOULD have changed the provider to fallback:
      //   - at least one table has a tag match (tag_hit)
      //   - NO table has a full match (tag_hit && sec_hit)
      // Only then does "enforce vs don't enforce" causally map to "FB vs TAGE".
      u64 tag_only_prov = NT; // best tag-only provider (tag_hit, !sec_hit)
      bool any_full_match = false;
      for (u64 i = 0; i < NT; i++) {
        bool th = static_cast<u64>(train_tag_hit[i]);
        bool sh = static_cast<u64>(train_sec_hit[i]);
        if (th && sh) { any_full_match = true; }
        else if (th && tag_only_prov == NT) { tag_only_prov = i; }
      }
      // Fires only when sec-tag would have rejected ALL TAGE matches → FB
      bool sec_would_reject_all = (tag_only_prov < NT) && !any_full_match;
      bool should_update = false;
      bool benefit_incr_val = true;
      if (sec_would_reject_all && num_branch > 0) {
        bool tage_right = (static_cast<u64>(train_pred[tag_only_prov]) & 1) ==
                          (static_cast<u64>(actual_dir_g[0]) & 1);
        bool fb_right = (static_cast<u64>(train_fb) & 1) ==
                        (static_cast<u64>(actual_dir_g[0]) & 1);
        if (fb_right != tage_right) {
          should_update = true;
          benefit_incr_val = fb_right; // true=incr (FB>TAGE, enforcement helps)
                                       // false=decr (TAGE>FB, enforcement hurts)
        }
      }
      // Always read and compute (fanout-safe), conditionally write back
      val<1> b_incr{benefit_incr_val ? 1ULL : 0ULL};
      auto new_benefit =
          ta_update_ctr(val<BENEFIT_REG_WIDTH>{benefit_ctr}, b_incr);
      if (should_update)
        benefit_ctr = new_benefit.fo1();
#ifdef TAGE_MONITOR
      mon.record_benefit_update(
          sec_would_reject_all && static_cast<u64>(do_train),
          should_update,
          benefit_incr_val);
      mon.record_benefit_ctr(static_cast<u64>(benefit_ctr));
#endif
    }

#ifdef TAGE_MONITOR
    if constexpr (DECAY_ENABLE || EPOCH_ENABLE)
      mon.record_pressure(static_cast<u64>(val<ACC_WIDTH>{acc_ctr}),
                          static_cast<u64>(val<ALLOC_WIDTH>{alloc_ctr}));
#endif

    // Decay LFSRs replaced by std::rand() — no tick needed.

    // ---- Step 8: History update ----
    // true_block uses framework's mispredict signal (not our computed
    // correct_pred) to avoid timing bleed from old_pred reg reads
    true_block = ~mispredict | val<1>{branch_dir[num_branch - 1]} | line_end();
    // NOTE: @prakhar @claude audit
    true_block.fanout(
        hard<NT * 2 + 2 + (USE_GSHARE ? 1 : 0)>{}); // NT fold_idx + NT fold_tag
                                                    // apply_update muxes +
                                                    // gh.update + monitor +
                                                    // (fb_fold if gshare)

#ifdef TAGE_MONITOR
    if (static_cast<u64>(true_block))
      mon.record_true_block();
#endif

    // Compute new fold values OUTSIDE execute_if — runs in parallel with
    // true_block gate, so timing is max(fold_computation, true_block)
    // instead of additive (true_block + fold_computation).
    auto hist_input = [&]() {
      if constexpr (HIST_MODE == HistUpdate::PATH)
        return val<PATHBITS>{block_end_info.next_pc.fo1() >> 2};
      else if constexpr (HIST_MODE == HistUpdate::DIR)
        return val<1>{branch_dir[num_branch - 1]};
      else // BOTH: concat(direction, path)
        return concat(val<1>{branch_dir[num_branch - 1]},
                      val<PATHBITS>{block_end_info.next_pc.fo1() >> 2});
    }();

    // NOTE: @prakhar @claude audit
    hist_input.fanout(hard<NT * 2 + 1 + (USE_GSHARE ? 1 : 0)>{});
    gh.template fanout_per_bit<GH_FANOUT>();

    // Per-table: compute new fold values, apply with mux (no execute_if gate).
    // select(true_block, new, old) avoids the execute_if timing bleed —
    // both paths resolve in parallel, mux adds ~10ps constant overhead.
    static_loop<NT>([&]<u64 I>() {
      auto &t = table<I>();
      auto new_idx = t.fold_idx.compute_update(
          gh, hard<TableCfg::HIST_LEN[I]>{}, hist_input);
      auto new_tag = t.fold_tag.compute_update(
          gh, hard<TableCfg::HIST_LEN[I]>{}, hist_input);
      t.fold_idx.apply_update(new_idx.fo1(), true_block);
      t.fold_tag.apply_update(new_tag.fo1(), true_block);
#ifdef DEBUG_PRINT
      if constexpr (I == 0) {
        std::cerr << "=== update_cycle fold[0] ===\n";
        new_idx.print("  new_idx=", "\n", true, std::cerr);
        new_tag.print("  new_tag=", "\n", true, std::cerr);
        t.fold_idx.get().print("  fold_idx_after=", "\n", true, std::cerr);
        t.fold_tag.get().print("  fold_tag_after=", "\n", true, std::cerr);
      }
#endif
    });
    if constexpr (USE_GSHARE) {
      auto new_fb = fb_fold.compute_update(gh, hard<GS_HIST>{}, hist_input);
      fb_fold.apply_update(new_fb.fo1(), true_block);
    }
    gh.update(hist_input, true_block);
#ifdef DEBUG_PRINT
    true_block.print("  true_block=", "\n", true, std::cerr);
    hist_input.print("  hist_input=", "\n", true, std::cerr);
#endif
  }
};
