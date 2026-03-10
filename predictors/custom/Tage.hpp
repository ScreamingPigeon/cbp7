#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"
#include "../common.hpp"
#include "TageTable.hpp"

using namespace hcm;

// ============================================================================
// Constexpr Helpers
// ============================================================================

// Generate a uniform std::array of N elements, all set to V.
template <typename T, std::size_t N>
constexpr std::array<T, N> uniform_array(T v) {
  std::array<T, N> a{};
  for (std::size_t i = 0; i < N; i++)
    a[i] = v;
  return a;
}

// Constexpr power for geometric series (integer approximation).
constexpr double constexpr_pow(double base, double exp) {
  // Use repeated squaring with integer part, linear interpolation for fraction.
  // Good enough for geometric history length computation.
  if (exp == 0.0) return 1.0;
  if (base == 0.0) return 0.0;
  // exp(exp * ln(base)) via integer approximation
  // For our use case, just compute iteratively
  int int_exp = static_cast<int>(exp);
  double frac = exp - int_exp;
  double result = 1.0;
  for (int i = 0; i < int_exp; i++)
    result *= base;
  // Linear interpolation for fractional part: result * (1 + frac*(base-1))
  result *= (1.0 + frac * (base - 1.0));
  return result;
}

// Generate a geometric history length series of N elements from MIN to MAX.
// Longest history first (index 0 = MAX), matching geometric_folds convention.
template <std::size_t N>
constexpr std::array<u64, N> geometric_hist(u64 min_h, u64 max_h) {
  static_assert(N >= 2, "Need at least 2 tables for geometric series");
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

// Max element of a constexpr array.
template <typename T, std::size_t N>
constexpr T array_max(const std::array<T, N> &a) {
  T m = a[0];
  for (std::size_t i = 1; i < N; i++)
    m = (a[i] > m) ? a[i] : m;
  return m;
}

// Min element of a constexpr array.
template <typename T, std::size_t N>
constexpr T array_min(const std::array<T, N> &a) {
  T m = a[0];
  for (std::size_t i = 1; i < N; i++)
    m = (a[i] < m) ? a[i] : m;
  return m;
}

// ============================================================================
// Default Config Structs
// ============================================================================

// Per-table configuration. Arrays indexed by table number.
// Override by defining your own struct with the same interface.
struct DefaultTableConfig {
  static constexpr u64 NUM_TABLES = 8;
  static constexpr u64 MINHIST = 2;
  static constexpr u64 MAXHIST = 100;

  static constexpr auto TABLE_SIZE =
      uniform_array<u64, NUM_TABLES>(1 << 11); // 2048
  static constexpr auto TAG_WIDTH = uniform_array<u64, NUM_TABLES>(11);
  static constexpr auto CTR_WIDTH = uniform_array<u64, NUM_TABLES>(3);
  static constexpr auto U_WIDTH = uniform_array<u64, NUM_TABLES>(1);
  static constexpr auto HIST_LEN =
      geometric_hist<NUM_TABLES>(MINHIST, MAXHIST);
};

// Allocation policy configuration.
struct DefaultAllocConfig {
  static constexpr u64 MAX_ALLOC = 1;
  static constexpr bool USE_ALLOC_FILTER = false;
  static constexpr bool PROTECT_RECENT_ALLOC = false;
};

// ============================================================================
// Table Tuple Generation
// ============================================================================

// Generate a TageTable type from config arrays at index I, plus global params.
template <typename Cfg, u64 BR_P_ENTRY, u64 NUM_BANKS, bool USE_AHEAD,
          bool SHARED_TAG, bool SHARED_U, bool U_STOR_FF, u64 DECAY_CTR,
          typename ResetFn, bool USE_FF_CACHE, u64 I>
using TableAt =
    TageTable<Cfg::TABLE_SIZE[I], Cfg::HIST_LEN[I], Cfg::TAG_WIDTH[I],
              Cfg::CTR_WIDTH[I], Cfg::U_WIDTH[I], BR_P_ENTRY, NUM_BANKS,
              USE_AHEAD, SHARED_TAG, SHARED_U, U_STOR_FF, DECAY_CTR, ResetFn,
              USE_FF_CACHE>;

// Build a std::tuple of TageTable types from the config.
template <typename Cfg, u64 BR_P_ENTRY, u64 NUM_BANKS, bool USE_AHEAD,
          bool SHARED_TAG, bool SHARED_U, bool U_STOR_FF, u64 DECAY_CTR,
          typename ResetFn, bool USE_FF_CACHE, typename Seq>
struct MakeTableTuple;

template <typename Cfg, u64 BR_P_ENTRY, u64 NUM_BANKS, bool USE_AHEAD,
          bool SHARED_TAG, bool SHARED_U, bool U_STOR_FF, u64 DECAY_CTR,
          typename ResetFn, bool USE_FF_CACHE, u64... Is>
struct MakeTableTuple<Cfg, BR_P_ENTRY, NUM_BANKS, USE_AHEAD, SHARED_TAG,
                      SHARED_U, U_STOR_FF, DECAY_CTR, ResetFn, USE_FF_CACHE,
                      std::index_sequence<Is...>> {
  using type =
      std::tuple<TableAt<Cfg, BR_P_ENTRY, NUM_BANKS, USE_AHEAD, SHARED_TAG,
                         SHARED_U, U_STOR_FF, DECAY_CTR, ResetFn,
                         USE_FF_CACHE, Is>...>;
};

// ============================================================================
// Forward Declaration
// ============================================================================

template <bool USE_AHEAD_V, typename TableCfg, typename AllocCfg,
          // Global hardware params
          u64 FETCH_WIDTH_V, u64 BIMODAL_SIZE_V, u64 BR_P_ENTRY_V,
          u64 NUM_BANKS_V, bool SHARED_TAG_V, bool SHARED_U_V,
          bool U_STOR_FF_V, u64 DECAY_CTR_V, typename ResetFn_V,
          bool USE_FF_CACHE_V,
          // P1 params
          bool P1_USE_GSHARE_V, u64 P1_TABLE_SIZE_V, u64 P1_HIST_V,
          // Meta-prediction params
          bool USE_META_V, u64 METABITS_V, u64 METAPIPE_V,
          u64 META_TABLE_SIZE_V,
          // Path history params
          bool USE_PATH_HIST_V, u64 PATH_HIST_WIDTH_V, u64 PATH_BITS_V>
struct TageImpl;

// ============================================================================
// TageBase — shared computed constants and types for both specializations
// ============================================================================

template <typename TableCfg, typename AllocCfg, u64 FETCH_WIDTH_V,
          u64 BIMODAL_SIZE_V, u64 BR_P_ENTRY_V, u64 NUM_BANKS_V,
          bool USE_AHEAD_V, bool SHARED_TAG_V, bool SHARED_U_V,
          bool U_STOR_FF_V, u64 DECAY_CTR_V, typename ResetFn_V,
          bool USE_FF_CACHE_V, bool P1_USE_GSHARE_V, u64 P1_TABLE_SIZE_V,
          u64 P1_HIST_V, bool USE_META_V, u64 METABITS_V, u64 METAPIPE_V,
          u64 META_TABLE_SIZE_V, bool USE_PATH_HIST_V, u64 PATH_HIST_WIDTH_V,
          u64 PATH_BITS_V>
struct TageBase {
  // ---- Import config ----
  static constexpr u64 NUM_TABLES = TableCfg::NUM_TABLES;
  static constexpr u64 FETCH_WIDTH = FETCH_WIDTH_V;
  static constexpr u64 LOG_FETCH_WIDTH = clog2(FETCH_WIDTH);
  static constexpr u64 BR_P_ENTRY = BR_P_ENTRY_V;
  static constexpr u64 BIMODAL_SIZE = BIMODAL_SIZE_V;
  static constexpr u64 NUM_BANKS = NUM_BANKS_V;

  // ---- Per-table max widths (for uniform result registers) ----
  static constexpr u64 MAX_TAG_WIDTH = array_max(TableCfg::TAG_WIDTH);
  static constexpr u64 MAX_CTR_WIDTH = array_max(TableCfg::CTR_WIDTH);
  static constexpr u64 MAX_U_WIDTH = array_max(TableCfg::U_WIDTH);
  static constexpr u64 MAX_TABLE_SIZE = array_max(TableCfg::TABLE_SIZE);
  static constexpr u64 MAX_IDX_BITS = clog2(MAX_TABLE_SIZE);
  static constexpr u64 MAXHIST =
      array_max(TableCfg::HIST_LEN); // longest history

  // ---- Tag hashing constants ----
  // When BR_P_ENTRY=1, offset is encoded in tag:
  //   HTAGBITS = TAG_WIDTH - LOG_FETCH_WIDTH (hashed portion)
  //   full tag = concat(offset, htag)
  // When BR_P_ENTRY>1, full tag is hashed:
  //   HTAGBITS = TAG_WIDTH
  static constexpr u64 MAX_HTAGBITS =
      (BR_P_ENTRY == 1) ? (MAX_TAG_WIDTH - LOG_FETCH_WIDTH) : MAX_TAG_WIDTH;

  // ---- Bimodal constants ----
  static constexpr u64 LOG_BIMODAL_SIZE = clog2(BIMODAL_SIZE);
  static constexpr u64 BINDEX_BITS = LOG_BIMODAL_SIZE - LOG_FETCH_WIDTH;

  // ---- P1 constants ----
  static constexpr u64 P1_INDEX_BITS = clog2(P1_TABLE_SIZE_V) - LOG_FETCH_WIDTH;

  // ---- Meta constants ----
  static constexpr u64 METABITS = METABITS_V;
  static constexpr u64 METAPIPE = METAPIPE_V;
  static constexpr u64 META_TABLE_SIZE = META_TABLE_SIZE_V;

  // ---- Path history constants ----
  static constexpr u64 PATH_BITS = PATH_BITS_V;
  static constexpr u64 PATH_HIST_WIDTH = PATH_HIST_WIDTH_V;

  // ---- Table tuple type ----
  using Tables =
      typename MakeTableTuple<TableCfg, BR_P_ENTRY_V, NUM_BANKS_V,
                              USE_AHEAD_V, SHARED_TAG_V, SHARED_U_V,
                              U_STOR_FF_V, DECAY_CTR_V, ResetFn_V,
                              USE_FF_CACHE_V,
                              std::make_index_sequence<NUM_TABLES>>::type;

  // ---- Static asserts ----
  static_assert(NUM_TABLES > 0, "Need at least 1 table");
  static_assert(FETCH_WIDTH >= 1, "FETCH_WIDTH must be >= 1");
  static_assert(std::has_single_bit(FETCH_WIDTH),
                "FETCH_WIDTH must be power of 2");
  static_assert(BR_P_ENTRY == 1 || BR_P_ENTRY == FETCH_WIDTH,
                "BR_P_ENTRY must be 1 or FETCH_WIDTH");
  static_assert(BR_P_ENTRY == 1 || MAX_TAG_WIDTH > LOG_FETCH_WIDTH,
                "TAG_WIDTH must be > LOG_FETCH_WIDTH when BR_P_ENTRY=1");
  static_assert(BIMODAL_SIZE >= FETCH_WIDTH,
                "BIMODAL_SIZE must be >= FETCH_WIDTH");
  static_assert(P1_TABLE_SIZE_V >= FETCH_WIDTH,
                "P1_TABLE_SIZE must be >= FETCH_WIDTH");
};

// ============================================================================
// TageImpl — Direct (non-ahead) specialization
// ============================================================================

template <typename TableCfg, typename AllocCfg, u64 FETCH_WIDTH_V,
          u64 BIMODAL_SIZE_V, u64 BR_P_ENTRY_V, u64 NUM_BANKS_V,
          bool SHARED_TAG_V, bool SHARED_U_V, bool U_STOR_FF_V,
          u64 DECAY_CTR_V, typename ResetFn_V, bool USE_FF_CACHE_V,
          bool P1_USE_GSHARE_V, u64 P1_TABLE_SIZE_V, u64 P1_HIST_V,
          bool USE_META_V, u64 METABITS_V, u64 METAPIPE_V,
          u64 META_TABLE_SIZE_V, bool USE_PATH_HIST_V, u64 PATH_HIST_WIDTH_V,
          u64 PATH_BITS_V>
struct TageImpl<false, TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V,
                BR_P_ENTRY_V, NUM_BANKS_V, SHARED_TAG_V, SHARED_U_V,
                U_STOR_FF_V, DECAY_CTR_V, ResetFn_V, USE_FF_CACHE_V,
                P1_USE_GSHARE_V, P1_TABLE_SIZE_V, P1_HIST_V, USE_META_V,
                METABITS_V, METAPIPE_V, META_TABLE_SIZE_V, USE_PATH_HIST_V,
                PATH_HIST_WIDTH_V, PATH_BITS_V> : predictor {

  using Base = TageBase<TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V,
                        BR_P_ENTRY_V, NUM_BANKS_V, false, SHARED_TAG_V,
                        SHARED_U_V, U_STOR_FF_V, DECAY_CTR_V, ResetFn_V,
                        USE_FF_CACHE_V, P1_USE_GSHARE_V, P1_TABLE_SIZE_V,
                        P1_HIST_V, USE_META_V, METABITS_V, METAPIPE_V,
                        META_TABLE_SIZE_V, USE_PATH_HIST_V, PATH_HIST_WIDTH_V,
                        PATH_BITS_V>;

  static constexpr u64 NUM_TABLES = Base::NUM_TABLES;
  static constexpr u64 FETCH_WIDTH = Base::FETCH_WIDTH;
  static constexpr u64 LOG_FETCH_WIDTH = Base::LOG_FETCH_WIDTH;
  static constexpr u64 BR_P_ENTRY = Base::BR_P_ENTRY;
  static constexpr u64 MAX_TAG_WIDTH = Base::MAX_TAG_WIDTH;
  static constexpr u64 MAX_CTR_WIDTH = Base::MAX_CTR_WIDTH;
  static constexpr u64 MAX_U_WIDTH = Base::MAX_U_WIDTH;
  static constexpr u64 MAX_IDX_BITS = Base::MAX_IDX_BITS;
  static constexpr u64 MAX_HTAGBITS = Base::MAX_HTAGBITS;
  static constexpr u64 MAXHIST = Base::MAXHIST;
  static constexpr u64 BINDEX_BITS = Base::BINDEX_BITS;
  static constexpr u64 P1_INDEX_BITS = Base::P1_INDEX_BITS;

  // ---- Computed RAM sizes ----
  static constexpr u64 P1_ENTRIES = P1_TABLE_SIZE_V / FETCH_WIDTH;
  static constexpr u64 BIM_ENTRIES = BIMODAL_SIZE_V / FETCH_WIDTH;
  // Meta table size: at least 1 to keep ram<> valid even when conditional is false
  static constexpr u64 META_RAM_SIZE =
      (META_TABLE_SIZE_V > 0) ? META_TABLE_SIZE_V : 1;

  // ======== Table Instances ========
  typename Base::Tables tables;

  // ======== History State ========
  geometric_folds<NUM_TABLES, TableCfg::MINHIST, MAXHIST, MAX_IDX_BITS,
                  MAX_HTAGBITS>
      gfolds;
  reg<1> true_block = 1;

  // Path history (conditional)
  std::conditional_t<USE_PATH_HIST_V, reg<PATH_HIST_WIDTH_V>, EmptyMember>
      path_hist;

  // ======== P1 State ========
  // P1 uses its own short history (gshare mode only)
  std::conditional_t<P1_USE_GSHARE_V, reg<P1_HIST_V>, EmptyMember>
      global_history1;
  reg<P1_INDEX_BITS> index1;
  arr<reg<1>, FETCH_WIDTH> readp1; // P1 read results
  reg<FETCH_WIDTH> p1;             // P1 predictions (concatenated)

  // P1 tables
  hcm::ram<val<1>, P1_ENTRIES> table1_pred[FETCH_WIDTH]{"P1 pred"};

  // ======== Bimodal (P2 base) ========
  reg<BINDEX_BITS> bindex;
  arr<reg<1>, FETCH_WIDTH> readb; // bimodal read results
  hcm::ram<val<1>, BIM_ENTRIES> bim[FETCH_WIDTH]{"bpred"};

  // ======== P2 (TAGE) Result Registers ========
  // Uniform max-width to keep cross-table logic simple.
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> gindex;    // per-table indices
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> htag;      // per-table hashed tags
  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> readt;     // read tags
  arr<reg<1>, NUM_TABLES> readc;                 // read prediction bits
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> readh;     // read hysteresis
  arr<reg<1>, NUM_TABLES> readu;                 // read u-bits
  reg<NUM_TABLES> notumask;                      // inverted u-bits

  // ======== Match Registers ========
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match;  // all matches per offset
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match1; // longest match per offset
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match2; // second longest per offset
  arr<reg<1>, FETCH_WIDTH> pred1;               // primary prediction per offset
  arr<reg<1>, FETCH_WIDTH> pred2;               // alternate prediction per offset
  reg<FETCH_WIDTH> p2;                          // final P2 predictions

  // ======== Meta-prediction ========
  // USE_META=true: signed counter for alt-on-newly-allocated
  // META_TABLE_SIZE=0: global counter
  // META_TABLE_SIZE>0: PC-indexed table
  std::conditional_t<USE_META_V && META_TABLE_SIZE_V == 0,
                     arr<reg<METABITS_V, i64>, METAPIPE_V>, EmptyMember>
      meta_global;
  std::conditional_t<USE_META_V && (META_TABLE_SIZE_V > 0),
                     hcm::ram<val<METABITS_V, i64>, META_RAM_SIZE>,
                     EmptyMember>
      meta_table;
  std::conditional_t<USE_META_V, arr<reg<1>, FETCH_WIDTH>, EmptyMember>
      newly_alloc;

  // ======== U-bit Management ========
  // Epoch mode (U_STOR_FF=true): global reset counter
  static constexpr u64 UCTRBITS = 8;
  std::conditional_t<U_STOR_FF_V, reg<UCTRBITS>, EmptyMember> uctr;
  // SRAM decay mode: adaptive threshold (runtime register)
  // TODO: decay threshold register and allocation failure tracking

  // ======== UPDATE_ONLY Zone ========
  // Hysteresis tables don't affect prediction latency
  zone UPDATE_ONLY;
  hcm::ram<val<1>, P1_ENTRIES> table1_hyst[FETCH_WIDTH]{"P1 hyst"};
  hcm::ram<val<1>, BIM_ENTRIES> bhyst[FETCH_WIDTH]{"bhyst"};

  // ======== Block State (simulation artifacts) ========
  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<LOG_FETCH_WIDTH>, FETCH_WIDTH> branch_offset;
  arr<reg<1>, FETCH_WIDTH> branch_dir;
  reg<FETCH_WIDTH> block_entry; // one-hot vector

  // ======== Predictor Interface ========

  void new_block(val<64> inst_pc) {
    val<LOG_FETCH_WIDTH> offset = inst_pc.fo1() >> 2;
    block_entry = offset.fo1().decode().concat();
    block_size = 1;
  }

  val<1> predict1([[maybe_unused]] val<64> inst_pc) override {
    // TODO: Phase 3 — P1 gshare/bimodal predict
    new_block(inst_pc);
    return val<1>{0};
  }

  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) override {
    // TODO: Phase 3 — P1 reuse predict
    return val<1>{0};
  }

  val<1> predict2([[maybe_unused]] val<64> inst_pc) override {
    // TODO: Phase 2 — P2 TAGE predict
    return val<1>{0};
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) override {
    // TODO: Phase 2 — P2 TAGE reuse predict
    return val<1>{0};
  }

  void update_condbr([[maybe_unused]] val<64> branch_pc,
                     [[maybe_unused]] val<1> taken,
                     [[maybe_unused]] val<64> next_pc) override {
    // TODO: Phase 2 — record branch offset and direction
  }

  void update_cycle([[maybe_unused]] instruction_info &block_end_info) override {
    // TODO: Phase 2 — allocation, counter/u-bit update, history update
  }
};

// ============================================================================
// TageImpl — Ahead (1-branch-ahead) specialization
// ============================================================================

template <typename TableCfg, typename AllocCfg, u64 FETCH_WIDTH_V,
          u64 BIMODAL_SIZE_V, u64 BR_P_ENTRY_V, u64 NUM_BANKS_V,
          bool SHARED_TAG_V, bool SHARED_U_V, bool U_STOR_FF_V,
          u64 DECAY_CTR_V, typename ResetFn_V, bool USE_FF_CACHE_V,
          bool P1_USE_GSHARE_V, u64 P1_TABLE_SIZE_V, u64 P1_HIST_V,
          bool USE_META_V, u64 METABITS_V, u64 METAPIPE_V,
          u64 META_TABLE_SIZE_V, bool USE_PATH_HIST_V, u64 PATH_HIST_WIDTH_V,
          u64 PATH_BITS_V>
struct TageImpl<true, TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V,
                BR_P_ENTRY_V, NUM_BANKS_V, SHARED_TAG_V, SHARED_U_V,
                U_STOR_FF_V, DECAY_CTR_V, ResetFn_V, USE_FF_CACHE_V,
                P1_USE_GSHARE_V, P1_TABLE_SIZE_V, P1_HIST_V, USE_META_V,
                METABITS_V, METAPIPE_V, META_TABLE_SIZE_V, USE_PATH_HIST_V,
                PATH_HIST_WIDTH_V, PATH_BITS_V> : predictor {

  using Base = TageBase<TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V,
                        BR_P_ENTRY_V, NUM_BANKS_V, true, SHARED_TAG_V,
                        SHARED_U_V, U_STOR_FF_V, DECAY_CTR_V, ResetFn_V,
                        USE_FF_CACHE_V, P1_USE_GSHARE_V, P1_TABLE_SIZE_V,
                        P1_HIST_V, USE_META_V, METABITS_V, METAPIPE_V,
                        META_TABLE_SIZE_V, USE_PATH_HIST_V, PATH_HIST_WIDTH_V,
                        PATH_BITS_V>;

  static constexpr u64 NUM_TABLES = Base::NUM_TABLES;
  static constexpr u64 FETCH_WIDTH = Base::FETCH_WIDTH;
  static constexpr u64 LOG_FETCH_WIDTH = Base::LOG_FETCH_WIDTH;
  static constexpr u64 BR_P_ENTRY = Base::BR_P_ENTRY;
  static constexpr u64 MAX_TAG_WIDTH = Base::MAX_TAG_WIDTH;
  static constexpr u64 MAX_CTR_WIDTH = Base::MAX_CTR_WIDTH;
  static constexpr u64 MAX_U_WIDTH = Base::MAX_U_WIDTH;
  static constexpr u64 MAX_IDX_BITS = Base::MAX_IDX_BITS;
  static constexpr u64 MAX_HTAGBITS = Base::MAX_HTAGBITS;
  static constexpr u64 MAXHIST = Base::MAXHIST;
  static constexpr u64 BINDEX_BITS = Base::BINDEX_BITS;
  static constexpr u64 P1_INDEX_BITS = Base::P1_INDEX_BITS;

  // ---- Computed RAM sizes ----
  static constexpr u64 P1_ENTRIES = P1_TABLE_SIZE_V / FETCH_WIDTH;
  static constexpr u64 BIM_ENTRIES = BIMODAL_SIZE_V / FETCH_WIDTH;
  static constexpr u64 META_RAM_SIZE =
      (META_TABLE_SIZE_V > 0) ? META_TABLE_SIZE_V : 1;

  // ---- Ahead-specific constants ----
  static constexpr u64 PATHS = BR_P_ENTRY + 1;
  static constexpr u64 LOGPATHS = clog2(PATHS);
  static constexpr u64 LANES = FETCH_WIDTH;
  static constexpr u64 LOGLANES = (LANES > 1) ? clog2(LANES) : 1;

  // ======== Table Instances ========
  typename Base::Tables tables;

  // ======== History State ========
  geometric_folds<NUM_TABLES, TableCfg::MINHIST, MAXHIST, MAX_IDX_BITS,
                  MAX_HTAGBITS>
      gfolds;
  reg<1> true_block = 1;

  // Path history (conditional)
  std::conditional_t<USE_PATH_HIST_V, reg<PATH_HIST_WIDTH_V>, EmptyMember>
      path_hist;

  // ======== Ahead Pipeline Registers ========
  // [0] = current block, [1] = previous block
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> gindex[2]; // per-table indices
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> htag[2];   // per-table hashed tags
  reg<LOGPATHS> XB[2]; // address bits for path bank selection
  reg<LOGLANES> XL;    // address bits for lane selection
  reg<LOGPATHS> path;  // path taken out of previous block

  // Cached bank predictions: read in previous cycle, selected in current
  // TODO: size and structure depends on TageTable ahead banking layout

  // ======== P1 State ========
  std::conditional_t<P1_USE_GSHARE_V, reg<P1_HIST_V>, EmptyMember>
      global_history1;
  reg<P1_INDEX_BITS> index1[2]; // pipelined P1 index
  arr<reg<1>, FETCH_WIDTH> readp1;
  reg<FETCH_WIDTH> p1;

  // P1 tables
  hcm::ram<val<1>, P1_ENTRIES> table1_pred[FETCH_WIDTH]{"P1 pred"};

  // ======== Bimodal (P2 base) ========
  reg<BINDEX_BITS> bindex[2]; // pipelined bimodal index
  arr<reg<1>, FETCH_WIDTH> readb;
  hcm::ram<val<1>, BIM_ENTRIES> bim[FETCH_WIDTH]{"bpred"};

  // ======== P2 Result Registers ========
  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> readt;
  arr<reg<1>, NUM_TABLES> readc;
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> readh;
  arr<reg<1>, NUM_TABLES> readu;
  reg<NUM_TABLES> notumask;

  // ======== Match Registers ========
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match;
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match1;
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match2;
  arr<reg<1>, FETCH_WIDTH> pred1;
  arr<reg<1>, FETCH_WIDTH> pred2;
  reg<FETCH_WIDTH> p2;

  // ======== Meta-prediction ========
  std::conditional_t<USE_META_V && META_TABLE_SIZE_V == 0,
                     arr<reg<METABITS_V, i64>, METAPIPE_V>, EmptyMember>
      meta_global;
  std::conditional_t<USE_META_V && (META_TABLE_SIZE_V > 0),
                     hcm::ram<val<METABITS_V, i64>, META_RAM_SIZE>,
                     EmptyMember>
      meta_table;
  std::conditional_t<USE_META_V, arr<reg<1>, FETCH_WIDTH>, EmptyMember>
      newly_alloc;

  // ======== U-bit Management ========
  static constexpr u64 UCTRBITS = 8;
  std::conditional_t<U_STOR_FF_V, reg<UCTRBITS>, EmptyMember> uctr;

  // ======== UPDATE_ONLY Zone ========
  zone UPDATE_ONLY;
  hcm::ram<val<1>, P1_ENTRIES> table1_hyst[FETCH_WIDTH]{"P1 hyst"};
  hcm::ram<val<1>, BIM_ENTRIES> bhyst[FETCH_WIDTH]{"bhyst"};

  // ======== Block State ========
  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<LOG_FETCH_WIDTH>, FETCH_WIDTH> branch_offset;
  arr<reg<1>, FETCH_WIDTH> branch_dir;
  reg<FETCH_WIDTH> block_entry;

  // ======== Predictor Interface ========

  void new_block(val<64> inst_pc) {
    val<LOG_FETCH_WIDTH> offset = inst_pc.fo1() >> 2;
    block_entry = offset.fo1().decode().concat();
    block_size = 1;
  }

  val<1> predict1([[maybe_unused]] val<64> inst_pc) override {
    // TODO: Phase 5.5 — ahead P1 predict
    new_block(inst_pc);
    return val<1>{0};
  }

  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) override {
    // TODO: Phase 5.5 — ahead P1 reuse predict
    return val<1>{0};
  }

  val<1> predict2([[maybe_unused]] val<64> inst_pc) override {
    // TODO: Phase 5.5 — ahead P2 predict (select from cached bank reads)
    return val<1>{0};
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) override {
    // TODO: Phase 5.5 — ahead P2 reuse predict
    return val<1>{0};
  }

  void update_condbr([[maybe_unused]] val<64> branch_pc,
                     [[maybe_unused]] val<1> taken,
                     [[maybe_unused]] val<64> next_pc) override {
    // TODO: Phase 5.5 — ahead update_condbr
  }

  void update_cycle([[maybe_unused]] instruction_info &block_end_info) override {
    // TODO: Phase 5.5 — ahead update_cycle (write to index[1], path banking)
    // TODO: Phase 8 — loop predictor override
  }
};

// ============================================================================
// User-facing Alias
// ============================================================================

template <typename TableCfg = DefaultTableConfig,
          typename AllocCfg = DefaultAllocConfig,
          // Global hardware params
          u64 FETCH_WIDTH = 8, u64 BIMODAL_SIZE = 4096, u64 BR_P_ENTRY = 1,
          u64 NUM_BANKS = 1, bool USE_AHEAD = false, bool SHARED_TAG = true,
          bool SHARED_U = true, bool U_STOR_FF = false, u64 DECAY_CTR = 1024,
          typename ResetFn = DefaultResetFn, bool USE_FF_CACHE = false,
          // P1 params
          bool P1_USE_GSHARE = true, u64 P1_TABLE_SIZE = 16384, u64 P1_HIST = 6,
          // Meta-prediction
          bool USE_META = true, u64 METABITS = 4, u64 METAPIPE = 2,
          u64 META_TABLE_SIZE = 0,
          // Path history
          bool USE_PATH_HIST = false, u64 PATH_HIST_WIDTH = 27,
          u64 PATH_BITS = 6>
using Tage =
    TageImpl<USE_AHEAD, TableCfg, AllocCfg, FETCH_WIDTH, BIMODAL_SIZE,
             BR_P_ENTRY, NUM_BANKS, SHARED_TAG, SHARED_U, U_STOR_FF,
             DECAY_CTR, ResetFn, USE_FF_CACHE, P1_USE_GSHARE, P1_TABLE_SIZE,
             P1_HIST, USE_META, METABITS, METAPIPE, META_TABLE_SIZE,
             USE_PATH_HIST, PATH_HIST_WIDTH, PATH_BITS>;
