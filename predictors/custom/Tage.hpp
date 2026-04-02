// Suppress warn
#include <memory>
#pragma GCC diagnostic ignored "-Wunused-parameter"

#ifndef TAGE_HPP
#define TAGE_HPP

#include "../../cbp.hpp"
#include "../../harcom.hpp"
#include "TageTable.hpp"

using namespace hcm;

// Empty type for conditional compilation
struct Empty {};

/* ========== PREDICTOR-LEVEL PARAMETERS ========== */
/*
 * PRED_BLK_SIZE:     (u64) Size of fetch block
 *
 * TODO:
 * BASE_PRED:         (enum) Which base predictor to use for p0
 * BASE_PRED_SIZE:    (u64) Number of entries in base predictor
 * BASE_PRED_HYST:    (bool) Whether to use shared hysteresis
 *                    bits for base predictor
 * BASE_HYST_SIZE:    (u64) Size of shared hysteresis table
 *                    if enabled
 *
 * NUM_TABLES:        (u64) Number of TAGE_TABLES
 *                    (typically 7-14)
 * T_TABLE_HIST_LEN:  (const u64[]) Array specifying history
 *                    length for each table
 * T_TABLE_SIZE:      (const u64[]) Array specifying number of
 *                    entries for each table
 *
 * GLOBAL_HIST_LEN:   (u64) Maximum global history length
 *                    (max of all T_TABLE_HIST_LEN)
 * PATH_WIDTH:        (u64) size of PC folded into path
 * PATH_HIST_LEN:     (u64) Path history length
 *                    (typically 16-27 bits)
 * USE_PATH_HIST:     (bool) Whether to use path history in
 *                    indexing
 * SPLIT_KSPACE:      use a different history for kernel hist
 *
 * TODO: Probabilistic allocation failure decay: Keep a counter that tracks
 * allocation rate. If allocation failure is high, increase probability for
 * usefulness decay on tag mismatch, restore once it is better
 *
 * TODO:
 * USE_ALT_ON_NA:     (bool) Whether to use alternate prediction
 *                    on newly allocated entries
 * USE_ALT_COUNTERS:  (u64) Number of USE_ALT_ON_NA counters

 * TODO:
 * ALLOC_POLICY:      (enum) Allocation policy on misprediction
 * RESET_PERIOD:      (u64) Period for resetting useful counters
 * RESET_POLICY:      (enum) Policy for useful counter reset
 */

/* ========== TAGE TABLE-LEVEL PARAMETERS ========== */
/*
 * T_TAG_WIDTH:       (const u64[]) Tag width in bits for each
 *                    table
 * T_CTR_WIDTH:       (u64) Prediction counter width
 *
 * T_U_WIDTH:         (u64) Useful counter width
 *
 * DECAY_BITS:        (u64) bits to randomly gener
 *
 * TODO:T_HASH_FUNC:  (enum[]) Hash function type for each table
 *                    (XOR, folded-XOR, etc.)
 */

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
  static constexpr auto CTR_WIDTH = uniform_array<u64, NUM_TABLES>(1);
  static constexpr auto HYST_WIDTH = uniform_array<u64, NUM_TABLES>(2);
  static constexpr auto U_WIDTH = uniform_array<u64, NUM_TABLES>(1);
  static constexpr auto HIST_LEN = geometric_hist<NUM_TABLES>(MINHIST, MAXHIST);
};

// Constexpr helper: generate per-table sizes via a callable.
template <std::size_t N, typename Fn>
constexpr std::array<u64, N> generate_table_sizes(Fn fn) {
  std::array<u64, N> s{};
  for (std::size_t i = 0; i < N; i++)
    s[i] = fn(i, N);
  return s;
}

// Size lambda: geometric scaling — short-history tables larger, long-history
// smaller. Index 0 = longest history (smallest table), index N-1 = shortest
// history (largest table). SIZE_RATIO=1: uniform. SIZE_RATIO=2: ~2x range from
// smallest to largest.
template <u64 SIZE, u64 SIZE_RATIO>
constexpr auto size_fn = [](u64 i, u64 n) -> u64 {
  if constexpr (SIZE_RATIO <= 1)
    return SIZE;
  else {
    double t = double(i) / std::max(1.0, double(n - 1));
    // i=0 (longest history) → scale < 1 (smaller); i=N-1 (shortest) → scale > 1
    // (larger)
    double scale = constexpr_pow(double(SIZE_RATIO), t - 0.5);
    u64 sz = u64(SIZE * scale);
    // Round up to power of 2, minimum 64
    u64 result = 64;
    while (result < sz)
      result *= 2;
    return result;
  }
};

// Sweep-friendly table config: uniform values with optional geometric size
// scaling.
template <u64 N = 8, u64 SIZE = 2048, u64 TAG = 11, u64 CTR = 1, u64 HYST = 2,
          u64 U = 1, u64 MINH = 2, u64 MAXH = 100, u64 SIZE_RATIO = 1>
struct SweepTableConfig {
  static constexpr u64 NUM_TABLES = N;
  static constexpr u64 MINHIST = MINH;
  static constexpr u64 MAXHIST = MAXH;
  static constexpr auto TABLE_SIZE =
      generate_table_sizes<N>(size_fn<SIZE, SIZE_RATIO>);
  static constexpr auto TAG_WIDTH = uniform_array<u64, N>(TAG);
  static constexpr auto CTR_WIDTH = uniform_array<u64, N>(CTR);
  static constexpr auto HYST_WIDTH = uniform_array<u64, N>(HYST);
  static constexpr auto U_WIDTH = uniform_array<u64, N>(U);
  static constexpr auto HIST_LEN = geometric_hist<N>(MINH, MAXH);
};

// Allocation policy configuration.
struct DefaultAllocConfig {
  static constexpr u64 MAX_ALLOC = 1;
  static constexpr bool USE_ALLOC_FILTER = false;
  static constexpr bool PROTECT_RECENT_ALLOC = false;
};

// ============================================================================
// Decay Threshold Adaptation Policies
// ============================================================================
// Functors controlling how the probabilistic u-bit decay threshold adapts.
// apply() is called once per threshold_tick; returns the new threshold value.
// Signals: correct = correct prediction, uctrsat = uctr saturated,
//          misp = misprediction (= ~correct).

// Mode 0: Original — decrement on uctrsat only, increment on correct
struct DecayConservative {
  template <u64 W>
  static val<W> apply(val<W> &t, val<1> &correct, val<1> &uctrsat,
                      [[maybe_unused]] val<1> &misp) {
    return select(uctrsat, update_ctr(t, hard<0>{}),
                  select(correct, update_ctr(t, hard<1>{}), val<W>{t}));
  }
};

// Mode 1: Decrement on misprediction, increment on correct
struct DecayMild {
  template <u64 W>
  static val<W> apply(val<W> &t, val<1> &correct,
                      [[maybe_unused]] val<1> &uctrsat, val<1> &misp) {
    return select(misp, update_ctr(t, hard<0>{}),
                  select(correct, update_ctr(t, hard<1>{}), val<W>{t}));
  }
};

// Mode 2: Decrement on misprediction, no increment (aggressive)
struct DecayAggressive {
  template <u64 W>
  static val<W> apply(val<W> &t, [[maybe_unused]] val<1> &correct,
                      [[maybe_unused]] val<1> &uctrsat, val<1> &misp) {
    return select(misp, update_ctr(t, hard<0>{}), val<W>{t});
  }
};

// Mode 3: Decrement on misprediction OR uctrsat, increment on correct
struct DecayHybrid {
  template <u64 W>
  static val<W> apply(val<W> &t, val<1> &correct, val<1> &uctrsat,
                      val<1> &misp) {
    return select(misp | uctrsat, update_ctr(t, hard<0>{}),
                  select(correct, update_ctr(t, hard<1>{}), val<W>{t}));
  }
};

// Mode 4: Always decrement (threshold races to 0, maximum decay)
struct DecayMax {
  template <u64 W>
  static val<W> apply(val<W> &t, [[maybe_unused]] val<1> &correct,
                      [[maybe_unused]] val<1> &uctrsat,
                      [[maybe_unused]] val<1> &misp) {
    return update_ctr(t, hard<0>{});
  }
};

// ============================================================================
// Table Tuple Generation
// ============================================================================

// Generate a TageTable type from config arrays at index I, plus global params.
template <typename Cfg, u64 BR_P_ENTRY, u64 NUM_BANKS, bool USE_AHEAD,
          bool SHARED_TAG, bool SHARED_U, bool SHARED_HYS, bool U_STOR_FF,
          u64 DECAY_CTR, typename ResetFn, bool USE_FF_CACHE, std::size_t I>
using TableAt =
    TageTable<Cfg::TABLE_SIZE[I], Cfg::HIST_LEN[I], Cfg::TAG_WIDTH[I],
              Cfg::CTR_WIDTH[I], Cfg::HYST_WIDTH[I], Cfg::U_WIDTH[I],
              BR_P_ENTRY, NUM_BANKS, USE_AHEAD, SHARED_TAG, SHARED_U,
              SHARED_HYS, U_STOR_FF, DECAY_CTR, ResetFn, USE_FF_CACHE>;

// Build a std::tuple of TageTable types from the config.
template <typename Cfg, u64 BR_P_ENTRY, u64 NUM_BANKS, bool USE_AHEAD,
          bool SHARED_TAG, bool SHARED_U, bool SHARED_HYS, bool U_STOR_FF,
          u64 DECAY_CTR, typename ResetFn, bool USE_FF_CACHE, typename Seq>
struct MakeTableTuple;

template <typename Cfg, u64 BR_P_ENTRY, u64 NUM_BANKS, bool USE_AHEAD,
          bool SHARED_TAG, bool SHARED_U, bool SHARED_HYS, bool U_STOR_FF,
          u64 DECAY_CTR, typename ResetFn, bool USE_FF_CACHE, std::size_t... Is>
struct MakeTableTuple<Cfg, BR_P_ENTRY, NUM_BANKS, USE_AHEAD, SHARED_TAG,
                      SHARED_U, SHARED_HYS, U_STOR_FF, DECAY_CTR, ResetFn,
                      USE_FF_CACHE, std::index_sequence<Is...>> {
  using type = std::tuple<
      TableAt<Cfg, BR_P_ENTRY, NUM_BANKS, USE_AHEAD, SHARED_TAG, SHARED_U,
              SHARED_HYS, U_STOR_FF, DECAY_CTR, ResetFn, USE_FF_CACHE, Is>...>;
};

// ============================================================================
// Forward Declaration
// ============================================================================

template <
    // Tage Table Default params
    u64 T_TAG_WIDTH = 7, u64 T_U_WIDTH = 2, u64 T_CTR_WIDTH = 3,
    u64 T_HYS_WIDTH = 2, bool T_USE_HYS = true, u64 T_DECAY_CTR = 1024,
    // Base Predictor Patameters
    // TODO

// ============================================================================
// TageBase — shared computed constants and types for both specializations
// ============================================================================

template <
    typename TableCfg, typename AllocCfg, u64 FETCH_WIDTH_V, u64 BIMODAL_SIZE_V,
    u64 BR_P_ENTRY_V, u64 NUM_BANKS_V, bool USE_AHEAD_V, bool SHARED_TAG_V,
    bool SHARED_U_V, bool SHARED_HYS_V, bool U_STOR_FF_V, u64 DECAY_CTR_V,
    u64 DECAY_GRAN_V, typename ResetFn_V, bool USE_FF_CACHE_V,
    bool P1_USE_GSHARE_V, u64 P1_TABLE_SIZE_V, u64 P1_HIST_V, bool USE_META_V,
    u64 METABITS_V, u64 METAPIPE_V, u64 META_TABLE_SIZE_V, bool USE_PATH_HIST_V,
    u64 PATH_HIST_WIDTH_V, u64 PATH_BITS_V, typename Overrider_V>
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
  static constexpr u64 MAX_HYST_WIDTH = array_max(TableCfg::HYST_WIDTH);
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
  using Tables = typename MakeTableTuple<
      TableCfg, BR_P_ENTRY_V, NUM_BANKS_V, USE_AHEAD_V, SHARED_TAG_V,
      SHARED_U_V, SHARED_HYS_V, U_STOR_FF_V, DECAY_CTR_V, ResetFn_V,
      USE_FF_CACHE_V, std::make_index_sequence<NUM_TABLES>>::type;

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
// Flat inline implementation matching tage_tt structure for optimal
// EPI/latency. All prediction and update logic is inline — no member function
// decomposition.

template <
    typename TableCfg, typename AllocCfg, u64 FETCH_WIDTH_V, u64 BIMODAL_SIZE_V,
    u64 BR_P_ENTRY_V, u64 NUM_BANKS_V, bool SHARED_TAG_V, bool SHARED_U_V,
    bool SHARED_HYS_V, bool U_STOR_FF_V, u64 DECAY_CTR_V, u64 DECAY_GRAN_V,
    typename DecayPolicy_V, typename ResetFn_V, bool USE_FF_CACHE_V,
    bool P1_USE_GSHARE_V, u64 P1_TABLE_SIZE_V, u64 P1_HIST_V, bool USE_META_V,
    u64 METABITS_V, u64 METAPIPE_V, u64 META_TABLE_SIZE_V, bool USE_PATH_HIST_V,
    u64 PATH_HIST_WIDTH_V, u64 PATH_BITS_V, typename Overrider_V>
struct TageImpl<false, TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V,
                BR_P_ENTRY_V, NUM_BANKS_V, SHARED_TAG_V, SHARED_U_V,
                SHARED_HYS_V, U_STOR_FF_V, DECAY_CTR_V, DECAY_GRAN_V,
                DecayPolicy_V, ResetFn_V, USE_FF_CACHE_V, P1_USE_GSHARE_V,
                P1_TABLE_SIZE_V, P1_HIST_V, USE_META_V, METABITS_V, METAPIPE_V,
                META_TABLE_SIZE_V, USE_PATH_HIST_V, PATH_HIST_WIDTH_V,
                PATH_BITS_V, Overrider_V> : predictor {

  using Overrider = Overrider_V;
  static constexpr u64 OVR = Overrider::ENABLED ? 1 : 0;

  using Base =
      TageBase<TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V, BR_P_ENTRY_V,
               NUM_BANKS_V, false, SHARED_TAG_V, SHARED_U_V, SHARED_HYS_V,
               U_STOR_FF_V, DECAY_CTR_V, DECAY_GRAN_V, ResetFn_V,
               USE_FF_CACHE_V, P1_USE_GSHARE_V, P1_TABLE_SIZE_V, P1_HIST_V,
               USE_META_V, METABITS_V, METAPIPE_V, META_TABLE_SIZE_V,
               USE_PATH_HIST_V, PATH_HIST_WIDTH_V, PATH_BITS_V, Overrider_V>;

  // ======== Constants ========
  static constexpr u64 NUM_TABLES = Base::NUM_TABLES;
  static constexpr u64 FETCH_WIDTH = Base::FETCH_WIDTH;
  static constexpr u64 LOG_FETCH_WIDTH = Base::LOG_FETCH_WIDTH;
  static constexpr u64 MAX_TAG_WIDTH = Base::MAX_TAG_WIDTH;
  static constexpr u64 MAX_CTR_WIDTH = Base::MAX_CTR_WIDTH;
  static constexpr u64 MAX_HYST_WIDTH = Base::MAX_HYST_WIDTH;
  static constexpr u64 MAX_U_WIDTH = Base::MAX_U_WIDTH;
  static constexpr u64 MAX_IDX_BITS = Base::MAX_IDX_BITS;
  static constexpr u64 MAX_HTAGBITS = Base::MAX_HTAGBITS;
  static constexpr u64 MAXHIST = Base::MAXHIST;
  static constexpr u64 BINDEX_BITS = Base::BINDEX_BITS;
  static constexpr u64 P1_INDEX_BITS = Base::P1_INDEX_BITS;
  static constexpr u64 UCTRBITS = 8;
  static constexpr u64 PATHBITS = Base::PATH_BITS;
  static constexpr u64 LINEADDR_BITS = std::max(BINDEX_BITS, MAX_IDX_BITS);
  static constexpr u64 P1_ENTRIES = P1_TABLE_SIZE_V / FETCH_WIDTH;
  static constexpr u64 BIM_ENTRIES = BIMODAL_SIZE_V / FETCH_WIDTH;
  // Unhashed wordline bits: low ROW_DIRECT bits of gindex come straight from
  // lineaddr (no XOR with history). Only upper bits are hashed. Removes the
  // XOR gate delay from SRAM row-select critical path (EV8 technique).
  static constexpr u64 ROW_DIRECT = std::min(u64(6), MAX_IDX_BITS);

  // Truncate gindex to per-table IDX_BITS (needed when size_ratio > 1)
  template <u64 I> auto tidx(auto &gi) {
    using Table = std::tuple_element_t<I, typename Base::Tables>;
    return val<Table::IDX_BITS>{gi};
  }

  // ======== State ========
  geometric_folds<NUM_TABLES, TableCfg::MINHIST, MAXHIST, MAX_IDX_BITS,
                  MAX_HTAGBITS>
      gfolds;
  reg<1> true_block = 1;

  // P1
  std::conditional_t<P1_USE_GSHARE_V, reg<P1_HIST_V>, EmptyMember>
      global_history1;
  reg<P1_INDEX_BITS> index1;
  arr<reg<1>, FETCH_WIDTH> readp1;
  reg<FETCH_WIDTH> p1;

  // P2
  reg<BINDEX_BITS> bindex;
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> gindex;
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> htag;

  arr<reg<1>, FETCH_WIDTH> readb;
  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> readt;
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> readc;
  arr<reg<std::max(u64(1), MAX_HYST_WIDTH)>, NUM_TABLES> readh;
  arr<reg<MAX_U_WIDTH>, NUM_TABLES> readu;
  reg<NUM_TABLES> notumask;

  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match;
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match1;
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match2;

  arr<reg<1>, FETCH_WIDTH> pred1;
  arr<reg<1>, FETCH_WIDTH> pred2;
  reg<FETCH_WIDTH> p2;

  // Meta-prediction
  std::conditional_t<USE_META_V, arr<reg<METABITS_V, i64>, METAPIPE_V>,
                     EmptyMember>
      meta;
  std::conditional_t<USE_META_V, arr<reg<1>, FETCH_WIDTH>, EmptyMember>
      newly_alloc;

  // U-bit reset
  reg<UCTRBITS> uctr;

  // Probabilistic u-bit decay threshold (SRAM mode only, DECAY_CTR_V > 0)
  static constexpr bool USE_PROB_DECAY = (DECAY_CTR_V > 0);
  std::conditional_t<USE_PROB_DECAY, reg<DECAY_CTR_V == 0 ? 1 : DECAY_CTR_V>,
                     EmptyMember>
      decay_threshold;

  // Path history
  std::conditional_t<USE_PATH_HIST_V, reg<PATH_HIST_WIDTH_V>, EmptyMember>
      path_hist;

  // Simulation artifacts
  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<LOG_FETCH_WIDTH>, FETCH_WIDTH> branch_offset;
  arr<reg<1>, FETCH_WIDTH> branch_dir;
  reg<FETCH_WIDTH> block_entry;

  // P1 tables
  hcm::ram<val<1>, P1_ENTRIES> table1_pred[FETCH_WIDTH]{"P1 pred"};

  // P2 TAGE tables (TageTable storage — tuple for non-uniform params)
  typename Base::Tables tables;

  // P2 bimodal
  hcm::ram<val<1>, BIM_ENTRIES> bim[FETCH_WIDTH]{"bpred"};

  zone UPDATE_ONLY;
  hcm::ram<val<1>, P1_ENTRIES> table1_hyst[FETCH_WIDTH]{"P1 hyst"};
  hcm::ram<val<1>, BIM_ENTRIES> bhyst[FETCH_WIDTH]{"bhyst"};

  bool params_printed = false;
  void print_params(std::ostream &os) const {
    os << "=== TAGE Parameters ===\n";
    os << "NUM_TABLES=" << NUM_TABLES << "  FETCH_WIDTH=" << FETCH_WIDTH
       << "  BR_P_ENTRY=" << BR_P_ENTRY_V << "  NUM_BANKS=" << NUM_BANKS_V
       << "\n";
    os << "BIMODAL_SIZE=" << BIMODAL_SIZE_V << "  USE_AHEAD=false\n";
    os << "Table sizes: ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::TABLE_SIZE[i];
    os << "\nTag widths:  ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::TAG_WIDTH[i];
    os << "\nCtr widths:  ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::CTR_WIDTH[i];
    os << "\nHyst widths: ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::HYST_WIDTH[i];
    os << "\nU widths:    ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::U_WIDTH[i];
    os << "\nHist lens:   ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::HIST_LEN[i];
    os << "\n";
    os << "SHARED_TAG=" << SHARED_TAG_V << "  SHARED_U=" << SHARED_U_V
       << "  SHARED_HYS=" << SHARED_HYS_V << "  U_STOR_FF=" << U_STOR_FF_V
       << "\n";
    os << "DECAY_CTR=" << DECAY_CTR_V << "  DECAY_GRAN=" << DECAY_GRAN_V
       << "  USE_FF_CACHE=" << USE_FF_CACHE_V << "\n";
    os << "P1: USE_GSHARE=" << P1_USE_GSHARE_V
       << "  TABLE_SIZE=" << P1_TABLE_SIZE_V << "  HIST=" << P1_HIST_V << "\n";
    os << "META: USE=" << USE_META_V << "  BITS=" << METABITS_V
       << "  PIPE=" << METAPIPE_V << "  TABLE_SIZE=" << META_TABLE_SIZE_V
       << "\n";
    os << "PATH: USE=" << USE_PATH_HIST_V << "  WIDTH=" << PATH_HIST_WIDTH_V
       << "  BITS=" << PATH_BITS_V << "\n";
    os << "Overrider: " << (Overrider::ENABLED ? "enabled" : "disabled")
       << "\n";
    os << "\n";
  }

#ifdef TAGE_MONITOR
  TageMonitor<NUM_TABLES, Base::MAX_TABLE_SIZE, FETCH_WIDTH> mon;
  bool mon_init = false;
  void ensure_monitor_init() {
    if (!mon_init) {
      std::array<u64, NUM_TABLES> sizes{};
      for (u64 i = 0; i < NUM_TABLES; i++)
        sizes[i] = TableCfg::TABLE_SIZE[i];
      mon.set_table_sizes(sizes);
      mon.set_p1_size(P1_ENTRIES);
      mon_init = true;
    }
  }
  ~TageImpl() {
    if constexpr (Overrider::ENABLED) {
      if constexpr (Overrider::IS_LOOP) {
        mon.cum.loop_lookups = overrider.stats.lookups;
        mon.cum.loop_hits = overrider.stats.hits;
        mon.cum.loop_overrides = overrider.stats.overrides;
      }
      if constexpr (Overrider::IS_SC) {
        mon.cum.sc_lookups = overrider.stats.lookups;
        mon.cum.sc_overrides = overrider.stats.overrides;
        mon.cum.sc_updates = overrider.stats.update_writes;
        mon.cum.sc_correct = overrider.stats.sc_correct;
        mon.cum.sc_wrong = overrider.stats.sc_wrong;
      }
    }
    mon.print_summary(std::cerr);
    if constexpr (Overrider::ENABLED) {
      if constexpr (requires { overrider.print_stats(); }) {
        overrider.print_stats();
      }
    }
  }
#endif

  // ======== Overrider (conditional) ========
  u64 last_raw_pc = 0; // SW shadow of last inst_pc for cache keying
  std::conditional_t<Overrider::ENABLED, Overrider, EmptyMember> overrider;
  // (overrider lookup now in P2, no cross-cycle regs needed)
  std::conditional_t<Overrider::ENABLED, arr<reg<1>, FETCH_WIDTH>, EmptyMember>
      override_active;
  std::conditional_t<Overrider::ENABLED, reg<FETCH_WIDTH>, EmptyMember>
      p2_before_override;

  // ======== Predictor Interface ========

public:
  // Interface to simulator
  void new_block(val<64> inst_pc) {
    base_addr = 0; // TODO: right shift inst_pc by 64 - BLOCK_BITS to get block
                   // bits set in base_addr
    // TODO: Compute tag and index
  }
  val<1> predict(val<64> inst_pc) {
    // Compute newRead and return prediction for slot
    return val<1>(0);
  }
  val<1> reuse_predict(val<64> inst_pc) {
    // use cached read and return prediction for slot
    return val<1>(0);
  }

  void update_condbr(val<64> branch_pc, val<1> taken, val<64> next_pc) {
    // Update registers
  }
  void update_cycle(instruction_info &block_end_info) {
    // updateRAM
  }

private:
  // TageTableTuple
  template <std::size_t... Is>
  static auto make_tables_impl(std::index_sequence<Is...>) {
    return std::make_tuple(
        std::make_unique<TageTable<std::get<Is>(T_TABLE_HIST_SIZE),
                                   std::get<Is>(T_TABLE_HIST_LEN), T_TAG_WIDTH,
                                   T_CTR_WIDTH, T_U_WIDTH, PRED_BLK_SIZE,
                                   T_DECAY_CTR, EN_N_BLK_RD>>()...);
  }

  // Tuple of tables
  decltype(make_tables_impl(
      std::make_index_sequence<NUM_TABLES>{})) tage_tables;

  // State registers
  reg<GLOBAL_HIST_LEN> br_hist_reg;
  std::conditional_t<SPLIT_KSPACE, reg<GLOBAL_HIST_LEN>, Empty> br_khist_reg;
  std::conditional_t<USE_PATH_HIST, reg<PATH_HIST_LEN>, Empty> p_hist_reg;
  std::conditional_t<USE_PATH_HIST & SPLIT_KSPACE, reg<PATH_HIST_LEN>, Empty>
      p_khist_reg;

  reg<BLOCK_BITS> base_addr;

    // read tables — start RAM reads BEFORE tag hash (tag not needed until
    // compare)
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      readb[offset] = bim[offset].read(bindex);
    }
    readb.fanout(hard<2>{});
    static_loop<NUM_TABLES>([&]<u64 I>() {
      auto &t = std::get<I>(tables);
      readt[I] = t.tag_ram[0].read(tidx<I>(gindex[I]));
      readc[I] = t.pred_ram[0].read(tidx<I>(gindex[I]));
      if constexpr (MAX_HYST_WIDTH > 0) {
        readh[I] = t.hyst_ram[0].read(tidx<I>(gindex[I]));
      }
      if constexpr (U_STOR_FF_V) {
        readu[I] = t.u_ff[0].select(tidx<I>(gindex[I]));
      } else {
        readu[I] = t.u_ram[0].read(tidx<I>(gindex[I]));
      }
    });
    readt.fanout(hard<FETCH_WIDTH + 1>{});
    readc.fanout(hard<3>{});
    if constexpr (MAX_HYST_WIDTH > 0) {
      readh.fanout(hard<2>{});
    }
    readu.fanout(hard<2>{});
    notumask = ~readu.concat();
    notumask.fanout(hard<2>{});

    // compute hashed tags — parallel with RAM reads (not needed until compare)
    for (u64 i = 0; i < NUM_TABLES; i++) {
      htag[i] =
          val<MAX_HTAGBITS>{lineaddr}.reverse() ^ gfolds.template get<1>(i);
    }
    htag.fanout(hard<2>{});

    // gather prediction bits for each offset
    val<NUM_TABLES> gpreds = [&]() -> val<NUM_TABLES> {
      if constexpr (MAX_CTR_WIDTH == 1) {
        return readc.concat();
      } else {
        arr<val<1>, NUM_TABLES> gpred_split = [&](int i) -> val<1> {
          return readc[i] >> hard<MAX_CTR_WIDTH - 1>{};
        };
        return gpred_split.fo1().concat();
      }
    }();
    gpreds.fanout(hard<FETCH_WIDTH>{});
    arr<val<NUM_TABLES + 1>, FETCH_WIDTH> preds = [&](u64 offset) {
      return concat(readb[offset], gpreds);
    };
    preds.fanout(hard<2 * FETCH_WIDTH>{});

    // tag comparisons
    if constexpr (BR_P_ENTRY_V == 1) {
      // offset encoded in upper tag bits
      arr<val<1>, NUM_TABLES> htagcmp_split = [&](int i) {
        return val<MAX_HTAGBITS>{readt[i]} == htag[i];
      };
      val<NUM_TABLES> htagcmp = htagcmp_split.fo1().concat();
      htagcmp.fanout(hard<FETCH_WIDTH>{});

      static_loop<FETCH_WIDTH>([&]<u64 offset>() {
        arr<val<1>, NUM_TABLES> tagcmp = [&](int i) {
          return val<LOG_FETCH_WIDTH>{readt[i] >> MAX_HTAGBITS} ==
                 hard<offset>{};
        };
        match[offset] = concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp);
      });
    } else {
      // full tag match, same for all offsets
      arr<val<1>, NUM_TABLES> tagcmp_split = [&](int i) {
        return val<MAX_TAG_WIDTH>{readt[i]} == htag[i];
      };
      val<NUM_TABLES> tagcmp = tagcmp_split.fo1().concat();
      for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
        match[offset] = concat(val<1>{1}, tagcmp);
      }
    }
    match.fanout(hard<2>{});

    // for each offset, find longest match and select primary prediction
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      match1[offset] = match[offset].one_hot();
    }
    match1.fanout(hard<3 + OVR>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      pred1[offset] = (match1[offset] & preds[offset]) != hard<0>{};
    }
    pred1.fanout(hard<2 + OVR>{});

    // for each offset, find second longest match and select secondary
    // prediction
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      match2[offset] = (match[offset] ^ match1[offset]).one_hot();
    }
    match2.fanout(hard<2>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      pred2[offset] = (match2[offset] & preds[offset]) != hard<0>{};
    }
    pred2.fanout(hard<2 + OVR>{});

    if constexpr (USE_META_V) {
      meta.fanout(hard<2>{});
      arr<val<1>, NUM_TABLES> weakctr = [&](int i) {
        return readh[i] == hard<0>{};
      };
      val<NUM_TABLES> coldctr = notumask & weakctr.fo1().concat();
      coldctr.fanout(hard<FETCH_WIDTH>{});
      val<1> metasign = (meta[METAPIPE_V - 1] >= hard<0>{});
      metasign.fanout(hard<FETCH_WIDTH>{});
      for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
        newly_alloc[offset] = (match1[offset] & coldctr) != hard<0>{};
      }
      newly_alloc.fanout(hard<2 + OVR>{});
      arr<val<1>, FETCH_WIDTH> altsel = [&](u64 offset) {
        arr<val<1>, 3> inputs = {metasign, newly_alloc[offset],
                                 match2[offset] != hard<0>{}};
        return inputs.fo1().fold_and();
      };
      if constexpr (!Overrider::ENABLED) {
        p2 =
            arr<val<1>, FETCH_WIDTH>{[&](u64 offset) {
              return select(altsel[offset].fo1(), pred2[offset], pred1[offset]);
            }}.concat();
      } else {
        // Compute TAGE prediction for offset 0 to pass to SC overrider
        altsel.fanout(hard<2>{});
        val<1> tage_pred_0 = select(altsel[0], pred2[0], pred1[0]);
        auto ovr = overrider.lookup(inst_pc, tage_pred_0, newly_alloc[0]);
        // Build p2 with override baked in at offset 0 — single write
        p2 = arr<val<1>, FETCH_WIDTH>{[&](u64 offset) {
               val<1> tage_pred =
                   select(altsel[offset], pred2[offset], pred1[offset]);
               if (offset == 0)
                 return select(ovr.candidate, ovr.pred, tage_pred);
               return tage_pred;
             }}.concat();
        override_active[0] = ovr.candidate;
        for (u64 offset = 1; offset < FETCH_WIDTH; offset++)
          override_active[offset] = val<1>{0};
#ifdef TAGE_MONITOR
        mon.record_loop_override(0, static_cast<u64>(ovr.candidate) != 0);
#endif
      }
    } else {
      if constexpr (!Overrider::ENABLED) {
        p2 = pred1.concat();
      } else {
        val<1> no_low_conf = val<1>{0};
        auto ovr = overrider.lookup(inst_pc, pred1[0], no_low_conf);
        p2 = arr<val<1>, FETCH_WIDTH>{[&](u64 offset) {
               if (offset == 0)
                 return select(ovr.candidate, ovr.pred, pred1[offset]);
               return pred1[offset];
             }}.concat();
        override_active[0] = ovr.candidate;
        for (u64 offset = 1; offset < FETCH_WIDTH; offset++)
          override_active[offset] = val<1>{0};
#ifdef TAGE_MONITOR
        mon.record_loop_override(0, static_cast<u64>(ovr.candidate) != 0);
#endif
      }
    }

#ifdef TAGE_MONITOR
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      u64 m1 = static_cast<u64>(match1[offset]);
      u64 m2 = static_cast<u64>(match2[offset]);
      u64 m = static_cast<u64>(match[offset]);
      bool p1_bit = (static_cast<u64>(p1) >> offset) & 1;
      bool primary_bit = static_cast<u64>(pred1[offset]) != 0;
      bool final_bit = (static_cast<u64>(p2) >> offset) & 1;
      bool meta_switched = false;
      if constexpr (USE_META_V) {
        meta_switched = (primary_bit != final_bit);
      }
      mon.record_prediction(offset, m1, m2, meta_switched, p1_bit, final_bit);
      mon.record_tag_matches(m);
    }
#endif

    p2.fanout(hard<FETCH_WIDTH>{});
    val<1> taken = (block_entry & p2) != hard<0>{};
    taken.fanout(hard<2>{});
    reuse_prediction(~val<1>{block_entry >> (FETCH_WIDTH - 1)});
    return taken;
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) override {
    val<1> tage_taken = ((block_entry << block_size) & p2) != hard<0>{};
    val<1> taken = [&]() -> val<1> {
      if constexpr (Overrider::ENABLED && Overrider::REUSE_LOOKUP) {
        inst_pc.fanout(hard<2>{});
        val<1> no_low_conf = val<1>{0};
        auto ovr = overrider.lookup(inst_pc, tage_taken, no_low_conf);
        return select(ovr.candidate, ovr.pred, tage_taken);
      } else {
        return tage_taken;
      }
    }();
    taken.fanout(hard<2>{});
    reuse_prediction(~val<1>{block_entry >> (FETCH_WIDTH - 1 - block_size)});
    block_size++;
    return taken;
  }

  void update_condbr(val<64> branch_pc, val<1> taken,
                     [[maybe_unused]] val<64> next_pc) override {
    assert(num_branch < FETCH_WIDTH);
    branch_offset[num_branch] = branch_pc.fo1() >> 2;
    branch_dir[num_branch] = taken.fo1();
    if constexpr (Overrider::ENABLED) {
      overrider.save_branch(branch_pc, num_branch);
    }
    num_branch++;
  }

  void update_cycle(instruction_info &block_end_info) override {
    val<1> &mispredict = block_end_info.is_mispredict;
    val<64> &next_pc = block_end_info.next_pc;

    if (num_branch == 0) {
      // no conditional branch in this block
      val<1> line_end = block_entry >> (FETCH_WIDTH - block_size);
      val<1> actual_block = ~(true_block & line_end.fo1());
      actual_block.fanout(hard<MAXHIST + NUM_TABLES * 2 + 2>{});
      execute_if(actual_block, [&]() {
        next_pc.fanout(hard<2>{});
        if constexpr (P1_USE_GSHARE_V) {
          global_history1 =
              (global_history1 << 1) ^ val<P1_HIST_V>{next_pc >> 2};
        }
        gfolds.update(val<PATHBITS>{next_pc >> 2});
        if constexpr (USE_PATH_HIST_V) {
          path_hist =
              (path_hist << PATHBITS) ^ val<PATH_HIST_WIDTH_V>{next_pc >> 2};
        }
        true_block = 1;
      });
      return;
    }

    mispredict.fanout(hard<NUM_TABLES + 2>{});
    val<1> correct_pred = ~mispredict;
    correct_pred.fanout(hard<NUM_TABLES + 2>{});
    index1.fanout(hard<FETCH_WIDTH * 3>{});
    p2.fanout(hard<2>{});
    bindex.fanout(hard<FETCH_WIDTH * 3>{});
    gindex.fanout(hard<4>{});
    htag.fanout(hard<3>{});
    readb.fanout(hard<2>{});
    readt.fanout(hard<4>{});
    readc.fanout(hard<2>{});
    if constexpr (MAX_HYST_WIDTH > 0) {
      readh.fanout(hard<3>{});
    }
    match1.fanout(hard<3>{});
    match2.fanout(hard<2>{});
    pred1.fanout(hard<2>{});
    pred2.fanout(hard<2 + NUM_TABLES>{});
    branch_offset.fanout(hard<FETCH_WIDTH + NUM_TABLES + 1>{});
    branch_dir.fanout(hard<2>{});
    gfolds.fanout(hard<2>{});
    readu.fanout(hard<2>{});
    if constexpr (USE_META_V) {
      meta.fanout(hard<2>{});
    }
    val<LOG_FETCH_WIDTH> last_offset = branch_offset[num_branch - 1];
    last_offset.fanout(hard<4 * NUM_TABLES + 2>{});
    u64 update_valid = (u64(1) << num_branch) - 1;
    arr<val<FETCH_WIDTH>, FETCH_WIDTH> update_mask = [&](u64 offset) {
      arr<val<1>, FETCH_WIDTH> match_offset = [&](u64 i) {
        return branch_offset[i] == offset;
      };
      return match_offset.fo1().concat() & update_valid;
    };
    update_mask.fanout(hard<2>{});

    arr<val<1>, FETCH_WIDTH> is_branch = [&](u64 offset) {
      return update_mask[offset] != hard<0>{};
    };
    is_branch.fanout(hard<6>{});

    val<FETCH_WIDTH> branch_mask = is_branch.concat();

    val<FETCH_WIDTH> actualdirs = branch_dir.concat();
    actualdirs.fanout(hard<FETCH_WIDTH>{});

    arr<val<1>, FETCH_WIDTH> branch_taken = [&](u64 offset) {
      return (actualdirs & update_mask[offset]) != hard<0>{};
    };
    branch_taken.fanout(hard<3>{});

#ifdef TAGE_MONITOR
    {
      bool is_misp = static_cast<u64>(mispredict) != 0;
      for (u64 i = 0; i < num_branch; i++) {
        bool actual = static_cast<u64>(branch_dir[i]) != 0;
        mon.record_outcome(i, actual, is_misp && (i == num_branch - 1));
      }
      mon.record_block_advance(num_branch);
    }
#endif

    arr<val<NUM_TABLES + 1>, FETCH_WIDTH> actual_match1 = [&](u64 offset) {
      return select(is_branch[offset], match1[offset], val<NUM_TABLES + 1>{0});
    };
    actual_match1.fanout(hard<2>{});

    val<NUM_TABLES> primary_mask = actual_match1.fold_or();
    primary_mask.fanout(hard<2>{});
    arr<val<1>, NUM_TABLES> primary = primary_mask.make_array(val<1>{});
    primary.fanout(hard<3>{});

    arr<val<1>, FETCH_WIDTH> primary_wrong = [&](u64 offset) {
      return pred1[offset] != branch_taken[offset];
    };
    primary_wrong.fanout(hard<2>{});

    // select candidate entries for allocation
    val<NUM_TABLES> mispmask =
        mispredict.replicate(hard<NUM_TABLES>{}).concat();
    arr<val<1>, NUM_TABLES> last_tagcmp = [&](int i) {
      if constexpr (BR_P_ENTRY_V == 1) {
        return readt[i] == concat(last_offset, htag[i]);
      } else {
        return val<MAX_TAG_WIDTH>{readt[i]} == htag[i];
      }
    };
    val<NUM_TABLES + 1> last_match1 =
        last_tagcmp.fo1().append(1).concat().one_hot();
    last_match1.fanout(hard<2>{});
    val<NUM_TABLES> postmask =
        mispmask.fo1() & val<NUM_TABLES>(last_match1 - 1);
    postmask.fanout(hard<2>{});
    val<NUM_TABLES> candallocmask = postmask & notumask;
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
    allocate.fanout(hard<7>{});

#ifdef TAGE_MONITOR
    if (static_cast<u64>(mispredict)) {
      u64 amask = 0;
      for (u64 i = 0; i < NUM_TABLES; i++)
        amask |= (static_cast<u64>(allocate[i]) & 1) << i;
      mon.record_allocation(amask != 0, amask);
    }
#endif

    // associate a branch direction to each global table
    arr<val<1>, NUM_TABLES> bdir = [&](u64 i) {
      if constexpr (BR_P_ENTRY_V == 1) {
        val<LOG_FETCH_WIDTH> tag_offset = readt[i] >> MAX_HTAGBITS;
        val<LOG_FETCH_WIDTH> offset =
            select(allocate[i], last_offset, tag_offset.fo1());
        offset.fanout(hard<FETCH_WIDTH>{});
        arr<val<1>, FETCH_WIDTH> match_offset = [&](u64 j) {
          return branch_offset[j] == offset;
        };
        return (match_offset.fo1().concat() & update_valid & actualdirs) !=
               hard<0>{};
      } else {
        return branch_taken[0]; // simplified for multi-slot
      }
    };
    bdir.fanout(hard<2>{});

    // tell if global prediction is incorrect
    arr<val<1>, NUM_TABLES> badpred1 = [&](u64 i) -> val<1> {
      if constexpr (MAX_CTR_WIDTH == 1) {
        return readc[i] != bdir[i];
      } else {
        return val<1>{readc[i] >> hard<MAX_CTR_WIDTH - 1>{}} != bdir[i];
      }
    };
    badpred1.fanout(hard<3>{});

    // does primary differ from alt?
    arr<val<1>, NUM_TABLES> altdiffer = [&](u64 i) -> val<1> {
      auto pred_dir = [&]() -> val<1> {
        if constexpr (MAX_CTR_WIDTH == 1) {
          return readc[i];
        } else {
          return readc[i] >> hard<MAX_CTR_WIDTH - 1>{};
        }
      }();
      if constexpr (BR_P_ENTRY_V == 1) {
        val<LOG_FETCH_WIDTH> tag_offset = readt[i] >> MAX_HTAGBITS;
        return pred_dir != pred2.select(tag_offset.fo1());
      } else {
        return pred_dir != pred2[0];
      }
    };

    // is owning branch's prediction correct?
    arr<val<1>, NUM_TABLES> goodpred = [&](u64 i) {
      if constexpr (BR_P_ENTRY_V == 1) {
        val<LOG_FETCH_WIDTH> tag_offset = readt[i] >> MAX_HTAGBITS;
        return (tag_offset.fo1() != last_offset) | correct_pred;
      } else {
        return correct_pred;
      }
    };
    goodpred.fanout(hard<2>{});

    // do P1 and P2 agree?
    val<FETCH_WIDTH> disagree_mask = (p1 ^ p2) & branch_mask.fo1();
    disagree_mask.fanout(hard<2>{});
    arr<val<1>, FETCH_WIDTH> disagree = disagree_mask.make_array(val<1>{});
    disagree.fanout(hard<2>{});

    // read the P1 hysteresis if P1 and P2 disagree
    arr<val<1>, FETCH_WIDTH> p1_weak = [&](u64 offset) -> val<1> {
      return execute_if(disagree[offset],
                        [&]() { return ~table1_hyst[offset].read(index1); });
    };

    // read the bimodal hysteresis if bimodal caused a misprediction
    arr<val<1>, FETCH_WIDTH> b_weak = [&](u64 offset) -> val<1> {
      val<1> bim_primary = actual_match1[offset] >> NUM_TABLES;
      return execute_if(bim_primary.fo1() & primary_wrong[offset],
                        [&]() { return ~bhyst[offset].read(bindex); });
    };

    // determine which primary global predictions are incorrect with weak hyst
    arr<val<1>, NUM_TABLES> g_weak = [&](u64 i) -> val<1> {
      if constexpr (MAX_HYST_WIDTH > 0) {
        return primary[i] & badpred1[i] & (readh[i] == hard<0>{});
      } else {
        return primary[i] & badpred1[i];
      }
    };
    g_weak.fanout(hard<2>{});

    // SC: compute new weights BEFORE need_extra_cycle (no RAM writes)
    if constexpr (Overrider::ENABLED) {
      overrider.train_compute(branch_taken, is_branch, mispredict, correct_pred,
                              num_branch);
    }

    // need extra cycle for modifying prediction bits, TAGE allocation, and SC writes
    val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
    val<1> extra_cycle = [&]() -> val<1> {
      val<1> base = some_badpred1.fo1() | mispredict | (disagree_mask != hard<0>{});
      if constexpr (Overrider::ENABLED) {
        // SC always needs extra cycle when it has staged writes
        return base | val<1>{overrider.should_write ? 1u : 0u};
      } else {
        return base;
      }
    }();
    static constexpr u64 OVR_BUDGET = Overrider::ENABLED ? Overrider::EXTRA_CYCLE_BUDGET : 0;
    extra_cycle.fanout(hard<NUM_TABLES * 2 + 1 + OVR_BUDGET>{});
    need_extra_cycle(extra_cycle);

    // SC: write staged data to RAM. extra_cycle is guaranteed=1 when
    // should_write=true (we ORed it into the extra_cycle computation).
    if constexpr (Overrider::ENABLED) {
      if (overrider.should_write) {
        overrider.train_write();
      }
    }

    // update meta counter
    if constexpr (USE_META_V) {
      arr<val<1>, FETCH_WIDTH> altdiff = [&](u64 offset) {
        return (match2[offset] != hard<0>{}) & (pred2[offset] != pred1[offset]);
      };
      arr<val<2, i64>, FETCH_WIDTH> meta_incr = [&](u64 offset) -> val<2, i64> {
        val<1> update_meta =
            is_branch[offset] & altdiff[offset].fo1() & newly_alloc[offset];
        val<1> bad_pred2 = (pred2[offset] != branch_taken[offset]);
        return select(update_meta.fo1(), concat(bad_pred2.fo1(), val<1>{1}),
                      val<2>{0});
      };
      for (u64 i = METAPIPE_V - 1; i != 0; i--) {
        meta[i] = meta[i - 1];
      }
      auto newmeta = meta[0] + meta_incr.fo1().fold_add();
      newmeta.fanout(hard<3>{});
      using meta_t = valt<decltype(meta[0])>;
      meta[0] = select(newmeta > meta_t::maxval, meta_t{meta_t::maxval},
                       select(newmeta < meta_t::minval, meta_t{meta_t::minval},
                              meta_t{newmeta}));
    }

    // overwrite the tag in the allocated entry (mispredict)
    static_loop<NUM_TABLES>([&]<u64 I>() {
      execute_if(allocate[I], [&]() {
        auto &t = std::get<I>(tables);
        if constexpr (BR_P_ENTRY_V == 1) {
          t.tag_ram[0].write(tidx<I>(gindex[I]), concat(last_offset, htag[I]));
        } else {
          t.tag_ram[0].write(tidx<I>(gindex[I]), val<MAX_TAG_WIDTH>{htag[I]});
        }
      });
#ifdef TAGE_MONITOR
      if (static_cast<u64>(allocate[I])) {
        u64 idx = static_cast<u64>(tidx<I>(gindex[I]));
        mon.record_tag_write(I, idx);
        mon.record_tage_alias(I, idx, idx ^ (I << 16));
      }
#endif
    });

    // update the u bits
    arr<val<1>, NUM_TABLES> update_u = [&](u64 i) {
      return primary[i] & altdiffer[i].fo1();
    };
    // if all post entries have the u bit set, reset their u bits
    val<1> noalloc = (candallocmask == hard<0>{});
    val<NUM_TABLES> uclearmask =
        postmask & noalloc.fo1().replicate(hard<NUM_TABLES>{}).concat();
    arr<val<1>, NUM_TABLES> uclear = uclearmask.fo1().make_array(val<1>{});
    uclear.fanout(hard<2>{});
    if constexpr (USE_PROB_DECAY) {
      // Probabilistic decay: every touched entry has a chance of being forced
      // to u=0. decay_fire is independent of uclear — provides background aging
      // pressure.
      val<DECAY_CTR_V> lfsr = val<DECAY_CTR_V>{static_cast<u64>(std::rand())};
      val<1> decay_fire = (lfsr > val<DECAY_CTR_V>{decay_threshold});
#ifdef TAGE_MONITOR
      if (static_cast<u64>(decay_fire)) {
        mon.cum.decay_fire_count++;
        mon.win.decay_fire_count++;
      }
#endif
      decay_fire.fanout(hard<NUM_TABLES>{});
      static_loop<NUM_TABLES>([&]<u64 I>() {
        execute_if(
            update_u[I].fo1() | allocate[I] | uclear[I] | decay_fire, [&]() {
              auto &t = std::get<I>(tables);
              val<1> newu =
                  goodpred[I].fo1() & ~allocate[I] & ~uclear[I] & ~decay_fire;
#ifdef TAGE_MONITOR
              mon.record_u_write(I, static_cast<u64>(tidx<I>(gindex[I])),
                                 static_cast<u64>(newu) != 0);
#endif
              if constexpr (U_STOR_FF_V) {
                t.write_u_ff_arr(0, tidx<I>(gindex[I]), newu);
              } else {
                t.u_ram[0].write(tidx<I>(gindex[I]), newu.fo1(), extra_cycle);
              }
            });
      });
    } else {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        execute_if(update_u[I].fo1() | allocate[I] | uclear[I], [&]() {
          auto &t = std::get<I>(tables);
          val<1> newu = goodpred[I].fo1() & ~allocate[I] & ~uclear[I];
#ifdef TAGE_MONITOR
          mon.record_u_write(I, static_cast<u64>(tidx<I>(gindex[I])),
                             static_cast<u64>(newu) != 0);
#endif
          if constexpr (U_STOR_FF_V) {
            t.write_u_ff_arr(0, tidx<I>(gindex[I]), newu);
          } else {
            t.u_ram[0].write(tidx<I>(gindex[I]), newu.fo1(), extra_cycle);
          }
        });
      });
    }

    // update P1 prediction if P1 and P2 disagree and the hysteresis bit is weak
    auto p2_split = p2.make_array(val<1>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(p1_weak[offset].fo1(), [&]() {
        table1_pred[offset].write(index1, p2_split[offset].fo1());
      });
#ifdef TAGE_MONITOR
      if (static_cast<u64>(p1_weak[offset]))
        mon.record_p1_write(offset, static_cast<u64>(index1),
                            static_cast<u64>(index1) ^ (offset << 16));
#endif
    }
    // update P1 hysteresis
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(is_branch[offset], [&]() {
        table1_hyst[offset].write(index1, ~disagree[offset]);
      });
    }

    // update incorrect bimodal prediction if primary provider and hyst is weak
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(b_weak[offset].fo1(),
                 [&]() { bim[offset].write(bindex, branch_taken[offset]); });
    }
    // update bimodal hysteresis if bimodal is primary provider
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      val<1> bim_primary = match1[offset] >> NUM_TABLES;
      execute_if(is_branch[offset] & bim_primary.fo1(), [&]() {
        bhyst[offset].write(bindex, ~primary_wrong[offset]);
      });
    }

    // update incorrect global prediction if primary provider and hyst is weak;
    // initialize global prediction in the allocated entry
    static_loop<NUM_TABLES>([&]<u64 I>() {
      execute_if(g_weak[I].fo1() | allocate[I], [&]() {
        std::get<I>(tables).pred_ram[0].write(tidx<I>(gindex[I]), bdir[I]);
      });
    });
    // update global prediction hysteresis if primary provider or allocated
    // entry. Skip silent writes (counter already saturated in update
    // direction).
    if constexpr (MAX_HYST_WIDTH > 0) {
      static constexpr u64 HW = std::max(u64(1), MAX_HYST_WIDTH);
      static constexpr u64 HMAX = (u64(1) << HW) - 1;
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> would_change = allocate[I] |
                              (badpred1[I] & (readh[I] != hard<0>{})) |
                              (~badpred1[I] & (readh[I] != hard<HMAX>{}));
        execute_if((primary[I] | allocate[I]) & would_change, [&]() {
          auto newhyst = select(allocate[I], val<HW>{0},
                                update_ctr(readh[I], ~badpred1[I]));
          std::get<I>(tables).hyst_ram[0].write(tidx<I>(gindex[I]),
                                                newhyst.fo1(), extra_cycle);
        });
      });
    }

    // u-bit epoch reset / threshold adaptation
    uctr.fanout(hard<3>{});
    val<NUM_TABLES> allocmask1 = collamask1.reverse();
    allocmask1.fanout(hard<2>{});
    val<1> faralloc =
        (((last_match1 >> 3) | allocmask1).one_hot() ^ allocmask1) == hard<0>{};
    val<1> uctrsat = (uctr == hard<decltype(uctr)::maxval>{});
    uctrsat.fanout(hard<2>{});
    uctr = select(correct_pred, uctr,
                  select(uctrsat, val<decltype(uctr)::size>{0},
                         update_ctr(uctr, faralloc.fo1())));
#ifdef TAGE_MONITOR
    if (static_cast<u64>(uctrsat)) {
      mon.cum.uctr_saturation_count++;
      mon.win.uctr_saturation_count++;
      if constexpr (!USE_PROB_DECAY && !U_STOR_FF_V) {
        mon.cum.epoch_reset_count++;
        mon.win.epoch_reset_count++;
      }
      mon.record_epoch_reset();
    }
#endif
    if constexpr (USE_PROB_DECAY) {
      // Adaptive threshold: decrease on allocation pressure, increase on
      // correct prediction. DECAY_GRAN_V controls update rate: 0=every uctr
      // change, 6=every 64 increments.
      val<1> threshold_tick = [&]() -> val<1> {
        if constexpr (DECAY_GRAN_V == 0) {
          return ~correct_pred; // update every misprediction cycle
        } else {
          return (uctr & hard<(u64(1) << DECAY_GRAN_V) - 1>{}) == hard<0>{};
        }
      }();
      val<1> misp = ~correct_pred;
      decay_threshold =
          select(threshold_tick,
                 DecayPolicy_V::template apply<DECAY_CTR_V>(
                     decay_threshold, correct_pred, uctrsat, misp),
                 val<DECAY_CTR_V>{decay_threshold});
    } else if constexpr (!U_STOR_FF_V) {
      // Epoch reset: bulk-clear all SRAM u-bits when uctr saturates.
      // Not available for FF u-bits (double write conflict with
      // write_u_ff_arr).
      execute_if(uctrsat, [&]() {
        static_loop<NUM_TABLES>(
            [&]<u64 I>() { std::get<I>(tables).u_ram[0].reset(); });
      });
    }
    // Note: U_STOR_FF with DECAY_CTR=0 has no epoch reset — u-bits only
    // clear through uclear (allocation pressure). Consider using DECAY_CTR>0.

    // update global history
    val<1> line_end = block_entry >> (FETCH_WIDTH - block_size);
    true_block = correct_pred | branch_dir[num_branch - 1] | line_end.fo1();
    true_block.fanout(hard<MAXHIST + NUM_TABLES * 2 + 2>{});
    execute_if(true_block, [&]() {
      next_pc.fanout(hard<2>{});
      if constexpr (P1_USE_GSHARE_V) {
        global_history1 = (global_history1 << 1) ^ val<P1_HIST_V>{next_pc >> 2};
      }
      gfolds.update(val<PATHBITS>{next_pc >> 2});
      if constexpr (USE_PATH_HIST_V) {
        path_hist =
            (path_hist << PATHBITS) ^ val<PATH_HIST_WIDTH_V>{next_pc >> 2};
      }
    });

    num_branch = 0; // done
  }
};
#endif // TAGE_HPP
