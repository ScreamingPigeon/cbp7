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
  if (exp == 0.0)
    return 1.0;
  if (base == 0.0)
    return 0.0;
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
    u64 hl =
        static_cast<u64>(static_cast<double>(min_h) * constexpr_pow(ratio, e));
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
  static constexpr auto HIST_LEN = geometric_hist<NUM_TABLES>(MINHIST, MAXHIST);
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
  using type = std::tuple<
      TableAt<Cfg, BR_P_ENTRY, NUM_BANKS, USE_AHEAD, SHARED_TAG, SHARED_U,
              U_STOR_FF, DECAY_CTR, ResetFn, USE_FF_CACHE, Is>...>;
};

// ============================================================================
// Forward Declaration
// ============================================================================

template <bool USE_AHEAD_V, typename TableCfg, typename AllocCfg,
          // Global hardware params
          u64 FETCH_WIDTH_V, u64 BIMODAL_SIZE_V, u64 BR_P_ENTRY_V,
          u64 NUM_BANKS_V, bool SHARED_TAG_V, bool SHARED_U_V, bool U_STOR_FF_V,
          u64 DECAY_CTR_V, typename ResetFn_V, bool USE_FF_CACHE_V,
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
      typename MakeTableTuple<TableCfg, BR_P_ENTRY_V, NUM_BANKS_V, USE_AHEAD_V,
                              SHARED_TAG_V, SHARED_U_V, U_STOR_FF_V,
                              DECAY_CTR_V, ResetFn_V, USE_FF_CACHE_V,
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
// This is the predictor

template <typename TableCfg, typename AllocCfg, u64 FETCH_WIDTH_V,
          u64 BIMODAL_SIZE_V, u64 BR_P_ENTRY_V, u64 NUM_BANKS_V,
          bool SHARED_TAG_V, bool SHARED_U_V, bool U_STOR_FF_V, u64 DECAY_CTR_V,
          typename ResetFn_V, bool USE_FF_CACHE_V, bool P1_USE_GSHARE_V,
          u64 P1_TABLE_SIZE_V, u64 P1_HIST_V, bool USE_META_V, u64 METABITS_V,
          u64 METAPIPE_V, u64 META_TABLE_SIZE_V, bool USE_PATH_HIST_V,
          u64 PATH_HIST_WIDTH_V, u64 PATH_BITS_V>
struct TageImpl<false, TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V,
                BR_P_ENTRY_V, NUM_BANKS_V, SHARED_TAG_V, SHARED_U_V,
                U_STOR_FF_V, DECAY_CTR_V, ResetFn_V, USE_FF_CACHE_V,
                P1_USE_GSHARE_V, P1_TABLE_SIZE_V, P1_HIST_V, USE_META_V,
                METABITS_V, METAPIPE_V, META_TABLE_SIZE_V, USE_PATH_HIST_V,
                PATH_HIST_WIDTH_V, PATH_BITS_V> : predictor {

  using Base =
      TageBase<TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V, BR_P_ENTRY_V,
               NUM_BANKS_V, false, SHARED_TAG_V, SHARED_U_V, U_STOR_FF_V,
               DECAY_CTR_V, ResetFn_V, USE_FF_CACHE_V, P1_USE_GSHARE_V,
               P1_TABLE_SIZE_V, P1_HIST_V, USE_META_V, METABITS_V, METAPIPE_V,
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
  // Meta table size: at least 1 to keep ram<> valid even when conditional is
  // false
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
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> gindex; // per-table indices
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> htag;   // per-table hashed tags
  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> readt; // read tags
  arr<reg<1>, NUM_TABLES> readc;             // read prediction bits
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> readh; // read hysteresis
  arr<reg<1>, NUM_TABLES> readu;             // read u-bits
  reg<NUM_TABLES> notumask;                  // inverted u-bits

  // ======== Match Registers ========
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match;  // all matches per offset
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match1; // longest match per offset
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match2; // second longest per offset
  arr<reg<1>, FETCH_WIDTH> pred1;               // primary prediction per offset
  arr<reg<1>, FETCH_WIDTH> pred2; // alternate prediction per offset
  reg<FETCH_WIDTH> p2;            // final P2 predictions

  // ======== Meta-prediction ========
  // USE_META=true: signed counter for alt-on-newly-allocated
  // META_TABLE_SIZE=0: global counter
  // META_TABLE_SIZE>0: PC-indexed table
  std::conditional_t<USE_META_V && META_TABLE_SIZE_V == 0,
                     arr<reg<METABITS_V, i64>, METAPIPE_V>, EmptyMember>
      meta_global;
  std::conditional_t<USE_META_V && (META_TABLE_SIZE_V > 0),
                     hcm::ram<val<METABITS_V, i64>, META_RAM_SIZE>, EmptyMember>
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

  // ======== Derived Constants ========
  static constexpr u64 LINEADDR_BITS = std::max(BINDEX_BITS, MAX_IDX_BITS);
  static constexpr u64 PATHBITS = Base::PATH_BITS;

  // ======== Phase 2: Hash Functions (§2.1) ========

  // Compute bimodal index and per-table TAGE indexes.
  void compute_indexes(val<LINEADDR_BITS> lineaddr) {
    bindex = lineaddr;
    bindex.fanout(hard<FETCH_WIDTH>{});
    for (u64 i = 0; i < NUM_TABLES; i++) {
      gindex[i] = lineaddr ^ gfolds.template get<0>(i);
    }
    if constexpr (USE_PATH_HIST_V) {
      for (u64 i = 0; i < NUM_TABLES; i++) {
        gindex[i] = val<MAX_IDX_BITS>{gindex[i]} ^ val<MAX_IDX_BITS>{path_hist};
      }
    }
    gindex.fanout(hard<2>{});
  }

  // Compute per-table hashed tags (the hash portion, without offset).
  void compute_tags(val<LINEADDR_BITS> lineaddr) {
    for (u64 i = 0; i < NUM_TABLES; i++) {
      htag[i] =
          val<MAX_HTAGBITS>{lineaddr}.reverse() ^ gfolds.template get<1>(i);
    }
    htag.fanout(hard<2>{});
  }

  // ======== Phase 2: Table Read (§2.3) ========

  // Read all TageTable instances and bimodal tables, populating result
  // registers (readt, readc, readh, readu, notumask, readb).
  // Caller must set gindex, htag, bindex before calling (via §2.1).
  void read_tables() {
    // Fanout for multi-use registers

    // Read bimodal tables (one per offset in the fetch block)
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      readb[offset] = bim[offset].read(bindex);
    }
    readb.fanout(hard<2>{});

    // Read each TAGE table via static_loop (tables are heterogeneous)
    static_loop<NUM_TABLES>([&]<u64 I>() {
      constexpr u64 TW = TableCfg::TAG_WIDTH[I];
      constexpr u64 CW = TableCfg::CTR_WIDTH[I];

      // For BR_P_ENTRY=1, compare_tag is just htag zero-extended to TAG_WIDTH.
      // TageTable's internal hit_reg won't be fully correct (missing offset),
      // but it only affects SRAM u-bit decay (1/DECAY_CTR probability).
      // For BR_P_ENTRY>1, htag IS the full tag.
      val<TW> compare_tag = htag[I];
      std::get<I>(tables).read(gindex[I], compare_tag);

      // Extract results into uniform-width registers
      readt[I] = std::get<I>(tables).getTag();
      val<CW> ctr = std::get<I>(tables).getCounter(0);
      readc[I] = ctr >> hard<CW - 1>{}; // MSB = prediction
      readh[I] =
          ctr & hard<(1ULL << (CW - 1)) - 1>{}; // lower bits = hysteresis only
      readu[I] = std::get<I>(tables).getU();
    });
    readt.fanout(hard<FETCH_WIDTH + 1>{});
    readc.fanout(hard<3>{});
    readh.fanout(hard<2>{});
    readu.fanout(hard<2>{});

    // Fanout for downstream use
    notumask = ~readu.concat();
    notumask.fanout(hard<2>{});
  }

  // ======== Phase 2: Match Logic (§2.4) ========

  // Build match masks, find provider (longest match) and altpred (second
  // longest) for each offset. Populates match, match1, match2, pred1, pred2.
  // Caller must have called read_tables() first.
  void compute_matches() {
    // Gather prediction bits: per-table readc + bimodal readb per offset.
    // match vector is (NUM_TABLES+1) wide: bit 0..NUM_TABLES-1 = tables,
    // bit NUM_TABLES = bimodal (always matches as fallback).
    val<NUM_TABLES> gpreds = readc.concat();
    gpreds.fanout(hard<FETCH_WIDTH>{});
    arr<val<NUM_TABLES + 1>, FETCH_WIDTH> preds = [&](u64 offset) {
      return concat(readb[offset], gpreds);
    };
    preds.fanout(hard<2 * FETCH_WIDTH>{});

    if constexpr (BR_P_ENTRY == 1) {
      // BR_P_ENTRY=1: offset encoded in upper tag bits.
      // Tag stored = concat(offset, htag). readt has full stored tag.
      // Match requires both htag match AND offset match.
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
      // BR_P_ENTRY>1: full tag match, same for all offsets.
      arr<val<1>, NUM_TABLES> tagcmp_split = [&](int i) {
        return val<MAX_TAG_WIDTH>{readt[i]} == htag[i];
      };
      val<NUM_TABLES> tagcmp = tagcmp_split.fo1().concat();
      for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
        match[offset] = concat(val<1>{1}, tagcmp);
      }
    }

    // Longest match → provider prediction
    match.fanout(hard<2>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      match1[offset] = match[offset].one_hot();
    }
    match1.fanout(hard<3>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      pred1[offset] = (match1[offset] & preds[offset]) != hard<0>{};
    }
    pred1.fanout(hard<2>{});

    // Second longest match → alternate prediction
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      match2[offset] = (match[offset] ^ match1[offset]).one_hot();
    }
    match2.fanout(hard<2>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      pred2[offset] = (match2[offset] & preds[offset]) != hard<0>{};
    }
    pred2.fanout(hard<2>{});
  }

  // ======== Phase 2: Meta-prediction (§2.5) ========

  // Compute final p2 prediction per offset using alt-on-newly-allocated logic.
  // When USE_META: if the provider looks recently allocated (weak counter,
  // u=0) and the meta counter says alt is better, use pred2 instead of pred1.
  // Caller must have called compute_matches() first.
  void compute_meta_prediction() {
    if constexpr (!USE_META_V) {
      p2 = pred1.concat();
    } else {
      // Detect newly-allocated providers: weak hysteresis AND u=0
      if constexpr (META_TABLE_SIZE_V == 0) {
        meta_global.fanout(hard<2>{});
      }
      arr<val<1>, NUM_TABLES> weakctr = [&](int i) {
        return readh[i] == hard<0>{};
      };
      val<NUM_TABLES> coldctr = notumask & weakctr.fo1().concat();
      coldctr.fanout(hard<FETCH_WIDTH>{});

      // Read meta counter sign
      auto metasign = [&]() -> val<1> {
        if constexpr (META_TABLE_SIZE_V == 0) {
          // Global counter mode
          return (meta_global[Base::METAPIPE - 1] >= hard<0>{});
        } else {
          // PC-indexed table mode — TODO: read meta_table in predict2
          return val<1>{0};
        }
      }();
      metasign.fanout(hard<FETCH_WIDTH>{});

      for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
        newly_alloc[offset] = (match1[offset] & coldctr) != hard<0>{};
      }
      newly_alloc.fanout(hard<2>{});

      // altsel = metasign AND newly_alloc AND (altpred exists)
      arr<val<1>, FETCH_WIDTH> altsel = [&](u64 offset) {
        arr<val<1>, 3> inputs = {metasign, newly_alloc[offset],
                                 match2[offset] != hard<0>{}};
        return inputs.fo1().fold_and();
      };

      p2 = arr<val<1>, FETCH_WIDTH>{[&](u64 offset) {
             return select(altsel[offset], pred2[offset], pred1[offset]);
           }}.concat();
    }
  }

  // ======== Phase 2: Allocation (§2.8) ========

  // Results from allocation routine, consumed by update_cycle.
  struct AllocResult {
    arr<val<1>, NUM_TABLES> allocate; // per-table: allocate this entry?
    val<NUM_TABLES> postmask;        // tables with longer history than provider
    val<NUM_TABLES> candallocmask;   // postmask & u=0
    val<NUM_TABLES> collamask1;      // first candidate (reversed one_hot)
    val<NUM_TABLES + 1> last_match1; // provider match for last branch
  };

  // Select candidate entries for allocation on misprediction.
  // Finds tables with longer history than the mispredicting branch's provider
  // that have u=0, then selects one with randomization.
  AllocResult allocate_entries(val<1> mispredict,
                               val<LOG_FETCH_WIDTH> last_offset) {

    // Find provider for the last (mispredicting) branch by full tag comparison
    arr<val<1>, NUM_TABLES> last_tagcmp = [&](int i) {
      if constexpr (BR_P_ENTRY == 1) {
        return readt[i] == concat(last_offset, htag[i]);
      } else {
        return val<MAX_TAG_WIDTH>{readt[i]} == htag[i];
      }
    };
    val<NUM_TABLES + 1> last_match1 =
        last_tagcmp.fo1().append(1).concat().one_hot();

    // Post-provider mask: tables with longer history than provider
    // (one_hot - 1 gives all bits below the provider)
    val<NUM_TABLES> mispmask =
        mispredict.replicate(hard<NUM_TABLES>{}).concat();
    val<NUM_TABLES> postmask =
        mispmask.fo1() & val<NUM_TABLES>(last_match1 - 1);

    // Candidate = post-provider AND u=0
    val<NUM_TABLES> candallocmask = postmask & notumask;

    // Select one candidate with randomization.
    // Reverse so one_hot picks shortest-history candidate first,
    // then reverse back for table indexing.
    val<NUM_TABLES> collamask = candallocmask.reverse();
    collamask.fanout(hard<2>{});
    val<NUM_TABLES> collamask1 = collamask.one_hot();
    collamask1.fanout(hard<3>{});
    val<NUM_TABLES> collamask2 = (collamask ^ collamask1).one_hot();
    // 25% chance of picking the second candidate instead
    val<NUM_TABLES> collamask12 =
        select(val<2>{std::rand()} == hard<0>{}, collamask2.fo1(), collamask1);
    arr<val<1>, NUM_TABLES> allocate =
        collamask12.fo1().reverse().make_array(val<1>{});

    return {allocate, postmask, candallocmask, collamask1, last_match1};
  }

  // ======== Phase 2: Counter/U-bit Update (§2.9) ========

  // Update TAGE table entries: tag writes for allocation, counter updates
  // for primary providers, u-bit updates. One TageTable::write() per table.
  //
  // Per table, the write triggers when: allocate OR primary OR update_u OR
  // uclear. The new entry packs:
  //   tag:  allocated → full tag (concat offset+htag or htag), else keep readt
  //   ctr:  allocated → (bdir << (CW-1)) | 0  (pred=bdir, hyst=0)
  //         g_weak   → (bdir << (CW-1)) | 0  (flip pred, hyst stays 0)
  //         primary  → (readc << (CW-1)) | update_ctr(hyst, ~badpred)
  //         else     → keep existing
  //   u:    goodpred & ~allocate & ~uclear
  void update_tables(
      arr<val<1>, NUM_TABLES> &allocate, arr<val<1>, NUM_TABLES> &bdir,
      arr<val<1>, NUM_TABLES> &badpred1, arr<val<1>, NUM_TABLES> &goodpred,
      arr<val<1>, NUM_TABLES> &primary, arr<val<1>, NUM_TABLES> &altdiffer,
      arr<val<1>, NUM_TABLES> &uclear, arr<val<1>, NUM_TABLES> &g_weak,
      val<LOG_FETCH_WIDTH> last_offset,
      val<1> extra_cycle) {
    // u-bit update condition: primary provider where alt prediction differs
    arr<val<1>, NUM_TABLES> update_u = [&](u64 i) {
      return primary[i] & altdiffer[i].fo1();
    };

    static_loop<NUM_TABLES>([&]<u64 I>() {
      constexpr u64 CW = TableCfg::CTR_WIDTH[I];
      constexpr u64 TW = TableCfg::TAG_WIDTH[I];
      constexpr u64 UW = TableCfg::U_WIDTH[I];
      constexpr u64 HYST_W = CW - 1;

      // Gate all writes by extra_cycle to prevent same-cycle read+write
      // on TageTable bank_ram (read happens in predict2, write here).
      val<1> should_write =
          extra_cycle & (allocate[I] | primary[I] | update_u[I].fo1() | uclear[I]);

      execute_if(should_write, [&]() {
        // --- Tag ---
        val<TW> new_tag = [&]() -> val<TW> {
          if constexpr (BR_P_ENTRY == 1) {
            return select(allocate[I], val<TW>{concat(last_offset, htag[I])},
                          val<TW>{readt[I]});
          } else {
            return select(allocate[I], val<TW>{htag[I]}, val<TW>{readt[I]});
          }
        }();

        // --- Counter (pred + hyst) ---
        val<HYST_W> cur_hyst = readh[I]; // truncated to per-table width
        val<HYST_W> updated_hyst = update_ctr(cur_hyst, ~badpred1[I]);
        val<1> new_pred = select(g_weak[I] | allocate[I], bdir[I], readc[I]);
        val<HYST_W> new_hyst =
            select(allocate[I], val<HYST_W>{0}, updated_hyst);
        // Pack: pred is MSB, hyst is lower bits
        val<CW> new_ctr = concat(new_pred, new_hyst);

        // --- U-bit ---
        val<UW> new_u = val<UW>{goodpred[I].fo1() & ~allocate[I] & ~uclear[I]};

        std::get<I>(tables).write(gindex[I], 0, 0, new_tag, new_ctr, new_u);
      });
    });
  }

  // ======== Phase 2: Bimodal Update (§2.10) ========

  // Part 1: Read bimodal hysteresis (BEFORE need_extra_cycle, same cycle as predict2).
  // Returns b_weak: 1 iff bimodal is primary, prediction wrong, and hysteresis weak.
  arr<val<1>, FETCH_WIDTH> compute_b_weak(
      arr<val<NUM_TABLES + 1>, FETCH_WIDTH> &actual_match1,
      arr<val<1>, FETCH_WIDTH> &primary_wrong) {
    return [&](u64 offset) -> val<1> {
      val<1> bim_primary = actual_match1[offset] >> NUM_TABLES;
      return execute_if(bim_primary.fo1() & primary_wrong[offset], [&]() {
        return ~bhyst[offset].read(bindex);
      });
    };
  }

  // Part 2: Write bimodal tables (AFTER need_extra_cycle, in the extra cycle).
  // Gate bim writes with extra_cycle since bim is read in predict2 (same cycle).
  // bhyst is in UPDATE_ONLY zone and not read in predict2, so no gating needed.
  void write_bimodal(arr<val<1>, FETCH_WIDTH> &is_branch,
                     arr<val<1>, FETCH_WIDTH> &branch_taken,
                     arr<val<1>, FETCH_WIDTH> &primary_wrong,
                     arr<val<1>, FETCH_WIDTH> &b_weak,
                     val<1> extra_cycle) {
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      val<1> bim_primary = match1[offset] >> NUM_TABLES;

      // Flip prediction if weak and wrong (gated by extra_cycle for bim RAM)
      execute_if(extra_cycle & b_weak[offset],
                 [&]() { bim[offset].write(bindex, branch_taken[offset]); });

      // Update hysteresis if bimodal is primary (bhyst is UPDATE_ONLY zone)
      execute_if(is_branch[offset] & bim_primary, [&]() {
        bhyst[offset].write(bindex, ~primary_wrong[offset]);
      });
    }
  }

  // ======== Phase 2: Meta Update (§2.11) ========

  // Update the meta counter (alt-on-newly-allocated).
  // Increment when alt was correct for a newly-allocated provider,
  // decrement when alt was wrong.
  void update_meta(arr<val<1>, FETCH_WIDTH> &is_branch,
                   arr<val<1>, FETCH_WIDTH> &branch_taken) {
    if constexpr (USE_META_V) {
      // For each offset: should we update meta?
      arr<val<1>, FETCH_WIDTH> altdiff = [&](u64 offset) {
        return (match2[offset] != hard<0>{}) & (pred2[offset] != pred1[offset]);
      };

      // Compute signed increment per offset:
      // concat(bad_pred2, 1): bad_pred2=0 → +1 (primary better, trust primary)
      //                       bad_pred2=1 → -1 (alt better, trust alt)
      // 0 when no update needed.
      arr<val<2, i64>, FETCH_WIDTH> meta_incr = [&](u64 offset) -> val<2, i64> {
        val<1> do_update =
            is_branch[offset] & altdiff[offset].fo1() & newly_alloc[offset];
        val<1> bad_pred2 = (pred2[offset] != branch_taken[offset]);
        return select(do_update.fo1(), concat(bad_pred2.fo1(), val<1>{1}),
                      val<2>{0});
      };

      if constexpr (META_TABLE_SIZE_V == 0) {
        // Global counter: pipeline shift and saturating add
        for (u64 i = Base::METAPIPE - 1; i != 0; i--) {
          meta_global[i] = meta_global[i - 1];
        }
        auto newmeta = meta_global[0] + meta_incr.fo1().fold_add();
        newmeta.fanout(hard<3>{});
        using meta_t = valt<decltype(meta_global[0])>;
        meta_global[0] =
            select(newmeta > meta_t::maxval, meta_t{meta_t::maxval},
                   select(newmeta < meta_t::minval, meta_t{meta_t::minval},
                          meta_t{newmeta}));
      } else {
        // PC-indexed table mode — TODO: read/write meta_table
      }
    }
  }

  // ======== Phase 2: History Update (§2.2) ========

  // Update global history, fold registers, P1 history, and path history.
  // Called unconditionally (caller gates with execute_if when needed).
  void update_history(val<64> next_pc) {
    val<PATHBITS> branch_bits = next_pc >> 2;
    if constexpr (P1_USE_GSHARE_V) {
      global_history1 = (global_history1 << 1) ^ val<P1_HIST_V>{next_pc >> 2};
    }
    gfolds.update(branch_bits);
    if constexpr (USE_PATH_HIST_V) {
      path_hist =
          (path_hist << PATHBITS) ^ val<PATH_HIST_WIDTH_V>{next_pc >> 2};
    }
  }

  // ======== Predictor Interface ========

  void new_block(val<64> inst_pc) {
    val<LOG_FETCH_WIDTH> offset = inst_pc.fo1() >> 2;
    block_entry = offset.fo1().decode().concat();
    block_entry.fanout(hard<6 * FETCH_WIDTH>{});
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

  val<1> predict2(val<64> inst_pc) override {
    val<LINEADDR_BITS> lineaddr = inst_pc >> (LOG_FETCH_WIDTH + 2);
    lineaddr.fanout(hard<1 + NUM_TABLES * 2>{});
    gfolds.fanout(hard<2>{});

    compute_indexes(lineaddr);
    compute_tags(lineaddr);
    read_tables();
    compute_matches();
    compute_meta_prediction();

    p2.fanout(hard<FETCH_WIDTH>{});
    val<1> taken = (block_entry & p2) != hard<0>{};
    taken.fanout(hard<2>{});
    reuse_prediction(~val<1>{block_entry >> (FETCH_WIDTH - 1)});
    return taken;
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) override {
    val<1> taken = ((block_entry << block_size) & p2) != hard<0>{};
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
    num_branch++;
  }

  void update_cycle(instruction_info &block_end_info) {
    val<1> &mispredict = block_end_info.is_mispredict;
    val<64> &next_pc = block_end_info.next_pc;

    // --- No conditional branch: just update history if needed ---
    if (num_branch == 0) {
      val<1> line_end = block_entry >> (FETCH_WIDTH - block_size);
      val<1> actual_block = ~(true_block & line_end.fo1());
      actual_block.fanout(hard<MAXHIST + NUM_TABLES * 2 + 2>{});
      execute_if(actual_block, [&]() {
        next_pc.fanout(hard<3>{});
        update_history(next_pc);
        true_block = 1;
      });
      return;
    }

    // --- Fanout cached state for update logic ---
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
    readc.fanout(hard<3>{});
    readh.fanout(hard<2>{});
    match1.fanout(hard<3>{});
    match2.fanout(hard<2>{});
    pred1.fanout(hard<2>{});
    pred2.fanout(hard<2 + NUM_TABLES>{});
    branch_offset.fanout(hard<FETCH_WIDTH + NUM_TABLES + 1>{});
    branch_dir.fanout(hard<2>{});
    gfolds.fanout(hard<2>{});
    if constexpr (USE_META_V) {
      if constexpr (META_TABLE_SIZE_V == 0) {
        meta_global.fanout(hard<2>{});
      }
      newly_alloc.fanout(hard<2>{});
    }

    val<LOG_FETCH_WIDTH> last_offset = branch_offset[num_branch - 1];
    last_offset.fanout(hard<4 * NUM_TABLES + 2>{});

    // --- Build per-offset masks ---
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

    // --- Per-offset: actual provider match (zero if no branch at offset) ---
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

    // --- Allocation ---
    auto alloc = allocate_entries(mispredict, last_offset);

    // --- Associate branch direction to each table ---
    alloc.allocate.fanout(hard<7>{});
    arr<val<1>, NUM_TABLES> bdir = [&](u64 i) {
      if constexpr (BR_P_ENTRY == 1) {
        val<LOG_FETCH_WIDTH> tag_offset = readt[i] >> MAX_HTAGBITS;
        val<LOG_FETCH_WIDTH> offset =
            select(alloc.allocate[i], last_offset, tag_offset.fo1());
        offset.fanout(hard<FETCH_WIDTH>{});
        arr<val<1>, FETCH_WIDTH> match_offset = [&](u64 j) {
          return branch_offset[j] == offset;
        };
        return (match_offset.fo1().concat() & update_valid & actualdirs) !=
               hard<0>{};
      } else {
        // BR_P_ENTRY>1: direction from the branch that owns this table's slot
        return branch_taken[0]; // simplified; full multi-slot TODO
      }
    };
    bdir.fanout(hard<2>{});

    // --- Per-table: is prediction incorrect? ---
    arr<val<1>, NUM_TABLES> badpred1 = [&](u64 i) {
      return readc[i] != bdir[i];
    };
    badpred1.fanout(hard<3>{});

    // --- Per-table: does primary differ from alt? ---
    arr<val<1>, NUM_TABLES> altdiffer = [&](u64 i) {
      if constexpr (BR_P_ENTRY == 1) {
        val<LOG_FETCH_WIDTH> tag_offset = readt[i] >> MAX_HTAGBITS;
        return readc[i] != pred2.select(tag_offset.fo1());
      } else {
        return readc[i] != pred2[0]; // simplified for multi-slot
      }
    };

    // --- Per-table: is owning branch's prediction correct? ---
    arr<val<1>, NUM_TABLES> goodpred = [&](u64 i) {
      if constexpr (BR_P_ENTRY == 1) {
        val<LOG_FETCH_WIDTH> tag_offset = readt[i] >> MAX_HTAGBITS;
        return (tag_offset.fo1() != last_offset) | correct_pred;
      } else {
        return correct_pred; // simplified for multi-slot
      }
    };

    // --- Weak detection for primary global predictions ---
    arr<val<1>, NUM_TABLES> g_weak = [&](u64 i) -> val<1> {
      return primary[i] & badpred1[i] & (readh[i] == hard<0>{});
    };

    // --- P1/P2 disagreement (for P1 update) ---
    val<FETCH_WIDTH> disagree_mask = (p1 ^ p2) & branch_mask.fo1();
    disagree_mask.fanout(hard<2>{});
    arr<val<1>, FETCH_WIDTH> disagree = disagree_mask.make_array(val<1>{});
    disagree.fanout(hard<2>{});

    arr<val<1>, FETCH_WIDTH> p1_weak = [&](u64 offset) -> val<1> {
      return execute_if(disagree[offset],
                        [&]() { return ~table1_hyst[offset].read(index1); });
    };

    // --- Bimodal hysteresis read (before extra cycle, same cycle as predict2) ---
    auto b_weak = compute_b_weak(actual_match1, primary_wrong);

    // --- Extra cycle for prediction bit flips and allocation ---
    val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
    val<1> extra_cycle =
        some_badpred1.fo1() | mispredict | (disagree_mask != hard<0>{});
    extra_cycle.fanout(hard<NUM_TABLES * 2 + 1>{});
    need_extra_cycle(extra_cycle);

    // --- Meta update (§2.11) ---
    update_meta(is_branch, branch_taken);

    // --- U-bit clear: if no allocation candidates, clear u-bits on post
    // entries ---
    val<1> noalloc = (alloc.candallocmask == hard<0>{});
    val<NUM_TABLES> uclearmask =
        alloc.postmask & noalloc.fo1().replicate(hard<NUM_TABLES>{}).concat();
    arr<val<1>, NUM_TABLES> uclear = uclearmask.fo1().make_array(val<1>{});
    uclear.fanout(hard<2>{});

    // --- TAGE table writes (§2.9) ---
    update_tables(alloc.allocate, bdir, badpred1, goodpred, primary, altdiffer,
                  uclear, g_weak, last_offset, extra_cycle);

    // --- P1 update: flip prediction if P1/P2 disagree and P1 is weak ---
    auto p2_split = p2.make_array(val<1>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(p1_weak[offset].fo1(), [&]() {
        table1_pred[offset].write(index1, p2_split[offset].fo1());
      });
    }
    // P1 hysteresis: agreement strengthens, disagreement weakens
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(is_branch[offset], [&]() {
        table1_hyst[offset].write(index1, ~disagree[offset]);
      });
    }

    // --- Bimodal writes (§2.10, after extra cycle) ---
    write_bimodal(is_branch, branch_taken, primary_wrong, b_weak, extra_cycle);

    // --- U-bit epoch reset (U_STOR_FF mode) ---
    if constexpr (U_STOR_FF_V) {
      uctr.fanout(hard<3>{});
      val<NUM_TABLES> allocmask1 = alloc.collamask1.reverse();
      allocmask1.fanout(hard<2>{});
      val<1> faralloc = (((alloc.last_match1 >> 3) | allocmask1).one_hot() ^
                         allocmask1) == hard<0>{};
      val<1> uctrsat = (uctr == hard<decltype(uctr)::maxval>{});
      uctrsat.fanout(hard<2>{});
      uctr = select(correct_pred, uctr,
                    select(uctrsat, val<decltype(uctr)::size>{0},
                           update_ctr(uctr, faralloc.fo1())));
      execute_if(uctrsat, [&]() {
        static_loop<NUM_TABLES>([&]<u64 I>() {
          std::get<I>(tables).reset_u(val<DefaultResetFn::MODE_BITS>{0});
        });
      });
    }

    // --- Global history update ---
    val<1> line_end = block_entry >> (FETCH_WIDTH - block_size);
    true_block = correct_pred | branch_dir[num_branch - 1] | line_end.fo1();
    true_block.fanout(hard<MAXHIST + NUM_TABLES * 2 + 2>{});
    execute_if(true_block, [&]() {
      next_pc.fanout(hard<3>{});
      update_history(next_pc);
    });

    num_branch = 0;
  }
};

// ============================================================================
// TageImpl — Ahead (1-branch-ahead) specialization
// ============================================================================

template <typename TableCfg, typename AllocCfg, u64 FETCH_WIDTH_V,
          u64 BIMODAL_SIZE_V, u64 BR_P_ENTRY_V, u64 NUM_BANKS_V,
          bool SHARED_TAG_V, bool SHARED_U_V, bool U_STOR_FF_V, u64 DECAY_CTR_V,
          typename ResetFn_V, bool USE_FF_CACHE_V, bool P1_USE_GSHARE_V,
          u64 P1_TABLE_SIZE_V, u64 P1_HIST_V, bool USE_META_V, u64 METABITS_V,
          u64 METAPIPE_V, u64 META_TABLE_SIZE_V, bool USE_PATH_HIST_V,
          u64 PATH_HIST_WIDTH_V, u64 PATH_BITS_V>
struct TageImpl<true, TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V,
                BR_P_ENTRY_V, NUM_BANKS_V, SHARED_TAG_V, SHARED_U_V,
                U_STOR_FF_V, DECAY_CTR_V, ResetFn_V, USE_FF_CACHE_V,
                P1_USE_GSHARE_V, P1_TABLE_SIZE_V, P1_HIST_V, USE_META_V,
                METABITS_V, METAPIPE_V, META_TABLE_SIZE_V, USE_PATH_HIST_V,
                PATH_HIST_WIDTH_V, PATH_BITS_V> : predictor {

  using Base =
      TageBase<TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V, BR_P_ENTRY_V,
               NUM_BANKS_V, true, SHARED_TAG_V, SHARED_U_V, U_STOR_FF_V,
               DECAY_CTR_V, ResetFn_V, USE_FF_CACHE_V, P1_USE_GSHARE_V,
               P1_TABLE_SIZE_V, P1_HIST_V, USE_META_V, METABITS_V, METAPIPE_V,
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
                     hcm::ram<val<METABITS_V, i64>, META_RAM_SIZE>, EmptyMember>
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

  void
  update_cycle([[maybe_unused]] instruction_info &block_end_info) override {
    // TODO: Phase 5.5 — ahead update_cycle (write to index[1], path banking)
    // TODO: Phase 8 — loop predictor override
  }
};

// ============================================================================
// User-facing Alias
// ============================================================================

// Template parameter names prefixed with TAGE_ to avoid collisions with
// macros defined in other headers (e.g. tage.hpp's #define USE_META).
template <typename TableCfg = DefaultTableConfig,
          typename AllocCfg = DefaultAllocConfig,
          // Global hardware params
          u64 TAGE_FETCH_WIDTH = 8, u64 TAGE_BIMODAL_SIZE = 4096,
          u64 TAGE_BR_P_ENTRY = 1, u64 TAGE_NUM_BANKS = 1,
          bool TAGE_USE_AHEAD = false, bool TAGE_SHARED_TAG = true,
          bool TAGE_SHARED_U = true, bool TAGE_U_STOR_FF = false,
          u64 TAGE_DECAY_CTR = 1024, typename ResetFn = DefaultResetFn,
          bool TAGE_USE_FF_CACHE = false,
          // P1 params
          bool TAGE_P1_USE_GSHARE = true, u64 TAGE_P1_TABLE_SIZE = 16384,
          u64 TAGE_P1_HIST = 6,
          // Meta-prediction
          bool TAGE_USE_META = true, u64 TAGE_METABITS = 4,
          u64 TAGE_METAPIPE = 2, u64 TAGE_META_TABLE_SIZE = 0,
          // Path history
          bool TAGE_USE_PATH_HIST = false, u64 TAGE_PATH_HIST_WIDTH = 27,
          u64 TAGE_PATH_BITS = 6>
using Tage =
    TageImpl<TAGE_USE_AHEAD, TableCfg, AllocCfg, TAGE_FETCH_WIDTH,
             TAGE_BIMODAL_SIZE, TAGE_BR_P_ENTRY, TAGE_NUM_BANKS,
             TAGE_SHARED_TAG, TAGE_SHARED_U, TAGE_U_STOR_FF, TAGE_DECAY_CTR,
             ResetFn, TAGE_USE_FF_CACHE, TAGE_P1_USE_GSHARE, TAGE_P1_TABLE_SIZE,
             TAGE_P1_HIST, TAGE_USE_META, TAGE_METABITS, TAGE_METAPIPE,
             TAGE_META_TABLE_SIZE, TAGE_USE_PATH_HIST, TAGE_PATH_HIST_WIDTH,
             TAGE_PATH_BITS>;
