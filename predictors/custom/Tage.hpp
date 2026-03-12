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
  static constexpr auto CTR_WIDTH = uniform_array<u64, NUM_TABLES>(1);
  static constexpr auto HYST_WIDTH = uniform_array<u64, NUM_TABLES>(2);
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
          bool SHARED_TAG, bool SHARED_U, bool SHARED_HYS,
          bool U_STOR_FF, u64 DECAY_CTR,
          typename ResetFn, bool USE_FF_CACHE, u64 I>
using TableAt =
    TageTable<Cfg::TABLE_SIZE[I], Cfg::HIST_LEN[I], Cfg::TAG_WIDTH[I],
              Cfg::CTR_WIDTH[I], Cfg::HYST_WIDTH[I], Cfg::U_WIDTH[I],
              BR_P_ENTRY, NUM_BANKS,
              USE_AHEAD, SHARED_TAG, SHARED_U, SHARED_HYS,
              U_STOR_FF, DECAY_CTR, ResetFn, USE_FF_CACHE>;

// Build a std::tuple of TageTable types from the config.
template <typename Cfg, u64 BR_P_ENTRY, u64 NUM_BANKS, bool USE_AHEAD,
          bool SHARED_TAG, bool SHARED_U, bool SHARED_HYS,
          bool U_STOR_FF, u64 DECAY_CTR,
          typename ResetFn, bool USE_FF_CACHE, typename Seq>
struct MakeTableTuple;

template <typename Cfg, u64 BR_P_ENTRY, u64 NUM_BANKS, bool USE_AHEAD,
          bool SHARED_TAG, bool SHARED_U, bool SHARED_HYS,
          bool U_STOR_FF, u64 DECAY_CTR,
          typename ResetFn, bool USE_FF_CACHE, u64... Is>
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

template <bool USE_AHEAD_V, typename TableCfg, typename AllocCfg,
          // Global hardware params
          u64 FETCH_WIDTH_V, u64 BIMODAL_SIZE_V, u64 BR_P_ENTRY_V,
          u64 NUM_BANKS_V, bool SHARED_TAG_V, bool SHARED_U_V,
          bool SHARED_HYS_V, bool U_STOR_FF_V,
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
          bool SHARED_HYS_V, bool U_STOR_FF_V, u64 DECAY_CTR_V,
          typename ResetFn_V, bool USE_FF_CACHE_V,
          bool P1_USE_GSHARE_V, u64 P1_TABLE_SIZE_V,
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
  using Tables =
      typename MakeTableTuple<TableCfg, BR_P_ENTRY_V, NUM_BANKS_V, USE_AHEAD_V,
                              SHARED_TAG_V, SHARED_U_V, SHARED_HYS_V,
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
// Flat inline implementation matching tage_tt structure for optimal EPI/latency.
// All prediction and update logic is inline — no member function decomposition.

template <typename TableCfg, typename AllocCfg, u64 FETCH_WIDTH_V,
          u64 BIMODAL_SIZE_V, u64 BR_P_ENTRY_V, u64 NUM_BANKS_V,
          bool SHARED_TAG_V, bool SHARED_U_V, bool SHARED_HYS_V,
          bool U_STOR_FF_V, u64 DECAY_CTR_V,
          typename ResetFn_V, bool USE_FF_CACHE_V, bool P1_USE_GSHARE_V,
          u64 P1_TABLE_SIZE_V, u64 P1_HIST_V, bool USE_META_V, u64 METABITS_V,
          u64 METAPIPE_V, u64 META_TABLE_SIZE_V, bool USE_PATH_HIST_V,
          u64 PATH_HIST_WIDTH_V, u64 PATH_BITS_V>
struct TageImpl<false, TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V,
                BR_P_ENTRY_V, NUM_BANKS_V, SHARED_TAG_V, SHARED_U_V,
                SHARED_HYS_V, U_STOR_FF_V, DECAY_CTR_V, ResetFn_V,
                USE_FF_CACHE_V, P1_USE_GSHARE_V, P1_TABLE_SIZE_V, P1_HIST_V,
                USE_META_V, METABITS_V, METAPIPE_V, META_TABLE_SIZE_V,
                USE_PATH_HIST_V, PATH_HIST_WIDTH_V, PATH_BITS_V> : predictor {

  using Base =
      TageBase<TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V, BR_P_ENTRY_V,
               NUM_BANKS_V, false, SHARED_TAG_V, SHARED_U_V, SHARED_HYS_V,
               U_STOR_FF_V, DECAY_CTR_V, ResetFn_V, USE_FF_CACHE_V,
               P1_USE_GSHARE_V, P1_TABLE_SIZE_V, P1_HIST_V, USE_META_V,
               METABITS_V, METAPIPE_V, META_TABLE_SIZE_V, USE_PATH_HIST_V,
               PATH_HIST_WIDTH_V, PATH_BITS_V>;

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

  // ======== Table type ========
  using Table0 = std::tuple_element_t<0, typename Base::Tables>;

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
  std::conditional_t<USE_META_V,
      arr<reg<METABITS_V, i64>, METAPIPE_V>, EmptyMember> meta;
  std::conditional_t<USE_META_V,
      arr<reg<1>, FETCH_WIDTH>, EmptyMember> newly_alloc;

  // U-bit reset
  reg<UCTRBITS> uctr;

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

  // P2 TAGE tables (TageTable storage)
  Table0 table[NUM_TABLES];

  // P2 bimodal
  hcm::ram<val<1>, BIM_ENTRIES> bim[FETCH_WIDTH]{"bpred"};

  zone UPDATE_ONLY;
  hcm::ram<val<1>, P1_ENTRIES> table1_hyst[FETCH_WIDTH]{"P1 hyst"};
  hcm::ram<val<1>, BIM_ENTRIES> bhyst[FETCH_WIDTH]{"bhyst"};

  // ======== Predictor Interface ========

  void new_block(val<64> inst_pc) {
    val<LOG_FETCH_WIDTH> offset = inst_pc.fo1() >> 2;
    block_entry = offset.fo1().decode().concat();
    block_entry.fanout(hard<6 * FETCH_WIDTH>{});
    block_size = 1;
  }

  val<1> predict1(val<64> inst_pc) override {
    inst_pc.fanout(hard<2>{});
    new_block(inst_pc);
    val<std::max(P1_INDEX_BITS, P1_HIST_V)> lineaddr =
        inst_pc >> (LOG_FETCH_WIDTH + 2);
    lineaddr.fanout(hard<2>{});
    if constexpr (P1_USE_GSHARE_V) {
      if constexpr (P1_HIST_V <= P1_INDEX_BITS) {
        index1 = lineaddr ^
                 (val<P1_INDEX_BITS>{global_history1}
                      << (P1_INDEX_BITS - P1_HIST_V));
      } else {
        index1 = global_history1.make_array(val<P1_INDEX_BITS>{})
                     .append(lineaddr)
                     .fold_xor();
      }
    } else {
      index1 = lineaddr;
    }
    index1.fanout(hard<FETCH_WIDTH>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      readp1[offset] = table1_pred[offset].read(index1);
    }
    readp1.fanout(hard<2>{});
    p1 = readp1.concat();
    p1.fanout(hard<FETCH_WIDTH>{});
    return (block_entry & p1) != hard<0>{};
  }

  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) override {
    return ((block_entry << block_size) & p1) != hard<0>{};
  }

  val<1> predict2(val<64> inst_pc) override {
    val<LINEADDR_BITS> lineaddr = inst_pc >> (LOG_FETCH_WIDTH + 2);
    lineaddr.fanout(hard<1 + NUM_TABLES * 2>{});
    gfolds.fanout(hard<2>{});

    // compute indexes
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
    gindex.fanout(hard<4>{});

    // compute hashed tags
    for (u64 i = 0; i < NUM_TABLES; i++) {
      htag[i] =
          val<MAX_HTAGBITS>{lineaddr}.reverse() ^ gfolds.template get<1>(i);
    }
    htag.fanout(hard<2>{});

    // read tables
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      readb[offset] = bim[offset].read(bindex);
    }
    readb.fanout(hard<2>{});
    for (u64 i = 0; i < NUM_TABLES; i++) {
      readt[i] = table[i].tag_ram[0].read(gindex[i]);
      readc[i] = table[i].pred_ram[0].read(gindex[i]);
      if constexpr (MAX_HYST_WIDTH > 0) {
        readh[i] = table[i].hyst_ram[0].read(gindex[i]);
      }
      if constexpr (U_STOR_FF_V) {
        readu[i] = table[i].u_ff[0].select(gindex[i]);
      } else {
        readu[i] = table[i].u_ram[0].read(gindex[i]);
      }
    }
    readt.fanout(hard<FETCH_WIDTH + 1>{});
    readc.fanout(hard<3>{});
    if constexpr (MAX_HYST_WIDTH > 0) {
      readh.fanout(hard<2>{});
    }
    readu.fanout(hard<2>{});
    notumask = ~readu.concat();
    notumask.fanout(hard<2>{});

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
        match[offset] =
            concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp);
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
    match1.fanout(hard<3>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      pred1[offset] = (match1[offset] & preds[offset]) != hard<0>{};
    }
    pred1.fanout(hard<2>{});

    // for each offset, find second longest match and select secondary prediction
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      match2[offset] = (match[offset] ^ match1[offset]).one_hot();
    }
    match2.fanout(hard<2>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      pred2[offset] = (match2[offset] & preds[offset]) != hard<0>{};
    }
    pred2.fanout(hard<2>{});

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
      newly_alloc.fanout(hard<2>{});
      arr<val<1>, FETCH_WIDTH> altsel = [&](u64 offset) {
        arr<val<1>, 3> inputs = {metasign, newly_alloc[offset],
                                 match2[offset] != hard<0>{}};
        return inputs.fo1().fold_and();
      };
      p2 = arr<val<1>, FETCH_WIDTH>{[&](u64 offset) {
             return select(altsel[offset].fo1(), pred2[offset], pred1[offset]);
           }}.concat();
    } else {
      p2 = pred1.concat();
    }

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

    arr<val<NUM_TABLES + 1>, FETCH_WIDTH> actual_match1 = [&](u64 offset) {
      return select(is_branch[offset], match1[offset],
                    val<NUM_TABLES + 1>{0});
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
    val<NUM_TABLES> collamask12 = select(
        val<2>{std::rand()} == hard<0>{}, collamask2.fo1(), collamask1);
    arr<val<1>, NUM_TABLES> allocate =
        collamask12.fo1().reverse().make_array(val<1>{});
    allocate.fanout(hard<7>{});

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
        return (readc[i] >> hard<MAX_CTR_WIDTH - 1>{}) != bdir[i];
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

    // need extra cycle for modifying prediction bits and for TAGE allocation
    val<1> some_badpred1 =
        (primary_mask & badpred1.concat()) != hard<0>{};
    val<1> extra_cycle =
        some_badpred1.fo1() | mispredict | (disagree_mask != hard<0>{});
    extra_cycle.fanout(hard<NUM_TABLES * 2 + 1>{});
    need_extra_cycle(extra_cycle);

    // update meta counter
    if constexpr (USE_META_V) {
      arr<val<1>, FETCH_WIDTH> altdiff = [&](u64 offset) {
        return (match2[offset] != hard<0>{}) &
               (pred2[offset] != pred1[offset]);
      };
      arr<val<2, i64>, FETCH_WIDTH> meta_incr =
          [&](u64 offset) -> val<2, i64> {
        val<1> update_meta = is_branch[offset] & altdiff[offset].fo1() &
                             newly_alloc[offset];
        val<1> bad_pred2 = (pred2[offset] != branch_taken[offset]);
        return select(update_meta.fo1(),
                      concat(bad_pred2.fo1(), val<1>{1}), val<2>{0});
      };
      for (u64 i = METAPIPE_V - 1; i != 0; i--) {
        meta[i] = meta[i - 1];
      }
      auto newmeta = meta[0] + meta_incr.fo1().fold_add();
      newmeta.fanout(hard<3>{});
      using meta_t = valt<decltype(meta[0])>;
      meta[0] =
          select(newmeta > meta_t::maxval, meta_t{meta_t::maxval},
                 select(newmeta < meta_t::minval, meta_t{meta_t::minval},
                        meta_t{newmeta}));
    }

    // overwrite the tag in the allocated entry (mispredict)
    for (u64 i = 0; i < NUM_TABLES; i++) {
      execute_if(allocate[i], [&]() {
        if constexpr (BR_P_ENTRY_V == 1) {
          table[i].tag_ram[0].write(gindex[i],
                                    concat(last_offset, htag[i]));
        } else {
          table[i].tag_ram[0].write(gindex[i],
                                    val<MAX_TAG_WIDTH>{htag[i]});
        }
      });
    }

    // update the u bits
    arr<val<1>, NUM_TABLES> update_u = [&](u64 i) {
      return primary[i] & altdiffer[i].fo1();
    };
    // if all post entries have the u bit set, reset their u bits
    val<1> noalloc = (candallocmask == hard<0>{});
    val<NUM_TABLES> uclearmask =
        postmask & noalloc.fo1().replicate(hard<NUM_TABLES>{}).concat();
    arr<val<1>, NUM_TABLES> uclear =
        uclearmask.fo1().make_array(val<1>{});
    uclear.fanout(hard<2>{});
    for (u64 i = 0; i < NUM_TABLES; i++) {
      execute_if(update_u[i].fo1() | allocate[i] | uclear[i], [&]() {
        val<1> newu = goodpred[i].fo1() & ~allocate[i] & ~uclear[i];
        if constexpr (U_STOR_FF_V) {
          table[i].write_u_ff_arr(0, gindex[i], newu);
        } else {
          table[i].u_ram[0].write(gindex[i], newu.fo1(), extra_cycle);
        }
      });
    }

    // update P1 prediction if P1 and P2 disagree and the hysteresis bit is weak
    auto p2_split = p2.make_array(val<1>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(p1_weak[offset].fo1(), [&]() {
        table1_pred[offset].write(index1, p2_split[offset].fo1());
      });
    }
    // update P1 hysteresis
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(is_branch[offset], [&]() {
        table1_hyst[offset].write(index1, ~disagree[offset]);
      });
    }

    // update incorrect bimodal prediction if primary provider and hyst is weak
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(b_weak[offset].fo1(), [&]() {
        bim[offset].write(bindex, branch_taken[offset]);
      });
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
    for (u64 i = 0; i < NUM_TABLES; i++) {
      execute_if(g_weak[i].fo1() | allocate[i], [&]() {
        table[i].pred_ram[0].write(gindex[i], bdir[i]);
      });
    }
    // update global prediction hysteresis if primary provider or allocated entry
    if constexpr (MAX_HYST_WIDTH > 0) {
      for (u64 i = 0; i < NUM_TABLES; i++) {
        execute_if(primary[i] | allocate[i], [&]() {
          auto newhyst = select(allocate[i],
                                val<std::max(u64(1), MAX_HYST_WIDTH)>{0},
                                update_ctr(readh[i], ~badpred1[i]));
          table[i].hyst_ram[0].write(gindex[i], newhyst.fo1(), extra_cycle);
        });
      }
    }

    // u-bit epoch reset
    uctr.fanout(hard<3>{});
    val<NUM_TABLES> allocmask1 = collamask1.reverse();
    allocmask1.fanout(hard<2>{});
    val<1> faralloc =
        (((last_match1 >> 3) | allocmask1).one_hot() ^ allocmask1) ==
        hard<0>{};
    val<1> uctrsat = (uctr == hard<decltype(uctr)::maxval>{});
    uctrsat.fanout(hard<2>{});
    uctr = select(correct_pred, uctr,
                  select(uctrsat, val<decltype(uctr)::size>{0},
                         update_ctr(uctr, faralloc.fo1())));
    if constexpr (U_STOR_FF_V) {
      execute_if(uctrsat, [&]() {
        for (u64 i = 0; i < NUM_TABLES; i++)
          table[i].reset_u(val<ResetFn_V::MODE_BITS>{0});
      });
    } else {
      execute_if(uctrsat, [&]() {
        for (auto &t : table)
          t.u_ram[0].reset();
      });
    }

    // update global history
    val<1> line_end = block_entry >> (FETCH_WIDTH - block_size);
    true_block = correct_pred | branch_dir[num_branch - 1] | line_end.fo1();
    true_block.fanout(hard<MAXHIST + NUM_TABLES * 2 + 2>{});
    execute_if(true_block, [&]() {
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
    });

    num_branch = 0; // done
  }
};

// ============================================================================
// TageImpl — Ahead (1-branch-ahead) specialization
// ============================================================================

template <typename TableCfg, typename AllocCfg, u64 FETCH_WIDTH_V,
          u64 BIMODAL_SIZE_V, u64 BR_P_ENTRY_V, u64 NUM_BANKS_V,
          bool SHARED_TAG_V, bool SHARED_U_V, bool SHARED_HYS_V,
          bool U_STOR_FF_V, u64 DECAY_CTR_V,
          typename ResetFn_V, bool USE_FF_CACHE_V, bool P1_USE_GSHARE_V,
          u64 P1_TABLE_SIZE_V, u64 P1_HIST_V, bool USE_META_V, u64 METABITS_V,
          u64 METAPIPE_V, u64 META_TABLE_SIZE_V, bool USE_PATH_HIST_V,
          u64 PATH_HIST_WIDTH_V, u64 PATH_BITS_V>
struct TageImpl<true, TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V,
                BR_P_ENTRY_V, NUM_BANKS_V, SHARED_TAG_V, SHARED_U_V,
                SHARED_HYS_V, U_STOR_FF_V, DECAY_CTR_V, ResetFn_V,
                USE_FF_CACHE_V, P1_USE_GSHARE_V, P1_TABLE_SIZE_V, P1_HIST_V,
                USE_META_V, METABITS_V, METAPIPE_V, META_TABLE_SIZE_V,
                USE_PATH_HIST_V, PATH_HIST_WIDTH_V, PATH_BITS_V> : predictor {

  using Base =
      TageBase<TableCfg, AllocCfg, FETCH_WIDTH_V, BIMODAL_SIZE_V, BR_P_ENTRY_V,
               NUM_BANKS_V, true, SHARED_TAG_V, SHARED_U_V, SHARED_HYS_V,
               U_STOR_FF_V, DECAY_CTR_V, ResetFn_V, USE_FF_CACHE_V,
               P1_USE_GSHARE_V, P1_TABLE_SIZE_V, P1_HIST_V, USE_META_V,
               METABITS_V, METAPIPE_V, META_TABLE_SIZE_V, USE_PATH_HIST_V,
               PATH_HIST_WIDTH_V, PATH_BITS_V>;

  static constexpr u64 NUM_TABLES = Base::NUM_TABLES;
  static constexpr u64 FETCH_WIDTH = Base::FETCH_WIDTH;
  static constexpr u64 LOG_FETCH_WIDTH = Base::LOG_FETCH_WIDTH;
  static constexpr u64 BR_P_ENTRY = Base::BR_P_ENTRY;
  static constexpr u64 MAX_TAG_WIDTH = Base::MAX_TAG_WIDTH;
  static constexpr u64 MAX_CTR_WIDTH = Base::MAX_CTR_WIDTH;
  static constexpr u64 MAX_HYST_WIDTH = Base::MAX_HYST_WIDTH;
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
  reg<UCTRBITS> uctr;

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
          u64 TAGE_FETCH_WIDTH = 16, u64 TAGE_BIMODAL_SIZE = 4096,
          u64 TAGE_BR_P_ENTRY = 1, u64 TAGE_NUM_BANKS = 1,
          bool TAGE_USE_AHEAD = false, bool TAGE_SHARED_TAG = true,
          bool TAGE_SHARED_U = true, bool TAGE_SHARED_HYS = true,
          bool TAGE_U_STOR_FF = false,
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
             TAGE_SHARED_TAG, TAGE_SHARED_U, TAGE_SHARED_HYS,
             TAGE_U_STOR_FF, TAGE_DECAY_CTR,
             ResetFn, TAGE_USE_FF_CACHE, TAGE_P1_USE_GSHARE, TAGE_P1_TABLE_SIZE,
             TAGE_P1_HIST, TAGE_USE_META, TAGE_METABITS, TAGE_METAPIPE,
             TAGE_META_TABLE_SIZE, TAGE_USE_PATH_HIST, TAGE_PATH_HIST_WIDTH,
             TAGE_PATH_BITS>;
