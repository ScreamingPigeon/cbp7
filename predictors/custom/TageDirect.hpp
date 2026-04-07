#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"
#ifdef TAGE_MONITOR
#include "TageMonitor.hpp"
#endif

using namespace hcm;

// ============================================================================
// TageDirect — Branch-rank-indexed TAGE with separate LINEINST/N,
//              P1 gshare base predictor, and pluggable history folding.
//
// Key differences from Tage.hpp:
//   - Branch rank in tag (log2(LANES) bits) instead of instruction offset
//   - LINEINST and N are independent (LINEINST=line boundary, N=max branches)
//   - P1 gshare is the sole base predictor (no bimodal)
//   - Pluggable fold functors for history hashing
//   - MISPREDICT_ONLY_WRITE mode for reduced extra cycles
//   - Self-contained: no dependency on common.hpp or TageTable.hpp
// ============================================================================

// ============================================================================
// Constexpr Helpers
// ============================================================================

namespace td {

constexpr u64 clog2(u64 x) {
  u64 r = 0;
  u64 v = x - 1;
  while (v > 0) {
    v >>= 1;
    r++;
  }
  return r;
}

constexpr u64 next_pow2(u64 x) {
  u64 r = 1;
  while (r < x)
    r <<= 1;
  return r;
}

constexpr double constexpr_pow(double base, double exp) {
  if (exp == 0.0)
    return 1.0;
  if (base == 0.0)
    return 0.0;
  return std::exp(exp * std::log(base));
}

template <typename T, std::size_t N>
constexpr std::array<T, N> uniform_array(T v) {
  std::array<T, N> a{};
  for (std::size_t i = 0; i < N; i++)
    a[i] = v;
  return a;
}

template <typename T, std::size_t N>
constexpr T array_max(const std::array<T, N> &a) {
  T m = a[0];
  for (std::size_t i = 1; i < N; i++)
    m = (a[i] > m) ? a[i] : m;
  return m;
}

template <typename T, std::size_t N>
constexpr T array_min(const std::array<T, N> &a) {
  T m = a[0];
  for (std::size_t i = 1; i < N; i++)
    m = (a[i] < m) ? a[i] : m;
  return m;
}

// Saturating counter update (inspired by common.hpp update_ctr)
template <u64 N, typename T = u64>
val<N, T> update_ctr(val<N, T> &ctr, val<1> incr) {
  ctr.fanout(hard<3>{});
  auto maxv = val<N, T>::maxval;
  auto minv = val<N, T>::minval;
  return select(incr, select(ctr == maxv, val<N, T>{maxv}, val<N, T>{ctr + 1}),
                select(ctr == minv, val<N, T>{minv}, val<N, T>{ctr - 1}));
}

// ============================================================================
// History Series
// ============================================================================

enum class HistSeries { GEOMETRIC, QUADRATIC, SUPEREXP, ROS };

template <std::size_t N>
constexpr std::array<u64, N> geometric_hist(u64 min_h, u64 max_h) {
  static_assert(N >= 2);
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

template <std::size_t N>
constexpr std::array<u64, N> quadratic_hist(u64 minh, u64 maxh, u64 d = 2,
                                            u64 k = 1) {
  std::array<u64, N> h{};
  h[0] = minh;
  for (std::size_t n = 1; n < N; n++)
    h[n] = h[n - 1] + d * n + k;
  double raw_max = h[N - 1], raw_min = h[0];
  u64 prev = 0;
  for (std::size_t n = 0; n < N; n++) {
    double t = (raw_max > raw_min)
                   ? (double(h[n]) - raw_min) / (raw_max - raw_min)
                   : 0.0;
    u64 scaled = minh + u64(t * (maxh - minh));
    scaled = (scaled > prev + 1 || n == 0) ? scaled : prev + 1;
    h[n] = scaled;
    prev = scaled;
  }
  for (std::size_t i = 0; i < N / 2; i++)
    std::swap(h[i], h[N - 1 - i]);
  return h;
}

template <std::size_t N>
constexpr std::array<u64, N> superexp_hist(u64 h0, double f, double m) {
  std::array<u64, N> h{};
  h[0] = h0;
  for (std::size_t n = 1; n < N; n++) {
    double mult = f * double(n) + m;
    u64 next = u64(double(h[n - 1]) * mult + 0.5);
    h[n] = (next > h[n - 1] + 1) ? next : h[n - 1] + 1;
  }
  for (std::size_t i = 0; i < N / 2; i++)
    std::swap(h[i], h[N - 1 - i]);
  return h;
}

template <std::size_t N>
constexpr std::array<u64, N> ros_hist(u64 minh, u64 maxh, u64 d = 2, u64 k = 1,
                                      double f = 0.1, double m = 1.1,
                                      std::size_t t = 15) {
  std::array<u64, N> h{};
  h[0] = minh;
  for (std::size_t n = 1; n < N; n++) {
    if (n < t) {
      h[n] = h[n - 1] + d * n + k;
    } else {
      double mult = f * double(n - t + 1) + m;
      u64 next = u64(double(h[n - 1]) * mult + 0.5);
      h[n] = (next > h[n - 1] + 1) ? next : h[n - 1] + 1;
    }
  }
  double raw_max = h[N - 1], raw_min = h[0];
  u64 prev = 0;
  for (std::size_t n = 0; n < N; n++) {
    double t_val = (raw_max > raw_min)
                       ? (double(h[n]) - raw_min) / (raw_max - raw_min)
                       : 0.0;
    u64 scaled = minh + u64(t_val * (maxh - minh));
    scaled = (scaled > prev + 1 || n == 0) ? scaled : prev + 1;
    h[n] = scaled;
    prev = scaled;
  }
  for (std::size_t i = 0; i < N / 2; i++)
    std::swap(h[i], h[N - 1 - i]);
  return h;
}

// ============================================================================
// Tag Width Functors — (table_index, num_tables) -> tag_width
// ============================================================================

template <u64 TAG> struct UniformTag {
  constexpr u64 operator()(u64, u64) const { return TAG; }
};

template <u64 TAG_LONG, u64 TAG_SHORT> struct GradedTag {
  constexpr u64 operator()(u64 i, u64 n) const {
    if (n <= 1)
      return TAG_LONG;
    return TAG_LONG - (TAG_LONG - TAG_SHORT) * i / (n - 1);
  }
};

template <u64 TAG_HI, u64 TAG_LO, u64 SPLIT> struct StepTag {
  constexpr u64 operator()(u64 i, u64) const {
    return (i < SPLIT) ? TAG_HI : TAG_LO;
  }
};

template <u64 BASE, u64 SCALE> struct LogTag {
  constexpr u64 operator()(u64 i, u64 n) const {
    u64 level = (n > 1) ? (n - 1 - i) : 0;
    return BASE + SCALE * level / 4;
  }
};

template <u64 N, typename TagFn>
constexpr std::array<u64, N> generate_tag_widths(TagFn fn) {
  std::array<u64, N> t{};
  for (u64 i = 0; i < N; i++)
    t[i] = fn(i, N);
  return t;
}

// ============================================================================
// Size Functors — (table_index, num_tables) -> table_size
// ============================================================================

template <u64 SIZE> struct UniformSize {
  constexpr u64 operator()(u64, u64) const { return SIZE; }
};

template <u64 SIZE, u64 RATIO> struct GeoSize {
  constexpr u64 operator()(u64 i, u64 n) const {
    if (RATIO <= 1 || n <= 1)
      return SIZE;
    double t = double(i) / double(n - 1);
    double scale = constexpr_pow(double(RATIO), t - 0.5);
    u64 sz = u64(SIZE * scale);
    u64 result = 64;
    while (result < sz)
      result *= 2;
    return result;
  }
};

template <u64 S_HI, u64 S_LO, u64 SPLIT> struct StepSize {
  constexpr u64 operator()(u64 i, u64) const {
    return (i < SPLIT) ? S_HI : S_LO;
  }
};

template <u64 BASE_SIZE> struct SqrtHistSize {
  constexpr u64 operator()(u64 i, u64 n) const {
    if (n <= 1)
      return BASE_SIZE;
    double scale = constexpr_pow(double(n - 1 - i + 1), 0.5);
    u64 sz = u64(BASE_SIZE * scale);
    u64 result = 64;
    while (result < sz)
      result *= 2;
    return result;
  }
};

template <std::size_t N, typename SizeFn>
constexpr std::array<u64, N> generate_table_sizes(SizeFn fn) {
  std::array<u64, N> s{};
  for (std::size_t i = 0; i < N; i++)
    s[i] = fn(i, N);
  return s;
}

// ============================================================================
// Table Config
// ============================================================================

template <u64 N = 8, u64 SIZE = 2048, u64 TAG = 11, u64 CTR = 1, u64 HYST = 2,
          u64 U = 1, u64 MINH = 2, u64 MAXH = 100, u64 SIZE_RATIO = 1,
          HistSeries HIST = HistSeries::GEOMETRIC,
          typename TagFn = UniformTag<TAG>,
          typename SizeFn = GeoSize<SIZE, SIZE_RATIO>>
struct TDTableConfig {
  static constexpr u64 NUM_TABLES = N;
  static constexpr u64 MINHIST = MINH;
  static constexpr u64 MAXHIST = MAXH;
  static constexpr auto TABLE_SIZE = generate_table_sizes<N>(SizeFn{});
  static constexpr auto TAG_WIDTH = generate_tag_widths<N>(TagFn{});
  static constexpr auto CTR_WIDTH = uniform_array<u64, N>(CTR);
  static constexpr auto HYST_WIDTH = uniform_array<u64, N>(HYST);
  static constexpr auto U_WIDTH = uniform_array<u64, N>(U);
  static constexpr auto HIST_LEN = []() {
    if constexpr (HIST == HistSeries::GEOMETRIC)
      return geometric_hist<N>(MINH, MAXH);
    else if constexpr (HIST == HistSeries::QUADRATIC)
      return quadratic_hist<N>(MINH, MAXH);
    else if constexpr (HIST == HistSeries::SUPEREXP)
      return superexp_hist<N>(MINH, 0.1, 1.1);
    else if constexpr (HIST == HistSeries::ROS)
      return ros_hist<N>(MINH, MAXH);
  }();
};

// ============================================================================
// Allocation Policy Configs
// ============================================================================

struct TDDefaultAllocConfig {
  static constexpr u64 MAX_ALLOC = 1;
  static constexpr bool NON_CONSECUTIVE = false;
  static constexpr bool CONF_GATE = false;
  static constexpr u64 PROB_START = 0;
  static constexpr bool PARTIAL_UPDATE = false;
  static constexpr bool MISPREDICT_ONLY_WRITE = false;
};

struct TDMispredOnlyAllocConfig {
  static constexpr u64 MAX_ALLOC = 1;
  static constexpr bool NON_CONSECUTIVE = false;
  static constexpr bool CONF_GATE = false;
  static constexpr u64 PROB_START = 0;
  static constexpr bool PARTIAL_UPDATE = false;
  static constexpr bool MISPREDICT_ONLY_WRITE = true;
};

// ============================================================================
// Decay Threshold Adaptation Policies
// ============================================================================

struct TDDecayMild {
  template <u64 W>
  static val<W> apply(val<W> &t, val<1> &correct,
                      [[maybe_unused]] val<1> &uctrsat, val<1> &misp) {
    return select(misp, td::update_ctr(t, hard<0>{}),
                  select(correct, td::update_ctr(t, hard<1>{}), val<W>{t}));
  }
};

struct TDDecayAggressive {
  template <u64 W>
  static val<W> apply(val<W> &t, [[maybe_unused]] val<1> &correct,
                      [[maybe_unused]] val<1> &uctrsat, val<1> &misp) {
    return select(misp, td::update_ctr(t, hard<0>{}), val<W>{t});
  }
};

// ============================================================================
// History Fold Functors
// ============================================================================
// A fold functor maintains a reg<F> and provides:
//   get() -> val<F>
//   update(global_history, hist_len, new_bits)
//   fanout(fo)

// Standard XOR fold (same as Seznec folded_gh)
template <u64 F> struct XORFold {
  reg<F> folded;

  val<F> get() { return folded; }

  void fanout(auto fo) { folded.fanout(fo); }

  template <u64 MAXL> void update(auto &gh, u64 ghlen, auto in) {
    constexpr u64 inbits = std::min(F, u64(decltype(in)::size));
    val<inbits> input = in;
    auto f = folded.make_array(val<1>{});
    val<1> outbit = gh[ghlen - 1];
    u64 outpos = ghlen % F;
    arr<val<1>, F> ff = [&](u64 i) {
      if (i == 0)
        return (outpos == 0) ? f[F - 1].fo1() ^ outbit.fo1() : f[F - 1].fo1();
      else
        return (outpos == i) ? f[i - 1].fo1() ^ outbit.fo1() : f[i - 1].fo1();
    };
    auto x = input.fo1().make_array(val<1>{});
    arr<val<1>, F> y = [&](u64 i) {
      return (i < x.size) ? x[i].fo1() ^ ff[i].fo1() : ff[i].fo1();
    };
    folded = y.fo1().concat();
  }
};

// Rotate-XOR fold (better hash distribution)
template <u64 F> struct RotateXORFold {
  reg<F> folded;

  val<F> get() { return folded; }

  void fanout(auto fo) { folded.fanout(fo); }

  template <u64 MAXL> void update(auto &gh, u64 ghlen, auto in) {
    constexpr u64 inbits = std::min(F, u64(decltype(in)::size));
    val<inbits> input = in;
    // Rotate left by 1 then XOR in new bits and XOR out old bit
    auto f = folded.make_array(val<1>{});
    val<1> outbit = gh[ghlen - 1];
    u64 outpos = ghlen % F;
    // Rotate: shift each bit left by 1 (with wrap)
    arr<val<1>, F> rotated = [&](u64 i) {
      return (i == 0) ? f[F - 1].fo1() : f[i - 1].fo1();
    };
    // XOR out the evicted bit at its fold position
    arr<val<1>, F> after_out = [&](u64 i) {
      return (i == outpos) ? rotated[i].fo1() ^ outbit.fo1() : rotated[i].fo1();
    };
    // XOR in new bits
    auto x = input.fo1().make_array(val<1>{});
    arr<val<1>, F> y = [&](u64 i) {
      return (i < x.size) ? x[i].fo1() ^ after_out[i].fo1()
                          : after_out[i].fo1();
    };
    folded = y.fo1().concat();
  }
};

// ============================================================================
// Global History Register
// ============================================================================

// Global history stored as array of 1-bit regs (supports >64 bits)
template <u64 MAXLEN> struct TDGlobalHistory {
  arr<reg<1>, MAXLEN> hist;

  val<1> operator[](u64 pos) { return hist[pos]; }

  void update(auto in) {
    constexpr u64 inbits = decltype(in)::size;
    auto inbits_arr = in.fo1().make_array(val<1>{});
    // Shift left by inbits: move existing bits up, insert new bits at bottom
    for (u64 i = MAXLEN - 1; i >= inbits; i--)
      hist[i] = hist[i - inbits];
    for (u64 i = 0; i < inbits && i < MAXLEN; i++)
      hist[i] = inbits_arr[i];
  }

  void fanout(auto fo) {
    for (u64 i = 0; i < MAXLEN; i++)
      hist[i].fanout(fo);
  }
};

// ============================================================================
// History Folder — manages global history + per-table folds
// ============================================================================

template <u64 NUM_TABLES, u64 MAXHIST, u64 IDX_FOLD_W, u64 TAG_FOLD_W,
          template <u64> class FoldFn = XORFold>
struct TDHistoryFolder {
  TDGlobalHistory<MAXHIST> gh;

  // Per-table: one fold for index, one fold for tag
  std::array<FoldFn<IDX_FOLD_W>, NUM_TABLES> idx_folds;
  std::array<FoldFn<TAG_FOLD_W>, NUM_TABLES> tag_folds;

  // History lengths per table (set externally via init or constexpr)
  std::array<u64, NUM_TABLES> hist_len;

  void init(const std::array<u64, NUM_TABLES> &hl) { hist_len = hl; }

  val<IDX_FOLD_W> get_idx_fold(u64 i) { return idx_folds[i].get(); }
  val<TAG_FOLD_W> get_tag_fold(u64 i) { return tag_folds[i].get(); }

  void fanout(auto fo) {
    for (u64 i = 0; i < NUM_TABLES; i++) {
      idx_folds[i].fanout(fo);
      tag_folds[i].fanout(fo);
    }
  }

  void update(auto branchbits) {
    branchbits.fanout(hard<NUM_TABLES * 2 + 1>{});
    gh.fanout(hard<std::max(u64(2), NUM_TABLES * 2 + 1)>{});
    for (u64 i = 0; i < NUM_TABLES; i++) {
      idx_folds[i].template update<MAXHIST>(gh, hist_len[i], branchbits);
      tag_folds[i].template update<MAXHIST>(gh, hist_len[i], branchbits);
    }
    gh.update(branchbits);
  }
};

// ============================================================================
// TAGE Table Storage — simple struct with plain ram<> members
// ============================================================================

template <u64 TABLE_SIZE, u64 TAG_WIDTH, u64 CTR_WIDTH, u64 HYST_WIDTH,
          u64 U_WIDTH>
struct TDTable {
  static constexpr u64 IDX_BITS = td::clog2(TABLE_SIZE);
  static constexpr u64 table_size = TABLE_SIZE;
  static constexpr u64 tag_width = TAG_WIDTH;
  static constexpr u64 ctr_width = CTR_WIDTH;
  static constexpr u64 hyst_width = HYST_WIDTH;
  static constexpr u64 u_width = U_WIDTH;

  hcm::ram<val<TAG_WIDTH>, TABLE_SIZE> tag_ram{"td_tag"};
  hcm::ram<val<CTR_WIDTH>, TABLE_SIZE> pred_ram{"td_pred"};
  hcm::ram<val<std::max(u64(1), HYST_WIDTH)>, TABLE_SIZE> hyst_ram{"td_hyst"};
  hcm::ram<val<U_WIDTH>, TABLE_SIZE> u_ram{"td_u"};
};

// Generate a TDTable type from config arrays at index I
template <typename Cfg, u64 I>
using TDTableAt =
    TDTable<Cfg::TABLE_SIZE[I], Cfg::TAG_WIDTH[I], Cfg::CTR_WIDTH[I],
            Cfg::HYST_WIDTH[I], Cfg::U_WIDTH[I]>;

// Build a tuple of TDTable types
template <typename Cfg, typename Seq> struct TDMakeTableTuple;

template <typename Cfg, u64... Is>
struct TDMakeTableTuple<Cfg, std::index_sequence<Is...>> {
  using type = std::tuple<TDTableAt<Cfg, Is>...>;
};

} // namespace td

// ============================================================================
// TageDirectImpl
// ============================================================================

template <typename TableCfg, typename AllocCfg, u64 LINEINST_V, u64 N_V,
          bool SHARED_TAG_V, bool SHARED_U_V, bool SHARED_HYS_V,
          bool U_STOR_FF_V, u64 DECAY_CTR_V, u64 DECAY_GRAN_V,
          typename DecayPolicy_V, bool P1_USE_GSHARE_V, u64 P1_TABLE_SIZE_V,
          u64 P1_HIST_V, bool USE_META_V, u64 METABITS_V, u64 METAPIPE_V,
          bool USE_PATH_HIST_V, u64 PATH_HIST_WIDTH_V, u64 PATH_BITS_V,
          template <u64> class FoldFn_V>
struct TageDirectImpl : predictor {

  // ======== Constants ========
  static constexpr u64 NUM_TABLES = TableCfg::NUM_TABLES;
  static constexpr u64 N = N_V;
  static constexpr u64 LANES = td::next_pow2(N);
  static constexpr u64 LOG_LANES = td::clog2(LANES);
  static constexpr u64 LINEINST = LINEINST_V;
  static constexpr u64 LOG_LINEINST = td::clog2(LINEINST);

  static constexpr u64 MAX_TAG_WIDTH = td::array_max(TableCfg::TAG_WIDTH);
  static constexpr u64 MAX_CTR_WIDTH = td::array_max(TableCfg::CTR_WIDTH);
  static constexpr u64 MAX_HYST_WIDTH = td::array_max(TableCfg::HYST_WIDTH);
  static constexpr u64 MAX_U_WIDTH = td::array_max(TableCfg::U_WIDTH);
  static constexpr u64 MAX_TABLE_SIZE = td::array_max(TableCfg::TABLE_SIZE);
  static constexpr u64 MAX_IDX_BITS = td::clog2(MAX_TABLE_SIZE);
  static constexpr u64 MAXHIST = td::array_max(TableCfg::HIST_LEN);
  static constexpr u64 MAX_HTAGBITS = MAX_TAG_WIDTH - LOG_LANES;

  static constexpr u64 P1_ENTRIES = P1_TABLE_SIZE_V;
  static constexpr u64 P1_INDEX_BITS = td::clog2(P1_TABLE_SIZE_V);
  static constexpr u64 LINEADDR_BITS = MAX_IDX_BITS;
  static constexpr u64 UCTRBITS = 8;
  static constexpr u64 PATHBITS = PATH_BITS_V;

  static constexpr bool USE_PROB_DECAY = (DECAY_CTR_V > 0);

  // ---- Table tuple type ----
  using Tables =
      typename td::TDMakeTableTuple<TableCfg,
                                    std::make_index_sequence<NUM_TABLES>>::type;

  // Truncate gindex to per-table IDX_BITS
  template <u64 I> auto tidx(auto &gi) {
    using Table = std::tuple_element_t<I, Tables>;
    return val<Table::IDX_BITS>{gi};
  }

  // ======== Static asserts ========
  static_assert(NUM_TABLES > 0);
  static_assert(N >= 1);
  static_assert(LANES >= N);
  static_assert(std::has_single_bit(LINEINST));
  static_assert(MAX_TAG_WIDTH > LOG_LANES, "TAG_WIDTH must be > LOG_LANES");

  // ======== State ========

  // History
  td::TDHistoryFolder<NUM_TABLES, MAXHIST, MAX_IDX_BITS, MAX_HTAGBITS, FoldFn_V>
      gfolds;
  bool gfolds_inited = false;
  reg<1> true_block = 1;

  // P1 gshare (gshareN pattern)
  std::conditional_t<P1_USE_GSHARE_V, reg<P1_HIST_V>, EmptyMember>
      global_history1;
  reg<P1_INDEX_BITS> index1;
  reg<LANES> X; // lane scrambling
  reg<LANES> unordered_pred1;
  arr<reg<1>, LANES> pred; // extracted per-lane P1 predictions

  // P2 TAGE
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> gindex;
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> htag;

  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> readt;
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> readc;
  arr<reg<std::max(u64(1), MAX_HYST_WIDTH)>, NUM_TABLES> readh;
  arr<reg<MAX_U_WIDTH>, NUM_TABLES> readu;
  reg<NUM_TABLES> notumask;

  arr<reg<1>, NUM_TABLES> htagcmp_reg;
  arr<reg<1>, NUM_TABLES> last_tagcmp_reg;

  // Per-rank match/prediction
  arr<reg<NUM_TABLES + 1>, LANES> match;
  arr<reg<NUM_TABLES + 1>, LANES> match1;
  arr<reg<NUM_TABLES + 1>, LANES> match2;
  arr<reg<1>, LANES> pred1_tage;
  arr<reg<1>, LANES> pred2_tage;
  reg<LANES> p2;

  // Meta
  std::conditional_t<USE_META_V, arr<reg<METABITS_V, i64>, METAPIPE_V>,
                     EmptyMember>
      meta;
  std::conditional_t<USE_META_V, arr<reg<1>, LANES>, EmptyMember> newly_alloc;

  // U-bit reset
  reg<UCTRBITS> uctr;
  std::conditional_t<USE_PROB_DECAY, reg<DECAY_CTR_V == 0 ? 1 : DECAY_CTR_V>,
                     EmptyMember>
      decay_threshold;

  // Path history
  std::conditional_t<USE_PATH_HIST_V, reg<PATH_HIST_WIDTH_V>, EmptyMember>
      path_hist;

  // Block tracking (gshareN_ahead_best pattern: numeric offset)
  reg<LOG_LINEINST> block_entry;
  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<1>, LANES> branch_dir;
  reg<N + 1> rank;  // one-hot: rank of current branch in block (gshareN_ahead pattern)

  // P1 storage
  hcm::ram<val<LANES>, P1_ENTRIES> p1_pred{"P1 pred"};
  zone UPDATE_ONLY;
  hcm::ram<val<1>, P1_ENTRIES> p1_hyst[LANES]{"P1 hyst"};

  // TAGE tables
  Tables tables;

  bool params_printed = false;

  // ======== Helpers ========

  val<1> line_end() { return (block_entry + block_size) >= hard<LINEINST>{}; }

  val<1> last_pred() {
    assert(num_branch <= N);
    return rank >> (N - num_branch);
  }

  void ensure_gfolds_init() {
    if (!gfolds_inited) {
      std::array<u64, NUM_TABLES> hl;
      for (u64 i = 0; i < NUM_TABLES; i++)
        hl[i] = TableCfg::HIST_LEN[i];
      gfolds.init(hl);
      gfolds_inited = true;
    }
  }

  void print_params(std::ostream &os) const {
    os << "=== TageDirect Parameters ===\n";
    os << "NUM_TABLES=" << NUM_TABLES << "  LINEINST=" << LINEINST
       << "  N=" << N << "  LANES=" << LANES << "\n";
    os << "Table sizes: ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::TABLE_SIZE[i];
    os << "\nTag widths:  ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::TAG_WIDTH[i];
    os << "\nHist lens:   ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::HIST_LEN[i];
    os << "\nP1: TABLE_SIZE=" << P1_TABLE_SIZE_V << "  HIST=" << P1_HIST_V;
    os << "\nMETA: USE=" << USE_META_V << "  BITS=" << METABITS_V
       << "  PIPE=" << METAPIPE_V;
    os << "\nMISPREDICT_ONLY_WRITE=" << AllocCfg::MISPREDICT_ONLY_WRITE;
    os << "\n\n";
  }

  // ======== Predictor Interface ========

  val<1> predict1(val<64> inst_pc) override {
    ensure_gfolds_init();
    if (!params_printed) {
      print_params(std::cerr);
      params_printed = true;
    }

    inst_pc.fanout(hard<3>{});
    true_block.fanout(hard<4>{});

    // Block entry: numeric offset (gshareN_ahead_best pattern)
    block_entry = select(true_block, val<LOG_LINEINST>{inst_pc >> 2},
                         val<LOG_LINEINST>{block_entry + block_size});
    block_entry.fanout(hard<LANES + 2>{});

    // Lane scrambling
    rank = select(true_block, val<N + 1>{1}, rank << num_branch);
    rank.fanout(hard<N + 2>{});
    X = select(true_block, val<LOG_LANES>{inst_pc >> 2}.decode().concat(),
               X.rotate_left(num_branch));
    X.fanout(hard<LANES>{});

    // P1 gshare: read on true blocks only
    if constexpr (P1_USE_GSHARE_V) {
      global_history1.fanout(hard<2>{});
    }
    execute_if(true_block, [&]() {
      val<P1_INDEX_BITS> lineaddr = inst_pc >> 2;
      if constexpr (P1_USE_GSHARE_V) {
        if constexpr (P1_HIST_V <= P1_INDEX_BITS) {
          index1 = lineaddr ^ (val<P1_INDEX_BITS>{global_history1}
                               << (P1_INDEX_BITS - P1_HIST_V));
        } else {
          index1 = global_history1.make_array(val<P1_INDEX_BITS>{})
                       .append(lineaddr)
                       .fold_xor();
        }
      } else {
        index1 = lineaddr;
      }
      std::cerr << "P1: gshare p1_pred READ\n";
      unordered_pred1 = p1_pred.read(index1);
    });

    // Extract per-lane predictions
    unordered_pred1.fanout(hard<LANES>{});
    for (u64 i = 0; i < LANES; i++) {
      pred[i] = (unordered_pred1 & X.rotate_left(i)) != hard<0>{};
    }
    pred.fanout(hard<LANES * 2>{});

    block_size = 1;
    num_branch = 0;
    reuse_prediction(~line_end());
    return pred[0];
  }

  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) override {
    block_size++;
    reuse_prediction(~line_end());
    return pred[num_branch];
  }

  val<1> predict2([[maybe_unused]] val<64> inst_pc) override {
    // TAGE disabled — return not-taken
    std::cerr << "P2: STUB (not-taken)\n";
    p2 = hard<0>{};
    p2.fanout(hard<LANES>{});
    reuse_prediction(~val<1>{block_entry >> (LINEINST - 1)});
    return hard<0>{};
#if 0  // --- TAGE P2 disabled ---
    val<LINEADDR_BITS> lineaddr = inst_pc >> 2;
    lineaddr.fanout(hard<1 + NUM_TABLES * 2>{});
    gfolds.fanout(hard<2>{});

    // Compute TAGE indices
    for (u64 i = 0; i < NUM_TABLES; i++) {
      if constexpr (USE_PATH_HIST_V) {
        gindex[i] =
            lineaddr ^ gfolds.get_idx_fold(i) ^ val<MAX_IDX_BITS>{path_hist};
      } else {
        gindex[i] = lineaddr ^ gfolds.get_idx_fold(i);
      }
    }
    gindex.fanout(hard<4>{});

    // Read TAGE tables
    static_loop<NUM_TABLES>([&]<u64 I>() {
      auto &t = std::get<I>(tables);
      std::cerr << "P2: tage tag_ram[" << I << "] READ\n";
      readt[I] = t.tag_ram.read(tidx<I>(gindex[I]));
      std::cerr << "P2: tage pred_ram[" << I << "] READ\n";
      readc[I] = t.pred_ram.read(tidx<I>(gindex[I]));
      if constexpr (MAX_HYST_WIDTH > 0) {
        std::cerr << "P2: tage hyst_ram[" << I << "] READ\n";
        readh[I] = t.hyst_ram.read(tidx<I>(gindex[I]));
      }
      std::cerr << "P2: tage u_ram[" << I << "] READ\n";
      readu[I] = t.u_ram.read(tidx<I>(gindex[I]));
    });
    readt.fanout(hard<LANES + 1>{});
    readc.fanout(hard<3>{});
    if constexpr (MAX_HYST_WIDTH > 0)
      readh.fanout(hard<2>{});
    readu.fanout(hard<2>{});
    notumask = ~readu.concat();
    notumask.fanout(hard<2>{});

    // Compute hashed tags (parallel with RAM reads)
    for (u64 i = 0; i < NUM_TABLES; i++) {
      htag[i] = val<MAX_HTAGBITS>{lineaddr}.reverse() ^ gfolds.get_tag_fold(i);
    }
    htag.fanout(hard<2>{});

    // Gather prediction bits per table
    val<NUM_TABLES> gpreds = [&]() -> val<NUM_TABLES> {
      if constexpr (MAX_CTR_WIDTH == 1) {
        return readc.concat();
      } else {
        arr<val<1>, NUM_TABLES> gp = [&](int i) -> val<1> {
          return readc[i] >> hard<MAX_CTR_WIDTH - 1>{};
        };
        return gp.fo1().concat();
      }
    }();
    gpreds.fanout(hard<LANES>{});

    // Per-rank preds: P1 gshare as fallback (position NUM_TABLES = bimodal
    // slot)
    pred.fanout(hard<LANES>{});
    arr<val<NUM_TABLES + 1>, LANES> preds = [&](u64 r) {
      return concat(pred[r], gpreds);
    };
    preds.fanout(hard<2 * LANES>{});

    // Per-table htag comparison
    static_loop<NUM_TABLES>([&]<u64 I>() {
      using Table = std::tuple_element_t<I, Tables>;
      static constexpr u64 PER_HTAG = Table::tag_width - LOG_LANES;
      htagcmp_reg[I] = (val<PER_HTAG>{readt[I]} == val<PER_HTAG>{htag[I]});
    });
    htagcmp_reg.fanout(hard<LANES + 1>{});

    // Per-rank tag match
    static_loop<LANES>([&]<u64 R>() {
      arr<val<1>, NUM_TABLES> tagcmp = [&](int i) {
        return val<LOG_LANES>{readt[i] >> MAX_HTAGBITS} == hard<R>{};
      };
      match[R] =
          concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp_reg.concat());
    });
    match.fanout(hard<2>{});

    // Provider/alt selection per rank
    for (u64 r = 0; r < LANES; r++) {
      match1[r] = match[r].one_hot();
    }
    match1.fanout(hard<3>{});
    for (u64 r = 0; r < LANES; r++) {
      pred1_tage[r] = (match1[r] & preds[r]) != hard<0>{};
    }
    pred1_tage.fanout(hard<2>{});

    for (u64 r = 0; r < LANES; r++) {
      match2[r] = (match[r] ^ match1[r]).one_hot();
    }
    match2.fanout(hard<2>{});
    for (u64 r = 0; r < LANES; r++) {
      pred2_tage[r] = (match2[r] & preds[r]) != hard<0>{};
    }
    pred2_tage.fanout(hard<2>{});

    // Meta prediction
    if constexpr (USE_META_V) {
      meta.fanout(hard<2>{});
      arr<val<1>, NUM_TABLES> weakctr = [&](int i) {
        return readh[i] == hard<0>{};
      };
      val<NUM_TABLES> coldctr = notumask & weakctr.fo1().concat();
      coldctr.fanout(hard<LANES>{});
      val<1> metasign = (meta[METAPIPE_V - 1] >= hard<0>{});
      metasign.fanout(hard<LANES>{});
      for (u64 r = 0; r < LANES; r++) {
        newly_alloc[r] = (match1[r] & coldctr) != hard<0>{};
      }
      newly_alloc.fanout(hard<2>{});
      arr<val<1>, LANES> altsel = [&](u64 r) {
        arr<val<1>, 3> inputs = {metasign, newly_alloc[r],
                                 match2[r] != hard<0>{}};
        return inputs.fo1().fold_and();
      };
      p2 = arr<val<1>, LANES>{[&](u64 r) {
             return select(altsel[r].fo1(), pred2_tage[r], pred1_tage[r]);
           }}.concat();
    } else {
      p2 = pred1_tage.concat();
    }

    p2.fanout(hard<LANES>{});
    val<1> taken = p2 >> num_branch;
    taken.fanout(hard<2>{});
    reuse_prediction(~val<1>{block_entry >> (LINEINST - 1)});
    return taken;
#endif // --- TAGE P2 disabled ---
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) override {
    block_size++;
    reuse_prediction(~line_end());
    return hard<0>{};
  }

  void update_condbr([[maybe_unused]] val<64> branch_pc, val<1> taken,
                     [[maybe_unused]] val<64> next_pc) override {
    assert(num_branch < N);
    branch_dir[num_branch] = taken.fo1();
    num_branch++;
    reuse_prediction(~(line_end() | last_pred()));
  }

  void update_cycle(instruction_info &block_end_info) override {
    val<1> &mispredict = block_end_info.is_mispredict;
    val<64> &next_pc = block_end_info.next_pc;

    std::cerr << "UC: ENTER (num_branch=" << num_branch
              << " misp=" << static_cast<u64>(block_end_info.is_mispredict)
              << ")\n";
    if (num_branch == 0) {
      // No conditional branches — just update history
      std::cerr << "UC: EXIT (no branches)\n";
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
        true_block = 1;
      });
      return;
    }

    // ---- Gshare-only update (TAGE disabled) ----
    mispredict.fanout(hard<2>{});
    val<1> correct_pred = ~mispredict;
    correct_pred.fanout(hard<2>{});
    index1.fanout(hard<LANES * 3>{});
    branch_dir.fanout(hard<2>{});
    gfolds.fanout(hard<2>{});
    X.fanout(hard<LANES + 1>{});

    // ---- P1 gshare update (read / need_extra / write pattern) ----
    // Same structure as gshareN reference: reads in cycle 1, writes in cycle 2.

    // Combinational: compute which lanes were accessed by branches in this block.
    // access[i] = OR of all X.rotate_left(branch_rank) masks — 1 if lane i was used.
    arr<val<1>, LANES> access =
        arr<val<LANES>, LANES>{[&](u64 i) -> val<LANES> {
          return X.rotate_left(i) & val<LANES>{-(i < num_branch)};
        }}.fold_or()
            .make_array(val<1>{});

    // Combinational: identify the lane that holds the mispredicted branch.
    // misp_bank is a one-hot mask ANDed with mispredict — all zero on correct prediction.
    val<LANES> misp_bank = X.rotate_left(num_branch - 1)
        & mispredict.replicate(hard<LANES>{}).concat();
    arr<val<1>, LANES> mispredicted = misp_bank.fo1().make_array(val<1>{});
    mispredicted.fanout(hard<2>{});

    // Cycle 1 RAM read: read hysteresis for the mispredicted lane.
    // Gated by mispredicted[i] so only one lane's RAM is actually read.
    // On correct prediction, mispredicted is all-zero so no RAM is accessed.
    arr<val<1>, LANES> weak = [&](u64 i) {
      return execute_if(mispredicted[i], [&]() {
#ifdef CHEATING_MODE
        if (static_cast<u64>(mispredict) && static_cast<u64>(mispredicted[i]))
          std::cerr << "UC: gshare p1_hyst[" << i << "] READ\n";
#endif
        return p1_hyst[i].read(index1);
      });
    };

    // Cycle boundary: grant an extra cycle when mispredict=1.
    // Everything above this line is cycle 1 (reads).
    // Everything below this line is cycle 2 (writes).
    need_extra_cycle(mispredict);

    // Cycle 2 RAM write: flip the prediction bit if hysteresis was weak.
    // p1_pred was read in predict1 (cycle 1), so this write is safe in cycle 2.
    execute_if(mispredict, [&]() {
      arr<val<1>, LANES> stored = unordered_pred1.make_array(val<1>{});
      arr<val<1>, LANES> bundle = [&](u64 i) {
        return select(weak[i].fo1(), branch_dir[num_branch - 1],
                      stored[i].fo1());
      };
#ifdef CHEATING_MODE
      if (static_cast<u64>(mispredict))
        std::cerr << "UC: gshare p1_pred WRITE\n";
#endif
      p1_pred.write(index1, bundle.fo1().concat());
    });

    // Cycle 2 RAM write: update hysteresis for all accessed lanes.
    // p1_hyst[i] was read above in cycle 1, so this write is safe in cycle 2.
    for (u64 i = 0; i < LANES; i++) {
      execute_if(access[i].fo1(), [&]() {
#ifdef CHEATING_MODE
        if (static_cast<u64>(mispredict) && static_cast<u64>(access[i]))
          std::cerr << "UC: gshare p1_hyst[" << i << "] WRITE\n";
#endif
        p1_hyst[i].write(index1, mispredicted[i].fo1());
      });
    }

#if 0 // --- TAGE update disabled ---
    correct_pred.fanout(hard<NUM_TABLES + 2>{});
    index1.fanout(hard<LANES * 3>{});
    gindex.fanout(hard<4>{});
    htag.fanout(hard<3>{});
    readt.fanout(hard<4>{});
    readc.fanout(hard<2>{});
    if constexpr (MAX_HYST_WIDTH > 0)
      readh.fanout(hard<3>{});
    match1.fanout(hard<3>{});
    match2.fanout(hard<2>{});
    pred1_tage.fanout(hard<2>{});
    pred2_tage.fanout(hard<2 + NUM_TABLES>{});
    branch_dir.fanout(hard<2>{});
    gfolds.fanout(hard<2>{});
    readu.fanout(hard<2>{});
    X.fanout(hard<LANES + 1>{});
    if constexpr (USE_META_V)
      meta.fanout(hard<2>{});

    val<LOG_LANES> last_rank = val<LOG_LANES>{num_branch - 1};
    last_rank.fanout(hard<4 * NUM_TABLES + 2>{});

    // Determine which lanes have branches
    arr<val<1>, LANES> is_branch = [&](u64 r) -> val<1> {
      return val<1>{r < num_branch ? 1u : 0u};
    };
    is_branch.fanout(hard<4>{});

    arr<val<1>, LANES> branch_taken = [&](u64 r) -> val<1> {
      if (r < num_branch)
        return branch_dir[r];
      return val<1>{0};
    };
    branch_taken.fanout(hard<3>{});

    // Per-rank match restricted to actual branches
    arr<val<NUM_TABLES + 1>, LANES> actual_match1 = [&](u64 r) {
      return select(is_branch[r], match1[r], val<NUM_TABLES + 1>{0});
    };
    actual_match1.fanout(hard<2>{});

    val<NUM_TABLES> primary_mask = actual_match1.fold_or();
    primary_mask.fanout(hard<2>{});
    arr<val<1>, NUM_TABLES> primary = primary_mask.make_array(val<1>{});
    primary.fanout(hard<3>{});

    arr<val<1>, LANES> primary_wrong = [&](u64 r) {
      return pred1_tage[r] != branch_taken[r];
    };
    primary_wrong.fanout(hard<2>{});

    // ---- Allocation ----
    val<NUM_TABLES> mispmask =
        mispredict.replicate(hard<NUM_TABLES>{}).concat();

    // Per-table tag comparison for last branch (uses last_tagcmp_reg)
    static_loop<NUM_TABLES>([&]<u64 I>() {
      using Table = std::tuple_element_t<I, Tables>;
      static constexpr u64 PER_HTAG = Table::tag_width - LOG_LANES;
      last_tagcmp_reg[I] = (val<LOG_LANES>{readt[I] >> PER_HTAG} == last_rank) &
                           (val<PER_HTAG>{readt[I]} == val<PER_HTAG>{htag[I]});
    });
    val<NUM_TABLES + 1> last_match1 =
        last_tagcmp_reg.fo1().append(1).concat().one_hot();
    last_match1.fanout(hard<2>{});

    // Postmask: tables above the provider
    val<NUM_TABLES> postmask = [&]() -> val<NUM_TABLES> {
      if constexpr (AllocCfg::PROB_START > 0) {
        val<2> rstart = val<2>{static_cast<u64>(std::rand())};
        val<NUM_TABLES> base = mispmask & val<NUM_TABLES>(last_match1 - 1);
        val<NUM_TABLES> skip1 = base & val<NUM_TABLES>(base - 1);
        val<NUM_TABLES> skip2 = skip1 & val<NUM_TABLES>(skip1 - 1);
        return select(rstart == hard<0>{}, skip2,
                      select(rstart == hard<1>{}, skip1, base));
      } else {
        return mispmask.fo1() & val<NUM_TABLES>(last_match1 - 1);
      }
    }();
    postmask.fanout(hard<2>{});

    // Candidate allocation mask
    val<NUM_TABLES> candallocmask = [&]() -> val<NUM_TABLES> {
      if constexpr (AllocCfg::CONF_GATE) {
        arr<val<1>, NUM_TABLES> weak_entry = [&](u64 i) -> val<1> {
          if constexpr (MAX_CTR_WIDTH == 1)
            return val<1>{1};
          else {
            auto ctr = readc[i];
            return (ctr == hard<0>{}) |
                   (ctr == hard<(1u << (MAX_CTR_WIDTH - 1)) - 1>{}) |
                   (ctr == hard<(1u << (MAX_CTR_WIDTH - 1))>{});
          }
        };
        return postmask & notumask & weak_entry.fo1().concat();
      } else {
        return postmask & notumask;
      }
    }();
    candallocmask.fanout(hard<2>{});
    val<NUM_TABLES> collamask = candallocmask.reverse();
    collamask.fanout(hard<2>{});
    val<NUM_TABLES> collamask1 = collamask.one_hot();
    collamask1.fanout(hard<3>{});

    arr<val<1>, NUM_TABLES> allocate = [&]() -> arr<val<1>, NUM_TABLES> {
      if constexpr (AllocCfg::MAX_ALLOC >= 2) {
        val<NUM_TABLES> pick2 = [&]() -> val<NUM_TABLES> {
          val<NUM_TABLES> basic2 = (collamask ^ collamask1).one_hot();
          if constexpr (AllocCfg::NON_CONSECUTIVE) {
            val<NUM_TABLES> neighbors = (collamask1 << 1) | (collamask1 >> 1);
            val<NUM_TABLES> nc_mask = (collamask ^ collamask1) & ~neighbors;
            val<NUM_TABLES> nc_pick = nc_mask.reverse().one_hot();
            return select(nc_mask != hard<0>{}, nc_pick, basic2);
          } else {
            return basic2;
          }
        }();
        return (collamask1 | pick2).reverse().make_array(val<1>{});
      } else {
        val<NUM_TABLES> collamask2 = (collamask ^ collamask1).one_hot();
        val<NUM_TABLES> collamask12 = select(val<2>{std::rand()} == hard<0>{},
                                             collamask2.fo1(), collamask1);
        return collamask12.fo1().reverse().make_array(val<1>{});
      }
    }();
    allocate.fanout(hard<7>{});

    // ---- Branch direction per table ----
    arr<val<1>, NUM_TABLES> bdir = [&](u64 i) {
      val<LOG_LANES> stored_rank = readt[i] >> MAX_HTAGBITS;
      val<LOG_LANES> use_rank =
          select(allocate[i], last_rank, stored_rank.fo1());
      // Select branch_dir at the stored rank
      return branch_dir.select(use_rank);
    };
    bdir.fanout(hard<2>{});

    arr<val<1>, NUM_TABLES> badpred1 = [&](u64 i) -> val<1> {
      if constexpr (MAX_CTR_WIDTH == 1)
        return readc[i] != bdir[i];
      else
        return val<1>{readc[i] >> hard<MAX_CTR_WIDTH - 1>{}} != bdir[i];
    };
    badpred1.fanout(hard<3>{});

    arr<val<1>, NUM_TABLES> altdiffer = [&](u64 i) -> val<1> {
      auto pred_dir = [&]() -> val<1> {
        if constexpr (MAX_CTR_WIDTH == 1)
          return readc[i];
        else
          return readc[i] >> hard<MAX_CTR_WIDTH - 1>{};
      }();
      val<LOG_LANES> stored_rank = readt[i] >> MAX_HTAGBITS;
      return pred_dir != pred2_tage.select(stored_rank.fo1());
    };

    arr<val<1>, NUM_TABLES> goodpred = [&](u64 i) {
      val<LOG_LANES> stored_rank = readt[i] >> MAX_HTAGBITS;
      return (stored_rank.fo1() != last_rank) | correct_pred;
    };
    goodpred.fanout(hard<2>{});

    // ---- Hysteresis weakness ----
    arr<val<1>, NUM_TABLES> g_weak = [&](u64 i) -> val<1> {
      if constexpr (MAX_HYST_WIDTH > 0)
        return primary[i] & badpred1[i] & (readh[i] == hard<0>{});
      else
        return primary[i] & badpred1[i];
    };
    g_weak.fanout(hard<2>{});

    // ---- P1 disagreement ----
    val<LANES> p1_concat = pred.concat();
    val<LANES> disagree_mask = (p1_concat ^ p2) & is_branch.concat();
    disagree_mask.fanout(hard<2>{});

    // ---- Extra cycle ----
    if constexpr (AllocCfg::MISPREDICT_ONLY_WRITE) {
      // Only need extra cycle on misprediction
      need_extra_cycle(mispredict);
    } else {
      val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
      val<1> extra_cycle =
          some_badpred1.fo1() | mispredict | (disagree_mask != hard<0>{});
      extra_cycle.fanout(hard<NUM_TABLES * 2 + 1>{});
      need_extra_cycle(extra_cycle);
    }

    // ---- TAGE tag write (allocation) ----
    static_loop<NUM_TABLES>([&]<u64 I>() {
      execute_if(allocate[I], [&]() {
        std::cerr << "UC: tage tag_ram[" << I << "] WRITE (alloc)\n";
        std::get<I>(tables).tag_ram.write(tidx<I>(gindex[I]),
                                          concat(last_rank, htag[I]));
      });
    });

    // ---- U-bit update ----
    arr<val<1>, NUM_TABLES> update_u = [&](u64 i) {
      return primary[i] & altdiffer[i].fo1();
    };
    val<1> noalloc = (candallocmask == hard<0>{});
    val<NUM_TABLES> uclearmask =
        postmask & noalloc.fo1().replicate(hard<NUM_TABLES>{}).concat();
    arr<val<1>, NUM_TABLES> uclear = uclearmask.fo1().make_array(val<1>{});
    uclear.fanout(hard<2>{});

    if constexpr (AllocCfg::MISPREDICT_ONLY_WRITE) {
      // Only write u-bits on misprediction
      static_loop<NUM_TABLES>([&]<u64 I>() {
        execute_if(mispredict & (allocate[I] | uclear[I]), [&]() {
          std::cerr << "UC: tage u_ram[" << I << "] WRITE (misp_only)\n";
          val<1> newu =
              select(allocate[I], val<1>{1}, val<1>{0}); // init u=1 on alloc
          std::get<I>(tables).u_ram.write(tidx<I>(gindex[I]), newu.fo1());
        });
      });
    } else if constexpr (USE_PROB_DECAY) {
      val<DECAY_CTR_V> lfsr = val<DECAY_CTR_V>{static_cast<u64>(std::rand())};
      val<1> decay_fire = (lfsr > val<DECAY_CTR_V>{decay_threshold});
      decay_fire.fanout(hard<NUM_TABLES>{});
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> newu =
            goodpred[I].fo1() & ~allocate[I] & ~uclear[I] & ~decay_fire;
        val<1> u_changed = (val<1>{readu[I]} != newu);
        execute_if((update_u[I].fo1() | allocate[I] | uclear[I] | decay_fire) &
                       u_changed,
                   [&]() {
                     std::cerr << "UC: tage u_ram[" << I << "] WRITE (decay)\n";
                     std::get<I>(tables).u_ram.write(tidx<I>(gindex[I]),
                                                     newu.fo1());
                   });
      });
    } else {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> newu = goodpred[I].fo1() & ~allocate[I] & ~uclear[I];
        val<1> u_changed = (val<1>{readu[I]} != newu);
        execute_if(
            (update_u[I].fo1() | allocate[I] | uclear[I]) & u_changed, [&]() {
              std::cerr << "UC: tage u_ram[" << I << "] WRITE (default)\n";
              std::get<I>(tables).u_ram.write(tidx<I>(gindex[I]), newu.fo1());
            });
      });
    }

    // ---- P1 gshare update (gshareN pattern) ----
    // All P1 RAM accesses gated by mispredict to ensure they happen in extra
    // cycle (p1_pred was read in predict1, so writes must be in a different
    // cycle)
    execute_if(mispredict, [&]() {
      arr<val<1>, LANES> access =
          arr<val<LANES>, LANES>{[&](u64 i) -> val<LANES> {
            return X.rotate_left(i) & val<LANES>{-(i < num_branch)};
          }}.fold_or()
              .make_array(val<1>{});

      val<LANES> misp_bank = X.rotate_left(num_branch - 1);
      arr<val<1>, LANES> mispredicted = misp_bank.fo1().make_array(val<1>{});
      mispredicted.fanout(hard<2>{});

      // Read hysteresis for mispredicted lane
      arr<val<1>, LANES> weak = [&](u64 i) {
        return execute_if(mispredicted[i], [&]() {
          std::cerr << "UC: gshare p1_hyst[" << i << "] READ\n";
          return p1_hyst[i].read(index1);
        });
      };

      // Flip prediction if weak
      arr<val<1>, LANES> stored = unordered_pred1.make_array(val<1>{});
      arr<val<1>, LANES> bundle = [&](u64 i) {
        return select(weak[i].fo1(), branch_dir[num_branch - 1],
                      stored[i].fo1());
      };
      std::cerr << "UC: gshare p1_pred WRITE\n";
      p1_pred.write(index1, bundle.fo1().concat());

      // Update hysteresis for all accessed lanes
      // TODO: Only do this on mispredict
      for (u64 i = 0; i < LANES; i++) {
        execute_if(access[i].fo1(), [&]() {
          std::cerr << "UC: gshare p1_hyst[" << i << "] WRITE\n";
          p1_hyst[i].write(index1, mispredicted[i].fo1());
        });
      }
    });

    // ---- TAGE counter update ----
    if constexpr (AllocCfg::MISPREDICT_ONLY_WRITE) {
      // Only write counters on misprediction; initialize saturated on alloc
      static_loop<NUM_TABLES>([&]<u64 I>() {
        execute_if(mispredict & (g_weak[I].fo1() | allocate[I]), [&]() {
          std::cerr << "UC: tage pred_ram[" << I << "] WRITE (misp_only)\n";
          if constexpr (MAX_CTR_WIDTH == 1) {
            std::get<I>(tables).pred_ram.write(tidx<I>(gindex[I]), bdir[I]);
          } else {
            // Saturate toward branch direction on allocation
            auto init_ctr =
                select(bdir[I], val<MAX_CTR_WIDTH>{(1u << MAX_CTR_WIDTH) - 1},
                       val<MAX_CTR_WIDTH>{0});
            std::get<I>(tables).pred_ram.write(
                tidx<I>(gindex[I]),
                select(allocate[I], init_ctr, val<MAX_CTR_WIDTH>{bdir[I]}));
          }
        });
      });
    } else {
      // Standard: update on weak+wrong or allocate, with silent elimination
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> old_dir = [&]() -> val<1> {
          if constexpr (MAX_CTR_WIDTH == 1)
            return readc[I];
          else
            return readc[I] >> hard<MAX_CTR_WIDTH - 1>{};
        }();
        val<1> pred_changed = (old_dir != bdir[I]) | allocate[I];
        execute_if((g_weak[I].fo1() | allocate[I]) & pred_changed, [&]() {
          std::cerr << "UC: tage pred_ram[" << I << "] WRITE (standard)\n";
          std::get<I>(tables).pred_ram.write(tidx<I>(gindex[I]), bdir[I]);
        });
      });
    }

    // ---- Hysteresis update ----
    if constexpr (MAX_HYST_WIDTH > 0 && !AllocCfg::MISPREDICT_ONLY_WRITE) {
      static constexpr u64 HW = std::max(u64(1), MAX_HYST_WIDTH);
      static constexpr u64 HMAX = (u64(1) << HW) - 1;
      if constexpr (AllocCfg::PARTIAL_UPDATE)
        altdiffer.fanout(hard<2>{});
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> would_change = allocate[I] |
                              (badpred1[I] & (readh[I] != hard<0>{})) |
                              (~badpred1[I] & (readh[I] != hard<HMAX>{}));
        val<1> should_update = [&]() -> val<1> {
          if constexpr (AllocCfg::PARTIAL_UPDATE)
            return allocate[I] | badpred1[I] | altdiffer[I].fo1();
          else
            return primary[I] | allocate[I];
        }();
        execute_if(should_update & would_change, [&]() {
          std::cerr << "UC: tage hyst_ram[" << I << "] WRITE (standard)\n";
          auto newhyst = select(allocate[I], val<HW>{0},
                                td::update_ctr(readh[I], ~badpred1[I]));
          std::get<I>(tables).hyst_ram.write(tidx<I>(gindex[I]), newhyst.fo1());
        });
      });
    } else if constexpr (MAX_HYST_WIDTH > 0 &&
                         AllocCfg::MISPREDICT_ONLY_WRITE) {
      // Initialize hysteresis saturated on allocation
      static constexpr u64 HW = std::max(u64(1), MAX_HYST_WIDTH);
      static constexpr u64 HMAX = (u64(1) << HW) - 1;
      static_loop<NUM_TABLES>([&]<u64 I>() {
        execute_if(mispredict & allocate[I], [&]() {
          std::cerr << "UC: tage hyst_ram[" << I << "] WRITE (misp_only)\n";
          std::get<I>(tables).hyst_ram.write(tidx<I>(gindex[I]), val<HW>{HMAX});
        });
      });
    }

    // ---- Meta counter update ----
    if constexpr (USE_META_V) {
      arr<val<1>, LANES> altdiff = [&](u64 r) {
        return (match2[r] != hard<0>{}) & (pred2_tage[r] != pred1_tage[r]);
      };
      arr<val<2, i64>, LANES> meta_incr = [&](u64 r) -> val<2, i64> {
        val<1> update_meta = is_branch[r] & altdiff[r].fo1() & newly_alloc[r];
        val<1> bad_pred2 = (pred2_tage[r] != branch_taken[r]);
        return select(update_meta.fo1(), concat(bad_pred2.fo1(), val<1>{1}),
                      val<2>{0});
      };
      for (u64 i = METAPIPE_V - 1; i != 0; i--)
        meta[i] = meta[i - 1];
      auto newmeta = meta[0] + meta_incr.fo1().fold_add();
      newmeta.fanout(hard<3>{});
      using meta_t = valt<decltype(meta[0])>;
      meta[0] = select(newmeta > meta_t::maxval, meta_t{meta_t::maxval},
                       select(newmeta < meta_t::minval, meta_t{meta_t::minval},
                              meta_t{newmeta}));
    }

    // ---- Epoch / decay threshold ----
    uctr.fanout(hard<3>{});
    val<NUM_TABLES> allocmask1 = collamask1.reverse();
    allocmask1.fanout(hard<2>{});
    val<1> faralloc =
        (((last_match1 >> 3) | allocmask1).one_hot() ^ allocmask1) == hard<0>{};
    val<1> uctrsat = (uctr == hard<decltype(uctr)::maxval>{});
    uctrsat.fanout(hard<2>{});
    uctr = select(correct_pred, uctr,
                  select(uctrsat, val<decltype(uctr)::size>{0},
                         td::update_ctr(uctr, faralloc.fo1())));

    if constexpr (USE_PROB_DECAY) {
      val<1> threshold_tick = [&]() -> val<1> {
        if constexpr (DECAY_GRAN_V == 0)
          return ~correct_pred;
        else
          return (uctr & hard<(u64(1) << DECAY_GRAN_V) - 1>{}) == hard<0>{};
      }();
      val<1> misp = ~correct_pred;
      decay_threshold =
          select(threshold_tick,
                 DecayPolicy_V::template apply<DECAY_CTR_V>(
                     decay_threshold, correct_pred, uctrsat, misp),
                 val<DECAY_CTR_V>{decay_threshold});
    } else {
      execute_if(uctrsat, [&]() {
        static_loop<NUM_TABLES>([&]<u64 I>() {
          std::cerr << "UC: tage u_ram[" << I << "] RESET\n";
          std::get<I>(tables).u_ram.reset();
        });
      });
    }

#endif // --- TAGE update disabled ---

    // ---- History update ----
    true_block = arr<val<1>, 4>{~mispredict, branch_dir[num_branch - 1],
                                last_pred(), line_end()}
                     .fold_or();
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

    std::cerr << "UC: EXIT (full update)\n";
    num_branch = 0;
  }
};

// ============================================================================
// User-facing Alias
// ============================================================================

template <typename TableCfg = td::TDTableConfig<8, 512, 11, 1, 2, 1, 2, 100, 4>,
          typename AllocCfg = td::TDDefaultAllocConfig, u64 TD_LINEINST = 1024,
          u64 TD_N = 7, bool TD_SHARED_TAG = true, bool TD_SHARED_U = true,
          bool TD_SHARED_HYS = true, bool TD_U_STOR_FF = false,
          u64 TD_DECAY_CTR = 8, u64 TD_DECAY_GRAN = 2,
          typename TD_DECAY_POLICY = td::TDDecayMild,
          bool TD_P1_USE_GSHARE = true, u64 TD_P1_TABLE_SIZE = 4096,
          u64 TD_P1_HIST = 6, bool TD_USE_META = true, u64 TD_METABITS = 4,
          u64 TD_METAPIPE = 2, bool TD_USE_PATH_HIST = false,
          u64 TD_PATH_HIST_WIDTH = 27, u64 TD_PATH_BITS = 6,
          template <u64> class TD_FOLD_FN = td::XORFold>

using TageDirect =
    TageDirectImpl<TableCfg, AllocCfg, TD_LINEINST, TD_N, TD_SHARED_TAG,
                   TD_SHARED_U, TD_SHARED_HYS, TD_U_STOR_FF, TD_DECAY_CTR,
                   TD_DECAY_GRAN, TD_DECAY_POLICY, TD_P1_USE_GSHARE,
                   TD_P1_TABLE_SIZE, TD_P1_HIST, TD_USE_META, TD_METABITS,
                   TD_METAPIPE, TD_USE_PATH_HIST, TD_PATH_HIST_WIDTH,
                   TD_PATH_BITS, TD_FOLD_FN>;
