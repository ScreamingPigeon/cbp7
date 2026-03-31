#pragma once

#include "../../harcom.hpp"
#include "../common.hpp"
#include "TageTable.hpp"

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

// Constexpr power: base^exp via exp(exp * ln(base)).
// std::exp and std::log are constexpr in C++20 (GCC 12+).
constexpr double constexpr_pow(double base, double exp) {
  if (exp == 0.0)
    return 1.0;
  if (base == 0.0)
    return 0.0;
  if (exp < 0.0)
    return 1.0 / constexpr_pow(base, -exp);
  return std::exp(exp * std::log(base));
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

// Quadratic history series: h[n] = h[n-1] + d*(n-1) + k
// Produces: h1, h1+k, h1+d+k, h1+3d+k, h1+6d+k, ...
// Index 0 = longest history (reversed to match geometric_hist convention).
template <std::size_t N>
constexpr std::array<u64, N> quadratic_hist(u64 h1, u64 d, u64 k) {
  static_assert(N >= 1);
  std::array<u64, N> h{};
  h[0] = h1;
  for (std::size_t n = 1; n < N; n++)
    h[n] = h[n - 1] + d * n + k;
  // Reverse: index 0 = longest
  for (std::size_t i = 0; i < N / 2; i++)
    std::swap(h[i], h[N - 1 - i]);
  return h;
}

// Super-exponential history series: multiplier grows linearly with n.
// h[n] = round(h[n-1] * (f*(n-t+1) + m)) for n >= 0
// Where t is the starting index (0-based into this array).
// Produces accelerating growth.
template <std::size_t N>
constexpr std::array<u64, N> superexp_hist(u64 h0, double f, double m) {
  static_assert(N >= 1);
  std::array<u64, N> h{};
  h[0] = h0;
  for (std::size_t n = 1; n < N; n++) {
    double mult = f * static_cast<double>(n) + m;
    u64 next = static_cast<u64>(static_cast<double>(h[n - 1]) * mult + 0.5);
    h[n] = (next > h[n - 1] + 1) ? next : h[n - 1] + 1;
  }
  // Reverse: index 0 = longest
  for (std::size_t i = 0; i < N / 2; i++)
    std::swap(h[i], h[N - 1 - i]);
  return h;
}

// Combined quadratic-to-super-exponential (Ros formula).
// Quadratic phase for n < t, super-exponential for n >= t.
// Ros defaults: h1=2, d=2, k=1, f=0.1, m=1.1, t=15
template <std::size_t N>
constexpr std::array<u64, N> ros_hist(u64 h1 = 2, u64 d = 2, u64 k = 1,
                                      double f = 0.1, double m = 1.1,
                                      std::size_t t = 15) {
  static_assert(N >= 1);
  std::array<u64, N> h{};
  h[0] = h1;
  for (std::size_t n = 1; n < N; n++) {
    if (n < t) {
      h[n] = h[n - 1] + d * n + k;
    } else {
      double mult = f * static_cast<double>(n - t + 1) + m;
      u64 next = static_cast<u64>(static_cast<double>(h[n - 1]) * mult + 0.5);
      h[n] = (next > h[n - 1] + 1) ? next : h[n - 1] + 1;
    }
  }
  // Reverse: index 0 = longest
  for (std::size_t i = 0; i < N / 2; i++)
    std::swap(h[i], h[N - 1 - i]);
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
// Table Size Generation
// ============================================================================

// Generate per-table sizes via a callable.
template <std::size_t N, typename Fn>
constexpr std::array<u64, N> generate_table_sizes(Fn fn) {
  std::array<u64, N> s{};
  for (std::size_t i = 0; i < N; i++)
    s[i] = fn(i, N);
  return s;
}

// Geometric size scaling — short-history tables larger, long-history smaller.
// Index 0 = longest history (smallest table), index N-1 = shortest (largest).
// SIZE_RATIO=1: uniform. SIZE_RATIO=2: ~2x range.
template <u64 SIZE, u64 SIZE_RATIO>
constexpr auto size_fn = [](u64 i, u64 n) -> u64 {
  if constexpr (SIZE_RATIO <= 1)
    return SIZE;
  else {
    double t = double(i) / std::max(1.0, double(n - 1));
    double scale = constexpr_pow(double(SIZE_RATIO), t - 0.5);
    u64 sz = u64(SIZE * scale);
    u64 result = 64;
    while (result < sz)
      result *= 2;
    return result;
  }
};

// ============================================================================
// Config Structs
// ============================================================================

enum class HistSeries { GEOMETRIC, QUADRATIC, SUPEREXP, ROS };

// Per-table config. SweepTableConfig<> with no args gives defaults:
// 8 tables, 2048 entries, TAG=11, CTR=1, HYST=2, U=1, hist 2..100 geometric.
template <u64 N = 8, u64 SIZE = 2048, u64 TAG = 11, u64 CTR = 1, u64 HYST = 2,
          u64 U = 1, u64 MINH = 2, u64 MAXH = 100, u64 SIZE_RATIO = 1,
          HistSeries HIST = HistSeries::GEOMETRIC>
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
  static constexpr auto HIST_LEN = []() {
    if constexpr (HIST == HistSeries::GEOMETRIC)  return geometric_hist<N>(MINH, MAXH);
    else if constexpr (HIST == HistSeries::QUADRATIC) return quadratic_hist<N>(MINH, 2, 1);
    else if constexpr (HIST == HistSeries::SUPEREXP)  return superexp_hist<N>(MINH, 0.1, 1.1);
    else if constexpr (HIST == HistSeries::ROS)       return ros_hist<N>(MINH, 2, 1, 0.1, 1.1, 15);
  }();
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
          u64 DECAY_CTR, typename ResetFn, bool USE_FF_CACHE, u64 I>
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
          u64 DECAY_CTR, typename ResetFn, bool USE_FF_CACHE, u64... Is>
struct MakeTableTuple<Cfg, BR_P_ENTRY, NUM_BANKS, USE_AHEAD, SHARED_TAG,
                      SHARED_U, SHARED_HYS, U_STOR_FF, DECAY_CTR, ResetFn,
                      USE_FF_CACHE, std::index_sequence<Is...>> {
  using type = std::tuple<
      TableAt<Cfg, BR_P_ENTRY, NUM_BANKS, USE_AHEAD, SHARED_TAG, SHARED_U,
              SHARED_HYS, U_STOR_FF, DECAY_CTR, ResetFn, USE_FF_CACHE, Is>...>;
};
