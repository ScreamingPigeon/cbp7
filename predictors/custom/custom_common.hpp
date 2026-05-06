#ifndef CUSTOM_COMMON_H
#define CUSTOM_COMMON_H

#include "../../harcom.hpp"

using namespace hcm;

// ============================================================================
// update_ctr — saturating up-down counter
// ============================================================================

template <u64 N, typename T>
[[nodiscard]] val<N, T> ta_update_ctr(val<N, T> ctr, val<1> incr) {
  // NOTE: @prakhar @claude ensure hardcoded fanout is correct
  // 6 reads: ==maxval, ==minval, +1, -1, select(ctr), select(ctr)
  ctr.fanout(hard<6>{});
  val<N, T> incsat = select(ctr == hard<ctr.maxval>{}, ctr, val<N, T>{ctr + 1});
  val<N, T> decsat = select(ctr == hard<ctr.minval>{}, ctr, val<N, T>{ctr - 1});
  return select(incr.fo1(), incsat.fo1(), decsat.fo1());
}

// ============================================================================
// ta_global_history — shift register updated by XOR with branch bits
// ============================================================================

template <u64 N> struct ta_global_history {
  arr<reg<1>, N> h;

  void update(valtype auto in) {
    auto input = in.fo1().make_array(val<1>{});
    static_assert(input.size <= N);
    for (u64 i = N - 1; i >= input.size; i--)
      h[i] = h[i - 1];
    for (u64 i = input.size - 1; i >= 1; i--)
      h[i] = h[i - 1] ^ input[i].fo1();
    h[0] = input[0].fo1();
  }

  // Unconditional write with enable mux — avoids execute_if gate timing.
  // When enable=1: shift and XOR in new input (normal update).
  // When enable=0: hold current state (registers still written, mux selects old
  // value).
  void update(valtype auto in, val<1> enable) {
    auto input = in.fo1().make_array(val<1>{});
    static_assert(input.size <= N);
    enable.fanout(hard<N>{});
    for (u64 i = N - 1; i >= input.size; i--)
      h[i] = select(enable, h[i - 1], h[i]);
    for (u64 i = input.size - 1; i >= 1; i--)
      h[i] = select(enable, h[i - 1] ^ input[i].fo1(), h[i]);
    h[0] = select(enable, input[0].fo1(), h[0]);
  }

  val<1> &operator[](u64 i) { return h[i]; }
  void fanout(hardval auto fo) { h.fanout(fo); }

  // Per-bit fanout: set each bit's fanout from a constexpr array
  template <auto const &FO_ARRAY> void fanout_per_bit() {
    static_assert(FO_ARRAY.size() == N);
    static_loop<N>([&]<u64 I>() { h[I].fanout(hard<FO_ARRAY[I]>{}); });
  }
};

// ============================================================================
// ta_folded_gh — incremental folded history with compute/apply split
//
// compute_update() returns the new fold value WITHOUT writing the register.
// apply_update() writes a precomputed value into the register.
// update() does both (compute + apply) for convenience.
//
// The split allows computing new fold values in parallel with a gate signal
// (e.g. true_block) and only applying writes inside execute_if, changing
// timing from additive to max(fold_computation, gate_signal).
// ============================================================================

template <u64 F> struct ta_folded_gh {
  static_assert(F != 0);

  reg<F> folded;

  val<F> get() { return folded; }
  void fanout(hardval auto fo) { folded.fanout(fo); }

  // Compute new fold value without writing the register
  template <u64 MAXL>
  [[nodiscard]] inline val<F> compute_update(ta_global_history<MAXL> &gh,
                                             hardval auto ghlen,
                                             valtype auto in) {
    constexpr u64 inbits = std::min(F, std::min(in.size, ghlen.value));
    val<inbits> input = in.fo1();
    auto f = folded.make_array(val<1>{});
    static_assert(f.size == F);
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
    return y.fo1().concat();
  }

  // Write a precomputed fold value into the register
  void apply_update(val<F> new_val) { folded = new_val.fo1(); }

  // Unconditional write with enable mux — avoids execute_if gate timing.
  // When enable=1: write new_val. When enable=0: hold current value.
  void apply_update(val<F> new_val, val<1> enable) {
    folded.fanout(hard<2>{}); // hold-value read in select below (single use)
    folded = select(enable.fo1(), new_val.fo1(), folded);
  }

  // Combined compute + apply (convenience, same as original folded_gh::update)
  template <u64 MAXL>
  void update(ta_global_history<MAXL> &gh, hardval auto ghlen,
              valtype auto in) {
    apply_update(compute_update(gh, ghlen, in));
  }
};

// ============================================================================
// History update mode for geometric_folds_ex
// ============================================================================

enum class HistUpdate { PATH, DIR, BOTH };
enum class DecayMiss { TAG, SEC, TAG_OR_SEC, TAG_AND_SEC };
enum class DecayOp { DECREMENT, HALVE, CLEAR };
enum class SiblingPolicy { NONE, ALL };
// What happens to u-bit when provider mispredicts
enum class UMispPolicy { UNTOUCHED, ZERO, DECREMENT };
// Alloc-failure u-bit policy for tables above provider
enum class UClearPolicy { DISABLED, ZERO, DECREMENT };

// ============================================================================
// geometric_folds_ex — geometric_folds with templated update mode
//
// Extends geometric_folds with three update modes:
//   update<PATH>(path_bits)            — fold in path bits only
//   update<DIR>(dir_bit)               — fold in direction bit only
//   update<BOTH>(dir_bit, path_bits)   — fold in concat(dir, path)
// ============================================================================

template <u64 NH, u64 MINH, u64 MAXH, u64... FOLDS> struct geometric_folds_ex {
  static_assert(NH >= 2);
  static constexpr u64 NF = sizeof...(FOLDS);

  static constexpr auto HLEN = []() {
    std::array<u64, NH> hlen;
    u64 prevhl = 0;
    for (u64 i = 0; i < NH; i++) {
      u64 hl = MINH * mypow(f64(MAXH) / MINH, f64(i) / (NH - 1));
      hl = std::max(prevhl + 1, hl);
      hlen[NH - 1 - i] = hl;
      prevhl = hl;
    }
    return hlen;
  }();

  static_assert(HLEN[0] == MAXH);

  ta_global_history<MAXH> gh;
  std::array<std::tuple<ta_folded_gh<FOLDS>...>, NH> folds;

  template <u64 J = 0> auto get(u64 i) {
    if (i >= NH) {
      std::cerr << "geometric_folds_ex: out of bound access\n";
      std::terminate();
    }
    return std::get<J>(folds[i]).get();
  }

  void fanout(hardval auto fo) {
    for (u64 i = 0; i < NH; i++) {
      static_loop<NF>([&]<u64 J>() { std::get<J>(folds[i]).fanout(fo); });
    }
  }

  // Internal: common update logic — feed bits into all folds then gh
  void do_update(valtype auto bits) {
    bits.fanout(hard<NH * NF + 1>{});
    gh.fanout(hard<std::max(u64(2), NF + 1)>{});
    static_loop<NH>([&]<u64 I>() {
      static_loop<NF>([&]<u64 J>() {
        std::get<J>(folds[I]).update(gh, hard<HLEN[I]>{}, bits);
      });
    });
    gh.update(bits);
  }

  // update<PATH>(path_bits) — path only
  template <HistUpdate M>
  void update(valtype auto path_bits)
    requires(M == HistUpdate::PATH)
  {
    do_update(path_bits);
  }

  // update<DIR>(dir_bit) — direction only (1-bit)
  template <HistUpdate M>
  void update(val<1> dir_bit)
    requires(M == HistUpdate::DIR)
  {
    do_update(dir_bit);
  }

  // update<BOTH>(dir_bit, path_bits) — concat(dir, path)
  template <HistUpdate M>
  void update(val<1> dir_bit, valtype auto path_bits)
    requires(M == HistUpdate::BOTH)
  {
    do_update(concat(dir_bit, path_bits));
  }
};

// ============================================================================
// Constexpr Helpers
// ============================================================================

namespace ta {

template <typename T, std::size_t N>
constexpr std::array<T, N> uniform_array(T v) {
  std::array<T, N> a{};
  for (std::size_t i = 0; i < N; i++)
    a[i] = v;
  return a;
}

// Two-zone array: v_lo for indices [0, SPLIT), v_hi for [SPLIT, N).
template <typename T, std::size_t N>
constexpr std::array<T, N> split_array(T v_lo, T v_hi, std::size_t split) {
  std::array<T, N> a{};
  for (std::size_t i = 0; i < N; i++)
    a[i] = (i < split) ? v_lo : v_hi;
  return a;
}

// Graded array: linear interpolation from v_first (index 0) to v_last (index N-1).
// Useful for per-table LFSR widths where T0 (longest history) differs from T(N-1).
template <typename T, std::size_t N>
constexpr std::array<T, N> graded_array(T v_first, T v_last) {
  std::array<T, N> a{};
  using S = std::make_signed_t<T>;
  for (std::size_t i = 0; i < N; i++)
    a[i] = N > 1 ? T(S(v_first) + (S(v_last) - S(v_first)) * S(i) / S(N - 1))
                  : v_first;
  return a;
}

constexpr u64 clog2(u64 x) {
  u64 r = 0, v = x - 1;
  while (v > 0) {
    v >>= 1;
    r++;
  }
  return r;
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

// Compute per-bit gh fanout: for each bit position 0..MAXHIST-1,
// count how many tables need gh[HIST_LEN[I]-1] at that position.
// Each table reads the outgoing bit twice (fold_idx + fold_tag).
// Add 1 for the gh.update() shift read.
template <std::size_t MAXHIST, std::size_t NT, auto const &HIST_LEN>
constexpr std::array<u64, MAXHIST> gh_per_bit_fanout() {
  std::array<u64, MAXHIST> fo{};
  // Base: each bit is read by update() for shift (h[i-1]) and mux hold (h[i])
  for (std::size_t i = 0; i < MAXHIST; i++)
    fo[i] = 3;
  for (std::size_t t = 0; t < NT; t++) {
    u64 bit = HIST_LEN[t] - 1;
    fo[bit] += 2; // fold_idx + fold_tag read the outgoing bit
  }
  return fo;
}

constexpr double constexpr_pow(double base, double exp) {
  if (exp == 0.0)
    return 1.0;
  if (base == 0.0)
    return 0.0;
  return std::exp(exp * std::log(base));
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

template <u64 SIZE, u64 RATIO> struct InvGeoSize {
  constexpr u64 operator()(u64 i, u64 n) const {
    if (RATIO <= 1 || n <= 1)
      return SIZE;
    double t = double(n - 1 - i) / double(n - 1);
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

template <u64 TOTAL_BITS, u64 TAG, u64 CTR, u64 HYST, u64 U>
struct ConstBitsSize {
  constexpr u64 operator()(u64, u64) const {
    u64 entry_bits = TAG + CTR + HYST + U;
    u64 entries = TOTAL_BITS / entry_bits;
    u64 result = 64;
    while (result * 2 <= entries)
      result *= 2;
    return result;
  }
};

// Linear interpolation from S_HI (table 0) to S_LO (table N-1), pow2 rounded.
// Mirrors GradedTag.
template <u64 S_HI, u64 S_LO> struct GradedSize {
  constexpr u64 operator()(u64 i, u64 n) const {
    if (n <= 1)
      return S_HI;
    double t = double(i) / double(n - 1);
    u64 sz = u64(double(S_HI) * (1.0 - t) + double(S_LO) * t);
    u64 result = 64;
    while (result < sz)
      result *= 2;
    return result;
  }
};

// Logarithmic scaling: larger tables for shorter histories (higher index).
// BASE + SCALE * level / 4, pow2 rounded. Mirrors LogTag.
template <u64 BASE, u64 SCALE> struct LogSize {
  constexpr u64 operator()(u64 i, u64 n) const {
    u64 level = (n > 1) ? (n - 1 - i) : 0;
    u64 sz = BASE + SCALE * level / 4;
    u64 result = 64;
    while (result < sz)
      result *= 2;
    return result;
  }
};

// Three-tier step: S1 for [0, SPLIT1), S2 for [SPLIT1, SPLIT2), S3 for rest.
template <u64 S1, u64 S2, u64 S3, u64 SPLIT1, u64 SPLIT2> struct MultiStepSize {
  static_assert(SPLIT1 <= SPLIT2, "SPLIT1 must be <= SPLIT2");
  constexpr u64 operator()(u64 i, u64) const {
    if (i < SPLIT1)
      return S1;
    if (i < SPLIT2)
      return S2;
    return S3;
  }
};

// Fixed total entry budget distributed equally across tables, pow2 rounded.
// Each table gets TOTAL/N entries (clamped to >= MIN_SIZE).
template <u64 TOTAL, u64 MIN_SIZE = 64> struct BudgetSize {
  constexpr u64 operator()(u64, u64 n) const {
    u64 per = (n > 0) ? TOTAL / n : TOTAL;
    u64 result = 64;
    while (result * 2 <= per)
      result *= 2;
    return (result < MIN_SIZE) ? MIN_SIZE : result;
  }
};

// Fixed total entry budget weighted by inverse history length.
// Shorter-history tables (higher index) get more entries. Pow2 rounded.
// Requires access to HIST_LEN at TATableConfig level — use via generate overload.
template <u64 TOTAL, u64 MIN_SIZE = 64> struct InvHistBudgetSize {
  // Weight for table i = 1/sqrt(hist_len[i]).
  // Called from generate_table_sizes with hist_len array.
  template <std::size_t N>
  constexpr std::array<u64, N>
  generate(const std::array<u64, N> &hist_len) const {
    std::array<u64, N> sizes{};
    double total_weight = 0.0;
    std::array<double, N> w{};
    for (std::size_t i = 0; i < N; i++) {
      w[i] = 1.0 / constexpr_pow(double(hist_len[i]), 0.5);
      total_weight += w[i];
    }
    for (std::size_t i = 0; i < N; i++) {
      u64 per = u64(double(TOTAL) * w[i] / total_weight);
      u64 result = 64;
      while (result * 2 <= per)
        result *= 2;
      sizes[i] = (result < MIN_SIZE) ? MIN_SIZE : result;
    }
    return sizes;
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
// Decay helpers
// ============================================================================

// Generate a tuple of reg<W> with per-element widths from a constexpr array.
// Usage: TALfsrTuple<LFSR_WIDTHS, std::make_index_sequence<NT>>::type
template <auto const &Widths, typename Seq> struct TALfsrTuple;

template <auto const &Widths, u64... Is>
struct TALfsrTuple<Widths, std::index_sequence<Is...>> {
  using type = std::tuple<reg<std::max(u64(1), Widths[Is])>...>;
};

// Default threshold functor: threshold = alloc_ctr (ignore accuracy).
// Users can define custom functors with the same interface.
struct DefaultDecayThresh {
  // Returns a val<LW> threshold from the global counters.
  // I = table index (allows per-table differentiation).
  template <u64 I, u64 LW, u64 AW, u64 PW>
  static auto compute(val<AW> /*accuracy_ctr*/, val<PW> alloc_ctr) {
    return val<LW>{alloc_ctr};
  }
};

// Fixed threshold: compile-time constant probability independent of counters.
// P(decay|miss) ≈ THRESH / 2^LW.  E.g. THRESH=16, LW=8 → ~6.25% per miss.
template <u64 THRESH>
struct FixedDecayThresh {
  template <u64 I, u64 LW, u64 AW, u64 PW>
  static auto compute(val<AW> /*acc*/, val<PW> /*alloc*/) {
    constexpr u64 clamped = THRESH < (u64(1) << LW) ? THRESH : (u64(1) << LW) - 1;
    return val<LW>{hard<clamped>{}};
  }
};

// Graded fixed threshold: per-table constant.  T0 (longest history, index 0)
// gets THRESH_LO (decays rarely), T(N-1) gets THRESH_HI (decays more often).
// Linear interpolation between them.
template <u64 THRESH_LO, u64 THRESH_HI, u64 NT>
struct GradedDecayThresh {
  template <u64 I, u64 LW, u64 AW, u64 PW>
  static auto compute(val<AW> /*acc*/, val<PW> /*alloc*/) {
    constexpr u64 t = NT > 1 ? THRESH_LO + (THRESH_HI - THRESH_LO) * I / (NT - 1)
                              : THRESH_LO;
    constexpr u64 clamped = t < (u64(1) << LW) ? t : (u64(1) << LW) - 1;
    return val<LW>{hard<clamped>{}};
  }
};

// Accuracy-gated threshold: threshold = FIXED * (1 - acc_ctr/maxval).
// High accuracy (acc_ctr near max) → low threshold → rare decay.
// Low accuracy (acc_ctr near 0) → threshold near FIXED → frequent decay.
template <u64 FIXED>
struct AccGatedDecayThresh {
  template <u64 I, u64 LW, u64 AW, u64 PW>
  static auto compute(val<AW> acc_ctr, val<PW> /*alloc*/) {
    // Invert acc: low accuracy → high value
    val<AW> inv_acc = val<AW>{hard<acc_ctr.maxval>{}} - acc_ctr;
    // Scale: thresh = FIXED * inv_acc / maxval ≈ (FIXED * inv_acc) >> AW
    // Approximate with shift to avoid HW multiply
    constexpr u64 clamped = FIXED < (u64(1) << LW) ? FIXED : (u64(1) << LW) - 1;
    return val<LW>{(val<AW + LW>{inv_acc} * hard<clamped>{}) >> AW};
  }
};

// Pressure-gated graded threshold: only decay when alloc_ctr exceeds GATE.
// Below GATE → threshold=0 (no decay). Above GATE → graded per-table threshold.
// Addresses workloads with low table contention (e.g. namd) where decay is harmful.
template <u64 THRESH_LO, u64 THRESH_HI, u64 NT, u64 GATE>
struct PressGatedDecayThresh {
  template <u64 I, u64 LW, u64 AW, u64 PW>
  static auto compute(val<AW> /*acc*/, val<PW> alloc_ctr) {
    constexpr u64 t = NT > 1 ? THRESH_LO + (THRESH_HI - THRESH_LO) * I / (NT - 1)
                              : THRESH_LO;
    constexpr u64 clamped = t < (u64(1) << LW) ? t : (u64(1) << LW) - 1;
    val<1> gate = alloc_ctr.fo1() > hard<GATE>{};
    return select(gate, val<LW>{hard<clamped>{}}, val<LW>{hard<0>{}});
  }
};

// Pressure-scaled graded threshold: scales threshold linearly with alloc_ctr.
// thresh[I] = graded_base[I] * alloc_ctr / max_alloc_ctr
// High pressure → full threshold. Low pressure → near-zero threshold.
template <u64 THRESH_LO, u64 THRESH_HI, u64 NT>
struct PressScaledDecayThresh {
  template <u64 I, u64 LW, u64 AW, u64 PW>
  static auto compute(val<AW> /*acc*/, val<PW> alloc_ctr) {
    constexpr u64 t = NT > 1 ? THRESH_LO + (THRESH_HI - THRESH_LO) * I / (NT - 1)
                              : THRESH_LO;
    constexpr u64 clamped = t < (u64(1) << LW) ? t : (u64(1) << LW) - 1;
    return val<LW>{(val<PW + LW>{alloc_ctr} * hard<clamped>{}) >> PW};
  }
};

// Default epoch trigger: fire when alloc_ctr saturates.
struct DefaultEpochTrigger {
  template <u64 AW, u64 PW>
  static val<1> should_fire(val<AW> /*acc_ctr*/, val<PW> alloc_ctr) {
    return alloc_ctr.fo1() == hard<alloc_ctr.maxval>{};
  }
};

// ============================================================================
// SecTagPolicy — per-table sec-tag enforcement
// ============================================================================
// Controls whether the sec-tag match is required for each TAGE table.
// apply<I, NT>() returns true (at compile time) if table I requires sec-tag.
// For pressure-gated variants, apply_runtime<I, NT>(alloc_ctr) returns a
// val<1> gate signal (1 = enforce sec-tag, 0 = skip).
//
// The RUNTIME flag distinguishes compile-time-only policies (no extra HW)
// from runtime-gated policies (adds a comparator on critical path).

// SecTagAll: all tables check sec-tag (default, current behavior)
struct SecTagAll {
  static constexpr bool RUNTIME = false;
  template <u64 I, u64 NT> static constexpr bool apply() { return true; }
};

// SecTagNone: no tables check sec-tag
struct SecTagNone {
  static constexpr bool RUNTIME = false;
  template <u64 I, u64 NT> static constexpr bool apply() { return false; }
};

// SecTagFloor<F>: only tables with index >= F check sec-tag.
// T0..T(F-1) skip sec-tag (long-history tables — namd-style entries survive).
// T(F)..T(NT-1) still check (short-history tables — stale entries rejected).
template <u64 F>
struct SecTagFloor {
  static constexpr bool RUNTIME = false;
  template <u64 I, u64 NT> static constexpr bool apply() { return I >= F; }
};

// SecTagCeil<C>: only tables with index < C check sec-tag.
// T0..T(C-1) check (long history — more likely to be stale).
// T(C)..T(NT-1) skip (short history — less likely to be stale).
template <u64 C>
struct SecTagCeil {
  static constexpr bool RUNTIME = false;
  template <u64 I, u64 NT> static constexpr bool apply() { return I < C; }
};

// SecTagMask<MASK>: bitmask — bit I set = table I checks sec-tag.
template <u64 MASK>
struct SecTagMask {
  static constexpr bool RUNTIME = false;
  template <u64 I, u64 NT> static constexpr bool apply() {
    return (MASK >> I) & 1;
  }
};

// SecTagPressGated<GATE>: skip sec-tag for ALL tables when alloc_ctr > GATE.
// Under high pressure, entries are churning fast → sec-tag mismatch is common
// but the TAGE entry might still be useful. Under low pressure, sec-tag filter
// remains active to reject stale entries.
template <u64 GATE>
struct SecTagPressGated {
  static constexpr bool RUNTIME = true;
  template <u64 I, u64 NT> static constexpr bool apply() { return true; }
  template <u64 I, u64 NT, u64 AW, u64 PW, u64 BW>
  static val<1> gate(val<AW> /*acc*/, val<PW> alloc_ctr, val<BW> /*benefit*/) {
    return val<1>{alloc_ctr.fo1() <= hard<GATE>{}};
  }
};

// SecTagPressGatedFloor<F, GATE>: combine floor + pressure gating.
// Tables < F always skip sec-tag. Tables >= F check sec-tag only when
// alloc_ctr <= GATE (low pressure).
template <u64 F, u64 GATE>
struct SecTagPressGatedFloor {
  static constexpr bool RUNTIME = true;
  template <u64 I, u64 NT> static constexpr bool apply() { return I >= F; }
  template <u64 I, u64 NT, u64 AW, u64 PW, u64 BW>
  static val<1> gate(val<AW> /*acc*/, val<PW> alloc_ctr, val<BW> /*benefit*/) {
    if constexpr (I < F) return val<1>{hard<0>{}};
    else return val<1>{alloc_ctr.fo1() <= hard<GATE>{}};
  }
};

// SecTagAccGated<GATE>: skip sec-tag when accuracy is high (acc_ctr > GATE).
// When accuracy is high, TAGE entries are likely valid → sec-tag mismatch
// is noise, not signal. When accuracy drops, re-enable sec-tag filtering.
template <u64 GATE>
struct SecTagAccGated {
  static constexpr bool RUNTIME = true;
  template <u64 I, u64 NT> static constexpr bool apply() { return true; }
  template <u64 I, u64 NT, u64 AW, u64 PW, u64 BW>
  static val<1> gate(val<AW> acc_ctr, val<PW> /*alloc*/, val<BW> /*benefit*/) {
    return val<1>{acc_ctr.fo1() <= hard<GATE>{}};
  }
};

// SecTagAdaptive<WIDTH, THRESH>: benefit-tracking adaptive sec-tag policy.
// Maintains a saturating up/down counter updated at training time:
//   - When sec-tag rejects a TAGE match and fallback was correct but TAGE wrong
//     → increment (sec-tag helped)
//   - When sec-tag rejects a TAGE match and TAGE was correct but fallback wrong
//     → decrement (sec-tag hurt)
// At prediction time, enforce sec-tag when benefit_ctr >= THRESH.
template <u64 WIDTH, u64 THRESH>
struct SecTagAdaptive {
  static constexpr bool RUNTIME = true;
  static constexpr u64 BENEFIT_WIDTH = WIDTH;
  static constexpr u64 BENEFIT_THRESH = THRESH;
  template <u64 I, u64 NT> static constexpr bool apply() { return true; }
  template <u64 I, u64 NT, u64 AW, u64 PW, u64 BW>
  static val<1> gate(val<AW> /*acc*/, val<PW> /*alloc*/, val<BW> benefit) {
    return val<1>{benefit.fo1() >= hard<THRESH>{}};
  }
};

// ============================================================================
// Default sec_tag hash: extract bits [2+SEC_TAG_BITS-1 : 2] from PC.
struct DefaultSecTagHash {
  template <u64 SEC_TAG_BITS> static val<SEC_TAG_BITS> apply(val<64> pc) {
    // pc is a by-value copy (read_credit=0); use fo1() for the single read
    return val<SEC_TAG_BITS>{pc.fo1() >> 2};
  }
};

// 3-way XOR sec_tag hash for 3-bit sec_tag: PC[4:2] ^ PC[9:7] ^ PC[14:12]
struct Xor3SecTagHash3 {
  template <u64 SEC_TAG_BITS> static val<SEC_TAG_BITS> apply(val<64> pc) {
    static_assert(SEC_TAG_BITS == 3, "Xor3SecTagHash3 is tuned for 3-bit sec_tag");
    pc.fanout(hard<3>{});
    return val<3>{pc >> 2} ^ val<3>{pc >> 7} ^ val<3>{pc >> 12};
  }
};

// PC[5:2] ^ PC[11:8] ^ PC[16:13]  — best 4-bit hash from block analyzer
struct Xor3SecTagHash {
  template <u64 SEC_TAG_BITS> static val<SEC_TAG_BITS> apply(val<64> pc) {
    static_assert(SEC_TAG_BITS == 4, "Xor3SecTagHash is tuned for 4-bit sec_tag");
    pc.fanout(hard<3>{});
    return val<4>{pc >> 2} ^ val<4>{pc >> 8} ^ val<4>{pc >> 13};
  }
};

// 3-way XOR sec_tag hash for 5-bit sec_tag: PC[6:2] ^ PC[13:9] ^ PC[20:16]
struct Xor3SecTagHash5 {
  template <u64 SEC_TAG_BITS> static val<SEC_TAG_BITS> apply(val<64> pc) {
    static_assert(SEC_TAG_BITS == 5, "Xor3SecTagHash5 is tuned for 5-bit sec_tag");
    pc.fanout(hard<3>{});
    return val<5>{pc >> 2} ^ val<5>{pc >> 9} ^ val<5>{pc >> 16};
  }
};

// XOR-folded sec_tag hash: XOR two non-overlapping bit slices of PC.
struct XorSecTagHash {
  template <u64 SEC_TAG_BITS> static val<SEC_TAG_BITS> apply(val<64> pc) {
    // pc is a by-value copy (read_credit=0); two reads — declare fanout
    pc.fanout(hard<2>{});
    return val<SEC_TAG_BITS>{pc >> 2} ^
           val<SEC_TAG_BITS>{pc >> (2 + SEC_TAG_BITS)};
  }
};

} // namespace ta

// ============================================================================
// rwram — Banked RAM with configurable bank index bits
//
// Allows simultaneous read and write (to different banks). If both access
// the same bank, the write is buffered and flushed when the bank is free.
//
// Template params:
//   N = data width in bits
//   M = total entries (power of 2)
//   B = number of banks (power of 2, >= 2)
//   BANK_SHIFT = which bit position starts bank selection (default 0)
//                0 = low bits select bank, K = bits [K, K+log2(B))
// ============================================================================

template <u64 N, u64 M, u64 B, u64 BANK_SHIFT = 0> struct ta_rwram {
  static_assert(std::has_single_bit(M));
  static_assert(B >= 2 && B <= 64);
  static_assert(std::has_single_bit(B));
  static constexpr u64 A = std::bit_width(M - 1); // address bits
  static constexpr u64 E = M / B;                 // entries per bank
  static_assert(E > 1);
  static constexpr u64 L = std::bit_width(E - 1);         // local address bits
  static constexpr u64 BANK_BITS = std::bit_width(B - 1); // bank ID bits
  static_assert(A == L + BANK_BITS);
  static_assert(BANK_SHIFT + BANK_BITS <= A);

  hcm::ram<val<N>, E> bank[B];
  reg<B> read_bank;

  // buffered write state
  reg<B> write_bank;
  reg<L> write_localaddr;
  reg<N> write_data;

  // conflict tracking (simulation only, no hardware cost)
  const char *name_ = "";
  u64 stat_writes = 0;        // total write() calls
  u64 stat_buffered = 0;      // writes buffered (nc=0, read-write same cycle)
  u64 stat_lost = 0;          // buffered writes overwritten before flush
  u64 stat_flushed = 0;       // buffered writes successfully flushed
  u64 pending_bank_ = 0;      // which bank has the pending buffered write (bitmask)
  u64 prev_read_bank_ = 0;    // previous cycle's read bank (bitmask)
  u64 prev_read_addr_ = 0;    // previous cycle's read address
  u64 stat_reads = 0;         // total read() calls
  u64 stat_read_same_prev = 0; // read bank == previous read bank (consecutive correlation)
  // XOR histogram: how many bits of the address change between consecutive reads
  // idx_xor_hist[i] = count of cycles where popcount(addr XOR prev_addr) == i
  static constexpr u64 XOR_HIST_SIZE = A + 1;
  u64 idx_xor_hist[XOR_HIST_SIZE] = {};
  // Bank-bits XOR histogram: how many of the BANK_BITS change
  u64 bank_xor_hist[BANK_BITS + 1] = {};
  // Per-bit flip rate: which bit positions flip between consecutive reads
  u64 bit_flip_count[A] = {};

  ta_rwram(const char *label = "") : bank{label}, name_{label} {}

  // Split full address into (local_addr, bank_id).
  auto split_addr(val<A> addr) {
    if constexpr (BANK_SHIFT == 0) {
      return hcm::split<L, BANK_BITS>(addr.fo1()); // split reads once
    } else if constexpr (BANK_SHIFT + BANK_BITS == A) {
      addr.fanout(hard<2>{}); // truncate + shift
      return std::pair{val<L>{addr}, val<BANK_BITS>{addr >> BANK_SHIFT}};
    } else {
      addr.fanout(hard<3>{}); // lo + bankid(>>) + hi(>>)
      val<BANK_SHIFT> lo = addr;
      val<BANK_BITS> bankid = addr >> BANK_SHIFT;
      val<A - BANK_SHIFT - BANK_BITS> hi = addr >> (BANK_SHIFT + BANK_BITS);
      val<L> localaddr = concat(lo, hi);
      return std::pair{localaddr, bankid};
    }
  }

  val<N> read(val<A> addr) {
    auto [localaddr, bankid] = split_addr(addr.fo1().connect(bank[0]));
    localaddr.fanout(hard<B>{}); // B bank reads
    arr<val<1>, B> banksel = bankid.fo1().decode();
    // NOTE: @prakhar @claude ensure hardcoded fanout is correct
    banksel.fanout(hard<2>{}); // execute_if + concat
    arr<val<N>, B> data = [&](u64 i) -> val<N> {
      return execute_if(banksel[i], [&]() { return bank[i].read(localaddr); });
    };
    read_bank = banksel.concat();
#ifdef TAGE_MONITOR
    {
      u64 rb = static_cast<u64>(read_bank);
      u64 av = static_cast<u64>(addr);
      stat_reads++;
      if (stat_reads > 1) {
        if (rb == prev_read_bank_)
          stat_read_same_prev++;
        u64 xor_val = av ^ prev_read_addr_;
        idx_xor_hist[std::popcount(xor_val)]++;
        for (u64 bit = 0; bit < A; bit++)
          if ((xor_val >> bit) & 1) bit_flip_count[bit]++;
        u64 bank_bits_now = (av >> BANK_SHIFT) & ((1u << BANK_BITS) - 1);
        u64 bank_bits_prev = (prev_read_addr_ >> BANK_SHIFT) & ((1u << BANK_BITS) - 1);
        bank_xor_hist[std::popcount(bank_bits_now ^ bank_bits_prev)]++;
      }
      prev_read_bank_ = rb;
      prev_read_addr_ = av;
    }
#endif
    return data.fo1().fold_or();
  }

  void write(val<A> addr, val<N> data, val<1> nc) {
    // noconflict=1: no read this cycle, write immediately.
    // noconflict=0: buffer write, flush when bank is free.
    // --- conflict tracking (simulation-only, monitor builds) ---
#ifdef TAGE_MONITOR
    stat_writes++;
    {
      u64 nc_val = static_cast<u64>(nc);
      u64 read_bitmask = static_cast<u64>(read_bank);
      if (pending_bank_ != 0) {
        if ((pending_bank_ & read_bitmask) == 0) {
          stat_flushed++;
        } else {
          if (nc_val == 0) stat_lost++;
        }
      }
      if (nc_val == 0) {
        stat_buffered++;
        u64 addr_val = static_cast<u64>(addr);
        u64 new_bank;
        if constexpr (BANK_SHIFT == 0) {
          new_bank = 1u << ((addr_val & ((1u << BANK_BITS) - 1)));
        } else {
          new_bank = 1u << (((addr_val >> BANK_SHIFT) & ((1u << BANK_BITS) - 1)));
        }
        pending_bank_ = new_bank;
      } else {
        if (pending_bank_ != 0 && (pending_bank_ & read_bitmask) == 0)
          pending_bank_ = 0;
      }
    }
#endif
    // --- end conflict tracking ---
    // Connect narrow noconflict (1-bit) to bank[0] so operations happen
    // at the bank location, reducing wire distance.
    auto noconflict = nc.fo1().connect(bank[0]);
    auto [localaddr, bankid] = split_addr(addr.fo1().connect(bank[0]));
    localaddr.fanout(
        hard<B + 1>{});         // B selects in loop + buffered write_localaddr
    data.fanout(hard<B + 1>{}); // B bank writes + buffered write_data
    noconflict.fanout(hard<B + 2>{}); // B bank selects + mask + buffered gate
    val<B> banksel = bankid.fo1().decode().concat();
    // NOTE: @prakhar @claude ensure hardcoded fanout is correct
    banksel.fanout(hard<2>{}); // current_write AND + buffered write_bank
    val<B> noconflict_mask = noconflict.replicate(hard<B>{}).concat();
    // NOTE: @prakhar @claude ensure hardcoded fanout is correct
    noconflict_mask.fanout(
        hard<2>{}); // current_write AND + buffered ~noconflict
    val<B> current_write = banksel & noconflict_mask;
    // NOTE: @prakhar @claude ensure hardcoded fanout is correct
    current_write.fanout(hard<3>{}); // make_array + buffered OR + read_bank OR
    arr<val<1>, B> current_write_split = current_write.make_array(val<1>{});
    // NOTE: @prakhar @claude ensure hardcoded fanout is correct
    current_write_split.fanout(hard<3>{}); // execute_if cond + 2 selects
    write_bank.fanout(hard<2>{});          // make_array(1) + buffered_done &(1)
    read_bank.fanout(hard<2>{});           // make_array(1) + buffered_done |(1)
    write_localaddr.fanout(hard<B>{}); // select in each of B execute_if lambdas
    write_data.fanout(hard<B>{});      // select in each of B execute_if lambdas
    arr<val<1>, B> write_bank_split = write_bank.make_array(val<1>{});
    arr<val<1>, B> read_bank_split = read_bank.make_array(val<1>{});
    for (u64 i = 0; i < B; i++) {
      execute_if(current_write_split[i] |
                     (write_bank_split[i].fo1() & ~read_bank_split[i].fo1()),
                 [&]() {
                   val<L> a = select(current_write_split[i], localaddr,
                                     write_localaddr);
                   val<N> d = select(current_write_split[i], data, write_data);
                   bank[i].write(a.fo1(), d.fo1());
                 });
    }
    // buffer the current write if not done
    val<1> buffered_done =
        (write_bank & (current_write | read_bank)) == hard<0>{};
    val<1> buffered_gate = buffered_done.fo1() | ~noconflict;
    // NOTE: @prakhar @claude validate hardcoded fanout
    buffered_gate.fanout(
        hard<3>{}); // write_bank(1) + write_localaddr(1) + write_data(1)
    execute_if(buffered_gate, [&]() {
      write_bank = banksel & ~noconflict_mask;
      execute_if(~noconflict, [&]() {
        write_localaddr = localaddr;
        write_data = data;
      });
    });
  }

  void reset() {
    for (u64 i = 0; i < B; i++)
      bank[i].reset();
  }

  void print_conflict_stats(std::ostream &os = std::cerr) const {
    if (stat_writes == 0) return;
    os << "rwram[" << name_ << "] writes=" << stat_writes
       << " buffered=" << stat_buffered
       << " (" << (100.0 * stat_buffered / stat_writes) << "%)"
       << " lost=" << stat_lost
       << " (" << (stat_buffered > 0 ? 100.0 * stat_lost / stat_buffered : 0.0) << "% of buffered)"
       << " flushed=" << stat_flushed
       << " read_same_prev=" << stat_read_same_prev
       << " (" << (stat_reads > 1 ? 100.0 * stat_read_same_prev / (stat_reads - 1) : 0.0) << "%)"
       << "\n";
    // Print XOR histograms (only for pred RAMs to avoid spam)
    if (stat_reads > 1) {
      u64 total = stat_reads - 1;
      os << "  idx_xor_popcount[" << name_ << "]:";
      for (u64 i = 0; i < XOR_HIST_SIZE; i++)
        if (idx_xor_hist[i] > 0)
          os << " " << i << ":" << (100.0 * idx_xor_hist[i] / total) << "%";
      os << "\n";
      os << "  bank_xor_popcount[" << name_ << "]:";
      for (u64 i = 0; i <= BANK_BITS; i++)
        os << " " << i << ":" << (100.0 * bank_xor_hist[i] / total) << "%";
      os << "\n";
      os << "  bit_flip_rate[" << name_ << "]:";
      for (u64 i = 0; i < A; i++)
        os << " b" << i << ":" << (100.0 * bit_flip_count[i] / total) << "%";
      os << "\n";
    }
  }
};

// ============================================================================
// TATable — Per-table storage for ahead-pipelined TAGE
//
// Each table has its own compile-time TAG_WIDTH and TABLE_SIZE.
// Pred RAMs are per-bank (BANKS copies). Hyst/u are update-only (plain ram).
// Pipeline regs (ahead[0]/[1]) are included per table.
// ============================================================================

template <u64 TABLE_SIZE, u64 TAG_WIDTH,
          u64 HIST_LEN,  // per-table history length
          u64 CTR_WIDTH, // per-lane prediction counter
          u64 HYST_WIDTH, u64 U_WIDTH, u64 SEC_TAG_BITS,
          u64 N,                   // max branches per block (= lanes of pred)
          bool SHARED_HYS = false, // shared hyst: 2 entries share 1 counter
          u64 TAG_RAM_WIDTH = TAG_WIDTH, // RAM width for tag (uniform across tables)
          u64 RWRAM_BANK_SHIFT = 0, // which address bit selects rwram bank
          u64 RWRAM_BANKS = 2>     // number of rwram banks (power of 2)
struct TATable {
  static constexpr u64 IDX_BITS = ta::clog2(TABLE_SIZE);
  static constexpr u64 PRED_BITS = N * CTR_WIDTH;
  static constexpr u64 table_size = TABLE_SIZE;
  static constexpr u64 tag_width = TAG_WIDTH;
  static constexpr u64 hist_len = HIST_LEN;
  static constexpr u64 ctr_width = CTR_WIDTH;
  static constexpr u64 hyst_width = HYST_WIDTH;
  static constexpr u64 u_width = U_WIDTH;
  static constexpr u64 sec_tag_bits = SEC_TAG_BITS;
  static constexpr u64 tag_ram_width = TAG_RAM_WIDTH;

  // When SHARED_HYS=true, halve the hyst table: pairs of entries share one
  // counter
  static constexpr u64 HYST_SIZE = SHARED_HYS ? (TABLE_SIZE / 2) : TABLE_SIZE;
  static constexpr u64 HYST_IDX_BITS = ta::clog2(std::max(u64(2), HYST_SIZE));

  static_assert(TABLE_SIZE >= 2 && std::has_single_bit(TABLE_SIZE),
                "TABLE_SIZE must be a power of 2 >= 2");
  static_assert(TAG_RAM_WIDTH >= TAG_WIDTH,
                "TAG_RAM_WIDTH must be >= TAG_WIDTH");

  static constexpr u64 rwram_bank_shift = RWRAM_BANK_SHIFT;

  // ---- RAMs (predict-path critical) ----
  hcm::ram<val<TAG_RAM_WIDTH>, TABLE_SIZE> tag_ram{"ta_tag"};
  ta_rwram<PRED_BITS, TABLE_SIZE, RWRAM_BANKS, RWRAM_BANK_SHIFT> pred_ram{"ta_pred"};
  hcm::ram<val<SEC_TAG_BITS>, TABLE_SIZE> sec_ram{"ta_sec"};

  // ---- Per-table folded histories (predict-path, at sec_ram location) ----
  ta_folded_gh<IDX_BITS> fold_idx;
  ta_folded_gh<TAG_WIDTH> fold_tag;

  // ---- RAMs (update-only) ----
  hcm::zone ta_update_zone;
  ta_rwram<std::max(u64(1), HYST_WIDTH), HYST_SIZE, RWRAM_BANKS, RWRAM_BANK_SHIFT> hyst_ram{"ta_hyst"};
  ta_rwram<U_WIDTH, TABLE_SIZE, RWRAM_BANKS, RWRAM_BANK_SHIFT> u_ram{"ta_u"};
};

// Generate a TATable type from config arrays at index I
template <typename Cfg, u64 I, u64 CTR_WIDTH, u64 HYST_WIDTH, u64 U_WIDTH,
          u64 SEC_TAG_BITS, u64 N, bool SHARED_HYS = false,
          u64 TAG_RAM_WIDTH = ta::array_max(Cfg::TAG_WIDTH)>
using TATableAt =
    TATable<Cfg::TABLE_SIZE[I], Cfg::TAG_WIDTH[I], Cfg::HIST_LEN[I], CTR_WIDTH,
            HYST_WIDTH, U_WIDTH, SEC_TAG_BITS, N, SHARED_HYS, TAG_RAM_WIDTH,
            Cfg::RWRAM_BANK_SHIFT[I], Cfg::RWRAM_BANKS[I]>;

// Build a reversed index sequence: <N-1, N-2, ..., 1, 0>
template <typename Seq> struct ReverseSeq;
template <u64... Is> struct ReverseSeq<std::index_sequence<Is...>> {
  static constexpr u64 N = sizeof...(Is);
  using type = std::index_sequence<(N - 1 - Is)...>;
};

// Build a tuple of TATable types
// When REVERSE=true, RAMs are declared in reverse order (T(NT-1) first)
// to influence HARCOM floorplan placement.
template <typename Cfg, u64 CTR_WIDTH, u64 HYST_WIDTH, u64 U_WIDTH,
          u64 SEC_TAG_BITS, u64 N, bool SHARED_HYS, bool REVERSE, typename Seq>
struct TAMakeTableTuple;

template <typename Cfg, u64 CTR_WIDTH, u64 HYST_WIDTH, u64 U_WIDTH,
          u64 SEC_TAG_BITS, u64 N, bool SHARED_HYS, bool REVERSE, u64... Is>
struct TAMakeTableTuple<Cfg, CTR_WIDTH, HYST_WIDTH, U_WIDTH, SEC_TAG_BITS, N,
                        SHARED_HYS, REVERSE, std::index_sequence<Is...>> {
  // Storage indices: reversed when REVERSE=true so biggest-first tables
  // get their RAMs declared first, influencing floorplan placement.
  static constexpr u64 NT = sizeof...(Is);
  static constexpr auto storage_idx(u64 i) {
    return REVERSE ? (NT - 1 - i) : i;
  }

  using type = std::tuple<
      TATableAt<Cfg, (REVERSE ? NT - 1 - Is : Is), CTR_WIDTH, HYST_WIDTH,
                U_WIDTH, SEC_TAG_BITS, N, SHARED_HYS>...>;
};

// ============================================================================
// Allocation Policy Infrastructure (shared by TageAhead, TageDirect, etc.)
// ============================================================================

// Allocation trigger: what condition causes TAGE to attempt allocation
enum class AllocTrigger {
  MISPREDICT, // final misprediction (default, conservative)
  TAGE_WRONG, // TAGE provider was wrong (even if meta corrected it)
  ALWAYS,     // every update cycle (most aggressive)
};

// Allocation action: how to gate allocation when triggered
enum class AllocAction {
  STANDARD,  // allocate in tables above provider with u=0 (default)
  FILTERED,  // probabilistically throttle by accuracy counter
  THROTTLED, // probabilistically throttle by alloc pressure counter
};

// ---- Allocation Target Policies (functors) ----
//
// Uniform interface:
//   apply(collamask, alloc_pressure, acc_pressure, rng)
// where all args are val<> types. Functors ignore what they don't need.
//
// collamask = reversed candidate mask (LSB = closest to provider)
// rng       = val<8> from caller's LFSR (hardware randomness)
//
// x & (x-1) clears the lowest set bit = skip closest candidate.

struct ClosestTarget {
  static constexpr const char *name() { return "Closest"; }
  template <u64 NT>
  static val<NT> apply(val<NT> collamask, valtype auto, valtype auto, val<8>) {
    return collamask;
  }
};

namespace target_detail {
template <u64 NT> val<NT> clear_lsb(val<NT> x) {
  x.fanout(hard<2>{});
  return x & val<NT>(x - 1);
}

template <u64 SKIP, u64 NT> val<NT> skip_n(val<NT> x) {
  static_assert(SKIP <= 4, "skip_n: SKIP > 4 not supported");
  if constexpr (SKIP == 0)
    return x;
  else if constexpr (SKIP == 1)
    return clear_lsb(x);
  else {
    val<NT> s = skip_n<SKIP - 1, NT>(x);
    return clear_lsb(s);
  }
}
} // namespace target_detail

// Deterministically skip SKIP closest candidates, always allocate further out
template <u64 SKIP = 1> struct DeterministicSkipTarget {
  static_assert(SKIP <= 4, "DeterministicSkipTarget: SKIP > 4 not supported");
  static constexpr const char *name() { return "DetSkip"; }
  template <u64 NT>
  static val<NT> apply(val<NT> collamask, valtype auto, valtype auto, val<8>) {
    return target_detail::skip_n<SKIP, NT>(collamask);
  }
};

// Skip SKIP closest with static probability PROB/256
template <u64 SKIP = 1, u64 PROB_256 = 64> struct StaticSkipTarget {
  static_assert(SKIP <= 4, "StaticSkipTarget: SKIP > 4 not supported");
  static constexpr const char *name() { return "StaticSkip"; }
  template <u64 NT>
  static val<NT> apply(val<NT> collamask, valtype auto, valtype auto,
                       val<8> rng) {
    collamask.fanout(hard<2>{});
    val<NT> skipped = target_detail::skip_n<SKIP, NT>(collamask);
    return select(rng < hard<PROB_256>{}, skipped, collamask);
  }
};

// Skip probability scales with alloc pressure (high pressure = skip more)
template <u64 SKIP = 1> struct AllocPressureSkipTarget {
  static_assert(SKIP <= 4, "AllocPressureSkipTarget: SKIP > 4 not supported");
  static constexpr const char *name() { return "AllocPressSkip"; }
  template <u64 NT>
  static val<NT> apply(val<NT> collamask, valtype auto ap, valtype auto,
                       val<8> rng) {
    collamask.fanout(hard<2>{});
    val<NT> skipped = target_detail::skip_n<SKIP, NT>(collamask);
    return select(val<8>{ap} > rng, skipped, collamask);
  }
};

// Skip probability scales with accuracy pressure (low accuracy = skip more)
template <u64 SKIP = 1> struct AccuracyPressureSkipTarget {
  static_assert(SKIP <= 4,
                "AccuracyPressureSkipTarget: SKIP > 4 not supported");
  static constexpr const char *name() { return "AccPressSkip"; }
  template <u64 NT>
  static val<NT> apply(val<NT> collamask, valtype auto, valtype auto acp,
                       val<8> rng) {
    collamask.fanout(hard<2>{});
    val<NT> skipped = target_detail::skip_n<SKIP, NT>(collamask);
    return select(val<8>{acp} > rng, skipped, collamask);
  }
};

// ---- Allocation Config Structs ----

struct TADefaultAllocConfig {
  static constexpr u64 MAX_ALLOC = 1;
  static constexpr bool NON_CONSECUTIVE = false;
  static constexpr AllocTrigger ALLOC_TRIGGER = AllocTrigger::MISPREDICT;
  static constexpr AllocAction ALLOC_ACTION = AllocAction::STANDARD;
  using TARGET_POLICY = ClosestTarget;
};

struct TAAllocSkip1 : TADefaultAllocConfig {
  using TARGET_POLICY = StaticSkipTarget<1, 64>;
};

struct TAAllocDetSkip1 : TADefaultAllocConfig {
  using TARGET_POLICY = DeterministicSkipTarget<1>;
};

struct TAAlloc2 : TADefaultAllocConfig {
  static constexpr u64 MAX_ALLOC = 2;
};

struct TAAllocTageWrong : TADefaultAllocConfig {
  static constexpr AllocTrigger ALLOC_TRIGGER = AllocTrigger::TAGE_WRONG;
};

struct TAAllocFiltered : TADefaultAllocConfig {
  static constexpr AllocAction ALLOC_ACTION = AllocAction::FILTERED;
};

struct TAAllocThrottled : TADefaultAllocConfig {
  static constexpr AllocAction ALLOC_ACTION = AllocAction::THROTTLED;
};

struct TAAlloc2PressSkip : TADefaultAllocConfig {
  static constexpr u64 MAX_ALLOC = 2;
  using TARGET_POLICY = AllocPressureSkipTarget<1>;
};

struct TAAllocPressSkip : TADefaultAllocConfig {
  using TARGET_POLICY = AllocPressureSkipTarget<1>;
};

#endif // CUSTOM_COMMON_H
