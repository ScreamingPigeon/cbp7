#ifndef CUSTOM_COMMON_H
#define CUSTOM_COMMON_H

#include "../../harcom.hpp"
#include "../common.hpp"

using namespace hcm;

// ============================================================================
// History update mode for geometric_folds_ex
// ============================================================================

enum class HistUpdate { PATH, DIR, BOTH };

// ============================================================================
// geometric_folds_ex — geometric_folds with templated update mode
//
// Extends geometric_folds with three update modes:
//   update<PATH>(path_bits)            — fold in path bits only
//   update<DIR>(dir_bit)               — fold in direction bit only
//   update<BOTH>(dir_bit, path_bits)   — fold in concat(dir, path)
// ============================================================================

template<u64 NH, u64 MINH, u64 MAXH, u64... FOLDS>
struct geometric_folds_ex {
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

  global_history<MAXH> gh;
  std::array<std::tuple<folded_gh<FOLDS>...>, NH> folds;

  template<u64 J = 0>
  auto get(u64 i) {
    if (i >= NH) {
      std::cerr << "geometric_folds_ex: out of bound access\n";
      std::terminate();
    }
    return std::get<J>(folds[i]).get();
  }

  void fanout(hardval auto fo) {
    for (u64 i = 0; i < NH; i++) {
      static_loop<NF>([&]<u64 J>() {
        std::get<J>(folds[i]).fanout(fo);
      });
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
  template<HistUpdate M>
  void update(valtype auto path_bits)
    requires (M == HistUpdate::PATH)
  {
    do_update(path_bits);
  }

  // update<DIR>(dir_bit) — direction only (1-bit)
  template<HistUpdate M>
  void update(val<1> dir_bit)
    requires (M == HistUpdate::DIR)
  {
    do_update(dir_bit);
  }

  // update<BOTH>(dir_bit, path_bits) — concat(dir, path)
  template<HistUpdate M>
  void update(val<1> dir_bit, valtype auto path_bits)
    requires (M == HistUpdate::BOTH)
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
  for (std::size_t i = 0; i < N; i++) a[i] = v;
  return a;
}

constexpr u64 clog2(u64 x) {
  u64 r = 0, v = x - 1;
  while (v > 0) { v >>= 1; r++; }
  return r;
}

template <typename T, std::size_t N>
constexpr T array_max(const std::array<T, N> &a) {
  T m = a[0];
  for (std::size_t i = 1; i < N; i++) m = (a[i] > m) ? a[i] : m;
  return m;
}

template <typename T, std::size_t N>
constexpr T array_min(const std::array<T, N> &a) {
  T m = a[0];
  for (std::size_t i = 1; i < N; i++) m = (a[i] < m) ? a[i] : m;
  return m;
}

constexpr double constexpr_pow(double base, double exp) {
  if (exp == 0.0) return 1.0;
  if (base == 0.0) return 0.0;
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
    u64 hl = static_cast<u64>(static_cast<double>(min_h) * constexpr_pow(ratio, e));
    hl = (hl > prev + 1) ? hl : prev + 1;
    h[N - 1 - i] = hl;
    prev = hl;
  }
  return h;
}

template <std::size_t N>
constexpr std::array<u64, N> quadratic_hist(u64 minh, u64 maxh, u64 d = 2, u64 k = 1) {
  std::array<u64, N> h{};
  h[0] = minh;
  for (std::size_t n = 1; n < N; n++)
    h[n] = h[n - 1] + d * n + k;
  double raw_max = h[N - 1], raw_min = h[0];
  u64 prev = 0;
  for (std::size_t n = 0; n < N; n++) {
    double t = (raw_max > raw_min) ? (double(h[n]) - raw_min) / (raw_max - raw_min) : 0.0;
    u64 scaled = minh + u64(t * (maxh - minh));
    scaled = (scaled > prev + 1 || n == 0) ? scaled : prev + 1;
    h[n] = scaled; prev = scaled;
  }
  for (std::size_t i = 0; i < N / 2; i++) std::swap(h[i], h[N - 1 - i]);
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
  for (std::size_t i = 0; i < N / 2; i++) std::swap(h[i], h[N - 1 - i]);
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
    double t_val = (raw_max > raw_min) ? (double(h[n]) - raw_min) / (raw_max - raw_min) : 0.0;
    u64 scaled = minh + u64(t_val * (maxh - minh));
    scaled = (scaled > prev + 1 || n == 0) ? scaled : prev + 1;
    h[n] = scaled; prev = scaled;
  }
  for (std::size_t i = 0; i < N / 2; i++) std::swap(h[i], h[N - 1 - i]);
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
    if (n <= 1) return TAG_LONG;
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
  for (u64 i = 0; i < N; i++) t[i] = fn(i, N);
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
    if (RATIO <= 1 || n <= 1) return SIZE;
    double t = double(i) / double(n - 1);
    double scale = constexpr_pow(double(RATIO), t - 0.5);
    u64 sz = u64(SIZE * scale);
    u64 result = 64;
    while (result < sz) result *= 2;
    return result;
  }
};

template <u64 SIZE, u64 RATIO> struct InvGeoSize {
  constexpr u64 operator()(u64 i, u64 n) const {
    if (RATIO <= 1 || n <= 1) return SIZE;
    double t = double(n - 1 - i) / double(n - 1);
    double scale = constexpr_pow(double(RATIO), t - 0.5);
    u64 sz = u64(SIZE * scale);
    u64 result = 64;
    while (result < sz) result *= 2;
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
    if (n <= 1) return BASE_SIZE;
    double scale = constexpr_pow(double(n - 1 - i + 1), 0.5);
    u64 sz = u64(BASE_SIZE * scale);
    u64 result = 64;
    while (result < sz) result *= 2;
    return result;
  }
};

template <u64 TOTAL_BITS, u64 TAG, u64 CTR, u64 HYST, u64 U>
struct ConstBitsSize {
  constexpr u64 operator()(u64, u64) const {
    u64 entry_bits = TAG + CTR + HYST + U;
    u64 entries = TOTAL_BITS / entry_bits;
    u64 result = 64;
    while (result * 2 <= entries) result *= 2;
    return result;
  }
};

template <std::size_t N, typename SizeFn>
constexpr std::array<u64, N> generate_table_sizes(SizeFn fn) {
  std::array<u64, N> s{};
  for (std::size_t i = 0; i < N; i++) s[i] = fn(i, N);
  return s;
}

} // namespace ta

// ============================================================================
// TATable — Per-table storage for ahead-pipelined TAGE
//
// Each table has its own compile-time TAG_WIDTH and TABLE_SIZE.
// Pred RAMs are per-bank (BANKS copies). Hyst/u are update-only (plain ram).
// Pipeline regs (ahead[0]/[1]) are included per table.
// ============================================================================

template <u64 TABLE_SIZE,
          u64 TAG_WIDTH,
          u64 HIST_LEN,     // per-table history length
          u64 CTR_WIDTH,    // per-lane prediction counter
          u64 HYST_WIDTH,
          u64 U_WIDTH,
          u64 SEC_TAG_BITS,
          u64 N>            // max branches per block (= lanes of pred)
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

  static_assert(TABLE_SIZE >= 2 && std::has_single_bit(TABLE_SIZE),
                "TABLE_SIZE must be a power of 2 >= 2");

  // ---- RAMs ----
  hcm::ram<val<TAG_WIDTH>, TABLE_SIZE>    tag_ram{"ta_tag"};
  hcm::ram<val<PRED_BITS>, TABLE_SIZE>    pred_ram{"ta_pred"};
  hcm::ram<val<SEC_TAG_BITS>, TABLE_SIZE> sec_ram{"ta_sec"};
  hcm::ram<val<std::max(u64(1), HYST_WIDTH)>, TABLE_SIZE> hyst_ram{"ta_hyst"};
  hcm::ram<val<U_WIDTH>, TABLE_SIZE>      u_ram{"ta_u"};

  // ---- Per-table folded histories (fold into exact widths) ----
  folded_gh<IDX_BITS> fold_idx;
  folded_gh<TAG_WIDTH> fold_tag;

};

// Generate a TATable type from config arrays at index I
template <typename Cfg, u64 I, u64 CTR_WIDTH, u64 HYST_WIDTH, u64 U_WIDTH,
          u64 SEC_TAG_BITS, u64 N>
using TATableAt =
    TATable<Cfg::TABLE_SIZE[I], Cfg::TAG_WIDTH[I], Cfg::HIST_LEN[I],
            CTR_WIDTH, HYST_WIDTH, U_WIDTH, SEC_TAG_BITS, N>;

// Build a tuple of TATable types
template <typename Cfg, u64 CTR_WIDTH, u64 HYST_WIDTH, u64 U_WIDTH,
          u64 SEC_TAG_BITS, u64 N, typename Seq>
struct TAMakeTableTuple;

template <typename Cfg, u64 CTR_WIDTH, u64 HYST_WIDTH, u64 U_WIDTH,
          u64 SEC_TAG_BITS, u64 N, u64... Is>
struct TAMakeTableTuple<Cfg, CTR_WIDTH, HYST_WIDTH, U_WIDTH,
                        SEC_TAG_BITS, N, std::index_sequence<Is...>> {
  using type = std::tuple<TATableAt<Cfg, Is, CTR_WIDTH, HYST_WIDTH, U_WIDTH,
                                    SEC_TAG_BITS, N>...>;
};

#endif // CUSTOM_COMMON_H
