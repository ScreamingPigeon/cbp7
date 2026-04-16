#ifndef CUSTOM_COMMON_H
#define CUSTOM_COMMON_H

#include "../../harcom.hpp"

using namespace hcm;

// ============================================================================
// update_ctr — saturating up-down counter
// ============================================================================

template<u64 N, typename T>
[[nodiscard]] val<N,T> ta_update_ctr(val<N,T> ctr, val<1> incr)
{
    ctr.fanout(hard<6>{});
    val<N,T> incsat = select(ctr==hard<ctr.maxval>{},ctr,val<N,T>{ctr+1});
    val<N,T> decsat = select(ctr==hard<ctr.minval>{},ctr,val<N,T>{ctr-1});
    return select(incr.fo1(),incsat.fo1(),decsat.fo1());
}

// ============================================================================
// ta_global_history — shift register updated by XOR with branch bits
// ============================================================================

template<u64 N>
struct ta_global_history {
    arr<reg<1>,N> h;

    void update(valtype auto in)
    {
        auto input = in.fo1().make_array(val<1>{});
        static_assert(input.size<=N);
        for (u64 i=N-1; i>=input.size; i--) h[i] = h[i-1];
        for (u64 i=input.size-1; i>=1; i--) h[i] = h[i-1] ^ input[i].fo1();
        h[0] = input[0].fo1();
    }

    // Unconditional write with enable mux — avoids execute_if gate timing.
    // When enable=1: shift and XOR in new input (normal update).
    // When enable=0: hold current state (registers still written, mux selects old value).
    // Explicit fo1() on h[i] breaks the combinational feedback loop
    // (h[i] read → select → h[i] write) that inflates timing.
    void update(valtype auto in, val<1> enable)
    {
        auto input = in.fo1().make_array(val<1>{});
        static_assert(input.size<=N);
        enable.fanout(hard<N>{});
        for (u64 i=N-1; i>=input.size; i--) {
            val<1> old = h[i].fo1();
            h[i] = select(enable, h[i-1], old);
        }
        for (u64 i=input.size-1; i>=1; i--) {
            val<1> old = h[i].fo1();
            h[i] = select(enable, h[i-1] ^ input[i].fo1(), old);
        }
        val<1> old0 = h[0].fo1();
        h[0] = select(enable, input[0].fo1(), old0);
    }

    val<1>& operator[] (u64 i) { return h[i]; }
    void fanout(hardval auto fo) { h.fanout(fo); }

    // Per-bit fanout: set each bit's fanout from a constexpr array
    template <auto const &FO_ARRAY>
    void fanout_per_bit() {
      static_assert(FO_ARRAY.size() == N);
      static_loop<N>([&]<u64 I>() {
        h[I].fanout(hard<FO_ARRAY[I]>{});
      });
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

template<u64 F>
struct ta_folded_gh {
    static_assert(F!=0);

    reg<F> folded;

    val<F> get() { return folded; }
    void fanout(hardval auto fo) { folded.fanout(fo); }

    // Compute new fold value without writing the register
    template<u64 MAXL>
    [[nodiscard]] inline val<F> compute_update(ta_global_history<MAXL> &gh, hardval auto ghlen, valtype auto in)
    {
        constexpr u64 inbits = std::min(F,std::min(in.size,ghlen.value));
        val<inbits> input = in.fo1();
        auto f = folded.make_array(val<1>{});
        static_assert(f.size==F);
        val<1> outbit = gh[ghlen-1];
        u64 outpos = ghlen % F;
        arr<val<1>,F> ff = [&](u64 i){
            if (i==0) return (outpos==0)? f[F-1].fo1()^outbit.fo1() : f[F-1].fo1();
            else return (outpos==i)? f[i-1].fo1()^outbit.fo1() : f[i-1].fo1();
        };
        auto x = input.fo1().make_array(val<1>{});
        arr<val<1>,F> y = [&](u64 i){return (i<x.size)? x[i].fo1()^ff[i].fo1() : ff[i].fo1();};
        return y.fo1().concat();
    }

    // Write a precomputed fold value into the register
    void apply_update(val<F> new_val)
    {
        folded = new_val;
    }

    // Unconditional write with enable mux — avoids execute_if gate timing.
    // When enable=1: write new_val. When enable=0: hold current value.
    // Explicit fo1() on folded breaks the combinational feedback loop
    // (folded read → select → folded write) that inflates timing.
    void apply_update(val<F> new_val, val<1> enable)
    {
        val<F> old_val = folded.fo1();
        folded = select(enable, new_val, old_val);
    }

    // Combined compute + apply (convenience, same as original folded_gh::update)
    template<u64 MAXL>
    void update(ta_global_history<MAXL> &gh, hardval auto ghlen, valtype auto in)
    {
        apply_update(compute_update(gh, ghlen, in));
    }
};

// ============================================================================
// History update mode for geometric_folds_ex
// ============================================================================

enum class HistUpdate { PATH, DIR, BOTH };
enum class DecayMiss { TAG, SEC, TAG_OR_SEC, TAG_AND_SEC };
enum class DecayOp { DECREMENT, HALVE, CLEAR };
enum class UProvUpdate {
  SET_OR_CLEAR, // current: correct → max_u, wrong → 0
  SET_ON_CORRECT, // Tage.hpp style: correct → max_u, wrong → no write
  INC_DEC,        // saturating: correct → u+1, wrong → u-1
  INC_ONLY        // correct → u+1, wrong → no write
};

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

  ta_global_history<MAXH> gh;
  std::array<std::tuple<ta_folded_gh<FOLDS>...>, NH> folds;

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

// Compute per-bit gh fanout: for each bit position 0..MAXHIST-1,
// count how many tables need gh[HIST_LEN[I]-1] at that position.
// Each table reads the outgoing bit twice (fold_idx + fold_tag).
// Add 1 for the gh.update() shift read.
template <std::size_t MAXHIST, std::size_t NT, auto const &HIST_LEN>
constexpr std::array<u64, MAXHIST> gh_per_bit_fanout() {
  std::array<u64, MAXHIST> fo{};
  // Base: each bit is read by update() for shift (h[i-1]) and mux hold (h[i])
  for (std::size_t i = 0; i < MAXHIST; i++) fo[i] = 3;
  for (std::size_t t = 0; t < NT; t++) {
    u64 bit = HIST_LEN[t] - 1;
    fo[bit] += 2; // fold_idx + fold_tag read the outgoing bit
  }
  return fo;
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

// ============================================================================
// Decay helpers
// ============================================================================

// Generate a tuple of reg<W> with per-element widths from a constexpr array.
// Usage: TALfsrTuple<LFSR_WIDTHS, std::make_index_sequence<NT>>::type
template <auto const &Widths, typename Seq>
struct TALfsrTuple;

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

// Epoch trigger functors.
// Interface: should_fire<AW, PW, EW>(acc_ctr, alloc_ctr, epoch_ctr) → val<1>
// acc_ctr: accuracy counter (high = predicting well)
// alloc_ctr: allocation pressure counter (high = can't allocate)
// epoch_ctr: free-running counter incremented each update_cycle

// AllocSaturate: fire when alloc_ctr saturates (original behavior).
struct AllocSaturateEpoch {
  template <u64 AW, u64 PW, u64 EW>
  static val<1> should_fire(val<AW> /*acc_ctr*/, val<PW> alloc_ctr,
                             val<EW> /*epoch_ctr*/) {
    return alloc_ctr == hard<alloc_ctr.maxval>{};
  }
};
using DefaultEpochTrigger = AllocSaturateEpoch;

// FixedIntervalEpoch<PERIOD>: fire every PERIOD update_cycles.
// PERIOD must be a power of 2 (uses bit test for zero-gate-cost).
template <u64 PERIOD>
struct FixedIntervalEpoch {
  static_assert(PERIOD > 0 && (PERIOD & (PERIOD - 1)) == 0,
                "PERIOD must be a power of 2");
  template <u64 AW, u64 PW, u64 EW>
  static val<1> should_fire(val<AW> /*acc_ctr*/, val<PW> /*alloc_ctr*/,
                             val<EW> epoch_ctr) {
    // Fire when low bits are all zero (i.e. counter is a multiple of PERIOD)
    constexpr u64 MASK = PERIOD - 1;
    return (epoch_ctr & hard<MASK>{}) == hard<0>{};
  }
};

// AllocAccJointEpoch<ALLOC_THRESH, ACC_THRESH>: fire when alloc pressure
// exceeds ALLOC_THRESH AND accuracy is below ACC_THRESH.
// Avoids resetting when predictor is accurate despite high pressure.
template <u64 ALLOC_THRESH, u64 ACC_THRESH>
struct AllocAccJointEpoch {
  template <u64 AW, u64 PW, u64 EW>
  static val<1> should_fire(val<AW> acc_ctr, val<PW> alloc_ctr,
                             val<EW> /*epoch_ctr*/) {
    return (alloc_ctr >= hard<ALLOC_THRESH>{}) &
           (acc_ctr <= hard<ACC_THRESH>{});
  }
};

// CountdownEpoch<PERIOD, ACC_GATE>: fire every PERIOD update_cycles, but
// only if acc_ctr < ACC_GATE. Skips reset when predictor is doing well.
// PERIOD must be a power of 2.
template <u64 PERIOD, u64 ACC_GATE>
struct CountdownEpoch {
  static_assert(PERIOD > 0 && (PERIOD & (PERIOD - 1)) == 0,
                "PERIOD must be a power of 2");
  template <u64 AW, u64 PW, u64 EW>
  static val<1> should_fire(val<AW> acc_ctr, val<PW> /*alloc_ctr*/,
                             val<EW> epoch_ctr) {
    constexpr u64 MASK = PERIOD - 1;
    val<1> interval_hit = (epoch_ctr & hard<MASK>{}) == hard<0>{};
    return interval_hit & (acc_ctr < hard<ACC_GATE>{});
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

template <u64 N, u64 M, u64 B, u64 BANK_SHIFT = 0>
struct ta_rwram {
  static_assert(std::has_single_bit(M));
  static_assert(B >= 2 && B <= 64);
  static_assert(std::has_single_bit(B));
  static constexpr u64 A = std::bit_width(M - 1);         // address bits
  static constexpr u64 E = M / B;                          // entries per bank
  static_assert(E > 1);
  static constexpr u64 L = std::bit_width(E - 1);          // local address bits
  static constexpr u64 BANK_BITS = std::bit_width(B - 1);  // bank ID bits
  static_assert(A == L + BANK_BITS);
  static_assert(BANK_SHIFT + BANK_BITS <= A);

  hcm::ram<val<N>, E> bank[B];
  reg<B> read_bank;

  // buffered write state
  reg<B> write_bank;
  reg<L> write_localaddr;
  reg<N> write_data;

  ta_rwram(const char *label = "") : bank{label} {}

  // Split full address into (local_addr, bank_id).
  auto split_addr(val<A> addr) {
    if constexpr (BANK_SHIFT == 0) {
      return hcm::split<L, BANK_BITS>(addr);
    } else if constexpr (BANK_SHIFT + BANK_BITS == A) {
      return std::pair{val<L>{addr}, val<BANK_BITS>{addr >> BANK_SHIFT}};
    } else {
      val<BANK_SHIFT> lo = addr;
      val<BANK_BITS> bankid = addr >> BANK_SHIFT;
      val<A - BANK_SHIFT - BANK_BITS> hi = addr >> (BANK_SHIFT + BANK_BITS);
      val<L> localaddr = concat(lo, hi);
      return std::pair{localaddr, bankid};
    }
  }

  val<N> read(val<A> addr) {
    auto [localaddr, bankid] = split_addr(addr.fo1());
    localaddr.fanout(hard<B>{});
    arr<val<1>, B> banksel = bankid.fo1().decode();
    banksel.fanout(hard<2>{});
    arr<val<N>, B> data = [&](u64 i) -> val<N> {
      return execute_if(banksel[i], [&]() { return bank[i].read(localaddr); });
    };
    read_bank = banksel.concat();
    return data.fo1().fold_or();
  }

  void write(val<A> addr, val<N> data, val<1> noconflict) {
    // noconflict=1: no read this cycle, write immediately.
    // noconflict=0: buffer write, flush when bank is free.
    auto [localaddr, bankid] = split_addr(addr.fo1());
    data.fanout(hard<B + 1>{});
    noconflict.fanout(hard<B + 2>{});
    val<B> banksel = bankid.fo1().decode().concat();
    banksel.fanout(hard<2>{});
    val<B> noconflict_mask = noconflict.replicate(hard<B>{}).concat();
    noconflict_mask.fanout(hard<2>{});
    val<B> current_write = banksel & noconflict_mask;
    current_write.fanout(hard<3>{});
    arr<val<1>, B> current_write_split = current_write.make_array(val<1>{});
    current_write_split.fanout(hard<3>{});
    arr<val<1>, B> write_bank_split = write_bank.make_array(val<1>{});
    arr<val<1>, B> read_bank_split = read_bank.make_array(val<1>{});
    for (u64 i = 0; i < B; i++) {
      execute_if(
          current_write_split[i] |
              (write_bank_split[i].fo1() & ~read_bank_split[i].fo1()),
          [&]() {
            val<L> a =
                select(current_write_split[i], localaddr, write_localaddr);
            val<N> d = select(current_write_split[i], data, write_data);
            bank[i].write(a.fo1(), d.fo1());
          });
    }
    // buffer the current write if not done
    val<1> buffered_done =
        (write_bank & (current_write | read_bank)) == hard<0>{};
    execute_if(buffered_done.fo1() | ~noconflict, [&]() {
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
};

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
          u64 N,            // max branches per block (= lanes of pred)
          bool SHARED_HYS = false,  // shared hyst: 2 entries share 1 counter
          u64 HYST_BANKS = 1,       // hyst banks (1=shared, N=per-branch)
          u64 U_BANKS = 1>          // u banks (1=shared, N=per-branch)
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
  static constexpr u64 hyst_banks = HYST_BANKS;
  static constexpr u64 u_banks = U_BANKS;

  // When SHARED_HYS=true, halve the hyst table: pairs of entries share one counter
  static constexpr u64 HYST_SIZE = SHARED_HYS ? (TABLE_SIZE / 2) : TABLE_SIZE;
  static constexpr u64 HYST_IDX_BITS = ta::clog2(std::max(u64(2), HYST_SIZE));

  static_assert(TABLE_SIZE >= 2 && std::has_single_bit(TABLE_SIZE),
                "TABLE_SIZE must be a power of 2 >= 2");
  static_assert(HYST_BANKS >= 1 && std::has_single_bit(HYST_BANKS),
                "HYST_BANKS must be a power of 2");
  static_assert(U_BANKS >= 1 && std::has_single_bit(U_BANKS),
                "U_BANKS must be a power of 2");

  // ---- RAMs ----
  hcm::ram<val<TAG_WIDTH>, TABLE_SIZE>    tag_ram{"ta_tag"};
  ta_rwram<PRED_BITS, TABLE_SIZE, 2>      pred_ram{"ta_pred"};
  ta_rwram<std::max(u64(1), HYST_WIDTH), HYST_SIZE, 2> hyst_ram[HYST_BANKS];
  ta_rwram<U_WIDTH, TABLE_SIZE, 2>        u_ram[U_BANKS];

  // ---- Per-table folded histories (fold into exact widths) ----
  ta_folded_gh<IDX_BITS> fold_idx;
  ta_folded_gh<TAG_WIDTH> fold_tag;

};

// Generate a TATable type from config arrays at index I
template <typename Cfg, u64 I, u64 CTR_WIDTH, u64 HYST_WIDTH, u64 U_WIDTH,
          u64 SEC_TAG_BITS, u64 N, bool SHARED_HYS = false,
          u64 HYST_BANKS = 1, u64 U_BANKS = 1>
using TATableAt =
    TATable<Cfg::TABLE_SIZE[I], Cfg::TAG_WIDTH[I], Cfg::HIST_LEN[I],
            CTR_WIDTH, HYST_WIDTH, U_WIDTH, SEC_TAG_BITS, N, SHARED_HYS,
            HYST_BANKS, U_BANKS>;

// Build a tuple of TATable types
template <typename Cfg, u64 CTR_WIDTH, u64 HYST_WIDTH, u64 U_WIDTH,
          u64 SEC_TAG_BITS, u64 N, bool SHARED_HYS,
          u64 HYST_BANKS, u64 U_BANKS, typename Seq>
struct TAMakeTableTuple;

template <typename Cfg, u64 CTR_WIDTH, u64 HYST_WIDTH, u64 U_WIDTH,
          u64 SEC_TAG_BITS, u64 N, bool SHARED_HYS,
          u64 HYST_BANKS, u64 U_BANKS, u64... Is>
struct TAMakeTableTuple<Cfg, CTR_WIDTH, HYST_WIDTH, U_WIDTH,
                        SEC_TAG_BITS, N, SHARED_HYS,
                        HYST_BANKS, U_BANKS, std::index_sequence<Is...>> {
  using type = std::tuple<TATableAt<Cfg, Is, CTR_WIDTH, HYST_WIDTH, U_WIDTH,
                                    SEC_TAG_BITS, N, SHARED_HYS,
                                    HYST_BANKS, U_BANKS>...>;
};

// ============================================================================
// Allocation Policy Infrastructure (shared by TageAhead, TageDirect, etc.)
// ============================================================================

// Allocation trigger: what condition causes TAGE to attempt allocation
enum class AllocTrigger {
  MISPREDICT,   // final misprediction (default, conservative)
  TAGE_WRONG,   // TAGE provider was wrong (even if meta corrected it)
  ALWAYS,       // every update cycle (most aggressive)
};

// Allocation action: how to gate allocation when triggered
enum class AllocAction {
  STANDARD,     // allocate in tables above provider with u=0 (default)
  FILTERED,     // probabilistically throttle by accuracy counter
  THROTTLED,    // probabilistically throttle by alloc pressure counter
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
  static constexpr const char* name() { return "Closest"; }
  template <u64 NT>
  static val<NT> apply(val<NT> collamask, valtype auto, valtype auto, val<8>) {
    return collamask;
  }
};

namespace target_detail {
  template <u64 NT>
  val<NT> clear_lsb(val<NT> x) {
    x.fanout(hard<2>{});
    return x & val<NT>(x - 1);
  }

  template <u64 SKIP, u64 NT>
  val<NT> skip_n(val<NT> x) {
    static_assert(SKIP <= 4, "skip_n: SKIP > 4 not supported");
    if constexpr (SKIP == 0) return x;
    else if constexpr (SKIP == 1) return clear_lsb(x);
    else { val<NT> s = skip_n<SKIP - 1, NT>(x); return clear_lsb(s); }
  }
} // namespace target_detail

// Deterministically skip SKIP closest candidates, always allocate further out
template <u64 SKIP = 1>
struct DeterministicSkipTarget {
  static_assert(SKIP <= 4, "DeterministicSkipTarget: SKIP > 4 not supported");
  static constexpr const char* name() { return "DetSkip"; }
  template <u64 NT>
  static val<NT> apply(val<NT> collamask, valtype auto, valtype auto, val<8>) {
    return target_detail::skip_n<SKIP, NT>(collamask);
  }
};

// Skip SKIP closest with static probability PROB/256
template <u64 SKIP = 1, u64 PROB_256 = 64>
struct StaticSkipTarget {
  static_assert(SKIP <= 4, "StaticSkipTarget: SKIP > 4 not supported");
  static constexpr const char* name() { return "StaticSkip"; }
  template <u64 NT>
  static val<NT> apply(val<NT> collamask, valtype auto, valtype auto, val<8> rng) {
    collamask.fanout(hard<2>{});
    val<NT> skipped = target_detail::skip_n<SKIP, NT>(collamask);
    return select(rng < hard<PROB_256>{}, skipped, collamask);
  }
};

// Skip probability scales with alloc pressure (high pressure = skip more)
template <u64 SKIP = 1>
struct AllocPressureSkipTarget {
  static_assert(SKIP <= 4, "AllocPressureSkipTarget: SKIP > 4 not supported");
  static constexpr const char* name() { return "AllocPressSkip"; }
  template <u64 NT>
  static val<NT> apply(val<NT> collamask, valtype auto ap, valtype auto, val<8> rng) {
    collamask.fanout(hard<2>{});
    val<NT> skipped = target_detail::skip_n<SKIP, NT>(collamask);
    return select(val<8>{ap} > rng, skipped, collamask);
  }
};

// Skip probability scales with accuracy pressure (low accuracy = skip more)
template <u64 SKIP = 1>
struct AccuracyPressureSkipTarget {
  static_assert(SKIP <= 4, "AccuracyPressureSkipTarget: SKIP > 4 not supported");
  static constexpr const char* name() { return "AccPressSkip"; }
  template <u64 NT>
  static val<NT> apply(val<NT> collamask, valtype auto, valtype auto acp, val<8> rng) {
    collamask.fanout(hard<2>{});
    val<NT> skipped = target_detail::skip_n<SKIP, NT>(collamask);
    return select(val<8>{acp} > rng, skipped, collamask);
  }
};

// ---- Replacement Policy Functors (Technique 4) ----
// Determine whether an entry is replaceable during allocation.
// Signature: is_replaceable(u, hyst, alloc_p, acc_p) → val<1>
// u = val<UW>, hyst = val<HW>, alloc_p/acc_p = valtype auto (pressure counters)

// Original behavior: replace when u == 0 (ignore confidence)
struct ReplaceUZero {
  static constexpr const char* name() { return "UZero"; }
  template <u64 UW, u64 HW>
  static val<1> is_replaceable(val<UW> u, val<HW>, valtype auto, valtype auto) {
    return u == hard<0>{};
  }
};

// Replace when u == 0 AND hysteresis below threshold (low confidence)
template <u64 THRESH = 1>
struct ReplaceUZeroWeakConf {
  static constexpr const char* name() { return "UZeroWeak"; }
  template <u64 UW, u64 HW>
  static val<1> is_replaceable(val<UW> u, val<HW> hyst, valtype auto, valtype auto) {
    return (u == hard<0>{}) & (val<HW>{hyst} <= val<HW>{hard<THRESH>{}});
  }
};

// Replace when u == 0, but relax confidence gate under high alloc pressure
template <u64 THRESH = 1, u64 ALLOC_W = 4>
struct ReplacePressureAdaptive {
  static constexpr const char* name() { return "PressAdapt"; }
  template <u64 UW, u64 HW>
  static val<1> is_replaceable(val<UW> u, val<HW> hyst, val<ALLOC_W> ap, valtype auto) {
    val<1> u_zero = (u == hard<0>{});
    val<1> weak = (val<HW>{hyst} <= val<HW>{hard<THRESH>{}});
    // Under high pressure (alloc_ctr near saturation), allow any u==0 entry
    val<1> high_pressure = val<1>{ap >> hard<ALLOC_W - 1>{}};
    return u_zero & (weak | high_pressure);
  }
};

// ---- Alt Bank Promotion Functors (Technique 3) ----
// When provider wrong AND alt correct, optionally set alt's u-bit to protect it.
// Signature: should_promote(prov_wrong, alt_correct, alloc_p, acc_p, rng) → val<1>

struct AltPromoteOff {
  static constexpr const char* name() { return "AltPromOff"; }
  static val<1> should_promote(val<1>, val<1>, valtype auto, valtype auto, val<8>) {
    return val<1>{0};
  }
};

struct AltPromoteAlways {
  static constexpr const char* name() { return "AltPromAlways"; }
  static val<1> should_promote(val<1> prov_wrong, val<1> alt_correct, valtype auto, valtype auto, val<8>) {
    return prov_wrong & alt_correct;
  }
};

template <u64 PROB_256 = 128>
struct AltPromoteProb {
  static constexpr const char* name() { return "AltPromProb"; }
  static val<1> should_promote(val<1> prov_wrong, val<1> alt_correct, valtype auto, valtype auto, val<8> rng) {
    return prov_wrong & alt_correct & (rng < hard<PROB_256>{});
  }
};

// Promote more aggressively under high alloc pressure (entries being evicted fast)
struct AltPromotePressure {
  static constexpr const char* name() { return "AltPromPress"; }
  static val<1> should_promote(val<1> prov_wrong, val<1> alt_correct, valtype auto ap, valtype auto, val<8> rng) {
    return prov_wrong & alt_correct & (val<8>{ap} > rng);
  }
};

// ---- Allocation Config Structs ----

struct TADefaultAllocConfig {
  static constexpr u64 MAX_ALLOC = 1;
  static constexpr bool NON_CONSECUTIVE = false;
  static constexpr u64 ALLOC_DECAY_SHIFT = 1; // prob halved per extra slot (L-TAGE: 1)
  static constexpr AllocTrigger ALLOC_TRIGGER = AllocTrigger::MISPREDICT;
  static constexpr AllocAction  ALLOC_ACTION  = AllocAction::STANDARD;
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

struct TAAlloc3 : TADefaultAllocConfig {
  static constexpr u64 MAX_ALLOC = 3;
};

// L-TAGE style: up to 3 slots, prob decay 1/2 per extra slot
struct TAAllocLTAGE : TADefaultAllocConfig {
  static constexpr u64 MAX_ALLOC = 3;
  static constexpr u64 ALLOC_DECAY_SHIFT = 1;
};

// Aggressive: up to 3 slots, no decay (always allocate all 3 if available)
struct TAAlloc3NoDec : TADefaultAllocConfig {
  static constexpr u64 MAX_ALLOC = 3;
  static constexpr u64 ALLOC_DECAY_SHIFT = 0;
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

#endif // CUSTOM_COMMON_H
