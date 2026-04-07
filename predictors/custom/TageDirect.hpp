#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"
#ifdef TAGE_MONITOR
#include "TDMonitor.hpp"
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
// Banked RAM wrapper — rwram with configurable bank index bits
// ============================================================================
//
// Like common.hpp's rwram but lets you choose which address bits select the
// bank. BANK_SHIFT=0 uses the low bits (same as rwram). BANK_SHIFT=K uses
// bits [K, K+log2(B)) as the bank ID.
//
// Template params:
//   N = data width in bits
//   M = total entries (power of 2)
//   B = number of banks (power of 2, >= 2)
//   BANK_SHIFT = which bit position starts the bank selection (default 0)

template <u64 N, u64 M, u64 B, u64 BANK_SHIFT = 0>
struct td_rwram {
  static_assert(std::has_single_bit(M));
  static_assert(B >= 2 && B <= 64);
  static_assert(std::has_single_bit(B));
  static constexpr u64 A = std::bit_width(M - 1);     // address bits
  static constexpr u64 E = M / B;                      // entries per bank
  static_assert(E > 1);
  static constexpr u64 L = std::bit_width(E - 1);      // local address bits
  static constexpr u64 BANK_BITS = std::bit_width(B - 1); // bank ID bits
  static_assert(A == L + BANK_BITS);
  static_assert(BANK_SHIFT + BANK_BITS <= A);

  hcm::ram<val<N>, E> bank[B];
  reg<B> read_bank;

  // buffered write
  reg<B> write_bank;
  reg<L> write_localaddr;
  reg<N> write_data;

  td_rwram(const char *label = "") : bank{label} {}

  // Extract bank ID and local address from a full address.
  // With BANK_SHIFT=0: addr = [local_hi | bankid | local_lo] where local_lo
  // is empty. With BANK_SHIFT=K: bankid = addr[K +: BANK_BITS], local =
  // remaining bits concatenated.
  auto split_addr(val<A> addr) {
    if constexpr (BANK_SHIFT == 0) {
      // Low bits = bank ID (same as reference rwram)
      return hcm::split<L, BANK_BITS>(addr);
    } else if constexpr (BANK_SHIFT + BANK_BITS == A) {
      // High bits = bank ID
      return std::pair{val<L>{addr}, val<BANK_BITS>{addr >> BANK_SHIFT}};
    } else {
      // Middle bits = bank ID: extract and reconstruct local addr
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

#ifdef TD_VERBOSE
  void debug_read_info(const char *name, u64 table_idx, val<A> addr) {
    u64 a = static_cast<u64>(addr);
    u64 bid = (a >> BANK_SHIFT) & (B - 1);
    u64 wb = static_cast<u64>(write_bank);
    std::cerr << "  " << name << "[" << table_idx << "] addr=" << a
              << " bank=" << bid;
    if ((wb >> bid) & 1)
      std::cerr << " (pending-write in bank)";
    std::cerr << "\n";
  }

  // Debug: extract bank ID from address and report conflict with read_bank.
  void debug_write_info(const char *name, u64 table_idx, val<A> addr,
                        val<1> noconflict) {
    u64 a = static_cast<u64>(addr);
    u64 bid = (a >> BANK_SHIFT) & (B - 1);
    u64 rb = static_cast<u64>(read_bank);
    u64 wb = static_cast<u64>(write_bank);
    bool read_conflict = (rb >> bid) & 1;
    bool write_buffered = !static_cast<u64>(noconflict);
    std::cerr << "  " << name << "[" << table_idx << "] addr=" << a
              << " bank=" << bid;
    if (write_buffered)
      std::cerr << " BUFFERED";
    else
      std::cerr << " DIRECT";
    if (read_conflict && write_buffered)
      std::cerr << " (read-conflict)";
    if (wb)
      std::cerr << " (pending_wb=0x" << std::hex << wb << std::dec << ")";
    std::cerr << "\n";
  }
#endif
};

// ============================================================================
// TAGE Table Storage — tag/pred use plain ram<>, hyst/u use td_rwram
// ============================================================================
//
// tag_ram, pred_ram: plain ram<>. Writes are gated by conditions that imply
// extra_cycle (allocate ⊆ mispredict, g_weak ⊆ some_badpred1), so they
// never conflict with P2 reads.
//
// hyst_ram, u_ram: td_rwram with B banks. Can read+write same cycle via
// banking. Write takes a noconflict bit (= extra_cycle) following the
// reference tage.hpp pattern.

template <u64 TABLE_SIZE, u64 TAG_WIDTH, u64 CTR_WIDTH, u64 HYST_WIDTH,
          u64 U_WIDTH, u64 RWRAM_BANKS = 4, u64 RWRAM_BANK_SHIFT = 0>
struct TDTable {
  static constexpr u64 IDX_BITS = td::clog2(TABLE_SIZE);
  static constexpr u64 table_size = TABLE_SIZE;
  static constexpr u64 tag_width = TAG_WIDTH;
  static constexpr u64 ctr_width = CTR_WIDTH;
  static constexpr u64 hyst_width = HYST_WIDTH;
  static constexpr u64 u_width = U_WIDTH;

  static_assert(TABLE_SIZE >= 2 * RWRAM_BANKS,
                "TABLE_SIZE must be >= 2*RWRAM_BANKS for td_rwram");

  hcm::ram<val<TAG_WIDTH>, TABLE_SIZE> tag_ram{"td_tag"};
  hcm::ram<val<CTR_WIDTH>, TABLE_SIZE> pred_ram{"td_pred"};
  td_rwram<std::max(u64(1), HYST_WIDTH), TABLE_SIZE, RWRAM_BANKS,
           RWRAM_BANK_SHIFT>
      hyst_ram{"td_hyst"};
  td_rwram<U_WIDTH, TABLE_SIZE, RWRAM_BANKS, RWRAM_BANK_SHIFT>
      u_ram{"td_u"};
};

// Generate a TDTable type from config arrays at index I
template <typename Cfg, u64 I, u64 BANKS, u64 BSHIFT>
using TDTableAt =
    TDTable<Cfg::TABLE_SIZE[I], Cfg::TAG_WIDTH[I], Cfg::CTR_WIDTH[I],
            Cfg::HYST_WIDTH[I], Cfg::U_WIDTH[I], BANKS, BSHIFT>;

// Build a tuple of TDTable types
template <typename Cfg, u64 BANKS, u64 BSHIFT, typename Seq>
struct TDMakeTableTuple;

template <typename Cfg, u64 BANKS, u64 BSHIFT, u64... Is>
struct TDMakeTableTuple<Cfg, BANKS, BSHIFT, std::index_sequence<Is...>> {
  using type = std::tuple<TDTableAt<Cfg, Is, BANKS, BSHIFT>...>;
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
          template <u64> class FoldFn_V,
          u64 RWRAM_BANKS_V = 4, u64 RWRAM_BANK_SHIFT_V = 0>
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
      typename td::TDMakeTableTuple<TableCfg, RWRAM_BANKS_V, RWRAM_BANK_SHIFT_V,
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
  reg<N + 1>
      rank; // one-hot: rank of current branch in block (gshareN_ahead pattern)

  // P1 storage
  hcm::ram<val<LANES>, P1_ENTRIES> p1_pred{"P1 pred"};
  zone UPDATE_ONLY;
  hcm::ram<val<1>, P1_ENTRIES> p1_hyst[N]{"P1 hyst"};

  // TAGE tables
  Tables tables;

#ifdef TAGE_MONITOR
  TDMonitor<NUM_TABLES, LANES, P1_ENTRIES, RWRAM_BANKS_V> mon;
  // Shadow state for meta tracking (set in predict2, read in update_cycle)
  std::array<bool, LANES> mon_altsel{};     // meta chose alt for this rank
  std::array<bool, LANES> mon_meta_active{}; // meta override was active
  ~TageDirectImpl() { mon.print_summary(); }
#endif

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
#ifdef TD_VERBOSE
      std::cerr << "P1: gshare p1_pred READ\n";
#endif
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

  val<1> predict2(val<64> inst_pc) override {
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
#ifdef TD_VERBOSE
      std::cerr << "P2: tage tag_ram[" << I << "] READ\n";
#endif
      readt[I] = t.tag_ram.read(tidx<I>(gindex[I]));
#ifdef TD_VERBOSE
      std::cerr << "P2: tage pred_ram[" << I << "] READ\n";
#endif
      readc[I] = t.pred_ram.read(tidx<I>(gindex[I]));
      if constexpr (MAX_HYST_WIDTH > 0) {
#ifdef TD_VERBOSE
        t.hyst_ram.debug_read_info("hyst_ram", I, tidx<I>(gindex[I]));
#endif
        readh[I] = t.hyst_ram.read(tidx<I>(gindex[I]));
#ifdef TAGE_MONITOR
        mon.record_rwram_read(I);
#endif
      }
#ifdef TD_VERBOSE
      t.u_ram.debug_read_info("u_ram", I, tidx<I>(gindex[I]));
#endif
      readu[I] = t.u_ram.read(tidx<I>(gindex[I]));
#ifdef TAGE_MONITOR
      mon.record_rwram_read(I);
#endif
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
#ifdef TAGE_MONITOR
      for (u64 r = 0; r < LANES; r++) {
        bool has_alt = static_cast<u64>(match2[r] != hard<0>{});
        bool pri_ne_alt = static_cast<u64>(pred1_tage[r]) != static_cast<u64>(pred2_tage[r]);
        mon_altsel[r] = static_cast<u64>(altsel[r]);
        mon_meta_active[r] = has_alt && pri_ne_alt && static_cast<u64>(altsel[r]);
      }
#endif
    } else {
      p2 = pred1_tage.concat();
#ifdef TAGE_MONITOR
      for (u64 r = 0; r < LANES; r++) { mon_altsel[r] = false; mon_meta_active[r] = false; }
#endif
    }

    p2.fanout(hard<LANES>{});
    val<1> taken = p2 >> num_branch;
    taken.fanout(hard<2>{});
    reuse_prediction(~val<1>{block_entry >> (LINEINST - 1)});
    return taken;
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) override {
    val<1> taken = p2 >> num_branch;
    taken.fanout(hard<2>{});
    block_size++;
    reuse_prediction(~line_end());
    return taken;
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

#ifdef TD_VERBOSE
    std::cerr << "UC: ENTER (num_branch=" << num_branch
              << " misp=" << static_cast<u64>(block_end_info.is_mispredict)
              << ")\n";
#endif
    if (num_branch == 0) {
      // No conditional branches — just update history
#ifdef TD_VERBOSE
      std::cerr << "UC: EXIT (no branches)\n";
#endif
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

    // ================================================================
    // UPDATE CYCLE — Cycle 1: combinational logic + RAM reads
    // All code above need_extra_cycle() executes in the same cycle as
    // predict1/predict2. No RAM writes allowed here — only reads and
    // pure combinational logic on values already in regs.
    // ================================================================

    // ---- Fanouts for values read in predict2 (stored in regs) ----
    mispredict.fanout(hard<NUM_TABLES + 2>{});
    val<1> correct_pred = ~mispredict;
    correct_pred.fanout(hard<NUM_TABLES + 2>{});
    index1.fanout(hard<LANES * 3>{});
    gindex.fanout(hard<4>{});
    htag.fanout(hard<3>{});
    readt.fanout(hard<4>{});   // tag values read from TAGE tables in P2
    readc.fanout(hard<2>{});   // counter values read in P2
    if constexpr (MAX_HYST_WIDTH > 0)
      readh.fanout(hard<3>{}); // hysteresis values read in P2
    match1.fanout(hard<3>{});  // provider match vectors computed in P2
    match2.fanout(hard<2>{});  // alt match vectors computed in P2
    pred1_tage.fanout(hard<2>{});
    pred2_tage.fanout(hard<2 + NUM_TABLES>{});
    branch_dir.fanout(hard<2>{});
    gfolds.fanout(hard<2>{});
    readu.fanout(hard<2>{});   // u-bit values read in P2
    X.fanout(hard<LANES + 1>{});
    if constexpr (USE_META_V)
      meta.fanout(hard<2>{});

    // ---- TAGE combinational logic ----
    // All pure logic — computes allocation decisions, counter update
    // conditions, etc. from the reg values read in predict2.
    // No RAM access here.

    // last_rank: which lane held the last conditional branch in this block
    val<LOG_LANES> last_rank = val<LOG_LANES>{num_branch - 1};
    last_rank.fanout(hard<4 * NUM_TABLES + 2>{});

    // Per-lane: was this lane used by a branch in this block?
    arr<val<1>, LANES> is_branch = [&](u64 r) -> val<1> {
      return val<1>{r < num_branch ? 1u : 0u};
    };
    is_branch.fanout(hard<4>{});

    // Per-lane: actual branch direction (0 for unused lanes)
    arr<val<1>, LANES> branch_taken = [&](u64 r) -> val<1> {
      if (r < num_branch) return branch_dir[r];
      return val<1>{0};
    };
    branch_taken.fanout(hard<3>{});

    // Restrict match vectors to lanes that actually had branches
    arr<val<NUM_TABLES + 1>, LANES> actual_match1 = [&](u64 r) {
      return select(is_branch[r], match1[r], val<NUM_TABLES + 1>{0});
    };
    actual_match1.fanout(hard<2>{});

    // primary_mask: OR of all per-lane provider matches — tells which
    // TAGE tables are the provider for at least one branch
    val<NUM_TABLES> primary_mask = actual_match1.fold_or();
    primary_mask.fanout(hard<2>{});
    arr<val<1>, NUM_TABLES> primary = primary_mask.make_array(val<1>{});
    primary.fanout(hard<3>{});

    // ---- Allocation decision (combinational) ----
    // On misprediction, try to allocate entries in tables above the provider.
    // mispmask gates allocation to mispredictions only.
    val<NUM_TABLES> mispmask = mispredict.replicate(hard<NUM_TABLES>{}).concat();

    // Tag comparison for the last branch's rank specifically
    static_loop<NUM_TABLES>([&]<u64 I>() {
      using Table = std::tuple_element_t<I, Tables>;
      static constexpr u64 PER_HTAG = Table::tag_width - LOG_LANES;
      last_tagcmp_reg[I] = (val<LOG_LANES>{readt[I] >> PER_HTAG} == last_rank)
                         & (val<PER_HTAG>{readt[I]} == val<PER_HTAG>{htag[I]});
    });
    val<NUM_TABLES + 1> last_match1 =
        last_tagcmp_reg.fo1().append(1).concat().one_hot();
    last_match1.fanout(hard<2>{});

    // postmask: tables above the provider (candidates for allocation)
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

    // candallocmask: tables that are both above provider AND have u=0
    val<NUM_TABLES> candallocmask = [&]() -> val<NUM_TABLES> {
      if constexpr (AllocCfg::CONF_GATE) {
        arr<val<1>, NUM_TABLES> weak_entry = [&](u64 i) -> val<1> {
          if constexpr (MAX_CTR_WIDTH == 1) return val<1>{1};
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

    // allocate[i]: final per-table allocation decision (one-hot or two-hot)
    arr<val<1>, NUM_TABLES> allocate = [&]() -> arr<val<1>, NUM_TABLES> {
      if constexpr (AllocCfg::MAX_ALLOC >= 2) {
        val<NUM_TABLES> pick2 = [&]() -> val<NUM_TABLES> {
          val<NUM_TABLES> basic2 = (collamask ^ collamask1).one_hot();
          if constexpr (AllocCfg::NON_CONSECUTIVE) {
            val<NUM_TABLES> neighbors = (collamask1 << 1) | (collamask1 >> 1);
            val<NUM_TABLES> nc_mask = (collamask ^ collamask1) & ~neighbors;
            val<NUM_TABLES> nc_pick = nc_mask.reverse().one_hot();
            return select(nc_mask != hard<0>{}, nc_pick, basic2);
          } else { return basic2; }
        }();
        return (collamask1 | pick2).reverse().make_array(val<1>{});
      } else {
        val<NUM_TABLES> collamask2 = (collamask ^ collamask1).one_hot();
        val<NUM_TABLES> collamask12 =
            select(val<2>{std::rand()} == hard<0>{}, collamask2.fo1(), collamask1);
        return collamask12.fo1().reverse().make_array(val<1>{});
      }
    }();
    allocate.fanout(hard<7>{});

    // ---- Per-table update conditions (combinational) ----

    // bdir[i]: branch direction for the rank stored in table i's tag.
    // On allocation, uses last_rank instead.
    arr<val<1>, NUM_TABLES> bdir = [&](u64 i) {
      val<LOG_LANES> stored_rank = readt[i] >> MAX_HTAGBITS;
      val<LOG_LANES> use_rank = select(allocate[i], last_rank, stored_rank.fo1());
      return branch_dir.select(use_rank);
    };
    bdir.fanout(hard<2>{});

    // badpred1[i]: did table i's counter predict wrong?
    arr<val<1>, NUM_TABLES> badpred1 = [&](u64 i) -> val<1> {
      if constexpr (MAX_CTR_WIDTH == 1) return readc[i] != bdir[i];
      else return val<1>{readc[i] >> hard<MAX_CTR_WIDTH - 1>{}} != bdir[i];
    };
    badpred1.fanout(hard<3>{});

    // altdiffer[i]: does provider disagree with alt prediction?
    arr<val<1>, NUM_TABLES> altdiffer = [&](u64 i) -> val<1> {
      auto pred_dir = [&]() -> val<1> {
        if constexpr (MAX_CTR_WIDTH == 1) return readc[i];
        else return readc[i] >> hard<MAX_CTR_WIDTH - 1>{};
      }();
      val<LOG_LANES> stored_rank = readt[i] >> MAX_HTAGBITS;
      return pred_dir != pred2_tage.select(stored_rank.fo1());
    };

    // goodpred[i]: was this table's prediction correct for its stored rank?
    arr<val<1>, NUM_TABLES> goodpred = [&](u64 i) {
      val<LOG_LANES> stored_rank = readt[i] >> MAX_HTAGBITS;
      return (stored_rank.fo1() != last_rank) | correct_pred;
    };
    goodpred.fanout(hard<2>{});

    // g_weak[i]: provider with wrong prediction AND weak hysteresis — flip counter
    arr<val<1>, NUM_TABLES> g_weak = [&](u64 i) -> val<1> {
      if constexpr (MAX_HYST_WIDTH > 0)
        return primary[i] & badpred1[i] & (readh[i] == hard<0>{});
      else
        return primary[i] & badpred1[i];
    };
    g_weak.fanout(hard<2>{});

    // P1 vs P2 disagreement — need extra cycle to update P1
    val<LANES> p1_concat = pred.concat();
    val<LANES> disagree_mask = (p1_concat ^ p2) & is_branch.concat();
    disagree_mask.fanout(hard<2>{});

    // U-bit clear helpers (combinational)
    arr<val<1>, NUM_TABLES> update_u = [&](u64 i) {
      return primary[i] & altdiffer[i].fo1();
    };
    val<1> noalloc = (candallocmask == hard<0>{});
    val<NUM_TABLES> uclearmask =
        postmask & noalloc.fo1().replicate(hard<NUM_TABLES>{}).concat();
    arr<val<1>, NUM_TABLES> uclear = uclearmask.fo1().make_array(val<1>{});
    uclear.fanout(hard<2>{});

    // ---- P1 gshare + TAGE reads / need_extra / writes ----
    // Same structure as gshareN reference: reads in cycle 1, writes in cycle 2.

    // Combinational: compute which lanes were accessed by branches in this
    // block. access[i] = OR of all X.rotate_left(branch_rank) masks — 1 if lane
    // i was used.
    arr<val<1>, LANES> access =
        arr<val<LANES>, LANES>{[&](u64 i) -> val<LANES> {
          return X.rotate_left(i) & val<LANES>{-(i < num_branch)};
        }}.fold_or()
            .make_array(val<1>{});

    // Combinational: identify the lane that holds the mispredicted branch.
    // misp_bank is a one-hot mask ANDed with mispredict — all zero on correct
    // prediction.
    val<LANES> misp_bank = X.rotate_left(num_branch - 1) &
                           mispredict.replicate(hard<LANES>{}).concat();
    arr<val<1>, LANES> mispredicted = misp_bank.fo1().make_array(val<1>{});
    mispredicted.fanout(hard<2>{});

    // Cycle 1 RAM read: read hysteresis for the mispredicted lane.
    // Gated by mispredicted[i] so only one lane's RAM is actually read.
    // On correct prediction, mispredicted is all-zero so no RAM is accessed.
    arr<val<1>, LANES> weak = [&](u64 i) -> val<1> {
      if (i >= N) return val<1>{0};
      return execute_if(mispredicted[i], [&]() {
#ifdef TD_VERBOSE
        if (static_cast<u64>(mispredict) && static_cast<u64>(mispredicted[i]))
          std::cerr << "UC: gshare p1_hyst[" << i << "] READ\n";
#endif
        return p1_hyst[i].read(index1);
      });
    };

    // ================================================================
    // Cycle boundary: grant an extra cycle for RAM writes.
    // tag_ram/pred_ram use plain ram<> — need extra_cycle for writes.
    // hyst_ram/u_ram use td_rwram — pass extra_cycle as noconflict,
    // so they can buffer writes when extra_cycle=0.
    // Everything above = cycle 1. Everything below = cycle 2.
    // ================================================================
    val<1> extra_cycle = [&]() -> val<1> {
      if constexpr (AllocCfg::MISPREDICT_ONLY_WRITE) {
        return mispredict;
      } else {
        val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
        return some_badpred1.fo1() | mispredict;
      }
    }();
#ifdef TD_VERBOSE
    std::cerr << "UC: extra_cycle=" << static_cast<u64>(extra_cycle) << "\n";
#endif
    extra_cycle.fanout(hard<NUM_TABLES * 2 + 1>{});
    need_extra_cycle(extra_cycle);

#ifdef TAGE_MONITOR
    // Block stats
    mon.record_block(static_cast<u64>(block_entry), block_size, num_branch,
                     static_cast<u64>(extra_cycle));
    // Per-branch outcome + prediction tracking
    for (u64 r = 0; r < num_branch; r++) {
      bool actual = static_cast<u64>(branch_dir[r]);
      bool p1_pr = static_cast<u64>(pred[r]);
      bool p2_pr = static_cast<u64>((p2 >> r) & hard<1>{});
      bool misp = (r == num_branch - 1) ? static_cast<u64>(mispredict) : false;
      mon.record_prediction(r, static_cast<u64>(match1[r]),
                            static_cast<u64>(match2[r]),
                            mon_meta_active[r], mon_altsel[r],
                            p1_pr, p2_pr);
      mon.record_outcome(r, actual, misp);
    }
    // Tag match tracking
    for (u64 i = 0; i < NUM_TABLES; i++) {
      bool matched = static_cast<u64>(htagcmp_reg[i]);
      mon.record_tag_lookup(i, matched);
    }
    // Allocation
    if (static_cast<u64>(mispredict)) {
      u64 amask = 0;
      for (u64 i = 0; i < NUM_TABLES; i++)
        if (static_cast<u64>(allocate[i])) amask |= (u64(1) << i);
      mon.record_allocation(amask != 0, amask);
    }
#endif

    // ================================================================
    // Cycle 2: RAM writes — P1 gshare + TAGE tables
    // All writes are gated by execute_if so they only fire when needed.
    // ================================================================

    // ---- P1 gshare writes ----
    // Flip prediction bit if hysteresis was weak on misprediction.
    // p1_pred was read in predict1 (cycle 1), safe to write in cycle 2.
    execute_if(mispredict, [&]() {
      arr<val<1>, LANES> stored = unordered_pred1.make_array(val<1>{});
      arr<val<1>, LANES> bundle = [&](u64 i) {
        return select(weak[i].fo1(), branch_dir[num_branch - 1],
                      stored[i].fo1());
      };
#ifdef TD_VERBOSE
      if (static_cast<u64>(mispredict))
        std::cerr << "UC: gshare p1_pred WRITE\n";
#endif
      p1_pred.write(index1, bundle.fo1().concat());
    });

    // Update P1 hysteresis for all accessed lanes.
    // p1_hyst[i] was read in cycle 1, safe to write in cycle 2.
    for (u64 i = 0; i < N; i++) {
      execute_if(access[i].fo1(), [&]() {
#ifdef TD_VERBOSE
        if (static_cast<u64>(mispredict) && static_cast<u64>(access[i]))
          std::cerr << "UC: gshare p1_hyst[" << i << "] WRITE\n";
#endif
        p1_hyst[i].write(index1, mispredicted[i].fo1());
#ifdef TAGE_MONITOR
        mon.record_p1_write(i, static_cast<u64>(index1));
#endif
      });
    }

    // ---- TAGE tag write (allocation) ----
    // Write new tag = concat(last_rank, htag) into allocated table entries.
    // tag_ram was read in predict2 (cycle 1), safe to write in cycle 2.
    // Gate: allocate[I] ⊆ mispredict ⊆ extra_cycle, so no third arg needed.
    static_loop<NUM_TABLES>([&]<u64 I>() {
      execute_if(allocate[I], [&]() {
#ifdef TD_VERBOSE
        if (static_cast<u64>(allocate[I]))
          std::cerr << "UC: tage tag_ram[" << I << "] WRITE (alloc)\n";
#endif
        std::get<I>(tables).tag_ram.write(tidx<I>(gindex[I]),
                                          concat(last_rank, htag[I]));
      });
    });

    // ---- TAGE counter update ----
    // pred_ram was read in predict2 (cycle 1), safe to write in cycle 2.
    if constexpr (AllocCfg::MISPREDICT_ONLY_WRITE) {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        execute_if(mispredict & (g_weak[I].fo1() | allocate[I]), [&]() {
#ifdef TD_VERBOSE
          std::cerr << "UC: tage pred_ram[" << I << "] WRITE (misp_only)\n";
#endif
          if constexpr (MAX_CTR_WIDTH == 1) {
            std::get<I>(tables).pred_ram.write(tidx<I>(gindex[I]), bdir[I]);
          } else {
            auto init_ctr = select(bdir[I],
                val<MAX_CTR_WIDTH>{(1u << MAX_CTR_WIDTH) - 1},
                val<MAX_CTR_WIDTH>{0});
            std::get<I>(tables).pred_ram.write(tidx<I>(gindex[I]),
                select(allocate[I], init_ctr, val<MAX_CTR_WIDTH>{bdir[I]}));
          }
        });
      });
    } else {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> old_dir = [&]() -> val<1> {
          if constexpr (MAX_CTR_WIDTH == 1) return readc[I];
          else return readc[I] >> hard<MAX_CTR_WIDTH - 1>{};
        }();
        val<1> pred_changed = (old_dir != bdir[I]) | allocate[I];
        execute_if((g_weak[I].fo1() | allocate[I]) & pred_changed, [&]() {
#ifdef TD_VERBOSE
          std::cerr << "UC: tage pred_ram[" << I << "] WRITE (standard)\n";
#endif
          std::get<I>(tables).pred_ram.write(tidx<I>(gindex[I]), bdir[I]);
        });
      });
    }

    // ---- TAGE hysteresis update ----
    // hyst_ram was read in predict2 (cycle 1), safe to write in cycle 2.
    if constexpr (MAX_HYST_WIDTH > 0 && !AllocCfg::MISPREDICT_ONLY_WRITE) {
      static constexpr u64 HW = std::max(u64(1), MAX_HYST_WIDTH);
      static constexpr u64 HMAX = (u64(1) << HW) - 1;
      if constexpr (AllocCfg::PARTIAL_UPDATE) altdiffer.fanout(hard<2>{});
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
          auto newhyst = select(allocate[I], val<HW>{0},
                                td::update_ctr(readh[I], ~badpred1[I]));
#ifdef TD_VERBOSE
          std::get<I>(tables).hyst_ram.debug_write_info("hyst_ram", I, tidx<I>(gindex[I]), extra_cycle);
#endif
          std::get<I>(tables).hyst_ram.write(tidx<I>(gindex[I]), newhyst.fo1(), extra_cycle);
#ifdef TAGE_MONITOR
          { auto &r = std::get<I>(tables).hyst_ram;
            u64 bid = (static_cast<u64>(gindex[I]) >> RWRAM_BANK_SHIFT_V) & (RWRAM_BANKS_V - 1);
            mon.record_rwram_write(I, bid, static_cast<u64>(extra_cycle),
                                   static_cast<u64>(r.read_bank), static_cast<u64>(r.write_bank)); }
#endif
        });
      });
    } else if constexpr (MAX_HYST_WIDTH > 0 && AllocCfg::MISPREDICT_ONLY_WRITE) {
      static constexpr u64 HW = std::max(u64(1), MAX_HYST_WIDTH);
      static constexpr u64 HMAX = (u64(1) << HW) - 1;
      static_loop<NUM_TABLES>([&]<u64 I>() {
        execute_if(mispredict & allocate[I], [&]() {
#ifdef TD_VERBOSE
          std::get<I>(tables).hyst_ram.debug_write_info("hyst_ram", I, tidx<I>(gindex[I]), extra_cycle);
#endif
          std::get<I>(tables).hyst_ram.write(tidx<I>(gindex[I]), val<HW>{HMAX}, extra_cycle);
#ifdef TAGE_MONITOR
          { auto &r = std::get<I>(tables).hyst_ram;
            u64 bid = (static_cast<u64>(gindex[I]) >> RWRAM_BANK_SHIFT_V) & (RWRAM_BANKS_V - 1);
            mon.record_rwram_write(I, bid, static_cast<u64>(extra_cycle),
                                   static_cast<u64>(r.read_bank), static_cast<u64>(r.write_bank)); }
#endif
        });
      });
    }

    // ---- TAGE u-bit update ----
    // u_ram uses td_rwram: pass extra_cycle as noconflict for buffered writes.
    if constexpr (AllocCfg::MISPREDICT_ONLY_WRITE) {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        execute_if(mispredict & (allocate[I] | uclear[I]), [&]() {
          val<1> newu = select(allocate[I], val<1>{1}, val<1>{0});
#ifdef TD_VERBOSE
          std::get<I>(tables).u_ram.debug_write_info("u_ram", I, tidx<I>(gindex[I]), extra_cycle);
#endif
          std::get<I>(tables).u_ram.write(tidx<I>(gindex[I]), newu.fo1(), extra_cycle);
#ifdef TAGE_MONITOR
          mon.record_u_write(I, static_cast<u64>(newu));
          { auto &r = std::get<I>(tables).u_ram;
            u64 bid = (static_cast<u64>(gindex[I]) >> RWRAM_BANK_SHIFT_V) & (RWRAM_BANKS_V - 1);
            mon.record_rwram_write(I, bid, static_cast<u64>(extra_cycle),
                                   static_cast<u64>(r.read_bank), static_cast<u64>(r.write_bank)); }
#endif
        });
      });
    } else if constexpr (USE_PROB_DECAY) {
      val<DECAY_CTR_V> lfsr = val<DECAY_CTR_V>{static_cast<u64>(std::rand())};
      val<1> decay_fire = (lfsr > val<DECAY_CTR_V>{decay_threshold});
      decay_fire.fanout(hard<NUM_TABLES>{});
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> newu = goodpred[I].fo1() & ~allocate[I] & ~uclear[I] & ~decay_fire;
        val<1> u_changed = (val<1>{readu[I]} != newu);
        execute_if((update_u[I].fo1() | allocate[I] | uclear[I] | decay_fire) & u_changed, [&]() {
#ifdef TD_VERBOSE
          std::get<I>(tables).u_ram.debug_write_info("u_ram", I, tidx<I>(gindex[I]), extra_cycle);
#endif
          std::get<I>(tables).u_ram.write(tidx<I>(gindex[I]), newu.fo1(), extra_cycle);
#ifdef TAGE_MONITOR
          mon.record_u_write(I, static_cast<u64>(newu));
          { auto &r = std::get<I>(tables).u_ram;
            u64 bid = (static_cast<u64>(gindex[I]) >> RWRAM_BANK_SHIFT_V) & (RWRAM_BANKS_V - 1);
            mon.record_rwram_write(I, bid, static_cast<u64>(extra_cycle),
                                   static_cast<u64>(r.read_bank), static_cast<u64>(r.write_bank)); }
#endif
        });
      });
    } else {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> newu = goodpred[I].fo1() & ~allocate[I] & ~uclear[I];
        val<1> u_changed = (val<1>{readu[I]} != newu);
        execute_if((update_u[I].fo1() | allocate[I] | uclear[I]) & u_changed, [&]() {
#ifdef TD_VERBOSE
          std::get<I>(tables).u_ram.debug_write_info("u_ram", I, tidx<I>(gindex[I]), extra_cycle);
#endif
          std::get<I>(tables).u_ram.write(tidx<I>(gindex[I]), newu.fo1(), extra_cycle);
#ifdef TAGE_MONITOR
          mon.record_u_write(I, static_cast<u64>(newu));
          { auto &r = std::get<I>(tables).u_ram;
            u64 bid = (static_cast<u64>(gindex[I]) >> RWRAM_BANK_SHIFT_V) & (RWRAM_BANKS_V - 1);
            mon.record_rwram_write(I, bid, static_cast<u64>(extra_cycle),
                                   static_cast<u64>(r.read_bank), static_cast<u64>(r.write_bank)); }
#endif
        });
      });
    }

    // ---- Meta counter update (regs only, no RAM) ----
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
      for (u64 i = METAPIPE_V - 1; i != 0; i--) meta[i] = meta[i - 1];
      auto newmeta = meta[0] + meta_incr.fo1().fold_add();
      newmeta.fanout(hard<3>{});
      using meta_t = valt<decltype(meta[0])>;
      meta[0] = select(newmeta > meta_t::maxval, meta_t{meta_t::maxval},
                       select(newmeta < meta_t::minval, meta_t{meta_t::minval},
                              meta_t{newmeta}));
    }

    // ---- Epoch / decay threshold (regs only, no RAM except u_ram.reset) ----
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
        if constexpr (DECAY_GRAN_V == 0) return ~correct_pred;
        else return (uctr & hard<(u64(1) << DECAY_GRAN_V) - 1>{}) == hard<0>{};
      }();
      val<1> misp = ~correct_pred;
      decay_threshold = select(threshold_tick,
          DecayPolicy_V::template apply<DECAY_CTR_V>(
              decay_threshold, correct_pred, uctrsat, misp),
          val<DECAY_CTR_V>{decay_threshold});
    } else {
      // Periodic u-bit reset — u_ram.reset() is a bulk clear, safe in cycle 2.
      // uctrsat only fires on mispredictions, so extra_cycle is guaranteed.
      execute_if(uctrsat, [&]() {
#ifdef TAGE_MONITOR
        mon.record_epoch_reset();
#endif
        static_loop<NUM_TABLES>([&]<u64 I>() {
#ifdef TD_VERBOSE
          if (static_cast<u64>(uctrsat))
            std::cerr << "UC: tage u_ram[" << I << "] RESET\n";
#endif
          std::get<I>(tables).u_ram.reset();
        });
      });
    }

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

#ifdef TD_VERBOSE
    std::cerr << "UC: EXIT (full update)\n";
#endif
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
          template <u64> class TD_FOLD_FN = td::XORFold,
          u64 TD_RWRAM_BANKS = 4, u64 TD_RWRAM_BANK_SHIFT = 0>

using TageDirect =
    TageDirectImpl<TableCfg, AllocCfg, TD_LINEINST, TD_N, TD_SHARED_TAG,
                   TD_SHARED_U, TD_SHARED_HYS, TD_U_STOR_FF, TD_DECAY_CTR,
                   TD_DECAY_GRAN, TD_DECAY_POLICY, TD_P1_USE_GSHARE,
                   TD_P1_TABLE_SIZE, TD_P1_HIST, TD_USE_META, TD_METABITS,
                   TD_METAPIPE, TD_USE_PATH_HIST, TD_PATH_HIST_WIDTH,
                   TD_PATH_BITS, TD_FOLD_FN, TD_RWRAM_BANKS,
                   TD_RWRAM_BANK_SHIFT>;
