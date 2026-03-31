#pragma once
// SCOverrider — Statistical Corrector overrider for TAGE.
//
// PC-only bias tables + GEHL tables indexed by PC XOR folded global history.
// Sum precomputed in predict1 (no tage_pred dependency) → only final mux
// in predict2 depends on TAGE, keeping SC off the P2 critical path.
//
// Two-phase training: train_compute() before need_extra_cycle (no RAM writes),
// train_write() after (RAM writes on the extra cycle).

#include "../../harcom.hpp"
#include "../common.hpp"

using namespace hcm;

// Constexpr geometric history length series
template <u64 N, u64 MINH, u64 MAXH>
constexpr std::array<u64, N> sc_geometric_hists() {
  std::array<u64, N> h{};
  if constexpr (N == 1) {
    h[0] = MAXH;
  } else {
    for (u64 i = 0; i < N; i++) {
      double ratio = static_cast<double>(i) / (N - 1);
      double val = MINH * std::pow(static_cast<double>(MAXH) / MINH, ratio);
      h[i] = std::max(static_cast<u64>(val + 0.5), (i > 0) ? h[i-1] + 1 : MINH);
    }
    h[N-1] = MAXH;
  }
  return h;
}

constexpr u64 sc_clog2(u64 x) {
  if (x <= 1) return 1;
  u64 r = 0, v = x - 1;
  while (v > 0) { v >>= 1; r++; }
  return r;
}

template <
    u64 NUM_BIAS    = 3,
    u64 BIAS_LOG    = 8,       // 256 entries per bias table
    u64 SC_CTR_BITS = 6,       // signed counter width
    u64 FETCH_WIDTH = 16,
    u64 NUM_GEHL    = 4,
    u64 GEHL_LOG    = 9,       // 512 entries per GEHL table
    u64 GEHL_MINH   = 4,
    u64 GEHL_MAXH   = 64,
    i64 THRESH_INIT = 200,
    u64 THRESH_MARGIN = 4>
struct SCOverrider {

  static constexpr bool ENABLED = true;
  static constexpr u64 EXTRA_CYCLE_BUDGET = 0;

  // Derived constants
  static constexpr u64 BIAS_SIZE = 1u << BIAS_LOG;
  static constexpr u64 GEHL_SIZE = (NUM_GEHL > 0) ? (1u << GEHL_LOG) : 2;
  static constexpr u64 NUM_COMPONENTS = NUM_BIAS + NUM_GEHL;
  static constexpr u64 GEHL_N = (NUM_GEHL > 0) ? NUM_GEHL : 1;
  static constexpr auto GEHL_HISTS = sc_geometric_hists<GEHL_N, GEHL_MINH, GEHL_MAXH>();
  static constexpr u64 SC_GHIST_LEN = (NUM_GEHL > 0) ? GEHL_MAXH : 1;
  static constexpr u64 SUM_BITS = SC_CTR_BITS + sc_clog2(NUM_COMPONENTS > 0 ? NUM_COMPONENTS : 1);

  // ======== Storage ========

  // Bias: PC-only indexed. Sum precomputed in predict1.
  hcm::ram<val<SC_CTR_BITS, i64>, BIAS_SIZE> bias[NUM_BIAS]{"sc_bias"};
  arr<reg<SC_CTR_BITS, i64>, NUM_BIAS> bias_w;

  // GEHL
  hcm::ram<val<SC_CTR_BITS, i64>, GEHL_SIZE> gehl[GEHL_N]{"sc_gehl"};
  arr<reg<SC_CTR_BITS, i64>, GEHL_N> gehl_w;

  // SC global history
  reg<SC_GHIST_LEN> sc_ghist;

  // Precomputed in prefetch — available in predict2 with zero tage_pred dependency
  reg<1> pre_tage_likely_wrong;  // sum < 0
  reg<1> pre_sc_confident;       // |sum| > threshold

  // Saved state for training
  reg<BIAS_LOG> saved_bias_idx;
  arr<reg<GEHL_LOG>, GEHL_N> saved_gehl_idx;
  reg<SUM_BITS, i64> saved_sum;
  reg<1> saved_tage_pred;

  // Staged write data
  arr<reg<SC_CTR_BITS>, NUM_BIAS> staged_bias;
  arr<reg<SC_CTR_BITS>, GEHL_N> staged_gehl;

  // Threshold
  reg<SUM_BITS, i64> sc_threshold;
  bool thresh_initialized = false;

  // Control
  bool did_prefetch = false;
  bool first_lookup = true;
  bool should_write = false;

  // Stats
#ifdef TAGE_MONITOR
  struct LoopStats {
    u64 lookups = 0;
    u64 hits = 0;
    u64 overrides = 0;
    u64 alloc_writes = 0;
    u64 update_writes = 0;
    u64 sc_correct = 0;
    u64 sc_wrong = 0;
    u64 prefetch_conf_nonzero = 0;
    u64 prefetch_num_nonzero = 0;
  } stats;
#endif

  struct LookupResult { val<1> candidate; val<1> pred; };

  // History fold helper
  template <u64 H, u64 F>
  static val<F> fold_ghist(reg<SC_GHIST_LEN> &ghist) {
    static_assert(H <= SC_GHIST_LEN);
    static_assert(F > 0);
    val<H> hist_bits = val<H>{ghist};
    auto chunks = hist_bits.make_array(val<F>{});
    return chunks.fo1().fold_xor();
  }

  // ======== prefetch — predict1 ========
  // Reads all RAMs AND precomputes sum/confidence into regs.
  // No tage_pred dependency → everything available before predict2.
  void prefetch(val<64> &block_pc, [[maybe_unused]] u64 raw_pc) {
    if (!thresh_initialized) {
      sc_threshold = val<SUM_BITS, i64>{THRESH_INIT};
      thresh_initialized = true;
    }

    block_pc.fanout(hard<NUM_BIAS + NUM_GEHL + 2>{});

    // Bias reads (PC-only)
    val<BIAS_LOG> bidx = val<BIAS_LOG>{block_pc >> 2};
    bidx.fanout(hard<NUM_BIAS + 1>{});
    for (u64 i = 0; i < NUM_BIAS; i++)
      bias_w[i] = bias[i].read(bidx);
    saved_bias_idx = bidx;

    // GEHL reads
    if constexpr (NUM_GEHL > 0) {
      sc_ghist.fanout(hard<NUM_GEHL + 2>{});
      static_loop<NUM_GEHL>([&]<u64 I>() {
        constexpr u64 H = GEHL_HISTS[I];
        val<GEHL_LOG> pc_bits = val<GEHL_LOG>{block_pc >> 2};
        val<GEHL_LOG> fold = fold_ghist<H, GEHL_LOG>(sc_ghist);
        val<GEHL_LOG> gidx = pc_bits ^ fold;
        gidx.fanout(hard<2>{});
        gehl_w[I] = gehl[I].read(gidx);
        saved_gehl_idx[I] = gidx;
      });
    }

    // Precompute sum from all components (no tage_pred dependency)
    arr<val<SC_CTR_BITS, i64>, NUM_COMPONENTS> all_w = [&](u64 i) -> val<SC_CTR_BITS, i64> {
      if constexpr (NUM_GEHL > 0) {
        if (i < NUM_BIAS) return bias_w[i];
        return gehl_w[i - NUM_BIAS];
      } else {
        return bias_w[i];
      }
    };
    auto total = all_w.fo1().fold_add();
    total.fanout(hard<4>{});

    // Precompute confidence and direction into regs
    pre_tage_likely_wrong = (total < hard<0>{});
    sc_threshold.fanout(hard<4>{});
    pre_sc_confident = (total > sc_threshold) | ((total + sc_threshold) < hard<0>{});
    saved_sum = total;

    did_prefetch = true;
    first_lookup = true;
  }

  // ======== lookup — predict2 / reuse_predict2 ========
  // Pure combinational. Only the final select depends on tage_pred.
  LookupResult lookup([[maybe_unused]] val<64> &inst_pc,
                      val<1> &tage_pred,
                      [[maybe_unused]] val<1> &tage_low_conf) {
    // Pre-extracted from predict1 — zero additional latency
    pre_tage_likely_wrong.fanout(hard<2>{});
    pre_sc_confident.fanout(hard<2>{});
    val<1> candidate = val<1>{pre_sc_confident} & val<1>{pre_tage_likely_wrong};
    val<1> sc_pred = ~tage_pred;  // override = flip TAGE

    if (first_lookup) {
      saved_tage_pred = tage_pred;
      first_lookup = false;
    }

#ifdef TAGE_MONITOR
    stats.lookups++;
#ifdef CHEATING_MODE
    if (static_cast<u64>(candidate)) stats.overrides++;
#endif
#endif

    return {candidate, sc_pred};
  }

  // ======== save_branch — update_condbr ========
  void save_branch([[maybe_unused]] val<64> &branch_pc,
                   [[maybe_unused]] u64 idx) {}

  // ======== train_compute — BEFORE need_extra_cycle ========
  void train_compute(arr<val<1>, FETCH_WIDTH> &branch_taken,
                     [[maybe_unused]] arr<val<1>, FETCH_WIDTH> &is_branch,
                     val<1> &mispredict,
                     [[maybe_unused]] val<1> &correct_pred,
                     u64 num_branch) {
    if (num_branch == 0 || !did_prefetch) {
      did_prefetch = false;
      should_write = false;
      return;
    }

    val<1> actual = branch_taken[0];
    val<1> tage_correct = ~(actual ^ val<1>{saved_tage_pred});
    tage_correct.fanout(hard<NUM_COMPONENTS + 4>{});

    // Selective training
    saved_sum.fanout(hard<6>{});
    val<1> sc_said_wrong = (saved_sum < hard<0>{});
    val<1> sc_correct = ~(sc_said_wrong ^ (actual ^ val<1>{saved_tage_pred}));
    val<1> sc_wrong = ~sc_correct;
    sc_threshold.fanout(hard<5>{});
    val<1> below_thresh = (saved_sum <= sc_threshold) & ((saved_sum + sc_threshold) >= hard<0>{});
    mispredict.fanout(hard<3>{});
    val<1> should_train = sc_wrong | below_thresh;
    should_train.fanout(hard<NUM_BIAS + NUM_GEHL + 1>{});

    // Stage bias writes (mux old/new)
    for (u64 i = 0; i < NUM_BIAS; i++) {
      val<SC_CTR_BITS, i64> old_w = val<SC_CTR_BITS, i64>{bias_w[i]};
      val<SC_CTR_BITS, i64> new_w = update_ctr(old_w, tage_correct);
      staged_bias[i] = val<SC_CTR_BITS>{select(should_train, new_w, old_w)};
    }

    // Stage GEHL writes
    if constexpr (NUM_GEHL > 0) {
      for (u64 i = 0; i < NUM_GEHL; i++) {
        val<SC_CTR_BITS, i64> old_w = val<SC_CTR_BITS, i64>{gehl_w[i]};
        val<SC_CTR_BITS, i64> new_w = update_ctr(old_w, tage_correct);
        staged_gehl[i] = val<SC_CTR_BITS>{select(should_train, new_w, old_w)};
      }
      sc_ghist = (sc_ghist << 1) | val<SC_GHIST_LEN>{actual};
    }

    should_write = true;

#ifdef TAGE_MONITOR
#ifdef CHEATING_MODE
    if (should_write) {
      stats.update_writes++;
      bool sc_pred_v = (static_cast<i64>(saved_sum) < 0);
      bool tage_was_correct = !static_cast<u64>(mispredict);
      if (!tage_was_correct && sc_pred_v) stats.sc_correct++;
      if (tage_was_correct && sc_pred_v) stats.sc_wrong++;
    }
#endif
#endif

    did_prefetch = false;
  }

  // ======== train_write — AFTER need_extra_cycle ========
  // Called from C++ if (should_write) in Tage.hpp. Plain ram.write(), no execute_if.
  void train_write() {
    saved_bias_idx.fanout(hard<NUM_BIAS + 1>{});
    for (u64 i = 0; i < NUM_BIAS; i++)
      bias[i].write(val<BIAS_LOG>{saved_bias_idx}, val<SC_CTR_BITS, i64>{staged_bias[i]});
    if constexpr (NUM_GEHL > 0) {
      for (u64 i = 0; i < NUM_GEHL; i++)
        gehl[i].write(val<GEHL_LOG>{saved_gehl_idx[i]}, val<SC_CTR_BITS, i64>{staged_gehl[i]});
    }
    should_write = false;
  }

  void print_stats() const {
#ifdef TAGE_MONITOR
    std::cerr << "\n=== SC Overrider Stats ===\n"
              << "  Lookups: " << stats.lookups
              << "  Overrides: " << stats.overrides
              << "  Updates: " << stats.update_writes
              << "\n  SC correct (TAGE wrong): " << stats.sc_correct
              << "  SC wrong (TAGE right): " << stats.sc_wrong << "\n";
#endif
  }
};
