#pragma once
// SCOverrider — Statistical Corrector overrider (standalone).
//
// Perceptron-based corrector. Sums signed weights from bias tables
// indexed by PC. When sum magnitude exceeds threshold and disagrees
// with TAGE, SC overrides.
//
// Uses ram<> (not rwram). Reads in predict1, writes ONLY when extra_cycle fires.

#include "../../harcom.hpp"
#include "../common.hpp"

using namespace hcm;

template <
    u64 BIAS_LOG = 8,             // log2(bias table entries) = 256
    u64 SC_CTR_BITS = 6,         // 6-bit signed counters
    u64 FETCH_WIDTH = 16,
    i64 SC_CONF_THRESH = 4>      // confidence threshold: override when |sum| > this
struct SCOverrider {

  static constexpr bool ENABLED = true;
  static constexpr u64 BIAS_SIZE = 1u << BIAS_LOG;
  static constexpr u64 NUM_BIAS = 3;

  // SC bias tables
  hcm::ram<val<SC_CTR_BITS>, BIAS_SIZE> bias[NUM_BIAS]{"sc_bias"};
  arr<reg<SC_CTR_BITS>, NUM_BIAS> bias_w;
  u64 saved_bias_idx = 0;
  bool did_prefetch = false;

#ifdef TAGE_MONITOR
  // Named LoopStats for compatibility with TageMonitor destructor print
  struct LoopStats {
    u64 lookups = 0;         // total lookup calls
    u64 hits = 0;            // tag matches (N/A for SC — always "hits")
    u64 overrides = 0;       // confident overrides fired
    u64 alloc_writes = 0;    // N/A for SC
    u64 update_writes = 0;   // SC weight updates
    u64 sc_correct = 0;      // SC override was correct (TAGE was wrong)
    u64 sc_wrong = 0;        // SC override was wrong (TAGE was right)
    u64 prefetch_conf_nonzero = 0;
    u64 prefetch_num_nonzero = 0;
  } stats;
#endif

  struct LookupResult { val<1> candidate; val<1> pred; };

#ifdef TAGE_MONITOR
  ~SCOverrider() {
    std::cerr << "\n=== SC Overrider Stats ===\n"
              << "  Lookups: " << stats.lookups
              << "  Overrides: " << stats.overrides
              << "  Updates: " << stats.update_writes
              << "\n  SC correct (TAGE wrong): " << stats.sc_correct
              << "  SC wrong (TAGE right): " << stats.sc_wrong
              << "\n";
  }
#endif

  // ======== prefetch — called from predict1 ========
  void prefetch(val<64> &block_pc, u64 raw_pc) {
    val<BIAS_LOG> bidx = val<BIAS_LOG>{block_pc >> 2};
    bidx.fanout(hard<NUM_BIAS>{});
    for (u64 i = 0; i < NUM_BIAS; i++)
      bias_w[i] = bias[i].read(bidx);
    bias_w.fanout(hard<FETCH_WIDTH + 2>{});
    if (raw_pc != 0) saved_bias_idx = (raw_pc >> 2) & (BIAS_SIZE - 1);
    did_prefetch = true;
  }

  // ======== lookup — called from predict2 / reuse_predict2 ========
  LookupResult lookup([[maybe_unused]] val<64> &inst_pc) {
    // Sum bias weights
    arr<val<SC_CTR_BITS, i64>, NUM_BIAS> all_w = [&](u64 i) -> val<SC_CTR_BITS, i64> {
      return val<SC_CTR_BITS, i64>{bias_w[i]};
    };
    auto sum = all_w.fo1().fold_add();

    // Prediction = sign of sum
    val<1> sc_pred = (sum < hard<0>{});

    // Confidence = |sum| > threshold
    val<1> sc_conf = (sum > hard<SC_CONF_THRESH>{}) |
                     (sum < hard<-SC_CONF_THRESH>{});

#ifdef TAGE_MONITOR
    stats.lookups++;
    if (static_cast<u64>(sc_conf)) stats.overrides++;
#endif

    return {sc_conf, sc_pred};
  }

  // ======== save_branch — called from update_condbr ========
  void save_branch([[maybe_unused]] val<64> &branch_pc,
                   [[maybe_unused]] u64 idx) {}

  // ======== train — called from update_cycle ========
  void train(arr<val<1>, FETCH_WIDTH> &branch_taken,
             [[maybe_unused]] arr<val<1>, FETCH_WIDTH> &is_branch,
             [[maybe_unused]] val<1> &mispredict,
             [[maybe_unused]] val<1> &correct_pred,
             val<1> &extra_cycle,
             u64 num_branch) {
    if (num_branch == 0 || !did_prefetch) {
      did_prefetch = false;
      return;
    }

    val<1> actual = branch_taken[0];
    actual.fanout(hard<NUM_BIAS + 1>{});

    // Update bias weights: increment toward actual direction
    // Only when extra_cycle fires (different HARCOM cycle from prefetch read)
    for (u64 i = 0; i < NUM_BIAS; i++) {
      val<SC_CTR_BITS> new_w = update_ctr(bias_w[i], actual);
      execute_if(extra_cycle, [&]() {
        bias[i].write(val<BIAS_LOG>{saved_bias_idx}, new_w);
      });
    }

#ifdef TAGE_MONITOR
#ifdef CHEATING_MODE
    if (static_cast<u64>(extra_cycle)) {
      stats.update_writes++;
      i64 sum = 0;
      for (u64 i = 0; i < NUM_BIAS; i++)
        sum += static_cast<i64>(bias_w[i]);
      bool sc_pred_v = (sum < 0);
      bool actual_v = static_cast<u64>(actual) & 1;
      bool tage_correct = !static_cast<u64>(mispredict);
      bool sc_correct = (sc_pred_v == actual_v);
      if (sc_correct && !tage_correct) stats.sc_correct++;
      if (!sc_correct && tage_correct) stats.sc_wrong++;
    }
#endif
#endif

    did_prefetch = false;
  }
};
