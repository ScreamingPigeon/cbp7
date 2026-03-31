#pragma once
// TageOverrider — pluggable prediction override interface for TAGE.
//
// After TAGE computes its P2 prediction, an optional overrider can examine
// the prediction + confidence signals and decide to flip it. The overrider
// runs its own table lookup in PARALLEL with TAGE's RAM reads, then makes
// its override decision after TAGE signals are available.
//
// When Overrider = NoOverrider (default), all override code compiles away
// via if constexpr (Overrider::ENABLED). Zero cost.
//
// === Overrider 4-Phase Interface ===
//
// An overrider type must provide:
//
//   static constexpr bool ENABLED = true;
//   static constexpr u64 EXTRA_CYCLE_BUDGET = N;  // extra_cycle fanout consumers
//
//   struct LookupResult { val<1> candidate; val<1> pred; };
//
//   // Phase 1: Prefetch (called from predict1)
//   // Read RAM tables at line-level index, extract into regs.
//   void prefetch(val<64> &block_pc, u64 raw_pc);
//
//   // Phase 2: Lookup (called from predict2 AND reuse_predict2)
//   // PURE COMBINATIONAL — no RAM reads, no reg writes (except on first call).
//   // tage_pred: TAGE's prediction for this branch
//   // tage_low_conf: 1 if TAGE provider is newly allocated (low confidence)
//   // Returns {candidate, pred} where candidate=1 means override is active.
//   LookupResult lookup(val<64> &inst_pc, val<1> &tage_pred, val<1> &tage_low_conf);
//
//   // Phase 3: Save branch info (called from update_condbr)
//   // Save branch info for training (plain C++).
//   void save_branch(val<64> &branch_pc, u64 idx);
//
//   // Phase 4: Train (called from update_cycle)
//   // Update table weights, gated by execute_if(extra_cycle, ...).
//   void train(arr<val<1>, FETCH_WIDTH> &branch_taken,
//              arr<val<1>, FETCH_WIDTH> &is_branch,
//              val<1> &mispredict,
//              val<1> &correct_pred,
//              val<1> &extra_cycle,
//              u64 num_branch);
//
// === Critical Rules ===
//
// 1. Use ram<>, NOT rwram<>, for overrider tables.
//    Read in predict1 (cycle N), write in update_cycle (cycle N+1 via extra_cycle).
//
// 2. lookup() must be PURE COMBINATIONAL.
//    No RAM reads, no reg writes (except first-call saves guarded by bool).
//    Called from both predict2 AND reuse_predict2.
//
// 3. Pass HARCOM types by reference, never by value.
//    reg<N> by value CRASHES. val<N> by value adds ~50% latency.
//
// 4. No destructors on overrider structs.
//    HARCOM types can't be destroyed. Use print_stats() instead.
//
// 5. extra_cycle fanout: declare EXTRA_CYCLE_BUDGET to cover all
//    execute_if consumers in train().

struct NoOverrider {
  static constexpr bool ENABLED = false;
  static constexpr unsigned long long EXTRA_CYCLE_BUDGET = 0;
  // No methods needed — all call sites behind if constexpr (Overrider::ENABLED)
};
