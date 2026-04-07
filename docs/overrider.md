# TageOverrider — Pluggable Prediction Override

## Overview

The TAGE predictor supports an optional overrider component that can override P2 predictions. When no overrider is used (`NoOverrider`), all override code compiles away — zero cost.

Two overriders are implemented:
- **SCOverrider** — Statistical Corrector (perceptron-based, bias tables). Current default.
- **LoopPredictorA** — Loop predictor (detects fixed-iteration loops).

## Quick Start

```cpp
// SC overrider (default):
using MyPredictor = Tage<>;  // uses SCOverrider<8, 6, FETCH_WIDTH, 4>

// Loop predictor:
using MyPredictor = Tage<SweepTableConfig<>, DefaultAllocConfig,
    16, 4096, 1, 1, false, true, true, true, false, 8, 2, DecayMild,
    DefaultResetFn, false, true, 4096, 6, true, 4, 2, 0, false, 27, 6,
    LoopPredictorA<64, 10, 10, 3, 16, 2>>;

// No overrider:
// Change last template param to NoOverrider
```

## Architecture — 4-Phase Interface

```
predict1:     overrider.prefetch(block_pc, raw_pc)
                → Read tables at line-level index
                → Extract weights/tags into regs

predict2:     overrider.lookup(inst_pc)
reuse_p2:       → Pure COMBINATIONAL (no RAM access)
                → Tag match / sum computation on pre-extracted regs
                → Returns {candidate, pred}
                → TageImpl: p2 = select(candidate, pred, tage_pred)

update_condbr: overrider.save_branch(branch_pc, idx)
                → Save branch info for training (plain C++)

update_cycle:  overrider.train(branch_taken, is_branch, mispredict,
                               correct_pred, extra_cycle, num_branch)
                → Update table weights
                → ram<> writes gated by execute_if(extra_cycle, ...)
```

## Writing a New Overrider

```cpp
template <u64 TABLE_SIZE = 64, /* your params */, u64 FETCH_WIDTH = 16>
struct MyOverrider {
    static constexpr bool ENABLED = true;

    // --- Storage ---
    hcm::ram<val<BITS>, SIZE> my_table{"name"};  // use ram<>, NOT rwram<>
    arr<reg<BITS>, N> cached_weights;              // pre-extracted in prefetch

    // --- Stats (for TageMonitor compatibility) ---
#ifdef TAGE_MONITOR
    struct LoopStats {  // must be named LoopStats for TageImpl compatibility
        u64 lookups = 0;
        u64 hits = 0;
        u64 overrides = 0;
        u64 alloc_writes = 0;
        u64 update_writes = 0;
        u64 prefetch_conf_nonzero = 0;
        u64 prefetch_num_nonzero = 0;
        // Add your own fields here
    } stats;
#endif

    struct LookupResult { val<1> candidate; val<1> pred; };

    // Phase 1: prefetch — called from predict1
    void prefetch(val<64> &block_pc, u64 raw_pc) {
        // Read RAM tables → store results in regs
        // Save indices for training (plain C++)
    }

    // Phase 2: lookup — called from predict2 AND reuse_predict2
    LookupResult lookup(val<64> &inst_pc) {
        // PURE COMBINATIONAL — no RAM reads, no reg writes
        // Compute prediction from pre-extracted regs
        return {confident, prediction};
    }

    // Phase 3: save_branch — called from update_condbr
    void save_branch(val<64> &branch_pc, u64 idx) {
        // Save branch tag/info for training
    }

    // Phase 4: train — called from update_cycle
    void train(arr<val<1>, FETCH_WIDTH> &branch_taken,
               arr<val<1>, FETCH_WIDTH> &is_branch,
               val<1> &mispredict,
               val<1> &correct_pred,
               val<1> &extra_cycle,
               u64 num_branch) {
        // Update table weights, gated by extra_cycle
        execute_if(extra_cycle, [&]() {
            my_table.write(index, new_value);  // ram write in different cycle from read
        });
    }
};
```

## Critical Rules

### 1. Use ram<>, NOT rwram<>, for overrider tables
`ram<>` with `execute_if(extra_cycle, ...)` for writes is the proven pattern.
Read in predict1 (cycle N), write in update_cycle (cycle N+1 after `need_extra_cycle`).
`rwram<>` causes bank conflicts under high write pressure.

### 2. lookup() must be PURE COMBINATIONAL
No RAM reads, no reg writes. Only read from pre-extracted regs filled in prefetch.
This is called from both predict2 AND reuse_predict2 — must work in any cycle.

### 3. Pass HARCOM types by reference
```cpp
void prefetch(val<64> &pc, u64 raw_pc);  // ✓ val by reference
void prefetch(val<64> pc, u64 raw_pc);   // ✗ val by value = 50% latency penalty
```

### 4. Single p2 write — no modify-after-write
```cpp
// ✗ WRONG: writes p2 twice
p2 = tage_prediction;
p2 = select(override, ovr_pred, p2);

// ✓ CORRECT: single expression
p2 = arr<val<1>, FW>{[&](u64 offset) {
       val<1> tage = ...;
       if (offset == 0) return select(ovr.candidate, ovr.pred, tage);
       return tage;
     }}.concat();
```

### 5. val<> cannot be reassigned
Use IIFE, arr::fold, or select chains instead:
```cpp
arr<val<BITS, i64>, N> weights = [&](u64 i) { return ...; };
auto sum = weights.fo1().fold_add();  // ✓ clean accumulation
```

### 6. Fanout budget for extra_cycle
The overrider's `execute_if(extra_cycle, ...)` calls consume extra_cycle fanout.
TageImpl declares: `extra_cycle.fanout(hard<NUM_TABLES * 2 + 1 + OVR * 5>{})`.
If your overrider needs more, increase the `OVR * 5` multiplier.

### 7. CHEATING_MODE + ram<> writes require -DREAD_WRITE_RAM
In CHEATING_MODE (monitor build), `execute_if` evaluates differently, causing
`ram<>` read+write conflicts. Add `-DREAD_WRITE_RAM` to monitor build flags.
Normal builds work fine without it.

### 8. No destructors on overrider structs
HARCOM types (reg, ram, arr) cannot be destroyed. Don't add `~MyOverrider()`.
Print stats by having TageImpl copy `overrider.stats` fields to TageMonitor.

## Existing Overriders

### SCOverrider (predictors/custom/SCOverrider.hpp)
- 3 bias tables (ram<>, 256 entries × 6-bit signed) + 4 GEHL tables (ram<>, 512 entries × 6-bit)
- Prediction = sign(sum of weights). Confidence = |sum| > threshold.
- Template: `SCOverrider<NUM_BIAS=3, BIAS_LOG=8, SC_CTR_BITS=6, FETCH_WIDTH=16, NUM_GEHL=4, GEHL_LOG=9, GEHL_MINH=4, GEHL_MAXH=64, THRESH_INIT=200, THRESH_MARGIN=4>`
- On gcc: 4.31% mispred (vs 5.08% without SC). 11:1 correct/wrong ratio.

### LoopPredictorA (predictors/custom/LoopPredictorA.hpp)
- Line-level RAM index + branch-level tag matching
- 2-way set-associative, 64 entries, rwram storage
- On NAMD: 97% override accuracy, 3.6% misprediction reduction
- Template: `LoopPredictorA<TABLE_SIZE=64, TAG_BITS=10, ITER_BITS=10, CONF_BITS=3, FETCH_WIDTH=16, LWAYS=2>`

## Files

- `predictors/custom/TageOverrider.hpp` — NoOverrider + concept docs
- `predictors/custom/SCOverrider.hpp` — Statistical Corrector
- `predictors/custom/LoopPredictorA.hpp` — Loop predictor
- `predictors/custom/Tage.hpp` — Call sites (search for `Overrider::ENABLED`)
- `docs/research_todos.md` — SC improvement roadmap (GEHL, dynamic threshold, etc.)
