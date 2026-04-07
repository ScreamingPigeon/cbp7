# SC Overrider — Design & Performance Lessons

## Summary

A Statistical Corrector (SC) overrider was implemented following Seznec's TAGE-SC-L papers. The SC is structurally functional but does not improve VFS due to training instability and threshold tuning issues across workloads.

## Architecture

### Final Working Design
- **3 bias tables** (PC-only indexed, `ram<>`, 256 entries x 6-bit signed)
- **4 GEHL tables** (`ram<>`, 512 entries x 6-bit signed, history lengths 4/10/25/64)
- **Sum precomputed in predict1** (no tage_pred dependency on critical path)
- **Two-phase training**: `train_compute()` before `need_extra_cycle`, `train_write()` after
- **SC forces extra_cycle** via C++ bool OR'd into extra_cycle computation
- **Training direction**: `update_ctr(weight, tage_correct)` — weights learn whether TAGE was right/wrong
- **Override decision**: `candidate = (sum < 0) & (|sum| > threshold)`, prediction = `~tage_pred`

### Template Parameters (current defaults)
```
SCOverrider<NUM_BIAS=3, BIAS_LOG=8, SC_CTR_BITS=6, FETCH_WIDTH=16,
            NUM_GEHL=4, GEHL_LOG=9, GEHL_MINH=4, GEHL_MAXH=64,
            THRESH_INIT=200, THRESH_MARGIN=4>
```

## HARCOM Constraints Discovered

### 1. `execute_if` always evaluates the lambda body
HARCOM's `execute_if(cond, lambda)` evaluates the lambda for cost/transistor tracking regardless of `cond`. If the lambda contains `ram.write()`, the RAM access is counted even when `cond=0`. This causes "single RAM access per cycle" crashes when the same RAM was read earlier.

**Implication**: Cannot wrap `ram.write()` in `execute_if` to conditionally skip writes. TAGE's own writes work because their `execute_if` conditions (allocate, g_weak) inherently imply `extra_cycle=1` — the cycle has already advanced.

### 2. `rwram.write()` has the same issue
`rwram::write()` internally uses `execute_if` for bank selection → same problem. Even unconditional `rwram.write()` triggers the internal `bank[i].write()` inside `execute_if`, conflicting with earlier `bank[i].read()` from `rwram::read()`.

### 3. Working RAM write pattern
The only pattern that works for SC (reads in predict1, writes in update_cycle):
1. **Precompute write data** into staging regs (before `need_extra_cycle`)
2. **Force extra_cycle=1** by OR'ing a C++ bool into the extra_cycle computation
3. **Write unconditionally** from a C++ `if (should_write)` guard (not `execute_if`) after `need_extra_cycle`
4. Use **plain `ram<>`** (not `rwram`) — no internal `execute_if`

```cpp
// In Tage.hpp extra_cycle computation:
val<1> extra_cycle = base_extra_cycle | val<1>{overrider.should_write ? 1u : 0u};
need_extra_cycle(extra_cycle);

// After need_extra_cycle:
if (overrider.should_write) {
    overrider.train_write();  // plain ram.write(), no execute_if
}
```

### 4. val<> sign extension forbidden
HARCOM does not allow sign extension (e.g., `val<12, i64>{val<8, i64>}` fails). All arithmetic operands must match width or use `hard<>` constants. Two's complement negation (`-x`) widens by 1 bit via subtraction. Use `(x + threshold) < 0` instead of `x < -threshold` to avoid negation widening.

### 5. Centered sum issues
The paper's centered sum formula (`2*ctr + 1` per component + TAGE vote) requires careful width management. HARCOM addition always widens by 1 bit. Using `auto` throughout and only truncating when storing to regs avoids most issues.

## Performance Issues

### 1. P2 Latency (direction-indexed bias)
With bias tables indexed by `{PC, tage_pred}`, the `select(tage_pred, bias_t, bias_nt)` mux puts the entire `fold_add → threshold compare → override select` chain on the P2 critical path after TAGE's prediction arrives.

| Config | P2 Latency | L2 (ceil) |
|--------|-----------|-----------|
| Baseline (no SC) | 1.87 | 2 |
| SC with direction-indexed bias | 3.83 | 4 |
| SC with PC-only bias | 1.95 | 2 |

**Fix**: PC-only bias. Sum precomputed in predict1. Only the final `select(candidate, sc_pred, tage_pred)` gate depends on TAGE — adds ~40ps instead of ~600ps.

### 2. Training Instability — Weight Saturation
Initial approach: `update_ctr(weight, actual)` (train toward actual branch direction). All weights saturate to the same extreme (e.g., -32 for mostly-taken branches) because bias + GEHL tables all learn unconditional branch direction, not TAGE correction.

**Fix**: `update_ctr(weight, tage_correct)` — train toward whether TAGE was right/wrong. Weights center around 0 for well-predicted branches.

### 3. Threshold Runaway
Adaptive threshold (`increment on SC wrong, decrement on SC right`) runs away to max when SC is wrong more often than right (common during warmup). Once saturated, SC never overrides → threshold never decrements.

**Fix**: Only adapt threshold when `|sum| <= threshold` (sum near decision boundary). Or disable adaptive threshold and use fixed value.

### 4. Selective Training Gate
Seznec's formula: train when `misprediction OR |sum| <= threshold + margin`. With `tage_correct` training, `misprediction` refers to SC's own misprediction, not TAGE's. Using TAGE's misprediction signal causes over-training.

### 5. Extra Cycle Overhead
SC forces an extra HARCOM cycle every block with branches (for RAM writes). This adds ~770K extra cycles on gcc (879K vs 111K baseline), reducing IPC from ~4.3 to ~4.0.

### 6. Cross-Workload Regression
SC with THRESH=200 improves gcc (32K vs 37K mispredictions, -14%) but regresses many other traces. On 20-trace eval:
- **VFS: 0.562** (down from 0.857 baseline)
- Some traces severely hurt: dcapo-kafka 37.4 MPKI, sampleflow 37.0 MPKI

The `tage_correct` training approach causes weight saturation on biased workloads where TAGE is consistently right — weights all go positive, sum always exceeds threshold, SC never overrides (harmless but wasteful), OR when TAGE is consistently wrong for specific PCs, SC overrides correctly on those but incorrectly on others sharing the same hash.

## Key Metrics (gcc trace, 1M warmup, 5M measure)

| Config | Mispred | P2 lat | EPI | Extra cycles | VFS (20-trace) |
|--------|---------|--------|-----|-------------|----------------|
| Baseline (no SC) | 37,387 | 1.87 | ~1900 | 111K | 0.857 |
| SC direction-bias, GEHL=4, T=200 | 32,115 | 3.83 | 2046 | 879K | 0.573 |
| SC PC-only bias, GEHL=4, T=200 | 32,948 | 1.95 | 1837 | 880K | 0.562 |
| SC bias-only (no GEHL), T=80 | 267,314 | 3.21 | 1697 | 879K | — |

## Files

- `predictors/custom/SCOverrider.hpp` — SC implementation (current state: PC-only bias + 4 GEHL)
- `predictors/custom/Tage.hpp` — Integration: two-phase train, extra_cycle forcing, lookup signature
- `predictors/custom/TageOverrider.hpp` — Updated interface documentation

## Future Work

1. **Fix training**: The `tage_correct` approach doesn't generalize. Consider training toward `actual` but with proper centering (TAGE vote in sum) and dynamic PC-indexed thresholds.
2. **Reduce extra cycle overhead**: Only force extra_cycle when SC was actually confident and overrode, not every block.
3. **Dynamic threshold**: PC-indexed threshold table (32 entries) adapting per-branch, not global.
4. **Chooser**: Use TAGE confidence (newly_alloc, provider table depth) to gate SC overrides.
5. **History sharing**: Use TAGE's global history instead of separate SC history to avoid divergence.
