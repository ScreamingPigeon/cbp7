# HARCOM Timing Optimization Guide

How to systematically reduce P2 latency in HARCOM-based branch predictors.

## The Timing Model

Every `val` in HARCOM carries a timestamp. The timestamp of an operation's
result = max(input timestamps) + operation delay. Reg timestamps carry forward
from write time — they do NOT reset at cycle boundaries.

Three components contribute to P2 latency:

| Component | What it is | How to measure |
|-----------|-----------|----------------|
| **Logic** | Gate delays (AND, OR, XOR, MUX, comparators, one_hot, fold_or) | Build with `-DFREE_FANOUT -DFREE_WIRING` |
| **Fanout** | Buffer tree delays for driving multiple readers of a value | `baseline - FREE_FANOUT_result` |
| **Wiring** | Manhattan-distance wire delay between RAM locations | `baseline - FREE_WIRING_result` |

Measure all three before optimizing. Attack the largest contributor first.


## Analysis Steps

### 1. Measure the breakdown

```bash
# Baseline
make cbp && ./cbp trace.gz name 1000000 40000000

# Remove fanout cost
g++ ... -DFREE_FANOUT -DPREDICTOR='...' && ./cbp ...

# Remove wiring cost
g++ ... -DFREE_WIRING -DPREDICTOR='...' && ./cbp ...

# Remove both (pure logic floor)
g++ ... -DFREE_FANOUT -DFREE_WIRING -DPREDICTOR='...' && ./cbp ...
```

Compute:
- Logic = both_free
- Fanout only = FREE_WIRING - both_free
- Wiring only = FREE_FANOUT - both_free

### 2. Build and run the timing probe

```bash
make timing-probe
./build/timing-probe-* trace.gz name 1000000 40000000
```

This dumps per-signal absolute timestamps. Sort by timestamp to find the
critical path. Look at the **gaps** between consecutive signals — the largest
gap is the bottleneck.

### 3. Identify the critical path

For a TAGE ahead-pipelined predictor, the P2 critical path is typically:

```
predict1: fold.get() → XOR → idx → RAM reads → prefetch_* writes
          ↓ (transparent latch — same cycle)
update_cycle: prefetch_* reads → full_hits → match → one_hot →
              provider/alt extraction → use_alt → select → pred[] write
          ↓
predict2: return pred[num_branch]
```

Every named `val` on this path adds fanout delay. Every reg read on this
path adds transparent-latch delay if written in the same cycle.


## Fanout Optimization

Fanout is typically the dominant cost (40-55% of P2). Three tools control it:

### fo1() — zero-cost single read

```cpp
auto result = some_val.fo1();  // 0 FO2 delay, single read only
```

Rules:
- **Must be called on an lvalue** (named variable, not a temporary)
- The returned value must be consumed **exactly once**
- Use on every single-read intermediate on the critical path

Pattern — name the intermediate, then fo1():
```cpp
auto remainder = match ^ match1;
val<MATCH_BITS> match2 = remainder.fo1().one_hot();
```

### fanout(hard<N>) — optimized buffer tree for N reads

```cpp
some_val.fanout(hard<N>{});  // log2(N) FO2 stages instead of N × FO2
```

Rules:
- N must be the **exact** number of reads (over-declaring wastes stages)
- `fanout(hard<1>)` = same as fo1()
- `fanout(hard<2>)` = 1 FO2 stage = same as default (no benefit)
- **Only helps for N ≥ 3** (log2(3)=2 stages vs 3 individual reads)
- Count ALL reads: inside current function + after return (via move) + training saves

### Reducing fanout by restructuring

When a wide value (e.g., PRED_BITS) is read N times via bit-shift:
```cpp
// BAD: provider_pred read N times → fanout(hard<N+2>) = 4 stages for N=8
static_loop<N>([&]<u64 I>() {
    pred[I] = select(use_alt, val<1>{alt_pred >> I}, val<1>{provider_pred >> I});
});
```

Split into per-bit array first, then fo1() each bit:
```cpp
// GOOD: provider_pred read once (make_array) → fanout(hard<3>), each bit fo1()
arr<val<1>, PRED_BITS> pp_bits = provider_pred.make_array(val<1>{});
arr<val<1>, PRED_BITS> ap_bits = alt_pred.make_array(val<1>{});
static_loop<N>([&]<u64 I>() {
    pred[I] = select(use_alt, ap_bits[I].fo1(), pp_bits[I].fo1());
});
```

This converts one high-fanout read into one low-fanout read + N zero-cost reads.


## Transparent-Latch Bypass

When predict1() and update_cycle() run in the same HARCOM cycle, any reg
written in predict1 and read in update_cycle goes through a transparent latch
(adds delay). If update_cycle also shifts that value to another reg and reads
the second reg, that's **two** latch crossings.

### The pattern

```cpp
// predict1 writes prefetch_pred[I] (latch 1)
// update_cycle:
current_pred[I] = prefetch_pred[I];   // reads prefetch (latch 1), writes current
// ... resolution reads current_pred[I] (latch 2) ...
```

### The fix

Read `prefetch_*` directly in the resolution chain. Still shift for next cycle,
but the resolution doesn't go through the second latch:

```cpp
// Resolution reads prefetch_* directly (latch 1 only)
arr<val<PRED_BITS>, NT + 1> table_preds = [&](u64 i) -> val<PRED_BITS> {
    if (i < NT) return val<PRED_BITS>{prefetch_pred[i]};  // not current_pred
    return val<PRED_BITS>{prefetch_fb};
};
// Shift still happens (needed for next cycle's train_* save)
current_pred[I] = prefetch_pred[I];
```

Note: this may or may not help depending on whether the HARCOM model stacks
latch penalties. Measure before and after.


## Wiring Optimization

See `docs/floorplan.md` for zone, region, connect, and distribute techniques.
In practice, once fanout is optimized, wiring is often <10% of P2. Focus on
fanout first.


## Reference: gshareN_ahead_best Patterns

The `gshareN_ahead_best` predictor achieves excellent timing through discipline:

1. **fo1() on every single-read value** — register-to-register copies, lambda
   intermediates, RAM read results used once
2. **Precisely-sized fanout()** — every count matches exact consumer count
3. **No unnecessary reg hops** — values flow directly where needed
4. **No zone/region/connect/distribute** — clean fo1/fanout is sufficient

When optimizing a new predictor, audit every named value on the critical path:
- Read once → fo1()
- Read 3+ times → fanout(hard<N>) with exact N
- Read 2 times → no action needed (default FO2 = fanout(hard<2>))


## Optimization Checklist

1. Measure breakdown (FREE_FANOUT, FREE_WIRING, both)
2. Run timing probe to find critical path and gaps
3. For each named val on the critical path:
   - Count its reads (ctrl-F for variable name)
   - 1 read: add fo1()
   - 2 reads: leave as-is
   - 3+ reads: add fanout(hard<N>) with exact count
4. Look for high-fanout values that can be split (make_array + fo1 pattern)
5. Check for unnecessary reg hops (can you read an earlier reg directly?)
6. Re-measure. Repeat from step 2 if fanout is still dominant.
7. Only after fanout is minimized: consider floorplan (zone/region/connect)
