# LINEINST vs P2 Latency Tradeoff

## Summary

`LINEINST` (the max instructions per fetch block) directly controls P2 latency
through the `block_entry` register width. Larger LINEINST = more IPB but higher
P2 latency. **Only LINEINST <= 32 fits in ceil=2 P2 cycles.**

## Results (test_lineinst, 2026-04-13)

Test predictor: full TAGE with 8 tables, meta/altsel, gshare P1, block machinery.
Trace: gcc (502-gcc-all_16112), 1M instructions. Baseline (no blocking): P2=1.58.

| LINEINST | Blocks  | P2     | ceil | delta from baseline |
|----------|---------|--------|------|---------------------|
| 16       | 230,734 | 1.76   | 2    | +0.18               |
| 32       | 213,797 | 1.98   | 2    | +0.40               |
| 64       | 204,224 | 2.34   | 3    | +0.76               |
| 128      | 200,979 | 2.47   | 3    | +0.89               |
| 256      | 197,148 | 2.47   | 3    | +0.89 (saturated)   |
| 1024     | 196,461 | 2.47   | 3    | +0.89 (saturated)   |

P2 saturates at LI=128 — the longest block in the gcc trace is ~128 instructions.

## Root Cause: Excess Fanout Reads in reuse_predict2

The P2 inflation is caused by **HARCOM's fanout exhaustion penalty** on regs
read during `reuse_predict2()`, NOT by the `line_end()` computation (which
is never on the P2 critical path).

### Mechanism (harcom.hpp:3769-3796)

1. `predict2()` writes regs (`pred1_tage`, `pred2_tage`, `altsel_bits`, etc.)
   and calls `.fanout(hard<2>{})` giving each reg 2 read credits.
2. `predict2()` itself consumes 1 credit (return value). First `reuse_predict2()`
   consumes the 2nd. Credits now exhausted.
3. Every subsequent `reuse_predict2()` call triggers the **excess read penalty**:
   HARCOM adds a FO2 inverter delay per read (`set_time(time() + fo2inv.delay())`).
4. Longer blocks (larger LINEINST) → more `reuse_predict2()` calls → more excess
   reads → higher accumulated delay → higher max P2.

### Controlled experiment (test_lineinst.hpp, gcc trace, 1M instructions)

| Mode | Description                              | LI=1024 P2 |
|------|------------------------------------------|------------|
| 0    | `reuse_prediction(~line_end())`          | 2.47       |
| 1    | no reuse_prediction at all               | 1.58       |
| 2    | compute line_end(), sink it, pass hard<1> | 2.47      |
| 3    | `reuse_prediction(hard<1>{})`            | 2.47       |

MODE 1→3: calling reuse_prediction alone inflates P2 (even with zero-cost arg).
MODE 2→0: line_end() computation adds nothing to P2. The inflation is entirely
from excess reg reads in reuse_predict2.

### Why increasing fanout doesn't help

| REUSE_FO | P2 (LI=1024) |
|----------|--------------|
| 2        | 2.47         |
| 4        | 2.48         |
| 8        | 2.49         |
| 16       | 2.50         |
| 32       | 2.51         |
| 64       | 2.52         |

`fanout(hard<FO>{})` adds logarithmic upfront delay (repeater tree) that applies
to ALL blocks. Most blocks are short (~5 instructions), so the upfront cost
exceeds the saved excess-read penalty. No single fanout value is optimal.

## Approaches That Don't Help

1. **Increasing fanout on reuse-path regs**: Upfront repeater tree cost exceeds
   saved excess-read penalty (see table above).

2. **Changing line_end() circuit**: Algebraic rearrangement, ROM lookup, one-hot
   encoding — all irrelevant because line_end() is not on the P2 path.

3. **One-hot block_entry encoding**: Wider register has higher fanout cost for
   the regs that ARE on the P2 path, making things worse.

## HARCOM Free Operators (from harcom.pdf Table 3.4)

These have zero hardware cost and can be exploited:
- `<<` shift left with hard count
- `>>` unsigned shift right with hard count
- `&` and `|` with a hard value

All comparisons (`>=`, `==`, `!=`) always have hardware cost.

## Design Recommendation

For a two-level predictor (P1+P2) with TAGE:

- **LINEINST=32** (binary `reg<5>`): P2=1.98, ceil=2. Tight fit but viable.
- **LINEINST=16** (binary `reg<4>`): P2=1.76, ceil=2. Comfortable margin.
- **LINEINST=64+**: ceil=3 minimum due to excess fanout reads.

The ceil=2 limit comes from reuse_predict2 fanout exhaustion, NOT from
line_end() arithmetic. Potential architectural fixes:
- Pre-shift altsel_bits into per-rank single-bit regs (avoids re-reading wide reg)
- Store reuse_predict2 results in a shift register written once in predict2
- Accept ceil=3 and use LINEINST=64+ for better IPB

## Reference Predictors

- **tage.hpp**: LINEINST=16, one-hot `reg<16>`, P2=1.86
- **gshareN_ahead_best**: LINEINST=1024, binary `reg<10>`, single-level
  (predict2 returns `pred[num_branch]` from reg — no excess reads because
  predict2 is trivial and pred[] is re-fanned in predict1)
- **TageDirect**: LINEINST=1024, two-level with reuse_prediction → P2 inflated
  to ~2.47 due to excess fanout reads in reuse_predict2



# What if prakhar notes

what if we still advance cycle but skip re-reading rams?
