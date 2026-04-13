# LINEINST vs P2 Latency Tradeoff

## Summary

`LINEINST` (the max instructions per fetch block) directly controls P2 latency
through the `block_entry` register width. Larger LINEINST = more IPB but higher
P2 latency. **Only LINEINST <= 32 fits in ceil=2 P2 cycles.**

## Results (test10 series, 2026-04-13)

Test predictor: full TAGE with 8 tables, meta/altsel, gshare P1, block machinery.
Trace: int_11, 4000 instructions. Baseline (no blocking): P2=1.58.

| LINEINST | reg bits | IPB (N=16) | P2     | ceil | delta from baseline |
|----------|----------|------------|--------|------|---------------------|
| 16       | 4        | 6.8        | 1.76   | 2    | +0.18               |
| 32       | 5        | 8.6        | 1.98   | 2    | +0.40               |
| 64       | 6        | 9.8        | 2.40   | 3    | +0.82               |
| 128      | 7        | 10.8       | 2.98   | 3    | +1.40               |
| 256      | 8        | 11.4       | 3.12   | 4    | +1.54 (saturated)   |
| 512      | 9        | 11.7       | 3.12   | 4    | +1.54 (saturated)   |
| 1024     | 10       | 12.0       | 3.12   | 4    | +1.54 (saturated)   |

IPB values from ahead_block_analyzer with N=16. P2 saturates at 3.12 for
reg >= 8 bits (HARCOM cycle floor).

## Root Cause

HARCOM measures P2 as the **max timing across ALL cycle outputs**, including
`reuse_prediction()` calls in predict1/reuse_predict1 — not just predict2's
return path. The `line_end()` computation reads `block_entry` (a reg written
from `inst_pc >> 2` in predict1), and this timing chain inflates P2.

The bottleneck is the **reg read timing** (fanout tree delay), not the
downstream circuit. All three circuit approaches give identical P2 at the
same LINEINST:

| Method                          | LINEINST=32 P2 | LINEINST=64 P2 |
|---------------------------------|----------------|----------------|
| binary: `(entry+size) >= LI`   | 1.98           | 2.40           |
| cmp-only: `entry >= (LI-size)` | 1.98           | 2.40           |
| ROM + free shift                | 1.98           | 2.40           |
| one-hot + shift                 | 2.40*          | 2.40           |

*One-hot at 32 bits is worse because the 32-bit reg has higher fanout cost
than a 5-bit binary reg.

## Encoding: One-hot vs Binary

One-hot encoding does NOT help — it makes things worse or equal:
- LINEINST=64: one-hot (64-bit reg) = binary (6-bit reg) = P2=2.40
- LINEINST=32: one-hot (32-bit reg) = P2=2.40, binary (5-bit reg) = P2=1.98

The one-hot shift `block_entry >> (LI-1-block_size)` is free (hard shift count),
but the wider register has higher fanout cost that dominates.

Note: `reg<N>` max is N=64 (HARCOM limitation), so one-hot LINEINST=1024 is
impossible.

## Approaches That Don't Help

1. **Precomputed remaining reg** (test9): `remaining = LINEINST - block_entry`
   in predict1, then `block_size >= remaining` in line_end. Same P2 because
   the reg read timing traces back through the same inst_pc chain.

2. **Removing reuse_prediction from P2 path** (test9): Even removing it from
   predict2 AND reuse_predict2 doesn't help — predict1/reuse_predict1's
   line_end() still inflates P2 via HARCOM's cycle-wide max.

3. **Removing reuse_prediction from everywhere except update_condbr** (test9):
   Still inflated P2 because update_condbr's line_end() contributes too.

4. **Algebraic rearrangement** (test10 L6): Eliminating the adder by rewriting
   `(entry + size) >= LI` as `entry >= (LI - size)` — no improvement because
   HARCOM's comparator internally uses an adder anyway.

5. **ROM lookup** (test10 L8): `rom<val<LI>, LI>` indexed by block_entry,
   then free shift by block_size — same P2 as comparator.

## HARCOM Free Operators (from harcom.pdf Table 3.4)

These have zero hardware cost and can be exploited:
- `<<` shift left with hard count
- `>>` unsigned shift right with hard count
- `&` and `|` with a hard value

All comparisons (`>=`, `==`, `!=`) always have hardware cost.

## Design Recommendation

For a two-level predictor (P1+P2) with TAGE:

- **LINEINST=32** (binary `reg<5>`): P2=1.98, IPB=8.6. Tight fit in ceil=2
  but gives good throughput. Matches the sweet spot for N=8 (Br/Blk ≈ 1.1).
- **LINEINST=16** (binary `reg<4>`): P2=1.76, IPB=6.8. Comfortable ceil=2
  margin. This is what tage.hpp uses (one-hot, but binary is equivalent).
- **LINEINST=64+**: ceil=3 minimum. Only viable if P2 latency budget allows it.

The IPB ceiling of ~8.6 (LINEINST=32) or ~6.8 (LINEINST=16) at ceil=2 is a
fundamental design limitation of HARCOM's timing model for two-level predictors
that use `reuse_prediction()` with `line_end()`.

## Reference Predictors

- **tage.hpp**: LINEINST=16, one-hot `reg<16>`, P2=1.86
- **gshareN_ahead_best**: LINEINST=1024, binary `reg<10>`, but single-level
  (no reuse_prediction in P2 path), so line_end doesn't inflate P2
- **TageDirect** (before fix): LINEINST=1024, binary `reg<10>`, two-level
  with reuse_prediction everywhere → P2 inflated to ~3.12
