# Ahead TAGE Architecture

## Goal
Combine TAGE accuracy (MPKI ~5.9) with ahead pipelining (P2 < 1 cycle) and BPB>1 (low EPI). Target VFS > 0.95.

## Targets (from VFS Explorer)

| Metric | Target | Current TAGE | Priority |
|--------|--------|-------------|----------|
| P2 ceil | **1** | 2 | #1 — worth +0.077 VFS alone |
| MPKI | **≤ 5.75** | 5.9 | #2 — maintain, don't regress |
| Extra cycle ratio | **≤ 0.10** | ~0.10 | #3 — keep low, don't add overhead |
| EPI | **≤ 900 fJ** | 1211 | #4 — BPB>1 is the path here |

At these targets: **VFS ≈ 0.977**

## Reference Points

| Design | VFS | MPKI | EPI | P2 ceil |
|--------|-----|------|-----|---------|
| Current TAGE | 0.886 | 5.9 | 1211 | 2 |
| Current TAGE @ P2=1 | 0.963 | 5.9 | 1211 | 1 |
| Current TAGE @ P2=1 + EPI=900 | ~0.977 | 5.9 | 900 | 1 |
| gshareN_ahead_best | 0.922 | 8.8 | ~460 | 1 |

P2 latency is the single biggest lever. EPI reduction adds ~0.01. Accuracy is already good.

## Design Principles

1. **predict2 does ZERO RAM reads** — all reads happen in predict1 of the previous block (ahead pipelining)
2. **Branch-rank indexed (lanes), not instruction-offset indexed** — N=4 lanes covers >99% of blocks, vs FETCH_WIDTH=16 offset-indexed RAMs
3. **One TAGE entry covers multiple branches** — BPB>1 reduces aliasing pressure and per-branch read cost
4. **P1 = P2 = same TAGE-ahead prediction** — no separate gshare, no divergence penalty, no wasted P1 tables. Like gshareN_ahead_best, both levels return the same reg value.

## Simulator Interface Recap

```
predict1(PC)         — first instruction in block
reuse_predict1(PC)   — each subsequent instruction
predict2(PC)         — first instruction
reuse_predict2(PC)   — each subsequent instruction
update_condbr(PC, taken, next_pc)  — ONLY for conditional branches
update_cycle(info)   — end of block
```

- Predictor is called for EVERY instruction, not just branches
- Non-branch predictions are ignored by simulator
- num_branch (C++ int) tracks branch rank within block — free per forum §1
- reuse_prediction(1) extends the block; reuse_prediction(0) or taken branch ends it

## Timing Model (CRITICAL)

P1 and P2 are called with the **same starting timestamp**:
```cpp
p1_result = p.predict1({instruction.pc, time});
p2_result = p.predict2({instruction.pc, time});
```

Both latencies are measured from the same `time`. You CANNOT do work in predict1
and "use" it in predict2 to hide latency — they are independent timing chains
from the same origin. The only way to achieve P2 < 1 cycle is to have predict2's
output derived entirely from **reg values written in a previous cycle** — not from
any computation starting at the current `time`.

HARCOM timing: a `reg` read returns the value stored in a previous cycle. Its
timing is independent of the current `time`. So if predict2 reads from regs that
were written during a previous block's predict1, the latency from `time` is just
the combinational logic on those reg values (tag compare + mux).

## Ahead Pipeline

```
Block N-1 cycle:
  predict1(block_{N-1}):
    ├─ P1 gshare: read table → p1 (used for P1 prediction)
    ├─ Pipeline shift: ahead[0] → ahead[1]
    ├─ TAGE: compute indices for block_{N-1}
    ├─ TAGE: read all tables → ahead[0] (stored in regs)
    │   These reads are for block_{N-1}'s SUCCESSORS
    └─ Return: P1 prediction

  predict2(block_{N-1}):
    ├─ Use ahead[1] (regs written at block_{N-2}'s predict1)
    ├─ Tag comparison on reg values (COMBINATIONAL)
    ├─ Provider selection → prediction
    └─ Return: TAGE prediction

    P2 latency = tag compare + mux on reg values ≈ 0.5-0.8 cycles
    (No RAM reads, no dependence on current `time`)

Block N cycle:
  predict1(block_N):
    ├─ P1 gshare: read table → p1
    ├─ Pipeline shift: ahead[0] → ahead[1]
    │   (ahead[1] now contains block_{N-1}'s reads)
    ├─ TAGE: read at block_N's index → ahead[0]
    └─ Return: P1 prediction

  predict2(block_N):
    ├─ Use ahead[1] (regs from block_{N-1}'s predict1)
    ├─ Tag compare + mux (COMBINATIONAL from regs)
    └─ Return: TAGE prediction
```

Key insight: predict2 of block N uses data read during predict1 of block N-1.
Since those reads are stored in regs, predict2's timing starts from the reg
values (previous cycle), not from the current `time`. The only latency charged
to P2 is the combinational logic AFTER reading the regs.

This is exactly how gshareN_ahead_best achieves P2=P1: both predict1 and
predict2 return `pred[num_branch]` which is a reg written in a previous cycle's
predict1. The reg read has zero latency relative to current `time`.

## Consequence for the Ahead Ambiguity

When predict1(block_{N-1}) reads TAGE tables, it doesn't know what block_N
will be. So the reads must cover ALL possible successor blocks. This is resolved
by either:

1. **Banking**: read multiple banks simultaneously, select the right one in predict2
   based on how block_{N-1} actually exited (path known by predict2 time)
2. **Secondary tag**: store a small tag identifying which successor this entry was
   trained for. Check in predict2. Miss → fall back to bimodal.

The banking approach reads more data (higher EPI per predict1) but guarantees a
hit. The secondary tag approach reads less but has a miss rate requiring a
fallback read (bimodal) in predict2 — which adds latency to P2 on misses.

## gshareN_ahead_best Analysis

### RAM Structure
- `ctr_hi`: Single RAM, 16K entries, each entry = LANES × BANKS bits (e.g. 8×8=64 bits for N=7)
- `ctr_lo`: Per-lane hysteresis RAMs in UPDATE_ONLY zone (never read during prediction)
- LANES = next power of 2 ≥ N. Each lane = one branch rank. No aliasing within a block.
- BANKS = N+1 (rounded to power of 2). Each bank = predictions for one exit path of previous block.

### Pipeline
- predict1: reads ctr_hi at hash(PC, ghist) → stored in block_pred[0] for NEXT cycle
- predict1 also: shifts block_pred[0]→[1], selects bank from [1] using path^XB, reorders lanes
- predict2: returns pred[num_branch] — same reg, zero additional latency
- reuse_predict1/2: returns pred[num_branch] — just advances branch rank counter

### Bank Select
- `path`: set in update_cycle — encodes which branch was taken (exit path of prev block)
- `XB`: PC bits XORed with path for bank scrambling — spreads utilization evenly across banks
- `XL`: PC bits for lane scrambling — spreads branch rank mappings across lanes

### RAM Access Pattern
- predict1: ONE read of ctr_hi (for future block)
- update_cycle, no mispredict: writes ctr_lo only (UPDATE_ONLY zone, no conflict)
- update_cycle, mispredict: need_extra_cycle(1), then:
  - Read ctr_lo to check hysteresis (weak/strong)
  - Write ctr_hi to flip prediction (if weak)
  - Write ctr_lo to update hysteresis
- Read and write of ctr_hi never conflict: read is in predict1 cycle, write is in extra cycle

### Key Insight: Correct Predictions Cost Almost Nothing
- ctr_hi only written on misprediction (~1% of blocks)
- ctr_lo written every block but it's UPDATE_ONLY (no prediction reads)
- No extra cycles on correct predictions
- This is why EPI is so low

## TAGE Ahead Entry Format (N=7, 4 tables)

### Per-Table Entry (one bank)
```
[tag (11b)] [ctr_0..ctr_6 (7b)] [hyst_0..hyst_6 (7b)] [u (1b)]
= 26 bits per entry per bank
```
- Tag: shared across all N branches (one index = one tag)
- Counter: per-branch (each branch needs its own direction)
- Hysteresis: per-branch (per-entry too coarse — one mispredict destabilizes all)
- U-bit: shared per-entry (simplest, evicts all branches together)

### Banking
Each table banked N+1 = 8 ways (one per exit path of predecessor).
4 tables × 8 banks = 32 physical RAMs.
Each predict1 does 32 reads (all banks, all tables).
predict2 selects 1 bank per table based on actual path → 4 tag checks.

### Secondary Tag Alternative
Instead of 8 banks: 1 entry + 2-3 bit secondary tag identifying which exit path.
4 tables × 1 bank = 4 reads in predict1 (much lower EPI).
predict2 checks secondary tag. Miss → bimodal fallback (small RAM read in P2).
Trade: lower EPI but ~3% miss rate on secondary tag → bimodal accuracy on those blocks.

### RAM Conflict Strategy (same as gshare)
- predict1: reads TAGE tables (for future block)
- update_cycle correct: write hysteresis/u only (can use UPDATE_ONLY zone or rwram)
- update_cycle mispredict: need_extra_cycle, write tag + counter + hysteresis + u
- Key: prediction reads and update writes are separated by extra cycle on mispredict
- On correct predictions: minimize writes (silent update elimination) to avoid conflicts

## Open Questions

### Q1: How to index TAGE tables ahead?
Current TAGE indexes with: hash(PC_current, folded_history)
Ahead needs: index with PC of previous block, but predict for current block.

Options:
- **A)** Index with PC_prev, bank by exit path (like gshareN_ahead)
- **B)** Index with PC_prev, secondary tag encodes which successor
- **C)** Index with predicted next_PC (speculative)

gshareN_ahead uses option A with BANKS paths. For TAGE, option B (secondary tag from ahead-prediction paper) seems cleanest.

### Q2: How many lanes (N)?
- N=4 covers most blocks (rarely >4 conditional branches per block)
- N=7 matches gshareN_ahead_best
- More lanes = more prediction bits per entry = wider reads
- Need to check branch count distribution across our traces

### Q3: Banking for ahead?
gshareN_ahead uses BANKS=8 (one per exit path). For TAGE:
- 2 banks (fall-through vs taken) covers the common case
- More banks = more storage but better coverage
- Secondary tag (from ahead-prediction paper) can replace extra banks

Path bank is NOT a choice during allocation — it is determined by how the
predecessor block exited. Fall-through → bank 0, taken after branch K → bank K.
The allocator's only choice is which TABLE within that bank (longest table with
u=0 on the determined path bank).

### Q4: What history is available?
At predict1 time, the folded history reflects block_{-1}'s state (not yet updated with block_0's branches). This is fine — ahead predictions are trained with "block X was followed by block Y" relationship.

### Q5: Training / update?
- update_cycle knows the actual branch outcomes
- Need to write back to the entry that was read 1 cycle ago (index stored in pipeline reg)
- Path/bank determined by actual exit path (known at update time)

### Q5.5: Gshare→TAGE promotion on gshare misprediction
When TAGE has a secondary tag miss and gshare provides the prediction:
- Gshare correct → no action, pattern not worth TAGE capacity
- Gshare mispredicts → attempt to allocate a TAGE entry for this
  predecessor→successor pair, but ONLY if the candidate entry's usefulness
  is below a threshold. Don't evict useful TAGE entries just because gshare
  failed once. The secondary tag of the new entry is set to hash(actual next_pc).

### Q6: Fallback on secondary tag miss?
Small ahead-pipelined gshare as fallback. NOT bimodal — gshare captures
history correlations at the same storage cost.
- Read in predict1 alongside TAGE tables: index = hash(PC_N ^ short_dir_hist)
- Result stored in pipeline reg
- On secondary tag miss in predict2: use gshare pipeline reg (combinational mux)
- Tiny table (e.g. 1K-4K entries, 2-bit counters, no tags)
- Update on mispredict only (same policy as TAGE) — initialize saturated
- Used for ~3-5% of blocks (secondary tag miss rate)

### Q7: Per-lane vs per-entry counters/hysteresis/u-bits?
- **Per-lane prediction counter**: each branch gets its own direction bit (essential)
- **Shared hysteresis**: one hyst per entry, covers all lanes (saves storage)
- **Shared u-bit**: one useful bit per entry (saves storage, but evicts all branches together)
- **Shared tag**: one tag per entry (essential — one lookup per line)

## Entry Format (tentative, BPB=4)

```
[tag (11b)] [pred_0 (1b)] [pred_1 (1b)] [pred_2 (1b)] [pred_3 (1b)] [hyst (3b)] [u (1b)]

Total: 11 + 4 + 3 + 1 = 19 bits per entry
Current: 11 + 1 + 3 + 1 = 16 bits per entry (but need 4x entries for same coverage)

Net: 19 bits × E entries vs 16 bits × 4E entries → ~3x storage reduction for same coverage
```

## EPI Estimate

Current TAGE per block:
- 8 TAGE tables × 16 lanes = 128 RAM reads in predict2
- Plus bimodal, P1 pred/hyst = ~64 more reads
- Total: ~192 RAM reads per block

Ahead TAGE per block:
- 8 TAGE tables × 1 read = 8 RAM reads in predict1 (ahead)
- predict2: 0 RAM reads (combinational)
- P1: can also be ahead-pipelined
- Bimodal fallback: 1 read (small table)
- Total: ~9 RAM reads per block

Estimated EPI reduction: ~20x fewer reads → actual EPI maybe 3-5x lower (accounting for wider entries, update writes, etc.)

## Missing Functional Behavior

### Gshare fallback index
The fallback gshare must be ahead-pipelined. Index with `hash(PC_{N-1}, short_hist)` —
same ahead indirection as TAGE. Read in predict1 of block N-1, result in pipeline reg
for predict2 of block N. If indexed by PC_N directly, that's a RAM read in predict2
which breaks P2 < 1 cycle.

### Cold start (first block)
Pipeline regs are empty on the first block. No ahead data available.
Cold-start policy: predict not-taken (or static bias) until the pipeline fills.
This is 1 block of warmup — negligible accuracy impact.

### true_block handling
When a block ends due to misprediction of a not-taken branch (predicted taken but
actually not-taken), the next "block" continues the same line from where we left off.
The predecessor didn't change — so history should NOT be updated. gshare_ahead tracks
this with a `true_block` reg:
- true_block = 1 if the block ended normally (taken branch, line end, or correct pred)
- true_block = 0 if the block ended due to misprediction of a not-taken branch
- Only update history, shift pipeline regs, and issue new RAM reads when true_block = 1
- When true_block = 0, reuse the same pipeline regs (same predecessor, same reads)

### num_branch overflow
If a block has more than N conditional branches, end the block at N via
`reuse_prediction(0)`. Like gshare_ahead: `if (num_branch == N) reuse_prediction(0)`.
The block is split at the Nth branch; the remaining instructions form a new block
in the next cycle.

## Implementation Notes

### Do NOT use common.hpp
Implement all components from scratch. We want full control over:
- Global history register (shift width, what bits go in, per-branch vs per-block)
- Folded history (fold widths, incremental update logic)
- Saturating counter update
- Any RAM banking / rwram logic

common.hpp's `geometric_folds`, `global_history`, `folded_gh`, `rwram`, `update_ctr`
are all black boxes with hidden fanout/timing assumptions. For the ahead design we
need to understand and control every gate on the critical path.

### Secondary Tag Computation
- Computed in `update_cycle(block N)` when actual `next_pc` is known
- `secondary_tag = hash(next_pc[low bits])`
- Written into the TAGE entry at the index computed from PC_{N-1} + history
- During prediction of block N+1: check `hash(PC_{N+1}) == stored_secondary_tag`
- Match → use prediction from pipeline regs. Miss → gshare fallback.
- TAGE entry selection requires BOTH: primary_tag_hit AND secondary_tag_hit.
  Primary tag = right predecessor + history. Secondary tag = right successor.
  If primary hits but secondary misses → wrong successor → gshare fallback.

### Pipeline Reg Shift Ordering
HARCOM regs: one write per cycle, can't read and write same reg in same cycle.
In predict1, the shift must happen in this order:
1. Write ahead[1] from ahead[0]  (shift — reads ahead[0], writes ahead[1])
2. Write ahead[0] from RAM read  (overwrites ahead[0] with new data)

This works because step 1 reads ahead[0] before step 2 overwrites it.
Both are writes to different regs (ahead[1] and ahead[0]) so no conflict.
This is exactly what gshare_ahead does:
```
block_pred[1] = block_pred[0];          // shift
block_pred[0] = ctr_hi.read(index[0]);  // new read
```

### Simulator Cycle Model (CRITICAL)
predict1, predict2, update_condbr (×N branches), and update_cycle ALL run in
the SAME cycle. `panel.next_cycle()` advances AFTER update_cycle returns.
- RAM reads (predict1) and RAM writes (update_cycle) are in the same cycle
- `need_extra_cycle(1)` grants a second cycle for writes that conflict with reads
- Reg writes in update_cycle take effect next cycle — so predict1 of the next
  block sees updated dir_hist and path_hist.

### Update Policy: Write-on-Mispredict-Only
RAM read (predict1) and RAM write (update_cycle) are in the same cycle — different
indices but same physical RAM. HARCOM enforces one access per RAM per cycle, so
writes require `need_extra_cycle`.

Options:
- **Always update** (current TAGE): strengthen counters + set u-bits on correct
  predictions. Requires extra cycle every block. Extra_cycle ratio ≈ 1.0 (bad).
- **Write on mispredict only** (gshare_ahead style): only update counter/u/hyst
  when prediction was wrong. Zero extra cycles on correct predictions.
  Extra_cycle ratio ≈ misprediction rate ≈ 0.6% (excellent).

To make write-on-mispredict-only work for u-bits:
- Initialize new entries with **saturated u-bits** (useful by default)
- Decrement u on misprediction or when entry is useless
- No need to "set u when provider correct and alt wrong" — avoids a write

For counters/hysteresis:
- Initialize new entries at mid-confidence
- Only flip prediction / weaken hysteresis on misprediction
- Accept slower training for much lower extra_cycle and EPI
- gshare_ahead validates this approach works in practice

### History Register Design Decisions
Current TAGE: gfolds updates once per block with 6 bits of next_pc.
For ahead TAGE we need to decide:
- How many bits per update? (PATHBITS=6 in reference)
- Per-block or per-branch updates? (per-block in reference, but per-branch
  gives more information at the cost of more reg writes)
- What information? (next_pc bits, branch direction, or both)
- How does staleness interact with the secondary tag?

## Implementation Plan

1. Define the ahead TAGE entry format and storage
2. Implement own history register + folded history (no common.hpp)
3. Implement predict1: ahead reads + pipeline shift
4. Implement predict2: combinational tag match + bank/secondary tag select + provider selection
5. Implement update_cycle: write-back to stored indices, history update
6. Test on gcc trace, measure P2 latency and accuracy
7. Iterate on banking, lane count, tag width, secondary tag width
