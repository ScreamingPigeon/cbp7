# 2-Way Sec Tag Banking (Ahead History Disambiguation)

## Problem

TageAhead predicts block B during predict1(A), but the history used for
B's table lookup is stale — it's missing block A's outcome. The history
only updates in update_cycle(A), by which time B's prefetch is already
in registers. Result: TAGE tag match rate is low (~30% of branches
served by TAGE, 70% fall back to bimodal at 60% accuracy).

With N=7, a block has at most 8 possible exits (7 taken branches +
fall-through). Each exit produces a unique next_pc and thus a unique
sec tag. The stale index is always the same for a given block — it's
determined by the fold computed from history through block A-1. What
changes between different resolutions of A is which *entry* at that
index is correct, and that's exactly what the sec tag captures.

## Design: 2-Way Sec Tag Banking

### Core Insight

In ahead TAGE, the index into each table is always stale — computed
from history that doesn't include block A's outcome. This is inherent
to ahead pipelining and cannot be fixed by recomputing the index.

However, the sec tag already captures which history path (A's actual
exit) the entry was trained for. The problem is that each table entry
currently stores only ONE sec tag / prediction pair. If the entry was
trained under a different exit of A, the sec tag won't match and we
fall back to bimodal.

**Fix**: Store 2 ways per entry — each with its own sec tag, prediction
counter, and u-bit. Both ways are read at the same stale index in
predict1. In update_cycle, the actual sec tag selects the matching
way. This covers the 2 most common exit paths of the preceding block.

### What Gets 2-Way Banking

| Component | Banked? | Reason |
|-----------|---------|--------|
| pred (ctr) | **Yes** | Different predictions per history path |
| u-bit | **Yes** | Independent usefulness per way |
| sec_tag | **Yes** | Each way identifies its history path |
| tag | No | Same primary tag (same stale index) |
| hyst | No | Shared — trained post-resolution, way already selected |

### No BTB Needed

Unlike the previous design (see "Abandoned: Successor BTB" below),
this approach does NOT require a successor BTB or a second fold
computation. Both ways share the same stale index — the sec tag alone
disambiguates. This eliminates:
- Successor BTB RAM
- Second fold computation in predict1
- Second set of table RAMs
- All bank 1 prefetch/current register pipelines
- GH_FANOUT doubling

### Pipeline

```
predict1(A):
  t=25ps   inst_pc ready
  t=31ps   fold_idx/fold_tag ready (from regs)
  t=~80ps  RAM reads complete (same index for both ways)
           → prefetch_pred[way0], prefetch_pred[way1]
           → prefetch_sec[way0], prefetch_sec[way1]
           → prefetch_u[way0], prefetch_u[way1]
  t=300ps  cycle end — all prefetched into regs

update_cycle(A):
  t=0ps    current_pred[way0/1], current_sec[way0/1], current_u[way0/1] ready
  t=~10ps  sec tag comparison: actual_sec vs way0_sec, way1_sec
  t=~15ps  mux: select matching way → current_pred, current_u
  t=~67ps  full_hits (tag match, +52ps from mux)
  ...existing provider resolution chain...
  t=~249ps pred scatter (+234ps chain + 15ps way mux)

  Critical path: ~249ps (was 234ps, +15ps for way mux)
  Within 300ps cycle ✓
```

### Sec Tag Matching Logic

In predict1, both ways' sec tags are read into prefetch registers
(same index, wider read). In update_cycle:

```
actual_sec = hash(actual_next_pc)    // from block_end_info.next_pc
way0_match = (current_sec[way0] == actual_sec)
way1_match = (current_sec[way1] == actual_sec)

if way0_match → use way 0 pred/u
elif way1_match → use way 1 pred/u
else → bimodal fallback (neither way trained for this path)
```

When b1_match (either way matches), the per-table sec tag check can be
skipped entirely — the history is confirmed correct, no ahead aliasing.
Controlled by `BANK_SEC_FALLBACK` param for the fallback (neither-match)
case.

## Training

### Way Selection on Training

Training writes go to the matching way:
- If way 0's sec tag matches actual → train way 0
- If way 1's sec tag matches actual → train way 1
- If neither matches → allocate into a way (replace the less useful
  one, i.e. lower u-bit, or way 0 by default)

### Allocation

On TAGE allocation (misprediction), the new entry is written into the
way that matches the current sec tag. If neither way matches, replace
the way with lower u (or way 0 on tie). The other way is untouched.

### Hyst Training

Hyst is shared across ways (not banked). After the way mux selects
which prediction to use, hyst training proceeds as before using the
selected way's prediction. The hyst index uses the same stale index.

## Storage Cost

| Component        | Current          | With 2-Way Banking          |
|------------------|------------------|-----------------------------|
| TAGE pred RAM    | NT RAMs          | NT RAMs (2× width per entry)|
| TAGE tag RAM     | NT RAMs          | Unchanged                   |
| TAGE u RAM       | NT × U_BANKS     | NT × U_BANKS (2× width)     |
| TAGE sec tag RAM | NT RAMs          | NT RAMs (2× width per entry)|
| TAGE hyst RAM    | NT × HYST_BANKS  | Unchanged                   |
| Successor BTB    | —                | Not needed                  |

Storage overhead: ~2× for pred, u, sec_tag per entry. Tag and hyst
unchanged. Total overhead much less than full table replication.

## Graceful Degradation

- Way 0 miss + way 1 miss → bimodal fallback (unchanged from current)
- One way populated, other empty → behaves like current 1-way design
- Both ways populated for same sec tag → degenerate but harmless
  (both match, either gives correct prediction)
- Can only improve over current TageAhead, never regress

## Interaction with Existing Features

### Allocation Policy
Orthogonal to MAX_ALLOC / AllocCfg / target policy. The allocator
picks which TABLE. Way selection is within the table entry.

### Per-branch Banking (HYST_BANKS, U_BANKS)
Orthogonal. HYST_BANKS banks by branch offset within the block.
2-way sec tag banking banks by history path. Both dimensions coexist.
For u_ram: each U_BANK gets 2 ways.

### Epoch / U-bit Reset
Epoch resets clear both ways' u-bits.

### FARALLOC_DIST
Unchanged — measures table distance, independent of ways.

## Relationship to Ahead-Prediction Paper (Cai et al., ISCA '25)

The paper proposes 2^M predictions with MUX at consumption. This is a
practical 2-way version: 2 predictions per entry covering the 2 most
common exit paths of the preceding block. With avg 1.6 branches/block,
2 ways covers the dominant cases without needing a BTB or second fold.

## Implementation Plan

### Template Params
- `bool ENABLE_BANKING = false` — master switch
- `u64 BANK_SEC_TAG_BITS = 3` — per-way sec tag width
- `bool BANK_SEC_FALLBACK = true` — use sec_ram check when neither way matches

### Changes to pred_ram, u_ram, sec_ram
Widen to store 2 ways. Options:
1. Double the data width (pack 2 ways into one entry)
2. Use the existing banking dimension (add another bank level)
3. Separate RAM arrays for way 0 and way 1

### predict1
Read both ways from same index (already done — just widen the read).
Store into prefetch_pred[way0/1], prefetch_sec[way0/1], etc.

### update_cycle
Before provider resolution: sec tag comparison → way mux.
Rest of resolution chain unchanged.

### Training
Route pred/u writes to matching way. Allocation replaces lower-u way.

## Open TODOs

1. Implement 2-way banking in TageAhead.hpp
2. Decide storage approach (wider entries vs separate arrays)
3. Sweep BANK_SEC_TAG_BITS (1, 3, 6)
4. Measure way 0 vs way 1 hit rates

---

## Abandoned: Successor BTB + Dual-Fold Design

The original design used a successor BTB to predict block A's next_pc,
then computed a second fold in predict1 to index a full second set of
table RAMs. This was abandoned because:

1. **Unnecessary complexity**: The stale index is always the same for a
   given block. Different history resolutions don't change the index —
   they change which entry at that index is correct. The sec tag already
   captures this distinction.
2. **High cost**: Required full table replication (2× tag, pred, u RAMs),
   a BTB RAM, second fold computation, doubled GH_FANOUT, and a full
   second set of prefetch/current pipeline registers.
3. **BTB overhead**: BTB training policy, sizing, pipelining decisions
   all added complexity with uncertain payoff.

The 2-way sec tag approach achieves the same disambiguation with just
wider entries — no BTB, no second fold, no table replication.

### Timing Probe Results (preserved for reference)

Measured with `/tmp/btb_timing_test.cpp` on gcc trace, 4K instructions.

| BTB Size | Direct (slack) | Prefetched (slack) |
|----------|----------------|--------------------|
| 64       | 271ps (29ps)   | 215ps (85ps)       |
| 128      | 279ps (21ps)   | 215ps (85ps)       |
| 256      | 288ps (12ps)   | 216ps (84ps)       |
