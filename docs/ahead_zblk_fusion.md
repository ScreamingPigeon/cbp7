# Ahead TAGE: Zero-Block Fusion Datapath

## Motivation

17.3% of blocks on gcc have 0 conditional branches (zero-blocks). These are
unconditional jumps, calls, returns — blocks that consume a full predict1 cycle
of TAGE reads but produce no useful predictions. 98% of them exit via
non-sequential control flow (jump targets, not fall-through).

Zero-block fusion steals one path slot from the current block's ahead reads to
pre-fetch the zero-block's successor. The zero-block itself needs no TAGE
prediction (no conditional branches), so the stolen slot is used productively.

This only works in an **ahead TAGE** design, where predict1 reads are for a
future block and predict2 is purely combinational on regs. In direct TAGE,
serializing zblk_read → TAGE_index → TAGE_read → tag_compare in one cycle
would exceed P2=2.

## Empirical Data (gcc trace, 40M instructions, LINEINST=1024, N=7)

```
Zero-blocks:           575,278 (17.3% of all blocks)
Sequential successor:   11,276 ( 2.0% of zero-blocks)
Non-sequential:        564,002 (98.0% of zero-blocks)
Successor has branches: 466,499 (81.1% — fusable, TAGE read is productive)
Successor also zero:   108,778 (18.9% — chain, need recursive fusion or skip)
```

## Hardware Structures

### TAGE Tables

```
ram<tage_entry, TABLE_SIZE>  tables[NUM_TABLES][NUM_PATHS]
```

- `NUM_TABLES` tagged tables, each with `NUM_PATHS` banks (one per exit path
  of the predecessor block).
- Entry: `[tag | ctr_0..ctr_{N-1} | hyst | u]`
- Total: `NUM_TABLES × NUM_PATHS` physical RAMs.
- In predict1, ALL banks of all tables are read in parallel. In predict2,
  the secondary tag selects which path bank's data is valid.

### ZBLK Table

```
ram<zblk_entry, ZTABLE_SIZE>  zblk_table
```

- Small BTB-like table for zero-block detection.
- Entry format:

```
┌────────────┬────────────────┬───────────────────┐
│ tag (T bits)│ size (S bits) │ next_pc (P bits)  │
└────────────┴────────────────┴───────────────────┘
```

- **tag**: Partial block_pc for disambiguation. T = 8-10 bits.
- **size**: Number of instructions in the zero-block. S = 4-6 bits
  (max 64, covers 99%+ of zero-blocks).
- **next_pc**: Successor block's PC. P = ~20 bits (aligned, low bits dropped).
- **Total entry width**: ~32-36 bits.
- **Table size**: 256-1024 entries (1-4.5 KB). <5% of TAGE storage.

### Pipeline Registers

```
// Per-table, per-path: tag, counter, hysteresis, u-bit read from TAGE
reg<...> ahead[2][NUM_TABLES][NUM_PATHS]   // [0]=current reads, [1]=previous

// Secondary tags for path selection
reg<SEC_TAG_BITS> sec_tag[2][NUM_TABLES][NUM_PATHS]

// ZBLK state
reg<1> zblk_active          // current block is a fused zero-block
u64    zblk_remaining = 0   // C++ counter: instructions left in zero portion
```

## Datapath: predict1 (Two Parallel Paths)

```
predict1(block_pc):
  ╔══════════════════════════════════════════════════════════════════╗
  ║  PATH A — Serve current block's prediction (from pipeline regs) ║
  ╚══════════════════════════════════════════════════════════════════╝
  │
  ├─ Compute secondary_tag = hash(block_pc)
  ├─ Compare against NUM_PATHS secondary tags in ahead[1]
  │    (ahead[1] was written by PREVIOUS block's predict1)
  ├─ Secondary tag hit → select matching path's data → provider select
  ├─ Secondary tag miss → bimodal fallback
  └─ Return: P1 prediction
  
  ╔══════════════════════════════════════════════════════════════════╗
  ║  PATH B — Ahead reads for NEXT block (stored into regs)         ║
  ╚══════════════════════════════════════════════════════════════════╝
  │
  ├─ [PARALLEL] Pipeline shift: ahead[0] → ahead[1]
  ├─ [PARALLEL] Read ZBLK table at hash(block_pc)
  │
  ├─ ZBLK tag compare (combinational)
  │
  ├─ If ZBLK HIT:
  │     zblk_active = 1
  │     zblk_remaining = zblk.size
  │     // Override path 0 (fall-through) with zblk successor read
  │     // Zero-blocks always exit via fall-through (no cond branches taken),
  │     // so path 0 is the exit path — overriding it is exact, not a stolen slot
  │     tage_read_pc[0] = hash(zblk.next_pc, folded_hist)
  │     tage_read_pc[1..NUM_PATHS-1] = hash(block_pc, folded_hist)
  │     // path 0: ahead read for zblk successor (replaces fall-through)
  │     // path 1..NUM_PATHS-1: taken-exit successor coverage (unchanged)
  │
  ├─ If ZBLK MISS:
  │     zblk_active = 0
  │     tage_read_pc[0..NUM_PATHS-1] = hash(block_pc, folded_hist)
  │     // all NUM_PATHS slots for current block's successors (normal)
  │
  ├─ Read tables[t][p] at tage_read_pc[p] for all t, p
  └─ Store results → ahead[0]
```

### Path Slot Allocation

```
                    ZBLK MISS                         ZBLK HIT
            ┌───────────────────┐            ┌───────────────────┐
  path 0    │ fall-through succ │   path 0   │ ZBLK next_pc read │ ← overridden
  path 1    │ taken-1 successor │   path 1   │ taken-1 successor │
  ...       │ ...               │   ...      │ ...               │
  path P-2  │ taken P-2 succ   │   path P-2 │ taken P-2 succ    │
  path P-1  │ taken P-1 succ   │   path P-1 │ taken P-1 succ    │
            └───────────────────┘            └───────────────────┘

On ZBLK hit: path 0 (fall-through) is overridden with the zblk successor read.
No coverage loss for taken paths. The zero-block's exit IS the fall-through
path, so the override is semantically exact — not a stolen slot.
Penalty for false positive: fall-through coverage lost for 1 cycle.
```

## Datapath: predict2 / reuse_predict2

```
predict2(block_pc):
  ├─ If zblk_active && zblk_remaining > 0:
  │     zblk_remaining--
  │     reuse_prediction(1)           // extend block through zero portion
  │     return NOT_TAKEN              // no cond branches, ignored by sim
  │
  ├─ If zblk_active && zblk_remaining == 0:
  │     // Now serving successor block predictions
  │     // ahead[1] contains reads from previous predict1
  │     // Secondary tag match against branch_pc selects path
  │     // Provider selection across NUM_TABLES → prediction
  │     (same logic as normal ahead TAGE predict2)
  │
  └─ If !zblk_active:
        // Normal ahead TAGE predict2
        Secondary tag match on ahead[1] → provider select → prediction
```

**No additional P2 latency**: zblk_active is a reg. The mux is a single
select on a 1-bit signal — negligible.

## RAM Structure and Banking

### TAGE Table RAMs

Each TAGE table-path combination is a banked rwram:

```
rwram<{tag, ctr[LANES]}, TABLE_SIZE, RWRAM_BANKS>  tables[NUM_TABLES][NUM_PATHS]
```

- `RWRAM_BANKS`: parameterized (2, 4, 8, ...). Higher = fewer conflicts, more area.
- predict1 reads at `index_current`, update_cycle writes at `index_previous`.
- Conflict iff `index_current % RWRAM_BANKS == index_previous % RWRAM_BANKS`.
- Conflict probability ≈ `1 / RWRAM_BANKS` (indices are ~uniformly distributed).

| RWRAM_BANKS | Conflict rate | Extra cycle ratio (approx) |
|-------------|--------------|---------------------------|
| 2           | ~50%         | ~50%                      |
| 4           | ~25%         | ~25%                      |
| 8           | ~12.5%       | ~12.5%                    |
| 16          | ~6.25%       | ~6.25%                    |

### Separate Update-Only RAMs

U-bits and hysteresis are **never read during prediction** — only during
allocation (mispredict path, inside extra_cycle). So they live in separate
UPDATE_ONLY RAMs with no banking needed:

```
ram<u_entry, TABLE_SIZE>     u_ram[NUM_TABLES]       // UPDATE_ONLY zone
ram<hyst_entry, HYST_SIZE>   hyst_ram[NUM_TABLES]    // UPDATE_ONLY zone
```

These are written every cycle (u-bit set/decay, hyst update) with zero
conflict — predict1 never touches them.

### Shared Hysteresis (parameterized)

When `SHARED_HYS = true`, one hysteresis entry covers two consecutive TAGE
entries. This halves the hysteresis RAM size:

```
HYST_SIZE = SHARED_HYS ? TABLE_SIZE / 2 : TABLE_SIZE
HYST_IDX_BITS = IDX_BITS - (SHARED_HYS ? 1 : 0)
hyst_index = stored_index >> (SHARED_HYS ? 1 : 0)
```

- Two entries at index `2k` and `2k+1` share `hyst_ram[k]`.
- Saves storage at the cost of slight interference between neighbors.
- Same pattern as current TageDirect `SHARED_HYS` / `HYST_IDX_BITS`.
- `SHARED_HYS` is a template parameter — no runtime cost when false.

### ZBLK Table RAM

```
ram<zblk_entry, ZTABLE_SIZE>  zblk_table
```

Read in predict1, written in update_cycle. Since predict1 reads at
`hash(block_pc_current)` and update_cycle writes at `hash(block_pc_previous)`,
these are different indices. Can use rwram with banking, or accept the
occasional extra_cycle for zblk writes (rare — only when a new zero-block
is detected).

## Datapath: update_cycle

### Cycle 1 (same cycle as predict1/predict2)

Combinational logic + reads from UPDATE_ONLY RAMs. No TAGE table writes yet.

```
update_cycle(block_end_info):
  │
  ├─ Compute update signals (combinational on regs):
  │     provider_t, provider_p    // which table/path was the provider
  │     stored_index              // index used for the provider's read (from pipeline reg)
  │     new_ctr = update_ctr(old_ctr, correct)  // saturating increment/decrement
  │
  ├─ Silent write elimination:
  │     ctr_changed = (new_ctr != old_ctr)
  │     // old_ctr is already in pipeline regs (read during predict1 of prev block)
  │     // If counter is already saturated in the correct direction, skip the write
  │     // This avoids an rwram write (and potential bank conflict) when unnecessary
  │
  ├─ Bank conflict detection:
  │     write_bank = stored_index % RWRAM_BANKS
  │     read_bank  = index_current % RWRAM_BANKS   // from this cycle's predict1
  │     bank_conflict = (write_bank == read_bank) & ctr_changed
  │
  ├─ need_extra_cycle(mispredict | bank_conflict)
  │
  ├─ U-bit / hysteresis writes (UPDATE_ONLY, no conflict):
  │     u_ram[provider_t].write(stored_index, new_u)
  │     hyst_ram[provider_t].write(stored_index, new_hyst)
```

### Cycle 2 (extra cycle, when granted)

TAGE table writes happen here — separated from predict1 reads.

```
  ├─ execute_if(extra_cycle, [&]() {
  │
  │     // Counter write (provider)
  │     execute_if(ctr_changed, [&]() {
  │       tables[provider_t][provider_p].write(stored_index, new_entry)
  │     })
  │
  │     // Allocation (mispredict only)
  │     execute_if(mispredict, [&]() {
  │       // Find table above provider with u=0
  │       // Read u_ram to check usefulness (UPDATE_ONLY, safe)
  │       // Write new tag + counter to tables[alloc_t][alloc_p]
  │       // Write u_ram, hyst_ram for new entry
  │     })
  │
  │     // ZBLK table training
  │     execute_if(is_zero_block & ~zblk_active, [&]() {
  │       zblk_table.write(hash(block_pc), {tag, size, next_pc})
  │     })
  │
  │     // ZBLK entry invalidation
  │     execute_if(zblk_active & mispredict, [&]() {
  │       // Evict or update stale zblk entry
  │     })
  │  })
  │
  ├─ true_block computation:
  │     true_block = ~mispredict | last_condbr_dir | last_pred() | line_end()
  │
  └─ History update (gated by true_block):
        execute_if(true_block, [&]() { gfolds.update(...); })
```

## Silent Write Elimination

Counter writes are the primary source of bank conflicts. Many updates are
**silent** — the counter is already saturated and the update wouldn't change it:

- Correct prediction, counter already at max → saturating increment is a no-op
- Misprediction, counter already at min → saturating decrement is a no-op

By comparing `new_ctr` against `old_ctr` (both available in regs — old from
the ahead pipeline read, new from combinational update logic), we gate the
write:

```
ctr_changed = (new_ctr != old_ctr)
bank_conflict = (write_bank == read_bank) & ctr_changed
need_extra_cycle(mispredict | bank_conflict)
execute_if(ctr_changed & extra_cycle, [&]() {
    tables[provider_t][provider_p].write(stored_index, new_entry)
})
```

### Impact estimate

With 3-bit signed counters (range -4..3), a correctly-predicted branch with
a saturated counter (+3) generates a silent write — no actual update needed.
Well-trained entries (the majority) are saturated. Rough estimate:

- ~60-70% of correct predictions have saturated counters → silent
- Effective write rate: ~30-40% of blocks instead of ~100%
- Combined with RWRAM_BANKS=4: conflict rate ≈ 0.35 × 25% ≈ 9%
- With RWRAM_BANKS=8: ≈ 0.35 × 12.5% ≈ 4%

Silent write elimination + moderate banking can bring extra_cycle ratio
below 10% while preserving full TAGE update semantics.

## true_block Interaction

```
true_block = ~mispredict | last_condbr_dir | last_pred() | line_end()
```

- `zblk_active` must be cleared on false blocks:
  `zblk_active = select(true_block, zblk_hit, 0)`
- A false block cannot use zblk fusion — the predecessor hasn't changed,
  pipeline regs are stale, reusing them is the correct (degraded) behavior.
- Zero-block portion of a fused block never mispredicts (no cond branches).
  Mispredictions only occur in the successor portion → standard true_block.

## ZBLK Miss / Stale Entry Handling

| Scenario | Detection | Action |
|----------|-----------|--------|
| Alias (wrong entry) | reuse_predict1(pc) outside expected range | End block via reuse_prediction(0), flush entry |
| Size too large | Conditional branch before zblk_remaining=0 | Switch to TAGE predictions early, update size |
| Size too small | zblk_remaining=0 but no branches yet | Successor TAGE predictions used early (benign) |
| next_pc wrong | Misprediction at block end | Accuracy loss for 1 block, update entry |
| Entry missing | zblk table miss | Normal ahead path, no degradation |
| False positive | Block actually has branches | Lose 1 path slot for 1 cycle, self-corrects next visit |

**Every failure mode is self-correcting.** Worst case: one block of degraded
accuracy or one path slot wasted, then the entry is updated or evicted.

## False Positive Cost Analysis

On a ZBLK false positive (tag alias — block isn't actually a zero-block):
- Path 0 (fall-through) was used for a bogus `next_pc` read.
- Taken-path slots (1..P-1) are unaffected.
- Only fall-through successor coverage is lost for one cycle.
- If the block actually exits via fall-through, secondary tag will miss
  on path 0 → bimodal fallback for that one block.
- If the block exits via taken branch, paths 1..P-1 are intact → no loss.
- The false positive is detected (cond branch seen during "zero" portion)
  and the entry is flushed. One-time penalty.

## Chained Zero-Blocks (18.9% of zero-blocks)

When the successor of a zero-block is also a zero-block:
- The fused block would need to extend through both.
- Requires knowing the second successor's PC (not in the entry).

**Recommendation**: Don't chain. The second zero-block gets its own predict1.
18.9% × 17.3% = 3.3% of blocks — not worth the complexity.

## EPI Impact Estimate

Without fusion: zero-block + successor = 2 predict1 cycles.
- 2 × (NUM_TABLES × NUM_PATHS) TAGE reads + 2 × 1 zblk reads = 2 × 9 + 2 = 20

With fusion: fused into 1 predict1 cycle.
- 1 × (NUM_TABLES × NUM_PATHS) TAGE reads + 1 × 1 zblk read = 9 + 1 = 10
- Saving: 10 reads per fused pair.

Across all blocks: 17.3% are zero-blocks.
- Saving: 0.173 × 10 = 1.73 reads/block on average
- Cost: 1 zblk read/block (every block checks the table)
- Net: **~0.73 reads/block saved** (~8% EPI reduction)

## Storage Budget

| Component | Entries | Width | Total |
|-----------|---------|-------|-------|
| ZBLK table | 256 | 36 bits | 1.1 KB |
| ZBLK table | 512 | 36 bits | 2.3 KB |
| ZBLK table | 1024 | 36 bits | 4.5 KB |

## Open Questions

1. **Optimal ZTABLE_SIZE**: Need to measure zblk working set (unique hot
   zero-block PCs) to size the table. Too small → misses on cold entries.
   Too large → wasted storage.

2. **Can the zblk table serve double duty?** E.g., as a BTB for unconditional
   branches, or as a block-size predictor for fetch width optimization.

3. **Interaction with SC/Loop overriders**: Fused block spans two physical
   blocks. Overriders see one logical block. Their PC-based indexing may
   need adjustment for the successor portion.

4. **zblk_remaining as C++ counter**: Cannot be a HARCOM reg (reuse_predict
   is called multiple times per cycle). Must use `u64 zblk_remaining` like
   `num_branch` — free in hardware, simulator artifact.

5. **Secondary tag for the overridden path 0**: The ZBLK successor's reads go
   into path 0. The secondary tag for that slot is `hash(zblk.next_pc)` so
   that when the successor block's predict2 runs, it can find its data in
   the pipeline regs via normal secondary tag matching.

6. **History staleness**: When the fused block's successor portion runs
   predict2, the folded history hasn't been updated with the zero-block's
   path information (there's no path — no branches). Is the history used
   for indexing the stolen slot's read (at `next_pc`) stale? It shouldn't
   be — the zero-block contributes nothing to branch history, so the
   history at the successor is the same as at the zero-block entry.
