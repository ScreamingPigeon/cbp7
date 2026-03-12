# BPB>1 Support and FF_CACHE Integration

## Context

The custom Tage predictor currently only works properly with BR_P_ENTRY=1 (1 branch per entry, offset in tag). When BR_P_ENTRY>1, multiple branches share one tag and each has its own prediction slot in a packed pred_ram entry (BPB = FETCH_WIDTH / NUM_BANKS slots per bank). This packed storage is broken:

1. `readc[I]` is `reg<MAX_CTR_WIDTH>` (1-bit) — truncates the full `PRED_ENTRY_BITS`-wide packed entry
2. Lines 699, 727 have `// simplified for multi-slot` stubs
3. `pred_ram` writes overwrite the entire packed entry instead of a single slot
4. `hyst_ram` writes same issue

FF_CACHE (TageTable feature: `static_assert(!USE_FF_CACHE || BPB > 1)`) is designed for BPB>1: one RAM read returns BPB predictions cached in FFs. Subsequent branches reuse from FFs via `reuseRead()` instead of re-reading RAM. Currently unusable because BPB>1 slot extraction is broken.

**Non-regression constraint**: BR_P_ENTRY=1 (default config) must produce identical EPI, latency, and mispredictions. All changes to BR_P_ENTRY=1 code paths must be zero-cost.

## Current BR_P_ENTRY>1 Code Paths

The code already has `if constexpr (BR_P_ENTRY_V == 1) / else` branches in:
- **Tag comparison** (lines 468-493): BR_P_ENTRY=1 encodes offset in tag; BR_P_ENTRY>1 uses full tag match for all offsets — **working**
- **Tag write** (lines 809-815): BR_P_ENTRY=1 writes `concat(offset, htag)`; BR_P_ENTRY>1 writes full tag — **working**
- **bdir** (lines 688-701): BR_P_ENTRY>1 stub: `return branch_taken[0]` — **broken**
- **altdiffer** (lines 723-728): BR_P_ENTRY>1 stub: `return pred_dir != pred2[0]` — **broken**
- **goodpred** (lines 732-741): BR_P_ENTRY>1 stub: just `return correct_pred` — **working** (simplified but correct)

## Changes Needed

### Files to Modify

1. **`predictors/custom/Tage.hpp`** — Main changes in TageImpl direct specialization
2. **`tests/test_tage_compile.cpp`** — Add BR_P_ENTRY>1 test instantiations

### Step 1: Add Constants

Add to TageBase or TageImpl:
```cpp
static constexpr u64 BPB = FETCH_WIDTH / NUM_BANKS_V;
static constexpr u64 LOG_BPB = clog2(BPB);
// Per-bank packed entry widths (using max across tables)
static constexpr u64 MAX_PRED_ENTRY_BITS = BPB * MAX_CTR_WIDTH;
static constexpr u64 MAX_HYST_ENTRY_BITS = SHARED_HYS_V ? MAX_HYST_WIDTH : BPB * MAX_HYST_WIDTH;
```

### Step 2: Widen Result Registers (BR_P_ENTRY>1 only)

Use conditional widths so BR_P_ENTRY=1 doesn't change:
```cpp
static constexpr u64 READC_WIDTH =
    (BR_P_ENTRY_V == 1) ? MAX_CTR_WIDTH : MAX_PRED_ENTRY_BITS;
static constexpr u64 READH_WIDTH =
    (BR_P_ENTRY_V == 1) ? std::max(u64(1), MAX_HYST_WIDTH)
                         : std::max(u64(1), MAX_HYST_ENTRY_BITS);

arr<reg<READC_WIDTH>, NUM_TABLES> readc;
arr<reg<READH_WIDTH>, NUM_TABLES> readh;
```

BR_P_ENTRY=1: readc stays `reg<1>`, readh stays `reg<2>`. Zero EPI change.

### Step 3: Per-Slot Prediction Extraction in predict2

Currently, `gpreds` is the same NUM_TABLES-bit mask for all offsets. With BR_P_ENTRY>1, each offset needs its own prediction extracted from its slot in the packed entry.

```cpp
if constexpr (BR_P_ENTRY_V == 1) {
    // Current code unchanged: one prediction bit per table
    val<NUM_TABLES> gpreds = readc.concat();  // (for CTR_WIDTH=1)
    gpreds.fanout(hard<FETCH_WIDTH>{});
    arr<val<NUM_TABLES + 1>, FETCH_WIDTH> preds = [&](u64 offset) {
        return concat(readb[offset], gpreds);
    };
} else {
    // Per-offset extraction from packed entry
    arr<val<NUM_TABLES + 1>, FETCH_WIDTH> preds = [&](u64 offset) {
        arr<val<1>, NUM_TABLES> slot_pred = [&](int i) -> val<1> {
            // Extract prediction for slot 'offset' from packed entry
            if constexpr (MAX_CTR_WIDTH == 1) {
                return readc[i] >> offset;  // bit 'offset'
            } else {
                return readc[i] >> (offset * MAX_CTR_WIDTH + MAX_CTR_WIDTH - 1);  // MSB
            }
        };
        return concat(readb[offset], slot_pred.fo1().concat());
    };
}
```

Note: `offset` here is a compile-time `u64` from the lambda, not a `val<>`. With BR_P_ENTRY>1, all offsets share the same tag match, so `match[offset]` is the same for all. The per-offset differentiation is in the prediction bits.

### Step 4: Fix bdir (line 699)

For BR_P_ENTRY>1: `bdir[i]` needs to find the actual direction of the branch in this table's slot.

```cpp
if constexpr (BR_P_ENTRY_V == 1) {
    // ... existing code (offset from tag) ...
} else {
    // Each table entry covers all BPB slots.
    // The "direction" for this table is a vector of all branches' directions.
    // For update, we need the direction of branches that are in this entry.
    // Simplification: use fold of branch_taken weighted by is_branch
    // Or: associate each table with the branch directions for all slots
    return branch_taken[0]; // TODO: proper multi-slot direction
}
```

Actually, for BR_P_ENTRY>1 with full tag match: one entry covers all branches. Multiple branches can be correct or incorrect. The update logic needs to determine:
- Which SLOT to update in pred_ram/hyst_ram (the mispredicted branch's slot)
- The direction for that specific slot

The correct implementation:
```cpp
} else {
    // Multi-slot: for provider table's slot corresponding to this table's matched branch,
    // use that branch's actual direction
    // With one tag per entry, all BPB slots share the same gindex.
    // For the mispredicted branch at last_offset: that's the slot to update.
    // For other tables (non-primary): their bdir is for whichever branch they were provider for.
    // Simplified: use last branch's direction for update
    return branch_taken[num_branch - 1]; // last branch in block
}
```

### Step 5: Fix altdiffer (line 727)

```cpp
} else {
    // Multi-slot: compare this table's prediction for each slot against
    // the alt prediction for the same slot
    // For the provider table and the branch it provides for:
    // pred_dir is the prediction from readc slot, pred2 is alt prediction for same offset
    // Simplified: use slot 0 for now, or the last branch's slot
    return pred_dir != pred2[0]; // TODO: proper slot indexing
}
```

### Step 6: Slot-Aware Writes to pred_ram

Currently (line 874): `pred_ram[0].write(gindex[I], bdir[I])` — writes 1-bit to full entry.

For BR_P_ENTRY>1, need read-modify-write:
```cpp
if constexpr (BR_P_ENTRY_V == 1) {
    // Current: write full entry (1-bit truncation is fine)
    std::get<I>(tables).pred_ram[0].write(gindex[I], bdir[I]);
} else {
    // Read-modify-write: update only the slot for last_offset
    // readc[I] has the current packed entry (from predict2 cache)
    val<MAX_PRED_ENTRY_BITS> new_pred = readc[I];
    // Clear and set the slot for last_offset
    // ... bit manipulation to update slot at position (last_offset * CTR_WIDTH) ...
    std::get<I>(tables).pred_ram[0].write(gindex[I], new_pred);
}
```

For allocation: initialize all slots to 0 except the allocating branch's slot:
```cpp
if (allocate[I]) {
    // Set only the allocated branch's slot
    val<MAX_PRED_ENTRY_BITS> init_pred = bdir[I] << (last_offset * CTR_WIDTH);
    pred_ram[0].write(gindex[I], init_pred);
}
```

### Step 7: Slot-Aware Writes to hyst_ram

Same pattern as pred_ram. For SHARED_HYS=true: one hysteresis per entry (shared across slots) — no change needed. For SHARED_HYS=false: packed per-slot hysteresis, need slot-aware R-M-W.

### Step 8: FF_CACHE Integration

When USE_FF_CACHE=true and BPB>1:
- **predict2** (first branch): Read RAM → results go to TageTable staging regs (automatically via `ram::read()`)
- **reuse_predict2** (subsequent branches): Call `reuseRead()` on each table → get cached packed entry → extract new slot

Currently reuse_predict2 just returns from precomputed p2. With FF_CACHE:
```cpp
val<1> reuse_predict2(val<64> inst_pc) override {
    if constexpr (USE_FF_CACHE_V && BR_P_ENTRY_V > 1) {
        // Re-extract predictions for new branch slot from cached entries
        static_loop<NUM_TABLES>([&]<u64 I>() {
            auto &t = std::get<I>(tables);
            t.reuseRead(0, readt[I], readc[I], readh[I], readu[I]);
        });
        // Recompute tag matching + provider selection for new offset
        // (tag match is same — cached from predict2)
        // Extract slot for new branch_offset → update pred1/pred2/p2
    } else {
        // Current: return from precomputed p2
        val<1> taken = ((block_entry << block_size) & p2) != hard<0>{};
        reuse_prediction(~val<1>{block_entry >> (FETCH_WIDTH - 1 - block_size)});
        block_size++;
        return taken;
    }
}
```

Note: For BR_P_ENTRY=1 (even with FF_CACHE), reuse doesn't help because each branch has a different tag → needs separate RAM reads. FF_CACHE only saves energy for BR_P_ENTRY>1.

## Verification

1. **Non-regression (BR_P_ENTRY=1)**: Run default `Tage<>` — must match current metrics exactly:
   - EPI = 1903, P2 latency = 1.86, mispredictions = 26,164
   - No storage increase (readc/readh widths unchanged for BR_P_ENTRY=1)
2. **BR_P_ENTRY>1 compile**: Instantiate with BR_P_ENTRY=2 or FETCH_WIDTH, verify clean compile
3. **BR_P_ENTRY>1 functional**: Run benchmark with BR_P_ENTRY>1, verify no crashes
4. **FF_CACHE + BPB>1**: Run with USE_FF_CACHE=true, BR_P_ENTRY>1, compare EPI vs without FF_CACHE
5. **Compile tests**: All existing T1-T9 pass, add new tests for BR_P_ENTRY>1 configs

## Implementation Order

1. Add constants (Step 1) — zero risk
2. Conditional register widths (Step 2) — verify non-regression immediately
3. Per-slot prediction extraction (Step 3) — core logic change
4. Fix bdir/altdiffer stubs (Steps 4-5) — correctness for BR_P_ENTRY>1
5. Slot-aware writes (Steps 6-7) — write path correctness
6. FF_CACHE integration (Step 8) — energy optimization
7. Test throughout
