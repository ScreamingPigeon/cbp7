# TageTable Class Design Plan

## Summary
One class instance per TAGE table. The predictor holds N instances and handles all cross-table logic (match priority, allocation selection, u-bit clearing).

## Decided

### Ownership
- **Predictor** owns fold registers, index/tag computation, ahead pipeline buffering, path/XB management, cross-table logic
- **Table** is a banked RAM wrapper with register caching — no fold, no pipeline, no timing awareness
- Predictor passes pre-computed index, tag, and bank index to the table

### Parameters (confirmed)
| Parameter | Type | Description |
|-----------|------|-------------|
| `TABLE_SIZE` | u64 | Number of entries per bank |
| `TABLE_HIST` | u64 | History length (used by predictor for fold, not by table directly) |
| `TAG_WIDTH` | u64 | Tag width in bits |
| `CTR_WIDTH` | u64 | Prediction counter width |
| `U_WIDTH` | u64 | Useful counter width |
| `N` | u64 | Max branches per cycle = total predictions per block (replaces PRED_BLK_SIZE) |
| `NUM_BANKS` | u64 | Banks partitioning branch slots. BPB = N/NUM_BANKS counters per bank. |
| `USE_AHEAD` | bool | 1-ahead pipelining. Doubles physical banks (2 × NUM_BANKS). |
| `SHARED_TAG` | bool | Share one tag across all branch-slot banks (NUM_BANKS dimension). See design notes. |
| `SHARED_U` | bool | Share one u-bit across all branch-slot banks. Constraint: `!SHARED_U` requires `!SHARED_TAG`. See design notes. |
| `U_STORAGE` | enum | `SRAM` or `FF`. Mutually exclusive u-bit management strategies. |
| `DECAY_CTR` | u64 | Only when `U_STORAGE = SRAM`. Probabilistic per-access decay probability (1/DECAY_CTR). |
| `ResetFn` | typename | Only when `U_STORAGE = FF`. Functor `(val<U_WIDTH>, val<MODE_BITS>) -> val<U_WIDTH>`. Default provides reset/rshift/decrement. Predictor selects mode at runtime. |
| `USE_FF_CACHE` | bool | Cache SRAM reads in FFs for reuse across block. Constraint: `!USE_FF_CACHE \|\| BPB > 1`. See design notes. |

### Interface (confirmed)
- Read: predictor passes index → table reads all banks, caches in registers, returns results
- Write: predictor passes index + bank index + data → table writes to specific bank
- Reuse: predictor reads cached registers (no RAM access)
- U reset (FF mode): predictor calls `reset_u(mode)` → table applies ResetFn to all FF entries

### Result access (confirmed)
- Read takes index + comparison tag. Table compares tag against all banks internally.
- Accessors: `getHit(bank)`, `getTag(bank)`, `getCounter(bank, slot)`, `getU(bank)`
- Predictor controls fanout per-field
- For ahead: predictor buffers the tag and passes the correct (previous cycle's) tag at read time

## Design Notes

### SHARED_TAG (Q12)
`SHARED_TAG` only applies to the branch-slot dimension (`NUM_BANKS`). The ahead duplication (`USE_AHEAD`) always has independent state — the two ahead copies are different pipeline stages with different indexes entirely.

- **`SHARED_TAG = true`** (default): One tag covers all `NUM_BANKS` branch-slot banks. All banks hit or miss together. Less storage (`1 × TAG_WIDTH` vs `NUM_BANKS × TAG_WIDTH`). Natural fit when banks partition branch positions within a single fetch block — same context, one tag suffices. Allocation is coarser (replace all counters at once).
- **`SHARED_TAG = false`**: Each bank has its own tag. Banks can track independent contexts at the same index. Partial hits possible. More storage, more flexible allocation. Needed when banks might serve branches from different contexts that collide on the same index.

Expected usage: `SHARED_TAG = true` for most configurations, unless ahead prediction requires per-bank independence.

### SHARED_U (Q13)
`SHARED_U` follows the same logic as `SHARED_TAG` and must be enforced via `static_assert(!SHARED_U || SHARED_TAG)` — per-branch U is pointless with a shared tag since eviction is all-or-nothing anyway.

- **`SHARED_U = true`** (default): One u-bit per entry (or per bank if `SHARED_TAG = false`). Coarser — one useful branch protects the entire entry.
- **`SHARED_U = false`**: One u-bit per bank. Finer eviction info, but only actionable with per-bank tags.

Interaction with U storage modes:
- **FF mode**: `SHARED_U = true` → `TABLE_SIZE` FFs. `SHARED_U = false` → `TABLE_SIZE × NUM_BANKS` FFs. `ResetFn` applied to each.
- **SRAM mode**: `SHARED_U = true` → one u-field per entry. `SHARED_U = false` → one u-field per bank entry. Decay on tag miss applies per-bank when `SHARED_TAG = false`.

### USE_FF_CACHE (Q15)
FF caching of SRAM reads for reuse across a fetch block. When enabled, the first branch in a block reads SRAM and caches the entire entry (tag + counters + u) in FFs. Subsequent branches read from cached FFs — saving `BPB - 1` SRAM reads per block.

- **`USE_FF_CACHE = true`**: Requires `BPB > 1` (`static_assert`). FF cost: `NUM_BANKS × (TAG_WIDTH + BPB × CTR_WIDTH + U_WIDTH)` bits. Enables `reuseRead()` pattern.
- **`USE_FF_CACHE = false`**: Every branch slot requires its own SRAM read. Only valid option when `BPB = 1` (nothing to reuse). Can also be used with `BPB > 1` if FF area is a concern (trades energy for area).

When `BPB = 1`, caching is pointless — one read gives one prediction with no reuse opportunity.

### Allocation (Q11)
- No separate allocation method. Same write interface for both updates and allocation.
- Predictor owns allocation policy (which table, when, steal logic, initial counter values).
- Predictor casts initial values to table's `CTR_WIDTH`/`U_WIDTH` and calls the normal write.
- Table just packs and stores what it's given.

### Confidence & Hysteresis (Q16, Q17)
- No separate confidence or hysteresis bits in the table.
- `CTR_WIDTH > 1` embeds both: MSB = prediction, distance from threshold = confidence/hysteresis.
- Reference tage's 1-bit pred + 2-bit hyst is just a different encoding of a 3-bit counter.
