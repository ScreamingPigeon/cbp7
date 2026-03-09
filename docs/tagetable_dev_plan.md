# TageTable Development Plan

Development plan for the parameterized `TageTable` class. Builds on design decisions in `agent/tagetable_plan.md`.

## File Locations

- **Implementation**: `predictors/custom/TageTable.hpp` (rewrite existing file)
- **Test harness**: `tests/test_tagetable.cpp` (new)
- **Integration**: `predictors/tage.hpp` (modify to use TageTable)
- **Sweep scripts**: `scripts/sweep_tagetable.py` (new)
- **Design docs**: `agent/tagetable_plan.md`, `agent/tagetable_questions.md`

---

## Phase 1: Storage Primitives ✅ COMPLETE

**Goal**: Declare all HARCOM private members, parameterized by the template parameter list. No logic yet — just storage.

**Status**: All storage members declared. 14 parameter combinations compile with `-Werror`.

**Completed**:
- 1.1: Full template parameter list (13 params + ResetFn typename), all static_asserts, computed constants
- 1.2: SRAM — `bank_ram[PHYS_BANKS]`, conditional `shared_tag_ram`, conditional `shared_u_ram`. Used `hcm::ram` to avoid namespace ambiguity. `std::conditional_t<..., EmptyMember>` for conditional elimination.
- 1.3: FF — conditional `u_ff` arrays. Result/cache registers: `idx_reg`, `tag_reg`, `hit_reg`, `u_reg`, `ctr_regs` (depth = BPB when cached, 1 otherwise).
- 1.4: Compile test (`tests/test_tagetable_compile.cpp`) with 14 instantiations. Makefile target `test-tagetable` added.

**Note**: `U_STORAGE` enum replaced with `U_STOR_FF` bool per user preference.

---

## Phase 2: Accessors and Read/Write Interface — IN PROGRESS

**Goal**: Implement the read/write methods and accessors. Build a test harness.

### 2.1 Read path ✅
```
read(val<IDX_BITS> index, val<TAG_WIDTH> compare_tag, u64 stage = 0)
```
- Reads all bank RAMs at `index` for the given ahead stage
- Unpacks entries via `if constexpr` based on SHARED_TAG/SHARED_U/U_STOR_FF layout
- Extracts individual counters via `static_loop<CACHED_CTRS>` into `ctr_regs`
- Compares stored tag(s) against `compare_tag`, stores in `hit_reg`
- Reads u-bits from SRAM (packed or shared RAM) or FF arrays
- SRAM mode: applies probabilistic u-bit decay on tag miss (`1/DECAY_CTR`)

### 2.2 Reuse read ✅
```
reuseRead(u64 bank, val<clog2(BPB)> slot) -> val<CTR_WIDTH>
```
- Returns cached counter from `ctr_regs` via `arr.select()`
- `static_assert` enforces `USE_FF_CACHE=true`

### 2.3 Accessors ✅
- `getHit(bank=0)` — respects `SHARED_TAG` (shared returns `hit_reg[0]`)
- `getTag(bank=0)` — respects `SHARED_TAG`
- `getCounter(bank, slot=0)` — respects `USE_FF_CACHE` (cached returns by slot, uncached returns slot 0)
- `getU(bank=0)` — respects `SHARED_U`

### 2.4 Write path ✅
```
write(val<IDX_BITS> index, u64 bank, u64 stage, val<TAG_WIDTH> tag,
      val<CTR_BITS> ctr_bits, val<U_WIDTH> u)
```
- Packs entry via `pack_bank_entry()` using `if constexpr` for layout
- Writes to `bank_ram[stage * NUM_BANKS + bank]`
- Writes shared tag RAM / shared u RAM / u FFs as applicable

Also:
- `writeBack(bank, stage, tag, u)` — convenience write from cached `ctr_regs`
- `writeReg(bank, slot, new_ctr)` — update single counter in cache before writeBack

### 2.5 U-bit reset (FF mode) ✅
```
reset_u(val<ResetFn::MODE_BITS> mode)
```
- Iterates all u-bit FF arrays, applies `ResetFn::apply<U_WIDTH>(u, mode)`
- `static_assert` enforces `U_STOR_FF=true`

### 2.6 U-bit decay (SRAM mode) ✅
- `should_decay()` returns `val<1>` true with probability `1/DECAY_CTR` using `std::rand()`
- Applied inside `read()` on tag miss: saturating decrement of `u_reg`
- Shared vs per-bank decay handled via `if constexpr`

### 2.7 Test harness — TODO

Build a standalone test using the `harcom_superuser` pattern from `test_harcom.cpp`:

**Unit tests**:
- Write an entry, read it back, verify all accessors return correct values
- Test tag hit and miss
- Test `reuseRead()` returns cached data
- Test write-after-read (update flow)
- Test U-bit decay triggers probabilistically
- Test U-bit FF reset modes (reset, rshift, decrement)
- Test with different parameter combinations (compile-time instantiations)

**Predictor-emulation patterns**:
- Predict flow: `read()` → `getHit()` → `getCounter()` → cycle boundary → `write()` (update)
- Allocation flow: `read()` → miss → `write()` with fresh entry
- Block reuse flow: `read()` → `reuseRead()` × (BPB-1)
- Multi-bank flow: read all banks, check hits independently

**Makefile target**: `test-tagetable` ✅ (currently runs compile test; will be updated for runtime tests)

---

## Phase 3: HARCOM Timing Verification

**Goal**: Ensure the datapath meets HARCOM rules — no two reads or writes to the same RAM in one cycle, proper fanout declarations.

### 3.1 Timing audit
- Trace through the read path: count RAM accesses per cycle per bank
- Trace through the write path: count RAM accesses per cycle per bank
- Verify read and write never hit the same bank RAM in the same cycle
  - Without ahead (`USE_AHEAD = false`): register-cache ensures read and write are different cycles
  - With ahead (`USE_AHEAD = true`): reads go to current-cycle banks, writes go to previous-cycle banks — verify no collision
- Verify `need_extra_cycle` is used when conflicts are possible

### 3.2 Fanout audit
- Every `val` used more than once must have `.fanout(hard<N>{})` or `.fo1()`/`.fo2()`
- Check all accessors — returned values need appropriate fanout hints
- Check internal temporaries during unpack/compare

### 3.3 Cycle budget
- Document expected cycles per operation:
  - Read: 1 cycle (all banks in parallel)
  - Reuse read: 0 cycles (FF only)
  - Write: 1 cycle (one bank)
  - U-bit FF reset: 0 cycles (combinational on FFs)
- Verify with `panel.next_cycle()` in test harness

### 3.4 Deliverable
- All parameter combinations pass timing rules
- Document any `need_extra_cycle` conditions

---

## Phase 4: Integration with Reference TAGE

**Goal**: Drop TageTable into the reference `predictors/tage.hpp` and verify no accuracy regression.

### 4.1 Create integration branch
- Copy `predictors/tage.hpp` → `predictors/tage_tagetable.hpp`
- Replace the raw RAM arrays with TageTable instances

### 4.2 Mapping reference TAGE to TageTable

Reference tage uses per-table:
- `gtag[NUMG]` — `ram<val<TAGW>>` → TageTable tag storage
- `gpred[NUMG]` — `ram<val<1>>` → TageTable counter (CTR_WIDTH=1 or map 1-bit pred to counter MSB)
- `ghyst[NUMG]` — `rwram<2>` → embedded in CTR_WIDTH (use CTR_WIDTH=3 to encode pred+hyst)
- `ubit[NUMG]` — `rwram<1>` → TageTable u-bit (U_WIDTH=1)

Parameters for equivalence:
```cpp
TageTable<
    TABLE_SIZE = (1 << LOGG),
    TABLE_HIST = HLEN[i],        // per-table
    TAG_WIDTH = TAGW,
    CTR_WIDTH = 3,               // 1-bit pred + 2-bit hyst encoded as 3-bit counter
    U_WIDTH = 1,
    N = 1,                       // reference tage: 1 prediction per entry
    NUM_BANKS = 1,
    USE_AHEAD = false,
    SHARED_TAG = true,
    SHARED_U = true,
    U_STORAGE = FF,              // match reference tage's rwram u-bit behavior
    USE_FF_CACHE = false         // BPB = 1, no caching
>
```

### 4.3 Adapter logic
- The reference tage does cross-table match using offset-in-tag. The TageTable doesn't know about offsets — the predictor handles tag construction/comparison.
- The predictor builds `concat(offset, htag)` as the comparison tag passed to TageTable.
- Map `update_ctr()` (saturating counter) to CTR_WIDTH counter updates.
- Map `ghyst` read-modify-write to counter update through TageTable's register-cache flow.

### 4.4 Validation
- Run both versions on the same traces
- Compare: misprediction count, energy (fJ), latency (ps)
- **Pass criterion**: identical misprediction count. Energy/latency may differ slightly due to different RAM packing but should not be worse by more than ~5%.

### 4.5 Deliverable
- `predictors/tage_tagetable.hpp` passes validation against reference `tage.hpp`
- Document any discrepancies and root causes

---

## Phase 5: Parameter Sweep

**Goal**: Sweep TageTable parameters and document HARCOM cost/latency tradeoffs.

### 5.1 Sweep dimensions

| Parameter | Sweep values |
|-----------|-------------|
| `CTR_WIDTH` | 1, 2, 3, 4 |
| `U_WIDTH` | 1, 2 |
| `N` (BPB) | 1, 2, 4, 8 |
| `NUM_BANKS` | 1, 2, 4 |
| `USE_AHEAD` | false, true |
| `SHARED_TAG` | false, true |
| `U_STORAGE` | SRAM, FF |
| `USE_FF_CACHE` | false, true |
| `TABLE_SIZE` | 256, 512, 1024, 2048 |

Not all combinations are valid — filter by constraints. Estimated valid configs: ~100-200.

### 5.2 Metrics to collect
- **Area**: total SRAM bits + total FF bits (from HARCOM)
- **Read energy**: fJ per read (from HARCOM)
- **Write energy**: fJ per write
- **Read latency**: ps critical path
- **Misprediction rate**: on reference traces (using Phase 4 integration)
- **VFS score**: full pipeline score from `vfs.py`

### 5.3 Sweep script (`scripts/sweep_tagetable.py`)
- Generate `params.yaml` per configuration
- Build + run via `make cbp && make run-cbp`
- Collect CSV output
- Filter invalid parameter combinations
- Output: `results/tagetable_sweep.csv`

### 5.4 Analysis
- Pareto front: accuracy vs area, accuracy vs energy
- Identify sweet spots for competition entry
- Document findings in `docs/tagetable_sweep_results.md`

### 5.5 Deliverable
- Sweep script + results CSV
- Pareto analysis with recommended configurations
