# TageTable Development Plan

Development plan for the parameterized `TageTable` class. Builds on design decisions in `agent/tagetable_plan.md`.

## File Locations

- **Implementation**: `predictors/custom/TageTable.hpp` (rewrite existing file)
- **Test harness**: `tests/test_tagetable.cpp` (new)
- **Integration**: `predictors/tage.hpp` (modify to use TageTable)
- **Sweep scripts**: `scripts/sweep_tagetable.py` (new)
- **Design docs**: `agent/tagetable_plan.md`, `agent/tagetable_questions.md`

---

## Phase 1: Storage Primitives

**Goal**: Declare all HARCOM private members, parameterized by the template parameter list. No logic yet — just storage.

### 1.1 Template declaration and static_asserts
- Full template parameter list (see `agent/tagetable_plan.md` for confirmed parameters)
- Computed constants: `BPB = N / NUM_BANKS`, `PHYS_BANKS = NUM_BANKS * (USE_AHEAD ? 2 : 1)`, entry widths
- All constraints:
  - `N % NUM_BANKS == 0`
  - `!SHARED_U || SHARED_TAG` (per-branch U requires per-bank tags)
  - `!USE_FF_CACHE || BPB > 1` (caching only useful with multiple branches per bank)

### 1.2 SRAM storage
- **Counter RAM**: `ram<val<BPB * CTR_WIDTH>, TABLE_SIZE>` per physical bank (or wider if tag/U packed in)
- **Tag RAM**: depends on `SHARED_TAG`
  - `SHARED_TAG = true`: tag stored once per entry (separate small RAM or packed into one bank)
  - `SHARED_TAG = false`: tag packed per-bank entry
- **U-bit RAM** (when `U_STORAGE = SRAM`): packed into entry alongside counters
- Entry layout varies by parameter combination — use `if constexpr` to select packing

### 1.3 FF storage
- **U-bit FFs** (when `U_STORAGE = FF`):
  - `SHARED_U = true`: `arr<reg<U_WIDTH>, TABLE_SIZE>`
  - `SHARED_U = false`: `arr<arr<reg<U_WIDTH>, NUM_BANKS>, TABLE_SIZE>`
- **Cache FFs** (when `USE_FF_CACHE = true`):
  - Per-bank: `arr<reg<CTR_WIDTH>, BPB>` for cached counters
  - `reg<TAG_WIDTH>` for cached tag (per bank or shared depending on `SHARED_TAG`)
  - `reg<U_WIDTH>` for cached u-bit
  - `reg<clog2(TABLE_SIZE)>` for cached index
  - `reg<1>` for cached hit

### 1.4 Deliverable
- Compiles with all parameter combinations
- No methods yet — just verify storage instantiation doesn't break HARCOM

---

## Phase 2: Accessors and Read/Write Interface

**Goal**: Implement the read/write methods and accessors. Build a test harness.

### 2.1 Read path
```
read(val<clog2(TABLE_SIZE)> index, val<TAG_WIDTH> compare_tag)
```
- Read all banks at `index`
- Unpack entries: extract tag, counters, u-bit per bank
- Compare stored tag(s) against `compare_tag`, store hit result
- If `USE_FF_CACHE`: cache all fields in registers

### 2.2 Reuse read (when `USE_FF_CACHE = true`)
```
reuseRead(bank, slot) -> val<CTR_WIDTH>
```
- Return cached counter from FFs, no RAM access

### 2.3 Accessors
- `getHit(bank)` → `val<1>` (or just `getHit()` when `SHARED_TAG = true`)
- `getTag(bank)` → `val<TAG_WIDTH>`
- `getCounter(bank, slot)` → `val<CTR_WIDTH>`
- `getU(bank)` → `val<U_WIDTH>` (or `getU()` when `SHARED_U = true`)

### 2.4 Write path
```
write(val<clog2(TABLE_SIZE)> index, bank, val<TAG_WIDTH> tag,
      arr<val<CTR_WIDTH>, BPB> counters, val<U_WIDTH> u)
```
- Pack fields into entry format
- Write to specified bank at specified index
- For `U_STORAGE = FF`: write u-bit to FF array instead

### 2.5 U-bit reset (FF mode only)
```
reset_u(val<MODE_BITS> mode)
```
- Apply `ResetFn` to all FF u-bit entries with given mode

### 2.6 U-bit decay (SRAM mode only)
- Probabilistic decrement on tag miss inside read path
- Uses `std::rand()` with `1/DECAY_CTR` probability

### 2.7 Test harness (`tests/test_tagetable.cpp`)

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

**Add to Makefile**:
```makefile
test-tagetable: tests/test_tagetable.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lz
```

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
