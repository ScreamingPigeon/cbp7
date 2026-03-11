# Custom TAGE Predictor Development Plan

Development plan for `predictors/custom/Tage.hpp` — a fully parameterized TAGE predictor
built on `TageTable` instances. Design decisions documented in
`memory/tage_predictor_design.md`.

## File Locations

- **Implementation**: `predictors/custom/Tage.hpp` (rewrite existing stub)
- **TageTable**: `predictors/custom/TageTable.hpp` (modifications needed)
- **Tests**: `tests/test_tage.cpp` (new)
- **Reference**: `predictors/tage.hpp` (validation baseline)
- **Ahead reference**: `predictors/gshareN_ahead.hpp` (ahead pipeline reference)
- **Design docs**: `memory/tage_predictor_design.md`

## TageTable Modifications (Pre-Phase, as needed)

TageTable changes required by the new predictor. Can be done incrementally as phases need them.

- **LFSR replacement**: Replace `std::rand()` in `should_decay()` with a proper LFSR register
- **Dynamic decay threshold**: Change `DECAY_CTR` from template param to runtime `val<>` input
  on `read()`. Decay fires when LFSR > threshold (predictor controls aggressiveness)
- **HARCOM region per instance**: Constructor accepts region name, all internal storage
  declared within that region
- **Validation**: Existing `test_tagetable.cpp` must still pass after each modification

---

## Phase 1: Scaffolding and Primitives

**Goal**: Full predictor skeleton with all template parameters, config structs, constexpr
functions, private members, and TageTable instantiation. Both direct and ahead
specializations stubbed. Everything compiles, nothing runs.

### 1.1 Config structs

```cpp
// Per-table config — arrays indexed by table number
struct DefaultTableConfig {
    static constexpr u64 NUM_TABLES = 8;
    static constexpr u64 MINHIST = 2;
    static constexpr u64 MAXHIST = 100;

    static constexpr std::array<u64, NUM_TABLES> TABLE_SIZE = uniform<...>(2048);
    static constexpr std::array<u64, NUM_TABLES> TAG_WIDTH  = uniform<...>(11);
    static constexpr std::array<u64, NUM_TABLES> CTR_WIDTH  = uniform<...>(3);
    static constexpr std::array<u64, NUM_TABLES> U_WIDTH    = uniform<...>(1);
    static constexpr std::array<u64, NUM_TABLES> HIST_LEN   = geometric<...>();
};

// Allocation policy
struct DefaultAllocConfig {
    static constexpr u64 MAX_ALLOC = 1;
    static constexpr bool USE_ALLOC_FILTER = false;
    static constexpr bool PROTECT_RECENT_ALLOC = false;
};
```

### 1.2 Constexpr helpers
- `uniform<T, N, V>()` — generate constant array
- `geometric<N, MIN, MAX>()` — generate geometric history length series
- `constexpr_max(std::array)` — max element of array (for MAX_TAG_WIDTH etc.)
- Per-table `static constexpr` functions where needed (hist_len, table_size, etc.)

### 1.3 Computed constants
- `MAX_TAG_WIDTH`, `MAX_CTR_WIDTH`, `MAX_U_WIDTH` across all tables
- `LOG_FETCH_WIDTH = clog2(FETCH_WIDTH)`
- `HTAGBITS` = `MAX_TAG_WIDTH - LOG_FETCH_WIDTH` (when BR_P_ENTRY=1)
- `MAXHIST` = max of HIST_LEN array

### 1.4 Table tuple generation
- `TableAt<Cfg, GlobalParams, I>` type alias → TageTable with per-table params from
  config arrays at index I, plus global params (BR_P_ENTRY, NUM_BANKS, USE_AHEAD, etc.)
- `make_table_tuple_t<Cfg, GlobalParams>` → `std::tuple<TableAt<..., 0>, TableAt<..., 1>, ...>`

### 1.5 TageImpl partial specializations
```cpp
template <bool USE_AHEAD, typename TableCfg, typename AllocCfg, /* global params */>
struct TageImpl;

template <typename TableCfg, typename AllocCfg, ...>
struct TageImpl<false, TableCfg, AllocCfg, ...> : predictor {
    // Direct mode — all members declared, methods stubbed
};

template <typename TableCfg, typename AllocCfg, ...>
struct TageImpl<true, TableCfg, AllocCfg, ...> : predictor {
    // Ahead mode — all members declared (incl. pipeline regs), methods stubbed
};

template <typename TableCfg = DefaultTableConfig, ...>
using Tage = TageImpl<USE_AHEAD, TableCfg, ...>;
```

### 1.6 Private members (both specializations)
- **Table tuple**: `Tables tables;` (the generated tuple of TageTable instances)
- **Global history**: `reg<MAXHIST> global_history;`
- **Path history**: `std::conditional_t<USE_PATH_HIST, reg<PATH_HIST_WIDTH>, EmptyMember>`
- **Fold registers**: per-table index and tag folds
- **P2 result registers** (max-width uniform):
  - `arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> readt;`
  - `arr<reg<1>, NUM_TABLES> readc;` (prediction bits)
  - `arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> readh;` (hysteresis)
  - `arr<reg<1>, NUM_TABLES> readu;`
- **Match registers**: `match`, `match1`, `match2` per offset
- **P2 prediction registers**: `pred1`, `pred2`, `p2`
- **Meta-prediction**: `meta` counter/table (conditional on USE_META)
- **U-bit management**: `uctr` (epoch mode), decay threshold reg (SRAM mode)
- **Block state**: `block_entry`, `block_size`, `num_branch`, `branch_offset`, `branch_dir`
- **Bimodal**: `bim[FETCH_WIDTH]`, `bhyst[FETCH_WIDTH]` (UPDATE_ONLY zone)
- **P1**: gshare/bimodal tables (conditional on P1_USE_GSHARE), own HARCOM region

Additional for ahead specialization:
- **Pipeline registers**: `index[2]` per table, `XB[2]`, `XL`, `path`
- **Cached bank predictions**: `block_pred[2]`

### 1.7 Method stubs
All predictor interface methods stubbed with `return val<1>{0};` or empty bodies:
- `predict1`, `reuse_predict1`, `predict2`, `reuse_predict2`
- `update_condbr`, `update_cycle`

### 1.8 Compile test
- Makefile target: `test-tage-compile`
- Instantiate with DefaultTableConfig (uniform params)
- Instantiate with a custom config (non-uniform TABLE_SIZE, TAG_WIDTH)
- Both direct and ahead specializations
- Must compile clean with `-Werror`

---

## Phase 2: Accessors and Routines ✓ COMPLETE

**Goal**: Implement all dataflow as callable functions invoked by predict/update methods.
Build the predict2, reuse_predict2, and update_cycle flows from these building blocks.

**Status**: All sub-phases (2.1-2.13) implemented. Direct-mode Tage runs full traces.
Mispredictions match reference exactly when `-DREAD_WRITE_RAM` is used (208,030 on gcc trace).
Without it, 219,832 mispredictions due to TageTable write gating by `extra_cycle`
(HARCOM `ram<>` enforces single access per cycle; bank_ram is read in predict2 and
written in update_tables). See `memory/tage_phase2.md` for full details.

### 2.1 Hash functions
- `compute_index(table_idx, lineaddr, fold)` → per-table index
- `compute_tag(table_idx, lineaddr, fold)` → per-table tag
  - When BR_P_ENTRY=1: returns hashed tag bits (predictor prepends offset)
  - When BR_P_ENTRY>1: returns full tag
- Optionally mix path history (if constexpr USE_PATH_HIST)

### 2.2 History update routines
- `update_global_history(taken, next_pc)` — shift global history, update path history
- `update_folds()` — update all fold registers for all tables

### 2.3 Table read routine
- `read_tables(lineaddr)` — `static_loop<NUM_TABLES>`:
  - Compute index[i] and tag[i]
  - Call `table.read(index, tag)`
  - Extract pred, hyst, u from counter encoding into result registers
  - Store tags in `readt`

### 2.4 Match logic
- `compute_matches()` — build match masks from tag comparisons
  - BR_P_ENTRY=1: per-offset matching with offset-in-tag
  - BR_P_ENTRY>1: per-table hit matching
- `find_provider()` — longest match → `match1`, prediction → `pred1`
- `find_altpred()` — second longest match → `match2`, prediction → `pred2`

### 2.5 Meta-prediction
- `compute_meta_prediction()` — alt-on-newly-allocated logic
  - Global counter or PC-indexed table
  - Returns final p2 prediction per offset

### 2.6 predict2 flow
Wire together: `read_tables` → `compute_matches` → `find_provider` → `find_altpred`
→ `compute_meta_prediction` → block entry masking → return taken

### 2.7 reuse_predict2 flow
BR_P_ENTRY=1: reuse cached p2 with shifted block_entry
BR_P_ENTRY>1: call `reuseRead()` on tables, recompute prediction for next slot

### 2.8 Allocation routine
- `allocate_entries(mispredict, last_offset, ...)` — post-provider search, candidate
  selection, randomization, fallback u-clear

### 2.9 Counter/u-bit update routine
- `update_tables(...)` — combined TageTable write per table:
  tag (new on alloc, old otherwise), pred, hyst, u packed into single write()
- `update_u_bits(...)` — u-bit increment/clear logic, adaptive decay threshold

### 2.10 Bimodal update routine
- `update_bimodal(...)` — conditional pred/hyst updates

### 2.11 Meta update routine
- `update_meta(...)` — meta counter update based on primary vs alt correctness

### 2.12 update_cycle flow
Wire together: branch recording → allocation → counter/u-bit update → bimodal update
→ P1 update → meta update → history update → extra_cycle signaling

### 2.13 Functional test
- Instantiate direct-mode Tage with default params
- Run predict2 → update_cycle sequence on synthetic data
- Verify no crashes, no HARCOM assertion failures

---

## Phase 3: P1 Bimodal/Gshare Predictor

**Goal**: Implement the fast first-level predictor. Own HARCOM region.

### 3.1 P1 storage (from Phase 1 members)
- Gshare mode: P1 history register, pred table, hyst table (UPDATE_ONLY)
- Bimodal mode: pred table only
- `if constexpr (P1_USE_GSHARE)` selects

### 3.2 predict1 / reuse_predict1
- Compute P1 index, read table, cache results
- Block entry masking, reuse_prediction signaling

### 3.3 P1 update (within update_cycle)
- Update pred when P1/P2 disagree and hysteresis is weak
- Update hysteresis

### 3.4 Test
- P1 predictions are non-trivial (not all 0 or all 1)
- P1 + P2 integration: full predict1 → predict2 → update flow

---

## Phase 4: HARCOM Cost and Latency Validation

**Goal**: Verify hardware costs and access latencies. Test harness similar to
`tests/test_tagetable.cpp` and Makefile target.

### 4.1 Test harness (`tests/test_tage.cpp`)
- Makefile target: `test-tage`
- Instantiate Tage with known parameters
- Run a sequence of predict/update cycles
- Print per-region HARCOM costs (bits, transistors, area, power)
- Print per-table costs (using per-table HARCOM regions)
- Print access latencies: P1 predict, P2 predict, update

### 4.2 Per-operation energy
- Measure energy delta for: P1 predict, P2 predict, P2 reuse predict, update cycle
- Compare against reference tage for equivalent params

### 4.3 Cost summary
- Total storage bits, transistors, area, static power
- Breakdown by component: P1, bimodal, tage tables, meta-prediction

---

## Phase 5: Timing Access Validation

**Goal**: Test plan + accessor tests to ensure no HARCOM timing violations
(same-cycle read/write conflicts, missing fanout, etc.)

### 5.1 Test plan
- Document all read/write access patterns per cycle:
  - predict1 cycle: P1 table reads
  - predict2 cycle: bimodal reads, all TageTable reads, meta read
  - update_cycle: TageTable writes, bimodal writes, P1 writes, meta write, history update
- Identify potential conflicts (read and write to same RAM in same cycle)

### 5.2 Accessor tests
- Test each routine from Phase 2 in isolation
- Verify `next_cycle()` between predict and update phases
- Verify no "array bound exceeded" or double-write crashes
- Test with USE_AHEAD=false: read/write separation across cycles
- Test with various param combos: different NUM_TABLES, TAG_WIDTH, CTR_WIDTH

### 5.3 Cycle budget verification
- Confirm extra_cycle usage matches expectations
- Count extra_cycles on a trace, compare with reference

---

## Phase 5.5: Validate Ahead Prediction

**Goal**: Implement and validate the ahead (1-branch-ahead) specialization.

### 5.5.1 Implement ahead predict path
- Pipeline registers (index[2], XB[2], XL, path)
- Path banking: read all PATHS banks, select via `path ^ XB`
- Lane XOR for even storage utilization
- predict1/predict2 using previous block's cached reads

### 5.5.2 Implement ahead update path
- Write uses `index[1]` and `path ^ XB[1]`
- Reuse shared update routines from Phase 2 (allocation, u-bit, meta, etc.)

### 5.5.3 Factor shared logic
- Extract common update logic into free functions or CRTP base
- Both direct and ahead specializations use the same helpers

### 5.5.4 Test
- Ahead mode compiles and runs without crashes
- HARCOM timing: verify read (current stage) and write (previous stage) don't conflict
- Compare accuracy: ahead vs direct on same traces
- Compare latency: ahead should have lower P2 critical path

---

## Phase 6: Reference Equivalence Validation

**Goal**: Configure custom Tage to match reference `tage.hpp` parameters exactly
and verify identical behavior.

### 6.1 Parameter mapping
- Map reference tage defaults to custom Tage config:
  NUMG=8, LOGG=11, TAGW=11, GHIST=100, LOGB=12, LOGP1=14, GHIST1=6
- BR_P_ENTRY=1, NUM_BANKS=1, USE_AHEAD=false, U_STOR_FF=false

### 6.2 Trace validation
- Run both predictors on the same traces
- **Pass criterion**: identical misprediction count
- Energy/latency may differ slightly (different RAM packing) — within ~5%

### 6.3 Divergence debugging
- If mispredictions differ: add logging to both predictors, compare per-branch decisions
- Bisect to find first divergence point

---

## Phase 7: Performance Counters and Logging

**Goal**: Add runtime performance monitoring for parameter tuning.

### 7.1 Predictor-level counters
- Mispredictions (total, P1-only, P2-only)
- Allocation attempts, successes, failures
- U-bit decay events, epoch resets
- Provider table distribution (which table provides prediction, histogram)
- Meta-prediction: alt-used count, alt-correct count
- Extra cycle count

### 7.2 Per-table counters (in TageTable)
- Hits, misses, allocations, evictions
- Read/write counts
- U-bit distribution (how many entries have u=0, u=1, etc.)

### 7.3 Logging infrastructure
- Gated by `-DTAGE_MONITOR` compile flag (+ `-DCHEATING_MODE`)
- Periodic snapshots (every N branches)
- CSV output for analysis
- Per-instance logs using HARCOM region names

### 7.4 Summary stats
- End-of-trace summary printed to stderr
- Per-table accuracy contribution

---

## Phase 8: Loop Predictor (Deferred)

**Goal**: Add optional loop predictor component.

- TODO: Design loop predictor table (loop count, current iteration, confidence)
- TODO: `USE_LOOP_PRED` bool param, `LOOP_TABLE_SIZE` param
- TODO: Override TAGE prediction when loop predictor is confident
- TODO: Own HARCOM region
- TODO: Validate on loop-heavy traces
