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
- **Comparison tool**: `scripts/compare_predictors.sh` (side-by-side predictor comparison)

## TageTable Modifications (Pre-Phase, as needed)

TageTable changes required by the new predictor. Can be done incrementally as phases need them.

- **LFSR replacement**: Replace `std::rand()` in `should_decay()` with a proper LFSR register
- **Dynamic decay threshold**: Change `DECAY_CTR` from template param to runtime `val<>` input
  on `read()`. Decay fires when LFSR > threshold (predictor controls aggressiveness)
- **HARCOM region per instance**: Constructor accepts region name, all internal storage
  declared within that region
- **Direct-wire mode** ✓: When `USE_FF_CACHE=false`, staging registers are conditionally
  eliminated via `std::conditional_t<USE_FF_CACHE, StagingRegs, EmptyMember>`. `read()` uses
  templated output parameters to write directly into caller-provided registers (no intermediate
  staging). Cascaded index fanout eliminated — caller provides all fanout, TageTable skips inner
  fanout when `USE_FF_CACHE=false`.
- **Validation**: Existing `test_tagetable.cpp` must still pass after each modification

---

## Phase 1: Scaffolding and Primitives ✓ COMPLETE

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

**Post-Phase-2 optimizations** (documented in `memory/latency_optimization.md`):
- Precise fanout declarations: reduced P2 latency from 11.91 → 2.38
- TageTable direct-wire mode (USE_FF_CACHE=false): eliminated staging registers, P2 latency 2.38 → 2.34
- Cascaded index fanout eliminated
- Current gap to reference (2.34 vs 1.86) attributed to: packed counter shift+mask, `ram` vs `rwram`, P1 stub impact

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

## Phase 3: P1 Bimodal/Gshare Predictor ✓ COMPLETE

**Goal**: Implement the fast first-level predictor.

**Status**: P1 gshare implemented and working. P1 latency matches reference (0.937).
P1/P2 disagreement detection and hysteresis-based P1 update in update_cycle working.

---

## Phase 3.5: rwram Migration for TageTable ✓ COMPLETE

**Goal**: Split packed counter into separate RAMs, use rwram for hyst/u-bits.

**Status**: TageTable now has 4 separate RAMs: `tag_ram` (plain ram), `pred_ram` (plain ram),
`hyst_ram` (rwram), `u_ram` (rwram). Same-cycle read+write works without extra_cycle gating.
Mispredictions match reference without `-DREAD_WRITE_RAM`.

## Phase 3.6: Flatten TageImpl for EPI/Latency Parity ✓ COMPLETE

**Goal**: Eliminate 2x EPI gap (3826 vs 1868 fJ) caused by member function overhead.

**Status**: Rewrote Direct TageImpl as flat inline code matching tage_tt structure. All 8
member functions (read_tables, compute_matches, etc.) eliminated. AllocResult struct removed.
Bimodal write gating fixed to match reference (no extra_cycle guard).

**Results** (gcc trace, 1K warmup, 5M measure):
- P2 latency: 2.42 → 1.86 (matches reference exactly)
- EPI: 3826 → 1906 fJ (within 2% of reference 1868)
- Mispredictions: 26,402 → 26,164 (better than reference 26,870)

---

## Phase 3.7: Non-Uniform Table Parameters ✓ COMPLETE

**Goal**: Support per-table TABLE_SIZE, TAG_WIDTH, CTR_WIDTH, HYST_WIDTH, U_WIDTH.
Real TAGE designs use increasing table sizes with longer history lengths (more entries
to reduce aliasing for longer patterns).

**Status**: Replaced `Table0 table[NUM_TABLES]` (C array of uniform type) with
`typename Base::Tables tables` (std::tuple of per-table TageTable types). All 6 table
access loops converted from `for` to `static_loop<NUM_TABLES>` with `std::get<I>(tables)`.
Zero EPI/latency regression (EPI=1903, P2=1.86). All 9 compile test instantiations pass.

### 3.7.1 Replace Table0 array with tuple
- Change `Table0 table[NUM_TABLES]` back to `typename Base::Tables tables;`
- Access via `std::get<I>(tables)` inside `static_loop<NUM_TABLES>`
- All read/write loops become `static_loop` instead of `for`

### 3.7.2 Widen result registers to MAX_* widths
- `readt`: `arr<reg<MAX_TAG_WIDTH>, NUM_TABLES>` (already done)
- `readc`: `arr<reg<MAX_CTR_WIDTH>, NUM_TABLES>` (already done)
- `readh`: `arr<reg<max(1, MAX_HYST_WIDTH)>, NUM_TABLES>` (already done)
- `readu`: `arr<reg<MAX_U_WIDTH>, NUM_TABLES>` (already done)
- `gindex`: `arr<reg<MAX_IDX_BITS>, NUM_TABLES>` — index width varies with table size

### 3.7.3 Inline all table access via static_loop
- predict2 reads: `static_loop<NUM_TABLES>([&]<u64 I>() { auto &t = std::get<I>(tables); ... })`
- update_cycle writes: same pattern, 4 separate static_loops (tag, u, pred, hyst)
- CRITICAL: keep all logic inline (no helper functions that take val params)

### 3.7.4 EPI regression test
- After conversion, verify EPI stays ≤1961 (within 5% of reference)
- If EPI regresses: the cause is static_loop/tuple overhead, need to investigate
- The Phase 3.6 finding was that member functions caused the 2x EPI, not static_loop itself
  (static_loop was eliminated as a cause earlier). So this should be safe.

### 3.7.5 Functional test with non-uniform config
- CustomTableConfig with TABLE_SIZE = {512, 1024, 2048, 4096}
- Verify each table actually has different RAM sizes (check HARCOM storage accounting)
- Run trace and verify no crashes

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

## Phase 6: Reference Equivalence Validation ✓ SUBSTANTIALLY COMPLETE

**Goal**: Configure custom Tage to match reference `tage.hpp` parameters exactly
and verify identical behavior.

**Status**: Default `Tage<>` with default config achieves:
- EPI within 2% of reference (1906 vs 1868 fJ)
- P2 latency matches exactly (1.86)
- P1 latency matches (0.937)
- Mispredictions slightly better than reference (26,164 vs 26,870)

Remaining: multi-trace validation (`make compare-all`), investigate small misprediction
difference (likely from bimodal write gating change or minor fanout differences).

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

---

## Tooling (Available)

- **`scripts/compare_predictors.sh`**: Side-by-side comparison of two predictors.
  - Single trace: `make compare TRACE=path/to/trace.gz PRED_A='Tage<>' PRED_B='tage<>'`
  - All traces: `make compare-all TRACE_DIR=path/to/traces/`
  - Metrics: MPKI, IPC, VFS, EPI, latency, power, area, transistors, mispredictions
    (short/long/total), extra cycles, blocks, energy/instr, and all VFS components.
  - Directory mode: per-trace table + aggregate averages + overall VFS.

## Priority Order

1. **Phase 4-5 (validation)** — cost/timing verification, multi-trace comparison
2. **Phase 5.5 (ahead mode)** — pipeline optimization for lower P2 latency
3. **Phase 7 (monitoring)** — performance counters for parameter tuning
4. **Phase 8 (loop predictor)** — accuracy improvement

## Metrics Snapshot (gcc trace, 1K warmup, 5M measure)

| Metric | `Tage<>` (custom) | `tage<>` (reference) |
|--------|-------------------|---------------------|
| P2 mispredictions | 26,164 | 26,870 |
| P1 latency | 0.937 | 0.937 |
| P2 latency | 1.86 | 1.86 |
| Energy/instr (fJ) | 1,903 | 1,868 |

---

## Parameter Sweep Results (2026-03-13)

### Infrastructure

Three sweep scripts in `scripts/`:
- **`gaussian_sweep.py`** — broad Gaussian sampling around defaults
- **`tage_sweep.py`** — systematic tier-based sweep (has its own `TageConfig`, not using `sweep_common`)
- **`coordinate_descent.py`** — greedy local optimizer using `SweepConfig` from `sweep_common.py`

Coordinate descent: iteratively sweep one parameter at a time through its full domain, greedily accept the best value, repeat until no parameter improves in a full iteration. Checkpointing via `checkpoint.json`; eval cache backed by results CSV.

### Coordinate Descent (20 representative traces)

**Baseline `Tage<>`**: VFS=0.8763, MPKI=6.00, EPI=1182, P2=2 cycles
**Best config `5029c8a8`**: VFS=0.8870 (+1.2%), MPKI=5.85, EPI=1245, P2=2 cycles

138 configs evaluated across ~2.5 iterations (converged partway through iter 3).

#### Changes from baseline

| Param | Baseline | Best | Rationale |
|-------|----------|------|-----------|
| `tag_width` | 11 | 12 | Fewer aliasing collisions |
| `hyst_width` | 2 | 3 | More stable counters |
| `minhist` | 2 | 8 | Shortest table uses more history |
| `maxhist` | 100 | 130 | Longer history range |
| `p1_table_size` | 16384 | 32768 | 2× P1 gshare improves fast prediction |
| `p1_hist` | 6 | 10 | More P1 history bits |
| `metabits` | 4 | 2 | Less meta inertia, faster adaptation |
| `metapipe` | 2 | 1 | Single meta pipeline stage |

Unchanged: `num_tables` (8), `base_size` (2048), `size_ratio` (1.0), `bimodal_size` (8192), `decay_ctr` (2048), `ctr_width` (1), `u_width` (1).

#### Parameter sensitivity (VFS spread across swept domain)

| Sensitivity | Params | Spread |
|-------------|--------|--------|
| **Critical** | `maxhist` (0.260), `p1_table_size` (0.257), `tag_width` (0.165) | >0.10 |
| **Moderate** | `base_size` (0.075), `num_tables` (0.045), `bimodal_size` (0.035), `hyst_width` (0.034), `metabits` (0.026) | 0.02–0.08 |
| **Low** | `minhist` (0.013), `p1_hist` (0.007), `decay_ctr` (0.006), `u_width` (0.004), `metapipe` (0.003) | <0.02 |

#### Key findings

1. **P1 is the biggest lever** — 2× P1 table + longer P1 history gave repeated gains across iterations. P1 accuracy directly determines L1 latency and pipeline utilization.
2. **History range matters** — maxhist=130 and minhist=8 improved accuracy across diverse workloads.
3. **Tag width has a cliff** — below 10 bits, aliasing destroys accuracy. 12 is the sweet spot.
4. **Core table geometry already optimal** — num_tables, base_size, size_ratio unchanged from baseline.
5. **Meta predictor prefers small/fast** — metabits 4→2 and metapipe 2→1 both helped.
6. **Multi-trace shifts priorities** — single-trace gaussian sweep favored smaller tables (EPI savings); 20-trace descent kept base_size=2048 (accuracy matters more across diverse workloads).
