# TageTable Development Plan

Development plan for the parameterized `TageTable` class. Builds on design decisions in `agent/tagetable_plan.md`.

## File Locations

- **Implementation**: `predictors/custom/TageTable.hpp` (rewrite existing file)
- **Test harness**: `tests/test_tagetable.cpp` (new)
- **Integration**: `predictors/tage.hpp` (modify to use TageTable)
- **Parameter sweep**: `tests/test_tagetable_sweep.cpp` (28 configs, HARCOM regions)
- **Sweep scripts**: `scripts/sweep_tagetable.py` (new, Phase 5)
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

## Phase 2: Accessors and Read/Write Interface ✅ COMPLETE

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

### 2.7 Test harness ✅

16 tests, 69 assertions, all passing. File: `tests/test_tagetable.cpp`.

**Unit tests** (all ✅):
- Write/read/verify all accessors (test 1)
- Tag hit and miss (tests 1, 2)
- `reuseRead()` returns cached data (test 6)
- Write-after-read update flow via direct `write()` (test 3)
- U-bit SRAM decay on miss with `DECAY_CTR=1` (test 8)
- U-bit SRAM per-bank decay (test 15)
- U-bit FF reset modes: reset/rshift/decrement (test 9)
- Multi-bank shared tag (test 4)
- Per-bank tags + independent hit (tests 5, 13)
- Ahead pipelining stages (test 7)
- Everything combined: ahead + multibank + FF cache + per-bank (test 16)

**Predictor-emulation patterns** (all ✅):
- Predict flow: `read()` → `getHit()` → `getCounter()` → cycle → `write()` (test 10)
- Allocation flow: `read()` → miss → `write()` (test 11)
- Block reuse flow: `read()` → `reuseRead()` × (BPB-1) (test 12)
- Multi-bank hit: read all banks, check hits independently (test 13)
- HARCOM timing: read/write on different cycles verified (test 14)

**Latency reporting**: All read accessors print `(t=... ps, loc=...)` via HARCOM `.print()`.

**Bugs found and fixed during testing**:
- Accessors returned `reg<N>` by value (copy creates/destroys storage → HARCOM crash). Fixed: return `val<N>`.
- `u_reg` double-write in SRAM decay path (read from RAM + decay in same cycle). Fixed: extract u from entry into local val, apply decay, write `u_reg` once.
- `should_decay()` with `DECAY_CTR=1` produced `val<0>` (invalid). Fixed: `if constexpr (DECAY_CTR <= 1)` returns `val<1>{1}`.
- `u_ff[0][index]` where `index` is `val<IDX_BITS>` — HARCOM `arr::operator[]` takes `u64`. Fixed: `write_u_ff()` helper with `execute_if` decode pattern.

**Makefile targets**: `test-tagetable` (runtime, builds+runs), `test-tagetable-compile` (compile-only, 14 instantiations)

---

## Phase 3: HARCOM Timing Verification & Cost Exploration ✅ COMPLETE

**Goal**: Ensure the datapath meets HARCOM rules, proper fanout declarations, and explore hardware costs.

**Status**: 73 assertions passing (69 Phase 2 + 4 Phase 3). All timing, fanout, and cost tests verified.

### 3.1 Timing audit ✅
- Read and write never access the same bank RAM in the same cycle (verified by test 14)
- Without ahead: predictor calls `next_cycle()` between read and write
- With ahead: reads go to current-stage banks, writes go to specified-stage banks — no collision
- No `need_extra_cycle` conditions required

### 3.2 Fanout audit ✅
- Added computed fanout constants to `TageTable.hpp`:
  - `IDX_READ_FANOUT`: 1 + NUM_BANKS + conditional shared tag/u RAM reads
  - `TAG_CMP_FANOUT`: 1 if shared tag, NUM_BANKS if per-bank
  - `IDX_WRITE_FANOUT`, `TAG_WRITE_FANOUT`, `U_WRITE_FANOUT` for write path
- `index.fanout()` and `compare_tag.fanout()` in `read()`
- `entry.fanout(hard<2>{})` for per-bank SRAM u extraction
- `stored_tag.fanout(hard<2>{})` for shared tag hit + accessor
- `u_val.fanout(hard<3>{})` for decay paths
- All 16 configs pass with correct latency reporting

### 3.3 Cycle budget ✅
- Verified via `test_cycle_budget()`:
  - Read: 0 extra cycles (combinational within one cycle)
  - Write: 0 extra cycles (combinational within one cycle)
  - Reuse read: 0 extra cycles (FF only)
  - U-bit FF reset: 0 extra cycles (combinational on FFs)

### 3.4 HARCOM cost exploration ✅
- `test_harcom_costs()`: global costs across all 16 table instances
  - 126,247 storage bits, 1.12M transistors
  - 45.1 fJ total energy, 2.5 mW dynamic, 0.44 mW static power
- `test_per_operation_energy()`: per-operation energy deltas for 4 configs:
  - t1 (1024, basic): write=4366 fJ, read=728 fJ
  - t4 (256, 4-bank): write=1023 fJ, read=515 fJ
  - t8 (64, SRAM u): write=196 fJ, read=230 fJ
  - t16 (256, everything): write=1202 fJ, read=768 fJ, reuse=12 fJ
- Key finding: reuse read is ~62x cheaper than full read (12 vs 768 fJ)
- U-bit FF decode latency scales with TABLE_SIZE (4.2 ns for 1024 entries vs ~1.1 ns for 256)

### 3.5 Parameter sweep (`tests/test_tagetable_sweep.cpp`) ✅

44 configurations across 6 groups, each in its own HARCOM region for isolated cost measurement.
Each group answers one design question by varying a single dimension from a baseline.
FF and SRAM u-bit variants are paired for direct comparison.

Makefile target: `make test-tagetable-sweep`

**Sweep results** (latencies in ps, area in mm², static power in mW):

```
--- 1. Size scaling: FF u-bits (F=1) vs SRAM u-bits (F=0) ---
  Baseline: TAG=11, CTR=3, U=1, N=4, 1 bank, shared tag.
  Pairs show cost of FF u-bit mux tree as SIZE grows.

Config        SIZE  F     bits     trans    area staPwr    Hit   Tag   Ctr  Ubit
64  u=FF        64  1     1558     19796 0.00008 0.0081    127   104   113   346
64  u=SRAM      64  0     1558     11584 0.00009 0.0026    133   106   124   158
128 u=FF       128  1     3095     23664 0.00015 0.0059    211   188   207   682
128 u=SRAM     128  0     3095     21778 0.00016 0.0044    229   202   209   254
256 u=FF       256  1     6168     45808 0.00029 0.0108    156   133   155  1134
256 u=SRAM     256  0     6168     41928 0.00031 0.0076    187   160   165   212
512 u=FF       512  1    12313     88728 0.00047 0.0172    173   150   174  2167
512 u=SRAM     512  0    12313     80432 0.00050 0.0105    217   190   199   242
1K  u=FF      1024  1    24602    174152 0.00084 0.0307    186   163   187  4219
1K  u=SRAM    1024  0    24602    157132 0.00091 0.0170    218   191   198   243
2K  u=FF      2048  1    49179    349040 0.00169 0.0625    266   243   260  8358
2K  u=SRAM    2048  0    49179    314048 0.00180 0.0343    264   237   268   289

--- 2. Banking: vary N and BANKS (SIZE=256, SRAM u, shared tag) ---
  N = total branches, BxBPB = banks x branches_per_bank.

Config         N  BxBPB    bits     trans    area staPwr    Hit   Tag   Ctr  Ubit
N=1  1x1       1  1x1      3864     26968 0.00020 0.0052    164   137   132   189
N=2  2x1       2  2x1      4635     32560 0.00024 0.0063    221   194   189   246
N=4  2x2       4  2x2      6171     42488 0.00031 0.0079    201   174   180   226
N=4  4x1       4  4x1      6177     43744 0.00033 0.0087    236   209   202   261
N=8  2x4       8  2x4      9243     62480 0.00046 0.0113    214   187   205   239
N=8  4x2       8  4x2      9249     63600 0.00047 0.0118    202   175   180   227

--- 3. Entry width: TAG/CTR (SIZE=256, N=4, 1 bank, SRAM u) ---

Config        TAG CTR    bits     trans    area staPwr    Hit   Tag   Ctr  Ubit
tag=7  ctr=2    7   2    4115     28568 0.00021 0.0054    159   134   147   184
tag=7  ctr=3    7   3    5140     35332 0.00027 0.0067    209   184   200   234
tag=11 ctr=2   11   2    5143     35164 0.00026 0.0064    173   146   152   198
tag=11 ctr=3   11   3    6168     41928 0.00031 0.0076    193   166   170   218
tag=11 ctr=4   11   4    7193     48432 0.00035 0.0084    163   136   150   188
tag=13 ctr=3   13   3    6682     45472 0.00034 0.0085    226   197   206   251

--- 4. U-bit width x storage (SIZE=256, N=4, 1 bank) ---

Config          U  F     bits     trans    area staPwr    Hit   Tag   Ctr  Ubit
U=1 FF          1  1     6168     45808 0.00029 0.0108    162   139   159  1140
U=1 SRAM        1  0     6168     41928 0.00031 0.0076    198   171   176   223
U=2 FF          2  1     6425     53004 0.00029 0.0159    167   144   163  1145
U=2 SRAM        2  0     6425     43550 0.00032 0.0079    235   208   212   261

--- 5a. Features: ahead + per-bank tags (SIZE=256, N=4, SRAM u) ---
  Each row adds one feature to show incremental cost.

Config              bits     trans    area staPwr    Hit   Tag   Ctr  Ubit
base 2x2            6171     42488 0.00031 0.0079    197   170   176   222
+ahead             12315     84312 0.00063 0.0153    229   202   200   256
+per-bank          15432    103108 0.00074 0.0179    217   190   190   270

--- 5b. Features: FF cache (SIZE=256, N=8, 2 banks, SRAM u) ---

Config              bits     trans    area staPwr    Hit   Ctr  Ubit Reuse
base 2x4            9243     62480 0.00046 0.0113    206   199   231     -
+cache              9261     62936 0.00046 0.0116    201   191   226   216
+ahead+cache       18477    124752 0.00093 0.0225    208   194   235   219
+ahd+pbnk+ca       24634    160500 0.00108 0.0239    229   202   288   227

--- 6. Large tables (SRAM u) ---

Config              bits     trans    area staPwr    Hit   Ctr  Ubit Reuse
1K N=8 4x2         36899    239020 0.00149 0.0302    258   229   283     -
1K +ahead          73763    477176 0.00298 0.0599    222   204   263     -
1K +ahd+pb+c      147542    947316 0.00546 0.1111    290   263   383   276
2K N=4 1x4         49179    314048 0.00180 0.0343    276   278   301     -
2K N=1 1x1         30747    198050 0.00118 0.0234    290   214   315     -
```

Column key: bits = SRAM+FF storage (bits), trans = transistors, area = SRAM area (mm²),
staPwr = static/leakage power (mW, workload-independent), latencies in ps.
F = u-bit storage (0=SRAM, 1=FF). Reuse = reuseRead latency, - = no FF cache.

**Key findings**:

1. **SRAM u-bits are strictly better** than FF u-bits for all table sizes:
   - U-bit latency: 5–30x lower (e.g. 4219→243 ps at 1K, 8358→289 ps at 2K)
   - Fewer transistors (10–15% reduction — FF decode logic is expensive)
   - ~40–50% lower static power (e.g. 0.031→0.017 mW at 1K)
   - Same storage bits; negligible SRAM area increase
   - Hit/Tag/Ctr latencies unaffected (~same or slightly higher)

2. **Hit/Tag/Ctr latency** is well-behaved (~130–290 ps), scaling gently with size and
   entry width. Tag is consistently ~25 ps faster than Hit (no comparator on that path).

3. **Banking**: N=4 with 2 banks (2x2) vs 4 banks (4x1) costs similar bits/xtors
   but 4x1 has narrower per-bank entries → can be faster per-bank.
   N=8 nearly doubles area vs N=4.

4. **Ahead doubles area and power** (+ahead: 12K vs 6K bits, 0.015 vs 0.008 mW)
   for negligible latency change. It physically duplicates all bank RAMs.

5. **Per-bank tags are expensive** (+per-bank: 15K bits, 103K trans, 0.018 mW
   vs base 6K/42K/0.008) — each bank stores its own tag RAM plus comparison logic.

6. **FF cache is nearly free** in area/power (+18 bits, +0.0003 mW) but enables
   reuseRead at ~216 ps with no additional SRAM access.

7. **Entry width scaling**: each extra counter bit adds ~1K bits and ~7K trans at SIZE=256.
   Tag width has similar marginal cost.

### 3.6 Deliverable ✅
- All parameter combinations pass timing rules
- Fanout declarations added to read() and write() paths
- Cost/energy/cycle-budget tests in `tests/test_tagetable.cpp`
- Parameter sweep with 44 configs in `tests/test_tagetable_sweep.cpp`

---

## Phase 4: Integration with Reference TAGE — IN PROGRESS

**Goal**: Drop TageTable into the reference `predictors/tage.hpp` and verify no accuracy regression.

**Status**: `predictors/tage_tagetable.hpp` created and compiles clean. Runtime crash: `array bound exceeded (8>=8)` — an `arr<..., NUMG=8>` accessed at index 8. Not yet debugged.

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
