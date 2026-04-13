# TAGE Parameter Union

Union of all parameters across the three TAGE implementations:
- `predictors/tage.hpp` — reference TAGE (`tage<>`)
- `predictors/custom/Tage.hpp` — custom TAGE (`Tage<>` / `TageImpl`)
- `predictors/custom/TageDirect.hpp` — rank-indexed TAGE (`TageDirect<>` / `TageDirectImpl`)

Legend: **bold** = present, `—` = absent/not applicable, *(italic)* = hardcoded/derived

---

## 1. Structural / Geometry

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| `NUM_TABLES` / `NUMG` | **NUMG=8** | **via TableCfg (8)** | **via TableCfg (8)** | Number of TAGE tagged tables |
| `FETCH_WIDTH` / `LINEINST` | *LINEINST = 2^(LOGLB-2)* | **FETCH_WIDTH=16** | **LINEINST=16** | Instructions per fetch line |
| `LOGLB` | **6** | — | — | Log2(line bytes); derives LINEINST |
| `N` (max branches) | *=LINEINST* | *=FETCH_WIDTH* | **N=7** | Max conditional branches per block |
| `LANES` | *=LINEINST* | *=FETCH_WIDTH* | *=next_pow2(N)* | Parallel prediction lanes |
| `BR_P_ENTRY` | *=1* (offset in tag) | **1** | — | Branches per RAM entry (1 or FW) |
| `NUM_BANKS` | *=1* | **1** | — | Branch-slot banks per table |

## 2. Table Config (per-table arrays)

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| `TABLE_SIZE` | *uniform 2^LOGG* | **per-table via SizeFn** | **per-table via SizeFn** | Entries per table |
| `LOGG` | **11** | — | — | Log2(table size), uniform |
| `TAG_WIDTH` / `TAGW` | *uniform TAGW* | **per-table via TagFn** | **per-table via TagFn** | Tag width in bits |
| `CTR_WIDTH` | *=1* | **per-table (uniform)** | **per-table (uniform)** | Prediction counter width |
| `HYST_WIDTH` | *=2* | **per-table (uniform)** | **per-table (uniform)** | Hysteresis counter width |
| `U_WIDTH` | *=1* | **per-table (uniform)** | **per-table (uniform)** | Usefulness counter width |
| `HIST_LEN` | *geometric(MINHIST,GHIST)* | **per-table via HistSeries** | **per-table via HistSeries** | History length per table |
| `MINHIST` | *=2* | **via TableCfg (2)** | **via TableCfg (8)** | Minimum history length |
| `MAXHIST` / `GHIST` | **100** | **via TableCfg (100)** | **via TableCfg (100)** | Maximum history length |
| `SIZE_RATIO` | — | **via SweepTableConfig (4)** | **via TDTableConfig (1)** | Geometric size ratio across tables |
| `HistSeries` | *GEOMETRIC only* | **enum (GEO/QUAD/SUPEREXP/ROS)** | **enum (GEO/QUAD/SUPEREXP/ROS)** | History length series shape |
| `TagFn` | — | **functor (UniformTag default)** | **functor (GradedTag default)** | Per-table tag width functor |
| `SizeFn` | — | **functor (GeoSize default)** | **functor (GeoSize default)** | Per-table size functor |

### Available Tag Functors (Tage<>, TageDirect<>)
- `UniformTag<TAG>` — same width for all tables
- `GradedTag<TAG_LONG, TAG_SHORT>` — linear interpolation
- `StepTag<TAG_HI, TAG_LO, SPLIT>` — step function
- `LogTag<BASE, SCALE>` — logarithmic scaling

### Available Size Functors (Tage<>, TageDirect<>)
- `UniformSize<SIZE>` — same size for all
- `GeoSize<SIZE, RATIO>` — geometric scaling
- `InvGeoSize<SIZE, RATIO>` — inverse geometric (TageDirect only)
- `StepSize<S_HI, S_LO, SPLIT>` — step function
- `SqrtHistSize<BASE>` — sqrt of history-based
- `ConstBitsSize<TOTAL, TAG, CTR, HYST, U>` — fixed total bits (Tage<> only)

## 3. P1 (Base Predictor)

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| `P1_USE_GSHARE` | *=true* | **true** | **true** | Use gshare for P1 |
| `P1_TABLE_SIZE` / `LOGP1` | **2^LOGP1 (2^14)** | **4096** | **4096** | P1 table entries (total) |
| `P1_HIST` / `GHIST1` | **6** | **6** | **6** | P1 gshare history length |
| Bimodal | **yes** (LOGB=12) | **yes** (BIMODAL_SIZE=4096) | **no** | Separate bimodal base predictor |
| `BIMODAL_SIZE` / `LOGB` | **2^LOGB (2^12)** | **4096** | — | Bimodal table entries |
| P1 RAM count | *LINEINST* | *FETCH_WIDTH* | *N* | Number of parallel P1 RAMs |

## 4. Meta Predictor

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| `USE_META` | **#define** | **true** | **true** | Enable meta (alt-vs-primary) predictor |
| `METABITS` | *=4* | **4** | **4** | Meta counter width (signed) |
| `METAPIPE` | *=2* | **2** | **2** | Meta pipeline depth (array of counters) |
| `META_TABLE_SIZE` | — | **0** | — | Meta table size (0 = single global counter) |

## 5. Decay / U-bit Reset

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| `DECAY_CTR` | — | **8** | **0** | Probabilistic decay LFSR width (0=epoch reset) |
| `DECAY_GRAN` | — | **2** | **2** | Decay timing granularity (0=every misp, >0=1-in-2^N cycles) |
| `DecayPolicy` | — | **DecayMild** | **TDDecayMild** | Pluggable decay functor |
| `EPOCH_CTR_BITS` / `UCTRBITS` | *=8* | *=8* | **8** | Epoch counter width (saturate → reset u-bits) |
| `RESET_UBITS` | **#define** | *always* | *always* | Enable epoch-based u-bit reset |

### Available Decay Policies
- **Tage<>**: `DecayConservative`, `DecayMild`, `DecayAggressive`, `DecayHybrid`, `DecayMax`
- **TageDirect<>**: `TDDecayMild`, `TDDecayAggressive`

## 6. Allocation Config

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| `AllocCfg` | *hardcoded* | **DefaultAllocConfig** | **TDDefaultAllocConfig** | Allocation policy struct |
| `MAX_ALLOC` | *=1* | **1** | **1** | Max tables allocated per misprediction |
| `NON_CONSECUTIVE` | *=false* | **false** | **false** | Allow non-consecutive table allocation |
| `CONF_GATE` | — | **false** | **false** | Gate allocation on u-bit confidence |
| `PROB_START` | — | **0** | **0** | Probabilistic allocation start (0=disabled) |
| `PARTIAL_UPDATE` | — | **false** | **false** | Partial update policy |
| `MISPREDICT_ONLY_WRITE` | — | **false** | **false** | Only write tables on misprediction |
| `ALLOC_TRIGGER` | — | — | **MISPREDICT** | What triggers allocation attempt |
| `ALLOC_ACTION` | — | — | **STANDARD** | How allocation is gated |
| `UCTR_POLICY` | — | — | **ALWAYS_INC** | What drives the u-bit reset counter |
| `DISAGREE_EXTRA_CYCLE` | — | — | **true** | Claim extra cycle on P1/P2 disagree |
| `TARGET_POLICY` | — | — | **ClosestTarget** | Allocation target selection functor |
| `ALLOC_PRESSURE_BITS` | — | — | **0** | Allocation pressure register width (0=off) |
| `ACCURACY_PRESSURE_BITS` | — | — | **0** | Accuracy pressure register width (0=off) |

### Available AllocCfg Presets
- **Tage<>**: `DefaultAllocConfig`, `Alloc2Config`, `Alloc2NCConfig`, `AllocConfGateConfig`, `AllocProbStartConfig`, `PartialUpdateConfig`, `AllocFullConfig`
- **TageDirect<>**: `TDDefaultAllocConfig`, `TDMispredOnlyAllocConfig`

### Available Target Policies (TageDirect<> only)
- `ClosestTarget` — allocate in closest table above provider with u=0
- `DeterministicSkipTarget<SKIP>` — skip N closest candidates
- `StaticSkipTarget<SKIP, PROB>` — probabilistic skip
- `AllocPressureSkipTarget<SKIP>` — skip gated by alloc pressure register
- `AccuracyPressureSkipTarget<SKIP>` — skip gated by accuracy pressure

### Available AllocTrigger Values (TageDirect<> only)
- `MISPREDICT` — final P2 misprediction
- `BASE_WRONG` — P1 gshare was wrong
- `DISAGREE` — P1 and P2 disagree
- `TAGE_MISS` — no TAGE table match
- `TAGE_WRONG` — TAGE provider wrong (even if meta/P1 corrected)
- `ALWAYS` — every update cycle

### Available UctrPolicy Values (TageDirect<> only)
- `FARALLOC` — increment when alloc is far from provider
- `PENALTY_NA` — Seznec-style: inc on alloc failure, dec on success (every cycle)
- `NOALLOC` — inc on failure, dec on success (mispredict only)
- `ALWAYS_INC` — always increment on misprediction

## 7. Table Storage / Banking

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| Table impl | flat arrays | `TageTable` tuple | `TDTable` tuple | Per-table storage type |
| Tag RAM type | `ram<>` | `ram<>` | `ram<>` | Tag storage |
| Pred RAM type | `ram<>` | `ram<>` | `ram<>` | Prediction counter storage |
| Hyst RAM type | `rwram<>` (4 banks) | `rwram<>` (4 banks) | `td_rwram<>` | Hysteresis storage |
| U-bit RAM type | `rwram<>` (4 banks) | `rwram<>` or `ff<>` | `td_rwram<>` | U-bit storage |
| `SHARED_TAG` | — | **true** | — | Share tag RAM across branch-slot banks |
| `SHARED_U` | — | **true** | — | Share u-bits across branch-slot banks |
| `SHARED_HYS` | **false** | **true** | **false** | Share hysteresis (halve hyst table) |
| `U_STOR_FF` | — | **false** | — | U-bits in flip-flops vs SRAM |
| `USE_FF_CACHE` | — | **false** | — | Cache SRAM reads in FFs for block reuse |
| `RWRAM_BANKS` | *=4* | *=4 (in TageTable)* | **4** | RWRAM internal banking factor |
| `RWRAM_BANK_SHIFT` | — | — | **0** | RWRAM bank index shift |
| `ResetFn` | — | **DefaultResetFn** | — | U-bit reset functor (FF mode) |

## 8. History Folding

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| Fold impl | `geometric_folds` | `geometric_folds` | `TDHistoryFolder` | History folding mechanism |
| `FoldFn` | — | — | **XORFold** | Pluggable fold functor |
| `USE_DIR_HIST` | — | — | **false** | Concat direction bit with path bits in fold |

### Available Fold Functors (TageDirect<> only)
- `XORFold<F>` — standard XOR folding
- `RotateXORFold<F>` — rotate + XOR folding

## 9. Path History

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| `USE_PATH_HIST` | — | **false** | **false** | Enable separate path history register |
| `PATH_HIST_WIDTH` | — | **27** | **27** | Path history register width |
| `PATH_BITS` / `PATHBITS` | *=6* | **6** | **6** | Bits of PC mixed per update |

## 10. Overrider

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| `Overrider` | — | **SCOverrider<3,8,6,FW>** | — | Pluggable overrider type |
| Overrider in P2 | — | **yes** (offset 0) | — | Override baked into predict2 return |

### Available Overriders (Tage<> only)
- `NoOverrider` — disabled
- `SCOverrider<NUM_BIAS, BIAS_LOG, SC_CTR_BITS, FW, NUM_GEHL, GEHL_LOG, ...>` — Statistical Corrector
- `LoopPredictorA<...>` — Loop predictor

## 11. Monitor

| Parameter | tage<> | Tage<> | TageDirect<> | Description |
|---|---|---|---|---|
| Monitor type | inline struct | `TageMonitor.hpp` | `TDMonitor.hpp` | Performance monitoring |
| Enabled by | `#ifdef TAGE_MONITOR` | `#ifdef TAGE_MONITOR` | `#ifdef TAGE_MONITOR` | Compile-time flag |

---

## Summary: Unique Parameters per Implementation

### Only in `tage<>` (reference)
- `LOGLB` (log2 line bytes — others use LINEINST/FETCH_WIDTH directly)
- `LOGB` (log2 bimodal size — others use BIMODAL_SIZE)
- `LOGG` (log2 table size, uniform — others use per-table arrays)

### Only in `Tage<>` (custom)
- `BR_P_ENTRY`, `NUM_BANKS` (branch-slot banking)
- `SHARED_TAG`, `SHARED_U` (per-bank sharing)
- `U_STOR_FF`, `USE_FF_CACHE`, `ResetFn` (FF u-bit mode)
- `BIMODAL_SIZE` (parameterized bimodal)
- `META_TABLE_SIZE` (per-PC meta)
- `Overrider` (SC/Loop pluggable)
- `ConstBitsSize` functor

### Only in `TageDirect<>` (rank-indexed)
- `N` (decoupled from LINEINST)
- `AllocTrigger`, `AllocAction`, `UctrPolicy` enums
- `DISAGREE_EXTRA_CYCLE`
- `TARGET_POLICY` (5 allocation target functors)
- `ALLOC_PRESSURE_BITS`, `ACCURACY_PRESSURE_BITS`
- `FoldFn` (pluggable history fold: XORFold, RotateXORFold)
- `USE_DIR_HIST`
- `RWRAM_BANK_SHIFT`
- `EPOCH_CTR_BITS` (parameterized, vs hardcoded 8)
- `InvGeoSize` functor
