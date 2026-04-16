# TageAhead Parameter Reference

All template parameters for `TageAhead<>` and `TATableConfig<>`, their defaults, and intended effects.

## TATableConfig

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N` | 8 | Number of TAGE tables |
| `SIZE` | 1024 | Base table size (entries) |
| `TAG` | 10 | Base tag width (bits) |
| `MINH` | 10 | Minimum history length |
| `MAXH` | 100 | Maximum history length |
| `SIZE_RATIO` | 1 | Size scaling ratio between tables |
| `HIST` | GEOMETRIC | History length series: GEOMETRIC, QUADRATIC, SUPEREXP, ROS |
| `TagFn` | GradedTag<TAG,TAG-1> | Tag width functor (graded: shorter tags for longer-history tables) |
| `SizeFn` | GeoSize<SIZE,SIZE_RATIO> | Table size functor |

## TageAhead — Core

| Parameter | Default | Description |
|-----------|---------|-------------|
| `TableCfg` | TATableConfig<> | Table configuration (sizes, tags, history lengths) |
| `N` | 8 | Max conditional branches per block |
| `PATHBITS` | 6 | Bits of next_pc injected into global history |
| `SEC_TAG_BITS` | 3 | Secondary tag width for ahead-pipeline ambiguity resolution |
| `USE_SEC_TAG` | true | Enable secondary tag matching |
| `CTR_WIDTH` | 1 | Prediction counter width per lane. 1 = direction only; >1 = sign-magnitude with hyst forming the magnitude. Paper recommendation: 3 |
| `HYST_WIDTH` | 3 | Hysteresis (confidence) counter width, separate from pred. Combined with CTR_WIDTH forms the effective saturating counter |
| `U_WIDTH` | 2 | Usefulness counter width. 1 = binary useful/not. 2 = enables gradual u accumulation (needed for INC_DEC, DIP) |
| `FB_CAPACITY` | 16384 | Fallback table size (bimodal or gshare) |
| `USE_GSHARE` | false | Use gshare (PC^history) vs bimodal (PC-only) as fallback |
| `GS_HIST` | 6 | Gshare history length (only when USE_GSHARE=true) |
| `META_WIDTH` | 5 | Meta counter width (decides provider vs alt on weak predictions) |
| `META_CAPACITY` | 256 | Meta table entries (PC-indexed RAM, not global register) |
| `META_PIPE` | 2 | Meta pipeline depth |
| `LINEINST` | 1024 | Line size in instructions |
| `SHARED_HYS` | true | Share one hyst counter between 2 adjacent entries |
| `HIST_MODE` | PATH | What goes into history: PATH (pc bits), DIR (branch direction), BOTH |

## Allocation Policy (AllocCfg)

Controlled via `AllocCfg` struct. Default: `TADefaultAllocConfig`.

| Field | Default | Description |
|-------|---------|-------------|
| `MAX_ALLOC` | 1 | Max tables to allocate into per misprediction |
| `NON_CONSECUTIVE` | false | Allow non-consecutive table allocation |
| `ALLOC_TRIGGER` | MISPREDICT | When to attempt allocation |
| `ALLOC_ACTION` | STANDARD | Allocation action type |
| `TARGET_POLICY` | ClosestTarget | Which table(s) to target. Options: ClosestTarget, StaticSkipTarget<skip,prob>, DeterministicSkipTarget<skip> |

Predefined configs: `TADefaultAllocConfig`, `TAAllocSkip1`, `TAAllocDetSkip1`, `TAAlloc2`

## Pressure Counters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ACC_WIDTH` | 4 | Accuracy counter width. Tracks overall prediction accuracy. Feeds into epoch trigger and decay threshold |
| `ALLOC_WIDTH` | 4 | Allocation pressure counter width. Tracks how often allocation fails (no u=0 slot). Feeds into epoch trigger |

## Probabilistic U-bit Decay

Per-miss probabilistic u-bit decrement. Alternative to epoch-based bulk reset.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DECAY_ENABLE` | false | Enable probabilistic decay |
| `DECAY_MISS` | TAG_OR_SEC | What constitutes a "miss" for decay: TAG, SEC, TAG_OR_SEC, TAG_AND_SEC |
| `DECAY_OP` | DECREMENT | Decay operation: DECREMENT (u-1), HALVE (u>>1), CLEAR (u=0) |
| `DECAY_LFSR_WIDTHS` | uniform(8) | Per-table LFSR widths for random threshold |
| `DecayThreshFn` | DefaultDecayThresh | Functor computing decay probability threshold from pressure counters |

## Epoch-based U-bit Reset

Bulk reset of all u-bits when trigger fires. Opens allocation slots by clearing u-bit protection.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `EPOCH_ENABLE` | true | Enable epoch-based u-bit reset |
| `EpochTriggerFn` | DefaultEpochTrigger | Functor deciding when to fire epoch reset (see trigger options below) |
| `EPOCH_CTR_WIDTH` | 16 | Free-running counter width for interval-based triggers |
| `EPOCH_RESET_ACC` | false | Reset acc_ctr to 0 when epoch fires |
| `EPOCH_RESET_ALLOC` | true | Reset alloc_ctr to 0 when epoch fires |

### Epoch Trigger Functors

All have signature: `should_fire<AW, PW, EW>(acc_ctr, alloc_ctr, epoch_ctr) → val<1>`

| Functor | Fires when | Notes |
|---------|------------|-------|
| `AllocSaturateEpoch` (default) | `alloc_ctr == max` | Original behavior. Very aggressive with low alloc success rate — fires ~102K times on gcc trace, u never accumulates |
| `FixedIntervalEpoch<PERIOD>` | `epoch_ctr % PERIOD == 0` | PERIOD must be power of 2. Ignores pressure entirely. Predictable cadence |
| `AllocAccJointEpoch<ALLOC_T, ACC_T>` | `alloc_ctr >= ALLOC_T AND acc_ctr <= ACC_T` | Only resets when pressure is high AND accuracy is poor. Avoids resetting during accurate-but-pressured periods |
| `CountdownEpoch<PERIOD, ACC_GATE>` | `epoch_ctr % PERIOD == 0 AND acc_ctr < ACC_GATE` | Regular interval but skips reset when predictor is doing well. PERIOD must be power of 2 |

## Replacement Policy (Technique 4)

Controls which entries are considered replaceable during allocation.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ReplacePolicyFn` | ReplaceUZero | Functor: `is_replaceable(u, hyst, alloc_p, acc_p) → bool`. Options: ReplaceUZero (u==0), ReplaceUZeroWeakConf<THRESH> (u==0 AND hyst<THRESH), ReplacePressureAdaptive<THRESH,SHIFT> |

**ReplaceUZero**: Only replace entries with u=0. Current default, simple.

**ReplaceUZeroWeakConf<T>**: Also requires hyst < T. Protects high-confidence entries even when u=0.

**ReplacePressureAdaptive<T,S>**: Loosens replacement criteria when allocation pressure is high.

## Alt Bank Promotion (Technique 3)

When provider is wrong and alt is correct, optionally set u-bit on alt's entry to protect it.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `AltPromoteFn` | AltPromoteOff | Functor: `should_promote(prov_wrong, alt_correct, alloc_p, acc_p, rng) → bool`. Options: AltPromoteOff, AltPromoteAlways, AltPromoteProb<PROB_256>, AltPromotePressure |

## DIP-like Allocation (Technique 6)

Give newly allocated entries initial u>0 so they survive longer before eviction.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DIP_PROB_256` | 0 | Probability (0-256) of setting u>0 on allocation. 0=off, 256=always |
| `DIP_INIT_U` | 1 | Initial u value when DIP fires. With U_WIDTH=2, value of 1 gives one extra epoch of protection |

## Provider U-bit Update Policy

How the provider entry's u-bit is updated after a prediction.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `U_PROV_UPDATE` | SET_OR_CLEAR | Enum controlling provider u-bit update behavior |

| Mode | Correct | Wrong | Notes |
|------|---------|-------|-------|
| `SET_OR_CLEAR` | u = max | u = 0 | Current default. Aggressive — one wrong clears all protection |
| `SET_ON_CORRECT` | u = max | no write | Tage.hpp style. Epochs handle clearing. More protective |
| `INC_DEC` | u = sat(u+1) | u = sat(u-1) | Gradual. Needs multiple wrongs to become evictable. Best with U_WIDTH>=2 |
| `INC_ONLY` | u = sat(u+1) | no write | Most protective. u only grows; only epochs/decay can shrink it |

## Fallback Trains Toward P2

When USE_GSHARE is enabled, the gshare fallback can train toward the TAGE provider's prediction (like Tage.hpp's P1→P2 training). When TAGE is provider and gshare disagrees, gshare's prediction is overwritten with TAGE's prediction.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `FB_TRAIN_P2` | false | Enable gshare-toward-TAGE training. Requires USE_GSHARE=true |

## Far-Allocation Epoch Pressure

Tage.hpp-style alloc pressure: instead of incrementing alloc_ctr when no slot is found, increment when allocation lands far from the provider (indicating under-provisioning). Only updates on misprediction.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `FARALLOC_DIST` | 0 | 0=off (use default alloc pressure). When >0, alloc_ctr increments on far alloc (≥DIST tables from provider), decrements on close alloc. Tage.hpp uses 3 |

## Per-branch Hyst/U Banking

Decouple hyst and u tracking per branch within a block to avoid training contamination.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `HYST_BANKS` | 1 | Hyst banks per table. 1=shared (all branches share one hyst per entry), N=per-branch |
| `HYST_BANK_BIT` | 0 | Starting bit of branch rank for hyst bank selection |
| `U_BANKS` | 1 | U banks per table. Same as HYST_BANKS but for usefulness counters |
| `U_BANK_BIT` | 0 | Starting bit of branch rank for u bank selection |

**Experimental finding**: Banking showed negligible improvement (~2K mispredictions with 4 banks). With avg 1.6 branches/block, contamination isn't the bottleneck. The real issue is allocation/capacity.

## Convenience Typedefs

| Name | Key differences from default |
|------|------------------------------|
| `TageAhead_AltPromAlways` | AltPromoteAlways |
| `TageAhead_AltPromProb128` | AltPromoteProb<128> |
| `TageAhead_AltPromPress` | AltPromotePressure |
| `TageAhead_ReplWeakConf1` | ReplaceUZeroWeakConf<1> |
| `TageAhead_ReplPressAdapt` | ReplacePressureAdaptive<1,4> |
| `TageAhead_DIP64` | DIP_PROB_256=64 |
| `TageAhead_DIP128` | DIP_PROB_256=128 |
| `TageAhead_AllTechniques` | ReplaceUZeroWeakConf<1> + AltPromoteAlways + DIP_PROB_256=64 |
| `TageAhead_HystBank2` | HYST_BANKS=2 (bit 1) |
| `TageAhead_Bank4` | HYST_BANKS=4, U_BANKS=4 (bit 0) |
| `TageAhead_Bank8` | HYST_BANKS=8, U_BANKS=8 (fully per-branch) |
