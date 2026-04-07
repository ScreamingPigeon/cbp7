# TageDirect Parameter Sweep Guide

## Quick Start

```bash
# 1. Build monitor + profile
make clean-monitor clean-profile
make cbp-monitor
# 2. Run monitor on target trace
make cbp-monitor TRACE=traces/502-gcc-all_16112_trace.gz MONITOR_MEASURE=40000000
# 3. Run profile for latency/energy
make cbp-profile-analyze TRACE=traces/502-gcc-all_16112_trace.gz
# 4. Full 20-trace VFS
make quick-eval
```

## VFS Score

VFS = f(IPC, EPI, latency). Target: >0.886 (current best).

Constraints:
- **P1 latency** must ceil to 1 cycle
- **P2 latency** must ceil to 2 cycles
- **EPI** (energy per instruction) directly penalizes hardware cost
- Misprediction rate drives IPC

## What to Look For in Monitor Output

| Metric | Healthy | Problem |
|--------|---------|---------|
| TAGE provider % | >10% | <5% means TAGE tables aren't matching |
| T0-T3 accuracy | >P1 accuracy | Below P1 = net negative, those tables hurt |
| Tag match rate | >2% per table | <1% = too much aliasing or stale entries |
| Allocation success | >95% | Low = u-bits blocking eviction |
| Epoch resets | >0 (if DECAY_CTR=0) | 0 = u-bits never clear, tables stale |
| Decay fires | >0 (if DECAY_CTR>0) | 0 = prob decay not working |
| P1 accuracy | ~82-85% | Baseline; TAGE should beat this |
| Extra cycle % | <30% | High = too many writes, hurts IPC |
| Meta chose alt % | <50% | 100% = meta always overrides provider |

## Parameters to Vary

All parameters are in the `TageDirect` alias in `TageDirect.hpp`. Lines marked `// NOTE: @modify`.

### Tier 1: Highest Impact

| Parameter | Default | Range | Effect |
|-----------|---------|-------|--------|
| `BASE_SIZE` (TDTableConfig arg 2) | 2048 | 512-8192 | Table entry count. Bigger = fewer aliases, more EPI |
| `MAXH` (TDTableConfig arg 8) | 400 | 100-1000 | Max history length. Longer = better long-range correlation, slower warmup |
| `MINH` (TDTableConfig arg 7) | 16 | 2-32 | Min history. Higher = less overlap with P1 gshare |
| `NTABLES` (TDTableConfig arg 1) | 8 | 4-12 | More tables = finer history granularity, more area |
| `SIZE_RATIO` (TDTableConfig arg 9) | 4 | 1-8 | 1=uniform size, >1=geometric scaling |
| `SizeFn` (TDTableConfig arg 12) | `GeoSize` | `UniformSize`, `InvGeoSize`, `StepSize`, `SqrtHistSize` | How table sizes scale across history lengths |
| `P1_TABLE_SIZE` | 4096 | 1024-16384 | P1 gshare RAM size |
| `P1_HIST` | 6 | 4-12 | P1 gshare history bits |

### Tier 2: Moderate Impact

| Parameter | Default | Range | Effect |
|-----------|---------|-------|--------|
| `TAG` (TDTableConfig arg 3) | 11 | 8-14 | Tag width. Wider = fewer false matches, more storage |
| `TagFn` (TDTableConfig arg 11) | `UniformTag<11>` | `GradedTag`, `StepTag`, `LogTag` | Per-table tag width |
| `DECAY_CTR` | 0 | 0-12 | 0=epoch reset, >0=probabilistic decay width |
| `DECAY_GRAN` | 2 | 0-4 | How often decay threshold updates (0=every mispredict) |
| `DECAY_POLICY` | `TDDecayMild` | `TDDecayConservative`, `TDDecayAggressive`, `TDDecayHybrid`, `TDDecayMax` | Decay aggressiveness (prob decay only) |
| `UCTRBITS` | 7 | 5-10 | Epoch counter width. Smaller = more frequent resets |
| `METABITS` | 4 | 3-6 | Meta predictor counter width |
| `METAPIPE` | 2 | 1-3 | Meta pipeline depth |
| `AllocCfg` | `TDDefaultAllocConfig` | See below | Allocation policy |

### Tier 3: Fine Tuning

| Parameter | Default | Range | Effect |
|-----------|---------|-------|--------|
| `LINEINST` | 1024 | 64-2048 | Max instructions per line. Affects block boundaries |
| `N` | 7 | 4-8 | Max branches per block. LANES = next_pow2(N) |
| `SHARED_TAG/U/HYS` | true | true/false | Share tag/u/hyst storage. true saves area |
| `U_STOR_FF` | false | true/false | Store u-bits in flip-flops (fast but expensive) |
| `CTR` (TDTableConfig arg 4) | 1 | 1-3 | Counter width per entry |
| `HYST` (TDTableConfig arg 5) | 2 | 1-3 | Hysteresis width |
| `U` (TDTableConfig arg 6) | 1 | 1-2 | Usefulness bit width |
| `RWRAM_BANKS` | 4 | 2-8 | rwram bank count. More = fewer conflicts |
| `RWRAM_BANK_SHIFT` | 0 | 0-4 | Which address bits select bank |
| `PATH_HIST` | false | true/false | Enable path history XOR into index |

### Allocation Configs

| Config | Description |
|--------|-------------|
| `TDDefaultAllocConfig` | MAX_ALLOC=1, standard |
| `Alloc2Config` | MAX_ALLOC=2 |
| `Alloc2NCConfig` | MAX_ALLOC=2, non-consecutive |
| `AllocConfGateConfig` | Only allocate when provider confidence is low |
| `AllocProbStartConfig` | Randomize allocation start table |
| `PartialUpdateConfig` | Partial counter updates |
| `AllocFullConfig` | All options enabled |

## Iterative Process

### Step 1: Baseline

Run monitor + profile on gcc trace. Record misprediction rate, provider distribution, P1 accuracy, EPI, latency.

### Step 2: Diagnose

- **TAGE tables rarely provide (<5%)**: Increase table sizes, reduce TAG width, check epoch resets
- **Short tables (T0-T3) worse than P1**: Increase MINH to reduce overlap, or use InvGeoSize to give them more entries
- **Long tables (T6-T7) dominate allocation**: Expected. If accuracy is low, reduce MAXH
- **Epoch resets = 0**: DECAY_CTR>0 uses prob decay (check decay_fire). Set DECAY_CTR=0 for epoch reset
- **High extra cycle %**: Reduce writes (MISPREDICT_ONLY_WRITE), reduce tables
- **High EPI**: Reduce table sizes, enable SHARED_*, reduce TAG/CTR/HYST widths

### Step 3: Sweep One Parameter

Change ONE parameter. Rebuild monitor, re-run same trace. Compare metrics.

```bash
# Example: try InvGeoSize (bigger short-history tables)
# Edit TageDirect.hpp SizeFn, then:
make clean-monitor && make cbp-monitor TRACE=traces/502-gcc-all_16112_trace.gz MONITOR_MEASURE=40000000
```

### Step 4: Validate on 20 Traces

Once single-trace metrics improve, run full evaluation:

```bash
# Edit params.yaml: predictor_type: "TageDirect<>"
make predictor-config && make quick-eval
```

### Step 5: Compare

```bash
make quick-eval-all PRED_A='TageDirect<>' PRED_B='Tage<>'
```

## Example Sweep Sequence

1. Baseline: `TageDirect<>` default params → monitor + quick-eval
2. Table sizing: try `InvGeoSize<2048,4>` vs `UniformSize<2048>` vs `GeoSize<2048,4>`
3. History range: sweep MINH in {4, 8, 16, 32}, MAXH in {200, 400, 800}
4. Epoch tuning: UCTRBITS in {5, 6, 7, 8} (resets every 32/64/128/256 far-allocs)
5. P1 sizing: P1_TABLE_SIZE in {2048, 4096, 8192}, P1_HIST in {4, 6, 8}
6. Table count: NTABLES in {6, 8, 10}
7. Tag width: TAG in {9, 11, 13} or GradedTag
8. Fine tune: RWRAM_BANKS, BANK_SHIFT, allocation config
