# Tage Parameter Sweep

## Infrastructure

### SweepTableConfig (`Tage.hpp`)

Parameterized table config for sweep use. 9 template params:

```cpp
template <u64 N = 8, u64 SIZE = 2048, u64 TAG = 11,
          u64 CTR = 1, u64 HYST = 2, u64 U = 1,
          u64 MINH = 2, u64 MAXH = 100, u64 SIZE_RATIO = 1>
struct SweepTableConfig;
```

- Defaults match `DefaultTableConfig` exactly (verified by `static_assert`)
- `SIZE_RATIO > 1` enables geometric table size scaling (long-history tables smaller, short-history larger)
- Generates predictor strings like: `Tage<SweepTableConfig<8,1024,10,1,2,1,3,130,1>, DefaultAllocConfig, 16, 8192, ...>`

### Sweep Script (`scripts/tage_sweep.py`)

```bash
# Gaussian sweep: 100 configs around defaults, 8 parallel jobs
python3 scripts/tage_sweep.py \
  --gaussian 100 --seed 42 \
  --traces traces/502-gcc-all_16112_trace.gz \
  --warmup 1000000 --measure 40000000 \
  --jobs 8 --build-dir out/sweep/bin \
  --output out/sweep/results.csv --resume

# Report top 20
python3 scripts/tage_sweep.py --report out/sweep/results.csv --top 20
```

**Modes:**
- `--tier 0`: Baseline only
- `--tier 1`: One-at-a-time sweep (~50 configs)
- `--gaussian N`: Isotropic Gaussian sampling around defaults (perturbs multiple params simultaneously)
- `--report FILE`: Print leaderboard from existing CSV

**Features:**
- Parallel build AND parallel run via `ThreadPoolExecutor`
- `--resume` skips completed `(config_id, trace)` pairs
- VFS computation matches `compare_predictors.sh`
- Parses stdout (accuracy) + stderr (hw metrics from `-DVERBOSE`)

### Makefile Targets

```bash
make sweep                          # Run full sweep (tier 1, all traces)
make sweep SWEEP_TIER=0 SWEEP_MEASURE=10000  # Quick baseline check
make sweep-report                   # Print top-20 leaderboard
```

Variables: `SWEEP_TIER`, `SWEEP_JOBS`, `SWEEP_WARMUP`, `SWEEP_MEASURE`

## Results: Gaussian Sweep (2026-03-12)

**Trace**: gcc (502-gcc-all_16112), 1M warmup, 40M measure
**Baseline `Tage<>`**: VFS=0.7111, MPKI=9.033, EPI=1536, P2=1.86 cycles

101 configs generated, 79 built (22 compile failures), 64 ran (15 OOB runtime errors).

### Leaderboard (top 10)

| Rank | ID | Tables | Size | Tag | MPKI | EPI | P2lat | VFS | dVFS |
|------|----|--------|------|-----|------|-----|-------|-----|------|
| 1 | bfc1ff8e | 8 | 1024 | 10 | 9.10 | 1295 | 1.71 | 0.7181 | +0.0069 |
| 2 | f1b5229c | 8 | 2048 | 11 | 8.94 | 1455 | 1.87 | 0.7155 | +0.0043 |
| 3 | 0afd487f | 10 | 1024 | 10 | 9.01 | 1566 | 1.76 | 0.7148 | +0.0037 |
| 4 | 1e2f86f0 | 8 | 1024 | 10 | 9.20 | 1198 | 1.68 | 0.7146 | +0.0034 |
| 5 | 5e0263f6 | 6 | 2048 | 12 | 9.08 | 1360 | 1.82 | 0.7144 | +0.0032 |
| 6 | ebb47990 | 6 | 2048 | 11 | 9.00 | 1239 | 1.89 | 0.7133 | +0.0022 |
| 7 | f49e39cc | 6 | 2048 | 12 | 9.28 | 1219 | 1.84 | 0.7118 | +0.0007 |
| 8 | ff46e43d | 8 | 2048 | 11 | 8.86 | 1621 | 1.89 | 0.7113 | +0.0001 |
| 9* | 4849b6a4 | 8 | 2048 | 11 | 9.03 | 1536 | 1.86 | 0.7111 | baseline |
| 10 | 5f3787e8 | 6 | 1024 | 10 | 9.69 | 1046 | 1.64 | 0.7101 | -0.0011 |

### Winner: bfc1ff8e (VFS +0.98%)

Full params vs baseline:

| Param | Winner | Baseline | Change |
|-------|--------|----------|--------|
| table_size | 1024 | 2048 | Halved — big EPI savings |
| tag_width | 10 | 11 | Narrower tags |
| minhist | 3 | 2 | Slightly longer min |
| maxhist | 130 | 100 | Extended range |
| bimodal_size | 8192 | 4096 | 2x bigger bimodal |
| decay_ctr | 2048 | 1024 | Slower u-bit decay |
| p1_table_size | 32768 | 16384 | 2x bigger P1 gshare |
| p1_hist | 7 | 6 | Slightly longer P1 hist |

### Key Insights

1. **VFS rewards EPI heavily** — shrinking TAGE entries trades small accuracy for large EPI wins
2. **Invest in P1/bimodal** — bigger gshare + bimodal compensates for smaller TAGE tables
3. **maxhist 130** — extending history range helps slightly
4. **decay_ctr 2048** — default 1024 decays useful bits too fast
5. **Single-trace caveat** — needs multi-trace validation before adopting

### Known Issues

- **OOB errors**: Small bimodal/P1 tables cause index overflow (bimodal index not masked for non-power-of-2 or non-default sizes)
- **Compile failures (~22%)**: HARCOM `val<1>` instantiation errors with certain param combos (likely tag/counter width constraints)
- **SIZE_RATIO > 1**: Causes OOB at runtime (gindex wider than small table RAM)
