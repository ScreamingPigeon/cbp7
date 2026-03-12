#!/usr/bin/env python3
"""
Gaussian parameter sweep for Tage predictor on representative traces.

Generates configs by perturbing params around a baseline, builds, runs on
the quick_eval representative trace subset, and reports VFS scores.

Key feature: enforces geometric table sizing (short history = more entries).
Instead of SweepTableConfig's SIZE_RATIO, we directly generate per-table
TABLE_SIZE arrays using a GeometricTableConfig.
"""

import argparse
import csv
import hashlib
import math
import os
import random
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field, fields, asdict
from pathlib import Path
from typing import Optional


# ============================================================================
# Representative traces (same as scripts/quick_eval.sh)
# ============================================================================

REPR_TRACES = [
    "502-gcc-all_16112_trace.gz",
    "505-mcf-1_14364_trace.gz",
    "531-deepsjeng-1_53703_trace.gz",
    "548-exchange2-1_207156_trace.gz",
    "508-namd-1_126793_trace.gz",
    "554-roms-1_62613_trace.gz",
    "llvm-2.25945_0_trace.gz",
    "gcc-1.52513_0_trace.gz",
    "java16-specjbb-64k-ir1000.19011_0_trace.gz",
    "dcapo-kafka-jdk8.524_0_trace.gz",
    "web_74_trace.gz",
    "web_130_trace.gz",
    "nodejs-octane_3483_trace.gz",
    "nodejs-http2_4818_trace.gz",
    "rsbench-1.730_0_trace.gz",
    "sampleflow-1.127768_0_trace.gz",
    "zstd-1.19139_0_trace.gz",
    "gap-sssp-usroad.287_0_trace.gz",
    "python3-pyperf-dulwich.2000_0_trace.gz",
    "lua-3.25585_0_trace.gz",
]


# ============================================================================
# Geometric table sizing
# ============================================================================

def geometric_table_sizes(num_tables: int, base_size: int, size_ratio: float) -> list[int]:
    """Generate per-table sizes: short history (high index) gets more entries.

    Index 0 = longest history = smallest table.
    Index N-1 = shortest history = largest table.

    size_ratio: ratio of largest to smallest table.
        1.0 = uniform, 2.0 = 2x range, 4.0 = 4x range.
    """
    if num_tables <= 1 or size_ratio <= 1.0:
        return [base_size] * num_tables

    sizes = []
    for i in range(num_tables):
        # t goes 0..1 from longest to shortest history
        t = i / max(1, num_tables - 1)
        # Scale from 1/sqrt(ratio) to sqrt(ratio) centered on base_size
        scale = size_ratio ** (t - 0.5)
        sz = int(base_size * scale)
        # Round up to nearest power of 2, minimum 64
        result = 64
        while result < sz:
            result *= 2
        sizes.append(result)
    return sizes


def geometric_hist_lengths(num_tables: int, minhist: int, maxhist: int) -> list[int]:
    """Generate geometric history lengths. Index 0 = longest, index N-1 = shortest."""
    if num_tables <= 1:
        return [maxhist]
    lengths = []
    prev = 0
    ratio = maxhist / minhist
    for i in range(num_tables):
        e = i / (num_tables - 1)
        hl = int(minhist * (ratio ** e))
        hl = max(hl, prev + 1)
        lengths.append(hl)
        prev = hl
    return list(reversed(lengths))  # index 0 = longest


# ============================================================================
# Config dataclass
# ============================================================================

@dataclass(frozen=True)
class SweepConfig:
    # Table config
    num_tables: int = 8
    base_size: int = 2048       # base table size (center of geometric range)
    size_ratio: float = 1.0     # ratio of largest to smallest table (1=uniform)
    tag_width: int = 11
    ctr_width: int = 1
    hyst_width: int = 2
    u_width: int = 1
    minhist: int = 2
    maxhist: int = 100
    # Tage-level
    bimodal_size: int = 4096
    decay_ctr: int = 1024
    # P1
    p1_table_size: int = 16384
    p1_hist: int = 6
    # Meta
    metabits: int = 4
    metapipe: int = 2

    @property
    def config_id(self) -> str:
        return hashlib.md5(self.predictor_string.encode()).hexdigest()[:8]

    @property
    def table_sizes(self) -> list[int]:
        return geometric_table_sizes(self.num_tables, self.base_size, self.size_ratio)

    @property
    def hist_lengths(self) -> list[int]:
        return geometric_hist_lengths(self.num_tables, self.minhist, self.maxhist)

    @property
    def predictor_string(self) -> str:
        """Build a predictor string using a custom GeometricTableConfig."""
        sizes = self.table_sizes
        hists = self.hist_lengths

        # We encode per-table sizes into a custom config struct.
        # Use direct Tage template with a custom table config.
        # Problem: SweepTableConfig only supports uniform sizes with SIZE_RATIO.
        # Solution: use a lambda-generated config or encode sizes directly.
        #
        # Actually, the simplest approach: use SweepTableConfig with SIZE_RATIO=1
        # when uniform, and for non-uniform, we need per-table sizes.
        #
        # The C++ code already has `generate_table_sizes` with `size_fn`.
        # But SIZE_RATIO in SweepTableConfig maps to the same geometric scaling.
        # Let's verify the sizes match, and if so, just use SIZE_RATIO.
        #
        # For the sweep, we'll quantize size_ratio to integer values that
        # SweepTableConfig accepts.

        sr_int = max(1, round(self.size_ratio))

        tc = (f"SweepTableConfig<{self.num_tables},{self.base_size},"
              f"{self.tag_width},{self.ctr_width},{self.hyst_width},"
              f"{self.u_width},{self.minhist},{self.maxhist},{sr_int}>")

        parts = [tc, "DefaultAllocConfig"]
        parts.append("16")                      # FETCH_WIDTH
        parts.append(str(self.bimodal_size))
        parts.append("1")                       # BR_P_ENTRY
        parts.append("1")                       # NUM_BANKS
        parts.append("false")                   # USE_AHEAD
        parts.append("true")                    # SHARED_TAG
        parts.append("true")                    # SHARED_U
        parts.append("true")                    # SHARED_HYS
        parts.append("false")                   # U_STOR_FF
        parts.append(str(self.decay_ctr))
        parts.append("DefaultResetFn")
        parts.append("false")                   # USE_FF_CACHE
        parts.append("true")                    # P1_USE_GSHARE
        parts.append(str(self.p1_table_size))
        parts.append(str(self.p1_hist))
        parts.append("true")                    # USE_META
        parts.append(str(self.metabits))
        parts.append(str(self.metapipe))
        parts.append("0")                       # META_TABLE_SIZE
        parts.append("false")                   # USE_PATH_HIST
        parts.append("27")                      # PATH_HIST_WIDTH
        parts.append("6")                       # PATH_BITS

        return f"Tage<{','.join(parts)}>"


DEFAULTS = SweepConfig()


# ============================================================================
# Gaussian config generation
# ============================================================================

def generate_gaussian(n_samples: int, seed: int = 42) -> list[SweepConfig]:
    rng = random.Random(seed)

    # (field_name, choices, sigma_in_index_space)
    param_specs = [
        ("num_tables",    [4, 6, 8, 10, 12],                         1.0),
        ("base_size",     [512, 1024, 2048, 4096],                    1.0),
        ("size_ratio",    [1.0, 1.5, 2.0, 3.0, 4.0],                 1.0),
        ("tag_width",     [7, 8, 9, 10, 11, 12, 13],                 1.2),
        ("ctr_width",     [1, 2, 3],                                  0.5),
        ("hyst_width",    [1, 2, 3],                                  0.5),
        ("u_width",       [1, 2],                                     0.3),
        ("minhist",       [2, 3, 4, 6, 8],                            0.8),
        ("maxhist",       [50, 75, 100, 130, 160, 200, 300],          1.2),
        ("bimodal_size",  [2048, 4096, 8192, 16384],                  0.8),
        ("decay_ctr",     [256, 512, 1024, 2048, 4096],               0.8),
        ("p1_table_size", [4096, 8192, 16384, 32768, 65536],          0.8),
        ("p1_hist",       [4, 5, 6, 7, 8, 10],                        0.8),
        ("metabits",      [2, 3, 4, 6, 8],                            0.6),
        ("metapipe",      [1, 2, 3],                                  0.4),
    ]

    # Find default index for each param
    defaults_dict = {f.name: getattr(DEFAULTS, f.name) for f in fields(SweepConfig)}
    default_indices = {}
    for name, choices, _ in param_specs:
        dv = defaults_dict[name]
        # Find closest index
        best_idx = 0
        best_dist = abs(choices[0] - dv) if isinstance(dv, (int, float)) else 0
        for j, c in enumerate(choices):
            d = abs(c - dv)
            if d < best_dist:
                best_dist = d
                best_idx = j
        default_indices[name] = best_idx

    configs = []
    seen = set()
    attempts = 0

    while len(configs) < n_samples and attempts < n_samples * 30:
        attempts += 1
        overrides = {}

        for name, choices, sigma in param_specs:
            di = default_indices[name]
            sampled_idx = rng.gauss(di, sigma)
            idx = max(0, min(len(choices) - 1, round(sampled_idx)))
            val = choices[idx]
            default_val = defaults_dict[name]
            if val != default_val:
                overrides[name] = val

        # Validate constraints
        minh = overrides.get("minhist", DEFAULTS.minhist)
        maxh = overrides.get("maxhist", DEFAULTS.maxhist)
        if minh >= maxh:
            continue

        if not overrides:
            continue

        try:
            cfg = SweepConfig(**overrides)
        except TypeError:
            continue

        if cfg.config_id in seen:
            continue
        seen.add(cfg.config_id)

        # Print the geometric table sizes for verification
        configs.append(cfg)

    return configs


# ============================================================================
# Build & Run
# ============================================================================

CXX = os.environ.get("CXX", "g++")
COMMON_FLAGS = "-std=c++20 -O3"
WARN_FLAGS = ("-Wall -Wextra -pedantic -Wold-style-cast "
              "-Wno-deprecated-declarations -Wno-mismatched-tags")


def build_config(cfg: SweepConfig, build_dir: Path, extra_flags: str = "") -> tuple:
    binary = build_dir / f"cbp_{cfg.config_id}"
    if binary.exists():
        return (cfg, binary, 0.0, None)

    build_dir.mkdir(parents=True, exist_ok=True)
    pred = cfg.predictor_string
    cmd = (f'{CXX} {COMMON_FLAGS} {WARN_FLAGS} {extra_flags} '
           f'-Itrace_files -o {binary} cbp.cpp -lz '
           f"-DPREDICTOR='{pred}'")

    t0 = time.monotonic()
    try:
        subprocess.run(cmd, shell=True, check=True,
                       capture_output=True, text=True, timeout=300)
        return (cfg, binary, time.monotonic() - t0, None)
    except subprocess.CalledProcessError as e:
        return (cfg, binary, time.monotonic() - t0, e.stderr[-500:])
    except subprocess.TimeoutExpired:
        return (cfg, binary, 300.0, "TIMEOUT")


def run_trace(binary: Path, trace: Path, warmup: int, measure: int) -> tuple:
    trace_name = trace.stem.replace("_trace", "")
    cmd = [str(binary), str(trace), trace_name, str(warmup), str(measure)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=600, check=True)
    except subprocess.CalledProcessError as e:
        return (None, f"Run failed: {e.stderr[:200]}")
    except subprocess.TimeoutExpired:
        return (None, "TIMEOUT")

    # Parse CSV output (single line)
    line = result.stdout.strip()
    if not line:
        return (None, "Empty output")

    try:
        parts = line.split(',')
        data = {
            "trace": parts[0],
            "instructions": int(parts[1]),
            "branches": int(parts[2]),
            "cond_branches": int(parts[3]),
            "predictions": int(parts[4]),
            "extra_cycles": int(parts[5]),
            "disagreements": int(parts[6]),
            "disagree_end": int(parts[7]),
            "mispredictions": int(parts[8]),
            "p1_latency": float(parts[9]),
            "p2_latency": float(parts[10]),
            "epi": float(parts[11]),
        }
        return (data, None)
    except (IndexError, ValueError) as e:
        return (None, f"Parse error: {e}")


# ============================================================================
# VFS computation
# ============================================================================

MISP_PENALTY = 8


def compute_metrics(runs: list[dict]) -> dict:
    """Compute aggregate metrics across multiple traces, matching predictor_metrics.py."""
    # First pass: find max P1/P2 latency across all traces
    p1_latency = 0
    p2_latency = 0
    for r in runs:
        p1_latency = max(p1_latency, math.ceil(r["p1_latency"]))
        p2_latency = max(p2_latency, math.ceil(r["p2_latency"]))

    # Second pass: compute per-trace IPC, CPI, EPI
    count = 0
    sum_inv_ipc = 0.0
    sum_cpi = 0.0
    sum_epi = 0.0
    sum_mpi = 0.0
    sum_mpki = 0.0

    for r in runs:
        instr = r["instructions"]
        preds = r["predictions"]
        extra = r["extra_cycles"]
        diverge = r["disagreements"]
        divend = r["disagree_end"]
        misp = r["mispredictions"]
        epi = r["epi"]

        mpi = misp / instr
        mpki = misp / instr * 1000

        if p2_latency <= p1_latency:
            cycles = preds * max(1, p2_latency)
        else:
            cycles = (preds * max(1, p1_latency) +
                      diverge * p2_latency -
                      divend * max(1, p1_latency))
        cycles += extra

        ipc = instr / cycles if cycles > 0 else 0.001
        cpi = mpi * (MISP_PENALTY + p2_latency)

        count += 1
        sum_inv_ipc += 1.0 / ipc
        sum_cpi += cpi
        sum_epi += epi
        sum_mpi += mpi
        sum_mpki += mpki

    avg_ipc = count / sum_inv_ipc if sum_inv_ipc > 0 else 0
    avg_cpi = sum_cpi / count
    avg_epi = sum_epi / count
    avg_mpi = sum_mpi / count
    avg_mpki = sum_mpki / count

    # VFS score
    IPCcbp0 = 8
    CPIcbp0 = 0.0315
    EPIcbp0 = 1000
    ALPHA = 1.625
    BETA = 4 * ALPHA / (ALPHA - 1) ** 2
    GAMMA = 2 / (ALPHA - 1)
    cbp_energy_ratio = 0.05
    EPI0 = EPIcbp0 / cbp_energy_ratio

    WPI0 = IPCcbp0 * CPIcbp0
    WPI = avg_ipc * avg_cpi
    speedup = (avg_ipc / IPCcbp0) * (1 + WPI0) / (1 + WPI)
    LAMBDA = 1 / (1 + WPI0 / 2) - cbp_energy_ratio
    normalizedEPI = ((avg_epi / EPIcbp0) * cbp_energy_ratio +
                     LAMBDA * speedup ** GAMMA) * (1 + WPI / 2)
    vfs = speedup * ALPHA * (1 - 2 / (1 + math.sqrt(1 + BETA / (speedup * normalizedEPI))))

    return {
        "ipc": avg_ipc, "cpi": avg_cpi, "epi": avg_epi,
        "mpi": avg_mpi, "mpki": avg_mpki, "vfs": vfs,
        "p1_lat": p1_latency, "p2_lat": p2_latency,
        "n_traces": count,
    }


# ============================================================================
# CSV I/O
# ============================================================================

CONFIG_FIELDS = [f.name for f in fields(SweepConfig)]
AGG_FIELDS = ["config_id", "predictor_string", "table_sizes",
              "n_traces", "mpki", "ipc", "epi", "p1_lat", "p2_lat", "vfs"]
ALL_FIELDS = AGG_FIELDS + CONFIG_FIELDS


def load_done(path: Path) -> set[str]:
    done = set()
    if not path.exists():
        return done
    with open(path) as f:
        for row in csv.DictReader(f):
            done.add(row.get("config_id", ""))
    return done


def append_row(path: Path, row: dict):
    need_header = not path.exists() or path.stat().st_size == 0
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ALL_FIELDS, extrasaction="ignore")
        if need_header:
            writer.writeheader()
        writer.writerow(row)


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Gaussian Tage sweep on representative traces")
    parser.add_argument("-n", "--num-configs", type=int, default=60,
                        help="Number of Gaussian-sampled configs")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--trace-dir", type=Path, default=Path("./traces"))
    parser.add_argument("--warmup", type=int, default=1000000)
    parser.add_argument("--measure", type=int, default=40000000)
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 4)
    parser.add_argument("--build-dir", type=Path, default=Path("out/gsweep/bin"))
    parser.add_argument("-o", "--output", type=Path, default=Path("out/gsweep/results.csv"))
    parser.add_argument("--extra-flags", type=str, default="")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--report", type=Path, default=None,
                        help="Print report from existing CSV")
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()

    if args.report:
        print_report(args.report, args.top)
        return

    # Validate traces
    trace_paths = []
    for t in REPR_TRACES:
        p = args.trace_dir / t
        if p.exists():
            trace_paths.append(p)
        else:
            print(f"WARNING: missing trace: {p}", file=sys.stderr)
    if not trace_paths:
        print("ERROR: no traces found", file=sys.stderr)
        sys.exit(1)
    print(f"Using {len(trace_paths)}/{len(REPR_TRACES)} representative traces", file=sys.stderr)

    # Generate configs (baseline + gaussian)
    configs = [DEFAULTS] + generate_gaussian(args.num_configs, args.seed)
    print(f"Generated {len(configs)} configs (1 baseline + {len(configs)-1} gaussian)",
          file=sys.stderr)

    # Show a sample with geometric table sizes
    print(f"\nBaseline table sizes: {DEFAULTS.table_sizes}", file=sys.stderr)
    for cfg in configs[1:4]:
        print(f"  {cfg.config_id}: tables={cfg.num_tables}, base_size={cfg.base_size}, "
              f"ratio={cfg.size_ratio}, sizes={cfg.table_sizes}, "
              f"hists={cfg.hist_lengths}", file=sys.stderr)
    print("", file=sys.stderr)

    # Resume support
    done = load_done(args.output) if args.resume else set()
    if done:
        print(f"Resuming: {len(done)} configs already complete", file=sys.stderr)
        configs = [c for c in configs if c.config_id not in done]
        if not configs:
            print("All configs already complete. Use --report to view results.", file=sys.stderr)
            return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.build_dir.mkdir(parents=True, exist_ok=True)

    # Phase 1: Build
    print(f"=== Building {len(configs)} configs (jobs={args.jobs}) ===", file=sys.stderr)
    binaries = {}
    build_errors = 0

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(build_config, cfg, args.build_dir, args.extra_flags): cfg
                   for cfg in configs}
        for i, future in enumerate(as_completed(futures), 1):
            cfg, binary, elapsed, err = future.result()
            if err:
                print(f"  [{i}/{len(configs)}] FAIL {cfg.config_id}: {err[:80]}", file=sys.stderr)
                build_errors += 1
            else:
                binaries[cfg.config_id] = binary
                status = f"{elapsed:.1f}s" if elapsed > 0 else "cached"
                print(f"  [{i}/{len(configs)}] OK   {cfg.config_id} ({status})", file=sys.stderr)

    print(f"\nBuilt {len(binaries)}/{len(configs)} ({build_errors} errors)\n", file=sys.stderr)

    # Phase 2: Run all traces per config, then aggregate
    csv_lock = threading.Lock()
    results_list = []

    def run_config(cfg: SweepConfig):
        binary = binaries.get(cfg.config_id)
        if not binary:
            return None

        runs = []
        for trace in trace_paths:
            data, err = run_trace(binary, trace, args.warmup, args.measure)
            if err:
                print(f"    FAIL {cfg.config_id} {trace.name}: {err[:60]}", file=sys.stderr)
                continue
            runs.append(data)

        if not runs:
            return None

        metrics = compute_metrics(runs)

        row = {"config_id": cfg.config_id,
               "predictor_string": cfg.predictor_string,
               "table_sizes": str(cfg.table_sizes)}
        row.update({f.name: getattr(cfg, f.name) for f in fields(SweepConfig)})
        row.update({k: f"{v:.6f}" if isinstance(v, float) else v
                    for k, v in metrics.items()})

        with csv_lock:
            append_row(args.output, row)

        return (cfg, metrics)

    runnable = [c for c in configs if c.config_id in binaries]
    total = len(runnable)
    print(f"=== Running {total} configs × {len(trace_paths)} traces (jobs={args.jobs}) ===",
          file=sys.stderr)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_config, cfg): cfg for cfg in runnable}
        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            if result:
                cfg, metrics = result
                print(f"  [{i}/{total}] {cfg.config_id} "
                      f"MPKI={metrics['mpki']:.3f} EPI={metrics['epi']:.0f} "
                      f"P2={metrics['p2_lat']} VFS={metrics['vfs']:.4f} "
                      f"sizes={cfg.table_sizes}",
                      file=sys.stderr)
            else:
                print(f"  [{i}/{total}] SKIP (build or run failure)", file=sys.stderr)

    print(f"\nDone. Results: {args.output}", file=sys.stderr)
    print_report(args.output, args.top)


def print_report(csv_path: Path, top_n: int = 20):
    if not csv_path.exists():
        print(f"No results file: {csv_path}", file=sys.stderr)
        return

    with open(csv_path) as f:
        rows = list(csv.DictReader(f))

    if not rows:
        print("No results.")
        return

    # Sort by VFS descending
    for r in rows:
        try:
            r["_vfs"] = float(r["vfs"])
        except (ValueError, KeyError):
            r["_vfs"] = 0.0
    rows.sort(key=lambda r: r["_vfs"], reverse=True)

    # Find baseline
    baseline_vfs = None
    for r in rows:
        if r["config_id"] == DEFAULTS.config_id:
            baseline_vfs = r["_vfs"]
            break

    print(f"\n{'Rank':>4}  {'ID':>8}  {'Tables':>6}  {'BaseSize':>8}  {'Ratio':>5}  "
          f"{'Tag':>3}  {'MinH':>4}  {'MaxH':>4}  "
          f"{'MPKI':>7}  {'EPI':>7}  {'P2lat':>5}  {'VFS':>8}  {'dVFS':>7}  "
          f"{'Sizes':>30}")
    print("-" * 140)

    for i, r in enumerate(rows[:top_n]):
        dvfs = f"{r['_vfs'] - baseline_vfs:+.4f}" if baseline_vfs else "N/A"
        marker = " *" if r["config_id"] == DEFAULTS.config_id else ""
        print(f"{i+1:4d}  {r['config_id']:>8}  "
              f"{r.get('num_tables','?'):>6}  {r.get('base_size','?'):>8}  "
              f"{r.get('size_ratio','?'):>5}  "
              f"{r.get('tag_width','?'):>3}  "
              f"{r.get('minhist','?'):>4}  {r.get('maxhist','?'):>4}  "
              f"{float(r.get('mpki',0)):7.3f}  {float(r.get('epi',0)):7.0f}  "
              f"{r.get('p2_lat','?'):>5}  {r['_vfs']:8.4f}  {dvfs:>7}{marker}  "
              f"{r.get('table_sizes','?'):>30}")

    if baseline_vfs:
        print(f"\n  * = baseline (Tage<> defaults, config_id={DEFAULTS.config_id})")
    print(f"  Total configs: {len(rows)}, showing top {min(top_n, len(rows))}")


if __name__ == "__main__":
    main()
