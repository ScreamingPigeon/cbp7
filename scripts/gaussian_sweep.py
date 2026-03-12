#!/usr/bin/env python3
"""
Gaussian parameter sweep for Tage predictor on representative traces.

Generates configs by perturbing params around a baseline, builds, runs on
the quick_eval representative trace subset, and reports VFS scores.

Key feature: enforces geometric table sizing (short history = more entries).
"""

import argparse
import os
import random
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import fields
from pathlib import Path

from sweep_common import (
    SweepConfig, DEFAULTS, CONFIG_FIELDS,
    build_config, run_trace, compute_metrics,
    resolve_traces, append_csv_row, load_csv_rows,
)


# ============================================================================
# Gaussian config generation
# ============================================================================

def generate_gaussian(n_samples: int, seed: int = 42) -> list[SweepConfig]:
    rng = random.Random(seed)

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

    defaults_dict = {f.name: getattr(DEFAULTS, f.name) for f in fields(SweepConfig)}
    default_indices = {}
    for name, choices, _ in param_specs:
        dv = defaults_dict[name]
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
        configs.append(cfg)

    return configs


# ============================================================================
# CSV I/O
# ============================================================================

AGG_FIELDS = ["config_id", "predictor_string", "table_sizes",
              "n_traces", "mpki", "ipc", "epi", "p1_lat", "p2_lat", "vfs"]
ALL_FIELDS = AGG_FIELDS + CONFIG_FIELDS


def load_done(path: Path) -> set[str]:
    done = set()
    rows = load_csv_rows(path)
    for row in rows:
        done.add(row.get("config_id", ""))
    return done


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

    trace_paths = resolve_traces(args.trace_dir)
    if not trace_paths:
        print("ERROR: no traces found", file=sys.stderr)
        sys.exit(1)
    print(f"Using {len(trace_paths)} representative traces", file=sys.stderr)

    configs = [DEFAULTS] + generate_gaussian(args.num_configs, args.seed)
    print(f"Generated {len(configs)} configs (1 baseline + {len(configs)-1} gaussian)",
          file=sys.stderr)

    print(f"\nBaseline table sizes: {DEFAULTS.table_sizes}", file=sys.stderr)
    for cfg in configs[1:4]:
        print(f"  {cfg.config_id}: tables={cfg.num_tables}, base_size={cfg.base_size}, "
              f"ratio={cfg.size_ratio}, sizes={cfg.table_sizes}, "
              f"hists={cfg.hist_lengths}", file=sys.stderr)
    print("", file=sys.stderr)

    done = load_done(args.output) if args.resume else set()
    if done:
        print(f"Resuming: {len(done)} configs already complete", file=sys.stderr)
        configs = [c for c in configs if c.config_id not in done]
        if not configs:
            print("All configs already complete. Use --report to view results.", file=sys.stderr)
            return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.build_dir.mkdir(parents=True, exist_ok=True)

    # Build
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

    # Run
    csv_lock = threading.Lock()

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
            append_csv_row(args.output, row, ALL_FIELDS)

        return (cfg, metrics)

    runnable = [c for c in configs if c.config_id in binaries]
    total = len(runnable)
    print(f"=== Running {total} configs x {len(trace_paths)} traces (jobs={args.jobs}) ===",
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
    rows = load_csv_rows(csv_path)
    if not rows:
        print("No results.")
        return

    for r in rows:
        try:
            r["_vfs"] = float(r["vfs"])
        except (ValueError, KeyError):
            r["_vfs"] = 0.0
    rows.sort(key=lambda r: r["_vfs"], reverse=True)

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
