#!/usr/bin/env python3
"""
Coordinate descent optimizer for Tage branch predictor parameters.

Starts from a known config (defaults or best from prior sweep), then
iteratively improves one parameter at a time. Every evaluation is logged
to CSV for later visualization.

Usage:
    python3 scripts/coordinate_descent.py --trace-dir traces -j 8
    python3 scripts/coordinate_descent.py --start-from best --prior-csv out/gsweep/results.csv
    python3 scripts/coordinate_descent.py --report out/descent/results.csv
"""

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import fields
from datetime import datetime
from pathlib import Path

from sweep_common import (
    SweepConfig, DEFAULTS, CONFIG_FIELDS,
    build_config, run_trace, compute_metrics,
    resolve_traces, append_csv_row, load_csv_rows,
)


# ============================================================================
# Parameter domains and sweep order
# ============================================================================

PARAM_DOMAINS = {
    "num_tables":    [4, 6, 8, 10, 12],
    "base_size":     [256, 512, 1024, 2048, 4096, 8192],
    "size_ratio":    [1.0, 2.0, 4.0, 8.0, 16.0, 32.0],
    "tag_width":     [7, 8, 9, 10, 11, 12, 13, 15],
    "ctr_width":     [1, 2, 3],
    "hyst_width":    [1, 2, 3],
    "u_width":       [1, 2],
    "minhist":       [2, 3, 4, 6, 8],
    "maxhist":       [50, 75, 100, 130, 160, 200, 300, 400],
    "bimodal_size":  [1024, 2048, 4096, 8192, 16384],
    "decay_ctr":     [0, 6, 8, 10, 12],
    "decay_gran":    [0, 2, 4, 6],
    "p1_table_size": [4096, 8192, 16384, 32768, 65536],
    "p1_hist":       [4, 5, 6, 7, 8, 10],
    "metabits":      [2, 3, 4, 6, 8],
    "metapipe":      [1, 2, 3],
    # structural parameters
    "shared_tag":       [True, False],
    "shared_u":         [True, False],
    "shared_hys":       [True, False],
    "u_stor_ff":        [True, False],
    "use_path_hist":    [True, False],
    "path_hist_width":  [10, 16, 20, 27, 32],
    "path_bits":        [3, 4, 6, 8],
}

# Sweep order: most impactful parameters first (based on prior sweep insights),
# then never-evaluated params, then new structural params
PARAM_ORDER = [
    "maxhist", "minhist",
    # frozen for now: "num_tables", "tag_width",
    #                  "bimodal_size", "p1_table_size", "hyst_width"
    "decay_ctr", "decay_gran", "p1_hist",
    "base_size", "size_ratio", "ctr_width", "metabits", "metapipe", "u_width",
    # structural params
    "use_path_hist", "path_hist_width", "path_bits",
    "shared_tag", "shared_u", "shared_hys", "u_stor_ff",
]


# ============================================================================
# CSV schema
# ============================================================================

DESCENT_FIELDS = [
    "timestamp", "iteration", "step", "varied_param", "config_id",
    "is_improvement", "vfs", "mpki", "ipc", "epi", "p1_lat", "p2_lat",
    "n_traces", "build_time_s", "eval_time_s",
] + CONFIG_FIELDS + ["table_sizes", "predictor_string"]


# ============================================================================
# Evaluation cache
# ============================================================================

class EvalCache:
    """Cache of config_id -> metrics, backed by CSV."""

    def __init__(self, csv_path: Path):
        self.csv_path = csv_path
        self.cache: dict[str, dict] = {}
        self._load()

    def _load(self):
        rows = load_csv_rows(self.csv_path)
        for r in rows:
            cid = r.get("config_id", "")
            if cid:
                try:
                    self.cache[cid] = {
                        "vfs": float(r["vfs"]),
                        "mpki": float(r["mpki"]),
                        "ipc": float(r["ipc"]),
                        "epi": float(r["epi"]),
                        "p1_lat": r["p1_lat"],
                        "p2_lat": r["p2_lat"],
                        "n_traces": r["n_traces"],
                    }
                except (ValueError, KeyError):
                    pass

    def get(self, config_id: str) -> dict | None:
        return self.cache.get(config_id)

    def put(self, config_id: str, metrics: dict):
        self.cache[config_id] = metrics

    def __len__(self):
        return len(self.cache)


# ============================================================================
# Checkpoint
# ============================================================================

def save_checkpoint(path: Path, iteration: int, step: int,
                    config: SweepConfig, vfs: float, total_evals: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "iteration": iteration,
        "step": step,
        "current_config": config.to_dict(),
        "current_vfs": vfs,
        "total_evals": total_evals,
    }
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def load_checkpoint(path: Path) -> dict | None:
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


# ============================================================================
# Core: evaluate a single config
# ============================================================================

def evaluate_config(cfg: SweepConfig, build_dir: Path, trace_paths: list[Path],
                    warmup: int, measure: int, jobs: int,
                    extra_flags: str = "") -> tuple[dict | None, float, float]:
    """Build and run a config on all traces. Returns (metrics, build_time, eval_time)."""
    # Build
    _, binary, build_time, build_err = build_config(cfg, build_dir, extra_flags)
    if build_err:
        print(f"    BUILD FAIL {cfg.config_id}: {build_err[:80]}", file=sys.stderr)
        return None, build_time, 0.0

    # Run traces in parallel
    t0 = time.monotonic()
    runs = []

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(run_trace, binary, t, warmup, measure): t
                   for t in trace_paths}
        for future in as_completed(futures):
            trace = futures[future]
            data, err = future.result()
            if err:
                print(f"      TRACE FAIL {trace.name}: {err[:60]}", file=sys.stderr)
                continue
            runs.append(data)

    eval_time = time.monotonic() - t0

    if not runs:
        return None, build_time, eval_time

    metrics = compute_metrics(runs)
    return metrics, build_time, eval_time


# ============================================================================
# Coordinate descent
# ============================================================================

def coordinate_descent(args):
    # Resolve traces
    trace_paths = resolve_traces(args.trace_dir)
    if not trace_paths:
        print("ERROR: no traces found", file=sys.stderr)
        sys.exit(1)
    print(f"Using {len(trace_paths)} representative traces", file=sys.stderr)

    # Output paths
    csv_path = args.output
    checkpoint_path = csv_path.parent / "checkpoint.json"
    build_dir = csv_path.parent / "bin"
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    # Load cache from prior results
    cache = EvalCache(csv_path)
    # Also load from prior CSV if specified
    if args.prior_csv and args.prior_csv.exists():
        prior_cache = EvalCache(args.prior_csv)
        for cid, m in prior_cache.cache.items():
            if cid not in cache.cache:
                cache.cache[cid] = m
        print(f"Loaded {len(prior_cache)} entries from prior CSV", file=sys.stderr)

    # Determine starting config
    start_config = DEFAULTS
    start_vfs = None

    if args.start_from == "best" and args.prior_csv and args.prior_csv.exists():
        rows = load_csv_rows(args.prior_csv)
        best_row = None
        best_vfs = -1.0
        for r in rows:
            try:
                v = float(r["vfs"])
                if v > best_vfs:
                    best_vfs = v
                    best_row = r
            except (ValueError, KeyError):
                continue
        if best_row:
            overrides = {}
            for fname in CONFIG_FIELDS:
                if fname in best_row:
                    val = best_row[fname]
                    field_type = type(getattr(DEFAULTS, fname))
                    try:
                        overrides[fname] = field_type(val)
                    except ValueError:
                        pass
            start_config = SweepConfig(**overrides)
            start_vfs = best_vfs
            print(f"Starting from best prior config: {start_config.config_id} (VFS={best_vfs:.4f})",
                  file=sys.stderr)

    # Resume from checkpoint
    total_evals = 0
    start_iteration = 1
    start_step = 0

    if args.resume:
        ckpt = load_checkpoint(checkpoint_path)
        if ckpt:
            start_iteration = ckpt["iteration"]
            start_step = ckpt["step"] + 1  # resume from next step
            total_evals = ckpt["total_evals"]
            start_config = SweepConfig(**ckpt["current_config"])
            start_vfs = ckpt["current_vfs"]
            # If we finished all steps in an iteration, advance
            if start_step >= len(PARAM_ORDER):
                start_iteration += 1
                start_step = 0
            print(f"Resuming from iteration {start_iteration}, step {start_step} "
                  f"(config={start_config.config_id}, VFS={start_vfs:.4f}, "
                  f"evals={total_evals})", file=sys.stderr)

    # Evaluate starting config if needed
    current = start_config
    if start_vfs is not None:
        best_vfs = start_vfs
    else:
        cached = cache.get(current.config_id)
        if cached:
            best_vfs = cached["vfs"]
            print(f"Start config {current.config_id} cached: VFS={best_vfs:.4f}", file=sys.stderr)
        else:
            print(f"Evaluating start config {current.config_id}...", file=sys.stderr)
            metrics, bt, et = evaluate_config(
                current, build_dir, trace_paths, args.warmup, args.measure,
                args.jobs, args.extra_flags)
            if not metrics:
                print("ERROR: start config failed to evaluate", file=sys.stderr)
                sys.exit(1)
            best_vfs = metrics["vfs"]
            cache.put(current.config_id, metrics)
            total_evals += 1
            row = _make_row(0, 0, "init", current, metrics, True, bt, et)
            append_csv_row(csv_path, row, DESCENT_FIELDS)
            print(f"Start: VFS={best_vfs:.6f} MPKI={metrics['mpki']:.3f} "
                  f"EPI={metrics['epi']:.0f}", file=sys.stderr)

    print(f"\n{'='*70}", file=sys.stderr)
    print(f"Starting coordinate descent from VFS={best_vfs:.6f}", file=sys.stderr)
    print(f"Config: {current.config_id} {current.to_dict()}", file=sys.stderr)
    print(f"{'='*70}\n", file=sys.stderr)

    # Main descent loop
    for iteration in range(start_iteration, start_iteration + args.max_iterations):
        improved_this_iteration = False
        iter_start = time.monotonic()

        step_range_start = start_step if iteration == start_iteration else 0

        for step_idx in range(step_range_start, len(PARAM_ORDER)):
            param = PARAM_ORDER[step_idx]
            domain = PARAM_DOMAINS[param]
            current_val = getattr(current, param)

            print(f"\n--- Iteration {iteration}, Step {step_idx}: {param} "
                  f"(current={current_val}) ---", file=sys.stderr)

            # Generate candidates (skip current value)
            candidates = []
            for val in domain:
                if val == current_val:
                    continue
                candidate = current.replace(**{param: val})
                # Constraint: minhist < maxhist
                if candidate.minhist >= candidate.maxhist:
                    continue
                candidates.append((val, candidate))

            if not candidates:
                print(f"  No valid candidates for {param}", file=sys.stderr)
                save_checkpoint(checkpoint_path, iteration, step_idx,
                                current, best_vfs, total_evals)
                continue

            # Build all candidates in parallel
            print(f"  Building {len(candidates)} candidates...", file=sys.stderr)
            built = {}
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                futures = {pool.submit(build_config, cfg, build_dir, args.extra_flags): (val, cfg)
                           for val, cfg in candidates}
                for future in as_completed(futures):
                    val, cfg = futures[future]
                    _, binary, bt, err = future.result()
                    if err:
                        print(f"    BUILD FAIL {param}={val}: {err[:60]}", file=sys.stderr)
                    else:
                        built[val] = (cfg, binary, bt)

            # Evaluate candidates (check cache first, run rest in parallel)
            best_param_val = current_val
            best_param_vfs = best_vfs

            to_eval = []
            for val, (cfg, binary, bt) in built.items():
                cached = cache.get(cfg.config_id)
                if cached:
                    vfs = cached["vfs"]
                    print(f"    {param}={val} (cached) VFS={vfs:.6f} "
                          f"{'*** NEW BEST' if vfs > best_param_vfs else ''}", file=sys.stderr)
                    row = _make_row(iteration, step_idx, param, cfg, cached,
                                    vfs > best_param_vfs, bt, 0)
                    append_csv_row(csv_path, row, DESCENT_FIELDS)
                    if vfs > best_param_vfs:
                        best_param_vfs = vfs
                        best_param_val = val
                else:
                    to_eval.append((val, cfg, binary, bt))

            # Run uncached evaluations
            for val, cfg, binary, bt in to_eval:
                t0 = time.monotonic()
                runs = []
                with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                    futures = {pool.submit(run_trace, binary, t, args.warmup, args.measure): t
                               for t in trace_paths}
                    for future in as_completed(futures):
                        data, err = future.result()
                        if err:
                            continue
                        runs.append(data)

                et = time.monotonic() - t0
                total_evals += 1

                if not runs:
                    print(f"    {param}={val} EVAL FAIL (no trace results)", file=sys.stderr)
                    continue

                metrics = compute_metrics(runs)
                cache.put(cfg.config_id, metrics)
                vfs = metrics["vfs"]
                is_better = vfs > best_param_vfs

                row = _make_row(iteration, step_idx, param, cfg, metrics, is_better, bt, et)
                append_csv_row(csv_path, row, DESCENT_FIELDS)

                marker = "*** NEW BEST" if is_better else ""
                print(f"    {param}={val} VFS={vfs:.6f} MPKI={metrics['mpki']:.3f} "
                      f"EPI={metrics['epi']:.0f} {marker}", file=sys.stderr)

                if vfs > best_param_vfs:
                    best_param_vfs = vfs
                    best_param_val = val

            # Update current if improved
            if best_param_val != current_val:
                current = current.replace(**{param: best_param_val})
                best_vfs = best_param_vfs
                improved_this_iteration = True
                print(f"  >> Updated {param}: {current_val} -> {best_param_val} "
                      f"(VFS={best_vfs:.6f})", file=sys.stderr)
            else:
                print(f"  >> No improvement for {param} (staying at {current_val})",
                      file=sys.stderr)

            save_checkpoint(checkpoint_path, iteration, step_idx,
                            current, best_vfs, total_evals)

        iter_elapsed = time.monotonic() - iter_start
        print(f"\n{'='*70}", file=sys.stderr)
        print(f"Iteration {iteration} complete ({iter_elapsed:.0f}s, {total_evals} total evals)",
              file=sys.stderr)
        print(f"Current best: VFS={best_vfs:.6f} config={current.config_id}",
              file=sys.stderr)
        print(f"Config: {current.to_dict()}", file=sys.stderr)
        print(f"{'='*70}", file=sys.stderr)

        if not improved_this_iteration:
            print(f"\nConverged! No improvement in iteration {iteration}.", file=sys.stderr)
            break

    # Final summary
    print(f"\n{'='*70}", file=sys.stderr)
    print(f"FINAL RESULT: VFS={best_vfs:.6f}", file=sys.stderr)
    print(f"Config ID: {current.config_id}", file=sys.stderr)
    print(f"Parameters: {current.to_dict()}", file=sys.stderr)
    print(f"Table sizes: {current.table_sizes}", file=sys.stderr)
    print(f"Predictor: {current.predictor_string}", file=sys.stderr)
    print(f"Total evaluations: {total_evals}", file=sys.stderr)
    print(f"Results: {csv_path}", file=sys.stderr)
    print(f"{'='*70}", file=sys.stderr)


def _make_row(iteration, step, param, cfg, metrics, is_improvement, build_time, eval_time):
    row = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "iteration": iteration,
        "step": step,
        "varied_param": param,
        "config_id": cfg.config_id,
        "is_improvement": is_improvement,
        "build_time_s": f"{build_time:.1f}",
        "eval_time_s": f"{eval_time:.1f}",
        "table_sizes": str(cfg.table_sizes),
        "predictor_string": cfg.predictor_string,
    }
    # Metrics (could be dict from cache or from compute_metrics)
    for k in ["vfs", "mpki", "ipc", "epi"]:
        v = metrics.get(k, 0)
        row[k] = f"{float(v):.6f}" if isinstance(v, (int, float)) else v
    for k in ["p1_lat", "p2_lat", "n_traces"]:
        row[k] = metrics.get(k, "")
    # Config params
    row.update(cfg.to_dict())
    return row


# ============================================================================
# Report
# ============================================================================

def print_report(csv_path: Path, top_n: int = 30):
    rows = load_csv_rows(csv_path)
    if not rows:
        print("No results.", file=sys.stderr)
        return

    # Deduplicate by config_id, keeping best VFS
    by_config: dict[str, dict] = {}
    for r in rows:
        cid = r.get("config_id", "")
        try:
            vfs = float(r["vfs"])
        except (ValueError, KeyError):
            continue
        if cid not in by_config or vfs > float(by_config[cid]["vfs"]):
            by_config[cid] = r

    sorted_configs = sorted(by_config.values(), key=lambda r: float(r["vfs"]), reverse=True)

    # Find baseline
    baseline_vfs = None
    for r in sorted_configs:
        if r["config_id"] == DEFAULTS.config_id:
            baseline_vfs = float(r["vfs"])
            break

    # Print trajectory (all evaluations in order)
    print(f"\n=== Descent Trajectory ({len(rows)} evaluations) ===\n")
    print(f"{'#':>4}  {'Iter':>4}  {'Param':>14}  {'Value':>8}  "
          f"{'VFS':>8}  {'MPKI':>7}  {'EPI':>7}  {'Better':>6}  {'ID':>8}")
    print("-" * 95)

    for i, r in enumerate(rows[:80]):
        val = r.get(r.get("varied_param", ""), "?")
        marker = " *" if r.get("is_improvement", "").lower() == "true" else ""
        try:
            print(f"{i+1:4d}  {r.get('iteration','?'):>4}  {r.get('varied_param','?'):>14}  "
                  f"{val:>8}  {float(r.get('vfs',0)):8.6f}  "
                  f"{float(r.get('mpki',0)):7.3f}  {float(r.get('epi',0)):7.0f}  "
                  f"{marker:>6}  {r.get('config_id','?'):>8}")
        except ValueError:
            pass

    # Print leaderboard
    print(f"\n=== Leaderboard (top {min(top_n, len(sorted_configs))}) ===\n")
    print(f"{'Rank':>4}  {'ID':>8}  {'VFS':>8}  {'dVFS':>7}  {'MPKI':>7}  "
          f"{'EPI':>7}  {'Tables':>6}  {'Base':>6}  {'Ratio':>5}  "
          f"{'Tag':>3}  {'MinH':>4}  {'MaxH':>4}  {'Bimodal':>7}")
    print("-" * 110)

    for i, r in enumerate(sorted_configs[:top_n]):
        vfs = float(r["vfs"])
        dvfs = f"{vfs - baseline_vfs:+.4f}" if baseline_vfs else "N/A"
        marker = " *" if r["config_id"] == DEFAULTS.config_id else ""
        print(f"{i+1:4d}  {r['config_id']:>8}  {vfs:8.6f}  {dvfs:>7}  "
              f"{float(r.get('mpki',0)):7.3f}  {float(r.get('epi',0)):7.0f}  "
              f"{r.get('num_tables','?'):>6}  {r.get('base_size','?'):>6}  "
              f"{r.get('size_ratio','?'):>5}  {r.get('tag_width','?'):>3}  "
              f"{r.get('minhist','?'):>4}  {r.get('maxhist','?'):>4}  "
              f"{r.get('bimodal_size','?'):>7}{marker}")

    if baseline_vfs:
        print(f"\n  * = baseline defaults (config_id={DEFAULTS.config_id})")

    # Print best config details
    if sorted_configs:
        best = sorted_configs[0]
        print(f"\n=== Best Config ===")
        print(f"  Config ID: {best['config_id']}")
        print(f"  VFS: {float(best['vfs']):.6f}")
        for f in CONFIG_FIELDS:
            if f in best:
                print(f"  {f}: {best[f]}")
        print(f"  Table sizes: {best.get('table_sizes', '?')}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Coordinate descent optimizer for Tage predictor parameters")
    parser.add_argument("--start-from", choices=["defaults", "best"], default="defaults",
                        help="Start from defaults or best config from --prior-csv")
    parser.add_argument("--prior-csv", type=Path, default=None,
                        help="CSV from prior sweep to seed cache and find best start config")
    parser.add_argument("--trace-dir", type=Path, default=Path("./traces"))
    parser.add_argument("--warmup", type=int, default=1000000)
    parser.add_argument("--measure", type=int, default=40000000)
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 4)
    parser.add_argument("-o", "--output", type=Path, default=Path("out/descent/results.csv"))
    parser.add_argument("--extra-flags", type=str, default="")
    parser.add_argument("--max-iterations", type=int, default=5,
                        help="Max outer iterations before stopping")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from checkpoint")
    parser.add_argument("--report", type=Path, default=None,
                        help="Print report from existing CSV (no build/run)")
    parser.add_argument("--top", type=int, default=30)
    args = parser.parse_args()

    if args.report:
        print_report(args.report, args.top)
        return

    coordinate_descent(args)


if __name__ == "__main__":
    main()
