#!/usr/bin/env python3
"""Energy breakdown: logic vs fanout vs wiring vs RAM for one or two predictors.

Runs each predictor 4 times with different HARCOM flag combos, averages EPI
across traces, then decomposes:

  RAM+Logic = EPI(FREE_FANOUT + FREE_WIRING)
  Fanout    = EPI(FREE_WIRING) - RAM+Logic
  Wiring    = EPI(FREE_FANOUT) - RAM+Logic
  Other     = Baseline - RAM+Logic - Fanout - Wiring

Usage:
  python3 scripts/energy_breakdown.py <pred1> [pred2] [options]

Options:
  --traces DIR      Trace directory (default: traces)
  --jobs N          Parallel jobs (default: 8)
  --warmup N        Warmup instructions (default: 1000)
  --measure N       Measure instructions (default: 40000)
  --build-dir DIR   Build directory (default: build)
"""

import sys
import os
import subprocess
import shutil
import tempfile
import concurrent.futures
from pathlib import Path

# ── Column indices (cbp.hpp output) ─────────────────────────────────────────
COL_NINSTR = 1
COL_EPI    = 11

# ── Flag combos ─────────────────────────────────────────────────────────────
FLAG_COMBOS = [
    ("baseline",    ""),
    ("free_wiring", "-DFREE_WIRING"),
    ("free_fanout", "-DFREE_FANOUT"),
    ("free_both",   "-DFREE_WIRING -DFREE_FANOUT"),
]

def build(predictor, extra_flags, build_dir, jobs=8):
    """Build binary for given predictor + flags. Returns binary path."""
    cmd = [
        "make", "cbp",
        f"PREDICTOR={predictor}",
        f"EXTRA_COMMON_FLAGS={extra_flags}",
        f"BUILD_DIR={build_dir}",
        f"-j{jobs}",
        "--no-print-directory",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"BUILD FAILED: {predictor} {extra_flags!r}", file=sys.stderr)
        print(result.stderr[-2000:], file=sys.stderr)
        return None
    # Find the built binary (most recently modified cbp-* in build_dir)
    bins = sorted(Path(build_dir).glob("cbp-????????"),
                  key=lambda p: p.stat().st_mtime, reverse=True)
    return str(bins[0]) if bins else None

def run_trace(binary, trace_path, warmup, measure):
    """Run binary on one trace. Returns (ninstr, epi) or None."""
    result = subprocess.run(
        [binary, str(trace_path), "x", str(warmup), str(measure)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return None
    line = result.stdout.strip().split("\n")[0]
    parts = line.split(",")
    if len(parts) <= COL_EPI:
        return None
    try:
        return float(parts[COL_NINSTR]), float(parts[COL_EPI])
    except ValueError:
        return None

def avg_epi(binary, traces, warmup, measure, jobs):
    """Run binary on all traces in parallel. Returns instruction-weighted avg EPI."""
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(run_trace, binary, t, warmup, measure): t for t in traces}
        for fut in concurrent.futures.as_completed(futs):
            r = fut.result()
            if r:
                results.append(r)
    if not results:
        return None
    total_instr = sum(n for n, _ in results)
    total_energy = sum(n * e for n, e in results)
    return total_energy / total_instr

def breakdown(epi):
    """Compute component breakdown from 4 EPI values dict."""
    base      = epi["baseline"]
    fw        = epi["free_wiring"]
    ff        = epi["free_fanout"]
    fb        = epi["free_both"]
    ram_logic = fb
    fanout    = fw - fb
    wiring    = ff - fb
    other     = base - ram_logic - fanout - wiring
    return {
        "Baseline": base,
        "RAM+Logic": ram_logic,
        "Fanout":    fanout,
        "Wiring":    wiring,
        "Other":     other,
    }

def print_breakdown(name, bd):
    base = bd["Baseline"]
    print(f"\n{'─'*50}")
    print(f"  {name}")
    print(f"{'─'*50}")
    for k, v in bd.items():
        pct = v / base * 100 if base else 0
        bar = "█" * int(pct / 2)
        print(f"  {k:<12} {v:7.1f} fJ  ({pct:5.1f}%)  {bar}")

def print_comparison(name1, bd1, name2, bd2):
    keys = list(bd1.keys())
    col_w = max(len(name1), len(name2), 20)
    print(f"\n{'═'*70}")
    print(f"  {'Component':<12}  {name1:>{col_w}}  {name2:>{col_w}}  {'Delta':>10}  {'Delta%':>7}")
    print(f"{'─'*70}")
    for k in keys:
        v1, v2 = bd1[k], bd2[k]
        delta = v2 - v1
        dpct  = (delta / v1 * 100) if v1 else 0
        sign  = "+" if delta >= 0 else ""
        print(f"  {k:<12}  {v1:>{col_w}.1f}  {v2:>{col_w}.1f}  {sign}{delta:>9.1f}  {sign}{dpct:>6.1f}%")
    print(f"{'═'*70}")

# ── Main ─────────────────────────────────────────────────────────────────────
def parse_args():
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pred1")
    p.add_argument("pred2", nargs="?")
    p.add_argument("--traces",    default="traces")
    p.add_argument("--jobs",      type=int, default=8)
    p.add_argument("--warmup",    type=int, default=1000)
    p.add_argument("--measure",   type=int, default=40000)
    p.add_argument("--build-dir", default="build")
    return p.parse_args()

def run_for(predictor, args):
    print(f"\n[{predictor}] Building 4 flag combos...", flush=True)
    traces = sorted(Path(args.traces).glob("*_trace.gz"))
    if not traces:
        print(f"No traces found in {args.traces}", file=sys.stderr)
        sys.exit(1)

    epi_map = {}
    for combo_name, flags in FLAG_COMBOS:
        print(f"  [{predictor}] {combo_name:<14} building...", end=" ", flush=True)
        binary = build(predictor, flags, args.build_dir, args.jobs)
        if binary is None:
            sys.exit(1)
        print(f"running {len(traces)} traces...", end=" ", flush=True)
        epi = avg_epi(binary, traces, args.warmup, args.measure, args.jobs)
        if epi is None:
            print(f"ERROR: no results", file=sys.stderr)
            sys.exit(1)
        epi_map[combo_name] = epi
        print(f"EPI={epi:.1f} fJ")

    return breakdown(epi_map)

def main():
    args = parse_args()
    bd1 = run_for(args.pred1, args)
    print_breakdown(args.pred1, bd1)

    if args.pred2:
        bd2 = run_for(args.pred2, args)
        print_breakdown(args.pred2, bd2)
        print_comparison(args.pred1, bd1, args.pred2, bd2)

if __name__ == "__main__":
    main()
