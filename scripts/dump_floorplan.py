#!/usr/bin/env python3
"""Dump HARCOM floorplan.gv per predictor for downstream wire-length analysis.

For each predictor passed on the command line, builds a binary, runs it briefly
on one trace (just enough for HARCOM to emit floorplan.gv), and copies the .gv
to a named output path. The floorplan is determined at construction time, so a
single quick run is sufficient — no need for full simulation.

Usage:
  python3 scripts/dump_floorplan.py <Predictor1> [Predictor2] ... [options]

Options:
  --out-dir DIR   Directory to write {Predictor}.gv files (default: out/floorplans)
  --trace PATH    Trace to run on (default: first trace under traces/)
  --warmup N      Warmup instructions (default: 1000)
  --measure N     Measure instructions (default: 1000)
  --build-dir DIR Build directory (default: build/floorplan_dump)
  --jobs N        Parallel build jobs (default: nproc)
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


import hashlib


# floorplan.gv is only emitted when harcom is built with -DFLOORPLAN
# (see harcom.hpp:3636 and cbp.hpp:165). We pass it via EXTRA_COMMON_FLAGS.
FLOORPLAN_FLAG = "-DFLOORPLAN"


def pred_hash(predictor, extra_common=FLOORPLAN_FLAG, extra_cbp=""):
    """Reproduce the Makefile's PRED_HASH:
        echo 'PREDICTOR_TYPE|EXTRA_COMMON_FLAGS|EXTRA_CBP_FLAGS' | md5sum | cut -c1-8
    Note: shell `echo` appends a trailing newline; we must include it."""
    key = f"{predictor}|{extra_common}|{extra_cbp}\n"
    return hashlib.md5(key.encode()).hexdigest()[:8]


def build(predictor, build_dir, jobs):
    """Build cbp binary for predictor; returns binary path (resolved via PRED_HASH)."""
    Path(build_dir).mkdir(parents=True, exist_ok=True)
    cmd = [
        "make", "cbp",
        f"PREDICTOR_TYPE={predictor}",
        f"EXTRA_COMMON_FLAGS={FLOORPLAN_FLAG}",
        f"BUILD_DIR={build_dir}",
        f"-j{jobs}",
        "--no-print-directory",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[{predictor}] BUILD FAILED", file=sys.stderr)
        print(res.stderr[-2000:], file=sys.stderr)
        return None
    # Resolve the exact binary via Makefile's hash scheme, not mtime
    h = pred_hash(predictor)
    binary = Path(build_dir) / f"cbp-{h}"
    if not binary.exists():
        print(f"[{predictor}] expected binary {binary} not found", file=sys.stderr)
        return None
    return str(binary)


def dump_one(predictor, build_dir, trace, warmup, measure, out_dir, jobs):
    print(f"[{predictor}] building...", flush=True)
    binary = build(predictor, build_dir, jobs)
    if binary is None:
        return False

    # Wipe any stale floorplan.gv so we fail loudly if this run doesn't produce one
    gv = Path("floorplan.gv")
    if gv.exists():
        gv.unlink()

    name = Path(trace).name.replace("_trace.gz", "")
    print(f"[{predictor}] running ({warmup} warmup, {measure} measure)...", flush=True)
    res = subprocess.run([binary, str(trace), name, str(warmup), str(measure)],
                         capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[{predictor}] run failed: {res.stderr[-500:]}", file=sys.stderr)
        return False

    if not gv.exists():
        print(f"[{predictor}] WARNING: floorplan.gv not produced "
              f"(was -DFLOORPLAN actually compiled in?)", file=sys.stderr)
        return False

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    dest = Path(out_dir) / f"{predictor}.gv"
    shutil.copy2(gv, dest)
    size_kb = dest.stat().st_size / 1024
    print(f"[{predictor}] → {dest}  ({size_kb:.1f} KiB)")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("predictors", nargs="+", help="predictor type names")
    ap.add_argument("--out-dir",   default="out/floorplans")
    ap.add_argument("--trace",     default=None,
                    help="trace path; default: first traces/*_trace.gz")
    ap.add_argument("--warmup",    type=int, default=1000)
    ap.add_argument("--measure",   type=int, default=1000)
    ap.add_argument("--build-dir", default="build/floorplan_dump")
    ap.add_argument("--jobs",      type=int, default=os.cpu_count() or 8)
    args = ap.parse_args()

    if args.trace is None:
        traces = sorted(Path("traces").glob("*_trace.gz"))
        if not traces:
            print("No traces found under traces/; pass --trace explicitly",
                  file=sys.stderr)
            sys.exit(1)
        args.trace = str(traces[0])

    ok = 0
    for pred in args.predictors:
        if dump_one(pred, args.build_dir, args.trace,
                    args.warmup, args.measure, args.out_dir, args.jobs):
            ok += 1
        print()

    print(f"Done: {ok}/{len(args.predictors)} floorplans → {args.out_dir}/")


if __name__ == "__main__":
    main()
