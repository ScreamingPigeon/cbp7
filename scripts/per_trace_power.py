#!/usr/bin/env python3
"""Per-trace energy breakdown: logic/RAM vs fanout vs wiring vs other.

Builds 4 variants of the predictor with HARCOM flag combos, runs each on the
trace set, then decomposes EPI per trace into:

  RAM+Logic = EPI(FREE_WIRING + FREE_FANOUT)
  Fanout    = EPI(FREE_WIRING)         - RAM+Logic
  Wiring    = EPI(FREE_FANOUT)         - RAM+Logic
  Other     = EPI(baseline)            - RAM+Logic - Fanout - Wiring

Also computes EPPC (Energy Per Prediction Cycle) = energy / (npred + extra),
which normalizes out MPKI-driven block-stretching effects on EPI.

Also dumps the HARCOM panel info (-DDEBUG_PRINT) and floorplan PDF once.

Usage:
  python3 scripts/per_trace_power.py <predictor> [--quick | --full] [options]

Options:
  --predictor PRED  Predictor type (default: TageAheadHC_IR_M2)
  --quick           Run on 20 representative traces (default)
  --full            Run on all *_trace.gz in trace dir
  --traces DIR      Trace directory (default: traces)
  --out CSV         Output CSV path (default: out/per_trace_power.csv)
  --jobs N          Parallel jobs (default: nproc)
  --warmup N        Warmup instructions (default: 1000000)
  --measure N       Measure instructions (default: 40000000)
  --build-dir DIR   Build directory (default: build/ptp)
  --panel-out PATH  Panel info dump path (default: out/panel_info.txt)
  --floorplan-out PATH  Floorplan PDF path (default: out/floorplan.pdf)
"""

import argparse
import concurrent.futures
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

# CSV column indices from cbp.hpp output
COL_NAME    = 0
COL_NINSTR  = 1
COL_NPRED   = 4
COL_EXTRA   = 5
COL_MISP    = 8
COL_EPI     = 11

FLAG_COMBOS = [
    ("baseline",    ""),
    ("free_wiring", "-DFREE_WIRING"),
    ("free_fanout", "-DFREE_FANOUT"),
    ("free_both",   "-DFREE_WIRING -DFREE_FANOUT"),
]

# 20 representative traces (mirrors scripts/quick_eval.sh)
QUICK_TRACES = [
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


def build(predictor, extra_flags, build_dir, tag, jobs=8, debug_print=False):
    """Build binary; returns path. Uses a tagged output name so 4 builds coexist."""
    Path(build_dir).mkdir(parents=True, exist_ok=True)
    out = Path(build_dir) / f"cbp-{tag}"
    flags = extra_flags + (" -DDEBUG_PRINT" if debug_print else "")
    cmd = [
        "make", "cbp",
        f"PREDICTOR_TYPE={predictor}",
        f"EXTRA_COMMON_FLAGS={flags}",
        f"BUILD_DIR={build_dir}",
        f"-j{jobs}",
        "--no-print-directory",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"BUILD FAILED ({tag}): {predictor} flags={flags!r}", file=sys.stderr)
        print(res.stderr[-2000:], file=sys.stderr)
        sys.exit(1)
    # most-recent cbp-* in build_dir
    bins = sorted(Path(build_dir).glob("cbp-????????"),
                  key=lambda p: p.stat().st_mtime, reverse=True)
    if not bins:
        print(f"No binary produced for {tag}", file=sys.stderr)
        sys.exit(1)
    src = bins[0]
    shutil.copy2(src, out)
    return str(out)


def run_one(binary, trace_path, warmup, measure):
    """Run binary on one trace. Returns CSV-split list or None on failure."""
    name = Path(trace_path).name.replace("_trace.gz", "")
    res = subprocess.run(
        [binary, str(trace_path), name, str(warmup), str(measure)],
        capture_output=True, text=True
    )
    if res.returncode != 0:
        return None
    line = res.stdout.strip().split("\n")[0] if res.stdout else ""
    if not line:
        return None
    parts = line.split(",")
    if len(parts) <= COL_EPI:
        return None
    return parts


def run_all(binary, traces, warmup, measure, jobs, label):
    """Run binary on all traces in parallel. Returns dict[name] -> parts."""
    print(f"  [{label}] running {len(traces)} traces ({jobs} jobs)...", flush=True)
    out = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(run_one, binary, t, warmup, measure): t for t in traces}
        done = 0
        for fut in concurrent.futures.as_completed(futs):
            parts = fut.result()
            done += 1
            if parts is None:
                continue
            out[parts[COL_NAME]] = parts
    print(f"  [{label}] done: {len(out)}/{len(traces)} succeeded", flush=True)
    return out


def dump_panel_and_floorplan(predictor, build_dir, trace, warmup, measure,
                             panel_out, floorplan_out, jobs):
    """Build with -DDEBUG_PRINT, run once, capture stderr (panel info), copy
    floorplan.gv → PDF."""
    print(f"\n[panel] building -DDEBUG_PRINT variant...", flush=True)
    binary = build(predictor, "", build_dir, "dbg", jobs=jobs, debug_print=True)
    name = Path(trace).name.replace("_trace.gz", "")
    print(f"[panel] running on {name} to capture panel info + floorplan...", flush=True)
    # floorplan.gv is written to the cwd by harcom; run from project root
    res = subprocess.run(
        [binary, str(trace), name, str(warmup), str(measure)],
        capture_output=True, text=True
    )
    Path(panel_out).parent.mkdir(parents=True, exist_ok=True)
    with open(panel_out, "w") as f:
        f.write(f"# Predictor: {predictor}\n")
        f.write(f"# Trace: {trace}\n")
        f.write(f"# Warmup: {warmup}  Measure: {measure}\n")
        f.write("# stderr capture (panel.print + DEBUG_PRINT signal traces)\n\n")
        f.write(res.stderr)
    print(f"[panel] panel info → {panel_out}", flush=True)

    # Convert floorplan.gv → PDF via graphviz if available
    gv = Path("floorplan.gv")
    if gv.exists():
        Path(floorplan_out).parent.mkdir(parents=True, exist_ok=True)
        dot = shutil.which("dot")
        if dot is None:
            print(f"[panel] WARNING: 'dot' not in PATH; copying floorplan.gv → "
                  f"{floorplan_out.replace('.pdf', '.gv')}", file=sys.stderr)
            shutil.copy2(gv, floorplan_out.replace(".pdf", ".gv"))
        else:
            r = subprocess.run([dot, "-Tpdf", "floorplan.gv", "-o", floorplan_out],
                               capture_output=True, text=True)
            if r.returncode != 0:
                print(f"[panel] dot failed: {r.stderr}", file=sys.stderr)
            else:
                print(f"[panel] floorplan → {floorplan_out}", flush=True)
    else:
        print(f"[panel] WARNING: floorplan.gv not generated", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--predictor", default="TageAheadHC_IR_M2")
    grp = p.add_mutually_exclusive_group()
    grp.add_argument("--quick", action="store_true", help="20 representative traces (default)")
    grp.add_argument("--full",  action="store_true", help="all *_trace.gz")
    grp.add_argument("--trace", help="single trace path or basename (e.g. 502-gcc-all_16112_trace.gz)")
    p.add_argument("--traces",        default="traces")
    p.add_argument("--out",           default="out/per_trace_power.csv")
    p.add_argument("--jobs",          type=int, default=os.cpu_count() or 8)
    p.add_argument("--warmup",        type=int, default=1_000_000)
    p.add_argument("--measure",       type=int, default=40_000_000)
    p.add_argument("--build-dir",     default="build/ptp")
    p.add_argument("--panel-out",     default="out/panel_info.txt")
    p.add_argument("--floorplan-out", default="out/floorplan.pdf")
    args = p.parse_args()

    # Resolve trace list
    if args.trace:
        tp = Path(args.trace)
        if not tp.exists():
            tp = Path(args.traces) / args.trace
        if not tp.exists():
            print(f"Trace not found: {args.trace}", file=sys.stderr)
            sys.exit(1)
        traces = [tp]
        mode = f"single ({tp.name})"
    elif args.full:
        traces = sorted(Path(args.traces).glob("*_trace.gz"))
        mode = "full"
    else:
        traces = []
        for name in QUICK_TRACES:
            tp = Path(args.traces) / name
            if tp.exists():
                traces.append(tp)
            else:
                print(f"WARNING: missing trace {tp}", file=sys.stderr)
        mode = "quick"
    if not traces:
        print(f"No traces in {args.traces}", file=sys.stderr)
        sys.exit(1)

    print(f"=== per_trace_power: {args.predictor} | {mode} ({len(traces)} traces) ===")

    # 1. Panel dump + floorplan (once)
    dump_panel_and_floorplan(args.predictor, args.build_dir, traces[0],
                             args.warmup, args.measure,
                             args.panel_out, args.floorplan_out, args.jobs)

    # 2. Build 4 variants
    print(f"\n=== building 4 flag combos ===")
    binaries = {}
    for tag, flags in FLAG_COMBOS:
        print(f"[build] {tag:<12} flags={flags!r}", flush=True)
        binaries[tag] = build(args.predictor, flags, args.build_dir, tag, args.jobs)

    # 3. Run all 4 on the trace set
    print(f"\n=== running {len(traces)} traces × 4 variants ===")
    results = {}
    for tag, _ in FLAG_COMBOS:
        results[tag] = run_all(binaries[tag], traces, args.warmup, args.measure,
                               args.jobs, tag)

    # 4. Per-trace decomposition
    names = sorted(set(results["baseline"].keys())
                   & set(results["free_wiring"].keys())
                   & set(results["free_fanout"].keys())
                   & set(results["free_both"].keys()))
    if not names:
        print("ERROR: no traces completed across all 4 variants", file=sys.stderr)
        sys.exit(1)

    rows = []
    for n in names:
        b   = results["baseline"][n]
        fw  = results["free_wiring"][n]
        ff  = results["free_fanout"][n]
        fb  = results["free_both"][n]

        ninstr = float(b[COL_NINSTR])
        npred  = float(b[COL_NPRED])
        extra  = float(b[COL_EXTRA])
        misp   = float(b[COL_MISP])
        epi_b  = float(b[COL_EPI])
        epi_fw = float(fw[COL_EPI])
        epi_ff = float(ff[COL_EPI])
        epi_fb = float(fb[COL_EPI])

        ram_logic = epi_fb
        fanout    = epi_fw - epi_fb
        wiring    = epi_ff - epi_fb
        other     = epi_b - ram_logic - fanout - wiring

        cycles = npred + extra
        # convert EPI (fJ/instr) to energy (fJ), then to EPPC (fJ/cycle)
        if cycles > 0:
            eppc_base   = epi_b  * ninstr / cycles
            eppc_ram    = ram_logic * ninstr / cycles
            eppc_fanout = fanout * ninstr / cycles
            eppc_wiring = wiring * ninstr / cycles
            eppc_other  = other  * ninstr / cycles
        else:
            eppc_base = eppc_ram = eppc_fanout = eppc_wiring = eppc_other = 0.0

        mpki = 1000.0 * misp / ninstr if ninstr > 0 else 0.0
        rows.append({
            "trace": n,
            "ninstr": int(ninstr),
            "npred": int(npred),
            "extra": int(extra),
            "mpki": mpki,
            "epi_baseline": epi_b,
            "epi_ram_logic": ram_logic,
            "epi_fanout": fanout,
            "epi_wiring": wiring,
            "epi_other": other,
            "eppc_baseline": eppc_base,
            "eppc_ram_logic": eppc_ram,
            "eppc_fanout": eppc_fanout,
            "eppc_wiring": eppc_wiring,
            "eppc_other": eppc_other,
        })

    # 5. Write CSV
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            r2 = {**r}
            for k, v in r2.items():
                if isinstance(v, float):
                    r2[k] = f"{v:.3f}"
            w.writerow(r2)
    print(f"\nCSV → {args.out}  ({len(rows)} traces)")

    # 6. Aggregate (instruction-weighted)
    tot_inst   = sum(r["ninstr"] for r in rows)
    tot_cycles = sum(r["npred"] + r["extra"] for r in rows)
    def w_epi(key):
        return sum(r[key] * r["ninstr"] for r in rows) / tot_inst
    def w_eppc(key):
        return sum(r[key.replace("eppc", "epi")] * r["ninstr"]
                   for r in rows) / tot_cycles

    base  = w_epi("epi_baseline")
    ram   = w_epi("epi_ram_logic")
    fan   = w_epi("epi_fanout")
    wir   = w_epi("epi_wiring")
    oth   = w_epi("epi_other")

    eppc_base = w_eppc("eppc_baseline")
    eppc_ram  = w_eppc("eppc_ram_logic")
    eppc_fan  = w_eppc("eppc_fanout")
    eppc_wir  = w_eppc("eppc_wiring")
    eppc_oth  = w_eppc("eppc_other")

    print(f"\n{'='*60}")
    print(f"  Aggregate (instruction-weighted, {len(rows)} traces)")
    print(f"{'='*60}")
    print(f"  {'Component':<14} {'EPI (fJ/inst)':>16} {'EPPC (fJ/cyc)':>16}")
    print(f"  {'-'*14} {'-'*16} {'-'*16}")
    for label, e, ep in [
        ("Baseline",  base, eppc_base),
        ("RAM+Logic", ram,  eppc_ram),
        ("Fanout",    fan,  eppc_fan),
        ("Wiring",    wir,  eppc_wir),
        ("Other",     oth,  eppc_oth),
    ]:
        pct = 100.0 * e / base if base else 0.0
        bar = "█" * int(pct / 2)
        print(f"  {label:<14} {e:>11.1f}  ({pct:4.1f}%) {ep:>16.1f}  {bar}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
