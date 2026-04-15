#!/usr/bin/env python3
"""
param_sweep.py — Iterative parameter sweep for TageAhead predictor.

Each iteration:
  Phase 1: Compile + timing check all configs (parallel).
  Phase 2: Quick-eval on timing-passing configs (parallel).
  Auto-pause after each iteration; best config becomes new base.

Usage:
  python3 scripts/param_sweep.py --init                 # Create default sweep_config.json
  python3 scripts/param_sweep.py --generate             # Generate configs for review/prune
  python3 scripts/param_sweep.py --timing -j N          # Phase 1: compile + timing check only
  python3 scripts/param_sweep.py --eval -j N            # Phase 2: quick-eval on timing-passed configs
  python3 scripts/param_sweep.py --run -j N             # Both phases back-to-back
  python3 scripts/param_sweep.py --resume -j N          # Resume interrupted phase / next iteration
  python3 scripts/param_sweep.py --report [--top 20]    # Print results ranked by VFS

Workflow:
  1. --generate           → creates configs.json (prune if desired)
  2. --timing             → compile + timing check, prints pass/fail summary
     (prune configs.json: remove configs you don't want to eval)
  3. --eval               → quick-eval on remaining timing-passed configs
  4. Review results, then --resume to start next iteration with updated base
"""

import argparse
import hashlib
import itertools
import json
import math
import os
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

# ---------------------------------------------------------------------------
# Representative traces (same 20 as quick_eval.sh)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------
TIMING_WARMUP = 4000
TIMING_MEASURE = 4000
QUICK_WARMUP = 1000000
QUICK_MEASURE = 40000000

# ---------------------------------------------------------------------------
# Compilation
# ---------------------------------------------------------------------------
CXX = "g++"
CXXFLAGS = [
    "-std=c++20", "-O3",
    "-Wall", "-Wextra", "-pedantic", "-Wold-style-cast", "-Werror",
    "-Wno-deprecated-declarations", "-Wno-mismatched-tags",
    "-Itrace_files",
]

# ---------------------------------------------------------------------------
# Graceful shutdown
# ---------------------------------------------------------------------------
_shutdown = False


def _sigint_handler(sig, frame):
    global _shutdown
    if _shutdown:
        print("\nForce exit.", file=sys.stderr)
        sys.exit(1)
    print("\nShutting down after current tasks complete...", file=sys.stderr)
    _shutdown = True


signal.signal(signal.SIGINT, _sigint_handler)

# ---------------------------------------------------------------------------
# Minimum values for perturbed params (clamp base-delta)
# ---------------------------------------------------------------------------
MIN_VALUES = {
    "num_tables": 2,
    "table_size": 64,
    "tag_width": 1,
    "min_hist": 2,
    "max_hist": 20,
    "fb_capacity": 256,
    "gs_hist": 1,
    "ctr_width": 1,
    "hyst_width": 1,
    "meta_width": 1,
    "pathbits": 1,
}

# ===================================================================
# PREDICTOR_TYPE string generation
# ===================================================================


def _trace_name(trace_file):
    """502-gcc-all_16112_trace.gz -> 502-gcc-all_16112"""
    name = trace_file
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith("_trace"):
        name = name[:-6]
    return name


def _resolve_tag_fn(name, tag_width, num_tables):
    t = tag_width
    nt = num_tables
    fns = {
        "UniformTag": f"ta::UniformTag<{t}>",
        "GradedTag": f"ta::GradedTag<{t},{max(1, t - 1)}>",
        "StepTag": f"ta::StepTag<{t},{max(1, t - 2)},{nt // 2}>",
        "LogTag": f"ta::LogTag<{max(1, t - 2)},1>",
    }
    if name not in fns:
        raise ValueError(f"Unknown tag_fn: {name}")
    return fns[name]


def _resolve_size_fn(name, table_size, size_ratio, num_tables, cfg):
    s = table_size
    sr = size_ratio
    nt = num_tables
    fns = {
        "UniformSize": f"ta::UniformSize<{s}>",
        "GeoSize": f"ta::GeoSize<{s},{sr}>",
        "InvGeoSize": f"ta::InvGeoSize<{s},{sr}>",
        "StepSize": f"ta::StepSize<{s * 2},{max(64, s // 2)},{nt // 2}>",
        "SqrtHistSize": f"ta::SqrtHistSize<{s}>",
        "ConstBitsSize": _const_bits_size(s, cfg),
    }
    if name not in fns:
        raise ValueError(f"Unknown size_fn: {name}")
    return fns[name]


def _const_bits_size(table_size, cfg):
    t = cfg["tag_width"]
    c = cfg["ctr_width"]
    h = cfg["hyst_width"]
    u = cfg.get("u_width", 1)
    total = table_size * (t + c + h + u)
    return f"ta::ConstBitsSize<{total},{t},{c},{h},{u}>"


def generate_predictor_type(cfg):
    """Build the full TageAhead<...> C++ template string from a config dict."""
    nt = cfg["num_tables"]
    s = cfg["table_size"]
    t = cfg["tag_width"]
    minh = cfg["min_hist"]
    maxh = cfg["max_hist"]
    sr = cfg.get("size_ratio", 1)

    tag_fn = _resolve_tag_fn(cfg["tag_fn"], t, nt)
    size_fn = _resolve_size_fn(cfg["size_fn"], s, sr, nt, cfg)

    table_cfg = (
        f"TATableConfig<{nt},{s},{t},{minh},{maxh},{sr},"
        f"ta::HistSeries::{cfg['hist_series']},{tag_fn},{size_fn}>"
    )

    epoch = cfg["u_reset"] == "epoch"

    parts = [
        table_cfg,
        str(cfg.get("n", 8)),
        str(cfg.get("pathbits", 6)),
        str(cfg.get("sec_tag_bits", 3)),
        "true" if cfg.get("use_sec_tag", True) else "false",
        str(cfg["ctr_width"]),
        str(cfg["hyst_width"]),
        str(cfg.get("u_width", 1)),
        str(cfg["fb_capacity"]),
        "true" if cfg["use_gshare"] else "false",
        str(cfg["gs_hist"]),
        str(cfg["meta_width"]),
        str(cfg.get("meta_capacity", 256)),
        str(cfg.get("meta_pipe", 2)),
        str(cfg.get("lineinst", 1024)),
        "true" if cfg.get("shared_hys", True) else "false",
        f"HistUpdate::{cfg.get('hist_mode', 'PATH')}",
        cfg["alloc_cfg"],
        str(cfg.get("acc_width", 4)),
        str(cfg.get("alloc_width", 4)),
        # decay
        "false" if epoch else "true",
        f"DecayMiss::{cfg.get('decay_miss', 'TAG_OR_SEC')}",
        f"DecayOp::{cfg.get('decay_op', 'DECREMENT')}",
        f"ta::uniform_array<u64,{nt}>(8)",
        "ta::DefaultDecayThresh",
        # epoch
        "true" if epoch else "false",
        "ta::DefaultEpochTrigger",
        "false",
        "true",
    ]

    return f"TageAhead<{','.join(parts)}>"


# ===================================================================
# Config generation
# ===================================================================


def _normalize(cfg, base):
    """Set irrelevant params to base values so duplicates collapse."""
    c = dict(cfg)
    if not c["use_gshare"]:
        c["gs_hist"] = base["gs_hist"]
    if c["u_reset"] == "epoch":
        c["decay_miss"] = base.get("decay_miss", "TAG_OR_SEC")
        c["decay_op"] = base.get("decay_op", "DECREMENT")
    return c


def config_hash(cfg):
    """8-char MD5 of config (excluding internal _ fields)."""
    kvs = sorted((k, v) for k, v in cfg.items() if not k.startswith("_"))
    return hashlib.md5(json.dumps(kvs).encode()).hexdigest()[:8]


def generate_configs(sweep_cfg):
    """Generate all ALWAYS x PERTURBED configs, deduplicated."""
    base = sweep_cfg["base"]
    always = sweep_cfg["always"]
    perturbed = sweep_cfg["perturbed"]

    always_keys = sorted(always.keys())
    always_combos = list(itertools.product(*[always[k] for k in always_keys]))

    configs = {}  # hash -> config dict

    for p_name, p_spec in sorted(perturbed.items()):
        delta = p_spec["delta"]
        base_val = base[p_name]

        for offset in [-delta, 0, delta]:
            raw = base_val + offset
            p_val = max(raw, MIN_VALUES.get(p_name, 1))
            # skip if clamped value equals another offset's value (avoids dup work)
            if p_val == base_val and offset != 0:
                continue

            for combo in always_combos:
                cfg = dict(base)
                cfg[p_name] = p_val
                for i, ak in enumerate(always_keys):
                    cfg[ak] = combo[i]

                cfg = _normalize(cfg, base)
                h = config_hash(cfg)
                if h not in configs:
                    cfg["_hash"] = h
                    cfg["_perturbed"] = p_name
                    cfg["_offset"] = offset
                    configs[h] = cfg

    return list(configs.values())


# ===================================================================
# Compilation & timing
# ===================================================================


def _compile(cfg, bin_dir, extra_flags):
    """Compile a config. Returns (binary_path, ok, stderr)."""
    h = cfg["_hash"]
    binary = bin_dir / f"cbp-{h}"
    if binary.exists():
        return str(binary), True, ""

    pred = generate_predictor_type(cfg)
    cmd = (
        [CXX]
        + CXXFLAGS
        + (extra_flags.split() if extra_flags else [])
        + ["-o", str(binary), "cbp.cpp", "-lz", f"-DPREDICTOR={pred}"]
    )
    try:
        r = subprocess.run(
            cmd, capture_output=True, text=True, cwd=str(PROJECT_DIR), timeout=300
        )
    except subprocess.TimeoutExpired:
        return str(binary), False, "compile timeout"
    if r.returncode != 0:
        return str(binary), False, r.stderr[:1000]
    return str(binary), True, ""


def _timing_check(binary, trace_path, trace_name):
    """Run 4k/4k timing probe. Returns (status, p1, p2)."""
    cmd = [
        binary,
        str(trace_path),
        trace_name,
        str(TIMING_WARMUP),
        str(TIMING_MEASURE),
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return "timeout", 0.0, 0.0
    if r.returncode != 0:
        return "crash", 0.0, 0.0

    try:
        line = r.stdout.strip().split("\n")[0]
        fields = line.split(",")
        p1 = float(fields[9])
        p2 = float(fields[10])
    except (IndexError, ValueError):
        return "parse_error", 0.0, 0.0

    if math.ceil(p1) > 1 or math.ceil(p2) > 1:
        return "timing_fail", p1, p2
    return "pass", p1, p2


def _compile_and_time(cfg, bin_dir, trace_path, trace_name, extra_flags):
    """Combined compile + timing for one config (runs in thread pool)."""
    binary, ok, err = _compile(cfg, bin_dir, extra_flags)
    if not ok:
        return cfg["_hash"], "compile_fail", 0.0, 0.0, err
    status, p1, p2 = _timing_check(binary, trace_path, trace_name)
    return cfg["_hash"], status, p1, p2, ""


# ===================================================================
# Quick eval
# ===================================================================


def _run_trace(binary, trace_path, trace_name, out_path):
    """Run binary on one trace, write output file. Returns success bool."""
    cmd = [
        binary,
        str(trace_path),
        trace_name,
        str(QUICK_WARMUP),
        str(QUICK_MEASURE),
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return False
    if r.returncode == 0 and r.stdout.strip():
        Path(out_path).write_text(r.stdout)
        return True
    return False


def _compute_metrics(out_dir):
    """Run predictor_metrics.py + vfs.py on a directory of .out files."""
    r = subprocess.run(
        ["python3", "predictor_metrics.py", str(out_dir)],
        capture_output=True,
        text=True,
        cwd=str(PROJECT_DIR),
    )
    if r.returncode != 0:
        return None
    metrics_line = r.stdout.strip()
    if not metrics_line:
        return None

    r2 = subprocess.run(
        ["python3", "vfs.py", metrics_line],
        capture_output=True,
        text=True,
        cwd=str(PROJECT_DIR),
    )
    if r2.returncode != 0:
        return None

    vfs = float(r2.stdout.strip())
    f = metrics_line.split(",")
    return {
        "ipc": float(f[0]),
        "cpi": float(f[1]),
        "epi": float(f[2]),
        "mpi": float(f[3]),
        "dpi": float(f[4]),
        "ppi": float(f[5]),
        "p1_lat": int(f[6]),
        "p2_lat": int(f[7]),
        "vfs": vfs,
    }


# ===================================================================
# Checkpoint
# ===================================================================


def _load_checkpoint(path):
    if path.exists():
        return json.loads(path.read_text())
    return {"iteration": 0, "phase": "init", "timing": {}, "eval": {}}


def _save_checkpoint(ckpt, path):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(ckpt, indent=2))
    tmp.rename(path)


# ===================================================================
# Results CSV
# ===================================================================

CSV_COLS = [
    "iteration",
    "hash",
    "perturbed",
    "offset",
    "alloc_cfg",
    "use_gshare",
    "u_reset",
    "hist_series",
    "tag_fn",
    "size_fn",
    "num_tables",
    "table_size",
    "tag_width",
    "min_hist",
    "max_hist",
    "fb_capacity",
    "gs_hist",
    "ctr_width",
    "hyst_width",
    "meta_width",
    "pathbits",
    "status",
    "p1_lat",
    "p2_lat",
    "vfs",
    "ipc",
    "cpi",
    "epi",
    "mpi",
]


def _ensure_csv(path):
    if not path.exists():
        path.write_text(",".join(CSV_COLS) + "\n")


def _append_result(csv_path, iteration, cfg, status, p1, p2, metrics=None):
    _ensure_csv(csv_path)
    m = metrics or {}
    vals = [
        iteration,
        cfg["_hash"],
        cfg.get("_perturbed", ""),
        cfg.get("_offset", 0),
        cfg["alloc_cfg"],
        cfg["use_gshare"],
        cfg["u_reset"],
        cfg["hist_series"],
        cfg["tag_fn"],
        cfg["size_fn"],
        cfg["num_tables"],
        cfg["table_size"],
        cfg["tag_width"],
        cfg["min_hist"],
        cfg["max_hist"],
        cfg["fb_capacity"],
        cfg["gs_hist"],
        cfg["ctr_width"],
        cfg["hyst_width"],
        cfg["meta_width"],
        cfg.get("pathbits", 6),
        status,
        f"{p1:.4f}",
        f"{p2:.4f}",
        f"{m.get('vfs', 0):.6f}",
        f"{m.get('ipc', 0):.6f}",
        f"{m.get('cpi', 0):.6f}",
        f"{m.get('epi', 0):.6f}",
        f"{m.get('mpi', 0):.6f}",
    ]
    with open(csv_path, "a") as f:
        f.write(",".join(str(v) for v in vals) + "\n")


# ===================================================================
# Iteration runner — split into Phase 1 (timing) and Phase 2 (eval)
# ===================================================================


def _iter_setup(sweep_cfg, iteration, out_dir):
    """Common setup: dirs, checkpoint, configs. Returns (iter_dir, bin_dir,
    ckpt_path, results_path, configs, ckpt) or None on failure."""
    iter_dir = out_dir / f"iter_{iteration}"
    bin_dir = out_dir / "bin"
    iter_dir.mkdir(parents=True, exist_ok=True)
    bin_dir.mkdir(parents=True, exist_ok=True)

    ckpt_path = out_dir / "checkpoint.json"
    results_path = out_dir / "results.csv"
    configs_path = iter_dir / "configs.json"
    ckpt = _load_checkpoint(ckpt_path)

    if configs_path.exists():
        configs = json.loads(configs_path.read_text())
        print(f"Loaded {len(configs)} configs from {configs_path}")
    else:
        configs = generate_configs(sweep_cfg)
        configs_path.write_text(json.dumps(configs, indent=2))
        print(f"Generated {len(configs)} configs -> {configs_path}")

    if not configs:
        print("No configs to run.")
        return None
    return iter_dir, bin_dir, ckpt_path, results_path, configs, ckpt


def run_timing(sweep_cfg, iteration, jobs, out_dir, extra_flags=""):
    """Phase 1: compile + timing check. Sets phase to 'timing_done'."""
    global _shutdown

    setup = _iter_setup(sweep_cfg, iteration, out_dir)
    if setup is None:
        return False
    iter_dir, bin_dir, ckpt_path, results_path, configs, ckpt = setup

    trace_dir = Path(sweep_cfg.get("trace_dir", "./traces"))
    timing_trace = sweep_cfg.get("timing_trace", REPR_TRACES[0])
    timing_trace_path = trace_dir / timing_trace
    timing_trace_name = _trace_name(timing_trace)

    if not timing_trace_path.exists():
        print(f"ERROR: timing trace not found: {timing_trace_path}", file=sys.stderr)
        return False

    timing_todo = [c for c in configs if c["_hash"] not in ckpt["timing"]]
    total = len(configs)
    done = total - len(timing_todo)

    print(f"\n{'=' * 60}")
    print(f"Phase 1: Timing Check  ({total} configs, {len(timing_todo)} remaining, {jobs} jobs)")
    print(f"{'=' * 60}")

    if timing_todo:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {}
            for cfg in timing_todo:
                if _shutdown:
                    break
                f = pool.submit(
                    _compile_and_time,
                    cfg,
                    bin_dir,
                    timing_trace_path,
                    timing_trace_name,
                    extra_flags,
                )
                futures[f] = cfg

            for f in as_completed(futures):
                if _shutdown:
                    pool.shutdown(wait=True, cancel_futures=True)
                    _save_checkpoint(ckpt, ckpt_path)
                    print("Checkpoint saved. Resume with --resume.")
                    return False

                h, status, p1, p2, err = f.result()
                ckpt["timing"][h] = {"status": status, "p1": p1, "p2": p2}
                done += 1

                sym = {"pass": "OK", "compile_fail": "CE", "timing_fail": "TF",
                       "timeout": "TO", "crash": "CR", "parse_error": "PE"}.get(status, "??")
                print(f"  [{done}/{total}] {h} {sym}  p1={p1:.3f} p2={p2:.3f}")
                if err and status == "compile_fail":
                    first_err = err.strip().split("\n")[-1][:120]
                    print(f"           {first_err}")

                _save_checkpoint(ckpt, ckpt_path)

    # Summary
    passed = [c for c in configs if ckpt["timing"].get(c["_hash"], {}).get("status") == "pass"]
    failed_counts = {}
    for c in configs:
        s = ckpt["timing"].get(c["_hash"], {}).get("status", "unknown")
        if s != "pass":
            failed_counts[s] = failed_counts.get(s, 0) + 1

    print(f"\nTiming: {len(passed)} pass, {total - len(passed)} fail")
    if failed_counts:
        print(f"  Failures: {failed_counts}")

    # Log timing failures to CSV
    for cfg in configs:
        h = cfg["_hash"]
        t = ckpt["timing"].get(h, {})
        if t.get("status") != "pass":
            _append_result(
                results_path, iteration, cfg, t.get("status", "unknown"),
                t.get("p1", 0), t.get("p2", 0),
            )

    ckpt["phase"] = "timing_done"
    _save_checkpoint(ckpt, ckpt_path)

    if not passed:
        print("No configs passed timing.")
        return False

    configs_path = iter_dir / "configs.json"
    print(f"\nPrune {configs_path} to remove configs you don't want to eval, then:")
    print(f"  python3 {sys.argv[0]} --eval -j {jobs}")
    return True


def run_eval(sweep_cfg, iteration, jobs, out_dir, extra_flags=""):
    """Phase 2: quick-eval on timing-passed configs. Sets phase to 'done'."""
    global _shutdown

    setup = _iter_setup(sweep_cfg, iteration, out_dir)
    if setup is None:
        return None
    iter_dir, bin_dir, ckpt_path, results_path, configs, ckpt = setup

    trace_dir = Path(sweep_cfg.get("trace_dir", "./traces"))

    # Only eval configs that are still in configs.json AND passed timing
    passed_cfgs = [
        c for c in configs
        if ckpt["timing"].get(c["_hash"], {}).get("status") == "pass"
    ]
    eval_todo = [c for c in passed_cfgs if c["_hash"] not in ckpt["eval"]]
    eval_total = len(passed_cfgs)
    eval_done = eval_total - len(eval_todo)

    print(f"\n{'=' * 60}")
    print(f"Phase 2: Quick Eval  ({eval_total} configs, {len(eval_todo)} remaining, {jobs} jobs)")
    print(f"{'=' * 60}")

    if eval_total == 0:
        print("No timing-passed configs in configs.json. Nothing to eval.")
        return None

    for cfg in eval_todo:
        if _shutdown:
            _save_checkpoint(ckpt, ckpt_path)
            print("Checkpoint saved. Resume with --resume.")
            return None

        h = cfg["_hash"]
        binary = str(bin_dir / f"cbp-{h}")
        if not Path(binary).exists():
            print(f"  {h}: binary missing, recompile needed. Skipping.")
            continue

        qe_dir = iter_dir / h / "quick_eval"
        qe_dir.mkdir(parents=True, exist_ok=True)

        # Collect trace tasks that still need running
        trace_tasks = []
        for trace in REPR_TRACES:
            tp = trace_dir / trace
            if not tp.exists():
                continue
            tname = _trace_name(trace)
            out_path = qe_dir / f"{tname}.out"
            if not out_path.exists() or out_path.stat().st_size == 0:
                trace_tasks.append((binary, str(tp), tname, str(out_path)))

        if trace_tasks:
            with ThreadPoolExecutor(max_workers=jobs) as pool:
                futs = {pool.submit(_run_trace, *t): t for t in trace_tasks}
                for ft in as_completed(futs):
                    if _shutdown:
                        pool.shutdown(wait=True, cancel_futures=True)
                        break
                    ft.result()

        if _shutdown:
            _save_checkpoint(ckpt, ckpt_path)
            print("Checkpoint saved. Resume with --resume.")
            return None

        metrics = _compute_metrics(qe_dir)
        if metrics:
            ckpt["eval"][h] = metrics
            t = ckpt["timing"][h]
            _append_result(results_path, iteration, cfg, "pass", t["p1"], t["p2"], metrics)
        else:
            ckpt["eval"][h] = {"vfs": 0.0, "error": True}
            t = ckpt["timing"][h]
            _append_result(results_path, iteration, cfg, "eval_fail", t["p1"], t["p2"])

        eval_done += 1
        vfs = metrics["vfs"] if metrics else 0.0
        print(f"  [{eval_done}/{eval_total}] {h} VFS={vfs:.6f}")
        _save_checkpoint(ckpt, ckpt_path)

    # ------------------------------------------------------------------
    # Rank results
    # ------------------------------------------------------------------
    ranked = sorted(
        [(h, m) for h, m in ckpt["eval"].items() if m.get("vfs", 0) > 0],
        key=lambda x: x[1]["vfs"],
        reverse=True,
    )

    if not ranked:
        print("No configs produced valid VFS scores.")
        return None

    print(f"\n{'=' * 60}")
    print(f"Iteration {iteration} Results — Top 10")
    print(f"{'=' * 60}")
    print(f"  {'Hash':>8}  {'VFS':>10}  {'IPC':>8}  {'MPI':>10}  {'EPI':>8}  {'P1':>3}  {'P2':>3}")
    for h, m in ranked[:10]:
        print(
            f"  {h:>8}  {m['vfs']:10.6f}  {m['ipc']:8.4f}  {m['mpi']:10.6f}"
            f"  {m['epi']:8.2f}  {m['p1_lat']:>3}  {m['p2_lat']:>3}"
        )

    best_h, best_m = ranked[0]
    worst_h, worst_m = ranked[-1]
    best_cfg = next(c for c in configs if c["_hash"] == best_h)

    print(f"\nBest:  {best_h}  VFS={best_m['vfs']:.6f}")
    print(f"Worst: {worst_h}  VFS={worst_m['vfs']:.6f}")

    # Update base to best config
    new_base = {k: v for k, v in best_cfg.items() if not k.startswith("_")}
    sweep_cfg["base"] = new_base

    ckpt["phase"] = "done"
    ckpt["best_hash"] = best_h
    ckpt["best_vfs"] = best_m["vfs"]
    _save_checkpoint(ckpt, ckpt_path)

    return best_cfg


# ===================================================================
# Report
# ===================================================================


def do_report(out_dir, top_n=20):
    results_path = out_dir / "results.csv"
    if not results_path.exists():
        print(f"No results found at {results_path}")
        return

    import csv

    with open(results_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Filter to passing configs with VFS
    passing = [r for r in rows if r["status"] == "pass" and float(r["vfs"]) > 0]
    passing.sort(key=lambda r: float(r["vfs"]), reverse=True)

    if not passing:
        print("No passing configs found.")
        return

    print(f"\n{'=' * 80}")
    print(f"Top {min(top_n, len(passing))} configs by VFS")
    print(f"{'=' * 80}")
    print(
        f"  {'#':>3}  {'Iter':>4}  {'Hash':>8}  {'VFS':>10}  {'IPC':>8}  {'MPI':>10}"
        f"  {'EPI':>8}  {'Perturbed':>14}  {'Δ':>4}"
    )
    print(f"  {'---':>3}  {'----':>4}  {'--------':>8}  {'----------':>10}  {'--------':>8}"
          f"  {'----------':>10}  {'--------':>8}  {'--------------':>14}  {'----':>4}")

    for i, r in enumerate(passing[:top_n]):
        print(
            f"  {i + 1:>3}  {r['iteration']:>4}  {r['hash']:>8}"
            f"  {float(r['vfs']):10.6f}  {float(r['ipc']):8.4f}"
            f"  {float(r['mpi']):10.6f}  {float(r['epi']):8.2f}"
            f"  {r['perturbed']:>14}  {r['offset']:>4}"
        )

    # Summary by perturbed param
    print(f"\n{'=' * 80}")
    print("Best VFS by perturbed param")
    print(f"{'=' * 80}")
    by_param = {}
    for r in passing:
        p = r["perturbed"] or "(base)"
        if p not in by_param or float(r["vfs"]) > float(by_param[p]["vfs"]):
            by_param[p] = r
    for p in sorted(by_param, key=lambda p: float(by_param[p]["vfs"]), reverse=True):
        r = by_param[p]
        print(f"  {p:>14}: VFS={float(r['vfs']):10.6f}  hash={r['hash']}")

    # Failure summary
    failures = {}
    for r in rows:
        s = r["status"]
        if s != "pass":
            failures[s] = failures.get(s, 0) + 1
    if failures:
        print(f"\nFailures: {failures}")


# ===================================================================
# Init
# ===================================================================

DEFAULT_CONFIG = {
    "base": {
        "num_tables": 12,
        "table_size": 1024,
        "tag_width": 11,
        "min_hist": 8,
        "max_hist": 100,
        "size_ratio": 1,
        "n": 8,
        "pathbits": 6,
        "sec_tag_bits": 3,
        "use_sec_tag": True,
        "ctr_width": 1,
        "hyst_width": 2,
        "u_width": 1,
        "fb_capacity": 8192,
        "use_gshare": True,
        "gs_hist": 6,
        "meta_width": 4,
        "meta_capacity": 256,
        "meta_pipe": 2,
        "lineinst": 1024,
        "shared_hys": True,
        "hist_mode": "PATH",
        "alloc_cfg": "TADefaultAllocConfig",
        "acc_width": 4,
        "alloc_width": 4,
        "hist_series": "GEOMETRIC",
        "tag_fn": "GradedTag",
        "size_fn": "GeoSize",
        "u_reset": "epoch",
        "decay_miss": "TAG_OR_SEC",
        "decay_op": "DECREMENT",
    },
    "always": {
        "alloc_cfg": [
            "TADefaultAllocConfig",
            "TAAllocDetSkip1",
            "TAAlloc2",
            "TAAllocTageWrong",
            "TAAllocFiltered",
            "TAAllocThrottled",
        ],
        "use_gshare": [True, False],
        "u_reset": ["epoch", "decay"],
        "hist_series": ["GEOMETRIC", "QUADRATIC", "SUPEREXP", "ROS"],
        "tag_fn": ["UniformTag", "GradedTag", "StepTag", "LogTag"],
        "size_fn": ["UniformSize", "GeoSize", "InvGeoSize", "StepSize", "SqrtHistSize"],
    },
    "perturbed": {
        "num_tables": {"delta": 2},
        "table_size": {"delta": 512},
        "tag_width": {"delta": 2},
        "min_hist": {"delta": 2},
        "max_hist": {"delta": 20},
        "fb_capacity": {"delta": 4096},
        "gs_hist": {"delta": 4},
        "ctr_width": {"delta": 1},
        "hyst_width": {"delta": 1},
        "meta_width": {"delta": 1},
        "pathbits": {"delta": 2},
    },
    "timing_trace": "502-gcc-all_16112_trace.gz",
    "trace_dir": "./traces",
    "out_dir": "out/sweep",
    "extra_flags": "",
}


# ===================================================================
# CLI
# ===================================================================


def main():
    parser = argparse.ArgumentParser(
        description="TageAhead iterative parameter sweep",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config", default="sweep_config.json", help="Sweep config file"
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=os.cpu_count(), help="Parallel workers"
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--init", action="store_true", help="Create default config")
    group.add_argument("--generate", action="store_true", help="Generate configs only")
    group.add_argument("--timing", action="store_true", help="Phase 1: compile + timing check")
    group.add_argument("--eval", action="store_true", help="Phase 2: quick-eval (after timing)")
    group.add_argument("--run", action="store_true", help="Both phases back-to-back")
    group.add_argument("--resume", action="store_true", help="Resume interrupted phase / next iter")
    group.add_argument("--report", action="store_true", help="Print results report")

    parser.add_argument("--top", type=int, default=20, help="Top N for report")

    args = parser.parse_args()
    config_path = Path(args.config)

    # ---- init ----
    if args.init:
        if config_path.exists():
            print(f"{config_path} already exists. Delete it first to re-init.")
            return
        config_path.write_text(json.dumps(DEFAULT_CONFIG, indent=2))
        print(f"Created {config_path}")
        print(f"Edit it, then: python3 {sys.argv[0]} --generate")
        return

    # Load config (needed for all other modes)
    if not config_path.exists():
        print(f"Config not found: {config_path}")
        print(f"Run: python3 {sys.argv[0]} --init")
        return
    sweep_cfg = json.loads(config_path.read_text())
    out_dir = Path(sweep_cfg.get("out_dir", "out/sweep"))

    # ---- report ----
    if args.report:
        do_report(out_dir, args.top)
        return

    # Helper: resolve current iteration from checkpoint
    def _current_iteration(advance_if_done=False):
        ckpt_path = out_dir / "checkpoint.json"
        ckpt = _load_checkpoint(ckpt_path)
        it = ckpt.get("iteration", 0)
        if advance_if_done and ckpt.get("phase") == "done":
            it += 1
            ckpt = {"iteration": it, "phase": "init", "timing": {}, "eval": {}}
            _save_checkpoint(ckpt, ckpt_path)
        return it, ckpt

    extra = sweep_cfg.get("extra_flags", "")

    # ---- generate ----
    if args.generate:
        it, ckpt = _current_iteration(advance_if_done=True)
        iter_dir = out_dir / f"iter_{it}"
        iter_dir.mkdir(parents=True, exist_ok=True)
        configs_path = iter_dir / "configs.json"

        configs = generate_configs(sweep_cfg)
        configs_path.write_text(json.dumps(configs, indent=2))

        always = sweep_cfg["always"]
        always_count = 1
        for v in always.values():
            always_count *= len(v)

        print(f"Iteration {it}: {len(configs)} unique configs")
        print(f"  ALWAYS combos: {always_count}")
        print(f"  PERTURBED params: {len(sweep_cfg['perturbed'])}")
        print(f"  Saved to: {configs_path}")
        print(f"\nPrune {configs_path} if desired, then:")
        print(f"  python3 {sys.argv[0]} --timing -j {args.jobs}")
        return

    # ---- timing (Phase 1 only) ----
    if args.timing:
        it, ckpt = _current_iteration(advance_if_done=True)
        if ckpt.get("phase") == "timing_done":
            print("Timing already complete. Run --eval or prune configs.json first.")
            return
        run_timing(sweep_cfg, it, args.jobs, out_dir, extra)
        return

    # ---- eval (Phase 2 only) ----
    if args.eval:
        it, ckpt = _current_iteration()
        if ckpt.get("phase") not in ("timing_done", "eval"):
            print("Run --timing first (or --run for both phases).")
            return
        best = run_eval(sweep_cfg, it, args.jobs, out_dir, extra)
        if best:
            config_path.write_text(json.dumps(sweep_cfg, indent=2))
            print(f"\nBase config updated in {config_path}")
            print(f"Iteration {it} complete. Auto-paused.")
            print(f"To start next iteration: python3 {sys.argv[0]} --generate")
        return

    # ---- run (both phases back-to-back) ----
    if args.run:
        it, ckpt = _current_iteration(advance_if_done=True)
        if ckpt.get("phase") not in ("init", None, "done"):
            print(f"Sweep in progress (phase={ckpt.get('phase')}). Use --resume.")
            return

        ok = run_timing(sweep_cfg, it, args.jobs, out_dir, extra)
        if not ok:
            return
        best = run_eval(sweep_cfg, it, args.jobs, out_dir, extra)
        if best:
            config_path.write_text(json.dumps(sweep_cfg, indent=2))
            print(f"\nBase config updated in {config_path}")
            print(f"Iteration {it} complete. Auto-paused.")
            print(f"To continue: python3 {sys.argv[0]} --resume -j {args.jobs}")
        return

    # ---- resume (continue wherever we left off) ----
    if args.resume:
        it, ckpt = _current_iteration(advance_if_done=True)
        phase = ckpt.get("phase", "init")

        if phase in ("init", None):
            # Start from timing
            ok = run_timing(sweep_cfg, it, args.jobs, out_dir, extra)
            if not ok:
                return
            best = run_eval(sweep_cfg, it, args.jobs, out_dir, extra)
        elif phase == "timing_done":
            best = run_eval(sweep_cfg, it, args.jobs, out_dir, extra)
        elif phase == "eval":
            best = run_eval(sweep_cfg, it, args.jobs, out_dir, extra)
        else:
            # phase == "timing" (interrupted during timing)
            ok = run_timing(sweep_cfg, it, args.jobs, out_dir, extra)
            if not ok:
                return
            best = run_eval(sweep_cfg, it, args.jobs, out_dir, extra)

        if best:
            config_path.write_text(json.dumps(sweep_cfg, indent=2))
            print(f"\nBase config updated in {config_path}")
            print(f"Iteration {it} complete. Auto-paused.")
            print(f"To continue: python3 {sys.argv[0]} --resume -j {args.jobs}")
        return


if __name__ == "__main__":
    main()
