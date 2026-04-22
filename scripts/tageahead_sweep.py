#!/usr/bin/env python3
"""
tageahead_sweep.py — Streaming parameter sweep for TageAhead predictor.

Pipeline: generate configs → build 2 at a time → timing-check on multiple
traces → keep/discard → repeat. Survivors get full 20-trace eval.

Usage:
  # 1. Create/edit sweep_config.json with sweep dimensions
  python3 scripts/tageahead_sweep.py --init

  # 2. Generate configs, then stream build+timing gate (2 at a time)
  python3 scripts/tageahead_sweep.py --run -j2

  # 3. Full eval on timing-passed configs
  python3 scripts/tageahead_sweep.py --eval -j2

  # 4. Report
  python3 scripts/tageahead_sweep.py --report [--top 20]
"""

import argparse
import hashlib
import itertools
import json
import math
import os
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
TRACE_DIR = PROJECT_DIR / "traces"
BUILD_DIR = PROJECT_DIR / "build" / "sweep"
RESULTS_DIR = PROJECT_DIR / "sweep_results"
CONFIG_FILE = PROJECT_DIR / "sweep_config.json"

CXX = os.environ.get("CXX", "g++")
COMMON_FLAGS = "-std=c++20 -O3"
WARN_FLAGS = ("-Wall -Wextra -pedantic -Wold-style-cast -Werror "
              "-Wno-deprecated-declarations -Wno-mismatched-tags")

# ---------------------------------------------------------------------------
# Timing gate traces — diverse workloads to catch timing regressions
# ---------------------------------------------------------------------------
TIMING_TRACES = [
    "502-gcc-all_16112_trace.gz",
    "531-deepsjeng-1_53703_trace.gz",
    "web_74_trace.gz",
    "nodejs-octane_3483_trace.gz",
    "lua-3.25585_0_trace.gz",
]
TIMING_WARMUP = 10000
TIMING_MEASURE = 10000

# Full eval traces
EVAL_TRACES = [
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
EVAL_WARMUP = 1_000_000
EVAL_MEASURE = 40_000_000

# VFS constants
MISP_PENALTY = 8


# ============================================================================
# TageAhead config → template string
# ============================================================================

# All TageAhead template params in order.  typename params use string values.
PARAM_ORDER = [
    # TATableConfig params (composed into TableCfg typename)
    "tc_num_tables",    # TATableConfig::N
    "tc_size",          # TATableConfig::SIZE
    "tc_tag",           # TATableConfig::TAG
    "tc_minh",          # TATableConfig::MINH
    "tc_maxh",          # TATableConfig::MAXH
    "tc_size_ratio",    # TATableConfig::SIZE_RATIO
    "tc_hist",          # TATableConfig::HIST (string: "GEOMETRIC", "QUADRATIC", etc.)
    "tc_tag_fn",        # TATableConfig::TagFn (string, e.g. "GradedTag<8,7>")
    "tc_size_fn",       # TATableConfig::SizeFn (string, e.g. "GeoSize<1024,1>")
    # TageAhead template params (positions 1-32)
    "N",
    "PATHBITS",
    "SEC_TAG_BITS",
    "USE_SEC_TAG",
    "NUM_PATHS",
    "SecTagHashFn",     # typename string
    "CTR_WIDTH",
    "HYST_WIDTH",
    "U_WIDTH",
    "FB_CAPACITY",
    "USE_GSHARE",
    "GS_HIST",
    "META_WIDTH",
    "META_CAPACITY",
    "META_PIPE",
    "LINEINST",
    "SHARED_HYS",
    "HIST_MODE",        # string: "PATH", "DIR", "BOTH"
    "AllocCfg",         # typename string
    "ACC_WIDTH",
    "ALLOC_WIDTH",
    "DECAY_ENABLE",
    "DECAY_MISS",       # string: "TAG", "SEC", etc.
    "DECAY_OP",         # string: "DECREMENT", "HALVE", "CLEAR"
    "DECAY_LFSR_WIDTH", # single u64 (applied uniformly)
    "DecayThreshFn",    # typename string
    "EPOCH_ENABLE",
    "EpochTriggerFn",   # typename string
    "EPOCH_RESET_ACC",
    "EPOCH_RESET_ALLOC",
    "FB_RECONCILE",
    "FARALLOC_DIST",
]

DEFAULTS = {
    # TATableConfig
    "tc_num_tables": 7,
    "tc_size": 1024,
    "tc_tag": 11,
    "tc_minh": 10,
    "tc_maxh": 100,
    "tc_size_ratio": 1,
    "tc_hist": "GEOMETRIC",
    "tc_tag_fn": "__default__",   # GradedTag<TAG, TAG-1>
    "tc_size_fn": "__default__",  # GeoSize<SIZE, SIZE_RATIO>
    # TageAhead
    "N": 8,
    "PATHBITS": 6,
    "SEC_TAG_BITS": 3,
    "USE_SEC_TAG": True,
    "NUM_PATHS": 1,
    "SecTagHashFn": "ta::DefaultSecTagHash",
    "CTR_WIDTH": 1,
    "HYST_WIDTH": 2,
    "U_WIDTH": 1,
    "FB_CAPACITY": 8192,
    "USE_GSHARE": True,
    "GS_HIST": 6,
    "META_WIDTH": 6,
    "META_CAPACITY": 1024,
    "META_PIPE": 2,
    "LINEINST": 1024,
    "SHARED_HYS": True,
    "HIST_MODE": "PATH",
    "AllocCfg": "TADefaultAllocConfig",
    "ACC_WIDTH": 4,
    "ALLOC_WIDTH": 4,
    "DECAY_ENABLE": False,
    "DECAY_MISS": "TAG_OR_SEC",
    "DECAY_OP": "DECREMENT",
    "DECAY_LFSR_WIDTH": 8,
    "DecayThreshFn": "ta::DefaultDecayThresh",
    "EPOCH_ENABLE": True,
    "EpochTriggerFn": "ta::DefaultEpochTrigger",
    "EPOCH_RESET_ACC": False,
    "EPOCH_RESET_ALLOC": True,
    "FB_RECONCILE": False,
    "FARALLOC_DIST": 0,
}


def bool_str(v):
    return "true" if v else "false"


# ANSI colors
GREEN = "\033[32m"
RED = "\033[31m"
YELLOW = "\033[33m"
CYAN = "\033[36m"
DIM = "\033[2m"
BOLD = "\033[1m"
RESET = "\033[0m"


# ============================================================================
# Table geometry computation (mirrors custom_common.hpp)
# ============================================================================

def _round_up_pow2(x: int) -> int:
    r = 64
    while r < x:
        r *= 2
    return r


def compute_table_sizes(cfg: dict) -> list[int]:
    """Compute per-table sizes matching C++ SizeFn logic."""
    n = cfg["tc_num_tables"]
    size = cfg["tc_size"]
    ratio = cfg["tc_size_ratio"]
    fn = cfg.get("tc_size_fn", "__default__")

    sizes = []
    for i in range(n):
        if fn in ("__default__", "GeoSize"):
            if ratio <= 1 or n <= 1:
                sizes.append(size)
            else:
                t = i / (n - 1)
                scale = ratio ** (t - 0.5)
                sizes.append(_round_up_pow2(int(size * scale)))
        elif fn == "InvGeoSize":
            if ratio <= 1 or n <= 1:
                sizes.append(size)
            else:
                t = (n - 1 - i) / (n - 1)
                scale = ratio ** (t - 0.5)
                sizes.append(_round_up_pow2(int(size * scale)))
        elif fn == "SqrtHistSize":
            if n <= 1:
                sizes.append(size)
            else:
                scale = (n - 1 - i + 1) ** 0.5
                sizes.append(_round_up_pow2(int(size * scale)))
        elif fn == "UniformSize":
            sizes.append(size)
        else:
            sizes.append(size)  # fallback
    return sizes


def compute_hist_lengths(cfg: dict) -> list[int]:
    """Compute per-table history lengths matching C++ logic."""
    n = cfg["tc_num_tables"]
    minh = cfg["tc_minh"]
    maxh = cfg["tc_maxh"]
    hist = cfg.get("tc_hist", "GEOMETRIC")

    if hist == "GEOMETRIC":
        if n < 2:
            return [maxh]
        h = [0] * n
        prev = 0
        ratio = maxh / minh
        for i in range(n):
            e = i / (n - 1)
            hl = int(minh * (ratio ** e))
            hl = max(hl, prev + 1)
            h[n - 1 - i] = hl
            prev = hl
        return h
    elif hist == "QUADRATIC":
        d, k = 2, 1
        h = [0] * n
        h[0] = minh
        for i in range(1, n):
            h[i] = h[i - 1] + d * i + k
        raw_max, raw_min = h[-1], h[0]
        prev = 0
        for i in range(n):
            t = (h[i] - raw_min) / (raw_max - raw_min) if raw_max > raw_min else 0
            scaled = minh + int(t * (maxh - minh))
            scaled = max(scaled, prev + 1) if i > 0 else scaled
            h[i] = scaled
            prev = scaled
        h.reverse()
        return h
    elif hist == "ROS":
        d, k, f, m, t_split = 2, 1, 0.1, 1.1, 15
        h = [0] * n
        h[0] = minh
        for i in range(1, n):
            if i < t_split:
                h[i] = h[i - 1] + d * i + k
            else:
                mult = f * (i - t_split + 1) + m
                nxt = int(h[i - 1] * mult + 0.5)
                h[i] = max(nxt, h[i - 1] + 1)
        raw_max, raw_min = h[-1], h[0]
        prev = 0
        for i in range(n):
            t = (h[i] - raw_min) / (raw_max - raw_min) if raw_max > raw_min else 0
            scaled = minh + int(t * (maxh - minh))
            scaled = max(scaled, prev + 1) if i > 0 else scaled
            h[i] = scaled
            prev = scaled
        h.reverse()
        return h
    elif hist == "SUPEREXP":
        f, m_val = 0.1, 1.1
        h = [0] * n
        h[0] = minh
        for i in range(1, n):
            mult = f * i + m_val
            nxt = int(h[i - 1] * mult + 0.5)
            h[i] = max(nxt, h[i - 1] + 1)
        h.reverse()
        return h
    else:
        return list(range(minh, minh + n))


def compute_tag_widths(cfg: dict) -> list[int]:
    """Compute per-table tag widths matching C++ TagFn logic."""
    n = cfg["tc_num_tables"]
    tag = cfg["tc_tag"]
    tag_fn = cfg.get("tc_tag_fn", "__default__")

    if tag_fn == "__default__":
        # GradedTag<TAG, TAG-1>
        hi, lo = tag, max(1, tag - 1)
        if n <= 1:
            return [hi]
        return [hi - (hi - lo) * i // (n - 1) for i in range(n)]
    elif tag_fn == "LogTag":
        # LogTag<BASE, SCALE>: BASE + SCALE*(n-1-i)/4
        scale = max(1, tag // 4)
        return [tag + scale * (n - 1 - i) // 4 for i in range(n)]
    elif tag_fn == "StepTag":
        split = max(1, n // 2)
        lo = max(1, tag - 2)
        return [tag if i < split else lo for i in range(n)]
    else:
        return [tag] * n


def format_arrays(cfg: dict) -> str:
    """Format table sizes, history lengths, and tag widths for display."""
    sizes = compute_table_sizes(cfg)
    hists = compute_hist_lengths(cfg)
    tags = compute_tag_widths(cfg)
    return (f"  sizes={sizes}\n"
            f"  hists={hists}\n"
            f"  tags ={tags}")


def config_to_predictor_string(cfg: dict) -> str:
    """Convert a config dict to a TageAhead<...> template string."""
    # Build TATableConfig
    tag = cfg["tc_tag"]
    sr = cfg["tc_size_ratio"]
    size = cfg["tc_size"]

    tag_fn = cfg.get("tc_tag_fn", "__default__")
    if tag_fn == "__default__":
        tag_fn = f"ta::GradedTag<{tag},{max(1, tag-1)}>"
    elif tag_fn == "UniformTag":
        tag_fn = f"ta::UniformTag<{tag}>"
    elif tag_fn == "LogTag":
        # LogTag<BASE, SCALE>: BASE + SCALE*level/4, level = n-1-i
        # Use SCALE = tag//4 so range spans ~TAG to ~1.25*TAG
        scale = max(1, tag // 4)
        tag_fn = f"ta::LogTag<{tag},{scale}>"
    elif tag_fn == "StepTag":
        # StepTag<TAG_HI, TAG_LO, SPLIT>: first SPLIT tables get TAG_HI, rest get TAG_LO
        # Split at midpoint; hi=TAG, lo=max(1, TAG-2)
        n = cfg["tc_num_tables"]
        split = max(1, n // 2)
        lo = max(1, tag - 2)
        tag_fn = f"ta::StepTag<{tag},{lo},{split}>"

    size_fn = cfg.get("tc_size_fn", "__default__")
    if size_fn == "__default__" or size_fn == "GeoSize":
        size_fn = f"ta::GeoSize<{size},{sr}>"
    elif size_fn == "InvGeoSize":
        size_fn = f"ta::InvGeoSize<{size},{sr}>"
    elif size_fn == "UniformSize":
        size_fn = f"ta::UniformSize<{size}>"
    elif size_fn == "SqrtHistSize":
        size_fn = f"ta::SqrtHistSize<{size}>"
    # else: use as-is (raw typename string)

    hist = cfg["tc_hist"]
    tc = (f"TATableConfig<{cfg['tc_num_tables']},{size},{tag},"
          f"{cfg['tc_minh']},{cfg['tc_maxh']},{sr},"
          f"ta::HistSeries::{hist},{tag_fn},{size_fn}>")

    # Build rest of TageAhead params
    hist_mode = cfg["HIST_MODE"]
    decay_miss = cfg["DECAY_MISS"]
    decay_op = cfg["DECAY_OP"]
    lfsr_w = cfg["DECAY_LFSR_WIDTH"]
    nt = cfg["tc_num_tables"]

    parts = [
        tc,
        str(cfg["N"]),
        str(cfg["PATHBITS"]),
        str(cfg["SEC_TAG_BITS"]),
        bool_str(cfg["USE_SEC_TAG"]),
        str(cfg["NUM_PATHS"]),
        cfg["SecTagHashFn"],
        str(cfg["CTR_WIDTH"]),
        str(cfg["HYST_WIDTH"]),
        str(cfg["U_WIDTH"]),
        str(cfg["FB_CAPACITY"]),
        bool_str(cfg["USE_GSHARE"]),
        str(cfg["GS_HIST"]),
        str(cfg["META_WIDTH"]),
        str(cfg["META_CAPACITY"]),
        str(cfg["META_PIPE"]),
        str(cfg["LINEINST"]),
        bool_str(cfg["SHARED_HYS"]),
        f"HistUpdate::{hist_mode}",
        cfg["AllocCfg"],
        str(cfg["ACC_WIDTH"]),
        str(cfg["ALLOC_WIDTH"]),
        bool_str(cfg["DECAY_ENABLE"]),
        f"DecayMiss::{decay_miss}",
        f"DecayOp::{decay_op}",
        f"ta::uniform_array<u64,{nt}>({lfsr_w})",
        cfg["DecayThreshFn"],
        bool_str(cfg["EPOCH_ENABLE"]),
        cfg["EpochTriggerFn"],
        bool_str(cfg["EPOCH_RESET_ACC"]),
        bool_str(cfg["EPOCH_RESET_ALLOC"]),
        bool_str(cfg["FB_RECONCILE"]),
        str(cfg["FARALLOC_DIST"]),
    ]

    return f"TageAhead<{','.join(parts)}>"


def config_id(cfg: dict) -> str:
    pred = config_to_predictor_string(cfg)
    return hashlib.md5(pred.encode()).hexdigest()[:10]


# ============================================================================
# Sweep config generation
# ============================================================================

DEFAULT_SWEEP = {
    "base": dict(DEFAULTS),
    "sweep": {
        # Tier 1: Table geometry
        "tc_num_tables": [5, 7, 9],
        "tc_size": [512, 1024, 2048],
        "tc_minh": [4, 10, 20],
        "tc_maxh": [60, 100, 200],
    },
    "timing_traces": TIMING_TRACES,
    "timing_warmup": TIMING_WARMUP,
    "timing_measure": TIMING_MEASURE,
    "p2_max_cycles": 1.0,
}


# Size functors that use tc_size_ratio; all others ignore it
RATIO_USING_SIZE_FNS = {"__default__", "GeoSize", "InvGeoSize"}

def generate_configs(sweep_cfg: dict) -> list[dict]:
    """Generate cross-product of sweep dimensions over base config,
    filtering redundant combos (e.g. ratio>1 for functors that ignore it)."""
    base = dict(sweep_cfg["base"])
    sweep_dims = sweep_cfg["sweep"]

    if not sweep_dims:
        return [base]

    keys = list(sweep_dims.keys())
    value_lists = [sweep_dims[k] for k in keys]

    seen_ids = set()
    configs = []
    for combo in itertools.product(*value_lists):
        cfg = dict(base)
        for k, v in zip(keys, combo):
            cfg[k] = v

        # Normalize: functors that don't use ratio → canonicalize ratio to 1
        if cfg.get("tc_size_fn") not in RATIO_USING_SIZE_FNS:
            cfg["tc_size_ratio"] = 1

        cid = config_id(cfg)
        if cid in seen_ids:
            continue
        seen_ids.add(cid)
        configs.append(cfg)

    return configs


# ============================================================================
# Build & Run
# ============================================================================

def build_config(cfg: dict) -> tuple[dict, Path, float, str | None]:
    """Compile a single config. Returns (cfg, binary, build_time, error)."""
    cid = config_id(cfg)
    binary = BUILD_DIR / f"cbp_{cid}"

    if binary.exists():
        return (cfg, binary, 0.0, None)

    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    pred = config_to_predictor_string(cfg)
    cmd = (f"{CXX} {COMMON_FLAGS} {WARN_FLAGS} "
           f"-Itrace_files -o {binary} cbp.cpp -lz "
           f"-DPREDICTOR='{pred}'")

    t0 = time.monotonic()
    try:
        r = subprocess.run(cmd, shell=True, check=True,
                           capture_output=True, text=True, timeout=600,
                           cwd=str(PROJECT_DIR))
        return (cfg, binary, time.monotonic() - t0, None)
    except subprocess.CalledProcessError as e:
        return (cfg, binary, time.monotonic() - t0,
                e.stderr[-800:] if e.stderr else "compile error")
    except subprocess.TimeoutExpired:
        return (cfg, binary, 600.0, "COMPILE_TIMEOUT")


def run_trace(binary: Path, trace: Path, warmup: int, measure: int
              ) -> tuple[dict | None, str | None]:
    """Run binary on trace. Returns (parsed_data, error)."""
    trace_name = trace.stem.replace("_trace", "")
    cmd = [str(binary), str(trace), trace_name, str(warmup), str(measure)]

    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=1200, check=True)
    except subprocess.CalledProcessError as e:
        return (None, f"run failed: {e.stderr[:300]}")
    except subprocess.TimeoutExpired:
        return (None, "RUN_TIMEOUT")

    line = r.stdout.strip()
    if not line:
        return (None, "empty output")

    try:
        parts = line.split(",")
        return ({
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
        }, None)
    except (IndexError, ValueError) as e:
        return (None, f"parse error: {e}")


def timing_check(binary: Path, sweep_cfg: dict) -> tuple[float, list[dict]]:
    """Run timing gate on multiple traces.
    Returns (worst_p2_latency, list of run results)."""
    traces = sweep_cfg.get("timing_traces", TIMING_TRACES)
    warmup = sweep_cfg.get("timing_warmup", TIMING_WARMUP)
    measure = sweep_cfg.get("timing_measure", TIMING_MEASURE)

    worst_p2 = 0.0
    results = []
    good_runs = 0
    for t in traces:
        tp = TRACE_DIR / t
        if not tp.exists():
            continue
        data, err = run_trace(binary, tp, warmup, measure)
        if err:
            results.append({"trace": t, "error": err})
            continue
        worst_p2 = max(worst_p2, data["p2_latency"])
        results.append(data)
        good_runs += 1
    # If every run crashed, treat as timing failure with sentinel P2
    if good_runs == 0 and results:
        worst_p2 = 999.0
    return worst_p2, results


# ============================================================================
# VFS computation (mirrors predictor_metrics.py + vfs.py)
# ============================================================================

def compute_vfs(runs: list[dict]) -> dict:
    """Compute aggregate VFS from eval runs."""
    p1_lat = 0
    p2_lat = 0
    for r in runs:
        p1_lat = max(p1_lat, math.ceil(r["p1_latency"]))
        p2_lat = max(p2_lat, math.ceil(r["p2_latency"]))

    n = 0
    sum_inv_ipc = 0.0
    sum_cpi = 0.0
    sum_epi = 0.0
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
        mpki = mpi * 1000

        if p2_lat <= p1_lat:
            cycles = preds * max(1, p2_lat)
        else:
            cycles = (preds * max(1, p1_lat) +
                      diverge * p2_lat -
                      divend * max(1, p1_lat))
        cycles += extra
        ipc = instr / cycles if cycles > 0 else 0.001
        cpi = mpi * (MISP_PENALTY + p2_lat)

        n += 1
        sum_inv_ipc += 1.0 / ipc
        sum_cpi += cpi
        sum_epi += epi
        sum_mpki += mpki

    if n == 0:
        return {"vfs": 0, "mpki": 0, "epi": 0, "ipc": 0,
                "p1_lat": p1_lat, "p2_lat": p2_lat, "n_traces": 0}

    avg_ipc = n / sum_inv_ipc
    avg_cpi = sum_cpi / n
    avg_epi = sum_epi / n
    avg_mpki = sum_mpki / n

    IPCcbp0 = 8
    CPIcbp0 = 0.0315
    EPIcbp0 = 1000
    ALPHA = 1.625
    BETA = 4 * ALPHA / (ALPHA - 1) ** 2
    GAMMA = 2 / (ALPHA - 1)
    cbp_energy_ratio = 0.05

    WPI0 = IPCcbp0 * CPIcbp0
    WPI = avg_ipc * avg_cpi
    speedup = (avg_ipc / IPCcbp0) * (1 + WPI0) / (1 + WPI)
    LAMBDA = 1 / (1 + WPI0 / 2) - cbp_energy_ratio
    nEPI = ((avg_epi / EPIcbp0) * cbp_energy_ratio +
            LAMBDA * speedup ** GAMMA) * (1 + WPI / 2)
    vfs = speedup * ALPHA * (1 - 2 / (1 + math.sqrt(1 + BETA / (speedup * nEPI))))

    return {
        "vfs": round(vfs, 6), "mpki": round(avg_mpki, 3),
        "epi": round(avg_epi, 1), "ipc": round(avg_ipc, 4),
        "p1_lat": p1_lat, "p2_lat": p2_lat, "n_traces": n,
    }


# ============================================================================
# Persistence
# ============================================================================

def results_file() -> Path:
    return RESULTS_DIR / "results.json"


def load_results() -> dict:
    rf = results_file()
    if rf.exists():
        return json.loads(rf.read_text())
    return {"timing_passed": {}, "timing_failed": {}, "eval": {}}


def save_results(data: dict):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    results_file().write_text(json.dumps(data, indent=2))


def cfg_summary(cfg: dict, sweep_dims: dict) -> str:
    """Short summary showing only swept params."""
    parts = []
    for k in sorted(sweep_dims.keys()):
        parts.append(f"{k}={cfg[k]}")
    return " ".join(parts)


# ============================================================================
# Commands
# ============================================================================

def cmd_init(args):
    """Create default sweep_config.json."""
    if CONFIG_FILE.exists() and not args.force:
        print(f"{CONFIG_FILE} already exists. Use --force to overwrite.")
        return
    CONFIG_FILE.write_text(json.dumps(DEFAULT_SWEEP, indent=2))
    print(f"Created {CONFIG_FILE}")
    print(f"Edit sweep dimensions in 'sweep' key, then run --run")


def cmd_run(args):
    """Stream build + timing gate, 2 configs at a time."""
    sweep_cfg = json.loads(CONFIG_FILE.read_text())
    configs = generate_configs(sweep_cfg)
    p2_max = sweep_cfg.get("p2_max_cycles", 1.0)
    results = load_results()
    sweep_dims = sweep_cfg.get("sweep", {})

    # Filter out already-processed configs
    todo = []
    for cfg in configs:
        cid = config_id(cfg)
        if cid in results["timing_passed"] or cid in results["timing_failed"]:
            continue
        todo.append(cfg)

    total = len(configs)
    done = total - len(todo)
    print(f"Total configs: {total}, already done: {done}, remaining: {len(todo)}")
    print(f"Timing gate: P2 ≤ {p2_max} on {len(sweep_cfg.get('timing_traces', TIMING_TRACES))} traces")
    print(f"Building {args.jobs} at a time\n")

    passed = 0
    failed = 0

    # Process in batches of args.jobs
    for batch_start in range(0, len(todo), args.jobs):
        batch = todo[batch_start:batch_start + args.jobs]
        idx_base = done + batch_start

        # Phase 1: Build in parallel
        build_results = []
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futs = {pool.submit(build_config, cfg): cfg for cfg in batch}
            for fut in as_completed(futs):
                build_results.append(fut.result())

        # Phase 2: Timing check sequentially per config
        for cfg, binary, bt, build_err in build_results:
            cid = config_id(cfg)
            idx = idx_base + batch.index(cfg)
            summary = cfg_summary(cfg, sweep_dims)

            if build_err:
                print(f"{DIM}[{idx+1}/{total}]{RESET} {RED}FAIL compile{RESET}  {summary}")
                print(f"         {DIM}{build_err[:200]}{RESET}")
                results["timing_failed"][cid] = {
                    "config": cfg, "error": f"compile: {build_err[:200]}"}
                failed += 1
                save_results(results)
                continue

            bt_str = f"{bt:.1f}s" if bt > 0 else "cached"
            print(f"{DIM}[{idx+1}/{total}]{RESET} built ({bt_str})  {CYAN}{summary}{RESET}")
            print(f"{DIM}{format_arrays(cfg)}{RESET}")

            # Timing gate
            worst_p2, trace_runs = timing_check(binary, sweep_cfg)
            p2_ceil = math.ceil(worst_p2) if worst_p2 > 0 else 999

            if worst_p2 <= p2_max:
                print(f"         {GREEN}PASS{RESET}  P2={worst_p2:.4f} (ceil={p2_ceil})")
                results["timing_passed"][cid] = {
                    "config": cfg,
                    "worst_p2": worst_p2,
                    "timing_runs": trace_runs,
                }
                passed += 1
            else:
                print(f"         {RED}FAIL{RESET}  P2={worst_p2:.4f} (ceil={p2_ceil}) > {p2_max}")
                results["timing_failed"][cid] = {
                    "config": cfg,
                    "worst_p2": worst_p2,
                    "timing_runs": trace_runs,
                }
                failed += 1
                # Clean up binary to save disk
                if binary.exists():
                    binary.unlink()

            save_results(results)

    print(f"\nDone. {GREEN}Passed: {passed}{RESET}, {RED}Failed: {failed}{RESET}")
    print(f"Total passing: {BOLD}{len(results['timing_passed'])}{RESET}")


def cmd_eval(args):
    """Full 20-trace eval on timing-passed configs."""
    results = load_results()
    passed = results["timing_passed"]

    # Filter out already-evaluated
    todo = {cid: info for cid, info in passed.items()
            if cid not in results["eval"]}

    print(f"Timing-passed: {len(passed)}, already evaluated: {len(results['eval'])}, "
          f"remaining: {len(todo)}")

    available_traces = [t for t in EVAL_TRACES if (TRACE_DIR / t).exists()]
    print(f"Eval traces available: {len(available_traces)}/{len(EVAL_TRACES)}\n")

    for i, (cid, info) in enumerate(todo.items()):
        cfg = info["config"]
        binary = BUILD_DIR / f"cbp_{cid}"

        if not binary.exists():
            # Rebuild
            cfg_r, binary, bt, err = build_config(cfg)
            if err:
                print(f"[{i+1}/{len(todo)}] {cid} rebuild failed: {err[:200]}")
                continue

        sweep_dims = json.loads(CONFIG_FILE.read_text()).get("sweep", {})
        summary = cfg_summary(cfg, sweep_dims)
        print(f"[{i+1}/{len(todo)}] eval {cid}  {summary}")
        print(format_arrays(cfg))

        runs = []
        for j, t in enumerate(available_traces):
            tp = TRACE_DIR / t
            data, err = run_trace(binary, tp, EVAL_WARMUP, EVAL_MEASURE)
            if err:
                print(f"  trace {j+1}/{len(available_traces)} {t}: ERROR {err}")
                continue
            mpki = data["mispredictions"] / data["instructions"] * 1000
            print(f"  trace {j+1}/{len(available_traces)} {t}: "
                  f"MPKI={mpki:.2f} P2={data['p2_latency']:.4f} "
                  f"EPI={data['epi']:.0f}")
            runs.append(data)

        if runs:
            metrics = compute_vfs(runs)
            print(f"  → VFS={metrics['vfs']:.4f}  MPKI={metrics['mpki']:.2f}  "
                  f"EPI={metrics['epi']:.0f}  IPC={metrics['ipc']:.4f}")
            results["eval"][cid] = {
                "config": cfg,
                "metrics": metrics,
                "runs": runs,
            }
        else:
            print(f"  → no successful runs")

        save_results(results)

    print(f"\nTotal evaluated: {len(results['eval'])}")


def cmd_report(args):
    """Print results ranked by VFS."""
    results = load_results()
    sweep_cfg = json.loads(CONFIG_FILE.read_text()) if CONFIG_FILE.exists() else {}
    sweep_dims = sweep_cfg.get("sweep", {})

    if results["eval"]:
        entries = []
        for cid, info in results["eval"].items():
            m = info["metrics"]
            entries.append((m["vfs"], cid, info["config"], m))

        entries.sort(reverse=True)
        top = args.top if args.top else len(entries)

        print(f"{'#':>3}  {'VFS':>8}  {'MPKI':>7}  {'EPI':>6}  {'IPC':>7}  "
              f"{'P2':>4}  {'ID':>10}  Config")
        print("-" * 100)
        for i, (vfs, cid, cfg, m) in enumerate(entries[:top]):
            summary = cfg_summary(cfg, sweep_dims)
            print(f"{i+1:3}  {vfs:8.4f}  {m['mpki']:7.2f}  {m['epi']:6.0f}  "
                  f"{m['ipc']:7.4f}  {m['p2_lat']:4}  {cid:>10}  {summary}")
    else:
        print("No eval results yet.")

    tp = len(results["timing_passed"])
    tf = len(results["timing_failed"])
    ev = len(results["eval"])
    print(f"\nSummary: {tp} timing-passed, {tf} timing-failed, {ev} evaluated")


def cmd_status(args):
    """Show sweep status."""
    if not CONFIG_FILE.exists():
        print("No sweep_config.json. Run --init first.")
        return

    sweep_cfg = json.loads(CONFIG_FILE.read_text())
    configs = generate_configs(sweep_cfg)
    results = load_results()

    sweep_dims = sweep_cfg.get("sweep", {})
    dim_str = ", ".join(f"{k}({len(v)})" for k, v in sweep_dims.items())

    print(f"Sweep dimensions: {dim_str}")
    print(f"Total configs: {len(configs)}")
    print(f"Timing passed: {len(results['timing_passed'])}")
    print(f"Timing failed: {len(results['timing_failed'])}")
    print(f"Evaluated: {len(results['eval'])}")


def cmd_show(args):
    """Show the predictor string for a config ID or the default config."""
    if args.config_id:
        results = load_results()
        for bucket in [results["timing_passed"], results["timing_failed"], results["eval"]]:
            if args.config_id in bucket:
                cfg = bucket[args.config_id]["config"]
                print(config_to_predictor_string(cfg))
                return
        print(f"Config ID {args.config_id} not found in results.")
    else:
        print(config_to_predictor_string(DEFAULTS))


# ============================================================================
# Main
# ============================================================================

def main():
    p = argparse.ArgumentParser(
        description="TageAhead parameter sweep (streaming build + timing gate)")
    sub = p.add_subparsers(dest="cmd")

    sp_init = sub.add_parser("init", help="Create default sweep_config.json")
    sp_init.add_argument("--force", action="store_true")

    sp_run = sub.add_parser("run", help="Build + timing gate (streaming)")
    sp_run.add_argument("-j", "--jobs", type=int, default=2,
                        help="Configs to build in parallel (default: 2)")

    sp_eval = sub.add_parser("eval", help="Full eval on timing-passed configs")
    sp_eval.add_argument("-j", "--jobs", type=int, default=1)

    sp_report = sub.add_parser("report", help="Print results ranked by VFS")
    sp_report.add_argument("--top", type=int, default=0)

    sp_status = sub.add_parser("status", help="Show sweep status")

    sp_show = sub.add_parser("show", help="Show predictor string for config")
    sp_show.add_argument("config_id", nargs="?", default=None)

    args = p.parse_args()
    if not args.cmd:
        p.print_help()
        return

    {"init": cmd_init, "run": cmd_run, "eval": cmd_eval,
     "report": cmd_report, "status": cmd_status, "show": cmd_show
     }[args.cmd](args)


if __name__ == "__main__":
    main()
