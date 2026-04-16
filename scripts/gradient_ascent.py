#!/usr/bin/env python3
"""
Gradient ascent optimizer for TageDirect branch predictor parameters.

Minimizes MPKI under an EPI budget constraint. Each iteration perturbs
every parameter by ±step, evaluates all neighbors in parallel, and moves
to the best feasible neighbor. Perturbation halves on plateau.

Usage:
    python3 scripts/gradient_ascent.py --config configs/gradient_default.yaml \
        --trace-dir traces/ -j 8 --perturbation 2

    python3 scripts/gradient_ascent.py --report out/gradient/results.csv --top 20
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import yaml
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field, fields
from datetime import datetime
from pathlib import Path
from typing import Any

# Add scripts/ to path for sweep_common imports
sys.path.insert(0, str(Path(__file__).parent))
from sweep_common import (
    run_trace, compute_metrics, resolve_traces, REPR_TRACES,
    append_csv_row, load_csv_rows, CXX, COMMON_FLAGS, WARN_FLAGS,
    SweepConfig, CONFIG_FIELDS as SWEEP_CONFIG_FIELDS,
    build_config as build_sweep_config,
)


# ============================================================================
# TDConfig: maps YAML params -> TageDirect<...> predictor string
# ============================================================================

# ---- C++ enum/functor mappings ----

HIST_SERIES_MAP = {
    "geometric":  "td::HistSeries::GEOMETRIC",
    "quadratic":  "td::HistSeries::QUADRATIC",
    "superexp":   "td::HistSeries::SUPEREXP",
    "ros":        "td::HistSeries::ROS",
}

# TagFn: map name -> lambda(tag) returning C++ type string
TAG_FN_MAP = {
    "graded":  lambda tag: f"td::GradedTag<{tag},{max(tag-3, 4)}>",
    "uniform": lambda tag: f"td::UniformTag<{tag}>",
    "step":    lambda tag: f"td::StepTag<{tag},{max(tag-3, 4)},3>",
    "log":     lambda tag: f"td::LogTag<{max(tag-3, 4)},1>",
}

# SizeFn: map name -> lambda(size, size_ratio) returning C++ type string
SIZE_FN_MAP = {
    "geo":      lambda sz, sr: f"td::GeoSize<{sz},{sr}>",
    "uniform":  lambda sz, sr: f"td::UniformSize<{sz}>",
    "invgeo":   lambda sz, sr: f"td::InvGeoSize<{sz},{sr}>",
    "step":     lambda sz, sr: f"td::StepSize<{sz},{max(sz//2,64)},3>",
    "sqrthist": lambda sz, sr: f"td::SqrtHistSize<{sz}>",
}

DECAY_POLICY_MAP = {
    "mild":       "td::TDDecayMild",
    "aggressive": "td::TDDecayAggressive",
}

FOLD_FN_MAP = {
    "xor":       "td::XORFold",
    "rotatexor": "td::RotateXORFold",
}

ALLOC_TRIGGER_MAP = {
    "mispredict": "td::AllocTrigger::MISPREDICT",
    "base_wrong": "td::AllocTrigger::BASE_WRONG",
    "disagree":   "td::AllocTrigger::DISAGREE",
    "tage_miss":  "td::AllocTrigger::TAGE_MISS",
    "tage_wrong": "td::AllocTrigger::TAGE_WRONG",
}

UCTR_POLICY_MAP = {
    "always_inc": "td::UctrPolicy::ALWAYS_INC",
    "faralloc":   "td::UctrPolicy::FARALLOC",
    "noalloc":    "td::UctrPolicy::NOALLOC",
    "penalty_na": "td::UctrPolicy::PENALTY_NA",
}


@dataclass(frozen=True)
class TDConfig:
    # Numeric params
    num_tables: int = 4
    size: int = 2048
    tag: int = 11
    ctr: int = 1
    hyst: int = 2
    minh: int = 11
    maxh: int = 48
    size_ratio: int = 1
    p1_table_size: int = 4096
    p1_hist: int = 6
    metabits: int = 4
    metapipe: int = 2
    decay_ctr: int = 0
    epoch_ctr_bits: int = 8
    # Boolean params
    shared_hys: bool = False
    conf_gate: bool = False
    max_alloc: int = 1
    # Functor/enum params (string keys into maps above)
    hist_series: str = "geometric"
    tag_fn: str = "graded"
    size_fn: str = "geo"
    decay_policy: str = "mild"
    fold_fn: str = "xor"
    alloc_trigger: str = "mispredict"
    uctr_policy: str = "always_inc"

    @property
    def config_id(self) -> str:
        return hashlib.md5(self.predictor_string.encode()).hexdigest()[:8]

    @property
    def alloc_config_str(self) -> str:
        # Build inline struct when non-default alloc params are used
        parts = []
        if self.conf_gate:
            parts.append("static constexpr bool CONF_GATE = true;")
        if self.max_alloc >= 2:
            parts.append(f"static constexpr u64 MAX_ALLOC = {self.max_alloc};")
        if self.alloc_trigger != "mispredict":
            parts.append(f"static constexpr td::AllocTrigger ALLOC_TRIGGER = "
                         f"{ALLOC_TRIGGER_MAP[self.alloc_trigger]};")
        if self.uctr_policy != "always_inc":
            parts.append(f"static constexpr td::UctrPolicy UCTR_POLICY = "
                         f"{UCTR_POLICY_MAP[self.uctr_policy]};")

        if not parts:
            return "td::TDDefaultAllocConfig"

        # Use named configs for common combos, inline struct for custom
        if self.conf_gate and self.max_alloc >= 2 and self.alloc_trigger == "mispredict" and self.uctr_policy == "always_inc":
            return "td::TDAllocConfGate2"
        elif self.conf_gate and self.max_alloc == 1 and self.alloc_trigger == "mispredict" and self.uctr_policy == "always_inc":
            return "td::TDAllocConfGate"
        elif self.max_alloc >= 2 and not self.conf_gate and self.alloc_trigger == "mispredict" and self.uctr_policy == "always_inc":
            return "td::TDAlloc2Press"

        # Fallback: generate a struct that inherits from TDDefaultAllocConfig
        # C++ doesn't support inline struct in template args, so we use
        # the named configs and warn if we can't match
        # For now, pick the closest named config
        if self.conf_gate:
            return "td::TDAllocConfGate"
        if self.max_alloc >= 2:
            return "td::TDAlloc2Press"
        return "td::TDDefaultAllocConfig"

    @property
    def table_config_str(self) -> str:
        hist_cpp = HIST_SERIES_MAP[self.hist_series]
        tag_fn_cpp = TAG_FN_MAP[self.tag_fn](self.tag)
        size_fn_cpp = SIZE_FN_MAP[self.size_fn](self.size, self.size_ratio)
        return (f"td::TDTableConfig<{self.num_tables},{self.size},"
                f"{self.tag},{self.ctr},{self.hyst},1,"
                f"{self.minh},{self.maxh},{self.size_ratio},"
                f"{hist_cpp},{tag_fn_cpp},{size_fn_cpp}>")

    @property
    def predictor_string(self) -> str:
        b = lambda v: "true" if v else "false"
        decay_cpp = DECAY_POLICY_MAP[self.decay_policy]
        fold_cpp = FOLD_FN_MAP[self.fold_fn]
        parts = [
            self.table_config_str,
            self.alloc_config_str,
            "1024",                     # LINEINST
            "7",                        # N (max branches per block)
            str(self.decay_ctr),        # DECAY_CTR
            "2",                        # DECAY_GRAN
            decay_cpp,                  # DecayPolicy
            "true",                     # P1_USE_GSHARE
            str(self.p1_table_size),    # P1_TABLE_SIZE
            str(self.p1_hist),          # P1_HIST
            "true",                     # USE_META
            str(self.metabits),         # METABITS
            str(self.metapipe),         # METAPIPE
            "false",                    # USE_PATH_HIST
            "27",                       # PATH_HIST_WIDTH
            "6",                        # PATH_BITS
            fold_cpp,                   # FoldFn
            "4",                        # RWRAM_BANKS
            "0",                        # RWRAM_BANK_SHIFT
            str(self.epoch_ctr_bits),   # EPOCH_CTR_BITS
            b(self.shared_hys),         # SHARED_HYS
            "false",                    # USE_DIR_HIST
        ]
        return f"TageDirect<{','.join(parts)}>"

    def replace(self, **kwargs) -> 'TDConfig':
        d = {f.name: getattr(self, f.name) for f in fields(self)}
        d.update(kwargs)
        return TDConfig(**d)

    def to_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in fields(self)}


TD_CONFIG_FIELDS = [f.name for f in fields(TDConfig)]


# ============================================================================
# Parameter space: neighbor generation
# ============================================================================

class ParamSpace:
    """Generates ±step neighbors for a parameter given its spec."""

    def __init__(self, name: str, spec: dict):
        self.name = name
        if "values" in spec:
            self.mode = "discrete"
            self.values = spec["values"]
        else:
            self.mode = "range"
            self.min_val = spec["min"]
            self.max_val = spec["max"]
            self.step = spec["step"]

    def neighbors(self, current_val, perturbation: int = 1) -> list:
        """Return list of neighbor values (±step*perturbation from current)."""
        results = []
        if self.mode == "discrete":
            try:
                idx = self.values.index(current_val)
            except ValueError:
                return results
            for delta in [-perturbation, perturbation]:
                ni = idx + delta
                if 0 <= ni < len(self.values) and ni != idx:
                    results.append(self.values[ni])
        else:
            for delta in [-perturbation, perturbation]:
                nv = current_val + delta * self.step
                if self.min_val <= nv <= self.max_val and nv != current_val:
                    results.append(nv)
        return results


# ============================================================================
# Build helpers
# ============================================================================

def build_td_config(cfg: TDConfig, build_dir: Path, extra_flags: str = "",
                    monitor: bool = False) -> tuple:
    """Build a TageDirect binary. Returns (cfg, binary_path, build_time, error)."""
    prefix = "mon" if monitor else "cbp"
    binary = build_dir / f"{prefix}_{cfg.config_id}"
    if binary.exists():
        return (cfg, binary, 0.0, None)

    build_dir.mkdir(parents=True, exist_ok=True)
    pred = cfg.predictor_string
    monitor_flags = "-DTAGE_MONITOR -DCHEATING_MODE -DREAD_WRITE_RAM" if monitor else ""
    cmd = (f'{CXX} {COMMON_FLAGS} {WARN_FLAGS} {extra_flags} {monitor_flags} '
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


def run_monitor(binary: Path, trace: Path, warmup: int, measure: int,
                output_path: Path) -> str | None:
    """Run monitor binary, capture stderr to file. Returns error or None."""
    trace_name = trace.stem.replace("_trace", "")
    cmd = [str(binary), str(trace), trace_name, str(warmup), str(measure)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(result.stderr)
        return None
    except subprocess.TimeoutExpired:
        return "TIMEOUT"
    except Exception as e:
        return str(e)[:200]


# ============================================================================
# CSV schema
# ============================================================================

GRADIENT_BASE_FIELDS = [
    "timestamp", "iteration", "config_id", "varied_param", "is_improvement",
    "mpki", "epi", "ipc", "vfs", "p1_lat", "p2_lat", "n_traces",
    "build_time_s", "eval_time_s",
]

# Default (TageDirect) field order; Tage uses its own
GRADIENT_FIELDS = GRADIENT_BASE_FIELDS + TD_CONFIG_FIELDS + ["predictor_string"]
GRADIENT_FIELDS_TAGE = GRADIENT_BASE_FIELDS + SWEEP_CONFIG_FIELDS + ["predictor_string"]


def _gradient_fields(predictor: str = "tagedirect") -> list[str]:
    if predictor == "tage":
        return GRADIENT_FIELDS_TAGE
    return GRADIENT_FIELDS


# ============================================================================
# Eval cache (reuses coordinate_descent pattern)
# ============================================================================

class EvalCache:
    def __init__(self, csv_path: Path):
        self.csv_path = csv_path
        self.cache: dict[str, dict] = {}
        self._load()

    def _load(self):
        for r in load_csv_rows(self.csv_path):
            cid = r.get("config_id", "")
            if cid:
                try:
                    self.cache[cid] = {
                        "mpki": float(r["mpki"]),
                        "epi": float(r["epi"]),
                        "ipc": float(r["ipc"]),
                        "vfs": float(r["vfs"]),
                        "p1_lat": r["p1_lat"],
                        "p2_lat": r["p2_lat"],
                        "n_traces": r["n_traces"],
                    }
                except (ValueError, KeyError):
                    pass

    def get(self, cid: str) -> dict | None:
        return self.cache.get(cid)

    def put(self, cid: str, metrics: dict):
        self.cache[cid] = metrics


# ============================================================================
# Checkpoint
# ============================================================================

def save_checkpoint(path: Path, iteration: int, config: TDConfig,
                    mpki: float, perturbation: int, total_evals: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump({
            "iteration": iteration,
            "current_config": config.to_dict(),
            "current_mpki": mpki,
            "perturbation": perturbation,
            "total_evals": total_evals,
        }, f, indent=2)


def load_checkpoint(path: Path) -> dict | None:
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def save_best_iteration(saved_dir: Path, iteration: int, cfg, metrics: dict):
    """Save params and metrics for a new best iteration to saved/."""
    saved_dir.mkdir(parents=True, exist_ok=True)
    entry = {
        "iteration": iteration,
        "config_id": cfg.config_id,
        "predictor_string": cfg.predictor_string,
        "metrics": {k: float(v) if isinstance(v, (int, float)) else v
                    for k, v in metrics.items()},
        "params": cfg.to_dict(),
    }
    path = saved_dir / f"iter_{iteration:03d}_{cfg.config_id}.json"
    with open(path, "w") as f:
        json.dump(entry, f, indent=2)
    print(f"  Saved best to {path}", file=sys.stderr)


# ============================================================================
# Core evaluation
# ============================================================================

def _build_config(cfg, build_dir: Path, extra_flags: str = "",
                  monitor: bool = False, predictor: str = "tagedirect"):
    """Dispatch to the correct build function based on predictor type."""
    if predictor == "tage":
        if monitor:
            # build_sweep_config doesn't support monitor; build manually
            return _build_tage_monitor(cfg, build_dir, extra_flags)
        return build_sweep_config(cfg, build_dir, extra_flags)
    return build_td_config(cfg, build_dir, extra_flags, monitor=monitor)


def _build_tage_monitor(cfg, build_dir: Path, extra_flags: str = ""):
    """Build a Tage monitor binary."""
    prefix = "mon"
    config_id = cfg.config_id
    binary = build_dir / f"{prefix}_{config_id}"
    if binary.exists():
        return (cfg, binary, 0.0, None)

    build_dir.mkdir(parents=True, exist_ok=True)
    pred = cfg.predictor_string
    monitor_flags = "-DTAGE_MONITOR -DCHEATING_MODE -DREAD_WRITE_RAM"
    cmd = (f'{CXX} {COMMON_FLAGS} {WARN_FLAGS} {extra_flags} {monitor_flags} '
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


def evaluate_config(cfg, build_dir: Path, trace_paths: list[Path],
                    warmup: int, measure: int, jobs: int,
                    extra_flags: str = "",
                    predictor: str = "tagedirect") -> tuple[dict | None, float, float]:
    """Build + run standard binary on all traces. Returns (metrics, build_t, eval_t)."""
    _, binary, build_time, err = _build_config(cfg, build_dir, extra_flags,
                                                predictor=predictor)
    if err:
        print(f"    BUILD FAIL {cfg.config_id}: {err[:80]}", file=sys.stderr)
        return None, build_time, 0.0

    t0 = time.monotonic()
    runs = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(run_trace, binary, t, warmup, measure): t
                   for t in trace_paths}
        for fut in as_completed(futures):
            data, run_err = fut.result()
            if run_err:
                print(f"      TRACE FAIL {futures[fut].name}: {run_err[:60]}",
                      file=sys.stderr)
                continue
            runs.append(data)

    eval_time = time.monotonic() - t0
    if not runs:
        return None, build_time, eval_time
    return compute_metrics(runs), build_time, eval_time


def evaluate_monitor(cfg, build_dir: Path, monitor_dir: Path,
                     monitor_trace: Path, warmup: int, measure: int,
                     extra_flags: str = "",
                     predictor: str = "tagedirect") -> float:
    """Build + run monitor binary on one trace. Returns build_time."""
    _, binary, bt, err = _build_config(cfg, build_dir, extra_flags,
                                        monitor=True, predictor=predictor)
    if err:
        print(f"    MON BUILD FAIL {cfg.config_id}: {err[:80]}", file=sys.stderr)
        return bt
    out = monitor_dir / f"{cfg.config_id}.txt"
    mon_err = run_monitor(binary, monitor_trace, warmup, measure, out)
    if mon_err:
        print(f"    MON RUN FAIL {cfg.config_id}: {mon_err[:60]}", file=sys.stderr)
    return bt


def make_row(iteration, param, cfg, metrics, is_improvement, build_time, eval_time):
    row = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "iteration": iteration,
        "config_id": cfg.config_id,
        "varied_param": param,
        "is_improvement": is_improvement,
        "build_time_s": f"{build_time:.1f}",
        "eval_time_s": f"{eval_time:.1f}",
        "predictor_string": cfg.predictor_string,
    }
    for k in ["mpki", "epi", "ipc", "vfs"]:
        v = metrics.get(k, 0)
        row[k] = f"{float(v):.6f}" if isinstance(v, (int, float)) else v
    for k in ["p1_lat", "p2_lat", "n_traces"]:
        row[k] = metrics.get(k, "")
    row.update(cfg.to_dict())
    return row


# ============================================================================
# Gradient ascent
# ============================================================================

def gradient_ascent(args):
    # Load YAML config
    with open(args.config) as f:
        ycfg = yaml.safe_load(f)

    epi_budget = ycfg["epi_budget"]
    max_iterations = ycfg["max_iterations"]

    # Select config class based on predictor type
    predictor = args.predictor
    if predictor == "tage":
        ConfigClass = SweepConfig
        config_fields = SWEEP_CONFIG_FIELDS
    else:
        ConfigClass = TDConfig
        config_fields = TD_CONFIG_FIELDS

    # Build parameter spaces
    param_spaces = {}
    for name, spec in ycfg["parameters"].items():
        param_spaces[name] = ParamSpace(name, spec)

    # Starting config
    start = ycfg["starting_point"]
    current = ConfigClass(**{k: v for k, v in start.items() if k in config_fields})

    # Resolve traces
    if args.traces:
        trace_paths = [Path(t) for t in args.traces]
    else:
        trace_paths = resolve_traces(args.trace_dir)
    if not trace_paths:
        print("ERROR: no traces found", file=sys.stderr)
        sys.exit(1)
    print(f"Using {len(trace_paths)} traces", file=sys.stderr)

    # Find monitor trace (prefer gcc)
    monitor_trace = trace_paths[0]
    for t in trace_paths:
        if "gcc" in t.name.lower():
            monitor_trace = t
            break
    print(f"Monitor trace: {monitor_trace.name}", file=sys.stderr)

    # Output paths
    csv_path = args.output
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_path = csv_path.parent / "checkpoint.json"
    build_dir = csv_path.parent / "bin"
    monitor_dir = csv_path.parent / "monitors"
    monitor_dir.mkdir(parents=True, exist_ok=True)

    cache = EvalCache(csv_path)
    perturbation = args.perturbation
    total_evals = 0
    start_iteration = 1

    # Resume
    if args.resume:
        ckpt = load_checkpoint(checkpoint_path)
        if ckpt:
            current = ConfigClass(**{k: v for k, v in ckpt["current_config"].items()
                                     if k in config_fields})
            perturbation = ckpt.get("perturbation", args.perturbation)
            total_evals = ckpt.get("total_evals", 0)
            start_iteration = ckpt["iteration"] + 1
            print(f"Resuming from iteration {start_iteration}, "
                  f"config={current.config_id}, perturbation={perturbation}",
                  file=sys.stderr)

    # Evaluate starting config
    cached = cache.get(current.config_id)
    if cached:
        current_mpki = cached["mpki"]
        current_epi = cached["epi"]
        print(f"Start (cached): MPKI={current_mpki:.3f} EPI={current_epi:.0f}",
              file=sys.stderr)
    else:
        print(f"Evaluating start config {current.config_id}...", file=sys.stderr)
        metrics, bt, et = evaluate_config(
            current, build_dir, trace_paths,
            args.warmup, args.measure, args.jobs, args.extra_flags,
            predictor=predictor)
        if not metrics:
            print("ERROR: start config failed", file=sys.stderr)
            sys.exit(1)
        current_mpki = metrics["mpki"]
        current_epi = metrics["epi"]
        cache.put(current.config_id, metrics)
        total_evals += 1
        row = make_row(0, "init", current, metrics, True, bt, et)
        append_csv_row(csv_path, row, _gradient_fields(predictor))
        # Run monitor for starting config
        evaluate_monitor(current, build_dir, monitor_dir, monitor_trace,
                         args.warmup, args.measure, args.extra_flags,
                         predictor=predictor)
        print(f"Start: MPKI={current_mpki:.3f} EPI={current_epi:.0f}",
              file=sys.stderr)

    if current_epi > epi_budget:
        print(f"WARNING: starting config EPI={current_epi:.0f} exceeds budget {epi_budget}",
              file=sys.stderr)

    print(f"\n{'='*70}", file=sys.stderr)
    print(f"Gradient ascent: MPKI={current_mpki:.3f} EPI budget={epi_budget}",
          file=sys.stderr)
    print(f"Perturbation={perturbation}, params={len(param_spaces)}",
          file=sys.stderr)
    print(f"{'='*70}\n", file=sys.stderr)

    # Main loop
    for iteration in range(start_iteration, start_iteration + max_iterations):
        iter_start = time.monotonic()
        print(f"\n--- Iteration {iteration} (perturbation={perturbation}) ---",
              file=sys.stderr)

        # Generate all neighbors
        neighbors: list[tuple[str, Any, TDConfig]] = []  # (param, value, config)
        for pname, pspace in param_spaces.items():
            cur_val = getattr(current, pname)
            for nv in pspace.neighbors(cur_val, perturbation):
                candidate = current.replace(**{pname: nv})
                # Constraint: minh < maxh (field names differ by predictor)
                minh = getattr(candidate, 'minh', None) or getattr(candidate, 'minhist', 0)
                maxh = getattr(candidate, 'maxh', None) or getattr(candidate, 'maxhist', 999)
                if minh >= maxh:
                    continue
                # Constraint: table sizes must be power of 2
                # (GeoSize with size_ratio>1 might produce non-pow2, but
                # the C++ constexpr rounds up, so we're safe)
                neighbors.append((pname, nv, candidate))

        if not neighbors:
            print("  No valid neighbors, stopping.", file=sys.stderr)
            break

        # Deduplicate by config_id
        seen = {current.config_id}
        unique_neighbors = []
        for pname, nv, cfg in neighbors:
            if cfg.config_id not in seen:
                seen.add(cfg.config_id)
                unique_neighbors.append((pname, nv, cfg))
        neighbors = unique_neighbors

        print(f"  {len(neighbors)} unique neighbors to evaluate", file=sys.stderr)

        # Phase 1: Build all standard binaries in parallel
        print(f"  Building standard binaries...", file=sys.stderr)
        build_results = {}
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futs = {pool.submit(_build_config, cfg, build_dir, args.extra_flags,
                                predictor=predictor): (pname, nv, cfg)
                    for pname, nv, cfg in neighbors}
            for fut in as_completed(futs):
                pname, nv, cfg = futs[fut]
                _, binary, bt, err = fut.result()
                if err:
                    print(f"    BUILD FAIL {pname}={nv}: {err[:60]}", file=sys.stderr)
                else:
                    build_results[cfg.config_id] = (pname, nv, cfg, binary, bt)

        # Phase 2: Evaluate (check cache first, run uncached in parallel)
        results: list[tuple[str, Any, TDConfig, dict, float, float]] = []

        to_run = []
        for cid, (pname, nv, cfg, binary, bt) in build_results.items():
            cached = cache.get(cid)
            if cached:
                results.append((pname, nv, cfg, cached, bt, 0.0))
            else:
                to_run.append((pname, nv, cfg, binary, bt))

        if to_run:
            print(f"  Running {len(to_run)} configs × {len(trace_paths)} traces...",
                  file=sys.stderr)

            # Run all (config × trace) pairs in parallel
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                # Submit all trace runs grouped by config
                cfg_futures: dict[str, list] = {}
                meta_by_cid = {}
                for pname, nv, cfg, binary, bt in to_run:
                    cid = cfg.config_id
                    meta_by_cid[cid] = (pname, nv, cfg, bt)
                    cfg_futures[cid] = []
                    for trace in trace_paths:
                        fut = pool.submit(run_trace, binary, trace,
                                          args.warmup, args.measure)
                        cfg_futures[cid].append(fut)

                # Collect results per config
                for cid, futs_list in cfg_futures.items():
                    pname, nv, cfg, bt = meta_by_cid[cid]
                    runs = []
                    t0 = time.monotonic()
                    for fut in futs_list:
                        data, run_err = fut.result()
                        if data:
                            runs.append(data)
                    et = time.monotonic() - t0
                    total_evals += 1
                    if runs:
                        metrics = compute_metrics(runs)
                        cache.put(cid, metrics)
                        results.append((pname, nv, cfg, metrics, bt, et))
                    else:
                        print(f"    EVAL FAIL {pname}={nv}", file=sys.stderr)

        # Phase 3: Build and run monitor for all evaluated configs
        print(f"  Running monitors...", file=sys.stderr)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            mon_futs = []
            for pname, nv, cfg, metrics, bt, et in results:
                mon_out = monitor_dir / f"{cfg.config_id}.txt"
                if not mon_out.exists():
                    mon_futs.append(pool.submit(
                        evaluate_monitor, cfg, build_dir, monitor_dir,
                        monitor_trace, args.warmup, args.measure,
                        args.extra_flags, predictor=predictor))
            for fut in as_completed(mon_futs):
                fut.result()  # collect errors

        # Filter by EPI budget and find best
        best_neighbor = None
        best_mpki = current_mpki
        best_metrics = None
        best_param = None
        best_val = None

        for pname, nv, cfg, metrics, bt, et in results:
            mpki = metrics["mpki"]
            epi = metrics["epi"]
            is_feasible = epi <= epi_budget
            is_better = is_feasible and mpki < best_mpki

            # Log to CSV
            row = make_row(iteration, pname, cfg, metrics, is_better, bt, et)
            append_csv_row(csv_path, row, _gradient_fields(predictor))

            marker = ""
            if not is_feasible:
                marker = " [OVER BUDGET]"
            elif is_better:
                marker = " *** BEST"

            print(f"    {pname}={nv}: MPKI={mpki:.3f} EPI={epi:.0f}{marker}",
                  file=sys.stderr)

            if is_better:
                best_mpki = mpki
                best_metrics = metrics
                best_neighbor = cfg
                best_param = pname
                best_val = nv

        iter_time = time.monotonic() - iter_start

        if best_neighbor:
            current = best_neighbor
            current_mpki = best_mpki
            print(f"\n  >> Improved: {best_param}={best_val} "
                  f"MPKI={best_mpki:.3f} ({iter_time:.0f}s)", file=sys.stderr)
            # Save best iteration to saved/
            save_best_iteration(csv_path.parent / "saved", iteration,
                                current, best_metrics)
        elif perturbation > 1:
            perturbation = max(1, perturbation // 2)
            print(f"\n  >> No improvement, reducing perturbation to {perturbation} "
                  f"({iter_time:.0f}s)", file=sys.stderr)
        else:
            print(f"\n  >> Converged! No improving neighbor. ({iter_time:.0f}s)",
                  file=sys.stderr)
            save_checkpoint(checkpoint_path, iteration, current,
                            current_mpki, perturbation, total_evals)
            break

        save_checkpoint(checkpoint_path, iteration, current,
                        current_mpki, perturbation, total_evals)

    # Final summary
    print(f"\n{'='*70}", file=sys.stderr)
    print(f"FINAL: MPKI={current_mpki:.3f}", file=sys.stderr)
    print(f"Config: {current.config_id}", file=sys.stderr)
    print(f"Params: {current.to_dict()}", file=sys.stderr)
    print(f"Predictor: {current.predictor_string}", file=sys.stderr)
    print(f"Total evaluations: {total_evals}", file=sys.stderr)
    print(f"Results: {csv_path}", file=sys.stderr)
    print(f"{'='*70}", file=sys.stderr)


# ============================================================================
# Report
# ============================================================================

def print_report(csv_path: Path, top_n: int = 20, predictor: str = "tagedirect"):
    rows = load_csv_rows(csv_path)
    if not rows:
        print("No results.", file=sys.stderr)
        return

    # Deduplicate by config_id, keeping lowest MPKI
    by_config: dict[str, dict] = {}
    for r in rows:
        cid = r.get("config_id", "")
        try:
            mpki = float(r["mpki"])
        except (ValueError, KeyError):
            continue
        if cid not in by_config or mpki < float(by_config[cid]["mpki"]):
            by_config[cid] = r

    sorted_configs = sorted(by_config.values(), key=lambda r: float(r["mpki"]))

    # Trajectory
    print(f"\n=== Gradient Ascent Trajectory ({len(rows)} evaluations) ===\n")
    print(f"{'#':>4}  {'Iter':>4}  {'Param':>14}  {'MPKI':>7}  {'EPI':>7}  "
          f"{'VFS':>8}  {'Better':>6}  {'ID':>8}")
    print("-" * 80)

    for i, r in enumerate(rows[:80]):
        marker = " *" if str(r.get("is_improvement", "")).lower() == "true" else ""
        try:
            print(f"{i+1:4d}  {r.get('iteration','?'):>4}  "
                  f"{r.get('varied_param','?'):>14}  "
                  f"{float(r.get('mpki',0)):7.3f}  "
                  f"{float(r.get('epi',0)):7.0f}  "
                  f"{float(r.get('vfs',0)):8.6f}  "
                  f"{marker:>6}  {r.get('config_id','?'):>8}")
        except ValueError:
            pass

    # Leaderboard
    print(f"\n=== Leaderboard (top {min(top_n, len(sorted_configs))}) ===\n")
    print(f"{'Rank':>4}  {'ID':>8}  {'MPKI':>7}  {'EPI':>7}  {'VFS':>8}  "
          f"{'T':>2}  {'Size':>5}  {'Tag':>3}  {'MinH':>4}  {'MaxH':>4}  "
          f"{'P1sz':>5}  {'P1h':>3}  {'SR':>2}")
    print("-" * 90)

    for i, r in enumerate(sorted_configs[:top_n]):
        try:
            print(f"{i+1:4d}  {r['config_id']:>8}  "
                  f"{float(r.get('mpki',0)):7.3f}  "
                  f"{float(r.get('epi',0)):7.0f}  "
                  f"{float(r.get('vfs',0)):8.6f}  "
                  f"{r.get('num_tables','?'):>2}  "
                  f"{r.get('size','?'):>5}  "
                  f"{r.get('tag','?'):>3}  "
                  f"{r.get('minh','?'):>4}  "
                  f"{r.get('maxh','?'):>4}  "
                  f"{r.get('p1_table_size','?'):>5}  "
                  f"{r.get('p1_hist','?'):>3}  "
                  f"{r.get('size_ratio','?'):>2}")
        except ValueError:
            pass

    # Best config
    if sorted_configs:
        best = sorted_configs[0]
        print(f"\n=== Best Config ===")
        print(f"  Config ID: {best['config_id']}")
        print(f"  MPKI: {float(best['mpki']):.3f}  EPI: {float(best['epi']):.0f}  "
              f"VFS: {float(best['vfs']):.6f}")
        cfg_fields = SWEEP_CONFIG_FIELDS if predictor == "tage" else TD_CONFIG_FIELDS
        for f in cfg_fields:
            if f in best:
                print(f"  {f}: {best[f]}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Gradient ascent optimizer for TageDirect/Tage predictor")
    parser.add_argument("--config", type=Path,
                        default=Path("configs/gradient_default.yaml"))
    parser.add_argument("--predictor", choices=["tagedirect", "tage"],
                        default="tagedirect",
                        help="Predictor type: 'tagedirect' (default) or 'tage'")
    parser.add_argument("--trace-dir", type=Path, default=Path("./traces"))
    parser.add_argument("--traces", nargs="*", default=None,
                        help="Override trace list (paths to .gz files)")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 4)
    parser.add_argument("--perturbation", type=int, default=1)
    parser.add_argument("--warmup", type=int, default=1000000)
    parser.add_argument("--measure", type=int, default=40000000)
    parser.add_argument("-o", "--output", type=Path, default=None)
    parser.add_argument("--extra-flags", type=str, default="")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--report", type=Path, default=None,
                        help="Print report from existing CSV (no build/run)")
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()

    if args.output is None:
        args.output = Path(f"out/gradient_{args.predictor}/results.csv")

    if args.report:
        print_report(args.report, args.top, predictor=args.predictor)
        return

    gradient_ascent(args)


if __name__ == "__main__":
    main()
