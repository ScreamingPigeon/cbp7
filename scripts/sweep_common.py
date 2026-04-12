#!/usr/bin/env python3
"""Shared infrastructure for Tage parameter sweep scripts."""

import csv
import hashlib
import math
import os
import subprocess
import time
from dataclasses import dataclass, fields
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
        t = i / max(1, num_tables - 1)
        scale = size_ratio ** (t - 0.5)
        sz = int(base_size * scale)
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
    return list(reversed(lengths))


# ============================================================================
# Config dataclass
# ============================================================================

@dataclass(frozen=True)
class SweepConfig:
    num_tables: int = 8
    base_size: int = 2048
    size_ratio: float = 1.0
    tag_width: int = 11
    ctr_width: int = 1
    hyst_width: int = 2
    u_width: int = 1
    minhist: int = 2
    maxhist: int = 100
    bimodal_size: int = 4096
    decay_ctr: int = 0
    decay_gran: int = 0
    p1_table_size: int = 16384
    p1_hist: int = 6
    metabits: int = 4
    metapipe: int = 2
    # structural parameters
    shared_tag: bool = True
    shared_u: bool = True
    shared_hys: bool = True
    u_stor_ff: bool = False
    use_path_hist: bool = False
    path_hist_width: int = 27
    path_bits: int = 6

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
        sr_int = max(1, round(self.size_ratio))
        tc = (f"SweepTableConfig<{self.num_tables},{self.base_size},"
              f"{self.tag_width},{self.ctr_width},{self.hyst_width},"
              f"{self.u_width},{self.minhist},{self.maxhist},{sr_int}>")
        b = lambda v: "true" if v else "false"
        parts = [tc, "DefaultAllocConfig",
                 "16", str(self.bimodal_size), "1", "1",
                 "false",  # USE_AHEAD
                 b(self.shared_tag), b(self.shared_u), b(self.shared_hys),
                 b(self.u_stor_ff),
                 str(self.decay_ctr), str(self.decay_gran),
                 "DefaultResetFn", "false",  # USE_FF_CACHE
                 "true", str(self.p1_table_size), str(self.p1_hist),
                 "true", str(self.metabits), str(self.metapipe),
                 "0",  # META_TABLE_SIZE (not implemented)
                 b(self.use_path_hist),
                 str(self.path_hist_width), str(self.path_bits)]
        return f"Tage<{','.join(parts)}>"

    def replace(self, **kwargs) -> 'SweepConfig':
        """Return a new config with specified fields changed."""
        d = {f.name: getattr(self, f.name) for f in fields(self)}
        d.update(kwargs)
        return SweepConfig(**d)

    def to_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in fields(self)}


DEFAULTS = SweepConfig()
CONFIG_FIELDS = [f.name for f in fields(SweepConfig)]


# ============================================================================
# Build & Run
# ============================================================================

CXX = os.environ.get("CXX", "g++")
COMMON_FLAGS = "-std=c++20 -O3"
WARN_FLAGS = ("-Wall -Wextra -pedantic -Wold-style-cast "
              "-Wno-deprecated-declarations -Wno-mismatched-tags")


def build_config(cfg: SweepConfig, build_dir: Path, extra_flags: str = "") -> tuple:
    """Build a single config. Returns (cfg, binary_path, build_time, error_or_None)."""
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
    """Run a binary on a trace. Returns (parsed_data, error_or_None)."""
    trace_name = trace.stem.replace("_trace", "")
    cmd = [str(binary), str(trace), trace_name, str(warmup), str(measure)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=1800, check=True)
    except subprocess.CalledProcessError as e:
        return (None, f"Run failed: {e.stderr[:200]}")
    except subprocess.TimeoutExpired:
        return (None, "TIMEOUT")

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
    p1_latency = 0
    p2_latency = 0
    for r in runs:
        p1_latency = max(p1_latency, math.ceil(r["p1_latency"]))
        p2_latency = max(p2_latency, math.ceil(r["p2_latency"]))

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
# CSV helpers
# ============================================================================

def append_csv_row(path: Path, row: dict, fieldnames: list[str]):
    need_header = not path.exists() or path.stat().st_size == 0
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        if need_header:
            writer.writeheader()
        writer.writerow(row)


def load_csv_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def resolve_traces(trace_dir: Path) -> list[Path]:
    """Find available representative traces in trace_dir."""
    trace_paths = []
    missing = 0
    for t in REPR_TRACES:
        p = trace_dir / t
        if p.exists():
            trace_paths.append(p)
        else:
            missing += 1
    return trace_paths
