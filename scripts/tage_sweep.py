#!/usr/bin/env python3
"""Tage predictor parameter sweep: build, run traces, aggregate results."""

import argparse
import csv
import hashlib
import io
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, fields, asdict
from pathlib import Path
from typing import Optional


# ============================================================================
# Configuration dataclass
# ============================================================================

@dataclass(frozen=True)
class TageConfig:
    # Table config (SweepTableConfig params)
    num_tables: int = 8
    table_size: int = 2048
    tag_width: int = 11
    ctr_width: int = 1
    hyst_width: int = 2
    u_width: int = 1
    minhist: int = 2
    maxhist: int = 100
    size_ratio: int = 1
    # Tage-level params
    bimodal_size: int = 4096
    shared_tag: bool = True
    shared_u: bool = True
    shared_hys: bool = True
    u_stor_ff: bool = False
    decay_ctr: int = 1024
    # P1
    p1_use_gshare: bool = True
    p1_table_size: int = 16384
    p1_hist: int = 6
    # Meta
    use_meta: bool = True
    metabits: int = 4
    metapipe: int = 2
    meta_table_size: int = 0
    # Path history
    use_path_hist: bool = False
    path_hist_width: int = 27
    path_bits: int = 6
    # Alloc
    max_alloc: int = 1

    @property
    def config_id(self) -> str:
        return hashlib.md5(self.predictor_string.encode()).hexdigest()[:8]

    @property
    def predictor_string(self) -> str:
        tc = (f"SweepTableConfig<{self.num_tables},{self.table_size},"
              f"{self.tag_width},{self.ctr_width},{self.hyst_width},"
              f"{self.u_width},{self.minhist},{self.maxhist},{self.size_ratio}>")

        # AllocConfig: only override if non-default
        if self.max_alloc != 1:
            alloc = f"DefaultAllocConfig"  # We'll handle via template param
        else:
            alloc = "DefaultAllocConfig"

        parts = [tc, alloc]
        # Positional template params after TableCfg, AllocCfg
        parts.append(str(16))  # FETCH_WIDTH (fixed)
        parts.append(str(self.bimodal_size))
        parts.append(str(1))   # BR_P_ENTRY (fixed)
        parts.append(str(1))   # NUM_BANKS (fixed)
        parts.append(str("false"))  # USE_AHEAD (fixed)
        parts.append("true" if self.shared_tag else "false")
        parts.append("true" if self.shared_u else "false")
        parts.append("true" if self.shared_hys else "false")
        parts.append("true" if self.u_stor_ff else "false")
        parts.append(str(self.decay_ctr))
        parts.append("DefaultResetFn")
        parts.append("false")  # USE_FF_CACHE (fixed)
        # P1
        parts.append("true" if self.p1_use_gshare else "false")
        parts.append(str(self.p1_table_size))
        parts.append(str(self.p1_hist))
        # Meta
        parts.append("true" if self.use_meta else "false")
        parts.append(str(self.metabits))
        parts.append(str(self.metapipe))
        parts.append(str(self.meta_table_size))
        # Path history
        parts.append("true" if self.use_path_hist else "false")
        parts.append(str(self.path_hist_width))
        parts.append(str(self.path_bits))

        return f"Tage<{','.join(parts)}>"


# ============================================================================
# Config generation
# ============================================================================

DEFAULTS = TageConfig()


def make_variant(group: str, **overrides) -> TageConfig:
    """Create a config variant, skipping if identical to defaults."""
    cfg = TageConfig(**overrides)
    if cfg == DEFAULTS:
        return None
    return cfg


def generate_tier0() -> list[TageConfig]:
    return [DEFAULTS]


def generate_tier1() -> list[TageConfig]:
    configs = []

    def add(**kw):
        cfg = TageConfig(**kw)
        if cfg != DEFAULTS and cfg not in configs:
            configs.append(cfg)

    # A: NUM_TABLES
    for n in [4, 6, 10, 12]:
        add(num_tables=n)

    # B: TABLE_SIZE
    for s in [512, 1024, 4096]:
        add(table_size=s)

    # B2: SIZE_RATIO
    for r in [2, 4]:
        add(size_ratio=r)

    # C: TAG_WIDTH
    for t in [7, 9, 13]:
        add(tag_width=t)

    # D: CTR_WIDTH
    for c in [2, 3]:
        add(ctr_width=c)

    # E: HYST_WIDTH
    for h in [1, 3]:
        add(hyst_width=h)

    # F: History ranges
    for minh, maxh in [(2, 50), (2, 200), (4, 100), (4, 200), (2, 300), (8, 200), (2, 400)]:
        add(minhist=minh, maxhist=maxh)

    # G: BIMODAL_SIZE
    for b in [1024, 2048, 8192, 16384]:
        add(bimodal_size=b)

    # H: DECAY_CTR
    for d in [256, 512, 2048]:
        add(decay_ctr=d)

    # I: P1 config
    add(p1_use_gshare=False)
    for sz in [8192, 32768]:
        add(p1_table_size=sz)
    for h in [4, 8]:
        add(p1_hist=h)

    # J: Meta
    add(use_meta=False)
    for mb in [2, 8]:
        add(metabits=mb)
    for mp in [1, 3]:
        add(metapipe=mp)
    add(meta_table_size=4096)

    # K: Path history
    for phw, pb in [(16, 4), (27, 6), (32, 8)]:
        add(use_path_hist=True, path_hist_width=phw, path_bits=pb)

    # L: Sharing combos
    for st, su, sh in [(False, True, True), (True, False, True),
                        (True, True, False), (False, False, False)]:
        add(shared_tag=st, shared_u=su, shared_hys=sh)

    # M: U_STOR_FF
    add(u_stor_ff=True)

    # N: MAX_ALLOC
    add(max_alloc=2)

    return configs


def generate_gaussian(n_samples: int = 100, seed: int = 42) -> list[TageConfig]:
    """Generate configs by perturbing multiple params simultaneously around defaults.

    Each continuous param is sampled from a Gaussian centered on the default,
    with sigma ~20-30% of the reasonable range. Values are snapped to the
    nearest valid discrete value. Duplicates and configs matching defaults
    are filtered out.
    """
    import random
    rng = random.Random(seed)

    # (field_name, default, valid_values_or_range, sigma_in_index_space)
    # For "choices" params: sigma is in index space over the sorted valid set.
    # For power-of-2 params: we perturb the log2 and snap.
    param_specs = [
        # name, default, choices (sorted), sigma (in index units)
        ("num_tables",   8,    [4, 6, 8, 10, 12],           1.0),
        ("table_size",   2048, [256, 512, 1024, 2048, 4096, 8192], 1.2),
        ("size_ratio",   1,    [1, 2, 4],                   0.7),
        ("tag_width",    11,   [7, 8, 9, 10, 11, 12, 13, 15], 1.5),
        ("ctr_width",    1,    [1, 2, 3],                   0.6),
        ("hyst_width",   2,    [1, 2, 3],                   0.6),
        ("u_width",      1,    [1, 2],                      0.4),
        ("minhist",      2,    [2, 3, 4, 6, 8],             0.8),
        ("maxhist",      100,  [50, 75, 100, 130, 160, 200, 300, 400], 1.5),
        ("bimodal_size", 4096, [1024, 2048, 4096, 8192, 16384], 1.0),
        ("decay_ctr",    1024, [256, 512, 1024, 2048],      0.8),
        ("p1_table_size",16384,[4096, 8192, 16384, 32768],  0.8),
        ("p1_hist",      6,    [4, 5, 6, 7, 8, 10],        1.0),
        ("metabits",     4,    [2, 3, 4, 6, 8],             0.8),
        ("metapipe",     2,    [1, 2, 3],                   0.5),
    ]

    # Build index of default position for each param
    default_indices = {}
    for name, default, choices, _ in param_specs:
        default_indices[name] = choices.index(default)

    configs = []
    seen = set()
    attempts = 0
    max_attempts = n_samples * 20

    while len(configs) < n_samples and attempts < max_attempts:
        attempts += 1
        overrides = {}

        for name, default, choices, sigma in param_specs:
            di = default_indices[name]
            # Sample from Gaussian centered on default index
            sampled_idx = rng.gauss(di, sigma)
            # Clamp and round to nearest valid index
            idx = max(0, min(len(choices) - 1, round(sampled_idx)))
            val = choices[idx]
            if val != default:
                overrides[name] = val

        # Ensure minhist < maxhist
        minh = overrides.get("minhist", 2)
        maxh = overrides.get("maxhist", 100)
        if minh >= maxh:
            continue

        # Skip if identical to defaults
        if not overrides:
            continue

        try:
            cfg = TageConfig(**overrides)
        except TypeError:
            continue

        if cfg.config_id in seen:
            continue
        seen.add(cfg.config_id)
        configs.append(cfg)

    return configs


def generate_configs(tier: int, n_gaussian: int = 0, seed: int = 42) -> list[TageConfig]:
    configs = generate_tier0()
    if tier >= 1:
        configs.extend(generate_tier1())
    if n_gaussian > 0:
        configs.extend(generate_gaussian(n_gaussian, seed))
        # Deduplicate
        seen = set()
        unique = []
        for c in configs:
            if c.config_id not in seen:
                seen.add(c.config_id)
                unique.append(c)
        configs = unique
    return configs


# ============================================================================
# Build & Run
# ============================================================================

CXX = os.environ.get("CXX", "g++")
COMMON_FLAGS = "-std=c++20 -O3 -DVERBOSE"
WARN_FLAGS = ("-Wall -Wextra -pedantic -Wold-style-cast -Werror "
              "-Wno-deprecated-declarations -Wno-mismatched-tags")


def build_config(cfg: TageConfig, build_dir: Path, extra_flags: str = "") -> tuple[TageConfig, Path, float, Optional[str]]:
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
        elapsed = time.monotonic() - t0
        return (cfg, binary, elapsed, None)
    except subprocess.CalledProcessError as e:
        elapsed = time.monotonic() - t0
        return (cfg, binary, elapsed, e.stderr[:500])
    except subprocess.TimeoutExpired:
        return (cfg, binary, 300.0, "TIMEOUT")


def run_trace(binary: Path, trace: Path, warmup: int, measure: int) -> tuple[Optional[dict], Optional[str]]:
    """Run a binary on a trace. Returns (parsed_results, error_or_None)."""
    trace_name = trace.stem.replace("_trace", "")
    cmd = [str(binary), str(trace), trace_name, str(warmup), str(measure),
           "--format", "human"]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=600, check=True)
    except subprocess.CalledProcessError as e:
        return (None, f"Run failed: {e.stderr[:300]}")
    except subprocess.TimeoutExpired:
        return (None, "TIMEOUT")

    stdout = result.stdout
    stderr = result.stderr

    def parse_val(pattern: str, text: str, multiline: bool = False) -> str:
        flags = re.MULTILINE if multiline else 0
        m = re.search(pattern + r"\s*:\s*(\S+)", text, flags)
        return m.group(1) if m else ""

    data = {
        "trace": trace_name,
        "instructions": parse_val("instructions", stdout),
        "branches": parse_val(r"^branches", stdout, multiline=True),
        "cond_branches": parse_val("conditional branches", stdout),
        "predictions": parse_val("predictions", stdout),
        "extra_cycles": parse_val("extra_cycles", stdout),
        "disagreements": parse_val("short mispredictions", stdout),
        "disagree_end": parse_val("block-ending short misp", stdout),
        "mispredictions": parse_val(r"^mispredictions", stdout, multiline=True),
        "p1_latency": parse_val("p1 latency", stdout),
        "p2_latency": parse_val("p2 latency", stdout),
        "epi": parse_val("energy per instruction", stdout),
        # HW metrics from stderr
        "storage_bits": parse_val(r"storage \(bits\)", stderr),
        "transistors": parse_val("transistors", stderr),
        "sram_area_mm2": parse_val(r"SRAM area \(mm2\)", stderr),
        "dynamic_power_mW": parse_val(r"dynamic power \(mW\)", stderr),
        "static_power_mW": parse_val(r"static power \(mW\)", stderr),
    }

    return (data, None)


# ============================================================================
# VFS computation (matches compare_predictors.sh)
# ============================================================================

MISP_PENALTY = 8


def compute_vfs(instr: int, preds: int, extra: int,
                diverge: int, diverge_end: int, misp: int,
                p1_lat: float, p2_lat: float, epi: float) -> dict:
    """Compute derived metrics. Returns dict with mpki, ipc, cpi, vfs."""
    p1_latency = math.ceil(p1_lat)
    p2_latency = math.ceil(p2_lat)

    mpki = misp / instr * 1000
    mpi = misp / instr

    if p2_latency <= p1_latency:
        cycles = preds * max(1, p2_latency)
    else:
        cycles = (preds * max(1, p1_latency) +
                  diverge * p2_latency -
                  diverge_end * max(1, p1_latency))
    cycles += extra
    ipc = instr / cycles if cycles > 0 else 0

    cpi = mpi * (MISP_PENALTY + p2_latency)

    IPCcbp0 = 8
    CPIcbp0 = 0.0315
    EPIcbp0 = 1000
    ALPHA = 1.625
    BETA = 4 * ALPHA / (ALPHA - 1) ** 2
    GAMMA = 2 / (ALPHA - 1)
    cbp_energy_ratio = 0.05
    WPI0 = IPCcbp0 * CPIcbp0
    WPI = ipc * cpi
    speedup = (ipc / IPCcbp0) * (1 + WPI0) / (1 + WPI)
    LAMBDA = 1 / (1 + WPI0 / 2) - cbp_energy_ratio
    normalizedEPI = ((epi / EPIcbp0) * cbp_energy_ratio +
                     LAMBDA * speedup ** GAMMA) * (1 + WPI / 2)
    vfs = speedup * ALPHA * (1 - 2 / (1 + math.sqrt(1 + BETA / (speedup * normalizedEPI))))

    return {"mpki": mpki, "ipc": ipc, "cpi": cpi, "vfs": vfs}


# ============================================================================
# CSV I/O
# ============================================================================

CONFIG_FIELDS = [f.name for f in fields(TageConfig)]
RUN_FIELDS = [
    "trace", "instructions", "branches", "cond_branches", "predictions",
    "extra_cycles", "disagreements", "disagree_end", "mispredictions",
    "p1_latency", "p2_latency", "epi",
    "storage_bits", "transistors", "sram_area_mm2",
    "dynamic_power_mW", "static_power_mW",
]
DERIVED_FIELDS = ["mpki", "ipc", "cpi", "vfs"]
ALL_FIELDS = ["config_id"] + CONFIG_FIELDS + RUN_FIELDS + DERIVED_FIELDS


def load_existing_results(path: Path) -> set[tuple[str, str]]:
    """Load (config_id, trace) pairs from existing CSV."""
    done = set()
    if not path.exists():
        return done
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            done.add((row.get("config_id", ""), row.get("trace", "")))
    return done


def write_header(f):
    writer = csv.DictWriter(f, fieldnames=ALL_FIELDS)
    writer.writeheader()
    return writer


def append_row(path: Path, row: dict):
    need_header = not path.exists() or path.stat().st_size == 0
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ALL_FIELDS)
        if need_header:
            writer.writeheader()
        writer.writerow(row)


# ============================================================================
# Report
# ============================================================================

def print_report(csv_path: Path, top_n: int = 20):
    """Aggregate across traces and print leaderboard."""
    if not csv_path.exists():
        print(f"No results file: {csv_path}", file=sys.stderr)
        return

    with open(csv_path) as f:
        rows = list(csv.DictReader(f))

    if not rows:
        print("No results to report.", file=sys.stderr)
        return

    # Group by config_id
    from collections import defaultdict
    by_config: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_config[row["config_id"]].append(row)

    # Aggregate per config
    results = []
    for cid, runs in by_config.items():
        n = len(runs)
        try:
            avg_mpki = sum(float(r["mpki"]) for r in runs) / n
            avg_epi = sum(float(r["epi"]) for r in runs) / n
            # Harmonic mean IPC
            ipc_inv_sum = sum(1.0 / float(r["ipc"]) for r in runs if float(r["ipc"]) > 0)
            hmean_ipc = n / ipc_inv_sum if ipc_inv_sum > 0 else 0
            max_p1 = max(float(r["p1_latency"]) for r in runs)
            max_p2 = max(float(r["p2_latency"]) for r in runs)
            avg_vfs = sum(float(r["vfs"]) for r in runs) / n
        except (ValueError, ZeroDivisionError):
            continue

        # Grab config params from first run
        sample = runs[0]
        results.append({
            "config_id": cid,
            "traces": n,
            "mpki": avg_mpki,
            "ipc": hmean_ipc,
            "epi": avg_epi,
            "p1_lat": max_p1,
            "p2_lat": max_p2,
            "vfs": avg_vfs,
            "predictor": sample.get("num_tables", "?") + "T/" +
                         sample.get("table_size", "?") + "/" +
                         sample.get("tag_width", "?") + "b",
        })

    results.sort(key=lambda r: r["vfs"], reverse=True)

    # Find baseline
    baseline_vfs = None
    for r in results:
        if r["config_id"] == DEFAULTS.config_id:
            baseline_vfs = r["vfs"]
            break

    # Print
    print(f"\n{'Rank':>4}  {'ID':>8}  {'Config':>16}  {'MPKI':>7}  "
          f"{'IPC':>7}  {'EPI':>7}  {'P2lat':>6}  {'VFS':>8}  {'dVFS':>7}  {'Traces':>6}")
    print("-" * 95)
    for i, r in enumerate(results[:top_n]):
        dvfs = f"{r['vfs'] - baseline_vfs:+.4f}" if baseline_vfs else "N/A"
        marker = " *" if r["config_id"] == DEFAULTS.config_id else ""
        print(f"{i+1:4d}  {r['config_id']:>8}  {r['predictor']:>16}  "
              f"{r['mpki']:7.3f}  {r['ipc']:7.4f}  {r['epi']:7.0f}  "
              f"{r['p2_lat']:6.2f}  {r['vfs']:8.4f}  {dvfs:>7}{marker}")

    if baseline_vfs:
        print(f"\n  * = baseline (config_id={DEFAULTS.config_id})")
    print(f"\n  Total configs: {len(results)}, showing top {min(top_n, len(results))}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Tage predictor parameter sweep")
    parser.add_argument("--tier", type=int, default=0, choices=[0, 1],
                        help="Sweep tier (0=baseline only, 1=one-at-a-time)")
    parser.add_argument("--gaussian", type=int, default=0, metavar="N",
                        help="Add N Gaussian-sampled configs around defaults")
    parser.add_argument("--seed", type=int, default=42,
                        help="RNG seed for Gaussian sampling")
    parser.add_argument("--traces", nargs="+", type=Path,
                        help="Trace files to run")
    parser.add_argument("--warmup", type=int, default=1000)
    parser.add_argument("--measure", type=int, default=10000)
    parser.add_argument("--jobs", "-j", type=int, default=os.cpu_count() or 4)
    parser.add_argument("--build-dir", type=Path, default=Path("out/sweep/bin"))
    parser.add_argument("--output", "-o", type=Path, default=Path("out/sweep/results.csv"))
    parser.add_argument("--extra-flags", type=str, default="")
    parser.add_argument("--resume", action="store_true",
                        help="Skip already-completed (config_id, trace) pairs")
    parser.add_argument("--report", type=Path, default=None,
                        help="Print report from existing CSV (no build/run)")
    parser.add_argument("--top", type=int, default=20)

    args = parser.parse_args()

    if args.report:
        print_report(args.report, args.top)
        return

    if not args.traces:
        print("Error: --traces required (unless using --report)", file=sys.stderr)
        sys.exit(1)

    # Validate traces exist
    for t in args.traces:
        if not t.exists():
            print(f"Error: trace not found: {t}", file=sys.stderr)
            sys.exit(1)

    configs = generate_configs(args.tier, args.gaussian, args.seed)
    print(f"Generated {len(configs)} configs (tier {args.tier}, gaussian {args.gaussian})",
          file=sys.stderr)

    # Load existing results for resume
    done = load_existing_results(args.output) if args.resume else set()
    if done:
        print(f"Resuming: {len(done)} existing (config, trace) pairs", file=sys.stderr)

    # Ensure output dir exists
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.build_dir.mkdir(parents=True, exist_ok=True)

    # Phase 1: Build all configs in parallel
    print(f"\n=== Building {len(configs)} configs (jobs={args.jobs}) ===", file=sys.stderr)
    binaries: dict[str, Path] = {}  # config_id -> binary path
    build_errors = 0

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(build_config, cfg, args.build_dir, args.extra_flags): cfg
            for cfg in configs
        }
        for i, future in enumerate(as_completed(futures), 1):
            cfg, binary, elapsed, err = future.result()
            if err:
                print(f"  [{i}/{len(configs)}] FAIL {cfg.config_id}: {err[:80]}",
                      file=sys.stderr)
                build_errors += 1
            else:
                binaries[cfg.config_id] = binary
                if elapsed > 0:
                    print(f"  [{i}/{len(configs)}] OK   {cfg.config_id} ({elapsed:.1f}s)",
                          file=sys.stderr)
                else:
                    print(f"  [{i}/{len(configs)}] SKIP {cfg.config_id} (cached)",
                          file=sys.stderr)

    print(f"\nBuilt {len(binaries)}/{len(configs)} configs "
          f"({build_errors} errors)", file=sys.stderr)

    # Phase 2: Run traces in parallel
    import threading
    csv_lock = threading.Lock()

    # Build work items
    work = []
    skipped = 0
    for cfg in configs:
        if cfg.config_id not in binaries:
            continue
        binary = binaries[cfg.config_id]
        for trace in args.traces:
            trace_name = trace.stem.replace("_trace", "")
            if (cfg.config_id, trace_name) in done:
                skipped += 1
                continue
            work.append((cfg, binary, trace))

    total_runs = len(work)
    print(f"\n=== Running {total_runs} jobs ({skipped} skipped, jobs={args.jobs}) ===",
          file=sys.stderr)

    completed = 0
    errors = 0

    def run_one(item):
        cfg, binary, trace = item
        run_data, run_err = run_trace(binary, trace, args.warmup, args.measure)
        if run_err:
            return (cfg, trace, None, run_err)

        try:
            derived = compute_vfs(
                int(run_data["instructions"]),
                int(run_data["predictions"]),
                int(run_data["extra_cycles"]),
                int(run_data["disagreements"]),
                int(run_data["disagree_end"]),
                int(run_data["mispredictions"]),
                float(run_data["p1_latency"]),
                float(run_data["p2_latency"]),
                float(run_data["epi"]),
            )
        except (ValueError, ZeroDivisionError) as e:
            return (cfg, trace, None, f"VFS: {e}")

        row = {"config_id": cfg.config_id}
        row.update(asdict(cfg))
        row.update(run_data)
        row.update({k: f"{v:.6f}" for k, v in derived.items()})
        return (cfg, trace, row, None)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_one, item): item for item in work}
        for i, future in enumerate(as_completed(futures), 1):
            cfg, trace, row, err = future.result()
            if err:
                print(f"  [{i}/{total_runs}] FAIL {cfg.config_id} {trace.name}: {err[:60]}",
                      file=sys.stderr)
                errors += 1
            else:
                with csv_lock:
                    append_row(args.output, row)
                completed += 1
                misp = row["mispredictions"]
                vfs_str = row["vfs"]
                print(f"  [{i}/{total_runs}] OK   {cfg.config_id} {trace.name} "
                      f"(misp={misp}, VFS={vfs_str})", file=sys.stderr)

    print(f"\n=== Done: {completed} runs, {skipped} skipped, {errors} errors ===",
          file=sys.stderr)
    print(f"Results: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
