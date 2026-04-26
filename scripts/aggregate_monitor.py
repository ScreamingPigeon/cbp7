#!/usr/bin/env python3
"""Aggregate per-trace TAMonitor CSV files into a mean-per-window CSV and summary.

Usage:
    python3 scripts/aggregate_monitor.py <monitor_dir> <out_csv> <out_summary>

Reads all *_csv.txt files in <monitor_dir>, aligns by window number,
and produces:
  - <out_csv>: per-window mean CSV (same format as single-trace CSV, with
    added trace_count column showing how many traces contributed to each window)
  - <out_summary>: per-trace and aggregate statistics
"""

import sys
import os
import re
from collections import defaultdict
from pathlib import Path


def parse_monitor_csv(path):
    """Parse a single monitor CSV file. Returns (headers, rows) where rows is list of dicts."""
    headers = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# "):
                headers = line[2:].split(",")
                continue
            if line.startswith("#"):
                continue
            if line.startswith("=") or line.startswith("Instructions"):
                break
            if headers:
                vals = line.split(",")
                row = {}
                for i, h in enumerate(headers):
                    try:
                        row[h] = float(vals[i]) if i < len(vals) and vals[i] else 0.0
                    except ValueError:
                        row[h] = 0.0
                rows.append(row)
    return headers, rows


def parse_summary_mpki(path):
    """Extract cumulative MPKI from a summary file."""
    try:
        text = open(path).read()
        m = re.search(r'MPKI:\s+([0-9.]+)', text)
        return float(m.group(1)) if m else None
    except Exception:
        return None


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <monitor_dir> <out_csv> <out_summary>")
        sys.exit(1)

    monitor_dir = Path(sys.argv[1])
    out_csv = sys.argv[2]
    out_summary = sys.argv[3]

    # Find all per-trace CSV files (exclude aggregate)
    csv_files = sorted(monitor_dir.glob("*_csv.txt"))
    csv_files = [f for f in csv_files if "aggregate" not in f.name]

    if not csv_files:
        print(f"No *_csv.txt files found in {monitor_dir}")
        sys.exit(1)

    # Parse all traces
    all_headers = None
    traces = {}  # name -> rows
    for f in csv_files:
        name = f.stem.replace("_csv", "")
        headers, rows = parse_monitor_csv(f)
        if not headers or not rows:
            continue
        if all_headers is None:
            all_headers = headers
        traces[name] = rows

    if not traces:
        print("No valid data found")
        sys.exit(1)

    print(f"Loaded {len(traces)} traces")

    # Determine max windows across all traces
    max_windows = max(len(rows) for rows in traces.values())
    numeric_cols = [h for h in all_headers if h != "win"]

    # Compute per-window mean across traces
    agg_rows = []
    for w in range(max_windows):
        sums = defaultdict(float)
        count = 0
        for name, rows in traces.items():
            if w < len(rows):
                count += 1
                for h in numeric_cols:
                    sums[h] += rows[w].get(h, 0.0)
        if count == 0:
            continue
        row = {"win": w}
        for h in numeric_cols:
            row[h] = sums[h] / count
        row["trace_count"] = count
        agg_rows.append(row)

    # Write aggregate CSV
    out_headers = all_headers + ["trace_count"]
    with open(out_csv, "w") as f:
        f.write("# " + ",".join(out_headers) + "\n")
        for row in agg_rows:
            vals = []
            for h in out_headers:
                v = row.get(h, 0.0)
                if h == "win" or h == "trace_count":
                    vals.append(str(int(v)))
                elif h in ("br",):
                    vals.append(f"{v:.0f}")
                else:
                    vals.append(f"{v:.4f}")
            f.write(",".join(vals) + "\n")

    print(f"Wrote aggregate CSV: {out_csv} ({len(agg_rows)} windows)")

    # Write aggregate summary
    with open(out_summary, "w") as f:
        f.write("=" * 70 + "\n")
        f.write(f"TAMonitor Aggregate Summary — {len(traces)} traces\n")
        f.write("=" * 70 + "\n\n")

        # Per-trace table
        f.write("--- Per-Trace Summary ---\n")
        f.write(f"{'Trace':<45s} {'Windows':>8s} {'MPKI':>10s} {'misp%':>8s} {'alloc%':>8s}\n")
        f.write(f"{'-----':<45s} {'-------':>8s} {'----':>10s} {'-----':>8s} {'------':>8s}\n")

        trace_mpkis = []
        for name in sorted(traces.keys()):
            rows = traces[name]
            n_win = len(rows)

            # Compute trace-level averages from window data
            total_br = sum(r.get("br", 0) for r in rows)
            if total_br > 0:
                avg_misp_pct = sum(r.get("misp%", 0) * r.get("br", 1) for r in rows) / total_br
            else:
                avg_misp_pct = sum(r.get("misp%", 0) for r in rows) / max(1, n_win)

            avg_mpki = sum(r.get("MPKI", 0) for r in rows) / max(1, n_win)
            avg_alloc = sum(r.get("alloc_ok%", 0) for r in rows) / max(1, n_win)

            # Also try to get cumulative MPKI from summary file
            summary_path = monitor_dir / f"{name}_summary.txt"
            cum_mpki = parse_summary_mpki(summary_path)
            mpki_str = f"{cum_mpki:.3f}" if cum_mpki is not None else f"{avg_mpki:.3f}*"
            if cum_mpki is not None:
                trace_mpkis.append(cum_mpki)

            f.write(f"{name:<45s} {n_win:>8d} {mpki_str:>10s} {avg_misp_pct:>7.2f}% {avg_alloc:>7.1f}%\n")

        f.write("\n")

        # Aggregate stats
        if trace_mpkis:
            import statistics
            f.write("--- Aggregate Stats ---\n")
            f.write(f"  Traces:      {len(trace_mpkis)}\n")
            f.write(f"  Mean MPKI:   {statistics.mean(trace_mpkis):.3f}\n")
            f.write(f"  Median MPKI: {statistics.median(trace_mpkis):.3f}\n")
            f.write(f"  Stdev MPKI:  {statistics.stdev(trace_mpkis):.3f}\n" if len(trace_mpkis) > 1 else "")
            f.write(f"  Min MPKI:    {min(trace_mpkis):.3f}\n")
            f.write(f"  Max MPKI:    {max(trace_mpkis):.3f}\n")
            f.write("\n")

        # Per-column aggregate stats across all windows of all traces
        f.write("--- Per-Column Mean (across all windows, all traces) ---\n")
        interesting = ["MPKI", "misp%", "alloc_ok%", "acc_avg", "extra%",
                       "true_blk%", "coll%", "Tlast_zu%", "dir_misp%", "bnd_misp%",
                       "micro_p50_mpki", "micro_p95_mpki", "phase_delta_avg",
                       "jaccard", "pingpong", "cf_fb_only%", "cf_tage_only%",
                       "prov_no_tag%", "prov_sec_rej%", "prov_meta_alt%",
                       "prov_meta_pri%", "prov_no_meta%",
                       "cf_sec_fb_acc%", "cf_sec_tage_acc%"]
        for col in interesting:
            if col not in numeric_cols:
                continue
            all_vals = []
            for rows in traces.values():
                all_vals.extend(r.get(col, 0) for r in rows)
            if not all_vals:
                continue
            import numpy as np
            arr = np.array(all_vals)
            f.write(f"  {col:<20s}  mean={np.mean(arr):.3f}  std={np.std(arr):.3f}"
                    f"  p50={np.median(arr):.3f}  p95={np.percentile(arr, 95):.3f}\n")

    print(f"Wrote aggregate summary: {out_summary}")


if __name__ == "__main__":
    main()
