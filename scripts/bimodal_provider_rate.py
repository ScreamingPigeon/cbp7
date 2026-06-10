#!/usr/bin/env python3
"""Compute per-trace bimodal provider rate from TAMonitor summary files.

Usage:
  python3 scripts/bimodal_provider_rate.py <monitor_dir>
  python3 scripts/bimodal_provider_rate.py <monitor_dir> --out out.csv
  python3 scripts/bimodal_provider_rate.py <monitor_dir> --histogram

The metric is:
  Provider Distribution -> Bimodal Provided / total Branches

The script also validates that the Provider Distribution rows sum to the
summary's total Branches count for each trace.
"""

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path


INSTRUCTIONS_RE = re.compile(r"\bBranches:\s+(\d+)")
PROVIDER_ROW_RE = re.compile(r"^\s*(Bimodal|T\d+)\s*\|\s*(\d+)\b")


@dataclass
class TraceRate:
    trace: str
    bimodal_provided: int
    branches: int
    provider_sum: int

    @property
    def rate_pct(self):
        return 100.0 * self.bimodal_provided / self.branches


def parse_summary(path):
    trace = path.name.removesuffix("_summary.txt")
    branches = None
    provider_sum = 0
    bimodal_provided = None
    in_provider_distribution = False

    with path.open() as f:
        for line in f:
            if branches is None:
                match = INSTRUCTIONS_RE.search(line)
                if match:
                    branches = int(match.group(1))

            if line.startswith("Provider Distribution:"):
                in_provider_distribution = True
                continue

            if in_provider_distribution and not line.strip():
                in_provider_distribution = False
                continue

            if not in_provider_distribution:
                continue

            match = PROVIDER_ROW_RE.match(line)
            if not match:
                continue

            provider, provided = match.group(1), int(match.group(2))
            provider_sum += provided
            if provider == "Bimodal":
                bimodal_provided = provided

    if branches is None:
        raise ValueError(f"{path}: missing Branches field")
    if bimodal_provided is None:
        raise ValueError(f"{path}: missing Bimodal row in Provider Distribution")

    return TraceRate(trace, bimodal_provided, branches, provider_sum)


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "trace",
            "bimodal_provided",
            "branches",
            "bimodal_provider_rate_pct",
        ])
        for row in rows:
            writer.writerow([
                row.trace,
                row.bimodal_provided,
                row.branches,
                f"{row.rate_pct:.6f}",
            ])


def print_stats(rows):
    total_bimodal = sum(row.bimodal_provided for row in rows)
    total_branches = sum(row.branches for row in rows)
    mean_pct = sum(row.rate_pct for row in rows) / len(rows)
    weighted_pct = 100.0 * total_bimodal / total_branches
    min_row = min(rows, key=lambda row: row.rate_pct)
    max_row = max(rows, key=lambda row: row.rate_pct)

    print(f"traces={len(rows)}")
    print(f"mean_unweighted_pct={mean_pct:.6f}")
    print(f"weighted_pct={weighted_pct:.6f}")
    print(f"min={min_row.trace} {min_row.rate_pct:.6f}")
    print(f"max={max_row.trace} {max_row.rate_pct:.6f}")


def print_histogram(rows):
    bin_width = 10
    bins = [0] * (100 // bin_width)
    for row in rows:
        bin_index = min(int(row.rate_pct // bin_width), len(bins) - 1)
        bins[bin_index] += 1

    max_count = max(bins)
    bar_width = 40

    print()
    print("histogram_bimodal_provider_rate_pct:")
    for index, count in enumerate(bins):
        start = index * bin_width
        end = start + bin_width
        label = f"{start:3d}-{end:3d}%"
        if index == len(bins) - 1:
            label = f"{start:3d}-100%"
        bar_len = round(bar_width * count / max_count) if max_count else 0
        print(f"{label} {count:4d} {'#' * bar_len}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute bimodal provider rate from *_summary.txt files."
    )
    parser.add_argument("monitor_dir", type=Path)
    parser.add_argument(
        "--out",
        type=Path,
        help="Output CSV path. Default: <monitor_dir>/bimodal_provider_rate.csv",
    )
    parser.add_argument(
        "--allow-provider-sum-mismatch",
        action="store_true",
        help="Warn instead of failing if provider rows do not sum to Branches.",
    )
    parser.add_argument(
        "--histogram",
        action="store_true",
        help="Print a 10-point text histogram of per-trace bimodal provider rates.",
    )
    args = parser.parse_args()

    monitor_dir = args.monitor_dir
    if not monitor_dir.is_dir():
        print(f"error: not a directory: {monitor_dir}", file=sys.stderr)
        return 1

    summary_files = sorted(
        path for path in monitor_dir.glob("*_summary.txt")
        if path.name != "aggregate_summary.txt"
    )
    if not summary_files:
        print(f"error: no *_summary.txt files found in {monitor_dir}", file=sys.stderr)
        return 1

    rows = [parse_summary(path) for path in summary_files]
    mismatches = [
        row for row in rows
        if row.provider_sum != row.branches
    ]
    if mismatches:
        for row in mismatches:
            diff = row.branches - row.provider_sum
            print(
                f"provider sum mismatch: {row.trace}: "
                f"branches={row.branches} provider_sum={row.provider_sum} diff={diff}",
                file=sys.stderr,
            )
        if not args.allow_provider_sum_mismatch:
            return 1

    out = args.out or monitor_dir / "bimodal_provider_rate.csv"
    write_csv(out, rows)
    print(f"wrote {out}")
    print_stats(rows)
    if args.histogram:
        print_histogram(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
