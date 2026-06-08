#!/usr/bin/env python3
"""Plot per-trace power breakdown comparing two predictors.

Reads two CSVs produced by per_trace_power.py and renders either:
  - stacked bar (default): two stacked bars per trace, one per predictor
  - line:                  one line per component per predictor

Usage:
  python3 scripts/plot_per_trace_power.py A.csv B.csv [options]

Options:
  --mode {bar,line}    Plot style (default: bar)
  --metric {epi,epc}  Use EPI (fJ/inst) or EPC (fJ/cyc) (default: epi)
  --out PATH           Output image (default: out/per_trace_power.png)
  --label-a NAME       Legend label for A (default: derived from filename)
  --label-b NAME       Legend label for B (default: derived from filename)
  --sort {trace,delta} Sort traces by name or by B-A baseline delta (default: trace)
  --width INCHES       Figure width (default: 18)
  --height INCHES      Figure height (default: 6)
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COMPONENTS = ["ram_logic", "fanout", "wiring"]
COMPONENT_LABELS = {
    "ram_logic": "RAM+Logic",
    "fanout":    "Fanout",
    "wiring":    "Wiring",
}
COLORS = {
    "ram_logic": "#4C72B0",   # blue
    "fanout":    "#DD8452",   # orange
    "wiring":    "#55A467",   # green
}


def load_csv(path, metric):
    """Returns dict trace_name -> {component: value, baseline: value}."""
    pfx = "epi" if metric == "epi" else "epc"
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            t = row["trace"]
            out[t] = {
                "baseline": float(row[f"{pfx}_baseline"]),
                **{c: float(row[f"{pfx}_{c}"]) for c in COMPONENTS},
            }
    return out


def default_label(path):
    name = Path(path).stem
    # strip common prefixes/suffixes
    for s in ("ptp_", "_quick", "_full"):
        name = name.replace(s, "")
    return name


def plot_bar(traces, a, b, label_a, label_b, ax, ylabel):
    n = len(traces)
    x = np.arange(n)
    w = 0.4
    for offset, data, hatch in [(-w/2, a, ""), (w/2, b, "//")]:
        bottoms = np.zeros(n)
        for comp in COMPONENTS:
            vals = np.array([data[t][comp] for t in traces])
            ax.bar(x + offset, vals, w, bottom=bottoms,
                   color=COLORS[comp], hatch=hatch, edgecolor="black",
                   linewidth=0.3)
            bottoms += vals
    # Legend: components + predictor hatches
    comp_handles = [plt.Rectangle((0, 0), 1, 1, color=COLORS[c],
                                  label=COMPONENT_LABELS[c]) for c in COMPONENTS]
    pred_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black",
                      hatch="",   label=label_a),
        plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black",
                      hatch="//", label=label_b),
    ]
    ax.legend(handles=comp_handles + pred_handles, ncol=3,
              loc="upper right", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([t.split("_")[0].split(".")[0][:15] for t in traces],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.3)


def plot_line(traces, a, b, label_a, label_b, ax, ylabel):
    n = len(traces)
    x = np.arange(n)
    for comp in COMPONENTS:
        va = np.array([a[t][comp] for t in traces])
        vb = np.array([b[t][comp] for t in traces])
        ax.plot(x, va, "-o",  color=COLORS[comp], linewidth=1.5,
                markersize=4, label=f"{COMPONENT_LABELS[comp]} ({label_a})")
        ax.plot(x, vb, "--s", color=COLORS[comp], linewidth=1.5,
                markersize=4, label=f"{COMPONENT_LABELS[comp]} ({label_b})")
    # Baseline as faint reference
    ba = np.array([a[t]["baseline"] for t in traces])
    bb = np.array([b[t]["baseline"] for t in traces])
    ax.plot(x, ba, ":", color="black", alpha=0.5, label=f"Total {label_a}")
    ax.plot(x, bb, ":", color="gray",  alpha=0.7, label=f"Total {label_b}")
    ax.set_xticks(x)
    ax.set_xticklabels([t.split("_")[0].split(".")[0][:15] for t in traces],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.legend(ncol=2, fontsize=7, loc="upper right")
    ax.grid(axis="y", alpha=0.3)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv_a")
    ap.add_argument("csv_b")
    ap.add_argument("--mode",   choices=["bar", "line"], default="bar")
    ap.add_argument("--metric", choices=["epi", "epc"], default="epi")
    ap.add_argument("--out",    default="out/per_trace_power.png")
    ap.add_argument("--label-a")
    ap.add_argument("--label-b")
    ap.add_argument("--sort",   choices=["trace", "delta"], default="trace")
    ap.add_argument("--width",  type=float, default=18)
    ap.add_argument("--height", type=float, default=6)
    args = ap.parse_args()

    a = load_csv(args.csv_a, args.metric)
    b = load_csv(args.csv_b, args.metric)

    common = sorted(set(a.keys()) & set(b.keys()))
    if not common:
        print("No common traces between the two CSVs", file=sys.stderr)
        sys.exit(1)
    only_a = set(a.keys()) - set(b.keys())
    only_b = set(b.keys()) - set(a.keys())
    if only_a or only_b:
        print(f"Skipping non-overlapping traces: A-only={len(only_a)} B-only={len(only_b)}",
              file=sys.stderr)

    if args.sort == "delta":
        common.sort(key=lambda t: b[t]["baseline"] - a[t]["baseline"], reverse=True)
    traces = common

    label_a = args.label_a or default_label(args.csv_a)
    label_b = args.label_b or default_label(args.csv_b)

    units = "fJ/instruction" if args.metric == "epi" else "fJ/prediction-cycle"
    ylabel = f"Energy ({units})"

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    if args.mode == "bar":
        plot_bar(traces, a, b, label_a, label_b, ax, ylabel)
    else:
        plot_line(traces, a, b, label_a, label_b, ax, ylabel)
    ax.set_title("Power Breakdown")

    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=130)
    print(f"→ {args.out}  ({len(traces)} traces, mode={args.mode}, metric={args.metric})")

    # Terminal diff: per-category average across traces
    units = "fJ/inst" if args.metric == "epi" else "fJ/cyc"
    n = len(traces)
    print(f"\n  Mean over {n} traces ({units})")
    w_label = max(len(label_a), len(label_b), 12)
    print(f"  {'Component':<12} {label_a:>{w_label}} {label_b:>{w_label}} "
          f"{'Δ':>10} {'Δ%':>8}")
    print(f"  {'-'*12} {'-'*w_label} {'-'*w_label} {'-'*10} {'-'*8}")
    for key, label in [("baseline", "Baseline")] + \
                      [(c, COMPONENT_LABELS[c]) for c in COMPONENTS]:
        va = sum(a[t][key] for t in traces) / n
        vb = sum(b[t][key] for t in traces) / n
        d  = vb - va
        dp = (d / va * 100.0) if va else 0.0
        sign = "+" if d >= 0 else ""
        print(f"  {label:<12} {va:>{w_label}.1f} {vb:>{w_label}.1f} "
              f"{sign}{d:>9.1f} {sign}{dp:>7.1f}%")


if __name__ == "__main__":
    main()
