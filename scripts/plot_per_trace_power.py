#!/usr/bin/env python3
"""Plot per-trace power breakdown comparing N predictors (N >= 2).

Reads N CSVs produced by per_trace_power.py and renders either:
  - stacked bar (default): one stacked bar per (trace, predictor)
  - line:                  one line per (component, predictor)

Usage:
  python3 scripts/plot_per_trace_power.py A.csv B.csv [C.csv ...] [options]

Options:
  --mode {bar,line}     Plot style (default: bar)
  --metric {epi,epc}    Use EPI (fJ/inst) or EPC (fJ/cyc) (default: epi)
  --out PATH            Output image (default: out/per_trace_power.png)
  --labels NAME ...     Legend labels, one per CSV (default: derived from filename)
  --sort {trace,delta}  Sort by name, or by (last - first) baseline delta
  --width INCHES        Figure width (default: 18)
  --height INCHES       Figure height (default: 6)
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), (0, (5, 1)), (0, (1, 1))]
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]


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
    for s in ("ptp_", "_quick", "_full"):
        name = name.replace(s, "")
    return name


def predictor_tone(idx, count):
    """Positive values lighten toward white; negative values darken toward black."""
    if count <= 1:
        return 0.0
    return np.linspace(0.45, -0.25, count)[idx]


def apply_tone(color, tone):
    rgb = np.array(mcolors.to_rgb(color))
    if tone >= 0:
        rgb = rgb + (1.0 - rgb) * tone
    else:
        rgb = rgb * (1.0 + tone)
    return mcolors.to_hex(np.clip(rgb, 0.0, 1.0))


def component_color(comp, pred_idx, pred_count):
    return apply_tone(COLORS[comp], predictor_tone(pred_idx, pred_count))


def predictor_color(pred_idx, pred_count):
    return apply_tone("#666666", predictor_tone(pred_idx, pred_count))


def add_split_legends(ax, predictors, include_baseline=False):
    n = len(predictors)
    comp_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=COLORS[c], edgecolor="black",
                      linewidth=0.3, label=COMPONENT_LABELS[c])
        for c in COMPONENTS
    ]
    if include_baseline:
        comp_handles.append(Line2D([0], [0], color="black", linestyle=":",
                                   linewidth=1.5, label="Baseline total"))
    comp_legend = ax.legend(handles=comp_handles, title="Power category",
                            loc="upper left", fontsize=8, title_fontsize=8)
    ax.add_artist(comp_legend)

    pred_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=predictor_color(i, n),
                      edgecolor="black", linewidth=0.3, label=label)
        for i, (label, _) in enumerate(predictors)
    ]
    ax.legend(handles=pred_handles, title="Predictor",
              loc="upper right", fontsize=8, title_fontsize=8)


def plot_bar(traces, predictors, ax, ylabel):
    """predictors: list of (label, data_dict)."""
    n = len(traces)
    N = len(predictors)
    x = np.arange(n)
    # Total group width = 0.85 of slot; each bar = 0.85 / N
    group_w = 0.85
    w = group_w / N
    offsets = [(-group_w/2) + (i + 0.5) * w for i in range(N)]

    for i, ((label, data), offset) in enumerate(zip(predictors, offsets)):
        bottoms = np.zeros(n)
        for comp in COMPONENTS:
            vals = np.array([data[t][comp] for t in traces])
            ax.bar(x + offset, vals, w, bottom=bottoms,
                   color=component_color(comp, i, N), edgecolor="black",
                   linewidth=0.3)
            bottoms += vals
    add_split_legends(ax, predictors)
    ax.set_xticks(x)
    ax.set_xticklabels([t.split("_")[0].split(".")[0][:15] for t in traces],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.3)


def plot_line(traces, predictors, ax, ylabel):
    n = len(traces)
    N = len(predictors)
    x = np.arange(n)
    for i, (label, data) in enumerate(predictors):
        ls = LINESTYLES[i % len(LINESTYLES)]
        mk = MARKERS[i % len(MARKERS)]
        for comp in COMPONENTS:
            vals = np.array([data[t][comp] for t in traces])
            ax.plot(x, vals, ls, marker=mk, color=component_color(comp, i, N),
                    linewidth=1.5, markersize=4)
        base = np.array([data[t]["baseline"] for t in traces])
        ax.plot(x, base, ":", color="black",
                alpha=0.3 + 0.4 * (i / max(N - 1, 1)))

    ax.set_xticks(x)
    ax.set_xticklabels([t.split("_")[0].split(".")[0][:15] for t in traces],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    add_split_legends(ax, predictors, include_baseline=True)
    ax.grid(axis="y", alpha=0.3)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csvs", nargs="+",
                    help="2 or more CSVs produced by per_trace_power.py")
    ap.add_argument("--mode",   choices=["bar", "line"], default="bar")
    ap.add_argument("--metric", choices=["epi", "epc"], default="epi")
    ap.add_argument("--out",    default="out/per_trace_power.png")
    ap.add_argument("--labels", nargs="+", default=None,
                    help="legend labels (one per CSV); default: derived from filename")
    ap.add_argument("--sort",   choices=["trace", "delta"], default="trace")
    ap.add_argument("--width",  type=float, default=18)
    ap.add_argument("--height", type=float, default=6)
    args = ap.parse_args()

    if len(args.csvs) < 2:
        ap.error("need at least 2 CSVs")
    if args.labels and len(args.labels) != len(args.csvs):
        ap.error(f"--labels expects {len(args.csvs)} entries, got {len(args.labels)}")

    # Load all CSVs
    datasets = [load_csv(p, args.metric) for p in args.csvs]
    labels = (args.labels if args.labels
              else [default_label(p) for p in args.csvs])
    predictors = list(zip(labels, datasets))

    # Common traces = intersection
    common = set(datasets[0].keys())
    for d in datasets[1:]:
        common &= set(d.keys())
    if not common:
        print("No common traces across all CSVs", file=sys.stderr)
        sys.exit(1)
    common = sorted(common)
    skipped = sum(len(set(d.keys()) - set(common)) for d in datasets)
    if skipped:
        print(f"Skipping {skipped} non-overlapping (predictor, trace) entries",
              file=sys.stderr)

    if args.sort == "delta":
        first, last = datasets[0], datasets[-1]
        common.sort(key=lambda t: last[t]["baseline"] - first[t]["baseline"],
                    reverse=True)
    traces = common

    units = "fJ/instruction" if args.metric == "epi" else "fJ/prediction-cycle"
    ylabel = f"Energy ({units})"

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    if args.mode == "bar":
        plot_bar(traces, predictors, ax, ylabel)
    else:
        plot_line(traces, predictors, ax, ylabel)
    ax.set_title("Power Breakdown")

    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=130)
    print(f"→ {args.out}  ({len(traces)} traces, "
          f"mode={args.mode}, metric={args.metric}, N={len(predictors)})")

    # Terminal table: per-component mean across traces, one column per predictor
    units = "fJ/inst" if args.metric == "epi" else "fJ/cyc"
    n = len(traces)
    print(f"\n  Mean over {n} traces ({units})")
    w_label = max(max(len(l) for l in labels), 12)
    header = f"  {'Component':<12} " + " ".join(f"{l:>{w_label}}" for l in labels)
    print(header)
    print("  " + "-" * 12 + " " + " ".join("-" * w_label for _ in labels))
    mean_rows = {}
    for key, prn in [("baseline", "Baseline")] + \
                    [(c, COMPONENT_LABELS[c]) for c in COMPONENTS]:
        vals = [sum(d[t][key] for t in traces) / n for _, d in predictors]
        mean_rows[key] = vals
        print(f"  {prn:<12} " + " ".join(f"{v:>{w_label}.1f}" for v in vals))

    print(f"\n  Fraction of baseline mean")
    print(header)
    print("  " + "-" * 12 + " " + " ".join("-" * w_label for _ in labels))
    for key, prn in [("baseline", "Baseline")] + \
                    [(c, COMPONENT_LABELS[c]) for c in COMPONENTS]:
        vals = [
            mean_rows[key][i] / mean_rows["baseline"][i]
            if mean_rows["baseline"][i] else 0.0
            for i in range(len(predictors))
        ]
        print(f"  {prn:<12} " + " ".join(f"{v:>{w_label}.4f}" for v in vals))

    # Ratios vs first predictor
    if len(predictors) > 1:
        print(f"\n  Ratios vs {labels[0]}:")
        for (label, data), idx in [(p, i) for i, p in enumerate(predictors[1:], 1)]:
            row = f"    {label:<{w_label}}  "
            for key, prn in [("baseline", "Baseline")] + \
                            [(c, COMPONENT_LABELS[c]) for c in COMPONENTS]:
                v0 = sum(predictors[0][1][t][key] for t in traces) / n
                v  = sum(data[t][key] for t in traces) / n
                row += f"{prn[:4]} ×{v/v0 if v0 else 0:.2f}  "
            print(row)


if __name__ == "__main__":
    main()
