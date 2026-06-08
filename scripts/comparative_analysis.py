#!/usr/bin/env python3
"""Comparative analysis plots for predictor power decomposition.

Produces three figures:

  breakdown   Stacked-bar aggregate component breakdown for two predictors,
              extrapolated from per_trace_power.py 20-trace fractions to the
              full 168-trace EPI distribution.

  best-worst  Top-5 best + top-5 worst EPI deltas per trace (full eval),
              shown as per-trace stacked component estimates for 2 or 3
              predictors. Components use 20-trace fractions applied to each
              full-eval per-trace EPI.

  scaling     Three scatter panels of structural metrics vs wiring EPI:
              storage, Σ hub Manhattan, and #wires.

Inputs are inferred by predictor name. Expected file layout:

  out/ptp_{Predictor}_quick.csv            (per_trace_power.py output)
  out/ptp_{Predictor}_panel.txt            (HARCOM panel summary)
  out/floorplans/{Predictor}.gv            (HARCOM floorplan)
  full_out/{shortname}/*.out  OR
  out/full_{shortname}/*.out               (full 168-trace eval CSV)

Usage:
  python3 scripts/comparative_analysis.py [--mode {breakdown,best-worst,scaling,all}]
                                          [--out-dir DIR]
                                          [--name-a NAME] [--name-b NAME]
                                          [--best-worst-names NAME NAME [NAME]]
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np


# ── Predictor manifest ──────────────────────────────────────────────────────
# Hard-coded mapping from internal predictor name → display label and
# locations of available artifacts.
PREDICTORS = {
    "TageDefault": {
        "label":      "example TAGE",
        "ptp_csv":    "out/ptp_TageDefault_quick.csv",
        "panel":      "out/ptp_TageDefault_panel.txt",
        "floorplan":  "out/floorplans/TageDefault.gv",
        "full_eval":  "full_out/default_tage",
    },
    "TageDefaultRABT": {
        "label":      "tage<>@RABT (NUMG=15)",
        "ptp_csv":    "out/ptp_TageDefaultRABT_quick.csv",
        "panel":      "out/ptp_TageDefaultRABT_panel.txt",
        "floorplan":  "out/floorplans/TageDefaultRABT.gv",
        "full_eval":  "full_out/default_tage_rabt",
    },
    "TageAheadHC_IR": {
        "label":      "RABT (MINHIST=8)",
        "ptp_csv":    "out/ptp_TageAheadHC_IR_quick.csv",
        "panel":      "out/ptp_TageAheadHC_IR_panel.txt",
        "floorplan":  "out/floorplans/TageAheadHC_IR.gv",
        "full_eval":  "full_out/hc_ir_fo_fix",
    },
    "TageAheadHC_IR_M2": {
        "label":      "RABT M2",
        "ptp_csv":    "out/ptp_TageAheadHC_IR_M2_quick.csv",
        "panel":      "out/ptp_TageAheadHC_IR_M2_panel.txt",
        "floorplan":  "out/floorplans/TageAheadHC_IR_M2.gv",
        "full_eval":  "out/full_m2",
    },
}


# CSV column indices from cbp.hpp run output
COL_NAME, COL_NINSTR, COL_NPRED, COL_EXTRA, COL_MISP, COL_EPI = 0, 1, 4, 5, 8, 11


# ── Loaders ─────────────────────────────────────────────────────────────────

def load_ptp_csv(path):
    """Returns list of dicts (one per trace) from per_trace_power.py CSV."""
    if not Path(path).exists():
        return None
    return list(csv.DictReader(open(path)))


def component_fractions(ptp_rows):
    """Compute mean fraction of baseline EPI per component across traces."""
    if not ptp_rows:
        return None
    fracs = {"ram_logic": [], "fanout": [], "wiring": [], "other": []}
    for r in ptp_rows:
        b = float(r["epi_baseline"])
        if b <= 0:
            continue
        for k in fracs:
            fracs[k].append(float(r[f"epi_{k}"]) / b)
    return {k: (np.mean(v) if v else 0.0) for k, v in fracs.items()}


def component_fraction_by_trace(ptp_rows):
    """Returns trace_name -> component fractions of baseline EPI."""
    if not ptp_rows:
        return None
    out = {}
    for r in ptp_rows:
        b = float(r["epi_baseline"])
        if b <= 0:
            continue
        out[r["trace"]] = {
            k: float(r[f"epi_{k}"]) / b
            for k in ("ram_logic", "fanout", "wiring", "other")
        }
    return out


def predictor_tone(idx, count):
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


def toned_color(color, pred_idx, pred_count):
    return apply_tone(color, predictor_tone(pred_idx, pred_count))


def load_full_eval(directory):
    """Returns dict trace_name -> {ninstr, npred, extra, misp, epi}."""
    if directory is None or not Path(directory).is_dir():
        return None
    out = {}
    for f in Path(directory).iterdir():
        if not f.name.endswith(".out"):
            continue
        try:
            line = open(f).readline().strip()
        except Exception:
            continue
        parts = line.split(",")
        if len(parts) <= COL_EPI:
            continue
        try:
            out[parts[COL_NAME]] = {
                "ninstr": int(parts[COL_NINSTR]),
                "npred":  int(parts[COL_NPRED]),
                "extra":  int(parts[COL_EXTRA]),
                "misp":   int(parts[COL_MISP]),
                "epi":    float(parts[COL_EPI]),
            }
        except ValueError:
            continue
    return out


def parse_panel(path):
    """Extract panel.print() metrics: storage, transistors, sram_area, dyn_power, sta_power."""
    if not Path(path).exists():
        return {}
    text = open(path).read()
    def grab(pattern, conv=float):
        m = re.search(pattern, text)
        return conv(m.group(1)) if m else None
    return {
        "storage_bits": grab(r"storage \(bits\)\s*:\s*([\d.]+)", int),
        "transistors":  grab(r"transistors\s*:\s*([\d.]+)", int),
        "sram_mm2":     grab(r"SRAM area \(mm2\)\s*:\s*([\d.]+)"),
        "dyn_mW":       grab(r"dynamic power \(mW\)\s*:\s*([\d.]+)"),
        "sta_mW":       grab(r"static power \(mW\)\s*:\s*([\d.]+)"),
    }


def parse_floorplan_metrics(path):
    """Quick parse of .gv to get #RAMs, hub Σ Manhattan, n_wires, mean_hub_um."""
    if not Path(path).exists():
        return {}
    rect_re = re.compile(
        r'rect(\d+)\s*\[label="(\d+)(?:\\n([^"]*))?",'
        r'\s*width=([\d.eE+-]+),\s*height=([\d.eE+-]+),'
        r'\s*pos="([\d.eE+-]+),([\d.eE+-]+)!"')
    rams = []
    for line in open(path):
        m = rect_re.search(line)
        if not m:
            continue
        rams.append({
            "id":    int(m.group(1)),
            "x":     float(m.group(6)),
            "y":     float(m.group(7)),
        })
    if not rams:
        return {}
    hub = next((r for r in rams if r["id"] == 0), rams[0])
    dists = [abs(r["x"]-hub["x"]) + abs(r["y"]-hub["y"]) for r in rams if r["id"] != hub["id"]]
    return {
        "n_rams":      len(rams),
        "n_wires":     len(dists),
        "sum_hub_mm":  sum(dists) / 1e3,
        "mean_hub_um": (sum(dists) / len(dists)) if dists else 0.0,
    }


# ── Fig 1: Aggregate breakdown (extrapolated to 168) ────────────────────────

def fig_breakdown(name_a, name_b, out_path):
    pa, pb = PREDICTORS[name_a], PREDICTORS[name_b]
    rows_a = load_ptp_csv(pa["ptp_csv"])
    rows_b = load_ptp_csv(pb["ptp_csv"])
    full_a = load_full_eval(pa["full_eval"])
    full_b = load_full_eval(pb["full_eval"])
    if not (rows_a and rows_b and full_a and full_b):
        print(f"[breakdown] missing inputs for {name_a} or {name_b}", file=sys.stderr)
        return
    frac_a = component_fractions(rows_a)
    frac_b = component_fractions(rows_b)

    # Instruction-weighted baseline EPI across full 168
    def w_mean_epi(d):
        tot_e = sum(v["epi"] * v["ninstr"] for v in d.values())
        tot_n = sum(v["ninstr"] for v in d.values())
        return tot_e / tot_n if tot_n else 0.0
    base_a = w_mean_epi(full_a)
    base_b = w_mean_epi(full_b)

    # Apply fractions to full-eval baseline → extrapolated per-component
    comps = ["ram_logic", "fanout", "wiring"]
    comp_a = [base_a * frac_a[c] for c in comps]
    comp_b = [base_b * frac_b[c] for c in comps]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    colors = {"ram_logic": "#4C72B0", "fanout": "#DD8452", "wiring": "#55A467"}
    labels = {"ram_logic": "RAM + Logic", "fanout": "Fanout", "wiring": "Wiring"}
    x = [0, 1]
    width = 0.6
    bottoms = [0.0, 0.0]
    for c, va, vb in zip(comps, comp_a, comp_b):
        ax.bar(x, [va, vb], width, bottom=bottoms,
               color=colors[c], label=labels[c],
               edgecolor="black", linewidth=0.5)
        bottoms[0] += va
        bottoms[1] += vb
    # Total annotations
    for xi, total in zip(x, bottoms):
        ax.text(xi, total + total*0.02, f"{total:.0f}\n fJ/inst",
                ha="center", fontsize=10, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([pa["label"], pb["label"]], fontsize=12)
    ax.set_ylabel("EPI (fJ / correct-path instruction)", fontsize=11)
    ax.set_title(f"Aggregate EPI breakdown — 168-trace extrapolation\n"
                 f"(component fractions from 20-trace ptp, baseline from full eval)",
                 fontsize=11)
    ax.legend(loc="upper left")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"  → {out_path}")


# ── Fig 2: Best/worst per-trace component breakdowns ────────────────────────

def fig_best_worst(names, top_n, out_path):
    if len(names) < 2 or len(names) > 3:
        print("[best-worst] pass 2 or 3 predictors", file=sys.stderr)
        return
    preds = []
    for name in names:
        p = PREDICTORS[name]
        full = load_full_eval(p["full_eval"])
        ptp_rows = load_ptp_csv(p["ptp_csv"])
        trace_fracs = component_fraction_by_trace(ptp_rows)
        mean_fracs = component_fractions(ptp_rows)
        if not (full and trace_fracs and mean_fracs):
            print(f"[best-worst] missing inputs for {name}", file=sys.stderr)
            return
        preds.append({
            "name": name,
            "label": p["label"],
            "full": full,
            "trace_fracs": trace_fracs,
            "mean_fracs": mean_fracs,
        })

    common = set(preds[0]["full"].keys())
    for p in preds[1:]:
        common &= set(p["full"].keys())
    if not common:
        print("[best-worst] no common full-eval traces", file=sys.stderr)
        return

    # Δ = last - first (positive: last predictor uses more energy)
    first, last = preds[0], preds[-1]
    ranked = []
    for t in common:
        ranked.append({
            "trace": t,
            "delta": last["full"][t]["epi"] - first["full"][t]["epi"],
        })
    ranked.sort(key=lambda r: r["delta"])
    selected = ranked[:top_n] + ranked[-top_n:]
    selected_traces = [r["trace"] for r in selected]
    split_at = len(ranked[:top_n])

    comps = ["ram_logic", "fanout", "wiring"]
    colors = {"ram_logic": "#4C72B0", "fanout": "#DD8452", "wiring": "#55A467"}
    labels = {"ram_logic": "RAM + Logic", "fanout": "Fanout", "wiring": "Wiring"}

    def estimate_components(pred, trace):
        fracs = pred["trace_fracs"].get(trace, pred["mean_fracs"])
        total = pred["full"][trace]["epi"]
        vals = [total * fracs[c] for c in comps]
        frac_sum = sum(fracs[c] for c in comps)
        if frac_sum > 0:
            scale = total / sum(vals) if sum(vals) else 0.0
            vals = [v * scale for v in vals]
        return vals

    fig, ax = plt.subplots(figsize=(14, 6.5))
    x = np.arange(len(selected_traces))
    group_w = 0.82
    bar_w = group_w / len(preds)
    offsets = [(-group_w / 2) + (i + 0.5) * bar_w for i in range(len(preds))]

    for i, pred in enumerate(preds):
        bottoms = np.zeros(len(selected_traces))
        for comp in comps:
            vals = np.array([
                estimate_components(pred, t)[comps.index(comp)]
                for t in selected_traces
            ])
            ax.bar(x + offsets[i], vals, bar_w, bottom=bottoms,
                   color=toned_color(colors[comp], i, len(preds)),
                   edgecolor="black", linewidth=0.35,
                   label=labels[comp] if i == 0 else None)
            bottoms += vals
        for xi, total in zip(x + offsets[i], bottoms):
            ax.text(xi, total + max(bottoms) * 0.015, f"{total:.0f}",
                    ha="center", va="bottom", fontsize=7, rotation=90)

    if split_at and split_at < len(selected_traces):
        ax.axvline(split_at - 0.5, color="black", linewidth=0.8, alpha=0.55)
        ymax = ax.get_ylim()[1]
        ax.text((split_at - 1) / 2, ymax * 0.98,
                f"Top {top_n}: improved / increased least",
                ha="center", va="top", fontsize=9, fontweight="bold")
        ax.text(split_at + (len(selected_traces) - split_at - 1) / 2, ymax * 0.98,
                f"Top {top_n}: increased most",
                ha="center", va="top", fontsize=9, fontweight="bold")

    trace_labels = [t[:24] for t in selected_traces]
    ax.set_xticks(x)
    ax.set_xticklabels(trace_labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Estimated EPI breakdown (fJ / instruction)", fontsize=10)
    ax.set_title("Per-trace power breakdown for best/worst full-eval EPI deltas\n"
                 f"ranked by {last['label']} − {first['label']}",
                 fontsize=11)
    ax.grid(axis="y", alpha=0.3)

    comp_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=colors[c], edgecolor="black",
                      linewidth=0.35, label=labels[c])
        for c in comps
    ]
    comp_legend = ax.legend(handles=comp_handles, title="Power category",
                            loc="upper left", bbox_to_anchor=(1.01, 1.0),
                            fontsize=8, title_fontsize=8)
    ax.add_artist(comp_legend)

    pred_handles = [
        plt.Rectangle((0, 0), 1, 1,
                      facecolor=toned_color("#666666", i, len(preds)),
                      edgecolor="black",
                      linewidth=0.35, label=p["label"])
        for i, p in enumerate(preds)
    ]
    ax.legend(handles=pred_handles, title="Predictor",
              loc="upper left", bbox_to_anchor=(1.01, 0.78),
              fontsize=8, title_fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  → {out_path}")


# ── Fig 3: Wiring-mass growth — more wires drive more total length ──────────

def fig_scaling(names, out_path):
    rows = []
    for n in names:
        p = PREDICTORS[n]
        fp = parse_floorplan_metrics(p["floorplan"])
        if not fp:
            continue
        rows.append({
            "name":         n,
            "label":        p["label"],
            "n_wires":      fp["n_wires"],
            "sum_hub_mm":   fp["sum_hub_mm"],
            "mean_hub_um":  fp["mean_hub_um"],
        })
    if len(rows) < 2:
        print(f"[scaling] need ≥2 predictors with floorplan data", file=sys.stderr)
        return

    # Predictor color
    palette = plt.cm.tab10.colors
    pred_color = {r["label"]: palette[i % len(palette)] for i, r in enumerate(rows)}

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(13, 5.5))
    labels = [r["label"] for r in rows]
    x = np.arange(len(rows))
    colors = [pred_color[l] for l in labels]

    # Left: number of wires (the cause)
    n_wires = [r["n_wires"] for r in rows]
    ax_left.bar(x, n_wires, color=colors, edgecolor="black", linewidth=0.5)
    for xi, v in zip(x, n_wires):
        ax_left.text(xi, v + max(n_wires)*0.015, str(v),
                     ha="center", fontsize=10, fontweight="bold")
    ax_left.set_xticks(x)
    ax_left.set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
    ax_left.set_ylabel("# wires  (≈ #RAMs − 1)", fontsize=11)
    ax_left.set_title("Cause: banking + sec_tag = more RAMs = more wires",
                      fontsize=11)
    ax_left.grid(axis="y", alpha=0.3)

    # Right: Σ hub mm (the effect), annotated with mean wire length
    sum_mm = [r["sum_hub_mm"] for r in rows]
    mean_um = [r["mean_hub_um"] for r in rows]
    ax_right.bar(x, sum_mm, color=colors, edgecolor="black", linewidth=0.5)
    for xi, v, mu in zip(x, sum_mm, mean_um):
        ax_right.text(xi, v + max(sum_mm)*0.015,
                      f"{v:.1f} mm\n({mu:.0f} μm × {n_wires[xi]})",
                      ha="center", fontsize=9)
    ax_right.set_xticks(x)
    ax_right.set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
    ax_right.set_ylabel("Σ RAM → register-hub Manhattan (mm)", fontsize=11)
    ax_right.set_title("Effect: total wire length grows with wire count",
                       fontsize=11)
    ax_right.grid(axis="y", alpha=0.3)

    fig.suptitle("Why RABT pays more for wiring:  more banked RAMs → more wires → more total length",
                 fontsize=12, y=1.00)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  → {out_path}")


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["breakdown", "best-worst", "scaling", "all"],
                    default="all")
    ap.add_argument("--name-a", default="TageDefault",
                    help="primary predictor A (breakdown / best-worst)")
    ap.add_argument("--name-b", default="TageAheadHC_IR_M2",
                    help="primary predictor B (breakdown / best-worst)")
    ap.add_argument("--best-worst-names", nargs="+", default=None,
                    help="2 or 3 predictors for figure 2; first is baseline, "
                         "last is used for trace ranking")
    ap.add_argument("--scaling-names", nargs="+",
                    default=["TageDefault", "TageDefaultRABT",
                             "TageAheadHC_IR", "TageAheadHC_IR_M2"],
                    help="predictors to include in scaling fig")
    ap.add_argument("--top-n", type=int, default=5,
                    help="best/worst trace count per side")
    ap.add_argument("--out-dir", default="out/comparative")
    args = ap.parse_args()

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    if args.mode in ("breakdown", "all"):
        print("Fig 1: breakdown")
        fig_breakdown(args.name_a, args.name_b,
                      str(Path(args.out_dir) / "fig1_breakdown.png"))
    if args.mode in ("best-worst", "all"):
        print(f"Fig 2: best/worst top-{args.top_n} traces")
        best_worst_names = args.best_worst_names or [args.name_a, args.name_b]
        fig_best_worst(best_worst_names, args.top_n,
                       str(Path(args.out_dir) / "fig2_best_worst.png"))
    if args.mode in ("scaling", "all"):
        print("Fig 3: scaling")
        fig_scaling(args.scaling_names,
                    str(Path(args.out_dir) / "fig3_scaling.png"))


if __name__ == "__main__":
    main()
