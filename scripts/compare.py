#!/usr/bin/env python3
"""Bar chart comparing per-trace MPKI of two predictor configurations.

Usage: python3 scripts/compare.py <dir1> <name1> <dir2> <name2> [--numbered]
  e.g. python3 scripts/compare.py full_out/hc_ir_fo_fix "RABT (ahead)" full_out/td_hc2 "TD_HC_2 (direct)"

  --numbered    Replace x-axis trace names with 1..N (and print the legend to stdout).
                Useful when comparing many traces — readable labels are impossible.
"""

import argparse
import os
import matplotlib.pyplot as plt

# CSV columns from cbp.hpp (0-indexed):
# name,ninstr,nbranch,ncondbr,npred,extra_cycles,disagree,disagree_at_end,mispredictions,p1_lat,p2_lat,epi
COL_NAME, COL_NINSTR, COL_MISP = 0, 1, 8

def load_mpki(directory):
    """Load per-trace MPKI from .out files. Returns dict {trace_name: mpki}."""
    results = {}
    for fname in os.listdir(directory):
        if not fname.endswith(".out"):
            continue
        with open(os.path.join(directory, fname)) as f:
            line = f.readline().strip()
            if not line:
                continue
            parts = line.split(",")
            results[parts[COL_NAME]] = float(parts[COL_MISP]) / float(parts[COL_NINSTR]) * 1000
    return results

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("dir1")
ap.add_argument("name1")
ap.add_argument("dir2")
ap.add_argument("name2")
ap.add_argument("--numbered", action="store_true",
                help="replace x-axis labels with 1..N and print legend to stdout")
ap.add_argument("--pct-vs-baseline", action="store_true",
                help="additionally plot %% MPKI change per trace vs baseline MPKI (scatter)")
args = ap.parse_args()
dir1, name1, dir2, name2 = args.dir1, args.name1, args.dir2, args.name2

mpki1 = load_mpki(dir1)
mpki2 = load_mpki(dir2)

# Common traces, sorted by MPKI difference (pred2 - pred1)
common = sorted(set(mpki1) & set(mpki2), key=lambda t: mpki2[t] - mpki1[t])
if not common:
    print("No common traces found.")
    sys.exit(1)

m1 = [mpki1[t] for t in common]
m2 = [mpki2[t] for t in common]
avg1 = sum(m1) / len(m1)
avg2 = sum(m2) / len(m2)
diff = avg2 - avg1
sign = "+" if diff >= 0 else ""

fig, ax = plt.subplots(figsize=(max(12, len(common) * 0.15), 5))

x = range(len(common))
w = 0.4
ax.bar([i - w/2 for i in x], m1, w, label=f"{name1} ({avg1:.3f})", color="#2196F3", alpha=0.8)
ax.bar([i + w/2 for i in x], m2, w, label=f"{name2} ({avg2:.3f}, {sign}{diff:.3f})", color="#FF9800", alpha=0.8)
ax.set_ylabel("MPKI", fontsize=18)
ax.tick_params(labelsize=14)
leg = ax.legend(prop={"weight": "bold", "size": 16})
ax.grid(axis="y", alpha=0.3)

ax.set_xticks(list(x))
if args.numbered:
    labels = [str(i + 1) for i in x]
    ax.set_xticklabels(labels, rotation=0, ha="center", fontsize=10)
    ax.set_xlabel("Trace index", fontsize=16)
else:
    # Shorten: strip digits after last _ or ., keep first 15 chars
    short = [t.split("_")[0].split(".")[0][:15] for t in common]
    ax.set_xticklabels(short, rotation=45, ha="right", fontsize=18)

plt.tight_layout()
out = "compare_mpki.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved {out}")
print(f"  {len(common)} common traces")
print(f"  {name1}: {avg1:.3f} MPKI")
print(f"  {name2}: {avg2:.3f} MPKI")
print(f"  delta:  {sign}{diff:.3f} ({sign}{diff/avg1*100:.1f}%)")

if args.numbered:
    print(f"\nTrace index legend ({len(common)} traces, sort order: Δ ascending):")
    w = len(str(len(common)))
    for i, t in enumerate(common, 1):
        d = mpki2[t] - mpki1[t]
        s = "+" if d >= 0 else ""
        print(f"  {i:>{w}}. {t:<45}  {name1}={mpki1[t]:7.3f}  "
              f"{name2}={mpki2[t]:7.3f}  Δ={s}{d:7.3f}")

if args.pct_vs_baseline:
    # Scatter: x = baseline (name1) MPKI, y = % change ((mpki2 - mpki1) / mpki1 * 100)
    xs = [mpki1[t] for t in common]
    ys = [(mpki2[t] - mpki1[t]) / mpki1[t] * 100 if mpki1[t] > 0 else 0 for t in common]
    fig2, ax2 = plt.subplots(figsize=(9, 6))
    ax2.scatter(xs, ys, s=30, alpha=0.7, color="#FF9800", edgecolor="black", linewidth=0.4)
    ax2.axhline(0, color="black", linewidth=0.8, alpha=0.6)
    # Reference line: avg % change
    mean_pct = sum(ys) / len(ys)
    ax2.axhline(mean_pct, color="#2196F3", linewidth=1.2, linestyle="--",
                label=f"mean Δ% = {mean_pct:+.1f}%")
    ax2.set_xlabel(f"Baseline MPKI ({name1})", fontsize=14)
    ax2.set_ylabel(f"Change vs baseline (%)  ({name2} − {name1}) / {name1}", fontsize=14)
    ax2.set_xscale("log")
    ax2.grid(True, alpha=0.3, which="both")
    ax2.legend(fontsize=11)
    plt.tight_layout()
    out2 = "compare_pct_vs_baseline.png"
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    print(f"\nSaved {out2}")
    print(f"  mean Δ% across {len(common)} traces: {mean_pct:+.2f}%")
    n_better = sum(1 for y in ys if y < 0)
    n_worse  = sum(1 for y in ys if y > 0)
    print(f"  {name2} better on {n_better}/{len(common)} traces, worse on {n_worse}")
