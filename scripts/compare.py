#!/usr/bin/env python3
"""Bar chart comparing per-trace MPKI of two predictor configurations.

Usage: python3 scripts/compare.py <dir1> <name1> <dir2> <name2> [--numbered]
  e.g. python3 scripts/compare.py full_out/hc_ir_fo_fix "RABT (ahead)" full_out/td_hc2 "TD_HC_2 (direct)"

  --numbered    Replace x-axis trace names with 1..N (and print the legend to stdout).
                Useful when comparing many traces — readable labels are impossible.
  --pct-threshold PCT
                With --pct-vs-baseline, also report stats after dropping traces
                whose percent MPKI increase is greater than PCT.
  --delta-hist   Save a histogram of raw ΔMPKI = MPKI_B - MPKI_A.
"""

import argparse
import math
import os
import sys
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

def short_trace_label(name, max_len=24):
    """Compact trace labels without dropping generic IDs like int_117/web_99."""
    base = name[:-2] if name.endswith("_0") else name
    family = base.split("_", 1)[0]
    keep_numeric_suffix = {"int", "web", "fp", "infra", "compress"}
    parts = base.rsplit("_", 1)
    if (len(parts) == 2 and parts[1].isdigit()
            and family not in keep_numeric_suffix):
        base = parts[0]
    return base[:max_len]

def delta_hist_bins(deltas):
    """Fixed 1-MPKI bin width, integer edges that cover [min, max]."""
    if not deltas:
        return [0, 1]
    lo = math.floor(min(deltas))
    hi = math.ceil(max(deltas))
    if hi <= lo:
        hi = lo + 1
    return [lo + i for i in range(hi - lo + 1)]

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
ap.add_argument("--pct-threshold", type=float, default=None,
                help="with --pct-vs-baseline, drop traces with %% increase greater "
                     "than this threshold from an additional filtered summary")
ap.add_argument("--delta-hist", action="store_true",
                help="additionally plot a histogram of raw ΔMPKI = MPKI_B - MPKI_A")
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
    deltas = [mpki2[t] - mpki1[t] for t in common]
    robust_delta = sum(deltas) / len(deltas)
    robust_base = sum(xs) / len(xs)
    robust_pct = robust_delta / robust_base * 100 if robust_base else 0.0
    fig2, ax2 = plt.subplots(figsize=(9, 6))
    ax2.scatter(xs, ys, s=30, alpha=0.7, color="#FF9800", edgecolor="black", linewidth=0.4)
    ax2.axhline(0, color="black", linewidth=0.8, alpha=0.6)
    # Reference line: aggregate percent from mean MPKI delta, not mean per-trace percent.
    ax2.axhline(robust_pct, color="#2196F3", linewidth=1.2, linestyle="--",
                label=f"avg ΔMPKI / avg MPKI = {robust_pct:+.1f}%")
    ax2.set_xlabel(f"Baseline MPKI ({name1})", fontsize=14)
    ax2.set_ylabel(f"Change vs baseline (%)  ({name2} − {name1}) / {name1}", fontsize=14)
    ax2.set_xscale("log")
    ax2.grid(True, alpha=0.3, which="both")
    ax2.legend(fontsize=11)
    plt.tight_layout()
    out2 = "compare_pct_vs_baseline.png"
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    print(f"\nSaved {out2}")
    print(f"  avg baseline MPKI: {robust_base:.3f}")
    print(f"  avg ΔMPKI across {len(common)} traces: {robust_delta:+.3f}")
    print(f"  avg ΔMPKI / avg baseline MPKI: {robust_pct:+.2f}%")
    n_better = sum(1 for y in ys if y < 0)
    n_worse  = sum(1 for y in ys if y > 0)
    print(f"  {name2} better on {n_better}/{len(common)} traces, worse on {n_worse}")

    if args.pct_threshold is not None:
        kept = [(t, x, y) for t, x, y in zip(common, xs, ys)
                if y <= args.pct_threshold]
        dropped = [(t, x, y) for t, x, y in zip(common, xs, ys)
                   if y > args.pct_threshold]
        if not kept:
            print(f"  threshold Δ% <= {args.pct_threshold:g}: kept 0/{len(common)} traces")
        else:
            kept_ys = [y for _, _, y in kept]
            kept_base = sum(x for _, x, _ in kept) / len(kept)
            kept_delta = sum(mpki2[t] - mpki1[t] for t, _, _ in kept) / len(kept)
            kept_pct = kept_delta / kept_base * 100 if kept_base else 0.0
            kept_better = sum(1 for y in kept_ys if y < 0)
            kept_worse = sum(1 for y in kept_ys if y > 0)
            print(f"  threshold Δ% <= {args.pct_threshold:g}: "
                  f"kept {len(kept)}/{len(common)}, dropped {len(dropped)}")
            print(f"  filtered avg baseline MPKI: {kept_base:.3f}")
            print(f"  filtered avg ΔMPKI: {kept_delta:+.3f}")
            print(f"  filtered avg ΔMPKI / avg baseline MPKI: {kept_pct:+.2f}%")
            print(f"  filtered {name2} better on {kept_better}/{len(kept)} traces, "
                  f"worse on {kept_worse}")
            if dropped:
                worst = sorted(dropped, key=lambda r: r[2], reverse=True)[:10]
                print("  dropped worst traces:")
                for t, x, y in worst:
                    print(f"    {t:<45} baseline={x:7.3f} MPKI  Δ%={y:+8.2f}")

if args.delta_hist:
    deltas = [mpki2[t] - mpki1[t] for t in common]
    bins = delta_hist_bins(deltas)
    fig3, ax3 = plt.subplots(figsize=(9, 6))
    counts, bin_edges, _ = ax3.hist(deltas, bins=bins, color="#FF9800", alpha=0.75,
                                    edgecolor="black", linewidth=0.5)
    ax3.axvline(0, color="black", linewidth=1.0, alpha=0.7)
    ax3.axvline(diff, color="#2196F3", linewidth=1.5, linestyle="--",
                label=f"mean ΔMPKI = {diff:+.2f}")
    ax3.set_xlabel(f"ΔMPKI ({name2} vs {name1})", fontsize=18)
    ax3.set_ylabel("Trace count", fontsize=18)
    ax3.tick_params(labelsize=14)
    ax3.grid(axis="y", alpha=0.3)
    ax3.legend(prop={"weight": "bold", "size": 16})

    sparse_ann = []
    for i, c in enumerate(counts):
        lo, hi = bin_edges[i], bin_edges[i + 1]
        center = 0.5 * (lo + hi)
        if c > 0 and center >= 10:
            traces = [
                t for t, d in zip(common, deltas)
                if (lo <= d < hi) or (i == len(counts) - 1 and lo <= d <= hi)
            ]
            if traces:
                sparse_ann.append((i, center, traces))
    if sparse_ann:
        ymax = max(counts) if len(counts) else 1
        y_top = ymax * 1.15  # modest headroom; annotations should not push to the top
        ax3.set_ylim(top=y_top)
        # Keep any single sparse-bin label block short; long bins are summarized
        # so the stacked offsets remain readable.
        max_lines = 9
        # Sort by center x for left-to-right placement, then place each label as
        # low as possible above its own bar. Labels are raised only when their
        # estimated boxes would overlap a previously placed label.
        sparse_ann_sorted = sorted(sparse_ann, key=lambda t: t[1])
        fig3.canvas.draw()
        renderer = fig3.canvas.get_renderer()
        min_offset_pts = 2.0
        label_gap_pts = 3.0
        distinct_gap_pts = 0.5
        font_size = 12
        placed = []
        for bin_idx, center, traces in sparse_ann_sorted:
            label_lines = [short_trace_label(t) for t in traces]
            if len(label_lines) > max_lines:
                hidden = len(label_lines) - (max_lines - 1)
                label_lines = label_lines[:max_lines - 1] + [f"… +{hidden} more"]
            label = "\n".join(label_lines)
            bar_count = counts[bin_idx]
            y_tip = bar_count
            measure = ax3.text(0, 0, label, fontsize=font_size,
                               ha="center", va="bottom", alpha=0.0)
            bbox = measure.get_window_extent(renderer=renderer)
            measure.remove()
            x_px, y_tip_px = ax3.transData.transform((center, y_tip))
            label_h_px = bbox.height
            label_half_w_px = bbox.width / 2.0
            gap_px = label_gap_pts * fig3.dpi / 72.0
            y_offset = min_offset_pts
            while True:
                label_bottom_px = y_tip_px + y_offset * fig3.dpi / 72.0
                label_top_px = label_bottom_px + label_h_px
                raised = False
                for prev in placed:
                    px, pbottom, ptop, phalf_w, poffset = prev
                    horiz = abs(x_px - px) < (label_half_w_px + phalf_w + gap_px)
                    vert = label_bottom_px < ptop + gap_px and label_top_px > pbottom - gap_px
                    same_offset = abs(y_offset - poffset) < distinct_gap_pts
                    if horiz and vert:
                        label_bottom_px = ptop + gap_px
                        y_offset = (label_bottom_px - y_tip_px) * 72.0 / fig3.dpi
                        raised = True
                        break
                    if same_offset:
                        y_offset = poffset + distinct_gap_pts
                        raised = True
                        break
                if not raised:
                    break
            _, y_tip_px = ax3.transData.transform((center, y_tip))
            label_bottom_px = y_tip_px + y_offset * fig3.dpi / 72.0
            label_top_px = label_bottom_px + label_h_px
            placed.append((x_px, label_bottom_px, label_top_px, label_half_w_px, y_offset))
            ax3.annotate(label, xy=(center, y_tip), xytext=(0, y_offset),
                         textcoords="offset points", annotation_clip=False,
                         ha="center", va="bottom", fontsize=font_size,
                         clip_on=False,
                         arrowprops={"arrowstyle": "-", "color": "black",
                                     "linewidth": 0.5, "alpha": 0.6,
                                     "shrinkA": 0, "shrinkB": 0})

    plt.tight_layout()
    out3 = "compare_delta_mpki_hist.png"
    plt.savefig(out3, dpi=150, bbox_inches="tight")

    ds = sorted(deltas)
    def quantile(p):
        if not ds:
            return 0.0
        idx = round((len(ds) - 1) * p)
        return ds[int(idx)]

    buckets = [
        ("improved < -1", lambda d: d < -1),
        ("near-neutral [-1,+1]", lambda d: -1 <= d <= 1),
        ("small regression (+1,+5]", lambda d: 1 < d <= 5),
        ("large regression (+5,+10]", lambda d: 5 < d <= 10),
        ("severe regression > +10", lambda d: d > 10),
    ]

    print(f"\nSaved {out3}")
    print(f"  ΔMPKI distribution ({name2} - {name1})")
    print(f"  mean:   {diff:+.3f}")
    print(f"  median: {quantile(0.50):+.3f}")
    print(f"  p75:    {quantile(0.75):+.3f}")
    print(f"  p90:    {quantile(0.90):+.3f}")
    print(f"  p95:    {quantile(0.95):+.3f}")
    print(f"  max:    {max(ds):+.3f}")
    print("  buckets:")
    for label, pred in buckets:
        print(f"    {label:<26} {sum(1 for d in deltas if pred(d)):>4}/{len(deltas)}")
