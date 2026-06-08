#!/usr/bin/env python3
"""Compare structural wire-length metrics across predictor floorplans.

Parses one or more HARCOM floorplan.gv files (produced by dump_floorplan.py
or per_trace_power.py) and reports for each:

  - N_RAMs                      raw RAM bank count
  - Σ area (mm²)                total SRAM footprint
  - Bounding-box (mm²)          floorplan envelope
  - Σ all-pair Manhattan (mm)   structural total wire length under HARCOM's
                                inter-RAM connect[] matrix model
  - Σ nearest-neighbor (mm)     local routing pressure
  - Mean / median pair dist     distribution shape

Usage:
  python3 scripts/wire_length_analyzer.py <a.gv> <b.gv> ...
  python3 scripts/wire_length_analyzer.py out/floorplans/*.gv --csv out/wire_length.csv

Notes:
  - HARCOM positions are in micrometers; we report in mm.
  - All-pair Manhattan computed in O(N log N) via sort+prefix-sum (per axis).
  - Per-pair wire energy in HARCOM is proportional to length × activity.
    This script isolates length; activity comes from per_trace_power.py CSVs.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

# Matches: rect42 [label="42\nt2_pred", width=17.58, height=10.05, pos="8.79,71.51!", ...
# Also handles: rect42 [label="42", ...   (no \n suffix)
RECT_RE = re.compile(
    r'rect(\d+)\s*\[label="(\d+)(?:\\n([^"]*))?",'
    r'\s*width=([\d.eE+-]+),'
    r'\s*height=([\d.eE+-]+),'
    r'\s*pos="([\d.eE+-]+),([\d.eE+-]+)!"'
)


def parse_gv(path):
    """Return list of dicts: {id, label, x, y, w, h} per RAM rectangle."""
    rams = []
    with open(path) as f:
        for line in f:
            m = RECT_RE.search(line)
            if not m:
                continue
            ram_id, num, label, w, h, x, y = m.groups()
            rams.append({
                "id":    int(ram_id),
                "label": label or "",
                "x":     float(x),
                "y":     float(y),
                "w":     float(w),
                "h":     float(h),
            })
    return rams


def sum_all_pair_manhattan(coords):
    """O(N log N) sum of |xi-xj|+|yi-yj| over all unordered pairs."""
    def axis_sum(vals):
        v = sorted(vals)
        total = 0.0
        prefix = 0.0
        for i, vi in enumerate(v):
            total += vi * i - prefix
            prefix += vi
        return total
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    return axis_sum(xs) + axis_sum(ys)


def nn_manhattan(coords):
    """Σ over each RAM of nearest-neighbor Manhattan distance. O(N²) — fine for ~400 RAMs."""
    n = len(coords)
    total = 0.0
    for i in range(n):
        xi, yi = coords[i]
        best = float("inf")
        for j in range(n):
            if j == i:
                continue
            d = abs(coords[j][0] - xi) + abs(coords[j][1] - yi)
            if d < best:
                best = d
        total += best
    return total


def analyze(path):
    rams = parse_gv(path)
    if not rams:
        return None
    coords = [(r["x"], r["y"]) for r in rams]
    n = len(rams)

    # Areas in um² → convert to mm² (×1e-6)
    sum_area_um2 = sum(r["w"] * r["h"] for r in rams)
    # bbox: use full rect extent (center ± w/2 etc)
    xmins = [r["x"] - r["w"]/2 for r in rams]
    xmaxs = [r["x"] + r["w"]/2 for r in rams]
    ymins = [r["y"] - r["h"]/2 for r in rams]
    ymaxs = [r["y"] + r["h"]/2 for r in rams]
    bbox_w_um = max(xmaxs) - min(xmins)
    bbox_h_um = max(ymaxs) - min(ymins)
    bbox_um2  = bbox_w_um * bbox_h_um

    # Distances in um → mm (×1e-3)
    sum_pair_um = sum_all_pair_manhattan(coords)
    sum_nn_um   = nn_manhattan(coords)

    n_pairs = n * (n - 1) // 2
    mean_pair_um = sum_pair_um / n_pairs if n_pairs else 0.0

    return {
        "n_rams":       n,
        "sum_area_mm2": sum_area_um2 / 1e6,
        "bbox_mm2":     bbox_um2 / 1e6,
        "bbox_w_mm":    bbox_w_um / 1e3,
        "bbox_h_mm":    bbox_h_um / 1e3,
        "sum_pair_mm":  sum_pair_um / 1e3,
        "sum_nn_mm":    sum_nn_um / 1e3,
        "n_pairs":      n_pairs,
        "mean_pair_um": mean_pair_um,
    }


def predictor_name(path):
    return Path(path).stem


def fmt_row(name, m):
    return (f"  {name:<28} {m['n_rams']:>6}  {m['sum_area_mm2']:>10.4f}"
            f"  {m['bbox_mm2']:>10.4f}  {m['sum_pair_mm']:>12.1f}"
            f"  {m['sum_nn_mm']:>10.1f}  {m['mean_pair_um']:>10.1f}")


def _table_group(label):
    """Group RAMs by the leading table identifier in their label (e.g. 't0', 't1').
    Returns 'misc' if the label doesn't start with t<digits>."""
    if not label:
        return "misc"
    m = re.match(r"^(t\d+)", label)
    return m.group(1) if m else "misc"


def render_bundles(gv_path, out_path, mode="bundle", top_n=200):
    """Render the floorplan with wire bundles overlaid.
    mode:
        bundle  — aggregate by table-prefix; one line per table-pair, thickness ∝ Σ length
        top     — draw the top_n longest individual pairs
        all     — draw every pair (slow, opaque for large N)
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import matplotlib.cm as cm

    rams = parse_gv(gv_path)
    n = len(rams)
    if n == 0:
        print(f"WARNING: nothing to render in {gv_path}", file=sys.stderr)
        return

    # Color per table group
    groups = sorted({_table_group(r["label"]) for r in rams})
    palette = cm.get_cmap("tab20", max(len(groups), 4))
    color_map = {g: palette(i) for i, g in enumerate(groups)}

    fig, ax = plt.subplots(figsize=(12, 12))
    # Draw rectangles
    for r in rams:
        g = _table_group(r["label"])
        rect = Rectangle((r["x"] - r["w"]/2, r["y"] - r["h"]/2),
                         r["w"], r["h"],
                         facecolor=color_map[g], edgecolor="black",
                         linewidth=0.3, alpha=0.65)
        ax.add_patch(rect)

    # Wires
    if mode == "bundle":
        # Aggregate inter-group lengths and group centroids
        group_members = defaultdict(list)
        for r in rams:
            group_members[_table_group(r["label"])].append((r["x"], r["y"]))
        centroids = {g: (sum(x for x, _ in pts)/len(pts),
                         sum(y for _, y in pts)/len(pts))
                     for g, pts in group_members.items()}
        inter = defaultdict(float)  # (g1, g2) -> Σ Manhattan
        for i in range(n):
            gi = _table_group(rams[i]["label"])
            for j in range(i+1, n):
                gj = _table_group(rams[j]["label"])
                if gi == gj:
                    continue
                key = (gi, gj) if gi < gj else (gj, gi)
                inter[key] += (abs(rams[i]["x"]-rams[j]["x"])
                               + abs(rams[i]["y"]-rams[j]["y"]))
        if inter:
            max_len = max(inter.values())
            for (g1, g2), length in inter.items():
                x1, y1 = centroids[g1]
                x2, y2 = centroids[g2]
                lw = 0.5 + 5.0 * (length / max_len)
                ax.plot([x1, x2], [y1, y2], "-",
                        color="black", linewidth=lw, alpha=0.35)
        title_suffix = f"bundle-by-table ({len(inter)} bundles)"
    elif mode == "top":
        pairs = []
        for i in range(n):
            for j in range(i+1, n):
                d = (abs(rams[i]["x"]-rams[j]["x"])
                     + abs(rams[i]["y"]-rams[j]["y"]))
                pairs.append((d, i, j))
        pairs.sort(reverse=True)
        pairs = pairs[:top_n]
        max_d = pairs[0][0] if pairs else 1.0
        for d, i, j in pairs:
            lw = 0.3 + 1.5 * (d / max_d)
            alpha = 0.2 + 0.6 * (d / max_d)
            ax.plot([rams[i]["x"], rams[j]["x"]],
                    [rams[i]["y"], rams[j]["y"]],
                    "-", color="red", linewidth=lw, alpha=alpha)
        title_suffix = f"top-{top_n} longest pairs"
    elif mode == "all":
        for i in range(n):
            for j in range(i+1, n):
                ax.plot([rams[i]["x"], rams[j]["x"]],
                        [rams[i]["y"], rams[j]["y"]],
                        "-", color="black", linewidth=0.1, alpha=0.02)
        title_suffix = f"all {n*(n-1)//2} pairs"
    else:
        title_suffix = "(no wires)"

    # Legend for table groups
    handles = [Rectangle((0, 0), 1, 1, facecolor=color_map[g], edgecolor="black",
                         label=g) for g in groups]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0),
              fontsize=8, ncol=1, frameon=True)

    ax.set_aspect("equal")
    ax.set_xlabel("x (μm)")
    ax.set_ylabel("y (μm)")
    ax.set_title(f"{Path(gv_path).stem} — {title_suffix}")
    ax.autoscale_view()
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"  → {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gv_files", nargs="+")
    ap.add_argument("--csv", default=None, help="optional CSV output path")
    ap.add_argument("--render", action="store_true",
                    help="also render PNG floorplans with wire bundles overlaid")
    ap.add_argument("--render-mode", choices=["bundle", "top", "all"], default="bundle",
                    help="how to draw wires: bundle (inter-table aggregate, default),"
                         " top (longest N pairs), all (every pair, slow)")
    ap.add_argument("--top-n", type=int, default=200,
                    help="for --render-mode top: how many longest pairs to draw")
    ap.add_argument("--render-out-dir", default="out/floorplans",
                    help="directory for rendered PNGs (default: out/floorplans)")
    args = ap.parse_args()

    results = []
    for p in args.gv_files:
        m = analyze(p)
        if m is None:
            print(f"WARNING: no RAMs parsed from {p}", file=sys.stderr)
            continue
        results.append((predictor_name(p), p, m))

    if not results:
        print("No floorplans parsed", file=sys.stderr)
        sys.exit(1)

    # Table to stdout
    print("Wire-length structural comparison")
    print("=" * 110)
    print(f"  {'predictor':<28} {'N_RAMs':>6}  {'area mm²':>10}"
          f"  {'bbox mm²':>10}  {'Σ pair mm':>12}  {'Σ NN mm':>10}  {'mean μm':>10}")
    print("  " + "-" * 108)
    for name, _, m in results:
        print(fmt_row(name, m))
    print()

    # Pairwise ratios using the first row as baseline
    if len(results) > 1:
        base_name, _, base = results[0]
        print(f"  Ratios vs {base_name}:")
        for name, _, m in results[1:]:
            print(f"    {name:<28}  N_RAMs ×{m['n_rams']/base['n_rams']:.2f}"
                  f"  area ×{m['sum_area_mm2']/base['sum_area_mm2']:.2f}"
                  f"  Σ pair ×{m['sum_pair_mm']/base['sum_pair_mm']:.2f}"
                  f"  Σ NN ×{m['sum_nn_mm']/base['sum_nn_mm']:.2f}")

    if args.render:
        print("\nRendering wire-bundle overlays:")
        for _, p, _ in results:
            out_png = Path(args.render_out_dir) / (Path(p).stem + f"_wires_{args.render_mode}.png")
            render_bundles(p, str(out_png), mode=args.render_mode, top_n=args.top_n)

    if args.csv:
        Path(args.csv).parent.mkdir(parents=True, exist_ok=True)
        keys = ["predictor", "path"] + list(results[0][2].keys())
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for name, p, m in results:
                row = {"predictor": name, "path": p, **m}
                w.writerow(row)
        print(f"\nCSV → {args.csv}")


if __name__ == "__main__":
    main()
