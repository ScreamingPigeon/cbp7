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
    n = len(rams)

    # Areas in um² → mm² (×1e-6)
    sum_area_um2 = sum(r["w"] * r["h"] for r in rams)
    xmins = [r["x"] - r["w"]/2 for r in rams]
    xmaxs = [r["x"] + r["w"]/2 for r in rams]
    ymins = [r["y"] - r["h"]/2 for r in rams]
    ymaxs = [r["y"] + r["h"]/2 for r in rams]
    bbox_w_um = max(xmaxs) - min(xmins)
    bbox_h_um = max(ymaxs) - min(ymins)
    bbox_um2  = bbox_w_um * bbox_h_um

    # Per-RAM → register-hub (RAM 0) Manhattan distance.
    # Since all registers site at the first RAM (id=0), this captures the
    # structural cost of any register-to-RAM transport in HARCOM's model.
    hub = next((r for r in rams if r["id"] == 0), rams[0])
    hub_dists_um = [abs(r["x"] - hub["x"]) + abs(r["y"] - hub["y"])
                    for r in rams if not (r["id"] == hub["id"])]
    n_wires    = len(hub_dists_um)
    sum_hub_um = sum(hub_dists_um)
    mean_um    = sum_hub_um / n_wires if n_wires else 0.0
    sorted_d   = sorted(hub_dists_um)
    median_um  = sorted_d[n_wires // 2] if n_wires else 0.0
    max_um     = sorted_d[-1] if n_wires else 0.0
    min_um     = sorted_d[0] if n_wires else 0.0
    # 90th / 99th percentile
    def pct(p):
        idx = max(0, min(n_wires - 1, int(round(p/100.0 * (n_wires - 1)))))
        return sorted_d[idx] if n_wires else 0.0
    p90 = pct(90); p99 = pct(99)

    return {
        "n_rams":       n,
        "n_wires":      n_wires,
        "sum_area_mm2": sum_area_um2 / 1e6,
        "bbox_mm2":     bbox_um2 / 1e6,
        "bbox_w_mm":    bbox_w_um / 1e3,
        "bbox_h_mm":    bbox_h_um / 1e3,
        "sum_hub_mm":   sum_hub_um / 1e3,
        "mean_hub_um":  mean_um,
        "median_hub_um": median_um,
        "p90_hub_um":   p90,
        "p99_hub_um":   p99,
        "max_hub_um":   max_um,
        "min_hub_um":   min_um,
        "_hub_dists_um": hub_dists_um,  # for plotting
    }


def predictor_name(path):
    return Path(path).stem


# ── Site / register-placement analysis ──────────────────────────────────────

def load_sites(sites_path):
    """Parse sites.txt: returns list of (name, site_id)."""
    out = []
    with open(sites_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.rsplit(None, 1)
            if len(parts) != 2:
                continue
            try:
                out.append((parts[0], int(parts[1])))
            except ValueError:
                continue
    return out


# Map per-table register name → which RAM(s) it would naturally live near.
# Returns (table_index, role) or None if not a per-table register.
PER_TABLE_REG_RE = re.compile(r"^(\w+?)\[(\d+)\]")


# Family → ("rabt" | "tage", role-suffix-or-label)
#   - RABT: per-table label is 't<i>_<role>' (e.g. t3_pred, t3_tag, t3_sec)
#   - tage<>: per-table is array indexing into a single labeled head
#             (e.g. gindex[3] → 3rd entry of the "tags" array)
_FAMILIES = {
    # ---- RABT ----
    "prefetch_pred":     ("rabt", "pred"),
    "prefetch_hyst":     ("rabt", "hyst"),
    "prefetch_u":        ("rabt", "u"),
    "prefetch_tag_hit":  ("rabt", "tag"),
    "prefetch_ctag":     ("rabt", "tag"),
    "prefetch_idx":      ("rabt", "tag"),
    "prefetch_sec":      ("rabt", "sec"),
    "prefetch_group_id": ("rabt", "tag"),
    # ---- tage<> (per-table registers indexed by table id, NUMG-wide) ----
    "gindex":            ("tage", "tags"),
    "htag":              ("tage", "tags"),
    "readt":             ("tage", "tags"),
    "readc":             ("tage", "gpred"),
    "readh":             ("tage", "ghyst"),
    "readu":             ("tage", "u"),
}


def _owning_ram_role(reg_name):
    """Return (kind, table_idx, role_or_label) or None.
    kind: 'rabt' (RABT-style 't<i>_<role>' label lookup) or
          'tage' (tage<>-style — use labeled array head as proxy)."""
    m = PER_TABLE_REG_RE.match(reg_name)
    if not m:
        return None
    family, idx = m.group(1), int(m.group(2))
    fam = _FAMILIES.get(family)
    if fam is None:
        return None
    return (fam[0], idx, fam[1])


def _build_label_index(rams):
    """Return dict label_suffix -> (x, y, rect_id) for fast lookup."""
    out = {}
    for r in rams:
        if r["label"]:
            out[r["label"]] = (r["x"], r["y"], r["id"])
    return out


def analyze_sites(gv_path, sites_path):
    """Returns dict with:
        site_histogram: {site_id: count}
        hub_site:       most-common site_id (the register hub)
        relevant_rows:  list of dicts with current vs counterfactual distance per per-table register
        total_current_um:  Σ Manhattan from hub to natural-owner over relevant regs
        total_counterfactual_um:  Σ over relevant regs (≈ 0, since intra-table)
    """
    rams = parse_gv(gv_path)
    label_idx = _build_label_index(rams)
    sites = load_sites(sites_path)

    # Site histogram
    from collections import Counter
    hist = Counter(s for _, s in sites)
    hub_site = max(hist, key=hist.get) if hist else 0
    # Coordinate of the hub RAM
    hub_coord = next(((r["x"], r["y"]) for r in rams if r["id"] == hub_site), None)

    rows = []
    total_current = 0.0
    total_counter = 0.0
    for name, site in sites:
        owner = _owning_ram_role(name)
        if owner is None:
            continue
        kind, idx, role = owner
        if kind == "rabt":
            target_label = f"t{idx}_{role}"
        else:  # 'tage' — labeled array head as proxy for i-th member
            target_label = role
        if target_label not in label_idx:
            continue
        tx, ty, _ = label_idx[target_label]
        # Current: site → owner-RAM
        sx, sy = next(((r["x"], r["y"]) for r in rams if r["id"] == site), (None, None))
        if sx is None:
            continue
        d_current = abs(sx - tx) + abs(sy - ty)
        # Counterfactual: register sites at the owner RAM → 0 distance
        d_counter = 0.0
        total_current += d_current
        total_counter += d_counter
        rows.append({
            "register":  name,
            "site":      site,
            "target":    target_label,
            "d_current_um": d_current,
            "d_counter_um": d_counter,
        })

    return {
        "site_histogram":  dict(hist),
        "hub_site":        hub_site,
        "hub_coord":       hub_coord,
        "relevant_rows":   rows,
        "total_current_um": total_current,
        "total_counterfactual_um": total_counter,
    }


def fmt_row(name, m):
    return (f"  {name:<28} {m['n_wires']:>7}  {m['sum_area_mm2']:>9.4f}"
            f"  {m['sum_hub_mm']:>10.2f}  {m['mean_hub_um']:>9.1f}"
            f"  {m['median_hub_um']:>9.1f}  {m['p90_hub_um']:>8.1f}"
            f"  {m['p99_hub_um']:>8.1f}  {m['max_hub_um']:>8.1f}")


def _table_group(label):
    """Group RAMs by the leading table identifier in their label (e.g. 't0', 't1').
    Returns 'misc' if the label doesn't start with t<digits>."""
    if not label:
        return "misc"
    m = re.match(r"^(t\d+)", label)
    return m.group(1) if m else "misc"


def render_bundles(gv_path, out_path, mode="bundle", top_n=200, hub_info=None):
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
        # Starburst: draw a line from every RAM to the register hub.
        # If no hub info passed, fall back to RAM 0 by coordinate of first rect.
        if hub_info and hub_info.get("hub_coord"):
            hx, hy = hub_info["hub_coord"]
        else:
            hx, hy = rams[0]["x"], rams[0]["y"]

        # Color by Manhattan distance (longer = more saturated)
        max_d = max(abs(r["x"]-hx) + abs(r["y"]-hy) for r in rams) or 1.0
        n_wires = 0
        for r in rams:
            if r["x"] == hx and r["y"] == hy:
                continue  # skip the hub itself
            d = abs(r["x"]-hx) + abs(r["y"]-hy)
            frac = d / max_d
            ax.plot([hx, r["x"]], [hy, r["y"]], "-",
                    color=(1.0, 0.4 + 0.3*(1-frac), 0.2),  # orange→red by distance
                    linewidth=0.4 + 1.0*frac, alpha=0.4 + 0.3*frac,
                    zorder=2)
            n_wires += 1
        title_suffix = f"RAM ↔ register-hub wires ({n_wires})"
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

    # Highlight register-hub RAM
    if hub_info and hub_info.get("hub_coord"):
        hx, hy = hub_info["hub_coord"]
        n_regs = sum(hub_info["site_histogram"].values())
        ax.plot(hx, hy, "*", color="yellow", markersize=24,
                markeredgecolor="black", markeredgewidth=1.5,
                label=f"Register hub: RAM {hub_info['hub_site']} ({n_regs} regs)",
                zorder=10)

    # Legend for table groups + hub marker
    handles = [Rectangle((0, 0), 1, 1, facecolor=color_map[g], edgecolor="black",
                         label=g) for g in groups]
    if hub_info and hub_info.get("hub_coord"):
        from matplotlib.lines import Line2D
        n_regs = sum(hub_info["site_histogram"].values())
        handles.append(Line2D([0], [0], marker="*", color="w",
                              markerfacecolor="yellow",
                              markeredgecolor="black", markersize=14,
                              label=f"hub ({n_regs} regs)"))
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
    ap.add_argument("--sites", action="store_true",
                    help="load {name}.sites.txt next to each .gv; mark register hub"
                         " on the floorplan; report current vs counterfactual"
                         " distance for per-table prefetch_*[i] registers")
    ap.add_argument("--plot-dist", action="store_true",
                    help="overlay CDF + histogram of RAM→hub Manhattan distances "
                         "for all input predictors → out/floorplans/wire_dist.png")
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

    # Table to stdout: per-RAM → register-hub wire-length distribution
    print("Wire-length distribution (RAM → register-hub)")
    print("=" * 110)
    print(f"  {'predictor':<28} {'#wires':>7}  {'area mm²':>9}"
          f"  {'Σ hub mm':>10}  {'mean μm':>9}  {'median μm':>9}"
          f"  {'p90 μm':>8}  {'p99 μm':>8}  {'max μm':>8}")
    print("  " + "-" * 108)
    for name, _, m in results:
        print(fmt_row(name, m))
    print()

    # Pairwise ratios using the first row as baseline
    if len(results) > 1:
        base_name, _, base = results[0]
        print(f"  Ratios vs {base_name}:")
        for name, _, m in results[1:]:
            print(f"    {name:<28}  #wires ×{m['n_wires']/base['n_wires']:.2f}"
                  f"  area ×{m['sum_area_mm2']/base['sum_area_mm2']:.2f}"
                  f"  Σ hub ×{m['sum_hub_mm']/base['sum_hub_mm']:.2f}"
                  f"  mean ×{m['mean_hub_um']/base['mean_hub_um']:.2f}")

    # Sites analysis (per-predictor)
    sites_results = {}
    if args.sites:
        print("\nRegister-hub / per-table register distance analysis")
        print("=" * 110)
        for name, p, _ in results:
            sites_path = Path(p).with_suffix("").with_suffix(".sites.txt")
            if not sites_path.exists():
                sites_path = Path(p).parent / f"{name}.sites.txt"
            if not sites_path.exists():
                print(f"  {name}: no sites file at {sites_path}", file=sys.stderr)
                continue
            sa = analyze_sites(p, str(sites_path))
            sites_results[name] = sa
            hub_hist = sa["site_histogram"]
            total_regs = sum(hub_hist.values())
            n_sites = len(hub_hist)
            n_relevant = len(sa["relevant_rows"])
            print(f"\n  {name}")
            print(f"    Total named registers:    {total_regs}")
            print(f"    Unique site IDs occupied: {n_sites}")
            print(f"    Register-hub site:        RAM {sa['hub_site']}"
                  f"  (holds {hub_hist[sa['hub_site']]} regs)")
            if n_relevant:
                cur_mm = sa["total_current_um"] / 1e3
                cnt_mm = sa["total_counterfactual_um"] / 1e3
                print(f"    Per-table prefetch regs matched: {n_relevant}")
                print(f"    Σ current-site → owner Manhattan:   {cur_mm:>8.2f} mm")
                print(f"    Σ counterfactual (intra-table):     {cnt_mm:>8.2f} mm")
                print(f"    Wire-length recoverable by reorder: "
                      f"{(cur_mm - cnt_mm):>8.2f} mm  "
                      f"({100*(cur_mm-cnt_mm)/cur_mm if cur_mm else 0:.1f}% of current)")
            else:
                print(f"    No per-table prefetch_*[i] registers matched naming heuristic")

    if args.plot_dist:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        fig, (ax_cdf, ax_pdf) = plt.subplots(1, 2, figsize=(14, 5))
        palette = plt.cm.tab10.colors
        for i, (name, _, m) in enumerate(results):
            d = np.array(m["_hub_dists_um"])
            d_sorted = np.sort(d)
            cdf = np.arange(1, len(d_sorted)+1) / len(d_sorted)
            color = palette[i % len(palette)]
            ax_cdf.plot(d_sorted, cdf, "-", color=color, linewidth=2,
                        label=f"{name} (n={len(d)})")
            ax_pdf.hist(d, bins=40, color=color, alpha=0.45,
                        label=f"{name}", edgecolor=color, linewidth=1)

        ax_cdf.set_xlabel("RAM → register-hub Manhattan distance (μm)")
        ax_cdf.set_ylabel("CDF")
        ax_cdf.set_title("CDF of RAM → register-hub wire lengths")
        ax_cdf.grid(True, alpha=0.3)
        ax_cdf.legend(fontsize=9, loc="lower right")

        ax_pdf.set_xlabel("RAM → register-hub Manhattan distance (μm)")
        ax_pdf.set_ylabel("count")
        ax_pdf.set_title("Histogram of RAM → register-hub wire lengths")
        ax_pdf.grid(True, alpha=0.3)
        ax_pdf.legend(fontsize=9)

        fig.tight_layout()
        out_png = Path(args.render_out_dir) / "wire_dist.png"
        Path(args.render_out_dir).mkdir(parents=True, exist_ok=True)
        fig.savefig(out_png, dpi=130)
        print(f"\nDistribution plot → {out_png}")

    if args.render:
        print("\nRendering wire-bundle overlays:")
        for name, p, _ in results:
            out_png = Path(args.render_out_dir) / (Path(p).stem + f"_wires_{args.render_mode}.png")
            render_bundles(p, str(out_png), mode=args.render_mode, top_n=args.top_n,
                           hub_info=sites_results.get(name))

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
