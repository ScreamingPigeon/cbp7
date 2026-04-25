#!/usr/bin/env python3
"""Visualize TAMonitor windowed CSV output across multiple figure pages.

Usage:
    Single trace:   python scripts/monitor_vis.py <csv_path>
    Multi-trace:    python scripts/monitor_vis.py --dir <eval_monitor_dir>

Single-trace mode produces per-page PNGs from one CSV file.
Multi-trace mode (--dir) reads all *_csv.txt + aggregate_csv.txt from an
eval-monitor output directory and produces:
    - Per-page PNGs for the aggregate (mean across traces)
    - Overlay pages showing all traces on the same axes for key metrics
"""

import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def parse_csv(path):
    rows = []
    headers = None
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
            rows.append(line.split(","))

    if not headers or not rows:
        print(f"No data found in {path}")
        return None, None, None

    data = {h: [] for h in headers}
    for row in rows:
        for i, h in enumerate(headers):
            try:
                data[h].append(float(row[i]) if i < len(row) and row[i] else 0.0)
            except ValueError:
                data[h].append(0.0)

    wins = data.get("win", list(range(len(rows))))
    return headers, data, wins


def detect_columns(headers, data):
    table_cols = [h for h in headers if
                  (h.startswith("T") and h.endswith("%") and
                   h not in ("true_blk%",) and
                   not h.startswith("Tlast"))]
    fb_col = next((h for h in headers if h in ("gs%", "bim%")), None)
    prov_cols = []
    if fb_col and fb_col in data:
        prov_cols.append(fb_col)
    prov_cols.extend(reversed(table_cols))
    return table_cols, fb_col, prov_cols


def has(data, *cols):
    return all(c in data and any(v != 0 for v in data[c]) for c in cols)


def pct(num, den):
    return 100.0 * num / den if den > 0 else 0.0


def save_fig(fig, base, suffix):
    out = f"{base}_{suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Saved {out}")
    plt.close(fig)


# ============================================================================
# Page 1: Prediction Quality
# ============================================================================
def page_quality(data, wins, prov_cols, base, title):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle(f"{title} — Prediction Quality", fontsize=13, fontweight="bold")

    # (0,0) MPKI over time
    ax = axes[0, 0]
    if "MPKI" in data:
        ax.plot(wins, data["MPKI"], color="#d62728", linewidth=1.2, label="MPKI")
        ax.set_ylabel("MPKI", color="#d62728")
        ax.tick_params(axis="y", labelcolor="#d62728")
    if "misp%" in data:
        ax2 = ax.twinx()
        ax2.plot(wins, data["misp%"], color="#1f77b4", linewidth=1.0,
                 alpha=0.7, label="misp%")
        ax2.set_ylabel("Misprediction %", color="#1f77b4")
        ax2.tick_params(axis="y", labelcolor="#1f77b4")
    ax.set_title("MPKI & Misprediction Rate")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (0,1) Provider distribution stacked area
    ax = axes[0, 1]
    if prov_cols:
        stack_data = np.array([data[c] for c in prov_cols])
        colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(prov_cols)))
        ax.stackplot(wins, stack_data, labels=prov_cols, colors=colors, alpha=0.8)
        ax.set_ylabel("Provider Share %")
        ax.legend(loc="upper right", fontsize=6, ncol=3)
        ax.set_ylim(0, 105)
    ax.set_title("Provider Distribution")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (1,0) Misprediction type breakdown (dir vs boundary)
    ax = axes[1, 0]
    plotted = False
    if has(data, "dir_misp%"):
        ax.plot(wins, data["dir_misp%"], color="#d62728", linewidth=1.2,
                label="Direction misp %")
        plotted = True
    if has(data, "bnd_misp%"):
        ax.plot(wins, data["bnd_misp%"], color="#ff7f0e", linewidth=1.2,
                label="Boundary misp %")
        plotted = True
    if plotted:
        ax.set_ylabel("% of mispredictions")
        ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Misprediction Type: Direction vs Boundary")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (1,1) MPKI with moving average overlay
    ax = axes[1, 1]
    if "MPKI" in data:
        mpki = np.array(data["MPKI"])
        ax.fill_between(wins, mpki, alpha=0.2, color="#1f77b4")
        ax.plot(wins, mpki, color="#1f77b4", linewidth=0.6, alpha=0.5, label="MPKI")
        # Moving average (window=10)
        if len(mpki) >= 10:
            kernel = np.ones(10) / 10
            ma = np.convolve(mpki, kernel, mode="valid")
            ma_wins = wins[4:4 + len(ma)]
            ax.plot(ma_wins, ma, color="#d62728", linewidth=1.5, label="MPKI (10-win MA)")
        ax.set_ylabel("MPKI")
        ax.legend(loc="upper right", fontsize=8)
    ax.set_title("MPKI Trend (with moving average)")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_fig(fig, base, "page1_quality")


# ============================================================================
# Page 2: Allocation & Table Health
# ============================================================================
def page_alloc(data, wins, base, title):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle(f"{title} — Allocation & Table Health", fontsize=13, fontweight="bold")

    # (0,0) Allocation success + pressure
    ax = axes[0, 0]
    if "alloc_ok%" in data:
        ax.plot(wins, data["alloc_ok%"], color="#2ca02c", linewidth=1.2,
                label="alloc_ok%")
    ax.set_ylabel("Allocation Success %", color="#2ca02c")
    ax.tick_params(axis="y", labelcolor="#2ca02c")
    ax2 = ax.twinx()
    plotted_r = False
    if "acc_avg" in data:
        ax2.plot(wins, data["acc_avg"], color="#ff7f0e", linewidth=0.8,
                 alpha=0.7, label="acc_avg", linestyle="--")
        plotted_r = True
    if "alloc_avg" in data:
        ax2.plot(wins, data["alloc_avg"], color="#9467bd", linewidth=0.8,
                 alpha=0.7, label="alloc_avg", linestyle="--")
        plotted_r = True
    ax2.set_ylabel("Counter Avg")
    la = ax.get_legend_handles_labels()
    lb = ax2.get_legend_handles_labels()
    ax.legend(la[0] + lb[0], la[1] + lb[1], loc="upper right", fontsize=7)
    ax.set_title("Allocation Success & Pressure")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (0,1) Entry table health — Tlast metrics
    ax = axes[0, 1]
    if "Tlast_zu%" in data:
        ax.plot(wins, data["Tlast_zu%"], color="#d62728", linewidth=1.2,
                label="Zero-Use Eviction %")
    ax.set_ylabel("Zero-Use %")
    if "Tlast_evict" in data:
        ax2 = ax.twinx()
        ax2.plot(wins, data["Tlast_evict"], color="#1f77b4", linewidth=0.8,
                 alpha=0.6, label="Evictions")
        if "Tlast_avgp" in data:
            ax2.plot(wins, data["Tlast_avgp"], color="#2ca02c", linewidth=0.8,
                     alpha=0.6, label="Avg Preds/Entry")
        ax2.set_ylabel("Count")
        la = ax.get_legend_handles_labels()
        lb = ax2.get_legend_handles_labels()
        ax.legend(la[0] + lb[0], la[1] + lb[1], loc="upper right", fontsize=7)
    ax.set_title("Entry Table (Tlast) Health")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (1,0) Collision rate + entry utilization
    ax = axes[1, 0]
    if "coll%" in data:
        ax.plot(wins, data["coll%"], color="#e377c2", linewidth=1.2,
                label="Collision %")
        ax.set_ylabel("Collision %", color="#e377c2")
    if "Tlast_used" in data:
        ax2 = ax.twinx()
        ax2.plot(wins, data["Tlast_used"], color="#17becf", linewidth=1.0,
                 label="Entries Used (tag match)")
        ax2.set_ylabel("Entries w/ Tag Match")
        la = ax.get_legend_handles_labels()
        lb = ax2.get_legend_handles_labels()
        ax.legend(la[0] + lb[0], la[1] + lb[1], loc="upper right", fontsize=7)
    ax.set_title("Collision Rate & Entry Utilization")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (1,1) Phase detection — misp rate delta
    ax = axes[1, 1]
    if has(data, "phase_delta_avg"):
        ax.plot(wins, data["phase_delta_avg"], color="#9467bd", linewidth=1.2,
                label="Avg |Δmisp rate|")
        ax.fill_between(wins, data["phase_delta_avg"], alpha=0.15, color="#9467bd")
    if has(data, "phase_max_delta"):
        ax.plot(wins, data["phase_max_delta"], color="#d62728", linewidth=0.8,
                alpha=0.7, linestyle="--", label="Max |Δmisp rate|")
    ax.set_ylabel("Misp Rate Delta (%)")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Phase Detection (256-branch window misp rate delta)")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_fig(fig, base, "page2_alloc")


# ============================================================================
# Page 3: Pipeline & Block Structure
# ============================================================================
def page_pipeline(data, wins, base, title):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle(f"{title} — Pipeline & Block Structure", fontsize=13, fontweight="bold")

    # (0,0) Extra cycles + true block %
    ax = axes[0, 0]
    if "extra%" in data:
        ax.plot(wins, data["extra%"], color="#d62728", linewidth=1.2,
                label="Extra Cycles %")
    ax.set_ylabel("Extra Cycles %", color="#d62728")
    ax.tick_params(axis="y", labelcolor="#d62728")
    if "true_blk%" in data:
        ax2 = ax.twinx()
        ax2.plot(wins, data["true_blk%"], color="#2ca02c", linewidth=1.0,
                 alpha=0.7, label="True Block %")
        ax2.set_ylabel("True Block %", color="#2ca02c")
        ax2.tick_params(axis="y", labelcolor="#2ca02c")
        la = ax.get_legend_handles_labels()
        lb = ax2.get_legend_handles_labels()
        ax.legend(la[0] + lb[0], la[1] + lb[1], loc="upper right", fontsize=7)
    ax.set_title("Pipeline Efficiency")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (0,1) Block structure — i/blk, br/blk, branch volume
    ax = axes[0, 1]
    if "i/blk" in data:
        ax.plot(wins, data["i/blk"], color="#1f77b4", linewidth=1.2,
                label="Instr/Block")
    if "br/blk" in data:
        ax.plot(wins, data["br/blk"], color="#ff7f0e", linewidth=1.0,
                alpha=0.8, label="Branches/Block")
    ax.set_ylabel("Per-Block Avg")
    if "br" in data:
        ax2 = ax.twinx()
        ax2.fill_between(wins, data["br"], alpha=0.15, color="#9467bd",
                         label="Branch Volume")
        ax2.plot(wins, data["br"], color="#9467bd", linewidth=0.6, alpha=0.4)
        ax2.set_ylabel("Branches / Window", color="#9467bd")
        ax2.tick_params(axis="y", labelcolor="#9467bd")
        la = ax.get_legend_handles_labels()
        lb = ax2.get_legend_handles_labels()
        ax.legend(la[0] + lb[0], la[1] + lb[1], loc="upper right", fontsize=7)
    ax.set_title("Block Structure & Branch Volume")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (1,0) Stuck & hard PCs per window
    ax = axes[1, 0]
    width = max(0.8, (wins[-1] - wins[0]) / len(wins) * 0.8) if len(wins) > 1 else 0.8
    if "win_stuck" in data:
        ax.bar(wins, data["win_stuck"], width=width, alpha=0.7,
               color="#d62728", label="Stuck PCs (no alloc, has misp)")
    if "win_hard" in data:
        ax.bar(wins, data["win_hard"], width=width * 0.5, alpha=0.7,
               color="#ff7f0e", label="Hard PCs (<60% bias)")
    ax.set_ylabel("PC Count")
    ax.set_title("Per-Window Problem Branches")
    ax.legend(loc="upper right", fontsize=7)
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (1,1) True block % vs extra% scatter
    ax = axes[1, 1]
    if "true_blk%" in data and "extra%" in data:
        sc = ax.scatter(data["true_blk%"], data["extra%"],
                        c=wins, cmap="viridis", s=12, alpha=0.7)
        plt.colorbar(sc, ax=ax, label="Window #")
        ax.set_xlabel("True Block %")
        ax.set_ylabel("Extra Cycles %")
    ax.set_title("True Block % vs Extra Cycles % (colored by time)")
    ax.grid(True, alpha=0.2)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_fig(fig, base, "page3_pipeline")


# ============================================================================
# Page 4: Advanced Features (micro MPKI, Jaccard, ping-pong, counterfactual)
# ============================================================================
def page_advanced(data, wins, base, title):
    fig, axes = plt.subplots(3, 2, figsize=(16, 15))
    fig.suptitle(f"{title} — Advanced Analysis", fontsize=13, fontweight="bold")

    # (0,0) Micro sliding-window MPKI p50/p95
    ax = axes[0, 0]
    if has(data, "micro_p50_mpki"):
        ax.plot(wins, data["micro_p50_mpki"], color="#2ca02c", linewidth=1.2,
                label="Micro MPKI p50")
    if has(data, "micro_p95_mpki"):
        ax.plot(wins, data["micro_p95_mpki"], color="#d62728", linewidth=1.2,
                label="Micro MPKI p95")
        ax.fill_between(wins, data["micro_p50_mpki"], data["micro_p95_mpki"],
                        alpha=0.12, color="#d62728")
    if "MPKI" in data:
        ax.plot(wins, data["MPKI"], color="#1f77b4", linewidth=0.8,
                alpha=0.5, linestyle="--", label="Window MPKI")
    ax.set_ylabel("MPKI")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Micro Sliding-Window MPKI (1024-branch window)")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (0,1) Micro MPKI p95-p50 spread (burstiness)
    ax = axes[0, 1]
    if has(data, "micro_p50_mpki", "micro_p95_mpki"):
        spread = [p95 - p50 for p50, p95 in
                  zip(data["micro_p50_mpki"], data["micro_p95_mpki"])]
        ax.fill_between(wins, spread, alpha=0.3, color="#d62728")
        ax.plot(wins, spread, color="#d62728", linewidth=1.0, label="p95 - p50")
        ax.set_ylabel("MPKI spread")
        ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Misprediction Burstiness (p95 − p50 micro MPKI)")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (1,0) Jaccard working set overlap
    ax = axes[1, 0]
    if has(data, "jaccard"):
        jac = data["jaccard"]
        # Filter out -1 (first window has no prev)
        valid_wins = [w for w, j in zip(wins, jac) if j >= 0]
        valid_jac = [j for j in jac if j >= 0]
        ax.plot(valid_wins, valid_jac, color="#1f77b4", linewidth=1.2, label="Jaccard")
        ax.fill_between(valid_wins, valid_jac, alpha=0.15, color="#1f77b4")
        ax.axhline(y=0.5, color="#999999", linestyle=":", linewidth=0.8, alpha=0.6)
        ax.set_ylabel("Jaccard Similarity")
        ax.set_ylim(-0.05, 1.05)
        ax.legend(loc="lower right", fontsize=8)
    ax.set_title("PC Working Set Overlap (Jaccard between consecutive windows)")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (1,1) Ping-pong evictions
    ax = axes[1, 1]
    if has(data, "pingpong"):
        ax.bar(wins, data["pingpong"], width=max(0.8, (wins[-1] - wins[0]) / len(wins) * 0.8) if len(wins) > 1 else 0.8,
               alpha=0.7, color="#e377c2", label="Ping-pong evictions")
        ax.set_ylabel("Ping-pong Count")
        ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Ping-Pong Evictions (same entry evicted by same PC)")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (2,0) Counterfactual: FB-only vs TAGE-only
    ax = axes[2, 0]
    if has(data, "cf_fb_only%"):
        ax.plot(wins, data["cf_fb_only%"], color="#ff7f0e", linewidth=1.2,
                label="FB-only correct %")
    if has(data, "cf_tage_only%"):
        ax.plot(wins, data["cf_tage_only%"], color="#2ca02c", linewidth=1.2,
                label="TAGE-only correct %")
    if has(data, "cf_fb_only%", "cf_tage_only%"):
        delta = [t - f for f, t in zip(data["cf_fb_only%"], data["cf_tage_only%"])]
        ax.fill_between(wins, delta, alpha=0.15,
                        color="#2ca02c", where=[d >= 0 for d in delta])
        ax.fill_between(wins, delta, alpha=0.15,
                        color="#ff7f0e", where=[d < 0 for d in delta])
    ax.set_ylabel("% of total branches")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Counterfactual: Fallback-only vs TAGE-only Correct")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    # (2,1) Counterfactual TAGE advantage (tage_only% - fb_only%)
    ax = axes[2, 1]
    if has(data, "cf_fb_only%", "cf_tage_only%"):
        advantage = [t - f for f, t in zip(data["cf_fb_only%"], data["cf_tage_only%"])]
        colors_bar = ["#2ca02c" if a >= 0 else "#d62728" for a in advantage]
        width = max(0.8, (wins[-1] - wins[0]) / len(wins) * 0.8) if len(wins) > 1 else 0.8
        ax.bar(wins, advantage, width=width, color=colors_bar, alpha=0.7)
        ax.axhline(y=0, color="black", linewidth=0.8)
        ax.set_ylabel("TAGE advantage (pp)")
        # moving average
        if len(advantage) >= 10:
            kernel = np.ones(10) / 10
            ma = np.convolve(advantage, kernel, mode="valid")
            ma_wins = wins[4:4 + len(ma)]
            ax.plot(ma_wins, ma, color="black", linewidth=1.5, linestyle="--",
                    label="10-win MA")
            ax.legend(loc="upper right", fontsize=8)
    ax.set_title("TAGE Advantage over Fallback (per window)")
    ax.set_xlabel("Window")
    ax.grid(True, alpha=0.2)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, base, "page4_advanced")


# ============================================================================
# Page 5: Histograms & Summary Bar Charts
# ============================================================================
def page_histograms(data, wins, prov_cols, base, title):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle(f"{title} — Distributions & Summaries", fontsize=13, fontweight="bold")

    # (0,0) MPKI histogram
    ax = axes[0, 0]
    if "MPKI" in data:
        mpki_vals = data["MPKI"]
        ax.hist(mpki_vals, bins=30, color="#1f77b4", alpha=0.7, edgecolor="black",
                linewidth=0.5)
        ax.axvline(np.mean(mpki_vals), color="#d62728", linestyle="--",
                   linewidth=1.5, label=f"Mean: {np.mean(mpki_vals):.3f}")
        ax.axvline(np.median(mpki_vals), color="#2ca02c", linestyle="--",
                   linewidth=1.5, label=f"Median: {np.median(mpki_vals):.3f}")
        p95 = np.percentile(mpki_vals, 95)
        ax.axvline(p95, color="#ff7f0e", linestyle=":",
                   linewidth=1.2, label=f"P95: {p95:.3f}")
        ax.set_xlabel("MPKI")
        ax.set_ylabel("Window Count")
        ax.legend(fontsize=8)
    ax.set_title("MPKI Distribution Across Windows")
    ax.grid(True, alpha=0.2)

    # (0,1) Average provider share bar chart
    ax = axes[0, 1]
    if prov_cols:
        avg_shares = [np.mean(data[c]) for c in prov_cols]
        colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(prov_cols)))
        y_pos = np.arange(len(prov_cols))
        ax.barh(y_pos, avg_shares, color=colors, alpha=0.8, edgecolor="black",
                linewidth=0.3)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(prov_cols, fontsize=8)
        ax.set_xlabel("Average Provider Share %")
        for i, v in enumerate(avg_shares):
            if v > 1.0:
                ax.text(v + 0.3, i, f"{v:.1f}%", va="center", fontsize=7)
    ax.set_title("Mean Provider Distribution")
    ax.grid(True, alpha=0.2, axis="x")

    # (1,0) Alloc success histogram
    ax = axes[1, 0]
    if "alloc_ok%" in data:
        ax.hist(data["alloc_ok%"], bins=30, color="#2ca02c", alpha=0.7,
                edgecolor="black", linewidth=0.5)
        ax.axvline(np.mean(data["alloc_ok%"]), color="#d62728", linestyle="--",
                   linewidth=1.5, label=f"Mean: {np.mean(data['alloc_ok%']):.1f}%")
        ax.set_xlabel("Allocation Success %")
        ax.set_ylabel("Window Count")
        ax.legend(fontsize=8)
    ax.set_title("Allocation Success Distribution")
    ax.grid(True, alpha=0.2)

    # (1,1) Jaccard histogram
    ax = axes[1, 1]
    if has(data, "jaccard"):
        valid_jac = [j for j in data["jaccard"] if j >= 0]
        if valid_jac:
            ax.hist(valid_jac, bins=30, color="#1f77b4", alpha=0.7,
                    edgecolor="black", linewidth=0.5)
            ax.axvline(np.mean(valid_jac), color="#d62728", linestyle="--",
                       linewidth=1.5, label=f"Mean: {np.mean(valid_jac):.3f}")
            ax.set_xlabel("Jaccard Similarity")
            ax.set_ylabel("Window Count")
            ax.legend(fontsize=8)
    ax.set_title("PC Working Set Stability Distribution")
    ax.grid(True, alpha=0.2)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_fig(fig, base, "page5_histograms")


# ============================================================================
# Page 6: Correlation Analysis
# ============================================================================
def page_correlations(data, wins, base, title):
    # Collect interesting numeric columns that have variance
    candidates = ["MPKI", "misp%", "alloc_ok%", "acc_avg", "extra%",
                  "true_blk%", "coll%", "Tlast_zu%", "dir_misp%", "bnd_misp%",
                  "micro_p50_mpki", "micro_p95_mpki", "phase_delta_avg", "phase_max_delta",
                  "jaccard", "pingpong", "cf_fb_only%", "cf_tage_only%"]
    cols = [c for c in candidates if c in data and np.std(data[c]) > 1e-9]

    if len(cols) < 4:
        return  # not enough data for meaningful correlations

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f"{title} — Correlations with MPKI", fontsize=13, fontweight="bold")

    # Scatter MPKI vs most interesting metrics
    scatter_vs = [c for c in cols if c != "MPKI"][:6]
    for idx, col in enumerate(scatter_vs):
        ax = axes[idx // 3, idx % 3]
        ax.scatter(data[col], data["MPKI"], s=10, alpha=0.5, color="#1f77b4")
        # Trend line
        if len(data[col]) > 2:
            z = np.polyfit(data[col], data["MPKI"], 1)
            p = np.poly1d(z)
            x_line = np.linspace(min(data[col]), max(data[col]), 100)
            ax.plot(x_line, p(x_line), color="#d62728", linewidth=1.2, linestyle="--")
            r = np.corrcoef(data[col], data["MPKI"])[0, 1]
            ax.set_title(f"MPKI vs {col}  (r={r:.3f})", fontsize=10)
        else:
            ax.set_title(f"MPKI vs {col}", fontsize=10)
        ax.set_xlabel(col)
        ax.set_ylabel("MPKI")
        ax.grid(True, alpha=0.2)

    # Hide unused subplots
    for idx in range(len(scatter_vs), 6):
        axes[idx // 3, idx % 3].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_fig(fig, base, "page6_correlations")


# ============================================================================
# Multi-trace overlay pages
# ============================================================================

def load_all_traces(dirpath):
    """Load all per-trace CSVs from a directory. Returns dict of name->(headers,data,wins)."""
    traces = {}
    for f in sorted(Path(dirpath).glob("*_csv.txt")):
        if "aggregate" in f.name:
            continue
        name = f.stem.replace("_csv", "")
        headers, data, wins = parse_csv(str(f))
        if data is not None:
            traces[name] = (headers, data, wins)
    return traces


def _short_name(name, maxlen=20):
    """Shorten trace name for legend."""
    if len(name) <= maxlen:
        return name
    return name[:maxlen-1] + "…"


def overlay_timeseries(traces, col, ylabel, title, base, suffix, log_y=False):
    """Plot a single column from all traces on one axis."""
    fig, ax = plt.subplots(figsize=(16, 6))
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(traces), 1)))
    for idx, (name, (headers, data, wins)) in enumerate(sorted(traces.items())):
        if col not in data:
            continue
        ax.plot(wins, data[col], linewidth=0.9, alpha=0.7, color=cmap[idx % 20],
                label=_short_name(name))
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Window")
    ax.set_title(title)
    if log_y:
        ax.set_yscale("log")
    ax.legend(loc="upper right", fontsize=6, ncol=2)
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    save_fig(fig, base, suffix)


def page_overlay_quality(traces, base):
    """Overlay page: MPKI + misp% + alloc_ok% across all traces."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle("Multi-Trace Overlay — Quality Metrics", fontsize=13, fontweight="bold")
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(traces), 1)))

    for metric, ax, ylabel in [
        ("MPKI", axes[0, 0], "MPKI"),
        ("misp%", axes[0, 1], "Misprediction %"),
        ("alloc_ok%", axes[1, 0], "Allocation Success %"),
        ("extra%", axes[1, 1], "Extra Cycles %"),
    ]:
        for idx, (name, (h, data, wins)) in enumerate(sorted(traces.items())):
            if metric not in data:
                continue
            ax.plot(wins, data[metric], linewidth=0.9, alpha=0.7,
                    color=cmap[idx % 20], label=_short_name(name))
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Window")
        ax.set_title(metric)
        ax.grid(True, alpha=0.2)

    # One shared legend at bottom
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=6,
              bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    save_fig(fig, base, "overlay_quality")


def page_overlay_advanced(traces, base):
    """Overlay page: advanced metrics across all traces."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("Multi-Trace Overlay — Advanced Metrics", fontsize=13, fontweight="bold")
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(traces), 1)))

    metrics = [
        ("micro_p50_mpki", axes[0, 0], "Micro MPKI p50"),
        ("micro_p95_mpki", axes[0, 1], "Micro MPKI p95"),
        ("phase_delta_avg", axes[0, 2], "Phase Delta Avg (%)"),
        ("jaccard", axes[1, 0], "Jaccard Similarity"),
        ("pingpong", axes[1, 1], "Ping-pong Count"),
        ("coll%", axes[1, 2], "Collision %"),
    ]

    for metric, ax, ylabel in metrics:
        for idx, (name, (h, data, wins)) in enumerate(sorted(traces.items())):
            if metric not in data:
                continue
            vals = data[metric]
            if metric == "jaccard":
                # Filter out -1 startup values
                valid_w = [w for w, v in zip(wins, vals) if v >= 0]
                valid_v = [v for v in vals if v >= 0]
                ax.plot(valid_w, valid_v, linewidth=0.9, alpha=0.7,
                        color=cmap[idx % 20], label=_short_name(name))
            else:
                ax.plot(wins, vals, linewidth=0.9, alpha=0.7,
                        color=cmap[idx % 20], label=_short_name(name))
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Window")
        ax.set_title(metric)
        ax.grid(True, alpha=0.2)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=6,
              bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    save_fig(fig, base, "overlay_advanced")


def page_trace_comparison_bars(traces, base):
    """Bar chart comparing per-trace aggregate metrics."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("Per-Trace Comparison — Aggregate Metrics", fontsize=13, fontweight="bold")

    metrics = [
        ("MPKI", axes[0, 0], "MPKI"),
        ("misp%", axes[0, 1], "Misprediction %"),
        ("alloc_ok%", axes[0, 2], "Allocation Success %"),
        ("extra%", axes[1, 0], "Extra Cycles %"),
        ("coll%", axes[1, 1], "Collision %"),
        ("Tlast_zu%", axes[1, 2], "Zero-Use Eviction %"),
    ]

    names = sorted(traces.keys())
    short_names = [_short_name(n, 15) for n in names]
    x = np.arange(len(names))

    for metric, ax, ylabel in metrics:
        vals = []
        for name in names:
            h, data, wins = traces[name]
            if metric in data and data[metric]:
                vals.append(np.mean(data[metric]))
            else:
                vals.append(0)
        colors = plt.cm.tab20(np.linspace(0, 1, len(names)))
        ax.bar(x, vals, color=colors, alpha=0.8, edgecolor="black", linewidth=0.3)
        ax.set_xticks(x)
        ax.set_xticklabels(short_names, rotation=45, ha="right", fontsize=6)
        ax.set_ylabel(ylabel)
        ax.set_title(metric)
        ax.grid(True, alpha=0.2, axis="y")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_fig(fig, base, "trace_comparison")


def visualize_single(path):
    """Generate all pages for a single CSV file."""
    headers, data, wins = parse_csv(path)
    if data is None:
        return

    base = str(Path(path).with_suffix(""))
    title = Path(path).stem
    table_cols, fb_col, prov_cols = detect_columns(headers, data)

    print(f"Parsed {len(wins)} windows, {len(headers)} columns from {path}")
    print(f"Generating visualization pages...")

    page_quality(data, wins, prov_cols, base, title)
    page_alloc(data, wins, base, title)
    page_pipeline(data, wins, base, title)
    page_advanced(data, wins, base, title)
    page_histograms(data, wins, prov_cols, base, title)
    page_correlations(data, wins, base, title)


def visualize_dir(dirpath):
    """Generate per-trace pages, aggregate pages, and overlay pages."""
    dirpath = Path(dirpath)

    # Load all per-trace CSVs
    traces = load_all_traces(dirpath)
    if not traces:
        print(f"No per-trace CSVs found in {dirpath}")
        return

    print(f"Loaded {len(traces)} traces from {dirpath}")

    # Per-trace: generate full page set for each trace
    for name, (headers, data, wins) in sorted(traces.items()):
        trace_base = str(dirpath / name)
        title = name
        _, _, prov_cols = detect_columns(headers, data)
        print(f"  Generating pages for {name} ({len(wins)} windows)...")
        page_quality(data, wins, prov_cols, trace_base, title)
        page_alloc(data, wins, trace_base, title)
        page_pipeline(data, wins, trace_base, title)
        page_advanced(data, wins, trace_base, title)
        page_histograms(data, wins, prov_cols, trace_base, title)
        page_correlations(data, wins, trace_base, title)

    # Aggregate: if aggregate_csv.txt exists, generate full page set
    agg_path = dirpath / "aggregate_csv.txt"
    if agg_path.exists():
        print(f"  Generating pages for aggregate...")
        visualize_single(str(agg_path))

    # Overlay pages: all traces on same axes
    overlay_base = str(dirpath / "all_traces")
    print(f"  Generating overlay pages...")
    page_overlay_quality(traces, overlay_base)
    page_overlay_advanced(traces, overlay_base)
    page_trace_comparison_bars(traces, overlay_base)

    print("Done.")


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "--dir":
        visualize_dir(sys.argv[2])
    elif len(sys.argv) >= 2:
        visualize_single(sys.argv[1])
    else:
        visualize_single("out/monitor_csv.txt")


if __name__ == "__main__":
    main()
