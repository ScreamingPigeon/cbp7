#!/usr/bin/env python3
"""Visualize TAMonitor windowed CSV output as line graphs and histograms."""

import sys
import csv
import matplotlib.pyplot as plt
import numpy as np

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "out/monitor_csv.txt"

    # Parse CSV, skip comment lines
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
            # Stop at non-CSV summary section
            if line.startswith("=") or line.startswith("Instructions"):
                break
            rows.append(line.split(","))

    if not headers or not rows:
        print(f"No data found in {path}")
        return

    # Parse into columns
    data = {h: [] for h in headers}
    for row in rows:
        for i, h in enumerate(headers):
            try:
                data[h].append(float(row[i]) if i < len(row) and row[i] else 0.0)
            except ValueError:
                data[h].append(0.0)

    wins = data.get("win", list(range(len(rows))))
    nwin = len(wins)

    # Detect table columns (T0%, T1%, ... or gs%/bim%)
    table_cols = [h for h in headers if
                  (h.startswith("T") and h.endswith("%") and
                   h not in ("true_blk%",) and
                   not h.startswith("Tlast"))]
    fb_col = next((h for h in headers if h in ("gs%", "bim%")), None)

    # ---- Build figure layout ----
    # Row 0: MPKI + misp%                | Provider distribution stacked area
    # Row 1: Allocation + pressure       | Entry table health (lifetime)
    # Row 2: Pipeline (extra%, true_blk) | Block structure (i/blk, br/blk, br)
    # Row 3: Collision + utilization     | Per-window stuck/hard PCs
    # Row 4: MPKI histogram              | Average provider share bar

    fig = plt.figure(figsize=(16, 22))
    gs = fig.add_gridspec(5, 2, hspace=0.38, wspace=0.3)

    # ================================================================
    # (0,0) MPKI & Misprediction Rate
    # ================================================================
    ax1 = fig.add_subplot(gs[0, 0])
    if "MPKI" in data:
        ax1.plot(wins, data["MPKI"], color="#d62728", linewidth=1.2, label="MPKI")
        ax1.set_ylabel("MPKI", color="#d62728")
        ax1.tick_params(axis="y", labelcolor="#d62728")
    if "misp%" in data:
        ax1b = ax1.twinx()
        ax1b.plot(wins, data["misp%"], color="#1f77b4", linewidth=1.0,
                  alpha=0.7, label="misp%")
        ax1b.set_ylabel("Misprediction %", color="#1f77b4")
        ax1b.tick_params(axis="y", labelcolor="#1f77b4")
    ax1.set_title("Prediction Quality")
    ax1.grid(True, alpha=0.2)

    # ================================================================
    # (0,1) Provider Distribution — stacked area
    # ================================================================
    ax2 = fig.add_subplot(gs[0, 1])
    prov_cols = []
    if fb_col and fb_col in data:
        prov_cols.append(fb_col)
    # Add table cols in reverse order (longest history first)
    prov_cols.extend(reversed(table_cols))
    if prov_cols:
        stack_data = np.array([data[c] for c in prov_cols])
        colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(prov_cols)))
        ax2.stackplot(wins, stack_data, labels=prov_cols, colors=colors, alpha=0.8)
        ax2.set_ylabel("Provider Share %")
        ax2.set_title("Provider Distribution")
        ax2.legend(loc="upper right", fontsize=6, ncol=3)
        ax2.set_ylim(0, 105)
    ax2.grid(True, alpha=0.2)

    # ================================================================
    # (1,0) Allocation Success + Pressure Counters
    # ================================================================
    ax3 = fig.add_subplot(gs[1, 0])
    if "alloc_ok%" in data:
        ax3.plot(wins, data["alloc_ok%"], color="#2ca02c", linewidth=1.2,
                 label="alloc_ok%")
    ax3.set_ylabel("Allocation Success %", color="#2ca02c")
    ax3.tick_params(axis="y", labelcolor="#2ca02c")
    ax3.set_title("Allocation & Pressure")
    ax3b = ax3.twinx()
    if "acc_avg" in data:
        ax3b.plot(wins, data["acc_avg"], color="#ff7f0e", linewidth=0.8,
                  alpha=0.7, label="acc_avg", linestyle="--")
    if "alloc_avg" in data:
        ax3b.plot(wins, data["alloc_avg"], color="#9467bd", linewidth=0.8,
                  alpha=0.7, label="alloc_avg", linestyle="--")
    ax3b.set_ylabel("Counter Avg")
    lines3a = ax3.get_legend_handles_labels()
    lines3b = ax3b.get_legend_handles_labels()
    ax3.legend(lines3a[0] + lines3b[0], lines3a[1] + lines3b[1],
               loc="upper right", fontsize=7)
    ax3.grid(True, alpha=0.2)

    # ================================================================
    # (1,1) Entry Table Health — lifetime metrics
    # ================================================================
    ax4 = fig.add_subplot(gs[1, 1])
    if "Tlast_zu%" in data:
        ax4.plot(wins, data["Tlast_zu%"], color="#d62728", linewidth=1.2,
                 label="Zero-Use Eviction %")
    if "Tlast_evict" in data:
        ax4_b = ax4.twinx()
        ax4_b.plot(wins, data["Tlast_evict"], color="#1f77b4", linewidth=0.8,
                   alpha=0.6, label="Evictions")
        if "Tlast_avgp" in data:
            ax4_b.plot(wins, data["Tlast_avgp"], color="#2ca02c", linewidth=0.8,
                       alpha=0.6, label="Avg Preds/Entry")
        ax4_b.set_ylabel("Count")
        lines4a = ax4.get_legend_handles_labels()
        lines4b = ax4_b.get_legend_handles_labels()
        ax4.legend(lines4a[0] + lines4b[0], lines4a[1] + lines4b[1],
                   loc="upper right", fontsize=7)
    ax4.set_ylabel("Zero-Use %")
    ax4.set_title("Entry Table (Tlast) Health")
    ax4.grid(True, alpha=0.2)

    # ================================================================
    # (2,0) Pipeline Efficiency — extra cycles + true block %
    # ================================================================
    ax5 = fig.add_subplot(gs[2, 0])
    has_pipe = False
    if "extra%" in data:
        ax5.plot(wins, data["extra%"], color="#d62728", linewidth=1.2,
                 label="Extra Cycles %")
        has_pipe = True
    ax5.set_ylabel("Extra Cycles %", color="#d62728")
    ax5.tick_params(axis="y", labelcolor="#d62728")
    if "true_blk%" in data:
        ax5b = ax5.twinx()
        ax5b.plot(wins, data["true_blk%"], color="#2ca02c", linewidth=1.0,
                  alpha=0.7, label="True Block %")
        ax5b.set_ylabel("True Block %", color="#2ca02c")
        ax5b.tick_params(axis="y", labelcolor="#2ca02c")
        has_pipe = True
        lines5a = ax5.get_legend_handles_labels()
        lines5b = ax5b.get_legend_handles_labels()
        ax5.legend(lines5a[0] + lines5b[0], lines5a[1] + lines5b[1],
                   loc="upper right", fontsize=7)
    elif has_pipe:
        ax5.legend(loc="upper right", fontsize=7)
    ax5.set_title("Pipeline Efficiency")
    ax5.grid(True, alpha=0.2)

    # ================================================================
    # (2,1) Block Structure — i/blk, br/blk, branch volume
    # ================================================================
    ax6 = fig.add_subplot(gs[2, 1])
    has_blk = False
    if "i/blk" in data:
        ax6.plot(wins, data["i/blk"], color="#1f77b4", linewidth=1.2,
                 label="Instr/Block")
        has_blk = True
    if "br/blk" in data:
        ax6.plot(wins, data["br/blk"], color="#ff7f0e", linewidth=1.0,
                 alpha=0.8, label="Branches/Block")
        has_blk = True
    ax6.set_ylabel("Per-Block Avg")
    if "br" in data:
        ax6b = ax6.twinx()
        ax6b.fill_between(wins, data["br"], alpha=0.15, color="#9467bd",
                          label="Branch Volume")
        ax6b.plot(wins, data["br"], color="#9467bd", linewidth=0.6, alpha=0.4)
        ax6b.set_ylabel("Branches / Window", color="#9467bd")
        ax6b.tick_params(axis="y", labelcolor="#9467bd")
        has_blk = True
        lines6a = ax6.get_legend_handles_labels()
        lines6b = ax6b.get_legend_handles_labels()
        ax6.legend(lines6a[0] + lines6b[0], lines6a[1] + lines6b[1],
                   loc="upper right", fontsize=7)
    elif has_blk:
        ax6.legend(loc="upper right", fontsize=7)
    ax6.set_title("Block Structure & Branch Volume")
    ax6.grid(True, alpha=0.2)

    # ================================================================
    # (3,0) Collision Rate + Entry Utilization
    # ================================================================
    ax7 = fig.add_subplot(gs[3, 0])
    if "coll%" in data:
        ax7.plot(wins, data["coll%"], color="#e377c2", linewidth=1.2,
                 label="Collision %")
        ax7.set_ylabel("Collision %", color="#e377c2")
    if "Tlast_used" in data:
        ax7b = ax7.twinx()
        ax7b.plot(wins, data["Tlast_used"], color="#17becf", linewidth=1.0,
                  label="Entries Used")
        ax7b.set_ylabel("Entries w/ Tag Match")
        lines7a = ax7.get_legend_handles_labels()
        lines7b = ax7b.get_legend_handles_labels()
        ax7.legend(lines7a[0] + lines7b[0], lines7a[1] + lines7b[1],
                   loc="upper right", fontsize=7)
    ax7.set_title("Collision Rate & Entry Utilization")
    ax7.grid(True, alpha=0.2)

    # ================================================================
    # (3,1) Per-Window Stuck & Hard PCs
    # ================================================================
    ax8 = fig.add_subplot(gs[3, 1])
    if "win_stuck" in data:
        ax8.bar(wins, data["win_stuck"], width=0.8, alpha=0.7,
                color="#d62728", label="Stuck PCs (no alloc, has misp)")
    if "win_hard" in data:
        ax8.bar(wins, data["win_hard"], width=0.4, alpha=0.7,
                color="#ff7f0e", label="Hard PCs (<60% bias)")
    ax8.set_ylabel("PC Count")
    ax8.set_title("Per-Window Problem Branches")
    ax8.legend(loc="upper right", fontsize=7)
    ax8.grid(True, alpha=0.2)

    # ================================================================
    # (4,0) MPKI Histogram
    # ================================================================
    ax9 = fig.add_subplot(gs[4, 0])
    if "MPKI" in data:
        mpki_vals = data["MPKI"]
        ax9.hist(mpki_vals, bins=30, color="#1f77b4", alpha=0.7, edgecolor="black",
                 linewidth=0.5)
        ax9.axvline(np.mean(mpki_vals), color="#d62728", linestyle="--",
                    linewidth=1.5, label=f"Mean: {np.mean(mpki_vals):.3f}")
        ax9.axvline(np.median(mpki_vals), color="#2ca02c", linestyle="--",
                    linewidth=1.5, label=f"Median: {np.median(mpki_vals):.3f}")
        p95 = np.percentile(mpki_vals, 95)
        ax9.axvline(p95, color="#ff7f0e", linestyle=":",
                    linewidth=1.2, label=f"P95: {p95:.3f}")
        ax9.set_xlabel("MPKI")
        ax9.set_ylabel("Window Count")
        ax9.set_title("MPKI Distribution Across Windows")
        ax9.legend(fontsize=8)
    ax9.grid(True, alpha=0.2)

    # ================================================================
    # (4,1) Average Provider Share — horizontal bar
    # ================================================================
    ax10 = fig.add_subplot(gs[4, 1])
    if prov_cols:
        avg_shares = [np.mean(data[c]) for c in prov_cols]
        colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(prov_cols)))
        y_pos = np.arange(len(prov_cols))
        ax10.barh(y_pos, avg_shares, color=colors, alpha=0.8, edgecolor="black",
                 linewidth=0.3)
        ax10.set_yticks(y_pos)
        ax10.set_yticklabels(prov_cols, fontsize=8)
        ax10.set_xlabel("Average Provider Share %")
        ax10.set_title("Mean Provider Distribution")
        for i, v in enumerate(avg_shares):
            if v > 1.0:
                ax10.text(v + 0.3, i, f"{v:.1f}%", va="center", fontsize=7)
    ax10.grid(True, alpha=0.2, axis="x")

    # ---- Common x-label ----
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]:
        ax.set_xlabel("Window", fontsize=8)

    fig.suptitle(f"TAMonitor Dashboard — {path}", fontsize=13, fontweight="bold")

    out = path.rsplit(".", 1)[0] + ".png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved to {out}")
    plt.show()

if __name__ == "__main__":
    main()
