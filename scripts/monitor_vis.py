#!/usr/bin/env python3
"""Visualize TDMonitor windowed CSV output as line graphs."""

import sys
import csv
import matplotlib.pyplot as plt

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
            rows.append(line.split(","))

    if not headers or not rows:
        print(f"No data found in {path}")
        return

    # Parse into columns
    ncols = len(headers)
    data = {h: [] for h in headers}
    for row in rows:
        for i, h in enumerate(headers):
            try:
                data[h].append(float(row[i]) if i < len(row) and row[i] else 0.0)
            except ValueError:
                data[h].append(0.0)

    wins = data.get("win", list(range(len(rows))))

    # Group columns for subplots
    groups = [
        ("Misprediction Rate", [h for h in headers if h == "misp%"]),
        ("Extra Cycle Rate", [h for h in headers if h == "extra%"]),
        ("Block Structure", [h for h in headers if h in ("i/blk", "br/blk")]),
        ("Provider Distribution", [h for h in headers if h.endswith("%") and
                                    (h.startswith("T") or h.startswith("p1%"))
                                    and h not in ("misp%", "extra%", "p1acc%", "p1p2dis%")]),
        ("P1 Accuracy & Disagree", [h for h in headers if h in ("p1acc%", "p1p2dis%")]),
        ("Allocation", [h for h in headers if h == "alloc_ok%"]),
    ]

    # Filter out empty groups
    groups = [(name, cols) for name, cols in groups if cols]

    fig, axes = plt.subplots(len(groups), 1, figsize=(12, 3 * len(groups)),
                              sharex=True)
    if len(groups) == 1:
        axes = [axes]

    for ax, (name, cols) in zip(axes, groups):
        for col in cols:
            ax.plot(wins, data[col], label=col, linewidth=0.8)
        ax.set_ylabel(name)
        ax.legend(loc="upper right", fontsize=7)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Window")
    fig.suptitle(f"TDMonitor — {path}", fontsize=11)
    fig.tight_layout()

    out = path.rsplit(".", 1)[0] + ".png"
    fig.savefig(out, dpi=150)
    print(f"Saved to {out}")
    plt.show()

if __name__ == "__main__":
    main()
