#!/usr/bin/env python3
"""Bar chart comparing per-trace MPKI of two predictor configurations.

Usage: python3 scripts/compare.py <dir1> <name1> <dir2> <name2>
  e.g. python3 scripts/compare.py full_out/hc_ir_fo_fix "RABT (ahead)" full_out/td_hc2 "TD_HC_2 (direct)"
"""

import sys
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

if len(sys.argv) != 5:
    print(__doc__)
    sys.exit(1)

dir1, name1, dir2, name2 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

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
