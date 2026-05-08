#!/usr/bin/env python3
"""VFS vs MPKI for different predictor latency configurations.

CBP-6 reference MPKI numbers (168 traces, from docs/example_and_reference_predictor_results.csv):

P1 candidates (single-cycle, fast overriding predictor):
  bimodal     ~32 KB               MPKI = 10.94
  bimodalN    ~16 KB  (N=4 banks)  MPKI = 10.48
  gshareN     ~8 KB   (N=4, GH=8) MPKI =  8.38
  gshare      ~32 KB  (GH=12)     MPKI =  6.84

P2 candidates (complex, 2-3 cycle latency):
  RUNLTS              192 KB       MPKI =  3.20   (runlts.pdf, best published CBP-6)
  TAGE-SC-L+Bullseye  187 KB       MPKI =  3.40   (bullseye.pdf)
  TAGE-SC-L (tuned)   192 KB       MPKI =  3.41   (runlts.pdf baseline)
  TASQ-SC-L (Seznec)  192 KB       MPKI =  3.41   (tage-sc-l deepdive.pdf)
  TAGE-SC-L (ref)      64 KB       MPKI =  3.50   (CBP-6 official reference)
  TAGE (no SC)         192 KB      MPKI =  3.67   (runlts.pdf)
  tage (example)       —           MPKI =  4.89   (CBP-6 example predictor)
  64KB TAGE-SC-L       64 KB       MPKI =  3.75   (tage-sc-l deepdive.pdf, CBP-5 code on CBP-6 traces)

Other context (not directly comparable, different trace sets):
  CBP-4 256Kbit TAGE-SC-L          MPKI =  2.37   (tage-sc.pdf, CBP-4 traces)
  CBP-5 64KB TAGE-SC-L             MPKI =  3.99   (CBP2016 paper, CBP-5 traces)
  Samsung Exynos M6     562 KB     MPKI =  2.54   (exynos_isca2020.pdf, internal traces)
  Samsung Exynos M1      99 KB     MPKI =  3.62   (exynos_isca2020.pdf, internal traces)

Measured FoM from reference CSV (168 CBP-6 traces, predictor_metrics.py):
  predictor   P1  P2  IPC    CPI      EPI    MPKI   VFS
  bimodal      1   1  5.875  0.1042    595   11.57  0.7792
  bimodalN     1   1  8.468  0.0995     55   11.06  0.8937
  gshare       1   2  5.863  0.0743    670    7.43  0.8429
  gshareN      1   1  7.837  0.0817     67    9.07  0.9179
  tage         1   2  5.769  0.0549   1246    5.49  0.8732
  RABT         1   1  9.009  0.0807   2837    8.96  0.8970  (full_out/hc_ir_fo_fix)

Modeling notes:
  - "N" variants (bimodalN, gshareN) are ahead-pipelined with banked SRAM, achieving
    IPC ~8-9 by predicting multiple branches per cycle. Non-N variants top out at ~5.8.
  - The reference predictor (TAGE-SC-L 64KB) reports P1=P2=0 and EPI=0 because it's
    pure C++ without HARCOM timing — its VFS=0.9729 is not comparable.
  - compute_vfs() overestimates IPC by ~5.8% vs measured (6.1 vs 5.77 for tage) because
    it uses a global average while predictor_metrics.py uses harmonic mean across traces.
    Graph 1 uses actual measured VFS to avoid this. Graph 2 still uses the model.

HARCOM latency analysis (300ps clock):
  bimodal/bimodalN: P1=P2=1 — just PC shift + SRAM read + select (~240-260ps).
  gshare/gshareN:   P1=2, P2=1 — GHR XOR before SRAM read pushes P1 to ~310ps.
    (In CBP-6 CSV these show P1=1, P2=2 because the example code structures them
     with P1=bimodal, P2=gshare rather than gshare standalone.)
  TAGE:             P1=1 (bimodal), P2=2 — tag fold + SRAM reads + tag comparison +
    sequential priority encoder across 8 tables (~150ps) + meta mux (~120ps) = ~700ps.
  TAGE-SC-L:        P1=1 (bimodal), P2=3 — SC index computation depends on TAGE's
    LongestMatchPred output, so SC can't start until TAGE finishes. SC then does 6-8
    table reads, weighted sum accumulation, and threshold check = ~1 extra cycle.
    P2=2 is optimistic (SC pipelined away or omitted, just TAGE+L).

EPI modeling for TAGE-SC-L:
  Measured tage example EPI = 1246 fJ (TAGE: bimodal + 8 tagged table reads).
  Per-read cost ≈ 1246/9 ≈ 138 fJ. SC adds ~6-8 more table reads → +830-1100 fJ.
  Estimated TAGE-SC-L EPI ≈ 2000-2400 fJ. Range 1800-3000 used for whiskers.
  Note: "N" variant EPIs (55-67 fJ) are low because banked designs only count the
  accessed bank's energy; static leakage of idle banks is not captured.

IPC modeling for TAGE-SC-L:
  We use the measured tage IPC=5.769 directly (P2=2 case) rather than compute_vfs()
  to avoid the harmonic mean error. For P2=3, we add the extra divergence cost:
    cycles_per_instr_p2_3 = 1/5.769 + diverge_ratio = 0.17334 + 0.00577 = 0.17911
    IPC_p2_3 = 1/0.17911 = 5.583
  This is valid because the delta is computed from the aggregate divergence ratio
  (0.577% of instructions have P1!=P2), which is exact. The harmonic mean issue
  only affects the base IPC, not the delta.

CPI modeling:
  CPI = MPI * (misprediction_penalty + p2_latency), where misprediction_penalty=8.
  This matches predictor_metrics.py exactly. For TAGE-SC-L at MPKI=3.50:
    P2=2: CPI = 0.00350 * 10 = 0.0350
    P2=3: CPI = 0.00350 * 11 = 0.0385

TAGE-SC-L bimodal fallback:
  TAGE-SC-L uses the same bimodal base predictor as the tage example:
    - 32K entries (2^15), 3-bit saturating counters (~12KB within the 64KB budget)
    - ~50%% of branches get no TAGE match and fall back to bimodal
    - Bimodal standalone accuracy on those branches: ~56%% (from sweep_findings.md)
  The diverge_ratio for TAGE-SC-L should be slightly higher than tage's 0.577%
  because better P2 accuracy means more P1-wrong/P2-right disagreements.
  Estimated at ~0.72%%, but we use tage's measured 0.577%% since the difference is small.
"""

import math
import numpy as np
import matplotlib.pyplot as plt

# VFS constants (from vfs.py)
IPCcbp0 = 8
CPIcbp0 = 0.0315
EPIcbp0 = 1000
ALPHA = 1.625
BETA = 4 * ALPHA / (ALPHA - 1) ** 2
GAMMA = 2 / (ALPHA - 1)
cbp_energy_ratio = 0.05
P2_TO_EXEC_STAGES = 9


def compute_vfs(mpki, p1, p2, epi, ipb=4.5, extra_ratio=0.05, diverge_ratio=0.01, divend_frac=0.8):
    mpi = mpki / 1000
    diverge_at_end = diverge_ratio * divend_frac

    if p2 <= p1:
        cycles_per_instr = max(1, p2) / ipb
    else:
        mid_diverge = diverge_ratio - diverge_at_end
        cycles_per_instr = (max(1, p1) / ipb
                            + mid_diverge * p2
                            + diverge_at_end * (p2 - max(1, p1)))
    cycles_per_instr += extra_ratio / ipb

    ipc = 1.0 / cycles_per_instr
    cpi = mpi * (P2_TO_EXEC_STAGES + p2 - p1)

    WPI0 = IPCcbp0 * CPIcbp0
    WPI = ipc * cpi
    speedup = (ipc / IPCcbp0) * (1 + WPI0) / (1 + WPI)
    LAMBDA = 1 / (1 + WPI0 / 2) - cbp_energy_ratio
    normalizedEPI = ((epi / EPIcbp0) * cbp_energy_ratio + LAMBDA * speedup ** GAMMA) * (1 + WPI / 2)
    vfs = speedup * ALPHA * (1 - 2 / (1 + math.sqrt(1 + BETA / (speedup * normalizedEPI))))
    return vfs


# --- VFS from IPC/CPI/EPI ---
def vfs_from_fom(ipc, cpi, epi):
    WPI0 = IPCcbp0 * CPIcbp0
    WPI = ipc * cpi
    speedup = (ipc / IPCcbp0) * (1 + WPI0) / (1 + WPI)
    LAMBDA = 1 / (1 + WPI0 / 2) - cbp_energy_ratio
    normalizedEPI = ((epi / EPIcbp0) * cbp_energy_ratio + LAMBDA * speedup ** GAMMA) * (1 + WPI / 2)
    return speedup * ALPHA * (1 - 2 / (1 + math.sqrt(1 + BETA / (speedup * normalizedEPI))))

# --- Actual measured predictors (from reference CSV + RABT) ---
# (label, mpki, vfs, marker, color, bold)
measured = [
    ("bimodal (P2=1cyc)",           11.57, 0.7792, "o",  "#B0BEC5", False),
    ("bimodalN (P1=P2=1)",          11.06, 0.8937, "o",  "#90CAF9", False),
    ("gshare (P2=2cyc)",             7.43, 0.8429, "s",  "#FF9800", False),
    ("gshareN (P1=P2=1)",            9.07, 0.9179, "o",  "#2196F3", False),
    ("TAGE example (P2=2cyc)",       5.49, 0.8732, "D",  "#4CAF50", False),
    # Ours
    ("RABT (P1=P2=1)",              8.96, 0.8970, "*",  "#FF0000", True),
]

# --- TAGE-SC-L variants: modeled using measured tage IPC + analytical correction ---
# See docstring for full reasoning. Summary:
#   - IPC at P2=2: use tage's measured harmonic-mean IPC=5.769 directly
#   - IPC at P2=3: add diverge_ratio (0.577%) extra cycles → IPC=5.583
#   - CPI: MPI * (8 + P2), matches predictor_metrics.py
#   - EPI: 1800–3000 fJ (tage=1246 + SC overhead of ~6-8 extra table reads)
EPI_LO, EPI_MID, EPI_HI = 1800, 2400, 3000

TAGE_IPC_P2_2 = 5.769   # measured harmonic mean from CSV (168 traces)
TAGE_DIVERGE_RATIO = 0.00577  # measured: 0.577% of instructions have P1!=P2
TAGE_IPC_P2_3 = 1.0 / (1.0 / TAGE_IPC_P2_2 + TAGE_DIVERGE_RATIO)  # = 5.583

# (label, mpki, ipc, p2, color, marker)
tage_variants = [
    ("TAGE-SC-L 64KB (P2=2cyc)", 3.50, TAGE_IPC_P2_2, 2, "#AB47BC", "^"),
    ("TAGE-SC-L 64KB (P2=3cyc)", 3.50, TAGE_IPC_P2_3, 3, "#7B1FA2", "v"),
]

from matplotlib.lines import Line2D

fig, ax = plt.subplots(figsize=(10, 7))

# Plot measured points
for label, mpki_val, vfs, marker, color, bold in measured:
    sz = 200 if bold else 100
    lw = 1.5 if bold else 0.5
    ax.scatter([mpki_val], [vfs], marker=marker, s=sz, color=color, zorder=6 if bold else 5,
               edgecolors="black", linewidths=lw, label=label)

# Plot TAGE variant whiskers (using vfs_from_fom with measured/corrected IPC)
for label, mpki_val, ipc, p2, color, marker in tage_variants:
    cpi = (mpki_val / 1000) * (8 + p2)  # matches predictor_metrics.py
    vfs_lo = vfs_from_fom(ipc, cpi, EPI_HI)
    vfs_mid = vfs_from_fom(ipc, cpi, EPI_MID)
    vfs_hi = vfs_from_fom(ipc, cpi, EPI_LO)
    # Whisker line
    ax.plot([mpki_val, mpki_val], [vfs_lo, vfs_hi], color=color, linewidth=1.5, alpha=0.6, zorder=4)
    # Caps
    cap_w = 0.08
    ax.plot([mpki_val - cap_w, mpki_val + cap_w], [vfs_lo, vfs_lo], color=color, linewidth=1.5, alpha=0.6, zorder=4)
    ax.plot([mpki_val - cap_w, mpki_val + cap_w], [vfs_hi, vfs_hi], color=color, linewidth=1.5, alpha=0.6, zorder=4)
    # Mid dot
    ax.scatter([mpki_val], [vfs_mid], marker=marker, s=100, color=color, zorder=5,
               edgecolors="black", linewidths=0.5, label=label)

ax.set_xlabel("MPKI", fontsize=18)
ax.set_ylabel("VFS", fontsize=18)
ax.tick_params(labelsize=14)
ax.set_xlim(3, 12)
ax.set_ylim(0.75, 0.95)
ax.grid(True, alpha=0.3)
ax.axhline(y=0.9, color="gray", linestyle=":", alpha=0.5)
leg = ax.legend(fontsize=12, loc="lower left")
for text in leg.get_texts():
    if "RABT" in text.get_text():
        text.set_fontweight("bold")

plt.tight_layout()
plt.savefig("latency_vs_vfs.png", dpi=150, bbox_inches="tight")
print("Saved latency_vs_vfs.png")

# ---- Graph 2: Where to invest (combined) ----
mpki_range = np.linspace(3, 10, 500)
dm = 0.01
epi_lo = 800
epi_hi = 3000
epi_mid = 1500
lat_colors = ["#2196F3", "#FF9800", "#E91E63"]
latency_configs = [(1, 1, "P2=1"), (1, 2, "P2=2"), (1, 3, "P2=3")]

fig2, ax_d = plt.subplots(figsize=(11, 6))

# Left y-axis: -dVFS/dMPKI (marginal accuracy value)
for i, (p1, p2, label) in enumerate(latency_configs):
    dvfs_lo = [-(compute_vfs(m + dm, p1, p2, epi_lo) - compute_vfs(m - dm, p1, p2, epi_lo)) / (2 * dm)
               for m in mpki_range]
    dvfs_hi = [-(compute_vfs(m + dm, p1, p2, epi_hi) - compute_vfs(m - dm, p1, p2, epi_hi)) / (2 * dm)
               for m in mpki_range]
    dvfs_mid = [-(compute_vfs(m + dm, p1, p2, epi_mid) - compute_vfs(m - dm, p1, p2, epi_mid)) / (2 * dm)
                for m in mpki_range]
    ax_d.fill_between(mpki_range, dvfs_hi, dvfs_lo, color=lat_colors[i], alpha=0.1)
    ax_d.plot(mpki_range, dvfs_mid, color=lat_colors[i], linewidth=2, label=f"-dVFS/dMPKI ({label})")

ax_d.set_xlabel("MPKI", fontsize=22, fontweight="bold")
ax_d.set_ylabel("-dVFS/dMPKI", fontsize=18, fontweight="bold", color="#333")
ax_d.tick_params(labelsize=16)
ax_d.set_xlim(3, 8)
ax_d.set_ylim(0.014, 0.026)
ax_d.grid(True, alpha=0.3)

# Right y-axis: VFS gain from latency reduction
ax_g = ax_d.twinx()

gain_2to1_lo = [compute_vfs(m, 1, 1, epi_lo) - compute_vfs(m, 1, 2, epi_lo) for m in mpki_range]
gain_2to1_hi = [compute_vfs(m, 1, 1, epi_hi) - compute_vfs(m, 1, 2, epi_hi) for m in mpki_range]
gain_2to1_mid = [compute_vfs(m, 1, 1, epi_mid) - compute_vfs(m, 1, 2, epi_mid) for m in mpki_range]
ax_g.fill_between(mpki_range, gain_2to1_hi, gain_2to1_lo, color="#2196F3", alpha=0.06)
l1, = ax_g.plot(mpki_range, gain_2to1_mid, color="#2196F3", linewidth=1.5, linestyle="--", alpha=0.8,
                label="VFS gain: P2 2->1")

gain_3to2_lo = [compute_vfs(m, 1, 2, epi_lo) - compute_vfs(m, 1, 3, epi_lo) for m in mpki_range]
gain_3to2_hi = [compute_vfs(m, 1, 2, epi_hi) - compute_vfs(m, 1, 3, epi_hi) for m in mpki_range]
gain_3to2_mid = [compute_vfs(m, 1, 2, epi_mid) - compute_vfs(m, 1, 3, epi_mid) for m in mpki_range]
ax_g.fill_between(mpki_range, gain_3to2_hi, gain_3to2_lo, color="#FF9800", alpha=0.06)
l2, = ax_g.plot(mpki_range, gain_3to2_mid, color="#FF9800", linewidth=1.5, linestyle="--", alpha=0.8,
                label="VFS gain: P2 3->2")

ax_g.set_ylabel("VFS gain (1-cycle reduction)", fontsize=18, fontweight="bold", color="gray")
ax_g.set_ylim(0.015, 0.04)
ax_g.tick_params(axis='y', colors='gray', labelsize=16)

# combined legend
lines_d, labels_d = ax_d.get_legend_handles_labels()
ax_d.legend(lines_d + [l1, l2], labels_d + [l1.get_label(), l2.get_label()],
            prop={"weight": "bold", "size": 14}, loc="upper right")


plt.tight_layout()
plt.savefig("vfs_derivative.png", dpi=150, bbox_inches="tight")
print("Saved vfs_derivative.png")
plt.show()
