# TageAhead Sweep Findings

All results on 20 representative traces with `-DFREE_FANOUT -DCHEATING_MODE` (accuracy-only, no timing/wiring constraints). VFS computed via `predictor_metrics.py` + `vfs.py`.

## Best Configs Overall (by VFS)

| # | Config | VFS | IPC | MPKI | EPI | Key Change |
|---|--------|-----|-----|------|-----|------------|
| 1 | s3_mw4_mc2048 | 0.824 | 9.31 | 14.33 | 2555 | Meta 4-bit/2048 entries |
| 2 | s5fix variants | 0.822 | 9.30 | 14.42 | 2577 | SecTagAdaptive (no effect, same as S1+decay) |
| 3 | best_8k / s3_mw4_mc1024 | 0.822 | 9.30 | 14.42 | 2557 | Meta 4-bit/1024 |
| 4 | sweep2_s1 (S1 baseline) | 0.821 | 9.30 | 14.48 | 2569 | 14T, N=7, L=256, decay TAG_OR_SEC/thresh=8 |
| 5 | s4_floor4 | 0.818 | 9.29 | 14.58 | 2569 | SecTag skip T0-T3 |
| 6 | best_16k | 0.817 | 9.32 | 14.35 | 2920 | 16K bimodal (high EPI) |
| 7 | sweep2_s1_bim2x | 0.819 | 9.33 | 14.26 | 2903 | 16K bimodal variant |

**Best MPKI** (not VFS): blk_n1_l128 at 13.94 MPKI — but VFS=0.714 due to low IPC (5.67).

---

## Sweep 1: Baseline Tuning (Decay Policies)

Base: TageAhead 14 tables, GradedSize 512-2048, 11-bit tag, MINH=8, MAXH=200, N=7, LINEINST=256.

| Config | MPKI | VFS | Decay Policy |
|--------|------|-----|-------------|
| decay_baseline | 15.45 | 0.800 | Epoch only, no decay |
| decay_v1 | 14.95 | 0.810 | TAG_OR_SEC, DECREMENT, fixed thresh=16 |
| decay_v4 | 14.96 | 0.810 | TAG, DECREMENT, graded 4-32 |
| decay_v2 | 15.04 | 0.808 | TAG_OR_SEC, DECREMENT, graded 4-32 |
| decay_v3 | 15.10 | 0.807 | TAG, DECREMENT, fixed thresh=16 |
| decay_v5 | 15.10 | 0.807 | TAG_AND_SEC, DECREMENT |

**Finding:** LFSR decay with fixed threshold gives best results. TAG_OR_SEC slightly better than TAG-only. V1/V4 both win at ~-3.3% MPKI vs baseline.

**Micro:** Decay reduces Tlast_zu% (last-table zero-use), meaning better utilization of long-history tables.

---

## Sweep 2: Decay Threshold + Allocation Tuning

Base: S1 = V1 winner (TAG_OR_SEC, thresh=8).

| Config | MPKI | VFS | Change vs S1 |
|--------|------|-----|-------------|
| Sweep2_S1 | 14.48 | 0.821 | Baseline (thresh=8) |
| Sweep2_S2 | 14.85 | 0.813 | thresh=24 |
| Sweep2_S3 | 14.95 | 0.811 | thresh=32 |
| Sweep2_S4 | 15.09 | 0.808 | Graded 8-64, TAG only |
| Sweep2_S5 | 14.98 | 0.810 | Graded 8-48, TAG only |
| Sweep2_S6 | 17.11 | 0.769 | MAX_ALLOC=2 (hurts) |
| Sweep2_S7 | 15.33 | 0.804 | V1+MAX_ALLOC=2 |
| Sweep2_S8 | 15.44 | 0.802 | V4+MAX_ALLOC=2 |
| Sweep2_S9 | 14.61 | 0.818 | LFSR graded 10-7 |
| Sweep2_S10 | 14.72 | 0.815 | LFSR graded 12-6 |
| Sweep2_S11 | 14.61 | 0.818 | Pressure-gated thresh 4-32 |
| Sweep2_S12 | 14.68 | 0.816 | Pressure-gated thresh 4-32 (gate@256) |
| Sweep2_S13 | 14.84 | 0.812 | Pressure-scaled 8-64 |
| Sweep2_S1_Bim2x | 14.26 | 0.819 | 16K bimodal (best MPKI, high EPI) |
| Sweep2_S1_Prov | 14.48 | 0.821 | Provider analysis (same as S1) |

**Finding:** Lower threshold is better (thresh=8 > 16 > 24 > 32). MAX_ALLOC=2 hurts badly. Doubled bimodal gives best raw MPKI but EPI penalty offsets VFS gain.

---

## Sweep 3: Meta Table Tuning

Base: S1 + decay.

| MW | MC=512 | MC=1024 | MC=2048 | MC=4096 |
|----|--------|---------|---------|---------|
| 1 | 14.67 / 0.817 | 14.64 / 0.818 | 14.70 / 0.816 | 14.65 / 0.818 |
| 2 | 14.45 / 0.821 | 14.48 / 0.821 | 14.45 / 0.822 | 14.43 / 0.822 |
| **4** | 14.44 / 0.821 | **14.42 / 0.822** | **14.33 / 0.824** | 14.44 / 0.820 |

**Finding:** MW=4 (4-bit meta counter) with MC=2048 is best at VFS=0.824. Wider meta counters help distinguish provider vs alt quality. MC=4096 slightly worse (diluted learning). MW=1 is clearly worse.

---

## Sweep 4: Sec-Tag Policy

Base: S1 + decay.

| Policy | MPKI | VFS | Notes |
|--------|------|-----|-------|
| SecTagNone | 16.60 | 0.776 | No sec-tag (baseline) |
| Floor4 | 14.58 | 0.818 | Skip T0-T3 |
| Floor7 | 15.20 | 0.805 | Skip T0-T6 |
| Floor10 | 15.95 | 0.789 | Skip T0-T9 |
| Ceil4 | 16.51 | 0.778 | Only T0-T3 |
| Ceil7 | 16.21 | 0.784 | Only T0-T6 |
| Ceil10 | 15.33 | 0.803 | Only T0-T9 |
| Press256 | 16.55 | 0.777 | Pressure-gated @256 |
| Press512 | 16.58 | 0.776 | Pressure-gated @512 |
| Press768 | 16.46 | 0.779 | Pressure-gated @768 |
| Acc256 | 16.61 | 0.776 | Accuracy-gated @256 |
| Acc512 | 16.64 | 0.775 | Accuracy-gated @512 |
| Acc768 | 16.49 | 0.778 | Accuracy-gated @768 |
| PGF7_512 | 16.60 | 0.776 | Combo (no better than None) |

**Finding:** Sec-tag is critical (+1.8 MPKI without it). Floor4 is best — short-history tables (T0-T3) have high collision rates so sec-tag filtering helps most there. Runtime gating policies (Press, Acc) don't work — they gate on proxy signals, not actual sec-tag benefit.

**Micro (sec-tag counterfactual):**
- With sec-tag on: ~26% of tag-matching branches rejected by sec-tag, fall back to bimodal (56% accuracy)
- Tag-only TAGE accuracy on those branches: ~51% — barely better than fallback
- Sec-tag delta: +1.3% (TAGE slightly better than FB), meaning sec-tag rejection is net positive

---

## Sweep 5: SecTagAdaptive (Benefit Counter)

Base: S1 + decay + Floor4.

All variants (8-bit/10-bit, thresholds 64-768) produced **identical results** to the base config. The benefit counter never activates because Floor4 already handles the easy cases.

**Finding:** Adaptive sec-tag is a no-op when combined with Floor4. The static policy already captures the benefit.

---

## Sweep 6: Block Size (N x LINEINST)

Base: S1 table config, varying N (branches/block) and LINEINST (instructions/line).

### MPKI

| | L=4 | L=8 | L=16 | L=32 | L=64 | L=128 | L=256 |
|--|-----|-----|------|------|------|-------|-------|
| **N=1** | 32.55 | 20.46 | 17.20 | 14.99 | 14.19 | 13.94 | — |
| **N=2** | 26.66 | 19.37 | 17.77 | 16.11 | 15.41 | 15.22 | 15.22 |
| **N=4** | 26.52 | 20.29 | 17.96 | 16.48 | 15.89 | 15.63 | — |
| **N=7** | — | — | — | — | — | — | 14.48 (S1) |

### VFS

| | L=4 | L=8 | L=16 | L=32 | L=64 | L=128 | L=256 |
|--|-----|-----|------|------|------|-------|-------|
| **N=1** | 0.510 | 0.655 | 0.695 | 0.709 | 0.713 | 0.714 | — |
| **N=2** | 0.542 | 0.714 | 0.767 | 0.774 | 0.777 | 0.777 | 0.776 |
| **N=4** | 0.537 | 0.708 | 0.768 | 0.785 | 0.789 | 0.789 | — |
| **N=7** | — | — | — | — | — | — | 0.821 (S1) |

### IPC

| | L=4 | L=8 | L=16 | L=128 | L=256 |
|--|-----|-----|------|-------|-------|
| **N=1** | 3.10 | 4.53 | 5.32 | 5.67 | — |
| **N=2** | 3.40 | 5.43 | 7.02 | 7.77 | 7.77 |
| **N=4** | 3.41 | 5.56 | 7.54 | 8.84 | — |
| **N=7** | — | — | — | — | 9.30 (S1) |

**Key Findings:**

1. **LINEINST dominates MPKI.** L=4 → L=128 cuts MPKI by ~50% regardless of N. The mechanism is sec-tag rejection: small lines cause more aliasing → higher sec-tag rejection → more fallback to bimodal.

2. **N dominates IPC (and therefore VFS).** N=7 gets 9.3 IPC vs N=2's 7.8 ceiling. The throughput advantage of wider blocks overwhelms the accuracy loss.

3. **N>1 provides negligible accuracy benefit.** At L=16, N=2 and N=4 are within 0.8 MPKI of N=1. Most blocks contain only 1 branch regardless of N.

4. **N=2 L=128 saturates at MPKI=15.2, VFS=0.777.** L=256 gives identical results.

5. **Best VFS remains N=7 L=256 (S1).** The IPC advantage of wider blocks is the dominant factor.

**Micro (N=1 L=16 vs S1):**

| Metric | S1 (N=7 L=256) | N1_L16 | N1_L128 |
|--------|----------------|--------|---------|
| alloc_ok% | 16.8% | 7.5% | 8.0% |
| true_blk% | 91.9% | 96.5% | 96.3% |
| Tlast_zu% | 37.3% | 50.0% | 50.2% |
| prov_sec_rej% | ~28% | 25.6% | 23.5% |
| prov_no_tag% | ~29% | 49.2% | 50.2% |
| cf_fb_only% | 16.9% | 20.5% | 20.6% |
| cf_tage_only% | 22.1% | 14.7% | 15.5% |
| dir_misp% | 50.3% | 54.1% | 53.3% |
| bnd_misp% | 49.6% | 45.4% | 46.0% |

The N=1 configs have much higher prov_no_tag% (50% vs 29%) — half of branches get no TAGE match at all. This is because N=1 makes more distinct block PCs competing for table entries. S1's N=7 packs more branches per entry, achieving better table utilization (Tlast_zu% 37% vs 50%).

---

## Sweep 7: History Depth for N=2

Base: N=2 BLKSWEEP config, varying MAXH (with MINH scaled proportionally).

Motivation: N=2 sees fewer branches per block → same MAXH covers less execution history. Scaling MAXH up should recover history depth.

| | MAXH=200 | MAXH=400 | MAXH=600 |
|--|----------|----------|----------|
| **L=8** | 19.37 | 22.23 (+2.9) | 23.47 (+4.1) |
| **L=16** | 17.77 | 18.54 (+0.8) | 20.48 (+2.7) |
| **L=128** | 15.22 | 16.84 (+1.6) | 18.19 (+3.0) |
| **L=256** | 15.22 | 16.81 (+1.6) | 18.18 (+3.0) |

**Finding:** Longer history **hurts uniformly** (+1.6 to +4.1 MPKI regression). The tables (512-2048 entries) are too small to support longer histories — longer folded hashes increase aliasing without adding useful correlation. H=200 is the sweet spot for these table sizes. Matching N=7's effective history depth would require proportionally larger tables, not just longer histories.

---

## Summary of Key Architectural Insights

1. **VFS is IPC-dominated.** Accuracy (MPKI) matters less than throughput for the scoring formula. N=7 wins despite slightly worse MPKI.

2. **Sec-tag is essential** but static policies (Floor4) are sufficient. Adaptive/runtime gating adds no value.

3. **Meta table quality matters.** 4-bit counters with 2048 entries give the best provider vs alt discrimination.

4. **Decay is a net win** at low thresholds (thresh=8). Higher thresholds or aggressive decay hurt.

5. **Table size is the binding constraint** for history depth. You can't just scale MAXH without scaling table sizes.

6. **Bimodal fallback accuracy (~57%) is the floor.** ~50% of branches get no TAGE match, and these dominate mispredictions. SC/Loop predictor would target this population.

7. **Remaining opportunity:** SC (statistical corrector), loop predictor, and IMLI — all target the fallback-stuck branches that TAGE can't help with.
