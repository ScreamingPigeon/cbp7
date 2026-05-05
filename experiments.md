# Experiment Log

All experiments run on 20 representative traces via `make eval-monitor`.
Baseline: S3_MW4_MC2048 (14 tables, UniformTag<11>, MINHIST=8, MAXHIST=200,
FB=8K, META MW=4/MC=2048, decay FixedThresh<8>, SecTagAll).

## Sweep 6a: Decay Tuning (LFSR widths + thresholds)

| Config | Description | MPKI | VFS | vs baseline |
|--------|-------------|------|-----|-------------|
| S3_MW4_MC2048 | baseline: uniform LFSR=8, FixedThresh<8> | 13.069 | 0.8238 | -- |
| S6_G10_7 | graded LFSR 10->7, FixedThresh<8> | 13.093 | 0.8232 | -0.07% |
| S6_G12_6 | graded LFSR 12->6, FixedThresh<8> | 13.247 | 0.8198 | -0.48% |
| S6_G10_7_T16 | graded LFSR 10->7, FixedThresh<16> | 13.214 | 0.8204 | -0.41% |
| S6_G12_6_T16 | graded LFSR 12->6, FixedThresh<16> | 13.502 | 0.8142 | -1.16% |
| S6_GT8_64 | uniform LFSR=8, GradedThresh 8->64 | 13.741 | 0.8087 | -1.83% |
| S6_PG512 | uniform LFSR=8, PressGated 4->32@512 | 13.334 | 0.8177 | -0.74% |

**Verdict**: More aggressive/graded decay consistently hurts. Baseline is optimal.

## Sweep 6b: Allocation, Bimodal, Meta Sizing

| Config | Description | MPKI | VFS | vs baseline |
|--------|-------------|------|-----|-------------|
| S3_MW4_MC2048 | baseline | 13.069 | 0.8238 | -- |
| S6_Alloc2 | MAX_ALLOC=2 | 13.570 | 0.8122 | -1.40% |
| S6_FB16K | bimodal 16K | 13.032 | 0.8219 | -0.23% |
| S6_MC4096 | meta 4096 entries | 13.183 | 0.8196 | -0.51% |
| S6_All3 | alloc2 + FB16K + MC4096 | 13.559 | 0.8125 | -1.37% |

**Verdict**: MAX_ALLOC=2 pollutes. FB 16K slight MPKI win but VFS drops. MC4096 worse. Baseline wins.

## Sweep 7: Structural (adaptive sec-tag, tag widths, graded tags, history lengths)

| Config | Description | MPKI | VFS | vs baseline |
|--------|-------------|------|-----|-------------|
| S3_MW4_MC2048 | baseline | | | -- |
| S7_Adapt64 | SecTagAdaptive<8,64> | | | |
| S7_Adapt96 | SecTagAdaptive<8,96> | | | |
| S7_Adapt128 | SecTagAdaptive<8,128> | | | |
| S7_Tag12 | UniformTag<12> | | | |
| S7_GT13_9 | GradedTag<13,9> (T0:13, T13:9) | | | |
| S7_GT14_8 | GradedTag<14,8> (T0:14, T13:8) | | | |
| S7_GT12_10 | GradedTag<12,10> (T0:12, T13:10) | | | |
| S7_H6_300 | MINHIST=6, MAXHIST=300 | | | |
| S7_H5_150 | MINHIST=5, MAXHIST=150 | | | |
| S7_H10_250 | MINHIST=10, MAXHIST=250 | | | |

## BR_P_ENTRY: Per-Group Tag Encoding (2026-05-04)

Encode group_id in tag so branches in different groups get independent resolution
chains and providers. `NUM_GROUPS = ceil(N / BR_P_ENTRY)`. Contiguous grouping
(INTERLEAVED=false). Same table/size/tag config as baseline 1C.

### Timing (gcc trace, all configs meet 1-cycle)

| Config | BR_P_ENTRY | NUM_GROUPS | P1 Latency (cycles) |
|--------|-----------|------------|---------------------|
| 1C     | 7 (=N)    | 1          | 0.677               |
| BPE4   | 4         | 2          | 0.863               |
| BPE2   | 2         | 4          | 0.820               |
| BPE1   | 1         | 7          | 0.737               |

### Accuracy (20 traces, eval-monitor 40M instr)

| Trace | 1C | BPE4 | BPE2 | BPE1 |
|-------|------|------|------|------|
| 502-gcc | 652,742 | 652,814 | 620,352 | **584,626** |
| 505-mcf | 685,446 | **661,175** | 688,652 | 668,608 |
| 508-namd | 198,046 | 191,132 | 198,301 | **184,906** |
| 531-deepsjeng | 876,443 | 874,463 | 858,298 | **817,062** |
| 548-exchange2 | **527,595** | 492,623 | **469,314** | 611,536 |
| 554-roms | 15,049 | **12,719** | 12,784 | 13,680 |
| dcapo-kafka | **246,089** | 250,341 | 316,093 | 403,362 |
| gap-sssp | 1,057,438 | 1,078,686 | 1,083,152 | **1,033,592** |
| gcc-1 | 279,416 | 277,929 | 270,682 | **248,516** |
| java16 | 380,819 | 376,433 | 360,273 | **329,676** |
| llvm-2 | 254,813 | 249,482 | **244,629** | 247,276 |
| lua-3 | 35,677 | 42,963 | **25,577** | 35,881 |
| nodejs-http2 | 358,138 | 353,560 | 349,533 | **335,027** |
| nodejs-octane | 352,280 | 342,387 | 340,671 | **304,558** |
| python3-dulwich | 275,048 | 269,332 | 260,205 | **238,509** |
| rsbench | 613,990 | 609,539 | 598,599 | **586,369** |
| sampleflow | **99,889** | **74,557** | 100,213 | 135,795 |
| web_130 | 622,304 | 626,916 | 620,588 | **598,124** |
| web_74 | 457,590 | 464,210 | 449,798 | **419,462** |
| zstd | **180,312** | 232,218 | 202,694 | 209,835 |
| **Total** | **8,169,124** | **8,133,479 (-0.4%)** | **8,070,408 (-1.2%)** | **8,006,400 (-2.0%)** |

BPE1 wins 12/20 traces (-2.0% aggregate). BPE1 regressions:
- dcapo-kafka +63.9%, sampleflow +35.9%, exchange2 +15.9%, zstd +16.4%
- These may be workloads where intra-block branch correlation benefits from
  shared providers.

### Bank Balance (mcf trace, write imbalance ratio B0:B1)

| Table | 1C | BPE4 | BPE1 |
|-------|------|------|------|
| T0 | 78:1 | 36:1 | 1.2:1 |
| T1 | 11:1 | 97:1 | 1.1:1 |
| T2 | 1.1:1 | 1.5:1 | 1.0:1 |
| T12 | 5.1:1 | 6.2:1 | 4.3:1 |
| T13 | 4.0:1 | 5.9:1 | 2.9:1 |

BPE1 fixes low-table imbalance completely. High-table imbalance persists
(index hash bias, not group-related).

### rwram conflict stats (1C, mcf)

100% of writes buffered (read+write same cycle). 53% of buffered writes lost
(overwritten before flush). Applies to all configs — structural issue with
1-write-per-cycle ta_rwram.

## TageAheadHC_IR: Hand-coded BPE1 (2026-05-04)

Ported BPE1 independent resolution chains from template TageAhead to
hand-coded TageAheadHC. Each entry stores 1-bit pred (was 7-bit). 11-bit
tag split into 8-bit htag + 3-bit group_id. 7 parallel resolve_chains.

### Single-trace results (gcc, 100K warmup + 100K measure)

| Config | Mispred | P2 (ns) | EPI (fJ) |
|--------|---------|---------|----------|
| TageAheadHC (baseline, no FREE_FANOUT) | 670 | 0.89 | 5305 |
| TageAheadHC_IR v1 (no FREE_FANOUT) | 859 | 0.997 | 5298 |
| TageAheadHC_IR v2 optimized (no FREE_FANOUT) | 859 | 0.997 | 5265 |
| TageAheadHC_IR v1 (FREE_FANOUT) | 859 | 0.61 | 5241 |
| TageAheadHC_IR v2 optimized (FREE_FANOUT) | 859 | 0.61 | 5211 |

### v2 optimizations applied

1. **Precompute sec_match per table**: sec_tag_now fanout 99→15 (7 FO2 stages → 4).
   Removed 84 redundant 5-bit comparators (was NT×NUM_GROUPS=98, now NT=14).
2. **Precompute htag_sec per table**: prefetch_tag_hit fanout 8→2, prefetch_sec fanout 8→2.
   Combined htag_hit & sec_match once, then fan out to 7 group chains.
3. **Precompute fb_bits via make_array**: prefetch_fb fanout 8→2.
   Split N-wide fb into individual bits once, each group reads its bit via fo1().
4. **Simplify PRED_BITS=1 resolve_chain**: Removed redundant replicate(hard<1>).concat()
   (identity for 1-bit). Removed `!= hard<0>{}` on 1-bit XOR (identity).
5. **RAM reorder**: tag, fold_idx, fold_tag, sec, pred, zone, hyst, u
   (was: tag, fold, pred, sec). Tiny 1-bit pred_rams placed after sec for tighter packing.
6. **Fixed fo1/fanout conflicts**: `merged` variable (fo1+fanout on same var),
   outer `pp` from resolve_chain return (fo1 but read twice → fanout(hard<2>)).

### Floorplan layout experiments

| Layout | P2 (ns) | EPI (fJ) | Notes |
|--------|---------|----------|-------|
| tag, fold, **pred**, sec, zone, hyst, u (v1) | 0.997 | 5298 | Tiny pred breaks predict cluster |
| tag, fold, sec, **pred**, zone, hyst, u (v2) | 0.997 | 5265 | Pred still in predict cluster |
| tag, fold, sec, zone, **pred**, hyst, u (v3) | **0.867** | **4910** | Pred in update zone ★ |

**Key insight**: With 1-bit pred_rams (7× smaller than HC's 7-bit), moving
them out of the predict cluster lets tag/fold/sec pack tightly. The 1-bit
pred value travels cheaply from the update zone. -13% latency, -7.3% EPI.

### Notes

- P2 latency 0.867ns vs HC's 0.89ns: IR is now *faster* than HC baseline
  thanks to smaller pred_rams and tighter predict cluster packing.
- Logic floor (FREE_FANOUT) is 0.61ns for both HC and HC_IR.
- MPKI is higher on gcc 100K (859 vs 670) but this is expected — IR trades
  intra-block correlation for per-branch independence. Full 20-trace eval
  shows -2% aggregate MPKI improvement for BPE1.

## Sweep 9: BPE=1 + Tag Width Sweep (2026-05-05)

Motivation: with BPE=1, `NUM_GROUPS=7`, `GROUP_BITS=3` — 3 of the 11 tag bits
encode group_id, leaving only 8 effective disambiguation bits. Explored shorter
raw tag widths and graded schemes to match or improve on 11-bit uniform.

Base config: S1 (TA1C_BASE + decay FixedThresh<8> TAG_OR_SEC), BR_P_ENTRY=1.
All configs use GradedSize<512,2048>. Effective bits = raw − GROUP_BITS(3).

| Config | Tag scheme | Raw range | Eff range | MPKI |
|--------|-----------|-----------|-----------|------|
| BPE1_GT11_7 | GradedTag<11,7> | 11→7 | 8→4 | **12.357** ← best |
| BPE1_U8 | UniformTag<8> | 8 | 5 | 12.518 |
| BPE1_GT12_8 | GradedTag<12,8> | 12→8 | 9→5 | 12.518 |
| BPE1_GT14_8 | GradedTag<14,8> | 14→8 | 11→5 | 12.448 |
| BPE1_U9 | UniformTag<9> | 9 | 6 | 12.556 |
| BPE1_U11 | UniformTag<11> | 11 | 8 | 12.556 |
| BPE1_U10 | UniformTag<10> | 10 | 7 | 12.758 |
| BPE1_GT10_6 | GradedTag<10,6> | 10→6 | 7→3 | 12.792 |
| BPE1_GT13_9 | GradedTag<13,9> | 13→9 | 10→6 | 12.792 |

**Verdict**: GradedTag<11,7> wins (+0.2 MPKI over uniform U11). Longer tags on
small/short-history tables, shorter on large/long-history tables is the right
tradeoff. Surprisingly, U8 (eff 5) ties GT12_8 (eff 9→5) — effective bits
matter more than raw bits. U10 and wider graded configs hurt, likely over-
tagging the small tables wastes index bits without helping disambiguation.

## HC_IR NT=15/16: Extra Tables Experiment (2026-05-05)

Motivation: v3 floorplan had white space between small (T0-4, 1024-entry)
and large (T5-13, 2048-entry) table clusters. Attempted to fill gap with
1-2 extra 2048-entry tables while keeping geometric_hist spacing.

### NT=16 (2 extra tables)

| Placement | P2 (ns) | EPI (fJ) | Misps (gcc 100K) |
|-----------|---------|----------|------------------|
| T14/T15 at end (after T13) | 1.103 | 6355 | 914 |
| T14 between T4-T5, T15 between T7-fb | 1.103 | 6355 | 914 |

Floorplan grew an extra row vertically. Reordering RAM declarations didn't
help — total area is the same regardless of placement order.

### NT=15 (1 extra table)

T14 placed between T4 and T5 in declaration order.

| Config | P2 (ns) | EPI (fJ) | Misps (gcc 100K) |
|--------|---------|----------|------------------|
| NT=14 (v3 baseline) | 0.867 | 4910 | 859 |
| NT=15 | 0.867 | 5328 | 891 |

Timing held at 0.867ns. But accuracy got *worse* (859→891 misps) and EPI
increased +8.5%. The geometric_hist<15>(8,200) reshuffles all history lengths
vs geometric_hist<14>(8,200), losing tuned sweet spots.

### Quick-eval (20/20 traces, 1M warmup + 40M measure)

| Config | VFS | MPKI | EPI (fJ) | IPC_cbp | P2 (ns) |
|--------|-----|------|----------|---------|---------|
| NT=14 | **0.878** | 10.81 | 2192 | 9.537 | 0.867 |
| NT=15 | 0.878 | 10.71 | 2374 | 9.549 | 0.867 |

NT=15 trades -0.1 MPKI for +182 fJ EPI. VFS essentially identical (0.8783 vs 0.8776).

**Verdict**: Extra tables don't pay off. NT=16 blows timing (+27%). NT=15
holds timing and marginally improves MPKI but the extra area/energy cancels
it out in VFS. The geometric history redistribution also reshuffles all
history lengths which may lose tuned sweet spots.

## Full Eval: TageAheadHC vs TageAheadHC_IR (2026-05-05)

All 168 traces, 1M warmup + 40M measure. HC = baseline (NT=14, 7-bit pred,
shared resolution). HC_IR = NT=15, BPE1 independent resolution (1-bit pred,
7 parallel resolve chains, 8-bit htag + 3-bit group_id).

### Aggregate

| Config | VFS | MPKI | EPI (fJ) | IPC_cbp | Max P2 (ns) |
|--------|-----|------|----------|---------|-------------|
| TageAheadHC | 0.8815 | 9.146 | 2571 | 8.903 | 0.890 |
| TageAheadHC_IR | **0.8848** | **9.089** | **2386** | 8.874 | **0.867** |
| Delta | **+0.37%** | **-0.63%** | **-7.2%** | -0.33% | -2.6% |

IR wins 79/168 traces, loses 88/168 (1 tie). Net improvement driven by
-7.2% EPI and -2.6% timing rather than raw accuracy (nearly flat MPKI).

### Top IR wins (misprediction reduction)

| Trace | HC misps | IR misps | Delta |
|-------|----------|----------|-------|
| roms-1.109375_0 | 27,023 | 6,348 | -76.5% |
| sampleflow-1.127768_0 | 50,506 | 24,885 | -50.7% |
| int_87 | 149,328 | 82,918 | -44.5% |
| int_117 | 1,235,001 | 814,430 | -34.1% |
| lua-3.19933_0 | 26,602 | 20,675 | -22.3% |
| fp_28 | 75,127 | 59,943 | -20.2% |
| dcapo-h2o-ai-jdk11 | 278,712 | 247,854 | -11.1% |
| 548-exchange2 | 164,742 | 149,857 | -9.0% |

### Top IR losses (misprediction increase)

| Trace | HC misps | IR misps | Delta |
|-------|----------|----------|-------|
| ren-reactors-jdk17.15692_0 | 3,299 | 17,935 | +443.6% |
| 554-roms-1_62613 | 2,069 | 3,581 | +73.1% |
| infra_75 | 14,080 | 23,496 | +66.9% |
| zstd-1.19139_0 | 82,573 | 114,480 | +38.6% |
| fp_11 | 104,952 | 138,719 | +32.2% |
| compress_7 | 173,592 | 221,259 | +27.5% |
| brotli-1.6219_0 | 126,999 | 159,105 | +25.3% |
| gap-bc-coauthors | 228,735 | 281,664 | +23.1% |

### Monitor stats (HC_IR, 20-trace eval-monitor)

| Metric | Mean | Median | P95 |
|--------|------|--------|-----|
| MPKI | 15.50 | 13.99 | 41.60 |
| alloc_ok% | 17.4% | 14.6% | 51.8% |
| Tlast_zu% | 58.0% | 62.2% | 82.0% |
| coll% | 0.3% | 0.3% | 0.9% |
| true_blk% | 91.5% | 93.1% | 99.7% |

Bank imbalance (mcf trace, NT=15):
- T0: 69:1, T1: 4:1, T2: 2.4:1 (small tables, group-hash bias)
- T3-T9: ~1.0-1.2:1 (well balanced)
- T10-T14: 1.6-7.3:1 (high tables, index hash bias)

rwram conflicts: 100% writes buffered, ~52% lost (overwritten before flush).

### Notes

- HC_IR's VFS advantage (+0.37%) comes primarily from -7.2% EPI (smaller
  1-bit pred_rams pack tighter) and -2.6% P2 latency, not MPKI.
- Large IR regressions (ren-reactors +443%, roms +73%) are low-absolute-count
  traces where shared providers helped intra-block correlation.
- Bank imbalance in high tables (T10-T14) suggests index hash doesn't spread
  well for long-history tables. Low tables (T0-T2) have group-related bias
  (69:1 for T0) despite BPE1 — likely due to program structure.

## Sweep 9b: BPE1 8-bank Config (2026-05-05)

Moved from 2-bank to 8-bank ta_rwram with BANK_SHIFT=1.
Config: S9_TC_U11_8B (15 tables, UniformTag<11>, GradedSize<512,2048>, 8 banks).

### Monitor baseline (eval_monitor_mb2, 20 traces)

| Metric | Value |
|--------|-------|
| Mean MPKI | 11.432 |
| Tlast_zu% | 27.6% |
| prov_sec_rej% | 24% |

rwram write loss by table (BANK_SHIFT=1, 8 banks, namd trace):

| Table group | Capacity | Loss % | Theoretical |
|-------------|----------|--------|-------------|
| T0-T2 (512) | small | 12.7-13.8% | 12.5% |
| T3 (1024) | mid | 19.9% | 12.5% |
| T4-T6 (2048) | large | 25.7-34.7% | 12.5% |
| T7 (2048) | largest | 39.0% | 12.5% |

Long-history tables show 2-3× theoretical loss rate. Root cause: folded
history register (ta_folded_gh) shifts by only 1 bit per cycle, so
consecutive block indices are highly correlated → bank bits repeat.

## Sweep 11: PC_IDX_SHIFT=5 (2026-05-05)

Hypothesis: TageAhead uses `idx = fold XOR (pc >> 2)` — shifting PC by only
2 leaves intra-block bits in the index. Tage.hpp uses `pc >> (LOG_FETCH_WIDTH + 2)`
= `>> 5` (block-aligned). Try matching that to reduce bank conflicts.

Added `PC_IDX_SHIFT` template parameter (default 2). PCS5 config uses shift=5.
Also shifts tag hash from `>> 4` to `>> 7`.

### Results (6/20 traces completed before kill)

| Trace | Baseline (>>2) | PCS5 (>>5) | Delta |
|-------|---------------|------------|-------|
| 505-mcf | 17.651 | 18.572 | +0.92 |
| 508-namd | 5.013 | 5.061 | +0.05 |
| 548-exchange2 | 5.980 | 7.117 | +1.14 |
| 554-roms | 0.065 | 0.065 | 0.00 |
| gcc-1 | 14.021 | 15.358 | +1.34 |
| llvm-2 | 11.522 | 12.280 | +0.76 |

rwram loss rates essentially unchanged (namd: T6 34.7% vs 34.7% baseline).

**Verdict**: Strictly worse. Shifting PC by 5 discards 3 bits of PC entropy
from the index, increasing aliasing without reducing bank conflicts. The
fold register dominates the index for long-history tables, so changing which
PC bits contribute has no effect on bank correlation. Killed after 6 traces.

## Sweep 10: Bank Shift Decorrelation (2026-05-05)

Motivation: with 8 banks and BANK_SHIFT=1 (baseline), the bank select bits
`index[SHIFT+2:SHIFT]` are highly correlated between consecutive predict/train
cycles for long-history tables (fold register shifts only 1 bit/cycle).
Try varying BANK_SHIFT per table to pick less-correlated index bits.

Max shift per table (constrained by `SHIFT + ceil_log2(banks) <= min(IDX, HYST_IDX)`):
- T0 (512 entries, hyst=256): max shift = 5
- T1-T4 (1024, hyst=512): max shift = 6
- T5-T14 (2048, hyst=1024): max shift = 7

### Single-trace results (505-mcf, 40M instr)

| Config | Shift (T0→T14) | MPKI | T0 loss% | T7 loss% | T12-14 loss% |
|--------|----------------|------|----------|----------|-------------|
| Baseline | uniform 1 | 17.651 | 14.8 | 16.0 | 47.6 |
| BS4 | uniform 4 | 17.636 | 14.4 | 15.5 | 24.2 |
| GBS_531 | 5,5,4,4,3,3,3,2,2,2,1×5 | 17.500 | 17.4 | 16.0 | 47.7 |
| **GBS_542** | **5,5,5,4,4,4,4,3,3,3,2×5** | **17.356** | 17.5 | 16.0 | 34.0 |
| BS_MAX | 5,6,6,6,6,7×10 | 17.827 | 17.8 | 19.7 | 55.2 |
| GBS_510 | 5,4,3,3,2,2,1,1,1,0×6 | 18.053 | 17.8 | 15.7 | 56.9 |

Note: T0 = longest history (512 entries), T14 = shortest history (2048 entries).

Key observations:
- **GBS_542 best MPKI** (−0.3 vs baseline), reduces T12-14 loss 47.6→34.0%
- BS4 (uniform 4) best T12-14 loss (24.2%) but less MPKI gain
- High shift on T0 (long history) consistently hurts — T0 loss rises 14.8→17.5%
- Aggressive shift everywhere (BS_MAX) or shift→0 (GBS_510) both degrade
- Sweet spot: moderate shift on short-history tables, higher on long-history

### Full 20-trace eval (GBS_542 vs baseline shift=1)

| Trace | Baseline | GBS_542 | Delta |
|-------|----------|---------|-------|
| 502-gcc | 17.566 | 17.963 | +0.40 |
| 505-mcf | 17.651 | 17.356 | −0.30 |
| 508-namd | 5.013 | 5.096 | +0.08 |
| 531-deepsjeng | 14.804 | 15.232 | +0.43 |
| 548-exchange2 | 5.980 | 7.017 | +1.04 |
| 554-roms | 0.065 | 0.057 | −0.01 |
| dcapo-kafka | 6.382 | 6.348 | −0.03 |
| gap-sssp | 37.253 | 37.232 | −0.02 |
| gcc-1 | 14.021 | 14.508 | +0.49 |
| java16 | 6.716 | 6.897 | +0.18 |
| llvm-2 | 11.522 | 11.944 | +0.42 |
| lua-3 | 2.632 | 2.639 | +0.01 |
| nodejs-http2 | 9.738 | 9.948 | +0.21 |
| nodejs-octane | 7.361 | 7.506 | +0.15 |
| python3-dulwich | 5.420 | 5.674 | +0.25 |
| rsbench | 42.737 | 42.531 | −0.21 |
| sampleflow | 2.097 | 3.984 | +1.89 |
| web_130 | 6.821 | 7.113 | +0.29 |
| web_74 | 6.475 | 6.625 | +0.15 |
| zstd | 8.390 | 8.457 | +0.07 |
| **Mean** | **11.432** | **11.706** | **+0.27** |

5 wins / 15 losses. Single-trace mcf result was misleading.

**Verdict**: BANK_SHIFT decorrelation hurts overall. While it reduces rwram
write loss on short-history tables (T12-14: 47.6%→34.0%), it changes the
aliasing pattern in ways that degrade accuracy on most traces. The write loss
reduction does not compensate for the aliasing disruption. BANK_SHIFT=1
(baseline) remains optimal.

## FB Banking: Address-banked Fallback RAM (2026-05-05)

Split single 8K×7-bit fb RAM into 4 address banks × 7 per-group 1-bit RAMs
= 28 RAMs of 2K×1-bit. Goal: reduce active SRAM area per cycle, improve
floorplan packing. Bank select from low 2 bits of fb_idx (SHIFT=0).

### Template TageAhead (eval-monitor, 20 traces, FREE_FANOUT)

| Config | Description | MPKI | VFS | vs baseline |
|--------|-------------|------|-----|-------------|
| BPE1_8B | baseline (no FB banking) | 11.432 | — | — |
| FB4 | 4 addr banks, per-group 1-bit, SHIFT=0 | 11.452 | — | +0.18% MPKI |
| FB4_S6 | same but SHIFT=6 | — | — | (killed, same binary bug) |

Accuracy neutral (+0.18%). EPI +5% from arr::select MUX overhead (FREE_FANOUT).

### HC_IR Timing (gcc-1, 100K warmup + 100K measure, no FREE_FANOUT)

| Config | P2 (ns) | EPI (fJ) | Misps |
|--------|---------|----------|-------|
| HC_IR MB2 (no FB banking) | 0.860 | 3056 | 2198 |
| HC_IR MB2 + FB4 (arr::select MUX) | 1.040 | 3655 | 2198 |
| HC_IR MB2 + FB4 (OR-fold execute_if) | 1.037 | 3132 | 2198 |

MUX approach regressed P2 by +0.18ns (+21%). OR-fold (execute_if masking +
fold_or) replaced 2-stage MUX<4> with AND+OR<4>. Result: timing barely
changed (-0.003ns) but energy dropped -14% (3655→3132 fJ). The timing
bottleneck is the 28 RAM reads' wiring, not the selection logic. All 4 banks
still read every cycle (HARCOM cannot gate RAM reads).

**Verdict**: FAILED. FB banking in HC_IR adds +0.18ns P2 regression (+21%)
regardless of MUX approach (arr::select vs OR-fold). The 28 small RAM reads
(4 banks × 7 groups, all read every cycle — HARCOM cannot gate reads)
dominate the critical path vs the single 8K×7 RAM. Reverted to single RAM.
FB banking remains viable in template TageAhead (FREE_FANOUT mode) where
timing is unconstrained.
