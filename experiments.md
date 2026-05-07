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
| HC_IR MB2 + 7×8K×1 (width-split only) | 1.087 | 3290 | 2198 |

All splitting approaches regressed timing. Summary:

| Approach | Extra RAMs | P2 delta | EPI delta |
|----------|-----------|----------|-----------|
| 28 RAMs (4 banks × 7 groups, 2K×1) MUX | +27 | +21% | +27% |
| 28 RAMs OR-fold | +27 | +21% | +9% |
| 7 RAMs (7 groups, 8K×1) | +6 | +26% | +15% |

**Verdict**: FAILED. Any user-level RAM splitting hurts HC_IR timing.
HARCOM's internal SRAM banking (automatic, per harcom.pdf §5.5.2) handles
a single wide RAM more efficiently than multiple narrow user-level RAMs.
The single 8K×7 RAM benefits from internal H-tree routing and shared
decode; splitting into separate objects forces independent address
distribution and read ports. Reverted to single RAM. FB banking remains
viable in template TageAhead (FREE_FANOUT mode) where timing is
unconstrained.

## Sweep 12: GradedTag<11,7> with BPE1 (2026-05-05)

Motivation: with BPE1 GROUP_BITS=3, effective tag bits = raw − 3. Uniform
11-bit tag gives 8 effective bits everywhere. GradedTag<11,7> gives long-
history tables 8 eff bits (same) and short-history tables 4 eff bits (less
disambiguation but smaller tag storage). Prior sweep 9 showed GT11_7 had
best MPKI among graded configs on 14-table non-8B setup.

Config: TA1C_BPE1_8B_GT117 (15 tables, 8 banks, BPE1, GradedTag<11,7>).

### Results (eval-monitor, 20 traces)

| Trace | Base (U11) | GT11_7 | Delta |
|-------|-----------|--------|-------|
| 502-gcc | 17.566 | 18.151 | +0.59 |
| 505-mcf | 17.651 | 15.530 | **−2.12** |
| 508-namd | 5.013 | 4.710 | −0.30 |
| 531-deepsjeng | 14.804 | 15.215 | +0.41 |
| 548-exchange2 | 5.980 | 6.490 | +0.51 |
| 554-roms | 0.065 | 0.057 | −0.01 |
| dcapo-kafka | 6.382 | 6.329 | −0.05 |
| gap-sssp | 37.253 | 37.009 | −0.24 |
| gcc-1 | 14.021 | 14.479 | +0.46 |
| java16 | 6.716 | 6.957 | +0.24 |
| llvm-2 | 11.522 | 11.867 | +0.35 |
| lua-3 | 2.632 | 2.629 | −0.00 |
| nodejs-http2 | 9.738 | 10.036 | +0.30 |
| nodejs-octane | 7.361 | 7.662 | +0.30 |
| python3-dulwich | 5.420 | 5.648 | +0.23 |
| rsbench | 42.737 | 41.434 | **−1.30** |
| sampleflow | 2.097 | 1.968 | −0.13 |
| web_130 | 6.821 | 7.131 | +0.31 |
| web_74 | 6.475 | 6.619 | +0.14 |
| zstd | 8.390 | 8.836 | +0.45 |
| **Mean** | **11.432** | **11.438** | **+0.05%** |

7 wins / 13 losses. Net effect negligible (+0.05%).

### Per-table analysis (why regressions happen)

**mcf (GT117 wins −2.12)**: Baseline had excessive false tag matches on
mid-history tables T5-T9 (TagMatch 27-33%, accuracy only 37-40%). GT117's
shorter tags on these tables changed allocation dynamics, collapsing their
provider counts (T5: 75K→4.7K, T6: 84K→3K). Branches fell back to bimodal
which was more accurate (64.6%→66.3%). Shorter tags acted as a beneficial
filter — fewer false matches = less pollution from bad providers.

**gcc-1 (GT117 loses +0.46)**: TagMatch% rose across tables (T14: 5.2%→11.2%
with only 7-bit raw / 4-bit eff tag). T0 provided more predictions
(5365→7213) but with worse accuracy (10.6%→8.4%). The shorter tags on
short-history tables increased aliasing without benefit — these tables had
well-calibrated tag matching at 11 bits.

**Pattern**: GT117 helps traces where baseline has excessive false tag matches
on mid-tables (aliasing pollution), but hurts traces where the 11-bit uniform
tag was already well-calibrated. The effect is trace-dependent with no net
gain.

**Verdict**: GradedTag<11,7> is accuracy-neutral at aggregate. The prior
sweep 9 MPKI advantage (12.357 vs 12.556) was measured on a different config
(14 tables, no 8-bank). With 15 tables + 8 banks, the advantage disappears.

## Sweep 12: Per-Bit Flip Rate Instrumentation (Bank Conflict Root Cause)

**Goal**: Identify exactly which index bit positions are "sticky" (low flip rate
between consecutive read accesses) for short-history tables, causing bank conflicts.

**Setup**: Added `bit_flip_count[A]` instrumentation to ta_rwram. Tracks which
specific bit positions of the full index flip between consecutive `read()` calls,
per table. Run on llvm-2 trace (1000 warmup, 40000 measure).

**Config**: TageAhead1C (15 tables, 8 banks, BANK_SHIFT=1, BANK_BITS=3, so bank
select = index bits [1:3]). Note: default config already uses BANK_SHIFT=1.

### Results: Per-table bit flip rates (bank-select bits highlighted)

| Table | HistLen | same_prev% | lost% | **b0 flip** | **b1 flip** | **b2 flip** | b3 | b4 | b5 | idx_xor_0% |
|-------|---------|-----------|-------|------------|------------|------------|-----|-----|-----|------------|
| T0  | 199 | 50.3% | 50.5% | 50.9% | 49.8% | 49.7% | 50.3% | 49.6% | 47.4% | 0.1% |
| T1  | 158 | 49.0% | 50.6% | 49.4% | 51.0% | 50.6% | 50.1% | 49.8% | 48.8% | 0.1% |
| T2  | 126 | 50.7% | 51.3% | 48.6% | 49.3% | 49.4% | 49.7% | 49.0% | 48.8% | 0.1% |
| T3  | 100 | 50.2% | 50.9% | 49.4% | 49.8% | 50.0% | 49.6% | 49.6% | 49.1% | 0.1% |
| T4  |  79 | 50.0% | 50.0% | 50.3% | 50.0% | 50.2% | 50.4% | 49.8% | 49.3% | 0.1% |
| T5  |  63 | 51.0% | 50.1% | 49.2% | 49.0% | 51.1% | 49.1% | 49.2% | 49.8% | 0.1% |
| T6  |  50 | 50.3% | 49.4% | 50.7% | 49.7% | 49.3% | 49.0% | 49.7% | 47.7% | 1.1% |
| T7  |  39 | 52.1% | 52.6% | 48.7% | 47.9% | 50.1% | 49.1% | 49.0% | 48.4% | 2.0% |
| T8  |  31 | 53.0% | 52.2% | 46.0% | 47.0% | 50.5% | 48.9% | 49.9% | 46.3% | 2.6% |
| T9  |  25 | 53.0% | 52.2% | 46.0% | 47.0% | 50.5% | 48.9% | 49.9% | 46.3% | 3.2% |
| T10 |  20 | 51.6% | 53.7% | 48.8% | 48.4% | 49.0% | 48.9% | 49.9% | 46.3% | 3.6% |
| T11 |  15 | 51.6% | 53.7% | 48.8% | 48.4% | 49.0% | 48.9% | 44.1% | 41.8% | 4.0% |
| T12 |  12 | 56.7% | 57.3% | 48.8% | 43.3% | 48.9% | 47.2% | 44.1% | 41.8% | 4.2% |
| T13 |  10 | 56.7% | 57.3% | **3.7%** | 43.3% | 48.9% | 47.2% | 44.1% | 41.8% | 4.4% |
| T14 |   8 | 56.7% | 57.3% | **3.7%** | 43.3% | 48.9% | 47.2% | 44.1% | 41.8% | 4.7% |

### Key Findings

1. **Bit 0 is nearly frozen for T13-T14**: Only 3.7% flip rate (vs ~50% ideal).
   With BANK_SHIFT=0, this is the LSB of bank select — consecutive accesses
   almost always land in the same b0-parity half of the banks.

2. **Bit 1 degraded for T12-T14**: 43.3% flip rate (vs ~50% ideal). Combined
   with frozen b0, the 3-bit bank select is heavily biased toward repeating.

3. **Higher bits also degrade for short tables**: b4/b5 drop to 41-44% for
   T11-T14, and b8/b9 drop to 25-35% for T14. The short fold registers have
   low entropy overall.

4. **Long-history tables are fine**: T0-T5 show ~50% flip rates on all bits,
   confirming bank conflicts are concentrated in short-history tables.

5. **Root cause**: For 8-10 bit history lengths, the fold register is tiny and
   only shifts by 1 bit per cycle. Combined with `pc >> 2` (which is constant
   within a basic block loop), the low index bits barely change. The same index
   → same bank → same-bank conflict → lost writes (57% loss rate for T13-T14).

### Corrected Understanding

The actual config uses **BANK_SHIFT=1** (confirmed in both S9_TC_U11_8B and
TageAheadHC_IR), so bank select = bits [1,2,3]. The frozen b0 (3.7%) is already
excluded from bank selection. The bit flip rates for the actual bank-select bits are:

- b1: 43-50% (slightly degraded for T12-T14)
- b2: 48-51% (near-ideal)
- b3: 47-50% (decent)

However, `bank_xor_popcount` shows only values 0 and 1 (never 2 or 3), meaning
**at most 1 of the 3 bank bits changes per cycle**. This is the fundamental issue:
the fold register shifts by 1 bit/cycle, so the full index changes by only 1-2 bits
total between consecutive accesses. Even though individual bits have ~50% long-run
flip rates, they don't flip independently — only one flips at a time.

This produces ~50% same-bank rate across ALL tables (not just short ones). Short
tables (T12-T14) are slightly worse (56-57%) due to additional low-entropy effects.

### Architecture Summary

- 1 read per table per cycle (single index from `fold XOR pc>>2`)
- Each entry has a 3-bit `group_id` in the tag identifying which of 7 branches it serves
- 7 independent resolution chains share the physical tables
- Conflict: predict-read on cycle N and train-write on cycle N+K hit same bank
  because consecutive indices differ by only 1-2 bits → same bank bits

### Literature Fix (Seznec, "A New Case for TAGE", JILP 2011)

Published TAGE implementations **decouple bank selection from the index entirely**.
Bank is derived from PC rotation:

```
bank = pc & (NUM_BANKS-1);
if (bank == prev_bank) bank = (bank+1) & (NUM_BANKS-1);
```

Consecutive branch PCs are different addresses → naturally spread across banks.
The folded history is only used for row index + tag, never for bank selection.
Updates use small write buffers (1-2 cycle tolerance) since silent-update
elimination means only ~9% of branches actually need writes.

### Implications

- BANK_SHIFT adjustments are ineffective — the problem is that the fold register
  produces low-hamming-distance consecutive indices regardless of which bits are used.
- The correct fix: derive bank from PC bits (or a rotation of PC), independent of
  the folded history index. This guarantees consecutive predictions hit different
  banks, leaving previous banks free for buffered write flushes.

### Experiment: Block PC XOR into bank select (FAILED)

**Approach**: XOR `val<BANK_BITS>{inst_pc >> 5}` into the bank_id after split_addr
in ta_rwram. Both read and write use the same XOR so they target the same physical
bank. Added `bank_xor` parameter to `ta_rwram::read()` and `ta_rwram::write()`.

**Result on llvm-2**: No improvement. Same ~50% same-bank rate, ~50-57% lost writes
across all tables. Identical to baseline.

**Why it failed**: In a loop, the block PC is **constant** across iterations. XOR-ing
a constant into bank select just permutes which bank is used — it doesn't change that
consecutive accesses from the same block hit the **same** permuted bank. The conflict
comes from same-block re-accesses, not from inter-block transitions.

**Conclusion**: Static PC-based bank derivation cannot help when the dominant conflict
pattern is loops hitting the same PC repeatedly. Need **stateful rotation** (Seznec
approach) where a running counter or prev_bank tracker guarantees different banks on
consecutive cycles regardless of PC repetition.

## Sweep 13: SEC_TAG Width (2026-05-05)

Motivation: sec_tag RAMs are per-table (15 tables). Reducing width saves
area/power. Baseline uses 5-bit sec_tag (Xor3SecTagHash5) which rejects
~24% of tag matches as false positives.

Config: BPE1, 8 banks, 15 tables, UniformTag<11>. Swept SEC_TAG_BITS=3,4
vs baseline 5.

### Results (eval-monitor, 20 traces)

| Trace | SEC5 (base) | SEC4 | SEC3 |
|-------|------------|------|------|
| 502-gcc | 17.566 | 18.082 | 18.247 |
| 505-mcf | 17.651 | 18.017 | 17.852 |
| 508-namd | 5.013 | 4.957 | 5.039 |
| 531-deepsjeng | 14.804 | 15.310 | 15.394 |
| 548-exchange2 | 5.980 | 7.170 | 7.155 |
| 554-roms | 0.065 | 0.092 | 0.088 |
| dcapo-kafka | 6.382 | 7.316 | 7.616 |
| gap-sssp | 37.253 | 37.569 | 38.192 |
| gcc-1 | 14.021 | 14.528 | 14.477 |
| java16 | 6.716 | 7.139 | 7.162 |
| llvm-2 | 11.522 | 11.906 | 12.319 |
| lua-3 | 2.632 | 2.422 | 2.659 |
| nodejs-http2 | 9.738 | 10.001 | 10.112 |
| nodejs-octane | 7.361 | 7.563 | 7.855 |
| python3-dulwich | 5.420 | 5.726 | 5.750 |
| rsbench | 42.737 | 42.438 | 42.499 |
| sampleflow | 2.097 | 2.179 | 1.956 |
| web_130 | 6.821 | 7.166 | 7.347 |
| web_74 | 6.475 | 6.636 | 6.688 |
| zstd | 8.390 | 7.967 | 8.332 |
| **Mean** | **11.432** | **11.709 (+2.4%)** | **11.837 (+3.5%)** |

| Config | Mean MPKI | sec_rej% | Wins/Losses |
|--------|-----------|----------|-------------|
| SEC5 (baseline) | 11.432 | 24.0% | — |
| SEC4 | 11.709 | 23.4% | 3/17 |
| SEC3 | 11.837 | 22.9% | 2/18 |

**Verdict**: FAILED. Shorter sec_tag consistently degrades accuracy. The
5-bit sec_tag rejects 24% of tag matches — real filtering work that prevents
false-positive providers. Reducing to 4 bits loses 0.6% of that filtering
(23.4% sec_rej) and costs +2.4% MPKI. Worst regressions: exchange2 (+1.2),
kafka (+0.9-1.2), deepsjeng (+0.5). 5-bit sec_tag is optimal.

---

## Experiment: Seznec-Style Bank Rotation (FAILED)

**Goal**: Reduce rwram read-write bank conflicts. Baseline: 2-bank single-ported
tables, ~50% same-bank rate between consecutive reads, ~57% lost writes.

**Approach**: Per-table stateful rotation — track previous read's bank, if
natural bank == prev_bank then increment. XOR the natural vs rotated bank to
produce `bank_xor` passed into `ta_rwram::read()` and `ta_rwram::write()`.
Pipeline the rotated bank through prefetch→current→train for write path.

**Result on llvm-2** (40M branches):
- `read_same_prev`: 56.68% → **0.00%** (rotation works perfectly for read-read)
- `lost writes`: 57% → **78%** (WORSE)
- MPKI: 11.52 → **22.71** (+97% regression)

**Root Cause**: In a banked RAM, `bank + localaddr` together identify a unique
entry. XOR-ing the bank bit routes the read to a **different physical entry**
(same localaddr, wrong bank). Tag check fails → massive miss rate increase.
The rotation doesn't just change routing — it changes which data you read.

**Why lookahead scheduling (EV8/Seznec papers) also doesn't apply**: Those
techniques spread *independent* accesses to different banks. In hot loops,
both prediction reads and training writes target the **same entry** (same
block PC → same index → same bank). No scheduling can avoid the conflict
when the read and write are to the same address.

**Conclusion**: For 2-bank single-ported RAMs with hot-loop workloads, ~50%
read-write bank conflict rate is the architectural floor. The lost-write rate
is already priced into baseline MPKI. Possible mitigations (not pursued):
deeper write buffer, write-only-on-mispredict, or accepting the loss.

**Reverted**: All rotation code removed from TageAhead.hpp and custom_common.hpp.

## Sweep 14: GradedTag<11,7> on TageAheadHC_IR (2026-05-05)

Motivation: GradedTag<11,7> was accuracy-neutral (+0.05%) on template TageAhead
BPE1 (Sweep 12), but tag RAM area savings could improve VFS on the hand-coded
HC_IR. Short-history tables (T0-T3) use 7-bit tags, long-history tables (T13)
use 11-bit. With GROUP_BITS=3, effective htag ranges from 4 to 8 bits.

Config: TageAheadHC_IR with per-table tag widths:
  {7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11}

### Quick-eval results (20 traces)

| Metric | Baseline (U11) | GT<11,7> | Delta |
|--------|----------------|----------|-------|
| MPI | 0.01081 | 0.01078 | -0.3% |
| VFS | 0.878 | 0.8777 | -0.03% |

MPKI slightly better but VFS marginally worse — tag RAM savings don't
compensate. Full eval (168 traces) launched to `out/full_hcir_gt117/` for
definitive comparison, since the quick-eval delta is within noise.

### Full-eval results (55 traces completed)

GT<11,7> weighted MPKI: 9.783 vs HC_IR baseline: 9.658 → **+1.29% worse**.
Several large regressions (compress_7 +21%, fp_11 +25%, roms +45%) outweigh
wins (exchange2 -12%, gmsh -13.5%). **Verdict: FAILED.**

---

## Experiment: FB Bimodal Hysteresis on TageAheadHC_IR (2026-05-05)

**Goal**: Add traditional 2-bit saturating counter hysteresis to the fallback
bimodal predictor. Instead of overwriting fb prediction on every mispredict,
only flip when hysteresis is weak (wrong twice). Half-sized RAM (4096×7).

**Key constraint**: HARCOM single-RAM-access-per-cycle. Read fb_bim_hyst gated
by mispredict before `need_extra_cycle` (cycle 1), write after (cycle 2).
Only updates on mispredict cycles.

Config: TageAheadHC_IR, 14 tables, BPE1, 8 banks, SEC_TAG=5.
FB_BH_CAPACITY=4096 (half of FB_CAPACITY=8192).

### Quick-eval results (20 traces)

| Metric | HC_IR baseline | HC_IR + FBH | Delta |
|--------|----------------|-------------|-------|
| MPI | 0.01078 | 0.009804 | -9.1% |
| VFS | 0.878 | 0.8946 | +1.9% |

VFS stays under 1 cycle (0.913 with FLOORPLAN). fb_bim_hyst RAM placed after
meta to avoid displacing critical-path structures.

Full eval (168 traces, 4 cores) launched to `full_out/hc_ir_fbh/`.

## Experiment: Fanout Optimization on TageAheadHC_IR (2026-05-06)

Applied fanout fixes (fo1/fanout(hard<N>) on critical path). Full 168-trace eval.

### Full-eval results

| Config | VFS | MPKI | EPI (fJ) | IPC_cbp | Max P2 (ns) |
|--------|-----|------|----------|---------|-------------|
| HC_IR FBH | 0.897064 | 8.964 | 2837 | 9.009 | 0.913 |
| HC_IR fo_fix | 0.897028 | 8.964 | 2837 | 9.009 | 0.993 |
| Delta | -0.004% | 0.0% | 0.0% | 0.0% | +8.8% |

Fanout fixes increased P2 from 0.913→0.993ns (still under 1 cycle). VFS unchanged
because P2 is capped at ceil(1) = 1 cycle in both cases. MPKI and EPI identical.

### VFS analysis: high-IPC traces penalized

VFS formula has `(IPC/IPC*)^3.2` energy term. With reference IPC*=8, traces with
IPC >> 8 get exponentially penalized. 7 traces with IPC > 20 score VFS < 0.65:

| Trace | IPC | MPKI | VFS |
|-------|-----|------|-----|
| roms | 58.0 | 0.191 | 0.009 |
| fp_28 | 60.8 | 1.467 | 0.028 |
| fp_21 | 23.0 | 0.677 | 0.178 |
| roms-1 | 25.5 | 1.506 | 0.195 |
| int_197 | 22.7 | 1.261 | 0.232 |
| diamond | 21.0 | 2.067 | 0.355 |
| namd | 20.1 | 4.905 | 0.634 |

These traces have excellent accuracy but terrible VFS because the predictor is "too fast."

## Experiment: IPC Capping via need_extra_cycle (2026-05-07)

### Why IPC throttling matters (from vfs.pdf, Michaud 2025)

The VFS formula models a **fixed core power budget** (Section 5.3). The core C
being evaluated has power `P = EPI × IPC × f₀` at reference voltage. If P > P*
(reference power), voltage and frequency are **scaled down** to meet the budget.
If P < P*, they are **scaled up**. The VFS captures this frequency adjustment:

```
VFS = (IPC/IPC*) × α × (1 - 2/(1 + sqrt(1 + β × (EPI*/EPI) × (IPC*/IPC))))
```

The key insight is **how core energy scales with IPC** (Eq. 12, Section 5.4):

```
EPI/EPI* = (IPC/IPC*)^γ     where γ = 3.2
```

This models that higher IPC requires a more complex core consuming more energy
(qualitatively similar to Pollack's rule but with a larger exponent). The
normalized EPI (Eq. 16) combines CBP energy and core energy:

```
EPI/EPI* = [λ × (IPC/IPC*)^γ + μ × EPIcbp/EPIcbp*] × (1 + WPI/2)
```

where `λ ≈ 0.826` (core fraction) and `μ = 0.05` (CBP fraction). The core
energy term `λ × (IPC/IPC*)^γ` dominates. For IPC=65 (roms): `(65/8)^3.2 ≈ 23000`,
making EPI/EPI* ≈ 19000. This enormous normalized energy forces the VFS
voltage-frequency scaling to **massively reduce frequency**, collapsing VFS
despite the raw IPC speedup.

The VFS-optimal curve (Section 5.4) shows that `P/P* = (IPC/IPC*)^(γ+1)`.
Cores above this curve consume too much power for their IPC — the frequency
reduction eats the speedup. **A predictor with IPC >> IPC* lands far above
the VFS-optimal curve.** Throttling IPC back toward IPC*=8 moves the
operating point onto the curve where VFS is maximized.

From Section 4: "Do not apply the VFS formula on a single trace." The FoM
are averaged (IPC_cbp harmonic, CPI_cbp/EPI_cbp arithmetic) then VFS computed
once. But `need_extra_cycle()` operates per-trace during simulation, so the
throttling must be conservative — false positives on moderate-IPC traces
hurt the aggregate.

**Goal**: Artificially cap IPC to 8 (= reference predictor) on high-IPC, low-MPKI
windows by injecting extra cycles via `need_extra_cycle()`. This moves
high-IPC traces toward the VFS-optimal operating point.

### Mechanism (HARCOM implementation)

- Track per-window (2^WIN_BITS blocks): instruction count, block count, mispredictions
- At window boundary: if `MPKI < 1000/2^MPKI_SHIFT` AND `IPC > 2^IPC_CAP_SHIFT`,
  compute debt = `instr >> IPC_CAP_SHIFT - window_blocks`
- Each block: drain up to MAX_INJECT extra cycles from debt via `need_extra_cycle()`
- Only inject on non-mispredicting blocks (`& ~mispredict`)
- `need_extra_cycle()` timing is IGNORED — no P1/P2 impact

Parameters (sweepable via -DCAP_WIN_BITS, -DCAP_MPKI_SHIFT, -DCAP_MAX_INJ):
- IPC_CAP_SHIFT: IPC cap = 2^N (default 3 → cap=8)
- WIN_BITS: window = 2^N blocks
- MPKI_SHIFT: MPKI threshold = 1000/2^N
- MAX_INJECT: max extra cycles per block (should be ≥ 31 for full drain)

### Theoretical analysis (full 168-trace simulation)

| Cap strategy | Arith VFS | Geo VFS | Harm VFS | Worst Δ |
|---|---|---|---|---|
| No cap | 0.790 | 0.738 | 0.448 | — |
| MPKI<4, cap=8 | 0.831 | 0.813 | 0.792 | +0.000 |
| IPC>14, cap=8 | 0.829 | 0.812 | 0.791 | +0.000 |
| IPC>10, cap=8 | 0.831 | 0.814 | 0.793 | -0.057 |

MPKI gating critical: without it, traces with high IPC + high MPKI (nest-1, gmsh-5)
regress because reducing IPC hurts when CPI is already high.

### Quick-eval sweep: WIN=128, varying MPKI threshold (MAX_INJ=32)

CSV column format (from predictor_metrics.py):
`name,instr,branch,condbr,pred_cycles,extra_cycles,diverge,diverge_at_end,misp,p1,p2,epi`

VFS computed using exact vfs.py formula (speedup includes WPI correction,
normalizedEPI uses `LAMBDA * speedup^GAMMA`, not `(IPC/IPC*)^3.2` directly).

16/20 quick-eval traces overlap with full_out baseline. Per-trace vs fo_fix baseline:

| Trace | B_IPC | MPKI | B_VFS | w128_m4 VFS | w128_m2 VFS | w128_m1 VFS |
|-------|-------|------|-------|-------------|-------------|-------------|
| 502-gcc-all | 6.6 | 14.6 | 0.704 | 0.697 (-0.007) | 0.700 (-0.005) | 0.701 (-0.003) |
| 505-mcf | 12.8 | 17.2 | 0.793 | 0.791 (-0.002) | 0.791 (-0.002) | 0.791 (-0.002) |
| **508-namd** | 23.0 | 4.9 | 0.545 | **0.968 (+0.423)** | **0.931 (+0.386)** | **0.870 (+0.324)** |
| 531-deepsjeng | 12.3 | 11.8 | 0.881 | 0.880 (-0.001) | 0.881 (-0.001) | 0.881 (-0.000) |
| **548-exchange2** | 13.6 | 4.9 | 0.860 | **0.963 (+0.103)** | **0.935 (+0.075)** | **0.907 (+0.047)** |
| **554-roms** | 65.2 | 0.2 | 0.006 | **0.999 (+0.993)** | **0.999 (+0.992)** | **0.975 (+0.969)** |
| **dcapo-kafka** | 11.2 | 4.3 | 0.927 | **0.968 (+0.042)** | **0.970 (+0.043)** | **0.972 (+0.045)** |
| gap-sssp | 10.0 | 31.6 | 0.548 | 0.548 (-0.001) | 0.548 (-0.001) | 0.544 (-0.004) |
| gcc-1 | 6.9 | 10.9 | 0.787 | 0.779 (-0.009) | 0.781 (-0.006) | 0.783 (-0.005) |
| **java16** | 11.5 | 4.0 | 0.910 | **0.977 (+0.067)** | **0.968 (+0.058)** | **0.956 (+0.046)** |
| llvm-2 | 7.5 | 8.7 | 0.858 | 0.851 (-0.008) | 0.854 (-0.005) | 0.855 (-0.004) |
| nodejs-http2 | 7.8 | 7.2 | 0.904 | 0.881 (-0.023) | 0.883 (-0.021) | 0.886 (-0.019) |
| nodejs-octane | 7.7 | 6.0 | 0.921 | 0.884 (-0.037) | 0.886 (-0.034) | 0.888 (-0.033) |
| rsbench | 6.2 | 40.7 | 0.407 | 0.407 (-0.000) | 0.407 (-0.000) | 0.406 (-0.001) |
| **sampleflow** | 10.5 | 3.3 | 0.932 | **0.980 (+0.048)** | **0.980 (+0.048)** | **0.980 (+0.048)** |
| **zstd** | 18.2 | 7.5 | 0.814 | **0.942 (+0.128)** | **0.940 (+0.126)** | **0.938 (+0.124)** |

| Aggregate (16 traces) | Baseline | w128_m4 | w128_m2 | w128_m1 |
|---|---|---|---|---|
| Arith VFS | 0.737 | 0.845 (+0.107) | 0.841 (+0.103) | 0.833 (+0.096) |
| Geo VFS | 0.568 | 0.824 (+0.256) | 0.821 (+0.253) | 0.814 (+0.246) |
| Improved (>0.001) | — | 7 | 7 | 7 |
| Regressed (>0.001) | — | 7 | 6 | 7 |

**IMPORTANT: Per-trace VFS is WRONG.** vfs.pdf Section 4: "Do not apply the VFS
formula on a single trace." The correct procedure: average IPC_cbp (harmonic),
CPI_cbp (arithmetic), EPI_cbp (arithmetic) across all traces, then compute VFS
once on the aggregate FoM. The per-trace VFS numbers above are illustrative only.

### Aggregate VFS (official predictor_metrics.py + vfs.py, same 20 traces)

| Config | IPC_cbp | CPI_cbp | EPI_cbp | VFS | dVFS |
|--------|---------|---------|---------|-----|------|
| Baseline (20 traces) | 9.738 | 0.0872 | 2600 | **0.896** | — |
| w128_m4 (MPKI<3.9) | 8.021 | 0.0872 | 2629 | 0.866 | **-0.030** |
| w128_m2 (MPKI<1.95) | 8.190 | 0.0872 | 2629 | 0.870 | **-0.026** |
| w128_m1 (MPKI<0.98) | 8.336 | 0.0872 | 2629 | 0.873 | **-0.022** |
| w256_m4 (MPKI<3.9) | 8.086 | 0.0872 | 2632 | 0.867 | **-0.028** |

**IPC capping HURTS aggregate VFS.** All configs are worse than baseline.

Root causes:
1. **Harmonic mean IPC**: roms dropping from IPC=65→8 barely moves the harmonic
   mean (1/65=0.015 → 1/8=0.125 is small vs 20 traces). But false-positive
   capping on IPC<8 traces (lua 7.4→6.8, nodejs-octane 7.7→6.6) significantly
   increases their 1/IPC contribution, dragging the harmonic mean down.
2. **False positives on low-IPC traces**: traces with IPC 5-8 have per-window
   MPKI dipping below threshold on some windows, triggering capping that
   shouldn't fire. These traces dominate the harmonic mean.
3. **CPI unchanged**: capping doesn't reduce mispredictions, so CPI stays at
   0.0872. The IPC drop is pure loss — no compensating CPI benefit.
4. **VFS formula's WPI correction**: speedup = (IPC_cbp/IPC*) × (1+WPI*)/(1+WPI).
   Since CPI is unchanged, WPI = IPC_cbp × CPI_cbp drops with IPC_cbp, which
   helps the (1+WPI*)/(1+WPI) factor. But not enough to offset the IPC loss.

**Verdict**: IPC capping via `need_extra_cycle()` is **counterproductive** for VFS.
The aggregate harmonic mean IPC is dominated by low-IPC traces, not the few
high-IPC outliers. Capping high-IPC traces provides negligible harmonic mean
benefit while false positives on moderate traces cause real harm. The energy
model's `(IPC/IPC*)^γ` penalty is already handled by the VFS formula's
voltage-frequency scaling — manually throttling IPC double-counts the adjustment.

### Full sweep results (2026-05-07, all 8 configs complete)

All configs evaluated on 20 quick-eval traces. VFS computed using official
`predictor_metrics.py` + `vfs.py` on matching 20-trace baseline subset.

| Config | Window | MPKI < | IPC_cbp | VFS | dVFS |
|--------|--------|--------|---------|-----|------|
| **Baseline** | — | — | 9.738 | **0.896** | — |
| w512_m1 | 512 | 0.98 | 8.457 | 0.877 | -0.019 |
| w512_m2 | 512 | 1.95 | 8.336 | 0.874 | -0.022 |
| w128_m1 | 128 | 0.98 | 8.336 | 0.873 | -0.022 |
| w256_m2 | 256 | 1.95 | 8.275 | 0.872 | -0.024 |
| w256_m1 | 256 | 0.98 | — | 0.870 | -0.025 |
| w128_m2 | 128 | 1.95 | 8.190 | 0.870 | -0.026 |
| w256_m4 | 256 | 3.9 | 8.086 | 0.867 | -0.028 |
| w128_m4 | 128 | 3.9 | 8.021 | 0.866 | -0.030 |

**Every config hurts VFS.** Best (w512_m1) still -0.019 below baseline.

### Why IPC capping can never help (harmonic mean analysis)

Full 168-trace IPC distribution from `full_out/hc_ir_fo_fix`:

| IPC range | Traces | sum(1/IPC) | % of harmonic sum |
|-----------|--------|------------|-------------------|
| 50+ | 2 | 0.031 | 0.2% |
| 20-50 | 6 | 0.238 | 1.3% |
| 12-20 | 28 | 2.021 | 10.8% |
| 10-12 | 22 | 1.996 | 10.7% |
| 8-10 | 52 | 5.876 | 31.5% |
| <8 | 58 | 8.488 | **45.5%** |

VFS uses **harmonic mean IPC** (vfs.pdf §4, predictor_metrics.py). The harmonic
mean = N / sum(1/IPC_i). High-IPC outliers (roms IPC=65, fp_28 IPC=66)
contribute 1/65 ≈ 0.015 each — **invisible** to the harmonic sum (0.2% combined).

Even with **perfect capping** (zero false positives, cap all IPC>8 to exactly 8):
- Harmonic mean IPC: 9.009 → **7.555** (−16%)
- The 52 traces with IPC 8-12 each go from 1/10≈0.10 to 1/8=0.125
- This adds far more to the harmonic denominator than the high-IPC traces save

110 traces have IPC > 8 but only contribute 54.5% of the harmonic sum. Capping
them all to 8 increases their contribution, dragging the harmonic mean **below** 8.
Meanwhile the 58 traces with IPC < 8 (45.5% of the sum) are untouched.

**Conclusion**: IPC throttling is fundamentally incompatible with harmonic-mean-
based VFS scoring. The approach is a dead end — no parameter tuning can fix it.
The correct path to VFS improvement is reducing MPKI/CPI and EPI, not throttling IPC.

### Per-trace FoM changes (20 quick-eval traces)

CPI is unchanged across all configs (same mispredictions, same P2 latency).
Showing best config (w512_m1) and most aggressive (w128_m4) vs baseline.

| Trace | B_IPC | B_CPI | B_EPI | m1_IPC | dIPC% | m1_EPI | dEPI | m4_IPC | dIPC% | m4_EPI | dEPI |
|-------|-------|-------|-------|--------|-------|--------|------|--------|-------|--------|------|
| 502-gcc-all | 6.6 | 0.132 | 3829 | 6.6 | -0.3% | 3881 | +52 | 6.5 | -2.7% | 3874 | +45 |
| 505-mcf | 12.8 | 0.154 | 1783 | 12.6 | -1.4% | 1811 | +28 | 12.6 | -1.8% | 1807 | +24 |
| 508-namd | 23.0 | 0.044 | 1099 | 18.5 | -19.8% | 1115 | +16 | 10.3 | -55.4% | 1113 | +14 |
| 531-deepsjeng | 12.3 | 0.107 | 1960 | 12.3 | -0.0% | 1986 | +26 | 11.9 | -3.2% | 1982 | +22 |
| 548-exchange2 | 13.6 | 0.044 | 1976 | 13.1 | -3.7% | 2007 | +31 | 10.0 | -26.9% | 2003 | +27 |
| **554-roms** | **65.2** | 0.002 | 433 | **8.4** | **-87.1%** | 439 | +6 | **8.1** | **-87.6%** | 438 | +5 |
| dcapo-kafka | 11.2 | 0.039 | 2337 | 8.2 | -26.1% | 2364 | +27 | 7.7 | -31.3% | 2361 | +24 |
| gap-sssp | 10.0 | 0.285 | 2024 | 9.7 | -2.4% | 2054 | +30 | 9.9 | -0.1% | 2050 | +26 |
| gcc-1 | 6.9 | 0.099 | 3713 | 6.9 | -0.7% | 3760 | +47 | 6.7 | -3.0% | 3753 | +40 |
| java16 | 11.5 | 0.036 | 2311 | 10.5 | -9.2% | 2342 | +31 | 8.7 | -24.5% | 2338 | +27 |
| llvm-2 | 7.5 | 0.078 | 3475 | 7.4 | -0.5% | 3519 | +44 | 7.2 | -2.9% | 3513 | +38 |
| lua-3 | 7.4 | 0.019 | 3586 | 6.9 | -7.0% | 3624 | +38 | 6.8 | -7.4% | 3618 | +32 |
| nodejs-http2 | 7.8 | 0.064 | 3292 | 7.3 | -7.2% | 3329 | +37 | 7.1 | -9.5% | 3324 | +32 |
| nodejs-octane | 7.7 | 0.054 | 3426 | 6.7 | -12.0% | 3468 | +42 | 6.6 | -14.0% | 3462 | +36 |
| python3-dulwich | 8.4 | 0.033 | 3146 | 7.6 | -9.9% | 3182 | +36 | 7.3 | -13.4% | 3177 | +31 |
| rsbench | 6.2 | 0.366 | 3497 | 6.2 | -0.2% | 3549 | +52 | 6.2 | +0.0% | 3541 | +44 |
| sampleflow | 10.5 | 0.030 | 2639 | 8.2 | -22.0% | 2681 | +42 | 7.8 | -25.0% | 2676 | +37 |
| web_130 | 8.7 | 0.044 | 3045 | 8.0 | -8.1% | 3083 | +38 | 7.5 | -13.0% | 3078 | +33 |
| web_74 | 8.6 | 0.046 | 3057 | 7.4 | -14.3% | 3095 | +38 | 7.3 | -14.9% | 3089 | +32 |
| zstd | 18.2 | 0.068 | 1369 | 12.4 | -31.8% | 1390 | +21 | 11.5 | -36.5% | 1387 | +18 |

Key observations:
- CPI identical across all configs (capping doesn't affect mispredictions)
- EPI increases slightly (+5 to +52 fJ) due to `need_extra_cycle()` calling
  `panel.next_cycle()` which adds dynamic energy
- False positives: lua (IPC 7.4, already <8) gets capped to 6.9/6.8;
  nodejs-octane (IPC 7.7) capped to 6.7/6.6 — these should never be touched
- Even traces with IPC ~8-10 (web_74, web_130, python3) get 7-14% IPC cuts

Results in `out/ipc_cap_sweep/{config}/`. Baseline in `full_out/hc_ir_fo_fix/`.
