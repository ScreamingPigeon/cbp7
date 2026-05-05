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
