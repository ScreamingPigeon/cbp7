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
