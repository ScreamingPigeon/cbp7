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
