# TageAheadHC Floorplan Experiments

Baseline: S3_MW4_MC2048 (template) — EPI=5523, P2=1.02333, misps=12627 (1M instr, gcc trace)

## Experiment 1: Type-grouped layout
**Order:** all tag_rams, all pred_rams, all sec_rams, all fold_idx, all fold_tag, fb, zone, meta, all hyst_rams, all u_rams
**Result:** EPI=7399, P2=1.04, misps=13188
**Timing (T0 cycle 2):** fold_idx@loc55 → stored_tag 593ps → pred[0] 1100ps
**Analysis:** Folds at loc=55 (last sec_ram), tag_ram0 at loc=0 — massive wiring delay (90ps vs S3's 503ps stored_tag). Worst layout.

## Experiment 2: Per-table clusters with zones (TATable order)
**Order per table:** tag, pred, sec, fold_idx, fold_tag, zone, hyst, u
**Result:** EPI=4917, P2=0.94, misps=12627
**Timing (T0 cycle 2):** fold_idx@loc3 → stored_tag 499ps → pred[0] 827ps
**Analysis:** Folds co-located with their table's predict RAMs (after sec_ram). 11% EPI improvement over template, 8% latency improvement. Zone separates predict cluster from update RAMs (hyst/u).

## Experiment 3: Per-table clusters, no zones
**Order per table:** tag, pred, sec, fold_idx, fold_tag, hyst, u (no zone boundary)
**Result:** EPI=5295, P2=1.05, misps=12627
**Analysis:** Removing zones hurts — floorplanner can't separate predict-path from update-path RAMs, leading to worse packing. Zones are essential.

## Experiment 4: Hyst/U before folds
**Order per table:** tag, pred, sec, hyst, u, fold_idx, fold_tag, zone
**Result:** EPI=5102, P2=1.0, misps=12627
**Analysis:** Moving hyst/u into predict cluster before folds places folds at u_ram location — too far from tag_ram. Fold→tag wiring penalty outweighs hyst/u read locality gain.

## Experiment 5: Folds between pred and sec
**Order per table:** tag, pred, fold_idx, fold_tag, sec, zone, hyst, u
**Result:** EPI=4924, P2=0.94, misps=12627
**Analysis:** Folds at pred_ram location. Essentially same as exp 2. Fold location doesn't matter much as long as it's in the predict cluster.

## Experiment 6: Folds right after tag ★ CURRENT BEST
**Order per table:** tag, fold_idx, fold_tag, pred, sec, zone, hyst, u
**Result:** EPI=4888, P2=0.94, misps=12627
**Analysis:** Folds at tag_ram location — optimal since fold→tag comparison is tightest constraint. Slight EPI improvement over exp 2 (4888 vs 4917, -0.6%).

## Experiment 7: FB before all tables
**Order:** fb, then per-table clusters
**Result:** EPI=6210, P2=1.12, misps=12344
**Analysis:** Large fb_ctr (8K entries) pushes all tables away from each other. Much worse. FB should stay after tables.

## Experiment 8: All RAMs in cluster, zone at end
**Order per table:** tag, fold_idx, fold_tag, pred, sec, hyst, u, zone
**Result:** EPI=5077, P2=1.0, misps=12627
**Analysis:** Zone at end of cluster (after u_ram) doesn't provide predict/update separation benefit. Worse than zone between sec and hyst.

## Summary table

| Exp | Layout | EPI | P2 | vs S3 EPI |
|-----|--------|-----|-----|-----------|
| S3  | TATable (template) | 5523 | 1.023 | baseline |
| 1   | Type-grouped | 7399 | 1.04 | +34% |
| 2   | Per-table (TATable order) | 4917 | 0.94 | -11% |
| 3   | Per-table, no zones | 5295 | 1.05 | -4% |
| 4   | Hyst/u before folds | 5102 | 1.0 | -8% |
| 5   | Folds after pred | 4924 | 0.94 | -11% |
| **6** | **Folds after tag** | **4888** | **0.94** | **-11.5%** |
| 7   | FB before tables | 6210 | 1.12 | +12% |
| 8   | Zone at end | 5077 | 1.0 | -8% |

## Key insights
1. **Fold placement is critical.** Folds must be in the predict cluster, near tag/pred/sec RAMs. Best: immediately after tag_ram.
2. **Zones matter.** Separating predict-path (tag/pred/sec) from update-path (hyst/u) with a zone boundary improves floorplanner packing.
3. **FB is too large to move.** The 8K-entry fb_ctr dominates floor area; moving it disrupts other clusters.
4. **Within-cluster order matters less.** As long as folds are in the predict cluster, exact RAM ordering is ~1% EPI difference.
