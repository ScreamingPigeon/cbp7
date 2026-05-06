# TODO: Low-Hanging Optimizations (Tier 1)

## 1. History length series (quadratic→geometric hybrid)
- Source: RUNLTS S3.1.2, Ros Deep Dive S2.1
- Current: pure geometric_hist<14>(8, 200)
- Change: arithmetic progression for short histories (denser coverage), then geometric for long
- Expected: 0.5-1% MPKI reduction
- Effort: change HIST_LEN array only, no structural change

## 2. Unshuffle bit permutation on fold_idx
- Source: Alpha EV8 S7.1
- Current: idx = fold_idx ^ (PC >> 2), simple XOR
- Change: apply bit permutation (bit i → bit i^f) before XOR with PC
- Expected: reduced aliasing, near-zero logic cost
- Effort: rewiring only

## 3. Ctr-based allocation filter
- Source: CBP2016 S2.1
- Current: allocate when u==0 only
- Change: also require hyst is weak (|2*ctr+1| < threshold) before evicting
- Expected: 0.2-0.5% MPKI reduction
- Effort: add one condition to alloc candidate mask

## 4. H2P allocation Bloom filter
- Source: Bullseye (Behrendt 2025)
- Current: no per-PC allocation suppression
- Change: ~256-bit Bloom filter tracking chronic mispredictors; skip alloc for flagged PCs
- Expected: 0.5-1% MPKI reduction (prevents thrashing)
- Effort: ~32 bytes storage, minimal logic

## 5. Larger bimodal fallback
- Source: RUNLTS S3.1.1
- Current: FB_CAPACITY=8192
- Change: increase to 16K-32K entries
- Expected: reduced FB aliasing, especially for large-footprint workloads
- Effort: constant change, costs storage budget

## 6. Multi-entry allocation (MAX_ALLOC=2)
- Source: TAGE-SC S2.2, New Case for TAGE S3.2.1
- Current: MAX_ALLOC=1, one entry per mispredict
- Change: allocate 2 entries on non-consecutive tables
- Expected: faster warmup after phase changes
- Effort: change one_hot→two_hot in alloc mask, add non-consecutive constraint

## 7. Altpred update on correct (when provider weak)
- Source: L-TAGE S1.2, TAGE-SC S2.2
- Current: only provider updated; alt never trained
- Change: when provider is weak AND correct, also strengthen alt entry
- Expected: keeps alt entries warm as better backup predictions
- Effort: track alt table index through pipeline, add conditional write

---

# Bugs to Fix

## B1. Overly conservative pred_ram update (HC_IR line 1058)
- `do_pred_update = t_m1[I] & t_phw & table_wrong`
- Only updates provider when hyst is weak AND wrong
- Standard TAGE updates provider on any misprediction regardless of hyst strength
- Fix: remove `t_phw` gate (or make it configurable)

## B2. Missing `.fo1()` on tag_hit pipeline shift (HC_IR line 514)
- `current_tag_hit[I] = prefetch_tag_hit[I]` — no `.fo1()`
- Other shifts (current_idx, current_ctag) all use `.fo1()`
- Risk: unintended fanout propagation from prefetch stage

## B3. FB bimodal hyst one-way ratchet (HC_IR ~line 960, Template line 2249)
- `bh_gate = do_train & t_m1[NT] & mispredict & bh_changed`
- Hyst only updated on mispredict; can go strong→weak but never weak→strong
- Once weak, FB prediction bit is permanently flippable
- Fix: also update hyst on correct predictions (remove mispredict gate, or
  read fb_bh on all cycles not just mispredict)

## B4. fb_gate fanout over-declaration (Template line 2237)
- `fb_gate.fanout(hard<5>{})` assumes FB_BANKS=2
- For FB_BANKS=1 only 1 execute_if needed
- Harmless but misleading for audits
