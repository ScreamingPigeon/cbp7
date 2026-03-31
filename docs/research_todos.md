# Research TODOs — Branch Prediction Optimizations

Extracted from comprehensive paper survey (2025-03-24).
Papers: L-TAGE, CBP2016 TAGE-SC-L, Ros Deep Dive, RUNLTS, Alpha EV8,
Ahead Prediction, Alt-Fetch, Omnipredictor, BranchNet, Branch Precomputation, Virt BTB.

## Priority 1 — High Impact

### Statistical Corrector (SC) [Phase 3-4]
- Perceptron-based, 6-bit signed counters. Prediction = sign of sum of centered predictions.
- Components: 3 bias tables, 2 global GEHL, 2 backward-branch GEHL, IMLI, 3 local history sets.
- Dynamic threshold: PC-indexed table (32-64 entries, 8-bit) + global 12-bit threshold.
- Dynamic multiplicative factors: each GEHL component weighted 1 or 2 via monitoring counters.
- TAGE-SC selection: SC preferred by default. When SC low-confidence AND TAGE high-confidence, use monitoring counter.
- **Expected: 6-8% misprediction reduction.** [CBP2016, Deep Dive]

### History Length Series — Quadratic + Super-Exponential
- Replace geometric series with: h_n = h_{n-1} + d*(n-1) + k for n < t; h_n = round(h_{n-1} * (f*(n-t+1) + m)) for n >= t.
- Ros: h1=2, d=2, k=1, f=0.1, m=1.1, t=15. Eliminates need for set-associative tables.
- RUNLTS: 0, 6, 14, 24, 36, 50, 66, 84, 104, 126, 150, 178, 212, 252, 300, 358, 426, 506, 602, 776, 1078, 1606, 2554, 4316.
- **Expected: slight MPKI reduction, simplifies design (direct-mapped only).** [Deep Dive, RUNLTS]

### Multi-Entry Allocation (NNN=2)
- Allocate 2 entries per misprediction instead of 1. More aggressive table filling.
- Current: MAX_ALLOC=1 in DefaultAllocConfig.
- **Expected: noticeable improvement at larger predictor sizes.** [CBP2016, RUNLTS]

### u-bit Confidence Gating
- Only replace entries with u=0 AND |2*ctr+1| < 5.
- If u=0 but counter is high-confidence, decrease |2*ctr+1| instead of evicting.
- Nearly closes gap between 1-bit and 2-bit u.
- **Expected: better entry retention, fewer useful evictions.** [CBP2016]

### Skip TAGE Allocation on Loop Hits
- When loop predictor is selected as provider and prediction is correct, bypass TAGE allocation.
- Avoids wasting TAGE entries on loop-predictable branches.
- **Easy to implement in TageImpl.** [Deep Dive]

## Priority 2 — Medium Impact

### Alternate Prediction Improvements
- Extend alt selection: consider ctr=0 OR ctr=1 (not just ctr=0).
- 256 choosers indexed by: hit_bank_chunk(3b) | alt_confidence(2b) | alt_bank_hit(1b) | hit_conf(2b).
- Track second alternate prediction to update more confidence counters.
- **Expected: 0.3-0.5% improvement.** [Deep Dive]

### IMLI (Inner-Most Loop Iteration)
- Track iteration count of innermost loop. 1 table in SC.
- BrIMLI + TaIMLI from Seznec cookbook: 256 entries with 8-bit usefulness.
- **Expected: 0.2-1% MPKI reduction, >1 MPKI on some traces.** [CBP2016, RUNLTS]

### Loop Predictor: Once Selected = Final
- Don't let SC override loop prediction. Loop is highest-confidence component.
- Current: loop fights with TAGE confidence gate. Should be final when confident.
- **Expected: fewer incorrect overrides of correct loop predictions.** [Deep Dive]

### Larger Bimodal for Web/JS
- RUNLTS uses 128K-entry (20 KiB) bimodal — 5x larger than baseline.
- Handles JavaScript/web application footprints where many branches are bimodal-predictable.
- Current: 4096 entries.
- **Expected: workload-dependent, helps web/JS traces significantly.** [RUNLTS]

### Set u-bit of Alternate Bank
- When alternate correctly predicts and hit bank does not, set the alternate's u-bit.
- Protects useful alternate entries from eviction.
- **Expected: marginal improvement.** [Deep Dive]

### Thrashing Detection (RUNLTS)
- Repurpose rare counter combinations (ctr=0/-1 with u=1) as "newly allocated" markers.
- Track ratio of successful predictions to evictions. Dynamically throttle allocation rate.
- Zero extra storage.
- **Expected: helps at larger predictor sizes where thrashing is an issue.** [RUNLTS]

## Priority 3 — Latency / Architecture

### Ahead Prediction — Secondary Tags
- Add secondary tag field to each TAGE entry encoding missing history pattern.
- 1 bit of secondary tag = 2.2% performance. 5 bits = near-full benefit.
- Ahead distance of 5 covers 3-cycle TAGE latency 91.3% of the time.
- Energy scales linearly with ahead distance (vs exponentially for prior work).
- **Expected: 4.4% IPC improvement. 19.65KB additional storage.** [Ahead Prediction]

### Single-Cycle Override for Ahead
- Keep 2-bit counter in single-cycle BTB entry. 3-bit usefulness counter.
- Override ahead prediction when counter > 2.
- **Expected: 1% additional IPC on top of ahead prediction.** [Ahead Prediction]

### Bank Interleaving (EV8-style)
- 4-way bank interleaving for single-ported SRAM. Bank = PC XOR previous bank.
- Guarantees no conflict between consecutive fetch blocks.
- Reduces area 3.3x, energy 2x vs multi-ported.
- **Expected: significant area/energy savings.** [Alpha EV8, A New Case for TAGE]

### Eliminate Second Read at Retire
- On correct predictions, skip the re-read for update. Checkpoint provider + counter state at prediction.
- Combined with silent update elimination: 1.13 accesses per retired branch.
- **Expected: ~50% fewer predictor accesses. 2 MPPKI accuracy cost.** [A New Case for TAGE]

## Priority 4 — Experimental / Novel

### Register-Value-Aware SC (RUNLTS)
- 12-bit digest per register: leading_zeros(5b) + trailing_zeros(3b) + value[5:0](6b).
- 8 banks covering logical registers. Usefulness + prediction weight tables.
- Digests invalidated after 256 subsequent decoded instructions.
- Largest single-feature contribution in RUNLTS.
- **Expected: significant MPKI reduction. Requires register state visibility.** [RUNLTS]

### Call-Stack History SC Component
- 47x8 history entries with 3-bit pointer. 5 GEHL tables.
- Captures call-stack-correlated branch behavior.
- **Expected: helps polymorphic/OOP workloads.** [RUNLTS]

### Mini-BranchNet (CNN-based)
- Convolutions replaced with lookup tables, 3-4 bit quantized FC layers.
- 32KB iso-latency with 64KB TAGE-SC-L. Only 8-25 hard-to-predict static branches need CNN.
- **Expected: 9.6% MPKI reduction. 4-cycle latency.** [BranchNet]

### Branch Precomputation (TEA Thread)
- Lightweight thread that speculatively executes dependence chains of hard-to-predict branches.
- 8-wide fetch+rename, shares backend with main thread.
- **Expected: 10.1% IPC improvement. 3.5% area overhead.** [Branch Precomputation]

### Omnipredictor — Unified TAGE for branches + loads + indirects
- Reinterpret 3-bit TAGE counter based on instruction type.
- Branches: direction. Loads: store distance. Indirects: BTB pointer.
- All share same physical TAGE tables.
- **Expected: "free" MDP and indirect prediction using existing TAGE.** [Omnipredictor]

## Current Bugs / Issues

### Loop Predictor rwram Bank Conflicts [Task #1]
- rwram writes rarely land (28784/236K = 12%). Same-set read+write always hits same bank.
- Confidence counter stays at 0, overrides never fire.
- Fix: cache last-read entry in reg, only hit RAM on set change. Write to both reg and RAM.
- File: predictors/custom/LoopPredictor.hpp line 57.

### P2 Latency Gap (2.087 vs reference 1.86)
- ceil(2.087)=3 vs ceil(1.86)=2. Huge VFS impact.
- Overrider confidence gate + mux chain adds ~20ps (now using P1 regs, but gate still in P2).
- Extra fanout from OVR=1 on match1, pred1, pred2, newly_alloc.
- Options: reduce tag width (11->9), fewer tables (8->6), disable meta.


# Prakhar's TODOs:

  CBP2016 TAGE-SC-L (#8)                                                                                                     
  - Multi-GEHL SC — Your SC is bias-only. Adding 2-4 GEHL tables with geometric history lengths is the expected biggest SC
  improvement (6-8% mispred reduction). Already Priority 1 in your todos.                                                    
  - TAGE-SC chooser — 3 confidence cases with 2 monitoring counters. 0.7% reduction. Directly addresses your SC
  cross-workload regression problem.                                                                                         
  - Confidence-gated allocation (u-bit trick) — Only evict u=0 AND low-confidence entries. Nearly closes 1-bit vs 2-bit u    
  gap. Zero extra storage. Already in your research_todos.                                                               
  - Non-consecutive allocation — Allocated entries must be in non-consecutive tables. Reduces destructive interference.      
  Simple check in allocation loop.                                                                                     
                                                                                                                             
  Ros Deep Dive (#14)
  - Quadratic-to-super-exponential history series — Already planned (Priority 1). Eliminates need for set-associative tables,
   denser medium-history coverage.                                                                                           
  - Skip TAGE allocation on loop hits — Easy, already planned in overrider coordination.                                     
  - Improved alt chooser — 256 counters, considers confidence < 2 (not just 0). Medium effort, moderate gain.                
  - Improved TAGE-SC chooser — 16 counters indexed by TAGE confidence + SC sum bin. Better than CBP2016's simpler version.   
  - Set u-bit of alternate bank — Trivial change, minor gain.                                                                
  - Loop predictor is final — Don't let SC override loop. Already planned.                                                   
                                                                                                                             
  A New Case for TAGE (#6)                                                                                                   
  - TAGE counter value as SC confidence — Add (2*ctr+1) * 8 to SC sum. Makes SC less likely to override high-confidence TAGE 
  predictions. Directly fixes your SC over-correction problem.                                                               
  - Simplified 1-bit u with global reset — You already use this pattern (DECAY_CTR). Validates your approach.                
  - Multi-entry allocation (NNN=2-4) — Already planned Priority 1. Faster warming for larger predictors.                     
  - Eliminating silent updates — Skip writes when old==new value. >90% write elimination. Would significantly reduce EPI.    
  - Local Statistical Corrector (LSC) — 5 LGEHL tables with local history. >8% mispred reduction on top of TAGE. Different   
  from global SC — targets the ~10 static branches causing most dynamic mispredictions.                                      
                                                                                                                             
  TAGE-SC original (#15)                                                                                                     
  - GEHL-based SC (global history) — Same as CBP2016 but with implementation details. 5-6 bit counters (you already use      
  6-bit). Prediction = sign(sum).                                                                                            
  - Bias indexed by TAGE output — Removing TAGE direction from bias index increases mispred by 1.6%. Your current PC-only    
  bias is missing this. But you found direction-indexed bias adds 2 cycles to P2 — the precomputed-sum trick (sum without    
  TAGE vote, add TAGE vote at the end) should solve this.                                                                    
  - Dynamic threshold (PC-indexed) — 32-entry table + global counter. 0.1% benefit. Low priority but addresses your threshold
   runaway issue.                                                                                                            
                                                                                                                             
  Worth Investigating (medium effort, clear benefit)
                                                                                                                             
  Samsung Exynos (#9)                                                                                                        
  - Scaled Hashed Perceptron (SHP) — This is the practical perceptron design Jimenez pointed to. 8-16 tables with geometric  
  history, sign/magnitude weights, XOR hash. Your SC is already a lightweight version of this. The key insight: BTB BIAS     
  weight — a per-branch dedicated weight stored in the BTB entry, doubled and added to sum. Eliminates the worst aliasing
  problem. You could add a dedicated per-PC bias as a 4th bias table with wider indexing.                                    
  - Always-taken branch filtering — Skip SC updates for always-taken branches. Reduces aliasing/pollution. Easy to add.
  - Stochastic search for history lengths — Better than hand-tuned. Your coordinate descent sweep already does this    
  partially.                                                                                                                 
                                                                                                                             
  L-TAGE (#10)                                                                                                               
  - USE_ALT_ON_NA — You already have this (meta predictor). 4-bit counter, use alt when provider newly allocated + counter   
  negative.                                                                                                                  
  - Graceful u-bit reset (alternating bits) — Alternate clearing u0 then u1 every 512K branches. You have probabilistic decay
   instead, which may be better.                                                                                             
  - Tag computation overlapped with table read — Tag match spans during index computation. Zero added latency for complex tag
   hashes. Relevant if you want wider tags without P2 cost.                                                                  
                                                                                                                             
  RUNLTS (#13)    
  - IMLI components (BrIMLI + TaIMLI) — Already in your Priority 2 todos. Reduced MPKI by ≥1 on three traces. 20K bits       
  storage. Good complement to loop predictor.                                                                                
  - Thrashing detection + dynamic allocation throttling — Zero extra storage (repurposes rare counter/u-bit combos). Already 
  in your todos.                                                                                                             
  - Large bimodal (128K entries) — Your sweep already explored bimodal_size. Current 8192 may be optimal for your VFS budget.
  - History length series — Same as Ros, validates the approach.                                                             
                                                                                                                             
  SONet (#12)                                                                                                                
  - H2P identification + TAGE space reclamation — When overrider handles a branch, skip TAGE allocation. You already plan    
  this. Validates the approach.                                                                                              
  - Misprediction concentration — Top 8 branches cause ~65% of mispredictions. Suggests your SC should focus resources on few
   branches rather than being general.                                                                                       
  - Per-branch dedicated weights — Each H2P branch gets its own weight set (8 networks in CAM). Too complex for HARCOM, but  
  the principle of per-branch specialization could inform a wider-indexed SC bias table.                                   
                                                                                              
