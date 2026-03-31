  ---             
  1. 738300a118.pdf — BranchNet (Zangeneh et al., MICRO 2020)                                                                
                                                                                                                             
  Core thesis: TAGE is limited by noisy global history; CNNs can learn to ignore noise for hard branches.                    
                                                                                                                             
  ┌───────────────────┬──────────────────────────────────┬────────────────────────────────┬─────────────────────────────┐    
  │   Optimization    │           How it works           │            Benefit             │         Limitations         │    
  ├───────────────────┼──────────────────────────────────┼────────────────────────────────┼─────────────────────────────┤    
  │ Hybrid TAGE + CNN │ Profile top hard-to-predict      │ Mini-BranchNet (32KB): 9.6%    │ Requires offline training   │
  │  for hard         │ static branches offline, train   │ MPKI reduction over 64KB       │ (6-18h on 4 GPUs), ISA/OS   │    
  │ branches          │ per-branch CNN models            │ TAGE-SC-L                      │ support to load models      │    
  ├───────────────────┼──────────────────────────────────┼────────────────────────────────┼─────────────────────────────┤
  │ Geometric history │ 5 slices at different history    │ Captures both short- and       │ Longer histories are        │    
  │  slices for       │ lengths (37, 77, 152, 302, 603)  │ long-range correlations        │ noisier                     │    
  │ features          │                                  │                                │                             │
  ├───────────────────┼──────────────────────────────────┼────────────────────────────────┼─────────────────────────────┤    
  │ Sum-pooling to    │ Aggregate features over windows; │ Makes predictor resilient to   │ Discards fine-grained       │
  │ compress noisy    │  wider pooling for longer        │ history shifts; dramatically   │ positional info             │    
  │ history           │ histories                        │ reduces storage                │                             │
  ├───────────────────┼──────────────────────────────────┼────────────────────────────────┼─────────────────────────────┤    
  │ Table lookups     │ Pre-compute all convolution      │ Eliminates runtime arithmetic  │ Requires binarized          │
  │ replacing         │ outputs, store in 256-entry      │ entirely                       │ activations, reduces        │    
  │ convolutions      │ tables                           │                                │ expressiveness              │
  ├───────────────────┼──────────────────────────────────┼────────────────────────────────┼─────────────────────────────┤    
  │                   │                                  │ Minimal accuracy loss;         │                             │
  │ Quantization (3-4 │ Fixed-point arithmetic with very │ validates that branch          │ Some accuracy loss          │    
  │  bit weights)     │  low precision for FC layers     │ prediction doesn't need        │                             │
  │                   │                                  │ high-precision weights         │                             │    
  ├───────────────────┼──────────────────────────────────┼────────────────────────────────┼─────────────────────────────┤
  │ IMLI counters for │ Track innermost loop iteration   │ Helps predict branches         │ Not implemented in this     │
  │  nested loops     │ count as additional predictor    │ correlated with outer loop     │ paper (future work)         │    
  │                   │ input                            │ iterations                     │                             │
  ├───────────────────┼──────────────────────────────────┼────────────────────────────────┼─────────────────────────────┤    
  │ Bias-free history │ Remove always-taken/biased       │ Improves effective history     │ Paper argues offline        │
  │  filtering        │ branches from global history     │ quality                        │ learning is more powerful   │    
  │                   │ before use                       │                                │                             │
  └───────────────────┴──────────────────────────────────┴────────────────────────────────┴─────────────────────────────┘    
                  
  ---                                                                                                                        
  2. aclp p392.pdf
                                                                                                                             
  Not relevant — Conference proceedings (APPT 2025), not a branch prediction paper.
                                                                                                                             
  ---             
  3. ahead-prediction.pdf — Enabling Ahead Prediction (Cai, Deshmukh, Patt, ISCA '25)                                        
                                                                                                                             
  ┌───────────────────────┬───────────────────────────────────┬────────────────────────────┬────────────────────────────┐
  │     Optimization      │           How it works            │          Benefit           │        Limitations         │    
  ├───────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────────┤    
  │ Secondary tag for     │ Additional tag field identifying  │ Energy 1.5x (vs 14.6x      │ Entries with same ahead    │    
  │ ahead TAGE entries    │ which missing history pattern the │ prior work); 4.4% IPC      │ history always conflict in │    
  │                       │  counter belongs to               │ improvement                │  every table               │    
  ├───────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────────┤
  │ Missing history hash  │ Hash of skipped branches'         │ Decouples tag width from   │                            │    
  │ function              │ predicted targets via XOR+rotate, │ ahead distance; handles    │ —                          │    
  │                       │  independent of ahead distance    │ indirect branches          │                            │
  ├───────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────────┤    
  │ Single-cycle override │ 3-bit confidence counter; bimodal │                            │ Only helps branches simple │
  │  of ahead predictor   │  overrides ahead when confident   │ 1% IPC improvement         │  enough for single-cycle   │    
  │                       │                                   │                            │ predictor                  │
  ├───────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────────┤    
  │ Late prediction       │ Use single-cycle prediction       │ 2.4% IPC vs always         │ Early flushes when         │
  │ handling              │ initially; flush if ahead         │ stalling                   │ predictions disagree       │    
  │ (non-stalling)        │ prediction disagrees              │                            │                            │
  ├───────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────────┤    
  │ Bank interleaving to  │ Reduce from 21 to 14 tables read  │ No significant MPKI hurt   │ Below ~15 tables starts to │
  │ reduce tables read    │ per prediction                    │                            │  noticeably increase MPKI  │    
  ├───────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────────┤
  │                       │ <3 missing history patterns 97%   │ Enables entire ahead       │ gcc and high-MPKI          │    
  │ Key empirical finding │ of the time                       │ design — only need a       │ benchmarks exhibit more    │    
  │                       │                                   │ handful of patterns        │ patterns                   │
  └───────────────────────┴───────────────────────────────────┴────────────────────────────┴────────────────────────────┘    
                  
  ---
  4. Alpha-Ev8-seznec.pdf — Alpha EV8 Predictor (Seznec et al.)
                                                                                                                             
  ┌───────────────────────────┬────────────────────────────────────┬──────────────────────────────┬─────────────────────┐ 
  │       Optimization        │            How it works            │           Benefit            │     Limitations     │    
  ├───────────────────────────┼────────────────────────────────────┼──────────────────────────────┼─────────────────────┤ 
  │                           │ Don't strengthen counters when all │ Better accuracy + space      │ More complex update │    
  │ Partial update policy     │  components agree; limit writes on │ utilization than total       │  logic              │ 
  │                           │  misprediction                     │ update                       │                     │    
  ├───────────────────────────┼────────────────────────────────────┼──────────────────────────────┼─────────────────────┤ 
  │ Separate                  │ MSB (prediction) and LSB           │ Correct predictions: only 1  │ Hysteresis array    │    
  │ prediction/hysteresis     │ (hysteresis) in separate arrays    │ read + at most 1 write       │ suffers more        │    
  │ arrays                    │                                    │                              │ aliasing            │ 
  ├───────────────────────────┼────────────────────────────────────┼──────────────────────────────┼─────────────────────┤    
  │ Shared hysteresis bits    │ Two prediction entries share one   │ Saves significant storage    │ Rare interference   │
  │                           │ hysteresis entry                   │ (EV8: 144 Kbits saved)       │ possible            │    
  ├───────────────────────────┼────────────────────────────────────┼──────────────────────────────┼─────────────────────┤
  │                           │ Compress up to 8 branch outcomes   │ Same accuracy as full        │ Loses per-branch    │    
  │ Block compressed history  │ per fetch block into 1 history bit │ history; represents ~3x more │ granularity within  │    
  │ (lghist)                  │  via XOR with PC bit               │  branches in same register   │ fetch block         │
  │                           │                                    │ width                        │                     │    
  ├───────────────────────────┼────────────────────────────────────┼──────────────────────────────┼─────────────────────┤
  │                           │                                    │                              │ Small impact        │
  │ Path info in history      │ XOR fetch block addresses into     │ De-aliases otherwise aliased │ without index       │    
  │                           │ index functions                    │  history paths               │ function            │
  │                           │                                    │                              │ constraints         │    
  ├───────────────────────────┼────────────────────────────────────┼──────────────────────────────┼─────────────────────┤
  │ Conflict-free bank        │ 4-way banking guaranteeing         │ Eliminates multi-ported      │ Restricts 2 bits of │    
  │ interleaving              │ consecutive fetch blocks access    │ memory; saves area/power     │  index function     │
  │                           │ different banks                    │                              │                     │    
  ├───────────────────────────┼────────────────────────────────────┼──────────────────────────────┼─────────────────────┤
  │                           │ Use history lengths >>             │ "Almost always beneficial"   │ Diminishing returns │    
  │ Very long history lengths │ log2(table_entries)                │ for large predictors         │  for smaller        │
  │                           │                                    │                              │ predictors          │    
  ├───────────────────────────┼────────────────────────────────────┼──────────────────────────────┼─────────────────────┤
  │ Carefully engineered      │ Minimize aliasing within and       │ 352 Kbit tuned predictor     │ Manual tuning       │    
  │ index/hash functions      │ across tables; favor history bits  │ matches 512 Kbit             │ required            │
  │                           │ for column index                   │ unconstrained                │                     │    
  └───────────────────────────┴────────────────────────────────────┴──────────────────────────────┴─────────────────────┘

  ---
  5. alt-fetch.pdf — Alternate Path Fetch (Deshmukh, Cai, Patt, ISCA 2024)
                                                                                                                             
  ┌─────────────────────────┬──────────────────────────────────────┬───────────────────────────┬────────────────────────┐ 
  │      Optimization       │             How it works             │          Benefit          │      Limitations       │    
  ├─────────────────────────┼──────────────────────────────────────┼───────────────────────────┼────────────────────────┤ 
  │                         │ 4 mini-TAGE-SC-Ls as banks, PC       │ 2 banks reduces aliasing; │ 4+ banks hurts MPKI by │    
  │ Banking TAGE-SC-L       │ determines bank                      │  enables parallel access  │  ~0.1 due to capacity  │ 
  │                         │                                      │                           │ issues                 │    
  ├─────────────────────────┼──────────────────────────────────────┼───────────────────────────┼────────────────────────┤ 
  │ H2P table               │ 128-entry 3-bit saturating counter   │ 95.4% misprediction       │ 89.6% wastage (many    │    
  │ (hard-to-predict        │ table, periodic decrement every 20K  │ coverage                  │ marked H2P are         │ 
  │ identification)         │ instructions                         │                           │ actually correct)      │    
  ├─────────────────────────┼──────────────────────────────────────┼───────────────────────────┼────────────────────────┤ 
  │ TAGE confidence for H2P │ Unsaturated counter = low            │ 56.3% coverage, 74.5%     │ Lower coverage than    │ 
  │  selection              │ confidence; path-sensitive           │ wastage; complementary to │ H2P table              │    
  │                         │ (per-history, not just per-PC)       │  H2P table                │                        │ 
  ├─────────────────────────┼──────────────────────────────────────┼───────────────────────────┼────────────────────────┤    
  │ APF (Alternate Path     │ Shallow parallel pipeline fetching   │ 5% geomean speedup; 2%    │ Only helps high-MPKI   │ 
  │ Fetch)                  │ the alternate path of H2P branches   │ core area vs 20% for true │ workloads              │    
  │                         │                                      │  16-wide                  │                        │ 
  ├─────────────────────────┼──────────────────────────────────────┼───────────────────────────┼────────────────────────┤    
  │ TAGE confidence         │ Prioritize low-TAGE-confidence       │ Lower wastage than        │                        │
  │ priority in scheduling  │ branches over H2P-only branches for  │ H2P-only                  │ —                      │    
  │                         │ APF                                  │                           │                        │
  └─────────────────────────┴──────────────────────────────────────┴───────────────────────────┴────────────────────────┘    
                  
  ---
  6. A new case for TAGE.pdf — A New Case for TAGE (Seznec)
                                                                                                                             
  ┌───────────────────────┬───────────────────────────────────────┬────────────────────────────┬────────────────────────┐ 
  │     Optimization      │             How it works              │          Benefit           │      Limitations       │    
  ├───────────────────────┼───────────────────────────────────────┼────────────────────────────┼────────────────────────┤ 
  │ Eliminating silent    │ Compare old/new values before         │ Eliminates >90% of writes; │ Requires reading old   │    
  │ updates               │ writing; skip if identical            │  avg 2.17 effective writes │ value                  │ 
  │                       │                                       │  per misprediction         │                        │    
  ├───────────────────────┼───────────────────────────────────────┼────────────────────────────┼────────────────────────┤ 
  │                       │ Scenario [C]: re-read only on         │ ~2.6% accuracy loss vs     │                        │    
  │ Eliminating           │ mispredictions, not correct           │ oracle; nearly halves      │ Scenario [B] (never    │ 
  │ retire-time re-read   │ predictions                           │ energy for correct         │ re-read) is too lossy  │    
  │                       │                                       │ predictions                │                        │ 
  ├───────────────────────┼───────────────────────────────────────┼────────────────────────────┼────────────────────────┤ 
  │ 4-way                 │ Split each table into 4 banks;        │ 3.3x area reduction, ~2x   │ Marginal accuracy loss │    
  │ bank-interleaved      │ algorithm guarantees no conflicts for │ energy reduction vs 3-port │  (627 vs 625 MPPKI)    │ 
  │ single-port tables    │  3-cycle intervals                    │  arrays                    │                        │    
  ├───────────────────────┼───────────────────────────────────────┼────────────────────────────┼────────────────────────┤ 
  │                       │ Fully-associative table recording     │ Recovers 3/4 of            │ Gap varies             │    
  │ Immediate Update      │ inflight branch outcomes; subsequent  │ delayed-update             │ significantly by       │ 
  │ Mimicker (IUM)        │ branches hitting same TAGE entry get  │ mispredictions             │ benchmark              │    
  │                       │ executed outcome                      │                            │                        │ 
  ├───────────────────────┼───────────────────────────────────────┼────────────────────────────┼────────────────────────┤ 
  │ Multi-entry           │ Allocate up to 4 entries on           │ Faster warming for large   │ Only beneficial for    │
  │ allocation            │ non-consecutive tables per            │ predictors                 │ large predictors       │    
  │                       │ misprediction                         │                            │ (256K-512Kbit)         │
  ├───────────────────────┼───────────────────────────────────────┼────────────────────────────┼────────────────────────┤    
  │ Simplified 1-bit      │ Set u when provider correct AND       │ "More efficient than 2-bit │                        │
  │ u-bit with global     │ altpred incorrect; 8-bit counter      │  useful counters"          │ —                      │    
  │ reset                 │ triggers global reset                 │                            │                        │
  ├───────────────────────┼───────────────────────────────────────┼────────────────────────────┼────────────────────────┤    
  │ Local Statistical     │ 5 LGEHL tables indexed by local       │ 30 Kbit LSC reduces        │ Only ~10 static        │
  │ Corrector (LSC)       │ history (32-entry direct-mapped local │ misprediction >8% on top   │ branches cause most    │    
  │                       │  history table)                       │ of TAGE+IUM                │ dynamic mispredictions │
  ├───────────────────────┼───────────────────────────────────────┼────────────────────────────┼────────────────────────┤    
  │ TAGE counter value as │ Centered counter (2*ctr+1) * 8 added  │ SC less likely to override │                        │
  │  SC confidence        │ to SC sum                             │  high-confidence TAGE      │ —                      │    
  │                       │                                       │ predictions                │                        │
  ├───────────────────────┼───────────────────────────────────────┼────────────────────────────┼────────────────────────┤    
  │ Dynamic threshold for │ Runtime-adjusted threshold for SC     │ "Ensures SC use is         │ Marginal benefit       │
  │  SC                   │ override decisions                    │ beneficial"                │ (0.1%)                 │    
  ├───────────────────────┼───────────────────────────────────────┼────────────────────────────┼────────────────────────┤
  │ TAGE-LSC with         │ Combined design: no re-read on        │ 569 MPPKI at 512Kbits —    │                        │    
  │ single-port           │ correct, single-port,                 │ competitive with SOTA, far │ —                      │    
  │ bank-interleaved      │ bank-interleaved                      │  more HW-practical         │                        │
  └───────────────────────┴───────────────────────────────────────┴────────────────────────────┴────────────────────────┘    
                  
  ---
  7. branch-precomputation.pdf — TEA (Deshmukh, Cai, Patt, MICRO 2024)
                                                                                                                             
  ┌──────────────────────┬─────────────────────────────────┬─────────────────────────────────┬───────────────────────────┐
  │     Optimization     │          How it works           │             Benefit             │        Limitations        │   
  ├──────────────────────┼─────────────────────────────────┼─────────────────────────────────┼───────────────────────────┤
  │                      │ 256-entry 3-bit counter table,  │ Identifies frequently           │ Needs tuning of decay     │   
  │ H2P counter table    │ periodic decrement every 50K    │ mispredicting branches          │ period                    │
  │                      │ instructions                    │                                 │                           │   
  ├──────────────────────┼─────────────────────────────────┼─────────────────────────────────┼───────────────────────────┤
  │ Backward dataflow    │ Trace dependence chain backward │ Lightweight yet accurate        │ Limited by buffer size    │   
  │ walk                 │  through retired instruction    │ chains; ~500 cycle walk         │ (512 entries)             │
  │                      │ buffer                          │                                 │                           │   
  ├──────────────────────┼─────────────────────────────────┼─────────────────────────────────┼───────────────────────────┤   
  │ Iterative chain      │ Previously identified chain     │ Chains grow to span thousands   │ Requires multiple dynamic │
  │ extension            │ instructions initiate future    │ of instructions                 │  occurrences              │   
  │                      │ walks                           │                                 │                           │
  ├──────────────────────┼─────────────────────────────────┼─────────────────────────────────┼───────────────────────────┤   
  │ Bit-mask combining   │ Per-basic-block bit-mask OR'd   │ 99.3% precomputation accuracy;  │ Adds unnecessary          │
  │ across control flows │ across different paths          │ removing masks drops coverage   │ instructions from OR'ing  │   
  │                      │                                 │ from 76% to 39%                 │                           │
  ├──────────────────────┼─────────────────────────────────┼─────────────────────────────────┼───────────────────────────┤
  │ Synchronized         │ TEA thread branches inherit     │ Enables both direction and      │                           │
  │ timestamps for early │ main thread timestamps          │ target misprediction resolution │ —                         │   
  │  flushes             │                                 │                                 │                           │
  ├──────────────────────┼─────────────────────────────────┼─────────────────────────────────┼───────────────────────────┤   
  │                      │ Dedicated frontend +            │ 10.1% IPC improvement; 3.5%     │ Data-dependent branches   │
  │ Overall TEA system   │ partitioned backend for on-core │ area overhead                   │ only; 8.5% peak power     │   
  │                      │  precomputation                 │                                 │ increase                  │
  └──────────────────────┴─────────────────────────────────┴─────────────────────────────────┴───────────────────────────┘   
                  
  Key insight: Some branches are fundamentally unpredictable by history-based methods (TAGE, perceptron) — they depend on    
  data values loaded from memory. This sets a ceiling on TAGE improvements.
                                                                                                                             
  ---             
  8. CBP2016-TAGE-SC-L-again.pdf — TAGE-SC-L at CBP-5 (Seznec, 2016)
                                                                                                                             
  ┌──────────────────────────┬───────────────────────────────────┬────────────────────────────┬────────────────────────┐  
  │       Optimization       │           How it works            │          Benefit           │      Limitations       │     
  ├──────────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────┤     
  │                          │ Short/long history tables share   │ 2-3% misprediction         │                        │     
  │ Bank interleaving        │ bank groups; different tag widths │ reduction — single most    │ Adds design complexity │     
  │                          │  per group                        │ significant contributor    │                        │     
  ├──────────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────┤  
  │ Multi-GEHL Statistical   │ Multiple GEHL components indexed  │ ~6% (8KB), ~8% (64KB)      │ Complexity scales with │     
  │ Corrector                │ by different history types        │ misprediction reduction    │  budget                │  
  │                          │ (global, backward, IMLI, local)   │                            │                        │     
  ├──────────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────┤  
  │ TAGE-SC chooser          │ 3 cases based on TAGE confidence  │ 0.7% misprediction         │                        │ 
  │ (confidence-based)       │ vs SC confidence; 2 monitoring    │ reduction over             │ —                      │     
  │                          │ counters                          │ always-using-SC            │                        │ 
  ├──────────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────┤     
  │ Confidence-gated         │ Entries with u=0 but medium/high  │ Nearly fills gap between   │ Only relevant when u   │     
  │ allocation (u-bit trick) │ confidence protected from         │ 1-bit and 2-bit u counters │ constrained to 1 bit   │ 
  │                          │ replacement                       │                            │                        │     
  ├──────────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────┤ 
  │ Partial associativity    │ 2-way set-associativity only for  │ 0.6% misprediction         │ Applying to all tables │
  │ (2-way intermediate      │ intermediate history length       │ reduction                  │  is counterproductive  │     
  │ tables)                  │ tables                            │                            │                        │
  ├──────────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────┤     
  │ IMLI counter (innermost  │ 8-bit counter tracking innermost  │                            │                        │
  │ loop iteration)          │ loop iteration; indexes GEHL      │ 0.2% improvement           │ —                      │     
  │                          │ table                             │                            │                        │
  ├──────────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────┤     
  │ Global backward branch   │ Separate history tracking only    │ 0.3% improvement           │ Marginal on its own    │
  │ history                  │ backward branches                 │                            │                        │     
  ├──────────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────┤
  │ Dynamic multiplicative   │ Per-component runtime weight (1x  │ Part of overall SC tuning  │ Extra storage for      │     
  │ factors for GEHL         │ or 2x)                            │                            │ monitoring counters    │     
  ├──────────────────────────┼───────────────────────────────────┼────────────────────────────┼────────────────────────┤
  │ SC bias with extra TAGE  │ 3 bias tables indexed by PC +     │ 0.3% improvement when      │                        │     
  │ info                     │ TAGE direction + TAGE confidence  │ augmented with TAGE info   │ —                      │     
  │                          │ + hitting table number            │                            │                        │
  └──────────────────────────┴───────────────────────────────────┴────────────────────────────┴────────────────────────┘     
                  
  ---
  9. exynos_isca2020.pdf — Samsung Exynos M1-M6 (ISCA 2020)
                                                                                                                             
  ┌─────────────────────┬───────────────────────────────────┬─────────────────────────────────┬─────────────────────────┐
  │    Optimization     │           How it works            │             Benefit             │       Limitations       │    
  ├─────────────────────┼───────────────────────────────────┼─────────────────────────────────┼─────────────────────────┤
  │ Scaled Hashed       │ 8-16 tables of weights in         │ MPKI 3.62→2.54 across M1→M6     │ Aliasing between        │    
  │ Perceptron (SHP)    │ sign/magnitude; XOR hash of       │ (25.6% SPECint reduction)       │ branches sharing        │
  │                     │ GHIST+PHIST+PC                    │                                 │ entries                 │    
  ├─────────────────────┼───────────────────────────────────┼─────────────────────────────────┼─────────────────────────┤
  │ Always-taken branch │ Skip SHP updates for always-taken │ Reduces aliasing/pollution in   │ —                       │    
  │  filtering          │  branches                         │ SHP tables                      │                         │    
  ├─────────────────────┼───────────────────────────────────┼─────────────────────────────────┼─────────────────────────┤
  │ Zero-bubble         │ Graph algorithm learns branch     │ Zero-bubble prediction for      │                         │    
  │ graph-based uBTB    │ sequences; locks and drives       │ tight loops; SHP disabled for   │ Limited capacity        │    
  │                     │ pipeline at 0 bubbles             │ power savings                   │                         │
  ├─────────────────────┼───────────────────────────────────┼─────────────────────────────────┼─────────────────────────┤    
  │ ZAT/ZOT             │ Copy targets of                   │ 1 taken branch per cycle        │                         │
  │ (zero-bubble        │ always/often-taken branches into  │ throughput                      │ Increased mBTB storage  │    
  │ always/often taken) │ predecessor's BTB entry           │                                 │                         │
  ├─────────────────────┼───────────────────────────────────┼─────────────────────────────────┼─────────────────────────┤    
  │ Mispredict Recovery │ Pre-record most likely 3 fetch    │ Example: 14 instructions in 5   │ Only for pre-identified │
  │  Buffer (MRB)       │ addresses after low-confidence    │ cycles vs 9 cycles without MRB  │  low-confidence         │    
  │                     │ branch                            │                                 │ branches                │
  ├─────────────────────┼───────────────────────────────────┼─────────────────────────────────┼─────────────────────────┤    
  │ VPC indirect        │ Break indirect into sequence of   │ Leverages existing SHP          │ O(n) training; doesn't  │
  │ prediction          │ virtual conditional predictions   │ infrastructure                  │ scale to hundreds of    │
  │                     │ via SHP                           │                                 │ targets                 │    
  ├─────────────────────┼───────────────────────────────────┼─────────────────────────────────┼─────────────────────────┤
  │ BTB BIAS weight     │ Signed BIAS in BTB entry, doubled │ Every branch has at least one   │ Increases BTB entry     │    
  │ (per-branch)        │  and added to SHP sum             │ unshared weight, reducing       │ size                    │    
  │                     │                                   │ aliasing impact                 │                         │
  ├─────────────────────┼───────────────────────────────────┼─────────────────────────────────┼─────────────────────────┤    
  │ Stochastic search   │ Empirically search GHIST/PHIST    │ Better MPKI than hand-tuned     │ Requires representative │
  │ for history         │ intervals per table               │ intervals                       │  workloads              │    
  │ optimization        │                                   │                                 │                         │
  └─────────────────────┴───────────────────────────────────┴─────────────────────────────────┴─────────────────────────┘    
                  
  ---
  10. L-TAGE.pdf — L-TAGE (Seznec)
                                                                                                                             
  ┌──────────────────────┬──────────────────────────────────────┬───────────────────────┬────────────────────────────────┐
  │     Optimization     │             How it works             │        Benefit        │          Limitations           │   
  ├──────────────────────┼──────────────────────────────────────┼───────────────────────┼────────────────────────────────┤
  │                      │ 4-bit counter; use alternate         │ "More efficient on    │ Only a single counter;         │   
  │ USE_ALT_ON_NA        │ prediction when provider is newly    │ several applications" │ application-global             │
  │                      │ allocated + counter negative         │                       │                                │   
  ├──────────────────────┼──────────────────────────────────────┼───────────────────────┼────────────────────────────────┤
  │ Graceful u-bit reset │ Alternately reset u0 then u1 every   │ Prevents permanently  │ Reset period is a tuning       │   
  │  (alternating bits)  │ 512K branches; flip bit significance │ "useful" entries      │ parameter                      │   
  │                      │  after reset                         │                       │                                │
  ├──────────────────────┼──────────────────────────────────────┼───────────────────────┼────────────────────────────────┤   
  │ Probabilistic        │ Start search at table i+1 (prob      │ Avoids ping-pong      │                                │
  │ allocation starting  │ 1/2), i+2 (prob 1/4), i+3 (prob 1/4) │ allocation            │ Simple 2-bit PRNG              │   
  │ point                │                                      │                       │                                │
  ├──────────────────────┼──────────────────────────────────────┼───────────────────────┼────────────────────────────────┤   
  │ Tag width scaling    │ Narrower tags for shorter-history    │ Saves storage; more   │ Too-small tags lead to false   │
  │ with history length  │ tables (7 bits for T1,T2 up to 15    │ entries per table     │ matches                        │   
  │                      │ bits for T12)                        │                       │                                │
  ├──────────────────────┼──────────────────────────────────────┼───────────────────────┼────────────────────────────────┤   
  │ Kernel/user history  │ Separate history sets; user history  │ Prevents kernel       │ Requires kernel/user address   │
  │ separation           │ only updated on user branches        │ branch pollution      │ discrimination                 │   
  ├──────────────────────┼──────────────────────────────────────┼───────────────────────┼────────────────────────────────┤
  │ Loop predictor (256  │ Identifies constant-iteration loops; │ L-TAGE: 3.314 vs      │ "Essentially marginal"         │   
  │ entries, 4-way)      │  14-bit counters; age-based          │ 3.368 misp/KI without │ benefit; "probably not worth   │   
  │                      │ replacement                          │  loop                 │ the complexity"                │
  ├──────────────────────┼──────────────────────────────────────┼───────────────────────┼────────────────────────────────┤   
  │                      │ Read 2^(X-1) adjacent entries per    │ 8C-Ahead: only ~0.1   │ Requires reading 2^(X-1)       │
  │ Ahead pipelining     │ table; select based on intermediate  │ misp/KI loss vs       │ entries per table              │   
  │                      │ branch info                          │ non-ahead             │                                │
  ├──────────────────────┼──────────────────────────────────────┼───────────────────────┼────────────────────────────────┤   
  │ Simplified 3-entry   │ Single stage of 3-entry XOR gates    │ Only +0.03 misp/KI vs │ Slightly less effective        │
  │ XOR index hash       │ for indexing                         │  full hash functions  │ distribution                   │   
  ├──────────────────────┼──────────────────────────────────────┼───────────────────────┼────────────────────────────────┤
  │ Tag computation      │ Tag match spans during index         │ Allows complex tag    │                                │   
  │ overlapped with      │ computation and table read           │ hashes at no latency  │ —                              │   
  │ table read           │                                      │ cost                  │                                │
  └──────────────────────┴──────────────────────────────────────┴───────────────────────┴────────────────────────────────┘   
                  
  ---
  11. omnipredictor.pdf — Omnipredictor (Perais & Seznec)
                                                                                                                             
  ┌──────────────────────────┬────────────────────────────────────┬──────────────────────────┬──────────────────────────┐
  │       Optimization       │            How it works            │         Benefit          │       Limitations        │    
  ├──────────────────────────┼────────────────────────────────────┼──────────────────────────┼──────────────────────────┤
  │                          │ Same TAGE counter interpreted as   │ Near-zero storage        │ Counter stealing causes  │    
  │ Multi-purpose 3-bit      │ direction confidence, store        │ overhead for MDP +       │ marginal direction       │
  │ counter                  │ distance, or BTB pointer depending │ indirect prediction      │ accuracy loss            │    
  │                          │  on instruction type               │                          │                          │
  ├──────────────────────────┼────────────────────────────────────┼──────────────────────────┼──────────────────────────┤    
  │ Probabilistic u-bit      │ Reset u with probability 1/256     │ Faster eviction of stale │ Probability is a tuning  │    
  │ reset for transient      │ when prediction was "silently      │  entries that are never  │ parameter                │
  │ entries                  │ wrong"                             │ formally mispredicted    │                          │    
  ├──────────────────────────┼────────────────────────────────────┼──────────────────────────┼──────────────────────────┤
  │                          │ Counter selects which BTB entry to │ In 4/6 benchmarks,       │ A2 costs pipeline        │
  │ TAGE-IT-BTB (indirect    │  use as indirect target; A2        │ matches full ITTAGE at   │ bubble; can't match      │    
  │ via TAGE + BTB)          │ fallback with gshare-hash          │ zero storage overhead    │ large ITTAGE on some     │
  │                          │                                    │                          │ benchmarks               │    
  ├──────────────────────────┼────────────────────────────────────┼──────────────────────────┼──────────────────────────┤
  │ Block-level prediction   │ Read predictions for entire fetch  │ Maximizes prediction     │ Only applicable in full  │
  │ with per-instruction     │ block; interpret post-Decode by    │ bandwidth utilization    │ processor pipeline       │    
  │ interpretation           │ instruction type                   │                          │ context                  │
  └──────────────────────────┴────────────────────────────────────┴──────────────────────────┴──────────────────────────┘    
                  
  Key transferable insight: Probabilistic u-bit decay for entries whose predictions are consistently overridden (by SC or    
  loop predictor) but never register as formal mispredictions. These "silently wrong" entries consume space without
  contributing.                                                                                                              
                  
  ---
  12. online ml bp.pdf — SONet (Online ML Branch Predictor)
                                                                                                                             
  ┌──────────────────────┬─────────────────────────────────────┬─────────────────────────────┬───────────────────────────┐
  │     Optimization     │            How it works             │           Benefit           │        Limitations        │   
  ├──────────────────────┼─────────────────────────────────────┼─────────────────────────────┼───────────────────────────┤   
  │ H2P branch           │ 32-slot table with 4-bit primary +  │ 16KB SONet + 64KB           │ Only works when           │   
  │ identification +     │ 6-bit secondary counters; periodic  │ TAGE-SC-L: 1.8% MPKI        │ mispredictions            │   
  │ offloading           │ decay every 1K instructions         │ reduction                   │ concentrated in few       │   
  │                      │                                     │                             │ branches                  │
  ├──────────────────────┼─────────────────────────────────────┼─────────────────────────────┼───────────────────────────┤   
  │ Misprediction        │ Top 8 branches cause ~65% of all    │ Validates focusing          │ Not all workloads exhibit │   
  │ concentration        │ mispredictions in high-error traces │ overrider resources on few  │  this pattern             │
  │ observation          │                                     │ branches                    │                           │   
  ├──────────────────────┼─────────────────────────────────────┼─────────────────────────────┼───────────────────────────┤
  │ TAGE entry space     │ When overrider handles a branch,    │                             │ Indirect benefit for      │
  │ reclamation          │ skip TAGE allocation for it         │ Reduces TAGE table pressure │ other branches is         │   
  │                      │                                     │                             │ "minimal"                 │
  ├──────────────────────┼─────────────────────────────────────┼─────────────────────────────┼───────────────────────────┤   
  │ Combined global +    │ 47 global + 12 local history inputs │ More effective learning of  │ Sweet spot exists; too    │   
  │ local history inputs │  to neural network                  │ interaction patterns        │ much history introduces   │
  │                      │                                     │                             │ noise                     │   
  ├──────────────────────┼─────────────────────────────────────┼─────────────────────────────┼───────────────────────────┤
  │                      │                                     │ Avoids destructive          │                           │   
  │ Per-branch dedicated │ Each H2P branch gets its own weight │ interference; higher        │ Can only track 8          │
  │  network             │  set (8 networks in CAM)            │ accuracy than shared        │ branches; ~14KB total     │   
  │                      │                                     │ weights                     │                           │
  ├──────────────────────┼─────────────────────────────────────┼─────────────────────────────┼───────────────────────────┤
  │ Trust tag /          │ Trust tag determines when to use    │ Graceful fallback for       │ —                         │
  │ confidence mechanism │ network vs fall back to TAGE        │ undertrained networks       │                           │   
  └──────────────────────┴─────────────────────────────────────┴─────────────────────────────┴───────────────────────────┘
                                                                                                                             
  ---             
  13. runlts.pdf — RUNLTS (CBP 2025 Submission)
                                                                                                                             
  Optimization: Novel history length series                                                                              
  How it works: Second-order arithmetic (constant 2nd differences) → geometric transition                                    
  Benefit: "Near-optimal" for large-scale predictors; denser medium-history coverage                                         
  Limitations: Doesn't use skewed-associative tables                                                                         
  ────────────────────────────────────────                                                                                   
  Optimization: Large bimodal (128K entries)                                                                                 
  How it works: 5x bimodal growth despite only 3x total budget increase                                                  
  Benefit: Handles large instruction footprints (JS, server workloads)                                                       
  Limitations: 20 KiB of budget                                                                                              
  ────────────────────────────────────────                                                                               
  Optimization: Register-value correlation (sR component)                                                                    
  How it works: 12-bit digest of register values (leading zeros, trailing zeros, LSBs for int; sign+exponent for FP);        
    Tomasulo-like valid tracking; per-bank usefulness selection
  Benefit: Largest single MPKI reduction; improved every trace category                                                      
  Limitations: 6.6 KiB storage; requires execution→predictor path; implementation complexity
  ────────────────────────────────────────
  Optimization: Thrashing detection + dynamic allocation throttling                                                          
  How it works: Repurpose rare counter/u-bit combos as "newly allocated" markers; track success vs eviction ratio
  Benefit: Boosts accuracy on large-footprint traces                                                                         
  Limitations: Moderate benefit; careful tuning needed
  ────────────────────────────────────────
  Optimization: IMLI components (BrIMLI + TaIMLI)                                                                            
  How it works: Count branches/taken-branches in innermost loop; greater weight assigned to IMLI
  Benefit: Reduced MPKI by at least 1 in three traces                                                                        
  Limitations: 20K bits storage
  ────────────────────────────────────────
  Optimization: Call-stack-based history (sC)                                                                                
  How it works: History stored in return stack; on call, top copied to new entry
  Benefit: Boosts accuracy on large-footprint traces                                                                         
  Limitations: 49K bits; overlaps with standard GEHL functionality
  ────────────────────────────────────────
  Optimization: SC meta predictors                                                                                           
  How it works: Two meta predictors (Meta1, Meta2) with confidence matrix for TAGE vs SC selection
  Benefit: More nuanced SC/TAGE arbitration                                                                                  
  Limitations: Additional decision logic complexity
  ────────────────────────────────────────
  Optimization: Overall RUNLTS                                                                                               
  How it works:
  Benefit: MPKI: 3.197 (vs TAGE-SC-L baseline 3.408); CycWpPKI: 140.3 (vs 145.2)                                             
  Limitations: 191.75 KiB total

  ---
  14. tage-sc-l deepdive.pdf — Deep Dive Into TAGE-SC-L (Ros, CBP 2025)
                                                                                                                             
  Optimization: Quadratic-to-super-exponential history series                                                            
  How it works: Starts quadratic (h + d*(n-1) + k), transitions to generalized geometric with linearly increasing multipliers
  Benefit: MPKI 3.4276; simplifies design — only direct-mapped tables needed                                                 
  Limitations: Needs parameter tuning (h1, d, k, f, m, t)                                                                    
  ────────────────────────────────────────                                                                                   
  Optimization: Improved loop predictor confidence                                                                           
  How it works: Confidence = significant bits of (iterations * consecutive_correct_predictions); loop predictor prediction is
                                                                                                                             
    final (not overridden by SC)
  Benefit: Part of combined MPKI 3.4216                                                                                      
  Limitations: Low coverage inherent to loop predictors
  ────────────────────────────────────────                                                                                   
  Optimization: Improved alternate prediction chooser                                                                        
  How it works: 256 chooser counters; considers confidence < 2 (not just 0); tracks second alternate prediction
  Benefit: Part of combined MPKI 3.4216                                                                                      
  Limitations: Only 128 of 256 counters effectively used
  ────────────────────────────────────────
  Optimization: Improved TAGE vs SC chooser                                                                                  
  How it works: 16 saturating chooser counters indexed by TAGE confidence + SC sum bin
  Benefit: Part of combined MPKI 3.4216                                                                                      
  Limitations: Some counters unused
  ────────────────────────────────────────
  Optimization: Skip TAGE allocation on loop predictor hits                                                                  
  How it works: Bypass TAGE allocation when loop predictor provides prediction
  Benefit: Part of MPKI 3.4156                                                                                               
  Limitations: Only helps loop-heavy workloads
  ────────────────────────────────────────
  Optimization: Improved u-bit: set when alternate correct but hit bank wrong                                                
  How it works: Alternate bank u increased when it correctly predicts but hit bank doesn't
  Benefit: Part of MPKI 3.4156                                                                                               
  Limitations: Minor optimization
  ────────────────────────────────────────
  Optimization: Loop predictor early direction change                                                                        
  How it works: On wrong direction after allocation, assume opposite direction, reset iterations to 2
  Benefit: Part of final MPKI 3.4120                                                                                         
  Limitations: Minor edge-case fix
  ────────────────────────────────────────
  Optimization: Direct-mapped tables (simplification)                                                                        
  How it works: New history series gives equal weight to all lengths — no set-associativity needed
  Benefit: Simplifies design (no associative search, no PRNG for allocation)                                                 
  Limitations: Needs enough tables (44 in paper)

  ---
  15. tage-sc.pdf — TAGE with Statistical Corrector (Seznec)
                                                                                                                             
  ┌────────────────────────────┬─────────────────────────────────┬───────────────────────────────┬──────────────────────┐
  │        Optimization        │          How it works           │            Benefit            │     Limitations      │    
  ├────────────────────────────┼─────────────────────────────────┼───────────────────────────────┼──────────────────────┤
  │                            │ 4+ tables indexed by global     │                               │                      │    
  │ GEHL-based SC (global      │ history at geometric lengths;   │ ~6.8% misprediction reduction │ ~46 Kbits storage    │
  │ history)                   │ 5-6 bit counters; prediction =  │  (256Kbit)                    │                      │    
  │                            │ sign(sum)                       │                               │                      │
  ├────────────────────────────┼─────────────────────────────────┼───────────────────────────────┼──────────────────────┤    
  │                            │                                 │ Removing TAGE direction from  │                      │
  │ Bias indexed by TAGE       │ SC bias indexed by PC + TAGE    │ bias index increases          │ —                    │    
  │ output                     │ predicted direction             │ mispredictions by 1.6% — most │                      │
  │                            │                                 │  influential single SC table  │                      │
  ├────────────────────────────┼─────────────────────────────────┼───────────────────────────────┼──────────────────────┤    
  │                            │ History stored with return      │                               │ Marginal ROI; half   │
  │ Return-stack-associated    │ stack; captures pre-call        │ 0.7% misprediction reduction  │ the benefit          │    
  │ history                    │ correlations                    │                               │ recaptured by        │
  │                            │                                 │                               │ doubling global GEHL │
  ├────────────────────────────┼─────────────────────────────────┼───────────────────────────────┼──────────────────────┤
  │ Multiple local histories   │ 3 local history tables with     │ 2nd local: ~0.5%; 3rd with    │ Diminishing returns  │
  │ (different hashing)        │ different sizes and hash        │ different hash: similar       │ per table            │    
  │                            │ functions                       │                               │                      │
  ├────────────────────────────┼─────────────────────────────────┼───────────────────────────────┼──────────────────────┤    
  │                            │ 3 counters at 3 SC confidence   │                               │ Simple version       │
  │ Confidence-based TAGE/SC   │ levels; TAGE wins when high     │ 0.7% misprediction reduction  │ nearly matches       │    
  │ chooser                    │ confidence + SC low confidence  │                               │ complex GEHL+LGEHL   │
  │                            │                                 │                               │ chooser              │    
  ├────────────────────────────┼─────────────────────────────────┼───────────────────────────────┼──────────────────────┤
  │ Dynamic threshold for SC   │ PC-indexed threshold table (32  │                               │                      │
  │ updates                    │ entries) + global threshold     │ 0.1% misprediction decrease   │ Very small benefit   │    
  │                            │ counter                         │                               │                      │
  ├────────────────────────────┼─────────────────────────────────┼───────────────────────────────┼──────────────────────┤    
  │ Loop predictor (64         │ Constant-iteration loops;       │ ~1% reduction (256Kbit),      │ Diminishing returns  │
  │ entries, 4-way skewed)     │ age-based replacement; SLIM for │ ~0.4% (no-limit)              │ at larger budgets    │    
  │                            │  speculative iterations         │                               │                      │
  ├────────────────────────────┼─────────────────────────────────┼───────────────────────────────┼──────────────────────┤    
  │ Corrector Filter (32Kbit   │ 64-entry associative table;     │                               │ Only for small       │
  │ SC)                        │ inverts TAGE prediction on      │ ~1.5% reduction               │ budget (32Kbit)      │    
  │                            │ high-confidence hits            │                               │                      │
  ├────────────────────────────┼─────────────────────────────────┼───────────────────────────────┼──────────────────────┤    
  │ TAGE saturates at 18-20    │ Beyond 20 tables, invest in SC  │ Guides design tradeoffs       │ —                    │
  │ tables                     │ rather than more TAGE tables    │                               │                      │    
  ├────────────────────────────┼─────────────────────────────────┼───────────────────────────────┼──────────────────────┤
  │ Allocation not in          │ Allocated entries must be in    │ Reduces destructive           │ —                    │    
  │ consecutive tables         │ non-consecutive tagged tables   │ interference                  │                      │    
  └────────────────────────────┴─────────────────────────────────┴───────────────────────────────┴──────────────────────┘
                                                                                                                             
  ---  
