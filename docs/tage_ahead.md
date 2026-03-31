  P1: Direct gshare. Not banked, not ahead. Same as Tage.hpp. P1 latency = 0.957.
                                                                                                                                                  
  P2 TAGE: Predecessor-indexed, 1-block ahead.     
  - predict1(N) reads TAGE tables at hash(PC_N) — these entries store predictions for blocks that follow block N
  - predict2(N+1) uses those reads (shifted via pipeline regs)                                                                                    
  - Entries are trained during update_cycle(N+1): "block N+1 followed block N, prediction should be Y"
                                                                                                                                                  
  Banking: BANKS=2 (PATHS=2). Each TAGE RAM entry stores 2 banks:                                                                                 
  - Bank 0: prediction for fall-through successor                                                                                                 
  - Bank 1: prediction for taken-branch successor
  - Bank selected by path ^ XB (path set in update_cycle from actual branch outcome)                                                              
  - 2x TAGE RAM width, but P1/bimodal NOT banked (big EPI saving vs previous attempt)
                                                                                                                                                  
  Secondary tag: Small tag (2-3 bits) stored per entry to disambiguate when BANKS=2 isn't enough (rare 3rd+ path, or different branches mapping to
   same bank). Compared at predict2 time — mismatch → bimodal fallback.                                                                           
                                                                                                                                                  
  Bimodal: Direct (not ahead). Read in predict2 as fallback when secondary tag misses.                                                            
                                                                                                                                                  
  Pipeline:                                                                                                                                       
  predict1(N):  shift ahead_*[0]→[1], read TAGE at PC_N → ahead_*[0]                                                                              
  predict2(N):  use ahead_*[1] (read at PC_{N-1}), bank select via path^XB,                                                                       
                secondary tag check, fallback to bimodal on miss                                                                                  
  update_cycle(N): write to tage_indexes[1] (PC_{N-1}'s index), set path                                                                          
                                                                                                                                                  
  Target: predict2 < 1.0 cycle (combinational: bank mux + sec tag compare + bimodal fallback mux)                                                 
