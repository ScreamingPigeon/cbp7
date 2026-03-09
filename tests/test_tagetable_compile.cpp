// Phase 1 compile test: verify TageTable instantiates with various parameter
// combinations. No runtime logic — just checks that all storage members are
// declared correctly.

#include "../predictors/custom/TageTable.hpp"

// Default configuration
static TageTable<> t_default;

// Minimal: 1 branch, 1 bank, no ahead, no cache
static TageTable<64, 8, 7, 3, 1, 1, 1, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t_minimal;

// Multi-bank with shared tag/u
static TageTable<256, 32, 11, 3, 2, 8, 4, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t_multibank_shared;

// Multi-bank with per-bank tag/u (SHARED_TAG=false, SHARED_U=false)
static TageTable<256, 32, 11, 3, 2, 8, 4, false, false, false, true, 1024,
                 DefaultResetFn, false>
    t_multibank_perbank;

// Per-bank tag, shared u (SHARED_TAG=false, SHARED_U=true)
static TageTable<256, 32, 11, 3, 2, 8, 4, false, false, true, true, 1024,
                 DefaultResetFn, false>
    t_perbank_tag_shared_u;

// Ahead pipelining enabled
static TageTable<512, 64, 11, 3, 1, 4, 2, true, true, true, true, 1024,
                 DefaultResetFn, false>
    t_ahead;

// Ahead with per-bank tags
static TageTable<512, 64, 11, 3, 1, 4, 2, true, false, false, true, 1024,
                 DefaultResetFn, false>
    t_ahead_perbank;

// U-bits in SRAM (U_STOR_FF=false)
static TageTable<256, 32, 11, 3, 2, 4, 1, false, true, true, false, 512,
                 DefaultResetFn, false>
    t_u_sram;

// U-bits in SRAM, per-bank
static TageTable<256, 32, 11, 3, 2, 4, 4, false, false, false, false, 512,
                 DefaultResetFn, false>
    t_u_sram_perbank;

// FF cache enabled (BPB > 1)
static TageTable<256, 32, 11, 3, 1, 4, 1, false, true, true, true, 1024,
                 DefaultResetFn, true>
    t_ff_cache;

// FF cache with multi-bank
static TageTable<256, 32, 11, 3, 1, 8, 2, false, true, true, true, 1024,
                 DefaultResetFn, true>
    t_ff_cache_multibank;

// Large table
static TageTable<2048, 128, 12, 4, 2, 4, 1, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t_large;

// Wide counters
static TageTable<256, 32, 11, 4, 2, 4, 1, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t_wide_ctr;

// Everything enabled: ahead + multi-bank + FF cache + per-bank tags
static TageTable<256, 32, 11, 3, 1, 8, 2, true, false, false, true, 1024,
                 DefaultResetFn, true>
    t_everything;

int main() { return 0; }
