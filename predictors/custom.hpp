#pragma once

#ifndef CUSTOM_HPP
#define CUSTOM_HPP

#include "../cbp.hpp"
#include "../harcom.hpp"
#include "custom/Tage.hpp"

using namespace hcm;

// Custom is the competition-entry facade.
//
// The original scaffold in this file was a dead shell that returned zeros and
// duplicated TAGE-specific parameters in a second, drifting interface. The
// live predictor architecture is now the pluggable overrider TAGE core in
// predictors/custom/Tage.hpp, so Custom is intentionally a thin alias over
// that implementation.
//
// This keeps the user-facing "Custom<...>" entry point while preserving the
// CBP-NG-native structure:
//   - predict1(): block-start prefetch
//   - predict2()/reuse_predict2(): combinational lookup / override
//   - update_condbr(): per-branch checkpointing
//   - update_cycle(): staged training + extra-cycle writeback
template <
    typename TableCfg = SweepTableConfig<8, 512, 11, 1, 2, 1, 2, 100, 4>,
    typename AllocCfg = DefaultAllocConfig,
    u64 FETCH_WIDTH = 16,
    u64 BIMODAL_SIZE = 4096,
    u64 BR_P_ENTRY = 1,
    u64 NUM_BANKS = 1,
    bool USE_AHEAD = false,
    bool SHARED_TAG = true,
    bool SHARED_U = true,
    bool SHARED_HYS = true,
    bool U_STOR_FF = false,
    u64 DECAY_CTR = 8,
    u64 DECAY_GRAN = 2,
    typename DECAY_POLICY = DecayMild,
    typename ResetFn = DefaultResetFn,
    bool USE_FF_CACHE = false,
    bool P1_USE_GSHARE = true,
    u64 P1_TABLE_SIZE = 4096,
    u64 P1_HIST = 6,
    bool USE_META = true,
    u64 METABITS = 4,
    u64 METAPIPE = 2,
    u64 META_TABLE_SIZE = 0,
    bool USE_PATH_HIST = false,
    u64 PATH_HIST_WIDTH = 27,
    u64 PATH_BITS = 6,
    typename Overrider = SCOverrider<3, 8, 6, FETCH_WIDTH>>
using Custom =
    Tage<TableCfg, AllocCfg, FETCH_WIDTH, BIMODAL_SIZE, BR_P_ENTRY, NUM_BANKS,
         USE_AHEAD, SHARED_TAG, SHARED_U, SHARED_HYS, U_STOR_FF, DECAY_CTR,
         DECAY_GRAN, DECAY_POLICY, ResetFn, USE_FF_CACHE, P1_USE_GSHARE,
         P1_TABLE_SIZE, P1_HIST, USE_META, METABITS, METAPIPE,
         META_TABLE_SIZE, USE_PATH_HIST, PATH_HIST_WIDTH, PATH_BITS,
         Overrider>;

#endif // CUSTOM_HPP
