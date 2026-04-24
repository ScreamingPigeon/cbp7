#include "predictors/always_taken.hpp"
#include "predictors/bimodal.hpp"
#include "predictors/bimodalN.hpp"
#include "predictors/custom/Tage.hpp"
#include "predictors/custom/TageAhead.hpp"
#include "predictors/custom/TageDirect.hpp"
#include "predictors/custom/TageDirectBim.hpp"
#include "predictors/experiment_perceptron.hpp"
#include "predictors/gshare.hpp"
#include "predictors/gshareN.hpp"
#include "predictors/gshareN_ahead.hpp"
#include "predictors/gshareN_ahead_best.hpp"
#include "predictors/never_taken.hpp"
#include "predictors/perceptron.hpp"
#include "predictors/tage.hpp"
#include "predictors/tutorial/tutorial.hpp"

// Keep `perceptron<>` as a stable user-facing name.
template <auto... Args> using perceptron = experiment_perceptron<Args...>;

// ============================================================================
// Competition configs: 1-cycle and 2-cycle tracks
// ============================================================================

// 1-cycle config: 14 tables, 1024 uniform, P2 ≈ 0.99
using TageAhead1C =
    TageAhead<TATableConfig<14, 1024, 11, 4, 300, 1, ta::HistSeries::GEOMETRIC,
                            ta::UniformTag<11>, ta::UniformSize<1024>>,
              8, 6, 4, true, 1, ta::Xor3SecTagHash>;

// 2-cycle config: 28 tables, StepSize 4096/2048 split@24, P2 ≈ 1.91
// 106K entries, 24 tables at 4096 + 4 at 2048, MAXH=1000
// SiblingPolicy::NONE — sibling skip hurts 2C (3.3% MPKI regression)
using TageAhead2C =
    TageAhead<TATableConfig<28, 2048, 12, 4, 1000, 2, ta::HistSeries::GEOMETRIC,
                            ta::UniformTag<12>, ta::StepSize<4096, 2048, 24>>,
              8, 6, 4, true, 1, ta::Xor3SecTagHash, 1, 2, 1,
              UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT,
              4096, false, 6,
              4, 256, 2, 1024, true, HistUpdate::PATH, TADefaultAllocConfig,
              SiblingPolicy::NONE>;

// ============================================================================
// Decay sweep configs (all 2C base + SiblingPolicy::NONE)
// Common: 28 tables, StepSize<4096,2048,24>, TAG=12, MAXH=1000
// ============================================================================

// Helper: 2C base params up to SiblingPolicy::NONE, parameterized by U_WIDTH
#define TA2C_BASE_U(UW)                                                        \
  TATableConfig<28, 2048, 12, 4, 1000, 2, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<12>, ta::StepSize<4096, 2048, 24>>,             \
      8, 6, 4, true, 1, ta::Xor3SecTagHash, 1, 2, UW,                           \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 4096, false, 6, 4,       \
      256, 2, 1024, true, HistUpdate::PATH, TADefaultAllocConfig,              \
      SiblingPolicy::NONE, 0, 16, 16
#define TA2C_BASE TA2C_BASE_U(1)

// A: Baseline — epoch enabled, no decay (already TageAhead2C)

// B: No management — neither epoch nor decay
using TA2C_B = TageAhead<TA2C_BASE, false, DecayMiss::TAG,
                         DecayOp::DECREMENT, // DECAY_ENABLE=false
                         ta::uniform_array<u64, 28>(8), ta::DefaultDecayThresh,
                         false>; // EPOCH_ENABLE=false

// C: TAG_OR_SEC, LW=8 (default decay, most aggressive trigger)
using TA2C_C =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, 28>(8), ta::DefaultDecayThresh, false>;

// D: TAG only, LW=8 (less aggressive trigger)
using TA2C_D =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG, DecayOp::DECREMENT,
              ta::uniform_array<u64, 28>(8), ta::DefaultDecayThresh, false>;

// E: TAG_AND_SEC, LW=8 (least aggressive trigger)
using TA2C_E =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG_AND_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, 28>(8), ta::DefaultDecayThresh, false>;

// F: TAG_OR_SEC, LW=12 (lower fire probability)
using TA2C_F =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, 28>(12), ta::DefaultDecayThresh, false>;

// G: TAG_OR_SEC, LW=6 (higher fire probability)
using TA2C_G =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, 28>(6), ta::DefaultDecayThresh, false>;

// H: TAG_OR_SEC, per-table LW: 12 for T0-13, 6 for T14-27
using TA2C_H =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::split_array<u64, 28>(12, 6, 14), ta::DefaultDecayThresh,
              false>;

#ifdef PREDICTOR
using branch_predictor = PREDICTOR;
#else
using branch_predictor = TageAhead1C;
#endif
