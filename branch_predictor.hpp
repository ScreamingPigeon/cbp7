#include "predictors/always_taken.hpp"
#include "predictors/bimodal.hpp"
#include "predictors/bimodalN.hpp"
#include "predictors/custom/Tage.hpp"
#include "predictors/custom/TageDirect.hpp"
#include "predictors/custom/TageAhead.hpp"
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
using TageAhead1C = TageAhead<
    TATableConfig<14, 1024, 11, 4, 500, 1,
                  ta::HistSeries::GEOMETRIC,
                  ta::UniformTag<11>,
                  ta::UniformSize<1024>>>;

// 2-cycle config: 28 tables, StepSize 4096/2048 split@24, P2 ≈ 1.91
// 106K entries, 24 tables at 4096 + 4 at 2048, MAXH=1000
// SiblingPolicy::NONE — sibling skip hurts 2C (3.3% MPKI regression)
using TageAhead2C = TageAhead<
    TATableConfig<28, 2048, 12, 4, 1000, 2,
                  ta::HistSeries::GEOMETRIC,
                  ta::UniformTag<12>,
                  ta::StepSize<4096, 2048, 24>>,
    8, 6, 3, true, 1, ta::DefaultSecTagHash, 1, 2, 1, 4096, false, 6,
    4, 256, 2, 1024, true, HistUpdate::PATH, TADefaultAllocConfig,
    SiblingPolicy::NONE>;

// Test configs: sibling skip variants of 2C
using TageAhead2C_NoSib = TageAhead<
    TATableConfig<28, 2048, 12, 4, 1000, 2,
                  ta::HistSeries::GEOMETRIC,
                  ta::UniformTag<12>,
                  ta::StepSize<4096, 2048, 24>>,
    8, 6, 3, true, 1, ta::DefaultSecTagHash, 1, 2, 1, 4096, false, 6,
    4, 256, 2, 1024, true, HistUpdate::PATH, TADefaultAllocConfig,
    SiblingPolicy::NONE>;

using TageAhead2C_SibFloor14 = TageAhead<
    TATableConfig<28, 2048, 12, 4, 1000, 2,
                  ta::HistSeries::GEOMETRIC,
                  ta::UniformTag<12>,
                  ta::StepSize<4096, 2048, 24>>,
    8, 6, 3, true, 1, ta::DefaultSecTagHash, 1, 2, 1, 4096, false, 6,
    4, 256, 2, 1024, true, HistUpdate::PATH, TADefaultAllocConfig,
    SiblingPolicy::ALL, 14>;

#ifdef PREDICTOR
using branch_predictor = PREDICTOR;
#else
using branch_predictor = TageAhead1C;
#endif
