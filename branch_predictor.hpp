#include "predictors/always_taken.hpp"
#include "predictors/bimodal.hpp"
#include "predictors/bimodalN.hpp"
#include "predictors/custom.hpp"
#include "predictors/custom/Tage.hpp"
#include "predictors/custom/LoopPredictorA.hpp"
#include "predictors/custom/LoopPredictorC.hpp"
#include "predictors/experiment_perceptron.hpp"
#include "predictors/gshare.hpp"
#include "predictors/never_taken.hpp"
#include "predictors/perceptron.hpp"
#include "predictors/tage.hpp"
#include "predictors/bimodalN.hpp"
#include "predictors/gshareN.hpp"
#include "predictors/gshareN_ahead.hpp"
#include "predictors/tutorial/tutorial.hpp"
#include "predictors/experiment_perceptron.hpp"
#include "predictors/gshareN_loop.hpp"
// NOTE: TageAhead pulls in TageConfig utilities that collide with
// predictors/custom/Tage.hpp in this TU. Keep it out of the default
// branch predictor include set unless specifically wiring an ahead predictor.
// #include "predictors/custom/TageAhead.hpp"

// Keep `perceptron<>` as a stable user-facing name.
template <auto... Args>
using perceptron = experiment_perceptron<Args...>;

#ifdef PREDICTOR
using branch_predictor = PREDICTOR;
#else
// using branch_predictor = bimodal<>;
// using branch_predictor = gshare<>;
// using branch_predictor = tage<>;
// using branch_predictor = experiment_perceptron<>;
using branch_predictor = gshareN_loop<>;
#endif
