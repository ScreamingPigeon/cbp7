#include "predictors/always_taken.hpp"
#include "predictors/bimodal.hpp"
#include "predictors/bimodalN.hpp"
#include "predictors/custom/Tage.hpp"
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

#ifdef PREDICTOR
using branch_predictor = PREDICTOR;
#else
// using branch_predictor = bimodal<>;
// using branch_predictor = gshare<>;
// using branch_predictor = tage<>;
using branch_predictor = TageDirectBim<>;
// using branch_predictor = experiment_perceptron<>;
#endif
