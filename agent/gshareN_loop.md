# Agent Guide: Constrained Parameter Sweep for `gshareN_loop`

This document tells an agent how to tune the top-level template parameters of `gshareN_loop` in the CBP-NG framework to maximize **`P2_cond_acc`** while enforcing the hard constraint that **both `predict2` and `reuse_predict2` latency stay below 2 cycles**.

It is written for the CBP-NG workflow in the AmpereComputing `cbp-ng` repo, where predictors are built by selecting a `PREDICTOR` template instantiation and then run on traces for accuracy, latency, energy, and VFS-style metrics. CBP-NG’s README also recommends using the official training traces and explains that predictors are implemented in HARCOM, which models latency, energy, and area directly.

---

## 1. Goal

Find a parameter setting for:

```cpp
gshareN_loop<
    LOGG, GHIST, LOGN, LINEINST_,
    LOOP_LOGSETS, LOOP_TAGW, LOOP_ITERW,
    LOOP_CONFW, LOOP_AGEW, LOOP_WAYS,
    MIN_TRIP, USE_LONG_LOOP_CONF,
    LONG_LOOP_CONF_PRODUCT_THRESH,
    LOOP_CONF_USE,
    USE_CTRW, USE_CTR_INIT, USE_CTR_OVERRIDE_THRESH,
    USE_CTR_INC_STEP, USE_CTR_DEC_STEP,
    ALLOC_CONF_INIT, ALLOC_AGE_INIT,
    CONF_BUMP_STEP,
    INVALIDATE_ON_WRONG_CONF, INVALIDATE_ON_EXCEEDED_ITER,
    INVALIDATE_ON_BAD_EXIT_MATCH, INVALIDATE_ON_TOO_SHORT,
    REQUIRE_FIRST_BRANCH_MISPRED_FOR_ALLOC,
    FORCE_EXTRA_CYCLE_ON_LOOKUP
>
```

that maximizes:

- **Primary objective:** `P2_cond_acc`

subject to:

- **Hard constraint 1:** `predict2 LatMax < 2.0`
- **Hard constraint 2:** `reuse_predict2 LatMax < 2.0`

Recommended secondary tie-breakers:

1. Higher `P2_cond_acc`
2. Lower `predict2 LatAvg`
3. Lower `reuse_predict2 LatAvg`
4. Lower `EPI(from per-call deltas)` or lower average dynamic energy if accuracy is nearly tied

Why use `LatMax` as the hard constraint instead of `LatAvg`? Because the user said latency must be “below 2.” If the maximum observed latency is above 2, then the design is not safely under budget even if the average is. The CBP-NG tooling reports both average and max predictor latency in cycles, so the safe interpretation is to enforce the bound on the maximum. The README also makes clear that predictor latency is a first-class modeled output of the simulator.

---

## 2. What the current results already say

From the user’s current runs:

- `P2_cond_acc = 0.930311`
- `predict2 LatAvg = 1.62132`, `LatMax = 1.86000`
- `reuse_predict2 LatAvg = 1.60707`, `LatMax = 1.86000`

This means the current design is already **feasible** under the `< 2` latency budget, but only with around **0.14 cycles of headroom** on the max. That is not much. So the sweep should prefer settings that either:

- keep the loop predictor logic size similar but improve learning/use behavior, or
- reduce logic on the P2 critical path slightly before trying to add more complexity.

Do **not** start by aggressively increasing table complexity or adding more tag/way logic everywhere. The current latency margin is too small for reckless expansion.

---

## 3. Repo facts that matter to the sweep

The CBP-NG README says:

- predictors are compiled by selecting `-DPREDICTOR="..."` at build time, including template instantiations;
- one-trace simulation is done with `./cbp TRACE NAME WARMUP MEASURE`;
- `run_all` exists for many traces;
- the official training set contains **168 traces** and is the recommended development set;
- the simulator output includes P1/P2 latency and dynamic energy fields;
- `predictor_metrics.py` and `vfs.py` exist for downstream aggregate analysis.

So the agent should not choose parameters from a single trace alone. It can use one trace for fast filtering, but the final decision must be made on a multi-trace set, ideally the official training traces.

---

## 4. Parameter meanings and expected effects

This section explains each top-level parameter in practical terms and how it is likely to affect:

- P2 accuracy
- P2 latency
- energy/area
- stability/generalization

### 4.1 Base predictor structure

#### `LOGG`
Base gshare table log-size. Since `index_bits = LOGG - LOGN`, larger `LOGG` means more indexed state in the base predictor.

Expected effect:
- Larger `LOGG` can improve the base predictor’s discrimination and reduce aliasing.
- It can indirectly affect loop usefulness because the loop predictor’s override usefulness depends partly on disagreement vs base behavior.
- Usually modest effect on P2 latency unless the implementation makes this lookup or downstream fanout larger.
- Can increase SRAM area and energy.

Guidance:
- Sweep locally around the current value first.
- Suggested range: `13..17`

#### `GHIST`
Global history length used in gshare index formation.

Expected effect:
- Too small: misses long-range correlation.
- Too large: more aliasing/noise for a small table.
- Changes the base predictor behavior and therefore the loop predictor’s override opportunities.

Guidance:
- Strong interaction with `LOGG`.
- Suggested range: `4..16`

#### `LOGN`
Number of branch positions per prediction block as `N = 1 << LOGN`.

Expected effect:
- Larger `LOGN` means more per-block branch positions and wider handling logic.
- Likely increases fanout, matching, and selection complexity.
- Very likely to hurt `predict2` and `reuse_predict2` latency if increased too much.

Guidance:
- Treat this as a structural parameter, not a fine tuning knob.
- Only sweep `LOGN` if you are willing to re-check latency carefully.
- Suggested range: keep near current value; test `2`, maybe `1`, and `3` only if needed.

#### `LINEINST_`
Prediction block size in instructions.

Expected effect:
- Larger blocks can improve reuse rate and amortize first-lookup cost.
- But they also widen block bookkeeping, fanout, and branch-position logic.
- Large values may increase complexity more than they help.

Guidance:
- Treat as structural.
- Only sweep a small set of powers of two supported cleanly by the implementation.
- Suggested range: current value and one smaller alternative first.

---

### 4.2 Loop table geometry

#### `LOOP_LOGSETS`
Loop predictor set count as `2^LOOP_LOGSETS`.

Expected effect:
- More sets reduce destructive aliasing between loops.
- Usually one of the safest ways to improve loop accuracy, because it mostly changes SRAM capacity rather than deep combinational logic.
- Area and energy increase.
- Latency impact is often smaller than adding more matching logic, but still must be checked.

Guidance:
- Good early sweep variable.
- Suggested range: `3..7`

#### `LOOP_WAYS`
Associativity of the loop table.

Expected effect:
- More ways reduce conflict misses.
- But way matching is on the prediction path. More ways means more tag compares and more select logic.
- This can directly hurt `predict2`/`reuse_predict2` latency.

Guidance:
- Because the current design already has only ~0.14 cycle max-latency margin, increasing `LOOP_WAYS` is risky.
- Prefer to try `1` and `2` first; only test `4` if many feasible points remain and accuracy stalls.

#### `LOOP_TAGW`
Tag width for loop entries.

Expected effect:
- Larger tag width lowers false hits and wrong-loop matching.
- Can help accuracy if aliasing is significant.
- More compare bits can add energy and possibly a little latency.

Guidance:
- Usually a medium-sensitivity parameter.
- Suggested range: `8..14`

#### `LOOP_ITERW`
Stored trip-count / current-iteration counter width.

Expected effect:
- Too small: loops with larger trip counts saturate or learn poorly.
- Too large: more storage bits, maybe slightly more datapath cost.
- Usually a modest latency effect but a real storage/energy effect.

Guidance:
- Good moderate sweep parameter.
- Suggested range: `6..12`

---

### 4.3 Loop confidence and aging

#### `LOOP_CONFW`
Confidence counter width.

Expected effect:
- Larger confidence width gives more granularity and can reduce flip-flopping.
- Too large can slow adaptation.
- Small storage cost increase.

Guidance:
- Suggested range: `2..4`

#### `LOOP_AGEW`
Age counter width for victim selection / entry persistence.

Expected effect:
- Larger age width can retain useful loop entries longer.
- Too large can make replacement too sticky.
- Mostly affects replacement behavior, not critical-path logic very much.

Guidance:
- Good sweep candidate.
- Suggested range: `2..6`

#### `LOOP_CONF_USE`
Threshold for using the loop predictor based on confidence.

Expected effect:
- Lower threshold makes loop override more aggressive.
- Higher threshold makes it more conservative.
- This is often a high-impact behavioral knob with very little hardware cost change.

Guidance:
- Very good tuning parameter.
- Must satisfy `LOOP_CONF_USE <= (1<<LOOP_CONFW)-1`.
- Sweep all legal values for each chosen `LOOP_CONFW`.

#### `USE_LONG_LOOP_CONF`
Enables the alternate confidence rule using a product of confidence and trip count.

Expected effect:
- Helps trusted long loops become usable earlier or more reliably.
- Can improve accuracy on steady long-running loops.
- May slightly increase logic on the decision path.

Guidance:
- Binary toggle; easy to test.
- Keep if it improves feasible points; otherwise disable to simplify behavior.

#### `LONG_LOOP_CONF_PRODUCT_THRESH`
Threshold used when `USE_LONG_LOOP_CONF` is true.

Expected effect:
- Lower threshold makes long loops become “confident” sooner.
- Higher threshold makes the predictor more conservative.

Guidance:
- Strongly interacts with `LOOP_CONFW`, `LOOP_ITERW`, and `LOOP_CONF_USE`.
- Suggested range: powers of two or near-powers like `16, 32, 64, 96, 128, 192, 256`.

---

### 4.4 Loop-use meta predictor

#### `USE_CTRW`
Width of the meta counter controlling whether loop override is allowed.

Expected effect:
- Larger width makes the override-use decision more stable.
- Smaller width adapts faster.
- Mostly a behavioral parameter, usually cheap.

Guidance:
- Good sweep parameter.
- Suggested range: `3..6`

#### `USE_CTR_INIT`
Initial value of the use counter.

Expected effect:
- Positive values bias toward trusting the loop predictor early.
- Negative values bias toward the base predictor early.

Guidance:
- Important when the meta predictor otherwise learns too slowly.
- Suggested range: `[-4, 4]`, clipped to legal range implied by `USE_CTRW`.

#### `USE_CTR_OVERRIDE_THRESH`
Threshold above which the loop predictor may override P1.

Expected effect:
- Lower threshold means more loop overrides.
- Higher threshold means fewer overrides.
- Often one of the most accuracy-sensitive knobs once feasible latency is established.

Guidance:
- Sweep around zero and small positive values.
- Suggested range depends on `USE_CTRW`; typically `[-2, -1, 0, 1, 2, 3]` when legal.

#### `USE_CTR_INC_STEP`, `USE_CTR_DEC_STEP`
Step sizes for meta counter adaptation.

Expected effect:
- Larger steps adapt faster but can be noisy.
- Asymmetric steps let you prefer caution or aggressiveness.

Guidance:
- Test a small set like `{1,2}` for each.
- Often `inc=1, dec=1` or `inc=1, dec=2` are reasonable.

---

### 4.5 Allocation and confidence update policy

#### `ALLOC_CONF_INIT`
Initial confidence of a newly allocated loop entry.

Expected effect:
- Higher values make new entries trusted sooner.
- Risk: noisy entries get used too early.

Guidance:
- Usually low values are safer.
- Sweep `0..min(2, CONF_MAX)` first.

#### `ALLOC_AGE_INIT`
Initial age assigned to newly allocated loop entries.

Expected effect:
- Larger values keep new entries alive longer.
- Too large can preserve junk.

Guidance:
- Good replacement-policy knob.
- Sweep low, medium, max: e.g. `{0, AGE_MAX/2, AGE_MAX}`.

#### `CONF_BUMP_STEP`
How much confidence increases when confirming the learned trip behavior.

Expected effect:
- Larger step means faster convergence to confident use.
- Can also over-trust unstable loops.

Guidance:
- Sweep `{1,2}`; maybe `3` only if `LOOP_CONFW` is wide enough.

---

### 4.6 Loop entry invalidation / cleanup policies

These booleans can have very strong behavioral effects and are cheap to test.

#### `INVALIDATE_ON_WRONG_CONF`
Invalidate an entry if it provided a confident wrong prediction.

Expected effect:
- More aggressive cleanup of bad learned entries.
- May improve robustness but can overreact.

#### `INVALIDATE_ON_EXCEEDED_ITER`
Invalidate if current iteration exceeds stored trip count.

Expected effect:
- Good if exceeded-trip count means the entry is stale/wrong.
- Bad if many loops vary slightly around a nominal count.

#### `INVALIDATE_ON_BAD_EXIT_MATCH`
Invalidate on exit mismatch when expected trip count does not match.

Expected effect:
- More aggressive when loop-exit learning is wrong.

#### `INVALIDATE_ON_TOO_SHORT`
Invalidate when observed trip count is below `MIN_TRIP`.

Expected effect:
- Prevents filling the loop table with short loops that are not worth specializing.
- Usually a helpful regularization knob.

Guidance for all four:
- Binary parameters: easy to include in sweep.
- These often improve generalization more than adding capacity.

---

### 4.7 Allocation gating and latency control

#### `REQUIRE_FIRST_BRANCH_MISPRED_FOR_ALLOC`
Only allocate a loop entry when the block had a misprediction and there was exactly one conditional branch.

Expected effect:
- More conservative allocation.
- Reduces pollution.
- Can miss some useful loops in multi-branch blocks.

Guidance:
- Important behavioral knob.
- Test both values.

#### `FORCE_EXTRA_CYCLE_ON_LOOKUP`
Requests an extra cycle when lookup/mispredict conditions apply.

Expected effect:
- Helps model realistic pipelining / latency constraints.
- Disabling it may reduce reported extra-cycle behavior, but be careful: if the actual modeled path remains too long, it may not really help the “below 2 cycles” requirement in the intended way.

Guidance:
- Test both values, but do not use it as a cheat. Prefer settings that remain safe and accurate under the more realistic option.
- If disabling it helps latency but hurts robustness or realism, reject it.

---

### 4.8 `MIN_TRIP`
Minimum loop trip count considered worth learning.

Expected effect:
- Larger values avoid wasting entries on short loops.
- Too large may ignore useful short, regular loops.

Guidance:
- High-impact regularization knob.
- Suggested range: `2..8`

---

## 5. Which parameters matter most first

Do **not** sweep everything uniformly from the start. Use this priority order.

### Tier A: tune these first
These are high leverage and usually low risk for latency:

- `LOOP_LOGSETS`
- `LOOP_TAGW`
- `LOOP_ITERW`
- `LOOP_CONF_USE`
- `MIN_TRIP`
- `USE_LONG_LOOP_CONF`
- `LONG_LOOP_CONF_PRODUCT_THRESH`
- `USE_CTRW`
- `USE_CTR_INIT`
- `USE_CTR_OVERRIDE_THRESH`
- `USE_CTR_INC_STEP`
- `USE_CTR_DEC_STEP`
- `ALLOC_CONF_INIT`
- `ALLOC_AGE_INIT`
- `CONF_BUMP_STEP`
- the invalidation booleans
- `REQUIRE_FIRST_BRANCH_MISPRED_FOR_ALLOC`

### Tier B: tune these once Tier A stalls
These can affect latency more or interact more strongly:

- `LOOP_WAYS`
- `LOOP_CONFW`
- `LOOP_AGEW`
- `LOGG`
- `GHIST`

### Tier C: structural / risky
Change only when you already understand the design:

- `LOGN`
- `LINEINST_`
- `FORCE_EXTRA_CYCLE_ON_LOOKUP`

---

## 6. Constraint-aware optimization strategy

Use a **constrained surrogate-model workflow**, not just brute-force ranking on one trace.

### Recommended statistical model
Use **constrained Bayesian optimization** with two surrogate models:

1. **Accuracy model**
   - Target: mean `P2_cond_acc`
   - Suggested surrogate: Gaussian process regression, random forest regressor, or gradient-boosted regressor

2. **Feasibility model**
   - Target: probability that the point satisfies both:
     - `predict2 LatMax < 2.0`
     - `reuse_predict2 LatMax < 2.0`
   - Suggested surrogate: logistic regression, random forest classifier, or gradient-boosted classifier

Score candidate points using:

```text
acquisition(x) = P(feasible | x) * U_acc(x)
```

where `U_acc(x)` is an upper-confidence or expected-improvement style utility for accuracy.

If Bayesian optimization tooling is inconvenient, a strong simpler alternative is:

- fit a regression model for `P2_cond_acc`
- fit a classifier for feasibility
- rank by `predicted_accuracy * predicted_feasibility_probability`
- then evaluate the top candidates exactly

### Why this is better than naive sweeping
Because many parameters interact and the search space is mixed:

- integers
- booleans
- thresholds constrained by widths

A surrogate model will find which dimensions matter most and avoid wasting runs on obviously infeasible high-latency regions.

---

## 7. Exact sweep workflow

### Stage 0: baseline capture
Always log the current baseline first.

For the current predictor configuration, store:

- all template parameters
- `P1_cond_acc`
- `P2_cond_acc`
- `predict2 LatAvg`
- `predict2 LatMax`
- `reuse_predict2 LatAvg`
- `reuse_predict2 LatMax`
- dynamic energy
- total transistors / SRAM bits if available

This is the reference point for all future comparisons.

---

### Stage 1: fast feasibility screen on a small trace subset
Use the user’s existing fast commands or a small representative subset of training traces.

For each candidate:

1. Build predictor with the chosen template parameters.
2. Run a short analyze window to get latency metrics.
3. Reject immediately if:
   - `predict2 LatMax >= 2.0`, or
   - `reuse_predict2 LatMax >= 2.0`

Only feasible points move on.

Why this works:
- Latency is usually much cheaper to check than full long-run accuracy across many traces.
- It quickly removes configurations that are structurally too expensive.

---

### Stage 2: medium accuracy pass on multiple traces
For all feasible points from Stage 1:

- run `acc` mode on a **moderate multi-trace subset**
- aggregate mean `P2_cond_acc`
- also compute variance / standard deviation across traces

Recommended aggregate score for this stage:

```text
screen_score = mean(P2_cond_acc) - 0.25 * stddev(P2_cond_acc)
```

This slightly penalizes unstable designs that spike on one trace but generalize poorly.

If trace count is small, also compute a bootstrap 95% confidence interval for mean `P2_cond_acc`.

---

### Stage 3: fit surrogate models
Using all runs so far:

- fit feasibility model on latency pass/fail
- fit regression model on mean `P2_cond_acc`
- propose next candidates in promising regions

Practical recommendation:
- use random search for the initial 50–200 points
- then switch to model-guided candidate proposal

---

### Stage 4: full evaluation on the official training set
The CBP-NG README recommends the official 168 training traces for development. Use them for the final comparison, not just one or two traces.

For the top feasible candidates:

- run across the full training set using `run_all` or equivalent scripted iteration
- compute:
  - mean `P2_cond_acc`
  - median `P2_cond_acc`
  - stddev across traces
  - bootstrap CI of mean
  - fraction of traces satisfying latency constraint if latency is measured per trace

Final winner rule:

1. reject any candidate with any observed `predict2` or `reuse_predict2` `LatMax >= 2.0`
2. among remaining candidates, choose the one with highest mean `P2_cond_acc`
3. if two are within statistical noise, choose the lower-energy / lower-area / lower-latency one

---

## 8. Parameter domains to sample

Use these domains as the initial legal search space.

```text
LOGG                           ∈ {13,14,15,16,17}
GHIST                          ∈ {4,6,8,10,12,14,16}
LOGN                           ∈ {1,2,3}          # use 3 sparingly
LINEINST_                      ∈ {32,64}          # keep powers of two supported by code

LOOP_LOGSETS                   ∈ {3,4,5,6,7}
LOOP_TAGW                      ∈ {8,10,12,14}
LOOP_ITERW                     ∈ {6,8,10,12}
LOOP_CONFW                     ∈ {2,3,4}
LOOP_AGEW                      ∈ {2,3,4,5,6}
LOOP_WAYS                      ∈ {1,2,4}
MIN_TRIP                       ∈ {2,3,4,5,6,8}

USE_LONG_LOOP_CONF             ∈ {false,true}
LONG_LOOP_CONF_PRODUCT_THRESH  ∈ {16,32,64,96,128,192,256}

USE_CTRW                       ∈ {3,4,5,6}
USE_CTR_INIT                   ∈ legal small integers centered at 0
USE_CTR_OVERRIDE_THRESH        ∈ legal small integers centered at 0
USE_CTR_INC_STEP               ∈ {1,2}
USE_CTR_DEC_STEP               ∈ {1,2}

ALLOC_CONF_INIT                ∈ legal values 0..CONF_MAX, start small
ALLOC_AGE_INIT                 ∈ {0, mid, max}
CONF_BUMP_STEP                 ∈ {1,2}

INVALIDATE_ON_WRONG_CONF       ∈ {false,true}
INVALIDATE_ON_EXCEEDED_ITER    ∈ {false,true}
INVALIDATE_ON_BAD_EXIT_MATCH   ∈ {false,true}
INVALIDATE_ON_TOO_SHORT        ∈ {false,true}
REQUIRE_FIRST_BRANCH_MISPRED_FOR_ALLOC ∈ {false,true}
FORCE_EXTRA_CYCLE_ON_LOOKUP    ∈ {false,true}
```

Important legality rules the agent must enforce:

- `LINEINST_` must be a power of two.
- `LOOP_WAYS` must match what the code supports.
- `LOOP_CONF_USE <= (1<<LOOP_CONFW)-1`.
- `USE_CTR_INIT` and `USE_CTR_OVERRIDE_THRESH` must fit within signed `USE_CTRW` range.
- `ALLOC_CONF_INIT <= CONF_MAX`.
- `ALLOC_AGE_INIT <= AGE_MAX`.

---

## 9. Recommended initial experiment batches

### Batch A: behavioral tuning with low latency risk
Keep geometry mostly fixed and vary:

- `LOOP_CONF_USE`
- `MIN_TRIP`
- `USE_LONG_LOOP_CONF`
- `LONG_LOOP_CONF_PRODUCT_THRESH`
- `USE_CTR_INIT`
- `USE_CTR_OVERRIDE_THRESH`
- invalidation booleans
- `REQUIRE_FIRST_BRANCH_MISPRED_FOR_ALLOC`

This is the fastest way to see whether the current design is under-using or over-using the loop predictor.

### Batch B: safe capacity tuning
Then vary:

- `LOOP_LOGSETS`
- `LOOP_TAGW`
- `LOOP_ITERW`
- `LOOP_AGEW`
- `ALLOC_AGE_INIT`

This addresses aliasing and persistence without immediately exploding match logic.

### Batch C: higher-risk geometry
Only after A/B:

- `LOOP_WAYS`
- `LOGG`
- `GHIST`
- maybe `LOOP_CONFW`

### Batch D: structural changes
Only if everything above saturates:

- `LOGN`
- `LINEINST_`
- `FORCE_EXTRA_CYCLE_ON_LOOKUP`

---

## 10. Decision rule for the final answer

The agent should report the best configuration using this exact logic:

1. Gather all candidates with valid runs.
2. Keep only those with:
   - `predict2 LatMax < 2.0`
   - `reuse_predict2 LatMax < 2.0`
3. Compute per-candidate aggregate `P2_cond_acc` over the selected trace set.
4. Prefer the candidate with highest mean `P2_cond_acc`.
5. If the top two are within the confidence interval or within `1e-4` to `3e-4` absolute accuracy of each other, break ties by:
   - lower `predict2 LatMax`
   - lower `reuse_predict2 LatMax`
   - lower dynamic energy / EPI
   - lower area/transistor count

This avoids overfitting to tiny measured differences that may not be statistically meaningful.

---

## 11. What the agent should output after a sweep

For every run, store one row with:

- predictor template string
- every top-level parameter value
- compile success/failure
- trace name or experiment set name
- `P1_cond_acc`
- `P2_cond_acc`
- `predict2 LatAvg`
- `predict2 LatMax`
- `reuse_predict2 LatAvg`
- `reuse_predict2 LatMax`
- energy metric(s)
- transistors / SRAM bits if available
- feasible flag

For final reporting, include:

- best feasible config
- baseline config
- absolute improvement in `P2_cond_acc`
- relative improvement in misprediction rate if possible
- latency margin to the 2-cycle limit
- whether gains are consistent across traces or concentrated in a few traces

---

## 12. Strong practical advice

1. **Do not optimize on `gcc_test_trace.gz` alone.** Use it only for quick filtering.
2. **Do not chase tiny gains** like `+0.00002` in `P2_cond_acc` unless they are repeatable across many traces.
3. **Do not expand `LOOP_WAYS` first.** It is accuracy-friendly but directly risky for P2 latency.
4. **Start with behavior knobs before capacity knobs.** Your current design already meets latency, so smarter use of the loop predictor is likely the cheapest way to improve accuracy.
5. **Treat `FORCE_EXTRA_CYCLE_ON_LOOKUP` carefully.** It is easy to misuse as an apparent shortcut.
6. **Keep a feasible frontier.** At all times, maintain the set of points that are legal and under 2 cycles, then compare only within that frontier.

---

## 13. Minimal command pattern to automate

CBP-NG supports compiling a specific predictor instantiation with `-DPREDICTOR="..."` and running either one trace or many traces.

The agent should automate something like:

```bash
make clean
make cbp-profile CXXFLAGS='-DPREDICTOR="gshareN_loop<...>"'
make cbp-profile-analyze
make cbp-profile-acc
```

or whatever exact local Makefile wrapper exists in the working tree.

For many traces, use the repo’s multi-trace flow or an equivalent loop around the binary; the README documents `run_all` for this purpose.

---

## 14. Best first knobs for this specific predictor

Given the current measured point is already feasible but close to the 2-cycle ceiling, the first knobs I would sweep are:

1. `LOOP_CONF_USE`
2. `USE_CTR_OVERRIDE_THRESH`
3. `USE_CTR_INIT`
4. `MIN_TRIP`
5. `INVALIDATE_ON_WRONG_CONF`
6. `INVALIDATE_ON_EXCEEDED_ITER`
7. `INVALIDATE_ON_BAD_EXIT_MATCH`
8. `INVALIDATE_ON_TOO_SHORT`
9. `REQUIRE_FIRST_BRANCH_MISPRED_FOR_ALLOC`
10. `LOOP_LOGSETS`
11. `LOOP_TAGW`
12. `LONG_LOOP_CONF_PRODUCT_THRESH`
13. `LOOP_ITERW`
14. `ALLOC_AGE_INIT`

Only after that would I test:

- `LOOP_WAYS`
- `LOGG`
- `GHIST`
- `LOGN`

---

## 15. Summary instruction to the agent

Your job is to maximize **`P2_cond_acc`** for `gshareN_loop` under the hard requirement that both **`predict2 LatMax`** and **`reuse_predict2 LatMax`** remain **strictly below 2 cycles**.

Use a staged constrained search:

- fast latency screen
- medium multi-trace accuracy screen
- surrogate-model-guided refinement
- final full-set evaluation on the official training traces

Prioritize behavioral/meta-policy parameters before risky structural expansion. Use statistical aggregation across traces, not one-trace wins, to choose the final configuration. The CBP-NG repo explicitly provides the simulator interface, compilation flow, multi-trace running flow, and recommends the official 168 training traces for development.