# HARCOM Floorplan & Wiring Optimization Guide

How to use HARCOM's layout features (`zone`, `region`, `connect`, `distribute`)
to reduce timing and energy in branch predictors.

## Background: Why Layout Matters

HARCOM models physical placement of RAMs on silicon. Each RAM occupies area,
and every val/reg has a **location** tied to the nearest RAM. When an operation
reads an input from a different location, that input must travel across a wire
whose delay and energy grow with Manhattan distance.

Without explicit layout control, HARCOM uses defaults:
- Operations happen at the location of the **latest-arriving input**
- Each read of a value runs a separate wire from source to destination
- All RAMs land in a single default region

These defaults are often suboptimal: wide values travel long distances, fanout
creates redundant wires, and unrelated RAMs are placed near each other.

## Technique 1: `zone` — Separate Predict-Path from Update-Path RAMs

### What it does
A `zone` groups RAMs together in the floorplan. RAMs in the same zone are
placed adjacent. RAMs in different zones may be placed far apart.

### The pattern
Most predictors have two classes of RAM:
1. **Predict-path RAMs** — read every cycle in predict1/predict2 (latency-critical)
2. **Update-only RAMs** — read/written only on mispredicts in update_cycle (not latency-critical)

Separating them into different zones lets the floorplanner pack predict-path
RAMs tightly (minimizing wire delay on the critical path), without being
forced to interleave them with large update-only RAMs.

### Example (from reference gshare)
```cpp
ram<val<1>,(1<<index1_bits)> table1_hi[LINEINST];  // predict-path
ram<val<1>,(1<<index2_bits)> table2_hi[LINEINST];  // predict-path

zone UPDATE_ONLY;
ram<val<1>,(1<<index1_bits)> table1_lo[LINEINST];  // hysteresis (update only)
ram<val<1>,(1<<index2_bits)> table2_lo[LINEINST];  // hysteresis (update only)
```

### Where to apply in our codebase

**gshareN_ahead**: `ctr_hi` is read every cycle (predict-path), `ctr_lo[LANES]`
is read only on mispredict (update-only). Adding `zone UPDATE_ONLY;` before
`ctr_lo` declaration separates them.

**Custom TAGE** (`Tage.hpp` / `TageTable.hpp`): The tag and prediction RAMs are
on the predict path. The hysteresis (`ghyst`) and useful-bit (`ubit`) RAMs are
update-only. Similarly, the bimodal hysteresis (`bhyst`) is update-only.

**SCOverrider**: The bias and GEHL weight RAMs are read in predict1 (via
`prefetch`), but writes happen in `train` gated by `extra_cycle`. These are
read-write-same-cycle RAMs (`rwram`), so the zone benefit is less clear — but
if any tables are read-only during prediction, they should be separated.

### Impact
- **Energy**: Shorter wires between predict-path RAMs → less wire energy per cycle
- **Timing**: Predict-path RAMs packed closer → lower wire delay on critical path
- **Area**: No change (same total RAM, just different placement)


## Technique 2: `connect` — Move Operations to Reduce Wire Width

### What it does
`val.connect(ram)` forces the val's location to be at `ram`'s location. The
val itself is unchanged — only where the operation is computed changes.

### The problem it solves
Consider reading a 64-bit value from a small fast RAM (location A) and ANDing
it with a 1-bit value from a large slow RAM (location B). By default, the AND
happens at B (latest-arriving input), requiring 64 wires from A→B. If instead
we `connect` the 1-bit value to location A, the AND happens at A, and only 1
wire runs from B→A. The result is 64× less wire energy.

### The pattern
```
Wide value at location X  ──(many wires)──→  operation at location Y
Narrow value at location Y ──(few wires)──→  operation at location Y

// Fix: move narrow value to wide value's location
narrow_val = narrow_val.connect(wide_ram);
// Now operation happens at X, only narrow wire from Y→X
```

### Example (from reference perceptron)
```cpp
// global_history bit (1-bit, at GHR location) needs to reach
// weight table (WEIGHT_BITS wide, at per-branch RAM location).
// Without connect: WEIGHT_BITS wires from weight RAM to GHR location.
// With connect: 1 wire from GHR to weight RAM location.
auto ghbit = global_history[j - 1].connect(site[i].weights);
ghbit.fanout(hard<WEIGHT_BITS>{});
return site[i].branch_weights[j] ^
       ghbit.replicate(hard<WEIGHT_BITS>{}).concat();
```

### Where to apply in our codebase

**Custom TAGE predict2**: Tag comparison reads a wide tag from a RAM and
compares with a hash of PC+history. The hash is narrow relative to the tag
width. If the hash is computed far from the tag RAM, `connect` the hash to the
tag RAM to avoid transporting the full tag width across the chip.

**gshareN_ahead predict1**: `global_history` (GHIST bits wide, at GHR location)
is folded into an index that addresses `ctr_hi` RAM. The fold result should be
`connect`ed to `ctr_hi` if the GHR is physically far from the prediction RAM.

**SCOverrider**: The TAGE prediction sum (a few bits) is combined with SC
weights (read from multiple GEHL RAMs at different locations). `connect` the
TAGE sum to the first GEHL table's location to minimize transport of the wider
weight values.

### Impact
- **Energy**: Dramatic reduction when wide values are kept local and only narrow values travel
- **Timing**: Modest improvement — wire delay is proportional to distance, which may decrease
- **When NOT to use**: Don't connect if both values are similarly wide; the default (operation at latest-arriving input) already minimizes delay


## Technique 3: `distribute` — Tree Distribution for Multi-RAM Writes

### What it does
`val.distribute(ram_array)` creates copies of a val at each RAM in the array,
using a distribution tree. Without it, writing the same value to N RAMs sends
N separate point-to-point wires from the val's location.

### The problem it solves
When the same index or data value is written to multiple RAMs (e.g., same
index across N TAGE tables, or same direction bit to N hysteresis RAMs), the
default point-to-point wiring creates N long wires. A distribution tree
creates a binary fan-out tree, reducing total wire length from O(N×D) to
O(N×log(D)) where D is the distance.

### Example (from reference hashed_perceptron)
```cpp
// index2 is used to address NTABLES weight tables, each at a different location.
// distribute creates copies at each table's location.
auto dindex2 = index2[i].distribute(wtable[i]);
```

### Where to apply in our codebase

**Custom TAGE update_cycle**: The same `branch_pc` and `ghist` values are used
to compute indices for all NUMG tables. Each table's tag/pred/hyst/ubit RAM is
at a different location. The index should be `distribute`d across the table
RAMs rather than sent point-to-point.

**gshareN_ahead update_cycle**: `concat(index[1],path)` is sent to all LANES
`ctr_lo[i]` RAMs. Since LANES=8, distributing this address saves 7 redundant
long wires.

**SCOverrider train**: The branch PC hash is used to address all GEHL tables.
Distributing the hash across the GEHL RAMs reduces wiring.

### Impact
- **Energy**: Significant reduction for high fan-out writes (N ≥ 4)
- **Timing**: Marginal — distribution tree adds buffer delay but removes long-wire delay
- **When NOT to use**: When N is small (2-3 RAMs) the benefit is minimal


## Technique 4: `region` — Isolate Independent Subsystems

### What it does
A `region` groups RAMs and assigns them to an isolated portion of the
floorplan. RAMs in different regions are placed independently. Regions also
control where C++ literal vals are located (via `region.enter()`).

### The pattern
Use regions when a predictor has independent subsystems that rarely exchange
data. Each subsystem gets its own region, so its RAMs are packed together
without interference from the other subsystem's RAMs.

### Where to apply in our codebase

**TAGE + SCOverrider + LoopPredictor**: These three subsystems interact only at
the prediction mux level (a few bits). Their RAMs should be in separate regions:
```cpp
region TAGE_REGION;
// ... TAGE table RAMs ...
region SC_REGION;
// ... SC bias and GEHL RAMs ...
region LOOP_REGION;
// ... Loop predictor RAMs ...
```

This prevents the floorplanner from interleaving SC RAMs with TAGE RAMs,
which would increase wire distances within each subsystem.

**gshareN_ahead**: Only has one subsystem — regions don't help here.

### Impact
- **Energy**: Reduces intra-subsystem wire lengths
- **Timing**: Can improve if predict-path RAMs within a region are packed tighter
- **When NOT to use**: When subsystems share many RAMs or exchange wide values frequently


## Summary: What to Apply Where

| Predictor | `zone` | `connect` | `distribute` | `region` |
|-----------|--------|-----------|--------------|----------|
| gshareN_ahead | ctr_lo → UPDATE_ONLY | GHR→ctr_hi | index→ctr_lo[LANES] | N/A |
| Custom TAGE | ghyst,ubit,bhyst → UPDATE_ONLY | tag hash→tag RAM | index→N tables | Separate from SC/Loop |
| SCOverrider | N/A (all rwram) | TAGE sum→GEHL[0] | PC hash→GEHL tables | Separate from TAGE |
| LoopPredictor | N/A (single rwram) | N/A | N/A | Separate from TAGE/SC |

## Measurement

Use these compiler flags to quantify the impact:
- `-DFREE_WIRING` — removes all wiring costs. Compare EPI with/without to see total wiring overhead.
- `-DFREE_FANOUT` — removes all fanout delays. Compare P1 latency with/without.
- `-DFLOORPLAN` — generates `floorplan.gv` for visualization with `dot`.

The EPI (energy per instruction) metric in the contest score includes wiring
energy. Reducing it directly improves the score without changing prediction
accuracy.
