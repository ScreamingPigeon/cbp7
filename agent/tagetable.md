# TageTable Interface Documentation

## Overview

`TageTable` is a hardware-modeled TAGE (TAgged GEometric history length) prediction table component written using the HARCOM library. It represents a single table in a TAGE predictor with support for:
- Block-based predictions (multiple predictions per fetch block)
- Tag matching for table hits
- Useful counters (u-bits) for allocation policies
- Per-instruction prediction counters

Each `TageTable` instance manages one table in the TAGE hierarchy, indexed by a combination of PC and global history with a specific history length.

---

## Template Parameters

```cpp
template <u64 TABLE_SIZE = 64,      // Number of entries in the table
          u64 TABLE_HIST = 64,      // History length for this table
          u64 TAG_WIDTH = 7,        // Tag width in bits
          u64 CTR_WIDTH = 3,        // Prediction counter width
          u64 U_WIDTH = 2,          // Useful counter width
          u64 PRED_BLK_SIZE = 8,    // Fetch block size (predictions per entry)
          u64 DECAY_CTR = 1024,     // Decay period for useful counters
          bool EN_N_BLK_RD = true>  // Enable newRead return value
class TageTable;
```

### Parameter Descriptions

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| `TABLE_SIZE` | Number of entries in the RAM | 16-1024 (power of 2) |
| `TABLE_HIST` | History length used for indexing this table | Geometric series: 8, 16, 32, 64, 128, etc. |
| `TAG_WIDTH` | Number of bits for tag comparison | 7-12 bits |
| `CTR_WIDTH` | Width of each prediction counter | 2-4 bits (3 is typical) |
| `U_WIDTH` | Width of useful counter | 1-3 bits (2 is typical) |
| `PRED_BLK_SIZE` | Number of predictions per entry (fetch block size) | 4-16 instructions |
| `DECAY_CTR` | Number of accesses before decaying u-bits | 256-4096 |
| `EN_N_BLK_RD` | Enable return value from `newRead()` | true (standard), false (void) |

### Computed Constants

```cpp
static constexpr u64 BITS_PER_INST = CTR_WIDTH;
static constexpr u64 BITS_PER_ENTRY = TAG_WIDTH + PRED_BLK_SIZE * BITS_PER_INST + U_WIDTH;
```

**Entry Layout** (MSB to LSB):
```
| TAG (TAG_WIDTH) | CTR[PRED_BLK_SIZE-1] | ... | CTR[1] | CTR[0] (CTR_WIDTH) | U (U_WIDTH) |
     MSB                                                                            LSB
```

Built via `concat(u_reg, ctr_regs.concat(), tag_reg)` where:
- `ctr_regs.concat()` places ctr_regs[0] at LSB, ctr_regs[PRED_BLK_SIZE-1] at MSB
- Total width: TAG_WIDTH + (PRED_BLK_SIZE × CTR_WIDTH) + U_WIDTH bits

---

## Public Interface

### 1. `newRead()` — Start New Block Access

```cpp
auto newRead(val<clog2(TABLE_SIZE)> idx,
             val<TAG_WIDTH> tag,
             val<clog2(PRED_BLK_SIZE)> slot_idx)
```

**Purpose**: Initiates a new table lookup for a fetch block.

**Parameters**:
- `idx`: Index into the table (computed from PC ⊕ folded_history)
- `tag`: Tag to compare against stored tag (computed from PC ⊕ history)
- `slot_idx`: Position within the fetch block (0 to PRED_BLK_SIZE-1)

**Behavior**:
1. Reads the entry at `table_ram[idx]`
2. Splits the entry using `split<U_WIDTH, PRED_BLK_SIZE * CTR_WIDTH, TAG_WIDTH>(entry)`:
   - Extracts u-bit (LSB)
   - Extracts combined counter bits (middle)
   - Extracts tag (MSB)
3. Stores values in registers:
   - `idx_reg = idx`: Caches the index
   - `tag_reg = tag`: Caches the provided tag (for future `updateBlock` calls)
   - `u_reg = u_bits`: Caches the extracted useful counter
   - `hit = (tag_bits == tag)`: Stores tag match result
4. Extracts individual counters using `static_loop<PRED_BLK_SIZE>`:
   - Creates a val array `ctrs` for SRAM read results
   - Iterates over each counter position with compile-time index `I`
   - Shifts combined counter bits right by `I * CTR_WIDTH`
   - Auto-truncates to `CTR_WIDTH` bits and stores in:
     - `ctrs[I]`: Val array (for return from SRAM read)
     - `ctr_regs[I]`: Register array (for cached reuse)
   - Layout: index 0 at LSB, index PRED_BLK_SIZE-1 at MSB
5. Returns the selected counter **from the SRAM read result** with FO2 optimization

**Implementation Details**:
```cpp
// Read and split entry
val<BITS_PER_ENTRY> entry = table_ram.read(idx);
auto [u_bits, ctr_bits_combined, tag_bits] =
    split<U_WIDTH, PRED_BLK_SIZE * CTR_WIDTH, TAG_WIDTH>(entry);

// Store state in registers
idx_reg = idx;
tag_reg = tag;
u_reg = u_bits;

// Extract individual counters into both val array and register array
arr<val<CTR_WIDTH>, PRED_BLK_SIZE> ctrs;
static_loop<PRED_BLK_SIZE>([&]<int I>() {
  constexpr u64 shift_amt = I * CTR_WIDTH;
  val<PRED_BLK_SIZE * CTR_WIDTH> shifted = ctr_bits_combined >> hard<shift_amt>{};
  ctrs[I] = shifted;     // Store in val array for return
  ctr_regs[I] = shifted; // Cache in register for reuse
});

// Check tag match
hit = (tag_bits == tag);

// Return from SRAM read result (not register) with fo2
return ctrs.select(slot_idx).fo2();
```

**Return Value**: `val<CTR_WIDTH>` — The prediction counter for the requested slot from SRAM read, with FO2 fanout optimization

**Usage Pattern**: Call once per fetch block when accessing a new PC.

**RAM Cost**: 1 RAM read from `table_ram`

**HARCOM Notes**:
- Uses `static_loop` for compile-time iteration (avoids runtime loop overhead)
- `hard<shift_amt>{}` creates zero-cost unsigned right shifts
- Auto-truncation when assigning wider val to narrower reg/val
- Maintains both val array (ctrs) and register array (ctr_regs) for dual purpose:
  - Val array: Return from SRAM read result (transient values)
  - Register array: Cache for `reuseRead()` calls (persistent state)
- `select()` on array creates hardware mux for slot selection
- `.fo2()` declares fanout of 2, optimizing read delay for the return value
- Returns from SRAM read (vals), not cached registers, ensuring fresh data

---

### 2. `reuseRead()` — Reuse Cached Entry

```cpp
auto reuseRead(val<clog2(PRED_BLK_SIZE)> slot_idx)
```

**Purpose**: Retrieves a prediction from the currently cached entry without accessing RAM.

**Parameters**:
- `slot_idx`: Position within the fetch block

**Behavior**:
1. Returns `ctr_regs[slot_idx]` (no RAM access)
2. Uses the entry cached by the previous `newRead()` call

**Usage Pattern**: Call for subsequent instructions in the same fetch block (PC+1, PC+2, etc.)

**RAM Cost**: 0 (uses cached registers)

**Hardware Benefit**: Amortizes the cost of RAM access across multiple predictions per block.

---

### 3. `writeReg()` — Update Cached Counter

```cpp
auto writeReg(val<BITS_PER_INST> new_data,
              val<clog2(PRED_BLK_SIZE)> slot_idx)
```

**Purpose**: Updates a prediction counter in the cached entry (register state only, not yet written to RAM).

**Parameters**:
- `new_data`: New prediction counter value (typically saturating counter update)
- `slot_idx`: Which counter to update within the block

**Behavior**:
1. Writes `new_data` to `ctr_regs[slot_idx]`
2. Does NOT write to RAM yet (deferred until `updateBlock()`)

**Usage Pattern**: Call when updating prediction counters after branch resolution.

**RAM Cost**: 0 (register write only)

---

### 4. `updateBlock()` — Commit Entry to RAM

```cpp
void updateBlock(val<1> use_regs,
                 val<TAG_WIDTH> tag,
                 val<BITS_PER_ENTRY> new_entry)
```

**Purpose**: Writes an entry back to the table RAM, either from cached registers or a new entry.

**Parameters**:
- `use_regs`: Select source (1 = use cached registers, 0 = use `new_entry`)
- `tag`: Tag to write (used when `use_regs == 1`, merged with cached counters)
- `new_entry`: Complete entry to write (used when `use_regs == 0`)

**Behavior**:
```cpp
if (use_regs == 1):
    entry = concat(tag_reg, ctr_regs, u_reg)  // Build from cached state
else:
    entry = new_entry                          // Use provided entry
table_ram.write(idx_reg, entry)
```

**Usage Pattern**:
- **Update existing entry**: `use_regs=1`, `tag` is current tag
  - Used when updating counters for a table hit
  - Writes cached `ctr_regs` (modified by `writeReg()`) back to RAM
- **Allocate new entry**: `use_regs=0`, `new_entry` is fresh entry
  - Used when allocating a new tag on misprediction
  - Initializes tag, counters, and u-bit

**RAM Cost**: 1 RAM write to `table_ram`

---

### 5. `getUsefulness()` — Read Useful Counter

```cpp
auto getUsefulness()
```

**Purpose**: Returns the cached useful counter value from the most recent `newRead()`.

**Return Value**: `reg<U_WIDTH>` — The u-bit value for the currently cached entry

**Usage Pattern**: Call after `newRead()` to check u-bit for allocation decisions

**Example**:
```cpp
table->newRead(idx, tag, slot);
val<U_WIDTH> u = table->getUsefulness();
val<1> can_allocate = (u == 0);  // Only allocate if u-bit is 0
```

---

### 6. `getHit()` — Read Tag Match Status

```cpp
auto getHit()
```

**Purpose**: Returns whether the most recent `newRead()` resulted in a tag match.

**Return Value**: `reg<1>` — 1 if tag matched, 0 if miss

**Usage Pattern**: Call after `newRead()` to determine if table provides a valid prediction

**Example**:
```cpp
table->newRead(idx, tag, slot);
val<1> hit = table->getHit();
val<1> pred = select(hit, provider_pred, base_pred);  // Use table if hit
```

---

### 7. `setThreshold()` — Configure Decay Threshold (TODO)

```cpp
auto setThreshold()
```

**Purpose**: Allows the predictor to set a dynamic threshold for probabilistic u-bit decay.

**Status**: Not yet implemented (placeholder for LFSR-based decay feature)

**Planned Usage**:
```cpp
// Adjust decay aggressiveness based on allocation failure rate
val<4> new_threshold = compute_threshold(alloc_failure_rate);
table->setThreshold(new_threshold);
```

See TageTable.hpp:97-100 for TODO notes.

---

### 8. `decrementU()` — Probabilistic U-bit Decay (TODO)

```cpp
auto decrementU()
```

**Purpose**: Implements probabilistic u-bit decay based on LFSR random number generation.

**Status**: Not yet implemented (private method, invoked within `updateBlock()`)

**Planned Behavior**:
1. Generate random bits using LFSR
2. Compare against threshold
3. If below threshold, decrement u-bit on tag miss
4. Otherwise, leave u-bit unchanged

**Planned Implementation** (TageTable.hpp:117-122):
```cpp
auto decrementU() {
  // Generate random bits from LFSR
  val<4> rand = lfsr & hard<0xF>{};

  // Probabilistic decay: decrement if rand < threshold
  val<1> should_decay = (rand < decay_threshold);
  val<U_WIDTH> decremented = select(u_reg == 0, u_reg, u_reg - 1);
  return select(should_decay, decremented, u_reg);
}
```

See the "Hardware RNG for Probabilistic Decay" section below for LFSR details.

---

## Internal State

### Registers (Persistent across cycles)

```cpp
reg<1> hit;                              // Tag match flag
reg<TAG_WIDTH> tag_reg;                  // Cached tag
reg<clog2(TABLE_SIZE)> idx_reg;          // Cached index
reg<U_WIDTH> u_reg;                      // Cached useful counter
arr<reg<CTR_WIDTH>, PRED_BLK_SIZE> ctr_regs;  // Cached prediction counters
```

**Purpose**: Cache the most recently accessed entry to enable:
- Block-level reuse (multiple predictions from one RAM read)
- Read-modify-write pattern for counter updates

### RAM

```cpp
hcm::ram<val<BITS_PER_ENTRY>, TABLE_SIZE> table_ram;
```

**Storage**: Each entry contains:
- Tag (TAG_WIDTH bits)
- Prediction counters (PRED_BLK_SIZE × CTR_WIDTH bits)
- Useful counter (U_WIDTH bits)

---

## Usage Workflow

### Typical Prediction Flow (Tag Hit)

```cpp
// 1. New fetch block at PC
val<idx_bits> idx = compute_index(pc, history);
val<TAG_WIDTH> tag = compute_tag(pc, history);
val<slot_bits> slot0 = 0;

// 2. Read table and get first prediction
auto pred0 = tage_table.newRead(idx, tag, slot0);

// 3. Reuse cached entry for subsequent instructions (PC+1, PC+2, ...)
auto pred1 = tage_table.reuseRead(val<slot_bits>{1});
auto pred2 = tage_table.reuseRead(val<slot_bits>{2});

// 4. On branch resolution, update counter
val<CTR_WIDTH> updated_ctr = saturating_update(pred0, actual_taken);
tage_table.writeReg(updated_ctr, slot0);

// 5. Commit to RAM (use cached registers)
tage_table.updateBlock(val<1>{1}, tag, val<BITS_PER_ENTRY>{0});
```

### Allocation Flow (Tag Miss → Allocate New Entry)

```cpp
// 1. Compute new entry components
val<TAG_WIDTH> new_tag = compute_tag(pc, history);
val<CTR_WIDTH> init_ctr = val<CTR_WIDTH>{taken ? 4 : 3};  // Weakly taken/not-taken
val<U_WIDTH> init_u = 0;

// 2. Build complete entry
arr<val<CTR_WIDTH>, PRED_BLK_SIZE> init_ctrs = {init_ctr, 0, 0, 0, ...};
val<BITS_PER_ENTRY> new_entry = concat(init_u, init_ctrs.concat(), new_tag);

// 3. Write new entry (bypass cached registers)
tage_table.updateBlock(val<1>{0}, new_tag, new_entry);
```

---

## HARCOM-Specific Considerations

### RAM Access Constraints

- **One access per cycle**: Each `TageTable` can perform at most one RAM operation (read OR write) per cycle
- **Cycle advancement**: The predictor's superuser must call `panel.next_cycle()` between consecutive RAM accesses to the same table
- **Floorplan**: `panel.make_floorplan()` must be called after declaring all `TageTable` instances

### Block-Level Optimization

The block-based design amortizes RAM costs:
- **Without blocking**: 8 instructions → 8 RAM reads
- **With blocking**: 8 instructions → 1 RAM read + 7 register reads

This significantly reduces energy and improves performance for fetch blocks.

### Opaque Values

All HARCOM values (`val`, `reg`) are opaque:
- Cannot use `if (hit)` — use `select(hit, x, y)` instead
- Cannot read counter values as integers — use `select()` for logic
- All operations on HARCOM types accumulate hardware cost

### Fanout Management

Consider using `fanout()` or `fo1()` on frequently-read values:
```cpp
auto pred = tage_table.newRead(idx, tag, slot0);
pred.fanout(hard<4>{});  // If pred will be read 4 times
```

---

## Integration with TAGE Predictor

### Multi-Table Hierarchy

A typical TAGE predictor has 5-14 tables with geometric history lengths:

```cpp
TageTable<128, 8,   7, 3, 2, 8, 1024> table0;   // Shortest history
TageTable<64,  16,  7, 3, 2, 8, 1024> table1;
TageTable<32,  32,  7, 3, 2, 8, 1024> table2;
TageTable<16,  64,  8, 3, 2, 8, 1024> table3;
TageTable<16,  128, 9, 3, 2, 8, 1024> table4;   // Longest history
```

### Provider/Alternate Selection

1. Query all tables in parallel (or sequentially across cycles)
2. Select the **longest-history hitting table** as the provider
3. Select the **next-longest hitting table** as the alternate
4. Use provider's prediction if confident (high u-bit), else alternate

### Allocation Policy

On misprediction:
1. Identify non-hitting tables with history longer than provider
2. Allocate 1-3 entries (starting from shortest eligible table)
3. Use `updateBlock(use_regs=0, ...)` to write fresh entries

### Useful Counter Management

- Increment u-bit when provider correct and alternate wrong
- Decrement u-bit when provider wrong
- Periodically reset u-bits (every DECAY_CTR accesses) to free entries

---

## Performance Notes

### Area Cost
- Dominated by `table_ram` (TABLE_SIZE × BITS_PER_ENTRY bits of SRAM)
- Registers add negligible area (BITS_PER_ENTRY bits)

### Energy Cost
- RAM reads: ~50-200 fJ depending on table size
- RAM writes: ~1.5× read cost
- Register operations: <1 fJ

### Latency
- RAM access: ~150-300 ps depending on floorplan
- Register reads: ~10-20 ps
- Critical path: `newRead()` → tag comparison → mux selection

### Optimization Tips
1. **Balance table sizes**: Smaller tables for long histories (fewer unique indices)
2. **Tune PRED_BLK_SIZE**: Larger blocks amortize RAM cost but increase entry size
3. **Minimize tag width**: Use just enough bits to avoid excessive aliasing
4. **Use regions**: Group related tables in the same HARCOM region for lower wiring cost

---

## Example: Complete Table Interaction

```cpp
// Initialization (in predictor constructor)
TageTable<64, 32, 7, 3, 2, 8, 1024> t2;
panel.make_floorplan();

// Prediction phase
val<6> idx = (pc ^ folded_history).fo1();
val<7> tag = hash_tag(pc, history);
val<3> slot = 0;

val<1> hit;
val<3> pred;
execute_if(access_this_table, [&]() {
    pred = t2.newRead(idx, tag, slot);
    hit = t2.get_hit();  // Hypothetical accessor
});

// Update phase (next cycle)
panel.next_cycle();

val<3> new_ctr = select(taken,
    select(pred == val<3>::maxval, pred, pred + 1),  // Increment
    select(pred == 0, pred, pred - 1));              // Decrement

execute_if(hit, [&]() {
    t2.writeReg(new_ctr, slot);
    t2.updateBlock(val<1>{1}, tag, val<BITS_PER_ENTRY>{0});
});
```

---

## Hardware RNG for Probabilistic Decay (Planned)

### Overview

TageTable.hpp:19-22 contains a TODO for implementing parameterized LFSR-based dynamic threshold probabilistic decay. This feature will adaptively decay u-bits based on allocation failure rate.

### Planned Architecture

```cpp
template <u64 TABLE_SIZE = 64, ..., u64 LFSR_WIDTH = 16>
class TageTable {
private:
  // LFSR for random number generation
  reg<LFSR_WIDTH> lfsr;

  // Dynamic decay threshold (adjusted by predictor)
  reg<4> decay_threshold;  // 0-15 range

  // Allocation tracking
  reg<10> alloc_attempts;
  reg<10> alloc_failures;
};
```

### LFSR Implementation

**16-bit maximal-length LFSR** (recommended):

```cpp
void tick_lfsr() {
  // Taps at positions 16,14,13,11 for maximal period
  val<1> feedback = lfsr[15] ^ lfsr[13] ^ lfsr[12] ^ lfsr[10];
  lfsr = concat(feedback, lfsr >> hard<1>{});
}
```

**Hardware cost**: 16 flip-flops + 3 XOR gates (~10 fJ/tick)

### Probabilistic Decay Logic

```cpp
auto decrementU() {
  // Extract random bits from LFSR
  val<4> rand_bits = lfsr;  // Bottom 4 bits

  // Compare against dynamic threshold
  val<1> should_decay = (rand_bits < decay_threshold);

  // Conditionally decrement u-bit
  val<U_WIDTH> decremented = select(u_reg == 0, u_reg, u_reg - 1);
  return select(should_decay, decremented, u_reg);
}
```

### Adaptive Threshold Adjustment

The predictor adjusts `decay_threshold` based on allocation success rate:

```cpp
// In Custom predictor's update_cycle()
val<10> failure_rate = alloc_failures << hard<10>{} / alloc_attempts;
val<1> high_failure = (failure_rate > hard<512>{});  // >50%

// Increase threshold if many failures (more aggressive decay)
// Decrease threshold if few failures (less aggressive decay)
val<4> new_thresh = select(high_failure,
                           decay_threshold + 1,
                           select(decay_threshold == 0,
                                  val<4>{0},
                                  decay_threshold - 1));

table->setThreshold(new_thresh);
```

### Decay Probability Mapping

| Threshold | Probability | Use Case |
|-----------|-------------|----------|
| 0 | 1/16 (6.25%) | Low allocation pressure |
| 4 | 5/16 (31%) | Moderate pressure |
| 8 | 9/16 (56%) | High pressure |
| 15 | 16/16 (100%) | Critical pressure (always decay) |

### Integration Points

1. **newRead()**: Tick LFSR on every access
   ```cpp
   tick_lfsr();  // Advance RNG state
   ```

2. **updateBlock()**: Apply probabilistic decay on tag misses
   ```cpp
   val<1> is_miss = !hit;
   val<U_WIDTH> new_u = select(is_miss, decrementU(), u_reg);
   ```

3. **setThreshold()**: Allow predictor to adjust aggressiveness
   ```cpp
   void setThreshold(val<4> new_threshold) {
     decay_threshold = new_threshold;
   }
   ```

### Benefits

1. **Adaptive allocation**: Frees entries when tables are full
2. **Low hardware cost**: Single 16-bit LFSR shared across all accesses
3. **Tunable**: Predictor controls decay rate based on observed behavior
4. **Deterministic**: Reproducible for debugging (fixed LFSR seed)

### Alternative: PC-based Pseudo-Random

Zero-cost alternative using existing values:

```cpp
auto decrementU() {
  // XOR existing state for pseudo-randomness
  val<4> pseudo_rand = (idx_reg ^ tag_reg);
  val<1> should_decay = (pseudo_rand < decay_threshold);
  // ... same decay logic ...
}
```

**Tradeoff**: No area cost, but less random (correlated with address patterns).

---

## Future Extensions

Potential enhancements to the TageTable interface:

1. ✅ **Accessor methods**: `getHit()`, `getUsefulness()` implemented
2. **Conditional update**: `updateBlock()` wrapped in `execute_if()` logic
3. ⏳ **LFSR-based decay**: Probabilistic u-bit decay (TODO, see above)
4. **Statistical correction**: Additional state for SC predictor integration
5. **Loop predictor support**: Iteration counter fields
6. **Per-table regions**: HARCOM region hints for better floorplan

---

## Related Files

- **`predictors/custom/TageTable.hpp`** — Implementation
- **`predictors/custom.hpp`** — Custom predictor using TageTable
- **`harcom.hpp`** — HARCOM library (~6K lines)
- **`claude/harcom.md`** — HARCOM language reference
- **`docs/harcom.pdf`** — Full HARCOM specification
