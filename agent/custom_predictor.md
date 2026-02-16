# Custom TAGE Predictor Architecture

## Overview

The Custom predictor is a TAGE (TAgged GEometric history length) branch predictor implementation using HARCOM for hardware cost modeling. It features:

- **Multiple TAGE tables** with geometric history lengths
- **Block-based predictions** (fetch blocks of 4-16 instructions)
- **Optional path history** for improved correlation
- **Kernel/user split history** support
- **Probabilistic u-bit decay** (planned)
- **Two-level prediction** (P1 fast, P2 accurate)

## File Structure

```
predictors/
├── custom.hpp              # Top-level Custom predictor (stub, references Tage.hpp)
└── custom/
    ├── Tage.hpp           # Main TAGE predictor implementation
    ├── TageTable.hpp      # Single TAGE table component
    └── loggers/           # Debugging/logging utilities
```

**Note**: `custom.hpp` includes `custom/Tage.hpp`, which contains the full implementation.

---

## Template Parameters

### Table-Level Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `T_TAG_WIDTH` | u64 | 7 | Tag width in bits for table matching |
| `T_U_WIDTH` | u64 | 2 | Useful counter width (allocation policy) |
| `T_CTR_WIDTH` | u64 | 3 | Prediction counter width per instruction |
| `T_HYS_WIDTH` | u64 | 2 | Hysteresis bit width (if used) |
| `T_USE_HYS` | bool | true | Enable hysteresis bits |
| `T_DECAY_CTR` | u64 | 1024 | Period for u-bit decay |

### Predictor-Level Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `PRED_BLK_SIZE` | u64 | 8 | Fetch block size (instructions per block) |
| `NUM_TABLES` | u64 | 5 | Number of TAGE tables in hierarchy |
| `T_TABLE_HIST_LEN` | array | {8,16,32,64,128} | History length for each table |
| `T_TABLE_HIST_SIZE` | array | {128,64,32,16,16} | Number of entries per table |
| `GLOBAL_HIST_LEN` | u64 | 128 | Maximum global history length |
| `USE_PATH_HIST` | bool | false | Enable path history tracking |
| `PATH_HIST_LEN` | u64 | 16 | Path history register width |
| `SPLIT_KSPACE` | bool | false | Separate kernel/user histories |
| `EN_N_BLK_RD` | bool | false | Enable multiple block reads per table |

---

## Architecture

### Component Hierarchy

```
Custom (predictor)
├── Base Predictor (TODO)
│   └── Bimodal table for T0
└── TAGE Tables (tuple of 5-14 tables)
    ├── Table 0: Short history (8 branches)
    ├── Table 1: Medium history (16 branches)
    ├── Table 2: Medium-long (32 branches)
    ├── Table 3: Long history (64 branches)
    └── Table 4: Very long (128 branches)
```

### State Registers

```cpp
// Global branch history
reg<GLOBAL_HIST_LEN> br_hist_reg;

// Conditional: Kernel history (if SPLIT_KSPACE)
std::conditional_t<SPLIT_KSPACE, reg<GLOBAL_HIST_LEN>, Empty> br_khist_reg;

// Conditional: Path history (if USE_PATH_HIST)
std::conditional_t<USE_PATH_HIST, reg<PATH_HIST_LEN>, Empty> p_hist_reg;

// Conditional: Kernel path history (if both)
std::conditional_t<USE_PATH_HIST & SPLIT_KSPACE, reg<PATH_HIST_LEN>, Empty> p_khist_reg;

// Base address for block-aligned PC
reg<BLOCK_BITS> base_addr;
```

**BLOCK_BITS** = `64 - clog2(PRED_BLK_SIZE) - 2`
- Aligns PC to block boundaries
- Example: PRED_BLK_SIZE=8 → BLOCK_BITS=59

---

## TAGE Table Tuple

The predictor uses template metaprogramming to create a heterogeneous tuple of TageTable instances:

```cpp
template <std::size_t... Is>
static auto make_tables_impl(std::index_sequence<Is...>) {
  return std::make_tuple(
      std::make_unique<TageTable<
          std::get<Is>(T_TABLE_HIST_SIZE),   // Size varies per table
          std::get<Is>(T_TABLE_HIST_LEN),    // History length varies
          T_TAG_WIDTH,                        // Shared
          T_CTR_WIDTH,                        // Shared
          T_U_WIDTH,                          // Shared
          PRED_BLK_SIZE,                      // Shared
          T_DECAY_CTR,                        // Shared
          EN_N_BLK_RD                         // Shared
      >>()...
  );
}

// Instantiate tuple with NUM_TABLES elements
decltype(make_tables_impl(std::make_index_sequence<NUM_TABLES>{})) tage_tables;
```

**Key points:**
- Each table has different size and history length
- Uses `std::make_unique` for heap allocation (HARCOM types may have complex lifetimes)
- Tuple allows accessing specific tables at compile time: `std::get<0>(tage_tables)`

---

## Prediction Pipeline

### 1. `new_block(val<64> inst_pc)`

**Purpose**: Initialize state for a new fetch block

**Current Implementation** (Tage.hpp:102-106):
```cpp
void new_block(val<64> inst_pc) {
  base_addr = 0; // TODO: Extract block-aligned PC
  // TODO: Compute tag and index for all tables
  // TODO: Call newRead() on all tables
  // TODO: Select provider and alternate tables
}
```

**Planned Implementation**:
```cpp
void new_block(val<64> inst_pc) {
  // 1. Extract block address
  base_addr = inst_pc >> hard<64 - BLOCK_BITS>{};

  // 2. Hash PC + history for each table
  auto [idx0, tag0] = compute_hash<0>(inst_pc, br_hist_reg);
  auto [idx1, tag1] = compute_hash<1>(inst_pc, br_hist_reg);
  // ... for all NUM_TABLES

  // 3. Read all tables in parallel (different RAM objects)
  val<3> ctr0 = std::get<0>(tage_tables)->newRead(idx0, tag0, slot_0);
  val<3> ctr1 = std::get<1>(tage_tables)->newRead(idx1, tag1, slot_0);
  // ... for all tables

  // 4. Determine provider (longest-history hitting table)
  val<1> hit0 = std::get<0>(tage_tables)->getHit();
  val<1> hit1 = std::get<1>(tage_tables)->getHit();
  // ... select provider and alternate
}
```

### 2. `predict1(val<64> inst_pc)` - Fast Prediction

**Purpose**: Return quick prediction (low latency, lower accuracy)

**Options**:
- Return base predictor result
- Return first hitting table
- Return cached prediction from `new_block()`

### 3. `predict2(val<64> inst_pc)` - Accurate Prediction

**Purpose**: Return high-confidence prediction (higher latency, better accuracy)

**Implementation**:
- Use full TAGE provider selection
- Check u-bit confidence
- Use alternate on low confidence (USE_ALT_ON_NA)

### 4. `reuse_predict1/2(val<64> inst_pc)` - Block Reuse

**Purpose**: Predict for subsequent instructions in the block (PC+1, PC+2, ...)

**Implementation**:
```cpp
val<1> reuse_predict1(val<64> inst_pc) {
  // Compute slot index within block
  val<slot_bits> slot_idx = inst_pc - (base_addr << hard<64-BLOCK_BITS>{});

  // Reuse cached table data
  val<3> ctr = provider_table->reuseRead(slot_idx);
  return (ctr >= hard<4>{});  // MSB = direction
}
```

### 5. `update_condbr(val<64> branch_pc, val<1> taken, val<64> next_pc)`

**Purpose**: Train predictor on branch outcome

**Implementation**:
```cpp
void update_condbr(val<64> branch_pc, val<1> taken, val<64> next_pc) {
  // 1. Update global history
  br_hist_reg = (br_hist_reg << hard<1>{}) | taken;

  // 2. Update path history (if enabled)
  if constexpr (USE_PATH_HIST) {
    val<PATH_HIST_LEN> pc_hash = branch_pc;
    p_hist_reg = (p_hist_reg << hard<4>{}) ^ pc_hash;
  }

  // 3. Update provider table counter
  val<3> new_ctr = saturating_update(old_ctr, taken);
  provider_table->writeReg(new_ctr, slot_idx);

  // 4. Update u-bit
  // Increment if provider correct and alternate wrong
  // Decrement if provider wrong

  // 5. Write back to RAM
  provider_table->updateBlock(val<1>{1}, tag, val<BITS_PER_ENTRY>{0});

  // 6. Handle allocation (if misprediction)
  // Allocate 1-3 entries in longer-history tables
}
```

### 6. `update_cycle(instruction_info &block_end_info)`

**Purpose**: End-of-block cleanup

**Implementation**:
- Periodic u-bit reset (every DECAY_CTR accesses)
- Statistical correction updates (if enabled)
- Allocation failure tracking

---

## Hash Functions

### Index Computation

```cpp
template <size_t TableIdx>
val<idx_bits> compute_index(val<64> pc, val<GLOBAL_HIST_LEN> history) {
  constexpr u64 HIST_LEN = T_TABLE_HIST_LEN[TableIdx];
  constexpr u64 TABLE_SIZE = T_TABLE_HIST_SIZE[TableIdx];
  constexpr u64 idx_bits = clog2(TABLE_SIZE);

  // Fold history to match table size
  val<idx_bits> hist_folded = fold_history(history, HIST_LEN, idx_bits);
  val<idx_bits> pc_bits = pc >> hard<2>{};  // Skip 4-byte alignment

  if constexpr (USE_PATH_HIST) {
    val<idx_bits> path_bits = p_hist_reg;
    return pc_bits ^ hist_folded ^ path_bits;
  } else {
    return pc_bits ^ hist_folded;
  }
}
```

### Tag Computation

```cpp
template <size_t TableIdx>
val<TAG_WIDTH> compute_tag(val<64> pc, val<GLOBAL_HIST_LEN> history) {
  constexpr u64 idx_bits = clog2(T_TABLE_HIST_SIZE[TableIdx]);

  // Use different PC bits than index
  val<TAG_WIDTH> pc_tag = pc >> hard<idx_bits + 2>{};
  val<TAG_WIDTH> hist_tag = fold_history(history, T_TABLE_HIST_LEN[TableIdx], TAG_WIDTH);

  if constexpr (USE_PATH_HIST) {
    val<TAG_WIDTH> path_tag = p_hist_reg >> hard<4>{};
    return pc_tag ^ hist_tag ^ path_tag;
  } else {
    return pc_tag ^ hist_tag;
  }
}
```

---

## Provider Selection Algorithm

```
For each table from longest to shortest history:
  If table hits (tag match):
    Provider = this table
    Alternate = next longest hitting table (or base predictor)
    Break

If no tables hit:
  Provider = base predictor
  Alternate = none

Prediction:
  If provider is base predictor OR provider has low confidence (u-bit):
    Use alternate prediction (if available)
  Else:
    Use provider prediction
```

---

## Allocation Policy

On misprediction:

1. **Identify candidate tables**: Non-hitting tables with history > provider history
2. **Select N tables** (typically 1-3, starting from shortest eligible)
3. **Allocate entries**:
   ```cpp
   val<TAG_WIDTH> new_tag = compute_tag<TableIdx>(pc, history);
   val<CTR_WIDTH> init_ctr = taken ? val<3>{4} : val<3>{3};  // Weak taken/not-taken
   val<U_WIDTH> init_u = 0;  // Fresh entry

   // Build new entry
   val<BITS_PER_ENTRY> new_entry = build_entry(new_tag, init_ctr, init_u);
   table->updateBlock(val<1>{0}, new_tag, new_entry);
   ```

---

## Planned Enhancements

### 1. Probabilistic U-bit Decay (TageTable.hpp:19-22)

**Goal**: Adaptively decay u-bits based on allocation failure rate

**Implementation**:
- Add LFSR-based RNG to TageTable
- Track allocation attempts/failures
- Increase decay probability when allocation failure rate is high
- Decrease when allocation succeeds

See `tagetable.md` for LFSR implementation details.

### 2. USE_ALT_ON_NA

**Goal**: Use alternate prediction for newly allocated entries

**Implementation**:
- Track entry age (0 = newly allocated)
- Force alternate prediction if provider age < threshold
- Counters to bias toward provider/alternate

### 3. Loop Predictor

**Goal**: Handle short loops with specialized predictor

**Components**:
- Loop iteration counter
- Loop target detector
- Confidence counter

### 4. Statistical Correction (SC)

**Goal**: Correct TAGE mispredictions with additional tables

**Components**:
- Multiple small tables with different hash functions
- Perceptron-style weights
- Final correction vote

---

## Hardware Cost Estimates

### Per-Table Costs (example: 128 entries, 8-instruction blocks)

| Component | Bits | Notes |
|-----------|------|-------|
| RAM (128×39 bits) | 4,992 | TAG(7) + CTR(3×8) + U(2) per entry |
| Registers (cached) | 39 | hit, tag_reg, idx_reg, u_reg, ctr_regs[8] |
| **Total per table** | **~5,031 bits** | |

### Full Predictor (5 tables)

| Component | Bits | Energy (fJ) |
|-----------|------|-------------|
| TAGE tables (5×) | ~25,000 | ~500-1000/prediction |
| History registers | 128 | ~10/update |
| Path history | 16 | ~2/update |
| Base predictor (TODO) | ~2,048 | ~50/prediction |
| **Total** | **~27,200 bits** | **~600-1100 fJ** |

---

## Usage Example

```cpp
// Instantiate in params.yaml or compile command:
Custom<
  /* T_TAG_WIDTH */ 7,
  /* T_U_WIDTH */ 2,
  /* T_CTR_WIDTH */ 3,
  /* PRED_BLK_SIZE */ 8,
  /* NUM_TABLES */ 5,
  /* T_TABLE_HIST_LEN */ {8, 16, 32, 64, 128},
  /* T_TABLE_HIST_SIZE */ {128, 64, 32, 16, 16},
  /* GLOBAL_HIST_LEN */ 128,
  /* USE_PATH_HIST */ true,
  /* PATH_HIST_LEN */ 16
> predictor;
```

---

## TODO List

### High Priority
- [ ] Implement `new_block()` hashing and table reads
- [ ] Implement provider selection logic
- [ ] Implement `predict1/2()` and `reuse_predict1/2()`
- [ ] Implement `update_condbr()` with counter updates
- [ ] Add base predictor (bimodal table)

### Medium Priority
- [ ] Implement allocation policy
- [ ] Add u-bit update logic
- [ ] Implement periodic u-bit reset
- [ ] Add LFSR for probabilistic decay

### Low Priority
- [ ] Add USE_ALT_ON_NA support
- [ ] Implement loop predictor
- [ ] Add statistical correction (SC)
- [ ] Optimize hash functions with folded_history

---

## Related Documentation

- **tagetable.md** — TageTable interface and implementation details
- **harcom.md** — HARCOM hardware modeling reference
- **CLAUDE.md** — Project overview and build system
