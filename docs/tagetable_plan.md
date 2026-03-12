# TageTable Design Plan

## Overview

`TageTable` is a single tagged history table used as a component inside the
parameterizable TAGE predictor. It handles:

- SRAM storage of `(tag | counters | optional-hysteresis | optional-U)` per entry
- Ahead-pipeline banking for 1-cycle-ahead block prediction
- U-bit management: storage, probabilistic decay, and epoch-based reset
- A uniform `prime → resolveBlock → fetch → commit` interface that is
  transparent regardless of whether the table is banked/ahead or not

---

## 1. Template Parameters

```cpp
template <
    // ── Core geometry ────────────────────────────────────────────────────
    u64  TABLE_SIZE      = 64,       // number of RAM rows
    u64  TABLE_HIST      = 64,       // history length this table folds (informational)
    u64  TAG_WIDTH       = 7,        // tag bits per entry
    u64  CTR_WIDTH       = 3,        // prediction counter bits per slot
    u64  U_WIDTH         = 2,        // usefulness counter bits
    u64  PRED_BLK_SIZE   = 8,        // instruction slots per fetch block

    // ── Hysteresis ───────────────────────────────────────────────────────
    bool USE_HYS         = false,    // enable shared hysteresis bits
    u64  HYS_WIDTH       = 1,        // hysteresis bits per entry (open item A)

    // ── U-bit storage ────────────────────────────────────────────────────
    UStorage U_STORAGE   = UStorage::SRAM,
    //   SRAM: U packed inside RAM entry alongside tag+ctrs. Only valid for BANKS==1.
    //   FF:   U in separate arr<reg<U_WIDTH>, TABLE_SIZE>. Required for BANKS>1 (D1).

    UDecay   U_DECAY     = UDecay::RAND,
    //   NONE: no automatic decay; predictor manages U explicitly via writeU()
    //   RAND: probabilistic decay via std::rand() -- zero HARCOM modeled cost
    //   LFSR: probabilistic decay via hardware LFSR reg -- non-zero HARCOM cost

    u64      DECAY_CTR   = 1024,     // decay fires with probability 1/DECAY_CTR

    // ── Epoch reset ──────────────────────────────────────────────────────
    bool         USE_EPOCHS     = false,
    EpochScope   EPOCH_SCOPE    = EpochScope::GLOBAL,
    //   GLOBAL:    Tage drives epochs by calling epochAction() externally
    //   PER_TABLE: TageTable manages its own epoch counter internally

    EpochAction  EPOCH_ACTION   = EpochAction::ZERO,
    //   ZERO:      u_ff[:] = 0 (hard reset, all U cleared)
    //   SOFT_FLIP: flip epoch_bit register; U reads XOR'd with epoch_bit (U_WIDTH==1 only)
    //   PHASE:     alternating MSB/LSB clear across epoch boundaries (U_WIDTH==2 natural fit)
    //   DECREMENT: all entries decremented by 1 (expensive -- see open item B)
    //   HALF:      all entries right-shifted by 1 (expensive -- see open item B)

    EpochTrigger EPOCH_TRIGGER  = EpochTrigger::ALLOC_FAILURE,
    //   Used only when EPOCH_SCOPE == PER_TABLE. Governs what event increments
    //   the internal epoch counter:
    //   ALLOC_FAILURE: fires when allocation found no U==0 candidate
    //   MISPRED:       fires on misprediction
    //   ALLOCATION:    fires on every successful new-entry allocation
    //   CYCLE:         fires on every update_cycle() call

    u64          EPOCH_CTR_BITS = 8,   // width of PER_TABLE epoch counter
                                        // epoch period = 2^EPOCH_CTR_BITS trigger events

    // ── Ahead pipeline / banking ─────────────────────────────────────────
    bool IS_AHEAD        = false,    // enable 1-cycle-ahead pipeline staging
    u64  BANKS           = 1,        // number of outcome banks; must be power of 2
                                     // BANKS==1: unbanked (SEQ1 ahead or non-ahead)
                                     // BANKS>1:  banked ahead (BANKED mode)

    bool USE_LANE_SCRAMBLE = false,  // XOR slot_idx with XL bits of block PC
                                     // distributes counter aliasing across slots

    // ── Pipeline ─────────────────────────────────────────────────────────
    bool EN_N_BLK_RD     = false     // legacy fo2() pipeline stage on fetch return
                                     // (open item D: interaction with IS_AHEAD)
>
class TageTable { ... };
```

---

## 2. Enumerations

```cpp
enum class UStorage   { SRAM, FF };
enum class UDecay     { NONE, RAND, LFSR };
enum class EpochScope { GLOBAL, PER_TABLE };

enum class EpochAction {
    ZERO,       // hard reset: all U entries set to 0
    SOFT_FLIP,  // cheap: flip epoch_bit (U_WIDTH==1 only)
    PHASE,      // two-phase alternating MSB/LSB clear (U_WIDTH==2 natural fit)
    DECREMENT,  // expensive multi-cycle sweep (see open item B)
    HALF        // expensive multi-cycle sweep (see open item B)
};

enum class EpochTrigger {
    ALLOC_FAILURE,  // all allocation candidates had U > 0
    MISPRED,        // any misprediction
    ALLOCATION,     // successful new entry allocation
    CYCLE           // every update_cycle() call
};
```

---

## 3. Static Assertions

```cpp
// Banking constraints
static_assert(std::has_single_bit(BANKS), "BANKS must be power of 2");
static_assert(BANKS == 1 || IS_AHEAD,     "Banking requires IS_AHEAD");

// U storage constraints
static_assert(BANKS == 1 || U_STORAGE == UStorage::FF,
    "Banked tables must use FF U storage: U is shared per-index across banks (D1)");
static_assert(U_STORAGE == UStorage::SRAM || BANKS == 1 || !IS_AHEAD,
    "SRAM U storage only valid for unbanked tables");

// Epoch constraints
static_assert(U_DECAY != UDecay::NONE || USE_EPOCHS,
    "Without decay, epochs must be enabled to prevent U bits getting permanently stuck");
static_assert(EPOCH_ACTION != EpochAction::SOFT_FLIP || U_WIDTH == 1,
    "SOFT_FLIP epoch action only valid for 1-bit U");
static_assert(EPOCH_ACTION != EpochAction::PHASE || U_WIDTH == 2,
    "PHASE epoch action natural fit for 2-bit U only");

// Hysteresis constraints
static_assert(!USE_HYS || HYS_WIDTH >= 1);
```

---

## 4. Derived Constants

```cpp
static constexpr u64 IDX_BITS   = clog2(TABLE_SIZE);
static constexpr u64 LOGBANKS   = clog2(BANKS);           // 0 when BANKS==1
static constexpr u64 LOGLANES   = clog2(PRED_BLK_SIZE);

// Per-bank sub-entry bit layout (LSB to MSB):
//   [ TAG_WIDTH | HYS_WIDTH*(USE_HYS) | CTR_WIDTH*PRED_BLK_SIZE | U_WIDTH*(U_STORAGE==SRAM) ]
static constexpr u64 BITS_PER_ENTRY =
    TAG_WIDTH
    + (USE_HYS ? HYS_WIDTH : 0)
    + PRED_BLK_SIZE * CTR_WIDTH
    + (U_STORAGE == UStorage::SRAM ? U_WIDTH : 0);

// Full RAM word width:
//   Unbanked (BANKS==1): BITS_PER_ENTRY
//   Banked   (BANKS>1):  BANKS * BITS_PER_ENTRY  (WIDE_ENTRY, Approach A)
static constexpr u64 RAM_WORD_BITS = BANKS * BITS_PER_ENTRY;
```

**Entry layout note**: Counter slots are packed `ctr[0]` at LSB through `ctr[PRED_BLK_SIZE-1]`
at MSB. Tag is at the top of the sub-entry. U (when SRAM mode) is at the very bottom.
For banked entries, sub-entries are concatenated: `[bank_0 | bank_1 | ... | bank_BANKS-1]`.

---

## 5. Hardware State

```cpp
// ── Core RAM ─────────────────────────────────────────────────────────────
//   WIDE_ENTRY: one row stores all BANKS sub-entries side by side.
//   Single read returns entire row; bank is selected post-read in resolveBlock().
hcm::ram<val<RAM_WORD_BITS>, TABLE_SIZE> table_ram;

// ── Per-cycle pipeline registers ─────────────────────────────────────────
reg<IDX_BITS>                      idx_reg;   // index latched at prime()
reg<TAG_WIDTH>                     tag_reg;   // tag latched at prime()
reg<1>                             hit;       // set at resolveBlock()
reg<U_WIDTH>                       u_reg;     // cached U for current entry (always present)
arr<reg<CTR_WIDTH>, PRED_BLK_SIZE> ctr_regs;  // filled at resolveBlock()

// Hysteresis cache (only when USE_HYS)
std::conditional_t<USE_HYS,
    reg<HYS_WIDTH>, _detail::Empty>           hys_reg;

// ── Ahead staging pipeline (only when IS_AHEAD) ───────────────────────────
//   staged[0]: just read this cycle (current block's prime())
//   staged[1]: previous cycle's read (ready for resolveBlock() of current block)
std::conditional_t<IS_AHEAD,
    arr<reg<RAM_WORD_BITS>, 2>,
    _detail::Empty>                            staged;

// ── FF U-bit backing store (only when U_STORAGE == FF) ───────────────────
std::conditional_t<U_STORAGE == UStorage::FF,
    arr<reg<U_WIDTH>, TABLE_SIZE>,
    _detail::Empty>                            u_ff;

// ── LFSR decay register (only when U_DECAY == LFSR) ──────────────────────
//   Hardware shift register; has non-zero HARCOM flip-flop + logic cost.
//   (Open item C: alternative is a single global LFSR in Tage passed into commit())
std::conditional_t<U_DECAY == UDecay::LFSR,
    reg<clog2(DECAY_CTR)>,
    _detail::Empty>                            lfsr;

// ── Epoch hardware ────────────────────────────────────────────────────────
//   Internal epoch counter (only when USE_EPOCHS && PER_TABLE)
std::conditional_t<USE_EPOCHS && EPOCH_SCOPE == EpochScope::PER_TABLE,
    reg<EPOCH_CTR_BITS>,
    _detail::Empty>                            epoch_ctr;

//   Soft-flip epoch bit (only when EPOCH_ACTION == SOFT_FLIP)
std::conditional_t<USE_EPOCHS && EPOCH_ACTION == EpochAction::SOFT_FLIP,
    reg<1>,
    _detail::Empty>                            epoch_bit;

//   Phase tracker for alternating MSB/LSB clear (only when EPOCH_ACTION == PHASE)
std::conditional_t<USE_EPOCHS && EPOCH_ACTION == EpochAction::PHASE,
    reg<1>,
    _detail::Empty>                            phase;

// ── Configurable decay threshold (runtime-adjustable) ────────────────────
reg<clog2(DECAY_CTR)>                          decay_threshold;
```

---

## 6. Public Interface

The interface is uniform for all table configurations. Tage always calls the
same methods in the same order. Differences (ahead vs non-ahead, banked vs
unbanked) are handled entirely inside TageTable via `if constexpr`.

### 6.1 Call Order Per Block

```
new_block():
    prime(idx, tag)                        ← all tables

predict1() / first instruction:
    resolveBlock(bank_sel, XL)             ← all tables, once per block

predict1() / reuse_predict1() per slot:
    fetch(slot_idx)                        ← all tables, per instruction
    [writeCounter(new_val, slot_idx)]      ← optional, if Tage updates ctr_regs for P1

update_cycle():
    [writeCounter(new_val, slot_idx)]      ← for provider table counter update
    [writeU(new_u)]                        ← explicit U management (when U_DECAY==NONE)
    [writeHysteresis(new_hys)]             ← if USE_HYS
    commit(bank_id, tag, use_regs,
           new_entry, mispredict,
           alloc_failed)                   ← all tables

[epochAction()]                            ← GLOBAL scope only, called by Tage when
                                             its epoch counter fires
```

### 6.2 Method Signatures

```cpp
// ── Stage 1: RAM read ────────────────────────────────────────────────────
// Called at new_block() for ALL tables.
//
// Non-ahead: triggers RAM read immediately into staging registers.
// Ahead:     triggers RAM read into staged[0]; shifts staged[0] → staged[1].
//
// idx:  pre-computed RAM index (Tage computes this; for ahead tables it is
//       derived from the predecessor block's PC + history, Approach A)
// tag:  pre-computed tag (for ahead tables this is the successor block's tag,
//       checked later at resolveBlock() time)
void prime(val<IDX_BITS>  idx,
           val<TAG_WIDTH> tag);

// ── Stage 2a: Bank selection and tag check ───────────────────────────────
// Called once per block, at predict1() time, after path is resolved.
//
// Non-ahead (BANKS==1): applies XL lane scramble mapping if USE_LANE_SCRAMBLE;
//                       fills ctr_regs and hys_reg from staged registers;
//                       checks tag → sets hit.
// Ahead banked:         selects sub-entry from staged[1] at bank_sel;
//                       checks tag → sets hit; fills ctr_regs, u_reg, hys_reg.
//
// bank_sel: path XOR XB where XB = low LOGBANKS bits of predecessor block PC.
//           For unbanked tables (BANKS==1), pass 0.
// XL:       low LOGLANES bits of current block PC for lane scrambling.
//           Ignored when USE_LANE_SCRAMBLE == false.
void resolveBlock(val<LOGBANKS> bank_sel,
                  val<LOGLANES> XL);

// ── Stage 2b: Counter fetch ──────────────────────────────────────────────
// Called per instruction after resolveBlock(). Pure ctr_regs lookup; no RAM
// access. Cheap for both ahead and non-ahead tables.
val<CTR_WIDTH> fetch(val<LOGLANES> slot_idx);

// ── Register write helpers (called before commit()) ──────────────────────
// Update ctr_regs in-place. Tage computes new counter value externally
// (e.g. via update_ctr()) and writes it back before commit().
void writeCounter(val<CTR_WIDTH> new_val, val<LOGLANES> slot_idx);

// Explicitly set u_reg. Required when U_DECAY == NONE (Tage owns U management).
// For RAND/LFSR modes, Tage may still call this to increment U (e.g. when
// provider was correct and altpred disagrees).
void writeU(val<U_WIDTH> new_u);

// Update hysteresis register. Only meaningful when USE_HYS == true.
void writeHysteresis(val<HYS_WIDTH> new_hys);

// ── Update: write to RAM ─────────────────────────────────────────────────
// Called at update_cycle().
//
// bank_id:     which bank sub-entry to update (0 for unbanked tables).
// tag:         tag to write (may differ from tag_reg on allocation).
// use_regs:    1 = build entry from {tag, ctr_regs, hys_reg, u_reg} (normal update)
//              0 = write new_entry verbatim (allocation of fresh entry)
// new_entry:   full sub-entry bits, used only when use_regs==0.
//              Layout must match BITS_PER_ENTRY. Tage constructs this with
//              the appropriate initial counter value, tag, and U=0.
// mispredict:  forwarded to PER_TABLE epoch counter if EPOCH_TRIGGER==MISPRED.
// alloc_failed:forwarded to PER_TABLE epoch counter if EPOCH_TRIGGER==ALLOC_FAILURE.
//
// Internally:
//   - Applies U decay (RAND/LFSR) on tag miss before writing u_reg back.
//   - For banked tables: performs WIDE_ENTRY read-modify-write using staged[1]
//     to preserve unchanged bank sub-entries.
//   - For PER_TABLE epochs: increments epoch_ctr and fires epochAction() if saturated.
void commit(val<LOGBANKS>       bank_id,
            val<TAG_WIDTH>      tag,
            val<1>              use_regs,
            val<BITS_PER_ENTRY> new_entry,
            val<1>              mispredict,
            val<1>              alloc_failed);

// ── Epoch control ────────────────────────────────────────────────────────
// For GLOBAL scope: called externally by Tage when its epoch counter fires.
// For PER_TABLE scope: called internally from commit(); calling externally is a no-op.
//
// Behavior determined by EPOCH_ACTION:
//   ZERO:      u_ff.reset() or zero U bits in SRAM entries
//   SOFT_FLIP: epoch_bit = ~epoch_bit  (1-bit U only; cheap single FF flip)
//   PHASE:     clear MSB or LSB of all u_ff entries per current phase; flip phase
//   DECREMENT: multi-cycle background sweep (see open item B)
//   HALF:      multi-cycle background sweep (see open item B)
void epochAction();

// ── Runtime configuration ────────────────────────────────────────────────
// Adjust decay firing probability at runtime (threshold compared against LFSR/rand).
void setDecayThreshold(val<clog2(DECAY_CTR)> new_threshold);

// ── Accessors ────────────────────────────────────────────────────────────
val<1>      getHit();           // tag match result from resolveBlock()
val<U_WIDTH> getUsefulness();   // u_reg from resolveBlock()

// Hysteresis accessor (only meaningful when USE_HYS == true)
val<HYS_WIDTH> getHysteresis(); // hys_reg from resolveBlock()

// ── Hardware geometry ────────────────────────────────────────────────────
// Returns critical-path latency through this table in picoseconds.
// Extracted from table_ram SRAM model + MUX overhead.
static constexpr u64 critical_path_ps();
```

---

## 7. Behavioral Notes

### 7.1 Non-ahead tables (`IS_AHEAD == false`)

`prime()` triggers the RAM read and immediately fills staging registers.
`resolveBlock()` checks the tag and fills `ctr_regs` from the just-read data.
`fetch()` is a pure register read. All three happen within the same block's
processing (same pipeline stage). No extra flip-flop overhead beyond what
`ctr_regs`, `tag_reg`, and `idx_reg` already represent.

### 7.2 Ahead unbanked tables (`IS_AHEAD == true`, `BANKS == 1`)

SEQ1 mode. `prime()` at `new_block(B0)` reads the RAM entry that will be
used to predict B1 (fall-through). The read is pipelined through `staged[]`
and arrives at `resolveBlock()` during B1's processing. Tag check happens
then. On a `pref_valid` miss (tag mismatch or base address mismatch), the
predictor falls back per `AHEAD_MISS_POLICY`.

### 7.3 Ahead banked tables (`IS_AHEAD == true`, `BANKS > 1`)

BANKED mode. Predecessor-indexed storage (Approach A): the RAM row index is
derived from the predecessor block B0's PC + history. All BANKS outcome
sub-entries are stored in the same row (WIDE_ENTRY). `prime()` reads the
entire row for block B0 into `staged[0]`. At `resolveBlock()` for block B1:
`staged[1]` (B0's read, shifted over one cycle) is available; `bank_sel =
path XOR XB` selects the sub-entry corresponding to B0's actual path out.
Tag check uses the stored tag against the successor block's expected tag.

`commit()` performs a read-modify-write: `staged[1]` provides the unchanged
bank sub-entries; only `bank_id`'s sub-entry is replaced.

### 7.4 U-bit decay in `commit()`

Decay fires when `hit == 0` (tag miss on the selected bank):

- `RAND`: `val<DECAY_BITS>{static_cast<unsigned>(std::rand())} == 0`
- `LFSR`: advance LFSR register, fire if `lfsr == 0`
- `NONE`: no automatic decay; u_reg reflects only what Tage wrote via `writeU()`

For `SOFT_FLIP` epoch action, all U reads and writes must XOR with `epoch_bit`.
This is handled internally in `getUsefulness()`, `writeU()`, and `commit()`.

### 7.5 Lane scrambling

When `USE_LANE_SCRAMBLE == true`, `resolveBlock(bank_sel, XL)` maps
`ctr_regs[i XOR XL]` into the natural slot `i` during the fill step.
`fetch(slot_idx)` always returns `ctr_regs[slot_idx]` (after the XL mapping
has been applied). `writeCounter(val, slot_idx)` writes to the natural
`ctr_regs[slot_idx]` directly; the XL mapping is purely a read-path concern.

---

## 8. Open Items

### A — Hysteresis storage location

**Options:**
1. Packed inside the RAM entry alongside tag + counters (current plan).
   - Read at `resolveBlock()` time into `hys_reg`; written at `commit()` time.
   - Hysteresis participates in WIDE_ENTRY layout; goes through the same
     ahead pipeline as counters.
2. Separate `UPDATE_ONLY` zone RAM of depth `TABLE_SIZE`.
   - Not read during prediction; only read/written at update time.
   - Saves read-path width; does not participate in `staged[]`.
   - Consistent with how `predictors/tage.hpp` structures hysteresis.

**Recommendation:** Option 2 (separate UPDATE_ONLY zone) is more consistent
with the reference implementation and cheaper on the predict critical path.
The `BITS_PER_ENTRY` formula above should exclude `HYS_WIDTH` if Option 2 is
chosen, with a separate `hcm::ram<val<HYS_WIDTH>, TABLE_SIZE>` in an
`UPDATE_ONLY` zone. Requires a decision before implementation.

### B — DECREMENT / HALF epoch actions

Multi-cycle sweeps over all `TABLE_SIZE` entries are expensive in hardware and
not straightforward to model in HARCOM. Two paths forward:

1. Skip these actions for now; keep `ZERO`, `SOFT_FLIP`, and `PHASE`.
2. Model as a background FSM: a `reg<IDX_BITS> sweep_ptr` and a `reg<1>
   sweeping` flag; one entry decremented/halved per `update_cycle()` call.

**Recommendation:** Skip in initial implementation. Add `DECREMENT` and `HALF`
only if experimental results suggest they outperform `PHASE`.

### C — LFSR scope: per-table vs global

**Per-table LFSR** (current plan): `reg<DECAY_BITS>` inside each TageTable.
Self-contained. Multiple tables have independent random streams. Cost:
`DECAY_BITS × NUM_TABLES` flip-flops.

**Global LFSR** (alternative): single `reg<DECAY_BITS>` in Tage, output
passed into `commit()` as an extra argument `val<DECAY_BITS> rng_bits`.
Cost: `DECAY_BITS` flip-flops total. Tables with same `DECAY_CTR` share one
stream (correlated decay across tables; acceptable for rate-control purposes).

**Recommendation:** Decide based on area budget. Suggest starting with
per-table for simplicity; add global LFSR pass-through as an opt-in parameter
`bool USE_GLOBAL_LFSR = false` with a corresponding `commit()` overload.

### D — `EN_N_BLK_RD` interaction with `IS_AHEAD`

`EN_N_BLK_RD` currently adds a `.fo2()` pipeline register on the `fetch()`
return value in the reference implementation. With `IS_AHEAD == true`, the
data has already been pipelined through `staged[]`. It is unclear whether the
`.fo2()` stage is still needed, redundant, or conflicting.

**Action required:** Trace the HARCOM fanout model through both paths and
determine if `EN_N_BLK_RD` should be: (a) ignored when `IS_AHEAD`, (b)
applied on top of the staged data, or (c) replaced entirely by the staged
pipeline. Add a `static_assert` to enforce a consistent choice.

---

## 9. Future Scaffolding (Not Implemented, Parameterized for Extension)

The following items are not in scope for the initial implementation but should
be reflected in parameter names or enum placeholders to avoid later interface
breaks.

| Item | Where | Notes |
|---|---|---|
| `HashFunc` enum (`FOLD_XOR`, ...) | Tage-level | TageTable receives pre-computed idx/tag; hash is external. Scaffold enum in Tage template params. |
| `USE_PATH_HIST`, `PATH_HIST_LEN` | Tage-level | Affects index computation in Tage only; TageTable is unaffected. |
| `BASE_PRED` enum | Tage-level | P0 base predictor type; no impact on TageTable interface. |
| `USE_ALT_ON_NA` | Tage-level | Meta counter preferring altpred on newly-allocated entries. |
| Multi-port RAM | TageTable | `READ_PORTS`, `WRITE_PORTS` params; default 1. Scaffold as reserved template params. |
| `USE_GLOBAL_LFSR` | TageTable / Tage | Per open item C above. |
| Statistical correction (SC) | Tage-level | Separate tables; no TageTable interface change. |
| Loop predictor | Tage-level | Separate component; no TageTable interface change. |
| Tag compression | TageTable | `bool TAG_COMPRESS = false`; affects BITS_PER_ENTRY but not interface. Reserve param slot. |
| Per-entry metadata | TageTable | `META_BITS = 0`; extra bits per entry for future use (e.g. confidence). Reserve field in entry layout. |
| `SPLIT_KSPACE` | Tage-level | Separate kernel history; explicitly out of scope but noted. |

---

## 10. Design Decisions Record

| Decision | Choice | Rationale |
|---|---|---|
| D1: U storage with banking | Shared per-index (FF array) | U tracks table-level usefulness, not per-path. Simpler; no U bits in wide entry. |
| D2: Banked storage impl | WIDE_ENTRY (single wide RAM row) | All banks share same predecessor-derived row index (Approach A); one RAM read returns all banks. Consistent with `gshareN_ahead`. |
| D3: Lane scrambling | Parameterized via `USE_LANE_SCRAMBLE` | Matches `gshareN_ahead` XL scrambling; zero cost when disabled. |
| D4: Interface transparency | Uniform `prime/resolveBlock/fetch/commit` | Tage calls identical methods on all tables. Ahead vs non-ahead distinction is internal to TageTable and in compile-time index selection in Tage (`AHEAD_TABLE_MASK`). |
| Banking index | Predecessor-indexed (Approach A) | Single row index for all banks; avoids BANKS separate index computations per table per block. |
| `USE_LSFR_DECAY` removal | Replaced by `UStorage` + `UDecay` | Original parameter conflated two independent concerns (storage location and decay mechanism). |
| Epoch scope | Both GLOBAL and PER_TABLE supported | Some tables benefit from coordinated reset; others benefit from autonomous rate control. |
| Epoch triggered | GLOBAL counter in Tage; internal counter for PER_TABLE | GLOBAL tables have no hardware epoch cost in TageTable; PER_TABLE tables are self-contained. |
