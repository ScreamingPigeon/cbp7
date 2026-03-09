#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"

using namespace hcm;

// Constexpr helper (kept for backward compatibility with Tage.hpp)
constexpr size_t clog2(std::uint64_t n) {
  int r = 0;
  std::uint64_t v = 1;
  while (v < n) {
    v <<= 1;
    ++r;
  }
  return r;
}

// Empty type for conditional member elimination
struct EmptyMember {
  EmptyMember() = default;
  EmptyMember(const char *) {}
};

// Default u-bit reset functor for FF mode
// Mode 0: reset to zero
// Mode 1: right shift by 1
// Mode 2: saturating decrement by 1
struct DefaultResetFn {
  static constexpr u64 MODE_BITS = 2;

  template <u64 W> static val<W> apply(val<W> u, val<MODE_BITS> mode) {
    auto shifted = u >> hard<1>{};
    auto decremented =
        select(u == hard<0>{}, val<W>{0}, val<W>{u - 1});
    return select(mode == hard<0>{}, val<W>{0},
                  select(mode == hard<1>{}, shifted, decremented));
  }
};

// ============================================================================
// TageTable: parameterized TAGE tagged table component
//
// One instance per table in the TAGE hierarchy. The predictor holds N
// instances and handles all cross-table logic (match priority, allocation,
// u-bit epoch triggers). The table is a banked RAM wrapper with optional
// register caching — no fold, no pipeline, no timing awareness.
//
// The predictor owns: fold registers, index/tag computation, ahead pipeline
// buffering, path/XB management, and all cross-table logic. It passes
// pre-computed index, tag, and bank index to the table.
//
// See agent/tagetable_plan.md for full design rationale.
// ============================================================================

template <
    u64 TABLE_SIZE = 1024,   // entries per bank (power of 2, >= 2)
    u64 TABLE_HIST = 64,     // history length (not used by table; predictor
                             //   uses for fold computation)
    u64 TAG_WIDTH = 11,      // tag width in bits
    u64 CTR_WIDTH = 3,       // prediction counter width (embeds hysteresis)
    u64 U_WIDTH = 1,         // useful counter width
    u64 N = 4,               // max branches per cycle = predictions per block
    u64 NUM_BANKS = 1,       // branch-slot banks; BPB = N / NUM_BANKS
    bool USE_AHEAD = false,  // 1-ahead pipelining; doubles physical banks
    bool SHARED_TAG = true,  // share one tag across branch-slot banks
    bool SHARED_U = true,    // share one u-bit across branch-slot banks
    bool U_STOR_FF = true,   // true = u-bits in flip-flops, false = u-bits in SRAM
    u64 DECAY_CTR = 1024,    // probabilistic decay period (U_STOR_FF=false only)
    typename ResetFn = DefaultResetFn, // u-bit reset functor (U_STOR_FF=true only)
    bool USE_FF_CACHE = false // cache SRAM reads in FFs for block reuse
    >
class TageTable {

  // ======== Computed Constants ========

public:
  static constexpr u64 bpb = N / NUM_BANKS;

  static constexpr u64 table_size = TABLE_SIZE;
  static constexpr u64 table_hist = TABLE_HIST;
  static constexpr u64 tag_width = TAG_WIDTH;
  static constexpr u64 ctr_width = CTR_WIDTH;
  static constexpr u64 u_width = U_WIDTH;
  static constexpr u64 n_branches = N;
  static constexpr u64 num_banks = NUM_BANKS;
  static constexpr bool use_ahead = USE_AHEAD;
  static constexpr bool shared_tag = SHARED_TAG;
  static constexpr bool shared_u = SHARED_U;
  static constexpr bool use_ff_cache = USE_FF_CACHE;

private:
  static constexpr u64 BPB = bpb;
  static constexpr u64 AHEAD_FACTOR = USE_AHEAD ? 2 : 1;
  static constexpr u64 PHYS_BANKS = NUM_BANKS * AHEAD_FACTOR;
  static constexpr u64 IDX_BITS = clog2(TABLE_SIZE);

  // Per-bank SRAM entry composition
  // Counters are always present. Tag and U are optionally packed per-bank.
  static constexpr u64 CTR_BITS = BPB * CTR_WIDTH;
  static constexpr u64 BANK_TAG_BITS = SHARED_TAG ? 0 : TAG_WIDTH;
  static constexpr u64 BANK_U_BITS =
      (!SHARED_U && !U_STOR_FF) ? U_WIDTH : 0;
  static constexpr u64 BANK_ENTRY_WIDTH =
      CTR_BITS + BANK_TAG_BITS + BANK_U_BITS;

  // Result register counts (depends on sharing mode)
  static constexpr u64 TAG_REG_COUNT = SHARED_TAG ? 1 : NUM_BANKS;
  static constexpr u64 U_REG_COUNT = SHARED_U ? 1 : NUM_BANKS;
  static constexpr u64 HIT_REG_COUNT = SHARED_TAG ? 1 : NUM_BANKS;

  // Counter cache depth per bank
  static constexpr u64 CACHED_CTRS = USE_FF_CACHE ? BPB : 1;

  // U-bit FF array count
  static constexpr u64 U_FF_ARRAYS = SHARED_U ? 1 : NUM_BANKS;

  // ======== Static Constraints ========

  static_assert(TABLE_SIZE >= 2,
                "TABLE_SIZE must be at least 2");
  static_assert(std::has_single_bit(TABLE_SIZE),
                "TABLE_SIZE must be power of 2");
  static_assert(N > 0,
                "Must predict at least one branch");
  static_assert(NUM_BANKS > 0,
                "Must have at least one bank");
  static_assert(N % NUM_BANKS == 0,
                "N must be divisible by NUM_BANKS");
  static_assert(!SHARED_TAG || SHARED_U,
                "Per-bank U (SHARED_U=false) requires per-bank tags "
                "(SHARED_TAG=false)");
  static_assert(!USE_FF_CACHE || BPB > 1,
                "FF caching requires BPB > 1 (multiple branches per bank)");
  static_assert(CTR_WIDTH > 0,
                "Counter width must be positive");
  static_assert(TAG_WIDTH > 0,
                "Tag width must be positive");
  static_assert(U_WIDTH > 0,
                "U-bit width must be positive");

  // ======== SRAM Storage ========

  // Per-bank entry RAM
  // Entry layout (MSB to LSB) depends on SHARED_TAG and U_STOR_FF:
  //   SHARED_TAG=true,  SHARED_U=true  or FF: [ctr[BPB-1]|...|ctr[0]]
  //   SHARED_TAG=true,  SHARED_U=false, SRAM: [ctr[BPB-1]|...|ctr[0]|u]
  //   SHARED_TAG=false, SHARED_U=true  or FF: [tag|ctr[BPB-1]|...|ctr[0]]
  //   SHARED_TAG=false, SHARED_U=false, SRAM: [tag|ctr[BPB-1]|...|ctr[0]|u]
  //
  // PHYS_BANKS = NUM_BANKS * AHEAD_FACTOR
  // Ahead stage s uses banks [s*NUM_BANKS .. (s+1)*NUM_BANKS - 1]
  hcm::ram<val<BANK_ENTRY_WIDTH>, TABLE_SIZE> bank_ram[PHYS_BANKS]{"bank"};

  // Shared tag RAM: one per ahead stage (only when SHARED_TAG=true)
  // When SHARED_TAG=false, tags are packed in bank_ram entries.
  std::conditional_t<SHARED_TAG,
      hcm::ram<val<TAG_WIDTH>, TABLE_SIZE>,
      EmptyMember> shared_tag_ram[AHEAD_FACTOR]{"stag"};

  // Shared U-bit RAM: one per ahead stage
  // Only when SHARED_U=true AND u-bits in SRAM.
  // When SHARED_U=false, u packed in bank_ram. When U_STOR_FF, u in FFs.
  std::conditional_t<(SHARED_U && !U_STOR_FF),
      hcm::ram<val<U_WIDTH>, TABLE_SIZE>,
      EmptyMember> shared_u_ram[AHEAD_FACTOR]{"su"};

  // ======== Flip-Flop Storage ========

  // U-bit flip-flops (only when U_STOR_FF=true)
  // SHARED_U=true:  1 array of TABLE_SIZE entries
  // SHARED_U=false: NUM_BANKS arrays of TABLE_SIZE entries each
  std::conditional_t<U_STOR_FF,
      arr<reg<U_WIDTH>, TABLE_SIZE>,
      EmptyMember> u_ff[U_FF_ARRAYS];

  // ======== Result / Cache Registers ========
  // These hold the results of the most recent read for accessor methods
  // and for write-back. When USE_FF_CACHE=true, counter regs also serve
  // as a block-level cache for reuseRead().

  // Cached index for write-back
  reg<IDX_BITS> idx_reg;

  // Tag result: stored tags from read (for accessor and write-back)
  // SHARED_TAG=true:  1 register
  // SHARED_TAG=false: NUM_BANKS registers
  reg<TAG_WIDTH> tag_reg[TAG_REG_COUNT];

  // Hit result: tag comparison results
  // SHARED_TAG=true:  1 register (all banks share hit/miss)
  // SHARED_TAG=false: NUM_BANKS registers (per-bank hit/miss)
  reg<1> hit_reg[HIT_REG_COUNT];

  // U-bit result
  // SHARED_U=true:  1 register
  // SHARED_U=false: NUM_BANKS registers
  reg<U_WIDTH> u_reg[U_REG_COUNT];

  // Counter result / cache registers
  // USE_FF_CACHE=true:  BPB counters per bank (full block cache for reuse)
  // USE_FF_CACHE=false: 1 counter per bank (current read result only)
  arr<reg<CTR_WIDTH>, CACHED_CTRS> ctr_regs[NUM_BANKS];

public:
  TageTable() = default;
  ~TageTable() = default;

  // Phase 2: read(), reuseRead(), write(), accessors
  // Phase 2: reset_u() for FF mode
};
