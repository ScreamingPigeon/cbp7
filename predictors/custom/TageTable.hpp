#pragma once

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

  // ======== Private Helpers ========

  // Unpack a bank entry into its component fields.
  // Returns counter bits. Writes tag_reg for this bank if applicable.
  // For per-bank SRAM u-bits: extracts u bits into u_val_out instead of
  // writing u_reg (caller applies decay before writing u_reg once).
  // bank_local: logical bank index within a stage [0..NUM_BANKS)
  // Unpack a bank entry into its component fields.
  // Returns counter bits. Writes tag_reg for this bank if applicable.
  // For per-bank SRAM u (!SHARED_U && !U_STOR_FF): writes u_reg directly
  // if apply_decay is false; if apply_decay is true, writes u_reg with
  // probabilistic decay applied based on hit_reg.
  val<CTR_BITS> unpack_bank_entry(val<BANK_ENTRY_WIDTH> entry,
                                  u64 bank_local) {
    if constexpr (!SHARED_TAG && !SHARED_U && !U_STOR_FF) {
      // [tag | ctrs | u]
      auto [u_bits, ctr_bits, tag_bits] =
          split<BANK_U_BITS, CTR_BITS, BANK_TAG_BITS>(entry);
      tag_reg[bank_local] = tag_bits;
      // u_reg written after tag comparison (see read())
      return ctr_bits;
    } else if constexpr (!SHARED_TAG && (SHARED_U || U_STOR_FF)) {
      // [tag | ctrs]
      auto [ctr_bits, tag_bits] = split<CTR_BITS, BANK_TAG_BITS>(entry);
      tag_reg[bank_local] = tag_bits;
      return ctr_bits;
    } else if constexpr (SHARED_TAG && !SHARED_U && !U_STOR_FF) {
      // [ctrs | u]
      auto [u_bits, ctr_bits] = split<BANK_U_BITS, CTR_BITS>(entry);
      // u_reg written after tag comparison (see read())
      return ctr_bits;
    } else {
      // [ctrs] only
      return entry;
    }
  }

  // Extract per-bank u bits from a bank entry (only for !SHARED_U && !U_STOR_FF).
  val<U_WIDTH> extract_u_from_entry(val<BANK_ENTRY_WIDTH> entry) {
    static_assert(!SHARED_U && !U_STOR_FF,
                  "extract_u_from_entry only for per-bank SRAM u");
    if constexpr (!SHARED_TAG) {
      // [tag | ctrs | u] — u is LSB
      auto [u_bits, rest] = split<BANK_U_BITS, CTR_BITS + BANK_TAG_BITS>(entry);
      return u_bits;
    } else {
      // [ctrs | u] — u is LSB
      auto [u_bits, rest] = split<BANK_U_BITS, CTR_BITS>(entry);
      return u_bits;
    }
  }

  // Extract individual counters from combined counter bits and store in regs.
  void extract_counters(val<CTR_BITS> ctr_bits, u64 bank_local) {
    static_loop<CACHED_CTRS>([&]<u64 I>() {
      constexpr u64 shift_amt = I * CTR_WIDTH;
      ctr_regs[bank_local][I] = ctr_bits >> hard<shift_amt>{};
    });
  }

  // Pack counter regs back into combined counter bits for a bank.
  val<CTR_BITS> pack_counters(u64 bank_local) {
    return ctr_regs[bank_local].concat();
  }

  // Pack a full bank entry from components.
  val<BANK_ENTRY_WIDTH> pack_bank_entry(val<CTR_BITS> ctr_bits,
                                        val<TAG_WIDTH> tag,
                                        val<U_WIDTH> u) {
    if constexpr (!SHARED_TAG && !SHARED_U && !U_STOR_FF) {
      return concat(u, ctr_bits, tag);
    } else if constexpr (!SHARED_TAG && (SHARED_U || U_STOR_FF)) {
      return concat(ctr_bits, tag);
    } else if constexpr (SHARED_TAG && !SHARED_U && !U_STOR_FF) {
      return concat(u, ctr_bits);
    } else {
      return ctr_bits;
    }
  }

  // Probabilistic u-bit decay (SRAM mode only).
  // Returns 1 with probability 1/DECAY_CTR.
  val<1> should_decay() {
    if constexpr (DECAY_CTR <= 1) {
      return val<1>{1}; // always decay
    } else {
      constexpr u64 DECAY_BITS = clog2(DECAY_CTR);
      val<DECAY_BITS> rng{static_cast<unsigned>(std::rand())};
      return rng == hard<0>{};
    }
  }

  // Write a single entry in a u-bit FF array at a dynamic val index.
  // Models a decoder driving each FF's write-enable.
  void write_u_ff(u64 arr_idx, val<IDX_BITS> index, val<U_WIDTH> u) {
    for (u64 i = 0; i < TABLE_SIZE; i++) {
      execute_if(index == val<IDX_BITS>{static_cast<unsigned>(i)},
                 [&]() { u_ff[arr_idx][i] = u; });
    }
  }

public:
  TageTable() = default;
  ~TageTable() = default;

  // ======== Read Interface ========

  // Compute index fanout: idx_reg + bank_ram reads + optional shared RAMs/FFs
  static constexpr u64 IDX_READ_FANOUT =
      1 // idx_reg
      + NUM_BANKS // bank_ram reads
      + (SHARED_TAG ? 1 : 0) // shared_tag_ram
      + (SHARED_U && U_STOR_FF ? 1 : 0) // u_ff select (shared)
      + (SHARED_U && !U_STOR_FF ? 1 : 0) // shared_u_ram
      + (!SHARED_U && U_STOR_FF ? NUM_BANKS : 0); // u_ff selects (per-bank)

  // Compare tag fanout: 1 comparison per hit_reg written
  static constexpr u64 TAG_CMP_FANOUT = SHARED_TAG ? 1 : NUM_BANKS;

  // Read all banks at index, compare tags, cache results.
  // stage: ahead stage index (0 when USE_AHEAD=false)
  void read(val<IDX_BITS> index, val<TAG_WIDTH> compare_tag,
            u64 stage = 0) {
    // Fanout declarations for multi-use vals
    index.fanout(hard<IDX_READ_FANOUT>{});
    if constexpr (TAG_CMP_FANOUT > 1) {
      compare_tag.fanout(hard<TAG_CMP_FANOUT>{});
    }

    idx_reg = index;

    // Read shared tag if applicable
    if constexpr (SHARED_TAG) {
      val<TAG_WIDTH> stored_tag = shared_tag_ram[stage].read(index);
      stored_tag.fanout(hard<2>{}); // tag_reg write + comparison
      tag_reg[0] = stored_tag;
      hit_reg[0] = (stored_tag == compare_tag);
    }

    // Read shared u from FFs if applicable (no decay for FFs)
    if constexpr (SHARED_U && U_STOR_FF) {
      u_reg[0] = u_ff[0].select(index);
    }

    // Read each bank: unpack entries, extract counters and tags.
    u64 bank_base = stage * NUM_BANKS;
    for (u64 b = 0; b < NUM_BANKS; b++) {
      val<BANK_ENTRY_WIDTH> entry = bank_ram[bank_base + b].read(index);

      val<CTR_BITS> ctr_bits = unpack_bank_entry(entry, b);
      extract_counters(ctr_bits, b);

      // Per-bank tag comparison (when !SHARED_TAG)
      if constexpr (!SHARED_TAG) {
        hit_reg[b] = (tag_reg[b] == compare_tag);
      }

      // Per-bank u from FFs (when !SHARED_U && U_STOR_FF)
      if constexpr (!SHARED_U && U_STOR_FF) {
        u_reg[b] = u_ff[b].select(index);
      }

      // Per-bank SRAM u: extract raw value (decay deferred to update path)
      if constexpr (!SHARED_U && !U_STOR_FF) {
        u_reg[b] = extract_u_from_entry(entry);
      }
    }

    // Shared SRAM u-bit: read raw value (decay deferred to update path).
    if constexpr (SHARED_U && !U_STOR_FF) {
      u_reg[0] = shared_u_ram[stage].read(index);
    }
  }

  // ======== Reuse Read (FF Cache) ========

  // Return cached counter for a given bank and slot.
  // Only valid when USE_FF_CACHE=true and a prior read() was done.
  auto reuseRead(u64 bank, val<clog2(BPB)> slot) {
    static_assert(USE_FF_CACHE, "reuseRead requires USE_FF_CACHE=true");
    return ctr_regs[bank].select(slot);
  }

  // ======== Accessors ========
  // Return val<> (not reg<>) to avoid creating/destroying temporary registers.

  val<1> getHit(u64 bank = 0) {
    if constexpr (SHARED_TAG) {
      return val<1>{hit_reg[0]};
    } else {
      return val<1>{hit_reg[bank]};
    }
  }

  val<TAG_WIDTH> getTag(u64 bank = 0) {
    if constexpr (SHARED_TAG) {
      return val<TAG_WIDTH>{tag_reg[0]};
    } else {
      return val<TAG_WIDTH>{tag_reg[bank]};
    }
  }

  val<CTR_WIDTH> getCounter(u64 bank, u64 slot = 0) {
    if constexpr (USE_FF_CACHE) {
      return val<CTR_WIDTH>{ctr_regs[bank][slot]};
    } else {
      return val<CTR_WIDTH>{ctr_regs[bank][0]};
    }
  }

  val<U_WIDTH> getU(u64 bank = 0) {
    if constexpr (SHARED_U) {
      return val<U_WIDTH>{u_reg[0]};
    } else {
      return val<U_WIDTH>{u_reg[bank]};
    }
  }

  // ======== Write Interface ========

  // Compute write fanout for index, tag, u
  static constexpr u64 IDX_WRITE_FANOUT =
      1 // bank_ram write
      + (SHARED_TAG ? 1 : 0) // shared_tag_ram write
      + (SHARED_U && !U_STOR_FF ? 1 : 0) // shared_u_ram write
      + (U_STOR_FF ? 1 : 0); // write_u_ff decode comparisons

  static constexpr u64 TAG_WRITE_FANOUT =
      1 // pack_bank_entry
      + (SHARED_TAG ? 1 : 0); // shared_tag_ram write

  static constexpr u64 U_WRITE_FANOUT =
      1 // pack_bank_entry or u_ff write
      + ((SHARED_U && !U_STOR_FF) || U_STOR_FF ? 1 : 0); // separate u storage

  // Write a single bank entry back to SRAM.
  // stage: ahead stage index (0 when USE_AHEAD=false)
  void write(val<IDX_BITS> index, u64 bank, u64 stage,
             val<TAG_WIDTH> tag,
             val<CTR_BITS> ctr_bits,
             val<U_WIDTH> u) {
    // Fanout declarations
    if constexpr (IDX_WRITE_FANOUT > 1) {
      index.fanout(hard<IDX_WRITE_FANOUT>{});
    }
    if constexpr (TAG_WRITE_FANOUT > 1) {
      tag.fanout(hard<TAG_WRITE_FANOUT>{});
    }
    if constexpr (U_WRITE_FANOUT > 1) {
      u.fanout(hard<U_WRITE_FANOUT>{});
    }

    u64 phys_bank = stage * NUM_BANKS + bank;
    val<BANK_ENTRY_WIDTH> entry = pack_bank_entry(ctr_bits, tag, u);
    bank_ram[phys_bank].write(index, entry);

    // Write shared tag RAM if applicable
    if constexpr (SHARED_TAG) {
      shared_tag_ram[stage].write(index, tag);
    }

    // Write u-bit
    if constexpr (SHARED_U && !U_STOR_FF) {
      shared_u_ram[stage].write(index, u);
    } else if constexpr (SHARED_U && U_STOR_FF) {
      write_u_ff(0, index, u);
    } else if constexpr (!SHARED_U && U_STOR_FF) {
      write_u_ff(bank, index, u);
    }
    // When !SHARED_U && !U_STOR_FF, u is already packed in bank entry
  }

  // Convenience: write back from cached registers (for counter updates).
  // Uses the index from the last read().
  void writeBack(u64 bank, u64 stage, val<TAG_WIDTH> tag, val<U_WIDTH> u) {
    val<CTR_BITS> ctr_bits = pack_counters(bank);
    write(idx_reg, bank, stage, tag, ctr_bits, u);
  }

  // Update a single counter in the cache registers.
  // Call before writeBack() to modify a counter value.
  void writeReg(u64 bank, u64 slot, val<CTR_WIDTH> new_ctr) {
    ctr_regs[bank][slot] = new_ctr;
  }

  // ======== U-bit Reset (FF mode only) ========

  // Apply ResetFn to all u-bit flip-flops with the given mode.
  void reset_u(val<ResetFn::MODE_BITS> mode) {
    static_assert(U_STOR_FF, "reset_u requires U_STOR_FF=true");
    for (u64 a = 0; a < U_FF_ARRAYS; a++) {
      for (u64 i = 0; i < TABLE_SIZE; i++) {
        u_ff[a][i] = ResetFn::template apply<U_WIDTH>(u_ff[a][i], mode);
      }
    }
  }
};
