#pragma once

#include <cstdlib>

#include "../../cbp.hpp"
#include "../../harcom.hpp"

using namespace hcm;

// Note: std::rand() is used for probabilistic decay instead of LSFR
// (see agent/rng.md). This has zero modeled cost in HARCOM.

// Constexpr helper
constexpr size_t clog2(std::uint64_t n) {
  int r = 0;
  std::uint64_t v = 1;
  while (v < n) {
    v <<= 1;
    ++r;
  }
  return r;
}

// Forward declaration for conditional Empty type (zero-cost placeholder)
namespace _tage_table_detail {
  struct Empty {};
}

template <u64 TABLE_SIZE = 64, u64 TABLE_HIST = 64, u64 TAG_WIDTH = 7,
          u64 CTR_WIDTH = 3, u64 U_WIDTH = 2, u64 PRED_BLK_SIZE = 8,
          u64 DECAY_CTR = 1024, bool EN_N_BLK_RD = true,
          bool USE_LSFR_DECAY = true>
class TageTable {

  static constexpr u64 BITS_PER_INST = CTR_WIDTH;
  // U is only packed in SRAM entry when USE_LSFR_DECAY = true
  static constexpr u64 BITS_PER_ENTRY =
      TAG_WIDTH + PRED_BLK_SIZE * BITS_PER_INST + (USE_LSFR_DECAY ? U_WIDTH : 0);

public:
  TageTable() {}
  ~TageTable() {}

  // Return critical path latency in picoseconds
  // This is the maximum latency across all public methods
  static constexpr u64 critical_path_ps() {
    // Extract RAM read latency from table_ram type
    // The latency is stored in static_ram::LATENCY
    using ram_type = hcm::ram<hcm::val<BITS_PER_ENTRY>, TABLE_SIZE>;
    constexpr u64 ram_read_latency = static_cast<u64>(ram_type::static_ram::LATENCY);

    // MUX selection overhead for arr.select()
    constexpr u64 mux_overhead = 20;  // ps

    return ram_read_latency + mux_overhead;
  }

  val<CTR_WIDTH> newRead(val<clog2(TABLE_SIZE)> idx, val<TAG_WIDTH> tag,
                         val<clog2(PRED_BLK_SIZE)> slot_idx) {
    // Read entry from RAM
    val<BITS_PER_ENTRY> entry = table_ram.read(idx);

    // Store index and tag
    idx_reg = idx;
    tag_reg = tag;

    // Read and cache U bits, extract counters (branch based on storage mode)
    if constexpr (USE_LSFR_DECAY) {
      // SRAM mode: U is packed at LSB of entry
      auto [u_bits, ctr_bits_combined, tag_bits] =
          split<U_WIDTH, PRED_BLK_SIZE * CTR_WIDTH, TAG_WIDTH>(entry);
      u_reg = u_bits;

      // Extract and store individual counters
      // Counter layout: ctr_regs[0] at LSB, ctr_regs[PRED_BLK_SIZE-1] at MSB
      arr<val<CTR_WIDTH>, PRED_BLK_SIZE> ctrs;
      static_loop<PRED_BLK_SIZE>([&]<int I>() {
        constexpr u64 shift_amt = I * CTR_WIDTH;
        val<CTR_WIDTH> shifted = (ctr_bits_combined >> hard<shift_amt>{});
        ctrs[I] = shifted; // Store in val array for return
        ctr_regs[I] = shifted; // Cache in register for reuse
      });

      // Check if tags match
      hit = (tag_bits == tag);

      // Return counter for requested slot with fo2 pipeline register
      if (EN_N_BLK_RD)
        return ctrs.select(slot_idx).fo2();
      else
        return ctrs.select(slot_idx);
    } else {
      // FF mode: U is in u_ff array, entry is just ctr_bits || tag (no U packed)
      auto [ctr_bits_combined, tag_bits] =
          split<PRED_BLK_SIZE * CTR_WIDTH, TAG_WIDTH>(entry);

      // Read U from FF array using MUX
      u_reg = u_ff.select(idx);

      // Extract and store individual counters
      arr<val<CTR_WIDTH>, PRED_BLK_SIZE> ctrs;
      static_loop<PRED_BLK_SIZE>([&]<int I>() {
        constexpr u64 shift_amt = I * CTR_WIDTH;
        val<CTR_WIDTH> shifted = (ctr_bits_combined >> hard<shift_amt>{});
        ctrs[I] = shifted;
        ctr_regs[I] = shifted;
      });

      // Check if tags match
      hit = (tag_bits == tag);

      // Return counter for requested slot with fo2 pipeline register
      if (EN_N_BLK_RD)
        return ctrs.select(slot_idx).fo2();
      else
        return ctrs.select(slot_idx);
    }
  }
  auto reuseRead(val<clog2(PRED_BLK_SIZE)> slot_idx) {
    return ctr_regs[slot_idx];
  }

  auto writeReg(val<BITS_PER_INST> new_data,
                val<clog2(PRED_BLK_SIZE)> slot_idx) {
    ctr_regs[slot_idx] = new_data;
  }

  void updateBlock(val<1> use_regs, val<TAG_WIDTH> tag,
                   val<BITS_PER_ENTRY> new_entry) {
    auto ctr_bits = ctr_regs.concat();

    if constexpr (USE_LSFR_DECAY) {
      // SRAM mode: U is packed in entry, apply probabilistic decay
      val<1> is_miss = ~hit;
      val<1> should_decay = is_miss & decrementU();

      // Decrement U only if it's non-zero, else keep it zero
      val<U_WIDTH> decremented_u = select(u_reg == hard<0>{},
                                          val<U_WIDTH>{0},
                                          val<U_WIDTH>{u_reg - 1});
      val<U_WIDTH> updated_u = select(should_decay, decremented_u, u_reg);

      // Build entry from registers: tag_reg || ctr_regs || updated_u
      val<BITS_PER_ENTRY> reg_entry = concat(updated_u, ctr_bits, tag_reg);

      // Select between reg_entry and new_entry based on use_regs
      val<BITS_PER_ENTRY> entry_to_write = select(use_regs, reg_entry, new_entry);

      // Write to RAM
      table_ram.write(idx_reg, entry_to_write);
    } else {
      // FF mode: U is in u_ff array, write u_reg to FF array
      u_ff[idx_reg] = u_reg;

      // Build entry without U: tag_reg || ctr_regs (no u_reg)
      val<BITS_PER_ENTRY> reg_entry = concat(ctr_bits, tag_reg);

      // Select between reg_entry and new_entry based on use_regs
      val<BITS_PER_ENTRY> entry_to_write = select(use_regs, reg_entry, new_entry);

      // Write to RAM
      table_ram.write(idx_reg, entry_to_write);
    }
  }

  auto getUsefulness() { return u_reg; }
  auto getHit() { return hit; }

  void setThreshold(val<clog2(DECAY_CTR)> new_threshold) {
    decay_threshold = new_threshold;
  }

private:
  // Instantiate the Registers for caching
  hcm::ram<val<BITS_PER_ENTRY>, TABLE_SIZE> table_ram;

  // Whether currently cached entry is a hit or not
  reg<1> hit;

  // Tag and index of currently cached entry
  reg<TAG_WIDTH> tag_reg;
  reg<clog2(TABLE_SIZE)> idx_reg;

  // U/CTR of currently cached entries
  // u_reg: always present — caches current entry's U for predictor access
  reg<U_WIDTH> u_reg;
  arr<reg<CTR_WIDTH>, PRED_BLK_SIZE> ctr_regs;

  // FF backing store for U bits (only instantiated when USE_LSFR_DECAY = false)
  // When USE_LSFR_DECAY = true: this is _tage_table_detail::Empty (zero cost)
  // When USE_LSFR_DECAY = false: this is arr<reg<U_WIDTH>, TABLE_SIZE> (flipflop array)
  std::conditional_t<USE_LSFR_DECAY, _tage_table_detail::Empty, arr<reg<U_WIDTH>, TABLE_SIZE>> u_ff;

  // Probabilistic decay threshold (configurable by predictor)
  reg<clog2(DECAY_CTR)> decay_threshold;

  val<1> decrementU() {
    // On tag miss, decrement U with probability 1/DECAY_CTR
    // Generate random bits and compare against threshold
    constexpr u64 DECAY_BITS = clog2(DECAY_CTR);
    val<DECAY_BITS> rng{static_cast<unsigned>(std::rand())};
    return rng == hard<0>{};  // true with probability 1/DECAY_CTR
  }
};
