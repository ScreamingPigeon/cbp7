#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"

using namespace hcm;

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

/* TODO:
 * I need to implement a parameterized LSFR to do
 * dynamic threshold probabilistic decay
 */

template <u64 TABLE_SIZE = 64, u64 TABLE_HIST = 64, u64 TAG_WIDTH = 7,
          u64 CTR_WIDTH = 3, u64 U_WIDTH = 2, u64 PRED_BLK_SIZE = 8,
          u64 DECAY_CTR = 1024, bool EN_N_BLK_RD = true>
class TageTable {

  static constexpr u64 BITS_PER_INST = CTR_WIDTH;
  static constexpr u64 BITS_PER_ENTRY =
      TAG_WIDTH + PRED_BLK_SIZE * BITS_PER_INST + U_WIDTH;

public:
  TageTable() {}
  ~TageTable() {}

  auto newRead(val<clog2(TABLE_SIZE)> idx, val<TAG_WIDTH> tag,
               val<clog2(PRED_BLK_SIZE)> slot_idx) {
    // Read entry from RAM
    val<BITS_PER_ENTRY> entry = table_ram.read(idx);

    // Split entry into components (LSB to MSB): u || counters || tag
    auto [u_bits, ctr_bits_combined, tag_bits] =
        split<U_WIDTH, PRED_BLK_SIZE * CTR_WIDTH, TAG_WIDTH>(entry);

    // Store index and tag
    idx_reg = idx;
    tag_reg = tag;

    // Store u-bit
    u_reg = u_bits;

    // Extract and store individual counters
    // Counter layout: ctr_regs[0] at LSB, ctr_regs[PRED_BLK_SIZE-1] at MSB
    arr<val<CTR_WIDTH>, PRED_BLK_SIZE> ctrs;
    static_loop<PRED_BLK_SIZE>([&]<int I>() {
      constexpr u64 shift_amt = I * CTR_WIDTH;
      val<PRED_BLK_SIZE * CTR_WIDTH> shifted =
          ctr_bits_combined >> hard<shift_amt>{};
      ctrs[I] = shifted; // Store in val array for return
                         // Hopefully, the wires for I != slot_idx are opted out
      ctr_regs[I] = shifted; // Cache in register for reuse
    });

    // Check if tags match
    hit = (tag_bits == tag);

    // Return from SRAM read result (not register) with fo2
    if (EN_N_BLK_RD)
      return ctrs.select(slot_idx).fo2();
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
    // Build entry from registers: tag_reg || ctr_regs || u_reg
    auto ctr_bits = ctr_regs.concat();
    val<BITS_PER_ENTRY> reg_entry = concat(u_reg, ctr_bits, tag_reg);

    // Select between reg_entry and new_entry based on use_regs
    val<BITS_PER_ENTRY> entry_to_write = select(use_regs, reg_entry, new_entry);

    // Write to RAM
    table_ram.write(idx_reg, entry_to_write);
  }

  auto getUsefulness() { return u_reg; }
  auto getHit() { return hit; }

  auto setThreshold() {
    // TODO: Let the predictor set the new
    // threshold for decrementing U
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
  reg<U_WIDTH> u_reg;
  arr<reg<CTR_WIDTH>, PRED_BLK_SIZE> ctr_regs;

  auto decrementU() {
    // TODO: If tag is a miss, compute probabilistic decay
    // if greater than threshold, return Ureg decrement.
    // Invoke inside updateBlock, so we knowreplace reg_entry
  };
};
