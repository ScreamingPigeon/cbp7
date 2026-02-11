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

template <u64 TABLE_SIZE = 64, u64 TABLE_HIST = 64, u64 TAG_WIDTH = 7,
          u64 CTR_WIDTH = 3, u64 U_WIDTH = 2, bool USE_HYS = true,
          u64 HYS_WIDTH = 2, u64 PRED_BLK_SIZE = 8, u64 DECAY_CTR = 1024>
class TageTable {

  static constexpr u64 BITS_PER_INST =
      CTR_WIDTH + U_WIDTH + ((USE_HYS) ? HYS_WIDTH : 0);
  static constexpr u64 BITS_PER_ENTRY =
      TAG_WIDTH + PRED_BLK_SIZE * BITS_PER_INST;

public:
  TageTable() {}
  ~TageTable() {}

  auto newRead(val<clog2(TABLE_SIZE)> idx, val<TAG_WIDTH> tag,
               val<clog2(PRED_BLK_SIZE)> slot_idx) {
    val<BITS_PER_ENTRY> ram_read = table_ram.read(idx);

    auto [read_tag, pred_data] =
        split<TAG_WIDTH, PRED_BLK_SIZE * BITS_PER_INST>(ram_read);

    val<1> is_hit = (read_tag == tag);

    // NOTE:
    // This is an interesting case. As of this commit,
    // I am fanning out is_hit and pred_data. To the cache_regs
    // as well as to the predictor to facilitate a same cycle
    // read. This means lower latency on newRead, but more
    // power and comb delay due to fanout.
    // An alternative could be to read through the cache_regs
    // on the next cycle via reuseRead() - less power, but 1 extra cycle of
    // latency.

    // Convert to array for per-slot access
    arr<val<BITS_PER_INST>, PRED_BLK_SIZE> pred_array =
        pred_data.template makearray<val<BITS_PER_INST>>();

    // Set regs
    hit = is_hit;
    tag_reg = tag;
    idx_reg = idx;
    executeif(is_hit, [&]() { cache_regs = pred_array; });

    // increment decay_ctr
    decay_ctr = decay_ctr + 1;

    // Return the selected slot using .select()
    return std::make_tuple(is_hit, pred_array.select(slot_idx));
  }
  auto reuseRead(val<clog2(PRED_BLK_SIZE)> slot_idx) {
    return cache_regs[slot_idx];
    decay_ctr = decay_ctr + 1;
  }

  auto writeReg(val<BITS_PER_INST> new_data,
                val<clog2(PRED_BLK_SIZE)> slot_idx) {
    cache_regs[slot_idx] = (new_data);
  }

  void updateRamFromRegs() {
    val<PRED_BLK_SIZE * BITS_PER_INST> regs_data = cache_regs.concat();
    // Combine tag and data into full entry
    val<BITS_PER_ENTRY> full_entry = concat(tag_reg, regs_data);
    // Write back to RAM
    table_ram.write(idx_reg, full_entry);
  }

private:
  // Instantiate the Registers for caching
  hcm::ram<val<BITS_PER_ENTRY>, TABLE_SIZE> table_ram;
  reg<1> hit;

  reg<clog2(DECAY_CTR)> decay_ctr = 0;
  reg<TAG_WIDTH> tag_reg;
  reg<clog2(TABLE_SIZE)> idx_reg;
  arr<reg<BITS_PER_INST>, PRED_BLK_SIZE> cache_regs;
  val<BITS_PER_INST> getSlotFromEntry(val<BITS_PER_ENTRY> ram_entry) {}

  void decayAllEntries() {
    // decays U on all entries
    // TODO:
  }
};
