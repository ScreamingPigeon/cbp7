#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"

using namespace hcm;

// Minimal ahead pipeline timing test — replicating gshareN_ahead_best pattern.

template <u64 LOGG = 14>
struct TageAheadTest : predictor {
  static constexpr u64 LOGLINEINST = 10;
  static constexpr u64 LINEINST = 1 << LOGLINEINST;
  static constexpr u64 N = 7;
  static constexpr u64 LOGLANES = std::bit_width(N - 1);
  static constexpr u64 LANES = 1 << LOGLANES;
  static constexpr u64 LOGBANKS = std::bit_width(N);
  static constexpr u64 BANKS = 1 << LOGBANKS;
  static constexpr u64 index_bits = LOGG - (LOGLANES + LOGBANKS);

  // Exact gshareN_ahead_best storage
  reg<index_bits> index[2];
  reg<LOGBANKS> path;
  reg<LOGBANKS> XB;
  reg<LANES> XL;

  arr<reg<LANES>, BANKS> block_pred[2];
  reg<LANES> unordered_pred;
  arr<reg<1>, LANES> pred;

  reg<1> true_block = 1;
  reg<1> last_condbr_dir = 1;
  reg<LOGLINEINST> block_entry;
  reg<N + 1> rank;

  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<1>, N> branch_dir;

  hcm::ram<arr<val<LANES>, BANKS>, (1 << index_bits)> ctr_hi;

  val<1> line_end() { return (block_entry + block_size) == hard<LINEINST>{}; }
  val<1> last_pred() { return rank >> (N - num_branch); }

  val<1> predict1(val<64> inst_pc) override {
    inst_pc.fanout(hard<4>{});
    true_block.fanout(hard<8 + BANKS * 2>{});

    block_entry = select(true_block, val<LOGLINEINST>{inst_pc >> 2},
                         val<LOGLINEINST>{block_entry + block_size});
    block_entry.fanout(hard<LINEINST + N + 1>{});

    rank = select(true_block, val<N + 1>{1}, rank << num_branch);
    rank.fanout(hard<N + 2>{});

    XL = select(true_block, val<LOGLANES>{inst_pc >> 6}.decode().concat(),
                XL.rotate_left(num_branch));
    XL.fanout(hard<LANES>{});

    execute_if(true_block, [&]() {
      index[1] = index[0];
      val<index_bits> pc_bits = inst_pc >> (LOGBANKS + 2);
      index[0] = pc_bits.fo1();
      block_pred[1] = block_pred[0].fo1();
      block_pred[0] = ctr_hi.read(index[0]);
      path = XB + num_branch + ~last_condbr_dir;
      unordered_pred = block_pred[1].select(path);
      unordered_pred.fanout(hard<LANES>{});
    });

    XB = select(true_block, val<LOGBANKS>{inst_pc >> 6},
                val<LOGBANKS>{XB.fo1() + num_branch});

    for (u64 i = 0; i < LANES; i++) {
      pred[i] = (unordered_pred & XL.rotate_left(i)) != hard<0>{};
    }
    pred.fanout(hard<LINEINST * 2>{});
    block_size = 1;
    num_branch = 0;
    reuse_prediction(~line_end());
    return pred[num_branch];
  }

  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) override {
    block_size++;
    reuse_prediction(~line_end());
    return pred[num_branch];
  }

  val<1> predict2([[maybe_unused]] val<64> inst_pc) override {
    return pred[num_branch];
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) override {
    return pred[num_branch];
  }

  void update_condbr([[maybe_unused]] val<64> branch_pc, val<1> taken,
                     [[maybe_unused]] val<64> next_pc) override {
    assert(num_branch < N);
    branch_dir[num_branch] = taken.fo1();
    num_branch++;
    reuse_prediction(~(line_end() | last_pred()));
  }

  void update_cycle(instruction_info &block_end_info) override {
    val<1> &mispredict = block_end_info.is_mispredict;

    if (num_branch == 0) {
      last_condbr_dir = 0;
      true_block = 1;
      return;
    }

    last_condbr_dir = branch_dir[num_branch - 1].fo1();

    need_extra_cycle(mispredict);

    // No writes — just update true_block
    true_block = arr<val<1>, 4>{~mispredict, last_condbr_dir, last_pred(), line_end()}
                     .fold_or();

    num_branch = 0;
  }
};
