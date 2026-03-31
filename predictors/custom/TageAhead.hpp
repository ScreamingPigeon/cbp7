#pragma once
#include "../../cbp.hpp"
#include "TageConfig.hpp"
#include <cassert>

#pragma GCC diagnostic ignored "-Wunused-parameter"

template <
    // Tage Params
    typename TableCfg = SweepTableConfig<>, u64 FETCH_WIDTH = 16,
    u64 BIMODAL_SIZE = 4096, u64 DECAY_CTR = 8, u64 DECAY_GRAN = 2,
    typename DecayPolicy = DecayAggressive, u64 UCTRBITS = 8, u64 PATHBITS = 4,

    // P1 Params
    u64 P1_TABLE_SIZE = 4096, u64 P1_HIST = 6,

    // Meta Prediction Params
    bool USE_META = true, u64 METABITS = 4, u64 METAPIPE = 2,

    // Ahead-specific
    u64 SEC_TAG_WIDTH = 4,
    u64 AHEAD_BANKS = 3,     // N: physical banks instantiated
    u64 AHEAD_BANK_BITS = 2> // log2(K): bank select width (K = 1 << AHEAD_BANK_BITS)

struct TageAhead : predictor {
  // Defining Constexprs
  static constexpr u64 NUM_TABLES = TableCfg::NUM_TABLES;
  static constexpr u64 LOG_FETCH_WIDTH = clog2(FETCH_WIDTH);
  static constexpr u64 AHEAD_K = (1 << AHEAD_BANK_BITS);
  static_assert(AHEAD_BANKS >= 1, "Need at least 1 bank");
  static_assert(AHEAD_BANKS <= AHEAD_K, "Cannot have more banks than select values");

  static constexpr u64 MAX_TAG_WIDTH = array_max(TableCfg::TAG_WIDTH);
  static constexpr u64 MAX_CTR_WIDTH = array_max(TableCfg::CTR_WIDTH);
  static constexpr u64 MAX_HYST_WIDTH = array_max(TableCfg::HYST_WIDTH);
  static constexpr u64 MAX_U_WIDTH = array_max(TableCfg::U_WIDTH);
  static constexpr u64 MAX_TABLE_SIZE = array_max(TableCfg::TABLE_SIZE);
  // This is for Tage BTW
  static constexpr u64 MAX_IDX_BITS = clog2(MAX_TABLE_SIZE);
  static constexpr u64 MINHIST = array_min(TableCfg::HIST_LEN);
  static constexpr u64 MAXHIST = array_max(TableCfg::HIST_LEN);

  // This is the portion of the Tag that is hashed
  // so if BPB > 1, the stored tag is actually [offset, hashed tag] where
  // offset is LOG_FETCH_WIDTH bits wide
  static constexpr u64 MAX_HTAGBITS = MAX_TAG_WIDTH - LOG_FETCH_WIDTH;

  static constexpr u64 LOG_BIMODAL_SIZE = clog2(BIMODAL_SIZE);
  static constexpr u64 BINDEX_BITS = LOG_BIMODAL_SIZE - LOG_FETCH_WIDTH;
  // this is the number of rows per BIMODAL ram.
  static constexpr u64 BIM_ENTRIES = BIMODAL_SIZE / FETCH_WIDTH;

  static constexpr u64 P1_INDEX_BITS = clog2(P1_TABLE_SIZE) - LOG_FETCH_WIDTH;
  static constexpr u64 P1_ENTRIES = P1_TABLE_SIZE / FETCH_WIDTH;

  // Bits needed to index ram tables
  static constexpr u64 LINEADDR_BITS = std::max(P1_INDEX_BITS, P1_HIST);

  // =========== Members ========================
  // Table tuple type
  // TODO: @Prakhar: Stop this hardcoding in future revision
  using Tables = typename MakeTableTuple<
      TableCfg, 1 /*BR_P_ENTRY*/, 1 /*NUM_BANKS*/, false /*USE_AHEAD*/,
      true /*SHARED_TAG*/, true /*SHARED_U*/, true /*SHARED_HYS*/,
      false /*U_STOR_FF*/, DECAY_CTR, DefaultResetFn, false /*USE_FF_CACHE*/,
      std::make_index_sequence<NUM_TABLES>>::type;

  // Truncate gindex to per-table IDX_BITS
  template <u64 I> auto tidx(auto &gi) {
    using Table = std::tuple_element_t<I, Tables>;
    return val<Table::IDX_BITS>{gi};
  }

  // Tracking History
  reg<P1_HIST> p1_hist_reg;
  geometric_folds<NUM_TABLES, MINHIST, MAXHIST, MAX_IDX_BITS, MAX_HTAGBITS>
      tage_folded_hists;
  reg<1> true_block = 1;

  // Gshare Primitives
  reg<P1_INDEX_BITS> p1_index1;
  arr<reg<1>, FETCH_WIDTH> readp1; // per-lane preds
  reg<FETCH_WIDTH> p1;
  hcm::ram<val<1>, P1_ENTRIES> p1_pred[FETCH_WIDTH];
  hcm::ram<val<1>, P1_ENTRIES> p1_hys[FETCH_WIDTH];

  // Tage Bimodal state
  reg<BINDEX_BITS> bindex;
  arr<reg<1>, FETCH_WIDTH> readb; // per-lane preds
  hcm::ram<val<1>, BIM_ENTRIES> bpred[FETCH_WIDTH];
  hcm::ram<val<1>, BIM_ENTRIES> bhyst[FETCH_WIDTH];

  // Tage Tables — one full set per bank (predecessor-indexed)
  // SIZE in TableCfg is per-bank. Total storage = SIZE * NUM_TABLES * AHEAD_BANKS.
  Tables tables[AHEAD_BANKS];

  // Ahead Pipeline regs — per-bank, 2-deep pipeline
  // [stage][bank]: stage 0 = current read, stage 1 = previous (ready for predict2)
  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> ahead_tag[2][AHEAD_BANKS];
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> ahead_pred[2][AHEAD_BANKS];
  arr<reg<std::max(u64(1), MAX_HYST_WIDTH)>, NUM_TABLES> ahead_hyst[2][AHEAD_BANKS];
  arr<reg<MAX_U_WIDTH>, NUM_TABLES> ahead_u[2][AHEAD_BANKS];

  // Shared across banks (same index for all banks of a predecessor)
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> tage_indexes[2];
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> ahead_htag[2];

  // Bank selection state
  reg<AHEAD_BANK_BITS> XB[2];  // address bits for bank scrambling, pipelined
  reg<AHEAD_BANK_BITS> path;   // path out of previous block (set in update_cycle)

  // Secondary tag — disambiguates when bank selection lands on a phantom bank
  reg<SEC_TAG_WIDTH> secondary_tag;

  // P2 result registers
  reg<NUM_TABLES> not_u_mask; // bitmask tracking usefulness =!=0
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH>
      match; // tag hits per offset, for N tables + bimodal, is a bitvector
  // Note this is one hot and not a log ctr to avoid making a mux chain.
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match1; // this is a one hot
  arr<reg<NUM_TABLES + 1>, FETCH_WIDTH> match2; // this is one hot

  arr<reg<1>, FETCH_WIDTH> tage_pred1; // primary pred
  arr<reg<1>, FETCH_WIDTH> tage_pred2; // alt pred

  reg<FETCH_WIDTH> p2; // final prediction reg

  // Meta prediction stuff
  arr<reg<METABITS>, METAPIPE> meta;    // signed counters
  arr<reg<1>, FETCH_WIDTH> newly_alloc; // per-offset flag was the provider
                                        // newly allocated?

  // U-bit management
  reg<UCTRBITS> epoch_ctr;
  std::conditional_t<(DECAY_CTR > 0), reg<DECAY_CTR == 0 ? 1 : DECAY_CTR>,
                     EmptyMember>
      decay_threshold; // adaptive decay threshold

  // Block-state (cpp based tracking)
  u64 num_branch = 0;
  u64 blk_size = 0;

  // Block-state (HARCOM)
  // keeps track of which instructions in the block are branches
  arr<reg<LOG_FETCH_WIDTH>, FETCH_WIDTH> br_offset;
  arr<reg<1>, FETCH_WIDTH> branch_dir;

  // so this jit marks the start of the block? and keeps getting left shifted
  // until we exit the block.
  reg<FETCH_WIDTH> blk_entry;

  //======== Simulator Inteface ========
  val<1> predict1(val<64> inst_pc) override {
    // ---- Block setup ----
    inst_pc.fanout(hard<4>{}); // P1 offset(1) + P1 lineaddr(1) + TAGE lineaddr(1) + XB(1)
    val<LOG_FETCH_WIDTH> offset = inst_pc >> 2;
    blk_entry = offset.fo1().decode().concat();
    blk_entry.fanout(hard<2 * FETCH_WIDTH + 4>{});
    blk_size = 1;

    // ---- P1 gshare ----
    val<LINEADDR_BITS> p1_lineaddr = inst_pc >> (LOG_FETCH_WIDTH + 2);
    if constexpr (P1_HIST <= P1_INDEX_BITS) {
      p1_index1 = p1_lineaddr.fo1() ^ (val<P1_INDEX_BITS>{p1_hist_reg}
                                       << (P1_INDEX_BITS - P1_HIST));
    } else {
      p1_index1 = p1_hist_reg.make_array(val<P1_INDEX_BITS>{})
                      .append(p1_lineaddr.fo1())
                      .fold_xor();
    }
    p1_index1.fanout(hard<FETCH_WIDTH>{});
    static_loop<FETCH_WIDTH>([&]<u64 l_offset>() {
      readp1[l_offset] = p1_pred[l_offset].read(p1_index1);
    });
    readp1.fanout(hard<2>{});
    p1 = readp1.concat();
    p1.fanout(hard<2 * FETCH_WIDTH + 1>{});

    // ---- Phase 2: ahead TAGE reads (predecessor-indexed, banked) ----

    // Pipeline shift: [0] → [1] for all banks
    for (u64 b = 0; b < AHEAD_BANKS; b++) {
      ahead_tag[0][b].fanout(hard<2>{});
      ahead_pred[0][b].fanout(hard<2>{});
      if constexpr (MAX_HYST_WIDTH > 0) {
        ahead_hyst[0][b].fanout(hard<2>{});
      }
      ahead_u[0][b].fanout(hard<2>{});
      for (u64 i = 0; i < NUM_TABLES; i++) {
        ahead_tag[1][b][i]  = ahead_tag[0][b][i];
        ahead_pred[1][b][i] = ahead_pred[0][b][i];
        if constexpr (MAX_HYST_WIDTH > 0) {
          ahead_hyst[1][b][i] = ahead_hyst[0][b][i];
        }
        ahead_u[1][b][i] = ahead_u[0][b][i];
      }
    }
    // Pipeline shift: shared index/tag state
    tage_indexes[0].fanout(hard<2>{});
    ahead_htag[0].fanout(hard<2>{});
    for (u64 i = 0; i < NUM_TABLES; i++) {
      tage_indexes[1][i] = tage_indexes[0][i];
      ahead_htag[1][i] = ahead_htag[0][i];
    }
    // Pipeline shift: bank selection
    XB[0].fanout(hard<2>{});
    XB[1] = XB[0];
    XB[1].fanout(hard<2>{}); // read in predict2 bank_sel + update_cycle

    // Compute XB[0] from current PC
    XB[0] = inst_pc >> 2;

    // Compute TAGE indices for current block (shared across banks)
    val<LINEADDR_BITS> tage_lineaddr = inst_pc >> (LOG_FETCH_WIDTH + 2);
    tage_lineaddr.fanout(hard<1 + NUM_TABLES * 2>{});
    tage_folded_hists.fanout(hard<2>{});
    for (u64 i = 0; i < NUM_TABLES; i++) {
      tage_indexes[0][i] = tage_lineaddr ^ tage_folded_hists.template get<0>(i);
    }
    tage_indexes[0].fanout(hard<AHEAD_BANKS * 4>{}); // each bank reads tag+pred+hyst+u

    // Compute hashed tags (parallel with RAM reads, shared across banks)
    for (u64 i = 0; i < NUM_TABLES; i++) {
      ahead_htag[0][i] = val<MAX_HTAGBITS>{tage_lineaddr}.reverse()
                          ^ tage_folded_hists.template get<1>(i);
    }

    // Read ALL banks into ahead_*[0][bank]
    for (u64 b = 0; b < AHEAD_BANKS; b++) {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        auto &t = std::get<I>(tables[b]);
        ahead_tag[0][b][I]  = t.tag_ram[0].read(tidx<I>(tage_indexes[0][I]));
        ahead_pred[0][b][I] = t.pred_ram[0].read(tidx<I>(tage_indexes[0][I]));
        if constexpr (MAX_HYST_WIDTH > 0) {
          ahead_hyst[0][b][I] = t.hyst_ram[0].read(tidx<I>(tage_indexes[0][I]));
        }
        ahead_u[0][b][I] = t.u_ram[0].read(tidx<I>(tage_indexes[0][I]));
      });
    }

    // TODO: compute + store secondary tag

    return (blk_entry & p1) != hard<0>{};
  };
  val<1> reuse_predict1(val<64> inst_pc) override {
    return ((blk_entry << blk_size) & p1) != hard<0>{};
  };
  // P2 proxies P1 for now
  val<1> predict2(val<64> inst_pc) override {
    val<1> taken = (blk_entry & p1) != hard<0>{};
    reuse_prediction(~val<1>{blk_entry >> (FETCH_WIDTH - 1)});
    return taken;
  };
  val<1> reuse_predict2(val<64> inst_pc) override {
    val<1> taken = ((blk_entry << blk_size) & p1) != hard<0>{};
    reuse_prediction(~val<1>{blk_entry >> (FETCH_WIDTH - 1 - blk_size)});
    blk_size++;
    return taken;
  };
  void update_condbr(val<64> branch_pc, val<1> taken,
                     val<64> next_pc) override {
    assert(num_branch < FETCH_WIDTH);
    br_offset[num_branch] = branch_pc.fo1() >> 2;
    branch_dir[num_branch] = taken.fo1();
    num_branch++;
  };
  void update_cycle(instruction_info &block_end_info) override {
    val<1> &mispredict = block_end_info.is_mispredict;
    val<64> &next_pc = block_end_info.next_pc;

    // No-branch early return: update history and exit
    if (num_branch == 0) {
      val<1> line_end = blk_entry >> (FETCH_WIDTH - blk_size);
      val<1> actual_block = ~(true_block & line_end.fo1());
      actual_block.fanout(hard<MAXHIST + NUM_TABLES * 2 + 2>{});
      execute_if(actual_block, [&]() {
        next_pc.fanout(hard<2>{});
        p1_hist_reg = (p1_hist_reg << 1) ^ val<P1_HIST>{next_pc >> 2};
        tage_folded_hists.update(val<PATHBITS>{next_pc >> 2});
        true_block = 1;
      });
      return;
    }

    // TODO Phase 1.2: P1 pred/hyst update when P1 disagrees with P2
    // Needs: disagree mask, p1_hys read, conditional writes

    // End-of-block history update
    val<1> correct_pred = ~mispredict;
    val<1> line_end = blk_entry >> (FETCH_WIDTH - blk_size);
    true_block = correct_pred | branch_dir[num_branch - 1] | line_end.fo1();
    true_block.fanout(hard<MAXHIST + NUM_TABLES * 2 + 2>{});
    execute_if(true_block, [&]() {
      next_pc.fanout(hard<2>{});
      p1_hist_reg = (p1_hist_reg << 1) ^ val<P1_HIST>{next_pc >> 2};
      tage_folded_hists.update(val<PATHBITS>{next_pc >> 2});
    });

    num_branch = 0;
  };

  // ==============
};
