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

    // Ahead-sepcific
    u64 SEC_TAG_WIDTH = 4>

struct TageAhead : predictor {
  // Defining Constexpres
  static constexpr u64 NUM_TABLES = TableCfg::NUM_TABLES;
  static constexpr u64 LOG_FETCH_WIDTH = clog2(FETCH_WIDTH);

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

  // Tage Tables
  Tables tables;

  // Ahead Pipeline specific regs
  // Note that the [2] is for 1 + N_AHEAD (1 in our case)
  // [0] is for current read, and [1] is for whatever was 1 ahead (i.e run
  // before)
  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> ahead_tag[2];
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> ahead_pred[2];
  arr<reg<std::max(u64(1), MAX_HYST_WIDTH)>, NUM_TABLES> ahead_hyst[2];
  arr<reg<MAX_U_WIDTH>, NUM_TABLES> ahead_u[2];
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> tage_indexes[2];
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> ahead_htag[2];
  reg<SEC_TAG_WIDTH> secondary_tag; // this is needed to determine when to use
                                    // bimodal on speculated next block

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
    // block setup
    // ===============
    inst_pc.fanout(hard<2>{});

    val<LOG_FETCH_WIDTH> offset =
        inst_pc >> 2; // get the offset for the start of the block
    blk_entry = offset.fo1().decode().concat(); // convert binary to one hot
    blk_entry.fanout(hard<FETCH_WIDTH + 4>{});  // setup fanout
    blk_size = 1;

    // p1 gshare logic
    // ===============

    // Take the instruction. Drop the Lowest
    // LOG FETCH_WIDTH bits, and
    // drop another 2 because instructions are 4 bytes wide
    val<LINEADDR_BITS> lineaddr = inst_pc >> (LOG_FETCH_WIDTH + 2);

    // TODO: @Prakhar. Add fanout for this lineaddr

    // compute the index
    if constexpr (P1_HIST <= P1_INDEX_BITS) {
      p1_index1 = lineaddr.fo1() ^ (val<P1_INDEX_BITS>{p1_hist_reg}
                                    << (P1_INDEX_BITS - P1_HIST));
    } else {
      p1_index1 = p1_hist_reg.make_array(val<P1_INDEX_BITS>{})
                      .append(lineaddr.fo1())
                      .fold_xor();
    }
    // TODO: @ Prakhar. Add p1_index1 fanout
    p1_index1.fanout(hard<FETCH_WIDTH>{});

    static_loop<FETCH_WIDTH>([&]<u64 l_offset>() {
      readp1[l_offset] = p1_pred[l_offset].read(
          p1_index1); // read prediction ram and load into regs
    });

    readp1.fanout(hard<2>{}); // read in p1 and in updat_cycle
    p1 = readp1.concat();
    p1.fanout(hard<FETCH_WIDTH + 1>{}); // reuse P1 calls and updat_cycle
                                        //
    // Return: "is any branch from the current instruction onward predicted
    // taken?" This looks wrong — if lane i is not-taken but lane i+3 is taken,
    // we return 1 for lane i. But the simulator only checks P1 accuracy for
    // conditional branches, and by the time we reach a branch at lane j via
    // reuse_predict1 calls, the block_entry mask has shifted to start at lane j
    // — so the prediction is correct for that branch. The != 0 reduction is
    // equivalent to a per-lane mux but avoids the mux tree latency (~0.18
    // cycles). All reference predictors use this pattern
    return (blk_entry & p1) != hard<0>{};
  };
  val<1> reuse_predict1(val<64> inst_pc) override {
    return ((blk_entry << blk_size) & p1) != hard<0>{};
  };
  // P2 proxies P1 for now — will be replaced with ahead TAGE in Phase 3
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
