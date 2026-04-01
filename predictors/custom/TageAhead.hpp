#pragma once
#include "../../cbp.hpp"
#include "TageConfig.hpp"
#include <cassert>

#pragma GCC diagnostic ignored "-Wunused-parameter"

template <
    // Tage Params
    typename TableCfg = SweepTableConfig<8, 512>, u64 FETCH_WIDTH = 16,
    u64 BIMODAL_SIZE = 4096, u64 DECAY_CTR = 8, u64 DECAY_GRAN = 2,
    typename DecayPolicy = DecayAggressive, u64 UCTRBITS = 8, u64 PATHBITS = 4,

    // P1 Params
    u64 P1_TABLE_SIZE = 4096, u64 P1_HIST = 6,

    // Meta Prediction Params
    bool USE_META = true, u64 METABITS = 4, u64 METAPIPE = 2,

    // Ahead-specific
    bool USE_SEC_TAG = false,     // enable secondary tag for additional disambiguation
    u64 SEC_TAG_WIDTH = 4,
    u64 AHEAD_BANKS = 3,          // N: physical banks instantiated
    u64 AHEAD_BANK_BITS = 2>      // log2(K): bank select width (K = 1 << AHEAD_BANK_BITS)

struct TageAhead : predictor {
  // Defining Constexprs
  static constexpr u64 NUM_TABLES = TableCfg::NUM_TABLES;
  static constexpr u64 LOG_FETCH_WIDTH = clog2(FETCH_WIDTH);
  static constexpr u64 AHEAD_K = (1 << AHEAD_BANK_BITS);
  static_assert(AHEAD_BANKS >= 1, "Need at least 1 bank");
  static_assert(AHEAD_BANKS <= AHEAD_K,
                "Cannot have more banks than select values");

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
  // One full set per bank (predecessor-indexed).
  // SIZE in TableCfg is per-bank. Total storage = SIZE * NUM_TABLES *
  // AHEAD_BANKS.
  Tables tables[AHEAD_BANKS];

  // Ahead Pipeline regs — per-bank, 2-deep pipeline
  // [stage][bank]: stage 0 = current read, stage 1 = previous (ready for
  // predict2)
  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> ahead_tag[2][AHEAD_BANKS];
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> ahead_pred[2][AHEAD_BANKS];
  arr<reg<MAX_HYST_WIDTH>, NUM_TABLES> ahead_hyst[2][AHEAD_BANKS];
  arr<reg<MAX_U_WIDTH>, NUM_TABLES> ahead_u[2][AHEAD_BANKS];

  // Precomputed in predict1: per-bank hashed tag compare result
  // Moves tag compare out of predict2's critical path.
  reg<NUM_TABLES> ahead_htagcmp[2][AHEAD_BANKS];

  // Shared across banks (same index for all banks of a predecessor)
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> tage_indexes[2];
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> ahead_htag[2];

  // Bank selection state
  reg<AHEAD_BANK_BITS> XB[2]; // address bits for bank scrambling, pipelined
  reg<AHEAD_BANK_BITS> path; // path out of previous block (set in update_cycle)

  // Secondary tag — disambiguates when bank selection lands on a phantom bank
  reg<SEC_TAG_WIDTH> secondary_tag;

  // P2 working registers — written in predict2, read in update_cycle
  // In predict2: use vals for prediction logic, write to these regs for
  // update_cycle
  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> readt;
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> readc;
  arr<reg<MAX_HYST_WIDTH>, NUM_TABLES> readh;
  arr<reg<MAX_U_WIDTH>, NUM_TABLES> readu;

  // P2 result registers
  reg<NUM_TABLES> not_u_mask; // bitmask tracking usefulness !=0
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
    inst_pc.fanout(hard<3>{}); // offset(1) + lineaddr(1) + XB(1)
    val<LOG_FETCH_WIDTH> offset = inst_pc >> 2;
    blk_entry = offset.fo1().decode().concat();
    blk_entry.fanout(hard<2 * FETCH_WIDTH + 4>{});
    blk_size = 1;

    // ---- P1 gshare ----
    val<LINEADDR_BITS> lineaddr = inst_pc >> (LOG_FETCH_WIDTH + 2);
    lineaddr.fanout(
        hard<1 + NUM_TABLES * 2>{}); // P1(1) + TAGE idx(N) + TAGE tag(N)
    if constexpr (P1_HIST <= P1_INDEX_BITS) {
      p1_index1 = lineaddr ^ (val<P1_INDEX_BITS>{p1_hist_reg}
                              << (P1_INDEX_BITS - P1_HIST));
    } else {
      p1_index1 = p1_hist_reg.make_array(val<P1_INDEX_BITS>{})
                      .append(lineaddr)
                      .fold_xor();
    }
    p1_index1.fanout(hard<FETCH_WIDTH * 4>{}); // predict1 reads(FW) + update_cycle reads(FW*3)
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
      ahead_htagcmp[0][b].fanout(hard<2>{});
      for (u64 i = 0; i < NUM_TABLES; i++) {
        ahead_tag[1][b][i] = ahead_tag[0][b][i];
        ahead_pred[1][b][i] = ahead_pred[0][b][i];
        if constexpr (MAX_HYST_WIDTH > 0) {
          ahead_hyst[1][b][i] = ahead_hyst[0][b][i];
        }
        ahead_u[1][b][i] = ahead_u[0][b][i];
      }
      ahead_htagcmp[1][b] = ahead_htagcmp[0][b];
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
    XB[0] = inst_pc >> 2;

    // Compute TAGE indices (shared lineaddr, already computed above)
    tage_folded_hists.fanout(hard<4>{}); // predict1(2) + update_cycle(2)
    for (u64 i = 0; i < NUM_TABLES; i++) {
      tage_indexes[0][i] = lineaddr ^ tage_folded_hists.template get<0>(i);
    }
    tage_indexes[0].fanout(
        hard<AHEAD_BANKS * 4>{}); // each bank reads tag+pred+hyst+u

    // Compute hashed tags (parallel with RAM reads, shared across banks)
    for (u64 i = 0; i < NUM_TABLES; i++) {
      ahead_htag[0][i] = val<MAX_HTAGBITS>{lineaddr}.reverse() ^
                         tage_folded_hists.template get<1>(i);
    }

    // Read ALL banks into ahead_*[0][bank]
    for (u64 b = 0; b < AHEAD_BANKS; b++) {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        auto &t = std::get<I>(tables[b]);
        ahead_tag[0][b][I] = t.tag_ram[0].read(tidx<I>(tage_indexes[0][I]));
        ahead_pred[0][b][I] = t.pred_ram[0].read(tidx<I>(tage_indexes[0][I]));
        if constexpr (MAX_HYST_WIDTH > 0) {
          ahead_hyst[0][b][I] = t.hyst_ram[0].read(tidx<I>(tage_indexes[0][I]));
        }
        ahead_u[0][b][I] = t.u_ram[0].read(tidx<I>(tage_indexes[0][I]));
      });
    }

    // ---- Precompute hashed tag compare per bank ----
    // ahead_htag[0] and ahead_tag[0][b] were just written by RAM reads above.
    // Comparing them here moves the tag compare out of predict2's critical path.
    ahead_htag[0].fanout(hard<AHEAD_BANKS * NUM_TABLES>{});
    for (u64 b = 0; b < AHEAD_BANKS; b++) {
      ahead_tag[0][b].fanout(hard<NUM_TABLES>{}); // one read per table for compare
      arr<val<1>, NUM_TABLES> cmp = [&](int i) -> val<1> {
        return val<MAX_HTAGBITS>{ahead_tag[0][b][i]} == ahead_htag[0][i];
      };
      ahead_htagcmp[0][b] = cmp.fo1().concat();
    }

    // TODO: compute + store secondary tag

    // OVERRIDER HOOK: overrider.prefetch(inst_pc, raw_pc)
    // Read overrider tables (SC/loop) at line-level index into regs.

    return (blk_entry & p1) != hard<0>{};
  };
  val<1> reuse_predict1(val<64> inst_pc) override {
    return ((blk_entry << blk_size) & p1) != hard<0>{};
  };
  val<1> predict2(val<64> inst_pc) override {
    // ---- Step 1: Bimodal read (direct, fallback for bank miss) ----
    bindex = inst_pc.fo1() >> (LOG_FETCH_WIDTH + 2);
    bindex.fanout(hard<FETCH_WIDTH * 4>{}); // predict2(FW) + update_cycle(FW*3)
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      readb[offset] = bpred[offset].read(bindex);
    }
    readb.fanout(hard<4>{}); // predict2(2) + update_cycle(2)

    // ---- Step 2: Bank select ----
    path.fanout(hard<4>{}); // predict2(2) + update_cycle(2)
    val<AHEAD_BANK_BITS> raw_bank_sel = path ^ XB[1];
    val<1> bank_hit = (raw_bank_sel < hard<AHEAD_BANKS>{});
    // Clamp to valid bank range — phantom banks map to 0 (result discarded via
    // bank_hit)
    val<AHEAD_BANK_BITS> bank_sel =
        select(bank_hit, raw_bank_sel, val<AHEAD_BANK_BITS>{0});
    bank_sel.fanout(hard<NUM_TABLES * 4 + 1>{});
    bank_hit.fanout(hard<FETCH_WIDTH>{});

    // ---- Step 3: Select correct bank from ahead_*[1][bank_sel] ----
    // Mux across AHEAD_BANKS using bank_sel. Results are vals for prediction
    // logic (same cycle), then written to regs for update_cycle (next cycle).

    for (u64 b = 0; b < AHEAD_BANKS; b++) {
      ahead_tag[1][b].fanout(hard<2>{});
      ahead_pred[1][b].fanout(hard<2>{});
      ahead_u[1][b].fanout(hard<2>{});
      if constexpr (MAX_HYST_WIDTH > 0) {
        ahead_hyst[1][b].fanout(hard<2>{});
      }
    }

    // Per-table bank mux → vals for prediction, regs for update_cycle
    arr<val<MAX_TAG_WIDTH>, NUM_TABLES> v_tag =
        [&](u64 i) -> val<MAX_TAG_WIDTH> {
      arr<val<MAX_TAG_WIDTH>, AHEAD_BANKS> opts =
          [&](u64 b) -> val<MAX_TAG_WIDTH> { return ahead_tag[1][b][i]; };
      return opts.fo1().select(bank_sel);
    };
    arr<val<MAX_CTR_WIDTH>, NUM_TABLES> v_pred =
        [&](u64 i) -> val<MAX_CTR_WIDTH> {
      arr<val<MAX_CTR_WIDTH>, AHEAD_BANKS> opts =
          [&](u64 b) -> val<MAX_CTR_WIDTH> { return ahead_pred[1][b][i]; };
      return opts.fo1().select(bank_sel);
    };
    arr<val<MAX_U_WIDTH>, NUM_TABLES> v_u = [&](u64 i) -> val<MAX_U_WIDTH> {
      arr<val<MAX_U_WIDTH>, AHEAD_BANKS> opts = [&](u64 b) -> val<MAX_U_WIDTH> {
        return ahead_u[1][b][i];
      };
      return opts.fo1().select(bank_sel);
    };

    // Write to regs (for update_cycle) and set fanout for val usage
    v_tag.fanout(hard<FETCH_WIDTH + 2>{}); // tag compare + reg write
    v_pred.fanout(hard<3>{});
    v_u.fanout(hard<2>{});
    for (u64 i = 0; i < NUM_TABLES; i++) {
      readt[i] = v_tag[i];
      readc[i] = v_pred[i];
      readu[i] = v_u[i];
    }
    arr<val<MAX_HYST_WIDTH>, NUM_TABLES> v_hyst =
        [&](u64 i) -> val<MAX_HYST_WIDTH> {
      arr<val<MAX_HYST_WIDTH>, AHEAD_BANKS> opts =
          [&](u64 b) -> val<MAX_HYST_WIDTH> { return ahead_hyst[1][b][i]; };
      return opts.fo1().select(bank_sel);
    };
    v_hyst.fanout(hard<2>{});
    for (u64 i = 0; i < NUM_TABLES; i++)
      readh[i] = v_hyst[i];

    val<NUM_TABLES> v_notumask = ~v_u.concat();
    v_notumask.fanout(hard<2>{});
    not_u_mask = v_notumask;

    // ---- Step 4: Tag compare + match + one_hot ----
    // Gather prediction bits (MSB of counter = direction prediction)
    val<NUM_TABLES> gpreds = [&]() -> val<NUM_TABLES> {
      if constexpr (MAX_CTR_WIDTH == 1) {
        return v_pred.concat();
      } else {
        arr<val<1>, NUM_TABLES> gp = [&](int i) -> val<1> {
          return v_pred[i] >> hard<MAX_CTR_WIDTH - 1>{};
        };
        return gp.fo1().concat();
      }
    }();
    gpreds.fanout(hard<FETCH_WIDTH>{});

    // Per-offset prediction vector: [bimodal, table0, table1, ..., tableN-1]
    arr<val<NUM_TABLES + 1>, FETCH_WIDTH> preds = [&](u64 offset) {
      return concat(readb[offset], gpreds);
    };
    preds.fanout(hard<2 * FETCH_WIDTH>{});

    // Tag compare: hashed portion — precomputed per bank in predict1.
    // Select the correct bank's precomputed result.
    ahead_htag[1].fanout(hard<3>{}); // update_cycle only now
    for (u64 b = 0; b < AHEAD_BANKS; b++)
      ahead_htagcmp[1][b].fanout(hard<2>{});
    arr<val<NUM_TABLES>, AHEAD_BANKS> htagcmp_opts = [&](u64 b) -> val<NUM_TABLES> {
      return ahead_htagcmp[1][b];
    };
    val<NUM_TABLES> htagcmp = htagcmp_opts.fo1().select(bank_sel);
    htagcmp.fanout(hard<FETCH_WIDTH>{});

    // BR_P_ENTRY=1: offset encoded in upper tag bits (still computed here)
    static_loop<FETCH_WIDTH>([&]<u64 offset>() {
      arr<val<1>, NUM_TABLES> tagcmp = [&](int i) {
        return val<LOG_FETCH_WIDTH>{v_tag[i] >> MAX_HTAGBITS} == hard<offset>{};
      };
      match[offset] = concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp);
    });
    match.fanout(hard<2>{});

    // Find longest match (provider) and second longest (alt)
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      match1[offset] = match[offset].one_hot();
    }
    match1.fanout(hard<6>{}); // predict2(3) + update_cycle(3)
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      tage_pred1[offset] = (match1[offset] & preds[offset]) != hard<0>{};
    }
    tage_pred1.fanout(hard<4>{}); // predict2(2) + update_cycle(2)

    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      match2[offset] = (match[offset] ^ match1[offset]).one_hot();
    }
    match2.fanout(hard<4>{}); // predict2(2) + update_cycle(2)
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      tage_pred2[offset] = (match2[offset] & preds[offset]) != hard<0>{};
    }
    tage_pred2.fanout(hard<4 + NUM_TABLES>{}); // predict2(2) + update_cycle(2+N)

    // ---- Step 5: Meta select ----
    // When provider is newly allocated (low confidence), meta counter decides
    // whether to trust primary (tage_pred1) or alternate (tage_pred2).
    if constexpr (bool(USE_META)) {
      meta.fanout(hard<4>{}); // predict2(2) + update_cycle(2)
      arr<val<1>, NUM_TABLES> weakctr = [&](int i) -> val<1> {
        return v_hyst[i] == hard<0>{};
      };
      val<NUM_TABLES> coldctr = v_notumask & weakctr.fo1().concat();
      coldctr.fanout(hard<FETCH_WIDTH>{});
      val<1> metasign = (meta[METAPIPE - 1] >= hard<0>{});
      metasign.fanout(hard<FETCH_WIDTH>{});
      for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
        newly_alloc[offset] = (match1[offset] & coldctr) != hard<0>{};
      }
      newly_alloc.fanout(hard<2>{});

      // Per-offset: if newly_alloc AND meta says use alt AND alt exists -> use
      // alt
      arr<val<1>, FETCH_WIDTH> altsel = [&](u64 offset) {
        arr<val<1>, 3> inputs = {metasign, newly_alloc[offset],
                                 match2[offset] != hard<0>{}};
        return inputs.fo1().fold_and();
      };

      // ---- Step 6: Final p2 = select(bank_hit, ahead_pred, bimodal) ----
      // OVERRIDER HOOK: auto ovr = overrider.lookup(inst_pc, tage_pred,
      // newly_alloc) Override decision is combinational on pre-extracted regs
      // from prefetch. Bake into p2: select(ovr.candidate, ovr.pred,
      // tage_or_bimodal)
      p2 = arr<val<1>, FETCH_WIDTH>{[&](u64 offset) {
             val<1> ahead = select(altsel[offset].fo1(), tage_pred2[offset],
                                   tage_pred1[offset]);
             return select(bank_hit, ahead, readb[offset]);
           }}.concat();
    } else {
      // No meta: primary prediction with bank_hit fallback
      p2 = arr<val<1>, FETCH_WIDTH>{[&](u64 offset) {
             return select(bank_hit, tage_pred1[offset], readb[offset]);
           }}.concat();
    }
    p2.fanout(hard<2 * FETCH_WIDTH + 3>{}); // predict2(2*FW+1) + update_cycle(2)
    val<1> taken = (blk_entry & p2) != hard<0>{};
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


    // ---- Fanout declarations ----
    mispredict.fanout(hard<NUM_TABLES + 2>{});
    val<1> correct_pred = ~mispredict;
    correct_pred.fanout(hard<NUM_TABLES + 2>{});
    // Fanouts for regs shared with predict2 are consolidated there.
    // Only declare fanout here for values unique to update_cycle.
    tage_indexes[1].fanout(hard<4>{});
    readt.fanout(hard<4>{});
    readc.fanout(hard<2>{});
    readh.fanout(hard<3>{});
    readu.fanout(hard<2>{});
    br_offset.fanout(hard<FETCH_WIDTH + NUM_TABLES + 1>{});
    branch_dir.fanout(hard<2 + AHEAD_BANK_BITS>{});


    // ---- Branch direction + update mask ----
    val<LOG_FETCH_WIDTH> last_offset = br_offset[num_branch - 1];
    last_offset.fanout(hard<4 * NUM_TABLES + 2>{});
    u64 update_valid = (u64(1) << num_branch) - 1;

    arr<val<FETCH_WIDTH>, FETCH_WIDTH> update_mask = [&](u64 offset) {
      arr<val<1>, FETCH_WIDTH> match_offset = [&](u64 i) {
        return br_offset[i] == offset;
      };
      return match_offset.fo1().concat() & update_valid;
    };
    update_mask.fanout(hard<2>{});

    arr<val<1>, FETCH_WIDTH> is_branch = [&](u64 offset) {
      return update_mask[offset] != hard<0>{};
    };
    is_branch.fanout(hard<6>{});

    val<FETCH_WIDTH> branch_mask = is_branch.concat();
    val<FETCH_WIDTH> actualdirs = branch_dir.concat();
    actualdirs.fanout(hard<FETCH_WIDTH>{});

    arr<val<1>, FETCH_WIDTH> branch_taken = [&](u64 offset) {
      return (actualdirs & update_mask[offset]) != hard<0>{};
    };
    branch_taken.fanout(hard<3>{});


    // ---- Provider identification ----
    arr<val<NUM_TABLES + 1>, FETCH_WIDTH> actual_match1 = [&](u64 offset) {
      return select(is_branch[offset], match1[offset], val<NUM_TABLES + 1>{0});
    };
    actual_match1.fanout(hard<2>{});

    val<NUM_TABLES> primary_mask = actual_match1.fold_or();
    primary_mask.fanout(hard<2>{});
    arr<val<1>, NUM_TABLES> primary = primary_mask.make_array(val<1>{});
    primary.fanout(hard<3>{});

    arr<val<1>, FETCH_WIDTH> primary_wrong = [&](u64 offset) {
      return tage_pred1[offset] != branch_taken[offset];
    };
    primary_wrong.fanout(hard<2>{});


    // ---- Allocation candidates ----
    val<NUM_TABLES> mispmask =
        mispredict.replicate(hard<NUM_TABLES>{}).concat();
    not_u_mask.fanout(hard<2>{});
    arr<val<1>, NUM_TABLES> last_tagcmp = [&](int i) {
      return readt[i] == concat(last_offset, ahead_htag[1][i]);
    };
    val<NUM_TABLES + 1> last_match1 =
        last_tagcmp.fo1().append(1).concat().one_hot();
    last_match1.fanout(hard<2>{});
    val<NUM_TABLES> postmask =
        mispmask.fo1() & val<NUM_TABLES>(last_match1 - 1);
    postmask.fanout(hard<2>{});
    val<NUM_TABLES> candallocmask = postmask & not_u_mask;
    candallocmask.fanout(hard<2>{});
    val<NUM_TABLES> collamask = candallocmask.reverse();
    collamask.fanout(hard<2>{});
    val<NUM_TABLES> collamask1 = collamask.one_hot();
    collamask1.fanout(hard<3>{});
    val<NUM_TABLES> collamask2 = (collamask ^ collamask1).one_hot();
    val<NUM_TABLES> collamask12 =
        select(val<2>{std::rand()} == hard<0>{}, collamask2.fo1(), collamask1);
    arr<val<1>, NUM_TABLES> allocate =
        collamask12.fo1().reverse().make_array(val<1>{});
    allocate.fanout(hard<7>{});


    // ---- Branch direction per table (BR_P_ENTRY=1) ----
    arr<val<1>, NUM_TABLES> bdir = [&](u64 i) {
      val<LOG_FETCH_WIDTH> tag_offset = readt[i] >> MAX_HTAGBITS;
      val<LOG_FETCH_WIDTH> offset =
          select(allocate[i], last_offset, tag_offset.fo1());
      offset.fanout(hard<FETCH_WIDTH>{});
      arr<val<1>, FETCH_WIDTH> match_offset = [&](u64 j) {
        return br_offset[j] == offset;
      };
      return (match_offset.fo1().concat() & update_valid & actualdirs) !=
             hard<0>{};
    };
    bdir.fanout(hard<2>{});

    // ---- Prediction correctness per table ----
    arr<val<1>, NUM_TABLES> badpred1 = [&](u64 i) -> val<1> {
      if constexpr (MAX_CTR_WIDTH == 1) {
        return readc[i] != bdir[i];
      } else {
        return val<1>{readc[i] >> hard<MAX_CTR_WIDTH - 1>{}} != bdir[i];
      }
    };
    badpred1.fanout(hard<3>{});

    // Does primary differ from alt?
    arr<val<1>, NUM_TABLES> altdiffer = [&](u64 i) -> val<1> {
      auto pred_dir = [&]() -> val<1> {
        if constexpr (MAX_CTR_WIDTH == 1) {
          return readc[i];
        } else {
          return readc[i] >> hard<MAX_CTR_WIDTH - 1>{};
        }
      }();
      val<LOG_FETCH_WIDTH> tag_offset = readt[i] >> MAX_HTAGBITS;
      return pred_dir != tage_pred2.select(tag_offset.fo1());
    };

    // Is owning branch's prediction correct?
    arr<val<1>, NUM_TABLES> goodpred = [&](u64 i) {
      val<LOG_FETCH_WIDTH> tag_offset = readt[i] >> MAX_HTAGBITS;
      return (tag_offset.fo1() != last_offset) | correct_pred;
    };
    goodpred.fanout(hard<2>{});


    // ---- Hysteresis weakness (for counter update gating) ----
    arr<val<1>, NUM_TABLES> g_weak = [&](u64 i) -> val<1> {
      return primary[i] & badpred1[i] & (readh[i] == hard<0>{});
    };
    g_weak.fanout(hard<2>{});

    // ---- P1 vs P2 disagreement ----
    val<FETCH_WIDTH> disagree_mask = (p1 ^ p2) & branch_mask.fo1();
    disagree_mask.fanout(hard<2>{});
    arr<val<1>, FETCH_WIDTH> disagree = disagree_mask.make_array(val<1>{});
    disagree.fanout(hard<2>{});

    // Read P1 hysteresis only if P1 and P2 disagree (silent: no read if agree)
    arr<val<1>, FETCH_WIDTH> p1_weak = [&](u64 offset) -> val<1> {
      return execute_if(disagree[offset],
                        [&]() { return ~p1_hys[offset].read(p1_index1); });
    };

    // Read bimodal hysteresis only if bimodal is provider and wrong
    arr<val<1>, FETCH_WIDTH> b_weak = [&](u64 offset) -> val<1> {
      val<1> bim_primary = actual_match1[offset] >> NUM_TABLES;
      return execute_if(bim_primary.fo1() & primary_wrong[offset],
                        [&]() { return ~bhyst[offset].read(bindex); });
    };

    // OVERRIDER HOOK: overrider.train_compute(...) before need_extra_cycle


    // ---- need_extra_cycle ----
    val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
    val<1> extra_cycle =
        some_badpred1.fo1() | mispredict | (disagree_mask != hard<0>{});
    extra_cycle.fanout(hard<NUM_TABLES * 2 + 1>{});
    need_extra_cycle(extra_cycle);

    // OVERRIDER HOOK: overrider.train_write() after need_extra_cycle

    // ---- Bank for writes (predecessor-indexed) ----
    val<AHEAD_BANK_BITS> write_bank_raw = path ^ XB[1];
    val<1> write_bank_valid = (write_bank_raw < hard<AHEAD_BANKS>{});
    val<AHEAD_BANK_BITS> write_bank =
        select(write_bank_valid, write_bank_raw, val<AHEAD_BANK_BITS>{0});
    write_bank.fanout(hard<NUM_TABLES * 4 + FETCH_WIDTH * 4>{});

    // ---- Meta counter update ----
    if constexpr (bool(USE_META)) {
      arr<val<1>, FETCH_WIDTH> altdiff = [&](u64 offset) {
        return (match2[offset] != hard<0>{}) &
               (tage_pred2[offset] != tage_pred1[offset]);
      };
      arr<val<2, i64>, FETCH_WIDTH> meta_incr = [&](u64 offset) -> val<2, i64> {
        val<1> update_meta =
            is_branch[offset] & altdiff[offset].fo1() & newly_alloc[offset];
        val<1> bad_pred2 = (tage_pred2[offset] != branch_taken[offset]);
        return select(update_meta.fo1(), concat(bad_pred2.fo1(), val<1>{1}),
                      val<2>{0});
      };
      for (u64 i = METAPIPE - 1; i != 0; i--)
        meta[i] = meta[i - 1];
      auto newmeta = meta[0] + meta_incr.fo1().fold_add();
      newmeta.fanout(hard<3>{});
      using meta_t = valt<decltype(meta[0])>;
      meta[0] = select(newmeta > meta_t::maxval, meta_t{meta_t::maxval},
                       select(newmeta < meta_t::minval, meta_t{meta_t::minval},
                              meta_t{newmeta}));
    }


    // ---- Tag write: allocation only (always dirty — new entry) ----
    static_loop<NUM_TABLES>([&]<u64 I>() {
      execute_if(allocate[I], [&]() {
        static_loop<AHEAD_BANKS>([&]<u64 b>() {
          execute_if(write_bank == hard<b>{}, [&]() {
            std::get<I>(tables[b]).tag_ram[0].write(
                tidx<I>(tage_indexes[1][I]),
                concat(last_offset, ahead_htag[1][I]));
          });
        });
      });
    });


    // ---- U-bit update: silent — only write when value changes ----
    arr<val<1>, NUM_TABLES> update_u = [&](u64 i) {
      return primary[i] & altdiffer[i].fo1();
    };
    val<1> noalloc = (candallocmask == hard<0>{});
    val<NUM_TABLES> uclearmask =
        postmask & noalloc.fo1().replicate(hard<NUM_TABLES>{}).concat();
    arr<val<1>, NUM_TABLES> uclear = uclearmask.fo1().make_array(val<1>{});
    uclear.fanout(hard<2>{});
    static_loop<NUM_TABLES>([&]<u64 I>() {
      val<1> newu = goodpred[I] & ~allocate[I] & ~uclear[I];
      newu.fanout(hard<AHEAD_BANKS + 1>{}); // inner loop + u_changed check
      // Silent: skip write if new u == old u
      val<1> u_changed = (newu != readu[I]) | allocate[I];
      execute_if((update_u[I] | allocate[I] | uclear[I]) & u_changed,
                 [&]() {
        static_loop<AHEAD_BANKS>([&]<u64 b>() {
          execute_if(write_bank == hard<b>{}, [&]() {
            std::get<I>(tables[b]).u_ram[0].write(tidx<I>(tage_indexes[1][I]),
                                                  newu, extra_cycle);
          });
        });
      });
    });


    // ---- P1 pred update: only when P1 disagrees with P2 and hyst is weak ----
    auto p2_split = p2.make_array(val<1>{});
    p2_split.fanout(hard<2>{});
    p1_weak.fanout(hard<2>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(p1_weak[offset], [&]() {
        p1_pred[offset].write(p1_index1, p2_split[offset]);
      });
    }
    // P1 hyst: write strong(1) if agree, weak(0) if disagree
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(is_branch[offset],
                 [&]() { p1_hys[offset].write(p1_index1, ~disagree[offset]); });
    }

    // ---- Bimodal pred: only write when wrong and hyst is weak ----
    b_weak.fanout(hard<2>{});
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      execute_if(b_weak[offset],
                 [&]() { bpred[offset].write(bindex, branch_taken[offset]); });
    }
    // Bimodal hyst: write if bimodal is provider
    for (u64 offset = 0; offset < FETCH_WIDTH; offset++) {
      val<1> bim_primary = match1[offset] >> NUM_TABLES;
      execute_if(is_branch[offset] & bim_primary.fo1(), [&]() {
        bhyst[offset].write(bindex, ~primary_wrong[offset]);
      });
    }

    // ---- TAGE counter: silent — only write when prediction was wrong + weak
    // hyst ---- g_weak already gates this: primary & badpred & weak_hyst
    // Allocation always writes (new entry with initial direction)
    static_loop<NUM_TABLES>([&]<u64 I>() {
      execute_if(g_weak[I] | allocate[I], [&]() {
        static_loop<AHEAD_BANKS>([&]<u64 b>() {
          execute_if(write_bank == hard<b>{}, [&]() {
            std::get<I>(tables[b]).pred_ram[0].write(
                tidx<I>(tage_indexes[1][I]), bdir[I]);
          });
        });
      });
    });

    // ---- TAGE hysteresis: silent — only write when value would change ----
    static constexpr u64 HMAX = (u64(1) << MAX_HYST_WIDTH) - 1;
    static_loop<NUM_TABLES>([&]<u64 I>() {
      // would_change: allocation OR (wrong pred & hyst not already 0)
      //               OR (correct pred & hyst not already max)
      val<1> would_change = allocate[I] |
                            (badpred1[I] & (readh[I] != hard<0>{})) |
                            (~badpred1[I] & (readh[I] != hard<HMAX>{}));
      execute_if((primary[I] | allocate[I]) & would_change, [&]() {
        auto newhyst = select(allocate[I], val<MAX_HYST_WIDTH>{0},
                              update_ctr(readh[I], ~badpred1[I]));
        newhyst.fanout(hard<AHEAD_BANKS>{});
        static_loop<AHEAD_BANKS>([&]<u64 b>() {
          execute_if(write_bank == hard<b>{}, [&]() {
            std::get<I>(tables[b]).hyst_ram[0].write(
                tidx<I>(tage_indexes[1][I]), newhyst, extra_cycle);
          });
        });
      });
    });

    // ---- U-bit epoch reset ----
    epoch_ctr.fanout(hard<3>{});
    val<NUM_TABLES> allocmask1 = collamask1.reverse();
    allocmask1.fanout(hard<2>{});
    val<1> faralloc =
        (((last_match1 >> 3) | allocmask1).one_hot() ^ allocmask1) == hard<0>{};
    val<1> uctrsat = (epoch_ctr == hard<decltype(epoch_ctr)::maxval>{});
    uctrsat.fanout(hard<2>{});
    epoch_ctr = select(correct_pred, epoch_ctr,
                       select(uctrsat, val<decltype(epoch_ctr)::size>{0},
                              update_ctr(epoch_ctr, faralloc.fo1())));
    // TODO: probabilistic decay (DECAY_CTR > 0) or epoch reset (DECAY_CTR == 0)

    // ---- Path update for next cycle's bank selection ----
    // path encodes which exit the block took:
    // 0 = fall-through (last branch not taken)
    // num_branch = Nth branch was taken
    path =
        num_branch &
        branch_dir[num_branch - 1].replicate(hard<AHEAD_BANK_BITS>{}).concat();
    // path fanout consolidated in predict2

    // ---- End-of-block history update ----
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
