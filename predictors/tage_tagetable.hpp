// TAGE predictor using TageTable for global table storage.
// Drop-in replacement for tage.hpp — same prediction logic, different storage backend.
// Pred(1 bit) + hyst(2 bits) encoded as 3-bit TageTable counter.

#define USE_META
#define RESET_UBITS

#include "../cbp.hpp"
#include "../harcom.hpp"
#include "common.hpp"
#include "custom/TageTable.hpp"

using namespace hcm;

template <u64 LOGLB = 6, u64 NUMG = 8, u64 LOGG = 11, u64 LOGB = 12,
          u64 TAGW = 11, u64 GHIST = 100, u64 LOGP1 = 14, u64 GHIST1 = 6>
struct tage_tt : predictor {
  // provides 2^(LOGLB-2) predictions per cycle
  // P2 is a TAGE, P1 is a gshare
  static_assert(LOGLB > 2);
  static_assert(NUMG > 0);
  static constexpr u64 MINHIST = 2;
  static constexpr u64 METABITS = 4;
  static constexpr u64 UCTRBITS = 8;
  static constexpr u64 PATHBITS = 6;
#ifdef USE_META
  static constexpr u64 METAPIPE = 2;
#endif
  static constexpr u64 LOGLINEINST = LOGLB - 2;
  static constexpr u64 LINEINST = 1 << LOGLINEINST;
  static_assert(LOGP1 > LOGLINEINST);
  static_assert(LOGB > LOGLINEINST);
  static constexpr u64 index1_bits = LOGP1 - LOGLINEINST;
  static constexpr u64 bindex_bits = LOGB - LOGLINEINST;
  static_assert(TAGW >
                LOGLINEINST); // the unhashed line offset is part of the tag
  static constexpr u64 HTAGBITS = TAGW - LOGLINEINST; // hashed tag bits

  // ======== TageTable type for global tables ========
  // Encodes 1-bit prediction + 2-bit hysteresis as 3-bit counter.
  // pred = counter[2] (MSB), hyst = counter[1:0].
  static constexpr u64 TT_CTR_WIDTH = 3;
  using GTable = TageTable<
      (1 << LOGG), // TABLE_SIZE
      0,           // TABLE_HIST (unused by table; predictor uses geometric_folds)
      TAGW,        // TAG_WIDTH
      TT_CTR_WIDTH,// CTR_WIDTH
      1,           // U_WIDTH
      1,           // N (1 prediction per entry)
      1,           // NUM_BANKS
      false,       // USE_AHEAD
      true,        // SHARED_TAG
      true,        // SHARED_U
      false,       // U_STOR_FF (SRAM mode — avoids FF write conflict with reset)
      1024,        // DECAY_CTR (probabilistic decay, replaces bulk reset)
      DefaultResetFn,
      false        // USE_FF_CACHE (BPB=1)
  >;

  geometric_folds<NUMG, MINHIST, GHIST, LOGG, HTAGBITS> gfolds;
  reg<1> true_block = 1;

  // for P1
  reg<GHIST1> global_history1;
  reg<index1_bits> index1;
  arr<reg<1>, LINEINST>
      readp1;       // prediction bits read from P1 table for each offset
  reg<LINEINST> p1; // P1 predictions

  // for P2
  reg<bindex_bits> bindex;       // bimodal table index
  arr<reg<LOGG>, NUMG> gindex;   // global tables indexes
  arr<reg<HTAGBITS>, NUMG> htag; // computed hashed tags

  arr<reg<1>, LINEINST> readb; // read bimodal prediction bit for each offset
  arr<reg<TAGW>, NUMG> readt;  // read tags (from TageTable)
  arr<reg<1>, NUMG> readc;     // read predictions (from TageTable counter MSB)
  arr<reg<2>, NUMG> readh;     // read hysteresis (from TageTable counter[1:0])
  arr<reg<1>, NUMG> readu;     // read u bits (from TageTable)
  reg<NUMG> notumask;          // read u bits, inverted

  arr<reg<NUMG + 1>, LINEINST> match;  // all matches for each offset
  arr<reg<NUMG + 1>, LINEINST> match1; // longest match for each offset
  arr<reg<NUMG + 1>, LINEINST> match2; // second longest match for each offset

  arr<reg<1>, LINEINST> pred1; // primary P2 prediction for each offset
  arr<reg<1>, LINEINST> pred2; // alternate P2 prediction for each offset
  reg<LINEINST> p2;            // final P2 predictions

#ifdef USE_META
  arr<reg<METABITS, i64>, METAPIPE> meta; // select between pred1 and pred2
  arr<reg<1>, LINEINST> newly_alloc;
#endif

#ifdef RESET_UBITS
  reg<UCTRBITS> uctr; // u bits counter (reset u bits when counter saturates)
#endif

  // simulation artifacts, hardware cost may not be real
  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<LOGLINEINST>, LINEINST> branch_offset;
  arr<reg<1>, LINEINST> branch_dir;
  reg<LINEINST> block_entry; // one-hot vector

  // P1 (gshare)
  ram<val<1>, (1 << index1_bits)> table1_pred[LINEINST]{
      "P1 pred"}; // P1 prediction bit

  // P2 (TAGE) — global tables replaced by TageTable
  GTable table[NUMG];

  // P2 — bimodal (unchanged)
  ram<val<1>, (1 << bindex_bits)> bim[LINEINST]{
      "bpred"}; // bimodal prediction bits

  zone UPDATE_ONLY;
  ram<val<1>, (1 << index1_bits)> table1_hyst[LINEINST]{
      "P1 hyst"}; // P1 hysteresis
  ram<val<1>, (1 << bindex_bits)> bhyst[LINEINST]{
      "bhyst"}; // bimodal hysteresis

  tage_tt() {
#ifdef TAGE_VERBOSE
    std::cerr << "TAGE_TT history lengths: ";
    for (u64 i = 0; i < NUMG; i++)
      std::cerr << gfolds.HLEN[i] << " ";
    std::cerr << std::endl;
    if (LOGG == HTAGBITS) {
      std::cerr << "WARNING: the tag function and index function are not "
                   "different enough\n";
    }
#endif
  }

  void new_block(val<64> inst_pc) {
    val<LOGLINEINST> offset = inst_pc.fo1() >> 2;
    block_entry = offset.fo1().decode().concat();
    block_entry.fanout(hard<6 * LINEINST>{});
    block_size = 1;
  }

  val<1> predict1([[maybe_unused]] val<64> inst_pc) {
    inst_pc.fanout(hard<2>{});
    new_block(inst_pc);
    val<std::max(index1_bits, GHIST1)> lineaddr = inst_pc >> LOGLB;
    lineaddr.fanout(hard<2>{});
    if constexpr (GHIST1 <= index1_bits) {
      index1 = lineaddr ^
               (val<index1_bits>{global_history1} << (index1_bits - GHIST1));
    } else {
      index1 = global_history1.make_array(val<index1_bits>{})
                   .append(lineaddr)
                   .fold_xor();
    }
    index1.fanout(hard<LINEINST>{});
    for (u64 offset = 0; offset < LINEINST; offset++) {
      readp1[offset] = table1_pred[offset].read(index1);
    }
    readp1.fanout(hard<2>{});
    p1 = readp1.concat();
    p1.fanout(hard<LINEINST>{});
    return (block_entry & p1) != hard<0>{};
  };

  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) {
    return ((block_entry << block_size) & p1) != hard<0>{};
  };

  val<1> predict2(val<64> inst_pc) {
    val<std::max(bindex_bits, LOGG)> lineaddr = inst_pc >> LOGLB;
    lineaddr.fanout(hard<1 + NUMG * 2>{});
    gfolds.fanout(hard<2>{});

    // compute indexes
    bindex = lineaddr;
    bindex.fanout(hard<LINEINST>{});
    for (u64 i = 0; i < NUMG; i++) {
      gindex[i] = lineaddr ^ gfolds.template get<0>(i);
    }
    gindex.fanout(hard<2>{}); // TageTable: 1 read + 1 write (was 4+4 with separate RAMs)

    // compute hashed tags
    for (u64 i = 0; i < NUMG; i++) {
      htag[i] = val<HTAGBITS>{lineaddr}.reverse() ^ gfolds.template get<1>(i);
    }
    htag.fanout(hard<2>{});

    std::cerr << "DBG: predict2 reading tables\n";
    // read tables — TageTable replaces 4 separate RAM reads
    for (u64 offset = 0; offset < LINEINST; offset++) {
      readb[offset] = bim[offset].read(bindex);
    }
    readb.fanout(hard<2>{});
    for (u64 i = 0; i < NUMG; i++) {
      // Read TageTable (comparison tag unused — predictor does its own matching)
      table[i].read(val<LOGG>{gindex[i]}, val<TAGW>{0}, 0);
      readt[i] = table[i].getTag();
      val<TT_CTR_WIDTH> ctr_val = table[i].getCounter(0, 0);
      auto [hyst_val, pred_val] = split<2, 1>(ctr_val);
      readc[i] = pred_val;
      readh[i] = hyst_val;
      readu[i] = table[i].getU();
    }
    readt.fanout(hard<LINEINST + 1>{});
    readc.fanout(hard<3>{});
    readh.fanout(hard<2>{});
    readu.fanout(hard<2>{});
    notumask = ~readu.concat();
    notumask.fanout(hard<2>{});

    // gather prediction bits for each offset
    val<NUMG> gpreds = readc.concat();
    gpreds.fanout(hard<LINEINST>{});
    arr<val<NUMG + 1>, LINEINST> preds = [&](u64 offset) {
      return concat(readb[offset], gpreds);
    };
    preds.fanout(hard<2 * LINEINST>{});

    // hashed tags comparisons
    arr<val<1>, NUMG> htagcmp_split = [&](int i) {
      return val<HTAGBITS>{readt[i]} == htag[i];
    };
    val<NUMG> htagcmp = htagcmp_split.fo1().concat();
    htagcmp.fanout(hard<LINEINST>{});

    // generate match mask for each offset
    static_loop<LINEINST>([&]<u64 offset>() {
      arr<val<1>, NUMG> tagcmp = [&](int i) {
        return val<LOGLINEINST>{readt[i] >> HTAGBITS} == hard<offset>{};
      };
      match[offset] =
          concat(val<1>{1}, tagcmp.fo1().concat() &
                                htagcmp); // bimodal is default when no match
    });
    match.fanout(hard<2>{});

    // for each offset, find longest match and select primary prediction
    for (u64 offset = 0; offset < LINEINST; offset++) {
      match1[offset] = match[offset].one_hot();
    }
    match1.fanout(hard<3>{});
    for (u64 offset = 0; offset < LINEINST; offset++) {
      pred1[offset] = (match1[offset] & preds[offset]) != hard<0>{};
    }
    pred1.fanout(hard<2>{});

    // for each offset, find second longest match and select secondary
    // prediction
    for (u64 offset = 0; offset < LINEINST; offset++) {
      match2[offset] = (match[offset] ^ match1[offset]).one_hot();
    }
    match2.fanout(hard<2>{});
    for (u64 offset = 0; offset < LINEINST; offset++) {
      pred2[offset] = (match2[offset] & preds[offset]) != hard<0>{};
    }
    pred2.fanout(hard<2>{});

#ifdef USE_META
    meta.fanout(hard<2>{});
    arr<val<1>, NUMG> weakctr = [&](int i) { return readh[i] == hard<0>{}; };
    val<NUMG> coldctr = notumask & weakctr.fo1().concat();
    coldctr.fanout(hard<LINEINST>{});
    val<1> metasign = (meta[METAPIPE - 1] >= hard<0>{});
    metasign.fanout(hard<LINEINST>{});
    for (u64 offset = 0; offset < LINEINST; offset++) {
      newly_alloc[offset] = (match1[offset] & coldctr) != hard<0>{};
    }
    newly_alloc.fanout(hard<2>{});
    arr<val<1>, LINEINST> altsel = [&](u64 offset) {
      arr<val<1>, 3> inputs = {metasign, newly_alloc[offset],
                               match2[offset] != hard<0> {}};
      return inputs.fo1().fold_and();
    };
    p2 = arr<val<1>, LINEINST>{[&](u64 offset) {
           return select(altsel[offset].fo1(), pred2[offset], pred1[offset]);
         }}.concat();
#else
    p2 = pred1.concat();
#endif
    p2.fanout(hard<LINEINST>{});
    val<1> taken = (block_entry & p2) != hard<0>{};
    taken.fanout(hard<2>{});
    reuse_prediction(~val<1>{block_entry >> (LINEINST - 1)});
    return taken;
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) {
    val<1> taken = ((block_entry << block_size) & p2) != hard<0>{};
    taken.fanout(hard<2>{});
    reuse_prediction(~val<1>{block_entry >> (LINEINST - 1 - block_size)});
    block_size++;
    return taken;
  }

  void update_condbr(val<64> branch_pc, val<1> taken,
                     [[maybe_unused]] val<64> next_pc) {
    assert(num_branch < LINEINST);
    branch_offset[num_branch] = branch_pc.fo1() >> 2;
    branch_dir[num_branch] = taken.fo1();
    num_branch++;
  }

    void update_cycle(instruction_info &block_end_info)
    {
        std::cerr << "DBG: update_cycle start, num_branch=" << num_branch << "\n";
        val<1> &mispredict = block_end_info.is_mispredict;
        val<64> &next_pc = block_end_info.next_pc;
        // updates for all conditional branches in the predicted block
        if (num_branch == 0) {
            // no conditional branch in this block
            val<1> line_end = block_entry >> (LINEINST-block_size);
            // update global history if previous block ended on a mispredicted not-taken branch
            // (we are still in the same line, this is the last chunk)
            // or if the block ends before the line boundary (unconditional jump)
            val<1> actual_block = ~(true_block & line_end.fo1());
            actual_block.fanout(hard<GHIST+NUMG*2+2>{});
            execute_if(actual_block, [&](){
                next_pc.fanout(hard<2>{});
                global_history1 = (global_history1 << 1) ^ val<GHIST1>{next_pc>>2};
                gfolds.update(val<PATHBITS>{next_pc>>2});
                true_block = 1;
            });
            return; // stop here
        }
        mispredict.fanout(hard<NUMG+2>{});
        val<1> correct_pred = ~mispredict;
        correct_pred.fanout(hard<NUMG+2>{});
        index1.fanout(hard<LINEINST*3>{});
        p2.fanout(hard<2>{});
        bindex.fanout(hard<LINEINST*3>{});
        gindex.fanout(hard<2>{}); // TageTable: 1 write per table (was 4 with separate RAMs)
        htag.fanout(hard<3>{});
        readb.fanout(hard<2>{});
        readt.fanout(hard<4>{});
        readc.fanout(hard<2>{});
        readh.fanout(hard<3>{});
        match1.fanout(hard<3>{});
        match2.fanout(hard<2>{});
        pred1.fanout(hard<2>{});
        pred2.fanout(hard<2+NUMG>{});
        branch_offset.fanout(hard<LINEINST+NUMG+1>{});
        branch_dir.fanout(hard<2>{});
        gfolds.fanout(hard<2>{});
        readu.fanout(hard<2>{});
#ifdef USE_META
        meta.fanout(hard<2>{});
#endif
        std::cerr << "DBG: before last_offset\n";
        val<LOGLINEINST> last_offset = branch_offset[num_branch-1];
        last_offset.fanout(hard<4*NUMG+2>{});

        u64 update_valid = (u64(1)<<num_branch)-1;
        arr<val<LINEINST>,LINEINST> update_mask = [&](u64 offset){
            arr<val<1>,LINEINST> match_offset = [&](u64 i){return branch_offset[i] == offset;};
            return match_offset.fo1().concat() & update_valid;
        };
        update_mask.fanout(hard<2>{});

        arr<val<1>,LINEINST> is_branch = [&](u64 offset){
            return update_mask[offset] != hard<0>{};
        };
        is_branch.fanout(hard<6>{});

        val<LINEINST> branch_mask = is_branch.concat();

        val<LINEINST> actualdirs = branch_dir.concat();
        actualdirs.fanout(hard<LINEINST>{});

        arr<val<1>,LINEINST> branch_taken = [&](u64 offset){
            return (actualdirs & update_mask[offset]) != hard<0>{};
        };
        branch_taken.fanout(hard<3>{});

        arr<val<NUMG+1>,LINEINST> actual_match1 = [&] (u64 offset) {
            return select(is_branch[offset],match1[offset],val<NUMG+1>{0});
        };
        actual_match1.fanout(hard<2>{});

        val<NUMG> primary_mask = actual_match1.fold_or();
        primary_mask.fanout(hard<2>{});
        arr<val<1>,NUMG> primary = primary_mask.make_array(val<1>{});
        primary.fanout(hard<3>{});

        arr<val<1>,LINEINST> primary_wrong = [&](u64 offset){
            return pred1[offset] != branch_taken[offset];
        };
        primary_wrong.fanout(hard<2>{});

        std::cerr << "DBG: before allocation\n";
        // select some candidate entries for allocation
        val<NUMG> mispmask = mispredict.replicate(hard<NUMG>{}).concat();
        std::cerr << "DBG: alloc A\n";
        arr<val<1>,NUMG> last_tagcmp = [&](int i){return readt[i] == concat(last_offset,htag[i]);};
        std::cerr << "DBG: alloc B\n";
        val<NUMG+1> last_match1 = last_tagcmp.fo1().append(1).concat().one_hot();
        std::cerr << "DBG: alloc C\n";
        last_match1.fanout(hard<2>{});
        val<NUMG> postmask = mispmask.fo1() & val<NUMG>(last_match1-1);
        std::cerr << "DBG: alloc D\n";
        postmask.fanout(hard<2>{});
        val<NUMG> candallocmask = postmask & notumask; // candidate post entries for allocation
        candallocmask.fanout(hard<2>{});
        std::cerr << "DBG: alloc E\n";
        // if multiple candidate entries, we select a single one, with some randomization
        val<NUMG> collamask = candallocmask.reverse();
        collamask.fanout(hard<2>{});
        val<NUMG> collamask1 = collamask.one_hot();
        collamask1.fanout(hard<3>{});
        std::cerr << "DBG: alloc F\n";
        val<NUMG> collamask2 = (collamask^collamask1).one_hot();
        val<NUMG> collamask12 = select(val<2>{std::rand()}==hard<0>{}, collamask2.fo1(), collamask1);
        arr<val<1>,NUMG> allocate = collamask12.fo1().reverse().make_array(val<1>{});
        allocate.fanout(hard<10>{}); // increased fanout for combined TageTable write
        std::cerr << "DBG: alloc G\n";

        std::cerr << "DBG: before bdir\n";
        // associate a branch direction to each global table
        arr<val<1>,NUMG> bdir = [&](u64 i) {
            std::cerr << "DBG: bdir i=" << i << "\n";
            val<LOGLINEINST> tag_offset = readt[i] >> HTAGBITS;
            val<LOGLINEINST> offset = select(allocate[i],last_offset,tag_offset.fo1());
            offset.fanout(hard<LINEINST>{});
            arr<val<1>,LINEINST> match_offset = [&](u64 j){return branch_offset[j] == offset;};
            return (match_offset.fo1().concat() & update_valid & actualdirs) != hard<0>{};
        };
        std::cerr << "DBG: after bdir\n";
        bdir.fanout(hard<2>{});

        // tell if global prediction is incorrect
        std::cerr << "DBG: before badpred1\n";
        arr<val<1>,NUMG> badpred1 = [&](u64 i){
            std::cerr << "DBG: badpred1 i=" << i << "\n";
            return readc[i] != bdir[i];
        };
        std::cerr << "DBG: badpred1 constructed, calling fanout\n";
        badpred1.fanout(hard<3>{});
        std::cerr << "DBG: badpred1 fanout done\n";

        // associate to each global table a bit telling if local prediction differs from secondary prediction
        std::cerr << "DBG: before altdiffer\n";
        arr<val<1>,NUMG> altdiffer = [&](u64 i){
            std::cerr << "DBG: altdiffer i=" << i << "\n";
            val<LOGLINEINST> tag_offset = readt[i] >> HTAGBITS;
            return readc[i] != pred2.select(tag_offset.fo1());
        };

        // associate to each global table a bit telling if prediction for owning branch is correct
        std::cerr << "DBG: before goodpred\n";
        arr<val<1>,NUMG> goodpred = [&](u64 i){
            std::cerr << "DBG: goodpred i=" << i << "\n";
            val<LOGLINEINST> tag_offset = readt[i] >> HTAGBITS;
            return (tag_offset.fo1() != last_offset) | correct_pred;
        };
        std::cerr << "DBG: before goodpred.fanout\n";
        goodpred.fanout(hard<2>{});
        std::cerr << "DBG: after goodpred.fanout\n";

        // do P1 and P2 agree?
        std::cerr << "DBG: p1 ^ p2\n";
        val<LINEINST> p1xp2 = p1 ^ p2;
        std::cerr << "DBG: branch_mask.fo1()\n";
        val<LINEINST> bm = branch_mask.fo1();
        std::cerr << "DBG: p1xp2 & bm\n";
        val<LINEINST> disagree_mask = p1xp2 & bm;
        std::cerr << "DBG: disagree_mask.fanout\n";
        disagree_mask.fanout(hard<2>{});
        std::cerr << "DBG: disagree_mask.make_array\n";
        arr<val<1>,LINEINST> disagree = disagree_mask.make_array(val<1>{});
        std::cerr << "DBG: disagree.fanout\n";
        disagree.fanout(hard<2>{});

        // read the P1 hysteresis if P1 and P2 disagree
        std::cerr << "DBG: before p1_weak\n";
        arr<val<1>,LINEINST> p1_weak = [&] (u64 offset) -> val<1> {
            // returns 1 iff disagreement and hysteresis is weak
            return execute_if(disagree[offset], [&](){
                return ~table1_hyst[offset].read(index1); // hyst=0 means weak
            });
        };

        // read the bimodal hysteresis if bimodal caused a misprediction
        std::cerr << "DBG: before b_weak\n";
        arr<val<1>,LINEINST> b_weak = [&] (u64 offset) -> val<1> {
            // returns 1 iff cause of misprediction and hysteresis is weak
            val<1> bim_primary = actual_match1[offset] >> NUMG;
            return execute_if(bim_primary.fo1() & primary_wrong[offset], [&](){
                return ~bhyst[offset].read(bindex); // hyst=0 means weak
            });
        };

        // determine which primary global predictions are incorrect with a weak hysteresis
        std::cerr << "DBG: before g_weak\n";
        arr<val<1>,NUMG> g_weak = [&] (u64 i) -> val<1> {
            // returns 1 iff incorrect primary prediction and hysteresis is weak
            return primary[i] & badpred1[i] & (readh[i]==hard<0>{});
        };
        g_weak.fanout(hard<2>{});

        // need extra cycle for modifying prediction bits and for TAGE allocation
        std::cerr << "DBG: before some_badpred1\n";
        val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
        val<1> extra_cycle = some_badpred1.fo1() | mispredict | (disagree_mask != hard<0>{});
        extra_cycle.fanout(hard<NUMG*2+1>{});
        need_extra_cycle(extra_cycle);

#ifdef USE_META
    // update meta counter
    arr<val<1>, LINEINST> altdiff = [&](u64 offset) {
      // for each offset, tell if primary and secondary predictions differ
      return (match2[offset] != hard<0>{}) & (pred2[offset] != pred1[offset]);
    };
    arr<val<2, i64>, LINEINST> meta_incr = [&](u64 offset) -> val<2, i64> {
      val<1> update_meta =
          is_branch[offset] & altdiff[offset].fo1() & newly_alloc[offset];
      val<1> bad_pred2 = (pred2[offset] != branch_taken[offset]);
      return select(update_meta.fo1(), concat(bad_pred2.fo1(), val<1>{1}),
                    val<2>{0});
    };
    for (u64 i = METAPIPE - 1; i != 0; i--) {
      meta[i] = meta[i - 1];
    }
    auto newmeta = meta[0] + meta_incr.fo1().fold_add();
    newmeta.fanout(hard<3>{});
    using meta_t = valt<decltype(meta[0])>;
    meta[0] = select(newmeta > meta_t::maxval, meta_t{meta_t::maxval},
                     select(newmeta < meta_t::minval, meta_t{meta_t::minval},
                            meta_t{newmeta}));
#endif

        std::cerr << "DBG: before combined write\n";
        // ======== Combined TageTable write ========
        // Reference tage writes tag, pred, hyst, u separately with different
        // conditions. With TageTable, we combine into a single write per table,
        // using selects to choose new vs old values for each field.

        // u-bit update logic
        arr<val<1>,NUMG> update_u = [&](u64 i){
            return primary[i] & altdiffer[i].fo1();
        };
        update_u.fanout(hard<2>{});
        // if all post entries have the u bit set, reset their u bits
        val<1> noalloc = (candallocmask == hard<0>{});
        val<NUMG> uclearmask = postmask & noalloc.fo1().replicate(hard<NUMG>{}).concat();
        arr<val<1>,NUMG> uclear = uclearmask.fo1().make_array(val<1>{});
        uclear.fanout(hard<3>{});

        for (u64 i = 0; i < NUMG; i++) {
            std::cerr << "DBG: combined write i=" << i << " need_write\n";
            val<1> need_write = allocate[i] | g_weak[i] | primary[i] |
                                update_u[i] | uclear[i];
            std::cerr << "DBG: combined write i=" << i << " execute_if\n";
            execute_if(need_write, [&]() {
                std::cerr << "DBG: write i=" << i << " wtag\n";
                // tag: new on allocate, old otherwise
                val<TAGW> wtag = select(allocate[i],
                    concat(last_offset, val<HTAGBITS>{htag[i]}),
                    val<TAGW>{readt[i]});

                std::cerr << "DBG: write i=" << i << " wpred\n";
                // pred: new on g_weak/allocate, old otherwise
                val<1> wpred = select(g_weak[i] | allocate[i],
                    val<1>{bdir[i]}, val<1>{readc[i]});

                std::cerr << "DBG: write i=" << i << " whyst\n";
                // hyst: new on primary/allocate, old otherwise
                val<2> whyst = select(primary[i] | allocate[i],
                    select(allocate[i], val<2>{0},
                           update_ctr(val<2>{readh[i]}, ~badpred1[i])),
                    val<2>{readh[i]});

                std::cerr << "DBG: write i=" << i << " wu\n";
                // u: new on update_u/allocate/uclear, old otherwise
                val<1> wu = select(update_u[i] | allocate[i] | uclear[i],
                    val<1>{goodpred[i] & ~allocate[i] & ~uclear[i]},
                    val<1>{readu[i]});

                std::cerr << "DBG: write i=" << i << " wctr+write\n";
                val<TT_CTR_WIDTH> wctr = concat(whyst, wpred);
                table[i].write(val<LOGG>{gindex[i]}, 0, 0, wtag, wctr, wu);
                std::cerr << "DBG: write i=" << i << " done\n";
            });
        }

        std::cerr << "DBG: before P1 update\n";
        // update P1 prediction if P1 and P2 disagree and the hysteresis bit is weak
        auto p2_split = p2.make_array(val<1>{});
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(p1_weak[offset].fo1(), [&](){
                // update with the P2 prediction, not with the actual branch direction
                table1_pred[offset].write(index1,p2_split[offset].fo1());
            });
        }
        // update P1 hysteresis
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(is_branch[offset],[&](){
                table1_hyst[offset].write(index1,~disagree[offset]);
            });
        }

        // update incorrect bimodal prediction if primary provider and hysteresis is weak
        for (u64 offset=0; offset<LINEINST; offset++) {
            execute_if(b_weak[offset].fo1(), [&](){
                bim[offset].write(bindex,branch_taken[offset]);
            });
        }
        // update bimodal hysteresis if bimodal is primary provider
        for (u64 offset=0; offset<LINEINST; offset++) {
            val<1> bim_primary = match1[offset] >> NUMG;
            execute_if(is_branch[offset] & bim_primary.fo1(), [&](){
                bhyst[offset].write(bindex,~primary_wrong[offset]);
            });
        }

        std::cerr << "DBG: before RESET_UBITS\n";
        // With U_STOR_FF=false (SRAM mode), u-bit decay is handled internally
        // by TageTable via DECAY_CTR probabilistic clearing on tag misses.
        // The bulk reset logic (uctr/uctrsat/reset_u) is not needed.

        std::cerr << "DBG: before history update\n";
        // update global history
        val<1> line_end = block_entry >> (LINEINST-block_size);
        true_block = correct_pred | branch_dir[num_branch-1] | line_end.fo1();
        true_block.fanout(hard<GHIST+NUMG*2+2>{});
        execute_if(true_block, [&](){
            next_pc.fanout(hard<2>{});
            global_history1 = (global_history1 << 1) ^ val<GHIST1>{next_pc>>2};
            gfolds.update(val<PATHBITS>{next_pc>>2});
        });

        num_branch = 0; // done
    }
};
