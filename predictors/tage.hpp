// this is a basic TAGE, not necessarily well optimized

#define USE_META
#define RESET_UBITS

#include "../cbp.hpp"
#include "../harcom.hpp"
#include "common.hpp"

#ifdef TAGE_MONITOR
#include <array>
#include <bitset>
#include <iostream>
#include <iomanip>
#include <unordered_set>
#endif

#if defined(TAGE_MONITOR) || defined(DUMP_SITES) || defined(TAGE_VERBOSE)
#include <fstream>
#include <iostream>
#endif

using namespace hcm;

template <u64 LOGLB = 6, u64 NUMG = 8, u64 LOGG = 11, u64 LOGB = 12,
          u64 TAGW = 11, u64 GHIST = 100, u64 LOGP1 = 14, u64 GHIST1 = 6,
          bool SHARED_HYS = false>
struct tage : predictor {
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
  arr<reg<TAGW>, NUMG> readt;  // read tags
  arr<reg<1>, NUMG> readc;     // read predictions
  arr<reg<2>, NUMG> readh;     // read hysteresis
  arr<reg<1>, NUMG> readu;     // read u bits
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

  // P2 (TAGE)
  ram<val<TAGW>, (1 << LOGG)> gtag[NUMG]{"tags"}; // tags
  ram<val<1>, (1 << LOGG)> gpred[NUMG]{"gpred"};  // predictions
  rwram<2, (1 << LOGG), 4> ghyst[NUMG]{"ghyst"};  // hysteresis
  rwram<1, (1 << LOGG), 4> ubit[NUMG]{"u"};       // "useful" bits
  ram<val<1>, (1 << bindex_bits)> bim[LINEINST]{
      "bpred"}; // bimodal prediction bits

  zone UPDATE_ONLY;
  ram<val<1>, (1 << index1_bits)> table1_hyst[LINEINST]{
      "P1 hyst"}; // P1 hysteresis
  ram<val<1>, (1 << bindex_bits)> bhyst[LINEINST]{
      "bhyst"}; // bimodal hysteresis

#ifdef TAGE_MONITOR
  struct TageMonitor {
    static constexpr u64 NT = NUMG;
    static constexpr u64 LI = LINEINST;
    static constexpr u64 WINDOW = 100000;

    struct Counters {
      u64 branches = 0, mispredictions = 0, blocks = 0, extra_cycles = 0;
      std::array<u64, NT + 1> provider_count{};   // NT = bimodal
      std::array<u64, NT + 1> provider_correct{};
      std::array<u64, NT> tag_lookups{}, tag_matches{};
      u64 p1_predictions = 0, p1_correct = 0;
      u64 p1p2_disagree = 0, p1_right_p2_wrong = 0, p2_right_p1_wrong = 0;
      u64 alloc_attempts = 0, alloc_success = 0, alloc_blocked = 0;
      std::array<u64, NT> alloc_per_table{};
      std::array<u64, NT + 1> alloc_from_provider{};
      std::array<std::array<u64, NT>, NT + 1> alloc_cascade{};
      u64 epoch_resets = 0, uctr_sum = 0, uctr_samples = 0;
      void reset() { *this = Counters{}; }
    };
    Counters cum, win;
    u64 window_num = 0;
    bool header_printed = false;
    std::unordered_set<u64> unique_branch_pcs;

    // Shadow state per offset (set in predict2/update_cycle)
    std::array<u64, LI> shadow_provider{};
    std::array<bool, LI> shadow_p1{}, shadow_p2{};

    static u64 decode_prov(u64 one_hot) {
      for (u64 i = 0; i <= NT; i++)
        if (one_hot & (u64(1) << i)) return i;
      return NT;
    }
    static double pct(u64 n, u64 d) { return d > 0 ? 100.0 * n / d : 0.0; }

    void record_prediction(u64 offset, u64 match1_val, bool p1, bool p2) {
      shadow_provider[offset] = decode_prov(match1_val);
      shadow_p1[offset] = p1;
      shadow_p2[offset] = p2;
    }

    void record_tag_lookup(u64 table, bool matched) {
      cum.tag_lookups[table]++; win.tag_lookups[table]++;
      if (matched) { cum.tag_matches[table]++; win.tag_matches[table]++; }
    }

    void record_outcome(u64 offset, bool actual, bool is_last_br, bool mispredict) {
      bool p2_ok = (shadow_p2[offset] == actual);
      bool p1_ok = (shadow_p1[offset] == actual);
      auto rec = [&](Counters &c) {
        c.branches++;
        if (is_last_br && mispredict) c.mispredictions++;
        u64 prov = shadow_provider[offset];
        c.provider_count[prov]++;
        if (p2_ok) c.provider_correct[prov]++;
        c.p1_predictions++;
        if (p1_ok) c.p1_correct++;
        if (shadow_p1[offset] != shadow_p2[offset]) {
          c.p1p2_disagree++;
          if (p1_ok && !p2_ok) c.p1_right_p2_wrong++;
          if (p2_ok && !p1_ok) c.p2_right_p1_wrong++;
        }
      };
      rec(cum); rec(win);
      if (win.branches >= WINDOW) { print_window(std::cerr); win.reset(); window_num++; }
    }

    void record_block(bool extra) {
      cum.blocks++; win.blocks++;
      if (extra) { cum.extra_cycles++; win.extra_cycles++; }
    }

    void record_alloc(u64 provider_idx, u64 alloc_mask) {
      auto rec = [&](Counters &c) {
        c.alloc_attempts++;
        if (alloc_mask) {
          c.alloc_success++;
          c.alloc_from_provider[provider_idx]++;
          for (u64 i = 0; i < NT; i++)
            if (alloc_mask & (u64(1) << i)) { c.alloc_per_table[i]++; c.alloc_cascade[provider_idx][i]++; }
        } else {
          c.alloc_blocked++;
        }
      };
      rec(cum); rec(win);
    }

    void record_uctr(u64 v) { cum.uctr_sum += v; cum.uctr_samples++; win.uctr_sum += v; win.uctr_samples++; }
    void record_epoch_reset() { cum.epoch_resets++; win.epoch_resets++; }
    void record_branch_pc(u64 pc) { unique_branch_pcs.insert(pc); }

    void print_window(std::ostream &os) {
      if (!header_printed) {
        os << "# win,br,misp%,extra%,bim%,";
        for (u64 i = 0; i < NT; i++) os << "T" << i << "%,";
        os << "p1acc%,p1p2dis%,alloc_ok%,uctr_avg\n";
        header_printed = true;
      }
      auto &w = win;
      os << std::fixed << std::setprecision(1);
      os << window_num << "," << w.branches << "," << pct(w.mispredictions, w.branches) << ","
         << pct(w.extra_cycles, w.blocks) << "," << pct(w.provider_count[NT], w.branches) << ",";
      for (u64 i = 0; i < NT; i++) os << pct(w.provider_count[i], w.branches) << ",";
      os << pct(w.p1_correct, w.p1_predictions) << ","
         << pct(w.p1p2_disagree, w.branches) << ","
         << pct(w.alloc_success, w.alloc_attempts) << ","
         << (w.uctr_samples > 0 ? double(w.uctr_sum) / w.uctr_samples : 0) << "\n";
    }

    void print_summary(std::ostream &os = std::cerr) const {
      const auto &c = cum;
      os << "\n=== Tage<> Monitor Summary ===\n";
      os << "Branches: " << c.branches << "  Mispredictions: " << c.mispredictions
         << " (" << std::fixed << std::setprecision(2) << pct(c.mispredictions, c.branches) << "%)\n";
      os << "Blocks: " << c.blocks << "  Extra cycles: " << c.extra_cycles
         << " (" << pct(c.extra_cycles, c.blocks) << "%)\n";
      os << "Unique branch PCs: " << unique_branch_pcs.size() << "\n\n";

      os << "Provider Distribution:\n";
      os << "  Table     | Provided  |  Correct  | Accuracy | TagMatch%\n";
      os << "  ----------+-----------+-----------+----------+----------\n";
      os << "  Bimodal   |" << std::setw(10) << c.provider_count[NT]
         << " |" << std::setw(10) << c.provider_correct[NT]
         << " |" << std::setw(7) << std::setprecision(1)
         << pct(c.provider_correct[NT], c.provider_count[NT]) << "% |      --\n";
      for (u64 i = 0; i < NT; i++) {
        os << "  T" << i << std::setw(8 - (i >= 10 ? 2 : 1)) << ""
           << "|" << std::setw(10) << c.provider_count[i]
           << " |" << std::setw(10) << c.provider_correct[i]
           << " |" << std::setw(7) << pct(c.provider_correct[i], c.provider_count[i]) << "%"
           << " |" << std::setw(7) << pct(c.tag_matches[i], c.tag_lookups[i]) << "%\n";
      }
      u64 tage_total = 0, tage_correct = 0;
      for (u64 i = 0; i < NT; i++) { tage_total += c.provider_count[i]; tage_correct += c.provider_correct[i]; }
      os << "  TAGE total: " << tage_total << " (" << pct(tage_total, c.branches) << "% of branches)"
         << "  Accuracy: " << pct(tage_correct, tage_total) << "%\n";

      os << "\nP1 (gshare):\n";
      os << "  Accuracy: " << c.p1_correct << "/" << c.p1_predictions
         << " (" << pct(c.p1_correct, c.p1_predictions) << "%)\n";
      os << "  Disagree with P2: " << c.p1p2_disagree
         << " (" << pct(c.p1p2_disagree, c.branches) << "%)\n";
      os << "  P1 right P2 wrong: " << c.p1_right_p2_wrong
         << "  P2 right P1 wrong: " << c.p2_right_p1_wrong << "\n";

      os << "\nAllocation:\n";
      os << "  Attempts: " << c.alloc_attempts
         << "  Success: " << c.alloc_success
         << " (" << pct(c.alloc_success, c.alloc_attempts) << "%)"
         << "  Blocked: " << c.alloc_blocked << "\n";
      os << "  Per table:";
      for (u64 i = 0; i < NT; i++) os << " T" << i << "=" << c.alloc_per_table[i];
      os << "\n";

      os << "\nAllocation Cascade (provider -> target):\n";
      os << "  Provider  | Allocs  |";
      for (u64 i = 0; i < NT; i++) os << std::setw(7) << "T" + std::to_string(i);
      os << "\n  ----------+---------+";
      for (u64 i = 0; i < NT; i++) os << "-------";
      os << "\n";
      for (u64 p = 0; p <= NT; p++) {
        if (c.alloc_from_provider[p] == 0) continue;
        if (p == NT) os << "  Bimodal   |";
        else os << "  T" << p << std::setw(8 - (p >= 10 ? 2 : 1)) << "" << "|";
        os << std::setw(8) << c.alloc_from_provider[p] << " |";
        for (u64 t = 0; t < NT; t++) os << std::setw(7) << c.alloc_cascade[p][t];
        os << "\n";
      }
      os << "\n  Epoch resets: " << c.epoch_resets
         << "  Avg uctr: " << (c.uctr_samples > 0 ? double(c.uctr_sum) / c.uctr_samples : 0) << "\n";
      os << "=== End Tage<> Monitor ===\n\n";
    }
  };

  TageMonitor mon;
  ~tage() { mon.print_summary(); }
#endif

  tage() {
#ifdef TAGE_VERBOSE
    std::cerr << "TAGE history lengths: ";
    for (u64 i = 0; i < NUMG; i++)
      std::cerr << gfolds.HLEN[i] << " ";
    std::cerr << std::endl;
    if (LOGG == HTAGBITS) {
      std::cerr << "WARNING: the tag function and index function are not "
                   "different enough\n";
    }
#endif
#ifdef DUMP_SITES
    std::ofstream f("sites.txt");
    f << "# register/value name → HARCOM site() (= RAM ID)\n";
    f << "true_block " << true_block.hcm_location() << "\n";
    f << "global_history1 " << global_history1.hcm_location() << "\n";
    f << "index1 " << index1.hcm_location() << "\n";
    f << "p1 " << p1.hcm_location() << "\n";
    f << "bindex " << bindex.hcm_location() << "\n";
    f << "notumask " << notumask.hcm_location() << "\n";
    f << "p2 " << p2.hcm_location() << "\n";
    f << "uctr " << uctr.hcm_location() << "\n";
    f << "block_entry " << block_entry.hcm_location() << "\n";
    for (u64 i = 0; i < LINEINST; i++) {
      f << "readp1[" << i << "] " << readp1[i].hcm_location() << "\n";
      f << "readb[" << i << "] " << readb[i].hcm_location() << "\n";
      f << "match[" << i << "] " << match[i].hcm_location() << "\n";
      f << "match1[" << i << "] " << match1[i].hcm_location() << "\n";
      f << "match2[" << i << "] " << match2[i].hcm_location() << "\n";
      f << "pred1[" << i << "] " << pred1[i].hcm_location() << "\n";
      f << "pred2[" << i << "] " << pred2[i].hcm_location() << "\n";
      f << "newly_alloc[" << i << "] " << newly_alloc[i].hcm_location() << "\n";
      f << "branch_offset[" << i << "] " << branch_offset[i].hcm_location() << "\n";
      f << "branch_dir[" << i << "] " << branch_dir[i].hcm_location() << "\n";
    }
    for (u64 i = 0; i < NUMG; i++) {
      f << "gindex[" << i << "] " << gindex[i].hcm_location() << "\n";
      f << "htag[" << i << "] " << htag[i].hcm_location() << "\n";
      f << "readt[" << i << "] " << readt[i].hcm_location() << "\n";
      f << "readc[" << i << "] " << readc[i].hcm_location() << "\n";
      f << "readh[" << i << "] " << readh[i].hcm_location() << "\n";
      f << "readu[" << i << "] " << readu[i].hcm_location() << "\n";
    }
    for (u64 i = 0; i < METAPIPE; i++)
      f << "meta[" << i << "] " << meta[i].hcm_location() << "\n";
    f.close();
    std::cerr << "[DUMP_SITES] wrote sites.txt\n";
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
    gindex.fanout(hard<4>{});

    // compute hashed tags
    for (u64 i = 0; i < NUMG; i++) {
      htag[i] = val<HTAGBITS>{lineaddr}.reverse() ^ gfolds.template get<1>(i);
    }
    htag.fanout(hard<2>{});

    // read tables
    for (u64 offset = 0; offset < LINEINST; offset++) {
      readb[offset] = bim[offset].read(bindex);
    }
    readb.fanout(hard<2>{});
    for (u64 i = 0; i < NUMG; i++) {
      readt[i] = gtag[i].read(gindex[i]);
      readc[i] = gpred[i].read(gindex[i]);
      readh[i] = ghyst[i].read(gindex[i]);
      readu[i] = ubit[i].read(gindex[i]);
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
#ifdef TAGE_MONITOR
    mon.record_branch_pc(static_cast<u64>(branch_pc));
#endif
    num_branch++;
  }

    void update_cycle(instruction_info &block_end_info)
    {
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
        gindex.fanout(hard<4>{});
        htag.fanout(hard<3>{});
        readb.fanout(hard<2>{});
        readt.fanout(hard<4>{});
        readc.fanout(hard<2>{});
        match1.fanout(hard<3>{});
        match2.fanout(hard<2>{});
        pred1.fanout(hard<2>{});
        pred2.fanout(hard<2+NUMG>{});
        branch_offset.fanout(hard<LINEINST+NUMG+1>{});
        branch_dir.fanout(hard<2>{});
        gfolds.fanout(hard<2>{});
#ifdef USE_META
        meta.fanout(hard<2>{});
#endif
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

        // select some candidate entries for allocation
        val<NUMG> mispmask = mispredict.replicate(hard<NUMG>{}).concat();
        arr<val<1>,NUMG> last_tagcmp = [&](int i){return readt[i] == concat(last_offset,htag[i]);};
        val<NUMG+1> last_match1 = last_tagcmp.fo1().append(1).concat().one_hot();
        last_match1.fanout(hard<2>{});
        val<NUMG> postmask = mispmask.fo1() & val<NUMG>(last_match1-1);
        postmask.fanout(hard<2>{});
        val<NUMG> candallocmask = postmask & notumask; // candidate post entries for allocation
        candallocmask.fanout(hard<2>{});
        // if multiple candidate entries, we select a single one, with some randomization
        val<NUMG> collamask = candallocmask.reverse();
        collamask.fanout(hard<2>{});
        val<NUMG> collamask1 = collamask.one_hot();
        collamask1.fanout(hard<3>{});
        val<NUMG> collamask2 = (collamask^collamask1).one_hot();
        val<NUMG> collamask12 = select(val<2>{std::rand()}==hard<0>{}, collamask2.fo1(), collamask1);
        arr<val<1>,NUMG> allocate = collamask12.fo1().reverse().make_array(val<1>{});
        allocate.fanout(hard<7>{});

        // associate a branch direction to each global table
        arr<val<1>,NUMG> bdir = [&](u64 i) {
            val<LOGLINEINST> tag_offset = readt[i] >> HTAGBITS;
            val<LOGLINEINST> offset = select(allocate[i],last_offset,tag_offset.fo1());
            offset.fanout(hard<LINEINST>{});
            arr<val<1>,LINEINST> match_offset = [&](u64 j){return branch_offset[j] == offset;};
            return (match_offset.fo1().concat() & update_valid & actualdirs) != hard<0>{};
        };
        bdir.fanout(hard<2>{});

        // tell if global prediction is incorrect
        arr<val<1>,NUMG> badpred1 = [&](u64 i){
            return readc[i] != bdir[i];
        };
        badpred1.fanout(hard<3>{});

        // associate to each global table a bit telling if local prediction differs from secondary prediction
        arr<val<1>,NUMG> altdiffer = [&](u64 i){
            val<LOGLINEINST> tag_offset = readt[i] >> HTAGBITS;
            return readc[i] != pred2.select(tag_offset.fo1());
        };

        // associate to each global table a bit telling if prediction for owning branch is correct
        arr<val<1>,NUMG> goodpred = [&](u64 i){
            val<LOGLINEINST> tag_offset = readt[i] >> HTAGBITS;
            return (tag_offset.fo1() != last_offset) | correct_pred;
        };

        // do P1 and P2 agree?
        val<LINEINST> disagree_mask = (p1 ^ p2) & branch_mask.fo1();
        disagree_mask.fanout(hard<2>{});
        arr<val<1>,LINEINST> disagree = disagree_mask.make_array(val<1>{});
        disagree.fanout(hard<2>{});

        // read the P1 hysteresis if P1 and P2 disagree
        arr<val<1>,LINEINST> p1_weak = [&] (u64 offset) -> val<1> {
            // returns 1 iff disagreement and hysteresis is weak
            return execute_if(disagree[offset], [&](){
                return ~table1_hyst[offset].read(index1); // hyst=0 means weak
            });
        };

        // read the bimodal hysteresis if bimodal caused a misprediction
        arr<val<1>,LINEINST> b_weak = [&] (u64 offset) -> val<1> {
            // returns 1 iff cause of misprediction and hysteresis is weak
            val<1> bim_primary = actual_match1[offset] >> NUMG;
            return execute_if(bim_primary.fo1() & primary_wrong[offset], [&](){
                return ~bhyst[offset].read(bindex); // hyst=0 means weak
            });
        };

        // determine which primary global predictions are incorrect with a weak hysteresis
        arr<val<1>,NUMG> g_weak = [&] (u64 i) -> val<1> {
            // returns 1 iff incorrect primary prediction and hysteresis is weak
            return primary[i] & badpred1[i] & (readh[i]==hard<0>{});
        };

        // need extra cycle for modifying prediction bits and for TAGE allocation
        val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
        val<1> extra_cycle = some_badpred1.fo1() | mispredict | (disagree_mask != hard<0>{});
        extra_cycle.fanout(hard<NUMG*2+1>{});
        need_extra_cycle(extra_cycle);

#ifdef TAGE_MONITOR
        mon.record_block(static_cast<u64>(extra_cycle));
        // Per-table tag match tracking (htag match, not full tag with offset)
        for (u64 i = 0; i < NUMG; i++) {
          bool matched = static_cast<u64>(val<HTAGBITS>{readt[i]} == htag[i]);
          mon.record_tag_lookup(i, matched);
        }
        // Per-offset predictions and outcomes
        for (u64 offset = 0; offset < LINEINST; offset++) {
          if (!static_cast<u64>(is_branch[offset])) continue;
          bool p1_pr = static_cast<u64>(readp1[offset]); // P1 gshare prediction for this offset
          bool p2_pr = static_cast<u64>((p2 >> offset) & hard<1>{});
          mon.record_prediction(offset, static_cast<u64>(match1[offset]), p1_pr, p2_pr);
          bool actual = static_cast<u64>(branch_taken[offset]);
          bool is_last = (offset == static_cast<u64>(last_offset));
          mon.record_outcome(offset, actual, is_last, static_cast<u64>(mispredict));
        }
        // Allocation tracking
        if (static_cast<u64>(mispredict)) {
          u64 lm1 = static_cast<u64>(last_match1);
          u64 prov_idx = NUMG; // default = bimodal
          for (u64 i = 0; i < NUMG; i++)
            if (lm1 & (u64(1) << i)) { prov_idx = i; break; }
          u64 amask = 0;
          for (u64 i = 0; i < NUMG; i++)
            if (static_cast<u64>(allocate[i])) amask |= (u64(1) << i);
          mon.record_alloc(prov_idx, amask);
        }
#endif

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

        // overwrite the tag in the allocated entry (mispredict)
        for (u64 i=0; i<NUMG; i++) {
            execute_if(allocate[i], [&](){gtag[i].write(gindex[i],concat(last_offset,htag[i]));});
        }

        // update the u bits
        arr<val<1>,NUMG> update_u = [&](u64 i){
            return primary[i] & altdiffer[i].fo1();
        };
        // if all post entries have the u bit set, reset their u bits
        val<1> noalloc = (candallocmask == hard<0>{});
        val<NUMG> uclearmask = postmask & noalloc.fo1().replicate(hard<NUMG>{}).concat();
        arr<val<1>,NUMG> uclear = uclearmask.fo1().make_array(val<1>{});
        uclear.fanout(hard<2>{});
        for (u64 i=0; i<NUMG; i++) {
            execute_if(update_u[i].fo1() | allocate[i] | uclear[i], [&]() {
                val<1> newu = goodpred[i].fo1() & ~allocate[i] & ~uclear[i];
                ubit[i].write(gindex[i],newu.fo1(),extra_cycle);
            });
        }

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

        // update incorrect global prediction if primary provider and the hysteresis is weak;
        // initialize global prediction in the allocated entry
        for (u64 i=0; i<NUMG; i++) {
            execute_if(g_weak[i].fo1() | allocate[i], [&](){
                gpred[i].write(gindex[i],bdir[i]);
            });
        }
        // update global prediction hysteresis if primary provider or allocated entry
        for (u64 i=0; i<NUMG; i++) {
            execute_if(primary[i] | allocate[i], [&](){
                // if allocated entry, set hysteresis to 0;
                // otherwise, increment hysteresis if correct pred, decrement if incorrect
                val<2> newhyst = select(allocate[i],val<2>{0},update_ctr(readh[i],~badpred1[i]));
                ghyst[i].write(gindex[i],newhyst.fo1(),extra_cycle);
            });
        }

#ifdef RESET_UBITS
    uctr.fanout(hard<3>{});
    val<NUMG> allocmask1 = collamask1.reverse();
    allocmask1.fanout(hard<2>{});
    val<1> faralloc =
        (((last_match1 >> 3) | allocmask1).one_hot() ^ allocmask1) == hard<0>{};
    val<1> uctrsat = (uctr == hard<decltype(uctr)::maxval>{});
    uctrsat.fanout(hard<2>{});
    uctr = select(correct_pred, uctr,
                  select(uctrsat, val<decltype(uctr)::size>{0},
                         update_ctr(uctr, faralloc.fo1())));
    execute_if(uctrsat, [&]() {
      for (auto &uram : ubit)
        uram.reset();
    });
#ifdef TAGE_MONITOR
    mon.record_uctr(static_cast<u64>(uctr));
    if (static_cast<u64>(uctrsat)) mon.record_epoch_reset();
#endif
#endif

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
