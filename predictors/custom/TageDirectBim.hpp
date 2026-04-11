#pragma once

#include "TageDirect.hpp"

// ============================================================================
// TageDirectBim — Offset-tagged TAGE with bimodal P1 base predictor.
//
// Fork of TageDirect: bimodal replaces gshare so TAGE tables become
// the primary provider, enabling the allocation cascade.
// TAGE tags use static PC-derived offset bits (inst_pc >> 2) instead of
// dynamic branch rank, so the same instruction always matches the same
// tag regardless of block boundaries.  Mirrors predictors/tage.hpp's
// concat(offset, htag) scheme.
// ============================================================================

// ============================================================================
// TageDirectBimImpl
// ============================================================================
template <typename TableCfg, typename AllocCfg, u64 LINEINST_V, u64 N_V,
          u64 DECAY_CTR_V, u64 DECAY_GRAN_V,
          typename DecayPolicy_V, u64 P1_TABLE_SIZE_V,
          bool USE_META_V, u64 METABITS_V, u64 METAPIPE_V,
          bool USE_PATH_HIST_V, u64 PATH_HIST_WIDTH_V, u64 PATH_BITS_V,
          template <u64> class FoldFn_V,
          u64 RWRAM_BANKS_V = 4, u64 RWRAM_BANK_SHIFT_V = 0,
          u64 EPOCH_CTR_BITS_V = 18,
          bool SHARED_HYS_V = false,
          bool USE_DIR_HIST_V = false>
struct TageDirectBimImpl : predictor {
  using AllocTrigger = td::AllocTrigger;
  using AllocAction  = td::AllocAction;
  using UctrPolicy   = td::UctrPolicy;

  // ======== Constants ========
  static constexpr u64 NUM_TABLES = TableCfg::NUM_TABLES;
  static constexpr u64 N = N_V;
  static constexpr u64 LANES = td::next_pow2(N);
  static constexpr u64 LOG_LANES = td::clog2(LANES);
  static constexpr u64 LINEINST = LINEINST_V;
  static constexpr u64 LOG_LINEINST = td::clog2(LINEINST);

  static constexpr u64 MAX_TAG_WIDTH = td::array_max(TableCfg::TAG_WIDTH);
  static constexpr u64 MAX_CTR_WIDTH = td::array_max(TableCfg::CTR_WIDTH);
  static constexpr u64 MAX_HYST_WIDTH = td::array_max(TableCfg::HYST_WIDTH);
  static constexpr u64 MAX_U_WIDTH = td::array_max(TableCfg::U_WIDTH);
  static constexpr u64 MAX_TABLE_SIZE = td::array_max(TableCfg::TABLE_SIZE);
  static constexpr u64 MAX_IDX_BITS = td::clog2(MAX_TABLE_SIZE);
  static constexpr u64 MAXHIST = td::array_max(TableCfg::HIST_LEN);
  static constexpr u64 MAX_HTAGBITS = MAX_TAG_WIDTH - LOG_LANES;

  static constexpr u64 P1_ENTRIES = P1_TABLE_SIZE_V;
  static constexpr u64 P1_INDEX_BITS = td::clog2(P1_TABLE_SIZE_V);
  static constexpr u64 LINEADDR_BITS = MAX_IDX_BITS;
  static constexpr u64 EPOCH_CTR_BITS = EPOCH_CTR_BITS_V;
  static constexpr u64 PATHBITS = PATH_BITS_V;

  static constexpr bool USE_PROB_DECAY = (DECAY_CTR_V > 0);

  // ---- Table tuple type ----
  using Tables =
      typename td::TDMakeTableTuple<TableCfg, RWRAM_BANKS_V, RWRAM_BANK_SHIFT_V,
                                    SHARED_HYS_V,
                                    std::make_index_sequence<NUM_TABLES>>::type;

  // Truncate gindex to per-table IDX_BITS
  template <u64 I> auto tidx(auto &gi) {
    using Table = std::tuple_element_t<I, Tables>;
    return val<Table::IDX_BITS>{gi};
  }

  // Truncate gindex to per-table HYST_IDX_BITS (= IDX_BITS or IDX_BITS-1 if SHARED_HYS)
  template <u64 I> auto hidx(auto &gi) {
    using Table = std::tuple_element_t<I, Tables>;
    return val<Table::HYST_IDX_BITS>{gi};
  }

  // ======== Static asserts ========
  static_assert(NUM_TABLES > 0);
  static_assert(N >= 1);
  static_assert(LANES >= N);
  static_assert(std::has_single_bit(LINEINST));
  static_assert(MAX_TAG_WIDTH > LOG_LANES, "TAG_WIDTH must be > LOG_LANES");

  // ======== State ========

  // History
  td::TDHistoryFolder<NUM_TABLES, MAXHIST, MAX_IDX_BITS, MAX_HTAGBITS, FoldFn_V,
                      TableCfg::HIST_LEN>
      gfolds;
  bool gfolds_inited = false;
  reg<1> true_block = 1;

  // P1 bimodal
  reg<P1_INDEX_BITS> index1;
  reg<LANES> unordered_pred1;
  arr<reg<1>, LANES> pred; // extracted per-lane P1 predictions

  // P2 TAGE
  arr<reg<MAX_IDX_BITS>, NUM_TABLES> gindex;
  arr<reg<MAX_HTAGBITS>, NUM_TABLES> htag;

  arr<reg<MAX_TAG_WIDTH>, NUM_TABLES> readt;
  arr<reg<MAX_CTR_WIDTH>, NUM_TABLES> readc;
  arr<reg<std::max(u64(1), MAX_HYST_WIDTH)>, NUM_TABLES> readh;
  arr<reg<MAX_U_WIDTH>, NUM_TABLES> readu;
  reg<NUM_TABLES> notumask;

  arr<reg<1>, NUM_TABLES> htagcmp_reg;
  arr<reg<1>, NUM_TABLES> last_tagcmp_reg;

  // Per-rank match/prediction
  arr<reg<NUM_TABLES + 1>, LANES> match;
  arr<reg<NUM_TABLES + 1>, LANES> match1;
  arr<reg<NUM_TABLES + 1>, LANES> match2;
  arr<reg<1>, LANES> pred1_tage;
  arr<reg<1>, LANES> pred2_tage;
  reg<LANES> p2;

  // Meta
  std::conditional_t<USE_META_V, arr<reg<METABITS_V, i64>, METAPIPE_V>,
                     EmptyMember>
      meta;
  std::conditional_t<USE_META_V, arr<reg<1>, LANES>, EmptyMember> newly_alloc;

  // U-bit reset
  reg<EPOCH_CTR_BITS> uctr;
  std::conditional_t<USE_PROB_DECAY, reg<DECAY_CTR_V == 0 ? 1 : DECAY_CTR_V>,
                     EmptyMember>
      decay_threshold;

  // Allocation pressure tracking
  static constexpr u64 ALLOC_PRESS_W = AllocCfg::ALLOC_PRESSURE_BITS;
  static constexpr u64 ACC_PRESS_W   = AllocCfg::ACCURACY_PRESSURE_BITS;
  std::conditional_t<(ALLOC_PRESS_W > 0),
      reg<ALLOC_PRESS_W == 0 ? 1 : ALLOC_PRESS_W>, EmptyMember> alloc_pressure;
  std::conditional_t<(ACC_PRESS_W > 0),
      reg<ACC_PRESS_W == 0 ? 1 : ACC_PRESS_W>, EmptyMember> accuracy_pressure;

  // Path history
  std::conditional_t<USE_PATH_HIST_V, reg<PATH_HIST_WIDTH_V>, EmptyMember>
      path_hist;

  // Block tracking (gshareN_ahead_best pattern: numeric offset)
  reg<LOG_LINEINST> block_entry;
  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<1>, LANES> branch_dir;     // direction per branch rank (0..num_branch-1)
  arr<reg<LOG_LANES>, LANES> branch_offset; // PC-derived offset per branch rank
  // is_branch computed from hardware in update_cycle (no C++ bitmask needed)
  reg<N + 1> rank;  // one-hot: rank of current branch (P1 lane scrambling)
  reg<LANES> X;     // lane scrambling register (P1 bimodal)

  // P1 storage
  hcm::ram<val<LANES>, P1_ENTRIES> p1_pred{"P1 pred"};
  zone UPDATE_ONLY;
  hcm::ram<val<1>, P1_ENTRIES> p1_hyst[LANES]{"P1 hyst"};

  // TAGE tables
  Tables tables;

#ifdef TAGE_MONITOR
  TDMonitor<NUM_TABLES, LANES, P1_ENTRIES, RWRAM_BANKS_V, MAX_TABLE_SIZE> mon;
  // Shadow state for meta tracking (set in predict2, read in update_cycle)
  std::array<bool, LANES> mon_altsel{};     // meta chose alt for this rank
  std::array<bool, LANES> mon_meta_active{}; // meta override was active
  ~TageDirectBimImpl() { mon.print_summary(); }
#endif

  bool params_printed = false;

  // ======== Helpers ========

  // End block at LINEINST boundary
  val<1> line_end() { return (block_entry + block_size) >= hard<LINEINST>{}; }

  val<1> last_pred() {
    return val<1>{num_branch >= N};
  }

  // XOR fold intra-line PC bits into LOG_LANES bits for stable lane assignment
  static u64 pc_hash_lane(u64 pc) {
    u64 lo = (pc >> 2) & (LANES - 1);
    u64 hi = (pc >> (2 + LOG_LANES)) & (LANES - 1);
    return lo ^ hi;
  }
  // pc_hash_lane(val<64>) not used — val→u64 cast requires cheating mode

  void ensure_gfolds_init() {
    if (!gfolds_inited) {
      std::array<u64, NUM_TABLES> hl;
      for (u64 i = 0; i < NUM_TABLES; i++)
        hl[i] = TableCfg::HIST_LEN[i];
      gfolds.init(hl);
      gfolds_inited = true;
    }
  }

  void print_params(std::ostream &os) const {
    os << "=== TageDirect Parameters ===\n";
    os << "NUM_TABLES=" << NUM_TABLES << "  LINEINST=" << LINEINST
       << "  N=" << N << "  LANES=" << LANES << "\n";
    os << "Table sizes: ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::TABLE_SIZE[i];
    os << "\nTag widths:  ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::TAG_WIDTH[i];
    os << "\nHist lens:   ";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << (i ? "," : "") << TableCfg::HIST_LEN[i];
    os << "\nP1: BIMODAL  TABLE_SIZE=" << P1_TABLE_SIZE_V;
    os << "\nMETA: USE=" << USE_META_V << "  BITS=" << METABITS_V
       << "  PIPE=" << METAPIPE_V;
    os << "\nMISPREDICT_ONLY_WRITE=" << AllocCfg::MISPREDICT_ONLY_WRITE;
    os << "\nAlloc: TRIGGER=" << td::to_string(AllocCfg::ALLOC_TRIGGER)
       << "  ACTION=" << td::to_string(AllocCfg::ALLOC_ACTION)
       << "  UCTR=" << td::to_string(AllocCfg::UCTR_POLICY);
    os << "\n  MAX_ALLOC=" << AllocCfg::MAX_ALLOC
       << "  NON_CONSECUTIVE=" << AllocCfg::NON_CONSECUTIVE
       << "  CONF_GATE=" << AllocCfg::CONF_GATE
       << "  PROB_START=" << AllocCfg::PROB_START
       << "  PARTIAL_UPDATE=" << AllocCfg::PARTIAL_UPDATE;
    os << "\n  ALLOC_PRESSURE_BITS=" << AllocCfg::ALLOC_PRESSURE_BITS
       << "  ACCURACY_PRESSURE_BITS=" << AllocCfg::ACCURACY_PRESSURE_BITS;
    os << "\n  DISAGREE_EXTRA_CYCLE=" << AllocCfg::DISAGREE_EXTRA_CYCLE;
    os << "\nSHARED_HYS=" << SHARED_HYS_V;
    os << "\n\n";
  }

  // ======== Predictor Interface ========

  val<1> predict1(val<64> inst_pc) override {
    ensure_gfolds_init();
    if (!params_printed) {
      print_params(std::cerr);
      params_printed = true;
    }
#ifdef TAGE_MONITOR
    mon.record_predict1();
#endif

    inst_pc.fanout(hard<3>{});
    true_block.fanout(hard<4>{});

    block_entry = select(true_block, val<LOG_LINEINST>{inst_pc >> 2},
                         val<LOG_LINEINST>{block_entry + block_size});
    block_entry.fanout(hard<LANES + 2>{});

    // Lane scrambling
    rank = select(true_block, val<N + 1>{1}, rank << num_branch);
    rank.fanout(hard<N + 2>{});
    X = select(true_block, val<LOG_LANES>{inst_pc >> 2}.decode().concat(),
               X.rotate_left(num_branch));
    X.fanout(hard<LANES>{});

    // P1 bimodal: read on true blocks only
    execute_if(true_block, [&]() {
      index1 = inst_pc >> 2;
      unordered_pred1 = p1_pred.read(index1);
    });

    // Extract per-lane predictions
    unordered_pred1.fanout(hard<LANES>{});
    for (u64 i = 0; i < LANES; i++) {
      pred[i] = (unordered_pred1 & X.rotate_left(i)) != hard<0>{};
    }
    pred.fanout(hard<LANES * 2>{});

    block_size = 1;
    num_branch = 0;
    reuse_prediction(~line_end());
    return pred[0];
  }

  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) override {
#ifdef TAGE_MONITOR
    mon.record_reuse_predict1();
#endif
    block_size++;
    reuse_prediction(~line_end());
    return pred[num_branch];
  }

  val<1> predict2(val<64> inst_pc) override {
#ifdef TAGE_MONITOR
    mon.record_predict2();
#endif
    // Line-level address: group LANES instructions per line (like tage.hpp's >> LOGLB)
    // Offset within line is stored in the tag via LOG_LANES bits of (inst_pc >> 2)
    val<LINEADDR_BITS> lineaddr = inst_pc >> (2 + LOG_LANES);
    lineaddr.fanout(hard<1 + NUM_TABLES * 2>{});
    gfolds.fanout(hard<2>{});

    // Compute TAGE indices
    for (u64 i = 0; i < NUM_TABLES; i++) {
      if constexpr (USE_PATH_HIST_V) {
        gindex[i] =
            lineaddr ^ gfolds.get_idx_fold(i) ^ val<MAX_IDX_BITS>{path_hist};
      } else {
        gindex[i] = lineaddr ^ gfolds.get_idx_fold(i);
      }
    }
    gindex.fanout(hard<4>{});

    // Read TAGE tables
    static_loop<NUM_TABLES>([&]<u64 I>() {
      auto &t = std::get<I>(tables);
#ifdef TD_VERBOSE
      std::cerr << "P2: tage tag_ram[" << I << "] READ\n";
#endif
      readt[I] = t.tag_ram.read(tidx<I>(gindex[I]));
#ifdef TD_VERBOSE
      std::cerr << "P2: tage pred_ram[" << I << "] READ\n";
#endif
      readc[I] = t.pred_ram.read(tidx<I>(gindex[I]));
      if constexpr (MAX_HYST_WIDTH > 0) {
#ifdef TD_VERBOSE
        t.hyst_ram.debug_read_info("hyst_ram", I, hidx<I>(gindex[I]));
#endif
        readh[I] = t.hyst_ram.read(hidx<I>(gindex[I]));
#ifdef TAGE_MONITOR
        mon.record_rwram_read(I);
#endif
      }
#ifdef TD_VERBOSE
      t.u_ram.debug_read_info("u_ram", I, tidx<I>(gindex[I]));
#endif
      readu[I] = t.u_ram.read(tidx<I>(gindex[I]));
#ifdef TAGE_MONITOR
      mon.record_rwram_read(I);
#endif
    });
    readt.fanout(hard<LANES + 1>{});
    readc.fanout(hard<3>{});
    if constexpr (MAX_HYST_WIDTH > 0)
      readh.fanout(hard<2>{});
    readu.fanout(hard<2>{});
    notumask = ~readu.concat();
    notumask.fanout(hard<2>{});

    // Compute hashed tags (parallel with RAM reads)
    for (u64 i = 0; i < NUM_TABLES; i++) {
      htag[i] = val<MAX_HTAGBITS>{lineaddr}.reverse() ^ gfolds.get_tag_fold(i);
    }
    htag.fanout(hard<2>{});

    // Gather prediction bits per table
    val<NUM_TABLES> gpreds = [&]() -> val<NUM_TABLES> {
      if constexpr (MAX_CTR_WIDTH == 1) {
        return readc.concat();
      } else {
        arr<val<1>, NUM_TABLES> gp = [&](int i) -> val<1> {
          return readc[i] >> hard<MAX_CTR_WIDTH - 1>{};
        };
        return gp.fo1().concat();
      }
    }();
    gpreds.fanout(hard<LANES>{});

    // Per-rank preds: P1 gshare as fallback (position NUM_TABLES = bimodal
    // slot)
    pred.fanout(hard<LANES>{});
    arr<val<NUM_TABLES + 1>, LANES> preds = [&](u64 r) {
      return concat(pred[r], gpreds);
    };
    preds.fanout(hard<2 * LANES>{});

    // Per-table htag comparison
    static_loop<NUM_TABLES>([&]<u64 I>() {
      using Table = std::tuple_element_t<I, Tables>;
      static constexpr u64 PER_HTAG = Table::tag_width - LOG_LANES;
      htagcmp_reg[I] = (val<PER_HTAG>{readt[I] >> LOG_LANES} == val<PER_HTAG>{htag[I]});
    });
    htagcmp_reg.fanout(hard<LANES + 1>{});

    // Per-rank tag match
    static_loop<LANES>([&]<u64 R>() {
      arr<val<1>, NUM_TABLES> tagcmp = [&](int i) {
        return val<LOG_LANES>{readt[i]} == hard<R>{};
      };
      match[R] =
          concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp_reg.concat());
    });
    match.fanout(hard<2>{});

    // Provider/alt selection per rank
    for (u64 r = 0; r < LANES; r++) {
      match1[r] = match[r].one_hot();
    }
    match1.fanout(hard<3>{});
    for (u64 r = 0; r < LANES; r++) {
      pred1_tage[r] = (match1[r] & preds[r]) != hard<0>{};
    }
    pred1_tage.fanout(hard<2>{});

    for (u64 r = 0; r < LANES; r++) {
      match2[r] = (match[r] ^ match1[r]).one_hot();
    }
    match2.fanout(hard<2>{});
    for (u64 r = 0; r < LANES; r++) {
      pred2_tage[r] = (match2[r] & preds[r]) != hard<0>{};
    }
    pred2_tage.fanout(hard<2>{});

    // Meta prediction
    if constexpr (USE_META_V) {
      meta.fanout(hard<2>{});
      arr<val<1>, NUM_TABLES> weakctr = [&](int i) {
        return readh[i] == hard<0>{};
      };
      val<NUM_TABLES> coldctr = notumask & weakctr.fo1().concat();
      coldctr.fanout(hard<LANES>{});
      val<1> metasign = (meta[METAPIPE_V - 1] >= hard<0>{});
      metasign.fanout(hard<LANES>{});
      for (u64 r = 0; r < LANES; r++) {
        newly_alloc[r] = (match1[r] & coldctr) != hard<0>{};
      }
      newly_alloc.fanout(hard<2>{});
      arr<val<1>, LANES> altsel = [&](u64 r) {
        arr<val<1>, 3> inputs = {metasign, newly_alloc[r],
                                 match2[r] != hard<0>{}};
        return inputs.fo1().fold_and();
      };
      p2 = arr<val<1>, LANES>{[&](u64 r) {
             return select(altsel[r].fo1(), pred2_tage[r], pred1_tage[r]);
           }}.concat();
#ifdef TAGE_MONITOR
      for (u64 r = 0; r < LANES; r++) {
        bool has_alt = static_cast<u64>(match2[r] != hard<0>{});
        bool pri_ne_alt = static_cast<u64>(pred1_tage[r]) != static_cast<u64>(pred2_tage[r]);
        mon_altsel[r] = static_cast<u64>(altsel[r]);
        bool is_newly_alloc = static_cast<u64>(newly_alloc[r]);
        mon_meta_active[r] = has_alt && pri_ne_alt && is_newly_alloc;
      }
#endif
    } else {
      p2 = pred1_tage.concat();
#ifdef TAGE_MONITOR
      for (u64 r = 0; r < LANES; r++) { mon_altsel[r] = false; mon_meta_active[r] = false; }
#endif
    }

    p2.fanout(hard<LANES>{});
    // Select prediction by static PC-derived offset (not dynamic rank)
    val<LOG_LANES> pc_offset = val<LOG_LANES>{inst_pc >> 2};
    val<LANES> offset_mask = pc_offset.fo1().decode().concat();
    val<1> taken = (offset_mask & p2) != hard<0>{};
    taken.fanout(hard<2>{});
    reuse_prediction(~line_end());
    return taken;
  }

  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) override {
#ifdef TAGE_MONITOR
    mon.record_reuse_predict2();
#endif
    // Select prediction by static PC-derived offset (not dynamic rank)
    val<LOG_LANES> pc_offset = val<LOG_LANES>{inst_pc >> 2};
    val<LANES> offset_mask = pc_offset.fo1().decode().concat();
    val<1> taken = (offset_mask & p2) != hard<0>{};
    taken.fanout(hard<2>{});
    block_size++;
    reuse_prediction(~line_end());
    return taken;
  }

  void update_condbr([[maybe_unused]] val<64> branch_pc, val<1> taken,
                     [[maybe_unused]] val<64> next_pc) override {
    assert(num_branch < N);
#ifdef TAGE_MONITOR
    mon.record_branch_pc(static_cast<u64>(branch_pc));
#endif
    branch_dir[num_branch] = taken.fo1();
    // Save static PC-derived offset for this branch (like reference tage.hpp)
    branch_offset[num_branch] = val<LOG_LANES>{branch_pc >> 2};
    num_branch++;
    reuse_prediction(~(line_end() | last_pred()));
  }

  void update_cycle(instruction_info &block_end_info) override {
    val<1> &mispredict = block_end_info.is_mispredict;
    val<64> &next_pc = block_end_info.next_pc;

#ifdef TD_VERBOSE
    std::cerr << "UC: ENTER (num_branch=" << num_branch
              << " misp=" << static_cast<u64>(block_end_info.is_mispredict)
              << ")\n";
#endif
    if (num_branch == 0) {
      // No conditional branches — just update history
#ifdef TD_VERBOSE
      std::cerr << "UC: EXIT (no branches)\n";
#endif
      execute_if(true_block, [&]() {
        next_pc.fanout(hard<2>{});
        if constexpr (USE_DIR_HIST_V)
          gfolds.update(concat(val<1>{0}, val<PATHBITS>{next_pc >> 2}));
        else
          gfolds.update(val<PATHBITS>{next_pc >> 2});
        if constexpr (USE_PATH_HIST_V) {
          path_hist =
              (path_hist << PATHBITS) ^ val<PATH_HIST_WIDTH_V>{next_pc >> 2};
        }
        true_block = 1;
      });
      return;
    }

    // ================================================================
    // UPDATE CYCLE — Cycle 1: combinational logic + RAM reads
    // All code above need_extra_cycle() executes in the same cycle as
    // predict1/predict2. No RAM writes allowed here — only reads and
    // pure combinational logic on values already in regs.
    // ================================================================

    // ---- Fanouts for values read in predict2 (stored in regs) ----
    mispredict.fanout(hard<NUM_TABLES + 2>{});
    val<1> correct_pred = ~mispredict;
    if constexpr (ACC_PRESS_W > 0)
      correct_pred.fanout(hard<NUM_TABLES + 3>{});
    else
      correct_pred.fanout(hard<NUM_TABLES + 2>{});
    index1.fanout(hard<LANES * 3>{});
    gindex.fanout(hard<4>{});
    htag.fanout(hard<3>{});
    readt.fanout(hard<4>{});   // tag values read from TAGE tables in P2
    readc.fanout(hard<2>{});   // counter values read in P2
    if constexpr (MAX_HYST_WIDTH > 0)
      readh.fanout(hard<3>{}); // hysteresis values read in P2
    match1.fanout(hard<3>{});  // provider match vectors computed in P2
    match2.fanout(hard<2>{});  // alt match vectors computed in P2
    pred1_tage.fanout(hard<2>{});
    pred2_tage.fanout(hard<2 + NUM_TABLES>{});
    branch_dir.fanout(hard<USE_DIR_HIST_V ? 3 : 2>{});
    gfolds.fanout(hard<2>{});
    readu.fanout(hard<2>{});   // u-bit values read in P2
    X.fanout(hard<LANES + 1>{});
    if constexpr (USE_META_V)
      meta.fanout(hard<2>{});

    // ---- TAGE combinational logic ----
    // All pure logic — computes allocation decisions, counter update
    // conditions, etc. from the reg values read in predict2.
    // No RAM access here.

    // last_offset: PC-derived offset of the last conditional branch
    branch_offset.fanout(hard<LANES + NUM_TABLES + 2>{});
    val<LOG_LANES> last_offset = val<LOG_LANES>{branch_offset[num_branch - 1]};
    last_offset.fanout(hard<4 * NUM_TABLES + 2>{});

    // Per-offset: was this offset used by a branch in this block?
    // Computed from hardware: check if any branch_offset[j] matches offset o.
    arr<val<1>, LANES> is_branch = [&](u64 o) -> val<1> {
      arr<val<1>, LANES> mo = [&](u64 j) -> val<1> {
        if (j >= num_branch) return val<1>{0};
        return val<LOG_LANES>{branch_offset[j]} == val<LOG_LANES>{o};
      };
      return mo.fo1().fold_or();
    };
    if constexpr (AllocCfg::ALLOC_TRIGGER == AllocTrigger::DISAGREE)
      is_branch.fanout(hard<5>{});
    else
      is_branch.fanout(hard<4>{});

    // Per-offset: actual branch direction (0 for offsets without branches).
    // Maps offset → direction using branch_offset[] matching (like reference).
    u64 update_valid = (u64(1) << num_branch) - 1;
    val<LANES> actualdirs = branch_dir.concat();
    actualdirs.fanout(hard<LANES + NUM_TABLES>{});
    arr<val<1>, LANES> branch_taken = [&](u64 o) -> val<1> {
      arr<val<1>, LANES> match_off = [&](u64 j) -> val<1> {
        if (j >= num_branch) return val<1>{0};
        return val<LOG_LANES>{branch_offset[j]} == val<LOG_LANES>{o};
      };
      return (match_off.fo1().concat() & val<LANES>{update_valid} & actualdirs) != hard<0>{};
    };
    branch_taken.fanout(hard<3>{});

    // Restrict match vectors to lanes that actually had branches
    arr<val<NUM_TABLES + 1>, LANES> actual_match1 = [&](u64 r) {
      return select(is_branch[r], match1[r], val<NUM_TABLES + 1>{0});
    };
    actual_match1.fanout(hard<2>{});

    // primary_mask: OR of all per-lane provider matches — tells which
    // TAGE tables are the provider for at least one branch
    val<NUM_TABLES> primary_mask = actual_match1.fold_or();
    primary_mask.fanout(hard<2>{});
    arr<val<1>, NUM_TABLES> primary = primary_mask.make_array(val<1>{});
    primary.fanout(hard<3>{});

    // ---- Allocation decision (combinational) ----

    // Tag comparison for the last branch's offset specifically
    // (computed first so last_match1 is available for TAGE_MISS trigger)
    static_loop<NUM_TABLES>([&]<u64 I>() {
      using Table = std::tuple_element_t<I, Tables>;
      static constexpr u64 PER_HTAG = Table::tag_width - LOG_LANES;
      last_tagcmp_reg[I] = (val<LOG_LANES>{readt[I]} == last_offset)
                         & (val<PER_HTAG>{readt[I] >> LOG_LANES} == val<PER_HTAG>{htag[I]});
    });
    val<NUM_TABLES + 1> last_match1 =
        last_tagcmp_reg.fo1().append(1).concat().one_hot();
    if constexpr (AllocCfg::ALLOC_TRIGGER == AllocTrigger::TAGE_MISS)
      last_match1.fanout(hard<3>{});
    else
      last_match1.fanout(hard<2>{});

    // Allocation trigger: what condition enables allocation attempt
    val<1> alloc_trigger = [&]() -> val<1> {
      if constexpr (AllocCfg::ALLOC_TRIGGER == AllocTrigger::MISPREDICT) {
        return mispredict;
      } else if constexpr (AllocCfg::ALLOC_TRIGGER == AllocTrigger::BASE_WRONG) {
        return val<1>{pred[num_branch - 1]} != branch_dir[num_branch - 1];
      } else if constexpr (AllocCfg::ALLOC_TRIGGER == AllocTrigger::DISAGREE) {
        // P1 is rank-indexed, P2 is offset-indexed — compare per branch
        arr<val<1>, LANES> disag = [&](u64 j) -> val<1> {
          if (j >= num_branch) return val<1>{0};
          val<LOG_LANES> off = val<LOG_LANES>{branch_offset[j]};
          val<LANES> mask = off.fo1().decode().concat();
          val<1> p2_j = (mask & p2) != hard<0>{};
          return pred[j] != p2_j;
        };
        return disag.fo1().fold_or();
      } else if constexpr (AllocCfg::ALLOC_TRIGGER == AllocTrigger::TAGE_MISS) {
        return val<1>{last_match1 >> NUM_TABLES};
      } else {
        static_assert(AllocCfg::ALLOC_TRIGGER == AllocTrigger::ALWAYS);
        return val<1>{1};
      }
    }();
    alloc_trigger.fanout(hard<5>{});
    val<NUM_TABLES> mispmask = alloc_trigger.replicate(hard<NUM_TABLES>{}).concat();

    // Allocation action: probabilistic gating of allocation attempts
    val<NUM_TABLES> gated_mispmask = [&]() -> val<NUM_TABLES> {
      if constexpr (AllocCfg::ALLOC_ACTION == AllocAction::STANDARD) {
        return mispmask;
      } else if constexpr (AllocCfg::ALLOC_ACTION == AllocAction::FILTERED) {
        static_assert(ACC_PRESS_W > 0, "FILTERED requires ACCURACY_PRESSURE_BITS > 0");
        val<ACC_PRESS_W> rv = val<ACC_PRESS_W>{static_cast<u64>(std::rand())};
        val<1> allow = (val<ACC_PRESS_W>{accuracy_pressure} >= rv);
        return mispmask & allow.replicate(hard<NUM_TABLES>{}).concat();
      } else if constexpr (AllocCfg::ALLOC_ACTION == AllocAction::THROTTLED) {
        static_assert(ALLOC_PRESS_W > 0, "THROTTLED requires ALLOC_PRESSURE_BITS > 0");
        val<ALLOC_PRESS_W> rv = val<ALLOC_PRESS_W>{static_cast<u64>(std::rand())};
        val<1> allow = (val<ALLOC_PRESS_W>{alloc_pressure} >= rv);
        return mispmask & allow.replicate(hard<NUM_TABLES>{}).concat();
      }
    }();

    // postmask: tables above the provider (candidates for allocation)
    val<NUM_TABLES> postmask = [&]() -> val<NUM_TABLES> {
      if constexpr (AllocCfg::PROB_START > 0) {
        val<2> rstart = val<2>{static_cast<u64>(std::rand())};
        val<NUM_TABLES> base = gated_mispmask & val<NUM_TABLES>(last_match1 - 1);
        val<NUM_TABLES> skip1 = base & val<NUM_TABLES>(base - 1);
        val<NUM_TABLES> skip2 = skip1 & val<NUM_TABLES>(skip1 - 1);
        return select(rstart == hard<0>{}, skip2,
                      select(rstart == hard<1>{}, skip1, base));
      } else {
        return gated_mispmask & val<NUM_TABLES>(last_match1 - 1);
      }
    }();
    if constexpr (AllocCfg::UCTR_POLICY == UctrPolicy::PENALTY_NA)
      postmask.fanout(hard<3>{});
    else
      postmask.fanout(hard<2>{});

    // candallocmask: tables that are both above provider AND have u=0
    val<NUM_TABLES> candallocmask = [&]() -> val<NUM_TABLES> {
      if constexpr (AllocCfg::CONF_GATE) {
        arr<val<1>, NUM_TABLES> weak_entry = [&](u64 i) -> val<1> {
          if constexpr (MAX_CTR_WIDTH == 1) return val<1>{1};
          else {
            auto ctr = readc[i];
            return (ctr == hard<0>{}) |
                   (ctr == hard<(1u << (MAX_CTR_WIDTH - 1)) - 1>{}) |
                   (ctr == hard<(1u << (MAX_CTR_WIDTH - 1))>{});
          }
        };
        return postmask & notumask & weak_entry.fo1().concat();
      } else {
        return postmask & notumask;
      }
    }();
    candallocmask.fanout(hard<2>{});
    val<NUM_TABLES> collamask = candallocmask.reverse();
    collamask.fanout(hard<2>{});
    val<NUM_TABLES> collamask1 = collamask.one_hot();
    collamask1.fanout(hard<3>{});

    // allocate[i]: final per-table allocation decision (one-hot or two-hot)
    arr<val<1>, NUM_TABLES> allocate = [&]() -> arr<val<1>, NUM_TABLES> {
      if constexpr (AllocCfg::MAX_ALLOC >= 2) {
        val<NUM_TABLES> pick2 = [&]() -> val<NUM_TABLES> {
          val<NUM_TABLES> basic2 = (collamask ^ collamask1).one_hot();
          if constexpr (AllocCfg::NON_CONSECUTIVE) {
            val<NUM_TABLES> neighbors = (collamask1 << 1) | (collamask1 >> 1);
            val<NUM_TABLES> nc_mask = (collamask ^ collamask1) & ~neighbors;
            val<NUM_TABLES> nc_pick = nc_mask.reverse().one_hot();
            return select(nc_mask != hard<0>{}, nc_pick, basic2);
          } else { return basic2; }
        }();
        return (collamask1 | pick2).reverse().make_array(val<1>{});
      } else {
        val<NUM_TABLES> collamask2 = (collamask ^ collamask1).one_hot();
        val<NUM_TABLES> collamask12 =
            select(val<2>{std::rand()} == hard<0>{}, collamask2.fo1(), collamask1);
        return collamask12.fo1().reverse().make_array(val<1>{});
      }
    }();
    allocate.fanout(hard<7>{});

    // ---- Per-table update conditions (combinational) ----

    // bdir[i]: branch direction for the offset stored in table i's tag.
    // On allocation, uses last_offset instead.  Follows reference pattern:
    // match branch_offset[j] against stored/allocated offset to find direction.
    arr<val<1>, NUM_TABLES> bdir = [&](u64 i) {
      val<LOG_LANES> stored_offset = val<LOG_LANES>{readt[i]};
      val<LOG_LANES> use_offset = select(allocate[i], last_offset, stored_offset.fo1());
      use_offset.fanout(hard<LANES>{});
      arr<val<1>, LANES> match_off = [&](u64 j) -> val<1> {
        if (j >= num_branch) return val<1>{0};
        return val<LOG_LANES>{branch_offset[j]} == use_offset;
      };
      return (match_off.fo1().concat() & val<LANES>{update_valid} & actualdirs) != hard<0>{};
    };
    bdir.fanout(hard<2>{});

    // badpred1[i]: did table i's counter predict wrong?
    arr<val<1>, NUM_TABLES> badpred1 = [&](u64 i) -> val<1> {
      if constexpr (MAX_CTR_WIDTH == 1) return readc[i] != bdir[i];
      else return val<1>{readc[i] >> hard<MAX_CTR_WIDTH - 1>{}} != bdir[i];
    };
    badpred1.fanout(hard<3>{});

    // altdiffer[i]: does provider disagree with alt prediction?
    arr<val<1>, NUM_TABLES> altdiffer = [&](u64 i) -> val<1> {
      auto pred_dir = [&]() -> val<1> {
        if constexpr (MAX_CTR_WIDTH == 1) return readc[i];
        else return readc[i] >> hard<MAX_CTR_WIDTH - 1>{};
      }();
      val<LOG_LANES> stored_offset = val<LOG_LANES>{readt[i]};
      return pred_dir != pred2_tage.select(stored_offset.fo1());
    };

    // goodpred[i]: was this table's prediction correct for its stored offset?
    arr<val<1>, NUM_TABLES> goodpred = [&](u64 i) {
      val<LOG_LANES> stored_offset = val<LOG_LANES>{readt[i]};
      return (stored_offset.fo1() != last_offset) | correct_pred;
    };
    goodpred.fanout(hard<2>{});

    // g_weak[i]: provider with wrong prediction AND weak hysteresis — flip counter
    arr<val<1>, NUM_TABLES> g_weak = [&](u64 i) -> val<1> {
      if constexpr (MAX_HYST_WIDTH > 0)
        return primary[i] & badpred1[i] & (readh[i] == hard<0>{});
      else
        return primary[i] & badpred1[i];
    };
    g_weak.fanout(hard<2>{});

    // P1 vs P2 disagreement — P1 is rank-indexed, P2 is offset-indexed.
    // Build per-branch disagreement: for each branch j at offset branch_offset[j],
    // compare pred[j] (P1, rank-indexed) with p2 bit at that offset.
    arr<val<1>, LANES> per_branch_disagree = [&](u64 o) -> val<1> {
      // Find if any branch maps to offset o, get its P1 prediction
      arr<val<1>, LANES> match_j = [&](u64 j) -> val<1> {
        if (j >= num_branch) return val<1>{0};
        return val<LOG_LANES>{branch_offset[j]} == val<LOG_LANES>{o};
      };
      match_j.fanout(hard<2>{});
      val<1> any_match = match_j.fold_or();
      val<1> p1_for_o = (match_j.fo1().concat() & pred.concat()) != hard<0>{};
      val<1> p2_for_o = (p2 >> o) & hard<1>{};
      return any_match & (p1_for_o != p2_for_o);
    };
    val<LANES> disagree_mask = per_branch_disagree.concat();
    disagree_mask.fanout(hard<2>{});

    // U-bit clear helpers (combinational)
    arr<val<1>, NUM_TABLES> update_u = [&](u64 i) {
      return primary[i] & altdiffer[i].fo1();
    };
    val<1> noalloc = (candallocmask == hard<0>{});
    // fanout: uclearmask + uctr_incr (NOALLOC/PENALTY_NA) + alloc_pressure
    if constexpr (AllocCfg::UCTR_POLICY == UctrPolicy::NOALLOC ||
                  AllocCfg::UCTR_POLICY == UctrPolicy::PENALTY_NA ||
                  ALLOC_PRESS_W > 0)
      noalloc.fanout(hard<3>{});
    else
      noalloc.fanout(hard<2>{});
    val<NUM_TABLES> uclearmask =
        postmask & noalloc.replicate(hard<NUM_TABLES>{}).concat();
    arr<val<1>, NUM_TABLES> uclear = uclearmask.fo1().make_array(val<1>{});
    uclear.fanout(hard<2>{});

    // ---- P1 gshare + TAGE reads / need_extra / writes ----
    // Same structure as gshareN reference: reads in cycle 1, writes in cycle 2.

    // Combinational: compute which lanes were accessed by branches in this
    // block. access[i] = OR of all X.rotate_left(branch_rank) masks — 1 if lane
    // i was used.
    arr<val<1>, LANES> access =
        arr<val<LANES>, LANES>{[&](u64 i) -> val<LANES> {
          return X.rotate_left(i) & val<LANES>{-(i < num_branch)};
        }}.fold_or()
            .make_array(val<1>{});

    // Combinational: identify the lane that holds the mispredicted branch.
    // misp_bank is a one-hot mask ANDed with mispredict — all zero on correct
    // prediction.
    val<LANES> misp_bank = X.rotate_left(num_branch - 1) &
                           mispredict.replicate(hard<LANES>{}).concat();
    arr<val<1>, LANES> mispredicted = misp_bank.fo1().make_array(val<1>{});
    mispredicted.fanout(hard<2>{});

    // Cycle 1 RAM read: read hysteresis for the mispredicted lane.
    // Gated by mispredicted[i] so only one lane's RAM is actually read.
    // On correct prediction, mispredicted is all-zero so no RAM is accessed.
    arr<val<1>, LANES> weak = [&](u64 i) -> val<1> {
      if (i >= N) return val<1>{0};
      return execute_if(mispredicted[i], [&]() {
#ifdef TD_VERBOSE
        if (static_cast<u64>(mispredict) && static_cast<u64>(mispredicted[i]))
          std::cerr << "UC: gshare p1_hyst[" << i << "] READ\n";
#endif
        return p1_hyst[i].read(index1);
      });
    };

    // ================================================================
    // Cycle boundary: grant an extra cycle for RAM writes.
    // tag_ram/pred_ram use plain ram<> — need extra_cycle for writes.
    // hyst_ram/u_ram use td_rwram — pass extra_cycle as noconflict,
    // so they can buffer writes when extra_cycle=0.
    // Everything above = cycle 1. Everything below = cycle 2.
    // ================================================================
    val<1> extra_cycle = [&]() -> val<1> {
      if constexpr (AllocCfg::MISPREDICT_ONLY_WRITE) {
        if constexpr (AllocCfg::ALLOC_TRIGGER != AllocTrigger::MISPREDICT)
          return mispredict | alloc_trigger;
        else
          return mispredict;
      } else {
        val<1> some_badpred1 = (primary_mask & badpred1.concat()) != hard<0>{};
        val<1> base = some_badpred1.fo1() | mispredict;
        if constexpr (AllocCfg::DISAGREE_EXTRA_CYCLE &&
                      AllocCfg::ALLOC_TRIGGER != AllocTrigger::MISPREDICT) {
          return base | (disagree_mask != hard<0>{}) | alloc_trigger;
        } else if constexpr (AllocCfg::DISAGREE_EXTRA_CYCLE) {
          return base | (disagree_mask != hard<0>{});
        } else if constexpr (AllocCfg::ALLOC_TRIGGER != AllocTrigger::MISPREDICT) {
          return base | alloc_trigger;
        } else {
          return base;
        }
      }
    }();
#ifdef TD_VERBOSE
    std::cerr << "UC: extra_cycle=" << static_cast<u64>(extra_cycle) << "\n";
#endif
    extra_cycle.fanout(hard<NUM_TABLES * 2 + 1>{});
    need_extra_cycle(extra_cycle);

#ifdef TAGE_MONITOR
    // Block stats
    mon.begin_update_cycle(static_cast<u64>(extra_cycle));
    mon.record_block(static_cast<u64>(block_entry), block_size, num_branch,
                     static_cast<u64>(extra_cycle));
    // Per-branch outcome + prediction tracking (p2 indexed by offset now)
    for (u64 r = 0; r < num_branch; r++) {
      bool actual = static_cast<u64>(branch_dir[r]);
      bool p1_pr = static_cast<u64>(pred[r]);
      u64 off_c = static_cast<u64>(branch_offset[r]) & (LANES - 1);
      bool p2_pr = (static_cast<u64>(p2) >> off_c) & 1;
      bool misp = (r == num_branch - 1) ? static_cast<u64>(mispredict) : false;
      mon.record_prediction(r, static_cast<u64>(match1[off_c]),
                            static_cast<u64>(match2[off_c]),
                            mon_meta_active[off_c], mon_altsel[off_c],
                            p1_pr, p2_pr);
      mon.record_outcome(r, actual, misp);
    }
    // Tag match tracking
    for (u64 i = 0; i < NUM_TABLES; i++) {
      bool matched = static_cast<u64>(htagcmp_reg[i]);
      mon.record_tag_lookup(i, matched);
    }
    // Allocation
    if (static_cast<u64>(alloc_trigger)) {
      u64 amask = 0;
      for (u64 i = 0; i < NUM_TABLES; i++) {
        if (static_cast<u64>(allocate[i])) {
          amask |= (u64(1) << i);
          mon.record_tage_write(i, static_cast<u64>(gindex[i]) % MAX_TABLE_SIZE);
        }
      }
      u64 pmask = static_cast<u64>(postmask);
      if (pmask == 0) {
        mon.record_alloc_blocked();
      }
      mon.record_allocation(amask != 0, amask);
      // Cascade: find provider from last_match1
      u64 lm1 = static_cast<u64>(last_match1);
      u64 prov_idx = NUM_TABLES; // default = P1
      for (u64 i = 0; i < NUM_TABLES; i++)
        if (lm1 & (u64(1) << i)) { prov_idx = i; break; }
      mon.record_alloc_cascade(prov_idx, amask);
    }
#endif

    // ================================================================
    // Cycle 2: RAM writes — P1 gshare + TAGE tables
    // All writes are gated by execute_if so they only fire when needed.
    // ================================================================

    // ---- P1 gshare writes ----
    // Flip prediction bit if hysteresis was weak on misprediction.
    // p1_pred was read in predict1 (cycle 1), safe to write in cycle 2.
    execute_if(mispredict, [&]() {
      arr<val<1>, LANES> stored = unordered_pred1.make_array(val<1>{});
      arr<val<1>, LANES> bundle = [&](u64 i) {
        return select(weak[i].fo1(), branch_dir[num_branch - 1],
                      stored[i].fo1());
      };
#ifdef TD_VERBOSE
      if (static_cast<u64>(mispredict))
        std::cerr << "UC: gshare p1_pred WRITE\n";
#endif
      p1_pred.write(index1, bundle.fo1().concat());
    });

    // Update P1 hysteresis for all accessed lanes.
    // p1_hyst[i] was read in cycle 1, safe to write in cycle 2.
    for (u64 i = 0; i < N; i++) {
      execute_if(access[i].fo1(), [&]() {
#ifdef TD_VERBOSE
        if (static_cast<u64>(mispredict) && static_cast<u64>(access[i]))
          std::cerr << "UC: gshare p1_hyst[" << i << "] WRITE\n";
#endif
        p1_hyst[i].write(index1, mispredicted[i].fo1());
#ifdef TAGE_MONITOR
        mon.record_p1_write(i, static_cast<u64>(index1));
#endif
      });
    }

    // ---- TAGE tag write (allocation) ----
    // Write new tag = concat(htag, last_offset) into allocated table entries.
    // tag_ram was read in predict2 (cycle 1), safe to write in cycle 2.
    // Gate: allocate[I] ⊆ mispredict ⊆ extra_cycle, so no third arg needed.
    static_loop<NUM_TABLES>([&]<u64 I>() {
      execute_if(allocate[I], [&]() {
#ifdef TD_VERBOSE
        if (static_cast<u64>(allocate[I]))
          std::cerr << "UC: tage tag_ram[" << I << "] WRITE (alloc)\n";
#endif
        std::get<I>(tables).tag_ram.write(tidx<I>(gindex[I]),
                                          concat(htag[I], last_offset));
#ifdef TAGE_MONITOR
        mon.mark_write();
#endif
      });
    });

    // ---- TAGE counter update ----
    // pred_ram was read in predict2 (cycle 1), safe to write in cycle 2.
    if constexpr (AllocCfg::MISPREDICT_ONLY_WRITE) {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        auto pred_gate = [&]() {
          if constexpr (AllocCfg::ALLOC_TRIGGER != AllocTrigger::MISPREDICT)
            return allocate[I] | (mispredict & g_weak[I].fo1());
          else
            return mispredict & (g_weak[I].fo1() | allocate[I]);
        }();
        execute_if(pred_gate, [&]() {
#ifdef TD_VERBOSE
          std::cerr << "UC: tage pred_ram[" << I << "] WRITE (misp_only)\n";
#endif
          if constexpr (MAX_CTR_WIDTH == 1) {
            std::get<I>(tables).pred_ram.write(tidx<I>(gindex[I]), bdir[I]);
          } else {
            auto init_ctr = select(bdir[I],
                val<MAX_CTR_WIDTH>{(1u << MAX_CTR_WIDTH) - 1},
                val<MAX_CTR_WIDTH>{0});
            std::get<I>(tables).pred_ram.write(tidx<I>(gindex[I]),
                select(allocate[I], init_ctr, val<MAX_CTR_WIDTH>{bdir[I]}));
          }
        });
      });
    } else {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> old_dir = [&]() -> val<1> {
          if constexpr (MAX_CTR_WIDTH == 1) return readc[I];
          else return readc[I] >> hard<MAX_CTR_WIDTH - 1>{};
        }();
        val<1> pred_changed = (old_dir != bdir[I]) | allocate[I];
        val<1> pred_eligible = g_weak[I].fo1() | allocate[I];
#ifdef TAGE_MONITOR
        if (static_cast<u64>(pred_eligible))
          mon.record_silent_pred(I, !static_cast<u64>(pred_changed));
#endif
        execute_if(pred_eligible & pred_changed, [&]() {
#ifdef TD_VERBOSE
          std::cerr << "UC: tage pred_ram[" << I << "] WRITE (standard)\n";
#endif
          std::get<I>(tables).pred_ram.write(tidx<I>(gindex[I]), bdir[I]);
#ifdef TAGE_MONITOR
          mon.mark_write();
#endif
        });
      });
    }

    // ---- TAGE hysteresis update ----
    // hyst_ram was read in predict2 (cycle 1), safe to write in cycle 2.
    if constexpr (MAX_HYST_WIDTH > 0 && !AllocCfg::MISPREDICT_ONLY_WRITE) {
      static constexpr u64 HW = std::max(u64(1), MAX_HYST_WIDTH);
      static constexpr u64 HMAX = (u64(1) << HW) - 1;
      if constexpr (AllocCfg::PARTIAL_UPDATE) altdiffer.fanout(hard<2>{});
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> would_change = allocate[I] |
                              (badpred1[I] & (readh[I] != hard<0>{})) |
                              (~badpred1[I] & (readh[I] != hard<HMAX>{}));
        val<1> should_update = [&]() -> val<1> {
          if constexpr (AllocCfg::PARTIAL_UPDATE)
            return allocate[I] | badpred1[I] | altdiffer[I].fo1();
          else
            return primary[I] | allocate[I];
        }();
#ifdef TAGE_MONITOR
        if (static_cast<u64>(should_update))
          mon.record_silent_hyst(I, !static_cast<u64>(would_change));
#endif
        execute_if(should_update & would_change, [&]() {
          auto newhyst = select(allocate[I], val<HW>{0},
                                td::update_ctr(readh[I], ~badpred1[I]));
#ifdef TD_VERBOSE
          std::get<I>(tables).hyst_ram.debug_write_info("hyst_ram", I, hidx<I>(gindex[I]), extra_cycle);
#endif
          std::get<I>(tables).hyst_ram.write(hidx<I>(gindex[I]), newhyst.fo1(), extra_cycle);
#ifdef TAGE_MONITOR
          { auto &r = std::get<I>(tables).hyst_ram;
            u64 bid = (static_cast<u64>(gindex[I]) >> RWRAM_BANK_SHIFT_V) & (RWRAM_BANKS_V - 1);
            mon.record_rwram_write(I, bid, static_cast<u64>(extra_cycle),
                                   static_cast<u64>(r.read_bank), static_cast<u64>(r.write_bank)); }
#endif
        });
      });
    } else if constexpr (MAX_HYST_WIDTH > 0 && AllocCfg::MISPREDICT_ONLY_WRITE) {
      static constexpr u64 HW = std::max(u64(1), MAX_HYST_WIDTH);
      static constexpr u64 HMAX = (u64(1) << HW) - 1;
      static_loop<NUM_TABLES>([&]<u64 I>() {
        auto hyst_gate = [&]() {
          if constexpr (AllocCfg::ALLOC_TRIGGER != AllocTrigger::MISPREDICT)
            return allocate[I];
          else
            return mispredict & allocate[I];
        }();
        execute_if(hyst_gate, [&]() {
#ifdef TD_VERBOSE
          std::get<I>(tables).hyst_ram.debug_write_info("hyst_ram", I, hidx<I>(gindex[I]), extra_cycle);
#endif
          std::get<I>(tables).hyst_ram.write(hidx<I>(gindex[I]), val<HW>{HMAX}, extra_cycle);
#ifdef TAGE_MONITOR
          { auto &r = std::get<I>(tables).hyst_ram;
            u64 bid = (static_cast<u64>(gindex[I]) >> RWRAM_BANK_SHIFT_V) & (RWRAM_BANKS_V - 1);
            mon.record_rwram_write(I, bid, static_cast<u64>(extra_cycle),
                                   static_cast<u64>(r.read_bank), static_cast<u64>(r.write_bank)); }
#endif
        });
      });
    }

    // ---- TAGE u-bit update ----
    // u_ram uses td_rwram: pass extra_cycle as noconflict for buffered writes.
    if constexpr (AllocCfg::MISPREDICT_ONLY_WRITE) {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        auto u_gate = [&]() {
          if constexpr (AllocCfg::ALLOC_TRIGGER != AllocTrigger::MISPREDICT)
            return allocate[I] | (mispredict & uclear[I]);
          else
            return mispredict & (allocate[I] | uclear[I]);
        }();
        execute_if(u_gate, [&]() {
          val<1> newu = select(allocate[I], val<1>{1}, val<1>{0});
#ifdef TD_VERBOSE
          std::get<I>(tables).u_ram.debug_write_info("u_ram", I, tidx<I>(gindex[I]), extra_cycle);
#endif
          std::get<I>(tables).u_ram.write(tidx<I>(gindex[I]), newu.fo1(), extra_cycle);
#ifdef TAGE_MONITOR
          mon.record_u_write(I, static_cast<u64>(newu));
          { auto &r = std::get<I>(tables).u_ram;
            u64 bid = (static_cast<u64>(gindex[I]) >> RWRAM_BANK_SHIFT_V) & (RWRAM_BANKS_V - 1);
            mon.record_rwram_write(I, bid, static_cast<u64>(extra_cycle),
                                   static_cast<u64>(r.read_bank), static_cast<u64>(r.write_bank)); }
#endif
        });
      });
    } else if constexpr (USE_PROB_DECAY) {
      val<DECAY_CTR_V> lfsr = val<DECAY_CTR_V>{static_cast<u64>(std::rand())};
      val<1> decay_fire = (lfsr > val<DECAY_CTR_V>{decay_threshold});
      decay_fire.fanout(hard<NUM_TABLES>{});
#ifdef TAGE_MONITOR
      if (static_cast<u64>(decay_fire)) mon.record_decay_fire();
#endif
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> newu = goodpred[I].fo1() & ~allocate[I] & ~uclear[I] & ~decay_fire;
        val<1> u_changed = (val<1>{readu[I]} != newu);
        val<1> u_eligible = update_u[I].fo1() | allocate[I] | uclear[I] | decay_fire;
#ifdef TAGE_MONITOR
        if (static_cast<u64>(u_eligible))
          mon.record_silent_u(I, !static_cast<u64>(u_changed));
#endif
        execute_if(u_eligible & u_changed, [&]() {
#ifdef TD_VERBOSE
          std::get<I>(tables).u_ram.debug_write_info("u_ram", I, tidx<I>(gindex[I]), extra_cycle);
#endif
          std::get<I>(tables).u_ram.write(tidx<I>(gindex[I]), newu.fo1(), extra_cycle);
#ifdef TAGE_MONITOR
          mon.record_u_write(I, static_cast<u64>(newu));
          { auto &r = std::get<I>(tables).u_ram;
            u64 bid = (static_cast<u64>(gindex[I]) >> RWRAM_BANK_SHIFT_V) & (RWRAM_BANKS_V - 1);
            mon.record_rwram_write(I, bid, static_cast<u64>(extra_cycle),
                                   static_cast<u64>(r.read_bank), static_cast<u64>(r.write_bank)); }
#endif
        });
      });
    } else {
      static_loop<NUM_TABLES>([&]<u64 I>() {
        val<1> newu = goodpred[I].fo1() & ~allocate[I] & ~uclear[I];
        val<1> u_changed = (val<1>{readu[I]} != newu);
        val<1> u_eligible_nd = update_u[I].fo1() | allocate[I] | uclear[I];
#ifdef TAGE_MONITOR
        if (static_cast<u64>(u_eligible_nd))
          mon.record_silent_u(I, !static_cast<u64>(u_changed));
#endif
        execute_if(u_eligible_nd & u_changed, [&]() {
#ifdef TD_VERBOSE
          std::get<I>(tables).u_ram.debug_write_info("u_ram", I, tidx<I>(gindex[I]), extra_cycle);
#endif
          std::get<I>(tables).u_ram.write(tidx<I>(gindex[I]), newu.fo1(), extra_cycle);
#ifdef TAGE_MONITOR
          mon.record_u_write(I, static_cast<u64>(newu));
          { auto &r = std::get<I>(tables).u_ram;
            u64 bid = (static_cast<u64>(gindex[I]) >> RWRAM_BANK_SHIFT_V) & (RWRAM_BANKS_V - 1);
            mon.record_rwram_write(I, bid, static_cast<u64>(extra_cycle),
                                   static_cast<u64>(r.read_bank), static_cast<u64>(r.write_bank)); }
#endif
        });
      });
    }

    // ---- Meta counter update (regs only, no RAM) ----
    if constexpr (USE_META_V) {
      arr<val<1>, LANES> altdiff = [&](u64 r) {
        return (match2[r] != hard<0>{}) & (pred2_tage[r] != pred1_tage[r]);
      };
      arr<val<2, i64>, LANES> meta_incr = [&](u64 r) -> val<2, i64> {
        val<1> update_meta = is_branch[r] & altdiff[r].fo1() & newly_alloc[r];
        val<1> bad_pred2 = (pred2_tage[r] != branch_taken[r]);
        return select(update_meta.fo1(), concat(bad_pred2.fo1(), val<1>{1}),
                      val<2>{0});
      };
      for (u64 i = METAPIPE_V - 1; i != 0; i--) meta[i] = meta[i - 1];
      auto newmeta = meta[0] + meta_incr.fo1().fold_add();
      newmeta.fanout(hard<3>{});
      using meta_t = valt<decltype(meta[0])>;
      meta[0] = select(newmeta > meta_t::maxval, meta_t{meta_t::maxval},
                       select(newmeta < meta_t::minval, meta_t{meta_t::minval},
                              meta_t{newmeta}));
    }

    // ---- Epoch / decay threshold (regs only, no RAM except u_ram.reset) ----
    uctr.fanout(hard<3>{});

    // UctrPolicy dispatch: compute uctr increment signal
    val<1> uctr_incr = [&]() -> val<1> {
      if constexpr (AllocCfg::UCTR_POLICY == UctrPolicy::FARALLOC) {
        val<NUM_TABLES> allocmask1 = collamask1.reverse();
        allocmask1.fanout(hard<2>{});
        return (((last_match1 >> 3) | allocmask1).one_hot() ^ allocmask1) == hard<0>{};
      } else if constexpr (AllocCfg::UCTR_POLICY == UctrPolicy::NOALLOC) {
        return noalloc;
      } else if constexpr (AllocCfg::UCTR_POLICY == UctrPolicy::ALWAYS_INC) {
        return val<1>{1};
      } else {
        return val<1>{0};  // PENALTY_NA: handled separately below
      }
    }();

    val<1> uctrsat = (uctr == hard<decltype(uctr)::maxval>{});
    uctrsat.fanout(hard<2>{});

    if constexpr (AllocCfg::UCTR_POLICY == UctrPolicy::PENALTY_NA) {
      // Seznec-style: update on every alloc attempt, not just mispredictions
      val<1> alloc_attempted = alloc_trigger & (postmask != hard<0>{});
      uctr = select(uctrsat, val<decltype(uctr)::size>{0},
                    select(alloc_attempted,
                           td::update_ctr(uctr, noalloc),
                           val<decltype(uctr)::size>{uctr}));
    } else {
      // FARALLOC, NOALLOC, ALWAYS_INC: only update on misprediction
      uctr = select(correct_pred, uctr,
                    select(uctrsat, val<decltype(uctr)::size>{0},
                           td::update_ctr(uctr, uctr_incr.fo1())));
    }

#ifdef TAGE_MONITOR
    mon.record_uctr(static_cast<u64>(uctr));
#endif

    // ---- Pressure register updates (regs only) ----
    if constexpr (ALLOC_PRESS_W > 0) {
      alloc_pressure = select(alloc_trigger,
          td::update_ctr(alloc_pressure, noalloc),
          val<ALLOC_PRESS_W>{alloc_pressure});
    }
    if constexpr (ACC_PRESS_W > 0) {
      accuracy_pressure = td::update_ctr(accuracy_pressure, ~correct_pred);
    }

    if constexpr (USE_PROB_DECAY) {
      val<1> threshold_tick = [&]() -> val<1> {
        if constexpr (DECAY_GRAN_V == 0) return ~correct_pred;
        else return (uctr & hard<(u64(1) << DECAY_GRAN_V) - 1>{}) == hard<0>{};
      }();
      val<1> misp = ~correct_pred;
      decay_threshold = select(threshold_tick,
          DecayPolicy_V::template apply<DECAY_CTR_V>(
              decay_threshold, correct_pred, uctrsat, misp),
          val<DECAY_CTR_V>{decay_threshold});
    } else {
      // Periodic u-bit reset — u_ram.reset() is a bulk clear, safe in cycle 2.
      // For PENALTY_NA, uctrsat can fire without mispredict, but alloc_trigger
      // is in extra_cycle so it's safe when allocation was attempted.
      execute_if(uctrsat & extra_cycle, [&]() {
#ifdef TAGE_MONITOR
        if (static_cast<u64>(uctrsat)) mon.record_epoch_reset();
#endif
        static_loop<NUM_TABLES>([&]<u64 I>() {
#ifdef TD_VERBOSE
          if (static_cast<u64>(uctrsat))
            std::cerr << "UC: tage u_ram[" << I << "] RESET\n";
#endif
          std::get<I>(tables).u_ram.reset();
        });
      });
    }

    // ---- History update ----
    true_block = arr<val<1>, 4>{~mispredict, branch_dir[num_branch - 1],
                                last_pred(), line_end()}
                     .fold_or();
    true_block.fanout(hard<MAXHIST + NUM_TABLES * 2 + 2>{});
    execute_if(true_block, [&]() {
      next_pc.fanout(hard<2>{});
      if constexpr (USE_DIR_HIST_V)
        gfolds.update(concat(branch_dir[num_branch - 1], val<PATHBITS>{next_pc >> 2}));
      else
        gfolds.update(val<PATHBITS>{next_pc >> 2});
      if constexpr (USE_PATH_HIST_V) {
        path_hist =
            (path_hist << PATHBITS) ^ val<PATH_HIST_WIDTH_V>{next_pc >> 2};
      }
    });

#ifdef TAGE_MONITOR
    mon.end_update_cycle();
#endif
#ifdef TD_VERBOSE
    std::cerr << "UC: EXIT (full update)\n";
#endif
    num_branch = 0;
  }
};

// ============================================================================
// User-facing Alias
// ============================================================================

template <typename TableCfg = td::TDTableConfig<>,
          typename AllocCfg = td::TDDefaultAllocConfig,
          u64 TD_LINEINST = 1024,
          u64 TD_N = 7,
          u64 TD_DECAY_CTR = 0, u64 TD_DECAY_GRAN = 2,
          typename TD_DECAY_POLICY = td::TDDecayMild,
          u64 TD_P1_TABLE_SIZE = 256,
          bool TD_USE_META = true, u64 TD_METABITS = 4,
          u64 TD_METAPIPE = 2,
          bool TD_USE_PATH_HIST = false,
          u64 TD_PATH_HIST_WIDTH = 27, u64 TD_PATH_BITS = 6,
          template <u64> class TD_FOLD_FN = td::XORFold,
          u64 TD_RWRAM_BANKS = 4, u64 TD_RWRAM_BANK_SHIFT = 0,
          u64 TD_EPOCH_CTR_BITS = 8,
          bool TD_SHARED_HYS = false,
          bool TD_USE_DIR_HIST = false>

using TageDirectBim =
    TageDirectBimImpl<TableCfg, AllocCfg, TD_LINEINST, TD_N, TD_DECAY_CTR,
                      TD_DECAY_GRAN, TD_DECAY_POLICY,
                      TD_P1_TABLE_SIZE, TD_USE_META, TD_METABITS,
                      TD_METAPIPE, TD_USE_PATH_HIST, TD_PATH_HIST_WIDTH,
                      TD_PATH_BITS, TD_FOLD_FN, TD_RWRAM_BANKS,
                      TD_RWRAM_BANK_SHIFT, TD_EPOCH_CTR_BITS, TD_SHARED_HYS,
                      TD_USE_DIR_HIST>;
