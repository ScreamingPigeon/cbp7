#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"
#include "custom_common.hpp"

using namespace hcm;

// ============================================================================
// TageAheadHC_IR — Hand-coded TageAhead with Independent Resolution (BPE1)
//
// Based on TageAheadHC (S3_MW4_MC2048) with per-branch independent
// resolution chains. Each TAGE entry stores 1-bit prediction for one
// branch. The 11-bit tag is split into 8-bit htag + 3-bit group_id.
// 7 independent resolve_chains run in parallel.
//
// Table sizes (T0 bumped from 512→1024):
//   T0-4:  1024 entries (IDX=10, HYST=512)
//   T5-13: 2048 entries (IDX=11, HYST=1024)
// ============================================================================

template <u64 LINEINST_V = 256>
struct TageAheadHC_IR_impl : predictor {

  // ======== Constants (hardcoded from S3_MW4_MC2048 + BPE1) ========
  static constexpr u64 NT = 15;           // number of tables
  static constexpr u64 N = 7;             // max conditional branches per block
  static constexpr u64 PATHBITS = 6;
  static constexpr u64 SEC_TAG_BITS = 5;
  static constexpr u64 CTR_WIDTH = 1;
  static constexpr u64 HYST_WIDTH = 2;
  static constexpr u64 U_WIDTH = 2;
  static constexpr u64 TAG_WIDTH = 11;    // UniformTag<11>
  static constexpr u64 MAXHIST = 200;
  static constexpr u64 LINEINST = LINEINST_V;
  static constexpr u64 LOGLINEINST = []() constexpr { u64 x = LINEINST_V, r = 0; while (x > 1) { x >>= 1; r++; } return r; }();
  static constexpr u64 FB_CAPACITY = 8192;
  static constexpr u64 FB_IDX_BITS = 13;  // clog2(8192)
  static constexpr u64 FB_BH_CAPACITY = FB_CAPACITY / 2;
  static constexpr u64 FB_BH_IDX_BITS = 12; // clog2(4096)
  static constexpr u64 META_WIDTH = 4;
  static constexpr u64 META_CAPACITY = 2048;
  static constexpr u64 META_IDX_BITS = 11; // clog2(2048)
  static constexpr u64 META_PIPE = 2;
  static constexpr u64 META_BANKS = 2;
  static constexpr u64 MATCH_BITS = NT + 1; // = 16
  static constexpr u64 ALLOC_PC_BITS = TAG_WIDTH + 2; // = 13
  static constexpr u64 ACC_WIDTH = 10;
  static constexpr u64 ALLOC_WIDTH = 10;
  static constexpr u64 MAX_IDX_BITS = 11; // clog2(2048)

  // ======== IR (Independent Resolution) constants ========
  static constexpr u64 NUM_GROUPS = N;    // = 7, one group per branch
  static constexpr u64 GROUP_BITS = 3;    // clog2(7) rounded up
  // Group-to-meta-bank mapping: i * META_BANKS / N
  // Groups 0-3 → bank 0, groups 4-6 → bank 1
  static constexpr std::array<u64, NUM_GROUPS> GROUP_TO_META_BANK = {
    0, 0, 0, 0, 1, 1, 1
  };
  static constexpr u64 PRED_BITS = 1;     // 1-bit prediction per entry
  static constexpr u64 FB_PRED_BITS = N * CTR_WIDTH; // = 7, fallback stays N-wide
  static constexpr u64 HTAG_WIDTH = TAG_WIDTH - GROUP_BITS; // = 8 (max, for pipeline regs)

  // GradedTag<11,7>: T0 (short hist) → 7, T13 (long hist) → 11
  // Formula: 7 + 4*i/13
  static constexpr std::array<u64, NT> PER_TABLE_TAG = {
    7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11
  };
  static constexpr std::array<u64, NT> PER_TABLE_HTAG = {
    4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8
  };

  // Per-table sizes: T0 bumped to 1024 (was GradedSize<512, 2048>)
  static constexpr std::array<u64, NT> TABLE_SIZE = {
    1024, 1024, 1024, 1024, 1024,
    2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048
  };
  // Per-table history lengths: geometric_hist<14>(8, 200)
  static constexpr std::array<u64, NT> HIST_LEN =
      ta::geometric_hist<NT>(8, 200);
  // Per-table IDX bits
  static constexpr std::array<u64, NT> IDX_BITS = {
    10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11
  };
  // Per-table hyst sizes (shared hyst: TABLE_SIZE/2)
  static constexpr std::array<u64, NT> HYST_SIZE = {
    512, 512, 512, 512, 512,
    1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024
  };

  // Per-bit GH fanout (USE_GSHARE=false, so no fb_fold contribution)
  static constexpr auto GH_FANOUT =
      ta::gh_per_bit_fanout<MAXHIST, NT, HIST_LEN>();

  // ======== Global History ========
  ta_global_history<MAXHIST> gh;

  // ======== Fallback (bimodal, ahead-pipelined) ========
  reg<FB_PRED_BITS> prefetch_fb;
  reg<FB_IDX_BITS> prefetch_fb_idx;

  // Piped PC for allocation tag recomputation
  reg<ALLOC_PC_BITS> prefetch_pc;

  // ======== Block tracking ========
  reg<1> true_block = 1;
  reg<1> last_condbr_dir = 1;
  reg<LOGLINEINST> block_entry;

  // Simulation artifacts (free in hardware)
  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<1>, N> branch_dir;

  // Secondary tag (computed in predict1 from inst_pc)
  reg<SEC_TAG_BITS> curr_sec_tag;

  // ======== Pipeline regs [NT] ========
  reg<1> prefetch_tag_hit[NT];       // htag-only match (lower HTAG_WIDTH bits)
  reg<PRED_BITS> prefetch_pred[NT];  // 1-bit prediction per entry
  reg<SEC_TAG_BITS> prefetch_sec[NT];
  reg<MAX_IDX_BITS> prefetch_idx[NT];
  reg<HYST_WIDTH> prefetch_hyst[NT];
  reg<U_WIDTH> prefetch_u[NT];
  reg<TAG_WIDTH> prefetch_ctag[NT];
  reg<GROUP_BITS> prefetch_group_id[NT]; // extracted from upper tag bits

  // Prediction regs (shared by predict1/2/reuse)
  arr<reg<1>, N> pred;

  // Meta pipeline (per-bank)
  reg<META_WIDTH, i64> meta_pipe[META_BANKS][META_PIPE];
  reg<META_IDX_BITS> meta_idx_pipe[META_BANKS][META_PIPE];

  // ====================================================================
  // RAMs — per-table clusters matching TATable declaration order:
  //   tag_ram, pred_ram, sec_ram, fold_idx, fold_tag, zone, hyst_ram, u_ram
  //
  // This places fold regs at their table's sec_ram location (co-located
  // with the predict-path RAMs they feed), and hyst/u in a separate
  // zone after each table's predict cluster.
  // ====================================================================

  // ============================================================
  // Per-table zone layout (IR-optimized):
  //   tag, fold_idx, fold_tag, sec, zone, pred, hyst, u
  // With 1-bit pred_rams (tiny), move them into the update zone
  // to let the predict cluster (tag/fold/sec) pack tightly.
  // pred_ram is read in predict1 (ahead pipeline) so the extra
  // wiring from the update zone is small — 1-bit value travels.
  // ============================================================

  // ---- Table 0 (1024 entries, IDX=10) ----
  hcm::ram<val<TAG_WIDTH>, 1024> tag_ram0{"t0_tag"};
  ta_folded_gh<10> fold_idx0;
  ta_folded_gh<TAG_WIDTH> fold_tag0;
  hcm::ram<val<SEC_TAG_BITS>, 1024> sec_ram0{"t0_sec"};
  hcm::zone zone0;
  ta_rwram<PRED_BITS, 1024, 8, 1> pred_ram0{"t0_pred"};
  ta_rwram<HYST_WIDTH, 512, 8, 1> hyst_ram0{"t0_hyst"};
  ta_rwram<U_WIDTH, 1024, 8, 1> u_ram0{"t0_u"};

  // ---- Table 1 (1024 entries, IDX=10) ----
  hcm::ram<val<TAG_WIDTH>, 1024> tag_ram1{"t1_tag"};
  ta_folded_gh<10> fold_idx1;
  ta_folded_gh<TAG_WIDTH> fold_tag1;
  hcm::ram<val<SEC_TAG_BITS>, 1024> sec_ram1{"t1_sec"};
  hcm::zone zone1;
  ta_rwram<PRED_BITS, 1024, 8, 1> pred_ram1{"t1_pred"};
  ta_rwram<HYST_WIDTH, 512, 8, 1> hyst_ram1{"t1_hyst"};
  ta_rwram<U_WIDTH, 1024, 8, 1> u_ram1{"t1_u"};

  // ---- Table 2 (1024 entries, IDX=10) ----
  hcm::ram<val<TAG_WIDTH>, 1024> tag_ram2{"t2_tag"};
  ta_folded_gh<10> fold_idx2;
  ta_folded_gh<TAG_WIDTH> fold_tag2;
  hcm::ram<val<SEC_TAG_BITS>, 1024> sec_ram2{"t2_sec"};
  hcm::zone zone2;
  ta_rwram<PRED_BITS, 1024, 8, 1> pred_ram2{"t2_pred"};
  ta_rwram<HYST_WIDTH, 512, 8, 1> hyst_ram2{"t2_hyst"};
  ta_rwram<U_WIDTH, 1024, 8, 1> u_ram2{"t2_u"};

  // ---- Table 3 (1024 entries, IDX=10) ----
  hcm::ram<val<TAG_WIDTH>, 1024> tag_ram3{"t3_tag"};
  ta_folded_gh<10> fold_idx3;
  ta_folded_gh<TAG_WIDTH> fold_tag3;
  hcm::ram<val<SEC_TAG_BITS>, 1024> sec_ram3{"t3_sec"};
  hcm::zone zone3;
  ta_rwram<PRED_BITS, 1024, 8, 1> pred_ram3{"t3_pred"};
  ta_rwram<HYST_WIDTH, 512, 8, 1> hyst_ram3{"t3_hyst"};
  ta_rwram<U_WIDTH, 1024, 8, 1> u_ram3{"t3_u"};

  // ---- Table 4 (1024 entries, IDX=10) ----
  hcm::ram<val<TAG_WIDTH>, 1024> tag_ram4{"t4_tag"};
  ta_folded_gh<10> fold_idx4;
  ta_folded_gh<TAG_WIDTH> fold_tag4;
  hcm::ram<val<SEC_TAG_BITS>, 1024> sec_ram4{"t4_sec"};
  hcm::zone zone4;
  ta_rwram<PRED_BITS, 1024, 8, 1> pred_ram4{"t4_pred"};
  ta_rwram<HYST_WIDTH, 512, 8, 1> hyst_ram4{"t4_hyst"};
  ta_rwram<U_WIDTH, 1024, 8, 1> u_ram4{"t4_u"};

  // ---- Table 14 (2048 entries, IDX=11) — placed here to fill gap ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram14{"t14_tag"};
  ta_folded_gh<11> fold_idx14;
  ta_folded_gh<TAG_WIDTH> fold_tag14;
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram14{"t14_sec"};
  hcm::zone zone14;
  ta_rwram<PRED_BITS, 2048, 8, 1> pred_ram14{"t14_pred"};
  ta_rwram<HYST_WIDTH, 1024, 8, 1> hyst_ram14{"t14_hyst"};
  ta_rwram<U_WIDTH, 2048, 8, 1> u_ram14{"t14_u"};

  // ---- Table 5 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram5{"t5_tag"};
  ta_folded_gh<11> fold_idx5;
  ta_folded_gh<TAG_WIDTH> fold_tag5;
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram5{"t5_sec"};
  hcm::zone zone5;
  ta_rwram<PRED_BITS, 2048, 8, 1> pred_ram5{"t5_pred"};
  ta_rwram<HYST_WIDTH, 1024, 8, 1> hyst_ram5{"t5_hyst"};
  ta_rwram<U_WIDTH, 2048, 8, 1> u_ram5{"t5_u"};

  // ---- Table 6 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram6{"t6_tag"};
  ta_folded_gh<11> fold_idx6;
  ta_folded_gh<TAG_WIDTH> fold_tag6;
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram6{"t6_sec"};
  hcm::zone zone6;
  ta_rwram<PRED_BITS, 2048, 8, 1> pred_ram6{"t6_pred"};
  ta_rwram<HYST_WIDTH, 1024, 8, 1> hyst_ram6{"t6_hyst"};
  ta_rwram<U_WIDTH, 2048, 8, 1> u_ram6{"t6_u"};

  // ---- Table 7 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram7{"t7_tag"};
  ta_folded_gh<11> fold_idx7;
  ta_folded_gh<TAG_WIDTH> fold_tag7;
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram7{"t7_sec"};
  hcm::zone zone7;
  ta_rwram<PRED_BITS, 2048, 8, 1> pred_ram7{"t7_pred"};
  ta_rwram<HYST_WIDTH, 1024, 8, 1> hyst_ram7{"t7_hyst"};
  ta_rwram<U_WIDTH, 2048, 8, 1> u_ram7{"t7_u"};

  // ---- Fallback (single 8K×7 RAM) ----
  hcm::ram<val<FB_PRED_BITS>, FB_CAPACITY> fb_ctr{"fb"};

  // ---- Table 8 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram8{"t8_tag"};
  ta_folded_gh<11> fold_idx8;
  ta_folded_gh<TAG_WIDTH> fold_tag8;
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram8{"t8_sec"};
  hcm::zone zone8;
  ta_rwram<PRED_BITS, 2048, 8, 1> pred_ram8{"t8_pred"};
  ta_rwram<HYST_WIDTH, 1024, 8, 1> hyst_ram8{"t8_hyst"};
  ta_rwram<U_WIDTH, 2048, 8, 1> u_ram8{"t8_u"};

  // ---- Table 9 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram9{"t9_tag"};
  ta_folded_gh<11> fold_idx9;
  ta_folded_gh<TAG_WIDTH> fold_tag9;
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram9{"t9_sec"};
  hcm::zone zone9;
  ta_rwram<PRED_BITS, 2048, 8, 1> pred_ram9{"t9_pred"};
  ta_rwram<HYST_WIDTH, 1024, 8, 1> hyst_ram9{"t9_hyst"};
  ta_rwram<U_WIDTH, 2048, 8, 1> u_ram9{"t9_u"};

  // ---- Table 10 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram10{"t10_tag"};
  ta_folded_gh<11> fold_idx10;
  ta_folded_gh<TAG_WIDTH> fold_tag10;
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram10{"t10_sec"};
  hcm::zone zone10;
  ta_rwram<PRED_BITS, 2048, 8, 1> pred_ram10{"t10_pred"};
  ta_rwram<HYST_WIDTH, 1024, 8, 1> hyst_ram10{"t10_hyst"};
  ta_rwram<U_WIDTH, 2048, 8, 1> u_ram10{"t10_u"};

  // ---- Table 11 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram11{"t11_tag"};
  ta_folded_gh<11> fold_idx11;
  ta_folded_gh<TAG_WIDTH> fold_tag11;
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram11{"t11_sec"};
  hcm::zone zone11;
  ta_rwram<PRED_BITS, 2048, 8, 1> pred_ram11{"t11_pred"};
  ta_rwram<HYST_WIDTH, 1024, 8, 1> hyst_ram11{"t11_hyst"};
  ta_rwram<U_WIDTH, 2048, 8, 1> u_ram11{"t11_u"};

  // ---- Table 12 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram12{"t12_tag"};
  ta_folded_gh<11> fold_idx12;
  ta_folded_gh<TAG_WIDTH> fold_tag12;
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram12{"t12_sec"};
  hcm::zone zone12;
  ta_rwram<PRED_BITS, 2048, 8, 1> pred_ram12{"t12_pred"};
  ta_rwram<HYST_WIDTH, 1024, 8, 1> hyst_ram12{"t12_hyst"};
  ta_rwram<U_WIDTH, 2048, 8, 1> u_ram12{"t12_u"};

  // ---- Table 13 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram13{"t13_tag"};
  ta_folded_gh<11> fold_idx13;
  ta_folded_gh<TAG_WIDTH> fold_tag13;
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram13{"t13_sec"};
  hcm::zone zone13;
  ta_rwram<PRED_BITS, 2048, 8, 1> pred_ram13{"t13_pred"};
  ta_rwram<HYST_WIDTH, 1024, 8, 1> hyst_ram13{"t13_hyst"};
  ta_rwram<U_WIDTH, 2048, 8, 1> u_ram13{"t13_u"};

  // ---- Meta (update-only, 2 banks) ----
  hcm::zone meta_zone;
  ta_rwram<META_WIDTH, META_CAPACITY, 8, 1> meta_ctr0{"meta0"};
  ta_rwram<META_WIDTH, META_CAPACITY, 8, 1> meta_ctr1{"meta1"};

  // ---- FB bimodal hysteresis (update-only, non-critical path) ----
  hcm::ram<val<FB_PRED_BITS>, FB_BH_CAPACITY> fb_bim_hyst{"fb_bim_hyst"};

  // ======== Training regs ========
  reg<MAX_IDX_BITS> train_idx[NT];
  reg<PRED_BITS> train_pred[NT];   // 1-bit per entry
  reg<HYST_WIDTH> train_hyst[NT];
  reg<U_WIDTH> train_u[NT];
  reg<FB_PRED_BITS> train_fb;
  reg<FB_IDX_BITS> train_fb_idx;
  reg<FB_PRED_BITS> fb_bh_read;
  reg<ALLOC_PC_BITS> train_pc;
  reg<TAG_WIDTH> train_ctag[NT];
  reg<SEC_TAG_BITS> train_sec_tag;
  reg<GROUP_BITS> train_group_id[NT]; // piped from current_group_id

  // Per-group resolution results (NUM_GROUPS=7 independent chains)
  reg<MATCH_BITS> train_match1[NUM_GROUPS];
  reg<PRED_BITS> train_provider_pred[NUM_GROUPS];
  reg<1> train_provider_weak[NUM_GROUPS];
  reg<1> train_altdiff[NUM_GROUPS];
  reg<1> train_valid = 0;

  // Piped tag/sec hit for decay + sibling
  reg<1> train_tag_hit[NT];  // htag-only match (piped)
  reg<1> train_sec_hit[NT];

  // Global pressure counters
  reg<ACC_WIDTH> acc_ctr;
  reg<ALLOC_WIDTH> alloc_ctr;

  // ======== Per-table RAM dispatch (compile-time index → named member) ========
#define HC_DISPATCH(name)                                                      \
  template <u64 I> auto &name##_at() {                                         \
    static_assert(I < NT);                                                     \
    if constexpr (I == 0) return name##0;                                      \
    else if constexpr (I == 1) return name##1;                                 \
    else if constexpr (I == 2) return name##2;                                 \
    else if constexpr (I == 3) return name##3;                                 \
    else if constexpr (I == 4) return name##4;                                 \
    else if constexpr (I == 5) return name##5;                                 \
    else if constexpr (I == 6) return name##6;                                 \
    else if constexpr (I == 7) return name##7;                                 \
    else if constexpr (I == 8) return name##8;                                 \
    else if constexpr (I == 9) return name##9;                                 \
    else if constexpr (I == 10) return name##10;                               \
    else if constexpr (I == 11) return name##11;                               \
    else if constexpr (I == 12) return name##12;                               \
    else if constexpr (I == 13) return name##13;                               \
    else return name##14;                                                      \
  }
  HC_DISPATCH(tag_ram)
  HC_DISPATCH(pred_ram)
  HC_DISPATCH(sec_ram)
  HC_DISPATCH(fold_idx)
  HC_DISPATCH(fold_tag)
  HC_DISPATCH(hyst_ram)
  HC_DISPATCH(u_ram)
#undef HC_DISPATCH

  // Per-table HYST_IDX_BITS
  static constexpr std::array<u64, NT> HYST_IDX_BITS = {
    9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
  };

  // ======== Helpers ========
  val<1> line_end() { return (block_entry + block_size) == hard<LINEINST>{}; }

  // ======== predict1 ========
  val<1> predict1([[maybe_unused]] val<64> inst_pc) {
    // inst_pc fanout: 2 per table (>>2, >>5) + fb(>>2) + prefetch_pc(>>2)
    //               + sec_tag(>>2,>>9,>>16) + meta_idx(>>2)
    inst_pc.fanout(hard<2 * NT + 6>{});

    // ================================================================
    // Stage 1: Meta pipeline shift (moved from update_cycle)
    //
    // Shift meta_pipe[] down, then read new meta counter value into [0].
    // meta_pipe[META_PIPE-1] (= [1]) holds the delayed meta value used
    // for provider-vs-alt selection in the resolution chain.
    //
    // Uses inst_pc for meta index (inst_pc(B) == next_pc(A)).
    // Shift FIRST, then write [0] — prevents stale-value bug.
    // ================================================================

    for (u64 mb = 0; mb < META_BANKS; mb++) {
      for (u64 i = META_PIPE - 1; i > 0; i--) {
        meta_pipe[mb][i] = meta_pipe[mb][i - 1].fo1();
        meta_idx_pipe[mb][i] = meta_idx_pipe[mb][i - 1].fo1();
      }
    }
    {
      auto meta_idx = val<META_IDX_BITS>{inst_pc >> 2};
      meta_idx.fanout(hard<2 + META_BANKS>{}); // 2 RAM reads + 2 idx_pipe writes
      meta_pipe[0][0] = meta_ctr0.read(meta_idx);
      meta_pipe[1][0] = meta_ctr1.read(meta_idx);
      meta_idx_pipe[0][0] = meta_idx;
      meta_idx_pipe[1][0] = meta_idx;
    }
    // meta_pipe[][1]: read twice — meta_use_alt here + old_meta in training
    meta_pipe[0][META_PIPE - 1].fanout(hard<2>{});
    meta_pipe[1][META_PIPE - 1].fanout(hard<2>{});
    // Per-bank meta_use_alt >= 0 means "trust alt over provider when provider is weak"
    arr<val<1>, META_BANKS> meta_use_alt = [&](u64 mb) -> val<1> {
      return val<META_WIDTH, i64>{meta_pipe[mb][META_PIPE - 1]} >= hard<0>{};
    };

    // ================================================================
    // Stage 2: Sec-tag + resolution (on OLD prefetch_* from previous predict1)
    //
    // inst_pc(B) == next_pc(A): the sec_tag hash is the same as what
    // update_cycle would have computed from block_end_info.next_pc.
    // Resolution reads prefetch_* (block A's entries) before they are
    // overwritten by this cycle's RAM reads.
    //
    // sec_tag_now also serves as train_sec_hit (same comparison).
    // ================================================================

    // Fanout for prefetch_* registers: resolution reads + shift fo1()
    for (u64 i = 0; i < NT; i++) {
      prefetch_sec[i].fanout(hard<2>{});       // sec_matches + shift
      prefetch_tag_hit[i].fanout(hard<2>{});   // htag_sec + shift
      prefetch_hyst[i].fanout(hard<2>{});      // weak_mask + shift
      prefetch_u[i].fanout(hard<2>{});         // weak_mask + shift
      prefetch_group_id[i].fanout(hard<1 + NUM_GROUPS>{}); // shift + 7 group== checks
      prefetch_pred[i].fanout(hard<1 + NUM_GROUPS>{});     // shift + 7 resolve chains
    }
    prefetch_fb.fanout(hard<2>{});  // fb_bits make_array + shift

    auto sec_tag_now = ta::Xor3SecTagHash5::apply<SEC_TAG_BITS>(inst_pc);
    // Fanout: NT sec_match + NT train_sec_hit + train_sec_tag write
    sec_tag_now.fanout(hard<2 * NT + 2>{});
    curr_sec_tag = sec_tag_now;

    // Precompute per-table sec_tag match
    arr<val<1>, NT> sec_matches = [&](u64 i) -> val<1> {
      return val<SEC_TAG_BITS>{prefetch_sec[i]} == val<SEC_TAG_BITS>{sec_tag_now};
    };

    // Precompute per-table combined htag_hit & sec_match
    arr<val<1>, NT> htag_sec = [&](u64 i) -> val<1> {
      return val<1>{prefetch_tag_hit[i]} & sec_matches[i].fo1();
    };
    htag_sec.fanout(hard<NUM_GROUPS>{}); // one read per group chain

    // Precompute per-group fb bits (split N-wide fb into individual bits)
    arr<val<1>, NUM_GROUPS> fb_bits =
        val<FB_PRED_BITS>{prefetch_fb}.make_array(val<1>{});

    // Precompute per-table weakness: entry is weak when hyst==0 AND u==0.
    val<NT> weak_mask = [&]() {
      arr<val<1>, NT> w = [&](u64 i) -> val<1> {
        return val<1>{val<HYST_WIDTH>{prefetch_hyst[i]} == hard<0>{}} &
               val<1>{val<U_WIDTH>{prefetch_u[i]} == hard<0>{}};
      };
      return w.fo1().concat();
    }();
    weak_mask.fanout(hard<NUM_GROUPS>{}); // one read per group chain

    // meta_use_alt: each bank serves its groups
    // Bank 0: groups 0-3 (4 reads), Bank 1: groups 4-6 (3 reads)
    meta_use_alt.fanout(hard<4>{}); // max 4 reads per bank (groups 0-3 for bank 0)

    // ---- resolve_chain lambda (1-bit version) ----
    auto resolve_chain = [&](arr<val<1>, NT> fh, val<PRED_BITS> fb_g,
                             u64 meta_bank = 0) {
      val<MATCH_BITS> match = concat(val<1>{1}, fh.fo1().concat());
      match.fanout(hard<3>{});

      val<1> ha = val<NT>{match} != hard<0>{};

      val<MATCH_BITS> match1 = match.one_hot();
      match1.fanout(hard<4>{}); // make_array + XOR + pw + train save

      val<MATCH_BITS> remainder = match ^ match1;
      val<MATCH_BITS> match2 = remainder.fo1().one_hot();

      // 1-bit table predictions + 1-bit fb
      fb_g.fanout(hard<2>{}); // min fanout; read once (when i==NT)
      arr<val<PRED_BITS>, NT + 1> table_preds = [&](u64 i) -> val<PRED_BITS> {
        if (i < NT)
          return val<PRED_BITS>{prefetch_pred[i]};
        return val<PRED_BITS>{fb_g};
      };
      for (u64 i = 0; i < NT + 1; i++) table_preds[i].fanout(hard<2>{});

      arr<val<1>, NT + 1> m1_bits = match1.make_array(val<1>{});
      match2.fanout(hard<2>{});
      arr<val<1>, NT + 1> m2_bits = match2.make_array(val<1>{});

      // 1-bit masking: m1/m2 bit directly ANDs with 1-bit pred
      arr<val<1>, NT + 1> pmask = [&](u64 i) {
        return m1_bits[i].fo1() & table_preds[i];
      };
      arr<val<1>, NT + 1> amask = [&](u64 i) {
        return m2_bits[i].fo1() & table_preds[i];
      };
      val<1> pp = pmask.fo1().fold_or();
      val<1> ap = amask.fo1().fold_or();
      pp.fanout(hard<2>{});
      ap.fanout(hard<2>{});

      val<1> pw = (val<NT>{match1} & weak_mask) != hard<0>{};
      pw.fanout(hard<2>{});

      val<1> ua = pw & meta_use_alt[meta_bank] & ha.fo1();

      // 1-bit: XOR is already the disagreement bit
      auto pp_xor_ap = pp ^ ap;
      val<1> ad = pp_xor_ap.fo1();

      return std::tuple{val<MATCH_BITS>{match1},
                        val<MATCH_BITS>{match2},
                        val<PRED_BITS>{pp},
                        val<PRED_BITS>{ap},
                        val<1>{pw},
                        ad.fo1(),
                        ua.fo1()};
    };

    // ---- Per-group resolution loop ----
    static_loop<NUM_GROUPS>([&]<u64 G>() {
      // Group hits: precomputed (htag & sec) & group_id match
      arr<val<1>, NT> group_hits = [&](u64 i) -> val<1> {
        val<1> group_ok =
            val<GROUP_BITS>{prefetch_group_id[i]} == hard<G>{};
        return htag_sec[i] & group_ok.fo1();
      };

      // Per-group fb bit (precomputed from make_array split)
      val<PRED_BITS> fb_g = fb_bits[G].fo1();

      auto [m1, m2, pp, ap, pw, ad, ua] =
          resolve_chain(group_hits.fo1(), fb_g.fo1(), GROUP_TO_META_BANK[G]);

      // pp is read twice: select + train save
      pp.fanout(hard<2>{});

      // Scatter: pred[G] = select(ua, ap, pp) — each group has 1 branch
      pred[G] = select(ua.fo1(), ap.fo1(), pp);

      // Save per-group resolution → train regs
      train_match1[G] = m1.fo1();
      train_provider_pred[G] = pp;
      train_provider_weak[G] = pw.fo1();
      train_altdiff[G] = ad.fo1();
    });

#ifdef DEBUG_PRINT
    for (u64 i = 0; i < N; i++)
      pred[i].print(("  pred[" + std::to_string(i) + "]=").c_str(), "\n",
                    true, std::cerr);
#endif

    // ================================================================
    // Stage 3: Pipeline shift (prefetch → train)
    //
    // Save OLD prefetch_* into train_* BEFORE RAM reads overwrite them.
    // train_sec_hit uses sec_tag_now (same as resolution sec_matches).
    // ================================================================

    train_sec_tag = sec_tag_now;
    static_loop<NT>([&]<u64 I>() {
      train_sec_hit[I] = (val<SEC_TAG_BITS>{prefetch_sec[I].fo1()} ==
                          val<SEC_TAG_BITS>{sec_tag_now});
      train_idx[I] = prefetch_idx[I].fo1();
      train_pred[I] = prefetch_pred[I].fo1();
      train_hyst[I] = prefetch_hyst[I].fo1();
      train_u[I] = prefetch_u[I].fo1();
      train_ctag[I] = prefetch_ctag[I].fo1();
      train_tag_hit[I] = prefetch_tag_hit[I].fo1();
      train_group_id[I] = prefetch_group_id[I].fo1();
    });
    train_fb = prefetch_fb.fo1();
    train_fb_idx = prefetch_fb_idx.fo1();
    train_pc = prefetch_pc.fo1();

    // ================================================================
    // Stage 4: Ahead RAM reads for next block
    // ================================================================

    static_loop<NT>([&]<u64 I>() {
      auto &fi = fold_idx_at<I>();
      auto &ft = fold_tag_at<I>();
      fi.fanout(hard<2>{}); // get() + compute_update
      ft.fanout(hard<2>{}); // get() + compute_update
      auto fold_idx_val = fi.get();
      auto idx = fold_idx_val.fo1() ^ val<IDX_BITS[I]>{inst_pc >> 2};
      idx.fanout(hard<6>{}); // 5 RAM reads + prefetch_idx write
      auto fold_tag_val = ft.get();
      auto computed_tag = fold_tag_val.fo1() ^ val<TAG_WIDTH>{inst_pc >> 5};
      computed_tag.fanout(hard<2>{}); // tag comparison + prefetch_ctag write

      auto stored_tag = tag_ram_at<I>().read(idx);
      stored_tag.fanout(hard<2>{}); // htag compare + group_id extract
      // Split tag: lower HTAG_WIDTH bits = htag, upper GROUP_BITS = group_id
      prefetch_tag_hit[I] =
          val<PER_TABLE_HTAG[I]>{stored_tag} == val<PER_TABLE_HTAG[I]>{computed_tag};
      prefetch_group_id[I] =
          val<GROUP_BITS>{stored_tag >> PER_TABLE_HTAG[I]};
      prefetch_ctag[I] = computed_tag;
      prefetch_pred[I] = pred_ram_at<I>().read(idx);
      prefetch_sec[I] = sec_ram_at<I>().read(idx);
      prefetch_idx[I] = idx;
      prefetch_hyst[I] = hyst_ram_at<I>().read(val<HYST_IDX_BITS[I]>{idx});
      prefetch_u[I] = u_ram_at<I>().read(idx);
#ifdef DEBUG_PRINT
      if constexpr (I == 0) {
        std::cerr << "\n=== predict1 table[0] ===\n";
        fold_idx_val.print("  fold_idx=", "\n", true, std::cerr);
        idx.print("  idx=", "\n", true, std::cerr);
        fold_tag_val.print("  fold_tag=", "\n", true, std::cerr);
        computed_tag.print("  ctag=", "\n", true, std::cerr);
        stored_tag.print("  stored_tag=", "\n", true, std::cerr);
        prefetch_tag_hit[I].print("  tag_hit=", "\n", true, std::cerr);
        prefetch_pred[I].print("  pred=", "\n", true, std::cerr);
        prefetch_sec[I].print("  sec=", "\n", true, std::cerr);
        prefetch_hyst[I].print("  hyst=", "\n", true, std::cerr);
      }
#endif
    });

    // Fallback ahead read (single 8K×7 RAM)
    auto fb_idx = val<FB_IDX_BITS>{inst_pc >> 2};
    fb_idx.fanout(hard<2>{}); // prefetch_fb_idx + fb read
    prefetch_fb_idx = fb_idx;
    prefetch_fb = fb_ctr.read(fb_idx);
    prefetch_pc = val<ALLOC_PC_BITS>{inst_pc >> 2};

    // Return prediction from resolution (just computed above)
    block_entry.fanout(hard<2 * LINEINST>{});
    pred.fanout(hard<2 * LINEINST>{});
    block_size = 1;
    num_branch = 0;
    reuse_prediction(~line_end());
    return pred[num_branch];
  }

  // ======== reuse_predict1 ========
  val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) {
    block_size++;
    reuse_prediction(~line_end());
    return pred[num_branch];
  }

  // ======== predict2 / reuse_predict2 ========
  val<1> predict2([[maybe_unused]] val<64> inst_pc) { return pred[num_branch]; }
  val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) { return pred[num_branch]; }

  // ======== update_condbr ========
  void update_condbr([[maybe_unused]] val<64> branch_pc, val<1> taken,
                     [[maybe_unused]] val<64> next_pc) {
    assert(num_branch < N);
    branch_dir[num_branch] = taken.fo1();
    num_branch++;
    reuse_prediction(~line_end() & val<1>{num_branch < N});
  }

  // ======== update_cycle ========
  void update_cycle([[maybe_unused]] instruction_info &block_end_info) {

    // ================================================================
    // Resolution + meta pipeline + shift all moved to predict1.
    // update_cycle only does: training + history updates.
    // ================================================================

    // Alloc-path: NT sec_ram writes use the sec_tag saved by predict1
    train_sec_tag.fanout(hard<NT>{});

    // ================================================================
    // Training setup
    //
    // Pre-read values saved by predict1's resolution loop.
    //
    // Fanout declarations for training regs:
    //   train_hyst[I]:    2 reads (t_hyst_weak + old_hyst in update)
    //   train_u[I]:       3 reads (u_zero + old_u base + old_u decay)
    //   train_tag_hit[I]: 2 reads (sibling skip + decay miss)
    //   train_sec_hit[I]: 2 reads (sibling skip + decay miss)
    //   train_idx[I]:     5 reads (pred/hyst/tag/sec/u RAM writes)
    // ================================================================

    // Pre-read values saved by predict1's per-group resolution loop.
    u64 train_group = num_branch > 0 ? num_branch - 1 : 0;
    arr<val<PRED_BITS>, NUM_GROUPS> old_provider_pred = [&](u64 g) -> val<PRED_BITS> {
      return val<PRED_BITS>{train_provider_pred[g].fo1()};
    };
    val<MATCH_BITS> pre_read_match1 = train_match1[train_group].fo1();
    arr<val<1>, NUM_GROUPS> pre_read_pw = [&](u64 g) -> val<1> {
      return train_provider_weak[g].fo1();
    };
    for (u64 g = 0; g < NUM_GROUPS; g++)
      pre_read_pw[g].fanout(hard<2>{}); // gate + dir
    arr<val<1>, NUM_GROUPS> pre_read_ad = [&](u64 g) -> val<1> {
      return train_altdiff[g].fo1();
    };
    for (u64 g = 0; g < NUM_GROUPS; g++)
      pre_read_ad[g].fanout(hard<3>{}); // t_ad + gate + dir

    branch_dir.fanout(hard<5>{}); // actual_dir_concat + per_group_wrong + entry_actual + true_block + last_condbr

    static_loop<NT>([&]<u64 I>() {
      train_hyst[I].fanout(hard<2>{});
      train_u[I].fanout(hard<3>{});       // u_zero + base old_u + decay old_u
      train_tag_hit[I].fanout(hard<2>{}); // sibling + decay
      train_sec_hit[I].fanout(hard<2>{}); // sibling + decay
      train_idx[I].fanout(hard<5>{});     // 5 RAM writes
    });

    // Use pre-read values from predict1's per-group resolution
    val<MATCH_BITS> t_match1 = pre_read_match1.fo1();
    t_match1.fanout(hard<2>{}); // make_array + alloc_base
    arr<val<1>, NT + 1> t_m1 = t_match1.make_array(val<1>{});
    static_loop<NT>([&]<u64 I>() {
      t_m1[I].fanout(hard<4>{});
    });
    // t_m1[NT]: fb_gate(1) + bh_gate(1)
    t_m1[NT].fanout(hard<2>{});

    // t_pw not needed — per-bank meta uses pre_read_pw[g] directly
    val<1> t_ad = pre_read_ad[train_group];

    // Hyst-only weakness: provider has hyst==0 (regardless of u).
    arr<val<1>, NT + 1> t_hyst_weak_arr = [&](u64 i) -> val<1> {
      if (i < NT)
        return t_m1[i] & val<1>{val<HYST_WIDTH>{train_hyst[i]} == hard<0>{}};
      return val<1>{0};
    };
    val<1> t_phw = t_hyst_weak_arr.fo1().fold_or();

    // Per-group train regs already saved in per-group loop above.
    // Guard: skip training until piped resolution regs have been populated
    val<1> do_train = train_valid.fo1();
    train_valid = 1;

    // ================================================================
    // Stage 5a: No conditional branches → update history, skip training
    //
    // When the block had no conditional branches (only unconditional),
    // we still need to update the folded histories and global history
    // with path bits, but there's nothing to train on.
    // ================================================================

    if (num_branch == 0) {
      val<PATHBITS> path_bits =
          val<PATHBITS>{block_end_info.next_pc.fo1() >> 2};
      path_bits.fanout(hard<NT * 2 + 1>{});
      gh.template fanout_per_bit<GH_FANOUT>();
      static_loop<NT>([&]<u64 I>() {
        auto &fi = fold_idx_at<I>();
        auto &ft = fold_tag_at<I>();
        fi.apply_update(fi.compute_update(
            gh, hard<HIST_LEN[I]>{}, path_bits));
        ft.apply_update(ft.compute_update(
            gh, hard<HIST_LEN[I]>{}, path_bits));
      });
      gh.update(path_bits);
      last_condbr_dir = 0;
      true_block = 1;
      return;
    }

    // ================================================================
    // Stage 6: Training — correctness signals
    //
    // mispredict: framework signal, gates extra_cycle + allocation
    // actual_dir: packed branch directions from update_condbr
    // any_provider_wrong: provider prediction != actual on any branch
    // ================================================================

    val<1> &mispredict = block_end_info.is_mispredict;
    // extra_cycle(1) + alloc_trigger(1) + fb_gate(1) + acc_ctr(1) + true_block(1) + bh_read(1) + bh_gate(1)
    mispredict.fanout(hard<7>{});

    // Read fb bimodal hyst BEFORE extra_cycle (cycle 1); write happens AFTER (cycle 2)
    train_fb_idx.fanout(hard<3>{});  // bh_read(1) + fb_write(1) + bh_write(1)
    fb_bh_read = execute_if(mispredict, [&]() {
      return fb_bim_hyst.read(val<FB_BH_IDX_BITS>{train_fb_idx});
    });

    need_extra_cycle(mispredict);

    // do_train gates all RAM writes: 4 per table + fb + bh + meta
    // fb(1) + bh(1) + meta_banks(2) + 4 per table
    do_train.fanout(hard<4 * NT + 2 + META_BANKS>{});

    // Full N-wide actual direction for fb path
    val<FB_PRED_BITS> actual_dir_full = arr<val<1>, N>{[&](u64 i) -> val<1> {
                                          return val<1>{branch_dir[i]};
                                        }}.concat();
    // fb_wrong(1) + fb_changed(1) + bh_write(1)
    actual_dir_full.fanout(hard<3>{});

    // Per-group provider wrong (used by per-bank meta training)
    arr<val<1>, NUM_GROUPS> per_group_wrong = [&](u64 g) -> val<1> {
      return (val<PRED_BITS>{old_provider_pred[g].fo1()} ^
              val<1>{branch_dir[g]}) != hard<0>{};
    };

    // t_pw no longer used directly (per-bank meta uses pre_read_pw per-group)
    t_phw.fanout(hard<NT>{});  // do_pred_update per table
    t_ad.fanout(hard<NT>{});   // u_correct per table

    // ================================================================
    // Stage 7: Allocation
    //
    // On mispredict, allocate in tables above provider with u==0.
    // AllocPressSkip: probabilistically skip the closest candidate
    // based on alloc pressure (alloc_ctr > rng → skip).
    // SiblingAll: skip entries with same primary tag but different
    // sec_tag (siblings) to avoid displacing related entries.
    // ================================================================

    // 7a. Trigger: MISPREDICT policy
    val<1> alloc_trigger = mispredict;
    val<NT> triggermask = alloc_trigger.fo1().replicate(hard<NT>{}).concat();

    // 7b. Action: STANDARD (no probabilistic gating)
    val<8> alloc_rng = val<8>{static_cast<u64>(std::rand()) & 0xFF};
    val<NT> gated_triggermask = triggermask.fo1();

    // 7c. Candidate mask: tables above provider (alloc_base) with u==0
    val<NT> alloc_base = val<NT>{t_match1 - 1};
    arr<val<1>, NT> u_zero = [&](u64 i) -> val<1> {
      return val<U_WIDTH>{train_u[i]} == hard<0>{};
    };
    val<NT> notumask = u_zero.fo1().concat();
    val<NT> postmask = alloc_base.fo1() & gated_triggermask.fo1();
    postmask.fanout(hard<2>{}); // candallocmask + uclearmask (U_CLEAR=DECREMENT)

    // 7d. Sibling skip (SiblingPolicy::ALL, FLOOR=0, group-aware)
    // Group-aware sibling: sibling if htag matches but
    // (group differs OR sec_tag differs).
    // not_sibling = ~htag_hit | (group_match & sec_hit)
    val<GROUP_BITS> alloc_gid_sib =
        val<GROUP_BITS>{u64(train_group)};
    alloc_gid_sib.fanout(hard<NT>{});
    val<NT> candallocmask = [&]() {
      val<NT> base = postmask & notumask.fo1();
      arr<val<1>, NT> not_sibling = [&](u64 i) -> val<1> {
        auto th = val<1>{train_tag_hit[i]};
        auto sh = val<1>{train_sec_hit[i]};
        auto gid = val<GROUP_BITS>{train_group_id[i].fo1()};
        val<1> group_match = gid.fo1() == val<GROUP_BITS>{alloc_gid_sib};
        return ~th.fo1() | (group_match.fo1() & sh.fo1());
      };
      val<NT> sibling_mask = not_sibling.fo1().concat();
      return base.fo1() & sibling_mask.fo1();
    }();
    candallocmask.fanout(hard<2>{}); // collamask + noalloc

    // Pressure counter fanout: collamask/target(1) + ta_update_ctr(1) = 2
    acc_ctr.fanout(hard<2>{});
    alloc_ctr.fanout(hard<2>{});

    // 7e. Target policy: AllocPressureSkipTarget<1>
    // Probabilistically skips the closest candidate when alloc pressure is high.
    // alloc_rng: single read, use fo1()
    val<NT> collamask =
        AllocPressureSkipTarget<1>::apply<NT>(
            candallocmask.reverse(), val<ALLOC_WIDTH>{alloc_ctr},
            val<ACC_WIDTH>{acc_ctr}, alloc_rng.fo1());

    // 7f. Final allocation: one-hot pick (MAX_ALLOC=1)
    arr<val<1>, NT> allocate =
        collamask.fo1().one_hot().reverse().make_array(val<1>{});
    // do_alloc(1) + alloc_target(1) + base_u_write(1) + decay(1)
    allocate.fanout(hard<4>{});

    val<NT> alloc_target = [&]() {
      arr<val<1>, NT> a = allocate;
      return a.fo1().concat();
    }();
    alloc_target.fanout(hard<2>{}); // != hard<0> check + any_alloc

    // 7g. uclear: on alloc failure, decrement u-bits above provider
    // UClearPolicy::DECREMENT — saturating decrement on failed alloc
    arr<val<1>, NT> uclear = [&]() -> arr<val<1>, NT> {
      val<1> noalloc = (candallocmask == hard<0>{});
      val<NT> uclearmask =
          postmask & noalloc.fo1().replicate(hard<NT>{}).concat();
      return uclearmask.fo1().make_array(val<1>{});
    }();

    // ================================================================
    // Stage 8: Fallback update
    //
    // When fallback is the provider AND we mispredicted, overwrite
    // the fb_ctr entry with actual_dir. FB_RECONCILE=false: no
    // agreement tracking or silent overwrites.
    // ================================================================

    // fb_write_data: only flip prediction where hyst is weak AND wrong
    train_fb.fanout(hard<3>{}); // fb_wrong(1) + fb_write_data(1) + fb_changed(1)
    fb_bh_read.fanout(hard<2>{}); // fb_bh_weak(1) + bh_changed(1)
    val<FB_PRED_BITS> fb_wrong = actual_dir_full ^ val<FB_PRED_BITS>{train_fb};
    fb_wrong.fanout(hard<2>{}); // fb_flip(1) + new_bh(1)
    val<FB_PRED_BITS> fb_bh_weak = ~val<FB_PRED_BITS>{fb_bh_read};
    val<FB_PRED_BITS> fb_flip = fb_wrong & fb_bh_weak.fo1();
    val<FB_PRED_BITS> fb_write_data = val<FB_PRED_BITS>{train_fb} ^ fb_flip.fo1();
    fb_write_data.fanout(hard<2>{}); // fb_changed compare + fb_ctr write

    val<1> fb_changed = fb_write_data != val<FB_PRED_BITS>{train_fb};
    val<1> fb_gate = do_train & t_m1[NT] & mispredict & fb_changed.fo1();
    fb_gate.fanout(hard<11>{}); // ta_rwram-style gate
    execute_if(fb_gate, [&]() {
      fb_ctr.write(val<FB_IDX_BITS>{train_fb_idx}, fb_write_data);
    });

    // Bimodal hysteresis update: correct→strong(1), wrong→weak(0)
    val<FB_PRED_BITS> new_bh = ~fb_wrong;
    new_bh.fanout(hard<2>{}); // bh_changed compare + fb_bim_hyst write
    val<1> bh_changed = new_bh != val<FB_PRED_BITS>{fb_bh_read};
    val<1> bh_gate = do_train & t_m1[NT] & mispredict & bh_changed.fo1();
    bh_gate.fanout(hard<11>{}); // ta_rwram-style gate
    execute_if(bh_gate, [&]() {
      fb_bim_hyst.write(val<FB_BH_IDX_BITS>{train_fb_idx}, new_bh);
    });

    // ================================================================
    // Stage 9: Meta counter update
    //
    // Meta tracks whether provider or alt is more accurate.
    // Updated when provider is newly-allocated (weak) AND provider/alt
    // disagree. Increment when provider was wrong (trust alt more),
    // decrement when provider was right (trust provider more).
    // ================================================================

    // Per-bank meta training with corrected gating:
    // bank_gate = OR_g(pw_g & ad_g) for groups in this bank
    // bank_wrong_dir = OR_g(pw_g & ad_g & wrong_g) — direction from gated groups
    // Bank 0: groups 0-3
    {
      arr<val<1>, 4> gate0_arr = [&](u64 j) -> val<1> {
        return pre_read_pw[j] & pre_read_ad[j];
      };
      arr<val<1>, 4> dir0_arr = [&](u64 j) -> val<1> {
        return pre_read_pw[j] & pre_read_ad[j] & per_group_wrong[j].fo1();
      };
      val<1> bank0_gate = gate0_arr.fo1().fold_or();
      val<1> bank0_dir = dir0_arr.fo1().fold_or();

      auto old_meta0 = val<META_WIDTH, i64>{meta_pipe[0][META_PIPE - 1]};
      old_meta0.fanout(hard<2>{});
      auto new_meta0 = ta_update_ctr(old_meta0, bank0_dir.fo1());
      new_meta0.fanout(hard<2>{});
      val<1> meta_gate0 = do_train & bank0_gate.fo1() & (new_meta0 != old_meta0);
      meta_gate0.fanout(hard<11>{});
      execute_if(meta_gate0, [&]() {
        meta_ctr0.write(val<META_IDX_BITS>{meta_idx_pipe[0][META_PIPE - 1].fo1()},
                        new_meta0, hard<0>{});
      });
    }
    // Bank 1: groups 4-6
    {
      arr<val<1>, 3> gate1_arr = [&](u64 j) -> val<1> {
        return pre_read_pw[j + 4] & pre_read_ad[j + 4];
      };
      arr<val<1>, 3> dir1_arr = [&](u64 j) -> val<1> {
        return pre_read_pw[j + 4] & pre_read_ad[j + 4] & per_group_wrong[j + 4].fo1();
      };
      val<1> bank1_gate = gate1_arr.fo1().fold_or();
      val<1> bank1_dir = dir1_arr.fo1().fold_or();

      auto old_meta1 = val<META_WIDTH, i64>{meta_pipe[1][META_PIPE - 1]};
      old_meta1.fanout(hard<2>{});
      auto new_meta1 = ta_update_ctr(old_meta1, bank1_dir.fo1());
      new_meta1.fanout(hard<2>{});
      val<1> meta_gate1 = do_train & bank1_gate.fo1() & (new_meta1 != old_meta1);
      meta_gate1.fanout(hard<11>{});
      execute_if(meta_gate1, [&]() {
        meta_ctr1.write(val<META_IDX_BITS>{meta_idx_pipe[1][META_PIPE - 1].fo1()},
                        new_meta1, hard<0>{});
      });
    }

    // ================================================================
    // Stage 10: Per-table merged writes
    //
    // For each table: allocation takes priority over training update.
    // pred_ram: alloc writes actual_dir, update writes actual_dir (same data)
    // hyst_ram: alloc writes 0, update writes ta_update_ctr(old, ~wrong)
    // tag_ram + sec_ram: alloc only (new entry setup)
    // u_ram: combined provider update + alloc(0) + uclear + decay
    //
    // Uses per-table wrong signal (this table's prediction vs actual)
    // rather than blanket any_provider_wrong.
    // ================================================================

    // entry_actual_dir: 1-bit actual direction for train_group's branch
    val<PRED_BITS> entry_actual_dir = val<1>{branch_dir[train_group]};
    // table_wrong×NT + pred_ram×NT
    entry_actual_dir.fanout(hard<2 * NT>{});

    train_pc.fanout(hard<NT>{});
    static_loop<NT>([&]<u64 I>() {
      val<1> do_alloc = allocate[I];
      // pred gate + hyst select + hyst gate + alloc gate
      do_alloc.fanout(hard<4>{});

      // Per-table wrong: does THIS table's stored pred disagree with actual?
      // For provider table: stored group == train_group (guaranteed by resolution).
      // For alloc target: new entry is for train_group.
      val<1> table_wrong =
          (val<PRED_BITS>{train_pred[I].fo1()} ^ entry_actual_dir) != hard<0>{};
      // do_pred_update + new_hyst + u_correct
      table_wrong.fanout(hard<3>{});

      // ---- pred_ram write ----
      val<1> do_pred_update = t_m1[I] & t_phw & table_wrong;
      val<1> gate_pred = do_train & (do_alloc | do_pred_update.fo1());
      gate_pred.fanout(hard<11>{}); // ta_rwram B=8: B+3=11
      execute_if(gate_pred, [&]() {
        pred_ram_at<I>().write(val<IDX_BITS[I]>{train_idx[I]},
                               entry_actual_dir, hard<0>{});
      });

      // ---- hyst_ram write ----
      // Alloc: write 0 (fresh entry, neutral hysteresis)
      // Update: saturating counter, increment on correct, decrement on wrong
      auto old_hyst = val<HYST_WIDTH>{train_hyst[I]};
      old_hyst.fanout(hard<2>{}); // ta_update_ctr + (new != old)
      auto new_hyst = ta_update_ctr(old_hyst, ~table_wrong);
      new_hyst.fanout(hard<2>{}); // select + (new != old)
      auto hyst_data = select(do_alloc, val<HYST_WIDTH>{0}, new_hyst);
      val<1> do_hyst_update = t_m1[I] & (new_hyst != old_hyst);
      val<1> gate_hyst = do_train & (do_alloc | do_hyst_update.fo1());
      gate_hyst.fanout(hard<11>{}); // ta_rwram B=8: B+3=11
      execute_if(gate_hyst, [&]() {
        hyst_ram_at<I>().write(val<HYST_IDX_BITS[I]>{train_idx[I]},
                               hyst_data.fo1(), hard<0>{});
      });

      // ---- tag_ram + sec_ram write (alloc only) ----
      // Prepend alloc group_id to htag for full TAG_WIDTH tag.
      val<1> gate_alloc = do_train & do_alloc;
      gate_alloc.fanout(hard<2>{}); // tag_ram + sec_ram
      execute_if(gate_alloc, [&]() {
        val<GROUP_BITS> alloc_gid = val<GROUP_BITS>{u64(train_group)};
        auto full_tag = concat(alloc_gid.fo1(), val<PER_TABLE_HTAG[I]>{train_ctag[I].fo1()});
        tag_ram_at<I>().write(val<IDX_BITS[I]>{train_idx[I]},
                              val<TAG_WIDTH>{full_tag.fo1()});
        sec_ram_at<I>().write(val<IDX_BITS[I]>{train_idx[I]},
                              val<SEC_TAG_BITS>{train_sec_tag});
      });

      // ---- u_ram write ----
      // Combined: provider correct → increment, uclear → decrement,
      // alloc → 0, decay → decrement (probabilistic).
      // U_MISP=UNTOUCHED: no modification on provider wrong.
      val<1> u_correct = t_m1[I] & t_ad & ~table_wrong;

      val<U_WIDTH> old_u = val<U_WIDTH>{train_u[I]};
      // u_inc(2) + u_dec(2) + decay_select(2) + merged_select(1) + merged_changed(1)
      old_u.fanout(hard<10>{});
      // Saturating increment: clamp at maxval
      auto u_inc =
          val<U_WIDTH>{old_u + val<U_WIDTH>{old_u != hard<old_u.maxval>{}}};
      // Saturating decrement: clamp at 0 (used by uclear DECREMENT + decay)
      auto u_dec =
          val<U_WIDTH>{old_u - val<U_WIDTH>{old_u != hard<0>{}}};

      // base_newu: flat OR-of-ANDs (mutually exclusive conditions)
      // u_correct → increment, uclear → decrement, alloc → 0 (vanishes in OR)
      val<U_WIDTH> base_newu = [&]() -> val<U_WIDTH> {
        auto uc_mask =
            val<U_WIDTH>{u_correct.fo1().replicate(hard<U_WIDTH>{}).concat()};
        auto inc_term = val<U_WIDTH>{uc_mask.fo1() & u_inc.fo1()};
        auto ucl_mask = val<U_WIDTH>{
            uclear[I].fo1().replicate(hard<U_WIDTH>{}).concat()};
        return inc_term.fo1() |
               val<U_WIDTH>{ucl_mask.fo1() & u_dec.fo1()};
      }();

      // ---- Probabilistic decay ----
      // TAG_OR_SEC: fire on tag miss OR sec miss (either rejection)
      // DECREMENT: saturating decrement of u on decay
      // FixedThresh<8>: fire when random < 8 (out of 256)
      // LFSR width = 8 for all tables (uniform)
      constexpr u64 LW = 8;
      val<1> tag_missed = ~val<1>{train_tag_hit[I]};
      val<1> sec_missed = ~val<1>{train_sec_hit[I]};
      val<1> decay_miss = tag_missed.fo1() | sec_missed.fo1(); // TAG_OR_SEC

      constexpr u64 THRESH = 8; // FixedDecayThresh<8>
      auto thresh = val<LW>{hard<THRESH>{}};
      val<LW> rng = val<LW>{static_cast<u64>(std::rand())};
      // Fire when: miss AND not allocating here AND random < threshold
      val<1> decay_fire =
          decay_miss.fo1() & ~allocate[I] & (rng.fo1() < thresh.fo1());

      // Decayed u: saturating decrement (DECREMENT policy)
      val<U_WIDTH> decayed_u =
          select(old_u == hard<0>{}, old_u, val<U_WIDTH>{old_u - 1});

      // Merge: base write takes priority, then decay, then hold
      base_newu.fanout(hard<2>{}); // != 0 check + select
      val<1> base_u_write = (base_newu != hard<0>{}) | allocate[I];
      base_u_write.fanout(hard<2>{}); // select + merged_write
      decay_fire.fanout(hard<2>{});   // select + merged_write
      val<U_WIDTH> merged =
          select(base_u_write, base_newu,
                 select(decay_fire, decayed_u.fo1(), old_u));
      val<1> merged_write = base_u_write | decay_fire;
      merged.fanout(hard<2>{}); // merged_changed + write
      val<1> merged_changed = merged != old_u;

      // Silent update elimination: only write when value actually changes
      val<1> gate_u = do_train & (merged_write.fo1() & merged_changed.fo1());
      gate_u.fanout(hard<11>{}); // u_ram ta_rwram B=8: B+3=11
      execute_if(gate_u, [&]() {
        u_ram_at<I>().write(val<IDX_BITS[I]>{train_idx[I]}, merged,
                            hard<0>{});
      });
    });

    // ================================================================
    // Stage 11: Pressure counter updates
    //
    // acc_ctr: accuracy pressure — increment on correct (~mispredict),
    //   decrement on mispredict. Used by AllocPressSkip to throttle
    //   allocation when accuracy is already good.
    // alloc_ctr: allocation pressure — increment when no alloc slot
    //   found (~any_alloc), decrement when allocation succeeds.
    //   High pressure → skip closest candidate (AllocPressSkipTarget).
    //
    // FARALLOC_DIST=0: no far-allocation extra decrement.
    // EPOCH_ENABLE=false: no bulk u-ram reset trigger.
    // ================================================================

    val<1> any_alloc = alloc_target != hard<0>{};
    auto new_acc = ta_update_ctr(val<ACC_WIDTH>{acc_ctr}, ~mispredict);
    auto new_alloc = ta_update_ctr(val<ALLOC_WIDTH>{alloc_ctr}, ~any_alloc.fo1());
    acc_ctr = new_acc.fo1();
    alloc_ctr = new_alloc.fo1();

    // ================================================================
    // Stage 12: History update
    //
    // true_block: was this a real (not speculative) block? Gates
    // fold/GH writes to prevent polluting history with wrong-path data.
    //   ~mispredict         → block was real, commit
    //   last_dir taken      → taken branch ends block at redirect point
    //   line_end            → block filled the fetch line, must commit
    //
    // hist_input: PATH mode — 6-bit next_pc path bits [7:2].
    // Fold updates run in parallel with true_block computation, then
    // select(true_block, new, old) adds only ~10ps mux overhead
    // (avoids execute_if timing bleed where gate chains serially).
    // ================================================================

    true_block = ~mispredict | val<1>{branch_dir[num_branch - 1]} | line_end();
    // NT fold_idx + NT fold_tag apply_update muxes + gh.update + last_condbr_dir
    true_block.fanout(hard<NT * 2 + 2>{});

    // PATH mode: 6-bit path from next_pc[7:2]
    auto hist_input = val<PATHBITS>{block_end_info.next_pc.fo1() >> 2};
    // NT fold_idx + NT fold_tag compute_update + gh.update
    hist_input.fanout(hard<NT * 2 + 1>{});
    gh.template fanout_per_bit<GH_FANOUT>();

    // Per-table: compute new fold values, apply with true_block mux.
    // Both new and old paths resolve in parallel; mux selects based
    // on whether this was a real block or speculative (wrong-path).
    static_loop<NT>([&]<u64 I>() {
      auto &fi = fold_idx_at<I>();
      auto &ft = fold_tag_at<I>();
      auto new_idx = fi.compute_update(gh, hard<HIST_LEN[I]>{}, hist_input);
      auto new_tag = ft.compute_update(gh, hard<HIST_LEN[I]>{}, hist_input);
      fi.apply_update(new_idx.fo1(), true_block);
      ft.apply_update(new_tag.fo1(), true_block);
    });
    gh.update(hist_input, true_block);
    last_condbr_dir = branch_dir[num_branch - 1];
  }
};

using TageAheadHC_IR    = TageAheadHC_IR_impl<256>;
using TageAheadHC_IR_32 = TageAheadHC_IR_impl<32>;
using TageAheadHC_IR_16 = TageAheadHC_IR_impl<16>;
