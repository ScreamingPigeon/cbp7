#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"
#include "custom_common.hpp"

using namespace hcm;

// ============================================================================
// TageAheadHC — Hand-coded TageAhead matching S3_MW4_MC2048
//
// Config: 14 tables, N=7, 6-bit path, 5-bit sec_tag (Xor3SecTagHash5),
//   CTR=1, HYST=2(shared), U=2, Bimodal 8K, Meta 4-bit/2048,
//   AllocPressSkip, SiblingAll, Decay TAG_OR_SEC/DECREMENT/FixedThresh8
//
// Table sizes (T0 bumped from 512→1024):
//   T0-4:  1024 entries (IDX=10, HYST=512)
//   T5-13: 2048 entries (IDX=11, HYST=1024)
// ============================================================================

struct TageAheadHC : predictor {

  // ======== Constants (hardcoded from S3_MW4_MC2048) ========
  static constexpr u64 NT = 14;           // number of tables
  static constexpr u64 N = 7;             // max conditional branches per block
  static constexpr u64 PATHBITS = 6;
  static constexpr u64 SEC_TAG_BITS = 5;
  static constexpr u64 CTR_WIDTH = 1;
  static constexpr u64 HYST_WIDTH = 2;
  static constexpr u64 U_WIDTH = 2;
  static constexpr u64 PRED_BITS = N * CTR_WIDTH; // = 7
  static constexpr u64 TAG_WIDTH = 11;    // UniformTag<11>
  static constexpr u64 MAXHIST = 200;
  static constexpr u64 LINEINST = 256;
  static constexpr u64 LOGLINEINST = 8;   // clog2(256)
  static constexpr u64 FB_CAPACITY = 8192;
  static constexpr u64 FB_IDX_BITS = 13;  // clog2(8192)
  static constexpr u64 META_WIDTH = 4;
  static constexpr u64 META_CAPACITY = 2048;
  static constexpr u64 META_IDX_BITS = 11; // clog2(2048)
  static constexpr u64 META_PIPE = 2;
  static constexpr u64 MATCH_BITS = NT + 1; // = 15
  static constexpr u64 ALLOC_PC_BITS = TAG_WIDTH + 2; // = 13
  static constexpr u64 ACC_WIDTH = 10;
  static constexpr u64 ALLOC_WIDTH = 10;
  static constexpr u64 MAX_IDX_BITS = 11; // clog2(2048)

  // Per-table sizes: T0 bumped to 1024 (was GradedSize<512, 2048>)
  static constexpr std::array<u64, NT> TABLE_SIZE = {
    1024, 1024, 1024, 1024, 1024,
    2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048
  };
  // Per-table history lengths: geometric_hist<14>(8, 200)
  static constexpr std::array<u64, NT> HIST_LEN =
      ta::geometric_hist<NT>(8, 200);
  // Per-table IDX bits
  static constexpr std::array<u64, NT> IDX_BITS = {
    10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11
  };
  // Per-table hyst sizes (shared hyst: TABLE_SIZE/2)
  static constexpr std::array<u64, NT> HYST_SIZE = {
    512, 512, 512, 512, 512,
    1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024
  };

  // Per-bit GH fanout (USE_GSHARE=false, so no fb_fold contribution)
  static constexpr auto GH_FANOUT =
      ta::gh_per_bit_fanout<MAXHIST, NT, HIST_LEN>();

  // ======== Global History ========
  ta_global_history<MAXHIST> gh;

  // ======== Fallback (bimodal, ahead-pipelined) ========
  reg<PRED_BITS> prefetch_fb;
  reg<PRED_BITS> current_fb;
  reg<FB_IDX_BITS> prefetch_fb_idx;
  reg<FB_IDX_BITS> current_fb_idx;

  // Piped PC for allocation tag recomputation
  reg<ALLOC_PC_BITS> prefetch_pc;
  reg<ALLOC_PC_BITS> current_pc;

  // ======== Block tracking ========
  reg<1> true_block = 1;
  reg<1> last_condbr_dir = 1;
  reg<LOGLINEINST> block_entry;

  // Simulation artifacts (free in hardware)
  u64 num_branch = 0;
  u64 block_size = 0;
  arr<reg<1>, N> branch_dir;

  // Secondary tag (precomputed in update_cycle from next_pc)
  reg<SEC_TAG_BITS> curr_sec_tag;

  // ======== Pipeline regs [NT] ========
  reg<TAG_WIDTH> prefetch_tag[NT];
  reg<1> prefetch_tag_hit[NT];
  reg<PRED_BITS> prefetch_pred[NT];
  reg<SEC_TAG_BITS> prefetch_sec[NT];
  reg<MAX_IDX_BITS> prefetch_idx[NT];
  reg<HYST_WIDTH> prefetch_hyst[NT];
  reg<U_WIDTH> prefetch_u[NT];
  reg<TAG_WIDTH> prefetch_ctag[NT];

  reg<TAG_WIDTH> current_tag[NT];
  reg<1> current_tag_hit[NT];
  reg<PRED_BITS> current_pred[NT];
  reg<SEC_TAG_BITS> current_sec[NT];
  reg<MAX_IDX_BITS> current_idx[NT];
  reg<HYST_WIDTH> current_hyst[NT];
  reg<U_WIDTH> current_u[NT];
  reg<TAG_WIDTH> current_ctag[NT];

  // Prediction regs (shared by predict1/2/reuse)
  arr<reg<1>, N> pred;

  // Meta pipeline
  reg<META_WIDTH, i64> meta_pipe[META_PIPE];
  reg<META_IDX_BITS> meta_idx_pipe[META_PIPE];

  // ====================================================================
  // RAMs — per-table clusters matching TATable declaration order:
  //   tag_ram, pred_ram, sec_ram, fold_idx, fold_tag, zone, hyst_ram, u_ram
  //
  // This places fold regs at their table's sec_ram location (co-located
  // with the predict-path RAMs they feed), and hyst/u in a separate
  // zone after each table's predict cluster.
  // ====================================================================

  // ============================================================
  // Per-table zone layout (Exp 6 — best):
  //   tag, fold_idx, fold_tag, pred, sec, zone, hyst, u
  // 14 zones: one per table separating predict/update paths.
  // ============================================================

  // ---- Table 0 (1024 entries, IDX=10) ----
  hcm::ram<val<TAG_WIDTH>, 1024> tag_ram0{"t0_tag"};
  ta_folded_gh<10> fold_idx0;
  ta_folded_gh<TAG_WIDTH> fold_tag0;
  ta_rwram<PRED_BITS, 1024, 2> pred_ram0{"t0_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 1024> sec_ram0{"t0_sec"};
  hcm::zone zone0;
  ta_rwram<HYST_WIDTH, 512, 2> hyst_ram0{"t0_hyst"};
  ta_rwram<U_WIDTH, 1024, 2> u_ram0{"t0_u"};

  // ---- Table 1 (1024 entries, IDX=10) ----
  hcm::ram<val<TAG_WIDTH>, 1024> tag_ram1{"t1_tag"};
  ta_folded_gh<10> fold_idx1;
  ta_folded_gh<TAG_WIDTH> fold_tag1;
  ta_rwram<PRED_BITS, 1024, 2> pred_ram1{"t1_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 1024> sec_ram1{"t1_sec"};
  hcm::zone zone1;
  ta_rwram<HYST_WIDTH, 512, 2> hyst_ram1{"t1_hyst"};
  ta_rwram<U_WIDTH, 1024, 2> u_ram1{"t1_u"};

  // ---- Table 2 (1024 entries, IDX=10) ----
  hcm::ram<val<TAG_WIDTH>, 1024> tag_ram2{"t2_tag"};
  ta_folded_gh<10> fold_idx2;
  ta_folded_gh<TAG_WIDTH> fold_tag2;
  ta_rwram<PRED_BITS, 1024, 2> pred_ram2{"t2_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 1024> sec_ram2{"t2_sec"};
  hcm::zone zone2;
  ta_rwram<HYST_WIDTH, 512, 2> hyst_ram2{"t2_hyst"};
  ta_rwram<U_WIDTH, 1024, 2> u_ram2{"t2_u"};

  // ---- Table 3 (1024 entries, IDX=10) ----
  hcm::ram<val<TAG_WIDTH>, 1024> tag_ram3{"t3_tag"};
  ta_folded_gh<10> fold_idx3;
  ta_folded_gh<TAG_WIDTH> fold_tag3;
  ta_rwram<PRED_BITS, 1024, 2> pred_ram3{"t3_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 1024> sec_ram3{"t3_sec"};
  hcm::zone zone3;
  ta_rwram<HYST_WIDTH, 512, 2> hyst_ram3{"t3_hyst"};
  ta_rwram<U_WIDTH, 1024, 2> u_ram3{"t3_u"};

  // ---- Table 4 (1024 entries, IDX=10) ----
  hcm::ram<val<TAG_WIDTH>, 1024> tag_ram4{"t4_tag"};
  ta_folded_gh<10> fold_idx4;
  ta_folded_gh<TAG_WIDTH> fold_tag4;
  ta_rwram<PRED_BITS, 1024, 2> pred_ram4{"t4_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 1024> sec_ram4{"t4_sec"};
  hcm::zone zone4;
  ta_rwram<HYST_WIDTH, 512, 2> hyst_ram4{"t4_hyst"};
  ta_rwram<U_WIDTH, 1024, 2> u_ram4{"t4_u"};

  // ---- Table 5 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram5{"t5_tag"};
  ta_folded_gh<11> fold_idx5;
  ta_folded_gh<TAG_WIDTH> fold_tag5;
  ta_rwram<PRED_BITS, 2048, 2> pred_ram5{"t5_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram5{"t5_sec"};
  hcm::zone zone5;
  ta_rwram<HYST_WIDTH, 1024, 2> hyst_ram5{"t5_hyst"};
  ta_rwram<U_WIDTH, 2048, 2> u_ram5{"t5_u"};

  // ---- Table 6 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram6{"t6_tag"};
  ta_folded_gh<11> fold_idx6;
  ta_folded_gh<TAG_WIDTH> fold_tag6;
  ta_rwram<PRED_BITS, 2048, 2> pred_ram6{"t6_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram6{"t6_sec"};
  hcm::zone zone6;
  ta_rwram<HYST_WIDTH, 1024, 2> hyst_ram6{"t6_hyst"};
  ta_rwram<U_WIDTH, 2048, 2> u_ram6{"t6_u"};

  // ---- Table 7 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram7{"t7_tag"};
  ta_folded_gh<11> fold_idx7;
  ta_folded_gh<TAG_WIDTH> fold_tag7;
  ta_rwram<PRED_BITS, 2048, 2> pred_ram7{"t7_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram7{"t7_sec"};
  hcm::zone zone7;
  ta_rwram<HYST_WIDTH, 1024, 2> hyst_ram7{"t7_hyst"};
  ta_rwram<U_WIDTH, 2048, 2> u_ram7{"t7_u"};

  // ---- Table 8 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram8{"t8_tag"};
  ta_folded_gh<11> fold_idx8;
  ta_folded_gh<TAG_WIDTH> fold_tag8;
  ta_rwram<PRED_BITS, 2048, 2> pred_ram8{"t8_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram8{"t8_sec"};
  hcm::zone zone8;
  ta_rwram<HYST_WIDTH, 1024, 2> hyst_ram8{"t8_hyst"};
  ta_rwram<U_WIDTH, 2048, 2> u_ram8{"t8_u"};

  // ---- Table 9 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram9{"t9_tag"};
  ta_folded_gh<11> fold_idx9;
  ta_folded_gh<TAG_WIDTH> fold_tag9;
  ta_rwram<PRED_BITS, 2048, 2> pred_ram9{"t9_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram9{"t9_sec"};
  hcm::zone zone9;
  ta_rwram<HYST_WIDTH, 1024, 2> hyst_ram9{"t9_hyst"};
  ta_rwram<U_WIDTH, 2048, 2> u_ram9{"t9_u"};

  // ---- Table 10 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram10{"t10_tag"};
  ta_folded_gh<11> fold_idx10;
  ta_folded_gh<TAG_WIDTH> fold_tag10;
  ta_rwram<PRED_BITS, 2048, 2> pred_ram10{"t10_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram10{"t10_sec"};
  hcm::zone zone10;
  ta_rwram<HYST_WIDTH, 1024, 2> hyst_ram10{"t10_hyst"};
  ta_rwram<U_WIDTH, 2048, 2> u_ram10{"t10_u"};

  // ---- Table 11 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram11{"t11_tag"};
  ta_folded_gh<11> fold_idx11;
  ta_folded_gh<TAG_WIDTH> fold_tag11;
  ta_rwram<PRED_BITS, 2048, 2> pred_ram11{"t11_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram11{"t11_sec"};
  hcm::zone zone11;
  ta_rwram<HYST_WIDTH, 1024, 2> hyst_ram11{"t11_hyst"};
  ta_rwram<U_WIDTH, 2048, 2> u_ram11{"t11_u"};

  // ---- Table 12 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram12{"t12_tag"};
  ta_folded_gh<11> fold_idx12;
  ta_folded_gh<TAG_WIDTH> fold_tag12;
  ta_rwram<PRED_BITS, 2048, 2> pred_ram12{"t12_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram12{"t12_sec"};
  hcm::zone zone12;
  ta_rwram<HYST_WIDTH, 1024, 2> hyst_ram12{"t12_hyst"};
  ta_rwram<U_WIDTH, 2048, 2> u_ram12{"t12_u"};

  // ---- Table 13 (2048 entries, IDX=11) ----
  hcm::ram<val<TAG_WIDTH>, 2048> tag_ram13{"t13_tag"};
  ta_folded_gh<11> fold_idx13;
  ta_folded_gh<TAG_WIDTH> fold_tag13;
  ta_rwram<PRED_BITS, 2048, 2> pred_ram13{"t13_pred"};
  hcm::ram<val<SEC_TAG_BITS>, 2048> sec_ram13{"t13_sec"};
  hcm::zone zone13;
  ta_rwram<HYST_WIDTH, 1024, 2> hyst_ram13{"t13_hyst"};
  ta_rwram<U_WIDTH, 2048, 2> u_ram13{"t13_u"};

  // ---- Fallback RAM ----
  hcm::ram<val<N>, FB_CAPACITY> fb_ctr{"fb"};

  // ---- Meta (update-only) ----
  hcm::zone meta_zone;
  ta_rwram<META_WIDTH, META_CAPACITY, 2> meta_ctr{"meta"};

  // ======== Training regs ========
  reg<MAX_IDX_BITS> train_idx[NT];
  reg<PRED_BITS> train_pred[NT];
  reg<HYST_WIDTH> train_hyst[NT];
  reg<U_WIDTH> train_u[NT];
  reg<PRED_BITS> train_fb;
  reg<FB_IDX_BITS> train_fb_idx;
  reg<ALLOC_PC_BITS> train_pc;
  reg<TAG_WIDTH> train_ctag[NT];
  reg<SEC_TAG_BITS> train_sec_tag;

  reg<MATCH_BITS> train_match1;
  reg<PRED_BITS> train_provider_pred;
  reg<1> train_provider_weak;
  reg<1> train_altdiff;
  reg<1> train_valid = 0;

  // Piped tag/sec hit for decay + sibling
  reg<1> train_tag_hit[NT];
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
    else return name##13;                                                      \
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
    9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10
  };

  // ======== Helpers ========
  val<1> line_end() { return (block_entry + block_size) == hard<LINEINST>{}; }

  // ======== predict1 ========
  val<1> predict1([[maybe_unused]] val<64> inst_pc) {
    // 2 reads per table (>>2, >>4) + fb (>>2) + prefetch_pc (>>2)
    inst_pc.fanout(hard<2 * NT + 2>{});

    // Ahead reads for next block (no true_block gate — see TageAhead comment)
    static_loop<NT>([&]<u64 I>() {
      auto &fi = fold_idx_at<I>();
      auto &ft = fold_tag_at<I>();
      fi.fanout(hard<2>{}); // get() + compute_update
      ft.fanout(hard<2>{}); // get() + compute_update
      auto fold_idx_val = fi.get();
      auto idx = fold_idx_val.fo1() ^ val<IDX_BITS[I]>{inst_pc >> 2};
      idx.fanout(hard<6>{}); // 5 RAM reads + prefetch_idx write
      auto fold_tag_val = ft.get();
      auto computed_tag = fold_tag_val.fo1() ^ val<TAG_WIDTH>{inst_pc >> 4};
      computed_tag.fanout(hard<2>{}); // tag comparison + prefetch_ctag write

      auto stored_tag = tag_ram_at<I>().read(idx);
      stored_tag.fanout(hard<2>{});
      prefetch_tag[I] = stored_tag;
      prefetch_tag_hit[I] =
          val<TAG_WIDTH>{stored_tag} == val<TAG_WIDTH>{computed_tag};
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

    // Fallback ahead read (bimodal: index = PC)
    auto fb_idx = val<FB_IDX_BITS>{inst_pc >> 2};
    fb_idx.fanout(hard<2>{}); // prefetch_fb_idx + fb_ctr.read
    prefetch_fb_idx = fb_idx;
    prefetch_fb = fb_ctr.read(fb_idx);
    prefetch_pc = val<ALLOC_PC_BITS>{inst_pc >> 2};

    // Return precomputed prediction from reg
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
    // Stage 1: Pipeline shift
    //
    // Save current_* → train_* (block A's data for training next cycle),
    // then shift prefetch_* → current_* (block B's ahead reads become
    // active for resolution this cycle).
    //
    // train_sec_hit is computed here by comparing current_sec against
    // curr_sec_tag (the sec-tag hash from the PREVIOUS cycle's next_pc).
    // This is correct: current_sec was prefetched using the predicted PC,
    // and curr_sec_tag was hashed from the actual next_pc of that block.
    // ================================================================

    // curr_sec_tag: NT reads for train_sec_hit comparisons + 1 for saving
    curr_sec_tag.fanout(hard<NT + 1>{});

    static_loop<NT>([&]<u64 I>() {
      train_idx[I] = current_idx[I].fo1();
      train_pred[I] = current_pred[I].fo1();
      train_hyst[I] = current_hyst[I].fo1();
      train_u[I] = current_u[I].fo1();
      train_ctag[I] = current_ctag[I].fo1();
      train_tag_hit[I] = current_tag_hit[I].fo1();
      // Sec-tag match: does stored sec_tag equal the hash from last cycle?
      train_sec_hit[I] = (val<SEC_TAG_BITS>{current_sec[I].fo1()} ==
                          val<SEC_TAG_BITS>{curr_sec_tag});
    });
    train_fb = current_fb.fo1();
    train_fb_idx = current_fb_idx.fo1();
    train_pc = current_pc.fo1();
    train_sec_tag = curr_sec_tag; // save for allocation writes

    // Fanout on prefetch_* regs: each read twice
    //   (1) shift into current_* and (2) resolution chain / weak_mask
    static_loop<NT>([&]<u64 I>() {
      prefetch_tag_hit[I].fanout(hard<2>{}); // shift + full_hits
      prefetch_pred[I].fanout(hard<2>{});    // shift + table_preds
      prefetch_hyst[I].fanout(hard<2>{});    // shift + weak_mask
      prefetch_u[I].fanout(hard<2>{});       // shift + weak_mask
      prefetch_sec[I].fanout(hard<2>{});     // shift + full_hits
    });
    prefetch_fb.fanout(hard<2>{}); // shift + table_preds[NT] in resolve_chain

    // Shift prefetch → current (unconditional)
    static_loop<NT>([&]<u64 I>() {
      current_tag[I] = prefetch_tag[I].fo1();
      current_tag_hit[I] = prefetch_tag_hit[I];
      current_pred[I] = prefetch_pred[I];
      current_sec[I] = prefetch_sec[I];
      current_idx[I] = prefetch_idx[I].fo1();
      current_hyst[I] = prefetch_hyst[I];
      current_u[I] = prefetch_u[I];
      current_ctag[I] = prefetch_ctag[I].fo1();
    });
    current_fb = prefetch_fb;
    current_fb_idx = prefetch_fb_idx.fo1();
    current_pc = prefetch_pc.fo1();

    // ================================================================
    // Stage 2: Sec-tag hash from next_pc
    //
    // Xor3SecTagHash5: val<5>{pc>>2} ^ val<5>{pc>>9} ^ val<5>{pc>>16}
    // This is the BOTTLENECK on the predict path (942ps in debug_print)
    // because next_pc arrives late from the pipeline.
    //
    // sec_tag_now is used for resolution (full_hits comparison) AND
    // stored into curr_sec_tag for next cycle's train_sec_hit.
    // ================================================================

    // next_pc fanout: sec_tag_now(1) + meta_idx(1) + hist path_bits(1)
    block_end_info.next_pc.fanout(hard<3>{});

    auto sec_tag_now = ta::Xor3SecTagHash5::apply<SEC_TAG_BITS>(
        block_end_info.next_pc);
    // Fanout: reg write(1) + NT full_hits comparisons
    sec_tag_now.fanout(hard<NT + 1>{});
    curr_sec_tag = sec_tag_now;
    // Alloc-path: NT sec_ram writes use the OLD sec_tag (saved before overwrite)
    train_sec_tag.fanout(hard<NT>{});
#ifdef DEBUG_PRINT
    std::cerr << "--- sec-tag path ---\n";
    sec_tag_now.print("  sec_tag_now=", "\n", true, std::cerr);
#endif

    // ================================================================
    // Stage 3: Meta pipeline
    //
    // Shift meta_pipe[] down, then read new meta counter value into [0].
    // meta_pipe[META_PIPE-1] (= [1]) holds the delayed meta value used
    // for provider-vs-alt selection in the resolution chain.
    //
    // Shift FIRST, then write [0] — prevents stale-value bug where [1]
    // would read the new RAM value instead of the properly delayed one.
    // ================================================================

    for (u64 i = META_PIPE - 1; i > 0; i--) {
      meta_pipe[i] = meta_pipe[i - 1].fo1();
      meta_idx_pipe[i] = meta_idx_pipe[i - 1].fo1();
    }
    {
      auto meta_idx = val<META_IDX_BITS>{block_end_info.next_pc >> 2};
      meta_idx.fanout(hard<2>{}); // meta_ctr.read + meta_idx_pipe[0] write
      meta_pipe[0] = meta_ctr.read(meta_idx);
      meta_idx_pipe[0] = meta_idx;
    }
    // meta_pipe[1]: read twice — meta_use_alt here + old_meta in training
    meta_pipe[META_PIPE - 1].fanout(hard<2>{});
    // meta_use_alt >= 0 means "trust alt over provider when provider is weak"
    val<1> meta_use_alt =
        val<META_WIDTH, i64>{meta_pipe[META_PIPE - 1]} >= hard<0>{};

    // ================================================================
    // Stage 4: Resolution chain
    //
    // Determines which TAGE table provides the prediction (provider)
    // and which is the alternate (alt). Single-chain with SecTagAll:
    // all tables require sec-tag match for a full hit.
    //
    // Critical path: sec_tag_now → full_hits → match → one_hot →
    //   provider/alt preds → use_alt select → scatter to pred regs
    //
    // The resolution reads prefetch_* directly (bypasses current_*
    // transparent-latch penalty — predict1 wrote prefetch, shift
    // copies to current, reading prefetch skips the second reg hop).
    // ================================================================

    // Save old prediction (block B) before scatter overwrites with B+1
    branch_dir.fanout(hard<3>{}); // true_block + hist_input + actual_dir

    // Precompute per-table weakness: entry is weak when hyst==0 AND u==0.
    // Done from prefetch_* to stay off the resolution critical path.
    val<NT> weak_mask = [&]() {
      arr<val<1>, NT> w = [&](u64 i) -> val<1> {
        return val<1>{val<HYST_WIDTH>{prefetch_hyst[i]} == hard<0>{}} &
               val<1>{val<U_WIDTH>{prefetch_u[i]} == hard<0>{}};
      };
      return w.fo1().concat();
    }();

    // full_hits: tag match AND sec-tag match (SecTagAll = enforce on all tables)
    // sec_tag_now is the timing bottleneck — depends on next_pc (905ps)
    arr<val<1>, NT> full_hits = [&](u64 i) {
      val<1> sec_match =
          val<SEC_TAG_BITS>{prefetch_sec[i]} ==
          val<SEC_TAG_BITS>{sec_tag_now};
      return val<1>{prefetch_tag_hit[i]} & sec_match;
    };

    // ---- resolve_chain lambda ----
    // Takes full_hits → builds match bitmask → finds provider (one_hot) →
    // selects provider/alt predictions → determines use_alt.
    //
    // match layout: bit 0 = fallback (always 1), bits 1..NT = tables
    // match1 = one_hot of match = provider (highest-priority match)
    // match2 = one_hot of (match ^ match1) = alt (second-highest)
    auto resolve_chain = [&](arr<val<1>, NT> fh) {
      // Concatenate fallback (always hit) with table hits
      val<MATCH_BITS> match = concat(val<1>{1}, fh.fo1().concat());
      match.fanout(hard<3>{}); // has_alt check + one_hot + remainder XOR

      // has_alt: any TAGE table matched? (excludes fallback bit 0)
      val<1> ha = val<NT>{match} != hard<0>{};

      // Provider = highest one_hot bit in match
      val<MATCH_BITS> match1 = match.one_hot();
      match1.fanout(hard<4>{}); // make_array + XOR(remainder) + pw + train save

      // Alt = highest one_hot bit in remainder (match with provider removed)
      val<MATCH_BITS> remainder = match ^ match1;
      val<MATCH_BITS> match2 = remainder.fo1().one_hot();

      // Gather per-table predictions (tables 0..NT-1 + fallback at index NT)
      arr<val<PRED_BITS>, NT + 1> table_preds = [&](u64 i) -> val<PRED_BITS> {
        if (i < NT)
          return val<PRED_BITS>{prefetch_pred[i]};
        return val<PRED_BITS>{prefetch_fb};
      };
      table_preds.fanout(hard<2>{}); // provider mask + alt mask

      // Extract provider and alt predictions via one_hot masking
      arr<val<1>, NT + 1> m1_bits = match1.make_array(val<1>{});
      match2.fanout(hard<2>{}); // make_array + return
      arr<val<1>, NT + 1> m2_bits = match2.make_array(val<1>{});

      // Provider pred = OR of (match1[i] & table_preds[i]) across all tables
      arr<val<PRED_BITS>, NT + 1> pmask = [&](u64 i) {
        return m1_bits[i].fo1().replicate(hard<PRED_BITS>{}).concat() &
               table_preds[i];
      };
      arr<val<PRED_BITS>, NT + 1> amask = [&](u64 i) {
        return m2_bits[i].fo1().replicate(hard<PRED_BITS>{}).concat() &
               table_preds[i];
      };
      val<PRED_BITS> pp = pmask.fo1().fold_or();
      val<PRED_BITS> ap = amask.fo1().fold_or();
      pp.fanout(hard<2>{}); // pp^ap + return
      ap.fanout(hard<2>{}); // pp^ap + return

      // Provider weakness: AND match1's table bits with weak_mask.
      // pw=1 means provider is newly-allocated (hyst==0, u==0).
      val<1> pw = (val<NT>{match1} & weak_mask.fo1()) != hard<0>{};
      pw.fanout(hard<2>{}); // use_alt computation + return

      // use_alt: trust alt over provider when provider is weak AND
      // meta says alt is better AND there IS a TAGE provider (ha).
      val<1> ua = pw & meta_use_alt.fo1() & ha.fo1();

      // altdiff: provider and alt predictions disagree on any branch
      auto pp_xor_ap = pp ^ ap;
      val<1> ad = pp_xor_ap.fo1() != hard<0>{};

      return std::tuple{val<MATCH_BITS>{match1},
                        val<MATCH_BITS>{match2},
                        val<PRED_BITS>{pp},
                        val<PRED_BITS>{ap},
                        val<1>{pw},
                        ad.fo1(),
                        ua.fo1()};
    };

    auto [match1, match2, provider_pred, alt_pred, provider_weak, altdiff,
          use_alt] = resolve_chain(full_hits.fo1());

#ifdef DEBUG_PRINT
    std::cerr << "--- resolve_chain output ---\n";
    block_end_info.next_pc.print("  next_pc=", "\n", true, std::cerr);
    match1.print("  match1=", "\n", true, std::cerr);
    provider_pred.print("  provider_pred=", "\n", true, std::cerr);
    alt_pred.print("  alt_pred=", "\n", true, std::cerr);
    provider_weak.print("  provider_weak=", "\n", true, std::cerr);
    use_alt.print("  use_alt=", "\n", true, std::cerr);
#endif
    // ---- Per-branch 1-bit scatter ----
    // Split wide provider/alt preds into per-bit arrays, then select
    // per-branch using 1-bit use_alt. This avoids PRED_BITS-wide
    // replication of use_alt (saves ~45ps ctrl-to-out delay).
    provider_pred.fanout(hard<2>{}); // make_array + train_provider_pred save
    use_alt.fanout(hard<N>{}); // one select per branch
    arr<val<1>, PRED_BITS> pp_bits = provider_pred.make_array(val<1>{});
    arr<val<1>, PRED_BITS> ap_bits = alt_pred.fo1().make_array(val<1>{});
    static_loop<N>([&]<u64 I>() {
      pred[I] = select(use_alt, ap_bits[I].fo1(), pp_bits[I].fo1());
    });

#ifdef DEBUG_PRINT
    for (u64 i = 0; i < N; i++)
      pred[i].print(("  pred[" + std::to_string(i) + "]=").c_str(), "\n",
                    true, std::cerr);
#endif

    // ================================================================
    // Stage 5: Training setup
    //
    // Read OLD piped resolution values (block A) BEFORE overwriting
    // with current resolution (block B). Training uses block A's
    // data while resolution just computed block B's.
    //
    // Fanout declarations for training regs:
    //   train_hyst[I]:    2 reads (t_hyst_weak + old_hyst in update)
    //   train_u[I]:       3 reads (u_zero + old_u base + old_u decay)
    //   train_tag_hit[I]: 2 reads (sibling skip + decay miss)
    //   train_sec_hit[I]: 2 reads (sibling skip + decay miss)
    //   train_idx[I]:     5 reads (pred/hyst/tag/sec/u RAM writes)
    // ================================================================

    static_loop<NT>([&]<u64 I>() {
      train_hyst[I].fanout(hard<2>{});
      train_u[I].fanout(hard<3>{});       // u_zero + base old_u + decay old_u
      train_tag_hit[I].fanout(hard<2>{}); // sibling + decay
      train_sec_hit[I].fanout(hard<2>{}); // sibling + decay
      train_idx[I].fanout(hard<5>{});     // 5 RAM writes
    });

    val<MATCH_BITS> t_match1 = train_match1.fo1();
    t_match1.fanout(hard<2>{}); // make_array + alloc_base
    arr<val<1>, NT + 1> t_m1 = t_match1.make_array(val<1>{});
    // t_m1[I] for I<NT: 4 reads each
    //   t_hyst_weak(1) + do_pred_update(1) + do_hyst_update(1) + base_u_write(1)
    static_loop<NT>([&]<u64 I>() {
      t_m1[I].fanout(hard<4>{});
    });
    // t_m1[NT]: 1 read (fb_gate). FB_RECONCILE=false, no extra fanout needed.

    val<PRED_BITS> t_pp = train_provider_pred.fo1();
    val<1> t_pw = train_provider_weak.fo1(); // newly-alloc weakness → meta gate
    val<1> t_ad = train_altdiff.fo1();       // provider/alt disagree → meta gate

    // Hyst-only weakness: provider has hyst==0 (regardless of u).
    // Gates pred/hyst counter updates — a useful entry (u>0) with
    // weak hyst that's wrong should still get its pred flipped.
    arr<val<1>, NT + 1> t_hyst_weak_arr = [&](u64 i) -> val<1> {
      if (i < NT)
        return t_m1[i] & val<1>{val<HYST_WIDTH>{train_hyst[i]} == hard<0>{}};
      return val<1>{0};
    };
    val<1> t_phw = t_hyst_weak_arr.fo1().fold_or();

    // Save current resolution → train regs (for NEXT cycle's training)
    train_match1 = match1.fo1();
    train_provider_pred = provider_pred;
    train_provider_weak = provider_weak.fo1();
    train_altdiff = altdiff.fo1();

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
    // extra_cycle(1) + alloc_trigger(1) + fb_gate(1) + acc_ctr(1) + true_block(1)
    mispredict.fanout(hard<5>{});
    need_extra_cycle(mispredict);

    // do_train gates all RAM writes: 4 per table + fb + meta
    do_train.fanout(hard<4 * NT + 2>{});

    // Pack all branch directions into a single PRED_BITS-wide value
    val<PRED_BITS> actual_dir = arr<val<1>, N>{[&](u64 i) -> val<1> {
                                  return val<1>{branch_dir[i]};
                                }}.concat();
    // any_prov_wrong(1) + table_wrong×NT + pred_ram×NT + fb_changed(1) + fb_ctr(1)
    actual_dir.fanout(hard<3 + 2 * NT>{});

    val<1> any_provider_wrong = (t_pp.fo1() ^ actual_dir) != hard<0>{};
    t_pw.fanout(hard<2>{});    // meta_gate + meta_update_dir
    t_phw.fanout(hard<NT>{});  // do_pred_update per table
    // u_correct×NT + meta_gate (U_MISP=UNTOUCHED, no u_wrong reads)
    t_ad.fanout(hard<NT + 1>{});

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

    // 7d. Sibling skip (SiblingPolicy::ALL, FLOOR=0)
    // Skip entries where primary tag matches but sec-tag doesn't (siblings).
    // Promotes allocation to the next higher table.
    val<NT> candallocmask = [&]() {
      val<NT> base = postmask & notumask.fo1();
      arr<val<1>, NT> not_sibling = [&](u64 i) -> val<1> {
        // Both sibling and decay read these (fanout(2) declared above)
        auto th = val<1>{train_tag_hit[i]};
        auto sh = val<1>{train_sec_hit[i]};
        // not_sibling = !(tag_hit & !sec_hit) = !tag_hit | sec_hit
        return ~(th.fo1() & ~sh.fo1());
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

    val<1> fb_changed = actual_dir != val<PRED_BITS>{train_fb.fo1()};
    val<1> fb_gate = do_train & t_m1[NT].fo1() & mispredict & fb_changed.fo1();
    fb_gate.fanout(hard<5>{}); // fb_ctr ta_rwram: B=2 banks + write_bank + localaddr + data
    execute_if(fb_gate, [&]() {
      fb_ctr.write(val<FB_IDX_BITS>{train_fb_idx.fo1()}, actual_dir);
    });

    // ================================================================
    // Stage 9: Meta counter update
    //
    // Meta tracks whether provider or alt is more accurate.
    // Updated when provider is newly-allocated (weak) AND provider/alt
    // disagree. Increment when provider was wrong (trust alt more),
    // decrement when provider was right (trust provider more).
    // ================================================================

    auto old_meta = val<META_WIDTH, i64>{meta_pipe[META_PIPE - 1]};
    old_meta.fanout(hard<2>{}); // ta_update_ctr + (new != old)
    auto new_meta = ta_update_ctr(old_meta, any_provider_wrong.fo1());
    new_meta.fanout(hard<2>{}); // (new != old) + meta_ctr.write
    // Gate: train valid AND provider weak AND provider/alt disagree AND value changed
    val<1> meta_gate = do_train & t_pw & t_ad & (new_meta != old_meta);
    meta_gate.fanout(hard<5>{}); // meta_ctr ta_rwram: B=2 + 3
    execute_if(meta_gate, [&]() {
      meta_ctr.write(val<META_IDX_BITS>{meta_idx_pipe[META_PIPE - 1].fo1()},
                     new_meta, hard<0>{});
    });

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

    train_pc.fanout(hard<NT>{});
    static_loop<NT>([&]<u64 I>() {
      val<1> do_alloc = allocate[I];
      // pred gate + hyst select + hyst gate + alloc gate
      do_alloc.fanout(hard<4>{});

      // Per-table wrong: does THIS table's stored pred disagree with actual?
      val<1> table_wrong =
          (val<PRED_BITS>{train_pred[I].fo1()} ^ actual_dir) != hard<0>{};
      // do_pred_update + new_hyst + u_correct
      table_wrong.fanout(hard<3>{});

      // ---- pred_ram write ----
      // Gated on hyst-only weakness (t_phw): a useful entry (u>0) with
      // weak hyst that's wrong should still get its prediction flipped.
      val<1> do_pred_update = t_m1[I] & t_phw & table_wrong;
      val<1> gate_pred = do_train & (do_alloc | do_pred_update.fo1());
      gate_pred.fanout(hard<5>{}); // ta_rwram B=2: B+3=5
      execute_if(gate_pred, [&]() {
        pred_ram_at<I>().write(val<IDX_BITS[I]>{train_idx[I]},
                               actual_dir, hard<0>{});
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
      gate_hyst.fanout(hard<5>{});
      execute_if(gate_hyst, [&]() {
        hyst_ram_at<I>().write(val<HYST_IDX_BITS[I]>{train_idx[I]},
                               hyst_data.fo1(), hard<0>{});
      });

      // ---- tag_ram + sec_ram write (alloc only) ----
      // New entry: write computed tag from predict1 time (piped through
      // train_ctag) and the saved sec_tag hash (train_sec_tag).
      val<1> gate_alloc = do_train & do_alloc;
      gate_alloc.fanout(hard<2>{}); // tag_ram + sec_ram
      execute_if(gate_alloc, [&]() {
        tag_ram_at<I>().write(val<IDX_BITS[I]>{train_idx[I]},
                              val<TAG_WIDTH>{train_ctag[I].fo1()});
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
      old_u.fanout(hard<8>{});
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
      merged.fanout(hard<2>{}); // merged_changed + return
      val<1> merged_changed = merged != old_u;

      // Silent update elimination: only write when value actually changes
      val<1> gate_u = do_train & (merged_write.fo1() & merged_changed.fo1());
      gate_u.fanout(hard<5>{}); // u_ram ta_rwram B=2: B+3=5
      execute_if(gate_u, [&]() {
        u_ram_at<I>().write(val<IDX_BITS[I]>{train_idx[I]}, merged.fo1(),
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
    auto new_alloc = ta_update_ctr(val<ALLOC_WIDTH>{alloc_ctr}, ~any_alloc);
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
    auto hist_input = val<PATHBITS>{block_end_info.next_pc >> 2};
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
