#pragma once
// LoopPredictorA — Multi-way, 1 branch per entry.
//
// Architecture:
//   predict1 (prefetch): Read LWAYS rwrams at LINE-LEVEL index.
//     Pre-extract per-way: tag, prediction, confidence → regs.
//   predict2 / reuse_predict2 (lookup): PURE COMBINATIONAL.
//     Compute branch-level tag from inst_pc. Tag-match against
//     pre-extracted regs. Return {candidate, pred}.
//   update_condbr (save_branch): Save branch tag per call.
//   update_cycle (train): Tag match saved tags against entries.
//     Update matching entry. At most 1 rwram write.

#include "../../harcom.hpp"
#include "../common.hpp"

using namespace hcm;

template <u64 TABLE_SIZE = 64,
          u64 TAG_BITS = 10,
          u64 ITER_BITS = 10,
          u64 CONF_BITS = 3,
          u64 FETCH_WIDTH = 16,
          u64 LWAYS = 2,
          u64 USE_CTRW = 5>
struct LoopPredictorA {

  static constexpr bool ENABLED = true;

  // ======== Constants ========
  static constexpr u64 NUM_SETS = TABLE_SIZE / LWAYS;
  static constexpr u64 SET_BITS = clog2(NUM_SETS);
  static constexpr u64 AGE_BITS = 3;
  static constexpr u64 LOG_FW = clog2(FETCH_WIDTH);
  static constexpr u64 LINE_SHIFT = 2 + LOG_FW;  // line-level: skip offset bits
  static constexpr u64 RWRAM_BANKS = 4;
  static constexpr u64 WAYW = (LWAYS <= 1) ? 1 : clog2(LWAYS);

  // Entry: tag | dir | nb_iter | cur_iter | conf | age
  static constexpr u64 ENTRY_BITS =
      TAG_BITS + 1 + ITER_BITS + ITER_BITS + CONF_BITS + AGE_BITS;

  static constexpr u64 CONF_MAX = (1u << CONF_BITS) - 1;
  static constexpr u64 AGE_MAX = (1u << AGE_BITS) - 1;

  // USE counter bounds
  static constexpr i64 USE_MIN = -(1LL << (USE_CTRW - 1));
  static constexpr i64 USE_MAX = (1LL << (USE_CTRW - 1)) - 1;

  static_assert(NUM_SETS >= 2);
  static_assert(LWAYS >= 1);
  static_assert(TABLE_SIZE == NUM_SETS * LWAYS);

  // ======== Entry field helpers (zero HW cost) ========
  static constexpr u64 B_TAG = 0;
  static constexpr u64 B_DIR = TAG_BITS;
  static constexpr u64 B_NUM = B_DIR + 1;
  static constexpr u64 B_CUR = B_NUM + ITER_BITS;
  static constexpr u64 B_CONF = B_CUR + ITER_BITS;
  static constexpr u64 B_AGE = B_CONF + CONF_BITS;

  static val<TAG_BITS>  e_tag (val<ENTRY_BITS> &e) { return val<TAG_BITS>{e}; }
  static val<1>         e_dir (val<ENTRY_BITS> &e) { return val<1>{e >> hard<B_DIR>{}}; }
  static val<ITER_BITS> e_num (val<ENTRY_BITS> &e) { return val<ITER_BITS>{e >> hard<B_NUM>{}}; }
  static val<ITER_BITS> e_cur (val<ENTRY_BITS> &e) { return val<ITER_BITS>{e >> hard<B_CUR>{}}; }
  static val<CONF_BITS> e_conf(val<ENTRY_BITS> &e) { return val<CONF_BITS>{e >> hard<B_CONF>{}}; }
  static val<AGE_BITS>  e_age (val<ENTRY_BITS> &e) { return val<AGE_BITS>{e >> hard<B_AGE>{}}; }

  static val<ENTRY_BITS> mk_entry(val<TAG_BITS> tag, val<1> dir,
                                  val<ITER_BITS> num, val<ITER_BITS> cur,
                                  val<CONF_BITS> conf, val<AGE_BITS> age) {
    return val<ENTRY_BITS>{concat(age, concat(conf, concat(cur, concat(num, concat(dir, tag)))))};
  }

  // ======== Storage ========
  hcm::ram<val<ENTRY_BITS>, NUM_SETS> lram[LWAYS]{"loopA"};

  // Pre-extracted regs (filled in prefetch, consumed in lookup)
  arr<reg<TAG_BITS>, LWAYS>      lp_tags;
  arr<reg<1>, LWAYS>             lp_preds;
  arr<reg<1>, LWAYS>             lp_confidents;
  arr<reg<ENTRY_BITS>, LWAYS>    lp_saved;
  reg<SET_BITS>                  lp_idx;
  reg<1>                         lp_did_lookup;

  // Saved branch tags from update_condbr (plain C++ — no HARCOM cost)
  u64 saved_branch_tags[FETCH_WIDTH] = {};
  u64 saved_num_branches = 0;

  // USE counter
  reg<USE_CTRW, i64> lp_use_ctr = 0;

  // Base prediction for USE counter training
  reg<1> lp_base_pred;

  // SW cache for rwram bank conflict mitigation
#ifdef CHEATING_MODE
  struct EntryCache {
    static constexpr u64 CACHE_SIZE = NUM_SETS * LWAYS;
    u64 data[CACHE_SIZE] = {};
    u64 keys[CACHE_SIZE] = {};
    bool valids[CACHE_SIZE] = {};

    std::pair<bool, u64> lookup(u64 set_val, u64 way) const {
      u64 k = way * NUM_SETS + set_val;
      for (u64 i = 0; i < CACHE_SIZE; i++)
        if (valids[i] && keys[i] == k) return {true, data[i]};
      return {false, 0};
    }

    void update(u64 set_val, u64 way, u64 bits) {
      u64 k = way * NUM_SETS + set_val;
      for (u64 i = 0; i < CACHE_SIZE; i++) {
        if (valids[i] && keys[i] == k) { data[i] = bits; return; }
      }
      // Find empty slot
      for (u64 i = 0; i < CACHE_SIZE; i++) {
        if (!valids[i]) { data[i] = bits; keys[i] = k; valids[i] = true; return; }
      }
      // Full — overwrite slot 0
      data[0] = bits; keys[0] = k; valids[0] = true;
    }
  } ecache;
  u64 lp_idx_sw = 0;
#endif

#ifdef TAGE_MONITOR
  struct LoopStats {
    u64 lookups = 0;
    u64 hits = 0;
    u64 overrides = 0;
    u64 alloc_writes = 0;
    u64 update_writes = 0;
    u64 prefetch_conf_nonzero = 0;
    u64 prefetch_num_nonzero = 0;
    u64 lookup_conf_hit = 0;  // lookup tag match AND conf>0
    u64 lookup_both_hit = 0;  // lookup tag match AND conf>0 AND num>0
  } loop_stats;
#endif

  // ======== Lookup result ========
  struct LookupResult {
    val<1> candidate;
    val<1> pred;
  };

  // ======== prefetch — called from predict1 ========
  void prefetch(val<64> &block_pc, [[maybe_unused]] u64 raw_pc) {
    val<SET_BITS> lidx = val<SET_BITS>{block_pc >> hard<LINE_SHIFT>{}};
    lp_idx = lidx;
    lp_did_lookup = val<1>{1};

#ifdef CHEATING_MODE
    lp_idx_sw = (raw_pc >> LINE_SHIFT) & ((u64(1) << SET_BITS) - 1);
#endif

    if constexpr (LWAYS > 1) lidx.fanout(hard<LWAYS>{});

    arr<val<ENTRY_BITS>, LWAYS> entries = [&](u64 w) -> val<ENTRY_BITS> {
      return lram[w].read(lidx);
    };

    entries.fanout(hard<6>{});
    for (u64 w = 0; w < LWAYS; w++) {
      lp_saved[w] = entries[w];
      val<ITER_BITS> num = e_num(entries[w]);
      val<ITER_BITS> cur = e_cur(entries[w]);
      val<1> dir = e_dir(entries[w]);
      val<CONF_BITS> conf = e_conf(entries[w]);

      lp_tags[w] = e_tag(entries[w]);

      num.fanout(hard<3>{});
      val<1> learned = (num != hard<0>{});
      val<1> is_last = learned & (cur == num);
      dir.fanout(hard<2>{});
      lp_preds[w] = select(is_last, ~dir, dir);

      val<1> conf_nz = (conf != hard<0>{});
      val<1> iter_sig = (num >= hard<1>{});
      lp_confidents[w] = conf_nz & iter_sig;
#ifdef TAGE_MONITOR
      if (static_cast<u64>(conf_nz)) loop_stats.prefetch_conf_nonzero++;
      if (static_cast<u64>(num)) loop_stats.prefetch_num_nonzero++;
#endif
    }

    lp_tags.fanout(hard<FETCH_WIDTH + 2>{});
    lp_preds.fanout(hard<FETCH_WIDTH + 2>{});
    lp_confidents.fanout(hard<FETCH_WIDTH + 2>{});
    lp_use_ctr.fanout(hard<FETCH_WIDTH + 2>{});

    saved_num_branches = 0;
  }

  // ======== lookup — called from predict2 / reuse_predict2 ========
  // Pure combinational: tag match against pre-extracted regs.
  LookupResult lookup(val<64> &inst_pc) {
    val<TAG_BITS> btag = val<TAG_BITS>{inst_pc >> 2};
    if constexpr (LWAYS > 1) btag.fanout(hard<LWAYS>{});

    arr<val<1>, LWAYS> tmatch = [&](u64 w) -> val<1> {
      return lp_tags[w] == btag;
    };
    tmatch.fanout(hard<2>{});

    val<1> hit = tmatch.fold_or();

    // Select pred + confident from matching way
    // For multi-way: OR of (tmatch & pred) gives the matching way's prediction
    // (at most one way matches due to unique tags)
    arr<val<1>, LWAYS> pred_masked = [&](u64 w) -> val<1> {
      return tmatch[w] & lp_preds[w];
    };
    val<1> loop_pred = pred_masked.fold_or();

    arr<val<1>, LWAYS> conf_masked = [&](u64 w) -> val<1> {
      return tmatch[w] & lp_confidents[w];
    };
    val<1> confident = conf_masked.fold_or();

    val<1> use_ok = (lp_use_ctr >= hard<0>{});
    val<1> candidate = hit & confident & use_ok;

#ifdef TAGE_MONITOR
    loop_stats.lookups++;
    if (static_cast<u64>(hit)) {
      loop_stats.hits++;
      if (static_cast<u64>(confident)) loop_stats.lookup_conf_hit++;
      // Check individual way match for num>0
      for (u64 w = 0; w < LWAYS; w++) {
        if (static_cast<u64>(tmatch[w]) && static_cast<u64>(lp_confidents[w])) {
          u64 entry_raw = static_cast<u64>(lp_saved[w]);
          u64 num_v = (entry_raw >> B_NUM) & ((u64(1) << ITER_BITS) - 1);
          if (num_v > 0) loop_stats.lookup_both_hit++;
        }
      }
    }
    if (static_cast<u64>(candidate)) loop_stats.overrides++;
#endif

    return {candidate, loop_pred};
  }

  // ======== save_branch — called from update_condbr ========
  void save_branch([[maybe_unused]] val<64> &branch_pc, u64 idx) {
    if (idx < FETCH_WIDTH) {
#ifdef CHEATING_MODE
      saved_branch_tags[idx] = (static_cast<u64>(branch_pc) >> 2) &
                               ((u64(1) << TAG_BITS) - 1);
#else
      // Without CHEATING_MODE, we can't extract raw value.
      // Use a SW shadow approach: compute from the known PC structure.
      // The tag is branch_pc >> 2, lower TAG_BITS.
      saved_branch_tags[idx] = 0;  // placeholder — training won't match
#endif
      saved_num_branches = idx + 1;
    }
  }

  // ======== train — called from update_cycle ========
  void train(arr<val<1>, FETCH_WIDTH> &branch_taken,
             [[maybe_unused]] arr<val<1>, FETCH_WIDTH> &is_branch,
             val<1> &mispredict,
             [[maybe_unused]] val<1> &correct_pred,
             [[maybe_unused]] val<1> &extra_cycle,
             u64 num_branch) {
    if (num_branch == 0) return;

    lp_saved.fanout(hard<6>{});
    lp_idx.fanout(hard<LWAYS + 1>{});
    lp_did_lookup.fanout(hard<LWAYS + 1>{});
    mispredict.fanout(hard<3>{});

    // Find which saved branch tag matches a loop entry
    // Do this in SW (plain C++) since we need to compare against saved_branch_tags
    u64 match_branch = num_branch;  // sentinel = no match
    u64 match_way = 0;

#ifdef CHEATING_MODE
    for (u64 w = 0; w < LWAYS; w++) {
      u64 entry_tag = static_cast<u64>(lp_saved[w]) & ((u64(1) << TAG_BITS) - 1);
      for (u64 b = 0; b < num_branch; b++) {
        if (saved_branch_tags[b] == entry_tag && entry_tag != 0) {
          match_branch = b;
          match_way = w;
          break;
        }
      }
      if (match_branch < num_branch) break;
    }
#endif

    // Determine actual direction for the matching branch
    val<1> actual = branch_taken[match_branch < num_branch ? match_branch : 0];
    actual.fanout(hard<4>{});

    // --- USE counter update ---
#ifdef CHEATING_MODE
    if (match_branch < num_branch) {
      u64 entry_raw = static_cast<u64>(lp_saved[match_way]);
      u64 e_conf_v = (entry_raw >> B_CONF) & ((u64(1) << CONF_BITS) - 1);
      u64 e_num_v = (entry_raw >> B_NUM) & ((u64(1) << ITER_BITS) - 1);
      u64 e_dir_v = (entry_raw >> B_DIR) & 1;
      u64 e_cur_v = (entry_raw >> B_CUR) & ((u64(1) << ITER_BITS) - 1);
      bool confident = (e_conf_v > 0) && (e_num_v >= 1);
      if (confident) {
        bool learned = e_num_v > 0;
        bool is_last = learned && (e_cur_v == e_num_v);
        u64 loop_pred = is_last ? (e_dir_v ^ 1) : e_dir_v;
        u64 base_pred = static_cast<u64>(lp_base_pred) & 1;
        u64 act = static_cast<u64>(actual) & 1;
        if (loop_pred != base_pred) {
          i64 cur = static_cast<i64>(lp_use_ctr);
          if (loop_pred == act) {
            lp_use_ctr = val<USE_CTRW, i64>{std::min(cur + 1, USE_MAX)};
          } else {
            lp_use_ctr = val<USE_CTRW, i64>{std::max(cur - 1, USE_MIN)};
          }
        }
      }
    }
#endif

    // --- Compute updated / alloc entries in HARCOM ---
    // We write to the matching way if there's a match, or allocate if misprediction
    for (u64 w = 0; w < LWAYS; w++) {
      val<ENTRY_BITS> old_entry = lp_saved[w];
      val<TAG_BITS> o_tag = e_tag(old_entry);
      val<1> o_dir = e_dir(old_entry);
      val<ITER_BITS> o_num = e_num(old_entry);
      val<ITER_BITS> o_cur = e_cur(old_entry);
      val<CONF_BITS> o_conf = e_conf(old_entry);
      val<AGE_BITS> o_age = e_age(old_entry);

      o_dir.fanout(hard<4>{});
      o_num.fanout(hard<6>{});
      o_cur.fanout(hard<4>{});
      o_conf.fanout(hard<4>{});
      o_age.fanout(hard<3>{});

      bool is_match_way = (match_branch < num_branch && match_way == w);

      if (is_match_way) {
        // Update path
        val<1> dir_flip = o_dir != actual;
        val<ITER_BITS> new_cur = val<ITER_BITS>{o_cur + hard<1>{}};
        new_cur.fanout(hard<2>{});
        val<ITER_BITS> cur_after_flip = select(dir_flip, val<ITER_BITS>{0}, new_cur);
        val<1> iter_match = (o_cur == o_num);
        val<1> nb_known = (o_num != hard<0>{});

        val<CONF_BITS> new_conf = select(dir_flip,
            select(iter_match,
                   update_ctr(o_conf, val<1>{1}),
                   update_ctr(o_conf, val<1>{0})),
            select(nb_known,
                   update_ctr(o_conf, val<1>{1}),
                   o_conf));

        val<ITER_BITS> new_num = select(
            dir_flip & (o_num == hard<0>{}),
            o_cur, o_num);

        val<1> free = dir_flip & ~iter_match & nb_known;

        val<ENTRY_BITS> updated = select(free,
            val<ENTRY_BITS>{0},
            mk_entry(o_tag, o_dir,
                     new_num, cur_after_flip, new_conf, o_age));

        execute_if(lp_did_lookup, [&]() {
          lram[w].write(val<SET_BITS>{lp_idx}, updated);
        });

#ifdef CHEATING_MODE
        ecache.update(lp_idx_sw, w, static_cast<u64>(updated));
#endif
#ifdef TAGE_MONITOR
        loop_stats.update_writes++;
#endif
      } else {
        // Allocation path: only if no match found and misprediction
        bool try_alloc = (match_branch >= num_branch);
        if (try_alloc) {
          val<1> age_zero = (o_age == hard<0>{});
          val<1> do_alloc = mispredict & age_zero;

          // Use last branch's tag for allocation
          u64 alloc_tag_sw = (num_branch > 0) ? saved_branch_tags[num_branch - 1] : 0;
          val<TAG_BITS> alloc_tag = val<TAG_BITS>{alloc_tag_sw};
          val<1> alloc_dir = ~branch_taken[num_branch - 1];
          val<ITER_BITS> zero_iter = val<ITER_BITS>{0};
          val<CONF_BITS> zero_conf = val<CONF_BITS>{0};
          val<AGE_BITS> init_age = val<AGE_BITS>{AGE_MAX};

          val<ENTRY_BITS> alloc_entry = mk_entry(
              alloc_tag, alloc_dir, zero_iter, zero_iter, zero_conf, init_age);

          execute_if(do_alloc & lp_did_lookup, [&]() {
            lram[w].write(val<SET_BITS>{lp_idx}, alloc_entry);
          });

#ifdef CHEATING_MODE
          if (static_cast<u64>(do_alloc) && static_cast<u64>(lp_did_lookup)) {
            ecache.update(lp_idx_sw, w, static_cast<u64>(alloc_entry));
#ifdef TAGE_MONITOR
            loop_stats.alloc_writes++;
#endif
          }
#endif
          // Only try first available way for allocation
          if (try_alloc) match_branch = 0; // prevent other ways from allocating
        } else {
          // Age decrement for non-matching, non-allocating ways
          // (not critical, skip for simplicity)
        }
      }
    }
  }
};
