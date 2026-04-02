#pragma once
// LoopPredictorC — Loop predictor overrider for TAGE.
//
// Key improvements over LoopPredictorA:
//
//   1. No CHEATING_MODE dependency.
//      First branch tag saved into a HARCOM reg on the first lookup(),
//      not a plain C++ array. Training is always functional.
//
//   2. USE/WITHLOOP counter updated fully in HARCOM.
//      execute_if(lp_did_lookup & disagrees, ...) in train_compute() — no #ifdef.
//
//   3. Paper-faithful defaults (L-TAGE CBP-2):
//      tag=14, iter=14, conf=2 bits (full-saturation threshold), age=8,
//      4-way, 64 sets (256 entries total), 7-bit WITHLOOP counter.
//
//   4. Confidence threshold: CONF_MAX (full saturation) by default.
//      Matches "3 consecutive consistent trips" from the L-TAGE paper.
//      Optional long-loop product gate (from gshareN_loop).
//
//   5. Age managed per L-TAGE paper:
//      - Initialized to ALLOC_AGE_INIT (255) on allocation.
//      - Decremented on non-victim ways when no free slot exists (is_age_dec).
//      - Incremented when confidently used as provider (was_provider).
//      - Reset to 0 on invalidation.
//
//   6. Priority MUX way selection (first_true_way + select_bit/select_entry)
//      instead of OR-masking fold — avoids ambiguity on multiple hits.
//
//   7. Full INVALIDATE_ON_* policy knobs (from gshareN_loop).
//
//   8. Uses ram<> per TageOverrider contract (read p1, write in train_write()).
//      Write data is staged explicitly to match Tage's overrider interface.
//
// Latency:
//   prefetch() — RAM read + pre-extraction into regs. All RAM access done here.
//   lookup()   — Pure combinational. Tag match + MUX. One first-call reg write
//                (lp_base_pred, guarded by lp_base_pred_saved bool).
//   save_branch() — no-op; update_condbr has already consumed branch_pc via fo1().
//   train_compute() — Combinational compute + stage writes into regs.
//   train_write()   — Commits staged writes to lram (called in extra cycle).

#include "../../harcom.hpp"
#include "../common.hpp"
#include "TageOverrider.hpp"

using namespace hcm;

template<
    // --- Table geometry ---
    u64 LOOP_LOGSETS = 6,      // sets = 2^LOOP_LOGSETS; x LOOP_WAYS = total entries
    u64 LOOP_TAGW    = 14,     // tag width  (L-TAGE paper: 14 bits)
    u64 LOOP_ITERW   = 14,     // iter width (L-TAGE paper: 14 bits)
    u64 LOOP_CONFW   = 2,      // conf width (L-TAGE paper:  2 bits)
    u64 LOOP_AGEW    = 8,      // age width  (L-TAGE paper:  8 bits, max 255)
    u64 LOOP_WAYS    = 4,      // associativity (L-TAGE paper: 4-way)
    u64 FETCH_WIDTH  = 16,

    // --- Confidence gate ---
    // LOOP_CONF_THRESH: default = CONF_MAX = full saturation.
    // Matches L-TAGE "3 consecutive consistent trips" with a 2-bit counter.
    u64  LOOP_CONF_THRESH   = ((1ULL << LOOP_CONFW) - 1),
    bool USE_LONG_LOOP_CONF = false,
    u64  LONG_LOOP_THRESH   = 128,

    // --- WITHLOOP counter (L-TAGE paper: 7-bit, override threshold = 0) ---
    u64  USE_CTRW         = 7,
    i64  USE_CTR_INIT     = 0,
    i64  USE_CTR_THRESH   = 0,
    i64  USE_CTR_INC_STEP = 1,
    i64  USE_CTR_DEC_STEP = 1,

    // --- Allocation / update policy ---
    u64  MIN_TRIP                 = 3,
    u64  ALLOC_AGE_INIT           = ((1ULL << LOOP_AGEW) - 1),  // paper: 255
    u64  ALLOC_CONF_INIT          = 0,
    u64  CONF_BUMP_STEP           = 1,
    bool INVALIDATE_ON_WRONG_CONF = true,
    bool INVALIDATE_ON_EXCEEDED   = true,
    bool INVALIDATE_ON_BAD_EXIT   = true,
    bool INVALIDATE_ON_TOO_SHORT  = true
>
struct LoopPredictorC {

    static constexpr bool ENABLED            = true;
    static constexpr bool IS_LOOP            = true;
    static constexpr bool IS_SC              = false;
    static constexpr u64  EXTRA_CYCLE_BUDGET = 0;
    static constexpr bool REUSE_LOOKUP       = false;

    // =========================================================
    // Derived constants
    // =========================================================
    static constexpr u64 LWAYS   = LOOP_WAYS;
    static_assert(LWAYS == 1 || LWAYS == 2 || LWAYS == 4);

    static constexpr u64 LSETS   = 1ULL << LOOP_LOGSETS;
    static constexpr u64 LLSETS  = LOOP_LOGSETS;

    static constexpr u64 LOG_FW     = clog2(FETCH_WIDTH);
    static constexpr u64 LINE_SHIFT = 2 + LOG_FW;

    static constexpr u64 CONF_MAX = (1ULL << LOOP_CONFW) - 1;
    static constexpr u64 AGE_MAX  = (1ULL << LOOP_AGEW)  - 1;
    static constexpr u64 WAYW     = (LWAYS <= 1) ? 1 : clog2(LWAYS);

    static_assert(LOOP_CONF_THRESH <= CONF_MAX);
    static_assert(ALLOC_CONF_INIT  <= CONF_MAX);
    static_assert(ALLOC_AGE_INIT   <= AGE_MAX);
    static_assert(CONF_BUMP_STEP   >= 1);
    static_assert(MIN_TRIP         >= 1);
    static_assert(USE_CTRW         >= 2);

    static constexpr i64 USE_MIN = -(1LL << (USE_CTRW - 1));
    static constexpr i64 USE_MAX =  (1LL << (USE_CTRW - 1)) - 1;
    static_assert(USE_CTR_INIT    >= USE_MIN && USE_CTR_INIT    <= USE_MAX);
    static_assert(USE_CTR_THRESH  >= USE_MIN && USE_CTR_THRESH  <= USE_MAX);
    static_assert(USE_CTR_INC_STEP >= 1);
    static_assert(USE_CTR_DEC_STEP >= 1);

    // =========================================================
    // Entry layout: [tag | dir | num | cur | conf | age]
    // =========================================================
    static constexpr u64 B_TAG  = 0;
    static constexpr u64 B_DIR  = LOOP_TAGW;
    static constexpr u64 B_NUM  = B_DIR  + 1;
    static constexpr u64 B_CUR  = B_NUM  + LOOP_ITERW;
    static constexpr u64 B_CONF = B_CUR  + LOOP_ITERW;
    static constexpr u64 B_AGE  = B_CONF + LOOP_CONFW;
    static constexpr u64 LENTRY = B_AGE  + LOOP_AGEW;

    // =========================================================
    // Storage — ram<> per TageOverrider contract
    // =========================================================
    ram<val<LENTRY>, LSETS> lram[LWAYS];

    // Pre-extracted regs: written by prefetch(), read by lookup()
    arr<reg<LOOP_TAGW>, LWAYS>   lp_tags;
    arr<reg<1>,         LWAYS>   lp_preds;
    arr<reg<1>,         LWAYS>   lp_confidents;
    arr<reg<LENTRY>,    LWAYS>   lp_saved;
    reg<LLSETS>                  lp_idx;
    reg<1>                       lp_did_lookup;

    // First branch tag — captured on the first lookup() call for the block.
    // update_condbr already consumes branch_pc via fo1(), so we cannot safely
    // re-read it again in save_branch().
    reg<LOOP_TAGW>               lp_branch_tag;

    // TAGE prediction for first branch — saved in lookup() first call.
    // Used in train() for WITHLOOP counter training.
    reg<1>                       lp_base_pred;
    bool                         lp_base_pred_saved = false;  // reset in prefetch()

    // WITHLOOP counter (L-TAGE paper: 7-bit saturating, threshold = 0)
    reg<USE_CTRW, i64>           lp_use_ctr = USE_CTR_INIT;

    // Two-phase training state (train_compute stages, train_write commits)
    arr<reg<LENTRY>, LWAYS>      staged_entries;
    arr<reg<1>,      LWAYS>      staged_write_mask;
    bool                         should_write = false;

#ifdef TAGE_MONITOR
    struct LoopStats {
        u64 lookups = 0;
        u64 hits = 0;
        u64 confident_hits = 0;
        u64 overrides = 0;
        u64 confident_correct = 0;
        u64 confident_wrong = 0;
        u64 use_blocked_confident = 0;
        u64 provider_correct = 0;
        u64 provider_wrong = 0;
        u64 disagree_updates = 0;
        u64 hit_updates = 0;
        u64 alloc_attempts = 0;
        u64 alloc_writes = 0;
        u64 age_dec_events = 0;
        u64 update_writes = 0;
        u64 sc_correct = 0;
        u64 sc_wrong = 0;
    } stats;
#endif

    // =========================================================
    // Entry field accessors (zero HW cost)
    // =========================================================
    static val<LOOP_TAGW>  e_tag (val<LENTRY>& e) { return val<LOOP_TAGW>{e}; }
    static val<1>          e_dir (val<LENTRY>& e) { return val<1>{e >> hard<B_DIR>{}}; }
    static val<LOOP_ITERW> e_num (val<LENTRY>& e) { return val<LOOP_ITERW>{e >> hard<B_NUM>{}}; }
    static val<LOOP_ITERW> e_cur (val<LENTRY>& e) { return val<LOOP_ITERW>{e >> hard<B_CUR>{}}; }
    static val<LOOP_CONFW> e_conf(val<LENTRY>& e) { return val<LOOP_CONFW>{e >> hard<B_CONF>{}}; }
    static val<LOOP_AGEW>  e_age (val<LENTRY>& e) { return val<LOOP_AGEW>{e >> hard<B_AGE>{}}; }

    static val<LENTRY> mk_entry(val<LOOP_TAGW>& tag, val<1>& dir,
                                val<LOOP_ITERW>& num, val<LOOP_ITERW>& cur,
                                val<LOOP_CONFW>& conf, val<LOOP_AGEW>& age) {
        return val<LENTRY>{concat(age, concat(conf, concat(cur, concat(num, concat(dir, tag)))))};
    }

    // =========================================================
    // Helpers
    // =========================================================

    // is_confident: full-saturation threshold by default (paper: 3 consistent trips).
    // Optional long-loop product gate from gshareN_loop.
    static val<1> is_confident(val<LOOP_CONFW>& conf, val<LOOP_ITERW>& num) {
        val<1> strong = (conf >= hard<LOOP_CONF_THRESH>{});
        if constexpr (USE_LONG_LOOP_CONF) {
            val<1> long_loop = ((conf * num) > hard<LONG_LOOP_THRESH>{});
            return strong | long_loop;
        } else {
            (void)num;
            return strong;
        }
    }

    // Priority MUX: first way where m[w]=1. No OR-masking ambiguity.
    static val<WAYW> first_true_way(arr<val<1>, LWAYS>& m) {
        if constexpr (LWAYS == 1) {
            return val<WAYW>{0};
        } else if constexpr (LWAYS == 2) {
            return select(m[0], val<WAYW>{0}, val<WAYW>{1});
        } else {
            return select(m[0], val<WAYW>{0},
                   select(m[1], val<WAYW>{1},
                   select(m[2], val<WAYW>{2},
                                val<WAYW>{3})));
        }
    }

    static val<1> way_is(val<WAYW>& way, u64 w) {
        return way == val<WAYW>{w};
    }

    static val<LENTRY> select_entry(arr<val<LENTRY>, LWAYS>& a, val<WAYW>& way) {
        if constexpr (LWAYS == 1) {
            (void)way; return a[0];
        } else if constexpr (LWAYS == 2) {
            return select(way_is(way, 0), a[0], a[1]);
        } else {
            return select(way_is(way, 0), a[0],
                   select(way_is(way, 1), a[1],
                   select(way_is(way, 2), a[2], a[3])));
        }
    }

    static val<1> select_bit(arr<val<1>, LWAYS>& a, val<WAYW>& way) {
        if constexpr (LWAYS == 1) {
            (void)way; return a[0];
        } else if constexpr (LWAYS == 2) {
            return select(way_is(way, 0), a[0], a[1]);
        } else {
            return select(way_is(way, 0), a[0],
                   select(way_is(way, 1), a[1],
                   select(way_is(way, 2), a[2], a[3])));
        }
    }

    static val<LOOP_CONFW> bump_conf(val<LOOP_CONFW>& c) {
        return select(c >= hard<CONF_MAX>{}, c,
                      val<LOOP_CONFW>{c + CONF_BUMP_STEP});
    }

    // Age helpers: saturating inc/dec
    static val<LOOP_AGEW> sat_age_inc(val<LOOP_AGEW>& a) {
        return select(a >= hard<AGE_MAX>{}, val<LOOP_AGEW>{AGE_MAX},
                                            val<LOOP_AGEW>{a + 1});
    }
    static val<LOOP_AGEW> sat_age_dec(val<LOOP_AGEW>& a) {
        return select(a > hard<0>{}, val<LOOP_AGEW>{a - 1},
                                     val<LOOP_AGEW>{0});
    }

    static val<USE_CTRW, i64> bump_use_up(val<USE_CTRW, i64>& c) {
        return select(c <= hard<USE_MAX - static_cast<i64>(USE_CTR_INC_STEP)>{},
                      val<USE_CTRW, i64>{c + USE_CTR_INC_STEP},
                      val<USE_CTRW, i64>{USE_MAX});
    }
    static val<USE_CTRW, i64> bump_use_down(val<USE_CTRW, i64>& c) {
        return select(c >= hard<USE_MIN + static_cast<i64>(USE_CTR_DEC_STEP)>{},
                      val<USE_CTRW, i64>{c - USE_CTR_DEC_STEP},
                      val<USE_CTRW, i64>{USE_MIN});
    }

    // =========================================================
    // LookupResult
    // =========================================================
    struct LookupResult {
        val<1> candidate;
        val<1> pred;
    };

    // =========================================================
    // prefetch — called from predict1
    //
    // Reads lram[*] at LINE-LEVEL index (same for all instructions in block).
    // Pre-extracts per-way tag, prediction, confidence into regs.
    // All RAM access for the prediction path is done here.
    // =========================================================
    void prefetch(val<64>& block_pc, [[maybe_unused]] u64 raw_pc) {
        lp_base_pred_saved = false;

        val<LLSETS> lidx = val<LLSETS>{block_pc >> hard<LINE_SHIFT>{}};
        lp_idx        = lidx;
        lp_did_lookup = val<1>{1};

        if constexpr (LWAYS > 1) lidx.fanout(hard<LWAYS>{});

        arr<val<LENTRY>, LWAYS> entries = [&](u64 w) -> val<LENTRY> {
            return lram[w].read(lidx);
        };

        entries.fanout(hard<6>{});
        for (u64 w = 0; w < LWAYS; w++) {
            lp_saved[w] = entries[w];

            val<LOOP_TAGW>  tag  = e_tag(entries[w]);
            val<1>          dir  = e_dir(entries[w]);
            val<LOOP_ITERW> num  = e_num(entries[w]);
            val<LOOP_ITERW> cur  = e_cur(entries[w]);
            val<LOOP_CONFW> conf = e_conf(entries[w]);

            lp_tags[w] = tag;

            num.fanout(hard<3>{});
            val<1> learned = (num != hard<0>{});
            val<1> is_last = learned & (cur == num);
            dir.fanout(hard<2>{});
            lp_preds[w]      = select(is_last, ~dir, dir);
            lp_confidents[w] = is_confident(conf, num);
        }

        lp_tags.fanout(hard<FETCH_WIDTH + 2>{});
        lp_preds.fanout(hard<FETCH_WIDTH + 2>{});
        lp_confidents.fanout(hard<FETCH_WIDTH + 2>{});
        lp_use_ctr.fanout(hard<FETCH_WIDTH + 2>{});
    }

    // =========================================================
    // lookup — called from predict2 / reuse_predict2
    //
    // Pure combinational: branch-level tag match against pre-extracted regs.
    // No RAM reads. Saves tage_pred on first call (for USE counter training).
    // =========================================================
    LookupResult lookup(val<64>& inst_pc, val<1>& tage_pred,
                        [[maybe_unused]] val<1>& tage_low_conf) {
        // First-call save: lp_base_pred is the first branch's TAGE prediction.
        // This is the counterfactual reference for WITHLOOP counter training.
        if (!lp_base_pred_saved) {
            lp_base_pred       = tage_pred;
            lp_branch_tag      = val<LOOP_TAGW>{inst_pc >> 2};
            lp_base_pred_saved = true;
        } else {
            // Match the hybrid gshareN_loop design: only the first branch in
            // the fetch block is eligible for loop override. This keeps
            // reuse_predict2 off the loop predictor's critical path.
            return {val<1>{0}, val<1>{0}};
        }

        // Branch-level tag: distinguishes branches within the same cache line
        val<LOOP_TAGW> btag = val<LOOP_TAGW>{inst_pc >> 2};
        if constexpr (LWAYS > 1) btag.fanout(hard<LWAYS>{});

        arr<val<1>, LWAYS> tmatch = [&](u64 w) -> val<1> {
            return lp_tags[w] == btag;
        };
        tmatch.fanout(hard<2>{});

        val<1> hit = tmatch.fold_or();

        // Lookup stays on the P2 critical path, so use masked OR selection
        // instead of priority-encode + MUX. Duplicated tags should be rare
        // because allocation only happens on miss.
        arr<val<1>, LWAYS> pred_masked = [&](u64 w) -> val<1> {
            return tmatch[w] & lp_preds[w];
        };
        val<1> loop_pred = pred_masked.fold_or();

        arr<val<1>, LWAYS> conf_masked = [&](u64 w) -> val<1> {
            return tmatch[w] & lp_confidents[w];
        };
        val<1> confident = conf_masked.fold_or();
        val<1> use_ok    = (lp_use_ctr >= hard<USE_CTR_THRESH>{});
        val<1> candidate = hit & confident & use_ok;

#ifdef TAGE_MONITOR
#ifdef CHEATING_MODE
        stats.lookups++;
        if (static_cast<u64>(hit) != 0) stats.hits++;
        if (static_cast<u64>(hit & confident) != 0) stats.confident_hits++;
        if (static_cast<u64>(candidate) != 0) stats.overrides++;
#endif
#endif

        return {candidate, loop_pred};
    }

    // =========================================================
    // save_branch — called from update_condbr
    //
    // No-op by design: Tage::update_condbr already consumes branch_pc via
    // branch_pc.fo1() before calling the overrider. The first branch tag is
    // captured in lookup() instead, when inst_pc is still safe to read.
    // =========================================================
    void save_branch(val<64>& branch_pc, u64 idx) {
        (void)branch_pc;
        (void)idx;
    }

    // =========================================================
    // train_compute — called from update_cycle, BEFORE need_extra_cycle
    //
    // Fully combinational: tag match, state machine update, USE counter update.
    // Stages per-way write data into HARCOM regs; sets should_write = true.
    // Per-way writes merge hit/alloc path with age-dec path into one write.
    // =========================================================
    void train_compute(arr<val<1>, FETCH_WIDTH>& branch_taken,
                       [[maybe_unused]] arr<val<1>, FETCH_WIDTH>& is_branch,
                       val<1>& mispredict,
                       [[maybe_unused]] val<1>& correct_pred,
                       u64 num_branch)
    {
        if (num_branch == 0) return;

        lp_saved.fanout(hard<6>{});
        mispredict.fanout(hard<2>{});

        val<1> first_dir = branch_taken[0];
        first_dir.fanout(hard<3>{});

        // --- Tag match against saved entries ---
        lp_branch_tag.fanout(hard<LWAYS + 6>{});
        arr<val<1>, LWAYS> umatch = [&](u64 w) -> val<1> {
            return (e_tag(lp_saved[w]) == lp_branch_tag);
        };
        umatch.fanout(hard<2>{});
        val<1>    uhit     = umatch.fold_or();
        uhit.fanout(hard<4>{});
        val<WAYW> uhit_way = first_true_way(umatch);

        arr<val<LENTRY>, LWAYS> saved_vals = [&](u64 w) -> val<LENTRY> {
            return lp_saved[w];
        };
        val<LENTRY> matched = select_entry(saved_vals, uhit_way);
        matched.fanout(hard<5>{});

        val<1>          m_dir  = e_dir(matched);
        val<LOOP_ITERW> m_num  = e_num(matched);
        val<LOOP_ITERW> m_cur  = e_cur(matched);
        val<LOOP_CONFW> m_conf = e_conf(matched);
        val<LOOP_AGEW>  m_age  = e_age(matched);

        m_dir.fanout(hard<4>{});
        m_num.fanout(hard<6>{});
        m_cur.fanout(hard<4>{});
        m_conf.fanout(hard<4>{});
        m_age.fanout(hard<2>{});

        // Recompute loop prediction (needed for was_correct and WITHLOOP training)
        val<1> m_learned    = (m_num != hard<0>{});
        val<1> m_is_last    = m_learned & (m_cur == m_num);
        val<1> loop_pred    = select(m_is_last, ~m_dir, m_dir);

        val<1> m_confident  = is_confident(m_conf, m_num);
        val<1> use_ok       = (lp_use_ctr >= hard<USE_CTR_THRESH>{});
        val<1> was_provider = uhit & m_confident & use_ok;
        val<1> was_correct  = (loop_pred == first_dir);
        was_correct.fanout(hard<2>{});

        // --- WITHLOOP counter update (L-TAGE paper: 7-bit WITHLOOP) ---
        // Train only when loop and base disagree: counterfactual signal.
        // Choosing loop would have changed the outcome, so we can learn from it.
        val<1> disagrees = uhit & m_confident & (loop_pred ^ val<1>{lp_base_pred});
        execute_if(val<1>{lp_did_lookup} & disagrees, [&](){
            val<USE_CTRW, i64> cur = lp_use_ctr;
            lp_use_ctr = select(was_correct,
                                bump_use_up(cur), bump_use_down(cur));
        });

        // --- Age: increment when confidently used as provider (L-TAGE paper) ---
        val<LOOP_AGEW> age_used = select(was_provider, sat_age_inc(m_age), m_age);

        // --- Build candidate write entries ---
        val<1> wrong_conf = INVALIDATE_ON_WRONG_CONF
                          ? (was_provider & ~was_correct) : val<1>{0};

        val<1> continuing   = (first_dir == m_dir);
        val<LOOP_ITERW> new_cur = val<LOOP_ITERW>{m_cur + 1};
        new_cur.fanout(hard<2>{});
        val<1> has_num    = (m_num != hard<0>{});
        val<1> exceeded   = has_num & (new_cur > m_num);
        val<1> trip_match = (m_cur == m_num);
        val<1> num_zero   = (m_num == hard<0>{});
        val<1> too_short  = (m_cur < hard<MIN_TRIP>{});

        val<1> zero_dir = val<1>{0};
        val<LOOP_ITERW> zero_iter = val<LOOP_ITERW>{0};
        val<LOOP_CONFW> zero_conf = val<LOOP_CONFW>{0};
        val<LOOP_AGEW> zero_age = val<LOOP_AGEW>{0};

        val<LENTRY> invalid_entry = mk_entry(
            lp_branch_tag, zero_dir, zero_iter, zero_iter, zero_conf, zero_age);

        // Continuing (same direction): increment current iteration count
        val<LENTRY> continued_entry = mk_entry(
            lp_branch_tag, m_dir, m_num, new_cur, m_conf, age_used);

        // Correct exit (cur==num): bump confidence, reset cur
        val<LOOP_CONFW> bumped = bump_conf(m_conf);
        val<LENTRY> confirm_entry = mk_entry(
            lp_branch_tag, m_dir, m_num, zero_iter, bumped, age_used);

        // First exit (num not yet known): learn trip count from cur
        val<LOOP_CONFW> alloc_conf_init = val<LOOP_CONFW>{ALLOC_CONF_INIT};
        val<LENTRY> learn_entry = mk_entry(
            lp_branch_tag, m_dir, m_cur, zero_iter, alloc_conf_init, m_age);

        val<LENTRY> exit_bad_entry   =
            INVALIDATE_ON_BAD_EXIT  ? invalid_entry : learn_entry;
        val<LENTRY> exit_short_entry =
            INVALIDATE_ON_TOO_SHORT ? invalid_entry : learn_entry;

        // Exit path (direction flipped)
        val<LENTRY> exit_entry = select(too_short, exit_short_entry,
            select(num_zero, learn_entry,
                select(trip_match, confirm_entry, exit_bad_entry)));

        // Continuing path: check for exceeding known trip count
        val<LENTRY> iter_entry = select(
            INVALIDATE_ON_EXCEEDED ? exceeded : val<1>{0},
            invalid_entry, continued_entry);

        val<LENTRY> hit_update = select(wrong_conf, invalid_entry,
            select(continuing, iter_entry, exit_entry));

        // Allocation entry: opposite direction, age=255, conf=0
        val<1> alloc_dir = ~first_dir;
        val<LOOP_AGEW> alloc_age_init = val<LOOP_AGEW>{ALLOC_AGE_INIT};
        val<LENTRY> alloc_entry = mk_entry(
            lp_branch_tag, alloc_dir, zero_iter, zero_iter, alloc_conf_init, alloc_age_init);

        // --- Victim selection: first way with age==0 (L-TAGE paper) ---
        arr<val<LOOP_AGEW>, LWAYS> ages = [&](u64 w) -> val<LOOP_AGEW> {
            return e_age(lp_saved[w]);
        };
        ages.fanout(hard<2>{});
        arr<val<1>, LWAYS> age_zero = [&](u64 w) -> val<1> {
            return (ages[w] == hard<0>{});
        };
        age_zero.fanout(hard<2>{});
        val<1>    found_victim = age_zero.fold_or();
        found_victim.fanout(hard<2>{});
        val<WAYW> victim_way   = first_true_way(age_zero);

        // is_age_dec: miss + no free slot → decrement all ways (L-TAGE paper)
        val<1> is_alloc   = ~uhit & mispredict & found_victim;
        val<1> is_age_dec = ~uhit & mispredict & ~found_victim;
        val<1> do_single  = uhit | is_alloc;
        do_single.fanout(hard<LWAYS>{});
        is_age_dec.fanout(hard<LWAYS>{});

        val<WAYW>   single_way  = select(uhit, uhit_way, victim_way);
        single_way.fanout(hard<LWAYS>{});
        val<LENTRY> single_data = select(uhit, hit_update, alloc_entry);
        single_data.fanout(hard<LWAYS>{});

        // Per-way: stage write data into regs for train_write().
        // Merge hit/alloc path with age-dec path into one write per way.
        for (u64 w = 0; w < LWAYS; w++) {
            val<LOOP_AGEW> old_age = ages[w];
            old_age.fanout(hard<2>{});

            val<LOOP_TAGW> aged_tag = e_tag(lp_saved[w]);
            val<1> aged_dir = e_dir(lp_saved[w]);
            val<LOOP_ITERW> aged_num = e_num(lp_saved[w]);
            val<LOOP_ITERW> aged_cur = e_cur(lp_saved[w]);
            val<LOOP_CONFW> aged_conf = e_conf(lp_saved[w]);
            val<LOOP_AGEW> aged_age = sat_age_dec(old_age);
            val<LENTRY> aged_entry = mk_entry(
                aged_tag, aged_dir, aged_num, aged_cur, aged_conf, aged_age);

            val<1> is_single_way = do_single & way_is(single_way, w);
            val<1> do_write      = is_single_way | is_age_dec;

            staged_entries[w]    = select(is_single_way, single_data, aged_entry);
            staged_write_mask[w] = do_write;
        }

#ifdef TAGE_MONITOR
#ifdef CHEATING_MODE
        if (static_cast<u64>(uhit & m_confident & was_correct) != 0) stats.confident_correct++;
        if (static_cast<u64>(uhit & m_confident & ~was_correct) != 0) stats.confident_wrong++;
        if (static_cast<u64>(uhit & m_confident & ~use_ok) != 0) stats.use_blocked_confident++;
        if (static_cast<u64>(uhit) != 0) stats.hit_updates++;
        if (static_cast<u64>(disagrees) != 0) stats.disagree_updates++;
        if (static_cast<u64>(was_provider & was_correct) != 0) stats.provider_correct++;
        if (static_cast<u64>(was_provider & ~was_correct) != 0) stats.provider_wrong++;
        if (static_cast<u64>((~uhit) & mispredict) != 0) stats.alloc_attempts++;
        if (static_cast<u64>(is_alloc) != 0) stats.alloc_writes++;
        if (static_cast<u64>(is_age_dec) != 0) stats.age_dec_events++;
#endif
#endif

        should_write = true;
    }

    // =========================================================
    // train_write — called from update_cycle, AFTER need_extra_cycle
    //
    // Commits staged writes to lram. Gated per-way by lp_did_lookup
    // and staged_write_mask so no spurious writes occur.
    // =========================================================
    void train_write() {
        if (!should_write) return;

        lp_did_lookup.fanout(hard<LWAYS>{});
        lp_idx.fanout(hard<LWAYS>{});

        for (u64 w = 0; w < LWAYS; w++) {
            execute_if(lp_did_lookup & staged_write_mask[w], [&](){
                lram[w].write(val<LLSETS>{lp_idx}, val<LENTRY>{staged_entries[w]});
            });
        }

#ifdef TAGE_MONITOR
#ifdef CHEATING_MODE
        for (u64 w = 0; w < LWAYS; w++) {
            if (static_cast<u64>(staged_write_mask[w]) != 0) {
                stats.update_writes++;
            }
        }
#endif
#endif

        should_write = false;
    }

    void print_stats() const {
#ifdef TAGE_MONITOR
        auto pct = [](u64 num, u64 den) -> double {
            return den ? (100.0 * static_cast<double>(num) / static_cast<double>(den)) : 0.0;
        };
        std::cerr << "\n=== Loop Overrider Detail ===\n"
                  << "  Lookups: " << stats.lookups
                  << "  Hits: " << stats.hits << " (" << pct(stats.hits, stats.lookups) << "%)\n"
                  << "  Confident hits: " << stats.confident_hits << " (" << pct(stats.confident_hits, stats.hits) << "% of hits)\n"
                  << "  Confident correct: " << stats.confident_correct
                  << "  Confident wrong: " << stats.confident_wrong;
        u64 confident_total = stats.confident_correct + stats.confident_wrong;
        std::cerr << " (" << pct(stats.confident_correct, confident_total) << "% correct)\n"
                  << "  Use-blocked confident hits: " << stats.use_blocked_confident << "\n"
                  << "  Overrides: " << stats.overrides << " (" << pct(stats.overrides, stats.lookups) << "% of lookups)\n"
                  << "  Provider correct: " << stats.provider_correct
                  << "  Provider wrong: " << stats.provider_wrong;
        u64 provider_total = stats.provider_correct + stats.provider_wrong;
        std::cerr << " (" << pct(stats.provider_correct, provider_total) << "% correct)\n"
                  << "  Hit updates: " << stats.hit_updates
                  << "  Disagree updates: " << stats.disagree_updates << "\n"
                  << "  Alloc attempts: " << stats.alloc_attempts
                  << "  Alloc writes: " << stats.alloc_writes
                  << "  Age-dec events: " << stats.age_dec_events << "\n"
                  << "  Committed writes: " << stats.update_writes << "\n";
#endif
    }
};
