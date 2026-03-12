#include "../cbp.hpp"
#include "../harcom.hpp"

using namespace hcm;

// =========================================================================
// gshareN_loop: gshare (p1) + hybrid loop predictor override (p2)
//
// Architecture (hybrid):
//   predict1: reads loop RAMs at LINE-LEVEL index (same for all instr in
//             a cache line). Pre-extracts per-way tags, predictions, and
//             confidence flags into reg arrays.
//   predict2 / reuse_predict2: PURE COMBINATIONAL (no reg writes).
//             Computes tag from actual inst_pc, tag-matches against the
//             pre-extracted tags, selects the matching way's pre-computed
//             prediction.
//   update_condbr: saves first branch's tag for update_cycle.
//   update_cycle: recomputes tag match from saved branch tag, updates the
//             matching entry in the loop RAMs.
//
// Key idea:
//   - RAM index is line-level, so predict1 can prefetch using block-start PC.
//   - Tag is branch-level, so predict2/update_cycle still distinguish
//     different branches within the same line.
//
// Default parameters preserve the current functionality/behavior.
// ==========================================================================

template<
    // ----------------------------
    // Base gshare parameters
    // ----------------------------
    u64 LOGG = 15,
    u64 GHIST = 8,
    u64 LOGN = 2,
    u64 LINEINST_ = 64,

    // ----------------------------
    // Loop table geometry
    // ----------------------------
    u64 LOOP_LOGSETS = 4,
    u64 LOOP_TAGW = 10,
    u64 LOOP_ITERW = 8,
    u64 LOOP_CONFW = 2,
    u64 LOOP_AGEW = 4,
    u64 LOOP_WAYS = 2,          // supported: 1, 2, 4

    // ----------------------------
    // Loop prediction policy
    // ----------------------------
    u64 MIN_TRIP = 4,
    bool USE_LONG_LOOP_CONF = true,
    u64 LONG_LOOP_CONF_PRODUCT_THRESH = 128,

    // Confidence threshold used by the strong-counter portion of
    // is_confident(). Default preserves current behavior because the
    // original code used conf == CONF_MAX, and the long-loop term is OR'ed in.
    u64 LOOP_CONF_USE = ((1ULL << LOOP_CONFW) - 1),

    // ----------------------------
    // Use counter policy
    // ----------------------------
    u64 USE_CTRW = 5,
    i64 USE_CTR_INIT = 0,
    i64 USE_CTR_OVERRIDE_THRESH = 0,
    i64 USE_CTR_INC_STEP = 1,
    i64 USE_CTR_DEC_STEP = 1,

    // ----------------------------
    // Allocation / update policy
    // ----------------------------
    u64 ALLOC_CONF_INIT = 0,
    u64 ALLOC_AGE_INIT = ((1ULL << LOOP_AGEW) - 1),
    u64 CONF_BUMP_STEP = 1,
    bool INVALIDATE_ON_WRONG_CONF = true,
    bool INVALIDATE_ON_EXCEEDED_ITER = true,
    bool INVALIDATE_ON_BAD_EXIT_MATCH = true,
    bool INVALIDATE_ON_TOO_SHORT = true,
    bool REQUIRE_FIRST_BRANCH_MISPRED_FOR_ALLOC = true,
    bool FORCE_EXTRA_CYCLE_ON_LOOKUP = true
>
struct gshareN_loop : predictor {

    // =====================================================================
    // Base gshare constants
    // =====================================================================
    static constexpr u64 N = 1 << LOGN;
    static_assert(LOGG > LOGN);
    static constexpr u64 index_bits = LOGG - LOGN;

    static constexpr u64 LINEINST = LINEINST_;
    static_assert(std::has_single_bit(LINEINST));
    static constexpr u64 LOGLINEINST = std::bit_width(LINEINST - 1);

    // =====================================================================
    // Loop predictor constants
    // =====================================================================
    static constexpr u64 LWAYS = LOOP_WAYS;
    static_assert(LWAYS == 1 || LWAYS == 2 || LWAYS == 4);

    static constexpr u64 LSETS  = 1 << LOOP_LOGSETS;
    static constexpr u64 LLSETS = LOOP_LOGSETS;

    static constexpr u64 B_TAG  = 0;
    static constexpr u64 B_DIR  = LOOP_TAGW;
    static constexpr u64 B_NUM  = B_DIR + 1;
    static constexpr u64 B_CUR  = B_NUM + LOOP_ITERW;
    static constexpr u64 B_CONF = B_CUR + LOOP_ITERW;
    static constexpr u64 B_AGE  = B_CONF + LOOP_CONFW;
    static constexpr u64 LENTRY = B_AGE + LOOP_AGEW;

    static constexpr u64 CONF_MAX = (1ULL << LOOP_CONFW) - 1;
    static constexpr u64 AGE_MAX  = (1ULL << LOOP_AGEW) - 1;

    static constexpr u64 LINE_SHIFT = 2 + LOGLINEINST;
    static constexpr u64 LOOP_WAYW = (LWAYS <= 1) ? 1 : std::bit_width(LWAYS - 1);

    static_assert(LOOP_CONF_USE <= CONF_MAX);
    static_assert(ALLOC_CONF_INIT <= CONF_MAX);
    static_assert(ALLOC_AGE_INIT <= AGE_MAX);
    static_assert(CONF_BUMP_STEP >= 1);
    static_assert(MIN_TRIP >= 1);
    static_assert(USE_CTRW >= 2);

    static constexpr i64 USE_CTR_MIN = -(1LL << (USE_CTRW - 1));
    static constexpr i64 USE_CTR_MAX =  (1LL << (USE_CTRW - 1)) - 1;

    static_assert(USE_CTR_INIT >= USE_CTR_MIN && USE_CTR_INIT <= USE_CTR_MAX);
    static_assert(USE_CTR_OVERRIDE_THRESH >= USE_CTR_MIN && USE_CTR_OVERRIDE_THRESH <= USE_CTR_MAX);
    static_assert(USE_CTR_INC_STEP >= 1);
    static_assert(USE_CTR_DEC_STEP >= 1);

    // =====================================================================
    // Base gshare state
    // =====================================================================
    reg<GHIST> global_history;
    reg<N> X;
    reg<index_bits> index;
    reg<N> unordered_pred;
    arr<reg<1>, N> pred;
    reg<1> true_block = 1;
    reg<LINEINST> block_entry;
    reg<N> rank;

    u64 num_branch = 0;
    u64 block_size = 0;
    arr<reg<1>, N> branch_dir;

    ram<val<N>, (1 << index_bits)> ctr_hi;
    ram<val<1>, (1 << index_bits)> ctr_lo[N];

    // =====================================================================
    // Loop predictor storage
    // =====================================================================
    ram<val<LENTRY>, LSETS> lram[LWAYS];

    // Saved from predict1 for p2/update_cycle
    arr<reg<LOOP_TAGW>, LWAYS> lp_tags;
    arr<reg<1>, LWAYS>         lp_preds;
    arr<reg<1>, LWAYS>         lp_confidents;
    arr<reg<LENTRY>, LWAYS>    lp_saved;

    reg<1> lp_did_lookup;
    reg<LLSETS> lp_idx;
    reg<1> lp_base_pred;

    // Saved from update_condbr for update_cycle
    reg<LOOP_TAGW> lp_branch_tag;

    // Use-gating counter
    reg<USE_CTRW, i64> lp_use_ctr = USE_CTR_INIT;

    // =====================================================================
    // Helpers (zero hardware cost)
    // =====================================================================
    static val<LOOP_TAGW>  e_tag (val<LENTRY> e) { return val<LOOP_TAGW>{e}; }
    static val<1>          e_dir (val<LENTRY> e) { return val<1>{e >> hard<B_DIR>{}}; }
    static val<LOOP_ITERW> e_num (val<LENTRY> e) { return val<LOOP_ITERW>{e >> hard<B_NUM>{}}; }
    static val<LOOP_ITERW> e_cur (val<LENTRY> e) { return val<LOOP_ITERW>{e >> hard<B_CUR>{}}; }
    static val<LOOP_CONFW> e_conf(val<LENTRY> e) { return val<LOOP_CONFW>{e >> hard<B_CONF>{}}; }
    static val<LOOP_AGEW>  e_age (val<LENTRY> e) { return val<LOOP_AGEW>{e >> hard<B_AGE>{}}; }

    static val<LENTRY> mk_entry(val<LOOP_TAGW> tag, val<1> dir,
                                val<LOOP_ITERW> num, val<LOOP_ITERW> cur,
                                val<LOOP_CONFW> conf, val<LOOP_AGEW> age)
    {
        return val<LENTRY>{
            concat(age, concat(conf, concat(cur, concat(num, concat(dir, tag)))))
        };
    }

    static val<1> is_confident(val<LOOP_CONFW> conf, val<LOOP_ITERW> num)
    {
        val<1> strong = (conf >= hard<LOOP_CONF_USE>{});
        if constexpr (USE_LONG_LOOP_CONF) {
            auto product = conf * num;
            val<1> long_loop = (product > hard<LONG_LOOP_CONF_PRODUCT_THRESH>{});
            return strong | long_loop;
        } else {
            (void)num;
            return strong;
        }
    }

    static val<LOOP_WAYW> first_true_way(arr<val<1>, LWAYS> m)
    {
        if constexpr (LWAYS == 1) {
            return val<LOOP_WAYW>{0};
        } else if constexpr (LWAYS == 2) {
            return select(m[0], val<LOOP_WAYW>{0}, val<LOOP_WAYW>{1});
        } else {
            return select(m[0], val<LOOP_WAYW>{0},
                   select(m[1], val<LOOP_WAYW>{1},
                   select(m[2], val<LOOP_WAYW>{2},
                                  val<LOOP_WAYW>{3})));
        }
    }

    static val<1> way_is(val<LOOP_WAYW> way, u64 w)
    {
        return way == val<LOOP_WAYW>{w};
    }

    static val<LENTRY> select_entry(arr<val<LENTRY>, LWAYS> a, val<LOOP_WAYW> way)
    {
        if constexpr (LWAYS == 1) {
            (void)way;
            return a[0];
        } else if constexpr (LWAYS == 2) {
            return select(way_is(way, 0), a[0], a[1]);
        } else {
            return select(way_is(way, 0), a[0],
                   select(way_is(way, 1), a[1],
                   select(way_is(way, 2), a[2], a[3])));
        }
    }

    static val<1> select_bit(arr<val<1>, LWAYS> a, val<LOOP_WAYW> way)
    {
        if constexpr (LWAYS == 1) {
            (void)way;
            return a[0];
        } else if constexpr (LWAYS == 2) {
            return select(way_is(way, 0), a[0], a[1]);
        } else {
            return select(way_is(way, 0), a[0],
                   select(way_is(way, 1), a[1],
                   select(way_is(way, 2), a[2], a[3])));
        }
    }

    static val<LOOP_WAYW> victim_way_from_age0(arr<val<1>, LWAYS> age_zero)
    {
        return first_true_way(age_zero);
    }

    static val<LOOP_CONFW> bump_conf(val<LOOP_CONFW> c)
    {
        return select(
            c >= hard<CONF_MAX>{},
            c,
            val<LOOP_CONFW>{c + CONF_BUMP_STEP}
        );
    }

    static val<USE_CTRW, i64> bump_use_ctr_up(val<USE_CTRW, i64> cur)
    {
        return select(cur <= hard<USE_CTR_MAX - static_cast<i64>(USE_CTR_INC_STEP)>{},
                      val<USE_CTRW, i64>{cur + USE_CTR_INC_STEP},
                      val<USE_CTRW, i64>{USE_CTR_MAX});
    }

    static val<USE_CTRW, i64> bump_use_ctr_down(val<USE_CTRW, i64> cur)
    {
        return select(cur >= hard<USE_CTR_MIN + static_cast<i64>(USE_CTR_DEC_STEP)>{},
                      val<USE_CTRW, i64>{cur - USE_CTR_DEC_STEP},
                      val<USE_CTRW, i64>{USE_CTR_MIN});
    }

    val<1> line_end() { return block_entry >> (LINEINST - block_size); }
    val<1> last_pred() { return rank.rotate_left(num_branch); }

    // =====================================================================
    // predict1: gshare + loop RAM read + pre-extraction
    // =====================================================================
    val<1> predict1([[maybe_unused]] val<64> inst_pc)
    {
        inst_pc.fanout(hard<5>{});
        global_history.fanout(hard<2>{});
        true_block.fanout(hard<6>{});

        // --- Base gshare ---
        block_entry = select(true_block,
                             val<LOGLINEINST>{inst_pc >> 2}.decode().concat(),
                             block_entry << block_size);
        block_entry.fanout(hard<LINEINST + N + 1>{});
        rank = select(true_block,
                      val<N>{1},
                      rank.rotate_left(num_branch));
        rank.fanout(hard<N + 1>{});
        X = select(true_block,
                   val<LOGN>{inst_pc >> 2}.decode().concat(),
                   X.rotate_left(num_branch));
        X.fanout(hard<N>{});
        execute_if(true_block, [&](){
            val<index_bits> pc_bits = inst_pc >> (LOGN + 2);
            if constexpr (GHIST <= index_bits) {
                index = pc_bits.fo1() ^ (val<index_bits>{global_history} << (index_bits - GHIST));
            } else {
                index = global_history.make_array(val<index_bits>{}).append(pc_bits.fo1()).fold_xor();
            }
            unordered_pred = ctr_hi.read(index);
        });
        unordered_pred.fanout(hard<N>{});
        for (u64 i = 0; i < N; i++) {
            pred[i] = (unordered_pred & X.rotate_left(i)) != hard<0>{};
        }
        pred.fanout(hard<LINEINST * 2>{});

        // --- Loop table pre-read ---
        lp_did_lookup = select(true_block, val<1>{1}, val<1>{0});

        execute_if(true_block, [&](){
            // line-level index; identical for all instructions in the line
            val<LLSETS> lidx = val<LLSETS>{inst_pc >> hard<LINE_SHIFT>{}};
            lp_idx = lidx;

            arr<val<LENTRY>, LWAYS> entries = [&](u64 w) -> val<LENTRY> {
                return lram[w].read(lidx);
            };

            for (u64 w = 0; w < LWAYS; w++) {
                lp_saved[w] = entries[w];
            }

            entries.fanout(hard<6>{});
            for (u64 w = 0; w < LWAYS; w++) {
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
                lp_preds[w] = select(is_last, ~dir, dir);

                lp_confidents[w] = is_confident(conf, num);
            }

            lp_base_pred = pred[0];
        });

        lp_tags.fanout(hard<LINEINST + 2>{});
        lp_preds.fanout(hard<LINEINST + 2>{});
        lp_confidents.fanout(hard<LINEINST + 2>{});
        lp_use_ctr.fanout(hard<LINEINST + 2>{});

        block_size = 1;
        num_branch = 0;

        reuse_prediction(~line_end());
        return pred[num_branch];
    }

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc)
    {
        block_size++;
        reuse_prediction(~line_end());
        return pred[num_branch];
    }

    // =====================================================================
    // p2 / reuse_p2: pure combinational override
    // =====================================================================
    val<1> loop_override(val<64> inst_pc)
    {
        val<1> base = pred[num_branch];
        if (num_branch != 0) return base;

        // branch-level tag; includes intra-line offset
        val<LOOP_TAGW> ltag = val<LOOP_TAGW>{inst_pc >> 2};

        ltag.fanout(hard<LWAYS>{});
        arr<val<1>, LWAYS> tmatch = [&](u64 w) -> val<1> {
            return (lp_tags[w] == ltag);
        };

        tmatch.fanout(hard<2>{});
        val<1> hit = tmatch.fold_or();
        val<LOOP_WAYW> hit_way = first_true_way(tmatch);

        val<1> loop_pred = select_bit(lp_preds, hit_way);
        val<1> confident = select_bit(lp_confidents, hit_way);

        val<1> use_nonneg = (lp_use_ctr >= hard<USE_CTR_OVERRIDE_THRESH>{});

        val<1> do_override = hit & confident & use_nonneg;
        return select(do_override, loop_pred, base);
    }

    val<1> predict2([[maybe_unused]] val<64> inst_pc)
    {
        inst_pc.fanout(hard<2>{});
        return loop_override(inst_pc);
    }

    val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc)
    {
        inst_pc.fanout(hard<2>{});
        return loop_override(inst_pc);
    }

    // =====================================================================
    // update_condbr: save first branch tag for update_cycle
    // =====================================================================
    void update_condbr([[maybe_unused]] val<64> branch_pc, val<1> taken,
                       [[maybe_unused]] val<64> next_pc)
    {
        assert(num_branch < N);
        if (num_branch == 0) {
            lp_branch_tag = val<LOOP_TAGW>{branch_pc >> 2};
        }
        branch_dir[num_branch] = taken.fo1();
        num_branch++;
        reuse_prediction(~(line_end() | last_pred()));
    }

    // =====================================================================
    // update_cycle: gshare + loop predictor update
    // =====================================================================
    void update_cycle([[maybe_unused]] instruction_info &block_end_info)
    {
        val<1>  &mispredict = block_end_info.is_mispredict;
        val<64> &next_pc    = block_end_info.next_pc;

        if (num_branch == 0) {
            global_history = (global_history << 1) ^ val<GHIST>{next_pc.fo1() >> 2};
            true_block = 1;
            return;
        }

        // --- Base gshare update ---
        static_assert(N <= 64);
        index.fanout(hard<N * 3>{});
        branch_dir.fanout(hard<N>{});
        mispredict.fanout(hard<N + 3 + 2>{});
        lp_did_lookup.fanout(hard<LWAYS + 3>{});

        arr<val<1>, N> access = arr<val<N>, N>{[&](u64 i)->val<N> {
            return X.rotate_left(i) & val<N>{-(i < num_branch)};
        }}.fold_or().make_array(val<1>{});

        val<N> misp_bank = X.rotate_left(num_branch - 1)
                         & mispredict.replicate(hard<N>{}).concat();
        arr<val<1>, N> mispredicted = misp_bank.fo1().make_array(val<1>{});
        mispredicted.fanout(hard<2>{});

        arr<val<1>, N> weak = [&](u64 i){
            return execute_if(mispredicted[i], [&](){
                return ctr_lo[i].read(index);
            });
        };

        if constexpr (FORCE_EXTRA_CYCLE_ON_LOOKUP) {
            need_extra_cycle(mispredict | lp_did_lookup);
        } else {
            need_extra_cycle(mispredict);
        }

        execute_if(mispredict, [&](){
            arr<val<1>, N> stored_pred = unordered_pred.make_array(val<1>{});
            arr<val<1>, N> bundle = [&](u64 i){
                return select(weak[i].fo1(),
                              branch_dir[num_branch - 1],
                              stored_pred[i].fo1());
            };
            ctr_hi.write(index, bundle.fo1().concat());
        });

        for (u64 i = 0; i < N; i++) {
            execute_if(access[i].fo1(), [&](){
                ctr_lo[i].write(index, mispredicted[i].fo1());
            });
        }

        true_block = arr<val<1>, 4>{
            ~mispredict, branch_dir[num_branch - 1], last_pred(), line_end()
        }.fold_or();

        execute_if(true_block, [&](){
            global_history = (global_history << 1)
                           ^ val<GHIST>{next_pc.fo1() >> 2};
        });

        // ==============================================================
        // Loop predictor update
        // ==============================================================
        val<1> first_dir = branch_dir[0];
        first_dir.fanout(hard<4>{});

        lp_saved.fanout(hard<6>{});
        lp_idx.fanout(hard<LWAYS + 1>{});

        // Re-derive tag match from saved branch tag
        lp_branch_tag.fanout(hard<LWAYS + 6>{});
        arr<val<1>, LWAYS> umatch = [&](u64 w) -> val<1> {
            return (e_tag(lp_saved[w]) == lp_branch_tag);
        };
        umatch.fanout(hard<2>{});
        val<1> uhit = umatch.fold_or();
        uhit.fanout(hard<4>{});
        val<LOOP_WAYW> uhit_way = first_true_way(umatch);

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
        m_age.fanout(hard<3>{});

        // Recompute provider / loop prediction
        val<1> m_confident = is_confident(m_conf, m_num);
        val<1> use_nonneg = (lp_use_ctr >= hard<USE_CTR_OVERRIDE_THRESH>{});
        val<1> was_provider = uhit & m_confident & use_nonneg;

        val<1> m_learned = (m_num != hard<0>{});
        val<1> m_is_last = m_learned & (m_cur == m_num);
        val<1> loop_pred = select(m_is_last, ~m_dir, m_dir);

        val<1> was_correct = (loop_pred == first_dir);
        was_correct.fanout(hard<2>{});

        // --- Use counter update ---
        val<1> disagrees = uhit & m_confident & (loop_pred ^ lp_base_pred.fo1());
        execute_if(lp_did_lookup & disagrees, [&](){
            val<USE_CTRW, i64> cur = lp_use_ctr;
            lp_use_ctr = select(was_correct, bump_use_ctr_up(cur), bump_use_ctr_down(cur));
        });

        // --- Pre-compute update entries ---
        val<1> first_mispred = REQUIRE_FIRST_BRANCH_MISPRED_FOR_ALLOC
                             ? (mispredict & val<1>{num_branch == 1})
                             : mispredict;
        first_mispred.fanout(hard<2>{});

        val<1> wrong_conf = INVALIDATE_ON_WRONG_CONF ? (was_provider & ~was_correct) : val<1>{0};
        val<1> continuing = (first_dir == m_dir);

        val<LOOP_ITERW> new_cur = val<LOOP_ITERW>{m_cur + 1};
        new_cur.fanout(hard<2>{});
        val<1> has_numiter = (m_num != hard<0>{});
        val<1> exceeded = has_numiter & (new_cur > m_num);

        val<LENTRY> continued_entry = mk_entry(
            lp_branch_tag, m_dir, m_num, new_cur, m_conf, m_age);

        val<1> numiter_zero = (m_num == hard<0>{});
        val<1> trip_matches = (m_cur == m_num);
        val<1> too_short = (m_cur < hard<MIN_TRIP>{});

        val<LENTRY> learn_entry = mk_entry(
            lp_branch_tag, m_dir, m_cur,
            val<LOOP_ITERW>{0}, val<LOOP_CONFW>{ALLOC_CONF_INIT}, m_age);

        val<LOOP_CONFW> bumped = bump_conf(m_conf);
        val<LENTRY> confirm_entry = mk_entry(
            lp_branch_tag, m_dir, m_num,
            val<LOOP_ITERW>{0}, bumped, m_age);

        val<LENTRY> invalid_entry = mk_entry(
            lp_branch_tag, val<1>{0}, val<LOOP_ITERW>{0},
            val<LOOP_ITERW>{0}, val<LOOP_CONFW>{0}, val<LOOP_AGEW>{0});

        val<LENTRY> exit_bad_entry =
            INVALIDATE_ON_BAD_EXIT_MATCH ? invalid_entry : learn_entry;

        val<LENTRY> exit_short_entry =
            INVALIDATE_ON_TOO_SHORT ? invalid_entry : learn_entry;

        val<LENTRY> exit_entry = select(too_short, exit_short_entry,
            select(numiter_zero, learn_entry,
                select(trip_matches, confirm_entry, exit_bad_entry)));

        val<LENTRY> iter_entry = select(
            INVALIDATE_ON_EXCEEDED_ITER ? exceeded : val<1>{0},
            invalid_entry,
            continued_entry
        );

        val<LENTRY> hit_update = select(wrong_conf, invalid_entry,
            select(continuing, iter_entry, exit_entry));

        // Allocation entry
        val<LENTRY> alloc_entry = mk_entry(
            lp_branch_tag, ~first_dir,
            val<LOOP_ITERW>{0}, val<LOOP_ITERW>{0},
            val<LOOP_CONFW>{ALLOC_CONF_INIT}, val<LOOP_AGEW>{ALLOC_AGE_INIT});

        // Victim finding
        arr<val<LOOP_AGEW>, LWAYS> ages = [&](u64 w) -> val<LOOP_AGEW> {
            return e_age(lp_saved[w]);
        };
        ages.fanout(hard<2>{});
        arr<val<1>, LWAYS> age_zero = [&](u64 w) -> val<1> {
            return (ages[w] == hard<0>{});
        };
        age_zero.fanout(hard<2>{});
        val<1> found_victim = age_zero.fold_or();
        found_victim.fanout(hard<2>{});
        val<LOOP_WAYW> victim_way = victim_way_from_age0(age_zero);

        // Merge write targets
        val<1> is_hit_write = uhit;
        is_hit_write.fanout(hard<3>{});
        val<1> is_alloc = ~uhit & first_mispred & found_victim;
        val<1> is_age_dec = ~uhit & first_mispred & ~found_victim;

        val<LOOP_WAYW> single_way = select(is_hit_write, uhit_way, victim_way);
        val<LENTRY> single_data = select(is_hit_write, hit_update, alloc_entry);
        val<1> do_single = is_hit_write | is_alloc;

        do_single.fanout(hard<LWAYS>{});
        is_age_dec.fanout(hard<LWAYS>{});
        single_data.fanout(hard<LWAYS>{});
        lp_idx.fanout(hard<LWAYS>{});

        // --- Write per-way RAMs ---
        for (u64 w = 0; w < LWAYS; w++) {
            val<LOOP_AGEW> old_age = ages[w];
            old_age.fanout(hard<2>{});
            val<LOOP_AGEW> dec_age = select(old_age > hard<0>{},
                val<LOOP_AGEW>{old_age - 1}, old_age);

            val<LENTRY> aged = mk_entry(
                e_tag(lp_saved[w]), e_dir(lp_saved[w]),
                e_num(lp_saved[w]), e_cur(lp_saved[w]),
                e_conf(lp_saved[w]), dec_age);

            val<1> do_write = (do_single & way_is(single_way, w)) | is_age_dec;
            val<1> use_single = do_single & way_is(single_way, w);
            val<LENTRY> wdata = select(use_single, single_data, aged);

            execute_if(lp_did_lookup & do_write, [&](){
                lram[w].write(lp_idx, wdata);
            });
        }
    }
};