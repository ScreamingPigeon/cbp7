#pragma once
// LoopPredictor — detects repeating branch patterns and overrides TAGE at loop exits.
//
// Implements the TageOverrider concept. Plugs into TageImpl via template parameter.
// Design follows TAGE-SC-L literature (Seznec CBP2025, Ros deep dive).
//
// The loop predictor detects branches that are taken N times then not-taken once
// (or vice versa). When confident, it overrides TAGE's prediction at the loop
// exit boundary where TAGE is typically wrong.

#include "../../harcom.hpp"
#include "../common.hpp"
#include "TageOverrider.hpp"

using namespace hcm;

template <u64 TABLE_SIZE = 64,      // total entries
          u64 ASSOC = 4,            // ways per set
          u64 TAG_BITS = 10,        // PC tag width
          u64 ITER_BITS = 10,       // iteration counter width
          u64 CONF_BITS = 3,        // confidence counter width (0..7)
          u64 FETCH_WIDTH = 16,     // from parent TAGE
          u64 BINDEX_BITS = 8>      // index bits from parent
struct LoopPredictor {

  static constexpr bool ENABLED = true;

  // ======== Derived constants ========
  static constexpr u64 NUM_SETS = TABLE_SIZE / ASSOC;
  static constexpr u64 SET_BITS = clog2(NUM_SETS);
  static constexpr u64 AGE_BITS = 3;
  static constexpr u64 CONF_THRESHOLD = 2;  // override when conf > this

  // Entry layout (LSB-first for split<>):
  //   [TAG_BITS] [ITER_BITS] [ITER_BITS] [CONF_BITS] [1 dir] [AGE_BITS]
  static constexpr u64 ENTRY_BITS =
      TAG_BITS + ITER_BITS + ITER_BITS + CONF_BITS + 1 + AGE_BITS;

  static constexpr u64 RWRAM_BANKS = 2;  // banking for rwram same-cycle read+write

  static_assert(NUM_SETS >= 2, "Need at least 2 sets");
  static_assert(ASSOC >= 1, "Need at least 1 way");
  static_assert(TABLE_SIZE == NUM_SETS * ASSOC, "TABLE_SIZE must be NUM_SETS * ASSOC");
  static_assert(NUM_SETS / RWRAM_BANKS > 1, "Need at least 2 entries per bank for rwram");

  // ======== HARCOM Storage ========

  // One rwram per way — supports same-cycle read+write via banking
  rwram<ENTRY_BITS, NUM_SETS, RWRAM_BANKS> loop_ram[ASSOC]{"loop"};

  // Lookup results (set in lookup, consumed in override_predict and train)
  arr<reg<ENTRY_BITS>, ASSOC> read_entries;
  reg<SET_BITS> lookup_set;
  reg<TAG_BITS> lookup_tag;
  reg<1> hit;             // any way matched
  reg<std::max(u64(1), clog2(ASSOC))> hit_way;  // which way matched
  reg<1> loop_has_candidate;  // confident match found (pre-computed in lookup)

  // Lookup result (returned as vals — no regs on critical path)
  // TageImpl applies the TAGE-confidence gate on top.
  struct LookupResult {
    val<1> candidate;  // confident tag match found
    val<1> pred;       // loop prediction (valid when candidate=1)
  };

  // ======== lookup() — parallel with TAGE reads ========
  // Called at START of predict2, concurrent with TAGE table reads.
  // Does RAM read + unpack + tag match + confidence + prediction.
  // Returns LookupResult (vals) — TageImpl applies the TAGE gate later.
  // regs, so the reg timing only reflects the final write, not a read-back.

  LookupResult lookup(val<64> &inst_pc, [[maybe_unused]] reg<BINDEX_BITS> &bindex) {
    // Hash PC
    val<SET_BITS> set_idx = val<SET_BITS>{inst_pc >> 2};
    lookup_set = set_idx;
    val<TAG_BITS> tag = val<TAG_BITS>{inst_pc >> (2 + SET_BITS)};
    tag.fanout(hard<ASSOC + 1>{});
    lookup_tag = tag;
    if constexpr (ASSOC > 1) {
      set_idx.fanout(hard<ASSOC>{});
    }

    // Per-way: single split, compute tag/conf/pred in one pass
    arr<val<3>, ASSOC> way_results = [&](u64 w) -> val<3> {
      val<ENTRY_BITS> entry = loop_ram[w].read(set_idx);
      read_entries[w] = entry;
      auto [e_tag, e_nb_iter, e_cur_iter, e_conf, e_dir, e_age] =
          split<TAG_BITS, ITER_BITS, ITER_BITS, CONF_BITS, 1, AGE_BITS>(entry);
      val<1> tag_match = (val<TAG_BITS>{e_tag} == tag);
      val<1> confident = (val<CONF_BITS>{e_conf} > hard<CONF_THRESHOLD>{});
      val<ITER_BITS> next_iter =
          val<ITER_BITS>{val<ITER_BITS>{e_cur_iter} + hard<1>{}};
      val<1> at_exit = (next_iter == val<ITER_BITS>{e_nb_iter});
      val<1> pred = select(at_exit, ~val<1>{e_dir}, val<1>{e_dir});
      return concat(tag_match, confident, pred);
    };

    // Compute candidate + pred — IIFE to avoid val reassignment
    auto result = [&]() -> LookupResult {
      if constexpr (ASSOC == 1) {
        // Direct-mapped: no priority encode, minimal gate chain
        way_results.fanout(hard<3>{});
        val<1> tag_match = val<1>{way_results[0]};
        val<1> candidate = tag_match & val<1>{way_results[0] >> 1};
        val<1> pred = val<1>{way_results[0] >> 2};
        hit = tag_match;
        hit_way = val<std::max(u64(1), clog2(ASSOC))>{0};
        return {candidate, pred};
      } else {
        // Multi-way: extract + priority encode
        way_results.fanout(hard<3>{});
        arr<val<1>, ASSOC> way_tag_match = [&](u64 w) -> val<1> {
          return val<1>{way_results[w]};
        };
        arr<val<1>, ASSOC> way_candidate = [&](u64 w) -> val<1> {
          return val<1>{way_results[w]} & val<1>{way_results[w] >> 1};
        };
        arr<val<1>, ASSOC> way_pred = [&](u64 w) -> val<1> {
          return val<1>{way_results[w] >> 2};
        };

        val<ASSOC> hit_mask = way_tag_match.fo1().concat();
        val<ASSOC> candidate_mask = way_candidate.fo1().concat();
        val<ASSOC> first_candidate = candidate_mask.one_hot();
        val<ASSOC> first_hit = hit_mask.one_hot();

        hit = (hit_mask != hard<0>{});
        if constexpr (ASSOC <= 2) {
          hit_way = val<std::max(u64(1), clog2(ASSOC))>{first_hit >> 1};
        } else {
          hit_way = encode(first_hit);
        }

        return {(candidate_mask != hard<0>{}),
                (first_candidate & way_pred.fo1().concat()) != hard<0>{}};
      }
    }();

    return result;
  }

  // ======== Phase 3: train() ========

  void train(arr<val<1>, FETCH_WIDTH> &branch_taken,
             arr<val<1>, FETCH_WIDTH> &is_branch,
             val<1> &mispredict,
             val<1> &correct_pred,
             arr<reg<1>, FETCH_WIDTH> &override_active,
             [[maybe_unused]] reg<FETCH_WIDTH> &p2_before_override) {
    // Only train on offset 0 (first branch in block) for now
    val<1> actual = branch_taken[0];
    val<1> is_br = is_branch[0];

    // Re-read the matching entry from saved lookup results
    // (read_entries was set during lookup, still valid)
    read_entries.fanout(hard<2>{});
    lookup_set.fanout(hard<ASSOC + 1>{});  // one per way write + one for alloc
    lookup_tag.fanout(hard<ASSOC + 1>{});
    hit.fanout(hard<ASSOC + 2>{});
    if constexpr (ASSOC > 1) {
      hit_way.fanout(hard<ASSOC>{});
    }

    // Combined update + allocation in a single loop (one RAM write per way max)
    val<1> should_alloc = ~val<1>{hit} & mispredict & is_br &
        (val<4>{static_cast<u64>(std::rand())} == hard<0>{});  // 1/16 probability

    static_loop<ASSOC>([&]<u64 w>() {
      auto [e_tag, e_nb_iter, e_cur_iter, e_conf, e_dir, e_age] =
          split<TAG_BITS, ITER_BITS, ITER_BITS, CONF_BITS, 1, AGE_BITS>(
              read_entries[w]);

      val<1> is_hit_way = val<1>{hit} &
          (val<std::max(u64(1), clog2(ASSOC))>{hit_way} == hard<w>{});

      // --- Update path (hit way) ---
      val<1> dir_flip = val<1>{e_dir} != actual;
      val<ITER_BITS> incremented_iter =
          val<ITER_BITS>{val<ITER_BITS>{e_cur_iter} + hard<1>{}};
      val<ITER_BITS> new_cur_iter =
          select(dir_flip, val<ITER_BITS>{0}, incremented_iter);
      val<1> iter_match =
          (val<ITER_BITS>{e_cur_iter} == val<ITER_BITS>{e_nb_iter});
      val<CONF_BITS> new_conf = select(dir_flip,
          select(iter_match,
                 update_ctr(val<CONF_BITS>{e_conf}, val<1>{1}),
                 update_ctr(val<CONF_BITS>{e_conf}, val<1>{0})),
          val<CONF_BITS>{e_conf});
      val<ITER_BITS> new_nb_iter = select(
          dir_flip & (val<ITER_BITS>{e_nb_iter} == hard<0>{}),
          val<ITER_BITS>{e_cur_iter},
          val<ITER_BITS>{e_nb_iter});
      val<1> free_entry = dir_flip & ~iter_match &
          (val<ITER_BITS>{e_nb_iter} != hard<0>{});
      val<AGE_BITS> new_age = select(
          val<1>{override_active[0]} & correct_pred,
          update_ctr(val<AGE_BITS>{e_age}, val<1>{1}),
          val<AGE_BITS>{e_age});
      val<ENTRY_BITS> updated_entry = select(free_entry,
          val<ENTRY_BITS>{0},
          concat(val<TAG_BITS>{e_tag}, new_nb_iter, new_cur_iter,
                 new_conf, select(dir_flip, actual, val<1>{e_dir}),
                 new_age));

      // --- Allocation path (miss way) ---
      val<1> age_zero = (val<AGE_BITS>{e_age} == hard<0>{});
      val<1> do_alloc = should_alloc & age_zero & ~is_hit_way;
      val<ENTRY_BITS> alloc_entry = concat(
          val<TAG_BITS>{lookup_tag},
          val<ITER_BITS>{0}, val<ITER_BITS>{0}, val<CONF_BITS>{0},
          ~actual, val<AGE_BITS>{0});

      // Single write: update if hit, allocate if miss, never both
      val<ENTRY_BITS> write_entry = select(is_hit_way,
          updated_entry, alloc_entry);
      val<1> do_write = (is_hit_way | do_alloc) & is_br;

      execute_if(do_write, [&]() {
        loop_ram[w].write(val<SET_BITS>{lookup_set}, write_entry, val<1>{0});
      });
    });
  }
};
