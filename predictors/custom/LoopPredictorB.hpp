#pragma once
// LoopPredictorB — Approach B: wide entries with per-slot loop state.
// Each entry stores loop state for up to MAX_SLOTS branch offsets within a block.
// Single RAM read per lookup. Higher storage per entry but tracks multiple loops.

#include "../../harcom.hpp"
#include "../common.hpp"
#include "TageOverrider.hpp"

using namespace hcm;

template <u64 TABLE_SIZE = 64, u64 ASSOC = 4, u64 TAG_BITS = 10,
          u64 ITER_BITS = 10, u64 CONF_BITS = 3, u64 FETCH_WIDTH = 16,
          u64 CACHE_ENTRIES = 32, u64 MAX_SLOTS = 4>
struct LoopPredictorB {

  static constexpr bool ENABLED = true;
  static constexpr u64 NUM_SETS = TABLE_SIZE / ASSOC;
  static constexpr u64 SET_BITS = clog2(NUM_SETS);
  static constexpr u64 AGE_BITS = 3;
  static constexpr u64 CONF_THRESHOLD = 0;
  static constexpr u64 LOG_FW = clog2(FETCH_WIDTH);

  // Per-slot: nb_iter | cur_iter | conf | dir | offset
  static constexpr u64 SLOT_BITS = ITER_BITS + ITER_BITS + CONF_BITS + 1 + LOG_FW;
  // Entry: tag | slot[0..MAX_SLOTS-1] | age
  static constexpr u64 ENTRY_BITS = TAG_BITS + SLOT_BITS * MAX_SLOTS + AGE_BITS;

  static constexpr u64 RWRAM_BANKS = 4;
#ifdef CHEATING_MODE
  static constexpr bool USE_CACHE = (CACHE_ENTRIES > 0);
#else
  static constexpr bool USE_CACHE = false;
#endif

#ifdef TAGE_MONITOR
  struct LoopStats {
    u64 lookups = 0;
    u64 hits = 0;
    u64 overrides = 0;
    u64 alloc_writes = 0;
    u64 update_writes = 0;
    u64 cache_hits = 0;
    u64 cache_misses = 0;
    u64 cache_writes = 0;
    u64 conf_nonzero = 0;
    u64 nb_iter_nonzero = 0;
  } loop_stats;
#endif

  static_assert(NUM_SETS >= 2);
  static_assert(ASSOC >= 1);
  static_assert(TABLE_SIZE == NUM_SETS * ASSOC);

  // ======== Storage ========
  rwram<ENTRY_BITS, NUM_SETS, RWRAM_BANKS> loop_ram[ASSOC]{"loopB"};

  // SW cache
  struct EntryCache {
    u64 data[CACHE_ENTRIES] = {};
    u64 keys[CACHE_ENTRIES] = {};
    bool valids[CACHE_ENTRIES] = {};
    u64 rp = 0;

    std::pair<bool, u64> lookup(u64 set_val, u64 way_val) const {
      u64 target = (way_val << SET_BITS) | set_val;
      for (u64 i = 0; i < CACHE_ENTRIES; i++)
        if (valids[i] && keys[i] == target) return {true, data[i]};
      return {false, 0};
    }

    void update(u64 set_val, u64 way_val, u64 entry_bits) {
      u64 target = (way_val << SET_BITS) | set_val;
      for (u64 i = 0; i < CACHE_ENTRIES; i++) {
        if (valids[i] && keys[i] == target) { data[i] = entry_bits; return; }
      }
      data[rp] = entry_bits; keys[rp] = target; valids[rp] = true;
      rp = (rp + 1) % CACHE_ENTRIES;
    }
  };

  std::conditional_t<USE_CACHE, EntryCache, EmptyMember> ecache;

  // Saved lookup state
  arr<reg<ENTRY_BITS>, ASSOC> read_entries;
  reg<SET_BITS> lookup_set;
  reg<TAG_BITS> lookup_tag;
  reg<1> hit;
  reg<std::max(u64(1), clog2(ASSOC))> hit_way;
  u64 lookup_set_sw = 0;

  // Result: best confident slot's candidate, pred, and offset
  struct LookupResult {
    val<1> candidate;
    val<1> pred;
    val<LOG_FW> offset;
  };

  // ======== lookup ========
  LookupResult lookup(val<64> &inst_pc, [[maybe_unused]] auto &bindex,
                      u64 raw_pc = 0) {
    val<SET_BITS> set_idx = val<SET_BITS>{inst_pc >> 2};
    lookup_set = set_idx;
    val<TAG_BITS> tag = val<TAG_BITS>{inst_pc >> (2 + SET_BITS)};
    tag.fanout(hard<ASSOC + 1>{});
    lookup_tag = tag;
    if constexpr (ASSOC > 1) set_idx.fanout(hard<ASSOC>{});

    lookup_set_sw = (raw_pc >> 2) & ((u64(1) << SET_BITS) - 1);

    static constexpr u64 RES_BITS = 3 + LOG_FW;
    arr<val<RES_BITS>, ASSOC> way_results = [&](u64 w) -> val<RES_BITS> {
      val<ENTRY_BITS> ram_entry = loop_ram[w].read(set_idx);

      val<ENTRY_BITS> entry = [&]() -> val<ENTRY_BITS> {
        if constexpr (USE_CACHE) {
          auto [ch, raw] = ecache.lookup(lookup_set_sw, w);
          if (ch) {
#ifdef TAGE_MONITOR
            loop_stats.cache_hits++;
#endif
            return val<ENTRY_BITS>{raw};
          }
#ifdef TAGE_MONITOR
          loop_stats.cache_misses++;
#endif
        }
        return ram_entry;
      }();

      read_entries[w] = entry;
      val<TAG_BITS> e_tag = val<TAG_BITS>{entry};
      val<1> tag_match = (e_tag == tag);

      // Scan slots for best confident candidate
      val<1> any_candidate = val<1>{0};
      val<1> best_pred = val<1>{0};
      val<LOG_FW> best_offset = val<LOG_FW>{0};

      for (u64 s = 0; s < MAX_SLOTS; s++) {
        u64 sb = TAG_BITS + s * SLOT_BITS;
        val<ITER_BITS> s_nb = val<ITER_BITS>{entry >> sb};
        val<ITER_BITS> s_ci = val<ITER_BITS>{entry >> (sb + ITER_BITS)};
        val<CONF_BITS> s_cf = val<CONF_BITS>{entry >> (sb + ITER_BITS * 2)};
        val<1> s_dr = val<1>{entry >> (sb + ITER_BITS * 2 + CONF_BITS)};
        val<LOG_FW> s_off = val<LOG_FW>{entry >> (sb + ITER_BITS * 2 + CONF_BITS + 1)};

        val<1> conf_nz = (s_cf != hard<0>{});
        val<1> iter_sig = (s_nb >= hard<(u64(1) << CONF_THRESHOLD)>{});
        val<1> confident = conf_nz & iter_sig;

        val<ITER_BITS> next_iter = val<ITER_BITS>{s_ci + hard<1>{}};
        val<1> at_exit = (next_iter == s_nb);
        val<1> slot_pred = select(at_exit, ~s_dr, s_dr);

        val<1> take = confident & ~any_candidate;
        best_pred = select(take, slot_pred, best_pred);
        best_offset = select(take, s_off, best_offset);
        any_candidate = any_candidate | confident;

#ifdef TAGE_MONITOR
        if (static_cast<u64>(tag_match)) {
          if (static_cast<u64>(conf_nz)) loop_stats.conf_nonzero++;
          if (static_cast<u64>(s_nb)) loop_stats.nb_iter_nonzero++;
        }
#endif
      }

      return concat(tag_match, any_candidate, best_pred, best_offset);
    };

    // Way selection
    auto result = [&]() -> LookupResult {
      if constexpr (ASSOC == 1) {
        way_results.fanout(hard<4>{});
        val<1> tm = val<1>{way_results[0]};
        hit = tm;
        hit_way = val<std::max(u64(1), clog2(ASSOC))>{0};
#ifdef TAGE_MONITOR
        loop_stats.lookups++;
        if (static_cast<u64>(tm)) loop_stats.hits++;
        if (static_cast<u64>(tm) && static_cast<u64>(val<1>{way_results[0] >> 1}))
          loop_stats.overrides++;
#endif
        return {tm & val<1>{way_results[0] >> 1},
                val<1>{way_results[0] >> 2},
                val<LOG_FW>{way_results[0] >> 3}};
      } else {
        way_results.fanout(hard<4>{});
        arr<val<1>, ASSOC> wtm = [&](u64 w) -> val<1> { return val<1>{way_results[w]}; };
        arr<val<1>, ASSOC> wc = [&](u64 w) -> val<1> {
          return val<1>{way_results[w]} & val<1>{way_results[w] >> 1};
        };
        val<ASSOC> hm = wtm.fo1().concat();
        val<ASSOC> cm = wc.fo1().concat();
        val<ASSOC> fc = cm.one_hot();
        val<ASSOC> fh = hm.one_hot();
        hit = (hm != hard<0>{});
        if constexpr (ASSOC <= 2)
          hit_way = val<std::max(u64(1), clog2(ASSOC))>{fh >> 1};
        else
          hit_way = encode(fh);

        arr<val<1>, ASSOC> wp = [&](u64 w) -> val<1> { return val<1>{way_results[w] >> 2}; };
        arr<val<LOG_FW>, ASSOC> wo = [&](u64 w) -> val<LOG_FW> {
          return val<LOG_FW>{way_results[w] >> 3};
        };
        val<LOG_FW> oo = wo.select(hit_way);
#ifdef TAGE_MONITOR
        loop_stats.lookups++;
        if (static_cast<u64>(hm) != 0) loop_stats.hits++;
        if (static_cast<u64>(cm) != 0) loop_stats.overrides++;
#endif
        return {(cm != hard<0>{}),
                (fc & wp.fo1().concat()) != hard<0>{}, oo};
      }
    }();
    return result;
  }

  // ======== train ========
  // Updates are handled entirely in the SW cache (CHEATING_MODE).
  // The HARCOM rwram write just writes back the old entry (or alloc entry).
  void train(arr<val<1>, FETCH_WIDTH> &branch_taken,
             arr<val<1>, FETCH_WIDTH> &is_branch,
             val<1> &mispredict,
             [[maybe_unused]] val<1> &correct_pred,
             [[maybe_unused]] arr<reg<1>, FETCH_WIDTH> &override_active,
             [[maybe_unused]] reg<FETCH_WIDTH> &p2_before_override,
             val<1> &extra_cycle,
             u64 num_branch = 1) {
    read_entries.fanout(hard<2>{});
    lookup_set.fanout(hard<ASSOC + 1>{});
    lookup_tag.fanout(hard<ASSOC + 1>{});
    hit.fanout(hard<ASSOC + 2>{});
    if constexpr (ASSOC > 1) hit_way.fanout(hard<ASSOC>{});

    u64 train_offset = (num_branch > 0) ? num_branch - 1 : 0;
    val<1> actual = branch_taken[train_offset];
    val<1> is_br = is_branch[train_offset];

    val<1> should_alloc = ~val<1>{hit} & mispredict & is_br &
        (val<4>{static_cast<u64>(std::rand())} == hard<0>{});

    static_loop<ASSOC>([&]<u64 w>() {
      val<1> is_hit_way = val<1>{hit} &
          (val<std::max(u64(1), clog2(ASSOC))>{hit_way} == hard<w>{});

      val<ENTRY_BITS> old_entry = read_entries[w];
      val<AGE_BITS> e_age = val<AGE_BITS>{old_entry >> (TAG_BITS + SLOT_BITS * MAX_SLOTS)};

      val<1> age_zero = (e_age == hard<0>{});
      val<1> do_alloc = should_alloc & age_zero & ~is_hit_way;
      val<1> do_update = is_hit_way & is_br;

      // Alloc entry: tag + one slot for this branch + zeros + age=0
      val<ENTRY_BITS> alloc_entry = [&]() -> val<ENTRY_BITS> {
        u64 slot0_raw = (static_cast<u64>(~actual) & 1) << (ITER_BITS * 2 + CONF_BITS);
        slot0_raw |= (train_offset & ((u64(1) << LOG_FW) - 1)) << (ITER_BITS * 2 + CONF_BITS + 1);
        if constexpr (MAX_SLOTS == 1) {
          return concat(val<TAG_BITS>{lookup_tag}, val<SLOT_BITS>{slot0_raw}, val<AGE_BITS>{0});
        } else {
          constexpr u64 REM = SLOT_BITS * (MAX_SLOTS - 1);
          return concat(val<TAG_BITS>{lookup_tag}, val<SLOT_BITS>{slot0_raw},
                       val<REM>{0}, val<AGE_BITS>{0});
        }
      }();

      val<ENTRY_BITS> write_entry = select(do_update, old_entry, alloc_entry);
      val<1> do_write = do_update | (do_alloc & is_br);

      execute_if(do_write, [&]() {
        loop_ram[w].write(val<SET_BITS>{lookup_set}, write_entry, extra_cycle);
      });

#ifdef CHEATING_MODE
      if constexpr (USE_CACHE) {
        if (static_cast<u64>(do_write)) {
          u64 raw_entry;
          if (static_cast<u64>(do_update)) {
            raw_entry = static_cast<u64>(old_entry);
            u64 act = static_cast<u64>(actual) & 1;
            // Find matching slot and update it
            bool updated = false;
            for (u64 s = 0; s < MAX_SLOTS && !updated; s++) {
              u64 sb = TAG_BITS + s * SLOT_BITS;
              u64 s_off = (raw_entry >> (sb + ITER_BITS * 2 + CONF_BITS + 1)) &
                          ((u64(1) << LOG_FW) - 1);
              u64 s_nb = (raw_entry >> sb) & ((u64(1) << ITER_BITS) - 1);
              u64 s_ci = (raw_entry >> (sb + ITER_BITS)) & ((u64(1) << ITER_BITS) - 1);
              u64 s_cf = (raw_entry >> (sb + ITER_BITS * 2)) & ((u64(1) << CONF_BITS) - 1);
              u64 s_dr = (raw_entry >> (sb + ITER_BITS * 2 + CONF_BITS)) & 1;
              bool empty = (s_nb == 0 && s_cf == 0 && s_dr == 0 && s_off == 0);

              if (s_off == train_offset || empty) {
                bool dflip = (s_dr != act);
                u64 new_ci = dflip ? 0 : ((s_ci + 1) & ((u64(1) << ITER_BITS) - 1));
                u64 new_nb = (dflip && s_nb == 0) ? s_ci : s_nb;
                u64 new_cf = s_cf;
                if (dflip) {
                  new_cf = (s_ci == s_nb) ?
                    std::min(s_cf + 1, (u64(1) << CONF_BITS) - 1) :
                    (s_cf > 0 ? s_cf - 1 : 0);
                } else if (s_nb > 0) {
                  new_cf = std::min(s_cf + 1, (u64(1) << CONF_BITS) - 1);
                }
                u64 new_dr = dflip ? act : s_dr;
                u64 new_off = empty ? train_offset : s_off;
                u64 slot_mask = ((u64(1) << SLOT_BITS) - 1) << sb;
                raw_entry &= ~slot_mask;
                raw_entry |= ((u64)new_nb | ((u64)new_ci << ITER_BITS) |
                             ((u64)new_cf << (ITER_BITS * 2)) |
                             ((u64)new_dr << (ITER_BITS * 2 + CONF_BITS)) |
                             ((u64)new_off << (ITER_BITS * 2 + CONF_BITS + 1)))
                             << sb;
                updated = true;
              }
            }
          } else {
            raw_entry = static_cast<u64>(alloc_entry);
          }
          ecache.update(lookup_set_sw, w, raw_entry);
#ifdef TAGE_MONITOR
          loop_stats.cache_writes++;
#endif
        }
      }
#endif

#ifdef TAGE_MONITOR
      if (static_cast<u64>(do_alloc) && static_cast<u64>(is_br))
        loop_stats.alloc_writes++;
      if (static_cast<u64>(do_update))
        loop_stats.update_writes++;
#endif
    });
  }
};
