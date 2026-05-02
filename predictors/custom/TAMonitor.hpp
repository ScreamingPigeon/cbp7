#pragma once
// TageAhead Performance Monitor — software-only instrumentation
// Gated by -DTAGE_MONITOR at compile time.
// All counters are plain C++ types (zero HARCOM cost).

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using u64 = uint64_t;
using i64 = int64_t;

template <u64 NUM_TABLES, u64 N, u64 MAX_TABLE_ENTRIES = 2048,
          bool USE_GSHARE = false, u64 FB_CAPACITY = 8192>
struct TAMonitor {

  static constexpr const char *FB_NAME = USE_GSHARE ? "Gshare" : "Bimodal";
  static constexpr const char *FB_NAME_PAD = USE_GSHARE ? "Gshare  " : "Bimodal ";

  static constexpr u64 WINDOW_SIZE = 100000; // branches per window

  static constexpr u64 MAX_BLOCK_INSTR = 64;
  static constexpr u64 MAX_BLOCK_BR = N + 1;
  static constexpr u64 MAX_RANK = N;

  // Allocation trajectory bins: 0-10, 10-50, 50-100, 100-300, 300+
  static constexpr u64 TRAJ_BINS = 5;
  static u64 traj_bin(u64 alloc_count) {
    if (alloc_count < 10) return 0;
    if (alloc_count < 50) return 1;
    if (alloc_count < 100) return 2;
    if (alloc_count < 300) return 3;
    return 4;
  }
  static const char *traj_label(u64 b) {
    static const char *labels[] = {"1-10", "10-50", "50-100", "100-300", "300+"};
    return labels[b];
  }

  // Entry age bins: 0-1K, 1K-10K, 10K-100K, 100K-1M, 1M+
  static constexpr u64 AGE_BINS = 5;
  static u64 age_bin(u64 age) {
    if (age < 1000) return 0;
    if (age < 10000) return 1;
    if (age < 100000) return 2;
    if (age < 1000000) return 3;
    return 4;
  }
  static const char *age_label(u64 b) {
    static const char *labels[] = {"0-1K", "1K-10K", "10K-100K", "100K-1M", "1M+"};
    return labels[b];
  }

  // Recurrence distance bins: 1-10, 10-100, 100-1K, 1K-10K, 10K+
  static constexpr u64 RECUR_BINS = 5;
  static u64 recur_bin(u64 dist) {
    if (dist < 10) return 0;
    if (dist < 100) return 1;
    if (dist < 1000) return 2;
    if (dist < 10000) return 3;
    return 4;
  }
  static const char *recur_label(u64 b) {
    static const char *labels[] = {"1-10", "10-100", "100-1K", "1K-10K", "10K+"};
    return labels[b];
  }

  // Alloc histogram bins: 1-10, 10-50, 50-100, 100-300, 300+
  static constexpr u64 ALLOC_HIST_BINS = 5;
  static u64 alloc_hist_bin(u64 ac) {
    if (ac < 10) return 0;
    if (ac < 50) return 1;
    if (ac < 100) return 2;
    if (ac < 300) return 3;
    return 4;
  }
  static const char *alloc_hist_label(u64 b) {
    static const char *labels[] = {"1-10", "10-50", "50-100", "100-300", "300+"};
    return labels[b];
  }

  // Micro-window MPKI bins: 0-5, 5-10, ..., 45-50, 50+
  static constexpr u64 MICRO_MPKI_BINS = 11;
  static constexpr u64 MICRO_WIN = 1024;

  // Feature 18: Phase detection via short-window misp rate delta
  // Two consecutive windows of PHASE_WIN branches. The absolute difference
  // in misprediction rate between them is the phase-change signal.
  static constexpr u64 PHASE_WIN = 256;
  // Delta-rate bins: [0,1%), [1,2%), [2,5%), [5,10%), [10,20%), 20%+
  static constexpr u64 PHASE_BINS = 6;
  static u64 phase_bin(double delta_pct) {
    if (delta_pct < 1.0) return 0;
    if (delta_pct < 2.0) return 1;
    if (delta_pct < 5.0) return 2;
    if (delta_pct < 10.0) return 3;
    if (delta_pct < 20.0) return 4;
    return 5;
  }
  static const char *phase_label(u64 b) {
    static const char *labels[] = {"<1%", "1-2%", "2-5%", "5-10%", "10-20%", "20%+"};
    return labels[b];
  }

  // Feature 5: Lifetime stats
  struct LifetimeStats {
    u64 evictions = 0;
    u64 total_preds = 0;
    u64 total_correct = 0;
    u64 zero_use = 0;
  };

  // Thrashing stats per table
  struct ThrashStats {
    static constexpr u64 LIFE_BUCKETS = 5;
    std::array<u64, LIFE_BUCKETS> life_hist{};
    std::array<u64, LIFE_BUCKETS> life_useful{};
    u64 useful_evict = 0;
    u64 useful_evict_preds = 0;
    u64 pingpong = 0;
    u64 repeat_alloc = 0;
    static u64 life_bucket(u64 lifetime) {
      if (lifetime < 10) return 0;
      if (lifetime < 100) return 1;
      if (lifetime < 1000) return 2;
      if (lifetime < 10000) return 3;
      return 4;
    }
    static const char *life_label(u64 b) {
      static const char *labels[] = {"<10", "10-100", "100-1K", "1K-10K", "10K+"};
      return labels[b];
    }
  };

  // Per-PC diagnostics
  struct PCDiag {
    u64 total = 0;
    u64 misp = 0;
    u64 taken = 0;
    u64 tage_prov = 0;
    u64 tage_correct = 0;
    u64 fb_prov = 0;
    u64 fb_correct = 0;
    u64 alloc_count = 0;
    // Feature 12: per block-position breakdown
    std::array<u64, MAX_RANK> pos_total{};
    std::array<u64, MAX_RANK> pos_misp{};
    // Feature 14: TAGE accuracy trajectory by alloc count
    std::array<u64, TRAJ_BINS> traj_tage_prov{};
    std::array<u64, TRAJ_BINS> traj_tage_correct{};
    // Feature 15: oscillation detection (consecutive correct streaks)
    u64 streak = 0;        // current consecutive correct count
    u64 streak_count = 0;  // number of completed streaks
    u64 streak_sum = 0;    // sum of streak lengths (for avg)
    u64 max_streak = 0;
    // Feature 16/24: counterfactual FB when TAGE provides
    u64 cf_tage_provides = 0;      // times TAGE was provider
    u64 cf_fb_would_correct = 0;   // of those, times FB would have been correct
    u64 cf_tage_correct_here = 0;  // of those, times TAGE was correct
  };

  // ======== Counters (cumulative + per-window) ========
  struct Counters {
    u64 branches = 0;
    u64 mispredictions = 0;
    u64 blocks = 0;
    u64 extra_cycles = 0;

    // Feature 10: misprediction type breakdown
    u64 misp_direction = 0;
    u64 misp_boundary = 0;

    // Pipeline call counts
    u64 predict1_calls = 0;
    u64 reuse_predict1_calls = 0;
    u64 predict2_calls = 0;
    u64 reuse_predict2_calls = 0;

    // Block structure
    u64 total_block_instr = 0;
    u64 total_block_branches = 0;
    std::array<u64, MAX_BLOCK_INSTR> block_instr_hist{};
    std::array<u64, MAX_BLOCK_BR> block_br_hist{};
    std::array<u64, 1024> entry_point_hist{};
    std::array<u64, 1024> exit_point_hist{};

    u64 true_block_count = 0;
    u64 train_skip_count = 0;

    // Provider distribution
    std::array<u64, NUM_TABLES + 1> provider_count{};
    std::array<u64, NUM_TABLES + 1> provider_correct{};
    std::array<u64, NUM_TABLES + 1> alt_count{};

    // Meta override
    u64 meta_override_count = 0;
    u64 meta_override_correct = 0;
    u64 meta_chose_alt = 0;
    u64 meta_chose_alt_correct = 0;
    u64 meta_chose_pri = 0;
    u64 meta_chose_pri_correct = 0;

    // Provider source breakdown (per branch):
    //   no_tag_match:   no TAGE table had a tag match at all
    //   sec_tag_reject: at least one tag match, but sec_tag filtered all out
    //   tage_meta_alt:  TAGE provider exists, meta chose alt (fallback)
    //   tage_meta_pri:  TAGE provider exists, meta chose primary (TAGE)
    //   tage_no_meta:   TAGE provider exists, meta not involved (provider strong)
    u64 prov_no_tag_match = 0;
    u64 prov_no_tag_match_correct = 0;
    u64 prov_sec_tag_reject = 0;
    u64 prov_sec_tag_reject_correct = 0;
    u64 prov_tage_meta_alt = 0;
    u64 prov_tage_meta_alt_correct = 0;
    u64 prov_tage_meta_pri = 0;
    u64 prov_tage_meta_pri_correct = 0;
    u64 prov_tage_no_meta = 0;
    u64 prov_tage_no_meta_correct = 0;

    // Counterfactual: for sec-tag-rejected branches, would the tag-only
    // provider have been correct? (i.e. what we lose by sec-tag filtering)
    u64 cf_sec_tag_only_total = 0;
    u64 cf_sec_tag_only_correct = 0;  // tag-only provider would have been right
    u64 cf_sec_fb_correct = 0;        // fallback was right (what we actually used)

    // Tag match rate
    std::array<u64, NUM_TABLES> tag_lookups{};
    std::array<u64, NUM_TABLES> tag_matches{};
    std::array<u64, NUM_TABLES> sec_matches{};
    std::array<u64, NUM_TABLES> full_matches{};

    // Allocation
    u64 alloc_attempts = 0;
    u64 alloc_success = 0;
    u64 alloc_fail = 0;
    u64 alloc_blocked = 0;
    u64 alloc_sibling_skip = 0;
    std::array<u64, NUM_TABLES> alloc_sibling_per_table{};
    std::array<u64, NUM_TABLES> alloc_per_table{};
    std::array<u64, NUM_TABLES + 1> alloc_from_provider{};
    std::array<std::array<u64, NUM_TABLES>, NUM_TABLES + 1> alloc_cascade{};

    // u-bit
    std::array<u64, NUM_TABLES> u_set_count{};
    std::array<u64, NUM_TABLES> u_clear_count{};
    u64 decay_fire_count = 0;
    u64 epoch_reset_count = 0;
    std::array<u64, NUM_TABLES> uclear_alloc_fail{};
    std::array<u64, NUM_TABLES> uclear_sibling{};
    std::array<u64, NUM_TABLES> uclear_epoch{};

    // Per-rank accuracy
    std::array<u64, N> rank_total{};
    std::array<u64, N> rank_correct{};

    // Pressure
    u64 acc_ctr_sum = 0;
    u64 alloc_ctr_sum = 0;
    u64 pressure_samples = 0;

    // Feature 5: Lifetime stats
    std::array<LifetimeStats, NUM_TABLES> lifetime_stats{};
    std::array<ThrashStats, NUM_TABLES> thrash{};

    // Feature 7: Tag collision
    u64 collision_checks = 0;
    u64 collision_hits = 0;
    std::array<u64, NUM_TABLES> per_table_coll_checks{};
    std::array<u64, NUM_TABLES> per_table_coll_hits{};

    // sec_tag consistency
    u64 sec_tag_checks = 0;
    u64 sec_tag_match_curr = 0;
    u64 sec_tag_match_now = 0;
    u64 sec_tag_match_both = 0;
    u64 sec_tag_match_neither = 0;

    // Feature 17: micro-window MPKI histogram
    std::array<u64, MICRO_MPKI_BINS> micro_mpki_hist{};

    // Feature 18: phase-change detection (short-window misp rate delta)
    std::array<u64, PHASE_BINS> phase_delta_hist{};
    u64 phase_samples = 0;
    double phase_delta_sum = 0;        // sum of |delta_rate|
    double phase_max_delta = 0;        // largest single delta seen

    // Feature 20: entry age at provider
    std::array<u64, AGE_BINS> age_prov_count{};
    std::array<u64, AGE_BINS> age_prov_correct{};

    // Feature 21: per-PC recurrence distance
    std::array<u64, RECUR_BINS> recurrence_hist{};

    // Feature 22: ping-pong evictions
    u64 pingpong_count = 0;

    // Feature 23: sec-tag mismatch breakdown
    u64 sec_mm_same_pc = 0;   // same PC, different path
    u64 sec_mm_diff_pc = 0;   // different PC entirely
    u64 sec_mm_invalid = 0;   // no valid entry

    // Feature 24: counterfactual fallback
    u64 cf_total = 0;         // all branches
    u64 cf_fb_correct = 0;    // FB would be correct
    u64 cf_tage_correct = 0;  // TAGE/actual was correct
    u64 cf_fb_only = 0;       // FB right, TAGE wrong
    u64 cf_tage_only = 0;     // TAGE right, FB wrong

    // Feature 25: sec-tag adaptive benefit counter diagnostics
    u64 ben_reject_all = 0;   // blocks where sec_would_reject_all fired (causal condition)
    u64 ben_incr = 0;         // benefit counter increments (FB > TAGE, enforcement helped)
    u64 ben_decr = 0;         // benefit counter decrements (TAGE > FB, enforcement hurt)
    u64 ben_tie = 0;          // causal condition fired but outcome tied (no update)
    u64 ben_ctr_sum = 0;      // sum of benefit_ctr value (for avg)
    u64 ben_ctr_samples = 0;  // number of samples

    void reset() { *this = Counters{}; }
  };

  Counters cum;
  Counters win;
  u64 window_num = 0;
  bool header_printed = false;

  // Per-table sizes
  std::array<u64, NUM_TABLES> table_sizes{};

  // Table occupancy
  std::array<std::bitset<MAX_TABLE_ENTRIES>, NUM_TABLES> tage_occupied{};
  std::array<u64, NUM_TABLES> tage_unique_entries{};
  std::bitset<FB_CAPACITY> fb_occupied{};
  u64 fb_unique_entries = 0;
  u64 fb_writes = 0;

  std::unordered_set<u64> unique_branch_pcs;

  // Shadow state
  std::array<u64, N> shadow_provider{};
  std::array<u64, N> shadow_alt{};
  std::array<bool, N> shadow_meta_overrode{};
  std::array<bool, N> shadow_meta_chose_alt{};
  std::array<bool, N> shadow_pred{};
  std::array<u64, N> shadow_pc{};
  // Feature 24: shadow fallback prediction per branch
  std::array<bool, N> shadow_fb_pred{};
  // Provider source breakdown shadow
  bool shadow_any_tag_hit = false;
  bool shadow_has_tage_provider = false;
  std::array<bool, N> shadow_tag_only_pred{};

  // Per-PC diagnostics
  std::unordered_map<u64, PCDiag> pc_diag;
  std::unordered_map<u64, PCDiag> win_pc_diag;

  // Feature 11: Per-PC ahead-context fanout
  std::unordered_map<u64, std::unordered_set<u64>> pc_block_contexts;

  // Feature 5: Entry lifetime state
  struct EntryLife {
    u64 alloc_time = 0;
    u64 pred_count = 0;
    u64 correct_count = 0;
    u64 stored_pc = 0;
    bool valid = false;
  };
  std::array<std::array<EntryLife, MAX_TABLE_ENTRIES>, NUM_TABLES> entry_life{};
  u64 global_branch_count = 0;

  // Feature 9: Entry utilization
  std::array<std::bitset<MAX_TABLE_ENTRIES>, NUM_TABLES> entry_ever_hit{};
  std::array<std::bitset<MAX_TABLE_ENTRIES>, NUM_TABLES> win_entry_ever_hit{};

  // Block PC pipeline
  u64 shadow_block_pc = 0;
  u64 current_block_pc = 0;
  u64 train_block_pc = 0;

  std::array<u64, N> current_pcs{};
  u64 current_num_branch = 0;
  std::array<u64, N> train_pcs{};
  u64 train_num_branch = 0;

  // Feature 17: micro sliding window for fine-grained MPKI
  std::array<uint8_t, MICRO_WIN> micro_buffer{}; // 1=mispredict, 0=correct
  u64 micro_idx = 0;
  u64 micro_misp_count = 0;
  bool micro_filled = false;

  // Feature 18: phase detection — non-overlapping short windows
  u64 phase_cur_misp = 0;   // misps in current PHASE_WIN window
  u64 phase_prev_misp = 0;  // misps in previous PHASE_WIN window
  u64 phase_win_pos = 0;    // position within current PHASE_WIN window
  bool phase_has_prev = false;

  // Feature 19: previous window PC set for Jaccard
  std::unordered_set<u64> prev_window_pcs;

  // Feature 21: per-PC last misprediction time
  std::unordered_map<u64, u64> pc_last_misp_time;

  // Feature 22: ping-pong detection — previous evicted PC per entry
  std::array<std::array<u64, MAX_TABLE_ENTRIES>, NUM_TABLES> prev_evicted_pc{};

  // ======== Recording methods ========

  void record_predict1() { cum.predict1_calls++; win.predict1_calls++; }
  void record_reuse_predict1() {
    cum.reuse_predict1_calls++;
    win.reuse_predict1_calls++;
  }
  void record_predict2() { cum.predict2_calls++; win.predict2_calls++; }
  void record_reuse_predict2() {
    cum.reuse_predict2_calls++;
    win.reuse_predict2_calls++;
  }

  void record_block(u64 block_entry, u64 block_size, u64 num_branch,
                     bool extra_cycle_fired) {
    auto record = [&](Counters &c) {
      c.blocks++;
      c.total_block_instr += block_size;
      c.total_block_branches += num_branch;
      if (extra_cycle_fired) c.extra_cycles++;
      if (block_size < MAX_BLOCK_INSTR)
        c.block_instr_hist[block_size]++;
      else
        c.block_instr_hist[MAX_BLOCK_INSTR - 1]++;
      if (num_branch < MAX_BLOCK_BR)
        c.block_br_hist[num_branch]++;
      if (block_entry < 1024)
        c.entry_point_hist[block_entry]++;
      u64 exit_pt = block_entry + block_size;
      if (exit_pt < 1024)
        c.exit_point_hist[exit_pt]++;
    };
    record(cum);
    record(win);
  }

  void record_true_block() {
    cum.true_block_count++;
    win.true_block_count++;
  }

  void record_train_skip() {
    cum.train_skip_count++;
    win.train_skip_count++;
  }

  void record_tag_lookup(u64 table, bool tag_hit, bool sec_hit) {
    auto record = [&](Counters &c) {
      c.tag_lookups[table]++;
      if (tag_hit) c.tag_matches[table]++;
      if (sec_hit) c.sec_matches[table]++;
      if (tag_hit && sec_hit) c.full_matches[table]++;
    };
    record(cum);
    record(win);
  }

  void record_prediction(u64 rank, u64 match1_val, u64 match2_val,
                          bool meta_overrode, bool meta_chose_alt,
                          bool pred_taken,
                          bool any_tag_hit = false,
                          bool has_tage_provider = false,
                          bool tag_only_pred = false) {
    u64 prov = decode_provider(match1_val);
    shadow_provider[rank] = prov;
    shadow_alt[rank] = decode_provider(match2_val);
    shadow_meta_overrode[rank] = meta_overrode;
    shadow_meta_chose_alt[rank] = meta_chose_alt;
    shadow_pred[rank] = pred_taken;
    shadow_any_tag_hit = any_tag_hit;
    shadow_has_tage_provider = has_tage_provider;
    shadow_tag_only_pred[rank] = tag_only_pred;
  }

  // Feature 24: record fallback prediction per branch (call from TageAhead)
  void record_fb_prediction(u64 rank, bool fb_taken) {
    if (rank < N) shadow_fb_pred[rank] = fb_taken;
  }

  void record_outcome(u64 rank, bool actual_taken, bool mispredict) {
    bool correct = (shadow_pred[rank] == actual_taken);
    bool fb_correct = (shadow_fb_pred[rank] == actual_taken);

    auto record = [&](Counters &c) {
      c.branches++;
      if (mispredict) {
        c.mispredictions++;
        if (!correct)
          c.misp_direction++;
        else
          c.misp_boundary++;
      }

      u64 prov = shadow_provider[rank];
      c.provider_count[prov]++;
      if (correct) c.provider_correct[prov]++;

      u64 alt = shadow_alt[rank];
      c.alt_count[alt]++;

      if (shadow_meta_overrode[rank]) {
        c.meta_override_count++;
        if (correct) c.meta_override_correct++;
        if (shadow_meta_chose_alt[rank]) {
          c.meta_chose_alt++;
          if (correct) c.meta_chose_alt_correct++;
        } else {
          c.meta_chose_pri++;
          if (correct) c.meta_chose_pri_correct++;
        }
      }

      // Provider source breakdown
      if (!shadow_any_tag_hit) {
        c.prov_no_tag_match++;
        if (correct) c.prov_no_tag_match_correct++;
      } else if (!shadow_has_tage_provider) {
        c.prov_sec_tag_reject++;
        if (correct) c.prov_sec_tag_reject_correct++;
        // Counterfactual: would tag-only provider have been correct?
        c.cf_sec_tag_only_total++;
        bool tag_only_correct = (shadow_tag_only_pred[rank] == actual_taken);
        if (tag_only_correct) c.cf_sec_tag_only_correct++;
        if (correct) c.cf_sec_fb_correct++; // FB was actually used and correct
      } else if (shadow_meta_overrode[rank]) {
        if (shadow_meta_chose_alt[rank]) {
          c.prov_tage_meta_alt++;
          if (correct) c.prov_tage_meta_alt_correct++;
        } else {
          c.prov_tage_meta_pri++;
          if (correct) c.prov_tage_meta_pri_correct++;
        }
      } else {
        c.prov_tage_no_meta++;
        if (correct) c.prov_tage_no_meta_correct++;
      }

      if (rank < N) {
        c.rank_total[rank]++;
        if (correct) c.rank_correct[rank]++;
      }

      // Feature 24: counterfactual fallback
      c.cf_total++;
      if (fb_correct) c.cf_fb_correct++;
      if (correct) c.cf_tage_correct++;
      if (fb_correct && !correct) c.cf_fb_only++;
      if (!fb_correct && correct) c.cf_tage_only++;

      // Feature 21: recurrence distance
      if (mispredict) {
        u64 pc = shadow_pc[rank];
        auto it = pc_last_misp_time.find(pc);
        if (it != pc_last_misp_time.end()) {
          u64 dist = global_branch_count - it->second;
          c.recurrence_hist[recur_bin(dist)]++;
        }
      }
    };
    record(cum);
    record(win);

    // Feature 21: update last misp time (non-windowed state)
    u64 pc = shadow_pc[rank];
    if (mispredict)
      pc_last_misp_time[pc] = global_branch_count;

    u64 prov = shadow_provider[rank];
    auto update_pc_diag = [&](PCDiag &d) {
      d.total++;
      if (mispredict) d.misp++;
      if (actual_taken) d.taken++;
      if (prov < NUM_TABLES) {
        d.tage_prov++;
        if (correct) d.tage_correct++;
        // Feature 14: trajectory — bin by current alloc count
        u64 tb = traj_bin(d.alloc_count);
        d.traj_tage_prov[tb]++;
        if (correct) d.traj_tage_correct[tb]++;
      } else {
        d.fb_prov++;
        if (correct) d.fb_correct++;
      }
      // Feature 12: per-position tracking
      if (rank < MAX_RANK) {
        d.pos_total[rank]++;
        if (mispredict) d.pos_misp[rank]++;
      }
      // Feature 15: oscillation — consecutive correct streaks
      if (correct) {
        d.streak++;
      } else {
        if (d.streak > 0) {
          d.streak_count++;
          d.streak_sum += d.streak;
          if (d.streak > d.max_streak) d.max_streak = d.streak;
        }
        d.streak = 0;
      }
      // Feature 16/24: counterfactual per-PC (when TAGE provides)
      if (prov < NUM_TABLES) {
        d.cf_tage_provides++;
        if (fb_correct) d.cf_fb_would_correct++;
        if (correct) d.cf_tage_correct_here++;
      }
    };
    update_pc_diag(pc_diag[pc]);
    update_pc_diag(win_pc_diag[pc]);

    // Feature 11: track block context for mispredicted PCs
    if (mispredict)
      pc_block_contexts[pc].insert(current_block_pc);

    // Feature 17: micro sliding window
    {
      if (micro_filled) {
        // Remove outgoing sample
        micro_misp_count -= micro_buffer[micro_idx];
      }
      micro_buffer[micro_idx] = mispredict ? 1 : 0;
      micro_misp_count += mispredict ? 1 : 0;
      micro_idx = (micro_idx + 1) % MICRO_WIN;
      if (!micro_filled && micro_idx == 0) micro_filled = true;
      if (micro_filled) {
        // Compute MPKI approximation: misps per 1K branches
        double micro_mpki = 1000.0 * micro_misp_count / MICRO_WIN;
        u64 bin = std::min(u64(micro_mpki / 5.0), MICRO_MPKI_BINS - 1);
        cum.micro_mpki_hist[bin]++;
        win.micro_mpki_hist[bin]++;
      }
    }

    // Feature 18: phase detection — non-overlapping PHASE_WIN-branch windows
    {
      phase_cur_misp += mispredict ? 1 : 0;
      phase_win_pos++;
      if (phase_win_pos >= PHASE_WIN) {
        if (phase_has_prev) {
          double cur_rate = 100.0 * phase_cur_misp / PHASE_WIN;
          double prev_rate = 100.0 * phase_prev_misp / PHASE_WIN;
          double delta = std::abs(cur_rate - prev_rate);
          u64 b = phase_bin(delta);
          cum.phase_delta_hist[b]++;
          win.phase_delta_hist[b]++;
          cum.phase_samples++;
          win.phase_samples++;
          cum.phase_delta_sum += delta;
          win.phase_delta_sum += delta;
          if (delta > cum.phase_max_delta) cum.phase_max_delta = delta;
          if (delta > win.phase_max_delta) win.phase_max_delta = delta;
        }
        phase_prev_misp = phase_cur_misp;
        phase_cur_misp = 0;
        phase_win_pos = 0;
        phase_has_prev = true;
      }
    }

    global_branch_count++;

    if (win.branches >= WINDOW_SIZE) {
      print_window(std::cerr);
      // Feature 19: save current window PCs for Jaccard
      prev_window_pcs.clear();
      for (auto &[wpc, wd] : win_pc_diag)
        prev_window_pcs.insert(wpc);
      win.reset();
      win_pc_diag.clear();
      for (auto &b : win_entry_ever_hit) b.reset();
      window_num++;
    }
  }

  void record_branch_pc(u64 pc) { unique_branch_pcs.insert(pc); }

  void record_branch_info(u64 rank, u64 pc, [[maybe_unused]] bool taken) {
    shadow_pc[rank] = pc;
    unique_branch_pcs.insert(pc);
  }

  // Allocation
  void record_allocation(bool success, u64 allocate_mask) {
    auto record = [&](Counters &c) {
      c.alloc_attempts++;
      if (success) {
        c.alloc_success++;
        for (u64 i = 0; i < NUM_TABLES; i++)
          if (allocate_mask & (u64(1) << i)) c.alloc_per_table[i]++;
      } else {
        c.alloc_fail++;
      }
    };
    record(cum);
    record(win);
  }

  void record_alloc_blocked() {
    cum.alloc_blocked++;
    win.alloc_blocked++;
  }

  void record_alloc_sibling_skip(u64 sibling_mask) {
    auto record = [&](Counters &c) {
      for (u64 i = 0; i < NUM_TABLES; i++)
        if (sibling_mask & (u64(1) << i)) {
          c.alloc_sibling_skip++;
          c.alloc_sibling_per_table[i]++;
        }
    };
    record(cum);
    record(win);
  }

  void record_alloc_cascade(u64 provider_idx, u64 allocate_mask) {
    auto record = [&](Counters &c) {
      if (allocate_mask == 0) return;
      c.alloc_from_provider[provider_idx]++;
      for (u64 i = 0; i < NUM_TABLES; i++)
        if (allocate_mask & (u64(1) << i))
          c.alloc_cascade[provider_idx][i]++;
    };
    record(cum);
    record(win);
  }

  u64 mask_index(u64 table, u64 index) const {
    return table_sizes[table] > 0 ? (index & (table_sizes[table] - 1)) : index;
  }

  void record_tage_write(u64 table, u64 index) {
    if (table >= NUM_TABLES) return;
    index = mask_index(table, index);
    if (index < MAX_TABLE_ENTRIES && !tage_occupied[table].test(index)) {
      tage_occupied[table].set(index);
      tage_unique_entries[table]++;
    }
  }

  void record_fb_write(u64 index) {
    fb_writes++;
    if (index < FB_CAPACITY && !fb_occupied.test(index)) {
      fb_occupied.set(index);
      fb_unique_entries++;
    }
  }

  void shift_block_pc(u64 num_branch) {
    train_block_pc = current_block_pc;
    train_pcs = current_pcs;
    train_num_branch = current_num_branch;
    current_block_pc = shadow_block_pc;
    for (u64 i = 0; i < num_branch && i < N; i++)
      current_pcs[i] = shadow_pc[i];
    current_num_branch = num_branch;
  }

  // Feature 5+20: Record provider entry usage + age
  void record_provider_entry(u64 table, u64 index, u64 num_branches,
                              u64 num_correct) {
    if (table >= NUM_TABLES) return;
    index = mask_index(table, index);
    if (index < MAX_TABLE_ENTRIES) {
      auto &e = entry_life[table][index];
      if (e.valid) {
        e.pred_count += num_branches;
        e.correct_count += num_correct;
        // Feature 20: entry age at provider
        u64 age = global_branch_count - e.alloc_time;
        u64 ab = age_bin(age);
        bool block_correct = (num_correct == num_branches);
        auto record_age = [&](Counters &c) {
          c.age_prov_count[ab]++;
          if (block_correct) c.age_prov_correct[ab]++;
        };
        record_age(cum);
        record_age(win);
      }
    }
  }

  // Feature 7+9: Tag collision check
  void record_collision_check(u64 table, u64 index, bool tag_hit,
                               u64 lookup_pc) {
    if (!tag_hit || table >= NUM_TABLES) return;
    index = mask_index(table, index);
    if (index >= MAX_TABLE_ENTRIES) return;
    entry_ever_hit[table].set(index);
    win_entry_ever_hit[table].set(index);
    auto record_coll = [&](Counters &c) {
      c.per_table_coll_checks[table]++;
      c.collision_checks++;
      auto &e = entry_life[table][index];
      if (e.valid && e.stored_pc != 0 && e.stored_pc != lookup_pc) {
        c.per_table_coll_hits[table]++;
        c.collision_hits++;
      }
    };
    record_coll(cum);
    record_coll(win);
  }

  // Feature 5+7+8+22: Record allocation (per table on write)
  void record_entry_alloc_diag(u64 table, u64 index) {
    if (table >= NUM_TABLES) return;
    index = mask_index(table, index);
    if (index >= MAX_TABLE_ENTRIES) return;
    auto &e = entry_life[table][index];
    // Eviction stats for old entry
    if (e.valid) {
      // Feature 22: ping-pong detection
      u64 evicted_pc = e.stored_pc;
      u64 prev_evict = prev_evicted_pc[table][index];
      if (prev_evict != 0 && prev_evict == train_block_pc && evicted_pc != 0) {
        cum.pingpong_count++;
        win.pingpong_count++;
      }
      prev_evicted_pc[table][index] = evicted_pc;

      auto record_evict = [&](Counters &c) {
        auto &ls = c.lifetime_stats[table];
        ls.evictions++;
        ls.total_preds += e.pred_count;
        ls.total_correct += e.correct_count;
        if (e.pred_count == 0) ls.zero_use++;
      };
      record_evict(cum);
      record_evict(win);
    }
    // Initialize new entry
    e.alloc_time = global_branch_count;
    e.pred_count = 0;
    e.correct_count = 0;
    e.stored_pc = train_block_pc;
    e.valid = true;
    // Feature 8: per-PC alloc tracking
    for (u64 i = 0; i < train_num_branch && i < N; i++) {
      pc_diag[train_pcs[i]].alloc_count++;
      win_pc_diag[train_pcs[i]].alloc_count++;
    }
  }

  void record_u_write(u64 table, bool new_u) {
    if (new_u) {
      cum.u_set_count[table]++;
      win.u_set_count[table]++;
    } else {
      cum.u_clear_count[table]++;
      win.u_clear_count[table]++;
    }
  }

  void record_decay_fire() {
    cum.decay_fire_count++;
    win.decay_fire_count++;
  }
  void record_epoch_reset() {
    cum.epoch_reset_count++;
    win.epoch_reset_count++;
  }

  void record_uclear_source(u64 table, u64 source) {
    auto record = [&](Counters &c) {
      if (source == 0) c.uclear_alloc_fail[table]++;
      else if (source == 1) c.uclear_sibling[table]++;
      else c.uclear_epoch[table]++;
    };
    record(cum);
    record(win);
  }

  void record_sec_tag_check(u64 stored, u64 curr, u64 now) {
    auto record = [&](Counters &c) {
      c.sec_tag_checks++;
      bool m_curr = (stored == curr);
      bool m_now = (stored == now);
      if (m_curr) c.sec_tag_match_curr++;
      if (m_now) c.sec_tag_match_now++;
      if (m_curr && m_now) c.sec_tag_match_both++;
      if (!m_curr && !m_now) c.sec_tag_match_neither++;
    };
    record(cum);
    record(win);
  }

  // Feature 23: sec-tag mismatch source breakdown
  // Call when sec_tag doesn't match (match-neither case).
  // stored_entry_pc = the PC stored in the entry, lookup_pc = current block PC.
  void record_sec_mismatch_source(u64 table, u64 index, u64 lookup_pc) {
    if (table >= NUM_TABLES) return;
    index = mask_index(table, index);
    if (index >= MAX_TABLE_ENTRIES) return;
    auto &e = entry_life[table][index];
    auto record = [&](Counters &c) {
      if (!e.valid) {
        c.sec_mm_invalid++;
      } else if (e.stored_pc == lookup_pc) {
        c.sec_mm_same_pc++;   // same PC, different ahead path
      } else {
        c.sec_mm_diff_pc++;   // different PC entirely (aliasing)
      }
    };
    record(cum);
    record(win);
  }

  void record_pressure(u64 acc_val, u64 alloc_val) {
    auto record = [&](Counters &c) {
      c.acc_ctr_sum += acc_val;
      c.alloc_ctr_sum += alloc_val;
      c.pressure_samples++;
    };
    record(cum);
    record(win);
  }

  // Feature 25: benefit counter update diagnostics.
  // fired=true when sec_would_reject_all was true (causal condition met).
  // updated=true when fb_right != tage_right (outcome was decisive).
  // incr=true when benefit counter was incremented (FB > TAGE).
  void record_benefit_update(bool fired, bool updated, bool incr) {
    auto record = [&](Counters &c) {
      if (fired) {
        c.ben_reject_all++;
        if (updated) {
          if (incr) c.ben_incr++;
          else      c.ben_decr++;
        } else {
          c.ben_tie++;
        }
      }
    };
    record(cum);
    record(win);
  }

  void record_benefit_ctr(u64 val) {
    auto record = [&](Counters &c) {
      c.ben_ctr_sum += val;
      c.ben_ctr_samples++;
    };
    record(cum);
    record(win);
  }

  // ======== Helpers ========

  static u64 decode_provider(u64 one_hot) {
    if (one_hot == 0) return NUM_TABLES;
    for (u64 i = 0; i <= NUM_TABLES; i++)
      if (one_hot & (u64(1) << i)) return i;
    return NUM_TABLES;
  }

  static double pct(u64 num, u64 den) {
    return den > 0 ? 100.0 * num / den : 0.0;
  }

  // ======== Window CSV output ========

  void print_window(std::ostream &os) {
    if (!header_printed) {
      os << "# win,br,misp%,MPKI,extra%,i/blk,br/blk,true_blk%,";
      os << (USE_GSHARE ? "gs%," : "bim%,");
      for (u64 i = 0; i < NUM_TABLES; i++)
        os << "T" << i << "%,";
      os << "alloc_ok%,acc_avg,alloc_avg,";
      os << "Tlast_evict,Tlast_zu%,Tlast_avgp,";
      os << "coll%,Tlast_used,win_stuck,win_hard,";
      os << "dir_misp%,bnd_misp%,";
      // New feature columns
      os << "micro_p50_mpki,micro_p95_mpki,";
      os << "phase_delta_avg,phase_max_delta,";
      os << "jaccard,";
      os << "pingpong,";
      os << "cf_fb_only%,cf_tage_only%,";
      os << "prov_no_tag%,prov_sec_rej%,prov_meta_alt%,prov_meta_pri%,prov_no_meta%,";
      os << "cf_sec_fb_acc%,cf_sec_tage_acc%,";
      os << "ben_rej%,ben_incr%,ben_decr%,ben_ctr_avg";
      os << "\n";
      header_printed = true;
    }
    auto &w = win;
    double win_mpki = w.total_block_instr > 0
                          ? 1000.0 * w.mispredictions / w.total_block_instr
                          : 0.0;
    os << std::fixed << std::setprecision(1);
    os << window_num << "," << w.branches << ","
       << pct(w.mispredictions, w.branches) << ","
       << std::setprecision(3) << win_mpki << "," << std::setprecision(1)
       << pct(w.extra_cycles, w.blocks) << ","
       << (w.blocks > 0 ? double(w.total_block_instr) / w.blocks : 0) << ","
       << (w.blocks > 0 ? double(w.total_block_branches) / w.blocks : 0) << ","
       << pct(w.true_block_count, w.blocks) << ",";
    os << pct(w.provider_count[NUM_TABLES], w.branches) << ",";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << pct(w.provider_count[i], w.branches) << ",";
    os << pct(w.alloc_success, w.alloc_attempts) << ","
       << (w.pressure_samples > 0 ? double(w.acc_ctr_sum) / w.pressure_samples : 0)
       << ","
       << (w.pressure_samples > 0 ? double(w.alloc_ctr_sum) / w.pressure_samples : 0)
       << ",";
    {
      u64 last = NUM_TABLES - 1;
      auto &ls = w.lifetime_stats[last];
      double avg_p = ls.evictions > 0 ? double(ls.total_preds) / ls.evictions : 0;
      os << ls.evictions << "," << pct(ls.zero_use, ls.evictions) << ","
         << std::setprecision(1) << avg_p << ",";
    }
    os << pct(w.collision_hits, w.collision_checks) << ",";
    os << win_entry_ever_hit[NUM_TABLES - 1].count() << ",";
    {
      u64 win_stuck = 0, win_hard = 0;
      for (auto &[pc, d] : win_pc_diag) {
        if (d.total < 3) continue;
        if (d.alloc_count == 0 && d.misp > 0) win_stuck++;
        double bias = std::max(pct(d.taken, d.total),
                               100.0 - pct(d.taken, d.total));
        if (bias < 60.0) win_hard++;
      }
      os << win_stuck << "," << win_hard;
    }
    os << "," << pct(w.misp_direction, w.mispredictions)
       << "," << pct(w.misp_boundary, w.mispredictions);

    // Feature 17: micro MPKI p50/p95
    {
      u64 total_samples = 0;
      for (u64 i = 0; i < MICRO_MPKI_BINS; i++) total_samples += w.micro_mpki_hist[i];
      double p50 = 0, p95 = 0;
      u64 running = 0;
      for (u64 i = 0; i < MICRO_MPKI_BINS; i++) {
        running += w.micro_mpki_hist[i];
        double bucket_mpki = (i + 0.5) * 5.0;
        if (total_samples > 0 && running * 2 >= total_samples && p50 == 0)
          p50 = bucket_mpki;
        if (total_samples > 0 && running * 20 >= total_samples * 19 && p95 == 0)
          p95 = bucket_mpki;
      }
      os << "," << std::setprecision(1) << p50 << "," << p95;
    }
    // Feature 18: phase delta avg and max
    os << "," << std::setprecision(2)
       << (w.phase_samples > 0 ? w.phase_delta_sum / w.phase_samples : 0)
       << "," << w.phase_max_delta;
    // Feature 19: Jaccard distance
    {
      double jaccard = -1.0;
      if (!prev_window_pcs.empty()) {
        u64 intersection = 0;
        for (auto &[pc, d] : win_pc_diag)
          if (prev_window_pcs.count(pc)) intersection++;
        u64 union_sz = prev_window_pcs.size() + win_pc_diag.size() - intersection;
        jaccard = union_sz > 0 ? double(intersection) / union_sz : 1.0;
      }
      os << "," << std::setprecision(3) << jaccard;
    }
    // Feature 22: pingpong
    os << "," << w.pingpong_count;
    // Feature 24: counterfactual
    os << "," << pct(w.cf_fb_only, w.cf_total)
       << "," << pct(w.cf_tage_only, w.cf_total);
    // Provider source breakdown
    u64 w_prov_total = w.prov_no_tag_match + w.prov_sec_tag_reject +
                       w.prov_tage_meta_alt + w.prov_tage_meta_pri +
                       w.prov_tage_no_meta;
    os << "," << pct(w.prov_no_tag_match, w_prov_total)
       << "," << pct(w.prov_sec_tag_reject, w_prov_total)
       << "," << pct(w.prov_tage_meta_alt, w_prov_total)
       << "," << pct(w.prov_tage_meta_pri, w_prov_total)
       << "," << pct(w.prov_tage_no_meta, w_prov_total);
    // Sec-tag counterfactual accuracy
    os << "," << pct(w.cf_sec_fb_correct, w.cf_sec_tag_only_total)
       << "," << pct(w.cf_sec_tag_only_correct, w.cf_sec_tag_only_total);
    // Feature 25: benefit counter diagnostics
    os << "," << pct(w.ben_reject_all, w.blocks)
       << "," << pct(w.ben_incr, w.ben_reject_all)
       << "," << pct(w.ben_decr, w.ben_reject_all)
       << "," << (w.ben_ctr_samples > 0
                      ? double(w.ben_ctr_sum) / w.ben_ctr_samples : 0.0);

    os << "\n";
  }

  // ======== End-of-trace summary ========

  void print_summary(std::ostream &os = std::cerr) const {
    const auto &c = cum;

    os << "\n=== TageAhead Monitor Summary ===\n";
    double mpki = c.total_block_instr > 0
                      ? 1000.0 * c.mispredictions / c.total_block_instr
                      : 0.0;
    os << "Instructions: " << c.total_block_instr
       << "  Branches: " << c.branches
       << "  Mispredictions: " << c.mispredictions << " ("
       << std::fixed << std::setprecision(2)
       << pct(c.mispredictions, c.branches) << "%)"
       << "  MPKI: " << std::setprecision(3) << mpki << "\n";
    os << "Blocks: " << c.blocks
       << "  Extra cycles: " << c.extra_cycles << " ("
       << std::setprecision(2) << pct(c.extra_cycles, c.blocks) << "%)\n";
    os << "True blocks: " << c.true_block_count << " ("
       << pct(c.true_block_count, c.blocks) << "%)\n";
    os << "Train skips (first cycle): " << c.train_skip_count << "\n";

    // Pipeline calls
    os << "\nPipeline Calls:\n";
    u64 total_p1 = c.predict1_calls + c.reuse_predict1_calls;
    os << "  predict1: " << c.predict1_calls
       << "  reuse: " << c.reuse_predict1_calls << "  (reuse: "
       << pct(c.reuse_predict1_calls, total_p1) << "%)\n";
    os << "  predict2: " << c.predict2_calls
       << "  reuse: " << c.reuse_predict2_calls << "\n";
    os << "  Avg instr/block (from P1): "
       << (c.predict1_calls > 0 ? double(total_p1) / c.predict1_calls : 0)
       << "\n";

    // Block structure
    os << "\nBlock Structure:\n";
    os << std::setprecision(1);
    os << "  Avg instr/block: "
       << (c.blocks > 0 ? double(c.total_block_instr) / c.blocks : 0) << "\n";
    os << "  Avg branches/block: "
       << (c.blocks > 0 ? double(c.total_block_branches) / c.blocks : 0) << "\n";
    os << "  Instr/block histogram: ";
    for (u64 i = 1; i < MAX_BLOCK_INSTR && i <= 16; i++)
      if (c.block_instr_hist[i] > 0)
        os << i << ":" << c.block_instr_hist[i] << " ";
    os << "\n  Branches/block histogram: ";
    for (u64 i = 0; i < MAX_BLOCK_BR; i++)
      if (c.block_br_hist[i] > 0)
        os << i << ":" << c.block_br_hist[i] << " ";
    os << "\n  Entry point top-5: ";
    print_top_n(os, c.entry_point_hist, 5);
    os << "\n  Exit point top-5: ";
    print_top_n(os, c.exit_point_hist, 5);
    os << "\n";

    // Provider distribution
    os << "\nProvider Distribution:\n";
    os << "  Table     | Provided  |  Correct  | Accuracy | TagMatch% "
          "| SecMatch% | FullMatch%\n";
    os << "  ----------+-----------+-----------+----------+-----------"
          "+-----------+-----------\n";
    os << "  " << FB_NAME_PAD << "  |" << std::setw(10) << c.provider_count[NUM_TABLES]
       << " |" << std::setw(10) << c.provider_correct[NUM_TABLES] << " |"
       << std::setw(7) << pct(c.provider_correct[NUM_TABLES], c.provider_count[NUM_TABLES])
       << "%" << " |      --" << " |      --" << " |      --\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      os << "  T" << i << std::setw(8 - (i >= 10 ? 2 : 1)) << "" << "|"
         << std::setw(10) << c.provider_count[i] << " |" << std::setw(10)
         << c.provider_correct[i] << " |" << std::setw(7)
         << pct(c.provider_correct[i], c.provider_count[i]) << "%"
         << " |" << std::setw(8) << pct(c.tag_matches[i], c.tag_lookups[i]) << "%"
         << " |" << std::setw(8) << pct(c.sec_matches[i], c.tag_lookups[i]) << "%"
         << " |" << std::setw(8) << pct(c.full_matches[i], c.tag_lookups[i]) << "%\n";
    }

    // TAGE-only
    u64 tage_total = 0, tage_correct = 0;
    for (u64 i = 0; i < NUM_TABLES; i++) {
      tage_total += c.provider_count[i];
      tage_correct += c.provider_correct[i];
    }
    os << "\nTAGE-only (" << tage_total << " branches, "
       << pct(tage_total, c.branches) << "% of all):\n";
    if (tage_total > 0) {
      os << "  Table  | Prov%  | Alt%   | ProvAcc\n";
      os << "  -------+--------+--------+--------\n";
      for (u64 i = 0; i < NUM_TABLES; i++) {
        os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
           << std::setw(6) << pct(c.provider_count[i], tage_total) << "%"
           << " |" << std::setw(6) << pct(c.alt_count[i], tage_total) << "%"
           << " |" << std::setw(6)
           << pct(c.provider_correct[i], c.provider_count[i]) << "%\n";
      }
      os << "  Overall TAGE accuracy: " << pct(tage_correct, tage_total) << "%\n";
    }

    // Meta
    os << "\nMeta Override:\n";
    os << "  Total: " << c.meta_override_count
       << "  Correct: " << c.meta_override_correct << " ("
       << pct(c.meta_override_correct, c.meta_override_count) << "%)\n";
    os << "  Chose alt:     " << c.meta_chose_alt
       << "  Correct: " << c.meta_chose_alt_correct << " ("
       << pct(c.meta_chose_alt_correct, c.meta_chose_alt) << "%)\n";
    os << "  Chose primary: " << c.meta_chose_pri
       << "  Correct: " << c.meta_chose_pri_correct << " ("
       << pct(c.meta_chose_pri_correct, c.meta_chose_pri) << "%)\n";

    // Provider source breakdown
    u64 prov_total = c.prov_no_tag_match + c.prov_sec_tag_reject +
                     c.prov_tage_meta_alt + c.prov_tage_meta_pri +
                     c.prov_tage_no_meta;
    os << "\nProvider Source Breakdown:\n";
    os << "  No tag match (FB only):    " << c.prov_no_tag_match << " ("
       << pct(c.prov_no_tag_match, prov_total) << "%)  acc="
       << pct(c.prov_no_tag_match_correct, c.prov_no_tag_match) << "%\n";
    os << "  Sec-tag rejected (FB):     " << c.prov_sec_tag_reject << " ("
       << pct(c.prov_sec_tag_reject, prov_total) << "%)  acc="
       << pct(c.prov_sec_tag_reject_correct, c.prov_sec_tag_reject) << "%\n";
    os << "  TAGE provider, meta→alt:   " << c.prov_tage_meta_alt << " ("
       << pct(c.prov_tage_meta_alt, prov_total) << "%)  acc="
       << pct(c.prov_tage_meta_alt_correct, c.prov_tage_meta_alt) << "%\n";
    os << "  TAGE provider, meta→pri:   " << c.prov_tage_meta_pri << " ("
       << pct(c.prov_tage_meta_pri, prov_total) << "%)  acc="
       << pct(c.prov_tage_meta_pri_correct, c.prov_tage_meta_pri) << "%\n";
    os << "  TAGE provider, no meta:    " << c.prov_tage_no_meta << " ("
       << pct(c.prov_tage_no_meta, prov_total) << "%)  acc="
       << pct(c.prov_tage_no_meta_correct, c.prov_tage_no_meta) << "%\n";

    // Counterfactual: sec-tag rejection cost
    if (c.cf_sec_tag_only_total > 0) {
      os << "\n  Sec-tag counterfactual (tag-hit but sec-tag rejected):\n";
      os << "    Branches: " << c.cf_sec_tag_only_total << "\n";
      os << "    FB accuracy (actual):      "
         << pct(c.cf_sec_fb_correct, c.cf_sec_tag_only_total) << "%\n";
      os << "    Tag-only TAGE accuracy:    "
         << pct(c.cf_sec_tag_only_correct, c.cf_sec_tag_only_total) << "%\n";
      double delta = static_cast<double>(c.cf_sec_tag_only_correct) -
                     static_cast<double>(c.cf_sec_fb_correct);
      os << "    Delta (TAGE - FB):         " << std::showpos
         << std::fixed << std::setprecision(0) << delta << std::noshowpos
         << " branches (" << std::showpos << std::setprecision(2)
         << (c.cf_sec_tag_only_total > 0
             ? 100.0 * delta / c.cf_sec_tag_only_total : 0.0)
         << "%" << std::noshowpos << ")\n";
    }

    // Allocation
    os << "\nAllocation:\n";
    os << "  Attempts: " << c.alloc_attempts
       << "  Success: " << c.alloc_success << " ("
       << pct(c.alloc_success, c.alloc_attempts) << "%)"
       << "  Fail: " << c.alloc_fail
       << "  Blocked: " << c.alloc_blocked
       << "  Sibling skips: " << c.alloc_sibling_skip << "\n";
    os << "  Per table:";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << " T" << i << "=" << c.alloc_per_table[i];
    os << "\n";
    if (c.alloc_sibling_skip > 0) {
      os << "  Sibling per table:";
      for (u64 i = 0; i < NUM_TABLES; i++)
        os << " T" << i << "=" << c.alloc_sibling_per_table[i];
      os << "\n";
    }
    os << "  Unique branch PCs: " << unique_branch_pcs.size() << "\n";

    // Allocation cascade
    os << "\nAllocation Cascade (provider -> target):\n";
    os << "  Provider  | Allocs  |";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << std::setw(7) << "T" + std::to_string(i);
    os << "\n  ----------+---------+";
    for (u64 i = 0; i < NUM_TABLES; i++) os << "-------";
    os << "\n";
    for (u64 p = 0; p <= NUM_TABLES; p++) {
      if (c.alloc_from_provider[p] == 0) continue;
      if (p == NUM_TABLES)
        os << "  " << FB_NAME_PAD << "  |";
      else
        os << "  T" << p << std::setw(8 - (p >= 10 ? 2 : 1)) << "" << "|";
      os << std::setw(8) << c.alloc_from_provider[p] << " |";
      for (u64 t = 0; t < NUM_TABLES; t++)
        os << std::setw(7) << c.alloc_cascade[p][t];
      os << "\n";
    }

    // Table occupancy
    os << "\nTAGE Table Occupancy:\n";
    os << "  Table  | Size  | Occupied |  Occ% | EverHit | Used% | Accuracy | Allocs\n";
    os << "  -------+-------+----------+-------+---------+-------+----------+-------\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      u64 tsize = table_sizes[i] > 0 ? table_sizes[i] : MAX_TABLE_ENTRIES;
      u64 ever_hit = entry_ever_hit[i].count();
      os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
         << std::setw(6) << tsize << " |" << std::setw(9) << tage_unique_entries[i]
         << " |" << std::setw(5) << pct(tage_unique_entries[i], tsize) << "%"
         << " |" << std::setw(8) << ever_hit << " |" << std::setw(5)
         << pct(ever_hit, tsize) << "%"
         << " |" << std::setw(7) << pct(c.provider_correct[i], c.provider_count[i]) << "%"
         << " |" << std::setw(7) << c.alloc_per_table[i] << "\n";
    }

    // Fallback occupancy
    os << "\n" << FB_NAME << " Fallback (" << FB_CAPACITY << " entries):\n";
    os << "  Writes: " << fb_writes << "  Occupied: " << fb_unique_entries
       << " (" << pct(fb_unique_entries, FB_CAPACITY) << "%)\n";
    os << "  Accuracy: " << c.provider_correct[NUM_TABLES] << "/"
       << c.provider_count[NUM_TABLES] << " ("
       << pct(c.provider_correct[NUM_TABLES], c.provider_count[NUM_TABLES]) << "%)\n";

    // u-bit
    os << "\nU-bit:\n";
    os << "  Table  | U set    | U clear  | Turnover\n";
    os << "  -------+----------+----------+----------\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      u64 total_u = c.u_set_count[i] + c.u_clear_count[i];
      os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
         << std::setw(9) << c.u_set_count[i] << " |" << std::setw(9)
         << c.u_clear_count[i] << " |" << std::setw(7)
         << pct(c.u_clear_count[i], total_u) << "%\n";
    }
    os << "  Decay fires: " << c.decay_fire_count
       << "  Epoch resets: " << c.epoch_reset_count << "\n";

    // u-clear source
    u64 total_uclear_af = 0, total_uclear_sib = 0, total_uclear_ep = 0;
    for (u64 i = 0; i < NUM_TABLES; i++) {
      total_uclear_af += c.uclear_alloc_fail[i];
      total_uclear_sib += c.uclear_sibling[i];
      total_uclear_ep += c.uclear_epoch[i];
    }
    if (total_uclear_af + total_uclear_sib + total_uclear_ep > 0) {
      os << "\n  U-clear Source Breakdown:\n";
      os << "  Table  | AllocFail | Sibling   | Epoch\n";
      os << "  -------+-----------+-----------+-----------\n";
      for (u64 i = 0; i < NUM_TABLES; i++) {
        if (c.uclear_alloc_fail[i] + c.uclear_sibling[i] + c.uclear_epoch[i] == 0)
          continue;
        os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
           << std::setw(10) << c.uclear_alloc_fail[i] << " |"
           << std::setw(10) << c.uclear_sibling[i] << " |"
           << std::setw(10) << c.uclear_epoch[i] << "\n";
      }
      os << "  Total  |" << std::setw(10) << total_uclear_af << " |"
         << std::setw(10) << total_uclear_sib << " |"
         << std::setw(10) << total_uclear_ep << "\n";
    }

    // Pressure
    if (c.pressure_samples > 0) {
      os << "\nPressure Counters:\n";
      os << "  Avg acc_ctr: " << double(c.acc_ctr_sum) / c.pressure_samples << "\n";
      os << "  Avg alloc_ctr: " << double(c.alloc_ctr_sum) / c.pressure_samples << "\n";
    }

    // sec_tag consistency
    if (c.sec_tag_checks > 0) {
      os << "\nSec-Tag Consistency (tag-hit entries only):\n";
      os << "  Checks: " << c.sec_tag_checks << "\n";
      os << "  Match curr (predict-time): " << c.sec_tag_match_curr
         << " (" << pct(c.sec_tag_match_curr, c.sec_tag_checks) << "%)\n";
      os << "  Match now  (update-time):  " << c.sec_tag_match_now
         << " (" << pct(c.sec_tag_match_now, c.sec_tag_checks) << "%)\n";
      os << "  Match both:  " << c.sec_tag_match_both
         << " (" << pct(c.sec_tag_match_both, c.sec_tag_checks) << "%)\n";
      os << "  Match neither: " << c.sec_tag_match_neither
         << " (" << pct(c.sec_tag_match_neither, c.sec_tag_checks) << "%)\n";
    }

    // ================================================================
    // Diagnostic Features
    // ================================================================

    // Feature 1: Top mispredicted PCs
    {
      os << "\n--- Feature 1: Top Mispredicted PCs ---\n";
      std::vector<std::pair<u64, const PCDiag *>> sorted;
      sorted.reserve(pc_diag.size());
      for (auto &[pc, d] : pc_diag)
        sorted.push_back({pc, &d});
      std::sort(sorted.begin(), sorted.end(),
                [](auto &a, auto &b) { return a.second->misp > b.second->misp; });
      os << "  PC               | Total   | Misp    | Misp%  | TAGE%  "
            "| TageAcc | FB Acc  | Bias%  | Allocs\n";
      os << "  -----------------+---------+---------+--------+--------"
            "+---------+---------+--------+-------\n";
      u64 cumul_misp = 0;
      u64 shown = 0;
      for (auto &[pc, d] : sorted) {
        if (shown >= 20 || d->misp == 0) break;
        cumul_misp += d->misp;
        os << "  0x" << std::hex << std::setw(13) << std::left << pc
           << std::dec << std::right
           << " |" << std::setw(8) << d->total
           << " |" << std::setw(8) << d->misp
           << " |" << std::setw(6) << pct(d->misp, d->total) << "%"
           << " |" << std::setw(6) << pct(d->tage_prov, d->total) << "%"
           << " |" << std::setw(7) << pct(d->tage_correct, d->tage_prov) << "%"
           << " |" << std::setw(7) << pct(d->fb_correct, d->fb_prov) << "%"
           << " |" << std::setw(6) << pct(d->taken, d->total) << "%"
           << " |" << std::setw(7) << d->alloc_count
           << "  [cumul " << pct(cumul_misp, c.mispredictions) << "%]\n";
        shown++;
      }
      os << "  Top-20 account for " << cumul_misp << "/"
         << c.mispredictions << " (" << pct(cumul_misp, c.mispredictions)
         << "%) of all mispredictions\n";
    }

    // Feature 2: Branch bias distribution
    {
      os << "\n--- Feature 2: Branch Bias Distribution ---\n";
      constexpr u64 BINS = 10;
      std::array<u64, BINS> bias_hist{};
      u64 easy_count = 0, hard_count = 0;
      for (auto &[pc, d] : pc_diag) {
        if (d.total < 5) continue;
        double bias = std::max(pct(d.taken, d.total), 100.0 - pct(d.taken, d.total));
        u64 bin = std::min(u64(bias / 10.0), BINS - 1);
        bias_hist[bin]++;
        if (bias >= 90.0) easy_count++;
        else if (bias < 60.0) hard_count++;
      }
      os << "  Bias (max(T%,NT%)) | PCs\n";
      os << "  -------------------+------\n";
      for (u64 i = 0; i < BINS; i++)
        os << "  " << std::setw(3) << i * 10 << "-" << std::setw(3)
           << (i + 1) * 10 << "%            |" << std::setw(6) << bias_hist[i] << "\n";
      os << "  Easy (>=90% biased): " << easy_count
         << "  Hard (<60% biased): " << hard_count << "\n";
    }

    // Feature 3: TAGE vs Fallback
    {
      os << "\n--- Feature 3: TAGE vs Fallback Per-PC ---\n";
      u64 tage_dom = 0, fb_dominant = 0, mixed = 0;
      u64 tage_dom_misp = 0, fb_dominant_misp = 0, mixed_misp = 0;
      for (auto &[pc, d] : pc_diag) {
        if (d.total < 5) continue;
        double tage_frac = pct(d.tage_prov, d.total);
        if (tage_frac >= 50.0) { tage_dom++; tage_dom_misp += d.misp; }
        else if (tage_frac <= 10.0) { fb_dominant++; fb_dominant_misp += d.misp; }
        else { mixed++; mixed_misp += d.misp; }
      }
      os << "  Category           | PCs    | Mispredictions\n";
      os << "  -------------------+--------+---------------\n";
      os << "  TAGE dominant >50% |" << std::setw(7) << tage_dom
         << " |" << std::setw(14) << tage_dom_misp << "\n";
      os << "  Fallback dom <=10% |" << std::setw(7) << fb_dominant
         << " |" << std::setw(14) << fb_dominant_misp << "\n";
      os << "  Mixed 10-50%       |" << std::setw(7) << mixed
         << " |" << std::setw(14) << mixed_misp << "\n";
    }

    // Feature 5: Entry lifetime
    {
      os << "\n--- Feature 5: Entry Lifetime ---\n";
      os << "  Table  | Evictions | AvgPreds | AvgAcc%  | ZeroUse | ZeroUse%\n";
      os << "  -------+-----------+----------+----------+---------+---------\n";
      for (u64 i = 0; i < NUM_TABLES; i++) {
        auto &ls = c.lifetime_stats[i];
        double avg_preds = ls.evictions > 0 ? double(ls.total_preds) / ls.evictions : 0;
        double avg_acc = ls.total_preds > 0 ? pct(ls.total_correct, ls.total_preds) : 0;
        os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
           << std::setw(10) << ls.evictions << " |" << std::setw(8)
           << std::setprecision(1) << avg_preds << " |" << std::setw(7)
           << std::setprecision(1) << avg_acc << "% |" << std::setw(8)
           << ls.zero_use << " |" << std::setw(7)
           << pct(ls.zero_use, ls.evictions) << "%\n";
      }
    }

    // Feature 7: Tag collision
    {
      os << "\n--- Feature 7: Tag Collision Rate ---\n";
      os << "  Total checks: " << c.collision_checks
         << "  Collisions: " << c.collision_hits << " ("
         << pct(c.collision_hits, c.collision_checks) << "%)\n";
      os << "  Table  | Checks   | Collisions | Rate\n";
      os << "  -------+----------+------------+------\n";
      for (u64 i = 0; i < NUM_TABLES; i++) {
        os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
           << std::setw(9) << c.per_table_coll_checks[i] << " |" << std::setw(11)
           << c.per_table_coll_hits[i] << " |" << std::setw(5)
           << pct(c.per_table_coll_hits[i], c.per_table_coll_checks[i]) << "%\n";
      }
    }

    // Feature 8: Fallback-stuck
    {
      os << "\n--- Feature 8: Fallback-Stuck Branches ---\n";
      u64 never_alloc = 0, alloc_but_fb = 0, alloc_and_tage = 0;
      u64 never_alloc_misp = 0, alloc_but_fb_misp = 0, alloc_tage_misp = 0;
      for (auto &[pc, d] : pc_diag) {
        if (d.total < 5) continue;
        if (d.alloc_count == 0) { never_alloc++; never_alloc_misp += d.misp; }
        else if (d.tage_prov * 2 < d.total) { alloc_but_fb++; alloc_but_fb_misp += d.misp; }
        else { alloc_and_tage++; alloc_tage_misp += d.misp; }
      }
      os << "  Category                     | PCs    | Mispred  | AvgAlloc\n";
      os << "  -----------------------------+--------+----------+---------\n";
      os << "  Never allocated to TAGE      |" << std::setw(7) << never_alloc
         << " |" << std::setw(9) << never_alloc_misp << " |      --\n";
      os << "  Allocated but FB dominates   |" << std::setw(7) << alloc_but_fb
         << " |" << std::setw(9) << alloc_but_fb_misp << " |";
      { u64 ta = 0, cnt = 0;
        for (auto &[pc, d] : pc_diag)
          if (d.total >= 5 && d.alloc_count > 0 && d.tage_prov * 2 < d.total) { ta += d.alloc_count; cnt++; }
        os << std::setw(8) << (cnt > 0 ? double(ta) / cnt : 0) << "\n";
      }
      os << "  Allocated and TAGE provides  |" << std::setw(7) << alloc_and_tage
         << " |" << std::setw(9) << alloc_tage_misp << " |";
      { u64 ta = 0, cnt = 0;
        for (auto &[pc, d] : pc_diag)
          if (d.total >= 5 && d.alloc_count > 0 && d.tage_prov * 2 >= d.total) { ta += d.alloc_count; cnt++; }
        os << std::setw(8) << (cnt > 0 ? double(ta) / cnt : 0) << "\n";
      }
    }

    // Per-rank accuracy
    {
      os << "\n--- Per-Rank Branch Accuracy ---\n";
      os << "  Rank | Branches  | Correct   | Accuracy\n";
      os << "  -----+-----------+-----------+---------\n";
      for (u64 r = 0; r < N; r++) {
        if (c.rank_total[r] == 0) break;
        os << "  " << std::setw(4) << r << " |" << std::setw(10) << c.rank_total[r]
           << " |" << std::setw(10) << c.rank_correct[r] << " |" << std::setw(7)
           << pct(c.rank_correct[r], c.rank_total[r]) << "%\n";
      }
    }

    // Feature 10: Misprediction type breakdown
    {
      os << "\n--- Feature 10: Misprediction Type Breakdown ---\n";
      os << "  Classifies each misprediction into two categories:\n";
      os << "    Direction:  the predicted taken/not-taken was wrong for the providing branch\n";
      os << "    Boundary:   the last branch in the block was predicted correctly, but an\n";
      os << "                earlier branch in the same fetch block was mispredicted, causing\n";
      os << "                the block to end at the wrong point\n\n";
      os << "  Total mispredictions: " << c.mispredictions << "\n";
      os << "  Direction (pred!=actual): " << c.misp_direction
         << " (" << pct(c.misp_direction, c.mispredictions) << "%)\n";
      os << "  Boundary  (block-end wrong): " << c.misp_boundary
         << " (" << pct(c.misp_boundary, c.mispredictions) << "%)\n";
    }

    // Feature 11: Per-PC ahead-context fanout
    {
      os << "\n--- Feature 11: Per-PC Ahead-Context Fanout ---\n";
      std::vector<std::pair<u64, const PCDiag *>> sorted;
      sorted.reserve(pc_diag.size());
      for (auto &[pc, d] : pc_diag)
        if (d.misp > 0) sorted.push_back({pc, &d});
      std::sort(sorted.begin(), sorted.end(),
                [](auto &a, auto &b) { return a.second->misp > b.second->misp; });
      os << "  PC               | Misp    | Contexts | Misp/Ctx | Diagnosis\n";
      os << "  -----------------+---------+----------+----------+----------\n";
      u64 shown = 0;
      for (auto &[pc, d] : sorted) {
        if (shown >= 30) break;
        auto it = pc_block_contexts.find(pc);
        u64 ctx = (it != pc_block_contexts.end()) ? it->second.size() : 0;
        double misp_per_ctx = ctx > 0 ? double(d->misp) / ctx : 0;
        const char *diag = (ctx <= 1) ? "single-context"
                           : (misp_per_ctx > 5.0) ? "hard-branch" : "context-fragmented";
        os << "  0x" << std::hex << std::setw(13) << std::left << pc
           << std::dec << std::right
           << " |" << std::setw(8) << d->misp << " |" << std::setw(9) << ctx
           << " |" << std::setw(8) << std::setprecision(1) << misp_per_ctx
           << " | " << diag << "\n";
        shown++;
      }
      u64 single_ctx = 0, fragmented = 0, hard = 0;
      u64 single_misp = 0, frag_misp = 0, hard_misp = 0;
      for (auto &[pc, d] : sorted) {
        auto it = pc_block_contexts.find(pc);
        u64 ctx = (it != pc_block_contexts.end()) ? it->second.size() : 0;
        double mpc = ctx > 0 ? double(d->misp) / ctx : 0;
        if (ctx <= 1) { single_ctx++; single_misp += d->misp; }
        else if (mpc > 5.0) { hard++; hard_misp += d->misp; }
        else { fragmented++; frag_misp += d->misp; }
      }
      os << "  Summary:\n";
      os << "    Single-context:     " << single_ctx << " PCs, " << single_misp << " mispredictions\n";
      os << "    Hard-branch (>5/ctx): " << hard << " PCs, " << hard_misp << " mispredictions\n";
      os << "    Context-fragmented: " << fragmented << " PCs, " << frag_misp << " mispredictions\n";
    }

    // Feature 12: Per-PC × block-position matrix
    {
      os << "\n--- Feature 12: Per-PC x Block-Position Misprediction ---\n";
      os << "  Cross-tabulates mispredictions by PC and rank (position within fetch block).\n";
      os << "  Columns:\n";
      os << "    PC:      branch program counter\n";
      os << "    Total:   total observations of this PC across all positions\n";
      os << "    Misp:    total mispredictions of this PC across all positions\n";
      os << "    R0..Rn:  per-position breakdown shown as misp/total; R0 is the first\n";
      os << "             branch in the fetch block, R1 the second, etc.\n";
      os << "    Conc?:   YES if >2/3 of this PC's mispredictions concentrate at a single\n";
      os << "             position (suggests the block boundary, not the branch itself,\n";
      os << "             is the root cause)\n\n";
      std::vector<std::pair<u64, const PCDiag *>> sorted;
      sorted.reserve(pc_diag.size());
      for (auto &[pc, d] : pc_diag)
        if (d.misp > 0) sorted.push_back({pc, &d});
      std::sort(sorted.begin(), sorted.end(),
                [](auto &a, auto &b) { return a.second->misp > b.second->misp; });
      u64 max_rank = 0;
      for (u64 r = 0; r < MAX_RANK; r++)
        if (c.rank_total[r] > 0) max_rank = r + 1;
      if (max_rank == 0) max_rank = 1;
      // Header
      os << "  PC               | Total   | Misp    |";
      for (u64 r = 0; r < max_rank; r++) {
        os << "  R" << std::setw(1) << r << "     |";
      }
      os << " Conc?\n";
      os << "  -----------------+---------+---------+";
      for (u64 r = 0; r < max_rank; r++) os << "---------+";
      os << "------\n";
      u64 shown = 0, position_concentrated = 0;
      for (auto &[pc, d] : sorted) {
        if (shown >= 20) break;
        os << "  0x" << std::hex << std::setw(13) << std::left << pc
           << std::dec << std::right
           << " |" << std::setw(8) << d->total << " |" << std::setw(8) << d->misp << " |";
        u64 max_pos_misp = 0, active_pos = 0;
        for (u64 r = 0; r < max_rank; r++) {
          if (d->pos_misp[r] > max_pos_misp) max_pos_misp = d->pos_misp[r];
          if (d->pos_total[r] > 0) active_pos++;
        }
        bool conc = (d->misp > 2 && active_pos > 1 && max_pos_misp > d->misp * 2 / 3);
        for (u64 r = 0; r < max_rank; r++) {
          if (d->pos_total[r] > 0) {
            os << std::setw(4) << d->pos_misp[r] << "/" << std::setw(4) << d->pos_total[r] << "|";
          } else {
            os << "    --   |";
          }
        }
        os << (conc ? " YES" : "  no") << "\n";
        if (conc) position_concentrated++;
        shown++;
      }
      os << "  Position-concentrated: " << position_concentrated
         << "/" << std::min(u64(20), u64(sorted.size())) << " of top mispredicted PCs\n";
      os << "\n  Aggregate per-position stats:\n";
      os << "  (Summed across all PCs — shows overall accuracy at each block rank)\n";
      os << "  Rank | Branches  | Misp      | MispRate\n";
      os << "  -----+-----------+-----------+---------\n";
      for (u64 r = 0; r < max_rank; r++) {
        u64 pt = 0, pm = 0;
        for (auto &[pc, d] : pc_diag) { pt += d.pos_total[r]; pm += d.pos_misp[r]; }
        if (pt == 0) continue;
        os << "  R" << std::setw(3) << std::left << r << std::right
           << " |" << std::setw(10) << pt
           << " |" << std::setw(10) << pm << " |" << std::setw(7) << pct(pm, pt) << "%\n";
      }
    }

    // Feature 13: Per-PC allocation outcome histogram
    {
      os << "\n--- Feature 13: Per-PC Allocation Outcome Histogram ---\n";
      std::array<u64, ALLOC_HIST_BINS> ah_pcs{}, ah_misp{};
      for (auto &[pc, d] : pc_diag) {
        if (d.alloc_count == 0) continue;
        u64 b = alloc_hist_bin(d.alloc_count);
        ah_pcs[b]++;
        ah_misp[b] += d.misp;
      }
      os << "  AllocCount | PCs     | Cumul Mispred\n";
      os << "  -----------+---------+--------------\n";
      for (u64 b = 0; b < ALLOC_HIST_BINS; b++) {
        os << "  " << std::setw(10) << alloc_hist_label(b) << " |"
           << std::setw(8) << ah_pcs[b] << " |" << std::setw(13) << ah_misp[b] << "\n";
      }
    }

    // Feature 14: TAGE accuracy trajectory per alloc count
    {
      os << "\n--- Feature 14: TAGE Accuracy Trajectory by Alloc Count ---\n";
      std::array<u64, TRAJ_BINS> total_prov{}, total_correct{};
      for (auto &[pc, d] : pc_diag) {
        for (u64 b = 0; b < TRAJ_BINS; b++) {
          total_prov[b] += d.traj_tage_prov[b];
          total_correct[b] += d.traj_tage_correct[b];
        }
      }
      os << "  AllocCount | TageProvided | TageAcc%\n";
      os << "  -----------+--------------+---------\n";
      for (u64 b = 0; b < TRAJ_BINS; b++) {
        os << "  " << std::setw(10) << traj_label(b) << " |"
           << std::setw(13) << total_prov[b] << " |" << std::setw(7)
           << pct(total_correct[b], total_prov[b]) << "%\n";
      }
    }

    // Feature 15: Per-PC oscillation detector
    {
      os << "\n--- Feature 15: Per-PC Oscillation Detector ---\n";
      // Focus on top mispredicted PCs
      std::vector<std::pair<u64, const PCDiag *>> sorted;
      sorted.reserve(pc_diag.size());
      for (auto &[pc, d] : pc_diag)
        if (d.misp > 0) sorted.push_back({pc, &d});
      std::sort(sorted.begin(), sorted.end(),
                [](auto &a, auto &b) { return a.second->misp > b.second->misp; });
      os << "  PC               | Total   | Misp    | AvgStreak | MaxStreak | Streaks | Diagnosis\n";
      os << "  -----------------+---------+---------+-----------+-----------+---------+----------\n";
      u64 shown = 0;
      for (auto &[pc, d] : sorted) {
        if (shown >= 20) break;
        double avg_streak = d->streak_count > 0 ? double(d->streak_sum) / d->streak_count : 0;
        const char *diag = (d->streak_count < 3) ? "too-few"
                           : (avg_streak < 3.0) ? "oscillating"
                           : (avg_streak < 10.0) ? "moderate"
                                                 : "stable-then-fail";
        os << "  0x" << std::hex << std::setw(13) << std::left << pc
           << std::dec << std::right
           << " |" << std::setw(8) << d->total
           << " |" << std::setw(8) << d->misp
           << " |" << std::setw(10) << std::setprecision(1) << avg_streak
           << " |" << std::setw(10) << d->max_streak
           << " |" << std::setw(8) << d->streak_count
           << " | " << diag << "\n";
        shown++;
      }
    }

    // Feature 16: Fallback accuracy for stuck PCs
    {
      os << "\n--- Feature 16: Fallback vs TAGE Per-PC (TAGE-provided branches) ---\n";
      u64 fb_better = 0, tage_better = 0, tie = 0;
      u64 fb_better_misp = 0, tage_better_misp = 0;
      for (auto &[pc, d] : pc_diag) {
        if (d.cf_tage_provides < 10) continue;
        double fb_acc = pct(d.cf_fb_would_correct, d.cf_tage_provides);
        double tage_acc = pct(d.cf_tage_correct_here, d.cf_tage_provides);
        if (fb_acc > tage_acc + 1.0) { fb_better++; fb_better_misp += d.misp; }
        else if (tage_acc > fb_acc + 1.0) { tage_better++; tage_better_misp += d.misp; }
        else { tie++; }
      }
      os << "  PCs where FB would be better: " << fb_better
         << " (mispredictions: " << fb_better_misp << ")\n";
      os << "  PCs where TAGE is better:     " << tage_better
         << " (mispredictions: " << tage_better_misp << ")\n";
      os << "  PCs tied (±1%):               " << tie << "\n";
      // Top 10 PCs where FB beats TAGE
      std::vector<std::tuple<double, u64, const PCDiag *>> fb_wins;
      for (auto &[pc, d] : pc_diag) {
        if (d.cf_tage_provides < 10) continue;
        double fb_acc = pct(d.cf_fb_would_correct, d.cf_tage_provides);
        double tage_acc = pct(d.cf_tage_correct_here, d.cf_tage_provides);
        if (fb_acc > tage_acc + 1.0)
          fb_wins.push_back({fb_acc - tage_acc, pc, &d});
      }
      std::sort(fb_wins.begin(), fb_wins.end(),
                [](auto &a, auto &b) { return std::get<0>(a) > std::get<0>(b); });
      if (!fb_wins.empty()) {
        os << "  Top PCs where FB beats TAGE:\n";
        os << "  PC               | TageN   | FBAcc%  | TageAcc%| Delta  | Misp\n";
        os << "  -----------------+---------+---------+---------+--------+------\n";
        u64 shown = 0;
        for (auto &[delta, pc, d] : fb_wins) {
          if (shown >= 10) break;
          os << "  0x" << std::hex << std::setw(13) << std::left << pc
             << std::dec << std::right
             << " |" << std::setw(8) << d->cf_tage_provides
             << " |" << std::setw(7) << pct(d->cf_fb_would_correct, d->cf_tage_provides) << "%"
             << " |" << std::setw(7) << pct(d->cf_tage_correct_here, d->cf_tage_provides) << "%"
             << " |" << std::setw(6) << std::setprecision(1) << delta << "%"
             << " |" << std::setw(6) << d->misp << "\n";
          shown++;
        }
      }
    }

    // Feature 17: Sliding-window MPKI distribution
    {
      os << "\n--- Feature 17: Micro-Window MPKI Distribution (1K-branch windows) ---\n";
      os << "  MPKI Range | Samples   | Fraction\n";
      os << "  -----------+-----------+---------\n";
      u64 total = 0;
      for (u64 i = 0; i < MICRO_MPKI_BINS; i++) total += c.micro_mpki_hist[i];
      for (u64 i = 0; i < MICRO_MPKI_BINS; i++) {
        os << "  ";
        if (i < MICRO_MPKI_BINS - 1)
          os << std::setw(3) << i * 5 << "-" << std::setw(3) << (i + 1) * 5;
        else
          os << "  50+  ";
        os << "    |" << std::setw(10) << c.micro_mpki_hist[i]
           << " |" << std::setw(7) << pct(c.micro_mpki_hist[i], total) << "%\n";
      }
    }

    // Feature 18: Phase detection — short-window misp rate delta
    {
      os << "\n--- Feature 18: Phase Detection (misp rate delta, " << PHASE_WIN << "-branch windows) ---\n";
      os << "  Compares misprediction rates between consecutive " << PHASE_WIN
         << "-branch windows.\n";
      os << "  A large delta signals a phase change in program behavior.\n\n";
      os << "  |Delta Rate| | Samples   | Fraction\n";
      os << "  ------------+-----------+---------\n";
      for (u64 i = 0; i < PHASE_BINS; i++) {
        os << "  " << std::setw(11) << phase_label(i) << " |"
           << std::setw(10) << c.phase_delta_hist[i] << " |" << std::setw(7)
           << pct(c.phase_delta_hist[i], c.phase_samples) << "%\n";
      }
      if (c.phase_samples > 0) {
        os << "  Mean |delta|: " << std::setprecision(2)
           << (c.phase_delta_sum / c.phase_samples) << "%";
        os << "  Max |delta|: " << c.phase_max_delta << "%\n";
      }
    }

    // Feature 20: Entry age at provider
    {
      os << "\n--- Feature 20: Entry Age at Provider ---\n";
      os << "  Age Range   | Providers | Accuracy\n";
      os << "  ------------+-----------+---------\n";
      for (u64 b = 0; b < AGE_BINS; b++) {
        os << "  " << std::setw(11) << age_label(b) << " |"
           << std::setw(10) << c.age_prov_count[b] << " |" << std::setw(7)
           << pct(c.age_prov_correct[b], c.age_prov_count[b]) << "%\n";
      }
    }

    // Feature 21: Per-PC recurrence distance
    {
      os << "\n--- Feature 21: Misprediction Recurrence Distance ---\n";
      os << "  Distance    | Count\n";
      os << "  ------------+------\n";
      for (u64 b = 0; b < RECUR_BINS; b++) {
        os << "  " << std::setw(11) << recur_label(b) << " |"
           << std::setw(10) << c.recurrence_hist[b] << "\n";
      }
    }

    // Feature 22: Allocation ping-pong
    {
      os << "\n--- Feature 22: Allocation Ping-Pong ---\n";
      os << "  Ping-pong evictions: " << c.pingpong_count << "\n";
    }

    // Feature 23: Sec-tag mismatch breakdown
    {
      u64 total_mm = c.sec_mm_same_pc + c.sec_mm_diff_pc + c.sec_mm_invalid;
      if (total_mm > 0) {
        os << "\n--- Feature 23: Sec-Tag Mismatch Breakdown ---\n";
        os << "  Total mismatches: " << total_mm << "\n";
        os << "  Same PC, different path: " << c.sec_mm_same_pc
           << " (" << pct(c.sec_mm_same_pc, total_mm) << "%)\n";
        os << "  Different PC (aliasing):  " << c.sec_mm_diff_pc
           << " (" << pct(c.sec_mm_diff_pc, total_mm) << "%)\n";
        os << "  Invalid/no entry:         " << c.sec_mm_invalid
           << " (" << pct(c.sec_mm_invalid, total_mm) << "%)\n";
      }
    }

    // Feature 24: Counterfactual fallback
    {
      os << "\n--- Feature 24: Counterfactual Fallback Analysis ---\n";
      os << "  Total branches: " << c.cf_total << "\n";
      os << "  FB would be correct:   " << c.cf_fb_correct
         << " (" << pct(c.cf_fb_correct, c.cf_total) << "%)\n";
      os << "  TAGE was correct:      " << c.cf_tage_correct
         << " (" << pct(c.cf_tage_correct, c.cf_total) << "%)\n";
      os << "  FB only (FB right, TAGE wrong): " << c.cf_fb_only
         << " (" << pct(c.cf_fb_only, c.cf_total) << "%)\n";
      os << "  TAGE only (TAGE right, FB wrong): " << c.cf_tage_only
         << " (" << pct(c.cf_tage_only, c.cf_total) << "%)\n";
      if (c.cf_total > 0) {
        double fb_mpki = c.total_block_instr > 0
                             ? 1000.0 * (c.cf_total - c.cf_fb_correct) / c.total_block_instr
                             : 0;
        os << "  Hypothetical FB-only MPKI: " << std::setprecision(3) << fb_mpki << "\n";
      }
    }

    // Feature 25: benefit counter diagnostics
    if (c.ben_ctr_samples > 0 || c.ben_reject_all > 0) {
      os << "\n--- Feature 25: Sec-Tag Adaptive Benefit Counter ---\n";
      os << "  sec_would_reject_all fired: " << c.ben_reject_all
         << " (" << pct(c.ben_reject_all, c.blocks) << "% of blocks)\n";
      u64 ben_total = c.ben_incr + c.ben_decr + c.ben_tie;
      os << "  Decisive outcomes: " << ben_total << "\n";
      os << "    Increments (FB > TAGE, enforce helped): " << c.ben_incr
         << " (" << pct(c.ben_incr, ben_total) << "%)\n";
      os << "    Decrements (TAGE > FB, enforce hurt):   " << c.ben_decr
         << " (" << pct(c.ben_decr, ben_total) << "%)\n";
      os << "    Ties (no update):                       " << c.ben_tie
         << " (" << pct(c.ben_tie, c.ben_reject_all) << "% of fired)\n";
      os << "  Avg benefit_ctr value: "
         << std::setprecision(1)
         << (c.ben_ctr_samples > 0
                 ? double(c.ben_ctr_sum) / c.ben_ctr_samples : 0.0) << "\n";
    }

    os << "\n=== End TageAhead Monitor ===\n\n";
  }

private:
  template <size_t SZ>
  static void print_top_n(std::ostream &os, const std::array<u64, SZ> &hist, u64 n) {
    std::array<std::pair<u64, u64>, 64> top{};
    u64 found = 0;
    for (u64 i = 0; i < SZ && found < 64; i++)
      if (hist[i] > 0) top[found++] = {hist[i], i};
    for (u64 i = 0; i < std::min(n, found); i++) {
      for (u64 j = i + 1; j < found; j++)
        if (top[j].first > top[i].first) std::swap(top[i], top[j]);
      os << top[i].second << ":" << top[i].first << " ";
    }
  }
};
