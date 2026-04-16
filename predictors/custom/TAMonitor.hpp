#pragma once
// TageAhead Performance Monitor — software-only instrumentation
// Gated by -DTAGE_MONITOR at compile time.
// All counters are plain C++ types (zero HARCOM cost).

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using u64 = uint64_t;

template <u64 NUM_TABLES, u64 N, u64 MAX_TABLE_ENTRIES = 2048,
          bool USE_GSHARE = false, u64 FB_CAPACITY = 8192,
          u64 HYST_BANKS = 1, u64 U_BANKS = 1, u64 NUM_PATHS = 1>
struct TAMonitor {

  static constexpr const char *FB_NAME = USE_GSHARE ? "Gshare" : "Bimodal";
  static constexpr const char *FB_NAME_PAD = USE_GSHARE ? "Gshare  " : "Bimodal ";

  static constexpr u64 WINDOW_SIZE = 100000; // branches per window

  static constexpr u64 MAX_BLOCK_INSTR = 64;
  static constexpr u64 MAX_BLOCK_BR = N + 1;

  // Feature 5: Lifetime stats (defined before Counters so it can be embedded)
  struct LifetimeStats {
    u64 evictions = 0;
    u64 total_preds = 0;
    u64 total_correct = 0;
    u64 zero_use = 0; // evicted before serving any prediction
  };

  // Feature 1+2+3+8: Per-PC diagnostics
  struct PCDiag {
    u64 total = 0;
    u64 misp = 0;
    u64 taken = 0;         // bias tracking (feature 2)
    u64 tage_prov = 0;     // TAGE was provider (feature 3)
    u64 tage_correct = 0;
    u64 fb_prov = 0;       // fallback was provider
    u64 fb_correct = 0;
    u64 alloc_count = 0;   // allocated to TAGE (feature 8)
  };

  // ======== Counters (cumulative + per-window) ========
  struct Counters {
    u64 branches = 0;
    u64 mispredictions = 0;
    u64 blocks = 0;
    u64 extra_cycles = 0;

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

    // True block rate (pipeline reuse)
    u64 true_block_count = 0;

    // Train valid skip (first-cycle guard)
    u64 train_skip_count = 0;

    // Provider distribution: 0..NUM_TABLES-1 = TAGE, NUM_TABLES = fallback
    std::array<u64, NUM_TABLES + 1> provider_count{};
    std::array<u64, NUM_TABLES + 1> provider_correct{};

    // Alt provider
    std::array<u64, NUM_TABLES + 1> alt_count{};

    // Meta override
    u64 meta_override_count = 0;
    u64 meta_override_correct = 0;
    u64 meta_chose_alt = 0;
    u64 meta_chose_alt_correct = 0;
    u64 meta_chose_pri = 0;
    u64 meta_chose_pri_correct = 0;

    // Tag match rate per table (primary + secondary)
    std::array<u64, NUM_TABLES> tag_lookups{};
    std::array<u64, NUM_TABLES> tag_matches{};
    std::array<u64, NUM_TABLES> sec_matches{};
    std::array<u64, NUM_TABLES> full_matches{};

    // Secondary tag rejection analysis: would the rejected entry have been correct?
    std::array<u64, NUM_TABLES> sec_reject_count{};      // tag hit, sec miss
    std::array<u64, NUM_TABLES> sec_reject_correct{};    // rejected but would have predicted correctly
    std::array<u64, NUM_TABLES> sec_reject_wrong{};      // rejected and would have predicted wrong

    // Multi-path sec tag banking (NUM_PATHS > 1)
    std::array<std::array<u64, NUM_PATHS>, NUM_TABLES> path_hit{};    // per-table, per-path sec match
    std::array<u64, NUM_TABLES> path_neither{};    // neither path matched (bimodal fallback)
    std::array<u64, NUM_TABLES> path_both{};       // both paths matched (degenerate)
    std::array<std::array<u64, NUM_PATHS>, NUM_TABLES> provider_from_path{};  // which path provided
    std::array<std::array<u64, NUM_PATHS>, NUM_TABLES> provider_path_correct{}; // correct per path
    std::array<std::array<u64, NUM_PATHS>, NUM_TABLES> alloc_to_path{};  // allocation target path

    // Allocation
    u64 alloc_attempts = 0;
    u64 alloc_success = 0;
    u64 alloc_fail = 0;
    u64 alloc_blocked = 0; // mispredict but fallback was provider (no candidate)
    std::array<u64, NUM_TABLES> alloc_per_table{};
    // Cascade: alloc_from_provider[i] = allocations when Ti was provider
    std::array<u64, NUM_TABLES + 1> alloc_from_provider{};
    std::array<std::array<u64, NUM_TABLES>, NUM_TABLES + 1> alloc_cascade{};

    // u-bit (per table, per bank)
    std::array<std::array<u64, U_BANKS>, NUM_TABLES> u_set_count{};
    std::array<std::array<u64, U_BANKS>, NUM_TABLES> u_clear_count{};
    u64 decay_fire_count = 0;
    u64 epoch_reset_count = 0;

    // hyst training (per table, per bank)
    std::array<std::array<u64, HYST_BANKS>, NUM_TABLES> hyst_inc_count{};
    std::array<std::array<u64, HYST_BANKS>, NUM_TABLES> hyst_dec_count{};

    // Pressure counters (sampled each update_cycle with branches)
    u64 acc_ctr_sum = 0;
    u64 alloc_ctr_sum = 0;
    u64 pressure_samples = 0;

    // Feature 5: Lifetime stats per table
    std::array<LifetimeStats, NUM_TABLES> lifetime_stats{};

    // Feature 7: Tag collision tracking
    u64 collision_checks = 0;
    u64 collision_hits = 0;
    std::array<u64, NUM_TABLES> per_table_coll_checks{};
    std::array<u64, NUM_TABLES> per_table_coll_hits{};

    void reset() { *this = Counters{}; }
  };

  Counters cum;
  Counters win;
  u64 window_num = 0;
  bool header_printed = false;

  // Table occupancy tracking (cumulative only — snapshots)
  std::array<std::bitset<MAX_TABLE_ENTRIES>, NUM_TABLES> tage_occupied{};
  std::array<u64, NUM_TABLES> tage_unique_entries{};

  // Fallback occupancy (gshare mode only)
  std::bitset<FB_CAPACITY> fb_occupied{};
  u64 fb_unique_entries = 0;
  u64 fb_writes = 0;

  // Unique branch PCs
  std::unordered_set<u64> unique_branch_pcs;

  // Shadow state (set during resolution, consumed during training)
  std::array<u64, N> shadow_provider{};
  std::array<u64, N> shadow_alt{};
  std::array<bool, N> shadow_meta_overrode{};
  std::array<bool, N> shadow_meta_chose_alt{};
  std::array<bool, N> shadow_pred{};
  std::array<u64, N> shadow_pc{};

  // Per-PC diagnostics (cumulative + per-window)
  std::unordered_map<u64, PCDiag> pc_diag;
  std::unordered_map<u64, PCDiag> win_pc_diag;

  // Feature 5: Entry lifetime state (per-entry, not windowed)
  struct EntryLife {
    u64 alloc_time = 0;
    u64 pred_count = 0;
    u64 correct_count = 0;
    u64 stored_pc = 0; // full PC for collision detection (feature 7)
    bool valid = false;
  };
  std::array<std::array<EntryLife, MAX_TABLE_ENTRIES>, NUM_TABLES> entry_life{};
  u64 global_branch_count = 0;

  // Feature 9: Entry utilization (cumulative + per-window)
  std::array<std::bitset<MAX_TABLE_ENTRIES>, NUM_TABLES> entry_ever_hit{};
  std::array<std::bitset<MAX_TABLE_ENTRIES>, NUM_TABLES> win_entry_ever_hit{};

  // Feature 10: Shadow u-values for histogram (per table, per entry, per bank)
  static constexpr u64 MAX_U_VAL = 16; // enough for U_WIDTH up to 4
  std::array<std::array<std::array<u64, U_BANKS>, MAX_TABLE_ENTRIES>, NUM_TABLES> u_shadow{};

  // Shadow per-table prediction state (for sec tag rejection analysis)
  // Saved at predict2, consumed at training (one cycle later)
  struct TablePredShadow {
    u64 pred_bits = 0;     // PRED_BITS-wide prediction from this table
    bool tag_hit = false;  // primary tag matched
    bool sec_hit = false;  // secondary tag matched
  };
  std::array<TablePredShadow, NUM_TABLES> shadow_table_pred{};
  std::array<TablePredShadow, NUM_TABLES> train_table_pred{};  // piped one cycle
  u64 shadow_num_branch = 0;
  u64 train_num_branch = 0;

  // Block PC pipeline (for collision/lifetime tracking)
  u64 shadow_block_pc = 0;   // set in predict1
  u64 current_block_pc = 0;  // shifted in update_cycle
  u64 train_block_pc = 0;    // shifted in update_cycle

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

  // Per-table prediction shadow (called at predict2 for each table)
  void record_table_pred(u64 table, u64 pred_bits, bool tag_hit, bool sec_hit) {
    shadow_table_pred[table] = {pred_bits, tag_hit, sec_hit};
  }
  void record_table_pred_num_branch(u64 nb) { shadow_num_branch = nb; }

  // Pipeline shift for sec rejection shadow (called at start of update_cycle)
  void shift_table_pred_shadow() {
    train_table_pred = shadow_table_pred;
    train_num_branch = shadow_num_branch;
  }

  // Evaluate sec tag rejections against actual direction (called during training)
  void eval_sec_rejections(u64 actual_dir_bits) {
    for (u64 t = 0; t < NUM_TABLES; t++) {
      auto &tp = train_table_pred[t];
      if (!tp.tag_hit || tp.sec_hit) continue; // only care about tag_hit && !sec_hit
      auto record = [&](Counters &c) {
        c.sec_reject_count[t]++;
        // Check per-branch: was this table's prediction correct for each branch?
        bool all_correct = true;
        for (u64 r = 0; r < train_num_branch; r++) {
          bool pred_r = (tp.pred_bits >> r) & 1;
          bool actual_r = (actual_dir_bits >> r) & 1;
          if (pred_r != actual_r) { all_correct = false; break; }
        }
        if (all_correct) c.sec_reject_correct[t]++;
        else c.sec_reject_wrong[t]++;
      };
      record(cum);
      record(win);
    }
  }

  // Tag lookups per table
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

  // Multi-path: per-table path hit distribution (called during path mux)
  void record_path_hits(u64 table, const bool path_hits[NUM_PATHS]) {
    auto record = [&](Counters &c) {
      u64 hit_count = 0;
      for (u64 p = 0; p < NUM_PATHS; p++) {
        if (path_hits[p]) { c.path_hit[table][p]++; hit_count++; }
      }
      if (hit_count == 0) c.path_neither[table]++;
      if (hit_count > 1) c.path_both[table]++;
    };
    record(cum);
    record(win);
  }

  // Multi-path: which path was the provider and was it correct
  void record_provider_path(u64 table, u64 path, bool correct) {
    auto record = [&](Counters &c) {
      c.provider_from_path[table][path]++;
      if (correct) c.provider_path_correct[table][path]++;
    };
    record(cum);
    record(win);
  }

  // Multi-path: allocation went to which path
  void record_alloc_path(u64 table, u64 path) {
    auto record = [&](Counters &c) {
      c.alloc_to_path[table][path]++;
    };
    record(cum);
    record(win);
  }

  // Provider/alt resolution — called per branch during resolution
  void record_prediction(u64 rank, u64 match1_val, u64 match2_val,
                          bool meta_overrode, bool meta_chose_alt,
                          bool pred_taken) {
    u64 prov = decode_provider(match1_val);
    shadow_provider[rank] = prov;
    shadow_alt[rank] = decode_provider(match2_val);
    shadow_meta_overrode[rank] = meta_overrode;
    shadow_meta_chose_alt[rank] = meta_chose_alt;
    shadow_pred[rank] = pred_taken;
  }

  void record_outcome(u64 rank, bool actual_taken, bool mispredict) {
    bool correct = (shadow_pred[rank] == actual_taken);
    auto record = [&](Counters &c) {
      c.branches++;
      if (mispredict) c.mispredictions++;

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
    };
    record(cum);
    record(win);

    // Per-PC diagnostics (both cumulative and window)
    u64 pc = shadow_pc[rank];
    u64 prov = shadow_provider[rank];
    auto update_pc_diag = [&](PCDiag &d) {
      d.total++;
      if (mispredict) d.misp++;
      if (actual_taken) d.taken++;
      if (prov < NUM_TABLES) {
        d.tage_prov++;
        if (correct) d.tage_correct++;
      } else {
        d.fb_prov++;
        if (correct) d.fb_correct++;
      }
    };
    update_pc_diag(pc_diag[pc]);
    update_pc_diag(win_pc_diag[pc]);
    global_branch_count++;

    if (win.branches >= WINDOW_SIZE) {
      print_window(std::cerr);
      win.reset();
      win_pc_diag.clear();
      for (auto &b : win_entry_ever_hit) b.reset();
      window_num++;
    }
  }

  void record_branch_pc(u64 pc) { unique_branch_pcs.insert(pc); }

  // Enhanced: save PC per rank for per-PC diagnostics
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

  // Table occupancy — called on allocation writes
  void record_tage_write(u64 table, u64 index) {
    if (table < NUM_TABLES && index < MAX_TABLE_ENTRIES &&
        !tage_occupied[table].test(index)) {
      tage_occupied[table].set(index);
      tage_unique_entries[table]++;
    }
  }

  // Fallback write — track occupancy
  void record_fb_write(u64 index) {
    fb_writes++;
    if (index < FB_CAPACITY && !fb_occupied.test(index)) {
      fb_occupied.set(index);
      fb_unique_entries++;
    }
  }

  // Block PC pipeline shift (call at start of update_cycle)
  void shift_block_pc() {
    train_block_pc = current_block_pc;
    current_block_pc = shadow_block_pc;
    shift_table_pred_shadow();
  }

  // Feature 5: Record provider entry usage (call once per block during resolution)
  void record_provider_entry(u64 table, u64 index, u64 num_branches,
                              u64 num_correct) {
    if (table < NUM_TABLES && index < MAX_TABLE_ENTRIES) {
      auto &e = entry_life[table][index];
      if (e.valid) {
        e.pred_count += num_branches;
        e.correct_count += num_correct;
      }
    }
  }

  // Feature 7+9: Tag collision check (call per table during resolution)
  void record_collision_check(u64 table, u64 index, bool tag_hit) {
    if (!tag_hit || index >= MAX_TABLE_ENTRIES) return;
    entry_ever_hit[table].set(index);
    win_entry_ever_hit[table].set(index);
    auto record_coll = [&](Counters &c) {
      c.per_table_coll_checks[table]++;
      c.collision_checks++;
      auto &e = entry_life[table][index];
      if (e.valid && e.stored_pc != 0 && e.stored_pc != current_block_pc) {
        c.per_table_coll_hits[table]++;
        c.collision_hits++;
      }
    };
    record_coll(cum);
    record_coll(win);
  }

  // Feature 5+7+8: Record allocation (call per table on allocation write)
  void record_entry_alloc_diag(u64 table, u64 index) {
    if (index >= MAX_TABLE_ENTRIES) return;
    auto &e = entry_life[table][index];
    // Eviction stats for old entry
    if (e.valid) {
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
    // Feature 8: per-PC alloc tracking (both cumulative and window)
    pc_diag[train_block_pc].alloc_count++;
    win_pc_diag[train_block_pc].alloc_count++;
  }

  // u-bit write (per bank) — tracks shadow u-values for histogram
  void record_u_write(u64 table, u64 bank, u64 index, u64 new_val) {
    if (new_val > 0) {
      cum.u_set_count[table][bank]++;
      win.u_set_count[table][bank]++;
    } else {
      cum.u_clear_count[table][bank]++;
      win.u_clear_count[table][bank]++;
    }
    if (index < MAX_TABLE_ENTRIES)
      u_shadow[table][index][bank] = new_val;
  }

  // hyst training (per bank)
  void record_hyst_update(u64 table, u64 bank, bool incremented) {
    if (incremented) {
      cum.hyst_inc_count[table][bank]++;
      win.hyst_inc_count[table][bank]++;
    } else {
      cum.hyst_dec_count[table][bank]++;
      win.hyst_dec_count[table][bank]++;
    }
  }

  void record_decay_fire() {
    cum.decay_fire_count++;
    win.decay_fire_count++;
  }
  void record_epoch_reset() {
    cum.epoch_reset_count++;
    win.epoch_reset_count++;
    // Clear shadow u-values (mirrors hardware epoch reset)
    for (auto &table : u_shadow)
      for (auto &entry : table)
        entry.fill(0);
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

  // ======== Helpers ========

  static u64 decode_provider(u64 one_hot) {
    if (one_hot == 0) return NUM_TABLES; // fallback (bimodal or gshare)
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
      // Feature 5: entry table lifetime
      os << "Tlast_evict,Tlast_zu%,Tlast_avgp,";
      // Feature 7: collision rate
      os << "coll%,";
      // Feature 9: entry table utilization
      os << "Tlast_used,";
      // Feature 8: never-allocated PCs (window)
      os << "win_stuck,win_hard";
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
    // provider distribution
    os << pct(w.provider_count[NUM_TABLES], w.branches) << ",";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << pct(w.provider_count[i], w.branches) << ",";
    os << pct(w.alloc_success, w.alloc_attempts) << ","
       << (w.pressure_samples > 0
               ? double(w.acc_ctr_sum) / w.pressure_samples
               : 0)
       << ","
       << (w.pressure_samples > 0
               ? double(w.alloc_ctr_sum) / w.pressure_samples
               : 0)
       << ",";
    // Feature 5: entry table (last table) lifetime
    {
      u64 last = NUM_TABLES - 1;
      auto &ls = w.lifetime_stats[last];
      double avg_p = ls.evictions > 0 ? double(ls.total_preds) / ls.evictions : 0;
      os << ls.evictions << ","
         << pct(ls.zero_use, ls.evictions) << ","
         << std::setprecision(1) << avg_p << ",";
    }
    // Feature 7: collision rate
    os << pct(w.collision_hits, w.collision_checks) << ",";
    // Feature 9: entry table utilization
    os << win_entry_ever_hit[NUM_TABLES - 1].count() << ",";
    // Feature 8: window stuck + hard PCs
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
       << (c.blocks > 0 ? double(c.total_block_branches) / c.blocks : 0)
       << "\n";
    os << "  Instr/block histogram: ";
    for (u64 i = 1; i < MAX_BLOCK_INSTR && i <= 16; i++) {
      if (c.block_instr_hist[i] > 0)
        os << i << ":" << c.block_instr_hist[i] << " ";
    }
    os << "\n";
    os << "  Branches/block histogram: ";
    for (u64 i = 0; i < MAX_BLOCK_BR; i++) {
      if (c.block_br_hist[i] > 0)
        os << i << ":" << c.block_br_hist[i] << " ";
    }
    os << "\n";
    os << "  Entry point top-5: ";
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
       << std::setw(7)
       << pct(c.provider_correct[NUM_TABLES], c.provider_count[NUM_TABLES])
       << "%" << " |      --" << " |      --" << " |      --\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      os << "  T" << i << std::setw(8 - (i >= 10 ? 2 : 1)) << "" << "|"
         << std::setw(10) << c.provider_count[i] << " |" << std::setw(10)
         << c.provider_correct[i] << " |" << std::setw(7)
         << pct(c.provider_correct[i], c.provider_count[i]) << "%"
         << " |" << std::setw(8)
         << pct(c.tag_matches[i], c.tag_lookups[i]) << "%"
         << " |" << std::setw(8)
         << pct(c.sec_matches[i], c.tag_lookups[i]) << "%"
         << " |" << std::setw(8)
         << pct(c.full_matches[i], c.tag_lookups[i]) << "%\n";
    }

    // TAGE-only provider/alt distribution
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
      os << "  Overall TAGE accuracy: " << pct(tage_correct, tage_total)
         << "%\n";
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

    // Secondary tag rejection analysis
    {
      u64 total_reject = 0, total_correct = 0, total_wrong = 0;
      for (u64 i = 0; i < NUM_TABLES; i++) {
        total_reject += c.sec_reject_count[i];
        total_correct += c.sec_reject_correct[i];
        total_wrong += c.sec_reject_wrong[i];
      }
      if (total_reject > 0) {
        os << "\nSec Tag Rejection Analysis (tag hit but sec miss — would pred have been correct?):\n";
        os << "  Table  | Rejected  | WouldCorrect | WouldWrong | Correct%\n";
        os << "  -------+-----------+--------------+------------+---------\n";
        for (u64 i = 0; i < NUM_TABLES; i++) {
          if (c.sec_reject_count[i] == 0) continue;
          os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
             << std::setw(10) << c.sec_reject_count[i] << " |"
             << std::setw(13) << c.sec_reject_correct[i] << " |"
             << std::setw(11) << c.sec_reject_wrong[i] << " |"
             << std::setw(7) << pct(c.sec_reject_correct[i], c.sec_reject_count[i]) << "%\n";
        }
        os << "  Total  |" << std::setw(10) << total_reject << " |"
           << std::setw(13) << total_correct << " |"
           << std::setw(11) << total_wrong << " |"
           << std::setw(7) << pct(total_correct, total_reject) << "%\n";
      }
    }

    // Multi-path sec tag banking
    if constexpr (NUM_PATHS > 1) {
      os << "\n--- Multi-Path Sec Tag Banking (NUM_PATHS=" << NUM_PATHS << ") ---\n";
      os << "  Table  |";
      for (u64 p = 0; p < NUM_PATHS; p++) os << "   P" << p << "hit ";
      os << "| Neither | Both    | P0prov  | P1prov  | P0acc   | P1acc   | P0alloc | P1alloc\n";
      os << "  -------+";
      for (u64 p = 0; p < NUM_PATHS; p++) os << "---------";
      os << "+---------+---------+---------+---------+---------+---------+---------+--------\n";
      for (u64 i = 0; i < NUM_TABLES; i++) {
        u64 lookups = c.tag_lookups[i];
        os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|";
        for (u64 p = 0; p < NUM_PATHS; p++)
          os << std::setw(7) << pct(c.path_hit[i][p], lookups) << "% ";
        os << "|" << std::setw(7) << pct(c.path_neither[i], lookups) << "%"
           << " |" << std::setw(7) << pct(c.path_both[i], lookups) << "%";
        // Provider path distribution + accuracy
        u64 total_prov = 0;
        for (u64 p = 0; p < NUM_PATHS; p++) total_prov += c.provider_from_path[i][p];
        for (u64 p = 0; p < std::min(NUM_PATHS, u64(2)); p++)
          os << " |" << std::setw(8) << c.provider_from_path[i][p];
        for (u64 p = 0; p < std::min(NUM_PATHS, u64(2)); p++)
          os << " |" << std::setw(6) << pct(c.provider_path_correct[i][p],
                                             c.provider_from_path[i][p]) << "%";
        for (u64 p = 0; p < std::min(NUM_PATHS, u64(2)); p++)
          os << " |" << std::setw(8) << c.alloc_to_path[i][p];
        os << "\n";
      }
      // Totals
      u64 t_lookups = 0;
      std::array<u64, NUM_PATHS> t_phit{}, t_pprov{}, t_pcorr{}, t_palloc{};
      u64 t_neither = 0, t_both = 0;
      for (u64 i = 0; i < NUM_TABLES; i++) {
        t_lookups += c.tag_lookups[i];
        t_neither += c.path_neither[i];
        t_both += c.path_both[i];
        for (u64 p = 0; p < NUM_PATHS; p++) {
          t_phit[p] += c.path_hit[i][p];
          t_pprov[p] += c.provider_from_path[i][p];
          t_pcorr[p] += c.provider_path_correct[i][p];
          t_palloc[p] += c.alloc_to_path[i][p];
        }
      }
      os << "  Total  |";
      for (u64 p = 0; p < NUM_PATHS; p++)
        os << std::setw(7) << pct(t_phit[p], t_lookups) << "% ";
      os << "|" << std::setw(7) << pct(t_neither, t_lookups) << "%"
         << " |" << std::setw(7) << pct(t_both, t_lookups) << "%";
      for (u64 p = 0; p < std::min(NUM_PATHS, u64(2)); p++)
        os << " |" << std::setw(8) << t_pprov[p];
      for (u64 p = 0; p < std::min(NUM_PATHS, u64(2)); p++)
        os << " |" << std::setw(6) << pct(t_pcorr[p], t_pprov[p]) << "%";
      for (u64 p = 0; p < std::min(NUM_PATHS, u64(2)); p++)
        os << " |" << std::setw(8) << t_palloc[p];
      os << "\n";
    }

    // Allocation
    os << "\nAllocation:\n";
    os << "  Attempts: " << c.alloc_attempts
       << "  Success: " << c.alloc_success << " ("
       << pct(c.alloc_success, c.alloc_attempts) << "%)"
       << "  Fail: " << c.alloc_fail
       << "  Blocked: " << c.alloc_blocked << "\n";
    os << "  Per table:";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << " T" << i << "=" << c.alloc_per_table[i];
    os << "\n";
    os << "  Unique branch PCs: " << unique_branch_pcs.size() << "\n";

    // Allocation cascade
    os << "\nAllocation Cascade (provider -> target):\n";
    os << "  Provider  | Allocs  |";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << std::setw(7) << "T" + std::to_string(i);
    os << "\n  ----------+---------+";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << "-------";
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
      u64 ever_hit = entry_ever_hit[i].count();
      os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
         << std::setw(6) << MAX_TABLE_ENTRIES << " |" << std::setw(9)
         << tage_unique_entries[i] << " |" << std::setw(5)
         << pct(tage_unique_entries[i], MAX_TABLE_ENTRIES) << "%"
         << " |" << std::setw(8) << ever_hit
         << " |" << std::setw(5)
         << pct(ever_hit, MAX_TABLE_ENTRIES) << "%"
         << " |" << std::setw(7)
         << pct(c.provider_correct[i], c.provider_count[i]) << "%"
         << " |" << std::setw(7) << c.alloc_per_table[i] << "\n";
    }

    // Fallback occupancy
    os << "\n" << FB_NAME << " Fallback (" << FB_CAPACITY << " entries):\n";
    os << "  Writes: " << fb_writes
       << "  Occupied: " << fb_unique_entries
       << " (" << pct(fb_unique_entries, FB_CAPACITY) << "%)\n";
    os << "  Accuracy: " << c.provider_correct[NUM_TABLES] << "/"
       << c.provider_count[NUM_TABLES] << " ("
       << pct(c.provider_correct[NUM_TABLES], c.provider_count[NUM_TABLES]) << "%)\n";

    // u-bit
    os << "\nU-bit:\n";
    os << "  Table  | U set    | U clear  | Turnover\n";
    os << "  -------+----------+----------+----------\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      u64 total_set = 0, total_clear = 0;
      for (u64 b = 0; b < U_BANKS; b++) {
        total_set += c.u_set_count[i][b];
        total_clear += c.u_clear_count[i][b];
      }
      u64 total_u = total_set + total_clear;
      os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
         << std::setw(9) << total_set << " |" << std::setw(9)
         << total_clear << " |" << std::setw(7)
         << pct(total_clear, total_u) << "%\n";
    }
    os << "  Decay fires: " << c.decay_fire_count
       << "  Epoch resets: " << c.epoch_reset_count << "\n";

    // Per-bank U-bit breakdown (only when banked)
    if constexpr (U_BANKS > 1) {
      os << "\nU-bit Per Bank:\n";
      os << "  Table  | Bank | U set    | U clear  | Turnover\n";
      os << "  -------+------+----------+----------+----------\n";
      for (u64 i = 0; i < NUM_TABLES; i++) {
        for (u64 b = 0; b < U_BANKS; b++) {
          u64 total_u = c.u_set_count[i][b] + c.u_clear_count[i][b];
          os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
             << std::setw(5) << b << " |" << std::setw(9) << c.u_set_count[i][b]
             << " |" << std::setw(9) << c.u_clear_count[i][b]
             << " |" << std::setw(7) << pct(c.u_clear_count[i][b], total_u) << "%\n";
        }
      }
    }

    // U-value distribution histogram (snapshot at end of trace)
    {
      // Find max u value present to size the histogram columns
      u64 max_seen = 0;
      for (u64 i = 0; i < NUM_TABLES; i++)
        for (u64 e = 0; e < MAX_TABLE_ENTRIES; e++)
          for (u64 b = 0; b < U_BANKS; b++)
            if (u_shadow[i][e][b] > max_seen) max_seen = u_shadow[i][e][b];
      u64 num_buckets = std::min(max_seen + 1, MAX_U_VAL);

      os << "\nU-value Distribution (end-of-trace snapshot):\n";
      os << "  Table  |";
      for (u64 v = 0; v < num_buckets; v++)
        os << " u=" << v << "   ";
      os << "\n  -------+";
      for (u64 v = 0; v < num_buckets; v++)
        os << "--------";
      os << "\n";
      for (u64 i = 0; i < NUM_TABLES; i++) {
        std::array<u64, MAX_U_VAL> hist{};
        u64 total = 0;
        for (u64 e = 0; e < MAX_TABLE_ENTRIES; e++) {
          for (u64 b = 0; b < U_BANKS; b++) {
            u64 v = u_shadow[i][e][b];
            if (v < MAX_U_VAL) hist[v]++;
            total++;
          }
        }
        os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|";
        for (u64 v = 0; v < num_buckets; v++)
          os << std::setw(5) << hist[v] << "(" << std::setw(2)
             << (total > 0 ? u64(100.0 * hist[v] / total + 0.5) : 0) << "%)";
        os << "\n";
      }
    }

    // Hyst training breakdown
    os << "\nHyst Training:\n";
    os << "  Table  | Hyst Inc | Hyst Dec | Inc%\n";
    os << "  -------+----------+----------+------\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      u64 total_inc = 0, total_dec = 0;
      for (u64 b = 0; b < HYST_BANKS; b++) {
        total_inc += c.hyst_inc_count[i][b];
        total_dec += c.hyst_dec_count[i][b];
      }
      u64 total_h = total_inc + total_dec;
      os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
         << std::setw(9) << total_inc << " |" << std::setw(9)
         << total_dec << " |" << std::setw(5)
         << pct(total_inc, total_h) << "%\n";
    }

    // Per-bank hyst breakdown (only when banked)
    if constexpr (HYST_BANKS > 1) {
      os << "\nHyst Training Per Bank:\n";
      os << "  Table  | Bank | Hyst Inc | Hyst Dec | Inc%\n";
      os << "  -------+------+----------+----------+------\n";
      for (u64 i = 0; i < NUM_TABLES; i++) {
        for (u64 b = 0; b < HYST_BANKS; b++) {
          u64 total_h = c.hyst_inc_count[i][b] + c.hyst_dec_count[i][b];
          os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
             << std::setw(5) << b << " |" << std::setw(9) << c.hyst_inc_count[i][b]
             << " |" << std::setw(9) << c.hyst_dec_count[i][b]
             << " |" << std::setw(5) << pct(c.hyst_inc_count[i][b], total_h) << "%\n";
        }
      }
    }

    // Pressure
    if (c.pressure_samples > 0) {
      os << "\nPressure Counters:\n";
      os << "  Avg acc_ctr: "
         << double(c.acc_ctr_sum) / c.pressure_samples << "\n";
      os << "  Avg alloc_ctr: "
         << double(c.alloc_ctr_sum) / c.pressure_samples << "\n";
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
      u64 easy_count = 0, hard_count = 0; // >90% biased = easy
      for (auto &[pc, d] : pc_diag) {
        if (d.total < 5) continue; // skip rare branches
        double bias = std::max(pct(d.taken, d.total),
                               100.0 - pct(d.taken, d.total));
        u64 bin = std::min(u64(bias / 10.0), BINS - 1);
        bias_hist[bin]++;
        if (bias >= 90.0) easy_count++;
        else if (bias < 60.0) hard_count++;
      }
      os << "  Bias (max(T%,NT%)) | PCs\n";
      os << "  -------------------+------\n";
      for (u64 i = 0; i < BINS; i++)
        os << "  " << std::setw(3) << i * 10 << "-" << std::setw(3)
           << (i + 1) * 10 << "%            |" << std::setw(6) << bias_hist[i]
           << "\n";
      os << "  Easy (>=90% biased): " << easy_count
         << "  Hard (<60% biased): " << hard_count << "\n";
    }

    // Feature 3: TAGE vs Fallback hit rate
    {
      os << "\n--- Feature 3: TAGE vs Fallback Per-PC ---\n";
      u64 tage_dom = 0, fb_dominant = 0, mixed = 0;
      u64 tage_dom_misp = 0, fb_dominant_misp = 0, mixed_misp = 0;
      for (auto &[pc, d] : pc_diag) {
        if (d.total < 5) continue;
        double tage_frac = pct(d.tage_prov, d.total);
        if (tage_frac >= 50.0) {
          tage_dom++;
          tage_dom_misp += d.misp;
        } else if (tage_frac <= 10.0) {
          fb_dominant++;
          fb_dominant_misp += d.misp;
        } else {
          mixed++;
          mixed_misp += d.misp;
        }
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
        double avg_preds = ls.evictions > 0
                               ? double(ls.total_preds) / ls.evictions
                               : 0;
        double avg_acc = ls.total_preds > 0
                             ? pct(ls.total_correct, ls.total_preds)
                             : 0;
        os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << "" << "|"
           << std::setw(10) << ls.evictions << " |" << std::setw(8)
           << std::setprecision(1) << avg_preds << " |" << std::setw(7)
           << std::setprecision(1) << avg_acc << "% |" << std::setw(8)
           << ls.zero_use << " |" << std::setw(7)
           << pct(ls.zero_use, ls.evictions) << "%\n";
      }
    }

    // Feature 7: Tag collision rate
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

    // Feature 8: Gshare-stuck branches
    {
      os << "\n--- Feature 8: Fallback-Stuck Branches ---\n";
      u64 never_alloc = 0, alloc_but_fb = 0, alloc_and_tage = 0;
      u64 never_alloc_misp = 0, alloc_but_fb_misp = 0, alloc_tage_misp = 0;
      for (auto &[pc, d] : pc_diag) {
        if (d.total < 5) continue;
        if (d.alloc_count == 0) {
          never_alloc++;
          never_alloc_misp += d.misp;
        } else if (d.tage_prov * 2 < d.total) {
          alloc_but_fb++;
          alloc_but_fb_misp += d.misp;
        } else {
          alloc_and_tage++;
          alloc_tage_misp += d.misp;
        }
      }
      os << "  Category                     | PCs    | Mispred  | AvgAlloc\n";
      os << "  -----------------------------+--------+----------+---------\n";
      os << "  Never allocated to TAGE      |" << std::setw(7) << never_alloc
         << " |" << std::setw(9) << never_alloc_misp << " |      --\n";
      os << "  Allocated but FB dominates   |" << std::setw(7) << alloc_but_fb
         << " |" << std::setw(9) << alloc_but_fb_misp << " |";
      {
        u64 total_alloc = 0, count = 0;
        for (auto &[pc, d] : pc_diag) {
          if (d.total >= 5 && d.alloc_count > 0 && d.tage_prov * 2 < d.total) {
            total_alloc += d.alloc_count;
            count++;
          }
        }
        os << std::setw(8) << (count > 0 ? double(total_alloc) / count : 0)
           << "\n";
      }
      os << "  Allocated and TAGE provides  |" << std::setw(7)
         << alloc_and_tage << " |" << std::setw(9) << alloc_tage_misp << " |";
      {
        u64 total_alloc = 0, count = 0;
        for (auto &[pc, d] : pc_diag) {
          if (d.total >= 5 && d.alloc_count > 0 &&
              d.tage_prov * 2 >= d.total) {
            total_alloc += d.alloc_count;
            count++;
          }
        }
        os << std::setw(8) << (count > 0 ? double(total_alloc) / count : 0)
           << "\n";
      }
    }

    os << "\n=== End TageAhead Monitor ===\n\n";
  }

private:
  template <size_t SZ>
  static void print_top_n(std::ostream &os,
                           const std::array<u64, SZ> &hist, u64 n) {
    std::array<std::pair<u64, u64>, 64> top{};
    u64 found = 0;
    for (u64 i = 0; i < SZ && found < 64; i++) {
      if (hist[i] > 0) top[found++] = {hist[i], i};
    }
    for (u64 i = 0; i < std::min(n, found); i++) {
      for (u64 j = i + 1; j < found; j++) {
        if (top[j].first > top[i].first) std::swap(top[i], top[j]);
      }
      os << top[i].second << ":" << top[i].first << " ";
    }
  }
};
