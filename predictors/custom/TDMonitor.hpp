#pragma once
// TageDirect Performance Monitor — software-only instrumentation
// Gated by -DTAGE_MONITOR -DCHEATING_MODE at compile time.
// All counters are plain C++ types (zero HARCOM cost).

#include <cstdint>
#include <array>
#include <bitset>
#include <iostream>
#include <iomanip>
#include <unordered_set>

using u64 = uint64_t;

template <u64 NUM_TABLES, u64 LANES, u64 P1_ENTRIES, u64 RWRAM_BANKS,
          u64 MAX_TABLE_ENTRIES = 2048>
struct TDMonitor {

  static constexpr u64 WINDOW_SIZE = 100000; // branches per window

  // Max values for histograms
  static constexpr u64 MAX_BLOCK_INSTR = 64; // histogram bins for instr/block
  static constexpr u64 MAX_BLOCK_BR = LANES + 1;

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
    u64 total_block_instr = 0;     // sum of instructions across all blocks
    u64 total_block_branches = 0;  // sum of branches across all blocks
    std::array<u64, MAX_BLOCK_INSTR> block_instr_hist{}; // histogram: instr/block
    std::array<u64, MAX_BLOCK_BR> block_br_hist{};       // histogram: branches/block
    std::array<u64, 1024> entry_point_hist{};   // block_entry distribution (capped at 1024)
    std::array<u64, 1024> exit_point_hist{};    // block_entry+block_size distribution

    // Provider distribution: slot 0..NUM_TABLES-1 = TAGE tables, NUM_TABLES = P1 fallback
    std::array<u64, NUM_TABLES + 1> provider_count{};
    std::array<u64, NUM_TABLES + 1> provider_correct{};

    // Alt provider
    std::array<u64, NUM_TABLES + 1> alt_count{};

    // Meta override
    u64 meta_override_count = 0;
    u64 meta_override_correct = 0;
    u64 meta_chose_alt = 0;       // meta selected alt over primary
    u64 meta_chose_alt_correct = 0;
    u64 meta_chose_pri = 0;       // meta selected primary over alt
    u64 meta_chose_pri_correct = 0;

    // Tag match rate per table
    std::array<u64, NUM_TABLES> tag_lookups{};
    std::array<u64, NUM_TABLES> tag_matches{};

    // P1 gshare
    u64 p1_predictions = 0;
    u64 p1_correct = 0;
    u64 p1_writes = 0;
    std::array<u64, LANES> p1_lane_predictions{};
    std::array<u64, LANES> p1_lane_correct{};

    // P1 vs P2 disagreement
    u64 p1p2_disagree = 0;
    u64 p1_right_p2_wrong = 0;
    u64 p2_right_p1_wrong = 0;

    // Allocation
    u64 alloc_attempts = 0;
    u64 alloc_success = 0;
    u64 alloc_fail = 0;
    u64 alloc_blocked = 0;  // mispredict but postmask==0 (P1 was provider)
    std::array<u64, NUM_TABLES> alloc_per_table{};
    // cascade: alloc_from_provider[i] = allocations when Ti was provider
    // alloc_from_provider[NUM_TABLES] = allocations when P1 was provider
    std::array<u64, NUM_TABLES + 1> alloc_from_provider{};
    std::array<std::array<u64, NUM_TABLES>, NUM_TABLES + 1> alloc_cascade{}; // [provider][target]

    // u-bit
    std::array<u64, NUM_TABLES> u_set_count{};
    std::array<u64, NUM_TABLES> u_clear_count{};
    u64 decay_fire_count = 0;
    u64 epoch_resets = 0;
    u64 uctr_sum = 0;
    u64 uctr_samples = 0;

    // rwram banking (hyst + u combined per table)
    std::array<u64, NUM_TABLES> rwram_reads{};
    std::array<u64, NUM_TABLES> rwram_writes{};
    std::array<u64, NUM_TABLES> rwram_writes_direct{};
    std::array<u64, NUM_TABLES> rwram_writes_buffered{};
    std::array<u64, NUM_TABLES> rwram_same_bank_conflict{};
    std::array<u64, NUM_TABLES> rwram_buffer_pending{};

    void reset() { *this = Counters{}; }
  };

  Counters cum;
  Counters win;
  u64 window_num = 0;
  bool header_printed = false;

  // P1 occupancy — 1 bit per entry per lane
  std::array<std::bitset<P1_ENTRIES>, LANES> p1_occupied{};
  std::array<u64, LANES> p1_unique_entries{};

  // TAGE table occupancy — 1 bit per entry per table
  std::array<std::bitset<MAX_TABLE_ENTRIES>, NUM_TABLES> tage_occupied{};
  std::array<u64, NUM_TABLES> tage_unique_entries{};

  // Unique branch PCs
  std::unordered_set<u64> unique_branch_pcs;

  // ======== Shadow state (set in predict2, consumed in update_cycle) ========
  static constexpr u64 MAX_FW = 64;
  std::array<u64, MAX_FW> shadow_provider{};
  std::array<u64, MAX_FW> shadow_alt{};
  std::array<bool, MAX_FW> shadow_meta_overrode{};
  std::array<bool, MAX_FW> shadow_meta_chose_alt{}; // true = meta picked alt, false = primary
  std::array<bool, MAX_FW> shadow_p1_pred{};
  std::array<bool, MAX_FW> shadow_p2_pred{};

  // ======== Recording methods ========

  void record_predict1() { cum.predict1_calls++; win.predict1_calls++; }
  void record_reuse_predict1() { cum.reuse_predict1_calls++; win.reuse_predict1_calls++; }
  void record_predict2() { cum.predict2_calls++; win.predict2_calls++; }
  void record_reuse_predict2() { cum.reuse_predict2_calls++; win.reuse_predict2_calls++; }

  // Block structure — called at end of block (update_cycle)
  void record_block(u64 block_entry, u64 block_size, u64 num_branch, bool extra_cycle_fired) {
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

      if (block_entry < 1024) c.entry_point_hist[block_entry]++;
      u64 exit_pt = block_entry + block_size;
      if (exit_pt < 1024) c.exit_point_hist[exit_pt]++;
    };
    record(cum);
    record(win);
  }

  // Prediction — called in predict2 per rank
  void record_prediction(u64 rank, u64 match1_val, u64 match2_val,
                          bool meta_overrode, bool meta_chose_alt,
                          bool p1_pred, bool p2_pred) {
    u64 prov = decode_provider(match1_val);
    shadow_provider[rank] = prov;
    shadow_alt[rank] = decode_provider(match2_val);
    shadow_meta_overrode[rank] = meta_overrode;
    shadow_meta_chose_alt[rank] = meta_chose_alt;
    shadow_p1_pred[rank] = p1_pred;
    shadow_p2_pred[rank] = p2_pred;
  }

  // Tag match — called per table in predict2
  void record_tag_lookup(u64 table, bool matched) {
    cum.tag_lookups[table]++; win.tag_lookups[table]++;
    if (matched) { cum.tag_matches[table]++; win.tag_matches[table]++; }
  }

  // Branch outcome — called per branch in update_cycle monitor block
  // mispredict: simulator-level (only true for block-ending misprediction)
  void record_outcome(u64 rank, bool actual_taken, bool mispredict) {
    // Per-branch correctness derived from actual vs predicted
    bool p2_ok = (shadow_p2_pred[rank] == actual_taken);
    auto record = [&](Counters &c) {
      c.branches++;
      if (mispredict) c.mispredictions++;

      u64 prov = shadow_provider[rank];
      c.provider_count[prov]++;
      if (p2_ok) c.provider_correct[prov]++;

      u64 alt = shadow_alt[rank];
      c.alt_count[alt]++;

      if (shadow_meta_overrode[rank]) {
        c.meta_override_count++;
        if (p2_ok) c.meta_override_correct++;
        if (shadow_meta_chose_alt[rank]) {
          c.meta_chose_alt++;
          if (p2_ok) c.meta_chose_alt_correct++;
        } else {
          c.meta_chose_pri++;
          if (p2_ok) c.meta_chose_pri_correct++;
        }
      }

      // P1
      c.p1_predictions++;
      bool p1_ok = (shadow_p1_pred[rank] == actual_taken);
      if (p1_ok) c.p1_correct++;
      if (rank < LANES) {
        c.p1_lane_predictions[rank]++;
        if (p1_ok) c.p1_lane_correct[rank]++;
      }

      // P1 vs P2
      if (shadow_p1_pred[rank] != shadow_p2_pred[rank]) {
        c.p1p2_disagree++;
        if (p1_ok && !p2_ok) c.p1_right_p2_wrong++;
        if (p2_ok && !p1_ok) c.p2_right_p1_wrong++;
      }
    };
    record(cum);
    record(win);

    if (win.branches >= WINDOW_SIZE) {
      print_window(std::cerr);
      win.reset();
      window_num++;
    }
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

  // Allocation blocked (mispredict/trigger fired but postmask==0)
  void record_alloc_blocked() {
    cum.alloc_blocked++; win.alloc_blocked++;
  }

  // Cascade: provider_idx = which table was provider (NUM_TABLES = P1),
  //          allocate_mask = which tables got allocated
  void record_alloc_cascade(u64 provider_idx, u64 allocate_mask) {
    auto record = [&](Counters &c) {
      if (allocate_mask == 0) return;
      c.alloc_from_provider[provider_idx]++;
      for (u64 i = 0; i < NUM_TABLES; i++)
        if (allocate_mask & (u64(1) << i)) c.alloc_cascade[provider_idx][i]++;
    };
    record(cum);
    record(win);
  }

  // Branch PC tracking
  void record_branch_pc(u64 pc) {
    unique_branch_pcs.insert(pc);
  }

  // u-bit write
  void record_u_write(u64 table, bool new_u) {
    if (new_u) { cum.u_set_count[table]++; win.u_set_count[table]++; }
    else { cum.u_clear_count[table]++; win.u_clear_count[table]++; }
  }

  void record_decay_fire() { cum.decay_fire_count++; win.decay_fire_count++; }
  void record_epoch_reset() { cum.epoch_resets++; win.epoch_resets++; }
  void record_uctr(u64 val) {
    cum.uctr_sum += val; cum.uctr_samples++;
    win.uctr_sum += val; win.uctr_samples++;
  }

  // P1 write — track per-lane occupancy
  void record_p1_write(u64 lane, u64 index) {
    cum.p1_writes++; win.p1_writes++;
    if (lane < LANES && index < P1_ENTRIES && !p1_occupied[lane].test(index)) {
      p1_occupied[lane].set(index);
      p1_unique_entries[lane]++;
    }
  }

  // TAGE table write — track occupancy
  void record_tage_write(u64 table, u64 index) {
    if (table < NUM_TABLES && index < MAX_TABLE_ENTRIES && !tage_occupied[table].test(index)) {
      tage_occupied[table].set(index);
      tage_unique_entries[table]++;
    }
  }

  // rwram banking
  void record_rwram_read(u64 table) {
    cum.rwram_reads[table]++; win.rwram_reads[table]++;
  }

  void record_rwram_write(u64 table, u64 bank_id, bool direct,
                          u64 read_bank_mask, u64 write_bank_mask) {
    cum.rwram_writes[table]++; win.rwram_writes[table]++;
    if (direct) {
      cum.rwram_writes_direct[table]++; win.rwram_writes_direct[table]++;
    } else {
      cum.rwram_writes_buffered[table]++; win.rwram_writes_buffered[table]++;
    }
    if ((read_bank_mask >> bank_id) & 1) {
      cum.rwram_same_bank_conflict[table]++; win.rwram_same_bank_conflict[table]++;
    }
    if (write_bank_mask != 0) {
      cum.rwram_buffer_pending[table]++; win.rwram_buffer_pending[table]++;
    }
  }

  // ======== Helpers ========

  static u64 decode_provider(u64 one_hot) {
    if (one_hot == 0) return NUM_TABLES; // P1 fallback
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
      os << "# win,br,misp%,extra%,i/blk,br/blk,";
      os << "p1%,";
      for (u64 i = 0; i < NUM_TABLES; i++) os << "T" << i << "%,";
      os << "p1acc%,p1p2dis%,alloc_ok%,uctr_avg";
      os << "\n";
      header_printed = true;
    }
    auto &w = win;
    os << std::fixed << std::setprecision(1);
    os << window_num << ","
       << w.branches << ","
       << pct(w.mispredictions, w.branches) << ","
       << pct(w.extra_cycles, w.blocks) << ","
       << (w.blocks > 0 ? double(w.total_block_instr) / w.blocks : 0) << ","
       << (w.blocks > 0 ? double(w.total_block_branches) / w.blocks : 0) << ",";
    // provider distribution
    os << pct(w.provider_count[NUM_TABLES], w.branches) << ",";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << pct(w.provider_count[i], w.branches) << ",";
    os << pct(w.p1_correct, w.p1_predictions) << ","
       << pct(w.p1p2_disagree, w.branches) << ","
       << pct(w.alloc_success, w.alloc_attempts) << ","
       << (w.uctr_samples > 0 ? double(w.uctr_sum) / w.uctr_samples : 0);
    os << "\n";
  }

  // ======== End-of-trace summary ========

  void print_summary(std::ostream &os = std::cerr) const {
    const auto &c = cum;

    os << "\n=== TageDirect Monitor Summary ===\n";
    os << "Branches: " << c.branches
       << "  Mispredictions: " << c.mispredictions
       << " (" << std::fixed << std::setprecision(2)
       << pct(c.mispredictions, c.branches) << "%)\n";
    os << "Blocks: " << c.blocks
       << "  Extra cycles: " << c.extra_cycles
       << " (" << pct(c.extra_cycles, c.blocks) << "%)\n";

    // Pipeline calls
    os << "\nPipeline Calls:\n";
    os << "  predict1: " << c.predict1_calls
       << "  reuse_predict1: " << c.reuse_predict1_calls
       << "  (reuse rate: " << pct(c.reuse_predict1_calls,
                                    c.predict1_calls + c.reuse_predict1_calls) << "%)\n";
    os << "  predict2: " << c.predict2_calls
       << "  reuse_predict2: " << c.reuse_predict2_calls
       << "  (reuse rate: " << pct(c.reuse_predict2_calls,
                                    c.predict2_calls + c.reuse_predict2_calls) << "%)\n";
    u64 total_p1 = c.predict1_calls + c.reuse_predict1_calls;
    os << "  Avg instr/block (from P1): "
       << (c.predict1_calls > 0 ? double(total_p1) / c.predict1_calls : 0) << "\n\n";

    // Block structure
    os << "Block Structure:\n";
    os << "  Avg instr/block: " << std::setprecision(1)
       << (c.blocks > 0 ? double(c.total_block_instr) / c.blocks : 0) << "\n";
    os << "  Avg branches/block: "
       << (c.blocks > 0 ? double(c.total_block_branches) / c.blocks : 0) << "\n";

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

    // Entry/exit point top-5
    os << "  Entry point top-5: ";
    print_top_n(os, c.entry_point_hist, 5);
    os << "\n  Exit point top-5: ";
    print_top_n(os, c.exit_point_hist, 5);
    os << "\n\n";

    // Provider distribution
    os << "Provider Distribution:\n";
    os << "  Table     | Provided  |  Correct  | Accuracy | TagMatch%\n";
    os << "  ----------+-----------+-----------+----------+----------\n";
    os << "  P1(gshr)  |" << std::setw(10) << c.provider_count[NUM_TABLES]
       << " |" << std::setw(10) << c.provider_correct[NUM_TABLES]
       << " |" << std::setw(7) << std::setprecision(1)
       << pct(c.provider_correct[NUM_TABLES], c.provider_count[NUM_TABLES]) << "%"
       << " |      --\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      os << "  T" << i << std::setw(8 - (i >= 10 ? 2 : 1)) << ""
         << "|" << std::setw(10) << c.provider_count[i]
         << " |" << std::setw(10) << c.provider_correct[i]
         << " |" << std::setw(7) << pct(c.provider_correct[i], c.provider_count[i]) << "%"
         << " |" << std::setw(7) << pct(c.tag_matches[i], c.tag_lookups[i]) << "%\n";
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
        os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << ""
           << "|" << std::setw(6) << pct(c.provider_count[i], tage_total) << "%"
           << " |" << std::setw(6) << pct(c.alt_count[i], tage_total) << "%"
           << " |" << std::setw(6) << pct(c.provider_correct[i], c.provider_count[i]) << "%\n";
      }
      os << "  Overall TAGE accuracy: " << pct(tage_correct, tage_total) << "%\n";
    }

    // Meta
    os << "\nMeta Override:\n";
    os << "  Total: " << c.meta_override_count
       << "  Correct: " << c.meta_override_correct
       << " (" << pct(c.meta_override_correct, c.meta_override_count) << "%)\n";
    os << "  Chose alt:     " << c.meta_chose_alt
       << "  Correct: " << c.meta_chose_alt_correct
       << " (" << pct(c.meta_chose_alt_correct, c.meta_chose_alt) << "%)\n";
    os << "  Chose primary: " << c.meta_chose_pri
       << "  Correct: " << c.meta_chose_pri_correct
       << " (" << pct(c.meta_chose_pri_correct, c.meta_chose_pri) << "%)\n";

    // P1
    os << "\nP1 Gshare:\n";
    os << "  Accuracy: " << c.p1_correct << "/" << c.p1_predictions
       << " (" << pct(c.p1_correct, c.p1_predictions) << "%)\n";
    os << "  Writes: " << c.p1_writes << "\n";
    os << "  Per-lane stats (" << P1_ENTRIES << " entries/lane):\n";
    os << "  Lane | Predictions |  Correct  | Accuracy | Occupied |  Occ%\n";
    os << "  -----+-------------+-----------+----------+----------+------\n";
    u64 total_occ = 0;
    for (u64 i = 0; i < LANES; i++) {
      os << "  " << std::setw(4) << i
         << " |" << std::setw(12) << c.p1_lane_predictions[i]
         << " |" << std::setw(10) << c.p1_lane_correct[i]
         << " |" << std::setw(7) << pct(c.p1_lane_correct[i], c.p1_lane_predictions[i]) << "%"
         << " |" << std::setw(9) << p1_unique_entries[i]
         << " |" << std::setw(5) << pct(p1_unique_entries[i], P1_ENTRIES) << "%\n";
      total_occ += p1_unique_entries[i];
    }
    os << "  Total occupancy: " << total_occ << "/" << P1_ENTRIES * LANES
       << " (" << pct(total_occ, P1_ENTRIES * LANES) << "%)\n";

    // P1 vs P2
    os << "\nP1 vs P2:\n";
    os << "  Disagree: " << c.p1p2_disagree
       << " (" << pct(c.p1p2_disagree, c.branches) << "%)\n";
    os << "  P1 right, P2 wrong: " << c.p1_right_p2_wrong
       << "  P2 right, P1 wrong: " << c.p2_right_p1_wrong << "\n";

    // Allocation
    os << "\nAllocation:\n";
    os << "  Attempts: " << c.alloc_attempts
       << "  Success: " << c.alloc_success
       << " (" << pct(c.alloc_success, c.alloc_attempts) << "%)"
       << "  Fail: " << c.alloc_fail
       << "  Blocked: " << c.alloc_blocked << "\n";
    os << "  Per table:";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << " T" << i << "=" << c.alloc_per_table[i];
    os << "\n";
    os << "  Unique branch PCs: " << unique_branch_pcs.size() << "\n";

    // Allocation cascade (provider → target)
    os << "\nAllocation Cascade (provider → target):\n";
    os << "  Provider  | Allocs  |";
    for (u64 i = 0; i < NUM_TABLES; i++) os << std::setw(7) << "T" + std::to_string(i);
    os << "\n  ----------+---------+";
    for (u64 i = 0; i < NUM_TABLES; i++) os << "-------";
    os << "\n";
    for (u64 p = 0; p <= NUM_TABLES; p++) {
      if (c.alloc_from_provider[p] == 0) continue;
      if (p == NUM_TABLES)
        os << "  P1(gshr)  |";
      else
        os << "  T" << p << std::setw(8 - (p >= 10 ? 2 : 1)) << "" << "|";
      os << std::setw(8) << c.alloc_from_provider[p] << " |";
      for (u64 t = 0; t < NUM_TABLES; t++)
        os << std::setw(7) << c.alloc_cascade[p][t];
      os << "\n";
    }

    // TAGE table occupancy
    os << "\nTAGE Table Occupancy:\n";
    os << "  Table  | Size  | Occupied |  Occ% | Accuracy | Allocs\n";
    os << "  -------+-------+----------+-------+----------+-------\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << ""
         << "|" << std::setw(6) << MAX_TABLE_ENTRIES
         << " |" << std::setw(9) << tage_unique_entries[i]
         << " |" << std::setw(5) << pct(tage_unique_entries[i], MAX_TABLE_ENTRIES) << "%"
         << " |" << std::setw(7) << pct(c.provider_correct[i], c.provider_count[i]) << "%"
         << " |" << std::setw(7) << c.alloc_per_table[i] << "\n";
    }

    // u-bit
    os << "\nU-bit Writes:\n";
    os << "  Table  | U set    | U clear  | Turnover\n";
    os << "  -------+----------+----------+----------\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      u64 total_u = c.u_set_count[i] + c.u_clear_count[i];
      os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << ""
         << "|" << std::setw(9) << c.u_set_count[i]
         << " |" << std::setw(9) << c.u_clear_count[i]
         << " |" << std::setw(7) << pct(c.u_clear_count[i], total_u) << "%\n";
    }
    os << "  Decay fires: " << c.decay_fire_count
       << "  Epoch resets: " << c.epoch_resets
       << "  Avg uctr: " << (c.uctr_samples > 0 ? double(c.uctr_sum) / c.uctr_samples : 0) << "\n";

    // rwram
    os << "\nRWRAM Bank Stats (per table, hyst+u combined):\n";
    os << "  Table  |   Reads  |  Writes  |  Direct  | Buffered | Conflict | BufPend\n";
    os << "  -------+----------+----------+----------+----------+----------+--------\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << ""
         << "|" << std::setw(9) << c.rwram_reads[i]
         << " |" << std::setw(9) << c.rwram_writes[i]
         << " |" << std::setw(9) << c.rwram_writes_direct[i]
         << " |" << std::setw(9) << c.rwram_writes_buffered[i]
         << " |" << std::setw(9) << c.rwram_same_bank_conflict[i]
         << " |" << std::setw(7) << c.rwram_buffer_pending[i] << "\n";
    }

    os << "=== End TageDirect Monitor ===\n\n";
  }

private:
  template <size_t SZ>
  static void print_top_n(std::ostream &os, const std::array<u64, SZ> &hist, u64 n) {
    // Find top-n entries by count
    std::array<std::pair<u64, u64>, 64> top{}; // {count, index}
    u64 found = 0;
    for (u64 i = 0; i < SZ && found < 64; i++) {
      if (hist[i] > 0) top[found++] = {hist[i], i};
    }
    // Simple selection sort for top-n
    for (u64 i = 0; i < std::min(n, found); i++) {
      for (u64 j = i + 1; j < found; j++) {
        if (top[j].first > top[i].first) std::swap(top[i], top[j]);
      }
      os << top[i].second << ":" << top[i].first << " ";
    }
  }
};
