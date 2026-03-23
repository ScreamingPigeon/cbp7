#pragma once
// TAGE Performance Monitor — software-only instrumentation
// Gated by -DTAGE_MONITOR -DCHEATING_MODE at compile time.
// All counters are plain C++ types (zero HARCOM cost).

#include <cstdint>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>

using u64 = uint64_t;

template <u64 NUM_TABLES, u64 MAX_TABLE_SIZE>
struct TageMonitor {

  static constexpr u64 WINDOW_SIZE = 100000; // branches per window

  // ======== Counter block (used for both cumulative and per-window) ========
  struct Counters {
    u64 branches = 0;
    u64 mispredictions = 0;
    std::array<u64, NUM_TABLES + 1> provider_count{};
    std::array<u64, NUM_TABLES + 1> provider_correct{};
    std::array<u64, NUM_TABLES + 1> alt_count{};
    u64 meta_override_count = 0;
    u64 meta_override_correct = 0;
    u64 alloc_attempts = 0;
    u64 alloc_success = 0;
    u64 alloc_fail = 0;
    std::array<u64, NUM_TABLES> alloc_per_table{};
    u64 p1_p2_disagree = 0;
    u64 p1_correct_p2_wrong = 0;
    u64 p2_correct_p1_wrong = 0;
    std::array<u64, NUM_TABLES> tag_lookups{};
    std::array<u64, NUM_TABLES> tag_matches{};
    std::array<u64, NUM_TABLES> u_set_count{};
    std::array<u64, NUM_TABLES> u_clear_count{};
    u64 decay_fire_count = 0;
    u64 epoch_reset_count = 0;

    void reset() { *this = Counters{}; }
  };

  Counters cum;  // cumulative (lifetime)
  Counters win;  // current window
  u64 window_num = 0;
  bool header_printed = false;

  // ======== Shadow memory (per-table, per-entry) ========
  // Tracks current state of each entry for occupancy/usefulness analysis.
  struct EntryState {
    bool occupied = false;  // has this entry ever been written (tag allocated)
    uint8_t u_value = 0;    // current u-bit value
  };
  std::array<std::array<EntryState, MAX_TABLE_SIZE>, NUM_TABLES> shadow_table{};
  std::array<u64, NUM_TABLES> table_size{};  // actual size per table (set at init)

  void set_table_sizes(const std::array<u64, NUM_TABLES> &sizes) {
    table_size = sizes;
  }

  // Called when a tag is written (allocation) to a table entry
  void record_tag_write(u64 table, u64 index) {
    if (index < MAX_TABLE_SIZE) {
      shadow_table[table][index].occupied = true;
    }
  }

  // Called when a u-bit is written
  void record_u_write(u64 table, u64 index, bool new_u) {
    if (index < MAX_TABLE_SIZE) {
      shadow_table[table][index].u_value = new_u ? 1 : 0;
    }
    if (new_u) { cum.u_set_count[table]++; win.u_set_count[table]++; }
    else { cum.u_clear_count[table]++; win.u_clear_count[table]++; }
  }

  // Called on epoch reset — clear all u-bits in shadow
  void record_epoch_reset() {
    for (u64 t = 0; t < NUM_TABLES; t++) {
      u64 sz = table_size[t] > 0 ? table_size[t] : MAX_TABLE_SIZE;
      for (u64 i = 0; i < sz; i++) {
        shadow_table[t][i].u_value = 0;
      }
    }
  }

  // Compute occupancy and avg usefulness from shadow state
  struct TableSnapshot {
    u64 occupied = 0;      // entries with tag written
    u64 total = 0;         // table size
    u64 u_nonzero = 0;     // entries with u > 0
    double avg_u = 0.0;    // avg u-value of occupied entries
  };

  TableSnapshot snapshot_table(u64 table) const {
    TableSnapshot s;
    s.total = table_size[table] > 0 ? table_size[table] : MAX_TABLE_SIZE;
    u64 u_sum = 0;
    for (u64 i = 0; i < s.total; i++) {
      auto &e = shadow_table[table][i];
      if (e.occupied) {
        s.occupied++;
        u_sum += e.u_value;
        if (e.u_value > 0) s.u_nonzero++;
      }
    }
    s.avg_u = s.occupied > 0 ? static_cast<double>(u_sum) / s.occupied : 0.0;
    return s;
  }

  // ======== Shadow prediction state (set in predict2, consumed in update_cycle) ========
  static constexpr u64 MAX_FW = 64;
  std::array<u64, MAX_FW> shadow_provider{};
  std::array<u64, MAX_FW> shadow_alt{};
  std::array<bool, MAX_FW> shadow_meta_overrode{};
  std::array<bool, MAX_FW> shadow_p1_pred{};
  std::array<bool, MAX_FW> shadow_p2_pred{};

  // ======== Helpers ========

  static u64 decode_provider(u64 one_hot) {
    if (one_hot == 0) return NUM_TABLES;
    for (u64 i = 0; i <= NUM_TABLES; i++) {
      if (one_hot & (u64(1) << i)) return i;
    }
    return 0;
  }

  void record_prediction(u64 offset, u64 match1_val, u64 match2_val,
                          bool meta_overrode, bool p1_pred, bool p2_pred) {
    shadow_provider[offset] = decode_provider(match1_val);
    shadow_alt[offset] = decode_provider(match2_val);
    shadow_meta_overrode[offset] = meta_overrode;
    shadow_p1_pred[offset] = p1_pred;
    shadow_p2_pred[offset] = p2_pred;
  }

  void record_tag_matches(u64 match_val) {
    for (u64 i = 0; i < NUM_TABLES; i++) {
      cum.tag_lookups[i]++;
      win.tag_lookups[i]++;
      if (match_val & (u64(1) << (i + 1))) {
        cum.tag_matches[i]++;
        win.tag_matches[i]++;
      }
    }
  }

  void record_outcome(u64 offset, bool actual_taken, bool mispredict) {
    auto record = [&](Counters &c) {
      c.branches++;
      if (mispredict) c.mispredictions++;
      u64 prov = shadow_provider[offset];
      if (prov <= NUM_TABLES) {
        c.provider_count[prov]++;
        if (!mispredict) c.provider_correct[prov]++;
      }
      u64 alt = shadow_alt[offset];
      if (alt <= NUM_TABLES) c.alt_count[alt]++;
      if (shadow_meta_overrode[offset]) {
        c.meta_override_count++;
        if (!mispredict) c.meta_override_correct++;
      }
      if (shadow_p1_pred[offset] != shadow_p2_pred[offset]) {
        c.p1_p2_disagree++;
        bool p1_correct = (shadow_p1_pred[offset] == actual_taken);
        bool p2_correct = (shadow_p2_pred[offset] == actual_taken);
        if (p1_correct && !p2_correct) c.p1_correct_p2_wrong++;
        if (p2_correct && !p1_correct) c.p2_correct_p1_wrong++;
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

  // ======== Window output (compact CSV per window) ========

  void print_window(std::ostream &os) {
    if (!header_printed) {
      os << "# window,branches,misp%,";
      os << "prov_bim%,";
      for (u64 i = 0; i < NUM_TABLES; i++) os << "prov_T" << i << "%,";
      os << "alloc_ok%,meta_ovr,p1p2_dis%,";
      // table snapshot columns
      for (u64 i = 0; i < NUM_TABLES; i++) os << "T" << i << "_occ%,";
      for (u64 i = 0; i < NUM_TABLES; i++) os << "T" << i << "_useful%,";
      for (u64 i = 0; i < NUM_TABLES; i++) os << "T" << i << "_avg_u,";
      os << "\n";
      header_printed = true;
    }
    auto pct = [](u64 num, u64 den) -> double {
      return den > 0 ? 100.0 * num / den : 0.0;
    };
    os << std::fixed << std::setprecision(1);
    os << window_num << "," << win.branches << ","
       << pct(win.mispredictions, win.branches) << ",";
    os << pct(win.provider_count[0], win.branches) << ",";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << pct(win.provider_count[i + 1], win.branches) << ",";
    os << pct(win.alloc_success, win.alloc_attempts) << ","
       << win.meta_override_count << ","
       << pct(win.p1_p2_disagree, win.branches) << ",";
    // table snapshots
    for (u64 i = 0; i < NUM_TABLES; i++) {
      auto s = snapshot_table(i);
      os << pct(s.occupied, s.total) << ",";
    }
    for (u64 i = 0; i < NUM_TABLES; i++) {
      auto s = snapshot_table(i);
      os << pct(s.u_nonzero, s.occupied) << ",";
    }
    for (u64 i = 0; i < NUM_TABLES; i++) {
      auto s = snapshot_table(i);
      os << std::setprecision(2) << s.avg_u << ",";
    }
    os << "\n";
  }

  // ======== End-of-trace summary ========

  void print_summary(std::ostream &os = std::cerr) const {
    auto pct = [](u64 num, u64 den) -> double {
      return den > 0 ? 100.0 * num / den : 0.0;
    };
    const auto &c = cum;

    os << "\n=== TAGE Monitor Summary ===\n";
    os << "Branches: " << c.branches
       << "  Mispredictions: " << c.mispredictions;
    if (c.branches > 0)
      os << " (" << std::fixed << std::setprecision(2)
         << pct(c.mispredictions, c.branches) << "%)";
    os << "\n\n";

    os << std::fixed << std::setprecision(1);
    os << "Provider Distribution:\n";
    os << "  Table     | Provided  |  Correct  | Accuracy | TagMatch%\n";
    os << "  ----------+-----------+-----------+----------+----------\n";
    os << "  Bimodal   |" << std::setw(10) << c.provider_count[0]
       << " |" << std::setw(10) << c.provider_correct[0]
       << " |" << std::setw(7) << pct(c.provider_correct[0], c.provider_count[0]) << "%"
       << " |      --\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      os << "  T" << i << std::setw(8 - (i >= 10 ? 2 : 1)) << ""
         << "|" << std::setw(10) << c.provider_count[i + 1]
         << " |" << std::setw(10) << c.provider_correct[i + 1]
         << " |" << std::setw(7) << pct(c.provider_correct[i + 1], c.provider_count[i + 1]) << "%"
         << " |" << std::setw(7) << pct(c.tag_matches[i], c.tag_lookups[i]) << "%\n";
    }

    os << "\nMeta Override: " << c.meta_override_count << " times";
    if (c.meta_override_count > 0)
      os << ", " << c.meta_override_correct << " correct ("
         << pct(c.meta_override_correct, c.meta_override_count) << "%)";
    os << "\n";

    os << "Allocation: " << c.alloc_attempts << " attempts, "
       << c.alloc_success << " success (" << pct(c.alloc_success, c.alloc_attempts) << "%), "
       << c.alloc_fail << " fail\n";
    os << "  Alloc/table:";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << " T" << i << "=" << c.alloc_per_table[i];
    os << "\n";

    os << "P1/P2 disagree: " << c.p1_p2_disagree
       << " (P1 right: " << c.p1_correct_p2_wrong
       << ", P2 right: " << c.p2_correct_p1_wrong << ")\n";

    // Table occupancy and usefulness
    os << "\nTable State (shadow memory):\n";
    os << "  Table  | Size  | Occupied |  Occ%  | U>0  | Useful% | Avg U\n";
    os << "  -------+-------+----------+--------+------+---------+------\n";
    for (u64 i = 0; i < NUM_TABLES; i++) {
      auto s = const_cast<TageMonitor*>(this)->snapshot_table(i);
      os << "  T" << i << std::setw(5 - (i >= 10 ? 1 : 0)) << ""
         << "|" << std::setw(6) << s.total
         << " |" << std::setw(9) << s.occupied
         << " |" << std::setw(6) << pct(s.occupied, s.total) << "%"
         << " |" << std::setw(5) << s.u_nonzero
         << " |" << std::setw(7) << pct(s.u_nonzero, s.occupied) << "%"
         << " |" << std::setw(5) << std::setprecision(2) << s.avg_u << "\n";
    }

    // U-bit write stats
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
       << "  Epoch resets: " << c.epoch_reset_count << "\n";
    os << "=== End TAGE Monitor ===\n\n";
  }
};
