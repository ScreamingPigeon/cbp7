#pragma once
// TageAhead Performance Monitor — software-only instrumentation
// Gated by -DTAGE_MONITOR at compile time.
// All counters are plain C++ types (zero HARCOM cost).

#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>

using u64 = uint64_t;

template <u64 NUM_TABLES, u64 N, u64 MAX_TABLE_ENTRIES = 2048>
struct TAMonitor {

  static constexpr u64 WINDOW_SIZE = 100000; // branches per window

  static constexpr u64 MAX_BLOCK_INSTR = 64;
  static constexpr u64 MAX_BLOCK_BR = N + 1;

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

    // True block rate (pipeline reuse)
    u64 true_block_count = 0;

    // Train valid skip (first-cycle guard)
    u64 train_skip_count = 0;

    // Provider distribution: 0..NUM_TABLES-1 = TAGE, NUM_TABLES = bimodal
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
    std::array<u64, NUM_TABLES> tag_matches{};    // primary tag match
    std::array<u64, NUM_TABLES> sec_matches{};    // secondary tag match
    std::array<u64, NUM_TABLES> full_matches{};   // both matched

    // Allocation
    u64 alloc_attempts = 0;
    u64 alloc_success = 0;
    u64 alloc_fail = 0; // noalloc (mispredict but no u=0 slot)
    std::array<u64, NUM_TABLES> alloc_per_table{};

    // u-bit
    std::array<u64, NUM_TABLES> u_set_count{};
    std::array<u64, NUM_TABLES> u_clear_count{};
    u64 decay_fire_count = 0;
    u64 epoch_reset_count = 0;

    // Pressure counters (sampled each update_cycle with branches)
    u64 acc_ctr_sum = 0;
    u64 alloc_ctr_sum = 0;
    u64 pressure_samples = 0;

    void reset() { *this = Counters{}; }
  };

  Counters cum;
  Counters win;
  u64 window_num = 0;
  bool header_printed = false;

  // Shadow state (set during resolution, consumed during training)
  // Since TageAhead does resolution in update_cycle, we set + consume in same call
  std::array<u64, N> shadow_provider{};
  std::array<u64, N> shadow_alt{};
  std::array<bool, N> shadow_meta_overrode{};
  std::array<bool, N> shadow_meta_chose_alt{};
  std::array<bool, N> shadow_pred{};

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

  void record_block(u64 block_size, u64 num_branch, bool extra_cycle_fired) {
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

  // u-bit write
  void record_u_write(u64 table, bool new_u) {
    if (new_u) { cum.u_set_count[table]++; win.u_set_count[table]++; }
    else { cum.u_clear_count[table]++; win.u_clear_count[table]++; }
  }

  void record_decay_fire() { cum.decay_fire_count++; win.decay_fire_count++; }
  void record_epoch_reset() {
    cum.epoch_reset_count++;
    win.epoch_reset_count++;
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
    if (one_hot == 0) return NUM_TABLES; // bimodal fallback
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
      os << "# win,br,misp%,extra%,i/blk,br/blk,true_blk%,";
      os << "bim%,";
      for (u64 i = 0; i < NUM_TABLES; i++) os << "T" << i << "%,";
      os << "alloc_ok%,acc_avg,alloc_avg";
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
               : 0);
    os << "\n";
  }

  // ======== End-of-trace summary ========

  void print_summary(std::ostream &os = std::cerr) const {
    const auto &c = cum;

    os << "\n=== TageAhead Monitor Summary ===\n";
    os << "Branches: " << c.branches
       << "  Mispredictions: " << c.mispredictions << " ("
       << std::fixed << std::setprecision(2)
       << pct(c.mispredictions, c.branches) << "%)\n";
    os << "Blocks: " << c.blocks
       << "  Extra cycles: " << c.extra_cycles
       << " (" << pct(c.extra_cycles, c.blocks) << "%)\n";
    os << "True blocks: " << c.true_block_count
       << " (" << pct(c.true_block_count, c.blocks) << "%)\n";
    os << "Train skips (first cycle): " << c.train_skip_count << "\n";

    // Pipeline calls
    os << "\nPipeline Calls:\n";
    os << "  predict1: " << c.predict1_calls
       << "  reuse: " << c.reuse_predict1_calls
       << "  (reuse: "
       << pct(c.reuse_predict1_calls,
              c.predict1_calls + c.reuse_predict1_calls)
       << "%)\n";
    os << "  predict2: " << c.predict2_calls
       << "  reuse: " << c.reuse_predict2_calls << "\n";

    // Block structure
    os << "\nBlock Structure:\n";
    os << "  Avg instr/block: " << std::setprecision(1)
       << (c.blocks > 0 ? double(c.total_block_instr) / c.blocks : 0) << "\n";
    os << "  Avg branches/block: "
       << (c.blocks > 0 ? double(c.total_block_branches) / c.blocks : 0)
       << "\n";
    os << "  Branches/block histogram: ";
    for (u64 i = 0; i < MAX_BLOCK_BR; i++) {
      if (c.block_br_hist[i] > 0)
        os << i << ":" << c.block_br_hist[i] << " ";
    }
    os << "\n";

    // Provider distribution
    os << "\nProvider Distribution:\n";
    os << "  Table     | Provided  |  Correct  | Accuracy | TagMatch% "
          "| SecMatch% | FullMatch%\n";
    os << "  ----------+-----------+-----------+----------+-----------"
          "+-----------+-----------\n";
    os << "  Bimodal   |" << std::setw(10) << c.provider_count[NUM_TABLES]
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

    // Allocation
    os << "\nAllocation:\n";
    os << "  Attempts: " << c.alloc_attempts
       << "  Success: " << c.alloc_success << " ("
       << pct(c.alloc_success, c.alloc_attempts) << "%)"
       << "  Fail: " << c.alloc_fail << "\n";
    os << "  Per table:";
    for (u64 i = 0; i < NUM_TABLES; i++)
      os << " T" << i << "=" << c.alloc_per_table[i];
    os << "\n";

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

    // Pressure
    if (c.pressure_samples > 0) {
      os << "\nPressure Counters:\n";
      os << "  Avg acc_ctr: "
         << double(c.acc_ctr_sum) / c.pressure_samples << "\n";
      os << "  Avg alloc_ctr: "
         << double(c.alloc_ctr_sum) / c.pressure_samples << "\n";
    }

    os << "\n";
  }
};
