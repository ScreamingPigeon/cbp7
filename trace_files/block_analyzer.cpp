// Block Analyzer for TageAhead banking decision
// Simulates block formation and analyzes predecessor→successor patterns.
//
// Usage: block_analyzer trace.gz [FETCH_WIDTH] [MAX_BANKS] [MAX_INSTRUCTIONS]
//
// Reports:
//   1. Block statistics (count, branches per block, taken/not-taken per block)
//   2. Entry/exit offset distributions
//   3. Successor count per predecessor
//   4. Successor frequency coverage (top-K coverage for K=1..MAX_BANKS)
//   5. Path distribution (fall-through vs taken-1 vs taken-2 etc.)
//   6. History staleness (predecessor changes for same block)
//   7. Predecessor→successor repeat rate
//   8. Collision rate per BANKS value
//   9. Hot/cold branch analysis
//  10. Misprediction model (per-PC bimodal prediction accuracy)

#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include "trace_reader.hpp"

struct BlockInfo {
  uint64_t start_pc;
  uint64_t start_offset; // lane within fetch line
  uint64_t exit_offset;  // lane of last instruction
  uint64_t successor_pc; // next block's PC
  uint64_t num_cond_branches;
  uint64_t num_taken;
  uint64_t num_not_taken;
  bool ends_with_taken; // block ended because of a taken branch
};

struct PredecessorEntry {
  // successor_pc → count
  std::unordered_map<uint64_t, uint64_t> successors;
  uint64_t total = 0;
};

struct BranchInfo {
  uint64_t total = 0;
  uint64_t taken = 0;
  uint64_t mispred_bimodal = 0; // mispredictions if using majority-vote predictor
};

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " trace.gz [FETCH_WIDTH=16] [MAX_BANKS=8] [MAX_INSTRUCTIONS=0]\n";
    return 1;
  }

  uint64_t FETCH_WIDTH = (argc > 2) ? std::stoull(argv[2]) : 16;
  uint64_t MAX_BANKS = (argc > 3) ? std::stoull(argv[3]) : 8;
  uint64_t MAX_INSTR = (argc > 4) ? std::stoull(argv[4]) : 0; // 0 = all

  uint64_t LOG_FETCH_WIDTH = 0;
  { uint64_t tmp = FETCH_WIDTH; while (tmp > 1) { tmp >>= 1; LOG_FETCH_WIDTH++; } }

  trace_reader reader(argv[1], "block_analysis");

  // ---- Tracking state ----
  uint64_t total_instructions = 0;
  uint64_t total_blocks = 0;

  // Block formation state
  uint64_t block_start_pc = 0;
  uint64_t block_start_offset = 0;
  uint64_t block_instr_count = 0;
  uint64_t block_cond_branches = 0;
  uint64_t block_taken = 0;
  uint64_t block_not_taken = 0;
  bool in_block = false;

  // Histograms
  std::map<uint64_t, uint64_t> branches_per_block_hist;
  std::map<uint64_t, uint64_t> taken_per_block_hist;
  std::map<uint64_t, uint64_t> entry_offset_hist;
  std::map<uint64_t, uint64_t> exit_offset_hist;
  std::map<uint64_t, uint64_t> block_size_hist;

  // Predecessor tracking
  // Key: line-level PC of predecessor block (PC >> (LOG_FW + 2))
  std::unordered_map<uint64_t, PredecessorEntry> predecessor_map;
  uint64_t last_block_line_pc = 0;
  bool has_last_block = false;

  // Path distribution: index = number of taken branches before block exit
  std::map<uint64_t, uint64_t> path_dist;

  // History staleness: for each block line PC, track which predecessors lead to it
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, uint64_t>>
      block_predecessors;

  // Repeat rate: consecutive same predecessor→successor
  uint64_t repeat_count = 0;
  uint64_t transition_count = 0;
  uint64_t last_pred_line_pc = 0;
  uint64_t last_succ_line_pc = 0;
  bool has_last_transition = false;

  // Per-branch stats
  std::unordered_map<uint64_t, BranchInfo> branch_stats;

  // Process trace
  instruction prev_inst;
  prev_inst.pc = 0;
  bool has_prev = false;

  auto finish_block = [&](uint64_t successor_pc) {
    total_blocks++;
    uint64_t exit_offset = (block_instr_count - 1 + block_start_offset) % FETCH_WIDTH;

    branches_per_block_hist[block_cond_branches]++;
    taken_per_block_hist[block_taken]++;
    entry_offset_hist[block_start_offset]++;
    exit_offset_hist[exit_offset]++;
    block_size_hist[block_instr_count]++;
    path_dist[block_taken]++;

    uint64_t block_line_pc = block_start_pc >> (LOG_FETCH_WIDTH + 2);
    uint64_t succ_line_pc = successor_pc >> (LOG_FETCH_WIDTH + 2);

    // Predecessor→successor tracking
    if (has_last_block) {
      predecessor_map[last_block_line_pc].successors[succ_line_pc]++;
      predecessor_map[last_block_line_pc].total++;

      // History staleness
      block_predecessors[succ_line_pc][last_block_line_pc]++;

      // Repeat rate
      transition_count++;
      if (has_last_transition && last_pred_line_pc == last_block_line_pc &&
          last_succ_line_pc == succ_line_pc) {
        repeat_count++;
      }
      last_pred_line_pc = last_block_line_pc;
      last_succ_line_pc = succ_line_pc;
      has_last_transition = true;
    }

    last_block_line_pc = block_line_pc;
    has_last_block = true;

    // Reset for next block
    in_block = false;
    block_instr_count = 0;
    block_cond_branches = 0;
    block_taken = 0;
    block_not_taken = 0;
  };

  try {
    while (MAX_INSTR == 0 || total_instructions < MAX_INSTR) {
      instruction inst = reader.next_instruction();
      total_instructions++;

      // Fix discontinuities (same as cbp.hpp)
      if (has_prev) {
        inst.taken_branch = inst.taken_branch || (inst.pc != prev_inst.pc + 4);
      }

      // Start new block if needed
      if (!in_block) {
        block_start_pc = inst.pc;
        block_start_offset = (inst.pc >> 2) & (FETCH_WIDTH - 1);
        block_instr_count = 0;
        in_block = true;
      }

      block_instr_count++;

      // Track conditional branches
      if (inst.inst_class == INST_CLASS::BR_COND) {
        block_cond_branches++;

        auto &bs = branch_stats[inst.pc];
        bs.total++;
        if (inst.taken_branch) {
          bs.taken++;
          block_taken++;
        } else {
          block_not_taken++;
        }
        // Bimodal misprediction model: predict majority direction
        bool predict_taken = (bs.taken * 2 >= bs.total);
        if (predict_taken != inst.taken_branch)
          bs.mispred_bimodal++;
      }

      // Block ends when:
      // 1. Taken branch (any type)
      // 2. Reached line boundary (FETCH_WIDTH instructions from start offset)
      // 3. End of trace (handled by catch)
      bool taken = inst.branch && inst.taken_branch;
      uint64_t current_offset =
          (block_start_offset + block_instr_count - 1) % FETCH_WIDTH;
      bool line_end = (current_offset == FETCH_WIDTH - 1);

      if (taken || line_end) {
        uint64_t successor_pc = inst.next_pc;
        if (!taken && line_end)
          successor_pc = inst.pc + 4; // fall-through
        finish_block(successor_pc);
      }

      prev_inst = inst;
      has_prev = true;
    }
  } catch (const out_of_instructions &) {
    if (in_block && has_prev) {
      finish_block(prev_inst.pc + 4);
    }
  }

  // ---- Report ----
  std::cout << "=== Block Analysis ===" << std::endl;
  std::cout << "FETCH_WIDTH: " << FETCH_WIDTH << std::endl;
  std::cout << "Total instructions: " << total_instructions << std::endl;
  std::cout << "Total blocks: " << total_blocks << std::endl;
  std::cout << "Avg instructions/block: " << std::fixed << std::setprecision(1)
            << double(total_instructions) / total_blocks << std::endl;
  std::cout << std::endl;

  // 1. Branches per block
  std::cout << "--- 1. Conditional branches per block ---" << std::endl;
  for (auto &[n, count] : branches_per_block_hist) {
    std::cout << "  " << n << " branches: " << count << " blocks ("
              << std::setprecision(1) << 100.0 * count / total_blocks << "%)"
              << std::endl;
  }
  std::cout << std::endl;

  // 2. Taken branches per block
  std::cout << "--- 2. Taken branches per block ---" << std::endl;
  for (auto &[n, count] : taken_per_block_hist) {
    std::cout << "  " << n << " taken: " << count << " blocks ("
              << std::setprecision(1) << 100.0 * count / total_blocks << "%)"
              << std::endl;
  }
  std::cout << std::endl;

  // 3. Entry offset distribution
  std::cout << "--- 3. Entry offset distribution ---" << std::endl;
  for (auto &[off, count] : entry_offset_hist) {
    std::cout << "  offset " << std::setw(2) << off << ": " << count << " ("
              << std::setprecision(1) << 100.0 * count / total_blocks << "%)"
              << std::endl;
  }
  std::cout << std::endl;

  // 4. Exit offset distribution
  std::cout << "--- 4. Exit offset distribution ---" << std::endl;
  for (auto &[off, count] : exit_offset_hist) {
    std::cout << "  offset " << std::setw(2) << off << ": " << count << " ("
              << std::setprecision(1) << 100.0 * count / total_blocks << "%)"
              << std::endl;
  }
  std::cout << std::endl;

  // 5. Block size distribution
  std::cout << "--- 5. Block size (instructions) distribution ---" << std::endl;
  for (auto &[sz, count] : block_size_hist) {
    std::cout << "  " << std::setw(3) << sz << " instrs: " << count << " ("
              << std::setprecision(1) << 100.0 * count / total_blocks << "%)"
              << std::endl;
  }
  std::cout << std::endl;

  // 6. Path distribution
  std::cout << "--- 6. Path distribution (taken count at block exit) ---"
            << std::endl;
  for (auto &[path, count] : path_dist) {
    std::cout << "  path " << path << " (";
    if (path == 0)
      std::cout << "fall-through";
    else
      std::cout << "taken after " << path << " branch(es)";
    std::cout << "): " << count << " (" << std::setprecision(1)
              << 100.0 * count / total_blocks << "%)" << std::endl;
  }
  std::cout << std::endl;

  // 7. Successor count per predecessor + coverage per BANKS
  std::cout << "--- 7. Successor count per predecessor ---" << std::endl;
  std::map<uint64_t, uint64_t> succ_count_hist;
  for (auto &[pred, entry] : predecessor_map) {
    succ_count_hist[entry.successors.size()]++;
  }
  for (auto &[n, count] : succ_count_hist) {
    std::cout << "  " << n << " successors: " << count << " predecessors ("
              << std::setprecision(1)
              << 100.0 * count / predecessor_map.size() << "%)" << std::endl;
  }
  std::cout << std::endl;

  // 8. Coverage and collision rate per BANKS value
  std::cout << "--- 8. Coverage by top-K successors (BANKS analysis) ---"
            << std::endl;
  for (uint64_t K = 1; K <= MAX_BANKS; K++) {
    uint64_t covered_transitions = 0;
    uint64_t total_transitions = 0;
    for (auto &[pred, entry] : predecessor_map) {
      // Sort successors by frequency
      std::vector<uint64_t> freq;
      for (auto &[succ, count] : entry.successors)
        freq.push_back(count);
      std::sort(freq.rbegin(), freq.rend());
      // Top-K coverage
      uint64_t covered = 0;
      for (uint64_t i = 0; i < std::min(K, (uint64_t)freq.size()); i++)
        covered += freq[i];
      covered_transitions += covered;
      total_transitions += entry.total;
    }
    double coverage = 100.0 * covered_transitions / total_transitions;
    double collision = 100.0 - coverage;
    std::cout << "  BANKS=" << K << ": coverage=" << std::setprecision(2)
              << coverage << "%, collision=" << std::setprecision(2)
              << collision << "%" << std::endl;
  }
  std::cout << std::endl;

  // 9. History staleness
  std::cout << "--- 9. History staleness (predecessors per block) ---"
            << std::endl;
  std::map<uint64_t, uint64_t> pred_count_hist;
  for (auto &[block, preds] : block_predecessors) {
    pred_count_hist[preds.size()]++;
  }
  for (auto &[n, count] : pred_count_hist) {
    std::cout << "  " << n << " predecessors: " << count << " blocks ("
              << std::setprecision(1)
              << 100.0 * count / block_predecessors.size() << "%)" << std::endl;
  }
  std::cout << std::endl;

  // 10. Predecessor→successor repeat rate
  std::cout << "--- 10. Repeat rate ---" << std::endl;
  std::cout << "  Total transitions: " << transition_count << std::endl;
  std::cout << "  Same as previous:  " << repeat_count << " ("
            << std::setprecision(1) << 100.0 * repeat_count / transition_count
            << "%)" << std::endl;
  std::cout << std::endl;

  // 11. Hot/cold branches (top 20)
  std::cout << "--- 11. Hot branches (top 20 by count) ---" << std::endl;
  std::vector<std::pair<uint64_t, BranchInfo>> sorted_branches(
      branch_stats.begin(), branch_stats.end());
  std::sort(sorted_branches.begin(), sorted_branches.end(),
            [](auto &a, auto &b) { return a.second.total > b.second.total; });
  uint64_t total_cond_branches = 0;
  uint64_t total_bimodal_mispred = 0;
  for (auto &[pc, info] : branch_stats) {
    total_cond_branches += info.total;
    total_bimodal_mispred += info.mispred_bimodal;
  }
  std::cout << "  Total conditional branches: " << total_cond_branches
            << std::endl;
  std::cout << "  Bimodal misprediction rate: " << std::setprecision(2)
            << 100.0 * total_bimodal_mispred / total_cond_branches << "%"
            << std::endl;
  std::cout << std::setw(18) << "PC" << std::setw(10) << "count"
            << std::setw(8) << "taken%" << std::setw(10) << "bim_misp%"
            << std::endl;
  for (uint64_t i = 0; i < std::min(uint64_t(20), (uint64_t)sorted_branches.size());
       i++) {
    auto &[pc, info] = sorted_branches[i];
    std::cout << "  0x" << std::hex << std::setw(14) << pc << std::dec
              << std::setw(10) << info.total << std::setw(7)
              << std::setprecision(1) << 100.0 * info.taken / info.total << "%"
              << std::setw(9) << std::setprecision(1)
              << 100.0 * info.mispred_bimodal / info.total << "%" << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
