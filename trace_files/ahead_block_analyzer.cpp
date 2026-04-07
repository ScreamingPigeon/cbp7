// Block Analyzer for Ahead TAGE Banking Decision
//
// Terminology:
//   LINEINST  = max instructions per aligned line (like a cache line boundary)
//   N         = max conditional branches predicted per block
//   BANKS     = number of ahead banks (one per possible exit path)
//   Lane      = branch rank within block (0 = first cond branch, 1 = second, ...)
//   Block     = prediction unit; ends at taken branch, line boundary, or N branches
//   Path      = exit path of a block: 0 = fall-through, k = taken after kth branch
//
// Usage: block_analyzer trace.gz [LINEINST=64] [N=7] [MAX_BANKS=8] [MAX_INSTRUCTIONS=0]
//
// Reports:
//   1. Block statistics (branches per block, block size)
//   2. Path distribution (how blocks exit)
//   3. Successor count per predecessor
//   4. Bank coverage: for K=1..MAX_BANKS, what % of transitions are captured
//   5. Secondary tag analysis: with K banks, how many unique successors remain per bank
//   6. Lane utilization (how many branches per block)
//   7. Hot/cold branch analysis

#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <cmath>
#include "trace_reader.hpp"

struct PredecessorEntry {
  // (path, successor_line_pc) → count
  // path encodes the exit path (0=fall-through, k=taken after kth branch)
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, uint64_t>> path_successors;
  // successor_line_pc → count (aggregated across all paths)
  std::unordered_map<uint64_t, uint64_t> all_successors;
  uint64_t total = 0;
};

struct BranchInfo {
  uint64_t total = 0;
  uint64_t taken = 0;
};

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " trace.gz [LINEINST=64] [N=7] [MAX_BANKS=8] [MAX_INSTRUCTIONS=0]\n";
    return 1;
  }

  uint64_t LINEINST = (argc > 2) ? std::stoull(argv[2]) : 64;
  uint64_t N = (argc > 3) ? std::stoull(argv[3]) : 7;
  uint64_t MAX_BANKS = (argc > 4) ? std::stoull(argv[4]) : 8;
  uint64_t MAX_INSTR = (argc > 5) ? std::stoull(argv[5]) : 0;

  uint64_t LOG_LINEINST = 0;
  { uint64_t tmp = LINEINST; while (tmp > 1) { tmp >>= 1; LOG_LINEINST++; } }

  trace_reader reader(argv[1], "block_analysis");

  // ---- State ----
  uint64_t total_instructions = 0;
  uint64_t total_blocks = 0;

  // Block formation
  uint64_t block_start_pc = 0;
  uint64_t block_start_offset = 0;
  uint64_t block_instr_count = 0;
  uint64_t block_cond_branches = 0;
  uint64_t block_path = 0; // exit path: 0=fall-through, k=taken after kth branch
  bool in_block = false;

  // Histograms
  std::map<uint64_t, uint64_t> branches_per_block_hist;
  std::map<uint64_t, uint64_t> block_size_hist;
  std::map<uint64_t, uint64_t> path_dist;
  std::map<uint64_t, uint64_t> lane_util_hist; // how many lanes used per block

  // Predecessor tracking (keyed by predecessor line PC)
  std::unordered_map<uint64_t, PredecessorEntry> predecessor_map;
  uint64_t last_block_line_pc = 0;
  uint64_t last_block_path = 0;
  bool has_last_block = false;

  // Per-branch stats
  std::unordered_map<uint64_t, BranchInfo> branch_stats;

  instruction prev_inst;
  prev_inst.pc = 0;
  bool has_prev = false;

  auto finish_block = [&](uint64_t successor_pc, bool taken_exit) {
    total_blocks++;

    branches_per_block_hist[block_cond_branches]++;
    block_size_hist[block_instr_count]++;
    lane_util_hist[std::min(block_cond_branches, N)]++;

    block_path = taken_exit ? block_cond_branches : 0;
    path_dist[block_path]++;

    uint64_t block_line_pc = block_start_pc >> (LOG_LINEINST + 2);
    uint64_t succ_line_pc = successor_pc >> (LOG_LINEINST + 2);

    // Predecessor→successor tracking
    if (has_last_block) {
      auto &entry = predecessor_map[last_block_line_pc];
      entry.path_successors[last_block_path][succ_line_pc]++;
      entry.all_successors[succ_line_pc]++;
      entry.total++;
    }

    last_block_line_pc = block_line_pc;
    last_block_path = block_path;
    has_last_block = true;

    in_block = false;
    block_instr_count = 0;
    block_cond_branches = 0;
    block_path = 0;
  };

  try {
    while (MAX_INSTR == 0 || total_instructions < MAX_INSTR) {
      instruction inst = reader.next_instruction();
      total_instructions++;

      if (has_prev) {
        inst.taken_branch = inst.taken_branch || (inst.pc != prev_inst.pc + 4);
      }

      if (!in_block) {
        block_start_pc = inst.pc;
        block_start_offset = (inst.pc >> 2) & (LINEINST - 1);
        block_instr_count = 0;
        in_block = true;
      }

      block_instr_count++;

      if (inst.inst_class == INST_CLASS::BR_COND) {
        block_cond_branches++;
        auto &bs = branch_stats[inst.pc];
        bs.total++;
        if (inst.taken_branch) bs.taken++;
      }

      bool taken = inst.branch && inst.taken_branch;
      uint64_t current_offset = block_start_offset + block_instr_count - 1;
      bool line_end = (current_offset >= LINEINST - 1);
      bool branch_limit = (block_cond_branches >= N && inst.inst_class == INST_CLASS::BR_COND);

      if (taken || line_end || branch_limit) {
        uint64_t successor_pc = taken ? inst.next_pc : (inst.pc + 4);
        finish_block(successor_pc, taken);
      }

      prev_inst = inst;
      has_prev = true;
    }
  } catch (const out_of_instructions &) {
    if (in_block && has_prev)
      finish_block(prev_inst.pc + 4, false);
  }

  // ---- Report ----
  std::cout << "=== Ahead TAGE Block Analysis ===" << std::endl;
  std::cout << "LINEINST=" << LINEINST << "  N=" << N
            << "  MAX_BANKS=" << MAX_BANKS << std::endl;
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

  // 2. Block size
  std::cout << "--- 2. Block size (instructions) ---" << std::endl;
  for (auto &[sz, count] : block_size_hist) {
    if (count > total_blocks / 200) // only show >0.5%
      std::cout << "  " << std::setw(3) << sz << " instrs: " << count << " ("
                << std::setprecision(1) << 100.0 * count / total_blocks << "%)"
                << std::endl;
  }
  std::cout << std::endl;

  // 3. Path distribution
  std::cout << "--- 3. Path distribution (exit path) ---" << std::endl;
  for (auto &[path, count] : path_dist) {
    std::cout << "  path " << path << " (";
    if (path == 0)
      std::cout << "fall-through/line-end";
    else
      std::cout << "taken after branch #" << path;
    std::cout << "): " << count << " (" << std::setprecision(1)
              << 100.0 * count / total_blocks << "%)" << std::endl;
  }
  std::cout << std::endl;

  // 4. Lane utilization
  std::cout << "--- 4. Lane utilization (cond branches per block, capped at N="
            << N << ") ---" << std::endl;
  uint64_t total_lanes_used = 0;
  for (auto &[lanes, count] : lane_util_hist) {
    total_lanes_used += lanes * count;
    std::cout << "  " << lanes << " lanes: " << count << " blocks ("
              << std::setprecision(1) << 100.0 * count / total_blocks << "%)"
              << std::endl;
  }
  std::cout << "  Avg lanes/block: " << std::setprecision(2)
            << double(total_lanes_used) / total_blocks << std::endl;
  std::cout << std::endl;

  // 5. Successor count per predecessor (total, ignoring path)
  std::cout << "--- 5. Unique successors per predecessor ---" << std::endl;
  std::map<uint64_t, uint64_t> succ_count_hist;
  for (auto &[pred, entry] : predecessor_map) {
    succ_count_hist[entry.all_successors.size()]++;
  }
  for (auto &[n, count] : succ_count_hist) {
    std::cout << "  " << n << " successors: " << count << " predecessors ("
              << std::setprecision(1)
              << 100.0 * count / predecessor_map.size() << "%)" << std::endl;
  }
  std::cout << std::endl;

  // 6. Bank coverage analysis (KEY metric)
  // For each K (number of banks), compute:
  //   a) If we use top-K most frequent successors per predecessor, what % of
  //      transitions are covered?
  //   b) If we partition by path (exit path → bank), what % covered?
  std::cout << "--- 6. Bank coverage analysis ---" << std::endl;
  std::cout << std::endl;

  // 6a. Top-K successor coverage (frequency-based, like gshare_ahead banking)
  std::cout << "  6a. Top-K successor coverage (frequency-based):" << std::endl;
  for (uint64_t K = 1; K <= MAX_BANKS; K++) {
    uint64_t covered = 0;
    uint64_t total = 0;
    for (auto &[pred, entry] : predecessor_map) {
      std::vector<uint64_t> freq;
      for (auto &[succ, count] : entry.all_successors)
        freq.push_back(count);
      std::sort(freq.rbegin(), freq.rend());
      uint64_t c = 0;
      for (uint64_t i = 0; i < std::min(K, (uint64_t)freq.size()); i++)
        c += freq[i];
      covered += c;
      total += entry.total;
    }
    std::cout << "    K=" << K << ": " << std::setprecision(2)
              << 100.0 * covered / total << "% covered" << std::endl;
  }
  std::cout << std::endl;

  // 6b. Path-based banking (bank = exit path, like our ahead TAGE design)
  // With B banks, map each path to a bank (path % B). Count coverage.
  std::cout << "  6b. Path-based banking (bank = path % B):" << std::endl;
  for (uint64_t B = 1; B <= MAX_BANKS; B++) {
    uint64_t covered = 0;
    uint64_t total = 0;
    for (auto &[pred, entry] : predecessor_map) {
      // For each bank, find the most frequent successor on paths mapping to that bank
      std::unordered_map<uint64_t, std::unordered_map<uint64_t, uint64_t>> bank_succs;
      for (auto &[path, succs] : entry.path_successors) {
        uint64_t bank = path % B;
        for (auto &[succ, count] : succs)
          bank_succs[bank][succ] += count;
      }
      // For each bank, take the most frequent successor
      for (auto &[bank, succs] : bank_succs) {
        uint64_t max_count = 0;
        for (auto &[succ, count] : succs)
          max_count = std::max(max_count, count);
        covered += max_count;
      }
      total += entry.total;
    }
    std::cout << "    B=" << B << ": " << std::setprecision(2)
              << 100.0 * covered / total << "% covered" << std::endl;
  }
  std::cout << std::endl;

  // 6c. Path-based banking + secondary tag (B banks, secondary tag disambiguates within bank)
  // With B banks and S-bit secondary tag, each bank can distinguish 2^S successors
  std::cout << "  6c. Path-based banking + secondary tag:" << std::endl;
  for (uint64_t B : {1, 2, 4}) {
    for (uint64_t sec_bits = 0; sec_bits <= 5; sec_bits++) {
      uint64_t max_per_bank = uint64_t(1) << sec_bits;
      uint64_t covered = 0;
      uint64_t total = 0;
      for (auto &[pred, entry] : predecessor_map) {
        std::unordered_map<uint64_t, std::unordered_map<uint64_t, uint64_t>> bank_succs;
        for (auto &[path, succs] : entry.path_successors) {
          uint64_t bank = path % B;
          for (auto &[succ, count] : succs)
            bank_succs[bank][succ] += count;
        }
        for (auto &[bank, succs] : bank_succs) {
          std::vector<uint64_t> freq;
          for (auto &[succ, count] : succs)
            freq.push_back(count);
          std::sort(freq.rbegin(), freq.rend());
          uint64_t c = 0;
          for (uint64_t i = 0; i < std::min(max_per_bank, (uint64_t)freq.size()); i++)
            c += freq[i];
          covered += c;
        }
        total += entry.total;
      }
      std::cout << "    B=" << B << " sec_tag=" << sec_bits << "b ("
                << max_per_bank << " entries/bank): "
                << std::setprecision(2) << 100.0 * covered / total
                << "% covered" << std::endl;
    }
  }
  std::cout << std::endl;

  // 7. Hot branches (top 20)
  std::cout << "--- 7. Hot branches (top 20) ---" << std::endl;
  uint64_t total_cond = 0;
  for (auto &[pc, info] : branch_stats) total_cond += info.total;
  std::vector<std::pair<uint64_t, BranchInfo>> sorted_br(
      branch_stats.begin(), branch_stats.end());
  std::sort(sorted_br.begin(), sorted_br.end(),
            [](auto &a, auto &b) { return a.second.total > b.second.total; });
  std::cout << "  Total conditional branches: " << total_cond << std::endl;
  std::cout << std::setw(18) << "PC" << std::setw(10) << "count"
            << std::setw(8) << "taken%" << std::endl;
  for (uint64_t i = 0; i < std::min(uint64_t(20), (uint64_t)sorted_br.size()); i++) {
    auto &[pc, info] = sorted_br[i];
    std::cout << "  0x" << std::hex << std::setw(14) << pc << std::dec
              << std::setw(10) << info.total << std::setw(7)
              << std::setprecision(1) << 100.0 * info.taken / info.total << "%"
              << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
