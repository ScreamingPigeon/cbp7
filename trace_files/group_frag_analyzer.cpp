// Group Fragmentation Analyzer
//
// Models the GROUP_BITS fragmentation cost in TageAheadHC_IR:
// each TAGE entry is tagged with a branch's ordinal position (group_id)
// among conditional branches in the fetch block. A branch that shifts
// position across different dynamic occurrences fragments its training
// signal across multiple group_id namespaces, requiring more capacity.
//
// Usage: group_frag_analyzer trace.gz [LINEINST=256] [N=7] [MAX_INSTR=0]
//
// Reports:
//   1. Group_id stability per branch PC
//   2. Capacity overhead distribution (how much extra capacity fragmentation needs)
//   3. Weighted aggregate capacity multiplier
//   4. Hot fragmented branches

#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <unordered_map>
#include <vector>
#include "trace_reader.hpp"

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " trace.gz [LINEINST=256] [N=7] [MAX_INSTR=0]\n";
    return 1;
  }

  uint64_t LINEINST  = (argc > 2) ? std::stoull(argv[2]) : 256;
  uint64_t N         = (argc > 3) ? std::stoull(argv[3]) : 7;
  uint64_t MAX_INSTR = (argc > 4) ? std::stoull(argv[4]) : 0;

  uint64_t LOG_LINEINST = 0;
  { uint64_t tmp = LINEINST; while (tmp > 1) { tmp >>= 1; LOG_LINEINST++; } }

  trace_reader reader(argv[1], "group_frag");

  // per branch PC: count of appearances at each group_id (0..N-1)
  std::unordered_map<uint64_t, std::vector<uint64_t>> pc_group_counts;

  // block formation state
  uint64_t total_instructions = 0;
  uint64_t total_blocks = 0;
  uint64_t block_start_pc = 0;
  uint64_t block_start_offset = 0;
  uint64_t block_cond_branches = 0;
  bool in_block = false;

  // current block's branch PCs in order (for assigning group_ids on close)
  std::vector<uint64_t> block_branch_pcs;

  instruction prev_inst;
  prev_inst.pc = 0;
  bool has_prev = false;

  auto finish_block = [&](bool /*taken_exit*/) {
    total_blocks++;
    // assign group_ids to the branches we saw
    for (uint64_t g = 0; g < block_branch_pcs.size(); g++) {
      uint64_t pc = block_branch_pcs[g];
      auto &counts = pc_group_counts[pc];
      if (counts.size() <= g) counts.resize(g + 1, 0);
      counts[g]++;
    }
    in_block = false;
    block_cond_branches = 0;
    block_branch_pcs.clear();
  };

  try {
    while (MAX_INSTR == 0 || total_instructions < MAX_INSTR) {
      instruction inst = reader.next_instruction();
      total_instructions++;

      if (has_prev)
        inst.taken_branch = inst.taken_branch || (inst.pc != prev_inst.pc + 4);

      if (!in_block) {
        block_start_pc = inst.pc;
        block_start_offset = (inst.pc >> 2) & (LINEINST - 1);
        in_block = true;
      }

      uint64_t current_offset = block_start_offset + (inst.pc - block_start_pc) / 4;

      if (inst.inst_class == INST_CLASS::BR_COND) {
        block_branch_pcs.push_back(inst.pc);
        block_cond_branches++;
      }

      bool taken = inst.branch && inst.taken_branch;
      bool line_end = (current_offset >= LINEINST - 1);
      bool branch_limit = (block_cond_branches >= N && inst.inst_class == INST_CLASS::BR_COND);

      if (taken || line_end || branch_limit)
        finish_block(taken);

      prev_inst = inst;
      has_prev = true;
    }
  } catch (const out_of_instructions &) {
    if (in_block && has_prev) finish_block(false);
  }

  // ── Analysis ─────────────────────────────────────────────────────────────

  // Per branch: total occurrences, dominant group count, unique group_ids
  struct BranchFrag {
    uint64_t pc;
    uint64_t total;
    uint64_t dominant;   // count at most-common group_id
    uint64_t unique_gids;
    double overhead;     // total / dominant — capacity multiplier
  };

  std::vector<BranchFrag> branches;
  uint64_t grand_total = 0;
  uint64_t grand_dominant = 0;

  for (auto &[pc, counts] : pc_group_counts) {
    uint64_t total = 0;
    uint64_t dominant = 0;
    uint64_t unique_gids = 0;
    for (uint64_t c : counts) {
      total += c;
      if (c > dominant) dominant = c;
      if (c > 0) unique_gids++;
    }
    double overhead = (dominant > 0) ? double(total) / dominant : 1.0;
    branches.push_back({pc, total, dominant, unique_gids, overhead});
    grand_total    += total;
    grand_dominant += dominant;
  }

  // sort by total occurrences desc
  std::sort(branches.begin(), branches.end(),
            [](auto &a, auto &b) { return a.total > b.total; });

  double aggregate_overhead = double(grand_total) / double(grand_dominant);

  std::cout << "=== Group Fragmentation Analysis ===" << std::endl;
  std::cout << "LINEINST=" << LINEINST << "  N=" << N << std::endl;
  std::cout << "Total instructions: " << total_instructions << std::endl;
  std::cout << "Total blocks: " << total_blocks << std::endl;
  std::cout << "Unique branch PCs: " << branches.size() << std::endl;
  std::cout << std::endl;

  // 1. Aggregate capacity overhead
  std::cout << "--- 1. Aggregate capacity overhead ---" << std::endl;
  std::cout << "  Total branch occurrences:    " << grand_total << std::endl;
  std::cout << "  Dominant-group occurrences:  " << grand_dominant << std::endl;
  std::cout << "  Non-dominant (fragmented):   " << (grand_total - grand_dominant)
            << " (" << std::fixed << std::setprecision(1)
            << 100.0 * (grand_total - grand_dominant) / grand_total << "%)" << std::endl;
  std::cout << "  Capacity multiplier:         " << std::setprecision(3)
            << aggregate_overhead << "x" << std::endl;
  std::cout << "    → HC_IR needs ~" << std::setprecision(1)
            << (aggregate_overhead - 1.0) * 100.0
            << "% more entries than a fixed-position scheme" << std::endl;
  std::cout << std::endl;

  // 2. Distribution by unique group_ids
  std::cout << "--- 2. Branches by number of unique group_ids ---" << std::endl;
  std::map<uint64_t, std::pair<uint64_t,uint64_t>> gid_hist; // unique_gids → (branch_count, occurrence_count)
  for (auto &b : branches) {
    gid_hist[b.unique_gids].first++;
    gid_hist[b.unique_gids].second += b.total;
  }
  std::cout << std::setw(12) << "unique_gids"
            << std::setw(12) << "branches"
            << std::setw(10) << "br%"
            << std::setw(14) << "occurrences"
            << std::setw(10) << "occ%" << std::endl;
  for (auto &[gids, p] : gid_hist) {
    std::cout << std::setw(12) << gids
              << std::setw(12) << p.first
              << std::setw(9)  << std::setprecision(1) << 100.0*p.first/branches.size() << "%"
              << std::setw(14) << p.second
              << std::setw(9)  << std::setprecision(1) << 100.0*p.second/grand_total << "%"
              << std::endl;
  }
  std::cout << std::endl;

  // 3. Overhead distribution buckets
  std::cout << "--- 3. Capacity overhead distribution (weighted by occurrences) ---" << std::endl;
  struct Bucket { double lo, hi; std::string label; uint64_t br_count = 0; uint64_t occ = 0; };
  std::vector<Bucket> buckets = {
    {1.0, 1.001, "1.00x (no fragmentation)"},
    {1.001, 1.05, "1.00–1.05x"},
    {1.05, 1.10, "1.05–1.10x"},
    {1.10, 1.25, "1.10–1.25x"},
    {1.25, 1.50, "1.25–1.50x"},
    {1.50, 2.00, "1.50–2.00x"},
    {2.00, 1e18, "2.00x+"},
  };
  for (auto &b : branches) {
    for (auto &bkt : buckets) {
      if (b.overhead >= bkt.lo && b.overhead < bkt.hi) {
        bkt.br_count++;
        bkt.occ += b.total;
        break;
      }
    }
  }
  std::cout << std::setw(28) << "overhead"
            << std::setw(10) << "branches"
            << std::setw(10) << "occ%" << std::endl;
  for (auto &bkt : buckets) {
    std::cout << std::setw(28) << bkt.label
              << std::setw(10) << bkt.br_count
              << std::setw(9)  << std::setprecision(1) << 100.0*bkt.occ/grand_total << "%"
              << std::endl;
  }
  std::cout << std::endl;

  // 4. Top 20 most fragmented hot branches (sorted by overhead * total)
  std::cout << "--- 4. Most impactful fragmented branches (top 20 by wasted occurrences) ---" << std::endl;
  std::vector<BranchFrag *> fragmented;
  for (auto &b : branches)
    if (b.unique_gids > 1) fragmented.push_back(&b);
  std::sort(fragmented.begin(), fragmented.end(),
            [](auto *a, auto *b) {
              return (a->total - a->dominant) > (b->total - b->dominant);
            });
  std::cout << std::setw(18) << "PC"
            << std::setw(12) << "total"
            << std::setw(12) << "dominant"
            << std::setw(10) << "wasted"
            << std::setw(10) << "overhead"
            << std::setw(10) << "gids" << std::endl;
  for (uint64_t i = 0; i < std::min(uint64_t(20), (uint64_t)fragmented.size()); i++) {
    auto *b = fragmented[i];
    std::cout << "  0x" << std::hex << std::setw(14) << b->pc << std::dec
              << std::setw(12) << b->total
              << std::setw(12) << b->dominant
              << std::setw(10) << (b->total - b->dominant)
              << std::setw(9)  << std::setprecision(2) << b->overhead << "x"
              << std::setw(10) << b->unique_gids << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
