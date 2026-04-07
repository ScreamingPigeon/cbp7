// Block Analyzer — sweep LINEINST and N from a single trace read.
//
// Usage: block_analyzer trace.gz [MAX_INSTRUCTIONS=0]
//
// Reads the trace once, records per-instruction events, then sweeps
// LINEINST={8,16,32,64,128,256,512,1024} x N={4,7,16,64} and reports
// total blocks, avg instructions/block, avg branches/block for each config.

#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdint>
#include "trace_reader.hpp"

struct Event {
  uint64_t pc;
  bool is_cond_branch;
  bool is_taken_branch; // any branch type, taken
  bool is_branch;
};

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " trace.gz [MAX_INSTRUCTIONS=0]\n";
    return 1;
  }

  uint64_t MAX_INSTR = (argc > 2) ? std::stoull(argv[2]) : 0;

  trace_reader reader(argv[1], "block_analysis");

  // ---- Pass 1: read trace into event log ----
  std::vector<Event> events;
  events.reserve(MAX_INSTR > 0 ? MAX_INSTR : 40000000);

  instruction prev_inst;
  prev_inst.pc = 0;
  bool has_prev = false;
  uint64_t total_instructions = 0;

  try {
    while (MAX_INSTR == 0 || total_instructions < MAX_INSTR) {
      instruction inst = reader.next_instruction();
      total_instructions++;

      if (has_prev) {
        inst.taken_branch = inst.taken_branch || (inst.pc != prev_inst.pc + 4);
      }

      events.push_back({
        inst.pc,
        inst.inst_class == INST_CLASS::BR_COND,
        inst.branch && inst.taken_branch,
        inst.branch
      });

      prev_inst = inst;
      has_prev = true;
    }
  } catch (const out_of_instructions &) {}

  std::cerr << "Read " << events.size() << " instructions" << std::endl;

  // ---- Pass 2: sweep LINEINST x N ----
  uint64_t lineinst_vals[] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
  uint64_t n_vals[] = {4, 7, 16, 64};

  std::cout << std::setw(8) << "LINEINST"
            << std::setw(5) << "N"
            << std::setw(12) << "TotalBlks"
            << std::setw(10) << "Inst/Blk"
            << std::setw(9) << "Br/Blk"
            << std::endl;
  std::cout << std::string(44, '-') << std::endl;

  for (uint64_t LINEINST : lineinst_vals) {
    for (uint64_t N : n_vals) {
      uint64_t total_blocks = 0;
      uint64_t total_cond = 0;

      // Block formation state
      uint64_t block_start_offset = 0;
      uint64_t block_instr_count = 0;
      uint64_t block_cond_branches = 0;
      bool in_block = false;

      for (auto &ev : events) {
        if (!in_block) {
          block_start_offset = (ev.pc >> 2) & (LINEINST - 1);
          block_instr_count = 0;
          block_cond_branches = 0;
          in_block = true;
        }

        block_instr_count++;
        if (ev.is_cond_branch) {
          block_cond_branches++;
          total_cond++;
        }

        bool taken = ev.is_taken_branch;
        uint64_t current_offset = block_start_offset + block_instr_count - 1;
        bool line_end = (current_offset >= LINEINST - 1);
        bool branch_limit = (block_cond_branches >= N && ev.is_cond_branch);

        if (taken || line_end || branch_limit) {
          total_blocks++;
          in_block = false;
        }
      }
      if (in_block) total_blocks++;

      double avg_ipb = double(events.size()) / total_blocks;
      double avg_bpb = double(total_cond) / total_blocks;

      std::cout << std::setw(8) << LINEINST
                << std::setw(5) << N
                << std::setw(12) << total_blocks
                << std::fixed << std::setprecision(1)
                << std::setw(10) << avg_ipb
                << std::setprecision(2)
                << std::setw(9) << avg_bpb
                << std::endl;
    }
    std::cout << std::endl;
  }

  return 0;
}
