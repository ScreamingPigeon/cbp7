#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "trace_reader.hpp"

static constexpr uint64_t WINDOW_SIZE = 1'000'000; // 1M instructions

struct BranchStats {
    uint64_t total = 0;
    uint64_t taken = 0;
    uint64_t flips = 0;
    bool last_outcome = false;
    bool has_last = false;
};

double entropy(uint64_t taken, uint64_t total) {
    if (total == 0) return 0.0;
    double p = (double)taken / total;
    if (p == 0.0 || p == 1.0) return 0.0;
    return -p * std::log2(p) - (1 - p) * std::log2(1 - p);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " trace.gz\n";
        return 1;
    }

    trace_reader reader(argv[1], "analysis");

    uint64_t total_insts = 0;
    uint64_t total_branches = 0;
    uint64_t taken_branches = 0;
    uint64_t backward_branches = 0;

    std::unordered_map<uint64_t, BranchStats> per_pc;
    std::unordered_map<uint64_t, uint64_t> pc_region_histogram;

    uint64_t window_inst = 0;
    uint64_t window_br = 0;
    uint64_t window_taken = 0;
    uint64_t window_id = 0;

    try {
        while (true) {
            instruction inst = reader.next_instruction();
            total_insts++;
            window_inst++;

            // PC region histogram (top 16 bits)
            uint64_t region = inst.pc >> 48;
            pc_region_histogram[region]++;

            if (inst.branch) {
                total_branches++;
                window_br++;

                if (inst.taken_branch) {
                    taken_branches++;
                    window_taken++;
                }

                if (inst.taken_branch && inst.next_pc < inst.pc)
                    backward_branches++;

                auto& stats = per_pc[inst.pc];
                stats.total++;
                if (inst.taken_branch)
                    stats.taken++;

                if (stats.has_last && stats.last_outcome != inst.taken_branch)
                    stats.flips++;

                stats.last_outcome = inst.taken_branch;
                stats.has_last = true;
            }

            // Window boundary
            if (window_inst >= WINDOW_SIZE) {
                double win_taken_rate = window_br ? 
                    (double)window_taken / window_br : 0.0;

                std::cout << "[Window " << window_id++ << "] "
                          << "Branches: " << window_br
                          << " TakenRate: " << std::fixed << std::setprecision(4)
                          << win_taken_rate
                          << " BranchDensity: "
                          << (double)window_br / window_inst
                          << "\n";

                window_inst = 0;
                window_br = 0;
                window_taken = 0;
            }
        }
    }
    catch (const out_of_instructions&) {
        // finished
    }

    std::cout << "\n==== GLOBAL STATS ====\n";
    std::cout << "Total instructions: " << total_insts << "\n";
    std::cout << "Total branches:     " << total_branches << "\n";
    std::cout << "Branch density:     "
              << (double)total_branches / total_insts << "\n";
    std::cout << "Taken rate:         "
              << (double)taken_branches / total_branches << "\n";
    std::cout << "Backward rate:      "
              << (double)backward_branches / total_branches << "\n";
    std::cout << "Unique branch PCs:  "
              << per_pc.size() << "\n";

    // Compute per-PC entropy and sort by frequency
    struct PCEntry {
        uint64_t pc;
        BranchStats stats;
        double ent;
    };

    std::vector<PCEntry> entries;

    for (auto& [pc, stats] : per_pc) {
        entries.push_back({
            pc,
            stats,
            entropy(stats.taken, stats.total)
        });
    }

    std::sort(entries.begin(), entries.end(),
              [](const PCEntry& a, const PCEntry& b) {
                  return a.stats.total > b.stats.total;
              });

    std::cout << "\n==== TOP 10 HOT BRANCHES ====\n";
    for (size_t i = 0; i < std::min<size_t>(10, entries.size()); i++) {
        auto& e = entries[i];
        std::cout << "PC: 0x" << std::hex << e.pc << std::dec
                  << " Count: " << e.stats.total
                  << " TakenRate: "
                  << (double)e.stats.taken / e.stats.total
                  << " Entropy: " << e.ent
                  << " FlipRate: "
                  << (double)e.stats.flips / e.stats.total
                  << "\n";
    }

    std::cout << "\n==== PC REGION HISTOGRAM (Top 8) ====\n";
    std::vector<std::pair<uint64_t, uint64_t>> regions(
        pc_region_histogram.begin(),
        pc_region_histogram.end());

    std::sort(regions.begin(), regions.end(),
              [](auto& a, auto& b) { return a.second > b.second; });

    for (size_t i = 0; i < std::min<size_t>(8, regions.size()); i++) {
        std::cout << "Region 0x"
                  << std::hex << regions[i].first
                  << std::dec << " Count: "
                  << regions[i].second << "\n";
    }

    return 0;
}