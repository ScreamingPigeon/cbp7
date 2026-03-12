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
    double p = static_cast<double>(taken) / total;
    if (p == 0.0 || p == 1.0) return 0.0;
    return -p * std::log2(p) - (1 - p) * std::log2(1 - p);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " trace.gz [region_shift]\n";
        std::cerr << "  region_shift: bits to shift PC right for region grouping\n";
        std::cerr << "    12 = 4KB pages (default, ARM v8)\n";
        std::cerr << "    21 = 2MB large pages\n";
        std::cerr << "    30 = 1GB sections\n";
        return 1;
    }

    // ARM v8 default: 12-bit shift = 4KB pages
    uint64_t region_shift = 12;
    if (argc > 2) {
        region_shift = std::stoull(argv[2]);
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
    std::unordered_map<uint64_t, uint64_t> window_region_histogram;

    // Memory access patterns
    std::unordered_map<uint64_t, uint64_t> jump_size_histogram;
    uint64_t last_pc = 0;
    bool has_last_pc = false;

    try {
        while (true) {
            instruction inst = reader.next_instruction();
            total_insts++;
            window_inst++;

            // PC region histogram (configurable shift)
            uint64_t region = inst.pc >> region_shift;
            pc_region_histogram[region]++;
            window_region_histogram[region]++;

            // Track PC jumps for memory access patterns
            if (has_last_pc) {
                uint64_t jump_size = (inst.pc > last_pc) ? (inst.pc - last_pc) : (last_pc - inst.pc);
                // Bucket jumps by magnitude (log scale)
                uint64_t bucket = 0;
                if (jump_size > 0) {
                    bucket = 63 - __builtin_clzll(jump_size);
                }
                jump_size_histogram[bucket]++;
            }
            last_pc = inst.pc;
            has_last_pc = true;

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
                    static_cast<double>(window_taken) / window_br : 0.0;

                std::cout << "[Window " << window_id << "] "
                          << "Branches: " << window_br
                          << " TakenRate: " << std::fixed << std::setprecision(4)
                          << win_taken_rate
                          << " BranchDensity: "
                          << static_cast<double>(window_br) / window_inst
                          << " | Top regions: ";

                // Show top 3 regions in this window
                std::vector<std::pair<uint64_t, uint64_t>> window_regions(
                    window_region_histogram.begin(),
                    window_region_histogram.end());
                std::sort(window_regions.begin(), window_regions.end(),
                         [](const auto& a, const auto& b) { return a.second > b.second; });

                for (size_t i = 0; i < std::min<size_t>(3, window_regions.size()); i++) {
                    std::cout << "0x" << std::hex << window_regions[i].first << std::dec
                              << "(" << window_regions[i].second << ") ";
                }
                std::cout << "\n";

                window_id++;
                window_inst = 0;
                window_br = 0;
                window_taken = 0;
                window_region_histogram.clear();
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
              << static_cast<double>(total_branches) / total_insts << "\n";
    std::cout << "Taken rate:         "
              << static_cast<double>(taken_branches) / total_branches << "\n";
    std::cout << "Backward rate:      "
              << static_cast<double>(backward_branches) / total_branches << "\n";
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
                  << static_cast<double>(e.stats.taken) / e.stats.total
                  << " Entropy: " << e.ent
                  << " FlipRate: "
                  << static_cast<double>(e.stats.flips) / e.stats.total
                  << "\n";
    }

    uint64_t region_size = 1ULL << region_shift;
    std::string region_unit = "B";
    uint64_t region_size_display = region_size;
    if (region_size >= (1ULL << 30)) {
        region_size_display >>= 30;
        region_unit = "GB";
    } else if (region_size >= (1ULL << 20)) {
        region_size_display >>= 20;
        region_unit = "MB";
    } else if (region_size >= (1ULL << 10)) {
        region_size_display >>= 10;
        region_unit = "KB";
    }

    std::cout << "\n==== MEMORY ACCESS PATTERNS ====\n";
    std::cout << "PC jump size distribution (bucket = log2 of bytes):\n";
    for (size_t i = 0; i <= 48; i++) {
        if (jump_size_histogram.count(i)) {
            uint64_t count = jump_size_histogram[i];
            uint64_t min_size = (i == 0) ? 0 : (1ULL << i);
            uint64_t max_size = (1ULL << (i + 1)) - 1;
            std::cout << "  [" << std::setw(2) << i << "] "
                      << std::setw(12) << min_size << "-"
                      << std::setw(12) << max_size << " bytes: "
                      << std::setw(10) << count << " ("
                      << std::fixed << std::setprecision(1)
                      << (100.0 * count / total_insts) << "%)\n";
        }
    }

    std::cout << "\n==== PC REGION HISTOGRAM (Top 16) ====\n";
    std::cout << "Region size: " << region_size_display << region_unit
              << " (shift=" << region_shift << ")\n";
    std::vector<std::pair<uint64_t, uint64_t>> regions(
        pc_region_histogram.begin(),
        pc_region_histogram.end());

    std::sort(regions.begin(), regions.end(),
              [](auto& a, auto& b) { return a.second > b.second; });

    for (size_t i = 0; i < std::min<size_t>(16, regions.size()); i++) {
        std::cout << "Region 0x"
                  << std::hex << regions[i].first
                  << std::dec << " Count: "
                  << regions[i].second << "\n";
    }

    return 0;
}