#pragma once
// GshareN-Ahead Performance Monitor — v3
// Gated by -DGSHARE_AHEAD_MONITOR -DCHEATING_MODE at compile time.
// All counters are plain C++ types; zero HARCOM cost.
//
// Fixes vs v2:
//   hysteresis: "should equal misp_weak" text corrected → "should equal mispreds"
//   hysteresis: confidence HIGH → MEDIUM (no same-cell repeat-miss evidence yet)
//   hysteresis: added per-cell consecutive-strong-miss streak tracking (CellStats)
//   hysteresis: ctr_hi_target_bit_flips noted as tautology; streak stats added
//   ctx_hash:   now includes index1 (XOR ghist ^ path ^ index1) for better aliasing proxy
//   dead-write: pending writes flushed as "unresolved" at summary time
//   section 10: accuracy per prev_num_branch bucket added (not just counts)

#include <cstdint>
#include <array>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <bit>
#include <cmath>

using u64 = uint64_t;

template <u64 LOGG, u64 GHIST, u64 N>
struct GshareAheadMonitor {

    static constexpr u64 LOGBANKS   = std::bit_width(N);
    static constexpr u64 LOGLANES   = std::bit_width(N > 1 ? N - 1 : 1ULL);
    static constexpr u64 LANES      = 1ULL << LOGLANES;
    static constexpr u64 BANKS      = 1ULL << LOGBANKS;
    static constexpr u64 INDEX_BITS = LOGG - LOGLANES - LOGBANKS;
    static constexpr u64 CTR_HI_ROWS = 1ULL << INDEX_BITS;
    static constexpr u64 INVALID_LANE = LANES; // sentinel for "no lane"

    static constexpr u64 WINDOW_SIZE = 100'000;

    // end_reason codes
    static constexpr u64 END_TAKEN      = 0;
    static constexpr u64 END_LINE_END   = 1;
    static constexpr u64 END_LAST_PRED  = 2;
    static constexpr u64 END_MISPREDICT = 3;
    static constexpr u64 END_ZERO_CBR   = 4; // branchless block — distinct from line_end

    // =========================================================================
    // Counter block
    // =========================================================================
    struct Counters {

        // --- 1. Core outcomes ---
        u64 blocks       = 0;
        u64 branches     = 0;
        u64 mispreds     = 0;
        u64 taken        = 0;
        u64 extra_cycles = 0;
        u64 fresh_blocks = 0;
        u64 cont_blocks  = 0;

        // --- 2. RAM activity ---
        u64 ctr_hi_reads             = 0;
        u64 ctr_hi_row_write_attempts = 0; // mispredict path entered (≤1 per mispredict)
        u64 ctr_hi_target_bit_flips   = 0; // direction bit actually changed
        u64 ctr_lo_reads             = 0;
        u64 ctr_lo_writes            = 0;

        // --- 3. Block formation ---
        std::array<u64, N + 2> num_cbr_hist{};
        u64 stop_zero_cbr   = 0;
        u64 stop_taken      = 0;
        u64 stop_line_end   = 0;
        u64 stop_last_pred  = 0;
        u64 stop_mispredict = 0;

        // --- 4. Accuracy breakdown ---
        u64 fresh_branches = 0;
        u64 fresh_correct  = 0;
        u64 cont_branches  = 0;
        u64 cont_correct   = 0;

        // --- 5. Rank stats (branch position within block) ---
        // These answer: do later-ranked branches predict worse?
        std::array<u64, N> rank_access{};
        std::array<u64, N> rank_misp{};

        // --- 6. Physical lane stats (from access_mask, XL-derived) ---
        // These answer: is the XL lane mapping balanced?
        std::array<u64, 16> phys_lane_access{};
        std::array<u64, 16> phys_lane_misp{};

        // --- 7. Path (bank) usage — split fresh vs served ---
        // fresh_path: incremented when a new path is selected for a fresh block
        // served_path: incremented for every block (fresh or continued) consuming predictions
        std::array<u64, 16> fresh_path_use{};
        std::array<u64, 16> fresh_path_misp{};
        std::array<u64, 16> served_path_use{};

        // --- 8. Hysteresis dynamics ---
        u64 misp_weak    = 0; // mispredict, phys lane was weak   → write attempted
        u64 misp_strong  = 0; // mispredict, phys lane was strong → no flip

        // --- 9. Global history update types ---
        u64 ghist_zero_cbr   = 0;
        u64 ghist_true_block = 0;
        u64 ghist_skipped    = 0;

        // --- 10. Continued-block root cause ---
        // split by how many branches the previous block had
        std::array<u64, N + 2> cont_by_prev_num_branch{};
        // accuracy per bucket: how many branches / how many correct in continued blocks
        std::array<u64, N + 2> cont_branches_by_prev{};
        std::array<u64, N + 2> cont_correct_by_prev{};

        void reset() { *this = Counters{}; }
    };

    Counters cum;
    Counters win;
    u64  window_num  = 0;
    bool hdr_printed = false;

    // =========================================================================
    // Per-block shadow state
    // =========================================================================
    bool shadow_is_fresh  = false;
    u64  shadow_path      = 0;
    u64  shadow_num_br    = 0;
    u64  shadow_index0    = 0; // index[0] from current fresh block — non-row-local PC context
    std::array<bool, N> shadow_pred{};
    std::array<bool, N> shadow_actual{};

    // State carried across blocks for continued-block root-cause analysis
    u64  prev_num_branch = 0;

    // Streak window accumulator (for CSV avg_streak column)
    u64  win_streak_sum    = 0;
    u64  win_flip_events   = 0;

    // =========================================================================
    // Shadow RAM: per-row hotspot and dead-write tracking for ctr_hi
    // =========================================================================
    struct RowStats {
        u64  reads          = 0;
        u64  write_attempts = 0; // how many times this row was the write target
        u64  mispreds       = 0;
        // dead-write tracking
        bool write_pending  = false; // last event was a write with no intervening read
        u64  dead_writes    = 0;     // writes overwritten before any read
        u64  live_writes    = 0;     // writes followed by ≥1 read before overwrite
        // context-collision tracking (proxy for destructive aliasing)
        bool ctx_valid      = false;
        u64  last_ctx_hash  = 0;     // hash(ghist, path) of last write
        u64  ctx_overwrites = 0;     // writes with a different ctx_hash than the prior write
    };
    std::unordered_map<u64, RowStats> row_shadow;

    // =========================================================================
    // Per-cell shadow: consecutive-strong-miss streak tracking for ctr_lo cells.
    // Key = (index1 << LOGBANKS) | path_val — the exact ctr_lo address.
    // Answers: how many times does the predictor miss on the same cell before it
    // finally gets a weak hit and is allowed to flip?
    // =========================================================================
    struct CellStats {
        u64 consecutive_strong = 0;    // current unbroken strong-miss run on this cell
        u64 total_strong       = 0;    // lifetime strong misses
        u64 total_weak         = 0;    // lifetime weak misses (flips)
        // Histogram of streak lengths that ended in a flip:
        //   streak_hist[k] = number of times a run of exactly k consecutive strong misses
        //   was terminated by a weak miss. Index capped at STREAK_HIST_SIZE-1.
        static constexpr u64 STREAK_HIST_SIZE = 17; // 0..15 exact, 16 = "≥16"
        std::array<u64, STREAK_HIST_SIZE> streak_hist{};
    };
    std::unordered_map<u64, CellStats> cell_shadow;

    // =========================================================================
    // Hook: on_predict1
    //   is_fresh  — true_block register value
    //   path_val  — new path (valid only when is_fresh; 0 otherwise)
    //   index0    — index[0] passed to ctr_hi.read (valid only when is_fresh)
    // =========================================================================
    void on_predict1(bool is_fresh, u64 path_val, u64 index0) {
        shadow_is_fresh = is_fresh;
        shadow_num_br   = 0;
        if (is_fresh) {
            shadow_path   = path_val;
            shadow_index0 = index0; // PC+ghist context of this block; non-local to index[1]
        }

        auto rec = [&](Counters &c) {
            c.blocks++;
            if (is_fresh) { c.fresh_blocks++; c.ctr_hi_reads++; }
            else            c.cont_blocks++;

            // Served-path usage (path that is active for this block's predictions)
            u64 pv = shadow_path < 16 ? shadow_path : 15u;
            c.served_path_use[pv]++;

            // Fresh-path usage (only when a new path is computed)
            if (is_fresh) {
                u64 fp = path_val < 16 ? path_val : 15u;
                c.fresh_path_use[fp]++;
            }

            // Continued-block root cause
            if (!is_fresh) {
                u64 pb = prev_num_branch < c.cont_by_prev_num_branch.size()
                       ? prev_num_branch : c.cont_by_prev_num_branch.size() - 1;
                c.cont_by_prev_num_branch[pb]++;
            }
        };
        rec(cum); rec(win);

        if (is_fresh) {
            auto &r = row_shadow[index0];
            // Dead-write accounting: if a write was pending and we're now reading,
            // that write was live (it got a read before overwrite)
            if (r.write_pending) { r.live_writes++; r.write_pending = false; }
            r.reads++;
        }
    }

    // =========================================================================
    // Hook: on_predict2 — capture prediction delivered at branch time
    //   rank     — num_branch at entry (0-indexed)
    //   pred_val — prediction bit returned
    // =========================================================================
    void on_predict2(u64 rank, bool pred_val) {
        if (rank < N) shadow_pred[rank] = pred_val;
    }

    // =========================================================================
    // Hook: on_update_condbr — record actual direction before num_branch increments
    //   rank  — num_branch at entry (0-indexed)
    //   taken — actual direction
    // =========================================================================
    void on_update_condbr(u64 rank, bool taken) {
        if (rank < N) shadow_actual[rank] = taken;
        shadow_num_br = rank + 1;

        auto rec = [&](Counters &c) {
            c.branches++;
            if (taken) c.taken++;
            if (rank < N) c.rank_access[rank]++;
        };
        rec(cum); rec(win);
    }

    // =========================================================================
    // Hook: on_ctr_lo_read — fires only on mispredict path, after weak[] computed
    //   phys_lane — decoded physical lane of mispredicted branch
    //   weak      — ctr_lo value at that lane (1=weak, 0=strong)
    //   index1    — index[1] (row address, same as passed to ctr_hi.write)
    //   path_val  — current path register value
    // =========================================================================
    void on_ctr_lo_read(u64 phys_lane, bool weak, u64 index1, u64 path_val) {
        auto rec = [&](Counters &c) {
            c.ctr_lo_reads++;
            if (weak) c.misp_weak++;
            else       c.misp_strong++;
        };
        rec(cum); rec(win);

        // Per-cell streak tracking.
        // Key includes phys_lane because ctr_lo is physically ctr_lo[LANES] —
        // each lane is a separate RAM. Omitting lane would merge distinct hardware cells.
        // Cell address = row * LANES * BANKS + lane * BANKS + path
        // i.e. (index1 * LANES + phys_lane) << LOGBANKS | path_val
        u64 cell_key = ((index1 * LANES + phys_lane) << LOGBANKS)
                       | (path_val & (BANKS - 1));
        auto &cs = cell_shadow[cell_key];

        // "Consecutive" means consecutive mispredict visits to this cell, not
        // all accesses. Correct predictions on this cell do not appear here, so
        // a streak of k means k strong-miss mispredicts in a row on this cell
        // with no weak-miss flip in between. Non-mispredict accesses are invisible.
        if (weak) {
            u64 streak = cs.consecutive_strong;
            u64 hist_idx = streak < CellStats::STREAK_HIST_SIZE - 1
                         ? streak : CellStats::STREAK_HIST_SIZE - 1;
            cs.streak_hist[hist_idx]++;
            cs.consecutive_strong = 0;
            cs.total_weak++;
            win_streak_sum  += streak;
            win_flip_events++;
        } else {
            cs.consecutive_strong++;
            cs.total_strong++;
        }
    }

    // =========================================================================
    // Hook: on_ctr_hi_write — called OUTSIDE execute_if, after the real write
    //   index1    — index[1] passed to ctr_hi.write (the row being updated)
    //   flipped   — true if the target direction bit actually changed value
    //   ctx_hash  — hash(ghist ^ path ^ shadow_index0) for context-collision tracking.
    //               shadow_index0 = index[0] of the current block (PC+ghist of the
    //               block being fetched), which is non-local to index[1], so it adds
    //               real discrimination within a row across different fetch contexts.
    // =========================================================================
    void on_ctr_hi_write(u64 index1, bool flipped, u64 ctx_hash) {
        auto rec = [&](Counters &c) {
            c.ctr_hi_row_write_attempts++;
            if (flipped) c.ctr_hi_target_bit_flips++;
        };
        rec(cum); rec(win);

        auto &r = row_shadow[index1];
        // Dead-write: if a write was already pending (no read in between), it died
        if (r.write_pending) r.dead_writes++;
        r.write_pending   = true;
        r.write_attempts++;

        // Context-collision tracking
        if (r.ctx_valid && r.last_ctx_hash != ctx_hash)
            r.ctx_overwrites++;
        r.last_ctx_hash = ctx_hash;
        r.ctx_valid     = true;
    }

    // =========================================================================
    // Hook: on_update_cycle — finalize block-level accounting
    //   mispredict    — is_mispredict
    //   num_branch    — branches in this block
    //   index1        — index[1] at cycle end (for mispred row tracking)
    //   zero_cbr      — num_branch == 0
    //   true_block_out— true_block value after all updates
    //   end_reason    — END_* constant
    //   num_lo_writes — number of ctr_lo writes (= num_branch for CBR blocks)
    //   access_mask   — bitmask of physical lanes accessed (from access[])
    //   misp_lane     — physical lane of mispredicted branch (INVALID_LANE if none)
    // =========================================================================
    void on_update_cycle(bool mispredict, u64 num_branch, u64 index1,
                         bool zero_cbr, bool true_block_out,
                         u64 end_reason, u64 num_lo_writes,
                         u64 access_mask, u64 misp_lane) {

        auto rec_block = [&](Counters &c) {
            u64 hi = num_branch < c.num_cbr_hist.size() ? num_branch
                                                        : c.num_cbr_hist.size() - 1;
            c.num_cbr_hist[hi]++;
            c.ctr_lo_writes += num_lo_writes;

            switch (end_reason) {
                case END_ZERO_CBR:   c.stop_zero_cbr++;   break;
                case END_TAKEN:      c.stop_taken++;       break;
                case END_LINE_END:   c.stop_line_end++;    break;
                case END_LAST_PRED:  c.stop_last_pred++;   break;
                case END_MISPREDICT: c.stop_mispredict++;  break;
            }

            if (mispredict) { c.mispreds++; c.extra_cycles++; }

            // Fresh-path mispred (attribute mispredict to the fresh path used)
            if (shadow_is_fresh) {
                u64 fp = shadow_path < 16 ? shadow_path : 15u;
                if (mispredict) c.fresh_path_misp[fp]++;
            }

            if (zero_cbr)            c.ghist_zero_cbr++;
            else if (true_block_out) c.ghist_true_block++;
            else                     c.ghist_skipped++;

            // Physical lane access / mispred from access_mask
            for (u64 l = 0; l < LANES && l < 16; l++) {
                if (access_mask & (u64(1) << l)) {
                    c.phys_lane_access[l]++;
                    if (mispredict && l == misp_lane)
                        c.phys_lane_misp[l]++;
                }
            }
        };
        rec_block(cum); rec_block(win);

        // Per-branch outcome recording
        for (u64 r = 0; r < num_branch && r < N; r++) {
            bool correct = (shadow_pred[r] == shadow_actual[r]);
            bool is_misp = mispredict && (r + 1 == num_branch);

            auto rec_br = [&](Counters &c) {
                if (shadow_is_fresh) {
                    c.fresh_branches++;
                    if (correct) c.fresh_correct++;
                } else {
                    c.cont_branches++;
                    if (correct) c.cont_correct++;
                    // accuracy per prev_num_branch bucket
                    u64 pb = prev_num_branch < c.cont_branches_by_prev.size()
                           ? prev_num_branch : c.cont_branches_by_prev.size() - 1;
                    c.cont_branches_by_prev[pb]++;
                    if (correct) c.cont_correct_by_prev[pb]++;
                }
                if (is_misp && r < N) c.rank_misp[r]++;
            };
            rec_br(cum); rec_br(win);
        }

        if (mispredict && index1 < CTR_HI_ROWS)
            row_shadow[index1].mispreds++;

        // Save state for next block's continued-block root-cause analysis
        prev_num_branch = num_branch;

        if (win.branches >= WINDOW_SIZE) {
            print_window(std::cerr);
            win.reset();
            win_streak_sum  = 0;
            win_flip_events = 0;
            window_num++;
        }
    }

    // =========================================================================
    // Window CSV output
    // =========================================================================
    void print_window(std::ostream &os) {
        if (!hdr_printed) {
            os << "# win,branches,misp%,xcyc%,"
               << "ctr_hi_rd,br/rd,ok/rd,"
               << "fresh_acc%,cont_acc%,"
               << "stop_0cbr%,stop_lend%,stop_lp%,stop_tk%,stop_mp%,"
               << "weak_miss%,avg_streak,";
            // avg_streak = avg consecutive strong misses before a flip this window.
            // Replaces the tautological flip% (which equals weak_miss% by design).
            for (u64 b = 0; b < BANKS && b < 16; b++) os << "fresh_path" << b << "%,";
            for (u64 l = 0; l < LANES && l < 16; l++) os << "plane" << l << "_mp%,";
            os << "\n";
            hdr_printed = true;
        }
        const auto &w = win;
        auto pct = [](u64 n, u64 d) { return d > 0 ? 100.0 * n / d : 0.0; };
        double br_rd = w.ctr_hi_reads ? (double)w.branches / w.ctr_hi_reads : 0.0;
        double ok_rd = w.ctr_hi_reads ? (double)(w.branches - w.mispreds) / w.ctr_hi_reads : 0.0;
        u64 total_fresh_path = 0;
        for (u64 b = 0; b < BANKS && b < 16; b++) total_fresh_path += w.fresh_path_use[b];

        os << std::fixed << std::setprecision(1);
        os << window_num << "," << w.branches << ","
           << pct(w.mispreds,     w.branches) << ","
           << pct(w.extra_cycles, w.branches) << ","
           << w.ctr_hi_reads << ","
           << br_rd << "," << ok_rd << ","
           << pct(w.fresh_correct, w.fresh_branches) << ","
           << pct(w.cont_correct,  w.cont_branches)  << ","
           << pct(w.stop_zero_cbr,   w.blocks) << ","
           << pct(w.stop_line_end,   w.blocks) << ","
           << pct(w.stop_last_pred,  w.blocks) << ","
           << pct(w.stop_taken,      w.blocks) << ","
           << pct(w.stop_mispredict, w.blocks) << ","
           << pct(w.misp_weak, w.mispreds) << ","
           << (win_flip_events > 0 ? (double)win_streak_sum / win_flip_events : 0.0) << ",";
        for (u64 b = 0; b < BANKS && b < 16; b++)
            os << pct(w.fresh_path_use[b], total_fresh_path) << ",";
        for (u64 l = 0; l < LANES && l < 16; l++)
            os << pct(w.phys_lane_misp[l], w.phys_lane_access[l] > 0
                                          ? w.phys_lane_access[l] : 1) << ",";
        os << "\n";
    }

    // =========================================================================
    // End-of-trace summary
    // =========================================================================
    void print_summary(std::ostream &os = std::cerr) const {
        auto pct = [](u64 n, u64 d) -> double { return d > 0 ? 100.0 * n / d : 0.0; };
        const auto &c = cum;

        double br_rd  = c.ctr_hi_reads ? (double)c.branches / c.ctr_hi_reads : 0.0;
        double ok_rd  = c.ctr_hi_reads ? (double)(c.branches - c.mispreds) / c.ctr_hi_reads : 0.0;
        double fresh_misp_pct = pct(c.fresh_branches - c.fresh_correct, c.fresh_branches);
        double cont_misp_pct  = pct(c.cont_branches  - c.cont_correct,  c.cont_branches);
        double stale_delta    = cont_misp_pct - fresh_misp_pct;

        os << "\n=== GshareN-Ahead Monitor v2"
           << " (LOGG=" << LOGG << " GHIST=" << GHIST << " N=" << N
           << " BANKS=" << BANKS << " LANES=" << LANES
           << " INDEX_BITS=" << INDEX_BITS << " CTR_HI_ROWS=" << CTR_HI_ROWS << ") ===\n";

        // -----------------------------------------------------------------
        // 1. Core outcomes
        // -----------------------------------------------------------------
        os << "\n[1. Core Outcomes]\n";
        os << std::fixed << std::setprecision(2);
        os << "  Blocks:          " << c.blocks
           << "  Fresh: " << c.fresh_blocks << " (" << pct(c.fresh_blocks, c.blocks) << "%)"
           << "  Continued: " << c.cont_blocks << " (" << pct(c.cont_blocks, c.blocks) << "%)\n";
        os << "  Branches:        " << c.branches
           << "  Taken: " << pct(c.taken, c.branches) << "%\n";
        os << "  Mispredictions:  " << c.mispreds
           << " (" << pct(c.mispreds, c.branches) << "%)\n";

        // -----------------------------------------------------------------
        // 2. RAM activity
        // -----------------------------------------------------------------
        os << "\n[2. RAM Activity (VFS cost proxy)]\n";
        os << "  ctr_hi reads:              " << c.ctr_hi_reads << "\n";
        os << "  ctr_hi row write attempts: " << c.ctr_hi_row_write_attempts
           << "  (should be ≤ mispreds=" << c.mispreds << ")\n";
        os << "  ctr_hi target bit flips:   " << c.ctr_hi_target_bit_flips
           << " (" << pct(c.ctr_hi_target_bit_flips, c.ctr_hi_row_write_attempts) << "% of write attempts)\n";
        os << "  ctr_lo reads:              " << c.ctr_lo_reads
           << "  writes: " << c.ctr_lo_writes << "\n";
        os << "  branches / ctr_hi_read:    " << br_rd
           << "  (bundle reuse depth)\n";
        os << "  correct_br / ctr_hi_read:  " << ok_rd
           << "  (direct VFS read-payoff proxy)\n";
        os << "  ctr_lo_writes / branch:    "
           << (c.branches ? (double)c.ctr_lo_writes / c.branches : 0.0) << "\n";

        // -----------------------------------------------------------------
        // 3. Block formation
        // -----------------------------------------------------------------
        os << "\n[3. Block Formation]\n";
        u64 total_blocks = c.blocks;
        os << "  Stop reasons:\n";
        os << "    zero_cbr  (branchless block):  " << pct(c.stop_zero_cbr,   total_blocks) << "%\n";
        os << "    line_end  (cache line limit):  " << pct(c.stop_line_end,   total_blocks) << "%\n";
        os << "    last_pred (ran out of slots):  " << pct(c.stop_last_pred,  total_blocks) << "%\n";
        os << "    taken     (taken branch exit): " << pct(c.stop_taken,      total_blocks) << "%\n";
        os << "    mispredict:                    " << pct(c.stop_mispredict, total_blocks) << "%\n";
        if (c.stop_last_pred > total_blocks / 10)
            os << "  ** last_pred > 10% of blocks: N=" << N << " may be too small **\n";
        os << "  Branches-per-block histogram:\n  Cnt: ";
        for (u64 i = 0; i <= N + 1; i++) os << std::setw(7) << i;
        os << "\n  Blk: ";
        for (u64 i = 0; i <= N + 1; i++) os << std::setw(7) << c.num_cbr_hist[i];
        os << "\n";

        // -----------------------------------------------------------------
        // 4. Accuracy breakdown
        // -----------------------------------------------------------------
        os << "\n[4. Accuracy Breakdown]  [confidence: HIGH]\n";
        os << "  Fresh blocks:     " << pct(c.fresh_correct, c.fresh_branches) << "% correct"
           << " (" << c.fresh_correct << "/" << c.fresh_branches << ")\n";
        os << "  Continued blocks: " << pct(c.cont_correct, c.cont_branches) << "% correct"
           << " (" << c.cont_correct << "/" << c.cont_branches << ")\n";
        os << "  Stale-reuse miss delta: +" << stale_delta << "pp on continued blocks";
        if (stale_delta > 1.0) os << "  ** stale continuation hurts — revisit true_block policy **";
        os << "\n";

        // -----------------------------------------------------------------
        // 5. Rank stats (branch position within block)
        // -----------------------------------------------------------------
        os << "\n[5. Branch-Rank Stats]  [confidence: HIGH]\n";
        os << "  Rank:  ";
        for (u64 r = 0; r < N; r++) os << std::setw(7) << r;
        os << "\n  Acc:   ";
        for (u64 r = 0; r < N; r++) os << std::setw(6) << pct(c.rank_access[r] - c.rank_misp[r], c.rank_access[r]) << "%";
        os << "\n  Misp:  ";
        for (u64 r = 0; r < N; r++) os << std::setw(7) << c.rank_misp[r];
        os << "\n  Misp%: ";
        for (u64 r = 0; r < N; r++) os << std::setw(6) << pct(c.rank_misp[r], c.rank_access[r]) << "%";
        os << "\n";

        // -----------------------------------------------------------------
        // 6. Physical lane utilization (XL-derived, from access_mask)
        // -----------------------------------------------------------------
        os << "\n[6. Physical Lane Utilization (XL mapping)]  [confidence: HIGH]\n";
        os << "  Lane | Access   | Acc%    | Mispr    | Misp%\n";
        os << "  -----+----------+---------+----------+------\n";
        for (u64 l = 0; l < LANES && l < 16; l++) {
            os << "  " << std::setw(4) << l
               << " |" << std::setw(9) << c.phys_lane_access[l]
               << " |" << std::setw(8) << pct(c.phys_lane_access[l], c.branches) << "%"
               << " |" << std::setw(9) << c.phys_lane_misp[l]
               << " |" << std::setw(6) << pct(c.phys_lane_misp[l], c.phys_lane_access[l]) << "%\n";
        }

        // -----------------------------------------------------------------
        // 7. Path (bank) utilization — fresh vs served
        // -----------------------------------------------------------------
        os << "\n[7. Path (Bank) Utilization]  [confidence: MEDIUM — served path conflates banks]\n";
        u64 total_fresh = 0, total_served = 0;
        for (u64 b = 0; b < BANKS && b < 16; b++) {
            total_fresh  += c.fresh_path_use[b];
            total_served += c.served_path_use[b];
        }
        double fresh_ent = 0.0, served_ent = 0.0;
        os << "  Bank | Fresh%  | FreshMisp% | Served%\n";
        os << "  -----+---------+------------+--------\n";
        for (u64 b = 0; b < BANKS && b < 16; b++) {
            double fp = total_fresh  > 0 ? (double)c.fresh_path_use[b]  / total_fresh  : 0.0;
            double sp = total_served > 0 ? (double)c.served_path_use[b] / total_served : 0.0;
            if (fp > 0.0) fresh_ent  -= fp * std::log2(fp);
            if (sp > 0.0) served_ent -= sp * std::log2(sp);
            os << "  " << std::setw(4) << b
               << " |" << std::setw(8) << pct(c.fresh_path_use[b], total_fresh) << "%"
               << " |" << std::setw(11) << pct(c.fresh_path_misp[b], c.fresh_path_use[b]) << "%"
               << " |" << std::setw(8) << pct(c.served_path_use[b], total_served) << "%\n";
        }
        double ideal_ent = std::log2((double)BANKS);
        os << std::setprecision(2);
        os << "  Fresh-path entropy:  " << fresh_ent  << " / " << ideal_ent << " (ideal)";
        if (fresh_ent < 0.75 * ideal_ent)
            os << "  ** skewed — bank distribution uneven; fewer banks or different XB may help **";
        os << "\n";
        os << "  Served-path entropy: " << served_ent << " / " << ideal_ent << " (ideal)\n";

        // -----------------------------------------------------------------
        // 8. Hysteresis dynamics
        // -----------------------------------------------------------------
        os << "\n[8. Hysteresis Dynamics]  [confidence: MEDIUM]\n";
        os << "  Note: raw weak/strong counts are trustworthy; the \"too sticky\" recommendation\n"
           << "  requires same-cell repeat-miss evidence — see streak stats below.\n";
        u64 total_misp_with_cbr = c.misp_weak + c.misp_strong;
        os << "  Mispredict + weak   (flip):    "
           << c.misp_weak << " (" << pct(c.misp_weak, total_misp_with_cbr) << "%)\n";
        os << "  Mispredict + strong (no flip): "
           << c.misp_strong << " (" << pct(c.misp_strong, total_misp_with_cbr) << "%)\n";
        os << "  ctr_hi row write attempts: " << c.ctr_hi_row_write_attempts
           << "  (= mispreds=" << c.mispreds << " — one write per mispredict, gated by nothing)\n";
        os << "  ctr_hi target bit flips:   " << c.ctr_hi_target_bit_flips
           << "  (= misp_weak=" << c.misp_weak << " — tautological: every weak miss flips by design)\n";
        if (total_misp_with_cbr > 0 && c.misp_strong > total_misp_with_cbr / 2)
            os << "  Raw signal: strong-miss rate > 50% — consistent with sticky hysteresis,\n"
               << "  but not conclusive without streak evidence (see below).\n";

        // Per-cell streak analysis — this is the counterfactual evidence
        {
            u64 total_cells    = cell_shadow.size();
            u64 total_streaks  = 0; // total flip events (streak endings)
            u64 streak_sum     = 0; // sum of streak lengths at flip time
            u64 max_streak     = 0;
            std::array<u64, CellStats::STREAK_HIST_SIZE> agg_hist{};
            u64 cells_with_long_streaks = 0; // cells with max streak ≥ 4
            for (auto &[key, cs] : cell_shadow) {
                u64 cell_max = 0;
                for (u64 k = 0; k < CellStats::STREAK_HIST_SIZE; k++) {
                    agg_hist[k] += cs.streak_hist[k];
                    total_streaks += cs.streak_hist[k];
                    streak_sum    += cs.streak_hist[k] * k;
                    if (cs.streak_hist[k] > 0 && k > cell_max) cell_max = k;
                }
                if (cell_max > max_streak) max_streak = cell_max;
                if (cell_max >= 4) cells_with_long_streaks++;
            }
            double avg_streak = total_streaks > 0 ? (double)streak_sum / total_streaks : 0.0;
            os << "  Per-cell streak analysis (consecutive strong misses before each flip):\n";
            os << "    Cells tracked: " << total_cells
               << "  Total flip events: " << total_streaks
               << "  Avg streak before flip: " << std::setprecision(2) << avg_streak << "\n";
            os << "    Max streak: " << max_streak
               << "  Cells with streak ≥ 4: " << cells_with_long_streaks
               << " (" << pct(cells_with_long_streaks, total_cells) << "% of cells)\n";
            os << "    Streak histogram (length → flip count):\n    ";
            for (u64 k = 0; k < CellStats::STREAK_HIST_SIZE; k++) {
                if (k == CellStats::STREAK_HIST_SIZE - 1)
                    os << "  ≥" << k << ":";
                else
                    os << "  " << std::setw(2) << k << ":";
                os << std::setw(6) << agg_hist[k];
            }
            os << "\n";
            if (avg_streak >= 2.0)
                os << "  ** Avg streak ≥ 2: cells typically miss multiple times before flipping —\n"
                   << "     hysteresis is likely too sticky. Medium-to-high confidence.\n";
            else if (avg_streak >= 1.0)
                os << "  Avg streak 1–2: some stickiness but not severe. Review per-cell distribution.\n";
            else
                os << "  Avg streak < 1: most flips happen after the first strong miss.\n"
                   << "  Hysteresis may be protecting against noise, not causing harm.\n";
        }

        // -----------------------------------------------------------------
        // 9. Global history update types
        // -----------------------------------------------------------------
        os << "\n[9. Global History Update Types]  [confidence: HIGH]\n";
        u64 total_gh = c.ghist_zero_cbr + c.ghist_true_block + c.ghist_skipped;
        os << "  Zero-CBR (branchless block, next_pc injected): "
           << c.ghist_zero_cbr   << " (" << pct(c.ghist_zero_cbr,   total_gh) << "%)\n";
        os << "  True-block (normal update):                    "
           << c.ghist_true_block << " (" << pct(c.ghist_true_block, total_gh) << "%)\n";
        os << "  Skipped (false block, true_block=0):           "
           << c.ghist_skipped   << " (" << pct(c.ghist_skipped,   total_gh) << "%)\n";
        if (c.ghist_zero_cbr > total_gh / 4)
            os << "  ** Zero-CBR injections > 25%: branchless-block PC noise"
               << " may pollute history **\n";

        // -----------------------------------------------------------------
        // 10. Continued-block root cause
        // -----------------------------------------------------------------
        os << "\n[10. Continued-Block Root Cause (by prev block's num_branch)]"
           << "  [confidence: HIGH]\n";
        if (c.cont_blocks > 0) {
            os << "  Note: continued blocks only arise when prev block's last CBR was NOT TAKEN,\n"
               << "  no mispredict, no last_pred, no line_end — so prev_num_branch drives variance.\n";
            os << "  Prev #CBR:  ";
            for (u64 i = 0; i <= N + 1; i++) os << std::setw(7) << i;
            os << "\n  Cont blks:  ";
            for (u64 i = 0; i <= N + 1; i++) os << std::setw(7) << c.cont_by_prev_num_branch[i];
            os << "\n  Cont%:      ";
            for (u64 i = 0; i <= N + 1; i++)
                os << std::setw(6) << pct(c.cont_by_prev_num_branch[i], c.cont_blocks) << "%";
            os << "\n";
            // Accuracy by bucket — the key question: which prev-block sizes produce bad predictions?
            os << "  Cont br:    ";
            for (u64 i = 0; i <= N + 1; i++) os << std::setw(7) << c.cont_branches_by_prev[i];
            os << "\n  Cont acc%:  ";
            for (u64 i = 0; i <= N + 1; i++)
                os << std::setw(6) << pct(c.cont_correct_by_prev[i], c.cont_branches_by_prev[i]) << "%";
            os << "\n  Cont misp%: ";
            for (u64 i = 0; i <= N + 1; i++) {
                u64 br = c.cont_branches_by_prev[i];
                u64 ok = c.cont_correct_by_prev[i];
                os << std::setw(6) << pct(br > ok ? br - ok : 0, br) << "%";
            }
            os << "\n";
        } else {
            os << "  No continued blocks observed.\n";
        }

        // -----------------------------------------------------------------
        // 11. ctr_hi Row Hotspots & Dead-Write Analysis
        // -----------------------------------------------------------------
        os << "\n[11. ctr_hi Row Hotspots & Dead-Write Analysis]"
           << "  [confidence: MEDIUM — ctx_hash is (ghist^path^index0_of_current_block)]\n";
        u64 rows_accessed  = row_shadow.size();
        u64 total_rw       = 0, total_rr = 0;
        u64 total_dead     = 0, total_live = 0, total_unresolved = 0;
        u64 total_ctx_owt  = 0, hot_rows   = 0;
        u64 max_writes     = 0;
        for (auto &[idx, r] : row_shadow) {
            total_rw       += r.write_attempts;
            total_rr       += r.reads;
            total_dead     += r.dead_writes;
            total_live     += r.live_writes;
            total_ctx_owt  += r.ctx_overwrites;
            // Rows that still have write_pending at end-of-trace are unresolved
            if (r.write_pending) total_unresolved++;
            if (r.write_attempts > max_writes) max_writes = r.write_attempts;
            if (r.write_attempts > 100) hot_rows++;
        }
        os << "  Rows ever accessed:   " << rows_accessed << " / " << CTR_HI_ROWS
           << " (" << pct(rows_accessed, CTR_HI_ROWS) << "%)\n";
        os << "  Avg writes/row:       "
           << (rows_accessed ? (double)total_rw / rows_accessed : 0.0) << "\n";
        os << "  Max writes/row:       " << max_writes
           << "  Hot rows (>100 writes): " << hot_rows << "\n";
        os << "  Dead writes:          " << total_dead
           << " (" << pct(total_dead, total_rw) << "% of write attempts)"
           << "  — overwritten before any read\n";
        os << "  Live writes:          " << total_live
           << " (" << pct(total_live, total_rw) << "% of write attempts)"
           << "  — consumed by ≥1 read before overwrite\n";
        os << "  Unresolved at EOT:    " << total_unresolved << " rows"
           << "  — pending write never read; neither dead nor live\n";
        if (total_rw > 0 && total_dead > total_rw / 3)
            os << "  ** Dead-write rate > 33%: many table updates are wasted churn;"
               << " consider write throttling or stronger history **\n";
        os << "  Ctx overwrites:       " << total_ctx_owt
           << " (" << pct(total_ctx_owt, total_rw) << "% of writes)"
           << "  — different (ghist^path^index0) context wrote same row\n";
        {
            std::vector<std::pair<u64,u64>> hot; // {write_attempts, index}
            for (auto &[idx, r] : row_shadow) hot.push_back({r.write_attempts, idx});
            std::sort(hot.rbegin(), hot.rend());
            u64 show = std::min<u64>(10, hot.size());
            if (show > 0) {
                os << "  Top " << show << " rows by write_attempts:\n";
                for (u64 i = 0; i < show; i++) {
                    u64 idx = hot[i].second;
                    auto &r = row_shadow.at(idx);
                    os << "    [" << std::setw(6) << idx << "]"
                       << "  writes="  << std::setw(5) << r.write_attempts
                       << "  reads="   << std::setw(5) << r.reads
                       << "  dead="    << std::setw(5) << r.dead_writes
                       << "  ctx_owt=" << std::setw(4) << r.ctx_overwrites
                       << "  mispreds="<< std::setw(5) << r.mispreds << "\n";
                }
            }
        }

        // -----------------------------------------------------------------
        // 12. VFS Proxy Summary
        // -----------------------------------------------------------------
        os << "\n[12. VFS Proxy Summary & Diagnoses]\n";
        os << std::fixed << std::setprecision(2);

        auto diag = [&](const char *metric, double val, const char *fmt,
                        bool flag, const char *advice, const char *conf) {
            os << "  " << std::left << std::setw(38) << metric
               << std::right << std::setw(7) << val << fmt;
            if (flag) os << "  [" << conf << "] → " << advice;
            os << "\n";
        };

        diag("branches/ctr_hi_read:", br_rd, "  ",
             br_rd < 2.0,
             "low reuse — bundle read cost amortized over <2 branches; tune block formation",
             "HIGH");
        diag("correct_branches/ctr_hi_read:", ok_rd, "  ",
             false, "", "HIGH");
        diag("strong-miss rate (%):",
             pct(c.misp_strong, total_misp_with_cbr), "%",
             total_misp_with_cbr > 0 && c.misp_strong > total_misp_with_cbr / 2,
             "consistent with sticky hysteresis — see streak stats in sec.8 for confidence",
             "MEDIUM");
        // Note: flip rate is a tautology (every weak miss flips by design).
        // The meaningful metric is avg_streak in sec.8.
        diag("flip rate of weak mispreds (tautological):",
             pct(c.ctr_hi_target_bit_flips, c.misp_weak > 0 ? c.misp_weak : 1), "%",
             false, "", "N/A");
        diag("last_pred stop rate (%):",
             pct(c.stop_last_pred, total_blocks), "%",
             c.stop_last_pred > total_blocks / 10,
             "N too small — increase N",
             "HIGH");
        diag("stale-reuse miss delta (pp):",
             stale_delta, "pp",
             stale_delta > 1.0,
             "revisit true_block policy or force fresh read more often",
             "HIGH");
        diag("dead-write rate (%):",
             pct(total_dead, total_rw), "%",
             total_rw > 0 && total_dead > total_rw / 3,
             "many updates wasted — consider write throttling",
             "MEDIUM");
        diag("ctx overwrite rate (%):",
             pct(total_ctx_owt, total_rw), "%",
             total_rw > 0 && total_ctx_owt > total_rw / 2,
             "high aliasing pressure — try more index bits or better hash",
             "MEDIUM");
        diag("fresh-path entropy:",
             fresh_ent, " bits",
             fresh_ent < 0.75 * ideal_ent,
             "bank selection skewed — fewer banks or different XB permutation",
             "MEDIUM");
        diag("zero-CBR ghist injections (%):",
             pct(c.ghist_zero_cbr, total_gh), "%",
             c.ghist_zero_cbr > total_gh / 4,
             "branchless-block PC noise may pollute history",
             "MEDIUM");

        os << "=== End GshareN-Ahead Monitor v3 ===\n\n";
    }
};
