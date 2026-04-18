// Timing probe: full update_cycle datapath breakdown.
// Build via: make timing-probe
// Run via:   make run-timing-probe

#include "trace_reader.hpp"
#include "harcom.hpp"
#include "cbp.hpp"
#include "branch_predictor.hpp"

#include <algorithm>
#include <vector>

branch_predictor pred_instance;

void harcom_superuser::timing_debug(predictor &p, uint64_t next_time) {
    auto &pp = static_cast<branch_predictor &>(p);
    auto t0 = int64_t(next_time);

    std::vector<std::pair<std::string, uint64_t>> signals;
    auto d = [&](const char *label, uint64_t t) {
        signals.emplace_back(label, t);
    };

    fprintf(stderr, "--- update_cycle datapath (cycle_end=%lu) ---\n", next_time);

    // Pipeline shift regs
    constexpr uint64_t NT = std::tuple_size_v<decltype(pp.tables)>;
    uint64_t worst_current = 0;
    for (uint64_t i = 0; i < NT; i++) {
        worst_current = std::max(worst_current, pp.current_pred[i].time());
        worst_current = std::max(worst_current, pp.current_tag_hit[i].time());
        worst_current = std::max(worst_current, pp.current_sec[i].time());
    }
    d("actual_dir", pp.dbg_actual_dir.time());
    d("alloc_base", pp.dbg_alloc_base.time());
    d("alloc_target (one_hot)", pp.dbg_alloc_target.time());
    d("alt_pred (fold_or)", pp.dbg_alt_pred.time());
    d("altdiff", pp.dbg_altdiff.time());
    d("any_provider_wrong", pp.dbg_any_prov_wrong.time());
    d("candallocmask", pp.dbg_candallocmask.time());
    d("curr_sec_tag (after hash)", pp.dbg_curr_sec_tag.time());
    d("curr_sec_tag (reg)", pp.curr_sec_tag.time());
    {
        uint64_t wdf = 0, wdm = 0;
        for (uint64_t i = 0; i < NT; i++) {
            wdf = std::max(wdf, pp.dbg_decay_fire[i].time());
            wdm = std::max(wdm, pp.dbg_decay_merged[i].time());
        }
        d("decay_fire (worst)", wdf);
        d("decay_merged (worst)", wdm);
    }
    d("epoch_fire", pp.dbg_epoch_fire.time());
    d("fb_pred", pp.dbg_fb_pred.time());
    d("fb_write (gate)", pp.dbg_fb_write.time());
    d("final_pred (select)", pp.dbg_final_pred.time());
    d("fold_apply[0] (update)", pp.dbg_fold_apply.time());
    d("fold_compute[0] (update)", pp.dbg_fold_compute.time());
    d("fold_early_write[0]", pp.dbg_fold_early_write.time());
    d("fold_idx[0] (predict1)", pp.dbg_fold_idx.time());
    d("fold_read_in_compute[0]", pp.dbg_fold_read_in_compute.time());
    d("fold_tag[0] (predict1)", pp.dbg_fold_tag.time());
    d("full_hits", pp.dbg_full_hits.time());
    d("gh[0] (update)", pp.dbg_gh_fanout.time());
    d("has_alt", pp.dbg_has_alt.time());
    d("hist_input (update)", pp.dbg_hist_input.time());
    {
        uint64_t wh = 0;
        for (uint64_t i = 0; i < NT; i++)
            wh = std::max(wh, pp.dbg_hyst_write[i].time());
        d("hyst_write (worst)", wh);
    }
    d("inst_pc (predict1)", pp.dbg_inst_pc.time());
    d("match (concat)", pp.dbg_match.time());
    d("match1 (one_hot)", pp.dbg_match1.time());
    d("match2 (one_hot)", pp.dbg_match2.time());
    d("meta_use_alt", pp.dbg_meta_use_alt.time());
    {
        constexpr uint64_t MP = sizeof(pp.meta_pipe) / sizeof(pp.meta_pipe[0]);
        uint64_t worst_meta = 0;
        for (uint64_t i = 0; i < MP; i++)
            worst_meta = std::max(worst_meta, pp.meta_pipe[i].time());
        d("meta_pipe (worst)", worst_meta);
    }
    d("meta_write (gate)", pp.dbg_meta_write.time());
    d("mispredict", pp.dbg_mispredict.time());
    d("next_pc (raw)", pp.dbg_next_pc.time());
    d("noalloc", pp.dbg_noalloc.time());
    d("notumask", pp.dbg_notumask.time());
    d("p1_return (predict1)", pp.dbg_p1_return.time());
    d("pipe_shift (current_*)", worst_current);
    {
        constexpr uint64_t PN = std::remove_reference_t<decltype(pp.pred)>::size;
        uint64_t worst_pred = 0;
        for (uint64_t i = 0; i < PN; i++)
            worst_pred = std::max(worst_pred, pp.pred[i].time());
        d("pred[I] (scatter)", worst_pred);
    }
    {
        uint64_t wp = 0;
        for (uint64_t i = 0; i < NT; i++)
            wp = std::max(wp, pp.dbg_pred_write[i].time());
        d("pred_write (worst)", wp);
    }
    d("provider_pred (fold_or)", pp.dbg_provider_pred.time());
    d("provider_weak", pp.dbg_provider_weak.time());
    {
        uint64_t wt = 0;
        for (uint64_t i = 0; i < NT; i++)
            wt = std::max(wt, pp.dbg_tag_write[i].time());
        d("tag_write (worst)", wt);
    }
    d("true_block", pp.dbg_true_block.time());
    d("uclearmask", pp.dbg_uclearmask.time());
    d("use_alt", pp.dbg_use_alt.time());
    {
        uint64_t wu = 0;
        for (uint64_t i = 0; i < NT; i++)
            wu = std::max(wu, pp.dbg_u_write[i].time());
        d("u_write (worst)", wu);
    }

    // Already alphabetical — just print
    for (auto &[label, t] : signals) {
        fprintf(stderr, "  %-30s t=%lu  delta=%+ld\n", label.c_str(), t, int64_t(t) - t0);
    }
    fprintf(stderr, "\n");
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <trace> <name> <warmup> <measure>\n";
        return 1;
    }
    trace_reader reader(argv[1], argv[2]);
    harcom_superuser sim(reader, true);
    sim.run(pred_instance, std::stoull(argv[3]), std::stoull(argv[4]));
}
