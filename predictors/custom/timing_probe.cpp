// Timing probe: full update_cycle datapath breakdown.
// Build via: make timing-probe
// Run via:   make run-timing-probe

#include "trace_reader.hpp"
#include "harcom.hpp"
#include "cbp.hpp"
#include "branch_predictor.hpp"

branch_predictor pred_instance;

void harcom_superuser::timing_debug(predictor &p, uint64_t next_time) {
    auto &pp = static_cast<branch_predictor &>(p);
    auto t0 = int64_t(next_time);

    auto d = [&](const char *label, uint64_t t) {
        fprintf(stderr, "  %-25s t=%lu  delta=%+ld\n", label, t, int64_t(t) - t0);
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
    d("pipe_shift (current_*)", worst_current);
    d("next_pc (raw)", pp.dbg_next_pc.time());
    d("curr_sec_tag (after hash)", pp.dbg_curr_sec_tag.time());
    d("curr_sec_tag (reg)", pp.curr_sec_tag.time());

    // Resolution chain
    d("full_hits", pp.dbg_full_hits.time());
    d("fb_pred", pp.dbg_fb_pred.time());
    d("match (concat)", pp.dbg_match.time());
    d("match1 (one_hot)", pp.dbg_match1.time());
    d("match2 (one_hot)", pp.dbg_match2.time());
    d("provider_pred (fold_or)", pp.dbg_provider_pred.time());
    d("alt_pred (fold_or)", pp.dbg_alt_pred.time());
    d("provider_weak", pp.dbg_provider_weak.time());
    d("has_alt", pp.dbg_has_alt.time());
    d("meta_use_alt", pp.dbg_meta_use_alt.time());
    d("use_alt", pp.dbg_use_alt.time());
    d("final_pred (select)", pp.dbg_final_pred.time());

    // Resolution gaps
    d("altdiff", pp.dbg_altdiff.time());
    d("actual_dir", pp.dbg_actual_dir.time());
    d("any_provider_wrong", pp.dbg_any_prov_wrong.time());

    // Allocation chain
    d("alloc_base", pp.dbg_alloc_base.time());
    d("notumask", pp.dbg_notumask.time());
    d("candallocmask", pp.dbg_candallocmask.time());
    d("alloc_target (one_hot)", pp.dbg_alloc_target.time());
    d("noalloc", pp.dbg_noalloc.time());
    d("uclearmask", pp.dbg_uclearmask.time());

    // Train write gates (worst across tables)
    d("fb_write (gate)", pp.dbg_fb_write.time());
    d("meta_write (gate)", pp.dbg_meta_write.time());
    {
        uint64_t wp = 0, wh = 0, wt = 0, wu = 0, wdf = 0, wdm = 0;
        for (uint64_t i = 0; i < NT; i++) {
            wp = std::max(wp, pp.dbg_pred_write[i].time());
            wh = std::max(wh, pp.dbg_hyst_write[i].time());
            wt = std::max(wt, pp.dbg_tag_write[i].time());
            wu = std::max(wu, pp.dbg_u_write[i].time());
            wdf = std::max(wdf, pp.dbg_decay_fire[i].time());
            wdm = std::max(wdm, pp.dbg_decay_merged[i].time());
        }
        d("pred_write (worst)", wp);
        d("hyst_write (worst)", wh);
        d("tag_write (worst)", wt);
        d("u_write (worst)", wu);
        d("decay_fire (worst)", wdf);
        d("decay_merged (worst)", wdm);
    }

    // Epoch
    d("epoch_fire", pp.dbg_epoch_fire.time());

    // Existing signals
    d("mispredict", pp.dbg_mispredict.time());
    d("true_block", pp.dbg_true_block.time());
    d("inst_pc (predict1)", pp.dbg_inst_pc.time());
    d("fold_idx[0] (predict1)", pp.dbg_fold_idx.time());
    d("fold_tag[0] (predict1)", pp.dbg_fold_tag.time());
    d("p1_return (predict1)", pp.dbg_p1_return.time());
    d("hist_input (update)", pp.dbg_hist_input.time());
    d("gh[0] (update)", pp.dbg_gh_fanout.time());
    d("fold_read_in_compute[0]", pp.dbg_fold_read_in_compute.time());
    d("fold_compute[0] (update)", pp.dbg_fold_compute.time());
    d("fold_early_write[0]", pp.dbg_fold_early_write.time());
    d("fold_apply[0] (update)", pp.dbg_fold_apply.time());

    // Final output
    uint64_t worst_pred = 0;
    constexpr uint64_t PN = std::remove_reference_t<decltype(pp.pred)>::size;
    for (uint64_t i = 0; i < PN; i++)
        worst_pred = std::max(worst_pred, pp.pred[i].time());
    d("pred[I] (scatter)", worst_pred);

    // Meta pipe
    constexpr uint64_t MP = sizeof(pp.meta_pipe) / sizeof(pp.meta_pipe[0]);
    uint64_t worst_meta = 0;
    for (uint64_t i = 0; i < MP; i++)
        worst_meta = std::max(worst_meta, pp.meta_pipe[i].time());
    d("meta_pipe (worst)", worst_meta);

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
