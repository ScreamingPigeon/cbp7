// bimodalB — RABT-style bimodal fallback, standalone.
//
// Indexing: single shared table, indexed by instruction-aligned PC.
//   - 1 pred table of (1 << LOGB) entries, each entry is N-bit (one bit per slot)
//   - 1 hyst table of (1 << LOGB_HYST) entries, also N-bit, HALF capacity
//     (matches RABT's FB_BH_CAPACITY = FB_CAPACITY / 2)
//   - Read entry at block start; pre-compute pred[i] for all i.
//   - Subsequent predict2/reuse_predict* just look up pred[num_branch].
//   - Bit selection by branch ordinal G (pure ordinal, no XOR).
//
// Mirrors RABT M2's fb_ctr / fb_bim_hyst (TageAheadHC_IR.hpp:247,315).
// Patterned on predictors/bimodalN.hpp (RAM read once, pre-cached per slot).

#include "../cbp.hpp"
#include "../harcom.hpp"

using namespace hcm;


template<u64 LOGB = 13, u64 N = 7>
struct bimodalB : predictor {
    static_assert(N > 0 && N <= 64);
    static_assert(LOGB >= 1);

    static constexpr u64 BBMAX     = 64;       // matches bimodalN convention
    static constexpr u64 LOGB_HYST = LOGB - 1; // half-capacity hyst (FB_BH_CAPACITY)

    reg<LOGB>      index;
    reg<LOGB_HYST> hyst_index;
    arr<reg<1>, N> read_pred;   // raw bits from pred RAM
    arr<reg<1>, N> pred;        // per-slot pre-computed predictions
    arr<reg<1>, N> bh_bits;     // hyst bits read in update_cycle

    u64 num_branch = 0;
    u64 bb_inst = 0;
    arr<reg<1>, N> branch_dir;

    ram<arr<val<1>, N>, (1 << LOGB)>      ctr_hi{"bimB_pred"};

    zone UPDATE_ONLY;
    ram<arr<val<1>, N>, (1 << LOGB_HYST)> ctr_lo{"bimB_hyst"};

    val<1> predict1([[maybe_unused]] val<64> inst_pc) {
        assert(num_branch == 0);
        assert(bb_inst == 0);
        inst_pc.fanout(hard<2>{});
        index      = inst_pc >> 2;                     // RABT-style index
        hyst_index = inst_pc >> 2;                     // low LOGB_HYST bits
        read_pred  = ctr_hi.read(index);               // ONE pred RAM read
        for (u64 i = 0; i < N; i++) pred[i] = read_pred[i];  // pure ordinal
        pred.fanout(hard<2 * BBMAX>{});
        bb_inst++;
        reuse_prediction(bb_inst < BBMAX);
        return pred[num_branch];
    }

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) {
        bb_inst++;
        reuse_prediction(bb_inst < BBMAX);
        return pred[num_branch];
    }

    val<1> predict2([[maybe_unused]] val<64> inst_pc) {
        return pred[num_branch];
    }

    val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) {
        return pred[num_branch];
    }

    void update_condbr([[maybe_unused]] val<64> branch_pc, val<1> taken,
                       [[maybe_unused]] val<64> next_pc) {
        assert(num_branch < N);
        branch_dir[num_branch] = taken.fo1();
        num_branch++;
        bb_inst = 0;
        if (num_branch == N) reuse_prediction(0);
    }

    void update_cycle(instruction_info &block_end_info) {
        val<1> &mispredict = block_end_info.is_mispredict;
        if (num_branch == 0) { bb_inst = 0; return; }

        index.fanout(hard<2>{});                   // 1 write + 1 fanout? actually 1 use (write)
        hyst_index.fanout(hard<2>{});              // 1 read + 1 write
        branch_dir.fanout(hard<N>{});              // each used in bundle and ~wrong
        mispredict.fanout(hard<3>{});              // need_extra_cycle + 2 execute_if

        // wrong[i] = (slot i was used) AND (pred[i] != actual)
        // for slot i not used, wrong[i] = 0 (treat as "correct")
        // → ~wrong[i] = 1 (strong) for unused or correct slots, 0 (weak) for wrong slot
        arr<val<1>, N> wrong = [&](u64 i) -> val<1> {
            return val<1>{i < num_branch} & (read_pred[i] ^ branch_dir[i]);
        };
        wrong.fanout(hard<2>{});                   // used in bundle and new_bh

        // === Cycle 1: hyst READ (unconditional — execute_if can't wrap arr), before extra-cycle ===
        bh_bits = ctr_lo.read(hyst_index);
        bh_bits.fanout(hard<N>{});                  // used per-slot in bundle

        // === Cycle boundary ===
        need_extra_cycle(mispredict);

        // === Cycle 2: WRITES, after extra-cycle boundary ===
        // pred update: flip slot i if wrong[i] AND hyst was weak (~bh_bits[i])
        execute_if(mispredict, [&]() {
            arr<val<1>, N> bundle = [&](u64 i) {
                val<1> flip = wrong[i] & ~bh_bits[i];
                return select(flip, ~read_pred[i], read_pred[i]);
            };
            ctr_hi.write(index, bundle.fo1());
        });
        // hyst update: write ~wrong for all N slots
        //   correct/unused (wrong=0) → strong (1)
        //   wrong (wrong=1) → weak (0)
        execute_if(mispredict, [&]() {
            arr<val<1>, N> new_bh = [&](u64 i) { return ~wrong[i]; };
            ctr_lo.write(hyst_index, new_bh.fo1());
        });

        num_branch = 0;
        bb_inst = 0;
    }
};
