// bimodalA — tage<>-style bimodal fallback, standalone (single level).
//
// Indexing: per-line-offset tables, indexed by block-aligned PC.
//   - LINEINST separate 1-bit RAMs, each (1<<index_bits) entries
//   - Read all LINEINST predictions at block start
//   - Bit selected by branch's fetch-line offset
//
// Mirrors tage<>'s `bim[offset]` / `bhyst[offset]` (predictors/tage.hpp:97-111).
// Update logic patterned on bimodal.hpp's P2 path (proven correct).

#include "../cbp.hpp"
#include "../harcom.hpp"

using namespace hcm;


template<u64 LOGLB = 6, u64 LOGB = 12>
struct bimodalA : predictor {
    static_assert(LOGLB > 2);
    static_assert(LOGB > LOGLB - 2);

    static constexpr u64 LOGLINEINST = LOGLB - 2;
    static constexpr u64 LINEINST    = 1 << LOGLINEINST;
    static constexpr u64 index_bits  = LOGB - LOGLINEINST;

    reg<index_bits> index;
    arr<reg<1>, LINEINST> ctr_hi;
    reg<LINEINST> preds;

    u64 num_branch = 0;
    u64 block_size = 0;
    arr<reg<LOGLINEINST>, LINEINST> branch_offset;
    arr<reg<1>, LINEINST> branch_dir;
    reg<LINEINST> block_entry;

    ram<val<1>, (1 << index_bits)> table_hi[LINEINST]{"bimA_pred"};

    zone UPDATE_ONLY;
    ram<val<1>, (1 << index_bits)> table_lo[LINEINST]{"bimA_hyst"};

    void new_block(val<64> inst_pc) {
        val<LOGLINEINST> offset = inst_pc.fo1() >> 2;
        block_entry = offset.fo1().decode().concat();
        block_entry.fanout(hard<LINEINST * 2>{});
        block_size = 1;
    }

    val<1> predict1(val<64> inst_pc) {
        inst_pc.fanout(hard<2>{});
        new_block(inst_pc);
        index = inst_pc >> LOGLB;
        index.fanout(hard<LINEINST>{});
        for (u64 i = 0; i < LINEINST; i++)
            ctr_hi[i] = table_hi[i].read(index);
        ctr_hi.fanout(hard<2>{});
        preds = ctr_hi.concat();
        preds.fanout(hard<LINEINST>{});
        return (block_entry & preds) != hard<0>{};
    }

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) {
        return ((block_entry << block_size) & preds) != hard<0>{};
    }

    val<1> predict2([[maybe_unused]] val<64> inst_pc) {
        val<1> taken = (block_entry & preds) != hard<0>{};
        taken.fanout(hard<2>{});
        reuse_prediction(~val<1>{block_entry >> (LINEINST - 1)});
        return taken;
    }

    val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc) {
        val<1> taken = ((block_entry << block_size) & preds) != hard<0>{};
        taken.fanout(hard<2>{});
        reuse_prediction(~val<1>{block_entry >> (LINEINST - 1 - block_size)});
        block_size++;
        return taken;
    }

    void update_condbr(val<64> branch_pc, val<1> taken,
                       [[maybe_unused]] val<64> next_pc) {
        assert(num_branch < LINEINST);
        branch_offset[num_branch] = branch_pc.fo1() >> 2;
        branch_dir[num_branch]    = taken.fo1();
        num_branch++;
    }

    void update_cycle([[maybe_unused]] instruction_info &block_end_info) {
        if (num_branch == 0) return;
        val<1> &mispredict = block_end_info.is_mispredict;
        mispredict.fanout(hard<LINEINST + 1>{});
        branch_offset.fanout(hard<LINEINST>{});
        branch_dir.fanout(hard<2>{});
        index.fanout(hard<LINEINST * 3>{});

        u64 update_valid = (u64(1) << num_branch) - 1;
        arr<val<LINEINST>, LINEINST> update_mask = [&](u64 offset) {
            arr<val<1>, LINEINST> match_offset = [&](u64 i) {
                return branch_offset[i] == offset;
            };
            return match_offset.fo1().concat() & update_valid;
        };
        update_mask.fanout(hard<2>{});
        arr<val<1>, LINEINST> is_branch = [&](u64 offset) {
            return update_mask[offset] != hard<0>{};
        };
        is_branch.fanout(hard<2>{});

        val<LINEINST> actualdirs = branch_dir.concat();
        actualdirs.fanout(hard<LINEINST>{});
        arr<val<1>, LINEINST> branch_taken = [&](u64 offset) {
            return (actualdirs & update_mask[offset]) != hard<0>{};
        };

        arr<val<1>, LINEINST> last_branch = branch_offset[num_branch - 1].decode();
        arr<val<1>, LINEINST> mispredicted = [&](u64 i) {
            return last_branch[i].fo1() & mispredict;
        };
        mispredicted.fanout(hard<2>{});

        // === Cycle 1: all hyst READS, before the extra-cycle boundary ===
        arr<val<1>, LINEINST> weak = [&](u64 i) -> val<1> {
            return execute_if(mispredicted[i], [&]() {
                return ~table_lo[i].read(index);
            });
        };

        // === Cycle boundary: extra cycle on mispredict ===
        need_extra_cycle(mispredict);

        // === Cycle 2: all WRITES, after the extra-cycle boundary ===
        for (u64 i = 0; i < LINEINST; i++) {
            execute_if(weak[i].fo1(), [&]() {
                table_hi[i].write(index, branch_taken[i].fo1());
            });
        }
        for (u64 i = 0; i < LINEINST; i++) {
            execute_if(is_branch[i], [&]() {
                table_lo[i].write(index, ~mispredicted[i]);
            });
        }

        num_branch = 0;
        block_size = 0;
    }
};
