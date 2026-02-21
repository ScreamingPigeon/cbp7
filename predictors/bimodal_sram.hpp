#include "../cbp.hpp"
#include "../harcom.hpp"

using namespace hcm;

//sram additions
template<u64 W, u64 ENTRIES, u64 MAXN=256, u64 MAXM=128>
struct sram_adapter {
    ram<val<W>, ENTRIES> mem; // functional storage used by simulator
    using model = sram_banked<ENTRIES, W, MAXN, MAXM>;

    template<typename I>
    val<W> read(I index) {
        // register storage with HARCOM panel for area accounting
        panel.update_storage(model::num_bits(), true);
        return mem.read(index);
    }

    template<typename I>
    void write(I index, val<W> v) {
        panel.update_storage(model::num_bits(), true);
        mem.write(index, v);
    }
};


template<u64 LOGLB=6, u64 LOGP2=17, u64 LOGP1=12, u64 P1_MAXN=256, u64 P1_MAXM=128, u64 P2_MAXN=256, u64 P2_MAXM=128>
struct bimodal_sram : predictor {
    static_assert(LOGLB>2);
    static constexpr u64 LOGLINEINST = LOGLB-2;
    static constexpr u64 LINEINST = 1<<LOGLINEINST;
    static_assert(LOGP1 > LOGLINEINST);
    static_assert(LOGP2 > LOGLINEINST);
    static constexpr u64 index1_bits = LOGP1-LOGLINEINST;
    static constexpr u64 index2_bits = LOGP2-LOGLINEINST;

    reg<index1_bits> index1;
    arr<reg<1>,LINEINST> ctr1_hi; // prediction bits read from P1 table
    reg<LINEINST> p1; // predictions
    reg<index2_bits> index2;
    arr<reg<1>,LINEINST> ctr2_hi; // prediction bits read from P2 table
    reg<LINEINST> p2; // predictions

    // simulation artifacts
    u64 num_branch = 0;
    u64 block_size = 0;
    arr<reg<LOGLINEINST>,LINEINST> branch_offset;
    arr<reg<1>,LINEINST> branch_dir;
    reg<LINEINST> block_entry;

    // Functional RAMs wrapped by sram_adapter for accounting
    sram_adapter<1,(1ull<<index1_bits),P1_MAXN,P1_MAXM> table1_hi[LINEINST]; // P1 prediction bits
    sram_adapter<1,(1ull<<index2_bits),P2_MAXN,P2_MAXM> table2_hi[LINEINST]; // P2 prediction bits

    zone UPDATE_ONLY;
    // hysteresis bits (0=weak, 1=strong)
    sram_adapter<1,(1ull<<index1_bits),P1_MAXN,P1_MAXM> table1_lo[LINEINST];
    sram_adapter<1,(1ull<<index2_bits),P2_MAXN,P2_MAXM> table2_lo[LINEINST];

    void new_block(val<64> inst_pc)
    {
        val<LOGLINEINST> offset = inst_pc.fo1() >> 2;
        block_entry = offset.fo1().decode().concat();
        block_entry.fanout(hard<LINEINST*2>{});
        block_size = 1;
    }

    val<1> predict1([[maybe_unused]] val<64> inst_pc)
    {
        inst_pc.fanout(hard<2>{});
        new_block(inst_pc);
        index1 = inst_pc >> LOGLB;
        index1.fanout(hard<LINEINST>{});
        for (u64 i=0; i<LINEINST; i++) {
            ctr1_hi[i] = table1_hi[i].read(index1);
        }
        ctr1_hi.fanout(hard<2>{});
        p1 = ctr1_hi.concat();
        p1.fanout(hard<LINEINST>{});
        return (block_entry & p1) != hard<0>{};
    };

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc)
    {
        return ((block_entry<<block_size) & p1) != hard<0>{};
    };

    val<1> predict2(val<64> inst_pc)
    {
        inst_pc.fanout(hard<2>{});
        index2 = inst_pc >> LOGLB;
        index2.fanout(hard<LINEINST>{});
        for (u64 i=0; i<LINEINST; i++) {
            ctr2_hi[i] = table2_hi[i].read(index2);
        }
        ctr2_hi.fanout(hard<2>{});
        p2 = ctr2_hi.concat();
        p2.fanout(hard<LINEINST>{});
        val<1> taken = (block_entry & p2) != hard<0>{};
        taken.fanout(hard<2>{});
        reuse_prediction(~val<1>{block_entry>>(LINEINST-1)});
        return taken;
    }

    val<1> reuse_predict2([[maybe_unused]] val<64> inst_pc)
    {
        val<1> taken = ((block_entry<<block_size) & p2) != hard<0>{};
        taken.fanout(hard<2>{});
        reuse_prediction(~val<1>{block_entry>>(LINEINST-1-block_size)});
        block_size++;
        return taken;
    }

    void update_condbr(val<64> branch_pc, val<1> taken, [[maybe_unused]] val<64> next_pc)
    {
        assert(num_branch < LINEINST);
        branch_offset[num_branch] = branch_pc.fo1() >> 2;
        branch_dir[num_branch] = taken.fo1();
        num_branch++;
    }

    void update_cycle([[maybe_unused]] instruction_info &block_end_info)
    {
        if (num_branch == 0) {
            return;
        }
        val<1> &mispredict = block_end_info.is_mispredict;
        mispredict.fanout(hard<LINEINST+1>{});
        branch_offset.fanout(hard<LINEINST>{});
        branch_dir.fanout(hard<2>{});
        index1.fanout(hard<LINEINST*3>{});
        index2.fanout(hard<LINEINST*3>{});
        u64 update_valid = (u64(1)<<num_branch)-1;
        arr<val<LINEINST>,LINEINST> update_mask = [&](u64 offset){
            arr<val<1>,LINEINST> match_offset = [&](u64 i){return branch_offset[i] == offset;};
            return match_offset.fo1().concat() & update_valid;
        };
        update_mask.fanout(hard<2>{});
        arr<val<1>,LINEINST> is_branch = [&](u64 offset){
            return update_mask[offset] != hard<0>{};
        };
        val<LINEINST> branch_mask = is_branch.fo1().concat();
        branch_mask.fanout(hard<3>{});
        val<LINEINST> actualdirs = branch_dir.concat();
        actualdirs.fanout(hard<LINEINST>{});
        arr<val<1>,LINEINST> branch_taken = [&](u64 offset){
            return (actualdirs & update_mask[offset]) != hard<0>{};
        };

        arr<val<1>,LINEINST> last_branch = branch_offset[num_branch-1].decode();
        arr<val<1>,LINEINST> mispredicted = [&](u64 i){
            return last_branch[i].fo1() & mispredict;
        };
        mispredicted.fanout(hard<2>{});

        val<LINEINST> disagree_mask = (p1 ^ p2) & branch_mask;
        disagree_mask.fanout(hard<2>{});
        arr<val<1>,LINEINST> disagree = disagree_mask.make_array(val<1>{});
        disagree.fanout(hard<2>{});

        arr<val<1>,LINEINST> p1_weak = [&] (u64 i) -> val<1> {
            return execute_if(disagree[i], [&](){
                return ~table1_lo[i].read(index1); // 0=weak
            });
        };

        arr<val<1>,LINEINST> p2_weak = [&] (u64 i) -> val<1> {
            return execute_if(mispredicted[i], [&](){
                return ~table2_lo[i].read(index2); // 0=weak
            });
        };

        need_extra_cycle(mispredict | (disagree_mask != hard<0>{}));

        execute_if(p1_weak.fo1().concat(), [&](u64 i){
            table1_hi[i].write(index1,ctr2_hi[i]);
        });
        execute_if(branch_mask,[&](u64 i){
            table1_lo[i].write(index1,~disagree[i]);
        });

        execute_if(p2_weak.fo1().concat(), [&](u64 i){
            table2_hi[i].write(index2,branch_taken[i].fo1());
        });
        execute_if(branch_mask,[&](u64 i){
            table2_lo[i].write(index2,~mispredicted[i]);
        });

        num_branch = 0;
    }
};
