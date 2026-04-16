// Test10: One-hot block_entry with LINEINST=64, bit-shift line_end
// L0 = binary block_entry, LINEINST=1024, arithmetic line_end (baseline ≈3.12)
// L1 = one-hot block_entry, LINEINST=64, bit-shift line_end
// L2 = binary block_entry, LINEINST=64, arithmetic line_end
// L3 = binary block_entry, LINEINST=16, arithmetic line_end
#pragma once

#ifndef TEST_LEVEL
#define TEST_LEVEL 0
#endif

#include "common.hpp"

using namespace hcm;

template <u64 TABLE_SIZE, u64 TAG_WIDTH>
struct MiniTable {
    static constexpr u64 IDX_BITS = std::bit_width(TABLE_SIZE - 1);
    static constexpr u64 table_size = TABLE_SIZE;
    static constexpr u64 tag_width = TAG_WIDTH;
    hcm::ram<val<TAG_WIDTH>, TABLE_SIZE> tag_ram{"tag"};
    hcm::ram<val<1>, TABLE_SIZE> pred_ram{"pred"};
    rwram<2, TABLE_SIZE, 4> hyst_ram{"hyst"};
    rwram<1, TABLE_SIZE, 4> u_ram{"u"};
};

using Tables = std::tuple<
    MiniTable<2048, 11>, MiniTable<2048, 11>, MiniTable<2048, 11>,
    MiniTable<2048, 10>, MiniTable<2048, 10>,
    MiniTable<2048, 9>,  MiniTable<2048, 9>,  MiniTable<2048, 8>
>;

struct test_t10 : predictor {
    static constexpr u64 N = 8;
    static constexpr u64 LOG_N = 3;
    static constexpr u64 IDX_BITS = 11;
    static constexpr u64 NUM_TABLES = 8;
    static constexpr u64 MAX_TAG_W = 11;
    static constexpr u64 MAX_HTAG_W = MAX_TAG_W - LOG_N;
    static constexpr u64 HYST_W = 2;
    static constexpr u64 U_W = 1;
    static constexpr u64 META_BITS = 4;
    static constexpr u64 META_PIPE = 2;

#if TEST_LEVEL == 0 || TEST_LEVEL == 9
    static constexpr u64 LINEINST = 1024;
    static constexpr u64 LOG_LINEINST = 10;
#elif TEST_LEVEL == 3
    static constexpr u64 LINEINST = 16;
    static constexpr u64 LOG_LINEINST = 4;
#elif TEST_LEVEL == 4 || TEST_LEVEL == 5 || TEST_LEVEL == 6 || TEST_LEVEL == 8
    static constexpr u64 LINEINST = 32;
    static constexpr u64 LOG_LINEINST = 5;
#elif TEST_LEVEL == 7
    static constexpr u64 LINEINST = 64;
    static constexpr u64 LOG_LINEINST = 6;
#elif TEST_LEVEL == 10
    static constexpr u64 LINEINST = 128;
    static constexpr u64 LOG_LINEINST = 7;
#elif TEST_LEVEL == 11
    static constexpr u64 LINEINST = 256;
    static constexpr u64 LOG_LINEINST = 8;
#elif TEST_LEVEL == 12
    static constexpr u64 LINEINST = 512;
    static constexpr u64 LOG_LINEINST = 9;
#else
    static constexpr u64 LINEINST = 64;
    static constexpr u64 LOG_LINEINST = 6;
#endif

    hcm::ram<val<1>, 256> p1_ram[N]{"p1"};
    arr<reg<1>, N> pred;

    geometric_folds<NUM_TABLES, 8, 100, IDX_BITS, MAX_HTAG_W> gfolds;
    Tables tables;

    arr<reg<IDX_BITS>, NUM_TABLES> gindex;
    arr<reg<MAX_TAG_W>, NUM_TABLES> readt;
    arr<reg<1>, NUM_TABLES> readc;
    arr<reg<HYST_W>, NUM_TABLES> readh;
    arr<reg<U_W>, NUM_TABLES> readu;
    reg<NUM_TABLES> notumask;
    arr<reg<MAX_HTAG_W>, NUM_TABLES> htag;

    arr<reg<NUM_TABLES + 1>, N> match_reg;
    arr<reg<NUM_TABLES + 1>, N> match1_reg;
    arr<reg<NUM_TABLES + 1>, N> match2_reg;
    arr<reg<1>, N> pred1_tage;
    arr<reg<1>, N> pred2_tage;
    arr<reg<1>, N> newly_alloc;
    arr<reg<META_BITS, i64>, META_PIPE> meta;

    reg<N> altsel_bits;
    arr<reg<1>, N> branch_dir;
    u64 num_branch = 0;
    u64 block_size = 0;

#if TEST_LEVEL == 1 || TEST_LEVEL == 5 || TEST_LEVEL == 9
    reg<LINEINST> block_entry;  // one-hot
#else
    reg<LOG_LINEINST> block_entry;  // binary (L0: 10-bit, L2: 6-bit)
#endif
    reg<N + 1> rank;
    reg<N> X;

#if TEST_LEVEL == 8
    // ROM: entry[i] bit j = 1 iff (i + j >= LINEINST)
    rom<val<LINEINST>, LINEINST> line_end_rom = [](u64 i) -> u64 {
        if (i == 0) return 0;
        return (~0ULL) << (LINEINST - i);
    };
#endif

    template <u64 I> auto tidx(auto &gi) {
        using Table = std::tuple_element_t<I, Tables>;
        return val<Table::IDX_BITS>{gi};
    }

#if TEST_LEVEL == 1 || TEST_LEVEL == 5 || TEST_LEVEL == 9
    val<1> line_end() { return block_entry >> (LINEINST - 1 - block_size); }
#elif TEST_LEVEL == 6 || TEST_LEVEL == 7
    // L6/L7: rearranged — no adder, just reg >= C++ constant
    val<1> line_end() {
        if (block_size >= LINEINST) return hard<1>{};
        return block_entry >= (LINEINST - block_size);
    }
#elif TEST_LEVEL == 8
    // L8: ROM lookup + free shift (>> with C++ count is free)
    val<1> line_end() {
        if (block_size >= LINEINST) return hard<1>{};
        return line_end_rom(block_entry) >> block_size;
    }
#else
    // L0, L2, L3, L4: binary arithmetic
    val<1> line_end() { return (block_entry + block_size) >= hard<LINEINST>{}; }
#endif
    val<1> last_pred() { return rank >> (N - num_branch); }

    val<1> predict1(val<64> inst_pc) override {
        inst_pc.fanout(hard<2>{});

#if TEST_LEVEL == 1 || TEST_LEVEL == 5 || TEST_LEVEL == 9
        block_entry = val<LOG_LINEINST>{inst_pc >> 2}.decode().concat();
        block_entry.fanout(hard<LINEINST>{});
#else
        // L0, L2, L3, L4: binary
        block_entry = val<LOG_LINEINST>{inst_pc >> 2};
        block_entry.fanout(hard<N + 2>{});
#endif
        rank = val<N + 1>{1};
        rank.fanout(hard<N + 2>{});
        X = val<LOG_N>{inst_pc >> 2}.decode().concat();
        X.fanout(hard<2>{});

        val<8> idx = inst_pc >> 5;
        idx.fanout(hard<N>{});
        for (u64 i = 0; i < N; i++) {
            pred[i] = p1_ram[i].read(idx);
        }
        pred.fanout(hard<2>{});
        block_size = 1;
        num_branch = 0;

        reuse_prediction(~line_end());
        return pred[0];
    }

    val<1> reuse_predict1([[maybe_unused]] val<64> inst_pc) override {
        block_size++;
        reuse_prediction(~line_end());
        return hard<0>{};
    }

    val<1> predict2(val<64> inst_pc) override {
        val<IDX_BITS> lineaddr = inst_pc >> 2;
        lineaddr.fanout(hard<1 + NUM_TABLES * 2>{});
        gfolds.fanout(hard<2>{});

        for (u64 i = 0; i < NUM_TABLES; i++) {
            gindex[i] = lineaddr ^ gfolds.template get<0>(i);
        }
        gindex.fanout(hard<4>{});
        for (u64 i = 0; i < NUM_TABLES; i++) {
            htag[i] = val<MAX_HTAG_W>{lineaddr}.reverse() ^ gfolds.template get<1>(i);
        }
        htag.fanout(hard<2>{});

        static_loop<NUM_TABLES>([&]<u64 I>() {
            auto &t = std::get<I>(tables);
            readt[I] = t.tag_ram.read(tidx<I>(gindex[I]));
            readc[I] = t.pred_ram.read(tidx<I>(gindex[I]));
            readh[I] = t.hyst_ram.read(tidx<I>(gindex[I]));
            readu[I] = t.u_ram.read(tidx<I>(gindex[I]));
        });
        readt.fanout(hard<N + 1>{});
        readc.fanout(hard<3>{});
        readh.fanout(hard<2>{});
        readu.fanout(hard<2>{});
        notumask = ~readu.concat();
        notumask.fanout(hard<2>{});

        arr<val<1>, NUM_TABLES> htagcmp_split = [&](int i) -> val<1> {
            return val<MAX_HTAG_W>{readt[i]} == val<MAX_HTAG_W>{htag[i]};
        };
        val<NUM_TABLES> htagcmp = htagcmp_split.fo1().concat();
        htagcmp.fanout(hard<N>{});

        val<NUM_TABLES> gpreds = readc.concat();
        gpreds.fanout(hard<N>{});

        pred.fanout(hard<N>{});
        arr<val<NUM_TABLES + 1>, N> preds = [&](u64 r) -> val<NUM_TABLES + 1> {
            return concat(val<1>{pred[r]}, gpreds);
        };
        preds.fanout(hard<2 * N>{});

        static_loop<N>([&]<u64 R>() {
            arr<val<1>, NUM_TABLES> tagcmp = [&](int i) {
                return val<LOG_N>{readt[i] >> MAX_HTAG_W} == hard<R>{};
            };
            match_reg[R] = concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp);
        });
        match_reg.fanout(hard<2>{});

        for (u64 r = 0; r < N; r++) {
            match1_reg[r] = match_reg[r].one_hot();
        }
        match1_reg.fanout(hard<3>{});
        for (u64 r = 0; r < N; r++) {
            pred1_tage[r] = (match1_reg[r] & preds[r]) != hard<0>{};
        }
        pred1_tage.fanout(hard<2>{});

        for (u64 r = 0; r < N; r++) {
            match2_reg[r] = (match_reg[r] ^ match1_reg[r]).one_hot();
        }
        match2_reg.fanout(hard<2>{});
        for (u64 r = 0; r < N; r++) {
            pred2_tage[r] = (match2_reg[r] & preds[r]) != hard<0>{};
        }
        pred2_tage.fanout(hard<2>{});

        // Meta/altsel
        meta.fanout(hard<2>{});
        arr<val<1>, NUM_TABLES> weakctr = [&](int i) -> val<1> {
            return readh[i] == hard<0>{};
        };
        val<NUM_TABLES> coldctr = notumask & weakctr.fo1().concat();
        coldctr.fanout(hard<N>{});
        val<1> metasign = (meta[META_PIPE - 1] >= hard<0>{});
        metasign.fanout(hard<N>{});
        for (u64 r = 0; r < N; r++) {
            newly_alloc[r] = (match1_reg[r] & coldctr) != hard<0>{};
        }
        newly_alloc.fanout(hard<2>{});

        // Write altsel_bits reg
        altsel_bits = arr<val<1>, N>{[&](u64 r) -> val<1> {
            return metasign & newly_alloc[r] & (match2_reg[r] != hard<0>{});
        }}.concat();

        // Return — no reuse_prediction in predict2
        val<1> alt = metasign & newly_alloc[num_branch] &
                     (match2_reg[num_branch] != hard<0>{});
        val<1> taken = select(alt.fo1(), pred2_tage[num_branch], pred1_tage[num_branch]);
        taken.fanout(hard<2>{});
        return taken;
    }

    val<1> reuse_predict2(val<64>) override {
        val<1> alt = altsel_bits >> num_branch;
        val<1> taken = select(alt.fo1(), pred2_tage[num_branch],
                              pred1_tage[num_branch]);
        taken.fanout(hard<2>{});
        // No reuse_prediction here — predict1/reuse_predict1 handles it
        return taken;
    }

    void update_condbr(val<64>, val<1> taken, val<64>) override {
        branch_dir[num_branch] = taken.fo1();
        num_branch++;
        reuse_prediction(~(line_end() | last_pred()));
    }

    void update_cycle(instruction_info &) override {
        gfolds.update(val<6>{hard<0>{}});
    }
};

using branch_predictor = test_t10;
