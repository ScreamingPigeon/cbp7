// TageTable parameter sweep: latency and area across configurations.
// Organized so each group answers one design question.
// FF/SRAM u-bit variants are paired for direct comparison.

#include <algorithm>
#include "../predictors/custom/TageTable.hpp"

using namespace hcm;

class harcom_superuser {
public:
  void next_cycle() { panel.next_cycle(); }
  auto get(valtype auto x) { return x.get(); }
  template <typename V> u64 time(V &x) { return x.timing; }
} hsu;

// ============================================================================
// Config declarations — each gets its own HARCOM region for isolated costs
// ============================================================================

#define DECL(id, SZ, HIST, TW, CW, UW, Nbr, NB, AH, ST, SU, UF, DC, CACHE)  \
  static region region_##id;                                                   \
  static TageTable<SZ, HIST, TW, CW, UW, Nbr, NB, AH, ST, SU, UF, DC,       \
                   DefaultResetFn, CACHE>                                      \
      sw_##id;

// Baseline: SIZE=256, TAG=11, CTR=3, U=1, N=4, BANKS=1, shared tag/u, no
// ahead, no cache.  Every group varies one dimension from this baseline.

// --- 1. Table size scaling (FF u-bits vs SRAM u-bits) ---
//        Vary SIZE.  Everything else = baseline.
DECL(s64f,   64,  32, 11, 3, 1, 4, 1, false, true, true, true,  1024, false)
DECL(s64s,   64,  32, 11, 3, 1, 4, 1, false, true, true, false,  512, false)
DECL(s128f, 128,  32, 11, 3, 1, 4, 1, false, true, true, true,  1024, false)
DECL(s128s, 128,  32, 11, 3, 1, 4, 1, false, true, true, false,  512, false)
DECL(s256f, 256,  32, 11, 3, 1, 4, 1, false, true, true, true,  1024, false)
DECL(s256s, 256,  32, 11, 3, 1, 4, 1, false, true, true, false,  512, false)
DECL(s512f, 512,  32, 11, 3, 1, 4, 1, false, true, true, true,  1024, false)
DECL(s512s, 512,  32, 11, 3, 1, 4, 1, false, true, true, false,  512, false)
DECL(s1Kf, 1024,  32, 11, 3, 1, 4, 1, false, true, true, true,  1024, false)
DECL(s1Ks, 1024,  32, 11, 3, 1, 4, 1, false, true, true, false,  512, false)
DECL(s2Kf, 2048,  32, 11, 3, 1, 4, 1, false, true, true, true,  1024, false)
DECL(s2Ks, 2048,  32, 11, 3, 1, 4, 1, false, true, true, false,  512, false)

// --- 2. Banking: vary N and NUM_BANKS (SIZE=256, SRAM u-bits) ---
//        Shows cost of increasing branch throughput.
DECL(n1_1x1, 256, 32, 11, 3, 1, 1, 1, false, true, true, false, 512, false)
DECL(n2_2x1, 256, 32, 11, 3, 1, 2, 2, false, true, true, false, 512, false)
DECL(n4_2x2, 256, 32, 11, 3, 1, 4, 2, false, true, true, false, 512, false)
DECL(n4_4x1, 256, 32, 11, 3, 1, 4, 4, false, true, true, false, 512, false)
DECL(n8_2x4, 256, 32, 11, 3, 1, 8, 2, false, true, true, false, 512, false)
DECL(n8_4x2, 256, 32, 11, 3, 1, 8, 4, false, true, true, false, 512, false)

// --- 3. Entry width: vary TAG and CTR (SIZE=256, N=4, 1 bank, SRAM u) ---
//        Shows cost of wider tags and counters.
DECL(t7c2,  256, 32,  7, 2, 1, 4, 1, false, true, true, false, 512, false)
DECL(t7c3,  256, 32,  7, 3, 1, 4, 1, false, true, true, false, 512, false)
DECL(t11c2, 256, 32, 11, 2, 1, 4, 1, false, true, true, false, 512, false)
DECL(t11c3, 256, 32, 11, 3, 1, 4, 1, false, true, true, false, 512, false)
DECL(t11c4, 256, 32, 11, 4, 1, 4, 1, false, true, true, false, 512, false)
DECL(t13c3, 256, 32, 13, 3, 1, 4, 1, false, true, true, false, 512, false)

// --- 4. U-bit width: U=1 vs U=2, FF vs SRAM (SIZE=256, N=4) ---
DECL(u1ff, 256, 32, 11, 3, 1, 4, 1, false, true, true, true,  1024, false)
DECL(u1sr, 256, 32, 11, 3, 1, 4, 1, false, true, true, false,  512, false)
DECL(u2ff, 256, 32, 11, 3, 2, 4, 1, false, true, true, true,  1024, false)
DECL(u2sr, 256, 32, 11, 3, 2, 4, 1, false, true, true, false,  512, false)

// --- 5. Features: ahead / per-bank tags / FF cache (SIZE=256, SRAM u) ---
//        Each adds one feature to the N=4, 2-bank baseline.
DECL(base,   256, 32, 11, 3, 1, 4, 2, false, true,  true,  false, 512, false)
DECL(ahead,  256, 32, 11, 3, 1, 4, 2, true,  true,  true,  false, 512, false)
DECL(pbank,  256, 32, 11, 3, 1, 4, 4, false, false, false, false, 512, false)
//        FF cache requires BPB>1, so use N=8, 2 banks (BPB=4).
DECL(cbase,  256, 32, 11, 3, 1, 8, 2, false, true,  true,  false, 512, false)
DECL(cache,  256, 32, 11, 3, 1, 8, 2, false, true,  true,  false, 512, true)
DECL(ah_ca,  256, 32, 11, 3, 1, 8, 2, true,  true,  true,  false, 512, true)
DECL(allsr,  256, 32, 11, 3, 1, 8, 2, true,  false, false, false, 512, true)

// --- 6. Realistic large configs (1K and 2K entries) ---
DECL(l1K,    1024, 64, 11, 3, 1, 8, 4, false, true,  true,  false, 512, false)
DECL(l1Kah,  1024, 64, 11, 3, 1, 8, 4, true,  true,  true,  false, 512, false)
DECL(l1Kall, 1024, 64, 11, 3, 1, 8, 4, true,  false, false, false, 512, true)
DECL(l2K,    2048, 32, 11, 3, 1, 4, 1, false, true,  true,  false, 512, false)
DECL(l2Kn1,  2048, 32, 11, 3, 1, 1, 1, false, true,  true,  false, 512, false)

// ============================================================================
// Measurement macros
// ============================================================================

#define HDR()                                                                  \
  printf("%-12s %5s %3s %3s %1s %2s %-5s %1s %1s %1s %1s "                    \
         "%7s %9s %7s %6s  %5s %5s %5s %5s %5s %5s\n",                        \
         "", "", "", "", "", "", "", "", "", "", "",                            \
         "", "", "mm\xC2\xB2", "mW",                                           \
         "ps", "ps", "ps", "ps", "ps", "ps");                                  \
  printf("%-12s %5s %3s %3s %1s %2s %-5s %1s %1s %1s %1s "                    \
         "%7s %9s %7s %6s  %5s %5s %5s %5s %5s %5s\n",                        \
         "Config", "SIZE", "TAG", "CTR", "U", " N", "BxBPB",                   \
         "A", "S", "F", "C",                                                   \
         "bits", "trans", "area", "staPwr",                                    \
         "Read", "Hit", "Tag", "Ctr", "Ubit", "Reuse");                       \
  printf("%-12s %5s %3s %3s %1s %2s %-5s %1s %1s %1s %1s "                    \
         "%7s %9s %7s %6s  %5s %5s %5s %5s %5s %5s\n",                        \
         "------------", "-----", "---", "---", "-", "--", "-----",            \
         "-", "-", "-", "-",                                                   \
         "-------", "---------", "-------", "------",                          \
         "-----", "-----", "-----", "-----", "-----", "-----");

#define MEAS(id, SZ, TW, CW, UW, Nbr, NB, AH, ST, UF, CACHE, NAME)           \
  {                                                                            \
    constexpr u64 IDX = clog2(SZ);                                             \
    constexpr u64 BPB = (Nbr) / (NB);                                         \
    constexpr u64 CBITS = (CW) * BPB;                                         \
    val<IDX> idx{1};                                                           \
    val<TW> tag{0x55};                                                         \
    val<CBITS> ctr{0};                                                         \
    val<UW> u{1};                                                              \
    sw_##id.write(idx, 0, 0, tag, ctr, u);                                     \
    hsu.next_cycle();                                                          \
    sw_##id.read(idx, tag, 0);                                                 \
    auto h = sw_##id.getHit(0);                                                \
    auto tg = sw_##id.getTag(0);                                               \
    auto c = sw_##id.getCounter(0, 0);                                         \
    auto uv = sw_##id.getU(0);                                                 \
    u64 ht = hsu.time(h);                                                      \
    u64 tt = hsu.time(tg);                                                     \
    u64 ct = hsu.time(c);                                                      \
    u64 ut = hsu.time(uv);                                                     \
    u64 rd = std::max({ht, tt, ct, ut});                                       \
    u64 bits = panel.storage(region_##id);                                     \
    u64 trans = panel.transistors(region_##id);                                \
    f64 area = panel.area_sram_mm2(region_##id);                               \
    f64 spwr = panel.sta_power_mW(region_##id);                                \
    printf("%-12s %5d %3d %3d %1d %2d %1dx%-3d %1d %1d %1d %1d "              \
           "%7lu %9lu %7.5f %6.4f  %5lu %5lu %5lu %5lu %5lu %5s\n",           \
           NAME, static_cast<int>(SZ), static_cast<int>(TW),                   \
           static_cast<int>(CW), static_cast<int>(UW),                         \
           static_cast<int>(Nbr), static_cast<int>(NB),                        \
           static_cast<int>(BPB), static_cast<int>(AH),                        \
           static_cast<int>(ST), static_cast<int>(UF),                         \
           static_cast<int>(CACHE), bits, trans, area, spwr,                   \
           rd, ht, tt, ct, ut, "-");                                           \
    hsu.next_cycle();                                                          \
  }

#define MEAS_CACHE(id, SZ, TW, CW, UW, Nbr, NB, AH, ST, UF, CACHE, NAME)     \
  {                                                                            \
    constexpr u64 IDX = clog2(SZ);                                             \
    constexpr u64 BPB = (Nbr) / (NB);                                         \
    constexpr u64 CBITS = (CW) * BPB;                                         \
    constexpr u64 SLOT_BITS = clog2(BPB);                                     \
    val<IDX> idx{1};                                                           \
    val<TW> tag{0x55};                                                         \
    val<CBITS> ctr{0};                                                         \
    val<UW> u{1};                                                              \
    sw_##id.write(idx, 0, 0, tag, ctr, u);                                     \
    hsu.next_cycle();                                                          \
    sw_##id.read(idx, tag, 0);                                                 \
    auto h = sw_##id.getHit(0);                                                \
    auto tg = sw_##id.getTag(0);                                               \
    auto c = sw_##id.getCounter(0, 0);                                         \
    auto uv = sw_##id.getU(0);                                                 \
    auto rv = sw_##id.reuseRead(0, val<SLOT_BITS>{0});                         \
    u64 ht = hsu.time(h);                                                      \
    u64 tt = hsu.time(tg);                                                     \
    u64 ct = hsu.time(c);                                                      \
    u64 ut = hsu.time(uv);                                                     \
    u64 rt = hsu.time(rv);                                                     \
    u64 rd = std::max({ht, tt, ct, ut});                                       \
    u64 bits = panel.storage(region_##id);                                     \
    u64 trans = panel.transistors(region_##id);                                \
    f64 area = panel.area_sram_mm2(region_##id);                               \
    f64 spwr = panel.sta_power_mW(region_##id);                                \
    printf("%-12s %5d %3d %3d %1d %2d %1dx%-3d %1d %1d %1d %1d "              \
           "%7lu %9lu %7.5f %6.4f  %5lu %5lu %5lu %5lu %5lu %5lu\n",          \
           NAME, static_cast<int>(SZ), static_cast<int>(TW),                   \
           static_cast<int>(CW), static_cast<int>(UW),                         \
           static_cast<int>(Nbr), static_cast<int>(NB),                        \
           static_cast<int>(BPB), static_cast<int>(AH),                        \
           static_cast<int>(ST), static_cast<int>(UF),                         \
           static_cast<int>(CACHE), bits, trans, area, spwr,                   \
           rd, ht, tt, ct, ut, rt);                                            \
    hsu.next_cycle();                                                          \
  }

#define SECTION(title)                                                         \
  printf("\n--- %s ---\n", title);

// ============================================================================

int main() {
  panel.make_floorplan();

  HDR()

  // ---- 1. Size scaling: FF u-bits vs SRAM u-bits (paired) ----
  // Baseline params held constant: TAG=11, CTR=3, U=1, N=4, 1 bank, shared.
  // Pairs show the cost of the FF u-bit mux tree as SIZE grows.
  SECTION("1. Size scaling: FF u-bits (F=1) vs SRAM u-bits (F=0)")
  MEAS(s64f,   64,  11, 3, 1, 4, 1, 0, 1, 1, 0, "64  u=FF")
  MEAS(s64s,   64,  11, 3, 1, 4, 1, 0, 1, 0, 0, "64  u=SRAM")
  MEAS(s128f, 128,  11, 3, 1, 4, 1, 0, 1, 1, 0, "128 u=FF")
  MEAS(s128s, 128,  11, 3, 1, 4, 1, 0, 1, 0, 0, "128 u=SRAM")
  MEAS(s256f, 256,  11, 3, 1, 4, 1, 0, 1, 1, 0, "256 u=FF")
  MEAS(s256s, 256,  11, 3, 1, 4, 1, 0, 1, 0, 0, "256 u=SRAM")
  MEAS(s512f, 512,  11, 3, 1, 4, 1, 0, 1, 1, 0, "512 u=FF")
  MEAS(s512s, 512,  11, 3, 1, 4, 1, 0, 1, 0, 0, "512 u=SRAM")
  MEAS(s1Kf, 1024,  11, 3, 1, 4, 1, 0, 1, 1, 0, "1K  u=FF")
  MEAS(s1Ks, 1024,  11, 3, 1, 4, 1, 0, 1, 0, 0, "1K  u=SRAM")
  MEAS(s2Kf, 2048,  11, 3, 1, 4, 1, 0, 1, 1, 0, "2K  u=FF")
  MEAS(s2Ks, 2048,  11, 3, 1, 4, 1, 0, 1, 0, 0, "2K  u=SRAM")

  // ---- 2. Banking: N and banks (SIZE=256, SRAM u) ----
  // Shows cost of increasing branch throughput.
  // N=total branches, BxBPB = BANKS x branches_per_bank.
  SECTION("2. Banking: vary N and BANKS (SIZE=256, SRAM u, shared tag)")
  MEAS(n1_1x1, 256, 11, 3, 1, 1, 1, 0, 1, 0, 0, "N=1  1x1")
  MEAS(n2_2x1, 256, 11, 3, 1, 2, 2, 0, 1, 0, 0, "N=2  2x1")
  MEAS(n4_2x2, 256, 11, 3, 1, 4, 2, 0, 1, 0, 0, "N=4  2x2")
  MEAS(n4_4x1, 256, 11, 3, 1, 4, 4, 0, 1, 0, 0, "N=4  4x1")
  MEAS(n8_2x4, 256, 11, 3, 1, 8, 2, 0, 1, 0, 0, "N=8  2x4")
  MEAS(n8_4x2, 256, 11, 3, 1, 8, 4, 0, 1, 0, 0, "N=8  4x2")

  // ---- 3. Entry width: tag and counter (SIZE=256, N=4, 1 bank, SRAM u) ----
  // Wider entries = more SRAM bits, higher read latency.
  SECTION("3. Entry width: TAG/CTR (SIZE=256, N=4, 1 bank, SRAM u)")
  MEAS(t7c2,  256,  7, 2, 1, 4, 1, 0, 1, 0, 0, "tag=7  ctr=2")
  MEAS(t7c3,  256,  7, 3, 1, 4, 1, 0, 1, 0, 0, "tag=7  ctr=3")
  MEAS(t11c2, 256, 11, 2, 1, 4, 1, 0, 1, 0, 0, "tag=11 ctr=2")
  MEAS(t11c3, 256, 11, 3, 1, 4, 1, 0, 1, 0, 0, "tag=11 ctr=3")
  MEAS(t11c4, 256, 11, 4, 1, 4, 1, 0, 1, 0, 0, "tag=11 ctr=4")
  MEAS(t13c3, 256, 13, 3, 1, 4, 1, 0, 1, 0, 0, "tag=13 ctr=3")

  // ---- 4. U-bit width and storage (SIZE=256, N=4) ----
  // Compares U=1 vs U=2 and FF vs SRAM storage for each.
  SECTION("4. U-bit width x storage (SIZE=256, N=4, 1 bank)")
  MEAS(u1ff, 256, 11, 3, 1, 4, 1, 0, 1, 1, 0, "U=1 FF")
  MEAS(u1sr, 256, 11, 3, 1, 4, 1, 0, 1, 0, 0, "U=1 SRAM")
  MEAS(u2ff, 256, 11, 3, 2, 4, 1, 0, 1, 1, 0, "U=2 FF")
  MEAS(u2sr, 256, 11, 3, 2, 4, 1, 0, 1, 0, 0, "U=2 SRAM")

  // ---- 5. Features: ahead / per-bank / cache (SIZE=256, SRAM u) ----
  // Each row adds one feature to show its incremental cost.
  // base = N=4, 2 banks, no ahead, shared tag, SRAM u.
  SECTION("5a. Features: ahead + per-bank tags (SIZE=256, N=4, SRAM u)")
  MEAS(base,  256, 11, 3, 1, 4, 2, 0, 1, 0, 0, "base 2x2")
  MEAS(ahead, 256, 11, 3, 1, 4, 2, 1, 1, 0, 0, "+ahead")
  MEAS(pbank, 256, 11, 3, 1, 4, 4, 0, 0, 0, 0, "+per-bank")

  SECTION("5b. Features: FF cache (SIZE=256, N=8, 2 banks, SRAM u)")
  MEAS(cbase,      256, 11, 3, 1, 8, 2, 0, 1, 0, 0, "base 2x4")
  MEAS_CACHE(cache, 256, 11, 3, 1, 8, 2, 0, 1, 0, 1, "+cache")
  MEAS_CACHE(ah_ca, 256, 11, 3, 1, 8, 2, 1, 1, 0, 1, "+ahead+cache")
  MEAS_CACHE(allsr, 256, 11, 3, 1, 8, 2, 1, 0, 0, 1, "+ahd+pbnk+ca")

  // ---- 6. Realistic large configs (1K/2K, SRAM u) ----
  SECTION("6. Large tables (SRAM u)")
  MEAS(l1K,          1024, 11, 3, 1, 8, 4, 0, 1, 0, 0, "1K N=8 4x2")
  MEAS(l1Kah,        1024, 11, 3, 1, 8, 4, 1, 1, 0, 0, "1K +ahead")
  MEAS_CACHE(l1Kall, 1024, 11, 3, 1, 8, 4, 1, 0, 0, 1, "1K +ahd+pb+c")
  MEAS(l2K,          2048, 11, 3, 1, 4, 1, 0, 1, 0, 0, "2K N=4 1x4")
  MEAS(l2Kn1,        2048, 11, 3, 1, 1, 1, 0, 1, 0, 0, "2K N=1 1x1")

  printf("\n");
  printf("Column legend:\n");
  printf("  SIZE   = table entries (rows)      TAG    = tag width (bits)       "
         "CTR = counter width (bits)\n");
  printf("  U      = u-bit width (bits)        N      = branches per block     "
         "BxBPB = banks x branches_per_bank\n");
  printf("  A      = ahead pipeline (0/1)      S      = shared tag (0/1)       "
         "F   = u-bit in FF (0=SRAM, 1=FF)\n");
  printf("  C      = FF cache (0/1)\n");
  printf("  bits   = SRAM+FF storage (bits)    trans  = transistors            "
         "area = SRAM area (mm2)\n");
  printf("  staPwr = static power (mW)         (leakage, independent of "
         "workload)\n");
  printf("  Hit    = getHit latency (ps)       Tag    = getTag latency (ps)    "
         "Ctr = getCounter latency (ps)\n");
  printf("  Ubit   = getU latency (ps)         Reuse  = reuseRead latency "
         "(ps), - = no FF cache\n");

  return 0;
}
