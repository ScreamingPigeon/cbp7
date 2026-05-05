#include "predictors/always_taken.hpp"
#include "predictors/bimodal.hpp"
#include "predictors/bimodalN.hpp"
#include "predictors/custom/Tage.hpp"
#include "predictors/custom/TageAhead.hpp"
#include "predictors/custom/TageAheadHC.hpp"
#include "predictors/custom/TageAheadHC_IR.hpp"
#include "predictors/custom/TageDirect.hpp"
// #include "predictors/custom/TageDirectBim.hpp"
#include "predictors/experiment_perceptron.hpp"
#include "predictors/gshare.hpp"
#include "predictors/gshareN.hpp"
#include "predictors/gshareN_ahead.hpp"
#include "predictors/gshareN_ahead_best.hpp"
#include "predictors/never_taken.hpp"
#include "predictors/perceptron.hpp"
#include "predictors/tage.hpp"
#include "predictors/tutorial/tutorial.hpp"

// Keep `perceptron<>` as a stable user-facing name.
template <auto... Args> using perceptron = experiment_perceptron<Args...>;

// ============================================================================
// Competition configs: 1-cycle and 2-cycle tracks
// ============================================================================

// 1-cycle config: see TageAhead1C typedef below TA1C_BASE macro

// ============================================================================
// Sweep helpers: 1C base macro, parameterized by AllocConfig
// ============================================================================
#define TA1C_BASE(ALLOC_CFG)                                                   \
  TATableConfig<15, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1,                               \
      7, false, /* BR_P_ENTRY=N, INTERLEAVED=false */                          \
      2, 2,                                                                    \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 2,     \
      1024, 2, 256, true, HistUpdate::PATH, ALLOC_CFG, SiblingPolicy::ALL,    \
      0, 10, 10

using TageAhead1C = TageAhead<TA1C_BASE(TAAllocPressSkip)>;

// BR_P_ENTRY sweep: test per-group tag encoding
#define TA1C_BPE(BPE, INTLV, ALLOC_CFG)                                       \
  TATableConfig<15, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1,                               \
      BPE, INTLV,                                                              \
      2, 2,                                                                    \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 2,     \
      1024, 2, 256, true, HistUpdate::PATH, ALLOC_CFG, SiblingPolicy::ALL,    \
      0, 10, 10

using TA1C_BPE4 = TageAhead<TA1C_BPE(4, false, TAAllocPressSkip)>;  // NUM_GROUPS=2
using TA1C_BPE2 = TageAhead<TA1C_BPE(2, false, TAAllocPressSkip)>;  // NUM_GROUPS=4
using TA1C_BPE1 = TageAhead<TA1C_BPE(1, false, TAAllocPressSkip)>;  // NUM_GROUPS=7
using TA1C_BPE4I = TageAhead<TA1C_BPE(4, true, TAAllocPressSkip)>;  // NUM_GROUPS=2, interleaved
using TA1C_BPE2I = TageAhead<TA1C_BPE(2, true, TAAllocPressSkip)>;  // NUM_GROUPS=4, interleaved
using TA1C_BPE1I = TageAhead<TA1C_BPE(1, true, TAAllocPressSkip)>;  // NUM_GROUPS=7, interleaved

// S1 base with variable SecTagPolicy (STP)
// S1 decay params + all defaults through REVERSE_TABLE_ORDER, then STP
#define S1_STP(STP)                                                            \
  TA1C_BASE(TAAllocPressSkip),                                                 \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      ta::uniform_array<u64, 15>(8), ta::FixedDecayThresh<8>, false,           \
      ta::DefaultEpochTrigger, false, true, false, 0, false, STP

// S1 base with variable META_WIDTH (MW) and META_CAPACITY (MC)
#define S1_META(MW, MC)                                                        \
  TATableConfig<15, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                         \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, MW,    \
      MC, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,                    \
      SiblingPolicy::ALL, 0, 10, 10,                                           \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      ta::uniform_array<u64, 15>(8), ta::FixedDecayThresh<8>, false

// S1 base with variable META_WIDTH/CAPACITY and doubled bimodal (16384)
#define S1_META_BIM2X(MW, MC)                                                  \
  TATableConfig<15, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                         \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 16384, false, 6, MW,   \
      MC, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,                    \
      SiblingPolicy::ALL, 0, 10, 10,                                           \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      ta::uniform_array<u64, 15>(8), ta::FixedDecayThresh<8>, false

// 1C base with doubled bimodal (16384 entries)
#define TA1C_BASE_BIM2X(ALLOC_CFG)                                             \
  TATableConfig<15, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                         \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 16384, false, 6, 2,    \
      1024, 2, 256, true, HistUpdate::PATH, ALLOC_CFG, SiblingPolicy::ALL,    \
      0, 10, 10

#if 0 // temporarily disabled for FLOORPLAN build
// ============================================================================
// Sweep 4: Sec-tag policy tuning on S1 base
// ============================================================================

// A. SecTagNone — disable sec-tag entirely (baseline counterfactual)
using S4_None = TageAhead<S1_STP(ta::SecTagNone)>;

// B. SecTagFloor — skip sec-tag for long-history tables (T0..T(F-1))
using S4_Floor4  = TageAhead<S1_STP(ta::SecTagFloor<4>)>;   // skip T0-T3
using S4_Floor7  = TageAhead<S1_STP(ta::SecTagFloor<7>)>;   // skip T0-T6 (half)
using S4_Floor10 = TageAhead<S1_STP(ta::SecTagFloor<10>)>;  // skip T0-T9

// C. SecTagCeil — skip sec-tag for short-history tables (T(C)..T13)
using S4_Ceil4  = TageAhead<S1_STP(ta::SecTagCeil<4>)>;   // only T0-T3 check
using S4_Ceil7  = TageAhead<S1_STP(ta::SecTagCeil<7>)>;   // only T0-T6 check
using S4_Ceil10 = TageAhead<S1_STP(ta::SecTagCeil<10>)>;  // only T0-T9 check

// D. SecTagPressGated — skip sec-tag under high allocation pressure
using S4_Press256  = TageAhead<S1_STP(ta::SecTagPressGated<256>)>;   // ~25%
using S4_Press512  = TageAhead<S1_STP(ta::SecTagPressGated<512>)>;   // ~50%
using S4_Press768  = TageAhead<S1_STP(ta::SecTagPressGated<768>)>;   // ~75%

// E. SecTagAccGated — skip sec-tag when accuracy is high
using S4_Acc256  = TageAhead<S1_STP(ta::SecTagAccGated<256>)>;   // skip when acc>256
using S4_Acc512  = TageAhead<S1_STP(ta::SecTagAccGated<512>)>;   // skip when acc>512
using S4_Acc768  = TageAhead<S1_STP(ta::SecTagAccGated<768>)>;   // skip when acc>768

// F. SecTagPressGatedFloor — skip T0-T6 always + rest pressure-gated
using S4_PGF7_512_Policy = ta::SecTagPressGatedFloor<7, 512>;
using S4_PGF7_512 = TageAhead<S1_STP(S4_PGF7_512_Policy)>;

// ============================================================================
// Sweep 5: SecTagAdaptive — benefit-tracking adaptive sec-tag policy
// ============================================================================
// 8-bit counter, vary threshold
using S5_A8_64_P   = ta::SecTagAdaptive<8, 64>;
using S5_A8_96_P   = ta::SecTagAdaptive<8, 96>;
using S5_A8_128_P  = ta::SecTagAdaptive<8, 128>;
using S5_A8_160_P  = ta::SecTagAdaptive<8, 160>;
using S5_A8_192_P  = ta::SecTagAdaptive<8, 192>;
using S5_A8_64     = TageAhead<S1_STP(S5_A8_64_P)>;
using S5_A8_96     = TageAhead<S1_STP(S5_A8_96_P)>;
using S5_A8_128    = TageAhead<S1_STP(S5_A8_128_P)>;
using S5_A8_160    = TageAhead<S1_STP(S5_A8_160_P)>;
using S5_A8_192    = TageAhead<S1_STP(S5_A8_192_P)>;
// 10-bit counter, vary threshold
using S5_A10_256_P = ta::SecTagAdaptive<10, 256>;
using S5_A10_512_P = ta::SecTagAdaptive<10, 512>;
using S5_A10_768_P = ta::SecTagAdaptive<10, 768>;
using S5_A10_256   = TageAhead<S1_STP(S5_A10_256_P)>;
using S5_A10_512   = TageAhead<S1_STP(S5_A10_512_P)>;
using S5_A10_768   = TageAhead<S1_STP(S5_A10_768_P)>;

// ============================================================================
// Sweep 3: Meta table tuning (META_WIDTH × META_CAPACITY) on S1 base
// ============================================================================
// MW=1
using S3_MW1_MC512  = TageAhead<S1_META(1, 512)>;
using S3_MW1_MC1024 = TageAhead<S1_META(1, 1024)>;
using S3_MW1_MC2048 = TageAhead<S1_META(1, 2048)>;
using S3_MW1_MC4096 = TageAhead<S1_META(1, 4096)>;
// MW=2
using S3_MW2_MC512  = TageAhead<S1_META(2, 512)>;
using S3_MW2_MC1024 = TageAhead<S1_META(2, 1024)>; // == Sweep2_S1
using S3_MW2_MC2048 = TageAhead<S1_META(2, 2048)>;
using S3_MW2_MC4096 = TageAhead<S1_META(2, 4096)>;
// MW=4
using S3_MW4_MC512  = TageAhead<S1_META(4, 512)>;
using S3_MW4_MC1024 = TageAhead<S1_META(4, 1024)>;
using S3_MW4_MC2048 = TageAhead<S1_META(4, 2048)>;
using S3_MW4_MC4096 = TageAhead<S1_META(4, 4096)>;

// ============================================================================
// Sweep 6: Graded LFSR widths on S3_MW4_MC2048 base
// ============================================================================
// Baseline (S3_MW4_MC2048): uniform LFSR=8, FixedThresh<8> → 3.1% decay/miss
// Graded LFSR: T0 (longest hist) wider LFSR (less decay), T13 narrower (more)

// LFSR width arrays for sweep 6
inline constexpr auto S6_LFSR_U8  = ta::uniform_array<u64, 14>(8);
inline constexpr auto S6_LFSR_G10_7 = ta::graded_array<u64, 14>(10, 7);
inline constexpr auto S6_LFSR_G12_6 = ta::graded_array<u64, 14>(12, 6);

// Helper: S3_MW4_MC2048 base + custom decay (LFSR via constexpr var, thresh via type)
#define S6_DECAY(LFSR_VAR, THRESH_FN)                                          \
  TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                         \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 4,     \
      2048, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,                  \
      SiblingPolicy::ALL, 0, 10, 10,                                           \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      LFSR_VAR, THRESH_FN, false

using S6_G10_7     = TageAhead<S6_DECAY(S6_LFSR_G10_7, ta::FixedDecayThresh<8>)>;
using S6_G12_6     = TageAhead<S6_DECAY(S6_LFSR_G12_6, ta::FixedDecayThresh<8>)>;
using S6_G10_7_T16 = TageAhead<S6_DECAY(S6_LFSR_G10_7, ta::FixedDecayThresh<16>)>;
using S6_G12_6_T16 = TageAhead<S6_DECAY(S6_LFSR_G12_6, ta::FixedDecayThresh<16>)>;
// GradedThresh and PressGated need type aliases to avoid comma-in-macro
using S6_GT8_64_Thresh = ta::GradedDecayThresh<8, 64, 14>;
using S6_PG512_Thresh  = ta::PressGatedDecayThresh<4, 32, 14, 512>;
using S6_GT8_64    = TageAhead<S6_DECAY(S6_LFSR_U8, S6_GT8_64_Thresh)>;
using S6_PG512     = TageAhead<S6_DECAY(S6_LFSR_U8, S6_PG512_Thresh)>;

// S3_MW4_MC2048 + MAX_ALLOC=2
using S6_Alloc2 = TageAhead<
  TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 4,
      2048, 2, 256, true, HistUpdate::PATH, TAAlloc2PressSkip,
      SiblingPolicy::ALL, 0, 10, 10,
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
      ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false>;

// S3_MW4_MC2048 + FB 16K
using S6_FB16K = TageAhead<
  TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 16384, false, 6, 4,
      2048, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,
      SiblingPolicy::ALL, 0, 10, 10,
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
      ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false>;

// S3_MW4_MC2048 + META 4096
using S6_MC4096 = TageAhead<
  TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 4,
      4096, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,
      SiblingPolicy::ALL, 0, 10, 10,
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
      ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false>;

// S3_MW4_MC2048 + ALL THREE (MAX_ALLOC=2 + FB 16K + META 4096)
using S6_All3 = TageAhead<
  TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 16384, false, 6, 4,
      4096, 2, 256, true, HistUpdate::PATH, TAAlloc2PressSkip,
      SiblingPolicy::ALL, 0, 10, 10,
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
      ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false>;

// ============================================================================
// Sweep 7: Structural accuracy improvements on S3_MW4_MC2048 base
// ============================================================================
// Helper: S3 base with custom TableCfg and SecTagPolicy, rest standard
#define S7_BASE(TABLE_CFG, STP)                                                \
  TABLE_CFG,                                                                   \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                         \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 4,     \
      2048, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,                  \
      SiblingPolicy::ALL, 0, 10, 10,                                           \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false,           \
      ta::DefaultEpochTrigger, false, true, false, 0, false, STP

// Pre-declare TableCfg and SecTagPolicy types to avoid comma-in-macro
using S7_TC11 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                              ta::UniformTag<11>, ta::GradedSize<512, 2048>>;
using S7_TC12 = TATableConfig<14, 1024, 12, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                              ta::UniformTag<12>, ta::GradedSize<512, 2048>>;
using S7_TC_GT13_9 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                              ta::GradedTag<13, 9>, ta::GradedSize<512, 2048>>;
using S7_TC_GT14_8 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                              ta::GradedTag<14, 8>, ta::GradedSize<512, 2048>>;
using S7_TC_GT12_10 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                              ta::GradedTag<12, 10>, ta::GradedSize<512, 2048>>;
// History length variants
using S7_TC_H6_300 = TATableConfig<14, 1024, 11, 6, 300, 1, ta::HistSeries::GEOMETRIC,
                              ta::UniformTag<11>, ta::GradedSize<512, 2048>>;
using S7_TC_H5_150 = TATableConfig<14, 1024, 11, 5, 150, 1, ta::HistSeries::GEOMETRIC,
                              ta::UniformTag<11>, ta::GradedSize<512, 2048>>;
using S7_TC_H10_250 = TATableConfig<14, 1024, 11, 10, 250, 1, ta::HistSeries::GEOMETRIC,
                              ta::UniformTag<11>, ta::GradedSize<512, 2048>>;

// SecTagPolicy aliases
using S7_STP_A96  = ta::SecTagAdaptive<8, 96>;
using S7_STP_A128 = ta::SecTagAdaptive<8, 128>;
using S7_STP_A64  = ta::SecTagAdaptive<8, 64>;

// --- A. Adaptive sec-tag (3 thresholds) ---
using S7_Adapt64   = TageAhead<S7_BASE(S7_TC11, S7_STP_A64)>;
using S7_Adapt96   = TageAhead<S7_BASE(S7_TC11, S7_STP_A96)>;
using S7_Adapt128  = TageAhead<S7_BASE(S7_TC11, S7_STP_A128)>;

// --- B. Tag width: uniform 12-bit ---
using S7_Tag12     = TageAhead<S7_BASE(S7_TC12, ta::SecTagAll)>;

// --- C. Graded tags ---
using S7_GT13_9    = TageAhead<S7_BASE(S7_TC_GT13_9, ta::SecTagAll)>;   // T0:13, T13:9
using S7_GT14_8    = TageAhead<S7_BASE(S7_TC_GT14_8, ta::SecTagAll)>;   // T0:14, T13:8
using S7_GT12_10   = TageAhead<S7_BASE(S7_TC_GT12_10, ta::SecTagAll)>;  // T0:12, T13:10

// --- D. History length tuning ---
using S7_H6_300    = TageAhead<S7_BASE(S7_TC_H6_300, ta::SecTagAll)>;   // min=6, max=300
using S7_H5_150    = TageAhead<S7_BASE(S7_TC_H5_150, ta::SecTagAll)>;   // min=5, max=150
using S7_H10_250   = TageAhead<S7_BASE(S7_TC_H10_250, ta::SecTagAll)>;  // min=10, max=250

// ============================================================================
// Sweep 8: rwram bank shift — which address bit selects the bank
// ============================================================================
// Baseline: BANK_SHIFT=0 (bit 0). ~50% of buffered writes lost.
// Higher shifts use higher-order bits which may have less correlation
// between consecutive accesses.
// Tables: T0-4 are 1024 (10-bit idx), T5-13 are 2048 (11-bit idx).
// Max safe uniform shift = 9 (fits 10-bit smallest idx).

// Bank shift arrays
inline constexpr auto S8_BS0 = ta::uniform_array<u64, 14>(0);
inline constexpr auto S8_BS1 = ta::uniform_array<u64, 14>(1);
inline constexpr auto S8_BS2 = ta::uniform_array<u64, 14>(2);
// Graded: small tables (T0-4) shift=X, large tables (T5-13) shift=Y
inline constexpr auto S8_BS_GRAD = ta::split_array<u64, 14>(2, 3, 5);
inline constexpr auto S8_BS_GRAD2 = ta::split_array<u64, 14>(1, 3, 5);  // lo=1, hi=3
inline constexpr auto S8_BS_GRAD3 = ta::split_array<u64, 14>(3, 5, 5);  // lo=3, hi=5

using S8_TC_BS0 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                ta::UniformTag<11>, ta::GradedSize<512, 2048>, S8_BS0>;
using S8_TC_BS1 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                ta::UniformTag<11>, ta::GradedSize<512, 2048>, S8_BS1>;
using S8_TC_BS2 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                ta::UniformTag<11>, ta::GradedSize<512, 2048>, S8_BS2>;
using S8_TC_GRAD = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                 ta::UniformTag<11>, ta::GradedSize<512, 2048>, S8_BS_GRAD>;
using S8_TC_GRAD2 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                  ta::UniformTag<11>, ta::GradedSize<512, 2048>, S8_BS_GRAD2>;
using S8_TC_GRAD3 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                  ta::UniformTag<11>, ta::GradedSize<512, 2048>, S8_BS_GRAD3>;

using S8_Shift0 = TageAhead<S7_BASE(S8_TC_BS0, ta::SecTagAll)>;
using S8_Shift1 = TageAhead<S7_BASE(S8_TC_BS1, ta::SecTagAll)>;
using S8_Shift2 = TageAhead<S7_BASE(S8_TC_BS2, ta::SecTagAll)>;
using S8_Graded = TageAhead<S7_BASE(S8_TC_GRAD, ta::SecTagAll)>;
using S8_Graded2 = TageAhead<S7_BASE(S8_TC_GRAD2, ta::SecTagAll)>;
using S8_Graded3 = TageAhead<S7_BASE(S8_TC_GRAD3, ta::SecTagAll)>;

// ---- B=4 banks (2-bit bank select) ----
inline constexpr auto S8_B4 = ta::uniform_array<u64, 14>(4);
// B=4, shift=0 (bits [1:0])
using S8_TC_4B_S0 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                   ta::UniformTag<11>, ta::GradedSize<512, 2048>,
                                   S8_BS0, S8_B4>;
// B=4, shift=1 (bits [2:1])
using S8_TC_4B_S1 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                   ta::UniformTag<11>, ta::GradedSize<512, 2048>,
                                   S8_BS1, S8_B4>;
using S8_4B_S0 = TageAhead<S7_BASE(S8_TC_4B_S0, ta::SecTagAll)>;
using S8_4B_S1 = TageAhead<S7_BASE(S8_TC_4B_S1, ta::SecTagAll)>;

// ---- B=8 banks (3-bit bank select) ----
inline constexpr auto S8_B8 = ta::uniform_array<u64, 14>(8);
using S8_TC_8B_S1 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                   ta::UniformTag<11>, ta::GradedSize<512, 2048>,
                                   S8_BS1, S8_B8>;
using S8_8B_S1 = TageAhead<S7_BASE(S8_TC_8B_S1, ta::SecTagAll)>;

// ============================================================================
// Best combined: MW4/MC1024 + 16K bimodal + LFSR decay
// ============================================================================
using Best_8K  = TageAhead<S1_META(4, 1024)>;         // == S3_MW4_MC1024
using Best_16K = TageAhead<S1_META_BIM2X(4, 1024)>;   // 16K bim + MW4/MC1024

// S1 with doubled bimodal: Fixed thresh=8, TAG_OR_SEC, FB_CAPACITY=16384
using Sweep2_S1_Bim2x = TageAhead<TA1C_BASE_BIM2X(TAAllocPressSkip),
    true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false>;

// ============================================================================
// Sweep 2: tuning V1/V4 winners, MAX_ALLOC=2, utilization-gated decay
// ============================================================================

// --- A. Threshold tuning around V1 (fixed, TAG_OR_SEC) ---
// S1: V1 but thresh=8 (~3.1%), half of V1
using Sweep2_S1 = TageAhead<TA1C_BASE(TAAllocPressSkip),
    true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false>;

// S2: V1 but thresh=24 (~9.4%), 1.5x of V1
using Sweep2_S2 = TageAhead<TA1C_BASE(TAAllocPressSkip),
    true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<24>, false>;

// S3: V1 but thresh=32 (~12.5%), 2x of V1
using Sweep2_S3 = TageAhead<TA1C_BASE(TAAllocPressSkip),
    true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<32>, false>;

// --- B. Threshold tuning around V4 (graded, TAG only) ---
// S4: V4 but graded 8→64 (wider range, TAG only)
using Sweep2_S4 = TageAhead<TA1C_BASE(TAAllocPressSkip),
    true, DecayMiss::TAG, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8), ta::GradedDecayThresh<8, 64, 14>, false>;

// S5: V4 but graded 8→48 (moderate range, TAG only)
using Sweep2_S5 = TageAhead<TA1C_BASE(TAAllocPressSkip),
    true, DecayMiss::TAG, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8), ta::GradedDecayThresh<8, 48, 14>, false>;

// --- C. MAX_ALLOC=2 (double allocation per misprediction) ---
// S6: baseline + MAX_ALLOC=2 (epoch, no decay)
using Sweep2_S6 = TageAhead<TA1C_BASE(TAAlloc2PressSkip)>;

// S7: V1 + MAX_ALLOC=2 (best fixed decay + double alloc)
using Sweep2_S7 = TageAhead<TA1C_BASE(TAAlloc2PressSkip),
    true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<16>, false>;

// S8: V4 + MAX_ALLOC=2 (best graded decay + double alloc)
using Sweep2_S8 = TageAhead<TA1C_BASE(TAAlloc2PressSkip),
    true, DecayMiss::TAG, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8), ta::GradedDecayThresh<4, 32, 14>, false>;

// --- D. Utilization-gated decay via per-table LFSR widths ---
// Longer-history tables (T0) get wider LFSR → lower P(decay),
// shorter-history tables (T13) get narrower LFSR → higher P(decay).
// Combined with fixed threshold so total P varies per table.

// S9: V1 thresh=16 + graded LFSR 10→7 (T0: 16/1024≈1.6%, T13: 16/128≈12.5%)
using Sweep2_S9 = TageAhead<TA1C_BASE(TAAllocPressSkip),
    true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
    ta::graded_array<u64, 14>(10, 7), ta::FixedDecayThresh<16>, false>;

// S10: V1 thresh=16 + graded LFSR 12→6 (T0: 16/4096≈0.4%, T13: 16/64≈25%)
using Sweep2_S10 = TageAhead<TA1C_BASE(TAAllocPressSkip),
    true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
    ta::graded_array<u64, 14>(12, 6), ta::FixedDecayThresh<16>, false>;

// --- E. Pressure-gated decay (only decay under allocation pressure) ---
// S11: graded 4→32, TAG_OR_SEC, gated at alloc_ctr > 512 (~50% pressure)
using Sweep2_S11 = TageAhead<TA1C_BASE(TAAllocPressSkip),
    true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8),
    ta::PressGatedDecayThresh<4, 32, 14, 512>, false>;

// S12: graded 4→32, TAG_OR_SEC, gated at alloc_ctr > 256 (~25% pressure)
using Sweep2_S12 = TageAhead<TA1C_BASE(TAAllocPressSkip),
    true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8),
    ta::PressGatedDecayThresh<4, 32, 14, 256>, false>;

// S13: pressure-scaled graded 8→64, TAG_OR_SEC (thresh proportional to pressure)
using Sweep2_S13 = TageAhead<TA1C_BASE(TAAllocPressSkip),
    true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
    ta::uniform_array<u64, 14>(8),
    ta::PressScaledDecayThresh<8, 64, 14>, false>;

// 2-cycle config: 28 tables, StepSize 4096/2048 split@24, P2 ≈ 1.91
// 106K entries, 24 tables at 4096 + 4 at 2048, MAXH=1000
// SiblingPolicy::NONE — sibling skip hurts 2C (3.3% MPKI regression)
using TageAhead2C =
    TageAhead<TATableConfig<28, 2048, 12, 4, 1000, 2, ta::HistSeries::GEOMETRIC,
                            ta::UniformTag<12>, ta::StepSize<4096, 2048, 24>>,
              8, 6, 4, true, 1, ta::Xor3SecTagHash, 1, 2, 1,
              UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 4096, false, 6,
              4, 256, 2, 1024, true, HistUpdate::PATH, TADefaultAllocConfig,
              SiblingPolicy::NONE>;

// ============================================================================
// Decay sweep configs (all 2C base + SiblingPolicy::NONE)
// Common: 28 tables, StepSize<4096,2048,24>, TAG=12, MAXH=1000
// ============================================================================

// Helper: 2C base params up to SiblingPolicy::NONE, parameterized by U_WIDTH
#define TA2C_BASE_U(UW)                                                        \
  TATableConfig<28, 2048, 12, 4, 1000, 2, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<12>, ta::StepSize<4096, 2048, 24>>,             \
      8, 6, 4, true, 1, ta::Xor3SecTagHash, 1, 2, UW, UMispPolicy::UNTOUCHED,  \
      UClearPolicy::DECREMENT, 4096, false, 6, 4, 256, 2, 1024, true,          \
      HistUpdate::PATH, TADefaultAllocConfig, SiblingPolicy::NONE, 0, 16, 16
#define TA2C_BASE TA2C_BASE_U(1)

// A: Baseline — epoch enabled, no decay (already TageAhead2C)

// B: No management — neither epoch nor decay
using TA2C_B = TageAhead<TA2C_BASE, false, DecayMiss::TAG,
                         DecayOp::DECREMENT, // DECAY_ENABLE=false
                         ta::uniform_array<u64, 28>(8), ta::DefaultDecayThresh,
                         false>; // EPOCH_ENABLE=false

// C: TAG_OR_SEC, LW=8 (default decay, most aggressive trigger)
using TA2C_C =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, 28>(8), ta::DefaultDecayThresh, false>;

// D: TAG only, LW=8 (less aggressive trigger)
using TA2C_D =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG, DecayOp::DECREMENT,
              ta::uniform_array<u64, 28>(8), ta::DefaultDecayThresh, false>;

// E: TAG_AND_SEC, LW=8 (least aggressive trigger)
using TA2C_E =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG_AND_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, 28>(8), ta::DefaultDecayThresh, false>;

// F: TAG_OR_SEC, LW=12 (lower fire probability)
using TA2C_F =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, 28>(12), ta::DefaultDecayThresh, false>;

// G: TAG_OR_SEC, LW=6 (higher fire probability)
using TA2C_G =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::uniform_array<u64, 28>(6), ta::DefaultDecayThresh, false>;

// H: TAG_OR_SEC, per-table LW: 12 for T0-13, 6 for T14-27
using TA2C_H =
    TageAhead<TA2C_BASE, true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,
              ta::split_array<u64, 28>(12, 6, 14), ta::DefaultDecayThresh,
              false>;

// ============================================================================
// Sweep 6: Block size sweep — N (branches/block) × LINEINST
// Base = TageAhead1C (best 1S config)
// ============================================================================
#define BLKSWEEP(N_VAL, LINE_VAL)                                              \
  TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      N_VAL, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                     \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 2,     \
      1024, 2, LINE_VAL, true, HistUpdate::PATH, TAAllocPressSkip,             \
      SiblingPolicy::ALL

using BLK_N1_L4  = TageAhead<BLKSWEEP(1, 4)>;
using BLK_N1_L8  = TageAhead<BLKSWEEP(1, 8)>;
using BLK_N1_L16 = TageAhead<BLKSWEEP(1, 16)>;
using BLK_N2_L4  = TageAhead<BLKSWEEP(2, 4)>;
using BLK_N2_L8  = TageAhead<BLKSWEEP(2, 8)>;
using BLK_N2_L16 = TageAhead<BLKSWEEP(2, 16)>;
using BLK_N4_L4  = TageAhead<BLKSWEEP(4, 4)>;
using BLK_N4_L8  = TageAhead<BLKSWEEP(4, 8)>;
using BLK_N4_L16 = TageAhead<BLKSWEEP(4, 16)>;
using BLK_N1_L32 = TageAhead<BLKSWEEP(1, 32)>;
using BLK_N1_L64 = TageAhead<BLKSWEEP(1, 64)>;
using BLK_N1_L128 = TageAhead<BLKSWEEP(1, 128)>;
using BLK_N2_L32  = TageAhead<BLKSWEEP(2, 32)>;
using BLK_N2_L64  = TageAhead<BLKSWEEP(2, 64)>;
using BLK_N2_L128 = TageAhead<BLKSWEEP(2, 128)>;
using BLK_N2_L256 = TageAhead<BLKSWEEP(2, 256)>;
using BLK_N4_L32  = TageAhead<BLKSWEEP(4, 32)>;
using BLK_N4_L64  = TageAhead<BLKSWEEP(4, 64)>;
using BLK_N4_L128 = TageAhead<BLKSWEEP(4, 128)>;

#define HISTSWEEP(N_VAL, LINE_VAL, MINH_VAL, MAXH_VAL)                        \
  TATableConfig<14, 1024, 11, MINH_VAL, MAXH_VAL, 1, ta::HistSeries::GEOMETRIC, \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      N_VAL, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                     \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 2,     \
      1024, 2, LINE_VAL, true, HistUpdate::PATH, TAAllocPressSkip,             \
      SiblingPolicy::ALL

// N2 history sweep: MAXH=200 (baseline), 400, 600. MINH scaled proportionally.
using HS_N2_L8_H200   = TageAhead<HISTSWEEP(2, 8, 8, 200)>;
using HS_N2_L8_H400   = TageAhead<HISTSWEEP(2, 8, 16, 400)>;
using HS_N2_L8_H600   = TageAhead<HISTSWEEP(2, 8, 24, 600)>;
using HS_N2_L16_H200  = TageAhead<HISTSWEEP(2, 16, 8, 200)>;
using HS_N2_L16_H400  = TageAhead<HISTSWEEP(2, 16, 16, 400)>;
using HS_N2_L16_H600  = TageAhead<HISTSWEEP(2, 16, 24, 600)>;
using HS_N2_L128_H200 = TageAhead<HISTSWEEP(2, 128, 8, 200)>;
using HS_N2_L128_H400 = TageAhead<HISTSWEEP(2, 128, 16, 400)>;
using HS_N2_L128_H600 = TageAhead<HISTSWEEP(2, 128, 24, 600)>;
using HS_N2_L256_H200 = TageAhead<HISTSWEEP(2, 256, 8, 200)>;
using HS_N2_L256_H400 = TageAhead<HISTSWEEP(2, 256, 16, 400)>;
using HS_N2_L256_H600 = TageAhead<HISTSWEEP(2, 256, 24, 600)>;

#endif // #if 0

// ============================================================================
// Sweep 9: S1 base + BPE=1 + graded tags
//
// S7_BASE is broken after BPE commit (missing BR_P_ENTRY/INTERLEAVED).
// S1_GT_BASE / S1_BPE1_BASE are correct replacements.
// ============================================================================

#define S1_GT_BASE(TABLE_CFG)                                                  \
  TABLE_CFG,                                                                   \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1,                               \
      7, false, /* BR_P_ENTRY=7, INTERLEAVED=false */                          \
      2, 2,                                                                    \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 2,     \
      1024, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,                  \
      SiblingPolicy::ALL, 0, 10, 10,                                           \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false,           \
      ta::DefaultEpochTrigger, false, true, false, 0, false, ta::SecTagAll

#define S1_BPE1_BASE(TABLE_CFG)                                                \
  TABLE_CFG,                                                                   \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1,                               \
      1, false, /* BR_P_ENTRY=1, INTERLEAVED=false */                          \
      2, 2,                                                                    \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 2,     \
      1024, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,                  \
      SiblingPolicy::ALL, 0, 10, 10,                                           \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false,           \
      ta::DefaultEpochTrigger, false, true, false, 0, false, ta::SecTagAll

// Uniform tag widths
using S9_TC_U8  = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                ta::UniformTag<8>,  ta::GradedSize<512, 2048>>;
using S9_TC_U9  = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                ta::UniformTag<9>,  ta::GradedSize<512, 2048>>;
using S9_TC_U10 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                ta::UniformTag<10>, ta::GradedSize<512, 2048>>;
using S9_TC_U11 = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                ta::UniformTag<11>, ta::GradedSize<512, 2048>>;

// Graded tag configs (raw widths; eff = raw - GROUP_BITS)
// BPE=7: GROUP_BITS=0, BPE=1: GROUP_BITS=3
using S9_TC_GT10_6  = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                    ta::GradedTag<10, 6>, ta::GradedSize<512, 2048>>;
using S9_TC_GT11_7  = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                    ta::GradedTag<11, 7>, ta::GradedSize<512, 2048>>;
using S9_TC_GT12_8  = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                    ta::GradedTag<12, 8>, ta::GradedSize<512, 2048>>;
using S9_TC_GT13_9  = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                    ta::GradedTag<13, 9>, ta::GradedSize<512, 2048>>;
using S9_TC_GT14_8  = TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                    ta::GradedTag<14, 8>, ta::GradedSize<512, 2048>>;

// BPE=7 baselines (GROUP_BITS=0, raw tag = effective tag)
using S9_GT_U8      = TageAhead<S1_GT_BASE(S9_TC_U8)>;      // uniform 8
using S9_GT_U9      = TageAhead<S1_GT_BASE(S9_TC_U9)>;      // uniform 9
using S9_GT_U10     = TageAhead<S1_GT_BASE(S9_TC_U10)>;     // uniform 10
using S9_GT_U11     = TageAhead<S1_GT_BASE(S9_TC_U11)>;     // uniform 11 (S1 baseline)
using S9_GT_GT10_6  = TageAhead<S1_GT_BASE(S9_TC_GT10_6)>;  // graded 10→6
using S9_GT_GT11_7  = TageAhead<S1_GT_BASE(S9_TC_GT11_7)>;  // graded 11→7
using S9_GT_GT12_8  = TageAhead<S1_GT_BASE(S9_TC_GT12_8)>;  // graded 12→8
using S9_GT_GT13_9  = TageAhead<S1_GT_BASE(S9_TC_GT13_9)>;  // graded 13→9
using S9_GT_GT14_8  = TageAhead<S1_GT_BASE(S9_TC_GT14_8)>;  // graded 14→8

// BPE=1 with 8 banks, 15 tables (matching HC_IR hardware config) for monitor runs
using S9_TC_U11_8B = TATableConfig<15, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                                   ta::UniformTag<11>, ta::GradedSize<512, 2048>,
                                   ta::uniform_array<u64, 15>(1),
                                   ta::uniform_array<u64, 15>(8)>;
#define S1_BPE1_8B_BASE(TABLE_CFG)                                             \
  TABLE_CFG,                                                                   \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1,                               \
      1, false, /* BR_P_ENTRY=1, INTERLEAVED=false */                          \
      2, 2,                                                                    \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 2,     \
      1024, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,                  \
      SiblingPolicy::ALL, 0, 10, 10,                                           \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      ta::uniform_array<u64, 15>(8), ta::FixedDecayThresh<8>, false,           \
      ta::DefaultEpochTrigger, false, true, false, 0, false, ta::SecTagAll

using TA1C_BPE1_8B = TageAhead<S1_BPE1_8B_BASE(S9_TC_U11_8B)>;
// META_BANKS variants: per-branch meta counters
using TA1C_BPE1_8B_MB4 = TageAhead<S1_BPE1_8B_BASE(S9_TC_U11_8B), 4>;
using TA1C_BPE1_8B_MB7 = TageAhead<S1_BPE1_8B_BASE(S9_TC_U11_8B), 7>;

// BPE=1 (NUM_GROUPS=7, GROUP_BITS=3, eff = raw - 3)
using S9_BPE1_U8      = TageAhead<S1_BPE1_BASE(S9_TC_U8)>;      // raw 8,  eff 5
using S9_BPE1_U9      = TageAhead<S1_BPE1_BASE(S9_TC_U9)>;      // raw 9,  eff 6
using S9_BPE1_U10     = TageAhead<S1_BPE1_BASE(S9_TC_U10)>;     // raw 10, eff 7
using S9_BPE1_U11     = TageAhead<S1_BPE1_BASE(S9_TC_U11)>;     // raw 11, eff 8
using S9_BPE1_GT10_6  = TageAhead<S1_BPE1_BASE(S9_TC_GT10_6)>;  // graded 10→6, eff 7→3
using S9_BPE1_GT11_7  = TageAhead<S1_BPE1_BASE(S9_TC_GT11_7)>;  // graded 11→7, eff 8→4
using S9_BPE1_GT12_8  = TageAhead<S1_BPE1_BASE(S9_TC_GT12_8)>;  // graded 12→8, eff 9→5
using S9_BPE1_GT13_9  = TageAhead<S1_BPE1_BASE(S9_TC_GT13_9)>;  // graded 13→9, eff 10→6
using S9_BPE1_GT14_8  = TageAhead<S1_BPE1_BASE(S9_TC_GT14_8)>;  // graded 14→8, eff 11→5

// ============================================================================
// TageDirect matching TageAheadHC structural params (no ahead pipeline)
// ============================================================================
using TD_HC_TC = td::TDTableConfig<14, 1024, 11, 1, 2, 2, 8, 200, 1,
                                    td::HistSeries::GEOMETRIC,
                                    td::UniformTag<11>,
                                    td::StepSize<1024, 2048, 5>>;
using TD_HC = TageDirect<TD_HC_TC,
                         td::TDDefaultAllocConfig,
                         256,   // LINEINST
                         7,     // N
                         8, 2,  // DECAY_CTR (8-bit LFSR), DECAY_GRAN
                         td::TDDecayMild,
                         true, 8192, // P1_USE_GSHARE, P1_TABLE_SIZE
                         6,     // P1_HIST
                         true, 4, // USE_META, METABITS
                         2,     // METAPIPE
                         false, // USE_PATH_HIST
                         27, 6, // PATH_HIST_WIDTH, PATH_BITS
                         td::XORFold,
                         8, 1,  // RWRAM_BANKS=8, BANK_SHIFT=1
                         8,     // EPOCH_CTR_BITS
                         true,  // SHARED_HYS
                         false>; // USE_DIR_HIST

#ifdef PREDICTOR
using branch_predictor = PREDICTOR;
#else
using branch_predictor = TageAheadHC;
#endif
