#include "predictors/always_taken.hpp"
#include "predictors/bimodal.hpp"
#include "predictors/bimodalN.hpp"
#include "predictors/custom/Tage.hpp"
#include "predictors/custom/TageAhead.hpp"
#include "predictors/custom/TageAheadHC.hpp"
#include "predictors/custom/TageDirect.hpp"
#include "predictors/custom/TageDirectBim.hpp"
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

// 1-cycle config: 14 tables, GradedSize 2048→512, P2 ≈ 0.99
using TageAhead1C =
    TageAhead<TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,
                            ta::UniformTag<11>, ta::GradedSize<512, 2048>>,
              7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,
              UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6,
              2, 1024, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,
              SiblingPolicy::ALL>;

// ============================================================================
// Sweep helpers: 1C base macro, parameterized by AllocConfig
// ============================================================================
#define TA1C_BASE(ALLOC_CFG)                                                   \
  TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                         \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, 2,     \
      1024, 2, 256, true, HistUpdate::PATH, ALLOC_CFG, SiblingPolicy::ALL,    \
      0, 10, 10

// S1 base with variable SecTagPolicy (STP)
// S1 decay params + all defaults through REVERSE_TABLE_ORDER, then STP
#define S1_STP(STP)                                                            \
  TA1C_BASE(TAAllocPressSkip),                                                 \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false,           \
      ta::DefaultEpochTrigger, false, true, false, 0, false, STP

// S1 base with variable META_WIDTH (MW) and META_CAPACITY (MC)
#define S1_META(MW, MC)                                                        \
  TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                         \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 8192, false, 6, MW,    \
      MC, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,                    \
      SiblingPolicy::ALL, 0, 10, 10,                                           \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false

// S1 base with variable META_WIDTH/CAPACITY and doubled bimodal (16384)
#define S1_META_BIM2X(MW, MC)                                                  \
  TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                         \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 16384, false, 6, MW,   \
      MC, 2, 256, true, HistUpdate::PATH, TAAllocPressSkip,                    \
      SiblingPolicy::ALL, 0, 10, 10,                                           \
      true, DecayMiss::TAG_OR_SEC, DecayOp::DECREMENT,                         \
      ta::uniform_array<u64, 14>(8), ta::FixedDecayThresh<8>, false

// 1C base with doubled bimodal (16384 entries)
#define TA1C_BASE_BIM2X(ALLOC_CFG)                                             \
  TATableConfig<14, 1024, 11, 8, 200, 1, ta::HistSeries::GEOMETRIC,           \
                ta::UniformTag<11>, ta::GradedSize<512, 2048>>,                \
      7, 6, 5, true, 1, ta::Xor3SecTagHash5, 1, 2, 2,                         \
      UMispPolicy::UNTOUCHED, UClearPolicy::DECREMENT, 16384, false, 6, 2,    \
      1024, 2, 256, true, HistUpdate::PATH, ALLOC_CFG, SiblingPolicy::ALL,    \
      0, 10, 10

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

#ifdef PREDICTOR
using branch_predictor = PREDICTOR;
#else
using branch_predictor = TageAhead1C;
#endif
