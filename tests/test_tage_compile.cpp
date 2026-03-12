// Phase 1 compile test: verify all Tage instantiations compile clean.
// No runtime behavior — just type instantiation.

#include "../predictors/custom/Tage.hpp"

// ---- Test 1: Default config (uniform params, direct mode) ----
using T1 = Tage<>;

// ---- Test 2: Non-uniform table sizes ----
struct CustomTableConfig {
  static constexpr u64 NUM_TABLES = 4;
  static constexpr u64 MINHIST = 4;
  static constexpr u64 MAXHIST = 64;

  static constexpr std::array<u64, NUM_TABLES> TABLE_SIZE = {512, 1024, 2048,
                                                              4096};
  static constexpr std::array<u64, NUM_TABLES> TAG_WIDTH = {9, 11, 13, 15};
  static constexpr std::array<u64, NUM_TABLES> CTR_WIDTH = {3, 3, 3, 3};
  static constexpr std::array<u64, NUM_TABLES> HYST_WIDTH = {2, 2, 2, 2};
  static constexpr std::array<u64, NUM_TABLES> U_WIDTH = {1, 1, 2, 2};
  static constexpr auto HIST_LEN =
      geometric_hist<NUM_TABLES>(MINHIST, MAXHIST);
};
using T2 = Tage<CustomTableConfig>;

// ---- Test 3: Ahead mode (direct defaults otherwise) ----
using T3 = Tage<DefaultTableConfig, DefaultAllocConfig,
                /*FETCH_WIDTH=*/8, /*BIMODAL_SIZE=*/4096,
                /*BR_P_ENTRY=*/1, /*NUM_BANKS=*/1,
                /*USE_AHEAD=*/true>;

// ---- Test 4: P1 bimodal mode (P1_USE_GSHARE=false) ----
using T4 = Tage<DefaultTableConfig, DefaultAllocConfig,
                /*FETCH_WIDTH=*/8, /*BIMODAL_SIZE=*/4096,
                /*BR_P_ENTRY=*/1, /*NUM_BANKS=*/1,
                /*USE_AHEAD=*/false, /*SHARED_TAG=*/true,
                /*SHARED_U=*/true, /*SHARED_HYS=*/true,
                /*U_STOR_FF=*/false,
                /*DECAY_CTR=*/1024, /*ResetFn=*/DefaultResetFn,
                /*USE_FF_CACHE=*/false,
                /*P1_USE_GSHARE=*/false>;

// ---- Test 5: Meta-prediction with PC-indexed table ----
using T5 = Tage<DefaultTableConfig, DefaultAllocConfig,
                /*FETCH_WIDTH=*/8, /*BIMODAL_SIZE=*/4096,
                /*BR_P_ENTRY=*/1, /*NUM_BANKS=*/1,
                /*USE_AHEAD=*/false, /*SHARED_TAG=*/true,
                /*SHARED_U=*/true, /*SHARED_HYS=*/true,
                /*U_STOR_FF=*/false,
                /*DECAY_CTR=*/1024, /*ResetFn=*/DefaultResetFn,
                /*USE_FF_CACHE=*/false,
                /*P1_USE_GSHARE=*/true, /*P1_TABLE_SIZE=*/16384,
                /*P1_HIST=*/6,
                /*USE_META=*/true, /*METABITS=*/4, /*METAPIPE=*/2,
                /*META_TABLE_SIZE=*/256>;

// ---- Test 6: No meta-prediction ----
using T6 = Tage<DefaultTableConfig, DefaultAllocConfig,
                /*FETCH_WIDTH=*/8, /*BIMODAL_SIZE=*/4096,
                /*BR_P_ENTRY=*/1, /*NUM_BANKS=*/1,
                /*USE_AHEAD=*/false, /*SHARED_TAG=*/true,
                /*SHARED_U=*/true, /*SHARED_HYS=*/true,
                /*U_STOR_FF=*/false,
                /*DECAY_CTR=*/1024, /*ResetFn=*/DefaultResetFn,
                /*USE_FF_CACHE=*/false,
                /*P1_USE_GSHARE=*/true, /*P1_TABLE_SIZE=*/16384,
                /*P1_HIST=*/6,
                /*USE_META=*/false>;

// ---- Test 7: U-bits in flip-flops (epoch mode) ----
using T7 = Tage<DefaultTableConfig, DefaultAllocConfig,
                /*FETCH_WIDTH=*/8, /*BIMODAL_SIZE=*/4096,
                /*BR_P_ENTRY=*/1, /*NUM_BANKS=*/1,
                /*USE_AHEAD=*/false, /*SHARED_TAG=*/true,
                /*SHARED_U=*/true, /*SHARED_HYS=*/true,
                /*U_STOR_FF=*/true>;

// ---- Test 8: Path history enabled ----
using T8 = Tage<DefaultTableConfig, DefaultAllocConfig,
                /*FETCH_WIDTH=*/8, /*BIMODAL_SIZE=*/4096,
                /*BR_P_ENTRY=*/1, /*NUM_BANKS=*/1,
                /*USE_AHEAD=*/false, /*SHARED_TAG=*/true,
                /*SHARED_U=*/true, /*SHARED_HYS=*/true,
                /*U_STOR_FF=*/false,
                /*DECAY_CTR=*/1024, /*ResetFn=*/DefaultResetFn,
                /*USE_FF_CACHE=*/false,
                /*P1_USE_GSHARE=*/true, /*P1_TABLE_SIZE=*/16384,
                /*P1_HIST=*/6,
                /*USE_META=*/true, /*METABITS=*/4, /*METAPIPE=*/2,
                /*META_TABLE_SIZE=*/0,
                /*USE_PATH_HIST=*/true, /*PATH_HIST_WIDTH=*/27>;

// ---- Test 9: Ahead + custom config ----
using T9 = Tage<CustomTableConfig, DefaultAllocConfig,
                /*FETCH_WIDTH=*/8, /*BIMODAL_SIZE=*/4096,
                /*BR_P_ENTRY=*/1, /*NUM_BANKS=*/1,
                /*USE_AHEAD=*/true>;

// ---- Test 10: SweepTableConfig with defaults (should match DefaultTableConfig) ----
using T10 = Tage<SweepTableConfig<>>;

// ---- Test 11: SweepTableConfig with SIZE_RATIO=2 ----
using T11 = Tage<SweepTableConfig<8, 2048, 11, 1, 2, 1, 2, 100, 2>>;

// ---- Test 12: SweepTableConfig with different table count ----
using T12 = Tage<SweepTableConfig<4, 1024, 9, 2, 2, 1, 4, 64>>;

// Force instantiation of all types
void force_instantiate() {
  // Direct mode types — just verify the struct can be constructed
  // (HARCOM needs simulator context to actually create, but the types must
  //  be well-formed)
  static_assert(T1::NUM_TABLES == 8);
  static_assert(T2::NUM_TABLES == 4);
  static_assert(T3::NUM_TABLES == 8);  // ahead
  static_assert(T4::NUM_TABLES == 8);  // bimodal P1
  static_assert(T5::NUM_TABLES == 8);  // PC-indexed meta
  static_assert(T6::NUM_TABLES == 8);  // no meta
  static_assert(T7::NUM_TABLES == 8);  // FF u-bits
  static_assert(T8::NUM_TABLES == 8);  // path history
  static_assert(T9::NUM_TABLES == 4);  // ahead + custom
  static_assert(T10::NUM_TABLES == 8); // sweep default
  static_assert(T11::NUM_TABLES == 8); // sweep ratio=2
  static_assert(T12::NUM_TABLES == 4); // sweep 4 tables

  // Verify SweepTableConfig<> matches DefaultTableConfig
  static_assert(SweepTableConfig<>::TABLE_SIZE == DefaultTableConfig::TABLE_SIZE);
  static_assert(SweepTableConfig<>::TAG_WIDTH == DefaultTableConfig::TAG_WIDTH);
  static_assert(SweepTableConfig<>::CTR_WIDTH == DefaultTableConfig::CTR_WIDTH);
  static_assert(SweepTableConfig<>::HYST_WIDTH == DefaultTableConfig::HYST_WIDTH);
  static_assert(SweepTableConfig<>::U_WIDTH == DefaultTableConfig::U_WIDTH);
  static_assert(SweepTableConfig<>::HIST_LEN == DefaultTableConfig::HIST_LEN);
}

int main() {
  std::cerr << "Phase 1 compile test: all " << 12
            << " instantiations compiled successfully.\n";
  return 0;
}
