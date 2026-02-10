

// Suppress warn
#pragma GCC diagnostic ignored "-Wunused-parameter"

#ifndef CUSTOM_HPP
#define CUSTOM_HPP

#include "../cbp.hpp"
#include "../harcom.hpp"
#include "custom/TageTable.hpp"

using namespace hcm;

/* ========== PREDICTOR-LEVEL PARAMETERS ========== */
/*
 * PRED_BLK_SIZE:     (u64) Size of fetch block
 *
 * TODO:
 * BASE_PRED:         (enum) Which base predictor to use for p0
 * BASE_PRED_SIZE:    (u64) Number of entries in base predictor
 * BASE_PRED_HYST:    (bool) Whether to use shared hysteresis
 *                    bits for base predictor
 * BASE_HYST_SIZE:    (u64) Size of shared hysteresis table
 *                    if enabled
 *
 * NUM_TABLES:        (u64) Number of TAGE_TABLES
 *                    (typically 7-14)
 * T_TABLE_HIST_LEN:  (const u64[]) Array specifying history
 *                    length for each table
 * T_TABLE_SIZE:      (const u64[]) Array specifying number of
 *                    entries for each table
 *
 * GLOBAL_HIST_LEN:   (u64) Maximum global history length
 *                    (max of all T_TABLE_HIST_LEN)
 * PATH_HIST_LEN:     (u64) Path history length
 *                    (typically 16-27 bits)
 * USE_PATH_HIST:     (bool) Whether to use path history in
 *                    indexing
 * SPLIT_KSPACE:      use a different history for kernel hist
 *
 * TODO:
 * USE_ALT_ON_NA:     (bool) Whether to use alternate prediction
 *                    on newly allocated entries
 * USE_ALT_COUNTERS:  (u64) Number of USE_ALT_ON_NA counters

 * TODO:
 * ALLOC_POLICY:      (enum) Allocation policy on misprediction
 * RESET_PERIOD:      (u64) Period for resetting useful counters
 * RESET_POLICY:      (enum) Policy for useful counter reset
 */

/* ========== TAGE TABLE-LEVEL PARAMETERS ========== */
/*
 * T_TAG_WIDTH:       (const u64[]) Tag width in bits for each
 *                    table
 * T_CTR_WIDTH:       (u64) Prediction counter width
 *
 * T_U_WIDTH:         (u64) Useful counter width
 *
 * T_HYS_WIDTH:       (u64) Hysterisis bit width
 * T_USE_HYS:         (bool) Use hysteresis btis
 *
 * TODO:T_HASH_FUNC:  (enum[]) Hash function type for each table
 *                    (XOR, folded-XOR, etc.)
 */

/* ========== OPTIONAL ENHANCEMENTS ========== */
/* TODO:
 * USE_LOOP_PRED:     (bool) Whether to include loop predictor
 * LOOP_PRED_SIZE:    (u64) Loop predictor table size
 * USE_STAT_CORR:     (bool) Whether to include SC
 * SC_NUM_TABLES:     (u64) Number of SC tables
 * SC_TABLE_SIZES:    (const u64[]) Sizes for SC tables
 */
template <
    // Tage Table Default params
    u64 T_TAG_WIDTH = 7, u64 T_CTR_WIDTH = 3, u64 T_HYS_WIDTH = 2,
    bool T_USE_HYS = true,
    // Base Predictor Patameters
    // TODO

    // Tage Predictor Default params
    u64 PRED_BLK_SIZE = 8, u64 NUM_TABLES = 5,
    std::array<u64, PRED_BLK_SIZE> T_TABLE_HIST_LEN = {8, 16, 32, 64, 128},
    std::array<u64, PRED_BLK_SIZE> T_TABLE_HIST_SIZE = {128, 64, 32, 16, 16},
    u64 GLOBAL_HIST_LEN = 128, bool USE_PATH_HIST = false,
    u64 PATH_HIST_LEN = 16, bool SPLIT_KSPACE = false>
class Custom : public predictor {
  /*
   * This is the top level class that implements *TEAM_NAME*'s Predictor
   * Please specify the parameters you need using the params.yaml file
   */
public:
  // Interface to simulator
  void new_block(val<64> inst_pc) {}
  val<1> predict1(val<64> inst_pc) { return val<1>(0); }
  val<1> reuse_predict1(val<64> inst_pc) { return val<1>(0); }
  val<1> predict2(val<64> inst_pc) { return val<1>(0); }
  val<1> reuse_predict2(val<64> inst_pc) { return val<1>(0); }
  void update_condbr(val<64> branch_pc, val<1> taken, val<64> next_pc) {}
  void update_cycle(instruction_info &block_end_info) {}

private:
  // TageTableTuple
  template <std::size_t... Is>
  static auto make_tables_impl(std::index_sequence<Is...>) {
    return std::tuple{
        TageTable<T_TABLE_HIST_SIZE[Is], T_TABLE_HIST_LEN[Is], T_TAG_WIDTH,
                  T_CTR_WIDTH, T_USE_HYS, T_HYS_WIDTH>{}...};
  }

  // Tuple of tables
  decltype(make_tables_impl(
      std::make_index_sequence<NUM_TABLES>{})) tage_tables;
};

#endif // CUSTOM_HPP
