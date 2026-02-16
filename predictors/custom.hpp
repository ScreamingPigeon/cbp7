

// Suppress warn
#include <memory>
#pragma GCC diagnostic ignored "-Wunused-parameter"

#ifndef CUSTOM_HPP
#define CUSTOM_HPP

#include "../cbp.hpp"
#include "../harcom.hpp"
#include "custom/Tage.hpp"
#include "custom/TageTable.hpp"

using namespace hcm;

// This instantiates TAGE AND our base predictor

template <
    // Tage Table Default params
    u64 T_TAG_WIDTH = 7, u64 T_U_WIDTH = 2, u64 T_CTR_WIDTH = 3,
    u64 T_HYS_WIDTH = 2, bool T_USE_HYS = true, u64 T_DECAY_CTR = 1024,
    // Base Predictor Patameters
    // TODO

    // Tage Predictor Default params
    u64 PRED_BLK_SIZE = 8, u64 NUM_TABLES = 5,
    std::array<u64, PRED_BLK_SIZE> T_TABLE_HIST_LEN = {8, 16, 32, 64, 128},
    std::array<u64, PRED_BLK_SIZE> T_TABLE_HIST_SIZE = {128, 64, 32, 16, 16},
    u64 GLOBAL_HIST_LEN = 64, bool USE_PATH_HIST = false,
    u64 PATH_HIST_LEN = 16, bool SPLIT_KSPACE = false, bool EN_N_BLK_RD = false>
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
  // Instantiate TAGE HERE
  Tage<T_TAG_WIDTH, T_U_WIDTH, T_CTR_WIDTH, T_HYS_WIDTH, T_USE_HYS, T_DECAY_CTR,
       PRED_BLK_SIZE, NUM_TABLES, T_TABLE_HIST_LEN, T_TABLE_HIST_SIZE,
       GLOBAL_HIST_LEN, USE_PATH_HIST, PATH_HIST_LEN, SPLIT_KSPACE, EN_N_BLK_RD> tage;
  // Instantiate Base Predictor here
};
#endif // CUSTOM_HPP
