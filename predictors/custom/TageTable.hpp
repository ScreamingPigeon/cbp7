#pragma once

#include "../../cbp.hpp"
#include "../../harcom.hpp"

template <u64 TABLE_SIZE = 64, u64 TABLE_HIST = 64, u64 TAG_WIDTH = 7,
          u64 CTR_WIDTH = 3, bool USE_HYS = true, u64 HYS_WIDTH = 2,
          u64 PRED_BLK_SIZE = 8>
class TageTable {
public:
  TageTable() {}
  ~TageTable() {}

private:
};
