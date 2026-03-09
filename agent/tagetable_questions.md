# TageTable Design Questions

## Parameters
- [x] 1. `PRED_BLK_SIZE = N` (same parameter) ✅
- [x] 2. `NUM_BANKS` partitions branches: `BPB = N / NUM_BANKS`. `USE_AHEAD` (bool) doubles banks for 1-ahead pipelining. Total physical banks = `NUM_BANKS * (USE_AHEAD ? 2 : 1)`. ✅
- [x] 3. Ahead banking (path) and branch-slot banking are orthogonal dimensions. `NUM_BANKS` handles branch slots, `USE_AHEAD` handles the path duplication. RMW comes free with ahead; without ahead, register-cache avoids R/W conflicts. ✅
- [x] 4. `DECAY_CTR` is a sub-parameter of `U_STORAGE = SRAM`. Only meaningful for probabilistic per-access decay. ✅
- [x] 5. `U_STORAGE = FF` uses a `ResetFn` functor `(val<U_WIDTH>, val<MODE_BITS>) -> val<U_WIDTH>`. Default provides reset/rshift/decrement. `U_STORAGE = SRAM` uses `DECAY_CTR`, no global reset. Mutually exclusive. ✅
- [x] 12. `SHARED_TAG` (bool) — applies to branch-slot banks only (NUM_BANKS dimension). Ahead duplication (×2) always has independent state. Default: shared (true). Per-bank needed when banks track independent contexts. ✅
- [x] 13. `SHARED_U` (bool) — follows same logic as `SHARED_TAG`. Per-branch U only useful with per-bank tags. Constraint: `SHARED_U = false` requires `SHARED_TAG = false` (static_assert). ✅
- [x] 14. Entry layout is an internal implementation detail. Table packs fields into SRAM based on parameter combination. Shared tag/U stored once (separate small RAM or designated bank). Predictor doesn't care about layout. ✅
- [x] 15. `USE_FF_CACHE` (bool) — cache SRAM reads in FFs for reuse across block. Only useful when BPB > 1. Constraint: `static_assert(!USE_FF_CACHE || BPB > 1)`. ✅
- [x] 16. No separate confidence bits. Confidence is derived from counter value (distance from threshold). Sufficient for TAGE. ✅
- [x] 17. No separate hysteresis bits. `CTR_WIDTH > 1` embeds hysteresis in counter saturation. Reference tage's 1-bit pred + 2-bit hyst is just a different encoding. ✅

## Interface
- [x] 6. Predictor manages the ahead pipeline. Predictor also owns the fold registers. Table is a banked RAM wrapper with register caching — no fold, no pipeline, no timing awareness. ✅
- [x] 7. Predictor manages path/XB. Predictor passes a bank index to the table. ✅
- [x] 8. Predictor passes bank index for writes. ✅
- [x] 9. Accessor methods (getHit(), getTag(), getCounter(slot), getU()). Matches HARCOM fanout control. ✅

## Scope
- [x] 10. Table takes a comparison tag on read, compares against all banks, exposes `getHit(bank)`. Predictor manages which tag to pass (handles ahead buffering). ✅
- [x] 11. Same write interface. No separate allocation method. Predictor owns allocation policy and initial values, casts to table's CTR_WIDTH/U_WIDTH. Table just stores what it's given. ✅
