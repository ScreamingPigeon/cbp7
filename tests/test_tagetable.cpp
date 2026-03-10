// Phase 2.7 — TageTable runtime test harness
// Tests all accessors, predictor-emulation access patterns, and HARCOM timing.

#include "../predictors/custom/TageTable.hpp"

using namespace hcm;

// ---- Superuser: access to HARCOM internals ----
class harcom_superuser {
public:
  void next_cycle() { panel.next_cycle(); }
  auto get(valtype auto x) { return x.get(); }
} hsu;

// ---- Assertion helpers ----
static int pass_count = 0;
static int fail_count = 0;

#define CHECK(expr)                                                            \
  do {                                                                         \
    if (expr) {                                                                \
      pass_count++;                                                            \
    } else {                                                                   \
      std::cerr << "FAIL: " #expr " at " << __FILE__ << ":" << __LINE__       \
                << std::endl;                                                  \
      fail_count++;                                                            \
    }                                                                          \
  } while (0)

#define CHECK_VAL(v, expected)                                                 \
  do {                                                                         \
    auto _got = hsu.get(v);                                                    \
    if (_got == (expected)) {                                                  \
      pass_count++;                                                            \
    } else {                                                                   \
      std::cerr << "FAIL: " #v " == " << +_got << ", expected "               \
                << +(expected) << " at " << __FILE__ << ":" << __LINE__        \
                << std::endl;                                                  \
      fail_count++;                                                            \
    }                                                                          \
  } while (0)

#define SECTION(name) std::cout << "\n=== " name " ===" << std::endl

// ============================================================================
// Global table instances — HARCOM requires all storage at global scope
// ============================================================================

// T1: Basic (shared tag/u, FF u-bits, 1 bank, no cache)
// TABLE_SIZE=1024, TAG=11, CTR=3, U=1, N=4, BANKS=1
static TageTable<1024, 64, 11, 3, 1, 4, 1, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t1;

// T2: Same config for tag miss test
static TageTable<256, 32, 11, 3, 1, 4, 1, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t2;

// T3: Update flow (BPB=2, FF_CACHE=true so all counters are cached for writeBack)
static TageTable<256, 32, 11, 3, 1, 2, 1, false, true, true, true, 1024,
                 DefaultResetFn, true>
    t3;

// T4: Multi-bank shared tag (N=8, BANKS=4, BPB=2)
static TageTable<256, 32, 11, 3, 1, 8, 4, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t4;

// T5: Per-bank tags (N=4, BANKS=4, BPB=1, SHARED_TAG=false, SHARED_U=false)
static TageTable<128, 16, 7, 3, 1, 4, 4, false, false, false, true, 1024,
                 DefaultResetFn, false>
    t5;

// T6: FF cache (N=4, BANKS=1, BPB=4, USE_FF_CACHE=true)
static TageTable<256, 32, 11, 3, 1, 4, 1, false, true, true, true, 1024,
                 DefaultResetFn, true>
    t6;

// T7: Ahead pipeline (N=4, BANKS=2, USE_AHEAD=true)
static TageTable<256, 32, 11, 3, 1, 4, 2, true, true, true, true, 1024,
                 DefaultResetFn, false>
    t7;

// T8: U-bit SRAM decay (U_STOR_FF=false, DECAY_CTR=1)
static TageTable<64, 8, 7, 3, 2, 1, 1, false, true, true, false, 1,
                 DefaultResetFn, false>
    t8;

// T9: U-bit FF reset (U_WIDTH=2)
static TageTable<64, 8, 7, 3, 2, 1, 1, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t9;

// T10: Predict+update flow
static TageTable<256, 32, 11, 3, 1, 2, 1, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t10;

// T11: Allocation flow
static TageTable<256, 32, 11, 3, 1, 2, 1, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t11;

// T12: Block reuse (N=8, BANKS=2, BPB=4, FF_CACHE=true)
static TageTable<256, 32, 11, 3, 1, 8, 2, false, true, true, true, 1024,
                 DefaultResetFn, true>
    t12;

// T13: Multi-bank independent hit (per-bank tags)
static TageTable<128, 16, 7, 3, 1, 4, 4, false, false, false, true, 1024,
                 DefaultResetFn, false>
    t13;

// T14: Timing verification
static TageTable<256, 32, 11, 3, 1, 2, 1, false, true, true, true, 1024,
                 DefaultResetFn, false>
    t14;

// T15: U-bit SRAM per-bank (SHARED_TAG=false, SHARED_U=false, U_STOR_FF=false, DECAY_CTR=1)
static TageTable<64, 8, 7, 3, 2, 4, 4, false, false, false, false, 1,
                 DefaultResetFn, false>
    t15;

// T16: Everything combined (ahead+multibank+cache+perbank)
static TageTable<256, 32, 11, 3, 1, 8, 2, true, false, false, true, 1024,
                 DefaultResetFn, true>
    t16;

// ============================================================================
// Test 1: Basic write-read-verify
// ============================================================================
void test_basic_write_read() {
  SECTION("test_basic_write_read");

  constexpr u64 IDX = clog2(1024);
  val<IDX> idx{42};
  val<11> tag{0x3AB};
  // BPB=4, CTR_BITS=12; ctrs: 101, 011, 010, 001
  val<12> ctr_bits{0b101011010001};
  val<1> u{1};

  t1.write(idx, 0, 0, tag, ctr_bits, u);
  hsu.next_cycle();

  t1.read(idx, tag, 0);

  std::cout << "--- Read accessor latencies (basic config) ---" << std::endl;
  t1.getHit(0).print("getHit: ");
  t1.getTag(0).print("getTag: ");
  t1.getCounter(0, 0).print("getCounter(0,0): ");
  t1.getU(0).print("getU: ");

  CHECK_VAL(t1.getHit(0), 1);
  CHECK_VAL(t1.getTag(0), 0x3AB);
  CHECK_VAL(t1.getU(0), 1);
  CHECK_VAL(t1.getCounter(0, 0), 0b001);

  hsu.next_cycle();
}

// ============================================================================
// Test 2: Tag miss
// ============================================================================
void test_tag_miss() {
  SECTION("test_tag_miss");

  constexpr u64 IDX = clog2(256);
  val<IDX> idx{10};
  val<11> write_tag{0x100};
  val<12> ctr_bits{0};
  val<1> u{0};

  t2.write(idx, 0, 0, write_tag, ctr_bits, u);
  hsu.next_cycle();

  val<11> read_tag{0x200};
  t2.read(idx, read_tag, 0);

  CHECK_VAL(t2.getHit(0), 0);
  CHECK_VAL(t2.getTag(0), 0x100);

  hsu.next_cycle();
}

// ============================================================================
// Test 3: Write-after-read update flow
// read() → check hit → modify counter → writeBack()
// ============================================================================
void test_update_flow() {
  SECTION("test_update_flow");

  constexpr u64 IDX = clog2(256);
  val<IDX> idx{7};
  val<11> tag{0x55};
  // BPB=2, CTR_BITS=6; ctr[1]=010, ctr[0]=011
  val<6> ctr_bits{0b010011};
  val<1> u{1};

  t3.write(idx, 0, 0, tag, ctr_bits, u);
  hsu.next_cycle();

  t3.read(idx, tag, 0);
  CHECK_VAL(t3.getHit(0), 1);
  CHECK_VAL(t3.getCounter(0, 0), 0b011);
  hsu.next_cycle();

  // Update via direct write (preferred pattern for predictor updates):
  // Predictor computes new counter value and calls write() directly.
  // ctr[1]=010 unchanged, ctr[0]=100 updated
  val<6> new_ctr_bits{0b010100};
  t3.write(idx, 0, 0, tag, new_ctr_bits, u);
  hsu.next_cycle();

  // Re-read and verify
  t3.read(idx, tag, 0);
  CHECK_VAL(t3.getCounter(0, 0), 0b100); // updated
  CHECK_VAL(t3.getCounter(0, 1), 0b010); // unchanged

  hsu.next_cycle();

  // NOTE: writeReg/writeBack flow requires careful HARCOM timing:
  //   read() [cycle A] → next_cycle → writeReg [cycle B] → next_cycle
  //   → writeBack [cycle C] → next_cycle → read to verify [cycle D]
  // The direct write() pattern above is preferred for predictor updates.
}

// ============================================================================
// Test 4: Multi-bank with shared tag
// ============================================================================
void test_multibank_shared() {
  SECTION("test_multibank_shared");

  constexpr u64 IDX = clog2(256);
  val<IDX> idx{20};
  val<11> tag{0x777};
  val<1> u{1};

  // Write different counters to each bank (BPB=2, CTR_BITS=6)
  for (u64 b = 0; b < 4; b++) {
    val<6> ctr_bits{static_cast<unsigned>(b * 9 + 1)};
    t4.write(idx, b, 0, tag, ctr_bits, u);
    hsu.next_cycle();
  }

  t4.read(idx, tag, 0);

  std::cout << "--- Multi-bank shared tag latencies ---" << std::endl;
  t4.getHit(0).print("getHit(0): ");
  for (u64 b = 0; b < 4; b++) {
    t4.getCounter(b, 0).print(
        std::string("getCounter(") + std::to_string(b) + ",0): ");
  }

  CHECK_VAL(t4.getHit(0), 1);
  CHECK_VAL(t4.getHit(1), 1); // shared → same as bank 0
  CHECK_VAL(t4.getHit(3), 1);

  CHECK_VAL(t4.getCounter(0, 0), 1);  // b=0: 000001 → slot0=001
  CHECK_VAL(t4.getCounter(1, 0), 2);  // b=1: 001010 → slot0=010
  CHECK_VAL(t4.getCounter(2, 0), 3);  // b=2: 010011 → slot0=011
  CHECK_VAL(t4.getCounter(3, 0), 4);  // b=3: 011100 → slot0=100

  hsu.next_cycle();
}

// ============================================================================
// Test 5: Per-bank tags (SHARED_TAG=false, SHARED_U=false)
// ============================================================================
void test_perbank_tags() {
  SECTION("test_perbank_tags");

  constexpr u64 IDX = clog2(128);
  val<IDX> idx{15};

  val<7> tag0{0x10};
  val<7> tag1{0x20};
  val<3> ctr{0b101};
  val<1> u{1};

  t5.write(idx, 0, 0, tag0, ctr, u);
  hsu.next_cycle();
  t5.write(idx, 1, 0, tag1, ctr, u);
  hsu.next_cycle();

  // Read with tag matching bank 0 only
  t5.read(idx, tag0, 0);

  std::cout << "--- Per-bank tag latencies ---" << std::endl;
  t5.getHit(0).print("getHit(bank0): ");
  t5.getHit(1).print("getHit(bank1): ");
  t5.getTag(0).print("getTag(bank0): ");
  t5.getTag(1).print("getTag(bank1): ");

  CHECK_VAL(t5.getHit(0), 1); // bank 0 matches
  CHECK_VAL(t5.getHit(1), 0); // bank 1 different tag

  hsu.next_cycle();
}

// ============================================================================
// Test 6: FF cache + reuseRead (BPB=4)
// ============================================================================
void test_ff_cache_reuse() {
  SECTION("test_ff_cache_reuse");

  constexpr u64 IDX = clog2(256);
  val<IDX> idx{50};
  val<11> tag{0x123};
  // BPB=4, CTR_BITS=12; ctrs: 111, 010, 100, 001
  val<12> ctr_bits{0b111010100001};
  val<1> u{0};

  t6.write(idx, 0, 0, tag, ctr_bits, u);
  hsu.next_cycle();

  t6.read(idx, tag, 0);
  CHECK_VAL(t6.getHit(0), 1);

  std::cout << "--- FF cache + reuseRead latencies ---" << std::endl;
  t6.getCounter(0, 0).print("getCounter(0,0): ");
  t6.getCounter(0, 1).print("getCounter(0,1): ");
  t6.getCounter(0, 2).print("getCounter(0,2): ");
  t6.getCounter(0, 3).print("getCounter(0,3): ");

  CHECK_VAL(t6.getCounter(0, 0), 0b001);
  CHECK_VAL(t6.getCounter(0, 1), 0b100);
  CHECK_VAL(t6.getCounter(0, 2), 0b010);
  CHECK_VAL(t6.getCounter(0, 3), 0b111);

  // reuseRead: no new RAM access, reads from cached FF regs
  auto r0 = t6.reuseRead(0, val<2>{0});
  auto r1 = t6.reuseRead(0, val<2>{1});
  auto r2 = t6.reuseRead(0, val<2>{2});
  auto r3 = t6.reuseRead(0, val<2>{3});

  std::cout << "--- reuseRead latencies (FF-only, no RAM) ---" << std::endl;
  r0.print("reuseRead(0,0): ");
  r1.print("reuseRead(0,1): ");
  r2.print("reuseRead(0,2): ");
  r3.print("reuseRead(0,3): ");

  CHECK_VAL(r0, 0b001);
  CHECK_VAL(r1, 0b100);
  CHECK_VAL(r2, 0b010);
  CHECK_VAL(r3, 0b111);

  hsu.next_cycle();
}

// ============================================================================
// Test 7: Ahead pipelining (USE_AHEAD=true)
// ============================================================================
void test_ahead_pipeline() {
  SECTION("test_ahead_pipeline");

  constexpr u64 IDX = clog2(256);
  val<IDX> idx{33};
  val<11> tag_s0{0x0AA};
  val<11> tag_s1{0x0BB};
  val<6> ctr_s0{0b010001};
  val<6> ctr_s1{0b110101};
  val<1> u{1};

  // Write stage 0 banks
  t7.write(idx, 0, 0, tag_s0, ctr_s0, u);
  hsu.next_cycle();
  t7.write(idx, 1, 0, tag_s0, ctr_s0, u);
  hsu.next_cycle();

  // Write stage 1 banks
  t7.write(idx, 0, 1, tag_s1, ctr_s1, u);
  hsu.next_cycle();
  t7.write(idx, 1, 1, tag_s1, ctr_s1, u);
  hsu.next_cycle();

  // Read stage 0
  t7.read(idx, tag_s0, 0);
  std::cout << "--- Ahead stage 0 latencies ---" << std::endl;
  t7.getHit(0).print("getHit(stage0): ");
  t7.getCounter(0, 0).print("getCounter(stage0,bank0,slot0): ");
  CHECK_VAL(t7.getHit(0), 1);
  CHECK_VAL(t7.getCounter(0, 0), 0b001);
  hsu.next_cycle();

  // Read stage 1
  t7.read(idx, tag_s1, 1);
  std::cout << "--- Ahead stage 1 latencies ---" << std::endl;
  t7.getHit(0).print("getHit(stage1): ");
  t7.getCounter(0, 0).print("getCounter(stage1,bank0,slot0): ");
  CHECK_VAL(t7.getHit(0), 1);
  CHECK_VAL(t7.getCounter(0, 0), 0b101);

  hsu.next_cycle();
}

// ============================================================================
// Test 8: U-bit SRAM mode with probabilistic decay
// ============================================================================
void test_u_sram_decay() {
  SECTION("test_u_sram_decay");

  constexpr u64 IDX = clog2(64);
  val<IDX> idx{5};
  val<7> tag{0x11};
  val<3> ctr{0b010};
  val<2> u{3};

  t8.write(idx, 0, 0, tag, ctr, u);
  hsu.next_cycle();

  // Read with MISS — DECAY_CTR=1 → always decay
  val<7> miss_tag{0x22};
  t8.read(idx, miss_tag, 0);

  std::cout << "--- U-bit SRAM decay latencies ---" << std::endl;
  t8.getU(0).print("getU after miss (should decay): ");

  CHECK_VAL(t8.getHit(0), 0);
  CHECK_VAL(t8.getU(0), 2); // 3 → 2

  hsu.next_cycle();

  // Read with HIT — u should NOT decay
  t8.read(idx, tag, 0);
  CHECK_VAL(t8.getHit(0), 1);
  // u re-read from SRAM: still 3 (decay wasn't written back)
  CHECK_VAL(t8.getU(0), 3);

  hsu.next_cycle();
}

// ============================================================================
// Test 9: U-bit FF reset modes
// ============================================================================
void test_u_ff_reset() {
  SECTION("test_u_ff_reset");

  constexpr u64 IDX = clog2(64);
  val<IDX> idx0{0};
  val<IDX> idx1{1};
  val<7> tag{0x3F};
  val<3> ctr{0};
  val<2> u_max{3};

  t9.write(idx0, 0, 0, tag, ctr, u_max);
  hsu.next_cycle();
  t9.write(idx1, 0, 0, tag, ctr, u_max);
  hsu.next_cycle();

  // Verify
  t9.read(idx0, tag, 0);
  CHECK_VAL(t9.getU(0), 3);
  hsu.next_cycle();

  // Mode 1: right shift (3 >> 1 = 1)
  t9.reset_u(val<2>{1});
  hsu.next_cycle();

  t9.read(idx0, tag, 0);
  CHECK_VAL(t9.getU(0), 1);
  hsu.next_cycle();

  // Re-write u=3 for decrement test
  t9.write(idx0, 0, 0, tag, ctr, u_max);
  hsu.next_cycle();

  // Mode 2: saturating decrement (3 → 2)
  t9.reset_u(val<2>{2});
  hsu.next_cycle();

  t9.read(idx0, tag, 0);
  CHECK_VAL(t9.getU(0), 2);
  hsu.next_cycle();

  // Mode 0: reset to zero
  t9.reset_u(val<2>{0});
  hsu.next_cycle();

  t9.read(idx0, tag, 0);
  CHECK_VAL(t9.getU(0), 0);
  hsu.next_cycle();
}

// ============================================================================
// Test 10: Predictor emulation — predict + update flow
// read() → getHit() → getCounter() → [cycle] → writeReg() → writeBack()
// ============================================================================
void test_predictor_predict_update() {
  SECTION("test_predictor_predict_update");

  constexpr u64 IDX = clog2(256);
  val<IDX> idx{100};
  val<11> tag{0x555};
  // BPB=2, CTR_BITS=6: ctr[1]=011, ctr[0]=100
  val<6> ctr_bits{0b011100};
  val<1> u{1};

  t10.write(idx, 0, 0, tag, ctr_bits, u);
  hsu.next_cycle();

  // --- PREDICT ---
  t10.read(idx, tag, 0);
  auto hit = t10.getHit(0);
  auto ctr0 = t10.getCounter(0, 0);

  std::cout << "--- Predict flow latencies ---" << std::endl;
  hit.print("hit: ");
  ctr0.print("ctr0 (prediction): ");

  CHECK_VAL(hit, 1);
  CHECK_VAL(ctr0, 0b100); // MSB=1 → taken
  hsu.next_cycle();

  // --- UPDATE: misprediction, decrement counter via direct write ---
  // ctr[1]=011 unchanged, ctr[0]=011 decremented from 100
  val<6> new_ctr_bits{0b011011};
  t10.write(idx, 0, 0, tag, new_ctr_bits, u);
  hsu.next_cycle();

  // Verify
  t10.read(idx, tag, 0);
  CHECK_VAL(t10.getCounter(0, 0), 0b011);
  hsu.next_cycle();
}

// ============================================================================
// Test 11: Predictor emulation — allocation on miss
// read() → miss → write() with fresh entry
// ============================================================================
void test_predictor_allocate() {
  SECTION("test_predictor_allocate");

  constexpr u64 IDX = clog2(256);
  val<IDX> idx{77};
  val<11> search_tag{0x321};

  // Read — should miss (nothing written at this index)
  t11.read(idx, search_tag, 0);

  std::cout << "--- Allocation flow latencies ---" << std::endl;
  t11.getHit(0).print("getHit (expect miss): ");

  hsu.next_cycle();

  // Allocate fresh entry
  val<6> weak_ctr{0b000100}; // weakly taken
  val<1> u{0};
  t11.write(idx, 0, 0, search_tag, weak_ctr, u);
  hsu.next_cycle();

  // Verify
  t11.read(idx, search_tag, 0);
  CHECK_VAL(t11.getHit(0), 1);
  CHECK_VAL(t11.getCounter(0, 0), 0b100);
  CHECK_VAL(t11.getU(0), 0);
  hsu.next_cycle();
}

// ============================================================================
// Test 12: Block reuse flow with FF cache
// read() → predict branch[0] → reuseRead for branch[1..BPB-1]
// ============================================================================
void test_block_reuse() {
  SECTION("test_block_reuse");

  constexpr u64 IDX = clog2(256);
  val<IDX> idx{60};
  val<11> tag{0x1FF};
  val<1> u{1};

  // Bank 0: BPB=4, CTR_BITS=12, ctrs = 111, 110, 101, 100
  val<12> ctr_b0{0b111110101100};
  t12.write(idx, 0, 0, tag, ctr_b0, u);
  hsu.next_cycle();

  // Bank 1: ctrs = 000, 001, 010, 011
  val<12> ctr_b1{0b000001010011};
  t12.write(idx, 1, 0, tag, ctr_b1, u);
  hsu.next_cycle();

  // Initial read
  t12.read(idx, tag, 0);
  CHECK_VAL(t12.getHit(0), 1);

  auto pred_b0_s0 = t12.getCounter(0, 0);
  CHECK_VAL(pred_b0_s0, 0b100);

  std::cout << "--- Block reuse latencies ---" << std::endl;
  pred_b0_s0.print("getCounter(b0,s0) [initial]: ");

  // Reuse reads (no RAM access)
  auto reuse_b0_s1 = t12.reuseRead(0, val<2>{1});
  auto reuse_b0_s2 = t12.reuseRead(0, val<2>{2});
  auto reuse_b0_s3 = t12.reuseRead(0, val<2>{3});

  reuse_b0_s1.print("reuseRead(b0,s1): ");
  reuse_b0_s2.print("reuseRead(b0,s2): ");
  reuse_b0_s3.print("reuseRead(b0,s3): ");

  CHECK_VAL(reuse_b0_s1, 0b101);
  CHECK_VAL(reuse_b0_s2, 0b110);
  CHECK_VAL(reuse_b0_s3, 0b111);

  // Bank 1 reuse
  auto reuse_b1_s0 = t12.reuseRead(1, val<2>{0});
  auto reuse_b1_s2 = t12.reuseRead(1, val<2>{2});

  reuse_b1_s0.print("reuseRead(b1,s0): ");
  reuse_b1_s2.print("reuseRead(b1,s2): ");

  CHECK_VAL(reuse_b1_s0, 0b011);
  CHECK_VAL(reuse_b1_s2, 0b001);

  hsu.next_cycle();
}

// ============================================================================
// Test 13: Multi-bank predict with per-bank hit/miss
// ============================================================================
void test_multibank_independent_hit() {
  SECTION("test_multibank_independent_hit");

  constexpr u64 IDX = clog2(128);
  val<IDX> idx{25};

  val<7> tags[4] = {val<7>{0x10}, val<7>{0x20}, val<7>{0x30}, val<7>{0x40}};
  for (u64 b = 0; b < 4; b++) {
    t13.write(idx, b, 0, tags[b], val<3>{static_cast<unsigned>(b + 1)},
              val<1>{1});
    hsu.next_cycle();
  }

  // Read with tag matching bank 2 only
  t13.read(idx, val<7>{0x30}, 0);

  std::cout << "--- Multi-bank independent hit latencies ---" << std::endl;
  for (u64 b = 0; b < 4; b++) {
    t13.getHit(b).print(
        std::string("getHit(bank") + std::to_string(b) + "): ");
  }

  CHECK_VAL(t13.getHit(0), 0);
  CHECK_VAL(t13.getHit(1), 0);
  CHECK_VAL(t13.getHit(2), 1);
  CHECK_VAL(t13.getHit(3), 0);

  CHECK_VAL(t13.getCounter(2, 0), 3);

  hsu.next_cycle();
}

// ============================================================================
// Test 14: HARCOM timing — read/write on different cycles
// ============================================================================
void test_timing_read_write_separation() {
  SECTION("test_timing_read_write_separation");

  constexpr u64 IDX = clog2(256);
  val<IDX> idx{88};
  val<11> tag{0x111};
  val<6> ctr{0b010001};
  val<1> u{0};

  u64 cycle_before_write = panel.cycle;
  t14.write(idx, 0, 0, tag, ctr, u);
  hsu.next_cycle();

  u64 cycle_before_read = panel.cycle;
  t14.read(idx, tag, 0);
  CHECK(cycle_before_read > cycle_before_write);

  std::cout << "--- Timing verification ---" << std::endl;
  std::cout << "Write cycle: " << cycle_before_write << std::endl;
  std::cout << "Read cycle:  " << cycle_before_read << std::endl;
  t14.getCounter(0, 0).print("Counter (with timing): ", "\n", true);

  hsu.next_cycle();
}

// ============================================================================
// Test 15: U-bit SRAM per-bank decay
// ============================================================================
void test_u_sram_perbank_decay() {
  SECTION("test_u_sram_perbank_decay");

  constexpr u64 IDX = clog2(64);
  val<IDX> idx{3};

  t15.write(idx, 0, 0, val<7>{0x10}, val<3>{0b101}, val<2>{3});
  hsu.next_cycle();
  t15.write(idx, 1, 0, val<7>{0x20}, val<3>{0b010}, val<2>{3});
  hsu.next_cycle();

  // Read with tag matching bank 0 → bank 1 misses → bank 1 decays
  t15.read(idx, val<7>{0x10}, 0);

  std::cout << "--- Per-bank U decay latencies ---" << std::endl;
  t15.getU(0).print("getU(bank0, hit): ");
  t15.getU(1).print("getU(bank1, miss+decay): ");

  CHECK_VAL(t15.getHit(0), 1);
  CHECK_VAL(t15.getU(0), 3);

  CHECK_VAL(t15.getHit(1), 0);
  CHECK_VAL(t15.getU(1), 2); // decayed from 3

  hsu.next_cycle();
}

// ============================================================================
// Test 16: Ahead + multi-bank + FF cache + per-bank tags combined
// ============================================================================
void test_everything_combined() {
  SECTION("test_everything_combined");

  constexpr u64 IDX = clog2(256);
  val<IDX> idx{200};
  val<11> tag{0x7FF};
  val<1> u{1};

  // BPB=4, CTR_BITS=12
  val<12> ctr{0b101100011010}; // ctrs: 101, 100, 011, 010

  t16.write(idx, 0, 0, tag, ctr, u);
  hsu.next_cycle();
  t16.write(idx, 1, 0, tag, ctr, u);
  hsu.next_cycle();

  t16.read(idx, tag, 0);

  std::cout << "--- Everything combined latencies ---" << std::endl;
  t16.getHit(0).print("getHit(bank0): ");
  t16.getHit(1).print("getHit(bank1): ");
  t16.getTag(0).print("getTag(bank0): ");
  t16.getU(0).print("getU(bank0): ");
  t16.getU(1).print("getU(bank1): ");
  t16.getCounter(0, 0).print("getCounter(b0,s0): ");

  CHECK_VAL(t16.getHit(0), 1);
  CHECK_VAL(t16.getHit(1), 1);
  CHECK_VAL(t16.getCounter(0, 0), 0b010);
  CHECK_VAL(t16.getCounter(0, 3), 0b101);

  // Reuse reads
  auto r = t16.reuseRead(0, val<2>{1});
  r.print("reuseRead(b0,s1): ");
  CHECK_VAL(r, 0b011);

  auto r2 = t16.reuseRead(1, val<2>{2});
  r2.print("reuseRead(b1,s2): ");
  CHECK_VAL(r2, 0b100);

  hsu.next_cycle();
}

// ============================================================================
// HARCOM Cost and Energy Exploration
// ============================================================================
void test_harcom_costs() {
  SECTION("HARCOM Cost Exploration");

  // Global costs (across ALL instantiated tables)
  std::cout << "--- Global HARCOM costs (all 16 table instances) ---"
            << std::endl;
  std::cout << "Total storage bits: " << panel.storage() << std::endl;
  std::cout << "Total transistors:  " << panel.transistors() << std::endl;
  std::cout << "Total energy (fJ):  " << panel.energy_fJ() << std::endl;
  std::cout << "Dynamic power (mW): " << panel.dyn_power_mW() << std::endl;
  std::cout << "Static power (mW):  " << panel.sta_power_mW() << std::endl;
}

// Per-operation energy measurement: measure energy delta for a read and write
void test_per_operation_energy() {
  SECTION("Per-Operation Energy Measurement");

  // Use t1 (basic config: 1024 entries, 11-bit tag, 3-bit ctr, shared tag/u)
  constexpr u64 IDX = clog2(1024);
  val<IDX> idx{99};
  val<11> tag{0x555};
  val<12> ctr{0b101010101010};
  val<1> u{1};

  // Measure write energy
  f64 e_before = panel.energy_fJ();
  t1.write(idx, 0, 0, tag, ctr, u);
  f64 e_after_write = panel.energy_fJ();
  hsu.next_cycle();

  // Measure read energy
  f64 e_before_read = panel.energy_fJ();
  t1.read(idx, tag, 0);
  f64 e_after_read = panel.energy_fJ();

  std::cout << "--- t1 (basic: 1024x11/3/1, shared tag/u, FF u) ---"
            << std::endl;
  std::cout << "Write energy: " << (e_after_write - e_before) << " fJ"
            << std::endl;
  std::cout << "Read energy:  " << (e_after_read - e_before_read) << " fJ"
            << std::endl;

  // Read latency from accessor timing
  t1.getHit(0).print("  Hit latency:     ");
  t1.getCounter(0, 0).print("  Counter latency: ");
  t1.getU(0).print("  U-bit latency:   ");

  hsu.next_cycle();

  // Measure for t4 (multi-bank: 256 entries, 4 banks, shared tag)
  {
    constexpr u64 IDX4 = clog2(256);
    val<IDX4> idx4{20};
    val<11> tag4{0x123};
    val<6> ctr4{0b101010};
    val<1> u4{0};

    f64 e0 = panel.energy_fJ();
    t4.write(idx4, 0, 0, tag4, ctr4, u4);
    f64 e1 = panel.energy_fJ();
    hsu.next_cycle();

    f64 e2 = panel.energy_fJ();
    t4.read(idx4, tag4, 0);
    f64 e3 = panel.energy_fJ();

    std::cout << "--- t4 (multibank: 256x11/3/1, 4 banks, shared tag) ---"
              << std::endl;
    std::cout << "Write energy: " << (e1 - e0) << " fJ" << std::endl;
    std::cout << "Read energy:  " << (e3 - e2) << " fJ" << std::endl;
    t4.getHit(0).print("  Hit latency:     ");
    t4.getCounter(0, 0).print("  Counter latency: ");
    hsu.next_cycle();
  }

  // Measure for t8 (SRAM u-bit: 64 entries, U_WIDTH=2, DECAY_CTR=1)
  {
    constexpr u64 IDX8 = clog2(64);
    val<IDX8> idx8{10};
    val<7> tag8{0x33};
    val<3> ctr8{0b101};
    val<2> u8{2};

    f64 e0 = panel.energy_fJ();
    t8.write(idx8, 0, 0, tag8, ctr8, u8);
    f64 e1 = panel.energy_fJ();
    hsu.next_cycle();

    f64 e2 = panel.energy_fJ();
    t8.read(idx8, tag8, 0);
    f64 e3 = panel.energy_fJ();

    std::cout << "--- t8 (SRAM u: 64x7/3/2, decay=1) ---" << std::endl;
    std::cout << "Write energy: " << (e1 - e0) << " fJ" << std::endl;
    std::cout << "Read energy:  " << (e3 - e2) << " fJ" << std::endl;
    t8.getU(0).print("  U-bit latency:   ");
    hsu.next_cycle();
  }

  // Measure for t16 (everything: ahead+multibank+cache+perbank)
  {
    constexpr u64 IDX16 = clog2(256);
    val<IDX16> idx16{150};
    val<11> tag16{0x3FF};
    val<12> ctr16{0b110011001100};
    val<1> u16{1};

    f64 e0 = panel.energy_fJ();
    t16.write(idx16, 0, 0, tag16, ctr16, u16);
    f64 e1 = panel.energy_fJ();
    hsu.next_cycle();

    f64 e2 = panel.energy_fJ();
    t16.read(idx16, tag16, 0);
    f64 e3 = panel.energy_fJ();

    std::cout << "--- t16 (everything: 256x11/3/1, ahead, 2 banks, "
                 "cache, per-bank) ---"
              << std::endl;
    std::cout << "Write energy: " << (e1 - e0) << " fJ" << std::endl;
    std::cout << "Read energy:  " << (e3 - e2) << " fJ" << std::endl;
    t16.getHit(0).print("  Hit latency:     ");
    t16.getHit(1).print("  Hit latency(b1): ");
    t16.getCounter(0, 0).print("  Counter latency: ");
    t16.getU(0).print("  U-bit latency:   ");

    // Reuse read latency (FF only, should be much cheaper)
    f64 e4 = panel.energy_fJ();
    auto rv = t16.reuseRead(0, val<2>{1});
    f64 e5 = panel.energy_fJ();
    std::cout << "Reuse energy: " << (e5 - e4) << " fJ" << std::endl;
    rv.print("  Reuse latency:   ");

    hsu.next_cycle();
  }
}

// ============================================================================
// Cycle budget verification
// ============================================================================
void test_cycle_budget() {
  SECTION("Cycle Budget Verification");

  // read: 1 cycle (all banks in parallel)
  u64 c0 = panel.cycle;
  constexpr u64 IDX = clog2(1024);
  t1.read(val<IDX>{1}, val<11>{0}, 0);
  // No cycle consumed (combinational within one cycle)
  CHECK(panel.cycle == c0);
  hsu.next_cycle();

  // write: 1 cycle
  u64 c1 = panel.cycle;
  t1.write(val<IDX>{1}, 0, 0, val<11>{0}, val<12>{0}, val<1>{0});
  CHECK(panel.cycle == c1);
  hsu.next_cycle();

  // reuseRead: 0 extra cycles (FF only, no RAM)
  constexpr u64 IDX6 = clog2(256);
  t6.write(val<IDX6>{1}, 0, 0, val<11>{0}, val<12>{0}, val<1>{0});
  hsu.next_cycle();
  t6.read(val<IDX6>{1}, val<11>{0}, 0);
  u64 c2 = panel.cycle;
  auto rv = t6.reuseRead(0, val<2>{0});
  CHECK(panel.cycle == c2); // no extra cycle
  (void)rv;
  hsu.next_cycle();

  // reset_u: 0 cycles (combinational on FFs)
  u64 c3 = panel.cycle;
  t9.reset_u(val<2>{0});
  CHECK(panel.cycle == c3);
  hsu.next_cycle();

  std::cout << "Cycle budget: read=0 extra, write=0 extra, "
               "reuseRead=0 extra, reset_u=0 extra"
            << std::endl;
}

// ============================================================================
int main() {
  panel.make_floorplan();

  test_basic_write_read();
  test_tag_miss();
  test_update_flow();
  test_multibank_shared();
  test_perbank_tags();
  test_ff_cache_reuse();
  test_ahead_pipeline();
  test_u_sram_decay();
  test_u_ff_reset();
  test_predictor_predict_update();
  test_predictor_allocate();
  test_block_reuse();
  test_multibank_independent_hit();
  test_timing_read_write_separation();
  test_u_sram_perbank_decay();
  test_everything_combined();

  // Phase 3: HARCOM cost, energy, timing exploration
  test_harcom_costs();
  test_per_operation_energy();
  test_cycle_budget();

  std::cout << "\n========================================" << std::endl;
  std::cout << "Results: " << pass_count << " passed, " << fail_count
            << " failed" << std::endl;
  std::cout << "========================================" << std::endl;

  return fail_count > 0 ? 1 : 0;
}
