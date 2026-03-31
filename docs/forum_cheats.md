# CBP-NG Forum Rulings & Tricks

Useful clarifications and "free" techniques confirmed by the organizers on the CBP-NG mailing list.

---

## 1. N-Style Predictors: Branch-Rank Steering Allowed

**Source**: CBP-NG mailing list, 2026-03-30 (Aaron + Pierre)

Using untimed C++ state like `num_branch` (updated in `update_condbr()`) to steer same-block predictions is **intentionally allowed**.

### Background

- The provided predictors (`bimodalN`, `gshareN`, `gshareN_ahead`, `perceptron`) all use C++ `num_branch` / `num_condbr` to track which conditional branch within a block is being predicted. This drives which prediction bit is returned and whether the block continues.
- A contestant (lpha) asked whether this violates NHC since branch-rank info originates from `update_condbr()` (at least one cycle after prediction) yet affects later same-block predictions without HARCOM timing cost.

### Organizer Response

**Aaron:**
> Our rationale for allowing the usage of branch-rank information:
> 1. The interface with the BTB is already intentionally "perfect"
> 2. We want to allow the exploration of predictors which predict based on the Nth branch after the entry point in addition to those which predict the Nth instruction after the entry point
> 3. We did not want to add additional unwanted complexity to the interface
> 4. We believe there are potential predictor/BTB organizations which will not result in significant additional cost for predictors of this type.

**Pierre:**
> The implicit assumption is that there exists an efficient BTB implementation that can work with whatever CBP contestants may come with under the championship's rules. [...] An N-style CBP probably requires an N-style level-1 BTB.

### Implications

- `num_branch` as C++ int in update_condbr / predict paths is valid
- BR_P_ENTRY>1 (multiple branches per entry, indexed by branch rank) is a legitimate direction
- No need to model branch-rank tracking in HARCOM -- cost assumed absorbed by the ideal BTB
- The Oracle BTB is "no assumptions" -- not an idealized traditional BTB

---

## 2. std::rand() for Random Number Generation

**Source**: CBP-NG organizer (Aaron), private correspondence

HARCOM has no built-in random number primitive. The approved approach is to use **`std::rand()`** directly, wrapped in a `val<N>{}` cast. Zero HARCOM cost.

### Organizer Quote

> "Using std::rand as the TAGE example predictor seems reasonable. Even if fully modeled in HARCOM, pseudo-random number generation shouldn't really impact the timing/latency of your predictor since you could pipeline it. Perhaps the main un-modeled cost is energy, and I think it unlikely that pseudo-random number generation would be used so frequently in a competitive predictor that it consumed a significant portion of its energy."
> -- Aaron (CBP-NG Organizer)

### Pattern

```cpp
#include <cstdlib>

// N-bit random value, zero HARCOM cost
val<N> rng{(unsigned)std::rand()};

// Example: probabilistic decay with probability 1/DECAY_CTR
constexpr u64 DECAY_BITS = clog2(DECAY_CTR);
val<DECAY_BITS> rng{(unsigned)std::rand()};
val<1> should_decay = (rng == hard<0>{});  // true ~1/DECAY_CTR of the time
```

### Use Cases

| Use case | Bits needed | Pattern |
|----------|-------------|---------|
| U-bit decay (1/1024) | 10 | `val<10>{(unsigned)std::rand()} == hard<0>{}` |
| U-bit decay (1/256) | 8 | `val<8>{(unsigned)std::rand()} == hard<0>{}` |
| Allocation tie-break | 2 | `select(val<2>{(unsigned)std::rand()} == hard<0>{}, candidate_b, candidate_a)` |

### Notes

- Cast to `unsigned` before wrapping to avoid signed/unsigned warnings
- No need to seed -- deterministic behavior aids debugging
- The reference TAGE predictor (`predictors/tage.hpp`) already uses this pattern
- This supersedes any LFSR-based designs in planning docs -- `std::rand()` is strictly better (zero cost)

---

## 3. VFS Scoring: P1 Latency Double-Count Fix

**Source**: CBP-NG mailing list, 2026-03-15 to 2026-03-19 (Kim Wen, Pierre, Aaron, lpha)

The original VFS formula double-counted P1 latency on P2 mispredictions (once in IPCcbp via `npred`, once in CPIcbp via the misprediction penalty). Fixed in commit [b6f179a](https://github.com/AmpereComputing/cbp-ng/commit/b6f179adb466dd1f0be75c51946c505d8f697877).

### Key Points

- **Mid-block P1/P2 divergence**: penalty = full P2 latency (block re-fetched from P2's target)
- **Block-end P1/P2 divergence**: penalty = P2 latency - P1 latency (would have paid P1 anyway)
- `misprediction_penalty` renamed to `p2_to_exec_stages`, value 8 → 9
- P1 latency now subtracted from CPIcbp misprediction penalty to avoid double-counting
- **Net effect: zero change for any predictor with P1 = 1 cycle** (all competitive predictors)
- Make sure to use the latest `predictor_metrics.py`

---

## 4. TAGE-SC-L Won't Win Unmodified; TAGE+Perceptron is the Target

**Source**: CBP-NG mailing list, 2026-02-24 to 2026-02-26 (Kim Wen, Aaron, Pierre, Daniel Jimenez)

### Organizer Signals

**Aaron** on TAGE-SC-L:
> To quote Andre's CBP2025 TAGE-SC paper, it "replicates most of the features that would prevent any reasonable direct hardware implementation". We believe it unlikely an unmodified re-implementation of CBP2025 TAGE-SC will achieve the highest VFS score for CBP-NG, despite its prediction accuracy.

**Pierre**:
> To obtain a good VFS, high prediction will not suffice, you will need to consider latency, throughput and energy as well.

Pierre pointed to Seznec's paper (https://hal.science/hal-04804900/) for hints on optimizing TAGE for latency/throughput/energy.

### Daniel Jimenez (Perceptron Inventor) on Hybrid Designs

> Two branch prediction techniques have enjoyed broad industry adoption in recent years: TAGE and perceptron. In practice, the highest-accuracy designs combine both. [...] AMD publicly says that they use TAGE and perceptron together.

Jimenez noted the original 2001 perceptron doesn't reflect modern designs. Pointed to **Samsung ISCA 2020** (https://ieeexplore.ieee.org/document/9138988) as a practical hashed multi-table perceptron with geometric history lengths -- much more implementation-friendly.

### Available Resources

- `predictors/perceptron.hpp` -- Aaron's basic HARCOM perceptron example (2001-style, released 2026-02-26)
- Samsung ISCA 2020 paper -- practical perceptron-based design for real hardware
- Seznec hal-04804900 -- hints on TAGE latency/energy optimization

### Implications

- Pure accuracy maximization (TAGE-SC-L style) will lose on VFS -- energy/latency matter
- **Hybrid TAGE + perceptron** is the expected competitive direction (industry standard)
- Our SCOverrider is already a lightweight perceptron-style component -- validates that approach
- The Samsung paper's hashed multi-table organization with geometric history lengths could inform a better SC design

---

## 5. execute_if Multi-Bit Masks Removed; Bit Replication Exploit Fixed

**Source**: CBP-NG mailing list, 2026-02-12 to 2026-03-03 (lpha, Pierre, Aaron)

### The Exploit

lpha discovered that `execute_if` could be abused for zero-latency data transfer. The trick: use a value as the `execute_if` condition to conditionally write registers, then detect what was written -- effectively replicating the condition bit without paying fanout/wiring costs. Combined with `fanout()` after `fo1()` granting free read credits, this could shave ~39ps off gshare latency.

### The Fix (Two Rounds)

**Round 1** (2026-02-17): Made `execute_if` predicate implicitly read by every reg/RAM write inside it. Added section 3.6.2 to the HARCOM manual. ~30% simulation slowdown for TAGE.

**Round 2** (2026-03-03): Round 1 was insufficient -- the predicate was read via a copy, still allowing low-latency replication. Pierre's final fix: **multi-bit masks in `execute_if` are now disabled**. All `execute_if` calls must use single-bit predicates. Multi-bit mask usage must be replaced with a loop of single-bit `execute_if` calls.

Updated HARCOM: https://github.com/PierreMichaud-Inria/cbp-ng

### Implications

- **All `execute_if` calls must use `val<1>` predicates** -- no multi-bit masks
- If you have `execute_if(multi_bit_mask, ...)`, replace with a loop: `for (u64 i = 0; i < N; i++) execute_if(mask_bit_i, ...)`
- The example predictors were updated accordingly
- `fanout()` after `fo1()` no longer grants free reads
- Make sure to use the latest `harcom.hpp`

### Pierre's Design Philosophy

> If you can duplicate a bit at no cost, you can replicate the bit as many times as you want with no cost, and this would break the consistency of the HARCOM model. The philosophy of HARCOM is that inefficient implementations are possible, but implementations with an unrealistically low cost should not be possible without conspicuously violating the rules of the HARCOM language.
