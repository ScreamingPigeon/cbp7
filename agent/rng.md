# Random Number Generation in the Custom Predictor

## Decision Summary

HARCOM has no built-in random number primitive. After discussion with the CBP7 organizers,
the approved approach is to use **`std::rand()`** directly, wrapped in a `val<N>{}` cast.

> "Using std::rand as the TAGE example predictor seems reasonable. Even if fully modeled in
> HARCOM, pseudo-random number generation shouldn't really impact the timing/latency of your
> predictor since you could pipeline it. Perhaps the main un-modeled cost is energy, and I
> think it unlikely that pseudo-random number generation would be used so frequently in a
> competitive predictor that it consumed a significant portion of its energy."
> — Aaron (CBP-NG Organizer)

This supersedes the earlier TODO in `TageTable.hpp` that planned an LFSR implementation.
LFSR adds hardware cost (flip-flops + XOR gates) that is not worth modeling given the
organizer's position.

---

## Pattern

```cpp
#include <cstdlib>

// N-bit random value, costs nothing in HARCOM's model
val<N> rng{(unsigned)std::rand()};
```

The reference TAGE predictor (`predictors/tage.hpp:380`) already uses this:

```cpp
val<NUMG> collamask12 =
    select(val<2>{std::rand()} == hard<0>{}, collamask2.fo1(), collamask1);
```

This uses 2 random bits to pick between two candidate allocation entries with roughly equal
probability.

---

## Use Cases in Our Predictor

### 1. Probabilistic U-bit Decay (`TageTable`)

The `DECAY_CTR` template parameter (default 1024) controls the decay probability.
On a tag miss during `updateBlock()`, decrement the u-bit with probability `1/DECAY_CTR`:

```cpp
#include <cstdlib>

// Inside TageTable::decrementU()
val<1> decrementU() {
    constexpr u64 DECAY_BITS = clog2(DECAY_CTR);
    val<DECAY_BITS> rng{(unsigned)std::rand()};
    return rng == hard<0>{};  // true with probability 1/DECAY_CTR
}
```

In `updateBlock()`, apply the decay on a tag miss:

```cpp
val<U_WIDTH> cur_u{u_reg};
val<U_WIDTH> dec_u = select(cur_u == hard<0>{}, cur_u, val<U_WIDTH>{cur_u - 1});
val<1> miss = val<1>{hit} == hard<0>{};
val<U_WIDTH> upd_u = select(miss & decrementU(), dec_u, cur_u);
// ... use upd_u instead of u_reg when building reg_entry
```

### 2. Allocation Candidate Selection (`Tage`)

When multiple tables are eligible for allocation on a misprediction, use random bits to
break ties rather than always picking the shortest or longest:

```cpp
// Flip between two candidate masks with ~50% probability
val<NUMG> chosen =
    select(val<2>{(unsigned)std::rand()} == hard<0>{}, candidate_b, candidate_a);
```

---

## Entropy Requirements

From the thread: 4–8 bits is sufficient for branch predictor use cases.

| Use case              | Bits needed | Notes                                  |
|-----------------------|-------------|----------------------------------------|
| Decay gate (1/1024)   | 10          | `clog2(DECAY_CTR)` bits, check == 0   |
| Decay gate (1/256)    | 8           | Smaller DECAY_CTR                     |
| Allocation tie-break  | 2           | Match reference TAGE pattern           |

---

## Notes

- Include `<cstdlib>` explicitly; do not rely on transitive includes from HARCOM headers.
- Cast to `unsigned` before wrapping: `val<N>{(unsigned)std::rand()}` avoids signed/unsigned
  conversion warnings.
- There is no need to seed `std::rand()` — deterministic behavior is fine and aids debugging.
- The LFSR sections in `tagetable.md` and the `TageTable.hpp` TODO comment are now outdated;
  `std::rand()` is the correct approach.
