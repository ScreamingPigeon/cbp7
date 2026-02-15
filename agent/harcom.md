# HARCOM Reference for Writing Hardware Component Models

HARCOM (**har**dware **com**plexity model) is a C++20 header-only library (`harcom.hpp`) for estimating the hardware complexity (delay, energy, area, transistors) of microarchitectural components directly inside a performance simulator. It is **not** an HDL. It replaces C++ integers with opaque types that track timing and cost. The full specification is in `docs/harcom.pdf`.

## Core Concept

Replace C++ integers with HARCOM types (`val`, `reg`, `ram`, etc.) in the part of the simulator you want to model. The rest of the simulator stays unchanged. Every operation on HARCOM types accumulates hardware cost (transistors, energy, delay). You cannot read HARCOM values as plain integers — this forces you to write code that maps to realistic hardware.

Namespace: `hcm`. Type aliases: `u64 = uint64_t`, `i64 = int64_t`, `f64 = double`, etc.

---

## Data Types

### `val<N>` and `val<N,T>` — Transient Values

A transient N-bit value with a timing (picoseconds) and a location. `T` is the underlying C++ integer type (default `u64`; use `i64` for signed). Max N = 64.

```cpp
val<8> x = 1;           // 8-bit unsigned, timing=0 (hardwired)
val<6,i64> y = -1;      // 6-bit signed
auto z = x + y;         // 9-bit val, timing = max(timing(x), timing(y)) + adder_delay
```

**Critical rules:**
- A `val` **cannot be modified** after creation. It is immutable.
- A `val`'s actual value is **private**. You cannot read it as a C++ int (compilation error unless `-DCHEATING_MODE`).
- Initializing from a C++ literal or variable creates a zero-timing val (hardwired constant).
- When assigning between different sizes: truncated if shorter, sign-extended if longer.

### `reg<N>` and `reg<N,T>` — Persistent State (Registers)

Derived from `val`. Represents a hardware register (flip-flop). Unlike `val`, a `reg` **can** be modified.

```cpp
reg<4> counter = 0;
counter = counter + 1;  // OK — one write per cycle
```

**Critical rules:**
- **All regs must have the same lifetime** — declare them all at the same scope (typically as struct members). A reg cannot be created after another reg has been destroyed.
- **A reg can be written at most once per clock cycle.** Call `panel.next_cycle()` (superuser only) between writes.
- Uninitialized regs default to zero.
- Declare regs as `static` variables or struct members to ensure correct lifetime.

### `arr<T,N>` — Fixed-Size Array

An array of N valtype objects (vals or regs).

```cpp
arr<val<3>,4> A = {1,3,0,2};
A[2].print();                    // subscript with C++ int
A.select(val<2>{1}).print();     // subscript with val (hardware mux cost)
```

**Key functions (hardware cost in bold):**
- **`select(idx_val)`** — mux read with a val index
- `concat()` — concatenate all elements into one val (treats elements as bit-vector chunks; index 0 = rightmost bits)
- **`fold_xor()`**, **`fold_or()`**, **`fold_and()`**, **`fold_add()`** etc. — reduce operations
- `shift_left(v)`, `shift_right(v)`, `make_array(v)`, `append(v)`, `truncate(hard<K>{})` — bit-vector manipulations
- **`fanout(hard<K>{})`**, **`fo1()`** — set fanout for delay optimization (see Fanout section)

### `hard<N>` — Compile-Time Constant

Represents a value known at design time (a hardware parameter). Zero cost, zero timing.

```cpp
val<8> y = x << hard<4>{};   // shift left by 4 — no hardware cost for the shift amount
auto r = x % hard<4>{};      // modulo requires hard divisor
```

You can substitute a C++ integer for `hard` in many operators, but some operators **require** `hard` (e.g., `%`, `>>` for unsigned right-shift, `&` and `|` with a constant).

### `ram<T,N>` — Random Access Memory (SRAM)

N-entry RAM storing type T (a val type or `arr<val<K>,M>`).

```cpp
ram<val<4>,1024> table;
ram<arr<val<64>,2>,256> wide_mem;
```

**Critical rules:**
- **All RAMs must have the same lifetime as regs.**
- **Only one RAM access (read OR write) per clock cycle** per RAM object. (Use `-DREAD_WRITE_RAM` to allow one read + one write.)
- Must call `panel.make_floorplan()` after declaring all RAMs, before using them.
- `read(addr)` returns data; `write(addr, data)` stores data.
- `reset()` zeros the entire RAM.
- RAM reads return data written by the most recent write whose timing <= the read's address/data timing. You cannot read data "from the future."

### `rom<T,N>` — Read-Only Memory (Lookup Table)

Initialized at creation, read with `operator()`. Modeled as combinational logic (decoder + OR trees).

```cpp
rom<val<3>,16> popcount = [](u64 i){ return std::popcount(i); };
popcount(val<4>{7}).print();  // reads the ROM
```

---

## Operators

All operators take valtypes (`val`/`reg`) as input and produce `val` as output.

| Operator | Notes |
|----------|-------|
| `==`, `!=`, `>`, `<`, `>=`, `<=` | Two vals of same size, or one val + one hard. Output: `val<1>`. |
| `&`, `\|`, `^` | Two vals, or one val + one hard. Output: same width as longest input. **`&` and `\|` with hard have zero cost.** |
| `~` | Bitwise NOT. Same width as input. |
| `<<` | One val + one hard shift count. Output: same width. |
| `>>` | **Unsigned** `>>` with hard: zero cost. Signed or val shift count: has cost. |
| `+`, `-` | Two vals or one val + one hard. Output: one bit longer than longest input. |
| `*` | Two vals or one val + one hard. Output: enough bits for the product (max 64). |
| `/` | Unsigned val dividend, hard divisor. |
| `%` | Val modulo hard divisor. |
| unary `-` | Negate. Same width as input. |

---

## Essential Functions

### `select(cond, x1, x0)` — 2-to-1 Mux

Returns `x1` if `cond` is true, `x0` otherwise. This is the primary way to do conditional logic (instead of `if/else` on opaque values).

```cpp
val<1> take = prediction;
val<64> target = select(take, branch_target, pc_plus_4);
```

### `execute_if(mask, F)` — Conditional Execution

Executes lambda `F(i)` for each bit of `mask` that is set. If `F` returns a val, `execute_if` returns an `arr` (zeros for unset mask bits). Essential for conditional RAM access.

```cpp
val<1> do_write = should_update;
execute_if(do_write, [&](){ table.write(idx, new_val); });
// conditional read:
val<4> result = execute_if(do_read, [&](){ return table.read(idx); });
```

When mask is zero, operations inside F consume no energy, but transistor count is still charged. Each RAM/reg access inside F still counts toward the one-access-per-cycle limit.

### `concat(a, b, ...)` — Bit Concatenation

Concatenates multiple unsigned vals into a single val. Rightmost argument becomes the least-significant bits.

```cpp
val<3> lo = 0b111;
val<4> hi = 0b0011;
val<7> combined = concat(lo, hi);  // 0b0011_111
```

### `split<sizes...>(v)` — Bit Splitting

Splits a val into multiple vals. Returns a structured binding.

```cpp
val<12> x = 0b000111000111;
auto [a, b, c] = split<4,3,5>(x);  // a=4 bits (rightmost), b=3 bits, c=5 bits (leftmost)
```

### `a_plus_bc(a, b, c)` — Fused Multiply-Add

Computes `a + b * c`.

### `fold(arr, op)` / `scan(arr, op)` — Tree Reduction / Prefix Sum

Applies a binary associative operation across an array in a tree structure.

---

## Fanout and `fo1()` — Managing Read Costs

**Every read of a named val/reg incurs delay and energy cost** (models a FO2 inverter chain). Multiple reads accumulate linearly, which can dominate timing.

### `fanout(hard<K>{})` — Declare Expected Reads

Makes the read delay grow **logarithmically** instead of linearly, by declaring the expected number of reads:

```cpp
val<4> x = some_computation;
x.fanout(hard<8>{});          // x will be read up to 8 times
auto a = x + 1;               // each read now has log delay
auto b = x & mask;
```

### `fo1()` — Read-Once Optimization

If a named val will be read exactly once, call `.fo1()` to avoid any read cost. After `fo1()`, reading the original variable again triggers an error.

```cpp
val<4> temp = x + y;
auto result = temp.fo1() & mask;  // temp is consumed; cannot read temp again
```

**Best practice:** Keep transient values unnamed (use `auto` with expressions) when possible. Name a val only when needed for readability, and use `fo1()` or `fanout()` to manage cost.

---

## RAM Usage Patterns

```cpp
// Declare all storage up front
ram<val<4>,1024> table;
ram<arr<val<8>,2>,512> wide_table;
panel.make_floorplan();     // REQUIRED after all RAM declarations

// Write
val<10> addr = compute_index;
val<4> data = compute_value;
table.write(addr, data);

// Must advance cycle before next access to same RAM
panel.next_cycle();         // superuser only

// Read
val<4> result = table.read(addr);

// Conditional access
execute_if(do_access, [&](){
  table.write(addr, data);
});
```

### RAM Arrays

Declare arrays of RAMs using C arrays:

```cpp
ram<val<4>,256> banks[8];
panel.make_floorplan();
```

### `distribute(ram_array)` — Write to Multiple RAMs

Efficiently sends a value to an array of RAM locations:

```cpp
val<64> data = 7;
arr<val<64>,16> d = data.distribute(mem);  // create copies at each RAM location
for (u64 i = 0; i < 16; i++)
  mem[i].write(0, d[i]);
```

### `connect(ram)` — Move Value to a RAM's Location

Forces an operation to happen at a specific RAM's location, reducing wiring cost:

```cpp
val<1> x1 = slow_ram.read(addr);
(x0.fo1() & x1.fo1().connect(fast_ram)).print();  // operation at fast_ram's location
```

---

## The `panel` Global Object

Tracks global hardware costs. Key functions:

| Function | Description |
|----------|-------------|
| `panel.make_floorplan()` | **Required** after declaring all RAMs. Generates spatial layout. |
| `panel.next_cycle()` | Advances clock cycle (**superuser only**). Enables re-writing regs and re-accessing RAMs. |
| `panel.print()` | Prints all costs: storage bits, transistors, SRAM area, energy, power. |
| `panel.energy_fJ()` | Returns total dynamic energy in femtojoules. |
| `panel.transistors()` | Returns total transistor count. |
| `panel.storage()` | Returns total storage in bits. |
| `panel.dyn_power_mW()` | Returns dynamic power in milliwatts. |
| `panel.sta_power_mW()` | Returns static/leakage power in milliwatts. |
| `panel.clock_cycle_ps` | Clock cycle period in picoseconds (default 300, modifiable by superuser). |

---

## The `harcom_superuser` Class

The superuser bridges HARCOM and the rest of the simulator. Defined in the global namespace. **The superuser can:**

- Call `panel.next_cycle()` to advance cycles
- Convert HARCOM values to/from C++ integers via private members: `x.get()`, `x.set_time(ps)`, `x.get_vt()` (returns `[value, timing]`)
- Assign to vals after construction

**The user (predictor code) cannot** access these private members. User code must follow the HARCOM language rules.

---

## The No-Hidden-Cost (NHC) Rule

Every statement's hardware cost must be accounted for. Guidelines:

- All data unknown at compile time **must** be a valtype (`val`/`reg`), not a C++ int.
- Do **not** use non-const C++ integers whose lifetime spans multiple clock cycles — use `reg` instead.
- Do **not** access private class members of HARCOM types.
- No type punning (`reinterpret_cast`, `union`, `memcpy`, `std::bit_cast` on HARCOM types).
- Do not put HARCOM types in `union`s or `std::vector`.
- `std::array`, `std::tuple` are compatible with HARCOM types.
- C++ loops, functions, classes, templates are all fine as long as NHC is not violated.
- Using C++ integers for loop counters and array indices is OK (they have no hardware cost — the loop is "unrolled" at compile time).

```cpp
// OK: loop index is a C++ int, body operates on HARCOM types
for (int i = 0; i < 3; i++)
  B[i] = A[i];  // equivalent to B[0]=A[0]; B[1]=A[1]; B[2]=A[2];
```

---

## Wiring and Regions

### Locations and Wiring Cost

Each RAM has a location in the floorplan. Regs and vals inherit locations from nearby RAMs. Reading a value from a different location incurs wiring delay/energy proportional to Manhattan distance.

### Regions

Group RAMs that interact frequently into the same region. RAMs in the same region are placed close together.

```cpp
region R1;
ram<val<4>,64> tbl_a;
ram<val<4>,64> tbl_b;
region R2;
ram<val<8>,1024> large_tbl;
panel.make_floorplan();
```

### Zones

Zones subdivide regions for finer floorplan control:

```cpp
zone Z1;
ram<val<8>,64> m1[8] {"m1[0]"};
zone Z2;
ram<val<4>,32> m2[8] {"m2[0]"};
panel.make_floorplan();
```

---

## Compiler Options

| Option | Effect |
|--------|--------|
| `-DFREE_FANOUT` | Removes all fanout delays (useful to see upper bound of fanout savings). |
| `-DCHECK_FANOUT` | Error if actual fanout exceeds declared fanout. |
| `-DCHEATING_MODE` | Allows converting val to C++ int (for `assert` debugging). |
| `-DFLOORPLAN` | `make_floorplan()` generates a `floorplan.gv` Graphviz file. |
| `-DFREE_WIRING` | Removes wiring costs (except inside RAMs). |
| `-DREAD_WRITE_RAM` | Allows one read + one write to the same RAM per cycle. |

Compile with `-Wall -Wextra -Werror` and `-std=c++20`.

---

## `static_loop` — Compile-Time Iteration

Iterates a lambda over a compile-time integer range. Useful for template metaprogramming with HARCOM types:

```cpp
static_loop<10>([]<int I>(){
  // I is a compile-time constant
});
```

---

## Printing and Debugging

```cpp
x.print("label=");         // decimal: "label=255 (t=44 ps, loc=0)"
x.printb("label=");        // binary:  "label=11111111 (t=44 ps, loc=0)"
x.print("x=", "\n", false); // no timing info
panel.print();              // global costs summary
```

---

## Common Patterns for Branch Predictors

### Saturating Counter Update
```cpp
// Use select() instead of if/else
reg<3> ctr = 4;
val<1> taken = ...;
val<3> incremented = select(ctr == val<3>::maxval, ctr, ctr + 1);
val<3> decremented = select(ctr == 0, ctr, ctr - 1);
ctr = select(taken, incremented, decremented);
```

### Hash/Index Computation
```cpp
val<10> idx = (pc ^ history).fo1();  // XOR folding, consume with fo1
val<8> tag = split<10,8>(concat(pc, history)).get<1>(); // NOT VALID — use positional split
auto [idx_part, tag_part] = split<10,8>(concat(pc_bits, hist_bits));
```

### Conditional Table Update
```cpp
execute_if(should_update, [&](){
  table.write(index, new_value);
});
```

### Reading Multiple Tables in One Cycle
Each `ram` object can only be accessed once per cycle. Use separate `ram` objects for structures that need concurrent reads:

```cpp
ram<val<4>,1024> table0;
ram<val<4>,1024> table1;
panel.make_floorplan();

// Both reads in the same cycle — OK, different RAM objects
val<4> r0 = table0.read(addr0);
val<4> r1 = table1.read(addr1);
```

---

## Summary of Key Constraints

1. `val` is immutable. `reg` is mutable (once per cycle).
2. All storage (`reg`, `ram`) must have the same lifetime — declare together.
3. One RAM access per cycle per RAM object (unless `-DREAD_WRITE_RAM`).
4. Call `panel.make_floorplan()` after all RAM declarations.
5. Use `select()` for conditional values, `execute_if()` for conditional side effects (especially RAM access).
6. Manage fanout with `fanout()` and `fo1()` to control delay.
7. Follow the NHC rule: all runtime data must be HARCOM types, not bare C++ ints.
8. The superuser handles the HARCOM/simulator boundary; user code stays within HARCOM constraints.
