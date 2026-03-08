# CBP Trace Analyzer Overview

This document explains how `trace_analyzer.cpp` uses `trace_reader.hpp` to analyze CBP trace files and describes all metrics printed by the analyzer.

The analyzer processes instruction traces and extracts statistics about branch behavior, instruction distribution, and control-flow structure.

---

# How the Analyzer Uses `trace_reader.hpp`

CBP traces store instruction-level execution information in a compressed binary format. The provided `trace_reader.hpp` file handles decoding this format.

The analyzer interacts with the trace reader as follows:

1. A `trace_reader` object is created using the trace path.

```
trace_reader reader(trace_path, "trace");
```

2. The analyzer repeatedly requests instructions using:

```
instruction inst = reader.next_instruction();
```

3. Each instruction returned contains:

- `pc` — program counter of the instruction
- `next_pc` — next program counter after execution
- `inst_class` — classification of instruction type
- `branch` — whether the instruction is a branch
- `taken_branch` — whether the branch was taken

4. The analyzer updates statistics based on this information.

5. When the end of the trace is reached, the reader throws an `out_of_instructions` exception which stops the analysis loop.

The analyzer itself does not interpret the binary trace format; it only consumes decoded instructions produced by `trace_reader`.

---

# Window-Based Statistics

The analyzer divides the execution stream into **fixed-size instruction windows** (for example 1M instructions per window). For each window it prints summary statistics.

Example:

```
[Window 0] Branches: 197805 TakenRate: 0.5924 BranchDensity: 0.1978
```

These metrics describe how branch behavior changes over time.

## Branches

Number of branch instructions executed within that window.

## Taken Rate

Fraction of branches that were taken in that window.

```
TakenRate = taken_branches / total_branches
```

This indicates whether branches tend to be taken or not taken during that portion of execution.

## Branch Density

Fraction of instructions that are branches within that window.

```
BranchDensity = branches / total_instructions
```

Higher density means the program is performing frequent control flow changes.

---

# Global Statistics

After processing the entire trace, the analyzer prints aggregate statistics.

Example:

```
==== GLOBAL STATS ====
Total instructions: 3005000
Total branches: 594405
Branch density: 0.1978
Taken rate: 0.6187
Backward rate: 0.3883
Unique branch PCs: 2083
```

## Total Instructions

Total number of instructions executed in the trace.

## Total Branches

Total number of branch instructions encountered.

## Branch Density

Fraction of instructions that are branches.

```
branch_density = total_branches / total_instructions
```

## Taken Rate

Overall probability that a branch is taken.

```
taken_rate = taken_branches / total_branches
```

## Backward Rate

Fraction of branches that jump **backward in the instruction stream**.

```
backward_branch = (next_pc < pc)
```

Backward branches are typically associated with **loops**.

High backward rates often indicate loop-heavy programs.

## Unique Branch PCs

Number of distinct program counter values that correspond to branch instructions.

This represents the number of unique branch sites in the program.

---

# Hot Branch Analysis

The analyzer identifies the most frequently executed branch instructions ("hot branches").

Example:

```
==== TOP 10 HOT BRANCHES ====
PC: 0xf02d88 Count: 35885 TakenRate: 0.8184 Entropy: 0.6835 FlipRate: 0.1428
```

Each entry represents a specific branch instruction address.

## PC

Program counter of the branch instruction.

## Count

Number of times this branch executed.

Hot branches are important because they dominate predictor performance.

## Taken Rate

Probability that this specific branch is taken.

```
taken_rate = taken_count / execution_count
```

Values close to 0 or 1 indicate highly predictable branches.

Values near 0.5 are difficult to predict.

---

# Entropy

Entropy measures the **predictability of branch outcomes**.

It is computed using binary entropy:

```
H = -(p log2 p + (1-p) log2 (1-p))
```

Where:

```
p = taken probability
```

Entropy ranges from:

| Value | Meaning |
|------|-------|
0 | perfectly predictable |
1 | maximally unpredictable |

Examples:

- `TakenRate ≈ 0 or 1` → low entropy
- `TakenRate ≈ 0.5` → high entropy

Entropy helps identify branches that are inherently difficult to predict.

---

# Flip Rate

Flip rate measures **how often branch direction changes between consecutive executions**.

```
FlipRate = number_of_direction_changes / executions
```

Example sequence:

```
T T T N T N N
```

Direction changes:

```
T→N
N→T
T→N
```

---

# Usage and Configuration

## Building

```bash
make trace-analyze      # Build the analyzer
make run-trace-analyze  # Build and run on default trace
```

## Running with Custom Parameters

The analyzer accepts an optional **region shift** parameter to control memory region granularity:

```bash
# Use Makefile with 4KB pages (default, ARM v8)
make run-trace-analyze TRACE=./traces/502-gcc-all_16112_trace.gz REGION_SHIFT=12

# 2MB large pages
make run-trace-analyze TRACE=./traces/502-gcc-all_16112_trace.gz REGION_SHIFT=21

# 1GB sections
make run-trace-analyze TRACE=./traces/502-gcc-all_16112_trace.gz REGION_SHIFT=30

# Direct binary (no output file)
./trace-analyze ./traces/502-gcc-all_16112_trace.gz 12
```

Output is saved to `out/trace_analysis.txt` when using the Makefile.

---

# Per-Window Region Distribution

Each 1M instruction window shows the **top 3 memory regions** (pages) that were active during that window.

Example:

```
[Window 0] Branches: 286905 TakenRate: 0.6829 BranchDensity: 0.2869 | Top regions: 0xdf1(493040) 0xdf0(108653) 0x1121(61908)
[Window 1] Branches: 276869 TakenRate: 0.6937 BranchDensity: 0.2769 | Top regions: 0xdf1(328262) 0xdf0(144707) 0x1121(105277)
[Window 8] Branches: 153590 TakenRate: 0.4569 BranchDensity: 0.1536 | Top regions: 0x10af(300391) 0x10a7(266996) 0x10eb(176219)
```

This reveals:
- **Phase behavior**: Different code regions dominate in different windows (phase changes)
- **Function boundaries**: Rapid region switches indicate calling patterns
- **Library activity**: Distinct memory regions may correspond to different libraries

**Interpretation example:**
- Windows 0-7 dominated by `0xdf1`, `0xdf0` → main function phase
- Window 8 switches to `0x10af`, `0x10a7` → different code section (library or other function)
- Windows 10+ show rapid changes → multiple function calls or loop phases

---

# Memory Access Patterns

The analyzer tracks **PC jump sizes** to reveal instruction locality and detect context switches.

Example output:

```
==== MEMORY ACCESS PATTERNS ====
PC jump size distribution (bucket = log2 of bytes):
  [ 2]            4-           7 bytes:   20866464 (86.8%)
  [ 3]            8-          15 bytes:     825956 (3.4%)
  [ 4]           16-          31 bytes:     217291 (0.9%)
  [ 7]          128-         255 bytes:     293256 (1.2%)
  [20]      1048576-     2097151 bytes:      23795 (0.1%)
  [21]      2097152-     4194303 bytes:     122633 (0.5%)
```

## Interpreting Jump Sizes

- **4-7 bytes (86%)**: Sequential instruction execution. Typical x86 instruction lengths.
- **8-63 bytes (3-5%)**: Moderate control flow (small jumps within functions)
- **64 bytes - 1MB**: Larger function calls or library jumps
- **1MB - 1GB**: Cross-function, cross-library jumps
- **> 1GB**: Rare; indicates potential context switches or kernel calls

## What This Reveals

- **High sequential execution** (>80% in 4-7B range) = tight code locality, good cache behavior
- **Few large jumps** = no heavy syscalls or context switches
- **Regular 2MB jumps** = pattern of large page allocations or multi-section code
- **Presence of 0xffff... regions** = kernel/stack activity detected

---

# PC Region Histogram

The analyzer groups all instructions into **memory regions** and counts occurrences, displayed with **configurable granularity**.

Example with 4KB pages (shift=12):

```
==== PC REGION HISTOGRAM (Top 16) ====
Region size: 4KB (shift=12)
Region 0xdf1 Count: 1515828
Region 0xdf0 Count: 1192845
Region 0x1121 Count: 974938
Region 0x10a7 Count: 901135
```

Example with 2MB pages (shift=21):

```
==== PC REGION HISTOGRAM (Top 16) ====
Region size: 2MB (shift=21)
Region 0x6 Count: 9603097
Region 0x8 Count: 5722933
Region 0x3 Count: 3424565
```

Example with 1GB sections (shift=30):

```
==== PC REGION HISTOGRAM (Top 16) ====
Region size: 1GB (shift=30)
Region 0x0 Count: 24010528
Region 0x3fffed7d6 Count: 18222
```

## Region Shift Values (ARM v8)

| Shift | Size | Region Count | Use Case |
|-------|------|--------------|----------|
| 12 | 4KB | 2^52 | Page-level analysis, fine-grained locality |
| 21 | 2MB | 2^43 | Large page (huge page) analysis |
| 30 | 1GB | 2^34 | Section-level, coarse memory layout |

**Default (12 bits)** follows ARM v8 4KB page architecture, enabling detailed memory locality analysis.
