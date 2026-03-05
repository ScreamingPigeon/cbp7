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

# PC Region Histogram

The analyzer groups branch PCs into **memory regions** and counts how many branches occur in each region.

Example:

```
==== PC REGION HISTOGRAM (Top 8) ====
Region 0x0 Count: 3001893
Region 0xffff Count: 3107
```

This provides a coarse view of where instructions originate in memory.
