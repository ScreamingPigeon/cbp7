# CBP Trace Analyzer Overview

This document describes how `trace_analyzer.cpp` uses the provided `trace_reader.hpp` to analyze CBP trace files and what metrics are produced by the program.

## How the Analyzer Uses `trace_reader.hpp`

The CBP traces are binary files (usually gzip-compressed) that contain instruction-level execution information. The file format itself is complex, so the provided `trace_reader.hpp` library handles all parsing and decoding.

The analyzer interacts with the trace through the `trace_reader` class.

Typical usage flow:

1. The program constructs a `trace_reader` object with the trace file path and a name for the trace.

   Example:

   ```cpp
   trace_reader reader(trace_path, "trace");
   ```

2. The analyzer repeatedly calls:

   ```cpp
   instruction inst = reader.next_instruction();
   ```

   Each call returns an `instruction` structure representing the next executed instruction in the trace.

3. The returned `instruction` struct contains the following fields used by the analyzer:

   - `pc`  
     The program counter of the instruction.

   - `next_pc`  
     The next program counter after the instruction executes. For branches that are taken, this will be the branch target.

   - `inst_class`  
     The instruction classification (ALU, LOAD, STORE, branch types, etc.).

   - `branch`  
     Boolean indicating whether the instruction is a branch.

   - `taken_branch`  
     Boolean indicating whether the branch was taken.

4. The analyzer updates counters and statistics for each instruction.

5. When the trace reaches the end of file, `trace_reader` throws an `out_of_instructions` exception, which terminates the analysis loop.

The analyzer itself does not interpret the binary trace format directly; it simply consumes decoded instructions provided by `trace_reader`.

---

## Metrics Produced by the Analyzer

### Total Instruction Count

The analyzer counts the total number of instructions processed from the trace.

This metric represents the full execution length of the workload contained in the trace.

---

### Total Branch Count

For every instruction where `inst.branch == true`, the analyzer increments the branch counter.

This measures how frequently control-flow instructions appear relative to total instructions.

---

### Taken and Not-Taken Branch Counts

For each branch instruction, the analyzer records whether it was taken or not.

Counters maintained:

- taken branches
- not-taken branches

These values describe the overall branch direction distribution.

---

### Global Branch Taken Rate

The global taken rate is computed as:

```
taken_branches / total_branches
```

This metric indicates whether branches in the program are generally biased toward taken or not-taken.

---

### Per-PC Branch Statistics

The analyzer keeps a table indexed by branch program counter (PC).  
For each branch PC, the analyzer tracks:

- number of times the branch executed
- number of taken outcomes
- number of not-taken outcomes

Example conceptual structure:

```
branch_stats[pc]:
    executions
    taken
    not_taken
```

This information helps identify how individual branch instructions behave.

---

### Hot Branch Identification

After processing the entire trace, the analyzer sorts branch PCs by execution count.

This reveals the most frequently executed branch instructions (often referred to as "hot branches"). These branches tend to dominate branch prediction performance because they appear most often in execution.

---

### Instruction Class Distribution

Using the `inst_class` field from the instruction structure, the analyzer counts how many instructions fall into each category, such as:

- ALU operations
- LOAD instructions
- STORE instructions
- BRANCH instructions
- CALL / RETURN instructions
- floating-point operations

This produces a high-level view of the workload's instruction mix.

---

## Typical Output

When the analysis finishes, the program prints summary statistics similar to the following:

```
Total Instructions: <value>
Total Branches: <value>
Taken Branches: <value>
Not Taken Branches: <value>
Global Taken Rate: <value>
```

It may also output:

- instruction class distribution
- most frequently executed branch PCs
- per-PC branch execution counts

---

## Purpose of the Analysis

These metrics provide a high-level characterization of the control-flow behavior of the program contained in the trace. This information is useful for:

- understanding branch behavior in the workload
- identifying frequently executed branches
- guiding branch predictor design and evaluation
- understanding instruction mix and program structure