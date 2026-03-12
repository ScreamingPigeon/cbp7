# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CBP-NG (Next-Generation Championship in Branch Prediction) simulator framework. The goal is to develop a high-performance branch predictor that maximizes a VFS (Voltage-Frequency-Scaled) score balancing prediction accuracy, energy efficiency, and area/complexity. Predictors are written using the HARCOM C++ library which models hardware costs through opaque types and operator overloading.

## Build Flow

The build is driven by `params.yaml` → `scripts/gen_predictor_config.py` → `build/predictor.mk` → `Makefile`. This YAML-first system (authored by Yuwei Sun) lets you tune predictor selection, template args, build type, sanitizers, and trace settings without touching the Makefile.

### Typical workflow

```bash
# 1. Edit params.yaml (predictor, template_args, build options, trace)
# 2. Generate build/predictor.mk from params.yaml
make predictor-config
# 3. Build the simulator (auto-includes build/predictor.mk)
make cbp
# 4. Run on a trace
./cbp ./gcc_test_trace.gz test 1000000 40000000
# Or via make (outputs to out/<trace_name>.out):
make run-cbp
```

### Other targets

```bash
make reference       # Build TAGE-SC-L reference predictor (accuracy baseline, not HARCOM-based)
make all             # Build both cbp and reference
make clean           # Remove cbp, reference, and build/predictor.mk
```

### Batch simulation and scoring

```bash
mkdir OUTDIR && ./run_all ./cbp ./traces OUTDIR
./predictor_metrics.py OUTDIR | ./vfs.py
```

### params.yaml options

`gen_predictor_config.py` reads `params.yaml` and writes `build/predictor.mk`, which overrides `PREDICTOR_TYPE`, `TRACE`, `TRACE_NAME`, `WARMUP`, `MEASURE`, `CBP_WARN_FLAGS`, `EXTRA_COMMON_FLAGS`, and `EXTRA_CBP_FLAGS` in the Makefile.

| Key | Effect |
|-----|--------|
| `predictor` / `template_args` | Assembled into `PREDICTOR_TYPE`, e.g. `Custom<>` or `tage<6,8,11,12,11,120,14,8>` |
| `predictor_type` | Overrides `predictor`+`template_args` with a raw string (escape hatch) |
| `build_type` | `release` (default, `-O3`), `debug` (`-O0 -g`), or `relwithdebinfo` (`-O2 -g`) |
| `warnings_as_errors` | Adds/removes `-Werror` from CBP warning flags (default `true`) |
| `native_arch` | Adds `-march=native` (default `false`) |
| `sanitize` | `none`, `address`, or `undefined` — injects `-fsanitize=...` flags |
| `defines` | List of macros, each becomes `-D<macro>` (e.g. `TAGE_VERBOSE`) |
| `extra_include_dirs` | List of paths, each becomes `-I<path>` |
| `trace` / `trace_name` / `warmup` / `measure` | Simulation parameters passed to run targets |
| `extra_common_flags` / `extra_cbp_flags` | Raw flags appended as-is |

CLI Makefile variable overrides (e.g. `make cbp CXX=clang++`) still take precedence over the generated `.mk` values.

### Bypassing params.yaml

Compile a specific predictor directly:
```bash
./compile cbp -DPREDICTOR="gshare<>"
```

## Build Requirements

- C++20: GCC 12+ or Clang 19+
- zlib (`-lz`) for gzip trace decompression
- Python 3 for utility scripts
- Optional: `nix-shell` via `shell.nix`

## Architecture

### Simulator Core

- **`cbp.cpp`** — Entry point. Parses args, runs simulation loop.
- **`cbp.hpp`** — Defines the `predictor` interface and `harcom_superuser` simulation engine.
- **`branch_predictor.hpp`** — Instantiates the predictor class selected by the `PREDICTOR` macro.
- **`trace_reader.hpp`** — Reads gzip-compressed branch trace files.
- **`harcom.hpp`** — HARCOM library (~6K lines). Provides `val<N>`, `reg<N>`, `arr<T,N>`, `ram<T,N>`, `rwram<N,M,B>`, `hard<V>` types that track energy, latency, and area.

### Two-Level Prediction Pipeline

Predictors implement two prediction levels for pipelining:
- **P1** (fast, lower accuracy): `predict1()` / `reuse_predict1()`
- **P2** (slower, higher accuracy): `predict2()` / `reuse_predict2()`
- **Training**: `update_condbr()` on misprediction, `update_cycle()` at block boundaries
- **Block control**: Predictors call `reuse_prediction(1)` to extend fetch blocks across multiple instructions

### Predictor Implementations (`predictors/`)

- **`common.hpp`** — Shared components: saturating counters, history folding (`folded_history`), banked RAM (`rwram`)
- **`bimodal.hpp`** — Simple PC-indexed prediction tables
- **`gshare.hpp`** — Global-history XOR predictor extending bimodal
- **`tage.hpp`** — Full TAGE predictor with geometric history tables, tag matching, meta-prediction, u-bit reset
- **`custom.hpp`** — Competition entry scaffold (the `Custom` class). Currently being developed on this branch.
- **`custom/TageTable.hpp`** — Template TAGE table structure (in progress)

### HARCOM Key Concepts

- **Opaque types**: `val<N>` values cannot be read directly; use `select()` (mux), `fold_xor()`, comparison operators instead of `if/else`
- **`hard<V>`**: Compile-time constants with zero latency cost
- **`reg<N>`**: Persistent state across cycles (writable once per cycle)
- **`ram`/`rwram`**: Hardware memory arrays tracked for area/energy; `rwram` supports banked read-modify-write
- **Automatic cost tracking**: All operations on HARCOM types accumulate energy and latency metrics transparently

### Scoring

Output CSV fields: trace_name, instructions, branches, conditional_branches, blocks, extra_cycles, P1/P2_divergences, divergences_at_block_end, P2_mispredictions, P1_latency, P2_latency, energy_per_instruction_fJ.

`predictor_metrics.py` aggregates IPC/CPI/EPI across traces; `vfs.py` computes the final VFS competition score.

## Current Development

The `prakhar` branch is developing a custom TAGE predictor in `predictors/custom.hpp` with parameterized `TageTable` components in `predictors/custom/TageTable.hpp`. The `Custom` class template takes extensive compile-time parameters for table sizes, history lengths, counter widths, etc.
