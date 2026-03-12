#!/usr/bin/env bash
# Compare two predictors side-by-side on one or more traces.
#
# Single trace:
#   scripts/compare_predictors.sh <pred_A> <pred_B> <trace_path> <warmup> <measure> [extra_flags]
#
# All traces in a directory:
#   scripts/compare_predictors.sh --dir <pred_A> <pred_B> <trace_dir> <warmup> <measure> [extra_flags]
#
# Examples:
#   scripts/compare_predictors.sh 'Tage<>' 'tage<>' ../traces/502-gcc-all_16112_trace.gz 1000000 40000000
#   scripts/compare_predictors.sh --dir 'Tage<>' 'tage<>' ../traces/ 1000000 40000000

set -euo pipefail

# --- Parse mode ---
DIR_MODE=false
if [[ "${1:-}" == "--dir" ]]; then
    DIR_MODE=true
    shift
fi

PRED_A="${1:?Usage: $0 [--dir] <pred_A> <pred_B> <trace_or_dir> <warmup> <measure> [extra_flags]}"
PRED_B="${2:?}"
TRACE_ARG="${3:?}"
WARMUP="${4:?}"
MEASURE="${5:?}"
EXTRA_FLAGS="${6:-}"

PYTHON="${PYTHON:-python3}"
CXX="${CXX:-g++}"
COMMON_FLAGS="-std=c++20 -O3 -DVERBOSE"
WARN_FLAGS="-Wall -Wextra -pedantic -Wold-style-cast -Wno-deprecated-declarations -Wno-mismatched-tags"
MISP_PENALTY=8
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# --- Build both predictors in parallel ---
echo "Building predictor A: $PRED_A" >&2
$CXX $COMMON_FLAGS $WARN_FLAGS $EXTRA_FLAGS -Itrace_files -o "$TMPDIR/cbp_a" cbp.cpp -lz -DPREDICTOR="$PRED_A" &
PID_A=$!

echo "Building predictor B: $PRED_B" >&2
$CXX $COMMON_FLAGS $WARN_FLAGS $EXTRA_FLAGS -Itrace_files -o "$TMPDIR/cbp_b" cbp.cpp -lz -DPREDICTOR="$PRED_B" &
PID_B=$!

wait $PID_A || { echo "Build failed for A: $PRED_A" >&2; exit 1; }
wait $PID_B || { echo "Build failed for B: $PRED_B" >&2; exit 1; }

# --- Derive trace name from filename ---
trace_name_from_path() {
    local base
    base="$(basename "$1")"
    # Strip _trace.gz or .gz suffix, then take the part before the first underscore-number
    base="${base%_trace.gz}"
    base="${base%.gz}"
    echo "$base"
}

# --- Parse human-readable output ---
parse_val() {
    grep -m1 "$1" "$2" | sed 's/.*: *//' | sed 's/ .*//'
}

# --- Run one trace, print comparison, output CSV line to fd 3 ---
run_one_trace() {
    local trace_path="$1"
    local trace_name
    trace_name="$(trace_name_from_path "$trace_path")"

    echo "Running $trace_name ..." >&2

    # Run predictor A
    "$TMPDIR/cbp_a" "$trace_path" "$trace_name" "$WARMUP" "$MEASURE" --format human \
        > "$TMPDIR/out_a.txt" 2> "$TMPDIR/hw_a.txt"

    # Run predictor B
    "$TMPDIR/cbp_b" "$trace_path" "$trace_name" "$WARMUP" "$MEASURE" --format human \
        > "$TMPDIR/out_b.txt" 2> "$TMPDIR/hw_b.txt"

    # Parse A
    local A_INSTR A_CONDBR A_PREDS A_EXTRA A_SHORT A_BLKSHORT A_MISPRED A_P1LAT A_P2LAT A_EPI
    local A_BITS A_XTORS A_AREA A_DYNPOW A_STAPOW
    A_INSTR=$(parse_val "instructions" "$TMPDIR/out_a.txt")
    A_CONDBR=$(parse_val "conditional" "$TMPDIR/out_a.txt")
    A_PREDS=$(parse_val "predictions" "$TMPDIR/out_a.txt")
    A_EXTRA=$(parse_val "extra_cycles" "$TMPDIR/out_a.txt")
    A_SHORT=$(parse_val "short mispredictions" "$TMPDIR/out_a.txt")
    A_BLKSHORT=$(parse_val "block-ending" "$TMPDIR/out_a.txt")
    A_MISPRED=$(parse_val "^mispredictions" "$TMPDIR/out_a.txt")
    A_P1LAT=$(parse_val "p1 latency" "$TMPDIR/out_a.txt")
    A_P2LAT=$(parse_val "p2 latency" "$TMPDIR/out_a.txt")
    A_EPI=$(parse_val "energy per" "$TMPDIR/out_a.txt")
    A_BITS=$(parse_val "storage" "$TMPDIR/hw_a.txt")
    A_XTORS=$(parse_val "transistors:" "$TMPDIR/hw_a.txt")
    A_AREA=$(parse_val "SRAM area" "$TMPDIR/hw_a.txt")
    A_DYNPOW=$(parse_val "dynamic power" "$TMPDIR/hw_a.txt")
    A_STAPOW=$(parse_val "static power" "$TMPDIR/hw_a.txt")

    # Parse B
    local B_INSTR B_CONDBR B_PREDS B_EXTRA B_SHORT B_BLKSHORT B_MISPRED B_P1LAT B_P2LAT B_EPI
    local B_BITS B_XTORS B_AREA B_DYNPOW B_STAPOW
    B_INSTR=$(parse_val "instructions" "$TMPDIR/out_b.txt")
    B_CONDBR=$(parse_val "conditional" "$TMPDIR/out_b.txt")
    B_PREDS=$(parse_val "predictions" "$TMPDIR/out_b.txt")
    B_EXTRA=$(parse_val "extra_cycles" "$TMPDIR/out_b.txt")
    B_SHORT=$(parse_val "short mispredictions" "$TMPDIR/out_b.txt")
    B_BLKSHORT=$(parse_val "block-ending" "$TMPDIR/out_b.txt")
    B_MISPRED=$(parse_val "^mispredictions" "$TMPDIR/out_b.txt")
    B_P1LAT=$(parse_val "p1 latency" "$TMPDIR/out_b.txt")
    B_P2LAT=$(parse_val "p2 latency" "$TMPDIR/out_b.txt")
    B_EPI=$(parse_val "energy per" "$TMPDIR/out_b.txt")
    B_BITS=$(parse_val "storage" "$TMPDIR/hw_b.txt")
    B_XTORS=$(parse_val "transistors:" "$TMPDIR/hw_b.txt")
    B_AREA=$(parse_val "SRAM area" "$TMPDIR/hw_b.txt")
    B_DYNPOW=$(parse_val "dynamic power" "$TMPDIR/hw_b.txt")
    B_STAPOW=$(parse_val "static power" "$TMPDIR/hw_b.txt")

    # Compute derived metrics via python
    local A_MPKI A_IPC A_CPI A_MPI A_VFS
    read A_MPKI A_IPC A_CPI A_MPI A_VFS < <(compute_vfs $A_INSTR $A_PREDS $A_EXTRA $A_SHORT $A_BLKSHORT $A_MISPRED $A_P1LAT $A_P2LAT $A_EPI)

    local B_MPKI B_IPC B_CPI B_MPI B_VFS
    read B_MPKI B_IPC B_CPI B_MPI B_VFS < <(compute_vfs $B_INSTR $B_PREDS $B_EXTRA $B_SHORT $B_BLKSHORT $B_MISPRED $B_P1LAT $B_P2LAT $B_EPI)

    local A_TOTAL=$((A_MISPRED + A_SHORT))
    local B_TOTAL=$((B_MISPRED + B_SHORT))

    # Print table
    local SEP="----------------------------+--------------------+--------------------"
    local FMT="%-28s| %18s | %18s\n"

    echo ""
    echo "=============================================================="
    echo "  Predictor Comparison: $trace_name"
    echo "=============================================================="
    printf "$FMT" "Metric" "$PRED_A" "$PRED_B"
    echo "$SEP"

    printf "$FMT" "Instructions" "$A_INSTR" "$B_INSTR"
    printf "$FMT" "Conditional branches" "$A_CONDBR" "$B_CONDBR"

    echo "$SEP"
    printf "$FMT" "Predictions (blocks)" "$A_PREDS" "$B_PREDS"
    printf "$FMT" "Extra cycles" "$A_EXTRA" "$B_EXTRA"
    printf "$FMT" "Short mispredictions" "$A_SHORT" "$B_SHORT"
    printf "$FMT" "Long mispredictions" "$A_MISPRED" "$B_MISPRED"
    printf "$FMT" "Total mispredictions" "$A_TOTAL" "$B_TOTAL"
    printf "$FMT" "MPKI" "$A_MPKI" "$B_MPKI"

    echo "$SEP"
    printf "$FMT" "P1 latency (cycles)" "$A_P1LAT" "$B_P1LAT"
    printf "$FMT" "P2 latency (cycles)" "$A_P2LAT" "$B_P2LAT"
    printf "$FMT" "IPC_cbp" "$A_IPC" "$B_IPC"
    printf "$FMT" "CPI_wrong" "$A_CPI" "$B_CPI"
    printf "$FMT" "MPI" "$A_MPI" "$B_MPI"
    printf "$FMT" "Energy/instr (fJ)" "$A_EPI" "$B_EPI"
    printf "$FMT" "VFS" "$A_VFS" "$B_VFS"

    echo "$SEP"
    printf "$FMT" "Storage (bits)" "$A_BITS" "$B_BITS"
    printf "$FMT" "Transistors" "$A_XTORS" "$B_XTORS"
    printf "$FMT" "SRAM area (mm2)" "$A_AREA" "$B_AREA"
    printf "$FMT" "Dynamic power (mW)" "$A_DYNPOW" "$B_DYNPOW"
    printf "$FMT" "Static power (mW)" "$A_STAPOW" "$B_STAPOW"
    echo "$SEP"

    # Emit CSV line for multi-trace aggregation (fd 3 if open, else skip)
    if { true >&3; } 2>/dev/null; then
        echo "$trace_name,$A_MPKI,$A_IPC,$A_CPI,$A_EPI,$A_VFS,$A_P1LAT,$A_P2LAT,$A_MISPRED,$A_SHORT,$B_MPKI,$B_IPC,$B_CPI,$B_EPI,$B_VFS,$B_P1LAT,$B_P2LAT,$B_MISPRED,$B_SHORT" >&3
    fi
}

# --- VFS computation (matches predictor_metrics.py + vfs.py) ---
compute_vfs() {
    $PYTHON -c "
import math
misp_penalty = $MISP_PENALTY
instr, preds, extra = $1, $2, $3
diverge, diverge_end, misp = $4, $5, $6
p1_lat, p2_lat, epi = $7, $8, $9

p1_latency = math.ceil(p1_lat)
p2_latency = math.ceil(p2_lat)

mpki = misp / instr * 1000
mpi = misp / instr

if p2_latency <= p1_latency:
    cycles = preds * max(1, p2_latency)
else:
    cycles = preds * max(1, p1_latency) + diverge * p2_latency - diverge_end * max(1, p1_latency)
cycles += extra
ipc = instr / cycles

cpi = mpi * (misp_penalty + p2_latency)

IPCcbp0 = 8; CPIcbp0 = 0.0315; EPIcbp0 = 1000
ALPHA = 1.625; BETA = 4*ALPHA / (ALPHA-1)**2; GAMMA = 2 / (ALPHA-1)
cbp_energy_ratio = 0.05
WPI0 = IPCcbp0 * CPIcbp0; WPI = ipc * cpi
speedup = (ipc/IPCcbp0) * (1+WPI0)/(1+WPI)
LAMBDA = 1/(1+WPI0/2) - cbp_energy_ratio
normalizedEPI = ((epi/EPIcbp0) * cbp_energy_ratio + LAMBDA * speedup**GAMMA) * (1+WPI/2)
vfs = speedup * ALPHA * (1-2/(1+math.sqrt(1+BETA/(speedup*normalizedEPI))))

print(f'{mpki:.3f} {ipc:.4f} {cpi:.6f} {mpi:.6f} {vfs:.6f}')
"
}

# --- Main ---
if [[ "$DIR_MODE" == false ]]; then
    # Single trace mode
    if [[ ! -f "$TRACE_ARG" ]]; then
        echo "Error: trace file not found: $TRACE_ARG" >&2
        exit 1
    fi
    run_one_trace "$TRACE_ARG"
else
    # Directory mode — run all *_trace.gz files
    if [[ ! -d "$TRACE_ARG" ]]; then
        echo "Error: trace directory not found: $TRACE_ARG" >&2
        exit 1
    fi

    TRACES=()
    for f in "$TRACE_ARG"/*_trace.gz; do
        [[ -f "$f" ]] && TRACES+=("$f")
    done

    if [[ ${#TRACES[@]} -eq 0 ]]; then
        echo "Error: no *_trace.gz files found in $TRACE_ARG" >&2
        exit 1
    fi

    echo "Found ${#TRACES[@]} traces in $TRACE_ARG" >&2

    CSV_FILE="$TMPDIR/results.csv"
    exec 3>"$CSV_FILE"

    for trace in "${TRACES[@]}"; do
        run_one_trace "$trace"
    done

    exec 3>&-

    # Print summary table
    echo ""
    echo ""
    echo "=========================================================================="
    echo "  Summary: $PRED_A  vs  $PRED_B  (${#TRACES[@]} traces)"
    echo "=========================================================================="

    $PYTHON -c "
import sys, math

lines = open('$CSV_FILE').read().strip().split('\n')
if not lines or lines == ['']:
    print('No results to summarize.')
    sys.exit(0)

# CSV columns:
# trace, a_mpki, a_ipc, a_cpi, a_epi, a_vfs, a_p1lat, a_p2lat, a_misp, a_short,
#        b_mpki, b_ipc, b_cpi, b_epi, b_vfs, b_p1lat, b_p2lat, b_misp, b_short

hdr_fmt = '%-30s| %8s %8s | %8s %8s | %8s %8s | %8s %8s'
row_fmt = '%-30s| %8s %8s | %8s %8s | %8s %8s | %8s %8s'
sep     = '-' * 30 + '+' + '-' * 19 + '+' + '-' * 19 + '+' + '-' * 19 + '+' + '-' * 19

print(hdr_fmt % ('Trace',
    'MPKI_A', 'MPKI_B',
    'IPC_A', 'IPC_B',
    'VFS_A', 'VFS_B',
    'EPI_A', 'EPI_B'))
print(sep)

# Accumulators for averages
n = 0
sum_a_mpki = sum_b_mpki = 0
sum_a_ipc_inv = sum_b_ipc_inv = 0  # harmonic mean
sum_a_cpi = sum_b_cpi = 0
sum_a_epi = sum_b_epi = 0
sum_a_vfs = sum_b_vfs = 0
max_a_p1 = max_b_p1 = 0
max_a_p2 = max_b_p2 = 0

for line in lines:
    cols = line.split(',')
    trace = cols[0]
    a_mpki, a_ipc, a_cpi, a_epi, a_vfs = float(cols[1]), float(cols[2]), float(cols[3]), float(cols[4]), float(cols[5])
    a_p1, a_p2 = float(cols[6]), float(cols[7])
    b_mpki, b_ipc, b_cpi, b_epi, b_vfs = float(cols[10]), float(cols[11]), float(cols[12]), float(cols[13]), float(cols[14])
    b_p1, b_p2 = float(cols[15]), float(cols[16])

    print(row_fmt % (trace[:30],
        f'{a_mpki:.3f}', f'{b_mpki:.3f}',
        f'{a_ipc:.2f}', f'{b_ipc:.2f}',
        f'{a_vfs:.4f}', f'{b_vfs:.4f}',
        f'{a_epi:.0f}', f'{b_epi:.0f}'))

    n += 1
    sum_a_mpki += a_mpki; sum_b_mpki += b_mpki
    sum_a_ipc_inv += 1/a_ipc; sum_b_ipc_inv += 1/b_ipc
    sum_a_cpi += a_cpi; sum_b_cpi += b_cpi
    sum_a_epi += a_epi; sum_b_epi += b_epi
    sum_a_vfs += a_vfs; sum_b_vfs += b_vfs
    max_a_p1 = max(max_a_p1, a_p1); max_b_p1 = max(max_b_p1, b_p1)
    max_a_p2 = max(max_a_p2, a_p2); max_b_p2 = max(max_b_p2, b_p2)

print(sep)

avg_a_mpki = sum_a_mpki / n; avg_b_mpki = sum_b_mpki / n
avg_a_ipc = n / sum_a_ipc_inv; avg_b_ipc = n / sum_b_ipc_inv  # harmonic mean
avg_a_cpi = sum_a_cpi / n; avg_b_cpi = sum_b_cpi / n
avg_a_epi = sum_a_epi / n; avg_b_epi = sum_b_epi / n
avg_a_vfs = sum_a_vfs / n; avg_b_vfs = sum_b_vfs / n

print(row_fmt % ('AVERAGE',
    f'{avg_a_mpki:.3f}', f'{avg_b_mpki:.3f}',
    f'{avg_a_ipc:.2f}', f'{avg_b_ipc:.2f}',
    f'{avg_a_vfs:.4f}', f'{avg_b_vfs:.4f}',
    f'{avg_a_epi:.0f}', f'{avg_b_epi:.0f}'))

print()
print('Averages (detailed):')
avg_fmt = '  %-22s  A: %-14s  B: %-14s'
print(avg_fmt % ('MPKI (arith)', f'{avg_a_mpki:.3f}', f'{avg_b_mpki:.3f}'))
print(avg_fmt % ('IPC_cbp (harmonic)', f'{avg_a_ipc:.4f}', f'{avg_b_ipc:.4f}'))
print(avg_fmt % ('CPI_wrong (arith)', f'{avg_a_cpi:.6f}', f'{avg_b_cpi:.6f}'))
print(avg_fmt % ('EPI (arith)', f'{avg_a_epi:.1f}', f'{avg_b_epi:.1f}'))
print(avg_fmt % ('VFS (arith)', f'{avg_a_vfs:.6f}', f'{avg_b_vfs:.6f}'))
print(avg_fmt % ('P1 latency (max)', f'{max_a_p1:.3f}', f'{max_b_p1:.3f}'))
print(avg_fmt % ('P2 latency (max)', f'{max_a_p2:.3f}', f'{max_b_p2:.3f}'))
print()

# VFS uses harmonic-mean IPC, arith-mean CPI, arith-mean EPI, max latencies
# (matching predictor_metrics.py aggregation)
def compute_vfs(ipc, cpi, epi):
    IPCcbp0 = 8; CPIcbp0 = 0.0315; EPIcbp0 = 1000
    ALPHA = 1.625; BETA = 4*ALPHA / (ALPHA-1)**2; GAMMA = 2 / (ALPHA-1)
    cbp_energy_ratio = 0.05
    WPI0 = IPCcbp0 * CPIcbp0; WPI = ipc * cpi
    speedup = (ipc/IPCcbp0) * (1+WPI0)/(1+WPI)
    LAMBDA = 1/(1+WPI0/2) - cbp_energy_ratio
    normalizedEPI = ((epi/EPIcbp0) * cbp_energy_ratio + LAMBDA * speedup**GAMMA) * (1+WPI/2)
    return speedup * ALPHA * (1-2/(1+math.sqrt(1+BETA/(speedup*normalizedEPI))))

agg_a_vfs = compute_vfs(avg_a_ipc, avg_a_cpi, avg_a_epi)
agg_b_vfs = compute_vfs(avg_b_ipc, avg_b_cpi, avg_b_epi)
print(avg_fmt % ('Aggregate VFS', f'{agg_a_vfs:.6f}', f'{agg_b_vfs:.6f}'))
print()
"
fi
