#!/bin/bash
#
# eval_monitor.sh — Run TAMonitor on the representative trace subset.
#
# Usage:
#   ./scripts/eval_monitor.sh <binary> [trace_dir] [out_dir] [jobs]
#
# Produces per-trace files:
#   <out_dir>/<trace>_csv.txt      — windowed CSV (stdout from monitor)
#   <out_dir>/<trace>_summary.txt  — end-of-trace summary (stderr from monitor)
#
# Then runs aggregate_monitor.py to produce:
#   <out_dir>/aggregate_csv.txt    — per-window mean across traces
#   <out_dir>/aggregate_summary.txt — combined summary stats

set -euo pipefail

EXEC="${1:?Usage: $0 <binary> [trace_dir] [out_dir] [jobs]}"
TRACE_DIR="${2:-./traces}"
OUT_DIR="${3:-out/eval_monitor}"
MAXJOBS="${4:-$(nproc)}"

WARMINST=1000000
SIMINST=40000000

# Same representative subset as quick_eval.sh
REPR_TRACES=(
    "502-gcc-all_16112_trace.gz"
    "505-mcf-1_14364_trace.gz"
    "531-deepsjeng-1_53703_trace.gz"
    "548-exchange2-1_207156_trace.gz"
    "508-namd-1_126793_trace.gz"
    "554-roms-1_62613_trace.gz"
    "llvm-2.25945_0_trace.gz"
    "gcc-1.52513_0_trace.gz"
    "java16-specjbb-64k-ir1000.19011_0_trace.gz"
    "dcapo-kafka-jdk8.524_0_trace.gz"
    "web_74_trace.gz"
    "web_130_trace.gz"
    "nodejs-octane_3483_trace.gz"
    "nodejs-http2_4818_trace.gz"
    "rsbench-1.730_0_trace.gz"
    "sampleflow-1.127768_0_trace.gz"
    "zstd-1.19139_0_trace.gz"
    "gap-sssp-usroad.287_0_trace.gz"
    "python3-pyperf-dulwich.2000_0_trace.gz"
    "lua-3.25585_0_trace.gz"
)

# Verify traces
MISSING=0
for t in "${REPR_TRACES[@]}"; do
    if [[ ! -f "${TRACE_DIR}/${t}" ]]; then
        echo "WARNING: missing trace: ${TRACE_DIR}/${t}" >&2
        MISSING=$((MISSING + 1))
    fi
done
if [[ $MISSING -gt 0 ]]; then
    echo "WARNING: ${MISSING} traces missing, running with available subset" >&2
fi

mkdir -p "${OUT_DIR}"

echo "Running monitor on ${#REPR_TRACES[@]} representative traces (${MAXJOBS} parallel jobs)..."

printf '%s\n' "${REPR_TRACES[@]}" | \
    xargs -P "${MAXJOBS}" -I{} bash -c '
        trace="$1"
        trace_dir="$2"
        out_dir="$3"
        exec_bin="$4"
        warmup="$5"
        measure="$6"
        tracefile="${trace_dir}/${trace}"
        if [[ ! -f "$tracefile" ]]; then exit 0; fi
        filename=$(basename "$trace")
        filename="${filename%.*}"
        name="${filename%_trace}"
        # Monitor: stdout=scoring line, stderr=windowed CSV then summary
        # Capture both, then split: CSV lines (# header + data) vs summary (starts with ===)
        "$exec_bin" "$tracefile" "$name" "$warmup" "$measure" \
            > "${out_dir}/${name}_score.txt" \
            2> "${out_dir}/${name}_raw.txt"
        # Split: lines before first "===" are windowed CSV, rest is summary
        awk "/^===/{found=1} !found{print}" "${out_dir}/${name}_raw.txt" > "${out_dir}/${name}_csv.txt"
        awk "/^===/{found=1} found{print}" "${out_dir}/${name}_raw.txt" > "${out_dir}/${name}_summary.txt"
        rm -f "${out_dir}/${name}_raw.txt"
        echo "  done: ${name}"
    ' _ {} "${TRACE_DIR}" "${OUT_DIR}" "${EXEC}" "${WARMINST}" "${SIMINST}"

echo ""

# Count successful runs
NRUNS=$(find "${OUT_DIR}" -name '*_csv.txt' -size +0c | wc -l)
if [[ $NRUNS -eq 0 ]]; then
    echo "ERROR: no successful runs" >&2
    exit 1
fi
echo "=== Monitor Eval: ${NRUNS}/${#REPR_TRACES[@]} traces completed ==="

# Per-trace MPKI summary
echo ""
printf "%-45s %10s\n" "Trace" "MPKI"
printf "%-45s %10s\n" "-----" "----"
for f in "${OUT_DIR}"/*_summary.txt; do
    [[ -s "$f" ]] || continue
    name=$(basename "$f" _summary.txt)
    mpki=$(grep -oP 'MPKI:\s+\K[0-9.]+' "$f" 2>/dev/null | head -1 || echo "?")
    printf "%-45s %10s\n" "$name" "$mpki"
done | sort

echo ""

# Aggregate across traces
echo "=== Generating aggregate CSV ==="
python3 scripts/aggregate_monitor.py "${OUT_DIR}" "${OUT_DIR}/aggregate_csv.txt" "${OUT_DIR}/aggregate_summary.txt"
echo ""
echo "Per-trace files:  ${OUT_DIR}/<trace>_{csv,summary}.txt"
echo "Aggregate CSV:    ${OUT_DIR}/aggregate_csv.txt"
echo "Aggregate summary: ${OUT_DIR}/aggregate_summary.txt"
