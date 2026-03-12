#!/bin/bash
#
# quick_eval.sh — Evaluate a predictor on a representative trace subset.
#
# Usage:
#   ./scripts/quick_eval.sh <binary> [trace_dir] [out_dir] [jobs]
#
# Example:
#   make cbp && ./scripts/quick_eval.sh ./cbp ./traces out/quick 8
#
# The subset covers all major workload categories from the 168 training traces:
#   SPEC CPU (int, fp, server), compilers, Java/JVM, web/browser, Node.js,
#   scientific/HPC, compression, graph analytics, scripting (Python/Lua/Julia),
#   media, databases, misc tools.
#
# After simulation, prints per-trace metrics and an aggregate VFS score.

set -euo pipefail

EXEC="${1:?Usage: $0 <binary> [trace_dir] [out_dir] [jobs]}"
TRACE_DIR="${2:-./traces}"
OUT_DIR="${3:-out/quick}"
MAXJOBS="${4:-$(nproc)}"

WARMINST=1000000
SIMINST=40000000

# --- Representative subset (20 traces, one per major category) ---
# Chosen to span: SPEC int, SPEC fp, SPEC server, compilers, Java/JVM,
# web/browser, Node.js, scientific, compression, graph, scripting,
# media, database, misc tools. One trace per family avoids redundancy
# while covering diverse branch behavior patterns.
REPR_TRACES=(
    # SPEC CPU 2017 (int + fp + server)
    "502-gcc-all_16112_trace.gz"          # compiler (int, high branch density)
    "505-mcf-1_14364_trace.gz"            # combinatorial optimization (int, pointer-chasing)
    "531-deepsjeng-1_53703_trace.gz"      # game tree search (int, deep recursion)
    "548-exchange2-1_207156_trace.gz"     # puzzle solving (int, heavy branching)
    "508-namd-1_126793_trace.gz"          # molecular dynamics (fp, regular loops)
    "554-roms-1_62613_trace.gz"           # ocean modeling (fp, scientific)

    # Compilers & tools
    "llvm-2.25945_0_trace.gz"            # LLVM compiler (large codebase, polymorphic)
    "gcc-1.52513_0_trace.gz"             # GCC compiler

    # Java / JVM server
    "java16-specjbb-64k-ir1000.19011_0_trace.gz"  # SPECjbb (JIT, GC, server)
    "dcapo-kafka-jdk8.524_0_trace.gz"    # Kafka (streaming, JVM)

    # Web / browser
    "web_74_trace.gz"                    # web browsing workload
    "web_130_trace.gz"                   # web browsing workload (different site)

    # Node.js / JavaScript
    "nodejs-octane_3483_trace.gz"        # V8 benchmark suite
    "nodejs-http2_4818_trace.gz"         # HTTP/2 server workload

    # Scientific / HPC
    "rsbench-1.730_0_trace.gz"           # Monte Carlo neutron transport
    "sampleflow-1.127768_0_trace.gz"     # statistical sampling

    # Compression
    "zstd-1.19139_0_trace.gz"            # Zstd compression

    # Graph analytics
    "gap-sssp-usroad.287_0_trace.gz"     # single-source shortest path

    # Scripting languages
    "python3-pyperf-dulwich.2000_0_trace.gz"  # Python (interpreter overhead)
    "lua-3.25585_0_trace.gz"             # Lua interpreter
)

# Verify all traces exist
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

# Run traces in parallel
echo "Running ${#REPR_TRACES[@]} representative traces (${MAXJOBS} parallel jobs)..."
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
        "$exec_bin" "$tracefile" "$name" "$warmup" "$measure" > "${out_dir}/${name}.out" 2>/dev/null
        echo "  done: ${name}"
    ' _ {} "${TRACE_DIR}" "${OUT_DIR}" "${EXEC}" "${WARMINST}" "${SIMINST}"

echo ""

# Count successful runs
NRUNS=$(find "${OUT_DIR}" -name '*.out' -size +0c | wc -l)
if [[ $NRUNS -eq 0 ]]; then
    echo "ERROR: no successful runs" >&2
    exit 1
fi
echo "=== Quick Eval: ${NRUNS}/${#REPR_TRACES[@]} traces completed ==="
echo ""

# Per-trace summary
printf "%-45s %8s %8s %8s %8s\n" "Trace" "MPKI" "EPI(fJ)" "P2lat" "IPC_cbp"
printf "%-45s %8s %8s %8s %8s\n" "-----" "----" "-------" "-----" "-------"

for f in "${OUT_DIR}"/*.out; do
    [[ -s "$f" ]] || continue
    line=$(head -1 "$f")
    IFS=',' read -r name instr branch condbr npred extra diverge divend misp p1lat p2lat epi <<< "$line"
    mpki=$(python3 -c "print(f'{1000*${misp}/${instr}:.3f}')")
    ipc=$(python3 -c "
import math
p2_latency = math.ceil(${p2lat})
p1_latency = math.ceil(${p1lat})
if p2_latency <= p1_latency:
    cycles = ${npred} * max(1, p2_latency)
else:
    cycles = ${npred} * max(1, p1_latency) + ${diverge} * p2_latency - ${divend} * max(1, p1_latency)
cycles += ${extra}
print(f'{${instr}/cycles:.3f}')
")
    printf "%-45s %8s %8s %8s %8s\n" "$name" "$mpki" "$epi" "$p2lat" "$ipc"
done | sort

echo ""

# Aggregate VFS score
METRICS=$(python3 predictor_metrics.py "${OUT_DIR}")
VFS=$(echo "$METRICS" | python3 vfs.py)

IFS=',' read -r avg_ipc avg_cpi avg_epi avg_mpi rest <<< "$METRICS"
echo "=== Aggregate (${NRUNS} traces) ==="
printf "  IPC_cbp:  %s\n" "$avg_ipc"
printf "  CPI_cbp:  %s\n" "$avg_cpi"
printf "  EPI(fJ):  %s\n" "$avg_epi"
printf "  MPI:      %s\n" "$avg_mpi"
printf "  VFS:      %s\n" "$VFS"
echo ""
echo "Full results in: ${OUT_DIR}/"
echo "To run ALL 168 traces: mkdir -p out/full && ./run_all ${EXEC} ${TRACE_DIR} out/full"
