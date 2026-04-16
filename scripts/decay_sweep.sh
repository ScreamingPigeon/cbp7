#!/bin/bash
set -e

# Decay sweep: compare different DECAY_CTR / policy configs
# All use minhist=4 (best VFS from iteration 1)

CXX="${CXX:-g++}"
BUILD=build
TRACE_DIR=traces
JOBS=4
OUTDIR=out/decay_sweep

COMMON="-std=c++20 -O3"
WARNS="-Wall -Wextra -pedantic -Wold-style-cast -Werror -Wno-deprecated-declarations -Wno-mismatched-tags"

# SweepTableConfig with minhist=4
TABLE="SweepTableConfig<8,2048,11,1,2,1,4,100,1>"
# Common prefix: TableCfg, AllocCfg, FW, BIM, BR, BANKS, STAG, SU, SHYS, UFF
PREFIX="$TABLE,DefaultAllocConfig,16,4096,1,1,true,true,true,false"

# Configs: vary DECAY_CTR, DECAY_GRAN, DECAY_POLICY
# After UFF: DECAY_CTR, DECAY_GRAN, POLICY — rest defaults
declare -a NAMES
declare -A CONFIGS

NAMES=(decay0_baseline decay8_mild_g2 decay8_aggressive_g0 decay4_mild_g0 decay8_hybrid_g2)
CONFIGS[decay0_baseline]="Tage<$PREFIX,0,0,DecayMild>"
CONFIGS[decay8_mild_g2]="Tage<$PREFIX,8,2,DecayMild>"
CONFIGS[decay8_aggressive_g0]="Tage<$PREFIX,8,0,DecayAggressive>"
CONFIGS[decay4_mild_g0]="Tage<$PREFIX,4,0,DecayMild>"
CONFIGS[decay8_hybrid_g2]="Tage<$PREFIX,8,2,DecayHybrid>"

mkdir -p $BUILD $OUTDIR

echo "=== Building configs ==="
for name in "${NAMES[@]}"; do
    pred="${CONFIGS[$name]}"
    echo "  Building $name..."
    $CXX $COMMON $WARNS -Itrace_files -o $BUILD/decay_$name cbp.cpp -lz -DPREDICTOR="$pred" &
done
wait
echo "  All builds done."

echo ""
echo "=== Running quick-eval (20 traces) ==="
for name in "${NAMES[@]}"; do
    echo ""
    echo "--- $name ---"
    echo "    ${CONFIGS[$name]}"
    scripts/quick_eval.sh ./$BUILD/decay_$name $TRACE_DIR $OUTDIR/$name $JOBS
done

echo ""
echo "=== Summary ==="
echo ""
printf "%-30s %10s %10s %10s %10s\n" "Config" "MPKI" "IPC" "VFS" "P2_lat"
printf "%-30s %10s %10s %10s %10s\n" "------------------------------" "----------" "----------" "----------" "----------"
for name in "${NAMES[@]}"; do
    csv="$OUTDIR/$name/results.csv"
    if [ -f "$csv" ]; then
        # results.csv: trace,mpki,epi,ipc,vfs,p1_latency,p2_latency
        # Get the AVERAGE line (last line)
        last=$(tail -1 "$csv")
        mpki=$(echo "$last" | awk -F, '{print $2}')
        ipc=$(echo "$last" | awk -F, '{print $4}')
        vfs=$(echo "$last" | awk -F, '{print $5}')
        p2=$(echo "$last" | awk -F, '{print $7}')
        printf "%-30s %10s %10s %10s %10s\n" "$name" "$mpki" "$ipc" "$vfs" "$p2"
    else
        printf "%-30s %10s\n" "$name" "NO RESULTS"
    fi
done
