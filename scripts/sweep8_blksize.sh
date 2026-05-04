#!/bin/bash
# Sweep 8: rwram bank shift — conflict stats + accuracy
# Quick run: 100k instructions per trace (conflict stats only need short runs)

set -euo pipefail

CONFIGS=(
  "S8_Shift0"   # baseline: bit 0
  "S8_Shift1"   # bit 1
  "S8_Shift2"   # bit 2
  "S8_Shift3"   # bit 3
  "S8_Shift5"   # bit 5
  "S8_Graded"   # T0-4: bit 2, T5-13: bit 3
)

WARMUP=1000
MEASURE=100000
TRACE_DIR="${1:-traces}"
JOBS="${2:-4}"

TRACES=(
  "502-gcc-all_16112_trace.gz"
  "531-deepsjeng-1_53703_trace.gz"
  "web_74_trace.gz"
  "nodejs-octane_3483_trace.gz"
  "llvm-2.25945_0_trace.gz"
)

echo "=== Sweep 8: Bank Shift (${#CONFIGS[@]} configs, ${#TRACES[@]} traces, ${MEASURE} instr) ==="
echo ""

for cfg in "${CONFIGS[@]}"; do
  echo "--- Building: ${cfg} ---"
  # Build monitor binary
  g++ -std=c++20 -O3 -Wall -Wextra -pedantic -Wno-deprecated-declarations -Wno-mismatched-tags \
    -DTAGE_MONITOR -DCHEATING_MODE -DFREE_FANOUT \
    -Itrace_files -o "build/sweep8_${cfg}" cbp.cpp -lz -DPREDICTOR="${cfg}" 2>&1 | grep -E "error:" | head -3
  if [ ! -f "build/sweep8_${cfg}" ]; then
    echo "  BUILD FAILED"
    continue
  fi

  outdir="out/sweep8_${cfg}"
  mkdir -p "${outdir}"

  for trace in "${TRACES[@]}"; do
    tracefile="${TRACE_DIR}/${trace}"
    [ -f "$tracefile" ] || continue
    name="${trace%_trace.gz}"
    "build/sweep8_${cfg}" "$tracefile" "$name" "$WARMUP" "$MEASURE" \
      > "${outdir}/${name}_score.txt" \
      2> "${outdir}/${name}_raw.txt"
    # Split raw into csv + summary
    awk '/^===/{found=1} !found{print}' "${outdir}/${name}_raw.txt" > "${outdir}/${name}_csv.txt"
    awk '/^===/{found=1} found{print}' "${outdir}/${name}_raw.txt" > "${outdir}/${name}_summary.txt"
    rm -f "${outdir}/${name}_raw.txt"
  done &
done
wait
echo ""

echo "=== Conflict Stats (per-table avg lost%) ==="
printf "%-15s" "Config"
for trace in "${TRACES[@]}"; do
  name="${trace%_trace.gz}"
  short=$(echo "$name" | cut -d- -f1-2 | head -c12)
  printf " %12s" "$short"
done
printf " %12s\n" "AVG_LOST%"

for cfg in "${CONFIGS[@]}"; do
  outdir="out/sweep8_${cfg}"
  printf "%-15s" "$cfg"
  total_lost=0
  count=0
  for trace in "${TRACES[@]}"; do
    name="${trace%_trace.gz}"
    sumfile="${outdir}/${name}_summary.txt"
    if [ -f "$sumfile" ]; then
      # Average lost% across all pred_ram lines (all 14 tables)
      avg=$(grep "rwram\[ta_pred\]" "$sumfile" | \
        grep -oP 'lost=\d+ \(\K[0-9.]+' | \
        awk '{s+=$1; n++} END {if(n>0) printf "%.1f", s/n; else print "N/A"}')
      printf " %11s%%" "$avg"
      total_lost=$(python3 -c "print(${total_lost} + ${avg})")
      count=$((count + 1))
    else
      printf " %12s" "N/A"
    fi
  done
  if [ $count -gt 0 ]; then
    overall=$(python3 -c "print(f'{${total_lost}/${count}:.1f}')")
    printf " %11s%%\n" "$overall"
  else
    printf " %12s\n" "N/A"
  fi
done
