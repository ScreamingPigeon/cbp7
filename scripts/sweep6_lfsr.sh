#!/bin/bash
# Sweep 6: Graded LFSR width / decay threshold variants
# Runs eval-monitor for each config, then compares VFS/MPKI

set -euo pipefail

CONFIGS=(
  "S3_MW4_MC2048"   # baseline: uniform LFSR=8, FixedThresh<8>
  "S6_G10_7"        # graded LFSR 10→7, FixedThresh<8>
  "S6_G12_6"        # graded LFSR 12→6, FixedThresh<8>
  "S6_G10_7_T16"    # graded LFSR 10→7, FixedThresh<16>
  "S6_G12_6_T16"    # graded LFSR 12→6, FixedThresh<16>
  "S6_GT8_64"       # uniform LFSR=8, GradedThresh 8→64
  "S6_PG512"        # uniform LFSR=8, PressGated 4→32@512
)

JOBS="${1:-4}"

echo "=== Sweep 6: LFSR/Decay Threshold (${#CONFIGS[@]} configs, ${JOBS} jobs) ==="
echo ""

for cfg in "${CONFIGS[@]}"; do
  outdir="out/sweep6_${cfg}"
  echo "--- Building and running: ${cfg} ---"
  make eval-monitor EVAL_MONITOR_PRED="${cfg}" EVAL_MONITOR_OUT="${outdir}" EVAL_MONITOR_JOBS="${JOBS}" -B 2>&1 | grep -E "done:|ERROR|Running"
  echo ""
done

echo ""
echo "=== Results ==="
printf "%-20s %10s %10s %10s\n" "Config" "IPC" "MPKI" "VFS"
printf "%-20s %10s %10s %10s\n" "------" "---" "----" "---"

for cfg in "${CONFIGS[@]}"; do
  outdir="out/sweep6_${cfg}"
  tmpdir="/tmp/sweep6_${cfg}"
  mkdir -p "${tmpdir}"
  for f in "${outdir}"/*_score.txt; do
    [ -s "$f" ] || continue
    bn=$(basename "$f" _score.txt)
    cp "$f" "${tmpdir}/${bn}.out"
  done
  metrics=$(python3 predictor_metrics.py "${tmpdir}" 2>/dev/null || echo "FAIL")
  if [ "$metrics" = "FAIL" ]; then
    printf "%-20s %10s %10s %10s\n" "$cfg" "FAIL" "FAIL" "FAIL"
  else
    vfs=$(echo "$metrics" | python3 vfs.py 2>/dev/null || echo "FAIL")
    IFS=',' read -r ipc cpi epi mpi rest <<< "$metrics"
    mpki=$(python3 -c "print(f'{float(\"${mpi}\")*1000:.3f}')")
    printf "%-20s %10s %10s %10s\n" "$cfg" "$ipc" "$mpki" "$vfs"
  fi
done
