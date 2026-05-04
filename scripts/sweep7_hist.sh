#!/bin/bash
# Sweep 7: Structural accuracy improvements
# Adaptive sec-tag, tag widths, graded tags, history lengths

set -euo pipefail

CONFIGS=(
  "S3_MW4_MC2048"   # baseline
  "S7_Adapt64"      # adaptive sec-tag thresh=64
  "S7_Adapt96"      # adaptive sec-tag thresh=96
  "S7_Adapt128"     # adaptive sec-tag thresh=128
  "S7_Tag12"        # uniform 12-bit tags
  "S7_GT13_9"       # graded tags 13->9
  "S7_GT14_8"       # graded tags 14->8
  "S7_GT12_10"      # graded tags 12->10
  "S7_H6_300"       # history 6-300
  "S7_H5_150"       # history 5-150
  "S7_H10_250"      # history 10-250
)

JOBS="${1:-4}"

echo "=== Sweep 7: Structural Accuracy (${#CONFIGS[@]} configs, ${JOBS} jobs) ==="
echo ""

for cfg in "${CONFIGS[@]}"; do
  outdir="out/sweep7_${cfg}"
  echo "--- Building and running: ${cfg} ---"
  make eval-monitor EVAL_MONITOR_PRED="${cfg}" EVAL_MONITOR_OUT="${outdir}" EVAL_MONITOR_JOBS="${JOBS}" -B 2>&1 | grep -E "done:|ERROR|Running"
  echo ""
done

echo ""
echo "=== Results ==="
printf "%-20s %10s %10s %10s\n" "Config" "IPC" "MPKI" "VFS"
printf "%-20s %10s %10s %10s\n" "------" "---" "----" "---"

for cfg in "${CONFIGS[@]}"; do
  outdir="out/sweep7_${cfg}"
  tmpdir="/tmp/sweep7_${cfg}"
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
