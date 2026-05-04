#!/usr/bin/env bash
set -uo pipefail

CONFIGS=(
  "BLK_N1_L4    blk_n1_l4"
  "BLK_N1_L8    blk_n1_l8"
  "BLK_N1_L16   blk_n1_l16"
  "BLK_N2_L4    blk_n2_l4"
  "BLK_N2_L8    blk_n2_l8"
  "BLK_N2_L16   blk_n2_l16"
  "BLK_N4_L4    blk_n4_l4"
  "BLK_N4_L8    blk_n4_l8"
  "BLK_N4_L16   blk_n4_l16"
)

for entry in "${CONFIGS[@]}"; do
  read -r pred outdir <<< "$entry"
  outpath="out/$outdir"

  if [[ -f "$outpath/aggregate_summary.txt" ]]; then
    mpki=$(grep 'Mean MPKI' "$outpath/aggregate_summary.txt" | awk '{print $NF}')
    echo "SKIP: $pred (already done) — Mean MPKI: ${mpki:-?}"
    continue
  fi

  echo "=== Building and running $pred → $outpath ==="
  rm -f build/cbp-monitor-eval
  if ! make eval-monitor EVAL_MONITOR_PRED="$pred" EVAL_MONITOR_OUT="$outpath" EVAL_MONITOR_JOBS=4 -j4 2>&1; then
    echo "  FAILED: $pred"
    continue
  fi

  mpki=$(grep 'Mean MPKI' "$outpath/aggregate_summary.txt" | awk '{print $NF}')
  echo "  DONE: $pred —   Mean MPKI:   ${mpki:-N/A}"
done

echo ""
echo "=== ALL DONE ==="
echo ""
echo "=== Summary ==="
printf "%-14s %s\n" "Config" "Mean MPKI"
printf "%-14s %s\n" "------" "---------"
for entry in "${CONFIGS[@]}"; do
  read -r pred outdir <<< "$entry"
  outpath="out/$outdir"
  if [[ -f "$outpath/aggregate_summary.txt" ]]; then
    mpki=$(grep 'Mean MPKI' "$outpath/aggregate_summary.txt" | awk '{print $NF}')
    printf "%-14s %s\n" "$pred" "${mpki:-N/A}"
  else
    printf "%-14s %s\n" "$pred" "MISSING"
  fi
done
