#!/usr/bin/env bash
set -uo pipefail

CONFIGS=(
  "BLK_N2_L32   blk_n2_l32"
  "BLK_N2_L64   blk_n2_l64"
  "BLK_N4_L32   blk_n4_l32"
  "BLK_N4_L64   blk_n4_l64"
  "BLK_N4_L128  blk_n4_l128"
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
