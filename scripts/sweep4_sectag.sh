#!/usr/bin/env bash
set -uo pipefail

CONFIGS=(
  "S4_None       s4_none"
  "S4_Floor4     s4_floor4"
  "S4_Floor7     s4_floor7"
  "S4_Floor10    s4_floor10"
  "S4_Ceil4      s4_ceil4"
  "S4_Ceil7      s4_ceil7"
  "S4_Ceil10     s4_ceil10"
  "S4_Press256   s4_press256"
  "S4_Press512   s4_press512"
  "S4_Press768   s4_press768"
  "S4_Acc256     s4_acc256"
  "S4_Acc512     s4_acc512"
  "S4_Acc768     s4_acc768"
  "S4_PGF7_512  s4_pgf7_512"
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
  if ! make eval-monitor EVAL_MONITOR_PRED="$pred" EVAL_MONITOR_OUT="$outpath" -j$(nproc) 2>&1; then
    echo "  FAILED: $pred"
    continue
  fi

  mpki=$(grep 'Mean MPKI' "$outpath/aggregate_summary.txt" | awk '{print $NF}')
  echo "  DONE: $pred —   Mean MPKI:   ${mpki:-N/A}"
done

echo ""
echo "=== ALL DONE ==="
