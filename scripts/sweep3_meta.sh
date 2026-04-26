#!/usr/bin/env bash
set -uo pipefail

CONFIGS=(
  "S3_MW1_MC512   s3_mw1_mc512"
  "S3_MW1_MC1024  s3_mw1_mc1024"
  "S3_MW1_MC2048  s3_mw1_mc2048"
  "S3_MW1_MC4096  s3_mw1_mc4096"
  "S3_MW2_MC512   s3_mw2_mc512"
  "S3_MW2_MC1024  s3_mw2_mc1024"
  "S3_MW2_MC2048  s3_mw2_mc2048"
  "S3_MW2_MC4096  s3_mw2_mc4096"
  "S3_MW4_MC512   s3_mw4_mc512"
  "S3_MW4_MC1024  s3_mw4_mc1024"
  "S3_MW4_MC2048  s3_mw4_mc2048"
  "S3_MW4_MC4096  s3_mw4_mc4096"
)

for entry in "${CONFIGS[@]}"; do
  read -r pred outdir <<< "$entry"
  outpath="out/$outdir"

  if [[ -f "$outpath/aggregate_summary.txt" ]]; then
    mpki=$(awk '/MPKI/ {print $2; exit}' "$outpath/aggregate_summary.txt" 2>/dev/null)
    echo "SKIP: $pred (already done) — Mean MPKI: ${mpki:-?}"
    continue
  fi

  echo "=== Building and running $pred → $outpath ==="
  rm -f build/cbp-monitor-eval
  if ! make eval-monitor EVAL_MONITOR_PRED="$pred" EVAL_MONITOR_OUT="$outpath" -j$(nproc) 2>&1; then
    echo "  FAILED: $pred"
    continue
  fi

  # Extract mean MPKI from aggregate summary
  mpki=$(awk '/MPKI/ {print $2; exit}' "$outpath/aggregate_summary.txt" 2>/dev/null)
  echo "  DONE: $pred —   Mean MPKI:   ${mpki:-N/A}"
done

echo ""
echo "=== ALL DONE ==="
