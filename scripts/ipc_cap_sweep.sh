#!/bin/bash
# IPC cap parameter sweep — builds separate binaries and runs quick_eval for each config.
set -euo pipefail

TRACE_DIR="${1:-./traces}"
BASE_OUT="out/ipc_cap_sweep"
JOBS="${2:-$(nproc)}"

# Config format: label WIN_BITS MPKI_SHIFT MAX_INJ
CONFIGS=(
  "w128_m4    7   8  32"
  "w128_m2    7   9  32"
  "w128_m1    7  10  32"
  "w256_m4    8   8  32"
  "w256_m2    8   9  32"
  "w256_m1    8  10  32"
  "w512_m2    9   9  32"
  "w512_m1    9  10  32"
)

mkdir -p "${BASE_OUT}"

for cfg in "${CONFIGS[@]}"; do
  read -r label wb ms mi <<< "$cfg"
  echo "=== Building: ${label} (WIN_BITS=${wb} MPKI_SHIFT=${ms} MAX_INJ=${mi}) ==="
  FLAGS="-DCAP_WIN_BITS=${wb} -DCAP_MPKI_SHIFT=${ms} -DCAP_MAX_INJ=${mi}"
  BIN="build/cbp-cap-${label}"
  make cbp EXTRA_CBP_FLAGS="${FLAGS}" 2>&1 | tail -1
  # Rename binary to unique name
  BUILT=$(ls -t build/cbp-* | head -1)
  cp "${BUILT}" "${BIN}"

  echo "=== Evaluating: ${label} ==="
  OUT_DIR="${BASE_OUT}/${label}"
  scripts/quick_eval.sh "./${BIN}" "${TRACE_DIR}" "${OUT_DIR}" "${JOBS}" 2>&1 | tee "${BASE_OUT}/${label}.log"
  echo ""
done

echo "=== All configs done ==="
echo ""
echo "Summary:"
printf "%-12s %10s %10s %10s %10s\n" "Config" "IPC" "CPI" "EPI" "VFS"
printf "%-12s %10s %10s %10s %10s\n" "------" "---" "---" "---" "---"
for cfg in "${CONFIGS[@]}"; do
  read -r label _ _ _ <<< "$cfg"
  OUT_DIR="${BASE_OUT}/${label}"
  if [[ -d "${OUT_DIR}" ]]; then
    METRICS=$(python3 predictor_metrics.py "${OUT_DIR}" 2>/dev/null || echo "0,0,0,0,0,0,0,0")
    VFS=$(echo "${METRICS}" | python3 vfs.py 2>/dev/null || echo "ERR")
    IFS=',' read -r ipc cpi epi mpi rest <<< "${METRICS}"
    printf "%-12s %10s %10s %10s %10s\n" "${label}" "${ipc}" "${cpi}" "${epi}" "${VFS}"
  fi
done
