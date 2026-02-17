#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-partial}"        # partial | complete
OUTDIR="${2:-benchmark_results}"
BOOTSTRAP="${3:-200}"       # used only in complete mode

echo "[1/3] Running local ADAM (${MODE})"
Rscript scripts/benchmark/run_adam_version.R \
  --engine local \
  --mode "${MODE}" \
  --outdir "${OUTDIR}" \
  --bootstrap "${BOOTSTRAP}"

echo "[2/3] Running Bioconductor ADAM (${MODE})"
Rscript scripts/benchmark/run_adam_version.R \
  --engine bioc \
  --mode "${MODE}" \
  --outdir "${OUTDIR}" \
  --bootstrap "${BOOTSTRAP}"

echo "[3/3] Comparing outputs"
Rscript scripts/benchmark/compare_adam_versions.R \
  --mode "${MODE}" \
  --outdir "${OUTDIR}" \
  --top-n 20

echo "Done. See ${OUTDIR}/overview_${MODE}.csv and ${OUTDIR}/metrics_${MODE}.csv"
