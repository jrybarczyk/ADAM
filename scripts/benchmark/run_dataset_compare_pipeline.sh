#!/usr/bin/env bash
set -euo pipefail

DATASET="${1:-airway}"      # airway | all | sce_mock
MODE="${2:-partial}"        # partial | complete
OUTDIR="${3:-benchmark_results}"
BOOTSTRAP="${4:-200}"
MAX_GENES="${5:-3000}"
PAIRING="${6:-first}"       # first | all

echo "[1/3] Running local ADAM (${DATASET}, ${MODE})"
Rscript scripts/benchmark/run_dataset_version.R \
  --engine local \
  --dataset "${DATASET}" \
  --mode "${MODE}" \
  --outdir "${OUTDIR}" \
  --bootstrap "${BOOTSTRAP}" \
  --max-genes "${MAX_GENES}" \
  --pairing "${PAIRING}"

echo "[2/3] Running Bioconductor ADAM (${DATASET}, ${MODE})"
Rscript scripts/benchmark/run_dataset_version.R \
  --engine bioc \
  --dataset "${DATASET}" \
  --mode "${MODE}" \
  --outdir "${OUTDIR}" \
  --bootstrap "${BOOTSTRAP}" \
  --max-genes "${MAX_GENES}" \
  --pairing "${PAIRING}"

echo "[3/3] Comparing outputs"
Rscript scripts/benchmark/compare_dataset_versions.R \
  --dataset "${DATASET}" \
  --mode "${MODE}" \
  --outdir "${OUTDIR}" \
  --top-n 20

echo "Done. See ${OUTDIR}/overview_${DATASET}_${MODE}.csv and ${OUTDIR}/metrics_${DATASET}_${MODE}.csv"
