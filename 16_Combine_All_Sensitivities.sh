#!/usr/bin/env bash
set -euo pipefail

mkdir -p logs

for variant in ANY_FEED MATTHEW_ACUTE DENSITY_PARTIAL; do
  tag="$(echo "${variant}" | tr '[:upper:]' '[:lower:]')"

  echo
  echo "============================================================"
  echo "COMBINING ${variant}"
  echo "============================================================"

  MODEL_VARIANT="${variant}" \
    Rscript 13_Combine_Sensitivity_LONG_Chains.R \
    2>&1 | tee "logs/combine_${tag}.log"
done

echo
echo "All three sensitivity variants combined."
