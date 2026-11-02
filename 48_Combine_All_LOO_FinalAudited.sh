#!/usr/bin/env bash
set -euo pipefail

sites=(Bitter FFRC Gaulin LeafNorth LeafSouth NAdder Noddy UCay WhiteB)

for site in "${sites[@]}"; do
  tag="$(echo "${site}" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]')"
  echo
  echo "============================================================"
  echo "COMBINING LOO: ${site}"
  echo "============================================================"
  OMIT_SITE="${site}" Rscript 46_Combine_LOO_FinalAudited_Chains.R 2>&1 | tee "logs/combine_loo_${tag}.log"
done

echo
echo "All nine leave-one-site-out fits combined."
