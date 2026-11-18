#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
MAX_PARALLEL="${MAX_PARALLEL:-9}"
NITER="${NITER:-150000}"
NBURN="${NBURN:-50000}"
THIN="${THIN:-10}"

mkdir -p logs

active_count() {
  tmux ls 2>/dev/null | grep '^iguana_' | wc -l || true
}

wait_for_slot() {
  while [ "$(active_count)" -ge "${MAX_PARALLEL}" ]; do sleep 30; done
}

start_if_missing() {
  outfile="$1"; session="$2"; cmd="$3"; logfile="$4"
  if [ -s "${outfile}" ]; then
    echo "SKIP complete: ${outfile}"
    return
  fi
  wait_for_slot
  if tmux has-session -t "${session}" 2>/dev/null; then
    echo "SKIP active: ${session}"
  else
    tmux new-session -d -s "${session}" "cd '${ROOT}' && ${cmd} 2>&1 | tee '${logfile}'"
    echo "START ${session} -> ${logfile}"
  fi
}

# Audited sensitivities
for variant in ANY_FEED DENSITY_PARTIAL MATTHEW_ACUTE; do
  case "${variant}" in
    ANY_FEED) tag="any_feed" ;;
    DENSITY_PARTIAL) tag="density_partial" ;;
    MATTHEW_ACUTE) tag="matthew_acute" ;;
  esac
  for chain in 1 2 3; do
    outfile="ModelRuns/Sensitivity_${variant}_AUDITED_LONG/Model_Sensitivity_${variant}_AUDITED_LONG_CHAIN${chain}.rds"
    start_if_missing "${outfile}" "iguana_${tag}_c${chain}" \
      "MODEL_VARIANT='${variant}' CHAIN_ID='${chain}' NITER='${NITER}' NBURN='${NBURN}' THIN='${THIN}' Rscript 12_Model_Sensitivity_LONG_1CHAIN.R" \
      "logs/${tag}_chain${chain}.log"
  done
done

# Targeted Leaf North + U Cay exclusion
for chain in 1 2 3; do
  outfile="ModelRuns/Sensitivity_Exclude_LeafNorth_UCay/Model_Sensitivity_Exclude_LeafNorth_UCay_CHAIN${chain}.rds"
  start_if_missing "${outfile}" "iguana_no_leafn_ucay_c${chain}" \
    "CHAIN_ID='${chain}' NITER='${NITER}' NBURN='${NBURN}' THIN='${THIN}' Rscript 41_Model_Sensitivity_Exclude_LeafNorth_UCay_LONG_1CHAIN.R" \
    "logs/sensitivity_no_leafn_ucay_chain${chain}.log"
done

# Leave-one-site-out fits
sites=(Bitter FFRC Gaulin LeafNorth LeafSouth NAdder Noddy UCay WhiteB)
for site in "${sites[@]}"; do
  tag="$(echo "${site}" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]')"
  for chain in 1 2 3; do
    outfile="ModelRuns/LOO_${tag}/Model_LOO_${tag}_CHAIN${chain}.rds"
    start_if_missing "${outfile}" "iguana_loo_${tag}_c${chain}" \
      "OMIT_SITE='${site}' CHAIN_ID='${chain}' NITER='${NITER}' NBURN='${NBURN}' THIN='${THIN}' Rscript 45_Model_LOO_FinalAudited_LONG_1CHAIN.R" \
      "logs/loo_${tag}_chain${chain}.log"
  done
done

echo
echo "============================================================"
echo "MISSING-ONLY ROBUSTNESS QUEUE SUBMITTED"
echo "MAX_PARALLEL=${MAX_PARALLEL}"
echo "Existing completed chain files were skipped."
echo "============================================================"
