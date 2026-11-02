#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
NITER="${NITER:-150000}"
NBURN="${NBURN:-50000}"
THIN="${THIN:-10}"

mkdir -p logs

sites=(Bitter FFRC Gaulin LeafNorth LeafSouth NAdder Noddy UCay WhiteB)

for site in "${sites[@]}"; do
  tag="$(echo "${site}" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]')"

  for chain in 1 2 3; do
    session="iguana_loo_${tag}_c${chain}"
    logfile="logs/loo_${tag}_chain${chain}.log"

    if tmux has-session -t "${session}" 2>/dev/null; then
      echo "${session}: already exists; not starting duplicate."
    else
      tmux new-session -d -s "${session}" \
        "cd '${ROOT}' && OMIT_SITE='${site}' CHAIN_ID='${chain}' NITER='${NITER}' NBURN='${NBURN}' THIN='${THIN}' Rscript 45_Model_LOO_FinalAudited_LONG_1CHAIN.R 2>&1 | tee '${logfile}'"
      echo "Started ${session} -> ${logfile}"
    fi
  done
done

echo
echo "============================================================"
echo "STARTED FINAL-AUDITED LEAVE-ONE-SITE-OUT BATCH"
echo "============================================================"
echo "Sites: ${sites[*]}"
echo "Chains per site: 3"
echo "Total LOO chains: 27"
echo "MCMC per chain: NITER=${NITER}, NBURN=${NBURN}, THIN=${THIN}"
echo
echo "Active LOO tmux sessions:"
tmux ls 2>/dev/null | grep '^iguana_loo_' || true
echo "============================================================"
