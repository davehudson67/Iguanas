#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
MAX_PARALLEL="${MAX_PARALLEL:-6}"
NITER="${NITER:-150000}"
NBURN="${NBURN:-50000}"
THIN="${THIN:-10}"

mkdir -p logs

wait_for_slot() {
  while true; do
    n=$(tmux ls 2>/dev/null | grep '^iguana_' | wc -l || true)
    if [ "${n}" -lt "${MAX_PARALLEL}" ]; then
      break
    fi
    sleep 30
  done
}

start_tmux_job() {
  session="$1"
  cmd="$2"
  logfile="$3"
  wait_for_slot
  if tmux has-session -t "${session}" 2>/dev/null; then
    echo "${session}: already exists; skipping."
  else
    tmux new-session -d -s "${session}" "cd '${ROOT}' && ${cmd} 2>&1 | tee '${logfile}'"
    echo "Started ${session} -> ${logfile}"
  fi
}

# 1) Audited sensitivity variants
for variant in ANY_FEED DENSITY_PARTIAL MATTHEW_ACUTE; do
  tag=$(echo "${variant}" | tr '[:upper:]' '[:lower:]')
  for chain in 1 2 3; do
    start_tmux_job "iguana_${tag}_c${chain}" \
      "VARIANT='${variant}' CHAIN_ID='${chain}' NITER='${NITER}' NBURN='${NBURN}' THIN='${THIN}' Rscript 12_Model_Sensitivity_LONG_1CHAIN.R" \
      "logs/${tag}_chain${chain}.log"
  done
done

# 2) Leaf North + U Cay exclusion
for chain in 1 2 3; do
  start_tmux_job "iguana_no_leafn_ucay_c${chain}" \
    "CHAIN_ID='${chain}' NITER='${NITER}' NBURN='${NBURN}' THIN='${THIN}' Rscript 41_Model_Sensitivity_Exclude_LeafNorth_UCay_LONG_1CHAIN.R" \
    "logs/sensitivity_no_leafn_ucay_chain${chain}.log"
done

# 3) Leave-one-site-out fits
sites=(Bitter FFRC Gaulin LeafNorth LeafSouth NAdder Noddy UCay WhiteB)
for site in "${sites[@]}"; do
  tag="$(echo "${site}" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]')"
  for chain in 1 2 3; do
    start_tmux_job "iguana_loo_${tag}_c${chain}" \
      "OMIT_SITE='${site}' CHAIN_ID='${chain}' NITER='${NITER}' NBURN='${NBURN}' THIN='${THIN}' Rscript 45_Model_LOO_FinalAudited_LONG_1CHAIN.R" \
      "logs/loo_${tag}_chain${chain}.log"
  done
done

echo
echo "All jobs have been submitted through the throttle."
echo "MAX_PARALLEL=${MAX_PARALLEL}"
