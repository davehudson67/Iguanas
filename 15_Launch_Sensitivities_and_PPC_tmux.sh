#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"

NITER="${NITER:-150000}"
NBURN="${NBURN:-50000}"
THIN="${THIN:-10}"
N_PPC="${N_PPC:-1000}"
PPC_CORES="${PPC_CORES:-4}"

mkdir -p logs

start_chain () {
  local variant="$1"
  local chain="$2"
  local short="$3"
  local session="iguana_${short}_c${chain}"
  local logfile="logs/${short}_chain${chain}.log"

  if tmux has-session -t "${session}" 2>/dev/null; then
    echo "${session}: already exists; not starting duplicate."
    return
  fi

  tmux new-session -d -s "${session}" \
    "cd '${ROOT}' && MODEL_VARIANT='${variant}' CHAIN_ID='${chain}' NITER='${NITER}' NBURN='${NBURN}' THIN='${THIN}' Rscript 12_Model_Sensitivity_LONG_1CHAIN.R 2>&1 | tee '${logfile}'"

  echo "Started ${session} -> ${logfile}"
}


for chain in 1 2 3; do
  start_chain "ANY_FEED" "${chain}" "anyfeed"
done


for chain in 1 2 3; do
  start_chain "MATTHEW_ACUTE" "${chain}" "matthewacute"
done


for chain in 1 2 3; do
  start_chain "DENSITY_PARTIAL" "${chain}" "densitypartial"
done


PPC_SESSION="iguana_ppc_primary"
PPC_LOG="logs/ppc_primary.log"

if tmux has-session -t "${PPC_SESSION}" 2>/dev/null; then
  echo "${PPC_SESSION}: already exists; not starting duplicate."
else
  tmux new-session -d -s "${PPC_SESSION}" \
    "cd '${ROOT}' && N_PPC='${N_PPC}' PPC_CORES='${PPC_CORES}' Rscript 14_PosteriorPredictiveChecks_Primary.R 2>&1 | tee '${PPC_LOG}'"

  echo "Started ${PPC_SESSION} -> ${PPC_LOG}"
fi


echo
echo "============================================================"
echo "STARTED SENSITIVITIES + PPC"
echo "============================================================"
echo "MCMC per chain: NITER=${NITER}, NBURN=${NBURN}, THIN=${THIN}"
echo "PPC: N_PPC=${N_PPC}, PPC_CORES=${PPC_CORES}"
echo
echo "Active Iguana tmux sessions:"
tmux ls 2>/dev/null | grep '^iguana_' || true
echo
echo "Useful status command:"
echo "  tmux ls | grep iguana_"
echo
echo "Tail a log, e.g.:"
echo "  tail -f logs/anyfeed_chain1.log"
echo "============================================================"
