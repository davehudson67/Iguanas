#!/usr/bin/env bash
set -euo pipefail
cd ~/Iguanas
mkdir -p logs

for chain in 1 2 3
do
  session="iguana_no_leafn_ucay_c${chain}"
  log="logs/sensitivity_no_leafn_ucay_chain${chain}.log"
  if tmux has-session -t "$session" 2>/dev/null
  then
    echo "$session already exists; not starting duplicate."
  else
    tmux new-session -d -s "$session" \
      "cd ~/Iguanas && CHAIN_ID=${chain} NITER=150000 NBURN=50000 THIN=10 Rscript 41_Model_Sensitivity_Exclude_LeafNorth_UCay_LONG_1CHAIN.R 2>&1 | tee ${log}"
    echo "Started $session -> $log"
  fi
done

tmux ls 2>/dev/null | grep "iguana_no_leafn_ucay" || true
