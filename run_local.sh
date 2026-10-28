#!/usr/bin/env bash
# Run the whole pipeline without Slurm (macOS, Linux, WSL, Git Bash).
#   bash run_local.sh            # everything
#   bash run_local.sh tasks 1 3  # only array indices 1..3
set -euo pipefail
cd "$(dirname "$0")"
source scripts/env.sh
total=$(( 2 * N_SEEDS ))
case "${1:-all}" in
  qc)    bash scripts/01_qc.sh ;;
  tasks) for i in $(seq "${2:-1}" "${3:-$total}"); do bash scripts/run_task.sh "$i"; done ;;
  agg)   bash scripts/03_aggregate.sh ;;
  all)   bash scripts/01_qc.sh
         for i in $(seq 1 "$total"); do bash scripts/run_task.sh "$i"; done
         bash scripts/03_aggregate.sh ;;
  *)     echo "usage: run_local.sh [all|qc|tasks FROM TO|agg]" >&2; exit 2 ;;
esac
