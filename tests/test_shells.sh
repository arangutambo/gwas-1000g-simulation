#!/bin/sh
# env.sh must resolve to the same paths and task mapping under bash and zsh.
cd "$(dirname "$0")/.." || exit 1
probe='source scripts/env.sh; source scripts/env.sh; echo "$PROJECT_ROOT|$THREADS|$(task_scenario 3)$(task_seed 3)|$(task_scenario 14)$(task_seed 14)|plink2=$P2_MODE"'
b=$(bash -c "$probe") || { echo "bash failed"; exit 1; }
if command -v zsh >/dev/null 2>&1; then
  z=$(zsh -c "$probe") || { echo "zsh failed"; exit 1; }
  [ "$b" = "$z" ] || { echo "MISMATCH: bash=[$b] zsh=[$z]"; exit 1; }
  echo "bash and zsh agree: $b"
else
  echo "zsh not installed; bash only: $b"
fi
