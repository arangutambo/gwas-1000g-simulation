#!/usr/bin/env bash
# Stage 3: combine per-task summaries and figures inputs.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$PROJECT_ROOT"
mkdir -p results/tables
Rscript R/aggregate.R
