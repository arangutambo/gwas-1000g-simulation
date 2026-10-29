#!/usr/bin/env bash
# Stage 4: LD score regression (Bulik-Sullivan et al. 2015) on the simulated scans.
#
# The report's own limitation says lambda cannot separate confounding from
# polygenicity. This stage tests whether the LD-score intercept can, using the
# simulation's known truth: scenario A has no environmental confounder,
# scenario B has one.
#
# Needs stages 1-2 to have run (results/qc and results/tasks/*/pheno.tsv).
# The LD score step takes about 45 minutes; it is skipped if its output exists.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$PROJECT_ROOT"

WORK="${1:-$PROJECT_ROOT/results/ldsc_work}"
mkdir -p "$WORK/scans" "$PROJECT_ROOT/results/tables"

echo "[ldsc] plink1 binary genotypes for the LD score reader"
plink2 --pfile "$QC/chr22_common" --make-bed --out "$WORK/chr22_common" > /dev/null

if [ -s "$WORK/ldscores.tsv" ]; then
  echo "[ldsc] reusing $WORK/ldscores.tsv"
else
  echo "[ldsc] computing LD scores (about 45 min)"
  Rscript R/ld_scores.R "$WORK/chr22_common" "$WORK/ldscores.tsv"
fi

# R/summarise_task.R deletes the full .glm.linear files because they are large,
# so regenerate them here from the saved phenotypes. Each scan takes under a second.
echo "[ldsc] re-running $((2 * N_SEEDS)) replicates x 2 models"
for S in A B; do
  for I in $(seq 1 "$N_SEEDS"); do
    PH="$TASKS/${S}_${I}/pheno.tsv"
    [ -s "$PH" ] || { echo "missing $PH - run stage 2 first" >&2; exit 1; }
    plink2 --pfile "$QC/chr22_common" --pheno "$PH" --pheno-name y \
      --glm allow-no-covars --out "$WORK/scans/${S}_${I}_naive" > /dev/null
    plink2 --pfile "$QC/chr22_common" --pheno "$PH" --pheno-name y \
      --covar "$QC/pcs.eigenvec" --covar-name PC1-PC10 --covar-variance-standardize \
      --glm hide-covar --out "$WORK/scans/${S}_${I}_pc" > /dev/null
  done
done

echo "[ldsc] regression"
Rscript R/ldsc.R "$WORK" "$PROJECT_ROOT/results"
echo "[ldsc] wrote results/tables/ldsc.csv and results/ldsc_summary.json"
