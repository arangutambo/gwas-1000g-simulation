#!/usr/bin/env bash
# Stage 2: one simulation replicate. Same entry point for a Slurm array task and a local run.
#   Slurm:  SLURM_ARRAY_TASK_ID is set.   Local:  bash scripts/run_task.sh 7
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd "$PROJECT_ROOT"

ID="${SLURM_ARRAY_TASK_ID:-${1:?usage: run_task.sh TASK_ID}}"
SCEN="$(task_scenario "$ID")"
SEED="$(task_seed "$ID")"
OUT="$TASKS/${SCEN}_${SEED}"
mkdir -p "$OUT"
echo "[task $ID] scenario=$SCEN seed=$SEED -> $OUT"

# 1. Choose causal variants (seeded, common SNPs only) and draw effect sizes.
Rscript R/pick_causal.R "$SEED" "$N_CAUSAL" "$OUT"

# 2. Pull causal dosages out of the pgen.
plink2 --pfile "$QC/chr22_common" --extract "$OUT/causal_ids.txt" --export A --out "$OUT/causal"

# 3. Build the phenotype for this scenario.
Rscript R/make_phenotype.R "$SCEN" "$SEED" "$OUT"

# 4. Association scans: naive, and adjusted for 10 PCs.
plink2 --pfile "$QC/chr22_common" --pheno "$OUT/pheno.tsv" --pheno-name y \
  --glm allow-no-covars --out "$OUT/naive"
plink2 --pfile "$QC/chr22_common" --pheno "$OUT/pheno.tsv" --pheno-name y \
  --covar "$QC/pcs.eigenvec" --covar-name PC1-PC10 --covar-variance-standardize \
  --glm hide-covar --out "$OUT/pcadj"

# 5. Score against the truth and estimate SNP heritability.
Rscript R/summarise_task.R "$SCEN" "$SEED" "$OUT"
