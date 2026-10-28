#!/bin/bash -l
#SBATCH --job-name=gwas_agg
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-0:20
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
bash scripts/03_aggregate.sh
