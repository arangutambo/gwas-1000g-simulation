#!/bin/bash -l
#SBATCH --job-name=gwas_task
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=0-1:00
#SBATCH --array=1-20%10
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# module load plink2 R
bash scripts/run_task.sh
