#!/bin/bash -l
#SBATCH --job-name=gwas_qc
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0-2:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# module load plink2 R   # site-specific; uncomment on a cluster that uses modules
bash scripts/01_qc.sh
