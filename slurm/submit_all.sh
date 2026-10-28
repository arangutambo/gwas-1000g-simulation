#!/usr/bin/env bash
# Submit QC -> array -> aggregate as a dependency chain. Run from the project root on a Slurm cluster.
set -euo pipefail
mkdir -p logs
qc=$(sbatch --parsable slurm/01_qc.slurm)
arr=$(sbatch --parsable --dependency=afterok:"$qc" slurm/02_task_array.slurm)
agg=$(sbatch --parsable --dependency=afterok:"$arr" slurm/03_aggregate.slurm)
echo "submitted: qc=$qc array=$arr aggregate=$agg"
