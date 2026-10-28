# Source this file from bash or zsh (macOS, Linux, WSL, Git Bash):  source scripts/env.sh
# It sets project paths and a plink2 wrapper. It does not change shell options,
# so it is safe to source from a script that has set -euo pipefail.

# Locate this file in either shell. The zsh-only expansion sits in a branch bash never takes.
if [ -n "${BASH_VERSION:-}" ]; then
  _env_src="${BASH_SOURCE[0]}"
elif [ -n "${ZSH_VERSION:-}" ]; then
  _env_src="${(%):-%x}"
else
  _env_src="$0"
fi
PROJECT_ROOT="$(cd "$(dirname "$_env_src")/.." && pwd)"
unset _env_src
export PROJECT_ROOT

export DATA="$PROJECT_ROOT/data"
export QC="$PROJECT_ROOT/results/qc"
export TASKS="$PROJECT_ROOT/results/tasks"
export VCF="${VCF:-$DATA/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz}"
export PANEL="${PANEL:-$DATA/integrated_call_samples_v3.20130502.ALL.panel}"

# Experiment size: 2 scenarios x N_SEEDS replicates = array indices 1..(2*N_SEEDS)
export N_SEEDS="${N_SEEDS:-10}"
export N_CAUSAL="${N_CAUSAL:-30}"

# Threads: Slurm's allocation if there is one, else all cores. getconf works on Linux and macOS.
THREADS="${SLURM_CPUS_PER_TASK:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
export THREADS

# plink2 lookup order: PLINK2 env var, a plink2 on PATH (cluster module or Windows plink2.exe),
# then the project's micromamba env. The lookup runs BEFORE the wrapper function is defined:
# once a function called plink2 exists, `command -v plink2` would find the function itself.
unset -f plink2 2>/dev/null || true
if [ -n "${PLINK2:-}" ]; then
  P2_MODE=bin; P2_BIN="$PLINK2"
elif _p2="$(command -v plink2 2>/dev/null)" && [ -n "$_p2" ]; then
  P2_MODE=bin; P2_BIN="$_p2"
elif command -v micromamba >/dev/null 2>&1 && micromamba run -n gwas plink2 --version >/dev/null 2>&1; then
  P2_MODE=mamba; P2_BIN=""
else
  P2_MODE=none; P2_BIN=""
fi
unset _p2
export P2_MODE P2_BIN

plink2() {
  case "$P2_MODE" in
    bin)   "$P2_BIN" --threads "$THREADS" "$@" ;;
    mamba) micromamba run -n gwas plink2 --threads "$THREADS" "$@" ;;
    *)     echo "plink2 not found: set PLINK2=/path/to/plink2 or create the env (micromamba create -f environment.yml)" >&2
           return 127 ;;
  esac
}

# Task id -> scenario + seed. Works for a Slurm array task or a local loop.
task_scenario() { if [ "$1" -le "$N_SEEDS" ]; then echo A; else echo B; fi; }
task_seed()     { echo $(( ($1 - 1) % N_SEEDS + 1 )); }
