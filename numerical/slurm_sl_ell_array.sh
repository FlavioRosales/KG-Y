#!/usr/bin/env bash
#SBATCH --job-name=sl_ell
#SBATCH --output=logs/sl_ell_%A_%a.out
#SBATCH --error=logs/sl_ell_%A_%a.err
#SBATCH --array=0-80
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00

set -euo pipefail

PREFIX=${PREFIX:-sl_part}
ELL=${SLURM_ARRAY_TASK_ID}

mkdir -p logs

./sl_solve_ell_parallel "$ELL" "$ELL" "$PREFIX"
