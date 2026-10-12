#!/bin/bash
#SBATCH --job-name=L4_SPIN
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).

set -euo pipefail
source "${CONDA_SETUP:-$HOME/miniconda3/etc/profile.d/conda.sh}"
conda activate "${EAT_ENV:-eat}"
cd "${SLURM_SUBMIT_DIR:?Submit from the experiment directory}"
[[ -s gotm.yaml && -s fabm.yaml && -x runrun.sh ]] || exit 1
[[ "${SLURM_NTASKS:?}" -eq 1 ]] || { echo 'Incorrect task count' >&2; exit 1; }
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1 GOTM_BIN=eat-gotm

echo "L4 job=${SLURM_JOB_ID} tasks=${SLURM_NTASKS} cwd=$(pwd)"
./runrun.sh
