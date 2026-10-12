#!/bin/bash
#SBATCH --job-name=L4_V3_500_T090_F010
#SBATCH --time=12:00:00
#SBATCH --ntasks=501
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --exclusive
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).

set -euo pipefail
source "${CONDA_SETUP:-$HOME/miniconda3/etc/profile.d/conda.sh}"
conda activate "${EAT_ENV:-eat}"
cd "${SLURM_SUBMIT_DIR:?Submit from the experiment directory}"
[[ -s gotm.yaml && -s fabm.yaml && -x runrun.sh ]] || exit 1
[[ "${SLURM_NTASKS:?}" -eq 501 ]] || { echo 'Incorrect task count' >&2; exit 1; }
source ./experiment.env
python preflight_l4.py

echo "L4 job=${SLURM_JOB_ID} tasks=${SLURM_NTASKS} cwd=$(pwd)"
./runrun.sh
