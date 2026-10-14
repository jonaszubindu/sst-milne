#!/bin/bash
#SBATCH --partition=aiub
#SBATCH --job-name=me_quicklook
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
##SBATCH --reservation=aiub-prime
#SBATCH --output=me_quicklook_%j.out
#
# Quicklook overviews + movies for inversion results, on a compute node
# (2900x2900 maps with tens of scans are too heavy for a login node).
# Only cube folders whose movie is missing or older than their newest scan
# are (re)made.
#
#   sbatch /path/to/sst-milne/slurm/me_quicklook.sbatch /path/to/me_out [/path/to/me_out_17 ...]
#
# Same environment variables as me_invert_array.sbatch (SST_MILNE_DIR,
# CONDA_ENV); needs ffmpeg in the environment (conda-forge ffmpeg).
set -eo pipefail
[[ $# -ge 1 ]] || { echo "usage: sbatch me_quicklook.sbatch OUT_DIR [OUT_DIR ...]"; exit 1; }

SST_MILNE_DIR=${SST_MILNE_DIR:-/unibe/scratch/zbindenj/sst-milne}
CONDA_ENV=${CONDA_ENV:-milne}

echo "Quicklooks on $SLURM_JOB_NODELIST, start $(date): $*"
cd "$SLURM_SUBMIT_DIR"
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
set -u

python -u "$SST_MILNE_DIR/me_quicklook.py" "$@" --movie --skip-existing
echo "End $(date)"
