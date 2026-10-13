#!/bin/bash
#SBATCH --partition=aiub
#SBATCH --job-name=build_pymilne
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
##SBATCH --reservation=aiub-prime
#SBATCH --output=build_pymilne_%j.out
#
# Build pyMilne ON A COMPUTE NODE: its setup.py compiles with -march=native,
# so a build on a different CPU type (e.g. the login node) can crash with
# "illegal instruction" on the nodes. Run once (and again if the inversion
# nodes are of another CPU generation).
#
#   sbatch /path/to/sst-milne/slurm/build_pymilne.sbatch
#
# Same environment variables as me_invert_array.sbatch. If the system gcc
# is older than 10, load a newer one first (e.g. add `module load GCC`).
set -eo pipefail

PYMILNE_DIR=${PYMILNE_DIR:-/unibe/scratch/zbindenj/pyMilne}
SST_MILNE_DIR=${SST_MILNE_DIR:-/unibe/scratch/zbindenj/sst-milne}
CONDA_ENV=${CONDA_ENV:-milne}

echo "Build on $SLURM_JOB_NODELIST, start $(date)"
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"
set -u

bash "$SST_MILNE_DIR/pymilne_build/build_linux.sh" "$PYMILNE_DIR"
echo "End $(date)"
