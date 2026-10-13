#!/bin/bash
#SBATCH --partition=aiub
#SBATCH --job-name=me_invert
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
##SBATCH --reservation=aiub-prime
#SBATCH --output=me_%A_%a.out
#
# Milne-Eddington inversion of SST/CRISP cubes as a Slurm job ARRAY: one
# cube per array task, each on a whole node. pyMilne parallelises with
# OpenMP threads inside one process (not MPI), so one cube cannot span
# nodes; several cubes run side by side instead. Scans within a cube stay
# sequential (each starts from the previous scan's result).
#
# Usage (from the folder where the logs should go):
#   find /path/to/SST_data -name 'nb_6302*_im.fits' | sort > cubes.txt
#   sbatch --array=0-$(( $(wc -l < cubes.txt) - 1 )) \
#       /path/to/sst-milne/slurm/me_invert_array.sbatch cubes.txt /path/to/me_out [me_invert options]
# e.g. options: --single-pass   --scans 0-9   --bin 2
# Throttle to N nodes at a time with --array=0-K%N.
#
# Paths (override with environment variables when submitting, e.g.
#   sbatch --export=ALL,PYMILNE_DIR=/x/pyMilne ...):
#   PYMILNE_DIR    pyMilne checkout, built with slurm/build_pymilne.sbatch
#   SST_MILNE_DIR  this repository
#   CONDA_ENV      conda environment with numpy, scipy, astropy, matplotlib
set -eo pipefail

LIST=${1:?usage: sbatch --array=0-N me_invert_array.sbatch cubes.txt OUT_DIR [options]}
OUT=${2:?output directory missing}
shift 2
EXTRA=("$@")

PYMILNE_DIR=${PYMILNE_DIR:-/unibe/scratch/zbindenj/pyMilne}
SST_MILNE_DIR=${SST_MILNE_DIR:-/unibe/scratch/zbindenj/sst-milne}
CONDA_ENV=${CONDA_ENV:-milne}

CUBE=$(sed -n "$(( SLURM_ARRAY_TASK_ID + 1 ))p" "$LIST")
[[ -n "$CUBE" && -f "$CUBE" ]] || { echo "ERROR: no cube for task $SLURM_ARRAY_TASK_ID in $LIST"; exit 1; }

echo "Job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} on $SLURM_JOB_NODELIST, start $(date)"
echo "Cube:   $CUBE"
echo "Output: $OUT"

cd "$SLURM_SUBMIT_DIR"
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

# All CPUs of the node go to pyMilne's OpenMP threads; numpy single-threaded.
NT=${SLURM_CPUS_ON_NODE:-$(nproc)}
export OMP_NUM_THREADS=$NT
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYMILNE_DIR
echo "Threads: $NT"
set -u

[[ -n $(ls "$PYMILNE_DIR"/pyMilne*.so 2>/dev/null) ]] || {
    echo "ERROR: pyMilne is not built in $PYMILNE_DIR (run slurm/build_pymilne.sbatch)"; exit 1; }

# pyMilne's progress counter (carriage returns) is filtered out of the log;
# the per-scan lines, warnings and errors stay.
python -u "$SST_MILNE_DIR/me_invert.py" "$CUBE" --out "$OUT" --per-cube-dirs \
    --nthreads "$NT" --quicklook "${EXTRA[@]}" 2>&1 \
    | stdbuf -oL tr '\r' '\n' | grep --line-buffered -v 'InvertMany'

echo "End $(date)"
