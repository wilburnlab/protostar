#!/bin/bash
# MS1 isotope-chromatogram extraction for the Counter-motivation figures.
#
# Drives `pipelines/30_extract_intermediates.py`: for the seed peptides (PROCAL +
# the prior notebooks' known-clear examples) across a dataset's non-ETD
# acquisitions, extracts each MS1 survey scan's M+0/M+1/M+2 isotope intensities +
# signed m/z error (every charge 1..4) over an RT window, with per-scan IIT/TIC.
# I/O-bound (streams each run's bundle-v3 peaks via protostar.bundles, one run at
# a time); the array shards the search-zip list so the tasks run concurrently.
# Needs v3 bundles (stage 10 after the 2026-11 reconversion): runs not yet
# reconverted are reported at the end, and a failed run makes the task exit 1. Resume-safe at shard
# granularity (a killed task just re-runs its [i::N] slice). Search zips extract
# to node-local $TMPDIR (not ESS) to avoid inode blow-up.
#
# Submit (Zolg2017, 16-way shard):
#   sbatch --account=<your-allocation> --array=0-15 \
#       --export=ALL,PROTOSTAR_DIR=$HOME/projects/protostar,DATASET=Zolg2017,N_SHARDS=16 \
#       pipelines/slurm/extract_intermediates.sbatch
# Match --array to N_SHARDS (0 .. N_SHARDS-1).
#
#SBATCH --job-name=protostar-extract-intermediates
# Set your OSC allocation — pass `--account=...` on submit (or uncomment + edit):
# #SBATCH --account=PASxxxx
#SBATCH --cluster=cardinal
#SBATCH --array=0-15
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=%x-%A_%a-%j.log

set -eo pipefail  # NOT -euo: conda activation trips on unbound vars

module load miniconda3/24.1.2-py310
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate constellation

# Pin the env's libstdc++ ahead of the spack/system one (module load shadows it →
# CXXABI_1.3.15 mismatch). The script self-heals too; belt and suspenders.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# --- knobs (override via --export) ---------------------------------------
PROTOSTAR_DIR="${PROTOSTAR_DIR:-$HOME/projects/protostar}"   # the OSC checkout
DATASET="${DATASET:-Zolg2017}"
N_SHARDS="${N_SHARDS:-16}"
PEPTIDES="${PEPTIDES:-seed}"     # procal | seed | all
RT_WINDOW="${RT_WINDOW:-60}"     # ± seconds around each PSM RT

idx="${SLURM_ARRAY_TASK_ID:-0}"

cd "$PROTOSTAR_DIR"
echo "extract_intermediates: dataset=$DATASET shard=$idx/$N_SHARDS peptides=$PEPTIDES host=$(hostname) tmpdir=${TMPDIR:-/tmp}"
python pipelines/30_extract_intermediates.py \
    --dataset "$DATASET" \
    --peptides "$PEPTIDES" \
    --rt-window "$RT_WINDOW" \
    --shard "$idx" --n-shards "$N_SHARDS"
echo "done shard $idx/$N_SHARDS"
