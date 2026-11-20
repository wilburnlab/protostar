#!/bin/bash
# Stage 10 — convert a dataset's .raw into Constellation parquet bundles (v3) on OSC.
#
# One pass writes proc/<dataset>/<stem>/ (centroid peaks) + <stem>.profile/ (MS1
# profile). Pre-v3 bundles in proc/<dataset>/{centroid,profile}/ are not touched.
#
# One extended Cardinal allocation per dataset, high worker count. Resume-safe:
# a wall-time kill just re-runs and skips finished bundles, so re-submitting is
# the recovery path. For very large datasets, the array-job escape hatch fans
# the work out — submit with `--array=0-(K-1)` and set N_SHARDS=K below.
#
# Submit (whole dataset):
#   sbatch --export=ALL,DATASET=Zolg2017 pipelines/slurm/convert.sbatch
# Submit (pilot: named files, no precursor-major copy):
#   sbatch --export=ALL,DATASET=Zolg2017,EXTRA_ARGS="--raw-name <stem> --no-precursor-major" \
#       pipelines/slurm/convert.sbatch
# Submit (sharded array of 20):
#   sbatch --array=0-19 --export=ALL,DATASET=Wilhelm2021,N_SHARDS=20 pipelines/slurm/convert.sbatch
#
#SBATCH --job-name=protostar-convert
# Set your OSC allocation — uncomment + edit, or pass `--account=...` on the sbatch command line
# (otherwise SLURM uses your default project):
# #SBATCH --account=PASxxxx
#SBATCH --cluster=cardinal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16   # convert is ESS-I/O-bound: 16 workers (see osc_runbook)
#SBATCH --time=48:00:00
#SBATCH --output=%x-%A_%a-%j.log

set -eo pipefail  # NOT -euo: conda activation trips on unbound vars

module load miniconda3/24.1.2-py310
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate constellation

# Pin the env's libstdc++ ahead of the spack/system one on the loader path, else
# sqlite3 -> libicui18n -> libstdc++ dies with a CXXABI_1.3.15 mismatch (the env's
# lib has it; `module load miniconda3` puts the older system lib first). Spawn
# workers inherit this, so the parallel converters resolve correctly too.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# --- knobs (override via --export) ---------------------------------------
PROTOSTAR_DIR="${PROTOSTAR_DIR:-$HOME/protostar}"   # repo on OSC home (rsynced)
DATASET="${DATASET:-all}"
EXTRA_ARGS="${EXTRA_ARGS:-}"                         # e.g. --raw-name X --no-precursor-major
WORKERS="${WORKERS:-${SLURM_CPUS_PER_TASK:-8}}"
# Optional: explicit Thermo CommonCore DLL dir (else auto-discovered in the env).
# export CONSTELLATION_THERMO_HOME=/path/to/thermo/dlls

# Array-job sharding: present only when submitted with --array.
SHARD_ARGS=()
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && -n "${N_SHARDS:-}" ]]; then
    SHARD_ARGS=(--shard "${SLURM_ARRAY_TASK_ID}" --n-shards "${N_SHARDS}")
fi

cd "$PROTOSTAR_DIR"
echo "convert: dataset=$DATASET workers=$WORKERS ${SHARD_ARGS[*]} $EXTRA_ARGS"
# shellcheck disable=SC2086  # EXTRA_ARGS is intentionally word-split
python pipelines/10_convert_raw.py \
    --dataset "$DATASET" \
    --workers "$WORKERS" \
    "${SHARD_ARGS[@]}" \
    $EXTRA_ARGS
