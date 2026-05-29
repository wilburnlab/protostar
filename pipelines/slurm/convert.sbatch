#!/bin/bash
# Stage 10 — convert a dataset's .raw into Constellation parquet bundles on OSC.
#
# One extended Cardinal allocation per dataset, high worker count. Resume-safe:
# a wall-time kill just re-runs and skips finished bundles, so re-submitting is
# the recovery path. For very large datasets, the array-job escape hatch fans
# the work out — submit with `--array=0-(K-1)` and set N_SHARDS=K below.
#
# Submit (centroid, whole dataset):
#   sbatch --export=ALL,DATASET=Zolg2017 pipelines/slurm/convert.sbatch
# Submit (profile pass):
#   sbatch --export=ALL,DATASET=Zolg2017,MODE_FLAG=--profile pipelines/slurm/convert.sbatch
# Submit (sharded array of 20):
#   sbatch --array=0-19 --export=ALL,DATASET=Wilhelm2021,N_SHARDS=20 pipelines/slurm/convert.sbatch
#
#SBATCH --job-name=protostar-convert
#SBATCH --account=<allocation>
#SBATCH --cluster=cardinal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=48:00:00
#SBATCH --output=%x-%A_%a-%j.log

set -eo pipefail  # NOT -euo: conda activation trips on unbound vars

module load miniconda3/24.1.2-py310
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate constellation

# --- knobs (override via --export) ---------------------------------------
PROTOSTAR_DIR="${PROTOSTAR_DIR:-$HOME/protostar}"   # repo on OSC home (rsynced)
DATASET="${DATASET:-all}"
MODE_FLAG="${MODE_FLAG:-}"                           # set to --profile for profile pass
WORKERS="${WORKERS:-${SLURM_CPUS_PER_TASK:-8}}"
# Optional: explicit Thermo CommonCore DLL dir (else auto-discovered in the env).
# export CONSTELLATION_THERMO_HOME=/path/to/thermo/dlls

# Array-job sharding: present only when submitted with --array.
SHARD_ARGS=()
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && -n "${N_SHARDS:-}" ]]; then
    SHARD_ARGS=(--shard "${SLURM_ARRAY_TASK_ID}" --n-shards "${N_SHARDS}")
fi

cd "$PROTOSTAR_DIR"
echo "convert: dataset=$DATASET mode='${MODE_FLAG:-centroid}' workers=$WORKERS ${SHARD_ARGS[*]}"
python pipelines/10_convert_raw.py \
    --dataset "$DATASET" \
    --workers "$WORKERS" \
    $MODE_FLAG \
    "${SHARD_ARGS[@]}"
