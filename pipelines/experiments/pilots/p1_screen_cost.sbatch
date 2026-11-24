#!/bin/bash
# Pilot P1 (extraction design): global MS1 sweep cost on Zolg2017 pool 1's four injections.
# Submit from the log directory (logs land in the submit dir):
#   sbatch --account=<allocation> \
#     --export=ALL,PROTOSTAR_DIR=<checkout>,DATA_ROOT=<data_root>,OUT=<data_root>/pilots/extraction/p1 \
#     <checkout>/pipelines/experiments/pilots/p1_screen_cost.sbatch
#SBATCH --job-name=protostar-p1-screen-cost
#SBATCH --cluster=cardinal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=8:00:00
#SBATCH --output=%x-%j.log
set -eo pipefail  # NOT -euo: conda activation trips on unbound vars
module load miniconda3/24.1.2-py310
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate constellation
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
: "${PROTOSTAR_DIR:?}" "${DATA_ROOT:?}" "${OUT:?}"
cd "$PROTOSTAR_DIR"
echo "host=$(hostname) protostar=$(git rev-parse --short HEAD) start=$(date -Is) tmpdir=${TMPDIR:-/tmp}"
python pipelines/experiments/pilots/p1_screen_cost.py --data-root "$DATA_ROOT" --out "$OUT" \
    --workers 4 8 16 ${EXTRA_ARGS:-}
echo "end=$(date -Is)"
