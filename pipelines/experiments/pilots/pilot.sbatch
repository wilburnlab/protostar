#!/bin/bash
# Run one extraction pilot (and optionally a second step) from pipelines/experiments/pilots/.
# Submit from the log directory (logs land in the submit dir), e.g. P2's sweep then its analysis:
#   sbatch --account=<allocation> --job-name=protostar-p2 --cpus-per-task=16 --mem=96G --time=6:00:00 \
#     --export=ALL,PROTOSTAR_DIR=<checkout>,STEP1="p2_screen_rows.py --data-root ... --out ...",STEP2="p2_analyze.py ..." \
#     <checkout>/pipelines/experiments/pilots/pilot.sbatch
#SBATCH --job-name=protostar-pilot
#SBATCH --cluster=cardinal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=%x-%j.log
set -eo pipefail  # NOT -euo: conda activation trips on unbound vars
module load miniconda3/24.1.2-py310
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate constellation
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
: "${PROTOSTAR_DIR:?}" "${STEP1:?}"
cd "$PROTOSTAR_DIR"
echo "host=$(hostname) protostar=$(git rev-parse --short HEAD) start=$(date -Is) tmpdir=${TMPDIR:-/tmp}"
# shellcheck disable=SC2086  # the steps are word-split on purpose
python pipelines/experiments/pilots/$STEP1
if [ -n "${STEP2:-}" ]; then
  echo "step 2 start=$(date -Is)"
  # shellcheck disable=SC2086
  python pipelines/experiments/pilots/$STEP2
fi
echo "end=$(date -Is)"
