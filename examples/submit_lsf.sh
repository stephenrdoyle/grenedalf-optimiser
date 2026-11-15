#!/usr/bin/env bash
# Submit the whole grenedalf-optimiser run as one LSF job (Sanger farm).
#   examples/submit_lsf.sh CONFIG
#
# Tools come from the conda environment in environment.yml (default name
# grenedalf-optimiser; set CONDA_ENV to use another), or from modules if that
# environment does not exist. Run the check step interactively first:
#   conda activate grenedalf-optimiser
#   bin/grenedalf_optimiser.sh CONFIG check
# The grid (about ten genome-wide grenedalf runs) is the slow part, hence the long queue.
# A failed run can be resumed step by step: bin/grenedalf_optimiser.sh CONFIG <step>

set -euo pipefail
CONFIG=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ENV_NAME=${CONDA_ENV:-grenedalf-optimiser}
mkdir -p logs

cat > logs/grenedalf_optimiser.job.sh <<JOB
#!/usr/bin/env bash
set -eo pipefail
if command -v conda > /dev/null && conda env list | awk '{print \$1}' | grep -qx "$ENV_NAME"; then
    source "\$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME"
    echo "using conda environment $ENV_NAME" >&2
else
    echo "conda environment $ENV_NAME not found: using modules" >&2
    module load grenedalf/0.6.3 mosdepth samtools genmap
fi
set -u
"$HERE/bin/grenedalf_optimiser.sh" "$CONFIG" all
JOB
chmod +x logs/grenedalf_optimiser.job.sh

bsub -q long -n 8 -R "span[hosts=1] select[mem>32000] rusage[mem=32000]" -M 32000 \
    -J grenedalf_optimiser -o logs/grenedalf_optimiser.%J.o -e logs/grenedalf_optimiser.%J.e \
    logs/grenedalf_optimiser.job.sh
