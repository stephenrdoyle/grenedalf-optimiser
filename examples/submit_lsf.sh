#!/usr/bin/env bash
# Submit the whole grenedalf-optimiser run as one LSF job (Sanger farm).
#   examples/submit_lsf.sh CONFIG
# Run the check step interactively first; it takes a minute or two.
# The grid (about ten genome-wide grenedalf runs) is the slow part, hence the long queue.
# Each step skips nothing it needs, so a failed run can be resumed step by step:
#   bin/grenedalf_optimiser.sh CONFIG <step>

set -euo pipefail
CONFIG=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mkdir -p logs

bsub -q long -n 8 -R "span[hosts=1] select[mem>32000] rusage[mem=32000]" -M 32000 \
    -J grenedalf_optimiser -o logs/grenedalf_optimiser.%J.o -e logs/grenedalf_optimiser.%J.e \
    "module load grenedalf/0.6.3 2>/dev/null || true
     module load mosdepth samtools genmap 2>/dev/null || true
     $HERE/bin/grenedalf_optimiser.sh $CONFIG all"
