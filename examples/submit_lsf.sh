#!/usr/bin/env bash
# Submit the whole grenedalf-optimiser run as one LSF job (Sanger farm).
#   examples/submit_lsf.sh CONFIG
#
# The tool activates its conda environment itself (CONDA_ENV in the config; create it
# once with: bin/grenedalf_optimiser.sh CONFIG setup). Run the check step first:
#   bin/grenedalf_optimiser.sh CONFIG check
# The grid (about ten genome-wide grenedalf runs) is the slow part, hence the long queue.
# A failed run can be resumed step by step: bin/grenedalf_optimiser.sh CONFIG <step>

set -euo pipefail
CONFIG=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mkdir -p logs

bsub -q long -n 8 -R "span[hosts=1] select[mem>32000] rusage[mem=32000]" -M 32000 \
    -J grenedalf_optimiser -o logs/grenedalf_optimiser.%J.o -e logs/grenedalf_optimiser.%J.e \
    "$HERE/bin/grenedalf_optimiser.sh" "$CONFIG" all
