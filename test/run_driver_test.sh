#!/usr/bin/env bash
# Runs bin/grenedalf_optimiser.sh end to end with stand-ins for the external tools
# (test/stub_tools), on the synthetic data from run_test.sh. Tests the driver's
# plumbing, not the tools.
#   test/run_driver_test.sh [OUTDIR]

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
D=${1:-$(mktemp -d)}
"$HERE/run_test.sh" "$D" > "$D.run_test.log" 2>&1 || { tail "$D.run_test.log"; exit 1; }

export SIM_DIR=$D
export PATH="$HERE/stub_tools:$PATH"
cat > "$D/config.sh" <<CONF
REF=$D/ref.fa
FAI=$D/ref.fa.fai
BAMS=($D/poolA_L3.bam $D/poolB_L3.bam $D/poolC_L3.bam $D/poolD_L3.bam $D/poolE_L3.bam $D/poolF_L3.bam)
OUT=$D/driver
POOL_SIZE=400
CHROM_REGEX=sim_
X_REGEX=chrX
MITO_REGEX=mito
TARGET_SNPS=20
LD_MAX_DIST=50000
CONF
for p in poolA_L3 poolB_L3 poolC_L3 poolD_L3 poolE_L3 poolF_L3; do touch "$D/$p.bam.bai"; echo bam > "$D/$p.bam"; done
"$HERE/../bin/grenedalf_optimiser.sh" "$D/config.sh" check
"$HERE/../bin/grenedalf_optimiser.sh" "$D/config.sh" all
for f in optimiser.mask.bed optimiser.depth_thresholds.tsv optimiser.total_mask.bed optimiser.ld.params optimiser.windows.params \
         optimiser.assess_runs.tsv optimiser.chrom_classes.tsv params.sh grenedalf_commands.sh optimiser_report.html grid/per_sample_masks/mask.sample_masks.tsv; do
    [[ -s "$D/driver/$f" ]] || { echo "FAIL: missing $f" >&2; exit 1; }
done
echo "driver test passed ($D/driver)"
