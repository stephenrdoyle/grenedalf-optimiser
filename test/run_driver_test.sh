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
SEX_CHROMS=sim_chrX
MITO_CHROMS=mito
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
# grenedalf's mask rule: per-sample runs must not fail, and no grid run may fail
[[ ! -s "$D/driver/grid/failed_runs.tsv" ]] || { echo "FAIL: grid runs failed: $(cat "$D/driver/grid/failed_runs.tsv")" >&2; exit 1; }
grep -q "^per_sample_masks" "$D/driver/optimiser.assess_runs.tsv" || { echo "FAIL: per_sample_masks run missing" >&2; exit 1; }
# Per-statistic judging: the check run is never recommended; a combination run is built and judged
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i} $1=="no_accessibility_mask" && $h["recommended"]=="True" {bad=1} END {exit bad}' \
    "$D/driver/optimiser.assess_runs.tsv" || { echo "FAIL: check run recommended" >&2; exit 1; }
# (the stand-in keeps more SNPs with per-sample masks, fitted bounds and min count 2, as on real data)
grep -q "^combined" "$D/driver/optimiser.assess_runs.tsv" || { echo "FAIL: no combination run" >&2; exit 1; }
awk -F'\t' '$1=="combined" && !($2==2 && $5=="per-sample" && $7=="fitted") {bad=1} END {exit bad}' "$D/driver/grid/runs.tsv" \
    || { echo "FAIL: combination is not min count 2 + per-sample + fitted: $(grep ^combined "$D/driver/grid/runs.tsv")" >&2; exit 1; }
grep -q "^combined" "$D/driver/grid/runs.tsv" || { echo "FAIL: combination run not recorded" >&2; exit 1; }
grep -q "data_metric" "$D/driver/optimiser.assess_runs.tsv" && grep -q "pair_snps" "$D/driver/optimiser.assess_runs.tsv" \
    || { echo "FAIL: data kept not per pair" >&2; exit 1; }
# No Tajima's D anywhere in the grid
! grep -q "tajimas_d" "$D/driver/optimiser.assess_stats.tsv" || { echo "FAIL: Tajima's D assessed" >&2; exit 1; }

# Cache: a second run with nothing changed reuses every slow result
"$HERE/../bin/grenedalf_optimiser.sh" "$D/config.sh" all > "$D/rerun.log" 2>&1
n_reused=$(grep -c "reusing" "$D/rerun.log" || true)
n_slow=$(grep -c -E "\] (mosdepth|profile|grenedalf|GenMap)" "$D/rerun.log" || true)
[[ $n_reused -gt 0 && $n_reused -eq $n_slow ]] || { echo "FAIL: rerun recomputed slow steps ($n_reused reused of $n_slow)" >&2; exit 1; }
# A setting that changes the depth masks: grenedalf runs recompute, mosdepth and GenMap are reused
echo "MAX_FOLD=2.5" >> "$D/config.sh"
"$HERE/../bin/grenedalf_optimiser.sh" "$D/config.sh" all > "$D/rerun2.log" 2>&1
grep -q "grenedalf recommended: reusing" "$D/rerun2.log" && { echo "FAIL: changed mask not recomputed" >&2; exit 1; }
grep -q "mosdepth poolA_L3 MAPQ>=0: reusing" "$D/rerun2.log" || { echo "FAIL: mosdepth recomputed for a mask change" >&2; exit 1; }
# Outputs from before the cache (no signatures): BAM-only steps adopted, the rest recomputed
find "$D/driver" -name "*signature" -delete
"$HERE/../bin/grenedalf_optimiser.sh" "$D/config.sh" all > "$D/rerun3.log" 2>&1
grep -q "mosdepth poolA_L3 MAPQ>=0: reusing" "$D/rerun3.log" || { echo "FAIL: pre-cache mosdepth bins not adopted" >&2; exit 1; }
grep -q "profile poolA_L3: reusing" "$D/rerun3.log" || { echo "FAIL: pre-cache profile not adopted" >&2; exit 1; }
grep -q "grenedalf recommended: reusing" "$D/rerun3.log" && { echo "FAIL: pre-cache grenedalf run adopted" >&2; exit 1; }

echo "driver test passed ($D/driver)"
