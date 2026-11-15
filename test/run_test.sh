#!/usr/bin/env bash
# Runs every grenedalf_optimiser.py step on synthetic data (test/simulate.py) and checks the
# built-in features are found. Needs only python3 with numpy and pandas.
#   test/run_test.sh [OUTDIR]

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OPTIMISER=(python3 "$HERE/../bin/grenedalf_optimiser.py")
SIM=(python3 "$HERE/simulate.py")
D=${1:-$(mktemp -d)}
mkdir -p "$D"
POOLS=(poolA_L3 poolB_L3 poolC_L3 poolD_L3 poolE_L3 poolF_L3)

fail() { echo "FAIL: $*" >&2; exit 1; }

"${SIM[@]}" base "$D"

depth=()
for p in "${POOLS[@]}"; do depth+=(--depth "$p" "$D/$p.mq0.regions.bed.gz" "$D/$p.mq30.regions.bed.gz"); done
"${OPTIMISER[@]}" mask --fai "$D/ref.fa.fai" --ref "$D/ref.fa" --softmask --genmap "$D/genmap.bedgraph" \
    --kmer 150 --chrom-regex sim_ "${depth[@]}" --out "$D/optimiser"

"${SIM[@]}" dist "$D"
a=()
for p in "${POOLS[@]}"; do a+=(--dist "$p" "$D/$p.mosdepth.region.dist.txt" --stats "$p" "$D/$p.stats"); done
"${OPTIMISER[@]}" thresholds "${a[@]}" --mask-summary "$D/optimiser.mask_summary.tsv" --out "$D/optimiser"

"${SIM[@]}" quantize "$D"
q=()
for p in "${POOLS[@]}"; do q+=(--quantized "$p" "$D/$p.quantized.bed.gz"); done
"${OPTIMISER[@]}" depthmask "${q[@]}" --thresholds "$D/optimiser.depth_thresholds.tsv" --chrom-class "$D/optimiser.chrom_class.tsv" \
    --fai "$D/ref.fa.fai" --base-mask "$D/optimiser.mask.bed" --mode per-sample --out "$D/persample" > /dev/null
"${OPTIMISER[@]}" depthmask "${q[@]}" --thresholds "$D/optimiser.depth_thresholds.tsv" --chrom-class "$D/optimiser.chrom_class.tsv" \
    --fai "$D/ref.fa.fai" --base-mask "$D/optimiser.mask.bed" --mode shared --out "$D/optimiser"

"${SIM[@]}" grenedalf "$D"
"${OPTIMISER[@]}" ld --frequency "$D/frequency.csv" --mask-bed "$D/optimiser.total_mask.bed" --min-count 3 \
    --min-depth 10 --max-dist 50000 --out "$D/optimiser"
"${OPTIMISER[@]}" windows --table "$D/clean.diversity.csv" --mask-bed "$D/optimiser.total_mask.bed" --ld "$D/optimiser.ld.params" \
    --out "$D/optimiser"

r=()
for p in "${POOLS[@]}"; do r+=(--depth "$p" "$D/$p.mq30.regions.bed.gz"); done
"${OPTIMISER[@]}" assess --run clean "$D/clean.diversity.csv" "$D/clean.fst.csv" \
    --run artefact "$D/artefact.diversity.csv" "$D/artefact.fst.csv" "${r[@]}" \
    --strict clean --min-snps 20 --out "$D/optimiser"

# --- checks
reasons=$(cut -f4 "$D/optimiser.mask_reasons.bed" | tr ',' '\n' | sort -u | tr '\n' ' ')
for want in excluded gap softmask low_mappability collapsed ambiguous; do
    [[ " $reasons " == *" $want "* ]] || fail "mask reason $want not found (got: $reasons)"
done
grep -q $'^sim_chr1\t100000\t110000\tcollapsed$' "$D/optimiser.mask_reasons.bed" || fail "collapsed repeat not at chr1:100-110 kb"
awk -F'\t' '$1=="sim_chr2" && $2>=150000 && $3<=155000 && $4=="poolB_L3"' "$D/optimiser.pool_specific_depth.tsv" | grep -q . \
    || fail "pool-specific deletion in poolB not reported"
grep -q $'^sim_chr2\t150000\t155000' "$D/optimiser.mask.bed" && fail "pool-specific deletion was masked for all pools"
awk -F'\t' 'NR>1 && $1=="poolA_L3" && $2=="X" && $4 >= 17 && $4 <= 23 {ok=1} END {exit !ok}' "$D/optimiser.depth_thresholds.tsv" \
    || fail "poolA X median not ~20"
grep -q "OK: grenedalf masked" "$D/optimiser.windows.txt" || fail "mask orientation check"
grep -q "RECOMMENDED_RUN=clean" "$D/optimiser.assess.params" || fail "assess should prefer the clean run"
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i} NR>1 && $1=="artefact" {exit !($h["artefact_vs_strict"] > 0.2)}' \
    "$D/optimiser.assess_runs.tsv" || fail "artefact run not scored worse than the clean run"
[[ $(wc -l < "$D/persample.sample_masks.tsv") -eq 6 ]] || fail "per-sample mask list"
# Simulated LD scale is 2 kb (correlation exp(-d/2000)): r2 excess should be gone by a few kb
source "$D/optimiser.ld.params"
[[ $LD_D10 != NA && $LD_D10 -ge 1000 && $LD_D10 -le 8000 ]] || fail "LD_D10 $LD_D10 not near the simulated 2 kb scale"
source "$D/optimiser.windows.params"
[[ $WINDOW -ge $LD_D10 ]] || fail "window $WINDOW narrower than the linkage scale $LD_D10"

echo "all checks passed ($D)"
