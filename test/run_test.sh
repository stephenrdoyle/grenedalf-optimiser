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
    --kmer 150 --exclude-chroms mito "${depth[@]}" --out "$D/optimiser"

"${SIM[@]}" dist "$D"
a=()
for p in "${POOLS[@]}"; do a+=(--dist "$p" "$D/$p.mosdepth.region.dist.txt" --stats "$p" "$D/$p.stats"); done
"${OPTIMISER[@]}" thresholds "${a[@]}" --x-chroms sim_chrX --mask-summary "$D/optimiser.mask_summary.tsv" --pool-sizes 400 --out "$D/optimiser"

c=()
for p in "${POOLS[@]}"; do c+=(--depth "$p" "$D/$p.mq30.regions.bed.gz" --depth-mq0 "$p" "$D/$p.mq0.regions.bed.gz"); done
"${OPTIMISER[@]}" chromclass "${c[@]}" --fai "$D/ref.fa.fai" --mask-bed "$D/optimiser.mask.bed" --x-chroms sim_chrX \
    --mito-chroms mito --pool-sizes 400 --out "$D/optimiser"

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

v=()
for p in "${POOLS[@]}"; do v+=(--depth "$p" "$D/$p.mq30.regions.bed.gz" --dist "$p" "$D/$p.mosdepth.region.dist.txt"); done
"${OPTIMISER[@]}" variance "${v[@]}" --mask-summary "$D/optimiser.mask_summary.tsv" --chrom-class "$D/optimiser.chrom_class.tsv" \
    --thresholds "$D/optimiser.depth_thresholds.tsv" --mask-bed "$D/optimiser.mask.bed" --out "$D/optimiser" > /dev/null
"${OPTIMISER[@]}" depthmask "${q[@]}" --thresholds "$D/optimiser.depth_thresholds.tsv" --chrom-class "$D/optimiser.chrom_class.tsv" \
    --fai "$D/ref.fa.fai" --base-mask "$D/optimiser.mask.bed" --mode shared --bound-method fitted --out "$D/fitted" > /dev/null

"${OPTIMISER[@]}" commands --prefix "$D/optimiser" --out-dir "$D" --fai "$D/ref.fa.fai" \
    --bam "${POOLS[@]/%/.bam}" --x-chroms sim_chrX --total-mask "$D/optimiser.total_mask.bed" --pool-size-arg 400 \
    --min-count 3 --min-depth 10 --max-depth 282 --window 6000 --stride 3000 --min-snps 50

# --- checks
reasons=$(cut -f4 "$D/optimiser.mask_reasons.bed" | tr ',' '\n' | sort -u | tr '\n' ' ')
for want in excluded gap softmask low_mappability collapsed ambiguous; do
    [[ " $reasons " == *" $want "* ]] || fail "mask reason $want not found (got: $reasons)"
done
grep -q $'^sim_chr1\t100000\t110000\tcollapsed$' "$D/optimiser.mask_reasons.bed" || fail "collapsed repeat not at chr1:100-110 kb"
awk -F'\t' '$1=="sim_chr2" && $2>=150000 && $3<=155000 && $4=="poolB_L3"' "$D/optimiser.pool_specific_depth.tsv" | grep -q . \
    || fail "pool-specific deletion in poolB not reported"
grep -q $'^sim_chr2\t150000\t155000' "$D/optimiser.mask.bed" && fail "pool-specific deletion was masked for all pools"
awk -F'\t' 'NR>1 && $1=="poolA_L3" && $2=="X" && $4 >= 26 && $4 <= 34 {ok=1} END {exit !ok}' "$D/optimiser.depth_thresholds.tsv" \
    || fail "poolA X median not ~30 (0.75 x 40)"
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

# X simulated at 0.6-0.8x autosome depth (mixed-sex pools): separate run, X pool size = 400 x ratio.
# +-15%: the simulated X is only 150 kb (150 bins, ~30 independent depth blocks); a real X has tens of thousands
# (300 for 50:50); mito always separate
source "$D/optimiser.chromclass.params"
[[ $X_SEPARATE == 1 && $MITO_SEPARATE == 1 ]] || fail "X / mito not flagged for separate runs"
awk -F'\t' 'BEGIN {r["poolA_L3"]=.75; r["poolB_L3"]=.75; r["poolC_L3"]=.7; r["poolD_L3"]=.8; r["poolE_L3"]=.75; r["poolF_L3"]=.6}
    {e = 400 * r[$1]; if ($2 < 0.85 * e || $2 > 1.15 * e) bad = 1} END {exit bad}' "$D/optimiser.x_pool_sizes.tsv" \
    || fail "X pool sizes not 400 x simulated X ratio: $(tr '\n' ' ' < "$D/optimiser.x_pool_sizes.tsv")"
# Known 50:50 design: X pool size 300 everywhere; poolF (simulated 0.6, i.e. 80% males) flagged as inconsistent
"${OPTIMISER[@]}" chromclass "${c[@]}" --fai "$D/ref.fa.fai" --mask-bed "$D/optimiser.mask.bed" --x-chroms sim_chrX \
    --mito-chroms mito --pool-sizes 400 --male-fraction 0.5 --out "$D/known" > /dev/null
awk -F'\t' '$2 != 300 {bad=1} END {exit bad}' "$D/known.x_pool_sizes.tsv" || fail "known 50:50 did not give X pool size 300"
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i} $2=="X" {m[$1]=$h["matches_known"]}
    END {exit !(m["poolA_L3"]=="True" && m["poolF_L3"]=="False")}' "$D/known.chrom_classes.tsv" \
    || fail "known sex ratio consistency check"
for run in autosomes X mito; do
    grep -q -- "--out-dir \$OUTDIR/$run " "$D/grenedalf_commands.sh" || fail "no $run run in grenedalf_commands.sh"
done
grep -q -- "--pool-sizes $D/optimiser.x_pool_sizes.tsv" "$D/grenedalf_commands.sh" || fail "X run without X pool sizes"
bash -n "$D/grenedalf_commands.sh" || fail "grenedalf_commands.sh is not valid bash"
# min count 3: theta pi / Watterson without Tajima's D, and a separate Tajima's D run at count 2
grep -q -- "--no-tajima-d" "$D/grenedalf_commands.sh" || fail "diversity at min count 3 without --no-tajima-d"
awk '/file-prefix autosomes_tajimas_d_/ {t=1} END {exit !t}' "$D/grenedalf_commands.sh" || fail "no separate Tajima's D run"
awk 'BEGIN {RS=""} /tajimas_d_/ && !/--filter-sample-min-count 2 / {bad=1} END {exit bad}' "$D/grenedalf_commands.sh" \
    || fail "Tajima's D run not at min count 2"

# Impossible X ratio (X depth above autosomes, as on reference-divergent real data): flagged,
# and 0.5 assumed rather than "no males"
mkdir -p "$D/hix"
for p in "${POOLS[@]}"; do
    gzip -dc "$D/$p.mq30.regions.bed.gz" | awk -F'\t' -v OFS='\t' '$1=="sim_chrX" {$4 = $4 * 1.6} 1' | gzip > "$D/hix/$p.mq30.regions.bed.gz"
done
h=()
for p in "${POOLS[@]}"; do h+=(--depth "$p" "$D/hix/$p.mq30.regions.bed.gz"); done
"${OPTIMISER[@]}" chromclass "${h[@]}" --fai "$D/ref.fa.fai" --mask-bed "$D/optimiser.mask.bed" --x-chroms sim_chrX \
    --mito-chroms mito --pool-sizes 400 --out "$D/hix/c" > /dev/null 2> "$D/hix/err"
source "$D/hix/c.chromclass.params"
[[ $X_RATIO_UNRELIABLE == 1 && $X_SEPARATE == 1 ]] || fail "impossible X ratio not flagged"
grep -q WARNING "$D/hix/err" || fail "no warning for an impossible X ratio"
awk -F'\t' '$2 != 300 {bad=1} END {exit bad}' "$D/hix/c.x_pool_sizes.tsv" || fail "impossible X ratio: 0.5 not assumed"
# Mito: per-pool bounds and per-pool masks
[[ $(wc -l < "$D/optimiser.mito_sample_masks.tsv") -eq 6 ]] || fail "no per-pool mito masks"
awk -F'\t' 'NR>1 {l[$6]=1} END {exit length(l) < 3}' "$D/optimiser.mito.tsv" || fail "mito bounds not per pool"
grep -q -- "--filter-mask-samples-bed-list $D/optimiser.mito_sample_masks.tsv" "$D/grenedalf_commands.sh" \
    || fail "mito run without per-pool masks"

# Coverage variance: simulated depth varies in 5 kb blocks shared by all pools (GC-like), so most
# variance persists at 1 kb, pools correlate, and poolB's deletion is a pool-specific low block
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i} NR>1 && $h["class"]=="A" && $h["share_1kb"] < 0.3 {bad=1} END {exit bad}' \
    "$D/optimiser.variance.tsv" || fail "variance: regional share at 1 kb unexpectedly low"
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i} NR>1 && $h["mean_corr_with_others"] < 0.5 {bad=1} END {exit bad}' \
    "$D/optimiser.variance.tsv" || fail "variance: pools' depth not correlated"
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i} NR>1 && $1=="poolB_L3" && $h["pool_specific"] > 0 {ok=1} END {exit !ok}' \
    "$D/optimiser.variance.tsv" || fail "variance: poolB deletion not pool-specific"
# Fitted bounds: within each pool's distribution, and a different mask from the fold bounds
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i} NR>1 && !($h["fit_lower"] < $h["median"] && $h["median"] < $h["fit_upper"]) {bad=1} END {exit bad}' \
    "$D/optimiser.depth_thresholds.tsv" || fail "fitted bounds do not bracket the median"
grep -q "BOUND_METHOD=fitted" "$D/fitted.depthmask.params" || fail "fitted bound method not recorded"
! cmp -s "$D/fitted.total_mask.bed" "$D/optimiser.total_mask.bed" || fail "fitted bounds gave the same mask as fold bounds"

# Nothing named: every contig autosomal, one nuclear run, no X or mito run
"${OPTIMISER[@]}" chromclass "${c[@]}" --fai "$D/ref.fa.fai" --mask-bed "$D/optimiser.mask.bed" --pool-sizes 400 \
    --out "$D/unnamed" > /dev/null
[[ $(awk -F'\t' 'NR>1 && $2!="A"' "$D/unnamed.chrom_classes.tsv" | wc -l) -eq 0 ]] || fail "unnamed contigs not all autosomal"
mkdir -p "$D/unnamed_cmds"
"${OPTIMISER[@]}" commands --prefix "$D/unnamed" --out-dir "$D/unnamed_cmds" --fai "$D/ref.fa.fai" \
    --bam "${POOLS[@]/%/.bam}" --total-mask "$D/optimiser.total_mask.bed" --pool-size-arg 400 \
    --min-count 2 --min-depth 10 --max-depth 282 --window 6000 --stride 3000 > /dev/null
grep -q -- "--out-dir \$OUTDIR/nuclear " "$D/unnamed_cmds/grenedalf_commands.sh" || fail "no single nuclear run when nothing is named"
! grep -q -- "--out-dir \$OUTDIR/\(X\|mito\) " "$D/unnamed_cmds/grenedalf_commands.sh" || fail "X or mito run without named contigs"

echo "all checks passed ($D)"
