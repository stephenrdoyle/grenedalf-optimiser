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
    --reference clean --min-snps 20 --out "$D/optimiser"

v=()
for p in "${POOLS[@]}"; do v+=(--depth "$p" "$D/$p.mq30.regions.bed.gz" --dist "$p" "$D/$p.mosdepth.region.dist.txt"); done
"${OPTIMISER[@]}" variance "${v[@]}" --mask-summary "$D/optimiser.mask_summary.tsv" --chrom-class "$D/optimiser.chrom_class.tsv" \
    --thresholds "$D/optimiser.depth_thresholds.tsv" --mask-bed "$D/optimiser.mask.bed" --out "$D/optimiser" > /dev/null
"${OPTIMISER[@]}" depthmask "${q[@]}" --thresholds "$D/optimiser.depth_thresholds.tsv" --chrom-class "$D/optimiser.chrom_class.tsv" \
    --fai "$D/ref.fa.fai" --base-mask "$D/optimiser.mask.bed" --mode shared --bound-method fitted --out "$D/fitted" > /dev/null

"${OPTIMISER[@]}" commands --prefix "$D/optimiser" --out-dir "$D" --fai "$D/ref.fa.fai" \
    --bam "${POOLS[@]/%/.bam}" --x-chroms sim_chrX --pool-size-arg 400 \
    --settings diversity "$D/persample.total_mask.bed" "$D/persample.sample_masks.tsv" 2 10 282 \
    --settings fst "$D/optimiser.total_mask.bed" - 3 10 282 \
    --window 6000 --stride 3000 --min-snps 50
# Smoothed depth masks: whole 100 bp blocks only
"${OPTIMISER[@]}" depthmask "${q[@]}" --thresholds "$D/optimiser.depth_thresholds.tsv" --chrom-class "$D/optimiser.chrom_class.tsv" \
    --fai "$D/ref.fa.fai" --mode per-sample --smooth 100 --out "$D/smooth" > /dev/null

# Relative lower bounds: masked where a pool is low against the others, not where all are low
rb=()
for p in "${POOLS[@]}"; do rb+=(--bins "$p" "$D/$p.mq30.regions.bed.gz"); done
"${OPTIMISER[@]}" depthmask "${q[@]}" "${rb[@]}" --thresholds "$D/optimiser.depth_thresholds.tsv" \
    --chrom-class "$D/optimiser.chrom_class.tsv" --fai "$D/ref.fa.fai" --base-mask "$D/optimiser.mask.bed" \
    --mode per-sample --bound-method relative --out "$D/relative" > /dev/null
# Callable fraction from a run's own masks, per pool and per pair (both pools' masks)
python3 - "$HERE/../bin" "$D" <<'PYTEST' || fail "per-pool callable fraction from masks"
import sys, pandas as pd
sys.path.insert(0, sys.argv[1])
from grenedalf_optimiser import mask_callable
d = sys.argv[2]
open(f"{d}/mc_total.bed", "w").write("c\t0\t10\n")
open(f"{d}/mc_a.bed", "w").write("c\t50\t100\n")
open(f"{d}/mc_b.bed", "w").write("c\t100\t150\n")
open(f"{d}/mc_list.tsv", "w").write(f"a\t{d}/mc_a.bed\nb\t{d}/mc_b.bed\n")
w = pd.DataFrame({"chrom": ["c", "c"], "start": [1, 101], "end": [100, 200]})
cal = mask_callable(w, f"{d}/mc_total.bed", f"{d}/mc_list.tsv", ["a", "b"], {"c": 200})
want = {"a": [0.4, 1.0], "b": [0.9, 0.5], "a:b": [0.4, 0.5]}
for k, v in want.items():
    assert [round(x, 3) for x in cal[k]] == v, (k, list(cal[k]), v)
PYTEST

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
for t in DIVERSITY FST; do
    grep -q "RECOMMENDED_RUN_$t=clean" "$D/optimiser.assess.params" || fail "assess should prefer the clean run for $t"
done
for t in diversity fst; do
    awk -F'\t' -v t="$t" 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i} NR>1 && $1=="artefact" {exit !($h["max_increase_" t] > 0.1)}' \
        "$D/optimiser.assess_runs.tsv" || fail "artefact run not scored worse than the clean run for $t"
done
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i} NR>1 {exit !($h["data_metric_diversity"]=="pool_snps" && $h["data_metric_fst"]=="pair_snps")}' \
    "$D/optimiser.assess_runs.tsv" || fail "data kept not per pool for diversity and per pair for fst"
# Different settings for diversity and fst reach their own commands
grep -A20 "grenedalf diversity" "$D/grenedalf_commands.sh" | grep -m1 -- "--filter-sample-min-count" | grep -q "count 2" \
    || fail "diversity command without its own min count"
grep -A20 "grenedalf fst" "$D/grenedalf_commands.sh" | grep -m1 -- "--filter-sample-min-count" | grep -q "count 3" \
    || fail "fst command without its own min count"
grep -A20 "grenedalf diversity" "$D/grenedalf_commands.sh" | grep -m1 -- "--filter-mask-samples-bed-list" | grep -q "persample" \
    || fail "diversity command without its per-pool masks"
awk '/^grenedalf fst/ {n++} n == 1' "$D/grenedalf_commands.sh" | awk '/file-prefix/ {print; exit} 1' \
    | grep -q -- "--filter-mask-samples-bed-list" && fail "autosome fst command given the diversity per-pool masks"
# (apart from the one-base pad at each chromosome end that every mask carries)
awk -F'\t' '($2 % 100) && $3 - $2 > 1 {bad=1} END {exit bad}' "$D/smooth.poolA_L3.depth_mask.bed" \
    || fail "smoothed mask not in 100 bp blocks"
[[ $(wc -l < "$D/smooth.poolA_L3.depth_mask.bed") -lt $(wc -l < "$D/persample.poolA_L3.depth_mask.bed") ]] \
    || fail "smoothing did not simplify the depth mask"
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
# No Tajima's D: every diversity command has --no-tajima-d, and there is no separate Tajima's D run
[[ $(grep -c "grenedalf diversity" "$D/grenedalf_commands.sh") -eq $(grep -c -- "--no-tajima-d" "$D/grenedalf_commands.sh") ]] \
    || fail "a diversity command computes Tajima's D"
! grep -q "tajimas_d" "$D/grenedalf_commands.sh" || fail "separate Tajima's D run in the commands"

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
# Relative bounds: poolB's deletion (chr2:150-155 kb) masked in poolB only
grep -q "BOUND_METHOD=relative" "$D/relative.depthmask.params" || fail "relative bound method not recorded"
awk -F'\t' '$1=="sim_chr2" && $2<=151000 && $3>=154000' "$D/relative.poolB_L3.depth_mask.bed" | grep -q . \
    || fail "relative bounds: poolB deletion not masked in poolB"
awk -F'\t' '$1=="sim_chr2" && $2<154000 && $3>151000 {s=($2>151000?$2:151000); e=($3<154000?$3:154000); n+=e-s}
    END {exit !(n < 1500)}' "$D/relative.poolA_L3.depth_mask.bed" || fail "relative bounds: poolB deletion masked in poolA"
! cmp -s "$D/fitted.total_mask.bed" "$D/optimiser.total_mask.bed" || fail "fitted bounds gave the same mask as fold bounds"

# Nothing named: every contig autosomal, one nuclear run, no X or mito run
"${OPTIMISER[@]}" chromclass "${c[@]}" --fai "$D/ref.fa.fai" --mask-bed "$D/optimiser.mask.bed" --pool-sizes 400 \
    --out "$D/unnamed" > /dev/null
[[ $(awk -F'\t' 'NR>1 && $2!="A"' "$D/unnamed.chrom_classes.tsv" | wc -l) -eq 0 ]] || fail "unnamed contigs not all autosomal"
mkdir -p "$D/unnamed_cmds"
"${OPTIMISER[@]}" commands --prefix "$D/unnamed" --out-dir "$D/unnamed_cmds" --fai "$D/ref.fa.fai" \
    --bam "${POOLS[@]/%/.bam}" --pool-size-arg 400 \
    --settings diversity "$D/optimiser.total_mask.bed" - 2 10 282 --settings fst "$D/optimiser.total_mask.bed" - 2 10 282 \
    --window 6000 --stride 3000 > /dev/null
grep -q -- "--out-dir \$OUTDIR/nuclear " "$D/unnamed_cmds/grenedalf_commands.sh" || fail "no single nuclear run when nothing is named"
! grep -q -- "--out-dir \$OUTDIR/\(X\|mito\) " "$D/unnamed_cmds/grenedalf_commands.sh" || fail "X or mito run without named contigs"

echo "all checks passed ($D)"
