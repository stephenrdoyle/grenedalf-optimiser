#!/usr/bin/env bash
# grenedalf_optimiser.sh CONFIG STEP
#
# Chooses grenedalf filter settings from the technical properties of pool-seq BAMs,
# before any biological signal is looked at. See README.md.
#
# Steps (run in order, or "all"):
#   mask       GenMap mappability + mosdepth depth bins -> accessibility mask
#   profile    per-pool depth distribution in the accessible genome, samtools stats
#   thresholds per-pool depth bounds, minimum allele count, quantize cut points
#   chromclass do X and mitochondrion depth differ from autosomes? (separate runs, pool sizes)
#   quantize   mosdepth --quantize per pool at those cut points
#   depthmask  depth masks for the recommended bounds
#   pilot      grenedalf with the recommended settings -> linkage scale, window width, mask check
#   grid       grenedalf over a grid of settings -> technical assessment
#   params     OUT/params.sh and OUT/grenedalf_commands.sh for the full grenedalf run
#
# pilot and grid use the whole genome unless REGION is set (e.g. one chromosome for a quick trial).
#   report     OUT/optimiser_report.html: every step's results in one self-contained page
#
# Needs on PATH (or set in the config): python3 with numpy and pandas, mosdepth,
# samtools, genmap (unless MAPPABILITY_BEDGRAPH is given) and grenedalf.

set -euo pipefail

if [[ $# -ne 2 ]]; then
    sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

CONFIG=$1
STEP=$2
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OPTIMISER=(python3 "$HERE/grenedalf_optimiser.py")

# Defaults; the config overrides them
THREADS=4
KMER=150
GENMAP_ERRORS=2
MAPPABILITY_BEDGRAPH=
SOFTMASK=0
BIN=1000
CHROM_REGEX=.
X_REGEX='chrX|_X_|^X$'
MITO_REGEX='mitochondri|^chrM$|^MT$|mtDNA'
REGION=               # empty: whole genome (recommended); or a grenedalf --filter-region for a quick trial
POOL_SIZES=           # per-pool sizes, POOL=N,POOL=N (overrides POOL_SIZE)
MALE_FRACTION=        # known male fraction (0.5 for 50:50), or POOL=m,...; empty = estimate from X depth
MAPQ=30
BASEQ=30
MIN_FOLD=0.5
MAX_FOLD=2.0
ABS_MIN_DEPTH=10
GRID_MIN_FOLDS=0.33,0.67
GRID_MAX_FOLDS=1.5,3.0
TARGET_SNPS=50
LD_REGION=            # chromosomes for the linkage scale (space separated); default: all analysed
LD_MAX_DIST=100000
PILOT_WINDOW=5000
FST_METHOD=unbiased-nei
TOLERANCE=0.05
GRENEDALF=grenedalf
MOSDEPTH=mosdepth
SAMTOOLS=samtools
GENMAP=genmap
# shellcheck source=/dev/null
source "$CONFIG"

: "${REF:?set REF in the config}" "${FAI:?set FAI}" "${OUT:?set OUT}"
[[ -n "${POOL_SIZE:-}${POOL_SIZES:-}" ]] || { echo "set POOL_SIZE (or POOL_SIZES) in the config" >&2; exit 1; }
[[ ${#BAMS[@]} -gt 0 ]] || { echo "set BAMS=(...) in the config" >&2; exit 1; }
REGION=${REGION:-${PILOT_REGION:-}}   # PILOT_REGION: older configs

# Pools are named as grenedalf names samples: the BAM file name without its extension
NAMES=()
for bam in "${BAMS[@]}"; do
    n=$(basename "$bam")
    NAMES+=("${n%.*}")
done

mkdir -p "$OUT"/{mappability,depth,profile,quantize,pilot,grid}
OUT=$(cd "$OUT" && pwd)
P="$OUT/optimiser"

# Pool sizes: POOL_SPEC for the python steps; POOL_SIZE_ARG for grenedalf --pool-sizes
# (a number, or a sample<TAB>size file when the pools differ)
if [[ -n "${POOL_SIZES:-}" ]]; then
    POOL_SPEC=$POOL_SIZES
    POOL_SIZE_ARG="$OUT/pool_sizes.tsv"
    tr ',' '\n' <<<"$POOL_SIZES" | tr '=' '\t' > "$POOL_SIZE_ARG"
else
    POOL_SPEC=$POOL_SIZE
    POOL_SIZE_ARG=$POOL_SIZE
fi

log() { echo "[$(date '+%F %T')] $*" >&2; }

load_params() {
    for f in "$P".mask.params "$P".thresholds.params "$P".chromclass.params "$P".depthmask.params \
             "$P".windows.params "$P".assess.params; do
        # shellcheck source=/dev/null
        [[ -f "$f" ]] && source "$f"
    done
    return 0
}

check_grenedalf() {
    local help
    help=$("$GRENEDALF" diversity --help 2>&1)
    for opt in --filter-mask-total-bed --filter-mask-sample-bed-list ${REGION:+--filter-region} \
               --filter-sample-min-count --filter-sample-min-read-depth --filter-sample-max-read-depth \
               --allow-file-overwriting; do
        grep -q -- "$opt" <<<"$help" || { echo "grenedalf diversity has no $opt; this tool needs grenedalf >= 0.6" >&2; exit 1; }
    done
}

# run_grenedalf DIR TOTAL_MASK SAMPLE_MASK_LIST MIN_COUNT MIN_DEPTH MAX_DEPTH WIDTH STRIDE
run_grenedalf() {
    local dir=$1 mask=$2 sample_masks=$3 mc=$4 mind=$5 maxd=$6 width=$7 stride=$8
    local args=(--reference-genome-fai "$FAI" --threads "$THREADS"
                --sam-min-map-qual "$MAPQ" --sam-min-base-qual "$BASEQ"
                --filter-mask-total-bed "$mask"
                --filter-sample-min-count "$mc"
                --filter-sample-min-read-depth "$mind" --filter-sample-max-read-depth "$maxd"
                --window-type interval --window-interval-width "$width" --window-interval-stride "$stride"
                --window-average-policy valid-loci
                --pool-sizes "$POOL_SIZE_ARG" --separator-char tab --out-dir "$dir" --allow-file-overwriting)
    if [[ -n "$REGION" ]]; then args+=(--filter-region "$REGION"); fi
    if [[ -n "$sample_masks" ]]; then args+=(--filter-mask-sample-bed-list "$sample_masks"); fi
    for bam in "${BAMS[@]}"; do args+=(--sam-path "$bam"); done
    mkdir -p "$dir"
    "$GRENEDALF" diversity "${args[@]}" > "$dir/diversity.log" 2>&1
    "$GRENEDALF" fst "${args[@]}" --method "$FST_METHOD" > "$dir/fst.log" 2>&1
}

step_mask() {
    local bedgraph=$MAPPABILITY_BEDGRAPH
    if [[ -z "$bedgraph" ]]; then
        bedgraph=$(ls "$OUT"/mappability/*.genmap.bedgraph 2>/dev/null | head -1 || true)
        if [[ -z "$bedgraph" ]]; then
            log "GenMap index and map (k=$KMER, e=$GENMAP_ERRORS)"
            rm -rf "$OUT/mappability/index"
            "$GENMAP" index -F "$REF" -I "$OUT/mappability/index"
            "$GENMAP" map -K "$KMER" -E "$GENMAP_ERRORS" -I "$OUT/mappability/index" \
                -O "$OUT/mappability" -bg -T "$THREADS"
            bedgraph=$(ls "$OUT"/mappability/*.genmap.bedgraph | head -1)
        fi
    fi
    local depth_args=()
    for i in "${!BAMS[@]}"; do
        n=${NAMES[$i]}
        for mq in 0 "$MAPQ"; do
            if [[ ! -s "$OUT/depth/$n.mq$mq.regions.bed.gz" ]]; then
                log "mosdepth $n MAPQ>=$mq"
                "$MOSDEPTH" -t "$THREADS" -n --by "$BIN" --mapq "$mq" "$OUT/depth/$n.mq$mq" "${BAMS[$i]}"
            fi
        done
        depth_args+=(--depth "$n" "$OUT/depth/$n.mq0.regions.bed.gz" "$OUT/depth/$n.mq$MAPQ.regions.bed.gz")
    done
    local soft=()
    [[ "$SOFTMASK" == 1 ]] && soft=(--softmask)
    "${OPTIMISER[@]}" mask --fai "$FAI" --ref "$REF" ${soft[@]+"${soft[@]}"} --genmap "$bedgraph" --kmer "$KMER" \
        --chrom-regex "$CHROM_REGEX" "${depth_args[@]}" --out "$P"
}

step_profile() {
    load_params
    for i in "${!BAMS[@]}"; do
        n=${NAMES[$i]}
        log "profile $n"
        "$MOSDEPTH" -t "$THREADS" -n --by "$ACCESSIBLE_BED" --mapq "$MAPQ" "$OUT/profile/$n" "${BAMS[$i]}"
        "$SAMTOOLS" stats -@ "$THREADS" "${BAMS[$i]}" > "$OUT/profile/$n.stats"
    done
}

step_thresholds() {
    local a=()
    for n in "${NAMES[@]}"; do
        a+=(--dist "$n" "$OUT/profile/$n.mosdepth.region.dist.txt" --stats "$n" "$OUT/profile/$n.stats")
    done
    "${OPTIMISER[@]}" thresholds "${a[@]}" --mask-summary "$P.mask_summary.tsv" --x-regex "$X_REGEX" \
        --min-fold "$MIN_FOLD" --max-fold "$MAX_FOLD" --abs-min-depth "$ABS_MIN_DEPTH" \
        --grid-min-folds "$GRID_MIN_FOLDS" --grid-max-folds "$GRID_MAX_FOLDS" \
        --base-qual "$BASEQ" --map-qual "$MAPQ" --pool-sizes "$POOL_SPEC" --out "$P"
}

step_chromclass() {
    load_params
    local a=()
    for n in "${NAMES[@]}"; do a+=(--depth "$n" "$OUT/depth/$n.mq$MAPQ.regions.bed.gz"); done
    "${OPTIMISER[@]}" chromclass "${a[@]}" --fai "$FAI" --mask-bed "$MASK_BED" --chrom-regex "$CHROM_REGEX" \
        --x-regex "$X_REGEX" --mito-regex "$MITO_REGEX" --pool-sizes "$POOL_SPEC" \
        ${MALE_FRACTION:+--male-fraction "$MALE_FRACTION"} \
        --min-fold "$MIN_FOLD" --max-fold "$MAX_FOLD" --abs-min-depth "$ABS_MIN_DEPTH" --out "$P"
}

step_quantize() {
    for i in "${!BAMS[@]}"; do
        n=${NAMES[$i]}
        spec=$(awk -v n="$n" '$1 == n {print $2}' "$P.quantize.tsv")
        log "mosdepth --quantize $spec $n"
        "$MOSDEPTH" -t "$THREADS" -n --quantize "$spec" --mapq "$MAPQ" "$OUT/quantize/$n" "${BAMS[$i]}"
    done
}

# depthmask PREFIX MIN_FOLD MAX_FOLD MODE BASE_MASK
depthmask() {
    local a=()
    for n in "${NAMES[@]}"; do a+=(--quantized "$n" "$OUT/quantize/$n.quantized.bed.gz"); done
    "${OPTIMISER[@]}" depthmask "${a[@]}" --thresholds "$P.depth_thresholds.tsv" --chrom-class "$P.chrom_class.tsv" \
        --fai "$FAI" --base-mask "$5" --min-fold "$2" --max-fold "$3" --mode "$4" --out "$1"
}

step_depthmask() {
    load_params
    depthmask "$P" "$MIN_FOLD" "$MAX_FOLD" shared "$MASK_BED"
}

step_pilot() {
    check_grenedalf
    load_params
    run_grenedalf "$OUT/pilot" "$TOTAL_MASK_BED" "$SAMPLE_MASK_LIST" "$MIN_COUNT" "$MIN_DEPTH" "$MAX_DEPTH" \
        "$PILOT_WINDOW" "$((PILOT_WINDOW / 2))"

    # Allele counts for the linkage scale. grenedalf frequency takes no sample filters,
    # so the ld step applies the mask, minimum count and depth bounds itself.
    if [[ ! -s "$OUT/pilot/frequency.csv.gz" ]]; then
        local sam=()
        for bam in "${BAMS[@]}"; do sam+=(--sam-path "$bam"); done
        "$GRENEDALF" frequency --reference-genome-fai "$FAI" --threads "$THREADS" \
            --sam-min-map-qual "$MAPQ" --sam-min-base-qual "$BASEQ" "${sam[@]}" \
            --write-sample-counts --separator-char tab --compress \
            --out-dir "$OUT/pilot" --allow-file-overwriting > "$OUT/pilot/frequency.log" 2>&1
    fi
    local ld_region=${LD_REGION:-${REGION%%:*}}
    # shellcheck disable=SC2086
    "${OPTIMISER[@]}" ld --frequency "$OUT/pilot/frequency.csv.gz" ${ld_region:+--region $ld_region} \
        --mask-bed "$TOTAL_MASK_BED" --min-count "$MIN_COUNT" --min-depth "$MIN_DEPTH" --max-depth "$MAX_DEPTH" \
        --max-dist "$LD_MAX_DIST" --out "$P"
    "${OPTIMISER[@]}" windows --table "$OUT/pilot/diversity.csv" --target-snps "$TARGET_SNPS" \
        --mask-bed "$TOTAL_MASK_BED" --ld "$P.ld.params" --out "$P"
}

step_grid() {
    check_grenedalf
    load_params
    local g="$OUT/grid" runs=() mc=$MIN_COUNT
    # Gaps and excluded chromosomes only: shows what the accessibility mask is doing
    awk '$4 ~ /excluded|gap/' "$P.mask_reasons.bed" | cut -f1-3 > "$g/minimal_mask.bed"

    # One setting changed at a time from the recommended one
    # Each run: name, min count, lower fold, upper fold, mask mode, accessibility mask
    runs+=("recommended $mc $MIN_FOLD $MAX_FOLD shared $MASK_BED")
    if [[ $mc -gt 1 ]]; then runs+=("min_count_$((mc - 1)) $((mc - 1)) $MIN_FOLD $MAX_FOLD shared $MASK_BED"); fi
    runs+=("min_count_$((mc + 1)) $((mc + 1)) $MIN_FOLD $MAX_FOLD shared $MASK_BED")
    for f in ${GRID_MIN_FOLDS//,/ }; do runs+=("min_fold_$f $mc $f $MAX_FOLD shared $MASK_BED"); done
    for f in ${GRID_MAX_FOLDS//,/ }; do runs+=("max_fold_$f $mc $MIN_FOLD $f shared $MASK_BED"); done
    runs+=("per_sample_masks $mc $MIN_FOLD $MAX_FOLD per-sample $MASK_BED")
    runs+=("no_accessibility_mask $mc $MIN_FOLD $MAX_FOLD shared $g/minimal_mask.bed")

    local assess=() spec r rmc rmin rmax rmode rbase
    for spec in "${runs[@]}"; do
        read -r r rmc rmin rmax rmode rbase <<<"$spec"
        log "grid run $r"
        mkdir -p "$g/$r"
        depthmask "$g/$r/mask" "$rmin" "$rmax" "$rmode" "$rbase" > "$g/$r/depthmask.txt"
        # shellcheck source=/dev/null
        ( source "$g/$r/mask.depthmask.params"
          run_grenedalf "$g/$r" "$TOTAL_MASK_BED" "$SAMPLE_MASK_LIST" "$rmc" "$MIN_DEPTH" "$MAX_DEPTH" \
              "$WINDOW" "$STRIDE" )
        assess+=(--run "$r" "$g/$r/diversity.csv" "$g/$r/fst.csv")
    done
    for n in "${NAMES[@]}"; do assess+=(--depth "$n" "$OUT/depth/$n.mq$MAPQ.regions.bed.gz"); done
    "${OPTIMISER[@]}" assess "${assess[@]}" --min-snps "$MIN_SNPS" --tolerance "$TOLERANCE" --out "$P"
}

step_params() {
    load_params
    local r=${RECOMMENDED_RUN:-recommended}
    # shellcheck source=/dev/null
    source "$OUT/grid/$r/mask.depthmask.params"
    local mc=$MIN_COUNT
    [[ $r == min_count_* ]] && mc=${r#min_count_}
    cat > "$OUT/params.sh" <<EOF
# grenedalf settings chosen by grenedalf_optimiser ($(date '+%F')), run: $r
# Technical criteria only; see $P.assess_runs.tsv and the other $P.* tables.
TOTAL_MASK_BED=$TOTAL_MASK_BED
SAMPLE_MASK_LIST=$SAMPLE_MASK_LIST
MIN_COUNT=$mc
MIN_DEPTH=$MIN_DEPTH
MAX_DEPTH=$MAX_DEPTH
POOL_SIZES=$POOL_SIZE_ARG
X_SEPARATE=${X_SEPARATE:-0}
X_POOL_SIZES=${X_POOL_SIZES:-}
MITO_SEPARATE=${MITO_SEPARATE:-0}
MAPQ=$MAPQ
BASEQ=$BASEQ
WINDOW=$WINDOW
STRIDE=$STRIDE
MIN_SNPS=$MIN_SNPS
# Ready-to-run commands, one run per chromosome class: $OUT/grenedalf_commands.sh
# grenedalf ... --filter-mask-total-bed \$TOTAL_MASK_BED \\
#     \${SAMPLE_MASK_LIST:+--filter-mask-sample-bed-list \$SAMPLE_MASK_LIST} \\
#     --filter-sample-min-count \$MIN_COUNT \\
#     --filter-sample-min-read-depth \$MIN_DEPTH --filter-sample-max-read-depth \$MAX_DEPTH \\
#     --sam-min-map-qual \$MAPQ --sam-min-base-qual \$BASEQ \\
#     --window-interval-width \$WINDOW --window-interval-stride \$STRIDE --window-average-policy valid-loci
# Drop windows with fewer than \$MIN_SNPS passed SNPs (total.passed) before interpreting them.
EOF
    log "wrote $OUT/params.sh"
    cat "$OUT/params.sh"
    "${OPTIMISER[@]}" commands --prefix "$P" --out-dir "$OUT" --fai "$FAI" --bam "${BAMS[@]}" --x-regex "$X_REGEX" \
        --total-mask "$TOTAL_MASK_BED" --sample-masks "$SAMPLE_MASK_LIST" --pool-size-arg "$POOL_SIZE_ARG" \
        --min-count "$mc" --min-depth "$MIN_DEPTH" --max-depth "$MAX_DEPTH" --window "$WINDOW" --stride "$STRIDE" \
        --min-snps "$MIN_SNPS" --map-qual "$MAPQ" --base-qual "$BASEQ" --fst-method "$FST_METHOD" --threads "$THREADS"
}

step_report() {
    "${OPTIMISER[@]}" report --prefix "$P" --out-dir "$OUT" --out "$OUT/optimiser_report.html" \
        ${REPORT_TITLE:+--title "$REPORT_TITLE"}
}

case "$STEP" in
    mask|profile|thresholds|chromclass|quantize|depthmask|pilot|grid|params|report) "step_$STEP" ;;
    all) for s in mask profile thresholds chromclass quantize depthmask pilot grid params report; do log "== $s"; "step_$s"; done ;;
    *) echo "unknown step: $STEP" >&2; exit 1 ;;
esac
