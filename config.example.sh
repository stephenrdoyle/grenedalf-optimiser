# grenedalf_optimiser config (sourced by bin/grenedalf_optimiser.sh). Copy and edit.
# Example values: H. contortus eprinomectin pools, WBPS18 reference.

REF_DIR=/nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/project_eprinomectin-resistance/01_references
REF=$REF_DIR/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa
FAI=$REF.fai
BAM_DIR=/path/to/bams
BAMS=("$BAM_DIR"/ARA_L3.bam "$BAM_DIR"/BET_L3.bam "$BAM_DIR"/BUN_L3.bam
      "$BAM_DIR"/CHI_L3.bam "$BAM_DIR"/LUC_L3.bam "$BAM_DIR"/MOU_L3.bam)
OUT=grenedalf_optimiser_out
THREADS=8

POOL_SIZE=400                         # haploid genomes per pool (200 diploid L3)
POOL_SIZES=                           # or per pool, if they differ: ARA_L3=400,BET_L3=380,...
MALE_FRACTION=                        # known fraction of males (0.5 for 50:50; or POOL=m,...); empty = estimate from X depth
SEX_CHROMS=hcontortus_chrX_Celeg_TT_arrow_pilon   # sex chromosome contig(s), comma separated; empty = none
MITO_CHROMS=mitochondrion                       # mitochondrial contig(s): own run; empty = none
# every other contig is treated as autosomal
REGION=                               # empty = whole genome (recommended); one chromosome for a quick trial

# Read filters (kept fixed; reported, not tuned)
MAPQ=30
BASEQ=30

# Mappability
KMER=150             # read length
GENMAP_ERRORS=2
MAPPABILITY_BEDGRAPH=   # reuse an existing GenMap bedGraph instead of running GenMap
SOFTMASK=0           # 1 = also mask soft-masked (lower-case) repeats; usually not needed when mappability is masked
BIN=1000             # mosdepth bin size for the consensus depth mask

# Starting depth bounds (x pool median) and the alternatives tried in the grid
MIN_FOLD=0.5
MAX_FOLD=2.0
ABS_MIN_DEPTH=10
GRID_MIN_FOLDS=0.33,0.67
GRID_MAX_FOLDS=1.5,3.0

# Windows
TARGET_SNPS=50       # SNPs per window for a precise window estimate
PILOT_WINDOW=5000
LD_MAX_DIST=100000

FST_METHOD=unbiased-nei
TOLERANCE=0.05       # artefact increase over the strictest run still accepted

# Tools: the conda environment from environment.yml (bin/grenedalf_optimiser.sh CONFIG setup)
# is activated by every step. CONDA_ENV=none uses whatever is on PATH, e.g. farm modules.
CONDA_ENV=grenedalf-optimiser
SCHEDULER=auto       # run: LSF when bsub is available, otherwise here; or lsf / local
LSF_QUEUE=long
LSF_MEM=32000        # MB; cores = THREADS
# Commands, if not the ones in the environment
GRENEDALF=grenedalf
MOSDEPTH=mosdepth
SAMTOOLS=samtools
GENMAP=genmap
