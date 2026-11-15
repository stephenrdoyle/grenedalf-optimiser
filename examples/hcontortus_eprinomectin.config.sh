# grenedalf-optimiser config: H. contortus eprinomectin L3 pools (WBPS18 reference), Sanger farm.
# Set BAM_DIR, then:  bin/grenedalf_optimiser.sh CONFIG run
# (creates the conda environment if needed, checks the inputs, submits one LSF job)

CONDA_ENV=grenedalf-optimiser         # created if missing; activated by every step
SCHEDULER=auto                        # LSF when bsub is available
LSF_QUEUE=long
LSF_MEM=32000

PROJECT=/nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/project_eprinomectin-resistance
REF=$PROJECT/01_references/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa
FAI=$REF.fai
BAM_DIR=/path/to/bams                 # <- the directory run_grenedalf.sh was run from
BAMS=("$BAM_DIR"/ARA_L3.bam "$BAM_DIR"/BET_L3.bam "$BAM_DIR"/BUN_L3.bam
      "$BAM_DIR"/CHI_L3.bam "$BAM_DIR"/LUC_L3.bam "$BAM_DIR"/MOU_L3.bam)
OUT=$PROJECT/grenedalf_optimiser_out
THREADS=8

POOL_SIZE=400                         # 200 diploid L3 per pool
MALE_FRACTION=                        # L3 sex ratio not known: estimate from X depth
CHROM_REGEX='hcontortus_chr[1-5X]_'
X_REGEX='chrX'
MITO_REGEX='mitochondrion'
REGION=                               # whole genome

KMER=150                              # read length: check (zcat a fastq | head -2 | tail -1 | wc -c)
MAPQ=30
BASEQ=30

REPORT_TITLE="H. contortus eprinomectin pools"
