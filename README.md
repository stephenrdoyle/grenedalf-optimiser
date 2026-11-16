# grenedalf-optimiser

Chooses grenedalf filter settings for pool-seq data from the technical properties of the
BAMs, **before any biological signal is looked at**.

Every threshold in a grenedalf run (minimum and maximum depth, minimum allele count,
masking, window size) changes what can be seen downstream. If the thresholds are set by
eye, or tuned until a result looks right, the analysis is fitted to its own output. This
tool never uses population labels, phenotypes, or the position of peaks. It judges
settings only on:

- **how much of the genome stays callable**: more data means more chance of discovery
- **whether every pool is called at the same sites**: uneven missingness between pools
  looks like differentiation
- **whether the statistics depend on technical covariates**: θπ, Watterson's θ and FST should not track read depth, depth imbalance between pools, or the callable
  fraction of a window. When they do, the filter is letting artefact through.

Each filter is loosened as far as it can go without adding such a dependence, compared
with the strictest setting tried.

## What it does

| Step | Tools | What it decides |
|---|---|---|
| `mask` | GenMap, mosdepth | **Accessibility mask**, the same for all pools: reference gaps; low mappability (most k-mers covering a base are not unique, so reads there get MAPQ 0 and the site looks "missing"); collapsed repeats (consensus depth > 2× the chromosome median: extra genomic copies stacked on one reference copy); ambiguous mapping (consensus: most reads below the MAPQ cut-off, which catches repeats GenMap cannot see because the reference lacks them); sequence absent from every pool. Optionally soft-masked repeats. Depth anomalies in **single** pools (deletions, CNVs) are written to `optimiser.pool_specific_depth.tsv` and **not masked**, because they may be biology. |
| `profile` | mosdepth, samtools | Each pool's depth distribution **inside the accessible genome only**, so repeats and gaps do not distort it. |
| `thresholds` | – | **Per-pool depth bounds**, autosomes and X separately (the X runs at lower depth when males are in the pool): fold × that pool's median; and, as an alternative the grid tests (`fitted_bounds`), the 0.5% and 99.5% points of a negative binomial fitted to the pool's main coverage, with a separate low-depth component, so the bounds follow each pool's actual spread. A single 30–300 for every pool removes a different fraction of the genome from a 60× pool than from a 150× pool. **Minimum allele count**: the smallest count that sequencing errors alone rarely reach at the pool's depth, given the base-quality cut-off. With the **pool size**: the effective number of chromosomes sampled at each pool's depth (n_eff = d·n / (d + n − 1)), and the lowest allele frequency the minimum count can see against one copy in the pool. |
| `chromclass` | – | **X and mitochondrion**: each pool's X depth against its autosomes, taken peak to peak (kernel density mode) over well-mapped bins, with a block-bootstrap 95% interval. Pools with males have X depth 1 − m/2 of autosome depth (m = male fraction), so the X is sampled from fewer copies and needs its own pool size, pool size × (1 − m/2). A ratio outside 0.5–1 is impossible for a sex ratio (on *H. contortus* field data the X had more depth than the autosomes, from reference divergence): then no pool's estimate is used, and the male fraction comes from `MALE_FRACTION`, or 0.5 is assumed, with a warning. The X gets its own grenedalf run whenever its pool size differs. The mitochondrion (haploid, many copies, one lineage per individual) always gets its own run, with per-pool depth bounds applied as per-pool masks, and a pool size of the number of individuals. |
| `variance` | – | **What kind of depth variation there is**: per-base spread against pure read-sampling noise, how much of the variance persists in 1 kb and 10 kb means (regional blocks vs scattered noise), and, for the sites a lower bound would cut, whether they lie in blocks and whether they are low in one pool only (population divergence, deletions) or in most pools (reference, technical). Also the 1 kb depth correlation between pools. |
| `quantize`, `depthmask` | mosdepth | Applies the per-pool bounds as BED masks. grenedalf has only one global depth filter, so per-pool bounds have to go in as masks. `shared` mode (default) merges them, so every pool is called at exactly the same sites. `per-sample` keeps more sites at the cost of uneven missingness; the grid tests both. |
| `pilot` | grenedalf | The whole genome with the recommended settings. **Linkage scale** (`ld`): how far apart two SNPs have to be before their allele frequencies vary independently across pools (an LD decay for pool-seq, no phasing needed). **Window width**: wide enough for `TARGET_SNPS` SNPs in 90% of windows, and no narrower than the linkage scale. Also checks that grenedalf applied the mask the right way round. |
| `grid` | grenedalf | The whole genome rerun with one setting changed at a time: minimum count ±1, lower and upper depth folds, fitted bounds, per-sample masks (and both together), and, as a check, no accessibility mask. `assess` compares each run with the reference run (the one keeping least data), **statistic by statistic**: a run is accepted only if none of θπ, θW or FST depends more on depth, depth imbalance or callable fraction than in the reference, beyond a tolerance. Then `combine` builds one more run from the best accepted choice of each setting, judged by the same rule. The accepted run keeping the most data (SNPs usable per pair of pools, counted the same way for every run) is recommended. |
| `params` | – | `OUT/params.sh`, and `OUT/grenedalf_commands.sh`: ready-to-run grenedalf commands with every chosen value filled in, one run per chromosome class (nuclear, or autosomes + X, plus mitochondrion). |
| `report` | – | `OUT/optimiser_report.html`: every step's results on one page (see below). |

### The linkage scale

Pool-seq reads are not phased, but linked SNPs share their frequency differences
between pools: if pool A carries more of a haplotype than pool B, every SNP on that
haplotype differs between A and B the same way. So the correlation (r², across pools)
of two SNPs' centred frequencies is high when they are close, and falls with distance
to a baseline for unlinked SNPs. That baseline is about 1/(pools − 1) plus sampling
noise, and the tool estimates it from SNP pairs on different chromosomes or far apart.
`LD_D50` and `LD_D10` are the distances at which the excess r² over baseline has fallen
to 50% and 10%. SNPs closer than `LD_D10` are not independent observations, so narrower
windows mostly repeat their neighbours. `width / LD_D10` is roughly the number of
independent segments in each window. That number, scaled to the whole genome, also
gives the effective number of independent tests downstream.

## Running

Requirements: python3 with numpy and pandas; mosdepth, samtools, GenMap and grenedalf
(0.6.3); bash. All of them are in the conda environment in `environment.yml` (Linux):

```bash
bin/grenedalf_optimiser.sh my_config.sh setup     # creates (or updates) it; run does this too if it is missing
```

Every step then activates the environment itself (`CONDA_ENV` in the config: a name or a
path; `none` uses whatever is on PATH, e.g. modules), so no `conda activate` is needed,
interactively or in a cluster job.

```bash
cp config.example.sh my_config.sh    # edit paths, pools, pool sizes
bin/grenedalf_optimiser.sh my_config.sh run
```

That is the only command you need, including after an update (`git pull`) or a config
change: run it again. `run` creates the conda environment if it is missing, runs the
`check` step (and stops there if anything is wrong), then runs the whole analysis,
submitted as one LSF job when `bsub` is available (queue `LSF_QUEUE`, memory `LSF_MEM`,
cores `THREADS`) or in the current shell otherwise (`SCHEDULER=local` forces this). Logs go
to `OUT/logs/`, the report to `OUT/optimiser_report.html`, the recommended commands to
`OUT/grenedalf_commands.sh`.

**Reruns are cheap.** The slow steps (GenMap, mosdepth, samtools stats, every grenedalf run)
record a signature of their inputs (settings, mask checksums, BAM size and date) and are
reused while it matches; the fast steps always rerun. So a rerun after an update only redoes
what the update affects. `FORCE=1 bin/grenedalf_optimiser.sh my_config.sh run` redoes
everything.

`check` on its own takes a minute or two and catches most setup problems: tools and
python packages, grenedalf's options, BAMs and their indexes, BAM contig names against the
reference, and that every contig named in `SEX_CHROMS` / `MITO_CHROMS` exists.

### Naming the sex chromosome and mitochondrion

Name them explicitly, comma separated if there are several contigs:

```bash
SEX_CHROMS=hcontortus_chrX_Celeg_TT_arrow_pilon
MITO_CHROMS=mitochondrion
```

Every other contig is treated as autosomal, and if neither is set, everything is. The
sex chromosome is the hemizygous one (X, or Z in ZW species): its depth and pool size are
tested against the autosomes. Do not list a Y or W here; it is absent from half the
individuals and does not fit that model. Mitochondrial contigs are left out of the
nuclear analysis and get their own run. `CHROM_REGEX` can optionally restrict the
analysis to matching contigs (default: all).

Individual steps exist for debugging or rebuilding one part (e.g. `bin/grenedalf_optimiser.sh my_config.sh report`); you never need them for a normal run. In order: `setup check mask profile thresholds chromclass variance quantize depthmask pilot grid params report`; each reads what the previous ones wrote in `OUT/`.

The `pilot` and `grid` steps run on the **whole genome** by default, so the settings
are judged on every chromosome rather than on one that may not be typical. The grid
is about ten genome-wide grenedalf diversity + FST runs, so it is the slow step
(hours, not minutes). For a quick trial, set `REGION` to one chromosome.

### How the grid affects the recommendation

The grid decides the **minimum count**, the **depth bounds** (lower and upper fold, or fitted
to each pool's coverage) and the **mask mode** (shared or per sample). It is a controlled
experiment: from a baseline, one setting changes at a time, and each run is scored on

- **data kept**: SNPs usable per pair of pools, the data FST works with, counted the same
  way whatever the mask mode (grenedalf's own totals count a site if any pool can use it,
  which would flatter per-sample masks);
- **artefact**: for θπ, θW and FST separately, how strongly the statistic follows read depth,
  depth imbalance between the pools of a pair, and the callable fraction of each window.
  None of these should matter.

The reference is the run keeping least data. A run is accepted only if **no statistic** gets
more dependent on any covariate than in the reference, beyond the tolerance (0.05), so an
improvement in θπ cannot hide a deterioration in FST. The best accepted choice of each setting
is then combined into one more run and judged the same way. Of all accepted runs, the one
keeping the most data is recommended, and `params.sh`, `grenedalf_commands.sh` and the report
use it. `no_accessibility_mask` is a check, never recommended: it shows what the mask does.

The grid does **not** change the accessibility mask, the window width (from SNP
density and the linkage scale), the X / mitochondrion decision, or the read filters
(MAPQ, base quality).

On the Sanger farm, `examples/hcontortus_eprinomectin.config.sh` is a filled-in config for the
eprinomectin pools (set `BAM_DIR`):

```bash
git clone https://github.com/stephenrdoyle/grenedalf-optimiser.git && cd grenedalf-optimiser
bin/grenedalf_optimiser.sh examples/hcontortus_eprinomectin.config.sh run
```

The memory goes on GenMap indexing and on the mask step, which holds one chromosome as
per-base arrays and the GenMap bedGraph.

## Outputs (`OUT/optimiser.*`)

| File | Read it for |
|---|---|
| `mask_summary.tsv` | Fraction of each chromosome masked, per reason |
| `mask_reasons.bed` | Where each masked block is and why |
| `pool_specific_depth.tsv` | Single-pool depth anomalies, left unmasked |
| `depth_thresholds.tsv` | Per pool and chromosome class: median, quantiles, bounds, fraction kept |
| `min_count.tsv` | Expected error-only sites per Mb for counts 1–8, per pool |
| `depthmask_summary.tsv` | Accessible bases lost to each pool's bounds, and to the shared mask |
| `ld_decay.tsv`, `ld.params` | r² against distance; LD_D50, LD_D10 |
| `windows.txt`, `window_widths.tsv` | SNP density, callable fraction, recommended width, mask check |
| `assess_runs.tsv` | Per grid run: callable bp, SNPs, median \|ρ\| of each statistic with depth, depth imbalance and callable fraction, artefact score, recommendation |
| `assess_stats.tsv` | The same per statistic column (each pool, each pair) |
| `OUT/params.sh` | Settings and masks for the full grenedalf run |
| `OUT/optimiser_report.html` | The report |

Downstream, drop windows with fewer than `MIN_SNPS` passed SNPs. Note that
`grenedalf frequency` applies no sample filters, so a script reading the frequency
table has to apply `TOTAL_MASK_BED`, `MIN_COUNT` and the depth bounds itself (the `ld`
step does this).

## Check before trusting the first run

These depend on grenedalf behaviour that could not be tested here (written against the
0.6.3 documentation):

- `check` stops the run if any grenedalf option the tool uses is missing from
  `grenedalf diversity --help` (written against grenedalf 0.6.3).
- `windows.txt` reports whether `total.masked` in the pilot table matches the BED. If it
  says *WRONG WAY ROUND*, grenedalf treats the BED as the regions to keep.
- The per-sample mask list is written as `sample<TAB>bed`. Check this against
  `grenedalf diversity --help` before using the `per_sample_masks` grid run.
- The male fraction estimated from X depth assumes the X maps like the autosomes. If the X
  differs in GC content or mappability, its depth (and so the estimate) shifts. When the
  sex ratio is known from the design, set `MALE_FRACTION` (for example 0.5): the known value
  then sets the X pool size, and the report flags pools whose X depth does not match it.
- Pool names are taken from the BAM file names without the extension, as grenedalf
  names samples. They must match the column prefixes in grenedalf's tables.

## Limits

- **Tajima's D is not computed**, in the grid or in `grenedalf_commands.sh` (`--no-tajima-d`):
  grenedalf computes it only with a minimum count of exactly 2, which would mean a different
  filter from the other statistics.
- The depth bounds use mosdepth depth (MAPQ filtered, not base-quality filtered), so they
  sit a few percent above grenedalf's own depth at the same site.
- mosdepth writes its depth distribution to two decimals. Medians are reliable, but the
  1% and 99% quantiles are coarse.
- Some dependence of diversity on depth is real: divergent haplotypes map less well and
  lose depth. The tool therefore compares runs with the strictest run, not with zero.
- The linkage scale assumes the pools differ in frequency at least somewhat. With
  near-identical pools there is no excess r² to measure, and the window width falls back
  to SNP density alone.

## Report

`OUT/optimiser_report.html` is a single self-contained file (inline SVG; no external
scripts, fonts or network access), so it can be copied off the cluster and opened in
any browser. It follows the light/dark setting, with a toggle. Sections:

- **Recommended settings**: headline numbers and the `params.sh` to use
- **Accessibility mask**: fraction of each chromosome masked, by reason, and the
  single-pool depth anomalies that were left in; the summary states how much more the shared
  depth mask removes than any single pool
- **Where the mask falls**: a karyoplot per chromosome (the ten longest analysed contigs),
  one row per mask reason plus the final mask, then one row per pool for sites outside
  that pool's depth bounds (long blocks in one pool: divergent haplotypes, deletions, copy-number changes)
- **Depth per pool**: each pool's depth curve with its kept range, and how much of the
  accessible genome each pool's bounds (and the shared mask) cost
- **Minimum allele count**: expected error-only sites per Mb at each count
- **Linkage scale**: r² against distance, with the baseline and decay distances
- **Window size**: windows reaching the SNP target by width, the linkage scale, and
  the mask orientation check
- **Filter grid**: data kept and artefact score for every run, and the median |ρ| of
  each statistic with each technical covariate

Hovering over any mark shows its values. Tables sit under each chart. An example built
from the synthetic test data is in [`docs/example_report.html`](docs/example_report.html).
To rebuild the report on its own:

```bash
bin/grenedalf_optimiser.sh my_config.sh report
# or: python3 bin/grenedalf_optimiser.py report --prefix OUT/optimiser --out report.html
```

## Test

```bash
test/run_test.sh
```

Simulates GenMap, mosdepth, samtools and grenedalf outputs, with known features built
in: a gap, a soft-masked block, a non-unique block, a collapsed repeat, an ambiguous-mapping
block, a deletion in one pool, X at half depth, a mitochondrion, a 2 kb linkage scale, and
a run where θπ follows depth. It then checks that each one is found. The shell driver
itself has not been run against the real tools yet.
