# grenedalf-tune

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
- **whether the statistics depend on technical covariates**: θπ, Watterson's θ, Tajima's
  D and FST should not track read depth, depth imbalance between pools, or the callable
  fraction of a window. When they do, the filter is letting artefact through.

Each filter is loosened as far as it can go without adding such a dependence, compared
with the strictest setting tried.

## What it does

| Step | Tools | What it decides |
|---|---|---|
| `mask` | GenMap, mosdepth | **Accessibility mask**, the same for all pools: reference gaps; low mappability (most k-mers covering a base are not unique, so reads there get MAPQ 0 and the site looks "missing"); collapsed repeats (consensus depth > 2× the chromosome median: extra genomic copies stacked on one reference copy); ambiguous mapping (consensus: most reads below the MAPQ cut-off, which catches repeats GenMap cannot see because the reference lacks them); sequence absent from every pool. Optionally soft-masked repeats. Depth anomalies in **single** pools (deletions, CNVs) are written to `tune.pool_specific_depth.tsv` and **not masked**, because they may be biology. |
| `profile` | mosdepth, samtools | Each pool's depth distribution **inside the accessible genome only**, so repeats and gaps do not distort it. |
| `thresholds` | – | **Per-pool depth bounds**, autosomes and X separately (the X runs at lower depth when males are in the pool): fold × that pool's median. A single 30–300 for every pool removes a different fraction of the genome from a 60× pool than from a 150× pool. **Minimum allele count**: the smallest count that sequencing errors alone rarely reach at the pool's depth, given the base-quality cut-off. |
| `quantize`, `depthmask` | mosdepth | Applies the per-pool bounds as BED masks. grenedalf has only one global depth filter, so per-pool bounds have to go in as masks. `shared` mode (default) merges them, so every pool is called at exactly the same sites. `per-sample` keeps more sites at the cost of uneven missingness; the grid tests both. |
| `pilot` | grenedalf | One chromosome with the recommended settings. **Linkage scale** (`ld`): how far apart two SNPs have to be before their allele frequencies vary independently across pools (an LD decay for pool-seq, no phasing needed). **Window width**: wide enough for `TARGET_SNPS` SNPs in 90% of windows, and no narrower than the linkage scale. Also checks that grenedalf applied the mask the right way round. |
| `grid` | grenedalf | Pilot chromosome rerun with one setting changed at a time: minimum count ±1, lower and upper depth folds, per-sample masks, and no accessibility mask. `assess` scores each run, and the most permissive run that adds no artefact is recommended. |
| `params` | – | `OUT/params.sh` for the genome-wide run. |

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
(≥ 0.6) on `PATH` or set in the config; bash (3.2 or later).

```bash
cp config.example.sh my_config.sh    # edit paths, pools, pilot chromosome
bin/grenedalf_tune.sh my_config.sh all
```

Or run one step at a time (`mask profile thresholds quantize depthmask pilot grid params`).
Each step reads what the previous ones wrote in `OUT/`. On the farm, for example:

```bash
bsub -n 8 -R "span[hosts=1] select[mem>16000] rusage[mem=16000]" -M 16000 \
    -o tune.%J.o -e tune.%J.e bin/grenedalf_tune.sh my_config.sh all
```

The memory goes on GenMap indexing and on the mask step, which holds one chromosome as
per-base arrays and the GenMap bedGraph.

## Outputs (`OUT/tune.*`)

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

Downstream, drop windows with fewer than `MIN_SNPS` passed SNPs. Note that
`grenedalf frequency` applies no sample filters, so a script reading the frequency
table has to apply `TOTAL_MASK_BED`, `MIN_COUNT` and the depth bounds itself (the `ld`
step does this).

## Check before trusting the first run

These depend on grenedalf behaviour that could not be tested here (written against the
0.6.3 documentation):

- `check_grenedalf` stops the run if `--filter-mask-total-bed`, `--filter-mask-sample-bed-list`
  or `--filter-region` is missing from `grenedalf diversity --help`.
- `windows.txt` reports whether `total.masked` in the pilot table matches the BED. If it
  says *WRONG WAY ROUND*, grenedalf treats the BED as the regions to keep.
- The per-sample mask list is written as `sample<TAB>bed`. Check this against
  `grenedalf diversity --help` before using the `per_sample_masks` grid run.
- Pool names are taken from the BAM file names without the extension, as grenedalf
  names samples. They must match the column prefixes in grenedalf's tables.

## Limits

- The depth bounds use mosdepth depth (MAPQ filtered, not base-quality filtered), so they
  sit a few percent above grenedalf's own depth at the same site.
- mosdepth writes its depth distribution to two decimals. Medians are reliable, but the
  1% and 99% quantiles are coarse.
- Some dependence of diversity on depth is real: divergent haplotypes map less well and
  lose depth. The tool therefore compares runs with the strictest run, not with zero.
- The linkage scale assumes the pools differ in frequency at least somewhat. With
  near-identical pools there is no excess r² to measure, and the window width falls back
  to SNP density alone.

## Test

```bash
test/run_test.sh
```

Simulates GenMap, mosdepth, samtools and grenedalf outputs, with known features built
in: a gap, a soft-masked block, a non-unique block, a collapsed repeat, an ambiguous-mapping
block, a deletion in one pool, X at half depth, a mitochondrion, a 2 kb linkage scale, and
a run where θπ follows depth. It then checks that each one is found. The shell driver
itself has not been run against the real tools yet.
