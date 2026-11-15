#!/usr/bin/env python3
"""
grenedalf_optimiser: choose grenedalf filter settings from the technical properties of
the pool-seq BAMs, before any biological signal is looked at.

Nothing here uses population labels, phenotypes, or the location of peaks. The
settings are judged only on technical grounds: how much of the genome stays
callable, whether every pool is callable at the same sites, and whether the
statistics (theta pi, Watterson's theta, Tajima's D, FST) depend on read depth or on
how much of a window could be called. A filter is loosened only as far as it can be
without adding such a dependence, so as little data as possible is thrown away.

Subcommands (run in this order; grenedalf_optimiser.sh runs the external tools between them):

  mask        accessibility mask shared by all pools: reference gaps, low mappability
              (GenMap), optionally soft-masked repeats, and regions whose depth is
              abnormal in the consensus of all pools (collapsed repeats: high depth;
              ambiguous mapping: most reads below the MAPQ cut-off; absent sequence:
              low depth). Depth anomalies in single pools are reported, not masked.
  thresholds  per-pool depth bounds (autosomes and X separately) from each pool's
              depth distribution inside the accessible genome, the minimum allele
              count implied by the base-quality cut-off, and the mosdepth --quantize
              cut points needed to build depth masks.
  depthmask   per-pool depth masks from mosdepth --quantize output; either merged
              into one mask (shared: every pool is called at the same sites) or kept
              per sample (grenedalf --filter-mask-sample-bed-list).
  windows     window width from SNP density in a pilot grenedalf table (the smallest
              width at which most windows hold enough SNPs for a precise estimate) and,
              with --ld, no narrower than the linkage scale from the ld step.
              Also checks that grenedalf applied the mask the right way round.
  chromclass  do the X and mitochondrion behave like autosomes (depth ratio, bootstrap
              interval)? If not, they get their own grenedalf runs and pool sizes.
  ld          linkage scale: how fast the correlation of allele frequencies across
              pools decays with distance between SNPs (an LD decay for pool-seq).
  assess      compares grenedalf runs over a grid of settings on technical criteria
              only, and names the most permissive setting that adds no artefact.
  commands    grenedalf_commands.sh: ready-to-run commands with the chosen settings.
  report      one HTML file with every step's results (self-contained, opens offline).
"""

import argparse
import gzip
import math
import os
import re
import sys

import numpy as np
import pandas as pd

# Mask reasons as bit flags, so a base can carry more than one
FLAGS = {
    "excluded": 1,         # chromosome not matching --chrom-regex (e.g. mitochondrion, unplaced)
    "gap": 2,              # N in the reference
    "softmask": 4,         # lower case in the reference (only with --softmask)
    "low_mappability": 8,  # most k-mers overlapping the base are not unique (GenMap)
    "collapsed": 16,       # consensus depth well above the chromosome median
    "ambiguous": 32,       # consensus: most reads have MAPQ below the cut-off
    "low_depth": 64,       # consensus depth well below the chromosome median
}

STAT_SUFFIXES = (".theta_pi", ".theta_watterson", ".tajimas_d", ".fst")


# ---------------------------------------------------------------------------
# Shared helpers

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def read_fai(path):
    return pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"],
                       dtype={"chrom": str, "length": np.int64})


def read_bed(path, ncol=3):
    names = ["chrom", "start", "end", "name"][:ncol]
    return pd.read_csv(path, sep="\t", header=None, usecols=range(ncol), names=names, comment="#",
                       dtype={"chrom": str, "start": np.int64, "end": np.int64})


def write_bed(df, path):
    df.to_csv(path, sep="\t", header=False, index=False)


def iter_fasta(path):
    """Yield (name, sequence as bytes) for each record."""
    name, chunks = None, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).encode()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).encode()


def runs(arr):
    """Start, end (exclusive) and value of each run of equal values."""
    if arr.size == 0:
        return np.array([], np.int64), np.array([], np.int64), arr[:0]
    change = np.flatnonzero(arr[1:] != arr[:-1]) + 1
    starts = np.concatenate(([0], change))
    ends = np.concatenate((change, [arr.size]))
    return starts, ends, arr[starts]


def merge_intervals(df):
    """Merge overlapping or adjacent intervals per chromosome."""
    if df.empty:
        return df[["chrom", "start", "end"]].copy()
    out = []
    for chrom, g in df.sort_values(["chrom", "start"]).groupby("chrom", sort=False):
        s, e = g["start"].to_numpy(), g["end"].to_numpy()
        emax = np.maximum.accumulate(e)
        idx = np.flatnonzero(np.concatenate(([True], s[1:] > emax[:-1])))
        out.append(pd.DataFrame({"chrom": chrom, "start": s[idx], "end": np.maximum.reduceat(e, idx)}))
    return pd.concat(out, ignore_index=True)


def bp(df):
    return int((df["end"] - df["start"]).sum()) if len(df) else 0


def bp_outside(a, b):
    """Bases of a not covered by b (both merged)."""
    return bp(merge_intervals(pd.concat([a, b]))) - bp(b)


def coverage_prefix(intervals, length):
    """Prefix sum of a 0/1 coverage track, so coverage of [s, e) is cs[e] - cs[s]."""
    d = np.zeros(length + 1, np.int64)
    if len(intervals):
        s = intervals["start"].clip(0, length).to_numpy()
        e = intervals["end"].clip(0, length).to_numpy()
        d += np.bincount(s, minlength=length + 1) - np.bincount(e, minlength=length + 1)
    cov = np.cumsum(d[:length]) > 0
    return np.concatenate(([0], np.cumsum(cov, dtype=np.int64)))


def read_mosdepth_bins(path, column):
    return pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3],
                       names=["chrom", "start", "end", column],
                       dtype={"chrom": str, "start": np.int64, "end": np.int64, column: np.float64})


def normalise_by_chrom(df, col):
    """Depth relative to the median of the chromosome's covered bins."""
    med = df[col].where(df[col] > 0).groupby(df["chrom"]).transform("median")
    return df[col] / med


def parse_pool_sizes(spec, names):
    """Pool size (haploid genomes) per pool, from "400" or "POOL=400,POOL2=380"."""
    if spec is None or str(spec).strip() == "":
        return {}
    spec = str(spec).strip()
    if "=" not in spec:
        return {n: float(spec) for n in names}
    sizes = dict((k.strip(), float(v)) for k, v in (x.split("=", 1) for x in spec.split(",") if x.strip()))
    missing = [n for n in names if n not in sizes]
    if missing:
        sys.exit(f"no pool size for {', '.join(missing)} in '{spec}'")
    return sizes


def write_params(path, params):
    with open(path, "w") as fh:
        for k, v in params.items():
            fh.write(f"{k}={v}\n")


# ---------------------------------------------------------------------------
# mask

def unique_fraction(length, gm, kmer, threshold, step=5_000_000):
    """Per base, the fraction of k-mers overlapping it that are unique.

    GenMap reports the mappability of the k-mer starting at each position. A read
    covering base j starts at one of j-k+1..j, so the base is reliably mapped when
    most of those k-mers are unique."""
    n_starts = max(length - kmer + 1, 1)
    uniq = gm[gm["value"] >= threshold]
    d = np.zeros(length + 1, np.int64)
    if len(uniq):
        s = uniq["start"].clip(0, length).to_numpy()
        e = uniq["end"].clip(0, length).to_numpy()
        d += np.bincount(s, minlength=length + 1) - np.bincount(e, minlength=length + 1)
    u = np.cumsum(d[:n_starts]) > 0
    cs = np.concatenate(([0], np.cumsum(u, dtype=np.int64)))
    frac = np.empty(length, np.float32)
    for a in range(0, length, step):
        j = np.arange(a, min(a + step, length))
        lo = np.clip(j - kmer + 1, 0, n_starts - 1)
        hi = np.clip(j, 0, n_starts - 1) + 1
        frac[a:a + j.size] = (cs[hi] - cs[lo]) / (hi - lo)
    return frac


def consensus_depth_flags(depth_args, max_norm, min_norm, min_ratio):
    """Bins flagged from the consensus (median across pools) of normalised depth."""
    wide = None
    pools = []
    for name, f_mq0, f_mq30 in depth_args:
        pools.append(name)
        a = read_mosdepth_bins(f_mq0, f"{name}.mq0")
        b = read_mosdepth_bins(f_mq30, f"{name}.mq30")
        ab = a.merge(b, on=["chrom", "start", "end"], how="inner", validate="one_to_one")
        if len(ab) != len(a):
            sys.exit(f"mosdepth bins differ between {f_mq0} and {f_mq30}")
        if wide is None:
            wide = ab
        else:
            n = len(wide)
            wide = wide.merge(ab, on=["chrom", "start", "end"], how="inner", validate="one_to_one")
            if len(wide) != n:
                sys.exit(f"mosdepth bins for {name} differ from the other pools (use the same --by)")
    norm0 = pd.DataFrame({p: normalise_by_chrom(wide, f"{p}.mq0") for p in pools})
    ratio = pd.DataFrame({p: (wide[f"{p}.mq30"] / wide[f"{p}.mq0"]).where(wide[f"{p}.mq0"] > 0)
                          for p in pools})
    wide["cons_norm_mq0"] = norm0.median(axis=1)
    wide["cons_mq_ratio"] = ratio.median(axis=1)
    wide["collapsed"] = wide["cons_norm_mq0"] > max_norm
    wide["low_depth"] = wide["cons_norm_mq0"] < min_norm
    wide["ambiguous"] = wide["cons_mq_ratio"] < min_ratio
    consensus = wide["collapsed"] | wide["low_depth"]
    specific = []
    for p in pools:
        odd = ((norm0[p] > max_norm) | (norm0[p] < min_norm)) & ~consensus
        if odd.any():
            specific.append(pd.DataFrame({
                "chrom": wide.loc[odd, "chrom"], "start": wide.loc[odd, "start"],
                "end": wide.loc[odd, "end"], "pool": p,
                "norm_depth": norm0.loc[odd, p].round(3),
                "consensus_norm_depth": wide.loc[odd, "cons_norm_mq0"].round(3)}))
    specific = pd.concat(specific, ignore_index=True) if specific else pd.DataFrame(
        columns=["chrom", "start", "end", "pool", "norm_depth", "consensus_norm_depth"])
    return wide, specific


def cmd_mask(args):
    fai = read_fai(args.fai)
    keep = re.compile(args.chrom_regex)

    genmap = {}
    if args.genmap:
        log(f"reading {args.genmap}")
        gm = pd.read_csv(args.genmap, sep="\t", header=None, usecols=[0, 1, 2, 3],
                         names=["chrom", "start", "end", "value"],
                         dtype={"chrom": "category", "start": np.int64, "end": np.int64, "value": np.float32})
        # GenMap may keep the whole FASTA header; match on its first word
        gm["chrom"] = gm["chrom"].cat.rename_categories(lambda c: c.split()[0])
        genmap = {str(c): g for c, g in gm.groupby("chrom", sort=False, observed=True)}
        del gm

    bins = None
    if args.depth:
        log("reading mosdepth bins")
        bins, specific = consensus_depth_flags(args.depth, args.max_norm_depth, args.min_norm_depth,
                                               args.min_mq_ratio)
        specific.to_csv(f"{args.out}.pool_specific_depth.tsv", sep="\t", index=False)
        bins_by_chrom = {c: g for c, g in bins.groupby("chrom", sort=False)}

    seqs = iter_fasta(args.ref) if args.ref else ((c, None) for c in fai["chrom"])
    lengths = dict(zip(fai["chrom"], fai["length"]))

    summary = []
    with open(f"{args.out}.mask.bed", "w") as f_mask, \
            open(f"{args.out}.mask_reasons.bed", "w") as f_reasons, \
            open(f"{args.out}.accessible.bed", "w") as f_acc, \
            open(f"{args.out}.excluded.bed", "w") as f_excl:
        for chrom, seq in seqs:
            if chrom not in lengths:
                sys.exit(f"{chrom} is in the FASTA but not in {args.fai}")
            length = int(lengths[chrom])
            flags = np.zeros(length, np.uint8)
            if not keep.search(chrom):
                flags[:] = FLAGS["excluded"]
            else:
                if seq is not None:
                    s = np.frombuffer(seq, np.uint8)
                    if s.size != length:
                        sys.exit(f"{chrom}: FASTA length {s.size} differs from the .fai ({length})")
                    flags[(s == ord("N")) | (s == ord("n"))] |= FLAGS["gap"]
                    if args.softmask:
                        flags[(s >= ord("a")) & (s <= ord("z")) & (s != ord("n"))] |= FLAGS["softmask"]
                if args.genmap:
                    if chrom not in genmap:
                        log(f"warning: {chrom} not in the GenMap output; mappability not masked there")
                    else:
                        frac = unique_fraction(length, genmap.pop(chrom), args.kmer, args.unique_threshold)
                        flags[frac < args.min_unique_frac] |= FLAGS["low_mappability"]
                if bins is not None and chrom in bins_by_chrom:
                    g = bins_by_chrom[chrom]
                    for reason in ("collapsed", "ambiguous", "low_depth"):
                        for s0, e0 in g.loc[g[reason], ["start", "end"]].itertuples(index=False):
                            flags[s0:min(e0, length)] |= FLAGS[reason]

            row = {"chrom": chrom, "length": length}
            for reason, bit in FLAGS.items():
                row[reason] = int(np.count_nonzero(flags & bit))
            row["masked"] = int(np.count_nonzero(flags))
            row["accessible"] = length - row["masked"]
            summary.append(row)

            starts, ends, values = runs(flags)
            m = values != 0
            pd.DataFrame({"c": chrom, "s": starts[m], "e": ends[m],
                          "r": [",".join(r for r, b in FLAGS.items() if v & b) for v in values[m]]}
                         ).to_csv(f_reasons, sep="\t", header=False, index=False)
            masked = flags != 0
            for out_fh, select in ((f_mask, True), (f_acc, False)):
                s2, e2, v2 = runs(masked)
                k = v2 == select
                pd.DataFrame({"c": chrom, "s": s2[k], "e": e2[k]}).to_csv(out_fh, sep="\t", header=False,
                                                                           index=False)
            if flags[0] == FLAGS["excluded"]:
                f_excl.write(f"{chrom}\t0\t{length}\n")
            log(f"{chrom}: {row['masked'] / length:.1%} masked")

    summary = pd.DataFrame(summary)
    total = summary.drop(columns="chrom").sum()
    total["chrom"] = "total"
    summary = pd.concat([summary, total.to_frame().T], ignore_index=True)
    for reason in list(FLAGS) + ["masked", "accessible"]:
        summary[f"{reason}_frac"] = (summary[reason].astype(float) / summary["length"].astype(float)).round(4)
    summary.to_csv(f"{args.out}.mask_summary.tsv", sep="\t", index=False)
    write_params(f"{args.out}.mask.params", {"MASK_BED": f"{args.out}.mask.bed",
                                             "EXCLUDED_BED": f"{args.out}.excluded.bed",
                                             "ACCESSIBLE_BED": f"{args.out}.accessible.bed"})
    log(f"accessible: {int(total['accessible']):,} of {int(total['length']):,} bp")


# ---------------------------------------------------------------------------
# thresholds

def read_dist(path):
    """mosdepth region.dist: per chromosome, P(depth >= d) -> probability of each depth."""
    d = pd.read_csv(path, sep="\t", header=None, names=["chrom", "depth", "cum"],
                    dtype={"chrom": str, "depth": np.int64, "cum": np.float64})
    pmfs = {}
    for chrom, g in d.groupby("chrom", sort=False):
        if chrom == "total":
            continue
        cum = g.set_index("depth")["cum"].sort_index()
        cum = cum.reindex(np.arange(cum.index.max() + 2)).bfill().fillna(0.0).to_numpy()
        pmf = np.clip(cum[:-1] - cum[1:], 0, None)
        if pmf.sum() > 0:
            pmfs[chrom] = pmf / pmf.sum()
    return pmfs


def pmf_quantile(pmf, q):
    return int(np.searchsorted(np.cumsum(pmf), q))


def depth_bounds(median, min_fold, max_fold, abs_min):
    lower = max(abs_min, int(math.floor(min_fold * median)))
    upper = max(lower, int(math.ceil(max_fold * median)))
    return lower, upper


def binom_sf(c, n, p):
    """P(X >= c), X ~ Binomial(n, p)."""
    if c <= 0:
        return 1.0
    if c > n:
        return 0.0
    logp, log1p = math.log(p), math.log1p(-p)
    terms = [math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) + k * logp + (n - k) * log1p
             for k in range(c, n + 1)]
    m = max(terms)
    return math.exp(m) * sum(math.exp(t - m) for t in terms)


def error_site_prob(c, depth, err):
    """Chance that sequencing errors alone give some non-reference base >= c reads."""
    return 1 - (1 - binom_sf(c, depth, err / 3)) ** 3


def read_samtools_stats(path):
    err, mapq = np.nan, {}
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "SN" and f[1].startswith("error rate"):
                err = float(f[2])
            elif f[0] == "MAPQ":
                mapq[int(f[1])] = int(f[2])
    return err, mapq


def cmd_thresholds(args):
    xre = re.compile(args.x_regex)
    summary = pd.read_csv(args.mask_summary, sep="\t", dtype={"chrom": str})
    weights = dict(zip(summary["chrom"], summary["accessible"]))
    min_folds = sorted(set(args.grid_min_folds) | {args.min_fold})
    max_folds = sorted(set(args.grid_max_folds) | {args.max_fold})

    rows, classes, quantize = [], {}, []
    for name, path in args.dist:
        pmfs = read_dist(path)
        by_class = {}
        for chrom, pmf in pmfs.items():
            cls = "X" if xre.search(chrom) else "A"
            classes[chrom] = cls
            w = float(weights.get(chrom, 0))
            if w <= 0:
                continue
            acc = by_class.setdefault(cls, [np.zeros(1), 0.0])
            n = max(acc[0].size, pmf.size)
            acc[0] = np.pad(acc[0], (0, n - acc[0].size)) + np.pad(pmf * w, (0, n - pmf.size))
            acc[1] += w
        cuts = set()
        for cls, (pmf, w) in sorted(by_class.items()):
            pmf = pmf / pmf.sum()
            median = pmf_quantile(pmf, 0.5)
            lower, upper = depth_bounds(median, args.min_fold, args.max_fold, args.abs_min_depth)
            for mf in min_folds:
                cuts.add(depth_bounds(median, mf, args.max_fold, args.abs_min_depth)[0])
            for xf in max_folds:
                cuts.add(depth_bounds(median, args.min_fold, xf, args.abs_min_depth)[1] + 1)
            cdf = np.cumsum(pmf)
            rows.append({
                "pool": name, "class": cls, "accessible_bp": int(w), "median": median,
                "q01": pmf_quantile(pmf, 0.01), "q05": pmf_quantile(pmf, 0.05),
                "q95": pmf_quantile(pmf, 0.95), "q99": pmf_quantile(pmf, 0.99),
                "frac_zero": round(float(pmf[0]), 4),
                "lower": lower, "upper": upper,
                "frac_below": round(float(cdf[lower - 1]) if lower > 0 else 0.0, 4),
                "frac_above": round(max(float(1 - cdf[upper]), 0.0) if upper < cdf.size else 0.0, 4),
                "abs_min_depth": args.abs_min_depth})
        cuts.discard(0)
        quantize.append({"pool": name, "spec": "0:" + ":".join(str(c) for c in sorted(cuts)) + ":"})

    thr = pd.DataFrame(rows)
    thr["frac_kept"] = (1 - thr["frac_below"] - thr["frac_above"]).round(4)

    # Minimum allele count: the smallest count that errors alone rarely reach at the
    # pool's typical autosomal depth, with error rate set by the base-quality cut-off.
    err = args.error_rate if args.error_rate else 10 ** (-args.base_qual / 10)
    stats = dict((n, read_samtools_stats(p)) for n, p in (args.stats or []))
    mc_rows = []
    for name in thr["pool"].unique():
        a = thr[(thr["pool"] == name) & (thr["class"] == "A")]
        depth = int(a["median"].iloc[0]) if len(a) else int(thr.loc[thr["pool"] == name, "median"].max())
        row = {"pool": name, "median_depth": depth, "error_rate": err}
        chosen = None
        for c in range(1, 9):
            p = error_site_prob(c, depth, err)
            row[f"err_sites_per_Mb_c{c}"] = round(p * 1e6, 2)
            if chosen is None and p <= args.site_alpha:
                chosen = c
        row["min_count"] = chosen if chosen is not None else 8
        if name in stats:
            e, mapq = stats[name]
            total = sum(mapq.values())
            row["samtools_error_rate"] = e
            row["frac_reads_mapq_ge_cutoff"] = (round(sum(v for q, v in mapq.items() if q >= args.map_qual) / total, 4)
                                                if total else np.nan)
        mc_rows.append(row)
    mc = pd.DataFrame(mc_rows)

    # With the pool size: how many chromosomes each pool effectively samples at its depth
    # (reads are a sample of a sample), and the lowest allele frequency the minimum count
    # can see, against the frequency of a single copy in the pool.
    sizes = parse_pool_sizes(args.pool_sizes, list(mc["pool"]))
    if sizes:
        global_mc = int(mc["min_count"].max())
        n = mc["pool"].map(sizes)
        d = mc["median_depth"].astype(float)
        mc["pool_size"] = n
        mc["n_eff"] = (d * n / (d + n - 1)).round(1)
        mc["min_detectable_freq"] = (global_mc / d).round(4)
        mc["single_copy_freq"] = (1 / n).round(4)
        thr["pool_size"] = thr["pool"].map(sizes)
        thr["n_eff_at_upper"] = (thr["upper"] * thr["pool_size"] / (thr["upper"] + thr["pool_size"] - 1)).round(1)
        thr["n_eff_at_median"] = (thr["median"] * thr["pool_size"] / (thr["median"] + thr["pool_size"] - 1)).round(1)

    thr.to_csv(f"{args.out}.depth_thresholds.tsv", sep="\t", index=False)
    mc.to_csv(f"{args.out}.min_count.tsv", sep="\t", index=False)
    pd.DataFrame(quantize).to_csv(f"{args.out}.quantize.tsv", sep="\t", index=False, header=False)
    pd.DataFrame(sorted(classes.items()), columns=["chrom", "class"]).to_csv(
        f"{args.out}.chrom_class.tsv", sep="\t", index=False)
    write_params(f"{args.out}.thresholds.params", {
        "MIN_COUNT": int(mc["min_count"].max()),
        # Global grenedalf bounds are only a backstop: the depth masks apply the per-pool bounds
        "MIN_DEPTH": int(thr["lower"].min()), "MAX_DEPTH": int(thr["upper"].max()),
        "MIN_FOLD": args.min_fold, "MAX_FOLD": args.max_fold})
    print(thr.to_string(index=False))
    print()
    print(mc[["pool", "median_depth", "min_count"] + [c for c in mc.columns if c.startswith("err_sites")]]
          .to_string(index=False))


# ---------------------------------------------------------------------------
# depthmask

def parse_bound(x):
    return np.inf if x in ("inf", "") else float(x)


def cmd_depthmask(args):
    thr = pd.read_csv(args.thresholds, sep="\t")
    classes = dict(pd.read_csv(args.chrom_class, sep="\t", dtype={"chrom": str}).itertuples(index=False))
    fai = read_fai(args.fai)
    base = (merge_intervals(read_bed(args.base_mask)) if args.base_mask
            else pd.DataFrame({"chrom": pd.Series(dtype=str), "start": pd.Series(dtype=np.int64),
                               "end": pd.Series(dtype=np.int64)}))
    genome_bp = int(fai.loc[fai["chrom"].isin(classes), "length"].sum())
    accessible_bp = genome_bp - bp(base[base["chrom"].isin(classes)])

    per_pool, rows = {}, []
    for name, path in args.quantized:
        t = thr[thr["pool"] == name].set_index("class")
        if t.empty:
            sys.exit(f"{name} is not in {args.thresholds}")
        bounds = {cls: depth_bounds(r["median"], args.min_fold, args.max_fold, int(r["abs_min_depth"]))
                  for cls, r in t.iterrows()}
        q = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "label"],
                        dtype={"chrom": str, "start": np.int64, "end": np.int64, "label": str})
        q = q[q["chrom"].isin(classes)]
        lohi = q["label"].str.split(":", expand=True)
        lo = lohi[0].map(parse_bound)
        hi = lohi[1].map(parse_bound)
        cls = q["chrom"].map(classes)
        lower = cls.map(lambda c: bounds.get(c, (np.inf, -np.inf))[0])
        upper = cls.map(lambda c: bounds.get(c, (np.inf, -np.inf))[1])
        inside = (lo >= lower) & (hi - 1 <= upper)
        straddle = ~inside & (hi - 1 >= lower) & (lo <= upper)
        if straddle.any():
            log(f"warning: {name}: {int(straddle.sum())} quantize bins straddle a bound and are masked; "
                "rerun the quantize step with the cut points from the thresholds step")
        out = merge_intervals(q.loc[~inside, ["chrom", "start", "end"]])
        per_pool[name] = out
        lost = bp_outside(out, base)
        row = {"pool": name}
        for c, (lo_c, hi_c) in sorted(bounds.items()):
            row[f"lower_{c}"], row[f"upper_{c}"] = lo_c, hi_c
        rows.append({**row, "accessible_bp": accessible_bp, "out_of_range_bp": lost,
                     "out_of_range_frac": round(lost / accessible_bp, 4) if accessible_bp else np.nan})

    # grenedalf's own (global) depth filters are set just wide enough not to undo the
    # per-pool bounds, which the masks already apply
    lowers = [v for r in rows for k, v in r.items() if k.startswith("lower_")]
    uppers = [v for r in rows for k, v in r.items() if k.startswith("upper_")]
    params = {"MIN_FOLD": args.min_fold, "MAX_FOLD": args.max_fold, "MASK_MODE": args.mode,
              "MIN_DEPTH": min(lowers), "MAX_DEPTH": max(uppers)}
    if args.mode == "shared":
        total = merge_intervals(pd.concat([base] + list(per_pool.values()), ignore_index=True))
        write_bed(total, f"{args.out}.total_mask.bed")
        shared = genome_bp - bp(total[total["chrom"].isin(classes)])
        rows.append({"pool": "shared", "accessible_bp": accessible_bp, "out_of_range_bp": accessible_bp - shared,
                     "out_of_range_frac": round(1 - shared / accessible_bp, 4) if accessible_bp else np.nan})
        params.update({"TOTAL_MASK_BED": f"{args.out}.total_mask.bed", "SAMPLE_MASK_LIST": ""})
    else:
        write_bed(base, f"{args.out}.total_mask.bed")
        with open(f"{args.out}.sample_masks.tsv", "w") as fh:
            for name, out in per_pool.items():
                write_bed(out, f"{args.out}.{name}.depth_mask.bed")
                fh.write(f"{name}\t{args.out}.{name}.depth_mask.bed\n")
        params.update({"TOTAL_MASK_BED": f"{args.out}.total_mask.bed",
                       "SAMPLE_MASK_LIST": f"{args.out}.sample_masks.tsv"})
    summary = pd.DataFrame(rows)
    summary.to_csv(f"{args.out}.depthmask_summary.tsv", sep="\t", index=False)
    write_params(f"{args.out}.depthmask.params", params)
    print(summary.to_string(index=False))


# ---------------------------------------------------------------------------
# windows

def read_grenedalf(path):
    with open_text(path) as fh:
        sep = "\t" if "\t" in fh.readline() else ","
    t = pd.read_csv(path, sep=sep, dtype={"chrom": str})
    t["width"] = t["end"] - t["start"] + 1
    t["valid"] = t["total.invariant"] + t["total.passed"]
    t["callable"] = t["valid"] / t["width"]
    return t


def cmd_windows(args):
    t = read_grenedalf(args.table)
    if args.chrom_regex:
        t = t[t["chrom"].str.contains(args.chrom_regex, regex=True)]
    t = t[t["valid"] > 0]
    width = int(t["width"].median())
    rate = (t["total.passed"] / t["width"]).to_numpy()   # SNPs per bp of window, callable or not
    density = t["total.passed"].sum() / t["valid"].sum()

    rows = []
    for w in range(args.round, args.max_width + 1, args.round):
        snps = rate * w
        rows.append({"width": w, "median_snps": round(float(np.median(snps)), 1),
                     "frac_windows_ge_target": round(float(np.mean(snps >= args.target_snps)), 4)})
    grid = pd.DataFrame(rows)
    ok = grid[grid["frac_windows_ge_target"] >= args.target_fraction]
    w_snps = int(ok["width"].iloc[0]) if len(ok) else int(grid["width"].iloc[-1])
    grid.to_csv(f"{args.out}.window_widths.tsv", sep="\t", index=False)

    # Linkage scale from the ld step: SNPs closer than LD_D10 still carry correlated
    # frequencies, so windows narrower than that are not independent of their neighbours.
    ld = {}
    if args.ld:
        with open(args.ld) as fh:
            ld = dict(line.strip().split("=", 1) for line in fh if "=" in line)
    d10 = float(ld["LD_D10"]) if ld.get("LD_D10", "NA") not in ("NA", "nan") else None
    rec = w_snps
    if d10:
        rec = max(w_snps, int(math.ceil(d10 / args.round) * args.round))

    lines = [f"pilot table: {args.table} ({len(t):,} windows of {width:,} bp)",
             f"SNPs per callable bp: {density:.4f}",
             f"callable fraction per window: median {t['callable'].median():.3f}, "
             f"10% of windows below {t['callable'].quantile(0.1):.3f}",
             f"width for >= {args.target_snps} SNPs in >= {args.target_fraction:.0%} of windows: {w_snps:,} bp"]
    if ld:
        lines.append(f"linkage scale (frequency correlation across pools): half decay {ld.get('LD_D50')} bp, "
                     f"90% decay {ld.get('LD_D10')} bp")
    if d10:
        lines.append(f"recommended width: {rec:,} bp, stride {rec // 2:,} bp "
                     f"({'set by linkage scale' if rec > w_snps else 'set by SNP density'}); "
                     f"about {rec / d10:.1f} independent segments per window")
    else:
        lines.append(f"recommended width: {rec:,} bp, stride {rec // 2:,} bp (SNP density only)")

    if args.mask_bed:
        # Did grenedalf mask what the BED says? total.masked should equal the BED overlap.
        mask = merge_intervals(read_bed(args.mask_bed))
        exp_total, obs_total = 0, int(t["total.masked"].sum())
        for chrom, g in t.groupby("chrom"):
            m = mask[mask["chrom"] == chrom]
            cs = coverage_prefix(m, int(g["end"].max()))
            exp_total += int((cs[g["end"].to_numpy()] - cs[g["start"].to_numpy() - 1]).sum())
        win_total = int(t["width"].sum())
        if exp_total == 0:
            verdict = "mask has no overlap with these windows"
        elif abs(obs_total - exp_total) <= 0.02 * exp_total:
            verdict = "OK: grenedalf masked the BED intervals"
        elif abs(obs_total - (win_total - exp_total)) <= 0.02 * (win_total - exp_total):
            verdict = "WRONG WAY ROUND: grenedalf masked everything except the BED (check --filter-mask-total-invert)"
        elif obs_total == 0:
            verdict = "NOT APPLIED: total.masked is 0"
        else:
            verdict = "MISMATCH: compare total.masked with the BED by hand"
        lines.append(f"mask check: expected {exp_total:,} masked bp, grenedalf reported {obs_total:,} -> {verdict}")

    write_params(f"{args.out}.windows.params", {"WINDOW": rec, "STRIDE": rec // 2, "MIN_SNPS": args.target_snps,
                                                "WINDOW_SNPS": w_snps})
    with open(f"{args.out}.windows.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


# ---------------------------------------------------------------------------
# chromclass

def block_bootstrap_ratio(x, a, n_boot, block, rng):
    """95% interval of median(x) / median(a), resampling blocks of adjacent bins."""
    def blocks(v):
        k = max(len(v) // block, 1)
        return [v[i * block:(i + 1) * block] for i in range(k)] if len(v) >= block else [v]
    bx, ba = blocks(x), blocks(a)
    out = np.empty(n_boot)
    for i in range(n_boot):
        sx = np.concatenate([bx[j] for j in rng.integers(0, len(bx), len(bx))])
        sa = np.concatenate([ba[j] for j in rng.integers(0, len(ba), len(ba))])
        out[i] = np.median(sx) / np.median(sa)
    return np.percentile(out, [2.5, 97.5])


def cmd_chromclass(args):
    """Do the X and the mitochondrion behave like autosomes? Depth ratios to autosomes per pool.

    X: in a pool with males (hemizygous), X depth is below autosome depth, and the X is
    sampled from fewer copies, so grenedalf needs a smaller X pool size: with a fraction m
    of males, X depth / autosome depth = 1 - m/2, and the X pool size is pool_size * (1 - m/2).
    A 50:50 pool gives 0.75. m is estimated from the depth ratio unless given with
    --male-fraction (known from the design), in which case the known value sets the X pool
    size and the estimate is reported as a check (consistent if its interval contains 1 - m/2).
    If any pool's ratio differs from 1 (95% interval excludes it, by more than
    --min-difference), the X is run separately with its own pool sizes.
    Mitochondrion: haploid, many copies per cell, one lineage per individual: always run
    separately, with its own depth bounds and a pool size of the number of individuals."""
    xre, mre, keep = re.compile(args.x_regex), re.compile(args.mito_regex), re.compile(args.chrom_regex)
    fai = read_fai(args.fai)
    mask = merge_intervals(read_bed(args.mask_bed)) if args.mask_bed else None
    sizes = parse_pool_sizes(args.pool_sizes, [n for n, _ in args.depth])
    known_m = parse_pool_sizes(args.male_fraction, [n for n, _ in args.depth])
    rng = np.random.default_rng(args.seed)

    def chrom_class(c):
        if mre.search(c):
            return "M"
        if not keep.search(c):
            return None
        return "X" if xre.search(c) else "A"

    rows, mito_bounds = [], {}
    for name, path in args.depth:
        b = read_mosdepth_bins(path, "depth")
        b["class"] = b["chrom"].map(chrom_class)
        b = b[b["class"].notna()].copy()
        if mask is not None:
            # A and X: bins at least half accessible; the mitochondrion is in the mask as excluded
            frac = np.zeros(len(b))
            for chrom, idx in b.groupby("chrom").indices.items():
                m = mask[mask["chrom"] == chrom]
                if m.empty:
                    continue
                g = b.iloc[idx]
                cs = coverage_prefix(m, int(g["end"].max()))
                frac[idx] = (cs[g["end"].to_numpy()] - cs[g["start"].to_numpy()]) / (g["end"] - g["start"]).to_numpy()
            b = b[(frac < 0.5) | (b["class"] == "M")]
        a = b.loc[b["class"] == "A", "depth"].to_numpy()
        if a.size == 0:
            sys.exit(f"{name}: no autosomal bins")
        med_a = float(np.median(a))
        n = sizes.get(name, np.nan)
        rows.append({"pool": name, "class": "A", "bins": a.size, "median_depth": med_a, "ratio_to_A": 1.0,
                     "ci_low": 1.0, "ci_high": 1.0, "differs": False, "pool_size": n})
        for cls in ("X", "M"):
            x = b.loc[b["class"] == cls, "depth"].to_numpy()
            if x.size == 0:
                continue
            ratio = float(np.median(x)) / med_a
            lo, hi = block_bootstrap_ratio(x, a, args.bootstrap, args.block, rng)
            differs = bool((hi < 1 or lo > 1) and abs(ratio - 1) > args.min_difference)
            row = {"pool": name, "class": cls, "bins": x.size, "median_depth": float(np.median(x)),
                   "ratio_to_A": round(ratio, 4), "ci_low": round(float(lo), 4), "ci_high": round(float(hi), 4),
                   "differs": differs}
            if cls == "X":
                m_est = float(np.clip(2 * (1 - ratio), 0, 1))
                row["male_fraction_est"] = round(m_est, 3)
                row["male_fraction_ci"] = f"{np.clip(2 * (1 - hi), 0, 1):.2f}-{np.clip(2 * (1 - lo), 0, 1):.2f}"
                if name in known_m:
                    m = known_m[name]
                    row["male_fraction_known"] = m
                    row["expected_ratio"] = round(1 - m / 2, 4)
                    row["matches_known"] = bool(lo <= 1 - m / 2 <= hi)
                else:
                    m = m_est
                row["male_fraction_used"] = round(m, 3)
                row["pool_size"] = round(n * (1 - m / 2)) if np.isfinite(n) else np.nan
            else:
                row["pool_size"] = round(n / 2) if np.isfinite(n) else np.nan
                mito_bounds[name] = depth_bounds(float(np.median(x)), args.min_fold, args.max_fold, args.abs_min_depth)
            rows.append(row)

    t = pd.DataFrame(rows)
    t.to_csv(f"{args.out}.chrom_classes.tsv", sep="\t", index=False)
    x = t[t["class"] == "X"]
    m = t[t["class"] == "M"]
    params = {"X_SEPARATE": int(bool(x["differs"].any())) if len(x) else 0,
              "MITO_CHROMS": ",".join(c for c in fai["chrom"] if mre.search(c))}
    if len(x) and x["pool_size"].notna().all():
        with open(f"{args.out}.x_pool_sizes.tsv", "w") as fh:
            for r in x.itertuples():
                fh.write(f"{r.pool}\t{int(r.pool_size)}\n")
        params["X_POOL_SIZES"] = f"{args.out}.x_pool_sizes.tsv"
    if len(m):
        params.update({"MITO_SEPARATE": 1, "MITO_MIN_DEPTH": min(v[0] for v in mito_bounds.values()),
                       "MITO_MAX_DEPTH": max(v[1] for v in mito_bounds.values())})
        if m["pool_size"].notna().all():
            with open(f"{args.out}.mito_pool_sizes.tsv", "w") as fh:
                for r in m.itertuples():
                    fh.write(f"{r.pool}\t{int(r.pool_size)}\n")
            params["MITO_POOL_SIZES"] = f"{args.out}.mito_pool_sizes.tsv"
    else:
        params["MITO_SEPARATE"] = 0
    write_params(f"{args.out}.chromclass.params", params)
    print(t.to_string(index=False))
    if len(x):
        print("\nX differs from autosomes in at least one pool: run it separately with its own pool sizes"
              if params["X_SEPARATE"] else "\nX behaves like the autosomes in every pool: analyse together")


# ---------------------------------------------------------------------------
# commands

def cmd_commands(args):
    """Ready-to-run grenedalf commands with the chosen settings, one run per chromosome class."""
    P = args.prefix
    cc = {}
    if os.path.exists(f"{P}.chromclass.params"):
        with open(f"{P}.chromclass.params") as fh:
            cc = dict(l.strip().split("=", 1) for l in fh if "=" in l)
    fai = read_fai(args.fai)
    xre = re.compile(args.x_regex)
    mito = [c for c in cc.get("MITO_CHROMS", "").split(",") if c]
    x_chroms = [c for c in fai["chrom"] if xre.search(c) and c not in mito]
    x_sep = cc.get("X_SEPARATE", "0") == "1" and x_chroms
    base = merge_intervals(read_bed(args.total_mask))
    whole = lambda chroms: fai[fai["chrom"].isin(chroms)].assign(start=0).rename(columns={"length": "end"})[
        ["chrom", "start", "end"]]

    runs = []
    if x_sep:
        # Autosome run: X added to the mask. X run: everything but X added to the mask.
        auto_mask = f"{P}.mask_autosomes.bed"
        write_bed(merge_intervals(pd.concat([base, whole(x_chroms)])), auto_mask)
        x_mask = f"{P}.mask_X.bed"
        write_bed(merge_intervals(pd.concat([base, whole([c for c in fai["chrom"] if c not in x_chroms])])), x_mask)
        runs.append(("autosomes", "autosomes (X masked; it gets its own pool sizes below)", auto_mask,
                     args.pool_size_arg, args.min_depth, args.max_depth, args.sample_masks))
        runs.append(("X", "X chromosome, with per-pool X pool sizes (fewer copies in pools with males)", x_mask,
                     cc.get("X_POOL_SIZES", args.pool_size_arg), args.min_depth, args.max_depth, args.sample_masks))
    else:
        runs.append(("nuclear", "autosomes" + (" and X (X depth matches the autosomes)" if x_chroms else ""),
                     args.total_mask, args.pool_size_arg, args.min_depth, args.max_depth, args.sample_masks))
    if mito and cc.get("MITO_SEPARATE", "0") == "1":
        m_mask = f"{P}.mask_mito.bed"
        write_bed(whole([c for c in fai["chrom"] if c not in mito]), m_mask)
        runs.append(("mito", "mitochondrion: haploid, pool size = individuals, its own depth bounds "
                     "(whole genome as one window)", m_mask, cc.get("MITO_POOL_SIZES", args.pool_size_arg),
                     cc.get("MITO_MIN_DEPTH", args.min_depth), cc.get("MITO_MAX_DEPTH", args.max_depth), ""))

    L = ["#!/usr/bin/env bash",
         f"# grenedalf commands recommended by grenedalf-optimiser ({pd.Timestamp.now():%Y-%m-%d}).",
         "# Settings chosen on technical criteria only; see optimiser_report.html for why.",
         "set -euo pipefail", "",
         f"OUTDIR={args.grenedalf_out}", f"THREADS={args.threads}", ""]
    for key, desc, mask_bed, pools, mind, maxd, smasks in runs:
        mito_run = key == "mito"
        L.append(f"# --- {desc}")
        for cmd in ("diversity", "fst"):
            c = [f"grenedalf {cmd}"]
            c += [f"--sam-path {b}" for b in args.bam]
            c += [f"--reference-genome-fai {args.fai}",
                  f"--sam-min-map-qual {args.map_qual} --sam-min-base-qual {args.base_qual}",
                  f"--filter-mask-total-bed {mask_bed}"]
            if smasks:
                c.append(f"--filter-mask-sample-bed-list {smasks}")
            c += [f"--filter-sample-min-count {args.min_count}",
                  f"--filter-sample-min-read-depth {mind} --filter-sample-max-read-depth {maxd}"]
            if mito_run:
                c.append("--window-type chromosomes")
            else:
                c.append(f"--window-type interval --window-interval-width {args.window} "
                         f"--window-interval-stride {args.stride}")
            c += ["--window-average-policy valid-loci", f"--pool-sizes {pools}"]
            if cmd == "fst":
                c.append(f"--method {args.fst_method}")
            c += ["--separator-char tab --threads $THREADS",
                  f"--out-dir $OUTDIR/{key} --file-prefix {key}_"]
            L.append(" \\\n    ".join(c))
            L.append("")
    if args.min_snps:
        L.append(f"# Downstream: drop windows with fewer than {args.min_snps} passed SNPs (total.passed).")
    out = os.path.join(args.out_dir, "grenedalf_commands.sh")
    with open(out, "w") as fh:
        fh.write("\n".join(L) + "\n")
    os.chmod(out, 0o755)
    print(f"wrote {out} ({len(runs)} run(s): {', '.join(r[0] for r in runs)})")


# ---------------------------------------------------------------------------
# ld

def read_frequency(path, region, min_count, min_depth, max_depth, mask_bed, chunksize=2_000_000):
    """SNPs from a grenedalf frequency table (--write-sample-counts), filtered as the
    diversity/fst runs are: masked positions dropped, counts below min_count set to zero,
    every pool within the depth bounds, still polymorphic. Returns chrom, pos, ALT freq."""
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n")
    sep = "\t" if "\t" in header else ","
    cols = header.split(sep)
    pat = re.compile(r"^(.+)\.(REF_CNT|ALT_CNT)$")
    ref, alt = {}, {}
    for c in cols:
        m = pat.match(c)
        if m and m.group(1) != "TOTAL":
            (ref if m.group(2) == "REF_CNT" else alt)[m.group(1)] = c
    samples = [x for x in ref if x in alt]
    if not samples:
        sys.exit(f"no <sample>.REF_CNT / <sample>.ALT_CNT columns in {path}")
    chrom_c = next(c for c in cols if c.lower() in ("chrom", "chromosome"))
    pos_c = next(c for c in cols if c.lower() in ("pos", "position"))
    rc, ac = [ref[x] for x in samples], [alt[x] for x in samples]
    mask = merge_intervals(read_bed(mask_bed)) if mask_bed else None

    out_c, out_p, out_f = [], [], []
    for chunk in pd.read_csv(path, sep=sep, usecols=[chrom_c, pos_c] + rc + ac, chunksize=chunksize,
                             dtype={chrom_c: str, pos_c: np.int64, **{c: np.int32 for c in rc + ac}}):
        if region:
            chunk = chunk[chunk[chrom_c].isin(region)]
        if chunk.empty:
            continue
        r, a = chunk[rc].to_numpy(), chunk[ac].to_numpy()
        r[r < min_count] = 0
        a[a < min_count] = 0
        d = r + a
        keep = ((d >= min_depth) & (d <= max_depth)).all(axis=1) & (r.sum(axis=1) > 0) & (a.sum(axis=1) > 0)
        ch, pos = chunk[chrom_c].to_numpy(), chunk[pos_c].to_numpy()
        if mask is not None:
            for c in np.unique(ch[keep]):
                m = mask[mask["chrom"] == c]
                if m.empty:
                    continue
                idx = np.flatnonzero(keep & (ch == c))
                p0 = pos[idx] - 1  # grenedalf positions are 1-based
                i = np.searchsorted(m["start"].to_numpy(), p0, side="right") - 1
                inside = (i >= 0) & (p0 < m["end"].to_numpy()[np.clip(i, 0, None)])
                keep[idx[inside]] = False
        out_c.append(ch[keep])
        out_p.append(pos[keep])
        out_f.append(a[keep] / d[keep])
    if not out_c:
        sys.exit("no SNPs left after filtering")
    return np.concatenate(out_c), np.concatenate(out_p), np.vstack(out_f), samples


def pair_r2(za, zb):
    """Squared correlation across pools of standardised frequency vectors."""
    return (np.einsum("ij,ij->i", za, zb) / za.shape[1]) ** 2


cmd_ld_doc = """
Pool-seq reads are not phased, but SNPs in linkage disequilibrium share their frequency
differences between pools: if pool A carries more of a haplotype than pool B, every SNP
on that haplotype differs between A and B the same way. So the correlation, across
pools, of two SNPs' centred frequencies is high when they are close and falls to a
baseline (set by the number of pools and by sampling noise) for unlinked SNPs. The
distance at which the excess over that baseline has decayed is the scale beyond which
SNPs carry independent information. It depends on recombination, LD and the pools'
sampling, not on any biological contrast between pools.

Outputs <out>.ld_decay.tsv (mean r2 per distance bin) and <out>.ld.params (LD_D50,
LD_D10: distances at which the excess r2 has fallen to 50% and 10% of its value in the
shortest bin)."""


def cmd_ld(args):
    chrom, pos, f, samples = read_frequency(args.frequency, set(args.region) if args.region else None,
                                            args.min_count, args.min_depth, args.max_depth, args.mask_bed)
    k = len(samples)
    if k < 3:
        sys.exit("need at least 3 pools for a correlation across pools")
    mean = f.mean(axis=1)
    keep = (np.minimum(mean, 1 - mean) >= args.min_maf) & (f.std(axis=1) > 0)
    chrom, pos, f = chrom[keep], pos[keep], f[keep]
    z = (f - f.mean(axis=1, keepdims=True)) / f.std(axis=1, keepdims=True)
    order = np.lexsort((pos, chrom))
    chrom, pos, z = chrom[order], pos[order], z[order]
    log(f"{len(pos):,} SNPs from {k} pools after filters")

    rng = np.random.default_rng(args.seed)
    edges = np.unique(np.round(np.logspace(np.log10(args.min_dist), np.log10(args.max_dist), args.bins + 1)))
    nb = edges.size - 1
    sums, sq, n = np.zeros(nb), np.zeros(nb), np.zeros(nb, np.int64)
    base_r2 = []
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        p, zc = pos[idx], z[idx]
        if p.size < 2:
            continue
        # For each anchor SNP and distance bin, one random partner SNP in that bin
        anchors = rng.choice(p.size, size=min(args.anchors, p.size), replace=False)
        for b in range(nb):
            lo = np.searchsorted(p, p[anchors] + edges[b], side="left")
            hi = np.searchsorted(p, p[anchors] + edges[b + 1], side="left")
            ok = hi > lo
            if not ok.any():
                continue
            a = anchors[ok]
            j = lo[ok] + (rng.random(ok.sum()) * (hi[ok] - lo[ok])).astype(np.int64)
            r2 = pair_r2(zc[a], zc[j])
            sums[b] += r2.sum()
            sq[b] += (r2 ** 2).sum()
            n[b] += r2.size
        # Baseline: pairs on the same chromosome far beyond the decay range
        a = rng.choice(p.size, size=args.anchors)
        b2 = rng.choice(p.size, size=args.anchors)
        sel = np.abs(p[a] - p[b2]) > 5 * args.max_dist
        if sel.any():
            base_r2.append(pair_r2(zc[a[sel]], zc[b2[sel]]))
    if np.unique(chrom).size > 1:
        # and pairs on different chromosomes
        a = rng.choice(len(pos), size=args.anchors)
        b2 = rng.choice(len(pos), size=args.anchors)
        sel = chrom[a] != chrom[b2]
        base_r2.append(pair_r2(z[a[sel]], z[b2[sel]]))
    base_r2 = np.concatenate(base_r2) if base_r2 else np.array([])
    if base_r2.size >= 100:
        baseline, baseline_src = float(base_r2.mean()), f"{base_r2.size:,} unlinked pairs"
    else:
        baseline, baseline_src = 1 / (k - 1), "expectation 1/(pools - 1); too few unlinked pairs"

    with np.errstate(invalid="ignore", divide="ignore"):
        mean_r2 = sums / n
        se = np.sqrt((sq / n - mean_r2 ** 2) / n)
    mid = np.sqrt(edges[:-1] * edges[1:])
    good = n >= args.min_pairs
    excess = mean_r2 - baseline
    e0 = excess[np.flatnonzero(good)[0]] if good.any() else np.nan
    rel = np.where(good, excess / e0, np.nan) if e0 > 0 else np.full(nb, np.nan)

    def decay_to(target):
        below = np.flatnonzero(good & (rel <= target))
        if not below.size:
            return np.nan
        i = below[0]
        prev = np.flatnonzero(good[:i])
        if not prev.size:
            return mid[i]
        prev = prev[-1]  # interpolate on log distance between the bins either side
        x0, x1, y0, y1 = np.log10(mid[prev]), np.log10(mid[i]), rel[prev], rel[i]
        return 10 ** (x0 + (y0 - target) * (x1 - x0) / (y0 - y1)) if y0 != y1 else mid[i]

    d50, d10 = decay_to(0.5), decay_to(0.1)
    tab = pd.DataFrame({"dist_from": edges[:-1].astype(int), "dist_to": edges[1:].astype(int),
                        "dist_mid": mid.round().astype(int), "pairs": n, "mean_r2": mean_r2.round(5),
                        "se": se.round(5), "excess_rel": np.round(rel, 4)})
    tab.to_csv(f"{args.out}.ld_decay.tsv", sep="\t", index=False)
    fmt = lambda v: "NA" if not np.isfinite(v) else str(int(round(v)))
    write_params(f"{args.out}.ld.params", {"LD_D50": fmt(d50), "LD_D10": fmt(d10),
                                           "LD_BASELINE_R2": round(baseline, 5), "LD_SNPS": len(pos)})
    print(tab.to_string(index=False))
    print(f"\nbaseline r2 {baseline:.4f} ({baseline_src})")
    print(f"excess r2 over baseline halves by {fmt(d50)} bp and falls to 10% by {fmt(d10)} bp")
    if not np.isfinite(e0) or e0 <= 0:
        print("no excess correlation at short distances: no linkage scale estimated")
    elif not np.isfinite(d10):
        print(f"excess r2 had not fallen to 10% by {args.max_dist:,} bp; rerun with a larger --max-dist")


# ---------------------------------------------------------------------------
# assess

def window_depth(windows, bins, col):
    """Mean of the bins whose midpoint lies in each window (grenedalf windows are 1-based, closed)."""
    out = np.full(len(windows), np.nan)
    for chrom, idx in windows.groupby("chrom").indices.items():
        b = bins[bins["chrom"] == chrom]
        if b.empty:
            continue
        mid = ((b["start"] + b["end"]) // 2).to_numpy()
        order = np.argsort(mid)
        mid, val = mid[order], b[col].to_numpy()[order]
        cs = np.concatenate(([0.0], np.cumsum(val)))
        w = windows.iloc[idx]
        lo = np.searchsorted(mid, w["start"].to_numpy() - 1, side="left")
        hi = np.searchsorted(mid, w["end"].to_numpy(), side="left")
        n = hi - lo
        with np.errstate(invalid="ignore", divide="ignore"):
            out[idx] = np.where(n > 0, (cs[hi] - cs[lo]) / n, np.nan)
    return out


def spearman(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 10:
        return np.nan
    return float(pd.Series(a[m]).rank().corr(pd.Series(b[m]).rank()))


def cmd_assess(args):
    depth = {}
    for name, path in args.depth or []:
        b = read_mosdepth_bins(path, "depth")
        b["norm"] = normalise_by_chrom(b, "depth")
        depth[name] = b

    run_rows, stat_rows = [], []
    for run, *tables in args.run:
        frames = [read_grenedalf(p) for p in tables if p != "-"]
        t = frames[0]
        for f in frames[1:]:
            stat_cols = [c for c in f.columns if c.endswith(STAT_SUFFIXES) or re.search(r":.+\.(passed|missing)$", c)]
            t = t.merge(f[["chrom", "start", "end"] + stat_cols], on=["chrom", "start", "end"], how="outer")
        # Windows overlap when stride < width; divide sums by width / stride to count each base once
        strides = t.groupby("chrom")["start"].apply(lambda x: np.diff(np.sort(x.to_numpy()))).explode()
        stride = float(strides.astype(float).median()) if strides.notna().any() else float(t["width"].median())
        step = max(float(t["width"].median()) / max(stride, 1.0), 1.0)
        windows_all = len(t)
        callable_bp = int(t["valid"].sum() / step)
        snps = int(t["total.passed"].sum() / step)
        # Statistics are only compared in windows with enough SNPs to be estimated precisely
        t = t[t["total.passed"] >= args.min_snps].reset_index(drop=True)
        wdepth = {p: window_depth(t, b, "norm") for p, b in depth.items()}
        stats = [c for c in t.columns if c.endswith(STAT_SUFFIXES)]

        for col in stats:
            x = t[col].to_numpy(dtype=float)
            row = {"run": run, "stat": col, "type": col.rsplit(".", 1)[1], "n_windows": int(np.isfinite(x).sum()),
                   "median": float(np.nanmedian(x)) if np.isfinite(x).any() else np.nan,
                   "rho_callable": spearman(x, t["callable"].to_numpy())}
            sample = col.rsplit(".", 1)[0]
            pools = sample.split(":")
            if all(p in wdepth for p in pools):
                d = np.vstack([wdepth[p] for p in pools])
                row["rho_depth"] = spearman(x, d.mean(axis=0))
                if len(pools) == 2:
                    with np.errstate(invalid="ignore", divide="ignore"):
                        row["rho_depth_imbalance"] = spearman(x, np.abs(np.log2(d[0] / d[1])))
            stat_rows.append(row)

        s = pd.DataFrame([r for r in stat_rows if r["run"] == run])
        pair_passed = [c for c in t.columns if re.search(r":.+\.passed$", c)]
        imbalance = np.nan
        if pair_passed:
            pp = t[pair_passed].to_numpy(dtype=float)
            with np.errstate(invalid="ignore", divide="ignore"):
                imbalance = float(np.nanmedian((pp.max(axis=1) - pp.min(axis=1)) / pp.mean(axis=1)))
        # Artefact score: for each statistic type (theta_pi, fst, ...) and technical
        # covariate, the median |rho| over pools or pairs; the score is the worst of these.
        row = {"run": run, "windows": windows_all, "windows_scored": len(t),
               "callable_bp": callable_bp, "snps": snps,
               "median_callable_frac": round(float(t["callable"].median()), 4) if len(t) else np.nan,
               "pair_callable_imbalance": round(imbalance, 4) if np.isfinite(imbalance) else np.nan}
        worst, worst_name = 0.0, ""
        for cov in ("rho_depth", "rho_callable", "rho_depth_imbalance"):
            if cov not in s:
                continue
            for typ, g in s.groupby("type"):
                v = g[cov].abs().dropna()
                if len(v):
                    v = v.median()
                    row[f"{typ}.{cov}"] = round(float(v), 4)
                    if v > worst:
                        worst, worst_name = float(v), f"{typ} vs {cov[4:]}"
        row["artefact"], row["worst"] = worst, worst_name
        run_rows.append(row)

    runs_df = pd.DataFrame(run_rows)
    stats_df = pd.DataFrame(stat_rows)

    # Some dependence on depth can be real (divergent haplotypes map less well), so each
    # run is judged against the strictest run rather than against zero.
    # "Permissive" = more callable bases, then more SNPs (min-count changes only the latter)
    by_data = runs_df.sort_values(["callable_bp", "snps"], kind="stable")
    strict = args.strict if args.strict else by_data["run"].iloc[0]
    ref = float(runs_df.loc[runs_df["run"] == strict, "artefact"].iloc[0])
    runs_df["artefact_vs_strict"] = (runs_df["artefact"] - ref).round(4)
    runs_df["acceptable"] = runs_df["artefact_vs_strict"] <= args.tolerance
    ok = runs_df[runs_df["acceptable"]]
    key = list(zip(ok["callable_bp"], ok["snps"]))
    best = ok["run"].iloc[key.index(max(key))] if key else strict  # ties: first run listed
    runs_df["recommended"] = runs_df["run"] == best
    for c in runs_df.columns:
        if runs_df[c].dtype == float:
            runs_df[c] = runs_df[c].round(4)

    runs_df.to_csv(f"{args.out}.assess_runs.tsv", sep="\t", index=False)
    stats_df.round(4).to_csv(f"{args.out}.assess_stats.tsv", sep="\t", index=False)
    print(runs_df.to_string(index=False))
    print(f"\nreference (strictest) run: {strict}; tolerance {args.tolerance}")
    print(f"most permissive run that adds no technical artefact: {best}")
    with open(f"{args.out}.assess.params", "w") as fh:
        fh.write(f"RECOMMENDED_RUN={best}\nTOLERANCE={args.tolerance}\n")


# ---------------------------------------------------------------------------

def floats(s):
    return [float(x) for x in s.split(",") if x]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("mask", help="accessibility mask shared by all pools")
    p.add_argument("--fai", required=True)
    p.add_argument("--ref", help="reference FASTA (for N gaps and --softmask)")
    p.add_argument("--softmask", action="store_true", help="also mask lower-case (repeat-masked) bases")
    p.add_argument("--genmap", help="GenMap mappability bedGraph (genmap map -bg)")
    p.add_argument("--kmer", type=int, default=150, help="k used for GenMap (read length) [150]")
    p.add_argument("--unique-threshold", type=float, default=1.0,
                   help="GenMap value counted as unique (1 = one occurrence) [1.0]")
    p.add_argument("--min-unique-frac", type=float, default=0.5,
                   help="mask a base when fewer than this fraction of the k-mers overlapping it are unique [0.5]")
    p.add_argument("--depth", nargs=3, action="append", metavar=("POOL", "MQ0_BINS", "MQX_BINS"),
                   help="mosdepth --by regions.bed.gz without and with the MAPQ cut-off; repeat per pool")
    p.add_argument("--max-norm-depth", type=float, default=2.0,
                   help="consensus depth / chromosome median above which a bin is a collapsed repeat [2.0]")
    p.add_argument("--min-norm-depth", type=float, default=0.25,
                   help="consensus depth / chromosome median below which a bin is masked as low depth [0.25]")
    p.add_argument("--min-mq-ratio", type=float, default=0.5,
                   help="consensus depth at MAPQ cut-off / depth at MAPQ 0 below which a bin is ambiguous [0.5]")
    p.add_argument("--chrom-regex", default=".", help="chromosomes to analyse; others are masked [all]")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_mask)

    p = sub.add_parser("thresholds", help="per-pool depth bounds and minimum allele count")
    p.add_argument("--dist", nargs=2, action="append", required=True, metavar=("POOL", "REGION_DIST"),
                   help="mosdepth .mosdepth.region.dist.txt from --by accessible.bed; repeat per pool")
    p.add_argument("--mask-summary", required=True, help="<out>.mask_summary.tsv from the mask step")
    p.add_argument("--stats", nargs=2, action="append", metavar=("POOL", "SAMTOOLS_STATS"))
    p.add_argument("--x-regex", default=r"chrX|_X_|^X$", help="X chromosome name pattern")
    p.add_argument("--min-fold", type=float, default=0.5, help="lower depth bound = fold x median [0.5]")
    p.add_argument("--max-fold", type=float, default=2.0, help="upper depth bound = fold x median [2.0]")
    p.add_argument("--abs-min-depth", type=int, default=10, help="lower bound never below this [10]")
    p.add_argument("--grid-min-folds", type=floats, default=[0.33, 0.67], help="lower folds tried in the grid")
    p.add_argument("--grid-max-folds", type=floats, default=[1.5, 3.0], help="upper folds tried in the grid")
    p.add_argument("--base-qual", type=int, default=30, help="grenedalf --sam-min-base-qual [30]")
    p.add_argument("--map-qual", type=int, default=30, help="grenedalf --sam-min-map-qual [30]")
    p.add_argument("--error-rate", type=float, help="per-base error rate [10^(-base-qual/10)]")
    p.add_argument("--pool-sizes", help="haploid genomes per pool: one number, or POOL=N,POOL=N")
    p.add_argument("--site-alpha", type=float, default=1e-4,
                   help="accepted chance per site that errors alone pass the minimum count [1e-4]")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_thresholds)

    p = sub.add_parser("depthmask", help="per-pool depth masks from mosdepth --quantize")
    p.add_argument("--quantized", nargs=2, action="append", required=True, metavar=("POOL", "QUANTIZED_BED"))
    p.add_argument("--thresholds", required=True)
    p.add_argument("--chrom-class", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--base-mask", help="accessibility mask to combine with")
    p.add_argument("--min-fold", type=float, default=0.5)
    p.add_argument("--max-fold", type=float, default=2.0)
    p.add_argument("--mode", choices=["shared", "per-sample"], default="shared",
                   help="shared: one mask, every pool called at the same sites; "
                        "per-sample: one mask per pool (more sites, uneven missingness) [shared]")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_depthmask)

    p = sub.add_parser("windows", help="window width from SNP density in a pilot grenedalf table")
    p.add_argument("--table", required=True, help="grenedalf diversity or fst table")
    p.add_argument("--target-snps", type=int, default=50, help="SNPs wanted per window [50]")
    p.add_argument("--target-fraction", type=float, default=0.9,
                   help="fraction of windows that should reach --target-snps [0.9]")
    p.add_argument("--round", type=int, default=500, help="width step [500]")
    p.add_argument("--max-width", type=int, default=100000)
    p.add_argument("--chrom-regex", help="windows used (default all)")
    p.add_argument("--mask-bed", help="total mask given to grenedalf; checks it was applied as intended")
    p.add_argument("--ld", help="<out>.ld.params from the ld step; the width is at least the linkage scale")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_windows)

    p = sub.add_parser("chromclass", help="test whether X and mitochondrion depth differ from autosomes")
    p.add_argument("--depth", nargs=2, action="append", required=True, metavar=("POOL", "MQX_BINS"),
                   help="mosdepth bins at the MAPQ cut-off; repeat per pool")
    p.add_argument("--fai", required=True)
    p.add_argument("--mask-bed", help="accessibility mask (bins mostly masked are ignored)")
    p.add_argument("--chrom-regex", default=".", help="analysed nuclear chromosomes")
    p.add_argument("--x-regex", default=r"chrX|_X_|^X$")
    p.add_argument("--mito-regex", default=r"mitochondri|^chrM$|^MT$|mtDNA")
    p.add_argument("--pool-sizes", help="haploid genomes per pool: one number, or POOL=N,POOL=N")
    p.add_argument("--male-fraction", help="known fraction of males per pool: one number, or POOL=m,POOL=m "
                   "(default: estimated from X depth)")
    p.add_argument("--min-difference", type=float, default=0.05,
                   help="smallest |ratio - 1| treated as a real difference [0.05]")
    p.add_argument("--min-fold", type=float, default=0.5)
    p.add_argument("--max-fold", type=float, default=2.0)
    p.add_argument("--abs-min-depth", type=int, default=10)
    p.add_argument("--bootstrap", type=int, default=500)
    p.add_argument("--block", type=int, default=50, help="bins per bootstrap block [50]")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_chromclass)

    p = sub.add_parser("commands", help="write grenedalf_commands.sh with the chosen settings")
    p.add_argument("--prefix", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--bam", nargs="+", required=True)
    p.add_argument("--x-regex", default=r"chrX|_X_|^X$")
    p.add_argument("--total-mask", required=True)
    p.add_argument("--sample-masks", default="")
    p.add_argument("--pool-size-arg", required=True, help="value for grenedalf --pool-sizes (number or file)")
    p.add_argument("--min-count", required=True)
    p.add_argument("--min-depth", required=True)
    p.add_argument("--max-depth", required=True)
    p.add_argument("--window", required=True)
    p.add_argument("--stride", required=True)
    p.add_argument("--min-snps", default="")
    p.add_argument("--map-qual", default=30)
    p.add_argument("--base-qual", default=30)
    p.add_argument("--fst-method", default="unbiased-nei")
    p.add_argument("--threads", default=4)
    p.add_argument("--grenedalf-out", default="grenedalf_out")
    p.set_defaults(func=cmd_commands)

    p = sub.add_parser("ld", help="linkage scale from the decay of between-pool frequency correlation",
                       description=cmd_ld_doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--frequency", required=True, help="grenedalf frequency table (--write-sample-counts)")
    p.add_argument("--region", nargs="+", help="chromosomes to use (default all)")
    p.add_argument("--mask-bed", help="total mask, as given to diversity/fst")
    p.add_argument("--min-count", type=int, default=2)
    p.add_argument("--min-depth", type=int, default=10)
    p.add_argument("--max-depth", type=int, default=100000)
    p.add_argument("--min-maf", type=float, default=0.05, help="minimum mean minor allele frequency [0.05]")
    p.add_argument("--min-dist", type=float, default=10)
    p.add_argument("--max-dist", type=float, default=100000)
    p.add_argument("--bins", type=int, default=30, help="log-spaced distance bins [30]")
    p.add_argument("--anchors", type=int, default=50000, help="anchor SNPs per chromosome [50000]")
    p.add_argument("--min-pairs", type=int, default=200, help="pairs needed for a bin to be used [200]")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_ld)

    p = sub.add_parser("assess", help="technical comparison of grenedalf runs over a grid of settings")
    p.add_argument("--run", nargs="+", action="append", required=True, metavar="NAME TABLE",
                   help="run name, then its grenedalf diversity and/or fst tables; repeat per run")
    p.add_argument("--depth", nargs=2, action="append", metavar=("POOL", "MQX_BINS"),
                   help="mosdepth bins at the MAPQ cut-off, named as the grenedalf samples")
    p.add_argument("--min-snps", type=int, default=50, help="windows with fewer SNPs are not scored [50]")
    p.add_argument("--strict", help="reference run (default: the run with the fewest callable bp)")
    p.add_argument("--tolerance", type=float, default=0.05,
                   help="artefact increase over the strict run still accepted [0.05]")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_assess)

    p = sub.add_parser("report", help="HTML report of all step outputs")
    p.add_argument("--prefix", required=True, help="output prefix of the other steps (e.g. OUT/optimiser)")
    p.add_argument("--out-dir", help="OUT directory, for params.sh and profile/*.dist (default: prefix directory)")
    p.add_argument("--dist", nargs=2, action="append", metavar=("POOL", "REGION_DIST"),
                   help="mosdepth region.dist per pool (default: OUT/profile/*.mosdepth.region.dist.txt)")
    p.add_argument("--title", help="report title")
    p.add_argument("--out", required=True, help="HTML file to write")
    p.set_defaults(func=lambda a: __import__("optimiser_report").build_report(a))

    args = ap.parse_args()
    if args.cmd == "assess" and any(len(r) < 2 for r in args.run):
        ap.error("--run needs a name and at least one table")
    args.func(args)


if __name__ == "__main__":
    main()
