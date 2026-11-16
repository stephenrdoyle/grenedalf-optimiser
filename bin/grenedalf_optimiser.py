#!/usr/bin/env python3
"""
grenedalf_optimiser: choose grenedalf filter settings from the technical properties of
the pool-seq BAMs, before any biological signal is looked at.

Nothing here uses population labels, phenotypes, or the location of peaks. The
settings are judged only on technical grounds: how much of the genome stays
callable, whether every pool is callable at the same sites, and whether the
statistics (theta pi, Watterson's theta, FST) depend on read depth or on
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
  variance    what kind of depth variation there is: spread against read-sampling noise,
              how much persists at 1-10 kb (blocks vs noise), and whether low depth is
              shared by the pools (technical) or specific to one (population divergence).
  depthmask   per-pool depth masks from mosdepth --quantize output; either merged
              into one mask (shared: every pool is called at the same sites) or kept
              per sample (grenedalf --filter-mask-samples-bed-list).
  windows     window width from SNP density in a pilot grenedalf table (the smallest
              width at which most windows hold enough SNPs for a precise estimate) and,
              with --ld, no narrower than the linkage scale from the ld step.
              Also checks that grenedalf applied the mask the right way round.
  chromclass  do the X and mitochondrion behave like autosomes (depth ratio, bootstrap
              interval)? If not, they get their own grenedalf runs and pool sizes.
  ld          linkage scale: how fast the correlation of allele frequencies across
              pools decays with distance between SNPs (an LD decay for pool-seq).
  combine     the best accepted choice of each setting, together, as one more grid run.
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
    "excluded": 1,         # not in the nuclear analysis: --exclude-chroms (mitochondrion) or not --chrom-regex
    "gap": 2,              # N in the reference
    "softmask": 4,         # lower case in the reference (only with --softmask)
    "low_mappability": 8,  # most k-mers overlapping the base are not unique (GenMap)
    "collapsed": 16,       # consensus depth well above the chromosome median
    "ambiguous": 32,       # consensus: most reads have MAPQ below the cut-off
    "low_depth": 64,       # consensus depth well below the chromosome median
}

# Tajima's D is not used: grenedalf computes it only with a minimum count of exactly 2
STAT_SUFFIXES = (".theta_pi", ".theta_watterson", ".fst")


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


def write_mask(df, path, fai):
    """Write a mask BED for grenedalf. grenedalf takes a mask's chromosome lengths from its
    last masked position per chromosome, and requires the total mask and every per-sample
    mask to agree exactly. So every mask also masks the last base of every reference
    chromosome: all masks then describe the same chromosomes, at the cost of one base each."""
    ends = pd.DataFrame({"chrom": fai["chrom"], "start": fai["length"] - 1, "end": fai["length"]})
    write_bed(merge_intervals(pd.concat([df[["chrom", "start", "end"]], ends], ignore_index=True)), path)


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


def name_list(spec):
    """Contig names from "a,b,c" (spaces allowed); empty or None -> empty set."""
    return {x.strip() for x in str(spec or "").split(",") if x.strip()}


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
    exclude = name_list(args.exclude_chroms)

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
            if chrom in exclude or not keep.search(chrom):
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
    total = {"chrom": "total", **summary.drop(columns="chrom").sum().to_dict()}
    summary = pd.concat([summary, pd.DataFrame([total])], ignore_index=True)
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


def nb_cdf(mu, size, kmax, lg_k1=None, lg_ks=None):
    """Negative binomial (mean mu, dispersion size) cdf at 0..kmax."""
    k = np.arange(kmax + 1, dtype=float)
    if lg_k1 is None:
        lg_k1 = np.vectorize(math.lgamma)(k + 1)
    if lg_ks is None:
        lg_ks = np.vectorize(math.lgamma)(k + size)
    p = size / (size + mu)
    logpmf = lg_ks - math.lgamma(size) - lg_k1 + size * math.log(p) + k * math.log1p(-p)
    return np.cumsum(np.exp(logpmf))


def fit_depth_distribution(pmf, from_peak=0.6, q_to=0.99):
    """Fit the main coverage component of a depth distribution: a negative binomial, with a
    fraction pi of sites in a separate low-depth component (reads that fail to map:
    divergent haplotypes, deletions). The fit uses depths from from_peak x the peak of the
    distribution (where most reads are) up to its 99th percentile, where the low component has ended and the
    cdf is pi + (1 - pi) * NB cdf; so the low tail neither drags the mean down nor
    inflates the spread. Returns (mu, size, pi)."""
    cdf = np.cumsum(pmf)
    k_all = np.arange(pmf.size)
    # Anchor on where most reads are (peak of depth x frequency): a dense low-depth tail can
    # outnumber the main peak in sites, but not in reads
    smooth = np.convolve(pmf * k_all, np.ones(5) / 5, mode="same")
    peak = int(np.argmax(smooth[1:]) + 1) if pmf.size > 2 else 1
    x0 = max(int(from_peak * peak), 1)
    x1 = max(int(np.searchsorted(cdf, q_to)), x0 + 2)
    xs = np.unique(np.linspace(x0, x1, 30).astype(int))
    target = cdf[xs]
    kmax = int(x1 * 3 + 20)
    k = np.arange(kmax + 1, dtype=float)
    lg = np.vectorize(math.lgamma)
    lg_k1 = lg(k + 1)
    pis = np.linspace(0, 0.7, 36)[:, None]
    best = (np.inf, float(peak), 10.0, 0.0)
    for size in np.logspace(-0.3, 3, 50):
        lg_ks = lg(k + size)
        for mu in np.linspace(0.5 * peak, 1.3 * peak, 46):
            c = nb_cdf(mu, size, kmax, lg_k1, lg_ks)[xs]
            loss = np.sum((pis + (1 - pis) * c[None, :] - target[None, :]) ** 2, axis=1)
            i = int(np.argmin(loss))
            if loss[i] < best[0]:
                best = (float(loss[i]), mu, size, float(pis[i, 0]))
    return best[1], best[2], best[3]


def nb_quantile(mu, size, q):
    kmax = int(mu * 10 + 50)
    return int(np.searchsorted(nb_cdf(mu, size, kmax), q))


def fitted_bounds(mu, size, q_low, q_high, abs_min):
    lower = max(abs_min, nb_quantile(mu, size, q_low))
    upper = max(lower, nb_quantile(mu, size, q_high))
    return lower, upper


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
    x_chroms = name_list(args.x_chroms)
    summary = pd.read_csv(args.mask_summary, sep="\t", dtype={"chrom": str})
    weights = dict(zip(summary["chrom"], summary["accessible"]))
    min_folds = sorted(set(args.grid_min_folds) | {args.min_fold})
    max_folds = sorted(set(args.grid_max_folds) | {args.max_fold})

    rows, classes, quantize = [], {}, []
    for name, path in args.dist:
        pmfs = read_dist(path)
        by_class = {}
        for chrom, pmf in pmfs.items():
            cls = "X" if chrom in x_chroms else "A"
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
            mu, size, pi_low = fit_depth_distribution(pmf)
            fl, fu = fitted_bounds(mu, size, args.fit_low_q, args.fit_high_q, args.abs_min_depth)
            cuts.update({fl, fu + 1})
            mean = float(np.sum(np.arange(pmf.size) * pmf))
            sd = float(np.sqrt(np.sum((np.arange(pmf.size) - mean) ** 2 * pmf)))
            rows.append({
                "pool": name, "class": cls, "accessible_bp": int(w), "median": median,
                "mean": round(mean, 2), "sd": round(sd, 2), "cv": round(sd / mean, 3) if mean else np.nan,
                "poisson_cv": round(1 / math.sqrt(mean), 3) if mean else np.nan,
                "fit_mu": round(mu, 2), "fit_size": round(size, 3), "fit_low_component": round(pi_low, 3),
                "fit_cv": round(math.sqrt(mu + mu * mu / size) / mu, 3),
                "fit_lower": fl, "fit_upper": fu,
                "fit_frac_below": round(float(cdf[fl - 1]) if fl > 0 else 0.0, 4),
                "fit_frac_above": round(max(float(1 - cdf[fu]), 0.0) if fu < cdf.size else 0.0, 4),
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
        if args.bound_method == "fitted":
            bounds = {cls: (int(r["fit_lower"]), int(r["fit_upper"])) for cls, r in t.iterrows()}
        else:
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
              "BOUND_METHOD": args.bound_method,
              "MIN_DEPTH": min(lowers), "MAX_DEPTH": max(uppers)}
    if args.mode == "shared":
        total = merge_intervals(pd.concat([base] + list(per_pool.values()), ignore_index=True))
        write_mask(total, f"{args.out}.total_mask.bed", fai)
        shared = genome_bp - bp(total[total["chrom"].isin(classes)])
        rows.append({"pool": "shared", "accessible_bp": accessible_bp, "out_of_range_bp": accessible_bp - shared,
                     "out_of_range_frac": round(1 - shared / accessible_bp, 4) if accessible_bp else np.nan})
        params.update({"TOTAL_MASK_BED": f"{args.out}.total_mask.bed", "SAMPLE_MASK_LIST": ""})
    else:
        write_mask(base, f"{args.out}.total_mask.bed", fai)
        with open(f"{args.out}.sample_masks.tsv", "w") as fh:
            for name, out in per_pool.items():
                write_mask(out, f"{args.out}.{name}.depth_mask.bed", fai)
                fh.write(f"{name}\t{args.out}.{name}.depth_mask.bed\n")
        params.update({"TOTAL_MASK_BED": f"{args.out}.total_mask.bed",
                       "SAMPLE_MASK_LIST": f"{args.out}.sample_masks.tsv"})
    summary = pd.DataFrame(rows)
    summary.to_csv(f"{args.out}.depthmask_summary.tsv", sep="\t", index=False)
    write_params(f"{args.out}.depthmask.params", params)
    print(summary.to_string(index=False))


# ---------------------------------------------------------------------------
# variance

def cmd_variance(args):
    """What kind of depth variation is there, and what kind of low depth is being cut?

    Per pool and chromosome class:
      spread       per-base mean, SD and CV, against the CV of pure read-sampling noise
                   (Poisson, 1/sqrt(mean)); var/mean = 1 for Poisson
      scale        share of the per-base variance that persists in 1 kb and 10 kb means.
                   Sampling noise averages away within a few read lengths; regional effects
                   (divergent haplotypes, GC, copy number) do not. High shares: depth varies
                   in blocks; low: it is scattered noise
      low depth    of the 1 kb bins below the pool's lower bound, the share in runs of at
                   least --block-bins bins (blocks, not speckle), and the share that no other
                   pool is low in (pool-specific: population divergence or deletion) or that
                   at least half the pools are low in (shared: reference or technical)
    Between pools: Spearman correlation of 1 kb depth (each pool normalised to its median)."""
    weights = {}
    if args.mask_summary:
        ms = pd.read_csv(args.mask_summary, sep="\t", dtype={"chrom": str})
        weights = dict(zip(ms["chrom"], ms["accessible"]))
    classes = dict(pd.read_csv(args.chrom_class, sep="\t", dtype={"chrom": str}).itertuples(index=False))
    thr = pd.read_csv(args.thresholds, sep="\t").set_index(["pool", "class"])
    dists = dict(args.dist or [])
    mask = merge_intervals(read_bed(args.mask_bed)) if args.mask_bed else None

    rows, norm, low = [], {}, {}
    for name, path in args.depth:
        b = read_mosdepth_bins(path, "d")
        b = b[b["chrom"].isin(classes)].copy()
        b["class"] = b["chrom"].map(classes)
        if mask is not None:   # bins at least half accessible
            frac = np.zeros(len(b))
            for chrom, idx in b.groupby("chrom").indices.items():
                m = mask[mask["chrom"] == chrom]
                if m.empty:
                    continue
                g = b.iloc[idx]
                cs = coverage_prefix(m, int(g["end"].max()))
                frac[idx] = (cs[g["end"].to_numpy()] - cs[g["start"].to_numpy()]) / (g["end"] - g["start"]).to_numpy()
            b = b[frac < 0.5]
        key = b["chrom"] + ":" + b["start"].astype(str)
        pmfs = read_dist(dists[name]) if name in dists else {}
        for cls, g in b.groupby("class"):
            # per-base spread from the mosdepth distribution (accessible genome)
            acc = np.zeros(1)
            for chrom, pmf in pmfs.items():
                w = float(weights.get(chrom, 0))
                if classes.get(chrom) != cls or w <= 0:
                    continue
                n = max(acc.size, pmf.size)
                acc = np.pad(acc, (0, n - acc.size)) + np.pad(pmf * w, (0, n - pmf.size))
            if acc.sum() > 0:
                acc = acc / acc.sum()
                k = np.arange(acc.size)
                mean_b = float(np.sum(k * acc)); var_b = float(np.sum((k - mean_b) ** 2 * acc))
            else:
                mean_b, var_b = float(g["d"].mean()), np.nan
            d1 = g["d"].to_numpy()
            g10 = g.assign(blk=g["start"] // 10_000).groupby(["chrom", "blk"])["d"].agg(["mean", "size"])
            d10 = g10.loc[g10["size"] >= 5, "mean"].to_numpy()
            lower = int(thr.loc[(name, cls), "lower"]) if (name, cls) in thr.index else np.nan
            is_low = d1 < lower
            # runs of consecutive low bins (adjacent starts on the same chromosome)
            gs = g.assign(low=is_low).sort_values(["chrom", "start"])
            brk = (gs["low"] != gs["low"].shift()) | (gs["chrom"] != gs["chrom"].shift()) | \
                  (gs["start"] - gs["end"].shift() != 0)
            run_id = brk.cumsum()
            run_len = gs.groupby(run_id)["low"].transform("size")
            in_block = (gs["low"] & (run_len >= args.block_bins)).sum()
            rows.append({
                "pool": name, "class": cls, "bins": len(g),
                "mean": round(mean_b, 2), "sd": round(math.sqrt(var_b), 2) if np.isfinite(var_b) else np.nan,
                "cv": round(math.sqrt(var_b) / mean_b, 3) if np.isfinite(var_b) and mean_b else np.nan,
                "poisson_cv": round(1 / math.sqrt(mean_b), 3) if mean_b else np.nan,
                "var_over_mean": round(var_b / mean_b, 2) if np.isfinite(var_b) and mean_b else np.nan,
                "share_1kb": round(float(np.var(d1)) / var_b, 3) if np.isfinite(var_b) and var_b else np.nan,
                "share_10kb": round(float(np.var(d10)) / var_b, 3) if np.isfinite(var_b) and var_b and d10.size else np.nan,
                "low_bins_frac": round(float(is_low.mean()), 4),
                "low_in_blocks": round(float(in_block / max(is_low.sum(), 1)), 3)})
            med = float(np.median(d1[d1 > 0])) if (d1 > 0).any() else np.nan
            norm.setdefault(name, []).append(pd.Series(d1 / med, index=key[g.index]))
            low.setdefault(name, []).append(pd.Series(is_low, index=key[g.index]))

    t = pd.DataFrame(rows)
    nmat = pd.DataFrame({p: pd.concat(v) for p, v in norm.items()}).dropna()
    lmat = pd.DataFrame({p: pd.concat(v) for p, v in low.items()}).fillna(False).astype(bool)
    corr = nmat.rank().corr()
    pools = list(lmat.columns)
    key_class = pd.Series(lmat.index.str.rsplit(":", n=1).str[0], index=lmat.index).map(classes)
    share = []
    for cls in sorted(key_class.dropna().unique()):
        lc = lmat[key_class == cls]
        n_low = lc.sum(axis=1)
        for p in pools:
            lp = lc[p]
            nl = int(lp.sum())
            share.append({"pool": p, "class": cls, "low_bins": nl,
                          "pool_specific": round(float((lp & (n_low == 1)).sum() / max(nl, 1)), 3),
                          "shared_half_or_more": round(float((lp & (n_low >= math.ceil(len(pools) / 2))).sum() / max(nl, 1)), 3),
                          "shared_all": round(float((lp & (n_low == len(pools))).sum() / max(nl, 1)), 3)})
    sh = pd.DataFrame(share)
    t = t.merge(sh, on=["pool", "class"], how="left")
    t["mean_corr_with_others"] = t["pool"].map(
        lambda p: round(float(corr.loc[p, [q for q in pools if q != p]].mean()), 3) if len(pools) > 1 else np.nan)
    t.to_csv(f"{args.out}.variance.tsv", sep="\t", index=False)
    corr.round(3).to_csv(f"{args.out}.depth_correlation.tsv", sep="\t")
    print(t.to_string(index=False))
    print("\n1 kb depth correlation between pools (Spearman)\n" + corr.round(3).to_string())


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
        inverted = win_total - exp_total
        diff = abs(obs_total - exp_total) / max(exp_total, 1)
        if exp_total == 0:
            verdict = "mask has no overlap with these windows"
        elif obs_total == 0:
            verdict = "NOT APPLIED: total.masked is 0"
        elif diff <= 0.05:
            verdict = (f"OK: grenedalf masked the BED intervals (within {100 * diff:.1f}%; small differences come "
                       "from how grenedalf counts masked positions at window and chromosome edges)")
        elif abs(obs_total - inverted) < abs(obs_total - exp_total):
            verdict = "WRONG WAY ROUND: grenedalf masked everything except the BED (check --filter-mask-total-bed-invert)"
        else:
            verdict = (f"MISMATCH: orientation is right but the counts differ by {100 * diff:.1f}%; "
                       "compare total.masked with the BED by hand")
        lines.append(f"mask check: expected {exp_total:,} masked bp, grenedalf reported {obs_total:,} -> {verdict}")

    write_params(f"{args.out}.windows.params", {"WINDOW": rec, "STRIDE": rec // 2, "MIN_SNPS": args.target_snps,
                                                "WINDOW_SNPS": w_snps})
    with open(f"{args.out}.windows.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


# ---------------------------------------------------------------------------
# chromclass

def density_mode(x, grid=1024):
    """Peak of a distribution: the maximum of a Gaussian kernel density estimate (robust
    Silverman bandwidth). Unlike the median, unaffected by a heavy tail on one side."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return float(np.median(x)) if x.size else np.nan
    lo, hi = np.percentile(x, [0.5, 99.5])
    q1, q3 = np.percentile(x, [25, 75])
    scale = min(np.std(x), (q3 - q1) / 1.349) or np.std(x) or 1.0
    bw = 0.9 * scale * x.size ** -0.2
    edges = np.linspace(lo, hi, grid + 1)
    counts, _ = np.histogram(x, bins=edges)
    step = edges[1] - edges[0]
    k = np.arange(-int(4 * bw / step) - 1, int(4 * bw / step) + 2) * step
    dens = np.convolve(counts, np.exp(-0.5 * (k / bw) ** 2), mode="same")
    return float((edges[:-1] + step / 2)[int(np.argmax(dens))])


def block_bootstrap_ratio(x, a, n_boot, block, rng, stat=np.median):
    """95% interval of stat(x) / stat(a), resampling blocks of adjacent bins."""
    def blocks(v):
        k = max(len(v) // block, 1)
        return [v[i * block:(i + 1) * block] for i in range(k)] if len(v) >= block else [v]
    bx, ba = blocks(x), blocks(a)
    out = np.empty(n_boot)
    for i in range(n_boot):
        sx = np.concatenate([bx[j] for j in rng.integers(0, len(bx), len(bx))])
        sa = np.concatenate([ba[j] for j in rng.integers(0, len(ba), len(ba))])
        out[i] = stat(sx) / stat(sa)
    return np.percentile(out, [2.5, 97.5])


def cmd_chromclass(args):
    """Do the X and the mitochondrion behave like autosomes? Depth ratios to autosomes per pool.

    X: in a pool with males (hemizygous), X depth is below autosome depth, and the X is
    sampled from fewer copies, so grenedalf needs a smaller X pool size: with a fraction m
    of males, X depth / autosome depth = 1 - m/2, and the X pool size is pool_size * (1 - m/2).
    A 50:50 pool gives 0.75. The ratio is only used when X and autosome depth have similar
    spread (robust CV within --max-spread-ratio): otherwise depth is not tracking copy number
    alike on both, as when divergence from the reference widens the autosomal distribution.

    The ratio is taken between the peaks (kernel density mode) of the two depth distributions,
    over well-mapped bins only (most reads pass the MAPQ cut-off). Medians are pulled down
    wherever reads fail to map, e.g. haplotypes divergent from the reference, which affects
    autosomes and X differently. A ratio outside 0.5-1 (beyond --min-difference) is
    impossible for a sex ratio and is reported as unreliable: the male fraction then comes
    from --male-fraction, or 0.5 is assumed. With --male-fraction given, the known value is
    always used and the estimate is reported as a check.

    Mitochondrion: haploid, many copies per cell, one lineage per individual: always run
    separately, with depth bounds per pool (copy number differs between pools) applied as
    per-pool masks, and a pool size of the number of individuals."""
    x_chroms, mito_chroms, keep = name_list(args.x_chroms), name_list(args.mito_chroms), re.compile(args.chrom_regex)
    fai = read_fai(args.fai)
    mask = merge_intervals(read_bed(args.mask_bed)) if args.mask_bed else None
    names = [n for n, _ in args.depth]
    sizes = parse_pool_sizes(args.pool_sizes, names)
    known_m = parse_pool_sizes(args.male_fraction, names)
    mq0 = dict(args.depth_mq0 or [])
    rng = np.random.default_rng(args.seed)

    def chrom_class(c):
        if c in mito_chroms:
            return "M"
        if not keep.search(c):
            return None
        return "X" if c in x_chroms else "A"

    rows, mito_rows, mito_masks = [], [], {}
    for name, path in args.depth:
        b = read_mosdepth_bins(path, "depth")
        if name in mq0:
            b0 = read_mosdepth_bins(mq0[name], "depth0")
            b = b.merge(b0, on=["chrom", "start", "end"], how="left")
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
        # Well-mapped bins: most reads pass the MAPQ cut-off (needs the MAPQ 0 bins)
        if "depth0" in b:
            with np.errstate(invalid="ignore", divide="ignore"):
                mapped = (b["depth"] / b["depth0"]).to_numpy()
            well = (mapped >= args.min_mapped_ratio) & (b["depth"].to_numpy() > 0)
        else:
            well = b["depth"].to_numpy() > 0
        a_all = b.loc[b["class"] == "A", "depth"].to_numpy()
        a = b.loc[(b["class"] == "A") & well, "depth"].to_numpy()
        if a.size == 0:
            sys.exit(f"{name}: no well-mapped autosomal bins")
        mode_a = density_mode(a)
        n = sizes.get(name, np.nan)
        rows.append({"pool": name, "class": "A", "bins": a.size, "median_depth": round(float(np.median(a_all)), 2),
                     "mode_depth": round(mode_a, 2), "ratio_to_A": 1.0, "ratio_median": 1.0, "ci_low": 1.0,
                     "ci_high": 1.0, "differs": False, "pool_size": n})

        x_all = b.loc[b["class"] == "X", "depth"].to_numpy()
        x = b.loc[(b["class"] == "X") & well, "depth"].to_numpy()
        if x.size:
            ratio = density_mode(x) / mode_a
            lo, hi = block_bootstrap_ratio(x, a, args.bootstrap, args.block, rng, stat=density_mode)
            # Depth only measures copy number if it behaves alike on X and autosomes: if one
            # is much more spread out (e.g. reference divergence on the autosomes), the ratio
            # of their peaks says nothing reliable about the sex ratio.
            rcv = lambda v: (np.subtract(*np.percentile(v, [75, 25])) / 1.349) / np.median(v) if v.size else np.nan
            cv_a, cv_x = float(rcv(a)), float(rcv(x))
            spread = max(cv_a, cv_x) / min(cv_a, cv_x) if min(cv_a, cv_x) > 0 else np.inf
            in_range = (0.5 - args.min_difference) <= ratio <= (1 + args.min_difference)
            plausible = bool(in_range and spread <= args.max_spread_ratio)
            m_est = float(np.clip(2 * (1 - ratio), 0, 1))
            row = {"pool": name, "class": "X", "bins": x.size, "median_depth": round(float(np.median(x_all)), 2),
                   "mode_depth": round(density_mode(x), 2), "ratio_to_A": round(ratio, 4),
                   "ratio_median": round(float(np.median(x_all) / np.median(a_all)), 4),
                   "ci_low": round(float(lo), 4), "ci_high": round(float(hi), 4),
                   "differs": bool((hi < 1 or lo > 1) and abs(ratio - 1) > args.min_difference),
                   "cv_A": round(cv_a, 3), "cv_X": round(cv_x, 3), "spread_ratio": round(spread, 2),
                   "plausible": plausible,
                   "unreliable_because": "" if plausible else
                   ("ratio outside 0.5-1" if not in_range else "") + (" and " if not in_range and spread > args.max_spread_ratio else "")
                   + (f"X and autosome depth spread differ {spread:.1f}x" if spread > args.max_spread_ratio else ""),
                   "male_fraction_est": round(m_est, 3) if plausible else np.nan,
                   "male_fraction_ci": f"{np.clip(2 * (1 - hi), 0, 1):.2f}-{np.clip(2 * (1 - lo), 0, 1):.2f}"
                   if plausible else "NA"}
            if name in known_m:
                m = known_m[name]
                row["male_fraction_known"] = m
                row["expected_ratio"] = round(1 - m / 2, 4)
                row["matches_known"] = bool(lo <= 1 - m / 2 <= hi) if plausible else np.nan
            row["_n"] = n
            rows.append(row)

        mt = b[b["class"] == "M"]
        if len(mt):
            med = float(np.median(mt["depth"]))
            lower, upper = depth_bounds(med, args.min_fold, args.max_fold, args.abs_min_depth)
            lo, hi = block_bootstrap_ratio(mt["depth"].to_numpy(), a, args.bootstrap, args.block, rng)
            out = mt[(mt["depth"] < lower) | (mt["depth"] > upper)]
            mito_masks[name] = out[["chrom", "start", "end"]]
            mito_rows.append({"pool": name, "median_depth": round(med, 1), "ratio_to_A": round(med / mode_a, 2),
                              "ci_low": round(float(lo), 2), "ci_high": round(float(hi), 2),
                              "lower": lower, "upper": upper,
                              "masked_bp": int((out["end"] - out["start"]).sum()),
                              "pool_size": round(n / 2) if np.isfinite(n) else np.nan})

    # Male fraction used for the X pool size. Whatever pushes one pool's X depth outside what a
    # sex ratio allows (reference divergence, GC, mapping) affects every pool, so if any
    # estimate is implausible, none is used: known values where given, else 0.5 for all.
    xrows = [r for r in rows if r["class"] == "X"]
    all_plausible = all(r["plausible"] for r in xrows)
    for r in xrows:
        n = r.pop("_n")
        if r["pool"] in known_m:
            m, source = known_m[r["pool"]], "known"
        elif all_plausible:
            m, source = r["male_fraction_est"], "estimated"
        else:
            m, source = 0.5, "assumed (X depth estimates implausible)"
        r["male_fraction_used"] = round(m, 3)
        r["male_fraction_source"] = source
        r["pool_size"] = round(n * (1 - m / 2)) if np.isfinite(n) else np.nan

    t = pd.DataFrame(rows)
    t.to_csv(f"{args.out}.chrom_classes.tsv", sep="\t", index=False)
    x = t[t["class"] == "X"]
    params = {"MITO_CHROMS": ",".join(c for c in fai["chrom"] if c in mito_chroms)}
    if len(x):
        # X gets its own run whenever its pool size differs from the autosomal one, or its depth does
        own_size = (x["male_fraction_used"] > 0.025).any()
        params["X_SEPARATE"] = int(bool(own_size or x["differs"].any()))
        params["X_RATIO_UNRELIABLE"] = int(not x["plausible"].all())
        if x["pool_size"].notna().all():
            with open(f"{args.out}.x_pool_sizes.tsv", "w") as fh:
                for r in x.itertuples():
                    fh.write(f"{r.pool}\t{int(r.pool_size)}\n")
            params["X_POOL_SIZES"] = f"{args.out}.x_pool_sizes.tsv"
    else:
        params["X_SEPARATE"] = 0
    if mito_rows:
        mt = pd.DataFrame(mito_rows)
        mt.to_csv(f"{args.out}.mito.tsv", sep="\t", index=False)
        with open(f"{args.out}.mito_sample_masks.tsv", "w") as fh:
            for name, m in mito_masks.items():
                path = f"{args.out}.mito_mask.{name}.bed"
                write_mask(m, path, fai)
                fh.write(f"{name}\t{path}\n")
        params.update({"MITO_SEPARATE": 1, "MITO_SAMPLE_MASKS": f"{args.out}.mito_sample_masks.tsv",
                       # grenedalf's global bounds only as a backstop: the per-pool masks apply each pool's range
                       "MITO_MIN_DEPTH": int(mt["lower"].min()), "MITO_MAX_DEPTH": int(mt["upper"].max())})
        if mt["pool_size"].notna().all():
            with open(f"{args.out}.mito_pool_sizes.tsv", "w") as fh:
                for r in mt.itertuples():
                    fh.write(f"{r.pool}\t{int(r.pool_size)}\n")
            params["MITO_POOL_SIZES"] = f"{args.out}.mito_pool_sizes.tsv"
    else:
        params["MITO_SEPARATE"] = 0
    write_params(f"{args.out}.chromclass.params", params)
    print(t.to_string(index=False))
    if mito_rows:
        print("\n" + pd.DataFrame(mito_rows).to_string(index=False))
    if len(x) and params["X_RATIO_UNRELIABLE"]:
        why = "; ".join(sorted(set(x.loc[~x["plausible"].astype(bool), "unreliable_because"])))
        print(f"\nWARNING: X / autosome depth cannot be read as a sex ratio here ({why}). The male fraction "
              "comes from --male-fraction, or 0.5 is assumed; set it from the study design.", file=sys.stderr)


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
    mito = [c for c in cc.get("MITO_CHROMS", "").split(",") if c]
    named_x = name_list(args.x_chroms)
    x_chroms = [c for c in fai["chrom"] if c in named_x and c not in mito]
    x_sep = cc.get("X_SEPARATE", "0") == "1" and x_chroms
    base = merge_intervals(read_bed(args.total_mask))
    whole = lambda chroms: fai[fai["chrom"].isin(chroms)].assign(start=0).rename(columns={"length": "end"})[
        ["chrom", "start", "end"]]

    runs = []
    if x_sep:
        # Autosome run: X added to the mask. X run: everything but X added to the mask.
        auto_mask = f"{P}.mask_autosomes.bed"
        write_mask(pd.concat([base, whole(x_chroms)]), auto_mask, fai)
        x_mask = f"{P}.mask_X.bed"
        write_mask(pd.concat([base, whole([c for c in fai["chrom"] if c not in x_chroms])]), x_mask, fai)
        runs.append(("autosomes", "autosomes (X masked; it gets its own pool sizes below)", auto_mask,
                     args.pool_size_arg, args.min_depth, args.max_depth, args.sample_masks))
        runs.append(("X", "X chromosome, with per-pool X pool sizes (fewer copies in pools with males)", x_mask,
                     cc.get("X_POOL_SIZES", args.pool_size_arg), args.min_depth, args.max_depth, args.sample_masks))
    else:
        runs.append(("nuclear", "autosomes" + (" and X (X depth matches the autosomes)" if x_chroms else ""),
                     args.total_mask, args.pool_size_arg, args.min_depth, args.max_depth, args.sample_masks))
    if mito and cc.get("MITO_SEPARATE", "0") == "1":
        m_mask = f"{P}.mask_mito.bed"
        write_mask(whole([c for c in fai["chrom"] if c not in mito]), m_mask, fai)
        runs.append(("mito", "mitochondrion: haploid, pool size = individuals, its own depth bounds "
                     "(whole genome as one window)", m_mask, cc.get("MITO_POOL_SIZES", args.pool_size_arg),
                     cc.get("MITO_MIN_DEPTH", args.min_depth), cc.get("MITO_MAX_DEPTH", args.max_depth),
                     cc.get("MITO_SAMPLE_MASKS", "")))

    L = ["#!/usr/bin/env bash",
         f"# grenedalf commands recommended by grenedalf-optimiser ({pd.Timestamp.now():%Y-%m-%d}).",
         "# Settings chosen on technical criteria only; see optimiser_report.html for why.",
         "set -euo pipefail", "",
         f"OUTDIR={args.grenedalf_out}", f"THREADS={args.threads}", ""]
    for key, desc, mask_bed, pools, mind, maxd, smasks in runs:
        mito_run = key == "mito"
        L.append(f"# --- {desc}")
        # Tajima's D is not computed (grenedalf allows it only with a minimum count of exactly 2)
        mc = str(args.min_count)
        jobs = [("diversity", mc, "", ["--no-tajima-d"]), ("fst", mc, "", [f"--method {args.fst_method}"])]
        for cmd, count, tag, extra in jobs:
            c = [f"grenedalf {cmd}"]
            c += [f"--sam-path {b}" for b in args.bam]
            c += [f"--reference-genome-fai {args.fai}",
                  f"--sam-min-map-qual {args.map_qual} --sam-min-base-qual {args.base_qual}",
                  f"--filter-mask-total-bed {mask_bed}"]
            if smasks:
                c.append(f"--filter-mask-samples-bed-list {smasks}")
            c += [f"--filter-sample-min-count {count}",
                  f"--filter-sample-min-read-depth {mind} --filter-sample-max-read-depth {maxd}"]
            if mito_run:
                c.append("--window-type chromosomes")
            else:
                c.append(f"--window-type interval --window-interval-width {args.window} "
                         f"--window-interval-stride {args.stride}")
            c += ["--window-average-policy valid-loci", f"--pool-sizes {pools}"] + extra
            c += ["--separator-char tab --threads $THREADS",
                  f"--out-dir $OUTDIR/{key} --file-prefix {key}_{tag}"]
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
        frames = [read_grenedalf(p) for p in tables if p != "-" and os.path.exists(p)]
        t = frames[0]
        for f in frames[1:]:
            stat_cols = [c for c in f.columns if (c.endswith(STAT_SUFFIXES) or re.search(r":.+\.(passed|missing)$", c))
                         and c not in t]
            t = t.merge(f[["chrom", "start", "end"] + stat_cols], on=["chrom", "start", "end"], how="outer")
        # Windows overlap when stride < width; divide sums by width / stride to count each base once
        strides = t.groupby("chrom")["start"].apply(lambda x: np.diff(np.sort(x.to_numpy()))).explode()
        stride = float(strides.astype(float).median()) if strides.notna().any() else float(t["width"].median())
        step = max(float(t["width"].median()) / max(stride, 1.0), 1.0)
        windows_all = len(t)
        callable_bp = int(t["valid"].sum() / step)
        snps = int(t["total.passed"].sum() / step)
        full = t
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
        # Data kept, like for like across mask modes. grenedalf's total.* columns count a site
        # if any pool can use it, which flatters per-pool masks; what FST can use is a pair's
        # own sites (<a>:<b>.passed), and what theta can use is a pool's own (<pool>.passed).
        pair_cols = [c for c in full.columns if re.search(r":.+\.passed$", c)]
        pool_cols = [c for c in full.columns if re.search(r"^[^:]+\.passed$", c) and not c.startswith("total.")]
        pair_snps = float(full[pair_cols].sum().mean() / step) if pair_cols else np.nan
        pool_snps = float(full[pool_cols].sum().mean() / step) if pool_cols else np.nan
        imbalance = np.nan
        if pair_cols:
            pp = t[[c for c in pair_cols if c in t]].to_numpy(dtype=float)
            with np.errstate(invalid="ignore", divide="ignore"):
                imbalance = float(np.nanmedian((pp.max(axis=1) - pp.min(axis=1)) / pp.mean(axis=1)))
        row = {"run": run, "windows": windows_all, "windows_scored": len(t),
               "callable_bp": callable_bp, "snps": snps,
               "pair_snps": round(pair_snps) if np.isfinite(pair_snps) else np.nan,
               "pool_snps": round(pool_snps) if np.isfinite(pool_snps) else np.nan,
               "median_callable_frac": round(float(t["callable"].median()), 4) if len(t) else np.nan,
               "pair_callable_imbalance": round(imbalance, 4) if np.isfinite(imbalance) else np.nan}
        # For each statistic type (theta_pi, theta_watterson, fst) and technical covariate:
        # the median |rho| over pools or pairs
        worst, worst_name = 0.0, ""
        for cov in ("rho_depth", "rho_callable", "rho_depth_imbalance"):
            if cov not in s:
                continue
            for typ, g in s.groupby("type"):
                v = g[cov].abs().dropna()
                if len(v):
                    v = float(v.median())
                    row[f"{typ}.{cov}"] = round(v, 4)
                    if v > worst:
                        worst, worst_name = v, f"{typ} vs {cov[4:]}"
        row["artefact"], row["worst"] = round(worst, 4), worst_name
        run_rows.append(row)

    runs_df = pd.DataFrame(run_rows)
    stats_df = pd.DataFrame(stat_rows)

    # Data kept: what FST can use (per pair) when available, else per pool, else all SNPs
    for metric in ("pair_snps", "pool_snps", "snps"):
        if metric in runs_df and runs_df[metric].notna().all():
            break
    runs_df["data_metric"] = metric
    runs_df["data_kept"] = runs_df[metric]

    # Judged per statistic against a reference run (the one keeping least data, unless
    # --strict): a run is accepted only if no statistic's dependence on any technical
    # covariate rises by more than --tolerance. Some dependence is real (divergent
    # haplotypes map less well), hence the comparison rather than an absolute threshold,
    # and one statistic's improvement cannot buy another's deterioration.
    checks = set(args.check_run or [])
    candidates = runs_df[~runs_df["run"].isin(checks)]
    strict = args.strict if args.strict else candidates.sort_values(["data_kept", "snps"], kind="stable")["run"].iloc[0]
    ref = runs_df.set_index("run").loc[strict]
    metric_cols = [c for c in runs_df.columns if ".rho_" in c]
    deltas = runs_df[["run"]].copy()
    for c in metric_cols:
        deltas[c] = (runs_df[c] - ref[c]).round(4)
    d = deltas.set_index("run")[metric_cols]
    runs_df["max_increase"] = d.max(axis=1, skipna=True).round(4).to_numpy()
    on = d.idxmax(axis=1, skipna=True).fillna("").str.replace(".rho_", " vs ").to_numpy()
    runs_df["max_increase_on"] = np.where(runs_df["max_increase"].fillna(0) > 0, on, "")
    runs_df["acceptable"] = (runs_df["max_increase"].fillna(0) <= args.tolerance) & ~runs_df["run"].isin(checks)
    runs_df["role"] = np.where(runs_df["run"].isin(checks), "check", np.where(runs_df["run"] == strict, "reference", "candidate"))
    ok = runs_df[runs_df["acceptable"]]
    key = list(zip(ok["data_kept"], ok["snps"]))
    best = ok["run"].iloc[key.index(max(key))] if key else strict  # ties: first run listed
    runs_df["recommended"] = runs_df["run"] == best

    runs_df.to_csv(f"{args.out}.assess_runs.tsv", sep="\t", index=False)
    deltas.to_csv(f"{args.out}.assess_deltas.tsv", sep="\t", index=False)
    stats_df.round(4).to_csv(f"{args.out}.assess_stats.tsv", sep="\t", index=False)
    show = ["run", "role", "data_kept", "callable_bp", "artefact", "max_increase", "max_increase_on", "acceptable",
            "recommended"]
    print(runs_df[show].to_string(index=False))
    print(f"\ndata kept measured as {metric}; reference run: {strict}; tolerance {args.tolerance} on every statistic")
    print(f"most permissive run that adds no technical artefact to any statistic: {best}")
    with open(f"{args.out}.assess.params", "w") as fh:
        fh.write(f"RECOMMENDED_RUN={best}\nTOLERANCE={args.tolerance}\nREFERENCE_RUN={strict}\n")


# ---------------------------------------------------------------------------
# combine

def cmd_combine(args):
    """The best accepted choice of each setting, together, as one more grid run.

    The grid changes one setting at a time from the baseline. For each setting (minimum
    count, lower depth bound, upper depth bound, bound method, mask mode) take the accepted
    run of that kind that keeps the most data, if it keeps more than the baseline, and
    combine them. Prints the combined run's spec, or nothing if it equals an existing run."""
    spec = pd.read_csv(args.runs, sep="\t", dtype=str)
    res = pd.read_csv(args.assess, sep="\t")
    m = spec.merge(res[["run", "acceptable", "data_kept"]], on="run", how="left")
    base = m[m["setting"] == "baseline"].iloc[0]
    combo = base.copy()
    chosen = []
    for setting, fields in (("count", ["mc"]), ("lower", ["min_fold"]), ("upper", ["max_fold"]),
                            ("method", ["method"]), ("mode", ["mode"])):
        c = m[(m["setting"] == setting) & (m["acceptable"].astype(str) == "True")
              & (m["data_kept"].astype(float) > float(base["data_kept"]))]
        if len(c):
            best = c.sort_values("data_kept").iloc[-1]
            for f in fields:
                combo[f] = best[f]
            chosen.append(best["run"])
    if combo["method"] == "fitted":   # fitted bounds replace the folds
        combo["min_fold"], combo["max_fold"] = base["min_fold"], base["max_fold"]
        fold_runs = set(m.loc[m["setting"].isin(["lower", "upper"]), "run"])
        chosen = [c for c in chosen if c not in fold_runs]
    fields = ["mc", "min_fold", "max_fold", "mode", "base", "method"]
    same = (m[fields].astype(str) == combo[fields].astype(str)).all(axis=1)
    if len(chosen) < 2 or same.any():
        log("combination: nothing new to combine" + (f" (best: {', '.join(chosen)})" if chosen else ""))
        return
    log(f"combination of {', '.join(chosen)}")
    print("combined " + " ".join(str(combo[f]) for f in fields))


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
    p.add_argument("--exclude-chroms", default="", help="contigs left out of the nuclear analysis, comma "
                   "separated (the mitochondrion: it gets its own run)")
    p.add_argument("--chrom-regex", default=".", help="optionally, only contigs matching this are analysed [all]")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_mask)

    p = sub.add_parser("thresholds", help="per-pool depth bounds and minimum allele count")
    p.add_argument("--dist", nargs=2, action="append", required=True, metavar=("POOL", "REGION_DIST"),
                   help="mosdepth .mosdepth.region.dist.txt from --by accessible.bed; repeat per pool")
    p.add_argument("--mask-summary", required=True, help="<out>.mask_summary.tsv from the mask step")
    p.add_argument("--stats", nargs=2, action="append", metavar=("POOL", "SAMTOOLS_STATS"))
    p.add_argument("--x-chroms", default="", help="sex chromosome (X or Z) contig names, comma separated; "
                   "all other contigs are autosomal")
    p.add_argument("--min-fold", type=float, default=0.5, help="lower depth bound = fold x median [0.5]")
    p.add_argument("--max-fold", type=float, default=2.0, help="upper depth bound = fold x median [2.0]")
    p.add_argument("--abs-min-depth", type=int, default=10, help="lower bound never below this [10]")
    p.add_argument("--fit-low-q", type=float, default=0.005,
                   help="fitted bounds: lower bound at this quantile of the fitted coverage distribution [0.005]")
    p.add_argument("--fit-high-q", type=float, default=0.995,
                   help="fitted bounds: upper bound at this quantile of the fitted coverage distribution [0.995]")
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
    p.add_argument("--bound-method", choices=["fold", "fitted"], default="fold",
                   help="fold: fold x median; fitted: quantiles of a negative binomial fitted to each pool's "
                        "coverage (bounds follow the pool's actual spread) [fold]")
    p.add_argument("--mode", choices=["shared", "per-sample"], default="shared",
                   help="shared: one mask, every pool called at the same sites; "
                        "per-sample: one mask per pool (more sites, uneven missingness) [shared]")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_depthmask)

    p = sub.add_parser("variance", help="kind of depth variation: spread, scale, sharing between pools")
    p.add_argument("--depth", nargs=2, action="append", required=True, metavar=("POOL", "MQX_BINS"),
                   help="mosdepth 1 kb bins at the MAPQ cut-off; repeat per pool")
    p.add_argument("--dist", nargs=2, action="append", metavar=("POOL", "REGION_DIST"),
                   help="mosdepth region.dist over the accessible genome (per-base spread)")
    p.add_argument("--mask-summary", help="<out>.mask_summary.tsv (accessible bp per chromosome)")
    p.add_argument("--chrom-class", required=True)
    p.add_argument("--thresholds", required=True, help="<out>.depth_thresholds.tsv (lower bounds)")
    p.add_argument("--mask-bed", help="accessibility mask; bins mostly masked are ignored")
    p.add_argument("--block-bins", type=int, default=5, help="run of low bins counted as a block [5]")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_variance)

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
    p.add_argument("--depth-mq0", nargs=2, action="append", metavar=("POOL", "MQ0_BINS"),
                   help="mosdepth bins without the MAPQ cut-off, to find well-mapped bins; repeat per pool")
    p.add_argument("--min-mapped-ratio", type=float, default=0.9,
                   help="well-mapped bin: depth at the MAPQ cut-off / depth at MAPQ 0 at least this [0.9]")
    p.add_argument("--fai", required=True)
    p.add_argument("--mask-bed", help="accessibility mask (bins mostly masked are ignored)")
    p.add_argument("--chrom-regex", default=".", help="optionally, only contigs matching this are analysed [all]")
    p.add_argument("--x-chroms", default="", help="sex chromosome (X or Z) contigs, comma separated")
    p.add_argument("--mito-chroms", default="", help="mitochondrial contigs, comma separated")
    p.add_argument("--pool-sizes", help="haploid genomes per pool: one number, or POOL=N,POOL=N")
    p.add_argument("--male-fraction", help="known fraction of males per pool: one number, or POOL=m,POOL=m "
                   "(default: estimated from X depth)")
    p.add_argument("--min-difference", type=float, default=0.05,
                   help="smallest |ratio - 1| treated as a real difference [0.05]")
    p.add_argument("--max-spread-ratio", type=float, default=1.5,
                   help="X/autosome depth ratio not used when their depth spreads (robust CV) differ more [1.5]")
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
    p.add_argument("--x-chroms", default="", help="sex chromosome contigs, comma separated")
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

    p = sub.add_parser("combine", help="combine the best accepted grid settings into one run")
    p.add_argument("--runs", required=True, help="grid/runs.tsv: run, setting and its values")
    p.add_argument("--assess", required=True, help="<out>.assess_runs.tsv")
    p.set_defaults(func=cmd_combine)

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
    p.add_argument("--strict", help="reference run (default: the run keeping the least data)")
    p.add_argument("--check-run", action="append", help="run evaluated but never recommended (e.g. no mask)")
    p.add_argument("--tolerance", type=float, default=0.05,
                   help="largest increase, on any statistic and covariate, over the reference run [0.05]")
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
