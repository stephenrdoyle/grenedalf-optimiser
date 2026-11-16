#!/usr/bin/env python3
"""
Synthetic stand-ins for the external tools' outputs (GenMap, mosdepth, samtools stats,
grenedalf), so grenedalf_optimiser.py can be tested without BAMs. Used by run_test.sh.

  simulate.py base DIR       reference, GenMap bedGraph, mosdepth bins, samtools stats,
                             and per-base depth (kept for the later phases)
  simulate.py dist DIR       mosdepth region.dist over DIR/optimiser.accessible.bed
  simulate.py quantize DIR   mosdepth --quantize output at the cut points in DIR/optimiser.quantize.tsv
  simulate.py grenedalf DIR  grenedalf diversity and fst tables for two runs: "clean" (statistics
                             independent of depth) and "artefact" (theta pi rises with depth)

Built in: an N gap and a soft-masked block (chr1), a non-unique block (chr1), a collapsed
repeat at 4x depth in every pool (chr1), an ambiguous-mapping block where most reads
have MAPQ 0 (chr2), a deletion in pool B only (chr2), X at 0.6-0.8x autosome depth (mixed-sex
pools; X_RATIO), and a mitochondrion
that should be excluded.
"""

import gzip
import sys

import numpy as np
import pandas as pd

CHROMS = {"sim_chr1": 300_000, "sim_chr2": 200_000, "sim_chrX": 150_000, "mito": 15_000}
POOLS = {"poolA_L3": 40, "poolB_L3": 80, "poolC_L3": 120, "poolD_L3": 60, "poolE_L3": 100, "poolF_L3": 140}
# X depth relative to autosomes = 1 - (male fraction)/2: mixed-sex pools, mostly near 50:50
X_RATIO = {"poolA_L3": 0.75, "poolB_L3": 0.75, "poolC_L3": 0.7, "poolD_L3": 0.8, "poolE_L3": 0.75, "poolF_L3": 0.6}
LD_SCALE = 2000  # bp; correlation of pool frequency deviations decays as exp(-d / LD_SCALE)
BIN = 1000
rng = np.random.default_rng(1)


def smooth_noise(n, sd, step=1000, span=5):
    """Smooth random log-scale variation along a chromosome."""
    k = n // step + 1
    x = np.convolve(rng.normal(0, sd, k + span), np.ones(span) / np.sqrt(span), "same")[:k]
    return np.repeat(x - x.mean(), step)[:n]  # centred: no chromosome-wide offset


def depth_arrays():
    """Per pool and chromosome: depth at MAPQ 0 and at MAPQ >= 30."""
    out = {}
    # Regional depth variation (GC-like), mostly shared by all pools, in 1 kb steps
    shared = {c: smooth_noise(n, 0.25) for c, n in CHROMS.items()}
    for p, mean in POOLS.items():
        for c, n in CHROMS.items():
            mu = np.full(n, mean * X_RATIO[p] if c == "sim_chrX" else mean, float)
            mu *= np.exp(shared[c] + smooth_noise(n, 0.08))
            if c == "mito":
                mu[:] = mean * 50
            if c == "sim_chr1":
                mu[20_000:22_000] = 0           # N gap
                mu[100_000:110_000] *= 4        # collapsed repeat
                mu[150_000:160_000] *= 1.0      # non-unique block: reads present, MAPQ low
            if c == "sim_chr2" and p == "poolB_L3":
                mu[150_000:155_000] = 0         # deletion in one pool
            r = 20.0                            # negative binomial overdispersion
            d0 = rng.negative_binomial(r, r / (r + np.maximum(mu, 1e-9))) * (mu > 0)
            frac30 = np.full(n, 0.97)
            if c == "sim_chr1":
                frac30[150_000:160_000] = 0.05
            if c == "sim_chr2":
                frac30[50_000:60_000] = 0.2     # ambiguous mapping not seen by GenMap
            d30 = rng.binomial(d0, frac30)
            out[(p, c)] = (d0.astype(np.int32), d30.astype(np.int32))
    return out


def base(d):
    with open(f"{d}/ref.fa", "w") as fa, open(f"{d}/ref.fa.fai", "w") as fai:
        for c, n in CHROMS.items():
            s = rng.choice(list("ACGT"), n)
            if c == "sim_chr1":
                s[20_000:22_000] = "N"
                s[60_000:63_000] = np.char.lower(s[60_000:63_000].astype(str))
            seq = "".join(s)
            fa.write(f">{c} synthetic\n")
            for i in range(0, n, 60):
                fa.write(seq[i:i + 60] + "\n")
            fai.write(f"{c}\t{n}\t0\t60\t61\n")
    with open(f"{d}/genmap.bedgraph", "w") as fh:
        for c, n in CHROMS.items():
            if c == "sim_chr1":
                fh.write(f"{c} synthetic\t0\t150000\t1\n{c} synthetic\t150000\t160000\t0.5\n"
                         f"{c} synthetic\t160000\t{n - 149}\t1\n")
            else:
                fh.write(f"{c} synthetic\t0\t{n - 149}\t1\n")
    depth = depth_arrays()
    np.savez_compressed(f"{d}/depth.npz", **{f"{p}|{c}|{q}": a for (p, c), pair in depth.items()
                                             for q, a in zip(("mq0", "mq30"), pair)})
    for p in POOLS:
        for qi, q in enumerate(("mq0", "mq30")):
            rows = []
            for c, n in CHROMS.items():
                a = depth[(p, c)][qi]
                for s in range(0, n, BIN):
                    rows.append((c, s, min(s + BIN, n), round(float(a[s:s + BIN].mean()), 2)))
            pd.DataFrame(rows).to_csv(f"{d}/{p}.{q}.regions.bed.gz", sep="\t", header=False, index=False)
        with open(f"{d}/{p}.stats", "w") as fh:
            fh.write("SN\terror rate:\t1.2e-02\t# mismatches / bases mapped (cigar)\n")
            for q, v in ((0, 500), (5, 200), (30, 1000), (60, 20000)):
                fh.write(f"MAPQ\t{q}\t{v}\n")


def load_depth(d):
    z = np.load(f"{d}/depth.npz")
    return {tuple(k.split("|")): z[k] for k in z.files}


def accessible(d):
    return pd.read_csv(f"{d}/optimiser.accessible.bed", sep="\t", header=None, names=["chrom", "start", "end"])


def dist(d):
    depth = load_depth(d)
    acc = accessible(d)
    for p in POOLS:
        lines, total = [], None
        for c, g in acc.groupby("chrom", sort=False):
            a = np.concatenate([depth[(p, c, "mq30")][s:e] for s, e in zip(g.start, g.end)])
            counts = np.bincount(a)
            total = counts if total is None else np.pad(total, (0, max(0, counts.size - total.size))) + \
                np.pad(counts, (0, max(0, total.size - counts.size)))
            cum = np.cumsum(counts[::-1])[::-1] / counts.sum()
            lines += [f"{c}\t{i}\t{cum[i]:.2f}" for i in range(cum.size - 1, -1, -1)]
        cum = np.cumsum(total[::-1])[::-1] / total.sum()
        lines += [f"total\t{i}\t{cum[i]:.2f}" for i in range(cum.size - 1, -1, -1)]
        with open(f"{d}/{p}.mosdepth.region.dist.txt", "w") as fh:
            fh.write("\n".join(lines) + "\n")


def quantize(d):
    depth = load_depth(d)
    spec = dict(pd.read_csv(f"{d}/optimiser.quantize.tsv", sep="\t", header=None).itertuples(index=False))
    for p in POOLS:
        cuts = [int(x) for x in spec[p].strip(":").split(":")]
        edges = np.array(cuts + [np.iinfo(np.int32).max])
        labels = [f"{lo}:{hi}" for lo, hi in zip(cuts, cuts[1:])] + [f"{cuts[-1]}:inf"]
        with gzip.open(f"{d}/{p}.quantized.bed.gz", "wt") as fh:
            for c, n in CHROMS.items():
                k = np.searchsorted(edges, depth[(p, c, "mq30")], side="right") - 1
                change = np.flatnonzero(k[1:] != k[:-1]) + 1
                starts = np.concatenate(([0], change))
                ends = np.concatenate((change, [n]))
                for s, e in zip(starts, ends):
                    fh.write(f"{c}\t{s}\t{e}\t{labels[k[s]]}\n")


def grenedalf(d):
    depth = load_depth(d)
    mask = pd.read_csv(f"{d}/optimiser.total_mask.bed", sep="\t", header=None, names=["chrom", "start", "end"])
    width, stride = 5000, 2500
    pools = list(POOLS)
    for run in ("clean", "artefact"):
        r = np.random.default_rng(7)
        div, fst = [], []
        for c in ("sim_chr1", "sim_chr2"):
            n = CHROMS[c]
            m = np.zeros(n, bool)
            for s, e in mask.loc[mask.chrom == c, ["start", "end"]].itertuples(index=False):
                m[s:e] = True
            for s in range(0, n - width + 1, stride):
                masked = int(m[s:s + width].sum())
                passed = int(r.poisson(0.02 * (width - masked)))
                inv = width - masked - passed - int(r.integers(0, 200))
                base_row = {"chrom": c, "start": s + 1, "end": s + width, "total.masked": masked,
                            "total.missing": 0, "total.empty": width - masked - passed - inv,
                            "total.numeric": 0, "total.invariant": inv, "total.passed": passed}
                row, frow = dict(base_row), dict(base_row)
                nd = {}
                for p in pools:
                    a = depth[(p, c, "mq30")][s:s + width]
                    nd[p] = a.mean() / POOLS[p]
                    pi = 0.01 * r.lognormal(0, 0.2)
                    if run == "artefact":
                        pi *= nd[p] ** 1.5
                    row[f"{p}.theta_pi"] = pi
                    row[f"{p}.theta_watterson"] = pi * 1.1
                    row[f"{p}.passed"] = passed
                    row[f"{p}.tajimas_d"] = r.normal(0, 0.5)
                for i, a in enumerate(pools):
                    for b in pools[i + 1:]:
                        f = abs(r.normal(0.02, 0.01))
                        if run == "artefact":
                            f += 0.05 * abs(np.log2(max(nd[a], 1e-3) / max(nd[b], 1e-3)))
                        frow[f"{a}:{b}.missing"] = 0
                        frow[f"{a}:{b}.numeric"] = 0
                        frow[f"{a}:{b}.passed"] = passed
                        frow[f"{a}:{b}.fst"] = f
                div.append(row)
                fst.append(frow)
        pd.DataFrame(div).to_csv(f"{d}/{run}.diversity.csv", sep="\t", index=False)
        pd.DataFrame(fst).to_csv(f"{d}/{run}.fst.csv", sep="\t", index=False)
    frequency(d, depth)


def frequency(d, depth):
    """grenedalf frequency table with linked SNPs: each pool's deviation from the
    mean frequency follows an Ornstein-Uhlenbeck process along the chromosome."""
    r = np.random.default_rng(11)
    pools = list(POOLS)
    rows = []
    for c in ("sim_chr1", "sim_chr2", "sim_chrX"):
        n = CHROMS[c]
        pos = np.sort(r.choice(n, size=int(0.02 * n), replace=False))
        rho = np.exp(-np.diff(pos, prepend=pos[0]) / LD_SCALE)
        z = np.zeros((pos.size, len(pools)))
        z[0] = r.normal(size=len(pools))
        for i in range(1, pos.size):
            z[i] = rho[i] * z[i - 1] + np.sqrt(1 - rho[i] ** 2) * r.normal(size=len(pools))
        logit = r.normal(0, 1.5, pos.size)[:, None] + 1.5 * z
        freq = 1 / (1 + np.exp(-logit))
        cols = {"chrom": c, "pos": pos + 1, "ref": "A", "alt": "C"}
        for j, p in enumerate(pools):
            dp = depth[(p, c, "mq30")][pos]
            alt = r.binomial(dp, freq[:, j])
            cols[f"{p}.REF_CNT"] = dp - alt
            cols[f"{p}.ALT_CNT"] = alt
        rows.append(pd.DataFrame(cols))
    pd.concat(rows).to_csv(f"{d}/frequency.csv", sep="\t", index=False)


if __name__ == "__main__":
    {"base": base, "dist": dist, "quantize": quantize, "grenedalf": grenedalf}[sys.argv[1]](sys.argv[2])
