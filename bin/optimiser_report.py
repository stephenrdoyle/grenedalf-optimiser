"""
HTML report for grenedalf_optimiser: one self-contained file (inline SVG, no external
scripts or fonts), so it opens offline from a cluster. Each section is drawn from the
step outputs that exist; missing steps are skipped.
"""

import html
import math
import os
import re

import numpy as np
import pandas as pd

# Colours: categorical slots 1-2, sequential blue ramp, status; light and dark (CSS vars)
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
       "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

CSS = """
.viz-root { color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --s1: #2a78d6; --s2: #eb6834; --s1-wash: rgba(42,120,214,0.10); --neutral: #b5b4ad;
  --good: #0ca30c; --critical: #d03b3b; --warning: #fab219; }
@media (prefers-color-scheme: dark) { :root:where(:not([data-theme="light"])) .viz-root { color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --s1: #3987e5; --s2: #d95926; --s1-wash: rgba(57,135,229,0.12); --neutral: #5f5e5a; } }
:root[data-theme="dark"] .viz-root { color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --s1: #3987e5; --s2: #d95926; --s1-wash: rgba(57,135,229,0.12); --neutral: #5f5e5a; }
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); }
.viz-root { font-family: system-ui, -apple-system, "Segoe UI", sans-serif; color: var(--ink);
  background: var(--page); line-height: 1.45; padding: 24px 16px 64px; }
main { max-width: 1040px; margin: 0 auto; }
h1 { font-size: 26px; margin: 0 0 4px; } h2 { font-size: 19px; margin: 0 0 4px; }
h3 { font-size: 14px; margin: 22px 0 6px; color: var(--ink2); font-weight: 600; }
p, li { color: var(--ink2); font-size: 14px; max-width: 75ch; }
.sub { color: var(--muted); font-size: 13px; margin: 0 0 20px; }
section { background: var(--surface); border: 1px solid var(--border); border-radius: 10px;
  padding: 20px; margin: 0 0 20px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; margin: 14px 0 6px; }
.tile { border: 1px solid var(--border); border-radius: 8px; padding: 12px 14px; }
.tile .label { font-size: 12px; color: var(--ink2); }
.tile .value { font-size: 24px; font-weight: 600; margin-top: 2px; }
.tile .note { font-size: 12px; color: var(--muted); }
.grid2 { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 16px; }
svg { display: block; width: 100%; height: auto; overflow: visible; }
svg text { font-family: inherit; font-variant-numeric: tabular-nums; }
.tick { fill: var(--muted); font-size: 11px; } .alabel { fill: var(--ink2); font-size: 12px; }
.ref { fill: var(--ink2); font-size: 11px; }
.gl { stroke: var(--grid); stroke-width: 1; } .ax { stroke: var(--axis); stroke-width: 1; }
.legend { display: flex; flex-wrap: wrap; gap: 14px; font-size: 12px; color: var(--ink2); margin: 4px 0 6px; }
.key { display: inline-block; width: 14px; height: 3px; border-radius: 2px; vertical-align: middle; margin-right: 6px; }
.dotkey { display: inline-block; width: 9px; height: 9px; border-radius: 50%; vertical-align: middle; margin-right: 6px; }
.tablewrap { overflow-x: auto; margin-top: 8px; }
table { border-collapse: collapse; font-size: 12px; font-variant-numeric: tabular-nums; }
th, td { padding: 4px 8px; text-align: right; border-bottom: 1px solid var(--grid); white-space: nowrap; }
th { color: var(--ink2); font-weight: 600; } td:first-child, th:first-child { text-align: left; }
tr.pick td { font-weight: 600; } td.cellpick { outline: 2px solid var(--ink); outline-offset: -2px; }
details { margin-top: 10px; } summary { cursor: pointer; font-size: 13px; color: var(--ink2); }
.status { display: inline-flex; align-items: center; gap: 6px; font-size: 13px; font-weight: 600; }
.status svg { width: 14px; height: 14px; display: inline; }
pre { font-size: 12px; background: var(--page); border: 1px solid var(--border); border-radius: 6px;
  padding: 10px; overflow-x: auto; color: var(--ink); }
.cmd { position: relative; } .cmd pre { max-height: 360px; padding-top: 34px; }
.copy { position: absolute; top: 8px; right: 8px; font: inherit; font-size: 12px; color: var(--ink2);
  background: var(--surface); border: 1px solid var(--border); border-radius: 6px; padding: 3px 10px; cursor: pointer; }
.theme { float: right; font: inherit; font-size: 12px; color: var(--ink2); background: var(--surface);
  border: 1px solid var(--border); border-radius: 6px; padding: 4px 10px; cursor: pointer; }
"""

THEME_JS = """
<script>
(function(){var r=document.documentElement,b=document.getElementById('theme');
try{var t=localStorage.getItem('optimiser-theme');if(t)r.setAttribute('data-theme',t);}catch(e){}
document.querySelectorAll('.copy').forEach(function(c){c.addEventListener('click',function(){
var t=document.getElementById(c.dataset.target).textContent;
(navigator.clipboard?navigator.clipboard.writeText(t):Promise.reject()).then(function(){c.textContent='Copied';
setTimeout(function(){c.textContent='Copy';},1500);},function(){c.textContent='Select and copy';});});});
b.addEventListener('click',function(){var d=r.getAttribute('data-theme')==='dark'||
(!r.getAttribute('data-theme')&&matchMedia('(prefers-color-scheme: dark)').matches);
var n=d?'light':'dark';r.setAttribute('data-theme',n);try{localStorage.setItem('optimiser-theme',n);}catch(e){}});})();
</script>
"""

ICON_OK = ('<svg viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="var(--good)"/>'
           '<path d="M4.5 8.2l2.3 2.3 4.7-4.9" stroke="#fff" stroke-width="2" fill="none"/></svg>')
ICON_BAD = ('<svg viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="var(--critical)"/>'
            '<path d="M5.2 5.2l5.6 5.6M10.8 5.2l-5.6 5.6" stroke="#fff" stroke-width="2"/></svg>')


ICON_INFO = ('<svg viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="var(--s1)"/>'
             '<path d="M8 7v4.5M8 4.6v.1" stroke="#fff" stroke-width="2" stroke-linecap="round"/></svg>')


def esc(x):
    return html.escape(str(x))


def compact(v):
    v = float(v)
    for div, suf in ((1e9, "G"), (1e6, "M"), (1e3, "K")):
        if abs(v) >= div:
            s = f"{v / div:.1f}".rstrip("0").rstrip(".")
            return s + suf
    return f"{v:g}"


def nice_ticks(lo, hi, n=5):
    if hi <= lo:
        hi = lo + 1
    step = (hi - lo) / n
    mag = 10 ** math.floor(math.log10(step))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= step), default=mag * 10)
    start = math.floor(lo / step) * step
    ticks = []
    t = start
    while t <= hi + step * 1e-9:
        if t >= lo - step * 1e-9:
            ticks.append(round(t, 10))
        t += step
    return ticks


def log_ticks(lo, hi):
    out = []
    for e in range(math.floor(math.log10(lo)), math.ceil(math.log10(hi)) + 1):
        for m in (1, 2, 5):
            v = m * 10 ** e
            if lo <= v <= hi:
                out.append(v)
    return out


class Plot:
    """A single-axis SVG plot. Coordinates in data space; x and/or y may be log10."""

    def __init__(self, xlim, ylim, xlog=False, ylog=False, w=480, h=260, ml=52, mr=16, mt=14, mb=40,
                 xlabel="", ylabel="", xfmt=compact, yfmt=compact):
        self.xlim, self.ylim, self.xlog, self.ylog = xlim, ylim, xlog, ylog
        self.w, self.h, self.ml, self.mr, self.mt, self.mb = w, h, ml, mr, mt, mb
        self.xlabel, self.ylabel, self.xfmt, self.yfmt = xlabel, ylabel, xfmt, yfmt
        self.parts = []

    def _t(self, v, lim, log):
        if log:
            v, lim = math.log10(max(v, 1e-12)), (math.log10(lim[0]), math.log10(lim[1]))
        return (v - lim[0]) / ((lim[1] - lim[0]) or 1)

    def sx(self, x):
        return self.ml + self._t(x, self.xlim, self.xlog) * (self.w - self.ml - self.mr)

    def sy(self, y):
        return self.h - self.mb - self._t(y, self.ylim, self.ylog) * (self.h - self.mt - self.mb)

    def frame(self):
        p = []
        yt = log_ticks(*self.ylim) if self.ylog else nice_ticks(*self.ylim)
        xt = log_ticks(*self.xlim) if self.xlog else nice_ticks(*self.xlim)
        for y in yt:
            Y = self.sy(y)
            p.append(f'<line class="gl" x1="{self.ml}" x2="{self.w - self.mr}" y1="{Y:.1f}" y2="{Y:.1f}"/>')
            p.append(f'<text class="tick" x="{self.ml - 6}" y="{Y + 4:.1f}" text-anchor="end">{esc(self.yfmt(y))}</text>')
        base = self.h - self.mb
        p.append(f'<line class="ax" x1="{self.ml}" x2="{self.w - self.mr}" y1="{base}" y2="{base}"/>')
        for x in xt:
            X = self.sx(x)
            p.append(f'<text class="tick" x="{X:.1f}" y="{base + 16}" text-anchor="middle">{esc(self.xfmt(x))}</text>')
        p.append(f'<text class="alabel" x="{(self.ml + self.w - self.mr) / 2:.1f}" y="{self.h - 6}" '
                 f'text-anchor="middle">{esc(self.xlabel)}</text>')
        p.append(f'<text class="alabel" transform="translate(12,{(self.mt + base) / 2:.1f}) rotate(-90)" '
                 f'text-anchor="middle">{esc(self.ylabel)}</text>')
        self.parts[:0] = p

    def band(self, x0, x1, label=""):
        X0, X1 = self.sx(max(x0, self.xlim[0])), self.sx(min(x1, self.xlim[1]))
        self.parts.append(f'<rect x="{X0:.1f}" y="{self.mt}" width="{max(X1 - X0, 0):.1f}" '
                          f'height="{self.h - self.mt - self.mb}" fill="var(--s1-wash)"><title>{esc(label)}</title></rect>')

    def vline(self, x, label, anchor="start", dy=0):
        if not (self.xlim[0] <= x <= self.xlim[1]):
            return
        X = self.sx(x)
        self.parts.append(f'<line x1="{X:.1f}" x2="{X:.1f}" y1="{self.mt}" y2="{self.h - self.mb}" '
                          f'stroke="var(--ink2)" stroke-width="1"/>')
        off = 4 if anchor == "start" else -4
        self.parts.append(f'<text class="ref" x="{X + off:.1f}" y="{self.mt - 10 - dy}" '
                          f'text-anchor="{anchor}">{esc(label)}</text>')

    def hline(self, y, label, anchor="end"):
        Y = self.sy(y)
        self.parts.append(f'<line x1="{self.ml}" x2="{self.w - self.mr}" y1="{Y:.1f}" y2="{Y:.1f}" '
                          f'stroke="var(--ink2)" stroke-width="1"/>')
        X = self.w - self.mr if anchor == "end" else self.ml + 6
        self.parts.append(f'<text class="ref" x="{X}" y="{Y - 5:.1f}" text-anchor="{anchor}">{esc(label)}</text>')

    def line(self, xs, ys, color, name, tips=None, dots=False):
        pts = [(self.sx(x), self.sy(y)) for x, y in zip(xs, ys) if np.isfinite(x) and np.isfinite(y)]
        if len(pts) > 1:
            d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)
            self.parts.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="2" '
                              f'stroke-linejoin="round" stroke-linecap="round"/>')
        for i, (x, y) in enumerate(pts):
            tip = esc(tips[i]) if tips else ""
            if dots:
                self.parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{color}" '
                                  f'stroke="var(--surface)" stroke-width="2"/>')
            if tip:  # hit target larger than the mark
                self.parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="9" fill="transparent"><title>{tip}</title></circle>')

    def dot(self, x, y, color, tip, label=None, r=5):
        X, Y = self.sx(x), self.sy(y)
        self.parts.append(f'<circle cx="{X:.1f}" cy="{Y:.1f}" r="{r}" fill="{color}" stroke="var(--surface)" stroke-width="2"/>')
        self.parts.append(f'<circle cx="{X:.1f}" cy="{Y:.1f}" r="{r + 6}" fill="transparent"><title>{esc(tip)}</title></circle>')
        if label:
            self.parts.append(f'<text class="ref" x="{X + r + 4:.1f}" y="{Y - r - 4:.1f}">{esc(label)}</text>')

    def svg(self, title):
        self.frame()
        return (f'<svg viewBox="0 0 {self.w} {self.h}" role="img" aria-label="{esc(title)}">'
                + "".join(self.parts) + "</svg>")


def hbars(labels, values, fmt, tips, title, w=480, highlight=None, colors=None, ref=None, xmax=None, ml=110):
    """Horizontal bars from one baseline, value at the tip; optional vertical reference line."""
    bar, gap = 18, 10
    h = len(labels) * (bar + gap) + (30 if ref else 24)
    top = 18 if ref else 6
    vmax = xmax or max(max(values), 1e-12)
    W = w - ml - 70
    p = [f'<line class="ax" x1="{ml}" x2="{ml}" y1="{top - 2}" y2="{h - 20}"/>']
    for i, (lab, v, tip) in enumerate(zip(labels, values, tips)):
        y = top + i * (bar + gap)
        bw = max(W * v / vmax, 1)
        r = min(4, bw / 2)
        col = colors[i] if colors else ("var(--s2)" if highlight and lab == highlight else "var(--s1)")
        # square at the baseline, rounded at the data end
        p.append(f'<path d="M{ml},{y} h{bw - r:.1f} a{r},{r} 0 0 1 {r},{r} v{bar - 2 * r} '
                 f'a{r},{r} 0 0 1 {-r},{r} h{-(bw - r):.1f} z" fill="{col}"><title>{esc(tip)}</title></path>')
        p.append(f'<text class="tick" x="{ml - 8}" y="{y + 13}" text-anchor="end">{esc(lab)}</text>')
        p.append(f'<text class="ref" x="{ml + bw + 6:.1f}" y="{y + 13}">{esc(fmt(v))}</text>')
    if ref:
        x = ml + W * ref[0] / vmax
        p.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{top - 4}" y2="{h - 20}" stroke="var(--ink2)" stroke-width="1"/>')
        p.append(f'<text class="ref" x="{x + 4:.1f}" y="{top - 6}">{esc(ref[1])}</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{esc(title)}">' + "".join(p) + "</svg>"


def seq_color(v, vmax):
    if not np.isfinite(v) or vmax <= 0:
        return "transparent", "var(--ink)"
    i = int(round(min(max(v / vmax, 0), 1) * (len(SEQ) - 1)))
    return SEQ[i], ("#ffffff" if i >= 6 else "#0b0b0b")


def heat_table(df, value_cols, fmt, vmax=None, pick=None, pick_row=None, caption="", scale=None):
    scale = scale or (lambda v: v)
    vals = np.vectorize(scale, otypes=[float])(df[value_cols].to_numpy(dtype=float)) if value_cols else np.zeros((0, 0))
    vmax = vmax if vmax is not None else (np.nanmax(vals) if np.isfinite(vals).any() else 1)
    out = ['<div class="tablewrap"><table>']
    if caption:
        out.append(f'<caption style="caption-side:bottom;text-align:left;color:var(--muted);'
                   f'font-size:12px;padding-top:6px">{esc(caption)}</caption>')
    out.append("<tr>" + "".join(f"<th>{esc(c)}</th>" for c in df.columns) + "</tr>")
    for ri, (_, row) in enumerate(df.iterrows()):
        cls = ' class="pick"' if pick_row is not None and pick_row(row) else ""
        cells = []
        for c in df.columns:
            v = row[c]
            if c in value_cols:
                bg, fg = seq_color(scale(float(v)) if pd.notna(v) else np.nan, vmax)
                extra = ' class="cellpick"' if pick and pick(row, c) else ""
                cells.append(f'<td{extra} style="background:{bg};color:{fg}">{esc(fmt(v)) if pd.notna(v) else ""}</td>')
            else:
                cells.append(f"<td>{esc(v)}</td>")
        out.append(f"<tr{cls}>" + "".join(cells) + "</tr>")
    out.append("</table></div>")
    return "".join(out)


def plain_table(df, max_rows=None):
    d = (df if max_rows is None else df.head(max_rows)).copy()
    for c in d.columns:
        v = pd.to_numeric(d[c], errors="coerce")
        if d[c].dtype.kind == "f" and v.notna().any() and (v.dropna() % 1 == 0).all():
            d[c] = v.map(lambda x: "–" if pd.isna(x) else f"{int(x):,}" if abs(x) >= 10000 else str(int(x)))
    return heat_table(d.astype(object).where(d.notna(), "–"), [], str)


def read_tsv(path, **kw):
    return pd.read_csv(path, sep="\t", **kw) if path and os.path.exists(path) and os.path.getsize(path) else None


def read_params(path):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return dict(l.strip().split("=", 1) for l in fh if "=" in l and not l.startswith("#"))


def bp(v):
    return f"{int(round(float(v))):,} bp"


def pct(v):
    return f"{100 * float(v):.1f}%"


def tile(label, value, note=""):
    return (f'<div class="tile"><div class="label">{esc(label)}</div><div class="value">{esc(value)}</div>'
            f'<div class="note">{esc(note)}</div></div>')


# ---------------------------------------------------------------------------

TARGETS = {"diversity": "diversity (θπ, θW per pool)", "fst": "FST (per pair of pools)"}
TARGET_TYPES = {"diversity": ("theta_pi", "theta_watterson"), "fst": ("fst",)}


def chosen_runs(P):
    """Grid run chosen for diversity and for FST (older outputs: one run for both)."""
    a = read_params(f"{P}.assess.params")
    fst = a.get("RECOMMENDED_RUN_FST") or a.get("RECOMMENDED_RUN")
    return {"diversity": a.get("RECOMMENDED_RUN_DIVERSITY") or fst, "fst": fst}


def run_file(out_dir, run, name):
    path = os.path.join(out_dir, "grid", run, name) if run and out_dir else None
    return path if path and os.path.exists(path) else None


def section_summary(P, out_dir):
    params = read_params(os.path.join(out_dir, "params.sh")) if out_dir else {}
    win = read_params(f"{P}.windows.params")
    ld = read_params(f"{P}.ld.params")
    mask = read_tsv(f"{P}.mask_summary.tsv", dtype={"chrom": str})
    # The grid may choose different depth folds, mask mode or minimum count for diversity and
    # FST: report the chosen runs' masks
    runs = chosen_runs(P)
    dms = {t: read_tsv(run_file(out_dir, r, "mask.depthmask_summary.tsv")) for t, r in runs.items()}
    modes = {t: read_params(run_file(out_dir, r, "mask.depthmask.params") or "").get("MASK_MODE", "shared")
             for t, r in runs.items()}
    dm = dms["fst"] if dms["fst"] is not None else read_tsv(f"{P}.depthmask_summary.tsv")
    dmp = {}
    for t, r in runs.items():
        f = run_file(out_dir, r, "mask.depthmask.params")
        dmp[t] = read_params(f) if f and os.path.exists(f) else {}
    tiles = []
    if mask is not None:
        tot = mask[mask["chrom"] == "total"].iloc[0]
        excl = int(tot["excluded"])
        analysed = int(tot["length"]) - excl
        acc = int(tot["accessible"])
        tiles.append(tile("Accessible genome", pct(acc / analysed) if analysed else "–",
                          f"{compact(acc)} of {compact(analysed)} bp analysed"))
    if dm is not None and (dm["pool"] == "shared").any():
        s = dm[dm["pool"] == "shared"].iloc[0]
        tiles.append(tile("Callable in every pool", pct(1 - s["out_of_range_frac"]),
                          "of the accessible genome, after the FST run's depth bounds"))
    mc_d, mc_f = params.get("DIVERSITY_MIN_COUNT"), params.get("FST_MIN_COUNT")
    if mc_d and mc_f and mc_d != mc_f:
        tiles.append(tile("Minimum allele count", f"{mc_d} / {mc_f}", "diversity / FST"))
    else:
        mc = mc_f or params.get("MIN_COUNT") or read_params(f"{P}.thresholds.params").get("MIN_COUNT")
        if mc:
            tiles.append(tile("Minimum allele count", mc, "reads supporting an allele"))
    if win:
        tiles.append(tile("Window", bp(win['WINDOW']), f"stride {bp(win['STRIDE'])}"))
    if ld.get("LD_D10", "NA") != "NA":
        tiles.append(tile("Linkage scale", bp(ld['LD_D10']), f"90% decay; half at {bp(ld['LD_D50'])}"))
    if runs["fst"] and runs["diversity"] == runs["fst"]:
        tiles.append(tile("Grid choice", runs["fst"].replace("_", " "), "for diversity and FST"))
    else:
        for t, lab in (("diversity", "Diversity settings"), ("fst", "FST settings")):
            if runs[t]:
                tiles.append(tile(lab, runs[t].replace("_", " "), "grid run, chosen for these statistics"))
    body = f'<div class="tiles">{"".join(tiles)}</div>'

    thr = read_tsv(f"{P}.depth_thresholds.tsv")
    if thr is not None:
        rows = []
        for pool, g in thr.groupby("pool", sort=False):
            g = g.set_index("class")
            r = {"pool": pool}
            if "pool_size" in g:
                r["pool size"] = int(g["pool_size"].iloc[0])
            def kept(d, c, t="fst"):
                ch = d.set_index("pool").loc[pool] if d is not None and pool in set(d["pool"]) else None
                lo, hi = g.loc[c, "lower"], g.loc[c, "upper"]
                if ch is not None and pd.notna(ch.get(f"lower_{c}", np.nan)):
                    lo, hi = int(ch[f"lower_{c}"]), int(ch[f"upper_{c}"])
                rel = (f", not below {dmp[t].get('MIN_FOLD')}× the other pools"
                       if dmp.get(t, {}).get("BOUND_METHOD") == "relative" else "")
                return f"{lo}–{hi}×{rel}"
            for c, lab in (("A", "autosomes"), ("X", "X")):
                if c in g.index:
                    r[f"{lab}: median"] = f"{g.loc[c, 'median']}×"
                    r[f"{lab}: kept"] = kept(dm, c)
                    if dms["diversity"] is not None and kept(dms["diversity"], c, "diversity") != r[f"{lab}: kept"]:
                        r[f"{lab}: kept"] += " (FST)"
                        r[f"{lab}: kept (diversity)"] = kept(dms["diversity"], c, "diversity")
            if "n_eff_at_median" in g and "A" in g.index:
                r["n_eff (autosomes)"] = g.loc["A", "n_eff_at_median"]
            rows.append(r)
        mp = read_params(os.path.join(out_dir, "params.sh")) if out_dir else {}
        lo_key = "FST_MIN_DEPTH" if "FST_MIN_DEPTH" in mp else "MIN_DEPTH"
        hi_key = lo_key.replace("MIN", "MAX")
        glob = (f" grenedalf's own <code>--filter-sample-min-read-depth</code> / <code>--filter-sample-max-read-depth</code> "
                f"are set to {mp[lo_key]} / {mp[hi_key]}, wide enough not to undo these."
                if lo_key in mp else "")
        cost = ""
        if dm is not None and (dm["pool"] == "shared").any():
            per = dm[dm["pool"] != "shared"]["out_of_range_frac"].astype(float)
            sh = float(dm.loc[dm["pool"] == "shared", "out_of_range_frac"].iloc[0])
            cost = (f" Each pool's own bounds remove {100 * per.min():.0f}–{100 * per.max():.0f}% of the accessible "
                    f"genome; the shared mask removes <strong>{100 * sh:.0f}%</strong> ({sh / per.max():.1f}× the "
                    "largest single pool), because the pools are out of range in different places. Most of that "
                    "is usually low depth; see <em>Where the mask falls</em> for whether it is long blocks in "
                    "single pools (divergent haplotypes, deletions) or scattered.")
        mode_txt = {"shared": "a site outside any pool's range is masked in every pool, so all pools are called "
                              "at the same sites",
                    "per-sample": "a site outside a pool's range is masked in that pool only"}
        how = "; ".join(f"for {t if t == 'diversity' else 'FST'}, {mode_txt.get(m, m)}" for t, m in modes.items()) \
            if len(set(modes.values())) > 1 else mode_txt.get(modes["fst"], "")
        body += ("<h3>Depth bounds per pool</h3>"
                 "<p>grenedalf has a single depth filter for all pools, so each pool's own range is applied as a "
                 f"mask: {how}.{glob} n_eff is the number of chromosomes effectively sampled at that depth from a pool of "
                 f"that size; depth well beyond the pool size adds little.{cost}</p>" + plain_table(pd.DataFrame(rows)))

    cmd_path = os.path.join(out_dir, "grenedalf_commands.sh") if out_dir else ""
    if cmd_path and os.path.exists(cmd_path):
        body += ("<h3>Recommended grenedalf commands</h3>"
                 "<p>Also saved as <code>grenedalf_commands.sh</code>. One run per chromosome class (see "
                 "<em>X chromosome and mitochondrion</em> below).</p>"
                 '<div class="cmd"><button class="copy" type="button" data-target="cmdtext">Copy</button>'
                 f'<pre id="cmdtext">{esc(open(cmd_path).read())}</pre></div>')
    if params:
        body += ("<details><summary>params.sh</summary><pre>"
                 + esc(open(os.path.join(out_dir, "params.sh")).read()) + "</pre></details>")
    return ("<section><h2>Recommended settings</h2><p>Chosen from technical properties of the data only: "
            "no population labels, phenotypes or peak positions were used.</p>" + body + "</section>")


def section_mask(P):
    m = read_tsv(f"{P}.mask_summary.tsv", dtype={"chrom": str})
    if m is None:
        return ""
    reasons = [r for r in ("gap", "softmask", "low_mappability", "collapsed", "ambiguous", "low_depth")
               if f"{r}_frac" in m and m[r].sum() > 0]
    kept = m[(m["excluded"] < m["length"]) | (m["chrom"] == "total")].copy()
    shown = kept[["chrom"] + [f"{r}_frac" for r in reasons] + ["masked_frac"]].copy()
    shown.columns = ["chromosome"] + [r.replace("_", " ") for r in reasons] + ["masked (any)"]
    excluded = m[(m["excluded"] == m["length"]) & (m["chrom"] != "total")]["chrom"].tolist()
    note = (f"Excluded entirely (not analysed): {', '.join(excluded[:8])}{' …' if len(excluded) > 8 else ''}."
            if excluded else "")
    specific = read_tsv(f"{P}.pool_specific_depth.tsv", dtype={"chrom": str})
    spec_html = ""
    if specific is not None and len(specific):
        counts = specific.groupby("pool").agg(bins=("start", "size"),
                                              low=("norm_depth", lambda x: int((x < 1).sum())),
                                              high=("norm_depth", lambda x: int((x >= 1).sum()))).reset_index()
        counts.columns = ["pool", "anomalous bins", "below range", "above range"]
        spec_html = ("<h3>Depth anomalies in single pools (reported, not masked)</h3>"
                     "<p>Bins where one pool's depth is outside the normal range while the consensus of all "
                     "pools is normal: candidate deletions or copy-number differences. They are left in, "
                     "because they may be biology. Full list: <code>pool_specific_depth.tsv</code>.</p>"
                     + plain_table(counts))
    return ("<section><h2>Accessibility mask</h2>"
            "<p>Fraction of each chromosome masked, by reason. A base can carry more than one reason, "
            "so rows do not sum to the total. Masked sites are reported by grenedalf as <em>masked</em>, "
            "not <em>missing</em>, and are the same for every pool.</p>"
            + heat_table(shown, list(shown.columns[1:]), pct, vmax=max(float(shown.iloc[:, 1:].max().max()), 1e-9),
                         pick_row=lambda r: r["chromosome"] == "total")
            + (f"<p class='sub'>{esc(note)}</p>" if note else "") + spec_html + "</section>")


MASK_ROWS = [("gap", "gap (N)"), ("softmask", "soft-masked"), ("low_mappability", "low mappability"),
             ("collapsed", "collapsed repeat"), ("ambiguous", "ambiguous mapping"), ("low_depth", "low depth (all pools)")]


def covered_upto(starts, ends, xs):
    """Bases covered by disjoint sorted intervals [starts, ends) to the left of each x."""
    if len(starts) == 0:
        return np.zeros(len(xs))
    cum = np.concatenate(([0], np.cumsum(ends - starts)))
    i = np.searchsorted(starts, xs, side="right")            # intervals starting before x
    prev_end = np.where(i > 0, ends[np.clip(i - 1, 0, None)], 0)
    prev_start = np.where(i > 0, starts[np.clip(i - 1, 0, None)], 0)
    partial = np.where(i > 0, np.clip(xs - prev_start, 0, prev_end - prev_start), 0)
    return cum[np.clip(i - 1, 0, None)] * (i > 0) + partial


def bin_fraction(iv, length, edges):
    """Fraction of each bin covered by the (merged) intervals of one chromosome."""
    s = iv["start"].to_numpy(); e = iv["end"].to_numpy()
    o = np.argsort(s)
    c = covered_upto(s[o], e[o], edges.astype(float))
    return np.diff(c) / np.diff(edges)


def section_karyoplot(P, out_dir, max_chroms=10, max_bins=300):
    """Where along each chromosome the mask falls: one panel per chromosome (the longest
    analysed ones), one row per mask reason, colour = fraction of each bin masked."""
    summary = read_tsv(f"{P}.mask_summary.tsv", dtype={"chrom": str})
    if summary is None or not os.path.exists(f"{P}.mask_reasons.bed"):
        return ""
    chroms = summary[(summary["chrom"] != "total") & (summary["excluded"] < summary["length"])]
    chroms = chroms.sort_values("length", ascending=False).head(max_chroms)
    if chroms.empty:
        return ""
    reasons_bed = pd.read_csv(f"{P}.mask_reasons.bed", sep="\t", header=None, names=["chrom", "start", "end", "why"],
                              dtype={"chrom": str})
    reasons_bed = reasons_bed[reasons_bed["chrom"].isin(chroms["chrom"])]
    # Final masks of the chosen grid runs: the total mask, plus each pool's own mask when the
    # run has per-pool masks (shown as the mean over pools)
    finals = {}
    for t, r in chosen_runs(P).items():
        tot = run_file(out_dir, r, "mask.total_mask.bed")
        if not tot or any(v[0] == r for v in finals.values()):
            continue
        sm = run_file(out_dir, r, "mask.sample_masks.tsv")
        pools = ([merge_bed(read_bed3(f)) for f in pd.read_csv(sm, sep="\t", header=None)[1]] if sm else None)
        finals[t] = (r, merge_bed(read_bed3(tot)), pools)
    if not finals and os.path.exists(f"{P}.total_mask.bed"):
        finals["fst"] = ("", merge_bed(read_bed3(f"{P}.total_mask.bed")), None)

    def final_label(t):
        what = "diversity and FST" if len(finals) == 1 else ("FST" if t == "fst" else "diversity")
        return f"final mask, {what}" + (" (pool mean)" if finals[t][2] is not None else "")

    def final_fraction(t, chrom, L, edges):
        _, total, pools = finals[t]
        tc = total[total["chrom"] == chrom]
        if pools is None:
            return bin_fraction(tc, L, edges)
        return np.mean([bin_fraction(merge_bed(pd.concat([tc, pm[pm["chrom"] == chrom]])), L, edges)
                        for pm in pools], axis=0)

    rows = [(k, lab) for k, lab in MASK_ROWS if (reasons_bed["why"].str.contains(k)).any()]
    rows.append(("any", "accessibility mask"))
    for t in finals:
        rows.append((f"final:{t}", final_label(t)))

    # Per pool: share of each bin (1 kb mosdepth bins) outside that pool's depth bounds
    pool_bins = {}
    thr = read_tsv(f"{P}.depth_thresholds.tsv")
    classes = read_tsv(f"{P}.chrom_class.tsv", dtype={"chrom": str})
    depth_dir = os.path.join(out_dir, "depth") if out_dir else ""
    if thr is not None and classes is not None and os.path.isdir(depth_dir):
        cls = dict(zip(classes["chrom"], classes["class"]))
        for pool, g in thr.groupby("pool", sort=False):
            files = [f for f in os.listdir(depth_dir) if f.startswith(pool + ".mq") and f.endswith(".regions.bed.gz")
                     and not f.startswith(pool + ".mq0.")]
            if not files:
                continue
            b = pd.read_csv(os.path.join(depth_dir, files[0]), sep="\t", header=None, names=["chrom", "start", "end", "d"],
                            dtype={"chrom": str})
            b = b[b["chrom"].isin(chroms["chrom"])]
            bounds = g.set_index("class")[["lower", "upper"]]
            c = b["chrom"].map(cls)
            lo = c.map(bounds["lower"]); hi = c.map(bounds["upper"])
            b["out"] = ((b["d"] < lo) | (b["d"] > hi)).astype(float)
            b["low"] = (b["d"] < lo).astype(float)
            pool_bins[pool] = b
    for pool in pool_bins:
        rows.append((f"pool:{pool}", f"{pool} out of range"))

    longest = int(chroms["length"].max())
    bin_size = max(10_000, int(math.ceil(longest / max_bins / 10_000) * 10_000))
    ml, mr, w = 168, 12, 900
    rowh, gap = 11, 1
    scale = (w - ml - mr) / longest
    panels, blocks = [], []
    block = 5_000_000 if longest > 20e6 else 1_000_000
    for r in chroms.itertuples():
        L = int(r.length)
        edges = np.append(np.arange(0, L, bin_size), L)
        sub = reasons_bed[reasons_bed["chrom"] == r.chrom]
        fr = {}
        for k, _ in rows:
            if k.startswith("pool:"):
                pb = pool_bins[k[5:]]
                pb = pb[pb["chrom"] == r.chrom]
                idx = np.minimum(pb["start"].to_numpy() // bin_size, len(edges) - 2)
                wgt = (pb["end"] - pb["start"]).to_numpy().astype(float)
                tot = np.bincount(idx, weights=wgt, minlength=len(edges) - 1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    fr[k] = np.nan_to_num(np.bincount(idx, weights=wgt * pb["out"].to_numpy(),
                                                      minlength=len(edges) - 1) / tot)
                continue
            if k.startswith("final:"):
                fr[k] = final_fraction(k[6:], r.chrom, L, edges)
                continue
            if k == "any":
                iv = sub
            else:
                iv = sub[np.array([k in w.split(",") for w in sub["why"]], dtype=bool)]
            fr[k] = bin_fraction(iv, L, edges)
        h = len(rows) * (rowh + gap) + 30 + (6 if pool_bins else 0)
        g = []
        for j, (k, lab) in enumerate(rows):
            y = 4 + j * (rowh + gap) + (6 if k.startswith("pool:") else 0)
            g.append(f'<text class="tick" x="{ml - 8}" y="{y + rowh - 2}" text-anchor="end">{esc(lab)}</text>')
            g.append(f'<rect x="{ml}" y="{y}" width="{L * scale:.1f}" height="{rowh}" fill="var(--page)"/>')
            v = fr[k]
            # one rect per run of bins with the same shade (keeps the file small)
            # shade proportional to the fraction masked; under 2% is not drawn
            steps = np.where(v >= 0.02, np.maximum(np.round(np.clip(v, 0, 1) * (len(SEQ) - 1)), 1), 0).astype(int)
            b = 0
            while b < len(v):
                e = b
                while e + 1 < len(v) and steps[e + 1] == steps[b]:
                    e += 1
                if steps[b] > 0:
                    x0, x1 = edges[b] * scale + ml, edges[e + 1] * scale + ml
                    tip = (f"{r.chrom}:{edges[b] / 1e6:.2f}-{edges[e + 1] / 1e6:.2f} Mb · {lab}: "
                           f"{100 * v[b:e + 1].mean():.0f}% masked")
                    g.append(f'<rect x="{x0:.1f}" y="{y}" width="{max(x1 - x0, 0.6):.1f}" height="{rowh}" '
                             f'fill="{SEQ[steps[b]]}"><title>{esc(tip)}</title></rect>')
                b = e + 1
        base = 4 + len(rows) * (rowh + gap) + 2 + (6 if pool_bins else 0)
        g.append(f'<line class="ax" x1="{ml}" x2="{ml + L * scale:.1f}" y1="{base}" y2="{base}"/>')
        step_mb = 5 if longest > 20e6 else (1 if longest > 4e6 else 0.1)
        for t in np.arange(0, L / 1e6 + 1e-9, step_mb):
            x = ml + t * 1e6 * scale
            g.append(f'<line class="ax" x1="{x:.1f}" x2="{x:.1f}" y1="{base}" y2="{base + 4}"/>'
                     f'<text class="tick" x="{x:.1f}" y="{base + 15}" text-anchor="middle">{t:g}</text>')
        w_bins = np.diff(edges) / L
        masked = [f"accessibility mask {100 * float((fr['any'] * w_bins).sum()):.1f}%"]
        for t in finals:
            masked.append(f"{final_label(t)} {100 * float((fr['final:' + t] * w_bins).sum()):.1f}%")
        # final masks in large blocks, for the table below
        bedges = np.append(np.arange(0, L, block), L)
        for t in finals:
            for b0, v in zip(bedges[:-1], final_fraction(t, r.chrom, L, bedges)):
                blocks.append({"target": t, "chrom": r.chrom, "block": int(b0 // block), "masked": float(v)})
        panels.append(f"<h3>{esc(r.chrom)} · {L / 1e6:.1f} Mb</h3><p class='sub'>{esc(' · '.join(masked))}</p>"
                      f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="mask along {esc(r.chrom)}">{"".join(g)}</svg>')
    key = "".join(f'<span style="display:inline-block;width:22px;height:10px;background:{c}"></span>'
                  for c in SEQ[1::2])
    legend = (f'<div class="legend"><span>fraction of each {bin_size / 1e3:g} kb bin masked: 0%&nbsp;{key}&nbsp;100%</span>'
              f'<span>x axis: Mb</span></div>')
    n_all = int(((summary["chrom"] != "total") & (summary["excluded"] < summary["length"])).sum())
    more = f" The {max_chroms} longest of {n_all} analysed contigs are shown." if n_all > max_chroms else ""
    region = region_table(pd.DataFrame(blocks), block, {t: final_label(t) for t in finals}) if blocks else ""
    return ("<section><h2>Where the mask falls</h2>"
            "<p>Each chromosome drawn to the same scale, with one row per mask reason, then the whole accessibility "
            "mask, then the final mask grenedalf uses for the chosen grid run (accessibility plus the depth bounds; "
            "with per-pool masks, the mean over pools of what each pool loses), then each pool's share of sites "
            "outside its own default depth bounds (from 1 kb depth bins). Long dark blocks in a single pool are "
            "regions that pool covers poorly or excessively (divergent haplotypes, deletions, copy-number changes); "
            "dark blocks in every pool are usually reference divergence; scattered speckle is noise. Darker = more "
            f"of the bin; hover for values.{more}</p>" + legend + "".join(panels) + region + "</section>")


def region_table(blocks, block, labels):
    """Final mask by large block, one row per chromosome: shows a region losing far more than
    the rest of the genome as a number, not just a shade."""
    out = []
    unit = block / 1e6
    for t, g in blocks.groupby("target", sort=False):
        wide = g.pivot(index="chrom", columns="block", values="masked").reindex(g["chrom"].unique())
        cols = [f"{b * unit:g}–{(b + 1) * unit:g}" for b in wide.columns]
        wide.columns = cols
        wide = wide.reset_index().rename(columns={"chrom": "chromosome"})
        med = float(np.nanmedian(g["masked"]))
        hot = g[g["masked"] > med + 0.15]
        note = (f" The genome-wide median block is {100 * med:.0f}% masked; "
                + (f"{len(hot)} block(s) lose more than 15 points over that: "
                   + ", ".join(f"{h.chrom} {h.block * unit:g}–{(h.block + 1) * unit:g} Mb ({100 * h.masked:.0f}%)"
                               for h in hot.sort_values("masked", ascending=False).head(8).itertuples())
                   + ("…" if len(hot) > 8 else "") + "." if len(hot) else "no block loses more than 15 points over that."))
        out.append(f"<h3>{esc(labels[t][0].upper() + labels[t][1:])}, by {unit:g} Mb block</h3>"
                   f"<p>Fraction of each block masked (Mb along the chromosome).{esc(note)}</p>"
                   + heat_table(wide, cols, lambda v: f"{100 * v:.0f}%", vmax=1.0))
    return "".join(out)


def read_bed3(path):
    return pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"],
                       dtype={"chrom": str})


def merge_bed(df):
    out = []
    for c, g in df.sort_values(["chrom", "start"]).groupby("chrom", sort=False):
        s, e = g["start"].to_numpy(), g["end"].to_numpy()
        emax = np.maximum.accumulate(e)
        idx = np.flatnonzero(np.concatenate(([True], s[1:] > emax[:-1])))
        out.append(pd.DataFrame({"chrom": c, "start": s[idx], "end": np.maximum.reduceat(e, idx)}))
    return pd.concat(out, ignore_index=True) if out else df


def depth_curves(P, dists):
    """Per pool: depth pmf for autosomes and X from mosdepth region.dist (weights from mask summary)."""
    from grenedalf_optimiser import read_dist  # same parser as the thresholds step
    summary = read_tsv(f"{P}.mask_summary.tsv", dtype={"chrom": str})
    classes = read_tsv(f"{P}.chrom_class.tsv", dtype={"chrom": str})
    if summary is None or classes is None:
        return {}
    weights = dict(zip(summary["chrom"], summary["accessible"]))
    cls = dict(zip(classes["chrom"], classes["class"]))
    out = {}
    for pool, path in dists:
        if not os.path.exists(path):
            continue
        acc = {}
        for chrom, pmf in read_dist(path).items():
            w = float(weights.get(chrom, 0))
            c = cls.get(chrom)
            if w <= 0 or c is None:
                continue
            a = acc.setdefault(c, np.zeros(1))
            n = max(a.size, pmf.size)
            acc[c] = np.pad(a, (0, n - a.size)) + np.pad(pmf * w, (0, n - pmf.size))
        out[pool] = {c: a / a.sum() for c, a in acc.items()}
    return out


def section_depth(P, dists):
    thr = read_tsv(f"{P}.depth_thresholds.tsv")
    if thr is None:
        return ""
    curves = depth_curves(P, dists) if dists else {}
    panels = []
    for pool, g in thr.groupby("pool", sort=False):
        g = g.set_index("class")
        xmax = float(g["upper"].max()) * 1.35
        cur = curves.get(pool, {})
        if not cur:
            continue
        # mosdepth reports P(depth >= d) to two decimals; plot that directly (differencing it is noisy)
        p = Plot((0, xmax), (0, 1), w=470, h=230, xlabel="read depth (MAPQ filtered)",
                 ylabel="sites with at least this depth", yfmt=lambda v: f"{100 * v:.0f}%")
        if "A" in g.index:
            p.band(g.loc["A", "lower"], g.loc["A", "upper"],
                   f"autosome bounds {g.loc['A', 'lower']}–{g.loc['A', 'upper']}")
        for c, color, name in (("A", "var(--s1)", "autosomes"), ("X", "var(--s2)", "X")):
            if c in cur:
                surv = (1 - np.cumsum(cur[c]) + cur[c])[: int(xmax) + 1]
                xs = np.arange(surv.size)
                if "fit_mu" in g and pd.notna(g.loc[c, "fit_mu"]):
                    from grenedalf_optimiser import nb_cdf
                    pi_ = float(g.loc[c, "fit_low_component"]) if "fit_low_component" in g else 0.0
                    fc = nb_cdf(float(g.loc[c, "fit_mu"]), float(g.loc[c, "fit_size"]), surv.size)
                    model = 1 - (pi_ + (1 - pi_) * np.concatenate(([0.0], fc[:-1])))[: surv.size]
                    p.line(xs, model, "var(--neutral)", "fitted model")
                lo, hi, med = g.loc[c, "lower"], g.loc[c, "upper"], g.loc[c, "median"]
                stepx = max(1, surv.size // 40)
                p.line(xs, surv, color, name,
                       tips=[f"{name}: {100 * surv[x]:.0f}% of sites at depth ≥ {x}; median {med}, kept {lo}–{hi}"
                             if x % stepx == 0 else "" for x in xs])
        a = g.loc["A"] if "A" in g.index else g.iloc[0]
        legend = ('<div class="legend"><span><span class="key" style="background:var(--s1)"></span>autosomes</span>'
                  + ('<span><span class="key" style="background:var(--s2)"></span>X</span>' if "X" in cur else "")
                  + '<span><span class="key" style="background:var(--s1-wash);height:10px"></span>kept (autosomes)</span>'
                  + ('<span><span class="key" style="background:var(--neutral)"></span>fitted coverage model</span>'
                     if "fit_mu" in g else "") + '</div>')
        panels.append(f"<div><h3>{esc(pool)} · median {a['median']}×, autosomes kept {a['lower']}–{a['upper']} "
                      f"({pct(a['frac_kept'])} of sites)</h3>{legend}{p.svg(f'{pool} depth distribution')}</div>")
    cols = [c for c in ("pool", "class", "median", "q05", "q95", "lower", "upper", "frac_below", "frac_above", "frac_kept")
            if c in thr]
    tbl = thr[cols].copy()
    for c in ("frac_below", "frac_above", "frac_kept"):
        if c in tbl:
            tbl[c] = tbl[c].map(pct)
    dm = read_tsv(f"{P}.depthmask_summary.tsv")
    cost = ""
    if dm is not None:
        cost = ("<h3>Accessible sites lost to the depth bounds</h3>"
                "<p>Per pool, and for the shared mask (a site outside the bounds in any pool is dropped in all, "
                "so every pool is called at the same sites).</p>"
                + hbars(dm["pool"].tolist(), dm["out_of_range_frac"].tolist(), pct,
                        [f"{r.pool}: {compact(r.out_of_range_bp)} bp ({pct(r.out_of_range_frac)})" for r in dm.itertuples()],
                        "sites lost to depth bounds", highlight="shared"))
    return ("<section><h2>Depth per pool</h2>"
            "<p>Each pool's depth inside the accessible genome, as the fraction of sites with at least a given "
            "depth; the shaded band is the range kept on autosomes, set from the pool's own median. The X is "
            "bounded separately because pools with males have lower X depth.</p>"
            + (f'<div class="grid2">{"".join(panels)}</div>' if panels else "")
            + f"<details><summary>Table</summary>{plain_table(tbl)}</details>" + cost + "</section>")


def section_min_count(P):
    mc = read_tsv(f"{P}.min_count.tsv")
    if mc is None:
        return ""
    cols = [c for c in mc.columns if c.startswith("err_sites_per_Mb_c")]
    t = mc[["pool", "median_depth"] + cols + ["min_count"]].copy()
    t.columns = ["pool", "median depth"] + [f"≥{c.rsplit('c', 1)[1]}" for c in cols] + ["error model"]
    vcols = [c for c in t.columns if c.startswith("≥")]
    chosen = dict(zip(t["pool"], t["error model"]))
    log1p = lambda v: math.log10(1 + v)
    mp = read_params(os.path.join(os.path.dirname(P), "params.sh"))
    d, f = mp.get("DIVERSITY_MIN_COUNT"), mp.get("FST_MIN_COUNT", mp.get("MIN_COUNT"))
    used = (f"{d} for diversity and {f} for FST" if d and f and d != f else (f or d))
    used_txt = (f" The count used is <strong>{used}</strong>: the error model's count is where the default "
                "run starts and is a floor: the filter grid tries one count above it, never below, because "
                "false SNPs from sequencing errors do not follow depth, so the grid cannot see them." if used else "")
    return ("<section><h2>Minimum allele count</h2>"
            f"<p>Expected sites per Mb where sequencing errors alone give at least <em>c</em> reads of one "
            f"non-reference base, at each pool's median depth and the error rate implied by the base-quality "
            f"cut-off ({mc['error_rate'].iloc[0]:g}); shading is on a log scale. The outlined count is the smallest "
            f"that keeps this below the accepted rate, and the largest across pools is the error model's count."
            f"{used_txt}</p>"
            + heat_table(t, vcols, lambda v: compact(round(v)) if v >= 10 else (f"{v:.1f}" if v >= 0.1 else
                                                 ("0" if v == 0 else "<0.1")), scale=log1p,
                         vmax=log1p(float(t[vcols].max().max())),
                         pick=lambda r, c: c == f"≥{chosen[r['pool']]}")
            + pool_size_note(mc)
            + "</section>")


def pool_size_note(mc):
    if "pool_size" not in mc:
        return ""
    t = mc[["pool", "pool_size", "median_depth", "n_eff", "min_detectable_freq", "single_copy_freq"]].copy()
    t["pool_size"] = t["pool_size"].astype(int)
    t.columns = ["pool", "pool size", "median depth", "n_eff", "lowest frequency counted", "one copy in pool"]
    low = (mc["min_detectable_freq"] > mc["single_copy_freq"]).any()
    return ("<h3>What the minimum count can see</h3><p>The lowest allele frequency reaching the minimum count "
            "at the median depth, against the frequency of a single copy in the pool."
            + (" Alleles present as one or a few copies in a pool fall below the minimum count; rare-variant "
               "statistics are affected more than common-variant ones." if low else "") + "</p>"
            + plain_table(t))


def section_variance(P):
    v = read_tsv(f"{P}.variance.tsv")
    thr = read_tsv(f"{P}.depth_thresholds.tsv")
    if v is None:
        return ""
    a = v[v["class"] == "A"] if (v["class"] == "A").any() else v
    ratio = float((a["cv"] / a["poisson_cv"]).median())
    s1, s10 = float(a["share_1kb"].median()), float(a["share_10kb"].median())
    blocks = float(a.loc[a["low_bins"] > 0, "low_in_blocks"].median()) if (a["low_bins"] > 0).any() else np.nan
    spec = float(a.loc[a["low_bins"] > 0, "pool_specific"].median()) if (a["low_bins"] > 0).any() else np.nan
    shared = float(a.loc[a["low_bins"] > 0, "shared_half_or_more"].median()) if (a["low_bins"] > 0).any() else np.nan
    corr_m = float(a["mean_corr_with_others"].median())
    li = [f"<li><strong>Spread:</strong> autosomal coverage varies {ratio:.1f}× more than read sampling alone "
          f"would (CV {a['cv'].median():.2f} against {a['poisson_cv'].median():.2f} for pure sampling noise). "
          "Fixed-fold bounds assume sampling-like spread; with this much more, 0.5× the median sits closer to the "
          "middle of the distribution than 2× does, so far more is cut below than above.</li>",
          f"<li><strong>Scale:</strong> {100 * s1:.0f}% of the per-base variance persists in 1 kb means and "
          f"{100 * s10:.0f}% in 10 kb means. "
          + ("Depth varies mostly in regional blocks (mapping, divergence, GC, copy number), not as scattered noise."
             if s1 >= 0.5 else "Much of the variation is local noise that averages out within 1 kb.") + "</li>"]
    if np.isfinite(blocks):
        li.append(f"<li><strong>Low depth:</strong> {100 * blocks:.0f}% of the 1 kb bins below a pool's lower bound "
                  f"lie in runs of at least 5 kb. Of a pool's low bins, {100 * spec:.0f}% are low in that pool only "
                  f"(population divergence, deletions: biology to report, not to mask everywhere) and "
                  f"{100 * shared:.0f}% are low in at least half the pools (reference or technical).</li>")
    li.append(f"<li><strong>Between pools:</strong> 1 kb depth correlates at {corr_m:.2f} (median), so most "
              "coverage structure is shared by all pools.</li>")
    cols = ["pool", "class", "mean", "cv", "poisson_cv", "share_1kb", "share_10kb", "low_bins_frac", "low_in_blocks",
            "pool_specific", "shared_half_or_more"]
    tv = v[[c for c in cols if c in v]].copy()
    for c in ("share_1kb", "share_10kb", "low_bins_frac", "low_in_blocks", "pool_specific", "shared_half_or_more"):
        if c in tv:
            tv[c] = tv[c].map(lambda x: "–" if pd.isna(x) else pct(x))
    tv.columns = [{"cv": "CV", "poisson_cv": "CV (sampling only)", "share_1kb": "variance at 1 kb",
                   "share_10kb": "at 10 kb", "low_bins_frac": "low 1 kb bins", "low_in_blocks": "low, in blocks",
                   "pool_specific": "low, this pool only", "shared_half_or_more": "low, ≥ half the pools"}.get(c, c)
                  for c in tv.columns]
    corr = read_tsv(f"{P}.depth_correlation.tsv")
    corr_html = ""
    if corr is not None:
        corr = corr.rename(columns={corr.columns[0]: "pool"})
        corr_html = ("<h3>1 kb depth correlation between pools (Spearman)</h3>"
                     + heat_table(corr, list(corr.columns[1:]), lambda x: f"{x:.2f}", vmax=1.0))
    bounds_html = ""
    if thr is not None and "fit_lower" in thr:
        tb = thr[["pool", "class", "median", "lower", "upper", "frac_below", "frac_above",
                  "fit_lower", "fit_upper", "fit_frac_below", "fit_frac_above", "fit_low_component"]].copy()
        for c in ("frac_below", "frac_above", "fit_frac_below", "fit_frac_above", "fit_low_component"):
            tb[c] = tb[c].map(pct)
        tb.columns = ["pool", "class", "median", "fold: from", "to", "cut below", "cut above",
                      "fitted: from", "to ", "cut below ", "cut above ", "low-depth component"]
        bounds_html = ("<h3>Fixed-fold against fitted bounds</h3><p>Fitted bounds are the 0.5% and 99.5% points of a "
                       "negative binomial fitted to each pool's main coverage, with a separate low-depth component "
                       "(sites whose reads largely fail to map), so the bounds follow the pool's own spread and cut "
                       "the same tail probability on both sides. The grid runs both (<code>fitted_bounds</code>) and the "
                       "artefact test decides.</p>" + plain_table(tb))
    return ("<section><h2>Coverage variance</h2><p>What kind of depth variation there is decides what a depth "
            "bound removes: sampling noise (harmless, and the estimators handle it), regional blocks shared by all "
            "pools (reference or technical: mask), or blocks in single pools (biology: report).</p>"
            "<ul>" + "".join(li) + "</ul>" + plain_table(tv) + corr_html + bounds_html + "</section>")


def section_chromclass(P):
    t = read_tsv(f"{P}.chrom_classes.tsv")
    if t is None:
        return ""
    cc = read_params(f"{P}.chromclass.params")
    x = t[t["class"] == "X"]
    mito = read_tsv(f"{P}.mito.tsv")
    parts = []
    if len(x):
        unreliable = cc.get("X_RATIO_UNRELIABLE") == "1"
        plaus = x["plausible"].astype(str) == "True" if "plausible" in x else pd.Series(True, index=x.index)
        # Dot with 95% interval per pool, on the X / autosome depth ratio (peak to peak)
        rowh, ml, w = 28, 100, 620
        h = len(x) * rowh + 52
        lo_v = min(0.4, float(x["ci_low"].min()) - 0.05)
        hi_v = max(1.15, float(x["ci_high"].max()) + 0.05)
        sx = lambda v: ml + (min(max(v, lo_v), hi_v) - lo_v) / (hi_v - lo_v) * (w - ml - 20)
        g = [f'<rect x="{sx(1.0):.1f}" y="18" width="{sx(hi_v) - sx(1.0):.1f}" height="{h - 44}" fill="var(--page)">'
             f'<title>X depth above autosome depth: not possible from copy number alone</title></rect>']
        for v in np.arange(np.ceil(lo_v * 10) / 10, hi_v + 1e-9, 0.1):
            g.append(f'<line class="gl" x1="{sx(v):.1f}" x2="{sx(v):.1f}" y1="22" y2="{h - 26}"/>'
                     f'<text class="tick" x="{sx(v):.1f}" y="{h - 10}" text-anchor="middle">{v:.1f}</text>')
        for v, lab, anchor in ((1.0, "no males", "end"), (0.75, "50:50", "middle"), (0.5, "all males", "start")):
            g.append(f'<line x1="{sx(v):.1f}" x2="{sx(v):.1f}" y1="18" y2="{h - 26}" stroke="var(--ink2)" stroke-width="1"/>'
                     f'<text class="ref" x="{sx(v) + {"end": -4, "start": 4, "middle": 0}[anchor]:.1f}" y="14" '
                     f'text-anchor="{anchor}">{lab}</text>')
        for i, (r, ok) in enumerate(zip(x.itertuples(), plaus)):
            y = 34 + i * rowh
            col = "var(--s1)" if ok else "var(--neutral)"
            g.append(f'<text class="tick" x="{ml - 8}" y="{y + 4}" text-anchor="end">{esc(r.pool)}</text>')
            g.append(f'<line x1="{sx(r.ci_low):.1f}" x2="{sx(r.ci_high):.1f}" y1="{y}" y2="{y}" stroke="{col}" '
                     f'stroke-width="2" stroke-linecap="round"/>')
            tip = (f"{r.pool}: X/autosome depth {r.ratio_to_A:.3f} (95% {r.ci_low:.3f}–{r.ci_high:.3f}; "
                   f"medians {getattr(r, 'ratio_median', np.nan):.3f}); male fraction used {r.male_fraction_used:.2f} "
                   f"({getattr(r, 'male_fraction_source', '')}); X pool size {r.pool_size:g}")
            g.append(f'<circle cx="{sx(r.ratio_to_A):.1f}" cy="{y}" r="5" fill="{col}" stroke="var(--surface)" stroke-width="2"/>'
                     f'<circle cx="{sx(r.ratio_to_A):.1f}" cy="{y}" r="11" fill="transparent"><title>{esc(tip)}</title></circle>')
        g.append(f'<text class="alabel" x="{(ml + w - 20) / 2:.0f}" y="{h + 6}" text-anchor="middle">'
                 f'X depth / autosome depth, peak to peak, well-mapped bins (95% block-bootstrap interval)</text>')
        svg = f'<svg viewBox="0 0 {w} {h + 12}" role="img" aria-label="X to autosome depth ratio">{"".join(g)}</svg>'
        legend = ('<div class="legend"><span><span class="dotkey" style="background:var(--s1)"></span>possible sex ratio (0.5–1)</span>'
                  '<span><span class="dotkey" style="background:var(--neutral)"></span>outside what a sex ratio can give</span></div>')
        verdict = ""
        if unreliable:
            src = x["male_fraction_source"].iloc[0] if "male_fraction_source" in x else ""
            bad = x.loc[~plaus]
            why = (bad["unreliable_because"].fillna("").astype(str) if "unreliable_because" in bad
                   else pd.Series("", index=bad.index))
            parts = []
            if why.str.contains("outside").any():
                parts.append(f"the X / autosome depth ratio is outside 0.5–1 in {int(why.str.contains('outside').sum())} "
                             f"of {len(x)} pools")
            sp = bad.loc[why.str.contains("spread"), "spread_ratio"].astype(float) if "spread_ratio" in bad else []
            if len(sp):
                rng = f"{sp.min():.1f}×" if round(sp.min(), 1) == round(sp.max(), 1) else f"{sp.min():.1f}–{sp.max():.1f}×"
                parts.append(f"X and autosome depth spreads differ {rng} ({len(sp)} of {len(x)} pools)")
            reasons = "; ".join(parts)
            verdict += ('<p class="status">' + ICON_BAD + "X depth does not reflect copy number in this data</p>"
                        f"<p>{esc(reasons[:1].upper() + reasons[1:]) + '. ' if reasons else ''}A sex ratio can only give X depth "
                        "between 0.5 and 1× the autosomes, and only if depth behaves alike on both; here it is shaped "
                        "by something else, most likely reference divergence (reads from haplotypes divergent from "
                        "the reference fail to map, more on the autosomes than the X, which also widens the autosomal "
                        "depth distribution). So no pool's estimate is used: the male fraction is "
                        + ("taken from <code>MALE_FRACTION</code>." if src == "known" else
                           "assumed to be 0.5. <strong>Set <code>MALE_FRACTION</code> from the study design.</strong>")
                        + "</p>")
        sep = cc.get("X_SEPARATE") == "1"
        verdict += (f'<p class="status">{ICON_INFO if sep else ICON_OK}'
                    + ("The X gets its own grenedalf run, with the X pool sizes below" if sep
                       else "The X can be analysed with the autosomes") + "</p>")
        cols = [c for c in ["pool", "ratio_to_A", "ci_low", "ci_high", "ratio_median", "cv_A", "cv_X", "male_fraction_est",
                            "male_fraction_known", "male_fraction_used", "male_fraction_source", "pool_size"] if c in x]
        tx = x[cols].copy()
        tx.columns = [{"ratio_to_A": "X / A depth (peaks)", "ci_low": "95% low", "ci_high": "95% high",
                       "ratio_median": "X / A (medians)", "cv_A": "spread (CV) A", "cv_X": "spread (CV) X",
                       "male_fraction_est": "male fraction (est.)",
                       "male_fraction_known": "male fraction (known)", "male_fraction_used": "used",
                       "male_fraction_source": "from", "pool_size": "X pool size"}.get(c, c) for c in cols]
        if "matches_known" in x and x["male_fraction_known"].notna().any() and not unreliable \
                and not x["matches_known"].fillna(True).astype(bool).all():
            verdict += ('<p class="status">' + ICON_BAD + "X depth does not match the known sex ratio in at least one "
                        "pool: check the pool composition, or X mappability</p>")
        parts.append("<h3>X chromosome</h3><p>Males carry one X, so in a pool with males the X has lower depth and "
                     "is sampled from fewer copies. With a fraction <em>m</em> of males, X depth / autosome depth = "
                     "1 − <em>m</em>/2 (0.75 for a 50:50 pool), and the X pool size is the pool size × (1 − <em>m</em>/2). "
                     "The ratio is taken between the peaks of the two depth distributions, over well-mapped bins, "
                     "because medians are pulled down wherever reads fail to map. <code>MALE_FRACTION</code>, when "
                     "set, is used instead, with the estimate kept as a check.</p>"
                     + legend + svg + verdict + plain_table(tx))
    if mito is not None and len(mito):
        tm = mito[[c for c in ["pool", "median_depth", "ratio_to_A", "lower", "upper", "masked_bp", "pool_size"]
                   if c in mito]].copy()
        tm.columns = [{"median_depth": "median depth", "ratio_to_A": "× autosome depth", "lower": "kept from",
                       "upper": "kept to", "masked_bp": "bp masked", "pool_size": "pool size (individuals)"}.get(c, c)
                      for c in tm.columns]
        parts.append("<h3>Mitochondrion</h3><p>Many copies per cell, haploid, and one lineage per individual: its "
                     "depth is far from the autosomes' and differs between pools (copy number varies), and its pool "
                     "size is the number of individuals, not genomes. It gets its own run, with each pool's depth "
                     "bounds applied as that pool's mask, and the whole genome as one window.</p>" + plain_table(tm))
    if not parts:
        return ""
    return "<section><h2>X chromosome and mitochondrion</h2>" + "".join(parts) + "</section>"


def section_ld(P):
    t = read_tsv(f"{P}.ld_decay.tsv")
    ld = read_params(f"{P}.ld.params")
    if t is None:
        return ""
    t = t[t["pairs"] > 0]
    base = float(ld.get("LD_BASELINE_R2", "nan"))
    ymax = float(np.nanmax(t["mean_r2"])) * 1.08
    ymin = min(float(np.nanmin(t["mean_r2"])), base) * 0.92
    p = Plot((t["dist_from"].min(), t["dist_to"].max()), (ymin, ymax), xlog=True, w=900, h=310, mt=26,
             xlabel="distance between SNPs (bp)", ylabel="mean r² across pools", yfmt=lambda v: f"{v:.2f}")
    if np.isfinite(base):
        p.hline(base, f"unlinked baseline {base:.3f}", anchor="start")
    for key, lab, anchor in (("LD_D50", "50% decay", "end"), ("LD_D10", "90% decay", "start")):
        if ld.get(key, "NA") != "NA":
            p.vline(float(ld[key]), f"{lab} {bp(ld[key])}", anchor=anchor)
    p.line(t["dist_mid"], t["mean_r2"], "var(--s1)", "mean r²", dots=True,
           tips=[f"{r.dist_from}–{r.dist_to} bp: r² {r.mean_r2:.3f} ± {r.se:.3f} ({r.pairs:,} pairs)" for r in t.itertuples()])
    return ("<section><h2>Linkage scale</h2>"
            "<p>Correlation, across pools, of the allele frequencies of SNP pairs at increasing distance. Linked SNPs "
            "share their frequency differences between pools; with distance this falls to the baseline of "
            "unlinked pairs. Beyond the 90% decay distance SNPs carry largely independent information, so windows "
            f"are not made narrower than that. Based on {compact(ld.get('LD_SNPS', 0))} SNPs.</p>"
            + p.svg("linkage disequilibrium decay")
            + f"<details><summary>Table</summary>{plain_table(t)}</details></section>")


def section_windows(P):
    t = read_tsv(f"{P}.window_widths.tsv")
    win = read_params(f"{P}.windows.params")
    ld = read_params(f"{P}.ld.params")
    txt = open(f"{P}.windows.txt").read() if os.path.exists(f"{P}.windows.txt") else ""
    if t is None:
        return ""
    status = ""
    m = re.search(r"mask check: (.*)", txt)
    if m:
        ok = "-> OK" in m.group(1)
        status = (f'<p class="status">{ICON_OK if ok else ICON_BAD}'
                  f'{"Mask applied as intended" if ok else "Mask check failed"}</p><p class="sub">{esc(m.group(1))}</p>')
    xmax = min(float(t["width"].max()), max(4 * float(win.get("WINDOW", t["width"].max())), 10000))
    t = t[t["width"] <= xmax]
    p = Plot((0, xmax), (0, 1), w=900, h=290, mt=40, xlabel="window width (bp)",
             ylabel=f"windows with ≥ {win.get('MIN_SNPS', '')} SNPs", yfmt=lambda v: f"{100 * v:.0f}%")
    p.line(t["width"], t["frac_windows_ge_target"], "var(--s1)", "windows reaching target",
           tips=[f"{r.width:,} bp: {pct(r.frac_windows_ge_target)} of windows, median {r.median_snps:g} SNPs"
                 for r in t.itertuples()])
    # Reference lines; one label per distinct position
    refs = {}
    if "WINDOW_SNPS" in win:
        refs.setdefault(float(win["WINDOW_SNPS"]), []).append("SNP target")
    if ld.get("LD_D10", "NA") != "NA":
        refs.setdefault(float(ld["LD_D10"]), []).append("linkage scale")
    chosen = float(win["WINDOW"]) if "WINDOW" in win else None
    for i, (x, names) in enumerate(sorted(refs.items())):
        if chosen is not None and x == chosen:
            names = ["chosen = " + " = ".join(names)]
        p.vline(x, f"{', '.join(names)} {bp(x)}", anchor="end" if i == 0 and len(refs) > 1 else "start", dy=14 * i)
    if chosen is not None and chosen not in refs:
        p.vline(chosen, f"chosen {bp(chosen)}", anchor="start", dy=14 * len(refs))
    lines = [l for l in txt.splitlines() if not l.startswith("mask check")]
    return ("<section><h2>Window size</h2>"
            "<p>Fraction of pilot windows that would hold enough SNPs for a precise estimate, by width. "
            "The chosen width is the larger of the SNP target and the linkage scale.</p>"
            + p.svg("window width") + status
            + "<ul>" + "".join(f"<li>{esc(l)}</li>" for l in lines) + "</ul></section>")


def section_grid(P):
    runs = read_tsv(f"{P}.assess_runs.tsv")
    stats = read_tsv(f"{P}.assess_stats.tsv")
    deltas = read_tsv(f"{P}.assess_deltas.tsv")
    targets = [t for t in TARGETS if runs is not None and f"data_kept_{t}" in runs]
    if not targets:
        return ""
    assess = read_params(f"{P}.assess.params")
    tol = float(assess.get("TOLERANCE", "nan"))
    ref_run = assess.get("REFERENCE_RUN", "")
    judged = [c for c in assess.get("JUDGED_COVARIATES", "").split(",") if c]
    cov_label = {"rho_rel_depth": "depth relative to the other pools", "rho_depth_imbalance": "depth imbalance",
                 "rho_callable": "callable fraction (each pool's own; for a pair, sites both can use)", "rho_depth": "depth"}
    metric_label = {"pair_snps": "SNPs usable per pair of pools (mean over pairs)",
                    "pool_snps": "SNPs usable per pool (mean over pools)", "snps": "SNPs"}
    role = runs["role"] if "role" in runs else pd.Series("candidate", index=runs.index)
    names = [r.replace("_", " ") + {"reference": " (reference)", "check": " (check)"}.get(ro, "")
             for r, ro in zip(runs["run"], role)]
    legend = ('<div class="legend"><span><span class="dotkey" style="background:var(--s2)"></span>chosen</span>'
              '<span><span class="dotkey" style="background:var(--s1)"></span>accepted (no statistic gets worse)</span>'
              '<span><span class="dotkey" style="background:var(--neutral)"></span>rejected, or a check</span></div>')
    nice = lambda c: c.replace(".rho_", " vs ").replace("_", " ")
    blocks = []
    for t in targets:
        rec_mask = runs[f"recommended_{t}"].astype(str) == "True"
        ok_mask = runs[f"acceptable_{t}"].astype(str) == "True"
        colors = ["var(--s2)" if r else ("var(--s1)" if a else "var(--neutral)") for r, a in zip(rec_mask, ok_mask)]
        kept = runs[f"data_kept_{t}"].astype(float)
        metric = str(runs[f"data_metric_{t}"].iloc[0])
        base = float(kept[runs["run"] == ref_run].iloc[0]) if (runs["run"] == ref_run).any() else float(kept.max())
        rel = lambda v: "" if v == base else f" ({100 * (v / base - 1):+.1f}%)"
        left = hbars(names, kept.tolist(), lambda v: compact(v) + rel(v),
                     [f"{r}: {compact(v)} ({metric})" for r, v in zip(runs["run"], kept)],
                     f"{t} data kept per grid run", colors=colors, ml=190)
        inc = runs[f"max_increase_{t}"].astype(float).fillna(0).clip(lower=0)
        on = runs[f"max_increase_on_{t}"].fillna("").astype(str)
        right = hbars(names, inc.tolist(), lambda v: f"+{v:.2f}" if v > 0 else "none",
                      [f"{r}: largest increase {v:+.3f} ({o})" for r, v, o in zip(runs["run"], inc, on)],
                      f"{t}: largest increase per grid run", colors=colors, ml=190,
                      xmax=max(float(inc.max()), tol if np.isfinite(tol) else 0, 0.01) * 1.15,
                      ref=(tol, f"accepted up to +{tol:.2f}") if np.isfinite(tol) else None)
        rec = runs.loc[rec_mask, "run"]
        rec = rec.iloc[0] if len(rec) else ref_run
        delta_html = ""
        if deltas is not None:
            dcols = [c for c in deltas.columns if ".rho_" in c and c.split(".rho_")[0] in TARGET_TYPES[t]
                     and (not judged or "rho_" + c.split(".rho_")[1] in judged)]
            if dcols:
                dt = deltas[["run"] + dcols].copy()
                dt.columns = ["run"] + [nice(c) for c in dcols]
                vmax = max(float(np.nanmax(np.clip(dt.iloc[:, 1:].to_numpy(dtype=float), 0, None))), tol, 0.05)
                delta_html = ("<p>Change in each statistic's dependence on each judged covariate (median |ρ| of the "
                              "run minus that of the reference); positive (shaded) = more dependence.</p>"
                              + heat_table(dt, list(dt.columns[1:]), lambda v: f"{v:+.2f}", vmax=vmax,
                                           scale=lambda v: max(v, 0.0), pick_row=lambda r, rec=rec: r["run"] == rec))
        combo = f"combined_{t}"
        combo_html = (f"<p>The <strong>{combo.replace('_', ' ')}</strong> run puts together the best accepted choice "
                      f"of each setting for {t if t == 'diversity' else 'FST'}, and is judged by the same rule.</p>"
                      if (runs["run"] == combo).any() else "")
        blocks.append(f"<h3>{esc(TARGETS[t][0].upper() + TARGETS[t][1:])}: <strong>{esc(rec.replace('_', ' '))}</strong></h3>"
                      f"<p>Data kept: {esc(metric_label.get(metric, metric).lower())}.</p>" + combo_html + legend
                      + f'<div class="grid2"><div>{left}</div><div>{right}</div></div>' + delta_html)
    # FST not judged: its dependence on the covariates is shown for information only
    if "fst" not in targets and deltas is not None:
        fcols = [c for c in deltas.columns if c.startswith("fst.rho_")
                 and (not judged or "rho_" + c.split(".rho_")[1] in judged)]
        if fcols:
            ft = deltas[["run"] + fcols].copy()
            ft.columns = ["run"] + [nice(c) for c in fcols]
            rec = chosen_runs(P)["fst"] or ref_run
            blocks.append(
                f"<h3>FST (per pair of pools): <strong>{esc(rec.replace('_', ' '))}</strong>, not judged</h3>"
                "<p>FST uses the settings chosen for diversity. It is not judged, because a real difference "
                "between pools in a region where every pool has low depth (divergent haplotypes that map poorly) "
                "makes FST follow depth and callable fraction across the genome just as an artefact would, so a "
                "single region can reject every setting that keeps it. The change in FST's dependence on each "
                "covariate is shown for information (JUDGE_FST=1 in the config judges it, with its own choice).</p>"
                + heat_table(ft, list(ft.columns[1:]), lambda v: f"{v:+.2f}",
                             vmax=max(float(np.nanmax(np.clip(ft.iloc[:, 1:].to_numpy(dtype=float), 0, None))),
                                      tol if np.isfinite(tol) else 0, 0.05),
                             scale=lambda v: max(v, 0.0), pick_row=lambda r, rec=rec: r["run"] == rec))

    cov_cols = [c for c in runs.columns if ".rho_" in c]
    heat = runs[["run"] + cov_cols].copy()
    heat.columns = ["run"] + [nice(c) for c in cov_cols]
    failed_html = ""
    failed_path = os.path.join(os.path.dirname(P), "grid", "failed_runs.tsv")
    if os.path.exists(failed_path) and os.path.getsize(failed_path):
        fr = pd.read_csv(failed_path, sep="\t", header=None, names=["run", "error"])
        failed_html = ('<p class="status">' + ICON_BAD + f"{len(fr)} grid run(s) failed and are left out of the "
                       "comparison</p>" + plain_table(fr))
    tbl_cols = ["run", "role", "pair_snps", "pool_snps", "callable_bp", "snps", "windows_scored"] + [
        f"{k}_{t}" for t in targets for k in ("max_increase", "max_increase_on", "acceptable", "recommended")]
    detail = ""
    if stats is not None:
        detail = ("<details><summary>Every statistic column (per pool / pair)</summary>"
                  + plain_table(stats) + "</details>")
    judged_txt = ", ".join(cov_label.get(c, c) for c in judged) or "the technical covariates"
    return ("<section><h2>Filter grid</h2>"
            "<p>The genome rerun with one setting changed at a time from the default settings. For each run and "
            "each statistic (θπ, θW per pool; FST per pair), the statistic is correlated across windows with "
            f"{esc(judged_txt)}; none of these should matter. Depth relative to the other pools is used rather "
            "than depth itself: where every pool has low depth (often divergent haplotypes that map poorly, which "
            "are also the most diverse regions) diversity drives depth, not the reverse, and judging on raw depth "
            "would penalise keeping those regions. Each run is compared with the default settings "
            f"(<strong>{esc(ref_run.replace('_', ' '))}</strong>), statistic by statistic: it is accepted if no "
            "statistic's dependence on any covariate rises by more than the tolerance, and of the accepted runs the "
            "one keeping the most data is chosen. "
            + ("Diversity and FST are separate grenedalf runs, so each gets its own choice.</p>" if "fst" in targets
               else "The choice is made on diversity, a within-pool measure, and FST uses the same settings.</p>")
            + failed_html + "".join(blocks)
            + "<details><summary>Median |ρ| of each run, every covariate (absolute depth reported, not judged)</summary>"
            + heat_table(heat, list(heat.columns[1:]), lambda v: f"{v:.2f}",
                         vmax=max(float(np.nanmax(heat.iloc[:, 1:].to_numpy(dtype=float))), 0.1)) + "</details>"
            + f"<details><summary>Run table</summary>{plain_table(runs[[c for c in tbl_cols if c in runs]])}</details>"
            + detail + "</section>")


def build_report(args):
    P = args.prefix
    out_dir = args.out_dir or os.path.dirname(P)
    dists = args.dist or []
    if not dists and out_dir and os.path.isdir(os.path.join(out_dir, "profile")):
        prof = os.path.join(out_dir, "profile")
        dists = [(f[: -len(".mosdepth.region.dist.txt")], os.path.join(prof, f)) for f in sorted(os.listdir(prof))
                 if f.endswith(".mosdepth.region.dist.txt")]
    sections = [section_summary(P, out_dir), section_mask(P), section_karyoplot(P, out_dir), section_depth(P, dists),
                section_variance(P), section_chromclass(P), section_min_count(P),
                section_ld(P), section_windows(P), section_grid(P)]
    title = args.title or "grenedalf optimiser report"
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1"><title>Grenedalf optimiser report</title>'
            f"<style>{CSS}</style></head><body><div class='viz-root'><main>"
            f'<button class="theme" id="theme" type="button">Light / dark</button>'
            f"<h1>{esc(title)}</h1><p class='sub'>{esc(os.path.abspath(P))}.* · "
            f"{pd.Timestamp.now():%Y-%m-%d %H:%M}</p>"
            + "".join(s for s in sections if s) + f"</main></div>{THEME_JS}</body></html>")
    with open(args.out, "w") as fh:
        fh.write(page)
    print(f"wrote {args.out}")
