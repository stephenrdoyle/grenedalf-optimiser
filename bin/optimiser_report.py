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
    d = df if max_rows is None else df.head(max_rows)
    return heat_table(d, [], str)


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

def section_summary(P, out_dir):
    params = read_params(os.path.join(out_dir, "params.sh")) if out_dir else {}
    win = read_params(f"{P}.windows.params")
    ld = read_params(f"{P}.ld.params")
    mask = read_tsv(f"{P}.mask_summary.tsv", dtype={"chrom": str})
    # The grid may choose different depth folds or mask mode: report the chosen run's masks
    run = read_params(f"{P}.assess.params").get("RECOMMENDED_RUN")
    run_dm = os.path.join(out_dir, "grid", run, "mask.depthmask_summary.tsv") if run and out_dir else None
    dm = read_tsv(run_dm) if run_dm and os.path.exists(run_dm) else read_tsv(f"{P}.depthmask_summary.tsv")
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
                          "of the accessible genome, after depth bounds"))
    mc = params.get("MIN_COUNT") or read_params(f"{P}.thresholds.params").get("MIN_COUNT")
    if mc:
        tiles.append(tile("Minimum allele count", mc, "reads supporting an allele"))
    if win:
        tiles.append(tile("Window", bp(win['WINDOW']), f"stride {bp(win['STRIDE'])}"))
    if ld.get("LD_D10", "NA") != "NA":
        tiles.append(tile("Linkage scale", bp(ld['LD_D10']), f"90% decay; half at {bp(ld['LD_D50'])}"))
    if run:
        tiles.append(tile("Grid choice", run.replace("_", " "), "most permissive without added artefact"))
    body = f'<div class="tiles">{"".join(tiles)}</div>'

    thr = read_tsv(f"{P}.depth_thresholds.tsv")
    if thr is not None:
        rows = []
        for pool, g in thr.groupby("pool", sort=False):
            g = g.set_index("class")
            r = {"pool": pool}
            if "pool_size" in g:
                r["pool size"] = int(g["pool_size"].iloc[0])
            chosen = dm.set_index("pool").loc[pool] if dm is not None and pool in set(dm["pool"]) else None
            for c, lab in (("A", "autosomes"), ("X", "X")):
                if c in g.index:
                    lo, hi = g.loc[c, "lower"], g.loc[c, "upper"]
                    if chosen is not None and pd.notna(chosen.get(f"lower_{c}", np.nan)):
                        lo, hi = int(chosen[f"lower_{c}"]), int(chosen[f"upper_{c}"])
                    r[f"{lab}: median"] = f"{g.loc[c, 'median']}×"
                    r[f"{lab}: kept"] = f"{lo}–{hi}×"
            if "n_eff_at_median" in g and "A" in g.index:
                r["n_eff (autosomes)"] = g.loc["A", "n_eff_at_median"]
            rows.append(r)
        mp = read_params(os.path.join(out_dir, "params.sh")) if out_dir else {}
        glob = (f" grenedalf's own <code>--filter-sample-min-read-depth</code> / <code>--filter-sample-max-read-depth</code> "
                f"are set to {mp['MIN_DEPTH']} / {mp['MAX_DEPTH']}, wide enough not to undo these."
                if "MIN_DEPTH" in mp else "")
        body += ("<h3>Depth bounds per pool</h3>"
                 "<p>grenedalf has a single depth filter for all pools, so each pool's own range is applied as a "
                 "mask: a site outside any pool's range is masked in every pool, so all pools are called at the same "
                 f"sites.{glob} n_eff is the number of chromosomes effectively sampled at that depth from a pool of "
                 "that size; depth well beyond the pool size adds little.</p>" + plain_table(pd.DataFrame(rows)))

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
                lo, hi, med = g.loc[c, "lower"], g.loc[c, "upper"], g.loc[c, "median"]
                stepx = max(1, surv.size // 40)
                p.line(xs, surv, color, name,
                       tips=[f"{name}: {100 * surv[x]:.0f}% of sites at depth ≥ {x}; median {med}, kept {lo}–{hi}"
                             if x % stepx == 0 else "" for x in xs])
        a = g.loc["A"] if "A" in g.index else g.iloc[0]
        legend = ('<div class="legend"><span><span class="key" style="background:var(--s1)"></span>autosomes</span>'
                  + ('<span><span class="key" style="background:var(--s2)"></span>X</span>' if "X" in cur else "")
                  + '<span><span class="key" style="background:var(--s1-wash);height:10px"></span>kept (autosomes)</span></div>')
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
    t.columns = ["pool", "median depth"] + [f"≥{c.rsplit('c', 1)[1]}" for c in cols] + ["chosen"]
    vcols = [c for c in t.columns if c.startswith("≥")]
    chosen = dict(zip(t["pool"], t["chosen"]))
    log1p = lambda v: math.log10(1 + v)
    return ("<section><h2>Minimum allele count</h2>"
            f"<p>Expected sites per Mb where sequencing errors alone give at least <em>c</em> reads of one "
            f"non-reference base, at each pool's median depth and the error rate implied by the base-quality "
            f"cut-off ({mc['error_rate'].iloc[0]:g}); shading is on a log scale. The outlined count is the smallest "
            f"that keeps this below the accepted rate; the run uses the largest across pools.</p>"
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


def section_chromclass(P):
    t = read_tsv(f"{P}.chrom_classes.tsv")
    if t is None:
        return ""
    cc = read_params(f"{P}.chromclass.params")
    x = t[t["class"] == "X"]
    m = t[t["class"] == "M"]
    parts = []
    if len(x):
        # Dot with 95% interval per pool, on the X / autosome depth ratio
        rowh, ml, w = 28, 100, 620
        h = len(x) * rowh + 52
        lo_v, hi_v = 0.4, 1.15
        sx = lambda v: ml + (min(max(v, lo_v), hi_v) - lo_v) / (hi_v - lo_v) * (w - ml - 20)
        g = []
        for v in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1):
            g.append(f'<line class="gl" x1="{sx(v):.1f}" x2="{sx(v):.1f}" y1="22" y2="{h - 26}"/>'
                     f'<text class="tick" x="{sx(v):.1f}" y="{h - 10}" text-anchor="middle">{v:.1f}</text>')
        for v, lab, anchor in ((1.0, "no males", "end"), (0.75, "50:50", "middle"), (0.5, "all males", "start")):
            g.append(f'<line x1="{sx(v):.1f}" x2="{sx(v):.1f}" y1="18" y2="{h - 26}" stroke="var(--ink2)" stroke-width="1"/>'
                     f'<text class="ref" x="{sx(v) + {"end": -4, "start": 4, "middle": 0}[anchor]:.1f}" y="14" '
                     f'text-anchor="{anchor}">{lab}</text>')
        for i, r in enumerate(x.itertuples()):
            y = 34 + i * rowh
            col = "var(--s2)" if r.differs else "var(--s1)"
            g.append(f'<text class="tick" x="{ml - 8}" y="{y + 4}" text-anchor="end">{esc(r.pool)}</text>')
            g.append(f'<line x1="{sx(r.ci_low):.1f}" x2="{sx(r.ci_high):.1f}" y1="{y}" y2="{y}" stroke="{col}" '
                     f'stroke-width="2" stroke-linecap="round"/>')
            tip = (f"{r.pool}: X/autosome depth {r.ratio_to_A:.3f} (95% {r.ci_low:.3f}–{r.ci_high:.3f}); "
                   f"male fraction est. {r.male_fraction_est:.2f}, used {r.male_fraction_used:.2f}; "
                   f"X pool size {r.pool_size:g}")
            g.append(f'<circle cx="{sx(r.ratio_to_A):.1f}" cy="{y}" r="5" fill="{col}" stroke="var(--surface)" stroke-width="2"/>'
                     f'<circle cx="{sx(r.ratio_to_A):.1f}" cy="{y}" r="11" fill="transparent"><title>{esc(tip)}</title></circle>')
        g.append(f'<text class="alabel" x="{(ml + w - 20) / 2:.0f}" y="{h + 6}" text-anchor="middle">'
                 f'median X depth / median autosome depth (95% block-bootstrap interval)</text>')
        svg = f'<svg viewBox="0 0 {w} {h + 12}" role="img" aria-label="X to autosome depth ratio">{"".join(g)}</svg>'
        legend = ('<div class="legend"><span><span class="dotkey" style="background:var(--s2)"></span>differs from autosomes</span>'
                  '<span><span class="dotkey" style="background:var(--s1)"></span>consistent with autosomes</span></div>')
        sep = cc.get("X_SEPARATE") == "1"
        verdict = (f'<p class="status">{ICON_INFO if sep else ICON_OK}'
                   + ("Analyse the X separately, with its own pool sizes" if sep else "The X can be analysed with the autosomes")
                   + "</p>")
        cols = ["pool", "ratio_to_A", "ci_low", "ci_high", "male_fraction_est", "male_fraction_ci"]
        names_ = ["pool", "X / autosome depth", "95% low", "95% high", "male fraction (est.)", "95% interval"]
        known = "male_fraction_known" in x and x["male_fraction_known"].notna().any()
        if known:
            cols += ["male_fraction_known", "matches_known"]
            names_ += ["male fraction (known)", "estimate consistent"]
        tx = x[cols + ["pool_size"]].copy()
        tx.columns = names_ + ["X pool size"]
        if known and not x["matches_known"].fillna(True).astype(bool).all():
            verdict += ('<p class="status">' + ICON_BAD + "X depth does not match the known sex ratio in at least one "
                        "pool: check the pool composition, or X mappability</p>")
        parts.append("<h3>X chromosome</h3><p>Males carry one X, so in a pool with males the X has lower depth and "
                     "is sampled from fewer copies. With a fraction <em>m</em> of males, X depth / autosome depth = "
                     "1 − <em>m</em>/2 (0.75 for a 50:50 pool), and the X pool size is the pool size × (1 − <em>m</em>/2). "
                     "The male fraction is estimated from the depth ratio, or taken from <code>MALE_FRACTION</code> when the "
                     "sex ratio is known, with the estimate kept as a check. If any pool "
                     "differs from 1 (interval excluding it, by more than 5%), the X gets its own grenedalf run with "
                     "these pool sizes; otherwise the autosomal pool size would be wrong for the X.</p>"
                     + legend + svg + verdict + plain_table(tx))
    if len(m):
        tm = m[["pool", "median_depth", "ratio_to_A", "pool_size"]].copy()
        tm.columns = ["pool", "median depth", "× autosome depth", "pool size (individuals)"]
        parts.append("<h3>Mitochondrion</h3><p>Many copies per cell, haploid, and one lineage per individual: its "
                     "depth is far from the autosomes', and its pool size is the number of individuals, not genomes. "
                     f"It gets its own run, with depth bounds {cc.get('MITO_MIN_DEPTH', '?')}–{cc.get('MITO_MAX_DEPTH', '?')}× "
                     "and the whole genome as one window.</p>" + plain_table(tm))
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
    if runs is None:
        return ""
    assess = read_params(f"{P}.assess.params")
    strict = float(runs["artefact"].iloc[0] - runs["artefact_vs_strict"].iloc[0])
    tol = float(assess.get("TOLERANCE", "nan"))
    rec_mask = runs["recommended"].astype(str) == "True"
    ok_mask = runs["acceptable"].astype(str) == "True"
    colors = ["var(--s2)" if r else ("var(--s1)" if a else "var(--neutral)") for r, a in zip(rec_mask, ok_mask)]
    names = [r.replace("_", " ") for r in runs["run"]]
    base_bp = float(runs.loc[rec_mask, "callable_bp"].iloc[0]) if rec_mask.any() else float(runs["callable_bp"].max())
    rel = lambda v: "" if v == base_bp else f" ({100 * (v / base_bp - 1):+.1f}%)"
    left = hbars(names, runs["callable_bp"].astype(float).tolist(), lambda v: compact(v) + rel(v),
                 [f"{r.run}: {r.callable_bp:,} callable bp, {r.snps:,} SNPs" for r in runs.itertuples()],
                 "callable bp per grid run", colors=colors, ml=150)
    limit = strict + tol if np.isfinite(tol) else None
    right = hbars(names, runs["artefact"].astype(float).tolist(), lambda v: f"{v:.2f}",
                  [f"{r.run}: {r.artefact:.3f} ({r.worst}), {r.artefact_vs_strict:+.3f} vs strictest" for r in runs.itertuples()],
                  "artefact score per grid run", colors=colors, ml=150,
                  xmax=max(float(runs["artefact"].max()), (limit or 0)) * 1.1,
                  ref=(limit, f"accepted up to {limit:.2f}") if limit is not None else None)
    legend = ('<div class="legend"><span><span class="dotkey" style="background:var(--s2)"></span>recommended</span>'
              '<span><span class="dotkey" style="background:var(--s1)"></span>accepted (no added artefact)</span>'
              '<span><span class="dotkey" style="background:var(--neutral)"></span>rejected</span></div>')
    charts = (f'<div class="grid2"><div><h3>Data kept: callable bp on the pilot region</h3>{left}</div>'
              f'<div><h3>Artefact score: worst median |ρ| with a technical covariate</h3>{right}</div></div>')
    cov_cols = [c for c in runs.columns if ".rho_" in c]
    heat = runs[["run"] + cov_cols].copy()
    heat.columns = ["run"] + [c.replace(".rho_", " vs ").replace("_", " ") for c in cov_cols]
    tbl_cols = ["run", "callable_bp", "snps", "windows_scored", "artefact", "worst", "artefact_vs_strict",
                "acceptable", "recommended"]
    rec = runs.loc[runs["recommended"].astype(str) == "True", "run"]
    rec = rec.iloc[0] if len(rec) else None
    detail = ""
    if stats is not None:
        detail = ("<details><summary>Every statistic column (per pool / pair)</summary>"
                  + plain_table(stats) + "</details>")
    return ("<section><h2>Filter grid</h2>"
            "<p>The pilot region rerun with one setting changed at a time. For each run, each statistic "
            "(θπ, θW, Tajima's D, FST) is correlated across windows with read depth, depth imbalance between the "
            "two pools of a pair, and the callable fraction of the window; none of these should matter. The artefact "
            "score is the worst median |ρ|. Runs no worse than the strictest run (within the tolerance) are accepted, "
            "and the one keeping the most data is recommended.</p>"
            + legend + charts
            + "<h3>Median |ρ| with each technical covariate</h3>"
            + heat_table(heat, list(heat.columns[1:]), lambda v: f"{v:.2f}", vmax=max(float(np.nanmax(heat.iloc[:, 1:].to_numpy(dtype=float))), 0.1),
                         pick_row=lambda r: r["run"] == rec)
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
    sections = [section_summary(P, out_dir), section_mask(P), section_depth(P, dists), section_chromclass(P), section_min_count(P),
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
