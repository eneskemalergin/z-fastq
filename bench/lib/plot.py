"""Benchmark figures, light and dark SVG, from the flat files bench.sh writes.

    plot.py REPORT_DIR FIGURE_DIR HEADLINE_TIER SCALING_SOURCES MAIN_FAMILIES [FAMILY:PARAM,...] [SECTIONS]

REPORT_DIR holds summary.tsv, style.tsv, cases.tsv, and edges.tsv. Every figure
is written twice per theme: an SVG for development and the PNG the README
shows (200 dpi on the theme's page color, 256 colors, no dithering). Figures:
the edge-case matrix, then for each main family (comma-separated): an
isocost-style overview against z-fastq; wall time, peak RSS, and minor page
faults by kind of input (one row per format, one bar group per declared shape
group); and throughput and peak RSS against input size for the scaling
sources. Each FAMILY:PARAM adds wall time and peak RSS against the value of
the case parameter PARAM. A tool's thread variants (style.tsv: cores above
1, `of` naming the tool) show in the overview as their mean marks only, joined
to the tool's in thread order.

SECTIONS draws related main families together, one column per family:
"slug:Title=family/Column,family/Column|..." gives per section SLUG-overview,
SLUG-speed, SLUG-memory, and SLUG-scaling, rows by input format; a family in a
section gets no figures of its own.
"""

import csv
import io
import math
import multiprocessing
import sys
import textwrap
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, FuncFormatter, LogLocator, NullFormatter, NullLocator  # noqa: E402
from PIL import Image  # noqa: E402

plt.rcParams.update({"font.family": "DejaVu Sans", "hatch.linewidth": 0.7, "svg.hashsalt": "z-fastq"})

# The z-fasta asset theme: GitHub's page colors, with Zig orange as the accent.
THEMES = {
    "light": dict(bg="#ffffff", ink="#1f2328", ink2="#59636e", muted="#818b98", grid="#d1d9e0", region="#fff4e0"),
    "dark": dict(bg="#0d1117", ink="#f0f6fc", ink2="#9198a1", muted="#656c76", grid="#3d444d", region="#2a1f0c"),
}
MARKERS = {"circle": "o", "square": "s", "diamond": "D", "triangle": "^", "triangle-down": "v", "plus": "P", "star": "*", "cross": "X"}
SUBJECT = "z-fastq"
# Reserved status colors (dataviz reference palette); every cell also carries its text.
# Strength: how much of the color a cell shows, per theme; stronger for worse verdicts.
VERDICTS = [("right", "#0ca30c", "right", (0.30, 0.42)), ("refused", "#fab219", "refused a valid input", (0.50, 0.45)),
            ("crash", "#ec835a", "crashed", (0.62, 0.55)), ("wrong", "#d03b3b", "wrong count, exit 0", (0.55, 0.50)),
            ("policy", None, "tool policy", (0, 0))]
METRICS = {
    "speed": dict(summary="wall_ns", summary_scale=1e-6, axis="wall time (ms)", title="Wall time"),
    "memory": dict(summary="rss_kib", summary_scale=1 / 1024, axis="peak RSS (MiB)", title="Peak memory"),
    "faults": dict(summary="minor_faults", summary_scale=1, axis="minor page faults", title="Minor page faults"),
}
CLASS_NOTE = "Solid: full structure check; hatched: partial check."


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def gmean(values):
    return math.exp(sum(math.log(v) for v in values) / len(values))


def t95(df):
    """Two-sided 95% quantile of Student's t."""
    table = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57, 6: 2.45, 7: 2.36, 8: 2.31, 9: 2.26, 10: 2.23,
             12: 2.18, 15: 2.13, 20: 2.09, 30: 2.04}
    return table.get(df) or table[max(k for k in table if k <= df)] if df <= 30 else 1.96


def mix(color, bg, weight):
    """`weight` of `color` over `bg`, as an opaque color (hatches stay crisp)."""
    return tuple(weight * a + (1 - weight) * b for a, b in zip(to_rgb(color), to_rgb(bg)))


def times(r):
    """A ratio as a short multiple: 0.62x, 1.4x, 24.5x, 310x."""
    if r >= 100:
        return f"{r:.0f}x"
    if r >= 0.995:
        return f"{r:.1f}x"
    return f"{r:.2g}x"


def count_axis(axis, unit="", dense=False):
    """Log axis with 1-2-5 ticks (1-1.5-2-3-5-7 when `dense`): 0.5, 2, 50, 2k, 1M."""
    def fmt(v, _):
        for div, suffix in ((1e6, "M"), (1e3, "k")):
            if v >= div:
                return f"{v / div:g}{suffix}{unit}"
        return f"{v:g}{unit}"
    axis.set_major_locator(LogLocator(base=10, subs=(1.0, 1.5, 2.0, 3.0, 5.0, 7.0) if dense else (1.0, 2.0, 5.0)))
    axis.set_major_formatter(FuncFormatter(fmt))
    axis.set_minor_formatter(NullFormatter())
    axis.set_minor_locator(NullLocator())


class Figures:
    def __init__(self, report, out, headline, scaling, mains, lines="", sections=""):
        self.report = Path(report)
        self.out = Path(out)
        self.out.mkdir(parents=True, exist_ok=True)
        self.headline = headline
        self.scaling = scaling.split()
        self.mains = mains.split(",")
        self.lines = [tuple(x.split(":", 1)) for x in lines.split(",") if x]
        self.sections = []
        for spec in filter(None, sections.split("|")):
            head, cols = spec.split("=", 1)
            slug, title = head.split(":", 1)
            self.sections.append((slug, title, [tuple(c.split("/", 1)) for c in cols.split(",")]))
        # The feature's phrases: the wrong-answer legend and the past tense of its job.
        self.words = dict(line.split("\t", 1) for line in (self.report / "words.tsv").read_text().splitlines())
        self.style = {r["lane"]: r for r in read(self.report / "style.tsv")}
        for name, r in self.style.items():
            r["cores"] = int(r.get("cores") or 1)
            r["of"] = r.get("of") or name
        self.order = sorted(self.style, key=lambda n: int(self.style[n]["order"]))
        self.cases = {(r["family"], r["case"]): r for r in read(self.report / "cases.tsv")}
        # Dots mark the input behind each value: a circle for a real run, one marker per synthetic shape.
        shapes = [i["source"] for i in sorted(self.cases.values(), key=lambda i: int(i["order"])) if i["kind"] == "shape"]
        self.dot = {src: m for src, m in zip(dict.fromkeys(shapes), ("s", "^", "D", "v", "P", "X"))}
        self.rows = {}
        for r in read(self.report / "summary.tsv"):
            self.rows.setdefault((r["family"], r["case"]), {})[r["lane"]] = r

    # ---- shared frame -------------------------------------------------------

    def frame(self, theme, nrows=1, ncols=1, width=8.0, height=5.0, grid="both", **kw):
        t = THEMES[theme]
        fig, axes = plt.subplots(nrows, ncols, figsize=(width, height), squeeze=False, **kw)
        fig.patch.set_facecolor(t["bg"])
        for ax in axes.flat:
            ax.set_facecolor(t["bg"])
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            for side in ("left", "bottom"):
                ax.spines[side].set_color(t["muted"])
                ax.spines[side].set_linewidth(0.8)
            ax.tick_params(colors=t["ink2"], labelsize=8.5, length=3, width=0.7)
            if grid:
                ax.grid(True, axis=grid, which="major", color=t["grid"], linewidth=0.6, alpha=0.8)
            ax.set_axisbelow(True)
        return fig, axes

    def finish(self, fig, stem, theme, title, subtitle, handles=None, ncol=None, keys=None):
        """Title and italic note centered on top, legend centered below, tight bounds."""
        t = THEMES[theme]
        fig.tight_layout(h_pad=1.6, w_pad=1.6)
        boxes = [ax.get_tightbbox() for ax in fig.axes]
        fh = fig.bbox.height
        top = max(b.y1 for b in boxes) / fh
        bottom = min(b.y0 for b in boxes) / fh
        inch = 1 / fig.get_figheight()
        lines = textwrap.wrap(subtitle, int(fig.get_figwidth() * 15))
        fig.text(0.5, top + 0.16 * inch, "\n".join(lines), ha="center", va="bottom", fontsize=9, color=t["ink2"],
                 style="italic", linespacing=1.45)
        fig.text(0.5, top + (0.26 + 0.2 * len(lines)) * inch, title, ha="center", va="bottom", fontsize=13,
                 color=t["ink"], fontweight="bold")
        if handles:
            fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, bottom - 0.06 * inch),
                       ncol=ncol or min(len(handles), 8), frameon=False, fontsize=9, labelcolor=t["ink"],
                       handlelength=1.8, handleheight=1.0, columnspacing=1.5, borderaxespad=0)
            if keys:
                fig.legend(handles=keys, loc="upper center", bbox_to_anchor=(0.5, bottom - 0.4 * inch), ncol=len(keys),
                           frameon=False, fontsize=8.5, labelcolor=t["ink2"], handlelength=1.0, columnspacing=1.3,
                           borderaxespad=0, handletextpad=0.4)
        path = self.out / f"{stem}-{theme}.svg"
        part = path.with_suffix(".svg.part")
        fig.savefig(part, format="svg", facecolor=t["bg"], bbox_inches="tight", pad_inches=0.2)
        part.replace(path)
        # The README's copy: a flat chart needs few colors, so a 256-color
        # palette without dithering keeps every edge and halves the file.
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=200, facecolor=t["bg"], bbox_inches="tight", pad_inches=0.2)
        png = path.with_suffix(".png")
        image = Image.open(buf).convert("RGB").quantize(256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
        image.save(png.with_suffix(".png.part"), format="png", optimize=True)
        png.with_suffix(".png.part").replace(png)
        plt.close(fig)

    def bar_style(self, lane, theme):
        st, bg = self.style[lane], THEMES[theme]["bg"]
        color, cls = st[theme], st["class"]
        if cls == "partial":
            return dict(facecolor=mix(color, bg, 0.30), edgecolor=color, hatch="////")
        return dict(facecolor=color, edgecolor=color, hatch=None)

    def dot_keys(self, theme, family=None):
        """Marker keys for the inputs a figure of FAMILY shows: real runs, then its shapes."""
        t = THEMES[theme]
        present = {i["source"] for i in self.cases.values() if family is None or i["family"] == family}
        keys = [("o", "real run")] + [(m, src) for src, m in self.dot.items() if src in present]
        return [Line2D([], [], marker=m, linestyle="none", markersize=5 if m == "o" else 5.5, markerfacecolor=t["bg"],
                       markeredgecolor=t["ink"], markeredgewidth=0.8, label=text) for m, text in keys]

    def bar_handles(self, lanes, theme):
        return [Patch(label=self.style[n]["label"], linewidth=0.8, **self.bar_style(n, theme)) for n in lanes]

    def mark(self, ax, lane, x, y, theme, size=9, label=True):
        """One marker per lane: class decides fill and size, the registry decides color and shape."""
        st, t = self.style[lane], THEMES[theme]
        color, cls = st[theme], st["class"]
        face = "none" if cls == "partial" else color
        edge = t["ink"] if cls == "subject" else color
        ax.plot(x, y, MARKERS[st["shape"]], markersize=size * (1.3 if cls == "subject" else 1),
                markerfacecolor=face, markeredgecolor=edge, markeredgewidth=1.0 if cls == "subject" else 1.5,
                linestyle="none", label=st["label"] if label else None, zorder=5 if cls == "subject" else 4)

    def mark_handles(self, lanes, theme):
        t = THEMES[theme]
        out = []
        for n in lanes:
            st = self.style[n]
            cls, color = st["class"], st[theme]
            out.append(Line2D([], [], marker=MARKERS[st["shape"]], linestyle="none", label=st["label"],
                              markersize=9 if cls == "subject" else 7.5,
                              markerfacecolor="none" if cls == "partial" else color,
                              markeredgecolor=t["ink"] if cls == "subject" else color, markeredgewidth=1.0 if cls == "subject" else 1.5))
        return out

    # ---- data ---------------------------------------------------------------

    def value(self, family, case, lane, metric):
        """One shape's value for a lane: the median of its runs."""
        spec = METRICS[metric]
        row = self.rows.get((family, case), {}).get(lane)
        if not row:
            return None
        v = float(row[spec["summary"]]) * spec["summary_scale"]
        return v if v > 0 else None

    def tier_cases(self, family, fmt, tier=None):
        tier = tier or self.headline
        out = [i for k, i in self.cases.items() if i["family"] == family and i["format"] == fmt and i["tier"] == tier
               and i["kind"] != "many" and k in self.rows and SUBJECT in self.rows[k]]
        return sorted(out, key=lambda i: int(i["order"]))

    def groups(self, family, fmt):
        """The headline shapes of one format, by declared group: [(label, [case, ...]), ...]."""
        out = {}
        for i in self.tier_cases(family, fmt):
            out.setdefault((int(i["group_order"]), i["group"]), []).append(i)
        return [(label, cases) for (_, label), cases in sorted(out.items())]

    def formats(self, family):
        return [f for f in ("plain", "gzip") if self.tier_cases(family, f)]

    # ---- grouped bars -------------------------------------------------------

    def group_bars(self, ax, groups, lanes, shape_values, theme):
        """One group per shape group, one bar per lane in a fixed slot. `shape_values(case, lane)`
        gives one value per shape; the bar is their geometric mean, the dots the shapes, the
        whisker their range (shapes differ by design, so an interval would not mean much). The
        box is the geometric mean of the lane's ratio to z-fastq over the shapes both ran."""
        t = THEMES[theme]
        n = len(lanes)
        width = 0.84 / n
        labels, lo_all, hi_all = [], math.inf, 0.0
        for gi, (_, cases) in enumerate(groups):
            spots = []
            for j, lane in enumerate(lanes):
                pairs = [(shape_values(c, lane), shape_values(c, SUBJECT)) for c in cases]
                vals = [v for v, _ in pairs if v]
                if not vals:
                    continue
                x = gi + (j - (n - 1) / 2) * width
                g, lo, hi = gmean(vals), min(vals), max(vals)
                ax.bar(x, g, width=width * 0.86, linewidth=0.8, zorder=2, **self.bar_style(lane, theme))
                if len(vals) > 1:
                    ax.errorbar(x, g, yerr=[[max(0.0, g - lo)], [max(0.0, hi - g)]], fmt="none", ecolor=t["ink"], elinewidth=0.8,
                                capsize=1.8, capthick=0.8, zorder=4)
                    xs = [x + width * 0.22 * (((k * 0.6180339) % 1) * 2 - 1) for k in range(len(vals))]
                    marks = [self.dot.get(c["source"], "o") for c, (v, _) in zip(cases, pairs) if v]
                    for m in dict.fromkeys(marks):
                        pick = [k for k, mk in enumerate(marks) if mk == m]
                        ax.scatter([xs[k] for k in pick], [vals[k] for k in pick], s=10 if m == "o" else 14, marker=m,
                                   facecolor=t["bg"], edgecolor=t["ink"], linewidths=0.7, alpha=0.95, zorder=5)
                ratios = [v / b for v, b in pairs if v and b]
                top = max(hi, max(vals))
                spots.append((x, top, lane))
                lo_all, hi_all = min(lo_all, lo), max(hi_all, top)
                if ratios:
                    labels.append((x, top, times(gmean(ratios)), lane))
            if any(s[2] == SUBJECT for s in spots):
                ref = gmean([v for v in (shape_values(c, SUBJECT) for c in cases) if v])
                ax.hlines(ref, spots[0][0] - width * 0.5, spots[-1][0] + width * 0.5, colors=self.style[SUBJECT][theme],
                          linestyles=(0, (4, 3)), linewidth=0.9, alpha=0.75, zorder=1.5)
        for x, top, text, lane in labels:
            color = self.style[lane][theme]
            ax.annotate(text, (x, top), xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
                        rotation=90, fontsize=6.8, fontweight="bold", color=color if lane == SUBJECT else t["ink"], zorder=6,
                        bbox=dict(boxstyle="round,pad=0.22", facecolor=t["bg"], edgecolor=color, linewidth=0.8))
        ax.set_yscale("log")
        if hi_all > 0:
            # Room above the tallest bar for its box: about 44 points of the axis height.
            ax_pts = ax.get_position().height * ax.figure.get_figheight() * 72
            frac = min(0.4, 44 / max(ax_pts, 60))
            a, b = math.log(lo_all / 1.6), math.log(hi_all)
            ax.set_ylim(math.exp(a), math.exp(a + (b - a) / (1 - frac)))
        ax.set_xlim(-0.5, len(groups) - 0.5)
        ax.set_xticks(range(len(groups)))
        synthetic = [all(c["kind"] == "shape" for c in cases) for _, cases in groups]
        ax.set_xticklabels([f"{label}\n{len(cases)} {('shape' if syn else 'real run') + ('' if len(cases) == 1 else 's')}"
                            for (label, cases), syn in zip(groups, synthetic)], fontsize=8.5, color=t["ink"], linespacing=1.3)
        ax.tick_params(axis="x", length=0)
        for gi in range(1, len(groups)):
            # A stronger rule where the pooled real runs end and the single synthetic shapes begin.
            edge = synthetic[gi] and not synthetic[gi - 1]
            ax.axvline(gi - 0.5, color=t["ink2"] if edge else t["grid"], linewidth=1.0 if edge else 0.6,
                       alpha=0.8 if edge else 0.6, zorder=0)
        count_axis(ax.yaxis)

    def bar_panels(self, family, lanes_for, shape_values, theme, share=True):
        """One row per format; returns the figure and the lanes shown."""
        fmts = [f for f in self.formats(family) if self.groups(family, f) and lanes_for(f)]
        if not fmts:
            return None, []
        panels = [(f, self.groups(family, f)) for f in fmts]
        shown = [n for n in self.order if any(n in lanes_for(f) for f in fmts)]
        widest = max(len(g) * (0.2 * len(lanes_for(f)) + 0.5) for f, g in panels)
        fig, axes = self.frame(theme, len(panels), 1, width=max(8.5, widest + 1.4), height=3.5 * len(panels) + 0.3, grid="y")
        for (fmt, groups), ax in zip(panels, axes[:, 0]):
            lanes = [n for n in self.order if n in lanes_for(fmt)]
            self.group_bars(ax, groups, lanes, lambda c, n, f=fmt: shape_values(f, c, n), theme)
            ax.set_title(f"{fmt} input", loc="left", fontsize=10, color=THEMES[theme]["ink"], fontweight="bold", pad=8)
        if share:
            lo = min(ax.get_ylim()[0] for ax in axes[:, 0])
            hi = max(ax.get_ylim()[1] for ax in axes[:, 0])
            for ax in axes[:, 0]:
                ax.set_ylim(lo, hi)
        return fig, shown

    def metric_bars(self, family, metric, theme):
        spec = METRICS[metric]

        def lanes_for(fmt):
            return {n for i in self.tier_cases(family, fmt) for n in self.rows[(family, i["case"])]
                    if self.value(family, i["case"], n, metric)}

        # Time rows share one scale so plain and gzip compare; memory and faults rows fit their own data.
        fig, shown = self.bar_panels(family, lanes_for, lambda f, c, n: self.value(family, c["case"], n, metric), theme,
                                              share=metric == "speed")
        if not fig:
            return False
        for ax in fig.axes:
            ax.set_ylabel(f"{spec['axis']}, log scale", fontsize=9, color=THEMES[theme]["ink2"])
        lower = "Lower is better. " if metric != "faults" else ""
        self.finish(fig, f"{family}-{metric}", theme, f"{spec['title']} by kind of input",
                    f"{lower}{self.headline.capitalize()} tier, about {self.words['size']} of FASTQ per input. Bars: geometric "
                    "mean over the inputs of a group (dots, each the median of its timed repeats; markers name the synthetic "
                    "shapes), whiskers their range. Boxes: multiple of "
                    f"z-fastq (dashed line). {CLASS_NOTE}", handles=self.bar_handles(shown, theme), keys=self.dot_keys(theme, family))
        return True

    # ---- overview -----------------------------------------------------------

    def overview(self, family, theme):
        fmts = self.formats(family)
        if not fmts:
            return False
        t = THEMES[theme]
        fig, axes = self.frame(theme, 1, len(fmts), width=5.2 * len(fmts) + 0.6, height=5.4)
        shown, nshapes = [], {}
        self.front_drawn = False
        for ax, fmt in zip(axes[0], fmts):
            data = {}
            cases = self.tier_cases(family, fmt)
            nshapes[fmt] = len(cases)
            for info in cases:
                ref_t = self.value(family, info["case"], SUBJECT, "speed")
                ref_m = self.value(family, info["case"], SUBJECT, "memory")
                for lane in self.rows[(family, info["case"])]:
                    lt, lm = self.value(family, info["case"], lane, "speed"), self.value(family, info["case"], lane, "memory")
                    if lt and lm and ref_t and ref_m:
                        data.setdefault(lane, []).append((lt / ref_t, lm / ref_m))
            self.isocost(ax, data, theme)
            ax.set_title(f"{fmt} input, {len(cases)} inputs", loc="left", fontsize=10, color=t["ink"], fontweight="bold", pad=8)
            shown += [n for n in data if n not in shown]
        self.finish(fig, f"{family}-overview", theme, "Time and memory relative to z-fastq",
                    f"Large marks: geometric mean over the {self.headline}-tier inputs; faint marks: single inputs; "
                    "whiskers: 95% interval of the mean; ellipse: where about half the inputs fall. Dashed curves: "
                    "equal cost C = (time + memory) / 2." + (" Solid line: tools nothing beats on both." if self.front_drawn else ""),
                    handles=self.mark_handles([n for n in self.order if n in shown], theme), ncol=len(shown))
        return True

    def isocost(self, ax, data, theme):
        """Per lane: faint marks for the shapes, the geometric mean as a large mark, whiskers for
        the 95% interval of that mean on each axis, and an ellipse holding about half the inputs.
        Equal-cost curves are linear cost (time ratio + RSS ratio) / 2, curved on log axes."""
        t = THEMES[theme]
        fig = ax.figure
        ax.set_xscale("log")
        ax.set_yscale("log")
        points, extent = {}, [(1.0, 1.0)]
        for lane in self.order:
            if lane not in data:
                continue
            lx = np.log([p[0] for p in data[lane]])
            ly = np.log([p[1] for p in data[lane]])
            n = len(lx)
            gx, gy = math.exp(lx.mean()), math.exp(ly.mean())
            points[lane] = (gx, gy)
            extent.append((gx, gy))
            st = self.style[lane]
            color = st[theme]
            variant = st["of"] != lane
            if lane != SUBJECT and not variant:
                ax.plot(np.exp(lx), np.exp(ly), MARKERS[st["shape"]], markersize=4.5, linestyle="none", alpha=0.35,
                        markerfacecolor=color, markeredgecolor="none", zorder=2)
                extent += list(zip(np.exp(lx), np.exp(ly)))
            if n >= 2 and lane != SUBJECT and not variant:
                q = t95(n - 1) / math.sqrt(n)
                hx, hy = q * lx.std(ddof=1), q * ly.std(ddof=1)
                ax.errorbar(gx, gy, xerr=[[gx - gx * math.exp(-hx)], [gx * math.exp(hx) - gx]],
                            yerr=[[gy - gy * math.exp(-hy)], [gy * math.exp(hy) - gy]], fmt="none", ecolor=color,
                            elinewidth=1.1, capsize=2.5, capthick=1.1, alpha=0.9, zorder=3)
            if n >= 4 and lane != SUBJECT and not variant:
                cov = np.cov(np.vstack([lx, ly]))
                vals, vecs = np.linalg.eigh(cov)
                if vals.min() > 1e-12:
                    # Radius sqrt(chi2(2, 0.5)) covers about half of a normal spread.
                    angle = np.linspace(0, 2 * math.pi, 120)
                    circle = np.vstack([np.cos(angle), np.sin(angle)]) * 1.1774
                    ring = vecs @ (np.sqrt(vals)[:, None] * circle)
                    ex, ey = np.exp(lx.mean() + ring[0]), np.exp(ly.mean() + ring[1])
                    ax.fill(ex, ey, facecolor=color, alpha=0.08, edgecolor="none", zorder=1.2)
                    ax.plot(ex, ey, color=color, linewidth=0.8, alpha=0.55, zorder=1.3)
                    extent += [(ex.min(), ey.min()), (ex.max(), ey.max())]
            self.mark(ax, lane, gx, gy, theme, label=False)
        # Each tool's path through its thread counts, from its own mark.
        for lane in self.order:
            path = sorted((self.style[n]["cores"], points[n]) for n in points if self.style[n]["of"] == lane)
            if len(path) > 1:
                ax.plot([q[1][0] for q in path], [q[1][1] for q in path], color=self.style[lane][theme], linewidth=1.0,
                        alpha=0.8, zorder=2.4)
        xs, ys = [p[0] for p in extent], [p[1] for p in extent]
        lo_x, hi_x = min(xs) / 1.35, max(xs) * 1.35
        lo_y, hi_y = min(ys) / 1.35, max(ys) * 1.35
        ax.set_xlim(lo_x, hi_x)
        ax.set_ylim(lo_y, hi_y)
        ax.fill_between([lo_x, 1], lo_y, 1, color=t["region"], zorder=0, linewidth=0)
        fixed = []
        # Name the tinted corner only where it has room for the words.
        if math.log(1 / lo_x) > 0.28 * math.log(hi_x / lo_x) and math.log(1 / lo_y) > 0.12 * math.log(hi_y / lo_y):
            ax.text(lo_x * 1.06, lo_y * 1.08, "faster and smaller\nthan z-fastq", fontsize=8, color=t["ink2"], style="italic",
                    va="bottom", ha="left")
            fixed.append(((lo_x * 1.06, lo_y * 1.08), 90, 24))
        for line in (ax.axvline, ax.axhline):
            line(1, color=self.style[SUBJECT][theme], linewidth=0.9, linestyle=(0, (4, 3)), alpha=0.7, zorder=1)
        # Equal-cost curves through readable levels between the cheapest and dearest corner.
        c_lo, c_hi = (lo_x + lo_y) / 2, (hi_x + hi_y) / 2
        steps = [m * 10 ** e for e in range(-2, 4) for m in (1, 1.5, 2, 3, 5, 7)]
        levels = [c for c in steps if c_lo < c < c_hi]
        if len(levels) > 6:
            pick = [levels[round(i * (len(levels) - 1) / 5)] for i in range(6)]
            levels = sorted(set(pick) | ({1} & set(levels)))
        grid_x = np.exp(np.linspace(math.log(lo_x), math.log(hi_x), 400))
        for c in levels:
            gy = 2 * c - grid_x
            keep = (gy > lo_y) & (gy < hi_y)
            if keep.sum() < 2:
                continue
            ax.plot(grid_x[keep], gy[keep], color=t["muted"], linewidth=0.7, linestyle=(0, (4, 3)), alpha=0.8, zorder=1)
            # Label at the curve's upper-left end, inside the axes with room for the text.
            room = keep & (grid_x > lo_x * 1.06) & (grid_x < hi_x / 1.45) & (gy > lo_y * 1.1) & (gy < hi_y / 1.1)
            if not room.any():
                continue
            i = np.flatnonzero(room)[0]
            lx, ly = grid_x[i], gy[i]
            fixed.append(((lx, ly - 0), 6 * len(f"C={c:g}") + 4, 10))
            ax.text(lx, ly, f"C={c:g}", fontsize=7, color=t["muted"], ha="left", va="center",
                    bbox=dict(boxstyle="round,pad=0.12", facecolor=t["bg"], edgecolor="none"), zorder=1.4)
        front = sorted(p for p in points.values() if not any(q[0] <= p[0] and q[1] <= p[1] and q != p for q in points.values()))
        if len(front) > 1:
            ax.plot([p[0] for p in front], [p[1] for p in front], color=t["ink2"], linewidth=1.0, zorder=2.5)
            self.front_drawn = True
        # Direct labels: the first of eight spots around the marker that overlaps no mark or label.
        scale = 72 / fig.dpi
        centers = {n: ax.transData.transform(p) * scale for n, p in points.items()}
        taken = [(x - 7, y - 7, x + 7, y + 7) for x, y in centers.values()]
        for (fx, fy), w, h in fixed:
            x, y = ax.transData.transform((fx, fy)) * scale
            taken.append((x, y - h / 2 if h == 10 else y, x + w, y + h / 2 if h == 10 else y + h))
        for lane in self.order:
            if lane not in points:
                continue
            text = self.style[lane]["label"].replace("_", " ")
            base = self.style[lane]["of"]
            # A thread variant is named by its option alone (-w 16); its tool's mark carries the name.
            if base != lane and base in self.style:
                text = text.replace(self.style[base]["label"].replace("_", " "), "", 1).strip()
            w, h = 5.4 * len(text) + 2, 11
            x, y = centers[lane]
            for dx, dy in ((9, 4), (9, -h - 4), (-w - 9, 4), (-w - 9, -h - 4), (-w / 2, 10), (-w / 2, -h - 10), (9, -h / 2), (-w - 9, -h / 2)):
                box = (x + dx, y + dy, x + dx + w, y + dy + h)
                if not any(box[0] < b[2] and b[0] < box[2] and box[1] < b[3] and b[1] < box[3] for b in taken):
                    break
            else:
                dx, dy = 9, 4
                box = (x + dx, y + dy, x + dx + w, y + dy + h)
            taken.append(box)
            ax.annotate(text, points[lane], xytext=(dx, dy), textcoords="offset points", fontsize=8.5, color=t["ink"],
                        va="bottom", ha="left", fontweight="bold" if lane == SUBJECT else "normal", zorder=6,
                        bbox=dict(boxstyle="square,pad=0.1", facecolor=t["bg"], edgecolor="none", alpha=0.85))
        count_axis(ax.xaxis, "x", dense=hi_x / lo_x < 12)
        count_axis(ax.yaxis, "x", dense=hi_y / lo_y < 12)
        ax.set_xlabel("wall time relative to z-fastq, log scale", color=t["ink2"], fontsize=9)
        ax.set_ylabel("peak RSS relative to z-fastq, log scale", color=t["ink2"], fontsize=9)

    # ---- sections: related families side by side ----------------------------

    def section_formats(self, cols):
        return [f for f in ("plain", "gzip") if any(self.tier_cases(fam, f) for fam, _ in cols)]

    def overview_data(self, family, fmt):
        """Per lane, (time, memory) over z-fastq's for each headline input of FAMILY in FMT."""
        data = {}
        cases = self.tier_cases(family, fmt)
        for info in cases:
            ref_t = self.value(family, info["case"], SUBJECT, "speed")
            ref_m = self.value(family, info["case"], SUBJECT, "memory")
            for lane in self.rows[(family, info["case"])]:
                lt, lm = self.value(family, info["case"], lane, "speed"), self.value(family, info["case"], lane, "memory")
                if lt and lm and ref_t and ref_m:
                    data.setdefault(lane, []).append((lt / ref_t, lm / ref_m))
        return data, len(cases)

    def section_overview(self, slug, title, cols, theme):
        fmts = self.section_formats(cols)
        if not fmts:
            return False
        t = THEMES[theme]
        fig, axes = self.frame(theme, len(fmts), len(cols), width=5.0 * len(cols) + 0.6, height=4.7 * len(fmts) + 0.4)
        shown = []
        self.front_drawn = False
        for r, fmt in enumerate(fmts):
            for c, (family, name) in enumerate(cols):
                ax = axes[r, c]
                data, n = self.overview_data(family, fmt)
                if not data:
                    ax.set_axis_off()
                    continue
                self.isocost(ax, data, theme)
                ax.set_title(f"{name}, {fmt} input, {n} input{'s' if n != 1 else ''}", loc="left", fontsize=10, color=t["ink"],
                             fontweight="bold", pad=8)
                shown += [x for x in data if x not in shown]
        self.finish(fig, f"{slug}-overview", theme, f"{title}: time and memory relative to z-fastq",
                    f"Large marks: geometric mean over the {self.headline}-tier inputs; faint marks: single inputs; whiskers: 95% "
                    "interval of the mean; ellipse: where about half the inputs fall. A tool's thread counts are joined by a line "
                    "from its one-thread mark. Dashed curves: equal cost C = (time + memory) / 2."
                    + (" Solid line: tools nothing beats on both." if self.front_drawn else ""),
                    handles=self.mark_handles([x for x in self.order if x in shown], theme), ncol=min(len(shown), 6))
        return True

    def section_bars(self, slug, title, cols, metric, theme):
        spec = METRICS[metric]
        fmts = self.section_formats(cols)
        panels = [(fmt, family, name) for fmt in fmts for family, name in cols]
        lanes_of = {}
        for fmt, family, _ in panels:
            lanes_of[(fmt, family)] = [x for x in self.order if any(self.value(family, i["case"], x, metric) for i in self.tier_cases(family, fmt))]
        if not any(lanes_of.values()):
            return False
        t = THEMES[theme]
        widest = max(len(self.groups(fam, fmt)) * (0.2 * len(lanes_of[(fmt, fam)]) + 0.5) for fmt, fam, _ in panels)
        fig, axes = self.frame(theme, len(fmts), len(cols), width=len(cols) * max(4.6, widest + 0.6) + 0.8,
                               height=3.4 * len(fmts) + 0.4, grid="y")
        shown = []
        for r, fmt in enumerate(fmts):
            for c, (family, name) in enumerate(cols):
                ax = axes[r, c]
                lanes = lanes_of[(fmt, family)]
                groups = self.groups(family, fmt)
                if not lanes or not groups:
                    ax.set_axis_off()
                    continue
                self.group_bars(ax, groups, lanes, lambda case, x, f=family: self.value(f, case["case"], x, metric), theme)
                ax.set_title(f"{name}, {fmt} input", loc="left", fontsize=10, color=t["ink"], fontweight="bold", pad=8)
                if c == 0:
                    ax.set_ylabel(f"{spec['axis']}, log scale", fontsize=9, color=t["ink2"])
                shown += [x for x in lanes if x not in shown]
        # One scale for every panel, so jobs and formats compare at a glance.
        live = [ax for ax in axes.flat if ax.axison]
        lo, hi = min(ax.get_ylim()[0] for ax in live), max(ax.get_ylim()[1] for ax in live)
        for ax in live:
            ax.set_ylim(lo, hi)
        self.finish(fig, f"{slug}-{metric}", theme, f"{title}: {spec['title'].lower()} by kind of input",
                    f"Lower is better. {self.headline.capitalize()} tier, about {self.words['size']} of FASTQ per input. Bars: "
                    "geometric mean over the inputs of a group (dots, each the median of its timed repeats), whiskers their range. "
                    f"Boxes: multiple of z-fastq (dashed line). {CLASS_NOTE}",
                    handles=self.bar_handles([x for x in self.order if x in shown], theme), ncol=min(len(shown), 6),
                    keys=self.dot_keys(theme, cols[0][0]))
        return True

    def section_scaling(self, slug, title, cols, theme):
        """Rows: formats; columns: per family, throughput then peak RSS against input size; one line per lane."""
        src = next((x for x in self.scaling if any(i["source"] == x and i["family"] == cols[0][0] for i in self.cases.values())), None)
        if not src:
            return False
        grid = {}
        for (fam, case), lanes in self.rows.items():
            info = self.cases.get((fam, case))
            if not info or info["source"] != src or fam not in [f for f, _ in cols]:
                continue
            mib = int(info["bytes"]) / 2**20
            for lane in lanes:
                tv, m = self.value(fam, case, lane, "speed"), self.value(fam, case, lane, "memory")
                if tv:
                    grid.setdefault((fam, info["format"], "speed"), {}).setdefault(lane, []).append((mib, mib / (tv / 1e3)))
                if m:
                    grid.setdefault((fam, info["format"], "memory"), {}).setdefault(lane, []).append((mib, m))
        fmts = [f for f in ("plain", "gzip") if any(len(v.get(SUBJECT, [])) > 1 for (_, ff, _), v in grid.items() if ff == f)]
        if not fmts:
            return False
        t = THEMES[theme]
        panels = [(fam, name, metric) for fam, name in cols for metric in ("speed", "memory")]
        fig, axes = self.frame(theme, len(fmts), len(panels), width=3.3 * len(panels) + 0.8, height=2.8 * len(fmts) + 0.4)
        sizes = [(1, "1 MiB"), (10, "10 MiB"), (100, "100 MiB"), (1024, "1 GiB")]
        xs = [q[0] for v in grid.values() for pts in v.values() for q in pts]
        sizes = [x for x in sizes if x[0] <= max(xs) * 1.5]
        shown = []
        for r, fmt in enumerate(fmts):
            for c, (fam, name, metric) in enumerate(panels):
                ax = axes[r, c]
                ax.set_xscale("log")
                ax.set_yscale("log")
                series = grid.get((fam, fmt, metric), {})
                for lane in self.order:
                    if len(series.get(lane, [])) < 2:
                        continue
                    pts = sorted(series[lane])
                    st = self.style[lane]
                    ax.plot([q[0] for q in pts], [q[1] for q in pts], color=st[theme], linestyle=(0, (4, 2)) if st["class"] == "partial" else "-",
                            linewidth=2.2 if st["class"] == "subject" else 1.2, zorder=3 if lane == SUBJECT else 2)
                    for q in pts:
                        self.mark(ax, lane, q[0], q[1], theme, size=4.5, label=False)
                    if lane not in shown:
                        shown.append(lane)
                ax.xaxis.set_major_locator(FixedLocator([v for v, _ in sizes]))
                ax.xaxis.set_major_formatter(FixedFormatter([n for _, n in sizes]))
                ax.xaxis.set_minor_locator(NullLocator())
                ax.tick_params(axis="x", labelsize=7.5)
                count_axis(ax.yaxis)
                ax.margins(y=0.12)
                if r == 0:
                    ax.set_title(f"{name}: {'throughput' if metric == 'speed' else 'peak memory'}", fontsize=9.5, color=t["ink"], pad=6,
                                 fontweight="bold")
                if r == len(fmts) - 1:
                    ax.set_xlabel("FASTQ size, log scale", fontsize=8.5, color=t["ink2"])
                unit = "MiB/s" if metric == "speed" else "peak RSS (MiB)"
                ax.set_ylabel(f"{fmt} input\n{unit}, log scale" if c == 0 else f"{unit}, log scale", fontsize=8.5, color=t["ink2"])
        label = next(i["label"] for i in self.cases.values() if i["source"] == src).split(": ", 1)[-1]
        self.finish(fig, f"{slug}-scaling", theme, f"{title}: throughput and peak memory by input size",
                    f"Throughput: higher is better, MiB of FASTQ per second of wall time; at 1 MiB, starting the process dominates. "
                    f"Peak memory: lower is better. The {src} shape ({label}); each point is the median of its runs; the tools at one "
                    "thread. Dashed: partial structure check.",
                    handles=self.mark_handles([x for x in self.order if x in shown], theme), ncol=min(len(shown), 6))
        return True

    # ---- scaling ------------------------------------------------------------

    def scaling_figure(self, family, metrics, theme):
        """Rows: formats; columns: each scaling source, once per metric (throughput,
        peak RSS); x: input size; one line per lane. With one source, both
        metrics share one figure side by side."""
        grid = {}
        for (fam, case), lanes in self.rows.items():
            info = self.cases.get((fam, case))
            if fam != family or not info or info["source"] not in self.scaling or info["kind"] == "many":
                continue
            mib = int(info["bytes"]) / 2**20
            for lane in lanes:
                t, m = self.value(family, case, lane, "speed"), self.value(family, case, lane, "memory")
                if t:
                    grid.setdefault((info["format"], info["source"], "speed"), {}).setdefault(lane, []).append((mib, mib / (t / 1e3)))
                if m:
                    grid.setdefault((info["format"], info["source"], "memory"), {}).setdefault(lane, []).append((mib, m))
        fmts = [f for f in ("plain", "gzip") if any(len(s.get(SUBJECT, [])) > 1 for (ff, _, _), s in grid.items() if ff == f)]
        sources = [s for s in self.scaling if any((f, s, metrics[0]) in grid for f in fmts)]
        if not fmts or not sources:
            return False
        cols = [(src, metric) for src in sources for metric in metrics]
        t = THEMES[theme]
        both = len(metrics) > 1
        fig, axes = self.frame(theme, len(fmts), len(cols), width=(3.6 if both else 2.5) * len(cols) + 0.8,
                               height=2.8 * len(fmts) + 0.4, sharex=True, sharey=False if both else True)
        sizes = [(1, "1 MiB"), (10, "10 MiB"), (100, "100 MiB"), (1024, "1 GiB"), (10240, "10 GiB")]
        xs = [p[0] for ser in grid.values() for pts in ser.values() for p in pts]
        sizes = [s for s in sizes if s[0] <= max(xs) * 1.5]
        shown = []
        for r, fmt in enumerate(fmts):
            for c, (src, metric) in enumerate(cols):
                ax = axes[r, c]
                ax.set_xscale("log")
                ax.set_yscale("log")
                series = grid.get((fmt, src, metric), {})
                for lane in self.order:
                    # A lane timed at one size only (headline-only) has no line to draw.
                    if len(series.get(lane, [])) < 2:
                        continue
                    pts = sorted(series[lane])
                    st = self.style[lane]
                    ls = (0, (4, 2)) if st["class"] == "partial" else "-"
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=st[theme], linestyle=ls,
                            linewidth=2.2 if st["class"] == "subject" else 1.2, zorder=3 if lane == SUBJECT else 2)
                    for p in pts:
                        self.mark(ax, lane, p[0], p[1], theme, size=4.5, label=False)
                    if lane not in shown:
                        shown.append(lane)
                ax.xaxis.set_major_locator(FixedLocator([v for v, _ in sizes]))
                ax.xaxis.set_major_formatter(FixedFormatter([n for _, n in sizes]))
                ax.xaxis.set_minor_locator(NullLocator())
                ax.set_xlim(min(0.65, min(xs) / 1.25), max(xs) * 1.4)
                ax.tick_params(axis="x", labelsize=7.5)
                count_axis(ax.yaxis)
                ax.margins(y=0.12)
                unit = "MiB/s" if metric == "speed" else "peak RSS (MiB)"
                if r == 0:
                    label = next(i["label"] for i in self.cases.values() if i["source"] == src)
                    head = ("Throughput" if metric == "speed" else "Peak memory") if both else textwrap.fill(label, 26)
                    ax.set_title(head, fontsize=9.5 if both else 9, color=t["ink"], pad=6, linespacing=1.2,
                                 fontweight="bold" if both else "normal")
                if r == len(fmts) - 1:
                    ax.set_xlabel("FASTQ size, log scale", fontsize=8.5, color=t["ink2"])
                if c == 0 or both:
                    ax.set_ylabel(f"{fmt} input\n{unit}, log scale" if c == 0 else f"{unit}, log scale", fontsize=9, color=t["ink2"])
        gz = " (gzip rows: the same reads compressed)" if "gzip" in fmts else ""
        what = (f"Each synthetic shape, 1 MiB to 1 GiB of FASTQ{gz}" if len(sources) > 1
                else "The {} shape ({}), 1 MiB to 1 GiB of FASTQ{}".format(sources[0], next(i['label'] for i in self.cases.values() if i['source'] == sources[0]).split(": ", 1)[-1], gz))
        if both:
            stem, title = f"{family}-scaling", "Throughput and peak memory by input size"
            note = (f"Throughput: higher is better, MiB of FASTQ {self.words['done']} per second of wall time; at 1 MiB, starting "
                    "the process dominates. Peak memory: lower is better; a flat line means memory does not grow with the input")
        elif metrics[0] == "speed":
            stem, title = f"{family}-scaling-speed", "Throughput by input size"
            note = (f"Higher is better. MiB of FASTQ {self.words['done']} per second of wall time; at 1 MiB, starting the process "
                    "dominates")
        else:
            stem, title = f"{family}-scaling-memory", "Peak memory by input size"
            note = ("Lower is better. A flat line means memory does not grow with the input; the kernel reports peak RSS in "
                    "256 KiB steps, so small wiggles are rounding")
        self.finish(fig, stem, theme, title, f"{note}. {what}. Each point is the median of its runs. Dashed: partial structure check.",
                    handles=self.mark_handles([n for n in self.order if n in shown], theme), ncol=len(shown))
        return True

    # ---- a case parameter --------------------------------------------------

    def param_figure(self, family, param, theme):
        """Columns: wall time and peak RSS; rows: formats; x: the case parameter; one line per lane."""
        grid = {}
        for (fam, case), lanes in self.rows.items():
            info = self.cases.get((fam, case))
            if fam != family or not info:
                continue
            values = dict(kv.split("=", 1) for kv in info.get("params", "-").split(";") if "=" in kv)
            if param not in values:
                continue
            x = float(values[param])
            for lane in lanes:
                t, m = self.value(family, case, lane, "speed"), self.value(family, case, lane, "memory")
                if t and m:
                    grid.setdefault(info["format"], {}).setdefault(lane, []).append((x, t / 1e3, m))
        fmts = [f for f in ("plain", "gzip") if f in grid]
        if not fmts:
            return False
        t = THEMES[theme]
        fig, axes = self.frame(theme, len(fmts), 2, width=9.0, height=3.4 * len(fmts) + 0.4, sharex=True)
        shown = []
        for r, fmt in enumerate(fmts):
            for c, (col, unit) in enumerate(((1, "wall time (s)"), (2, "peak RSS (MiB)"))):
                ax = axes[r, c]
                ax.set_xscale("log")
                ax.set_yscale("log")
                for lane in self.order:
                    pts = sorted(grid[fmt].get(lane, []))
                    if len(pts) < 2:
                        continue
                    st = self.style[lane]
                    ls = (0, (4, 2)) if st["class"] == "partial" else "-"
                    ax.plot([p[0] for p in pts], [p[col] for p in pts], color=st[theme], linestyle=ls,
                            linewidth=2.2 if st["class"] == "subject" else 1.2, zorder=3 if lane == SUBJECT else 2)
                    for p in pts:
                        self.mark(ax, lane, p[0], p[col], theme, size=4.5, label=False)
                    if lane not in shown:
                        shown.append(lane)
                count_axis(ax.xaxis)
                count_axis(ax.yaxis)
                ax.margins(y=0.12)
                if r == len(fmts) - 1:
                    ax.set_xlabel(f"{param.upper()}, log scale", fontsize=8.5, color=t["ink2"])
                ax.set_ylabel(f"{fmt} input\n{unit}, log scale" if c == 0 else f"{unit}, log scale", fontsize=9, color=t["ink2"])
        label = next(i["label"] for i in self.cases.values() if i["family"] == family)
        self.finish(fig, f"{family}-lines", theme, f"Time and peak memory against {param.upper()}",
                    f"Lower is better. {label}, one input; each point is the median of its runs. A line that rises with "
                    f"{param.upper()} holds what it keeps in memory. Dashed: partial structure check.",
                    handles=self.mark_handles([n for n in self.order if n in shown], theme), ncol=len(shown))
        return True

    # ---- edge cases ---------------------------------------------------------

    def edges(self, path, theme):
        rows = read(path)
        if not rows:
            return
        t = THEMES[theme]
        tools = list(dict.fromkeys(r["tool"] for r in rows))
        labels = {}
        for r in rows:
            labels[int(r["row"])] = (r["label"], r["expected"])
        nrows, ncols = len(labels), len(tools)
        fig, axes = self.frame(theme, width=2.1 + 3.0 + 1.15 * ncols, height=0.32 * nrows + 1.0, grid=None)
        ax = axes[0, 0]
        for side in ("left", "bottom"):
            ax.spines[side].set_visible(False)
        ax.set_xlim(-2.6, ncols)
        ax.set_ylim(nrows, -1.1)
        ax.set_xticks([])
        ax.set_yticks([])
        style = {v: (c, s[0 if theme == "light" else 1]) for v, c, _, s in VERDICTS}
        for r in rows:
            y, x = int(r["row"]), tools.index(r["tool"])
            color, strength = style[r["verdict"]]
            if color:
                ax.add_patch(Rectangle((x + 0.04, y + 0.06), 0.92, 0.88, facecolor=mix(color, t["bg"], strength), edgecolor="none"))
            # A long verdict wraps onto two lines in a smaller size, so it stays inside its cell.
            text, size = r["cell"], 8
            if len(text) > 16:
                text, size = textwrap.fill(text, 16), 6.8
            ax.text(x + 0.5, y + 0.5, text, ha="center", va="center", fontsize=size, color=t["ink"], linespacing=1.05,
                    fontweight="bold" if r["verdict"] in ("wrong", "crash") else "normal")
        prev = None
        for y, (label, expected) in sorted(labels.items()):
            group = "reject" if expected.startswith("reject") else "policy" if expected == "tool policy" else "value"
            if prev and group != prev:
                ax.axhline(y + 0.03, xmin=0, xmax=1, color=t["grid"], linewidth=0.8)
            prev = group
            ax.text(-2.55, y + 0.5, label, ha="left", va="center", fontsize=8.5, color=t["ink"])
            ax.text(-0.1, y + 0.5, expected, ha="right", va="center", fontsize=8, color=t["muted"])
        ax.text(-0.1, -0.55, "expected", ha="right", va="center", fontsize=8, color=t["muted"], style="italic")
        for x, tool in enumerate(tools):
            ax.text(x + 0.5, -0.55, tool, ha="center", va="center", fontsize=9, color=t["ink"], fontweight="bold" if tool == SUBJECT else "normal")
        handles = [Patch(facecolor=mix(c, t["bg"], s[0 if theme == "light" else 1]) if c else t["bg"], edgecolor=t["grid"] if not c else "none",
                         label=self.words["wrong"] if v == "wrong" else w) for v, c, w, s in VERDICTS]
        self.finish(fig, "edges", theme, "What each tool does on edge cases",
                    "Each cell is what the tool printed; color shows whether that was right.", handles=handles)

    def run(self):
        # The two themes draw in two processes; each figure is independent.
        with multiprocessing.get_context("fork").Pool(len(THEMES)) as pool:
            pool.map(self.draw, list(THEMES))

    def draw(self, theme):
        edges = self.report / "edges.tsv"
        if edges.exists():
            self.edges(edges, theme)
        grouped = {fam for _, _, cols in self.sections for fam, _ in cols}
        for slug, title, cols in self.sections:
            self.section_overview(slug, title, cols, theme)
            for metric in ("speed", "memory"):
                self.section_bars(slug, title, cols, metric, theme)
            self.section_scaling(slug, title, cols, theme)
        for family in self.mains:
            if family in grouped:
                continue
            if any(f == family for f, _ in self.rows):
                self.overview(family, theme)
                for metric in METRICS:
                    self.metric_bars(family, metric, theme)
                # One scaling source: throughput and memory side by side in one figure.
                if len([src for src in self.scaling if any(i["source"] == src and i["family"] == family for i in self.cases.values())]) == 1:
                    self.scaling_figure(family, ["speed", "memory"], theme)
                else:
                    for metric in ("speed", "memory"):
                        self.scaling_figure(family, [metric], theme)
        for family, param in self.lines:
            self.param_figure(family, param, theme)


if __name__ == "__main__":
    Figures(*sys.argv[1:8]).run()
