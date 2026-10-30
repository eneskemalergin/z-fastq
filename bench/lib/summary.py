"""Summary figures across the benchmark features, light and dark, from the tables
summary.sh writes; nothing here reads results or times anything.
    summary.py TABLES_DIR FIGURE_DIR
TABLES_DIR holds jobs.tsv, broken.tsv, threads.tsv, and style.tsv (one line per
tool label). Figures, one question each:
  01-jobs     time and peak memory per job: the fastest complete tool and the
              fastest partial one, as multiples of z-fastq's
  02-broken   of the malformed inputs each tool was given, how many it rejected
  03-threads  the jobs with thread counts: each tool's wall time by threads,
              against z-fastq on one core
  04-threads-memory  the same runs on time against memory, as multiples of
              z-fastq's: each tool's path through its thread counts, with
              equal-cost curves
The figures use plot.py's theme, frame, and file output.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import plot  # noqa: E402
from plot import MARKERS, SUBJECT, THEMES, VERDICTS, count_axis, mix, read, times  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator  # noqa: E402


class Summary:
    def __init__(self, tables, out):
        self.tables = Path(tables)
        self.out = Path(out)
        self.out.mkdir(parents=True, exist_ok=True)
        self.style = {r["tool"]: r for r in read(self.tables / "style.tsv")}
        self.jobs = read(self.tables / "jobs.tsv")
        self.broken = read(self.tables / "broken.tsv")
        self.threads = read(self.tables / "threads.tsv")
        # plot.py's frame and finish need only the output directory.
        self.fig = plot.Figures.__new__(plot.Figures)
        self.fig.out = self.out

    def accent(self, theme):
        return self.style[SUBJECT][theme]

    # ---- 01: each job against the fastest other tools -------------------------

    def jobs_figure(self, theme):
        t = THEMES[theme]
        rows = self.jobs
        n = len(rows)
        fig, axes = self.fig.frame(theme, 1, 2, width=11.5, height=0.42 * n + 1.6, grid="x")
        metrics = [("Wall time", "c_time", "p_time"), ("Peak memory", "c_mem", "p_mem")]
        gray = t["ink2"]
        for ax, (name, ck, pk) in zip(axes[0], metrics):
            ax.set_xscale("log")
            vals = [1.0]
            for y, r in enumerate(rows):
                for key, label_key, face, dy in ((ck, "complete", gray, -0.17), (pk, "partial", "none", 0.17)):
                    if r[key] == "-":
                        continue
                    v = float(r[key])
                    vals.append(v)
                    ax.plot([1, v], [y + dy, y + dy], color=mix(gray, t["bg"], 0.35), linewidth=1.0, zorder=1)
                    ax.plot(v, y + dy, "o", markersize=6.5, markerfacecolor=face, markeredgecolor=gray, markeredgewidth=1.3, zorder=3)
                    ax.annotate(f"{r[label_key]} {times(v)}", (max(v, 1), y + dy), xytext=(6, 0), textcoords="offset points",
                                ha="left", va="center", fontsize=7.6, color=t["ink"])
            ax.axvline(1, color=self.accent(theme), linewidth=1.6, linestyle=(0, (4, 3)), zorder=2)
            lo, hi = min(vals), max(vals)
            ax.set_xlim(min(0.6, lo / 1.4), hi * 6)
            ax.set_ylim(n - 0.5, -0.7)
            ax.set_yticks(range(n))
            ax.set_yticklabels([r["job"] for r in rows] if ax is axes[0, 0] else [], fontsize=9, color=t["ink"])
            ax.tick_params(axis="y", length=0)
            for y in range(1, n):
                ax.axhline(y - 0.5, color=t["grid"], linewidth=0.5, zorder=0)
            ax.text(1, -0.62, " z-fastq", color=self.accent(theme), fontsize=8, fontweight="bold", va="center", ha="left")
            ax.set_title(name, loc="left", fontsize=10, color=t["ink"], fontweight="bold", pad=8)
            ax.set_xlabel("the other tool's, as a multiple of z-fastq's, log scale", fontsize=8.5, color=t["ink2"])
            count_axis(ax.xaxis, "x")
        handles = [Line2D([], [], color=self.accent(theme), linewidth=1.6, linestyle=(0, (4, 3)), label="z-fastq (1x)"),
                   Line2D([], [], marker="o", linestyle="none", markersize=6.5, markerfacecolor=gray, markeredgecolor=gray,
                          label="fastest tool that does the whole job"),
                   Line2D([], [], marker="o", linestyle="none", markersize=6.5, markerfacecolor="none", markeredgecolor=gray,
                          markeredgewidth=1.3, label="fastest tool that does part of it")]
        self.fig.finish(fig, "01-jobs", theme, "Each job against the fastest other tools",
                        "Right of the dashed line, the other tool takes longer or holds more memory than z-fastq. Every tool on one "
                        "CPU, 1 GiB plain inputs, the geometric mean of the ratios within each batch; the memory panel shows the same "
                        "two tools as the time panel.", handles=handles, ncol=3)

    # ---- 02: broken input -----------------------------------------------------

    def broken_figure(self, theme):
        t = THEMES[theme]
        rows = self.broken
        n = len(rows)
        fig, axes = self.fig.frame(theme, 1, 1, width=9.5, height=0.36 * n + 1.5, grid="x")
        ax = axes[0, 0]
        colors = {v: mix(c, t["bg"], s[0 if theme == "light" else 1] + 0.25) for v, c, _, s in VERDICTS if c}
        parts = [("rejected", colors["right"], "rejected"), ("passed", colors["wrong"], "passed it on, exit 0"),
                 ("crashed", colors["crash"], "crashed")]
        for y, r in enumerate(rows):
            given, left = int(r["given"]), 0.0
            for key, color, _ in parts:
                share = int(r[key]) / given
                if share:
                    ax.barh(y, share, left=left, height=0.62, color=color, edgecolor="none", zorder=2)
                left += share
            ax.text(1.02, y, f"{r['rejected']} of {given}", va="center", ha="left", fontsize=8.2, color=t["ink"],
                    fontweight="bold" if r["tool"] == SUBJECT else "normal")
        ax.set_xlim(0, 1.16)
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_yticks(range(n))
        ax.set_yticklabels([r["tool"] for r in rows], fontsize=9, color=t["ink"])
        for lbl in ax.get_yticklabels():
            if lbl.get_text() == SUBJECT:
                lbl.set_fontweight("bold")
        ax.tick_params(axis="y", length=0)
        ax.xaxis.set_major_locator(FixedLocator([0, 0.25, 0.5, 0.75, 1]))
        ax.xaxis.set_major_formatter(FixedFormatter(["0", "25%", "50%", "75%", "100%"]))
        ax.set_xlabel("share of the malformed inputs the tool was given", fontsize=8.5, color=t["ink2"])
        self.fig.finish(fig, "02-broken", theme, "What each tool does with a broken record or pair",
                        "Every edge case across the feature benchmarks that a tool must reject: a truncated or malformed record, a "
                        "mismatched pair, a damaged gzip member. A tool counts only on the jobs it does; rows with no single right "
                        "answer are left out.", handles=[Patch(color=c, label=w) for _, c, w in parts], ncol=3)

    # ---- 03: threads ----------------------------------------------------------

    def threads_figure(self, theme):
        t = THEMES[theme]
        jobs = list(dict.fromkeys(r["job"] for r in self.threads))
        if not jobs:
            return
        fig, axes = self.fig.frame(theme, 1, len(jobs), width=5.2 * len(jobs) + 0.6, height=4.4)
        tools_shown, counts = [], set()
        for ax, job in zip(axes[0], jobs):
            rows = [r for r in self.threads if r["job"] == job]
            ax.set_xscale("log", base=2)
            ax.set_yscale("log")
            for tool in dict.fromkeys(r["tool"] for r in rows):
                pts = sorted((int(r["threads"]), float(r["wall_s"])) for r in rows if r["tool"] == tool)
                if tool == SUBJECT:
                    ax.axhline(pts[0][1], color=self.accent(theme), linewidth=1.6, linestyle=(0, (4, 3)), zorder=2)
                    continue
                st = self.style[tool]
                gray = t["ink2"]
                ax.plot([p[0] for p in pts], [p[1] for p in pts], color=gray, linewidth=1.3, zorder=3,
                        marker=MARKERS[st["shape"]], markersize=5.5, markerfacecolor="none", markeredgecolor=gray, markeredgewidth=1.3)
                ax.annotate(tool, pts[-1], xytext=(8, 0), textcoords="offset points", va="center", ha="left", fontsize=8, color=t["ink"])
                counts.update(p[0] for p in pts)
                if tool not in tools_shown:
                    tools_shown.append(tool)
            ax.set_title(job, loc="left", fontsize=10, color=t["ink"], fontweight="bold", pad=8)
            ax.set_xlabel("threads, log 2 scale", fontsize=8.5, color=t["ink2"])
            ax.set_ylabel("wall time (s), log scale", fontsize=8.5, color=t["ink2"])
            ax.margins(y=0.15)
            ax.set_xlim(0.8, 40)
            count_axis(ax.yaxis)
        ticks = sorted(counts)
        for ax in axes.flat:
            ax.xaxis.set_major_locator(FixedLocator(ticks))
            ax.xaxis.set_major_formatter(FixedFormatter([str(x) for x in ticks]))
            ax.xaxis.set_minor_locator(NullLocator())
        handles = [Line2D([], [], color=self.accent(theme), linewidth=1.6, linestyle=(0, (4, 3)), label="z-fastq, 1 thread")]
        for tool in tools_shown:
            st = self.style[tool]
            handles.append(Line2D([], [], color=t["ink2"], linewidth=1.3, marker=MARKERS[st["shape"]], markersize=5.5,
                                  markerfacecolor="none", markeredgecolor=t["ink2"], label=tool))
        self.fig.finish(fig, "03-threads", theme, "One core against many",
                        "Lower is better. The other tools at 1 to 16 threads on as many physical cores, z-fastq on one (dashed), in "
                        "the same batches; 1 GiB plain inputs, geometric means. Where a line drops below the dashed one, that tool "
                        "matches z-fastq's wall time, spending its threads to do it.", handles=handles, ncol=len(handles))

    # ---- 04: threads on time against memory -----------------------------------

    def threads_memory_figure(self, theme):
        t = THEMES[theme]
        jobs = list(dict.fromkeys(r["job"] for r in self.threads))
        if not jobs:
            return
        fig, axes = self.fig.frame(theme, 1, len(jobs), width=5.6 * len(jobs) + 0.6, height=5.2)
        shown = []
        for ax, job in zip(axes[0], jobs):
            rows = [r for r in self.threads if r["job"] == job]
            z = next(r for r in rows if r["tool"] == SUBJECT)
            zt, zm = float(z["wall_s"]), float(z["mib"])
            ax.set_xscale("log")
            ax.set_yscale("log")
            pts_all = [(1.0, 1.0)]
            for tool in dict.fromkeys(r["tool"] for r in rows):
                if tool == SUBJECT:
                    continue
                st = self.style[tool]
                pts = sorted((int(r["threads"]), float(r["wall_s"]) / zt, float(r["mib"]) / zm) for r in rows if r["tool"] == tool)
                pts_all += [(p[1], p[2]) for p in pts]
                gray = t["ink2"]
                ax.plot([p[1] for p in pts], [p[2] for p in pts], color=gray, linewidth=1.3, zorder=3,
                        marker=MARKERS[st["shape"]], markersize=5.5, markerfacecolor="none", markeredgecolor=gray, markeredgewidth=1.3)
                # Thread counts whose points nearly coincide share one label ("8, 16").
                groups = []
                for n, x, y in pts:
                    if groups and abs(math.log(x / groups[-1][1])) < 0.04 and abs(math.log(y / groups[-1][2])) < 0.04:
                        groups[-1][0].append(n)
                    else:
                        groups.append(([n], x, y))
                for ns, x, y in groups:
                    ax.annotate(", ".join(map(str, ns)), (x, y), xytext=(0, 7), textcoords="offset points", ha="center", va="bottom",
                                fontsize=7, color=t["ink2"])
                n, x, y = pts[0]
                ax.annotate(tool, (x, y), xytext=(8, -2), textcoords="offset points", ha="left", va="top", fontsize=8.5, color=t["ink"])
                if tool not in shown:
                    shown.append(tool)
            xs, ys = [p[0] for p in pts_all], [p[1] for p in pts_all]
            lo_x, hi_x = min(min(xs), 1) / 1.6, max(xs) * 1.6
            lo_y, hi_y = min(min(ys), 1) / 1.6, max(ys) * 2.2
            ax.set_xlim(lo_x, hi_x)
            ax.set_ylim(lo_y, hi_y)
            # The corner where a tool would be faster and smaller than z-fastq, then the 1x lines.
            ax.fill_between([lo_x, 1], lo_y, 1, color=t["region"], zorder=0, linewidth=0)
            for line in (ax.axvline, ax.axhline):
                line(1, color=self.accent(theme), linewidth=0.9, linestyle=(0, (4, 3)), alpha=0.7, zorder=1)
            # Equal-cost curves C = (time ratio + memory ratio) / 2, as in the feature overviews.
            c_lo, c_hi = (lo_x + lo_y) / 2, (hi_x + hi_y) / 2
            steps = [m * 10 ** e for e in range(-2, 5) for m in (1, 1.5, 2, 3, 5, 7)]
            levels = [c for c in steps if c_lo < c < c_hi]
            if len(levels) > 6:
                levels = sorted({levels[round(i * (len(levels) - 1) / 5)] for i in range(6)})
            grid = np.exp(np.linspace(math.log(lo_x), math.log(hi_x), 400))
            for c in levels:
                gy = 2 * c - grid
                keep = (gy > lo_y) & (gy < hi_y)
                if keep.sum() < 2:
                    continue
                ax.plot(grid[keep], gy[keep], color=t["muted"], linewidth=0.7, linestyle=(0, (4, 3)), alpha=0.8, zorder=1)
                room = keep & (grid > lo_x * 1.06) & (grid < hi_x / 1.45) & (gy > lo_y * 1.1) & (gy < hi_y / 1.1)
                if room.any():
                    i = np.flatnonzero(room)[0]
                    ax.text(grid[i], gy[i], f"C={c:g}", fontsize=7, color=t["muted"], ha="left", va="center",
                            bbox=dict(boxstyle="round,pad=0.12", facecolor=t["bg"], edgecolor="none"), zorder=1.4)
            ax.plot(1, 1, "o", markersize=11, markerfacecolor=self.accent(theme), markeredgecolor=t["ink"], markeredgewidth=1.0, zorder=5)
            ax.annotate("z-fastq, 1 thread", (1, 1), xytext=(9, -4), textcoords="offset points", ha="left", va="top", fontsize=8.5,
                        color=t["ink"], fontweight="bold")
            ax.set_title(job, loc="left", fontsize=10, color=t["ink"], fontweight="bold", pad=8)
            ax.set_xlabel("wall time relative to z-fastq, log scale", fontsize=8.5, color=t["ink2"])
            ax.set_ylabel("peak RSS relative to z-fastq, log scale", fontsize=8.5, color=t["ink2"])
            count_axis(ax.xaxis, "x", dense=hi_x / lo_x < 12)
            count_axis(ax.yaxis, "x", dense=hi_y / lo_y < 12)
        handles = [Line2D([], [], marker="o", linestyle="none", markersize=9, markerfacecolor=self.accent(theme),
                          markeredgecolor=t["ink"], label="z-fastq, 1 thread")]
        for tool in shown:
            st = self.style[tool]
            handles.append(Line2D([], [], color=t["ink2"], linewidth=1.3, marker=MARKERS[st["shape"]], markersize=5.5,
                                  markerfacecolor="none", markeredgecolor=t["ink2"], label=f"{tool}, by threads"))
        self.fig.finish(fig, "04-threads-memory", theme, "What the threads cost in memory",
                        "The runs of the threads figure on time against memory, each as a multiple of z-fastq's on one core; numbers "
                        "on a path are thread counts. Moving left, a tool gets faster with more threads; it stays high up because its "
                        "memory does not shrink with them. Dashed curves: equal cost C = (time + memory) / 2; the tinted corner is "
                        "faster and smaller than z-fastq.", handles=handles, ncol=len(handles))

    def run(self):
        for theme in THEMES:
            self.jobs_figure(theme)
            self.broken_figure(theme)
            self.threads_figure(theme)
            self.threads_memory_figure(theme)


if __name__ == "__main__":
    Summary(*sys.argv[1:3]).run()
