"""Wiki schematics, light and dark, in the bench figures' theme.

    tools/venv/bin/python tools/wiki_figures.py [OUT_DIR]

Writes each figure to OUT_DIR (default wiki/assets) as an SVG for editing and
the PNG the wiki shows (200 dpi on the theme's page color, 256 colors, no
dithering). Schematics: the parts of a FASTQ record, two-file and interleaved
pair layouts, the two sampling modes, adapter trimming by mate overlap and by
adapter search, and the 3' quality-tail window. Titles and legends belong to the
wiki page, so the figures carry none of their own. Bench figures are not
copied: the wiki links to them in the repository, so a bench rerun updates it.
"""

import io
import sys
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402
from PIL import Image  # noqa: E402

plt.rcParams.update({"font.family": "DejaVu Sans", "svg.hashsalt": "z-fastq"})

THEMES = {
    "light": dict(bg="#ffffff", ink="#1f2328", ink2="#59636e", muted="#818b98", grid="#d1d9e0",
                  insert="#8fb8de", adapter="#f7a41d", region="#fff4e0"),
    "dark": dict(bg="#0d1117", ink="#f0f6fc", ink2="#9198a1", muted="#656c76", grid="#3d444d",
                 insert="#3b6ea5", adapter="#f7a41d", region="#2a1f0c"),
}


MONO = "DejaVu Sans Mono"


def save(fig, out, stem, theme):
    bg = THEMES[theme]["bg"]
    fig.savefig(out / f"{stem}-{theme}.svg", format="svg", facecolor=bg, bbox_inches="tight", pad_inches=0.15)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, facecolor=bg, bbox_inches="tight", pad_inches=0.15)
    image = Image.open(buf).convert("RGB").quantize(256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    image.save(out / f"{stem}-{theme}.png", format="png", optimize=True)
    plt.close(fig)


def bar(ax, x0, x1, y, color, t, height=0.42):
    ax.add_patch(Rectangle((x0, y - height / 2), x1 - x0, height, facecolor=color, edgecolor=t["bg"], linewidth=1.2))


def arrow(ax, x0, x1, y, t):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=11, color=t["ink2"], linewidth=1.1))


def cut(ax, x, y0, y1, t):
    ax.plot([x, x], [y0, y1], color=t["ink"], linewidth=1.3, linestyle=(0, (3, 2)))


def label(ax, x, y, text, t, color=None, **kw):
    ax.text(x, y, text, color=color or t["ink"], fontsize=9, **kw)


def adapter_trimming(out, theme):
    t = THEMES[theme]
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 3.0), gridspec_kw=dict(hspace=0.12, height_ratios=(2.6, 2.2)))
    for ax in axes:
        ax.set_facecolor(t["bg"])
        ax.axis("off")

    ax = axes[0]
    ax.set_xlim(-52, 118)
    ax.set_ylim(-0.6, 2.0)
    label(ax, -52, 1.75, "Insert shorter than the reads: both mates are cut to the overlap", t, fontweight="bold")
    ax.add_patch(Rectangle((0, -0.45), 60, 2.0, facecolor=t["region"], edgecolor="none"))
    bar(ax, 0, 60, 1.0, t["insert"], t)
    bar(ax, 60, 100, 1.0, t["adapter"], t)
    arrow(ax, 0, 104, 1.0, t)
    label(ax, -50, 0.92, "mate 1", t, color=t["ink2"])
    bar(ax, -40, 0, 0.0, t["adapter"], t)
    bar(ax, 0, 60, 0.0, t["insert"], t)
    arrow(ax, 64, -44, 0.0, t)
    label(ax, 66, -0.08, "mate 2, reverse complement", t, color=t["ink2"])
    cut(ax, 0, -0.45, 1.55, t)
    cut(ax, 60, -0.45, 1.55, t)
    label(ax, 30, 1.32, "overlap", t, color=t["ink2"], ha="center")

    ax = axes[1]
    ax.set_xlim(-52, 118)
    ax.set_ylim(-0.6, 1.6)
    label(ax, -52, 1.35, "No overlap found: each mate is searched for its own adapter", t, fontweight="bold")
    bar(ax, 0, 74, 0.55, t["insert"], t)
    bar(ax, 74, 100, 0.55, t["adapter"], t)
    arrow(ax, 0, 104, 0.55, t)
    label(ax, -50, 0.47, "mate 1", t, color=t["ink2"])
    cut(ax, 74, 0.05, 1.05, t)
    label(ax, 76, 0.95, "cut where the adapter starts", t, color=t["ink2"])
    ax.add_patch(FancyArrowPatch((-4, -0.2), (40, -0.2), arrowstyle="-|>", mutation_scale=10,
                                 color=t["muted"], linewidth=1.0, linestyle=(0, (2, 2))))
    label(ax, 42, -0.28, "the search starts up to 4 bases before the read", t, color=t["muted"])

    handles = [Rectangle((0, 0), 1, 1, facecolor=t["insert"]), Rectangle((0, 0), 1, 1, facecolor=t["adapter"])]
    fig.legend(handles, ["insert", "adapter"], loc="lower center", ncol=2, frameon=False, fontsize=9,
               labelcolor=t["ink"], bbox_to_anchor=(0.5, -0.04))
    save(fig, out, "adapter-trimming", theme)


def quality_tail(out, theme):
    t = THEMES[theme]
    qualities = [36, 37, 36, 35, 37, 34, 36, 33, 35, 34, 32, 33, 30, 31, 28, 30, 27, 22, 25, 18, 14, 12, 9, 11, 6, 4, 5, 2]
    window, threshold = 4, 20
    fig, ax = plt.subplots(figsize=(7.4, 2.8))
    ax.set_facecolor(t["bg"])
    starts = range(len(qualities) - window, 0, -1)
    passing = next(s for s in starts if sum(qualities[s:s + window]) / window >= threshold)
    ax.axvspan(passing - 0.5, passing + window - 0.5, color=t["region"], zorder=0)
    for i, q in enumerate(qualities):
        kept = i <= passing
        ax.bar(i, q, width=0.8, color=t["insert"] if kept else t["muted"], zorder=2)
    ax.axhline(threshold, color=t["adapter"], linewidth=1.3, zorder=3)
    ax.text(len(qualities) - 0.4, threshold + 1, "Q = 20", color=t["adapter"], fontsize=9, ha="right", va="bottom")
    ax.axvline(passing + 0.5, color=t["ink"], linewidth=1.3, linestyle=(0, (3, 2)), zorder=3)
    ax.annotate("first window from the 3' end\nwith mean quality at least Q", xy=(passing + 1.5, 31), xytext=(passing + 4.5, 33),
                color=t["ink2"], fontsize=9, va="center", arrowprops=dict(arrowstyle="-", color=t["ink2"], linewidth=0.8))
    ax.annotate("cut after its first base", xy=(passing + 0.5, 6), xytext=(passing - 10, 6), color=t["ink"], fontsize=9,
                va="center", arrowprops=dict(arrowstyle="-|>", color=t["ink"], linewidth=0.8))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["grid"])
    ax.tick_params(colors=t["ink2"], labelsize=8)
    ax.set_xticks([0, len(qualities) - 1], ["5'", "3'"])
    ax.set_ylabel("Phred quality", color=t["ink2"], fontsize=9)
    ax.set_ylim(0, 40)
    save(fig, out, "quality-tail", theme)


def fastq_record(out, theme):
    t = THEMES[theme]
    fig, ax = plt.subplots(figsize=(7.4, 2.2))
    ax.set_facecolor(t["bg"])
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.6, 4.2)
    lines = [("@read42 1:N:0:ATCACG", "header: @, then the identifier"), ("ACGTNACGTTAGC", "sequence"),
             ("+", "plus line"), ("IIIIH#IIIGG?I", "quality: one byte per base")]
    for row, (text, note) in enumerate(lines):
        y = 3.4 - row
        ax.text(2, y, text, family=MONO, fontsize=11, color=t["ink"], va="center")
        ax.text(52, y, note, fontsize=9, color=t["ink2"], va="center")
    ax.plot([2, 2], [2.4, 0.4], alpha=0)
    ax.annotate("", xy=(1.2, 2.55), xytext=(1.2, 0.25),
                arrowprops=dict(arrowstyle="-", color=t["adapter"], linewidth=1.4))
    ax.text(52, -0.35, "sequence and quality have the same length", fontsize=9, color=t["adapter"], va="center")
    save(fig, out, "fastq-record", theme)


def pair_layouts(out, theme):
    t = THEMES[theme]
    fig, ax = plt.subplots(figsize=(7.4, 2.6))
    ax.set_facecolor(t["bg"])
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.4, 6.6)
    colors = [t["insert"], t["adapter"]]
    for column, (x, name) in enumerate([(4, "R1.fastq"), (22, "R2.fastq")]):
        ax.text(x + 6, 6.1, name, fontsize=9, color=t["ink"], ha="center", family=MONO)
        for i in range(3):
            ax.add_patch(Rectangle((x, 4.6 - i * 1.4), 12, 1.0, facecolor=colors[column], edgecolor=t["bg"]))
            ax.text(x + 6, 5.1 - i * 1.4, f"pair {i + 1}, R{column + 1}", fontsize=8, color=t["bg"], ha="center", va="center")
    ax.text(66, 6.1, "interleaved.fastq", fontsize=9, color=t["ink"], ha="center", family=MONO)
    for i in range(6):
        ax.add_patch(Rectangle((60, 5.2 - i * 0.95), 12, 0.8, facecolor=colors[i % 2], edgecolor=t["bg"]))
        ax.text(66, 5.6 - i * 0.95, f"pair {i // 2 + 1}, R{i % 2 + 1}", fontsize=8, color=t["bg"], ha="center", va="center")
    ax.add_patch(FancyArrowPatch((37, 4.0), (57, 4.0), arrowstyle="-|>", mutation_scale=12, color=t["ink2"], linewidth=1.2))
    ax.text(47, 4.35, "interleave", fontsize=9, color=t["ink2"], ha="center")
    ax.add_patch(FancyArrowPatch((57, 2.2), (37, 2.2), arrowstyle="-|>", mutation_scale=12, color=t["ink2"], linewidth=1.2))
    ax.text(47, 2.55, "deinterleave", fontsize=9, color=t["ink2"], ha="center")
    ax.text(78, 3.0, "the names of each\npair are checked\nbefore it is written", fontsize=8.5, color=t["muted"], va="center")
    save(fig, out, "pair-layouts", theme)


def sampling(out, theme):
    t = THEMES[theme]
    fig, ax = plt.subplots(figsize=(7.4, 2.2))
    ax.set_facecolor(t["bg"])
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0.7, 4.4)
    fraction = {1, 4, 5, 9, 13}
    count = {2, 4, 8, 11}
    rows = [(3.0, "--fraction 0.3", fraction, "each record is kept or not on its own:\nabout 30%, here 5 of 15"),
            (0.9, "--count 4", count, "exactly 4 records, chosen after\na first pass counts all 15")]
    for y, title, kept, note in rows:
        ax.text(0, y + 0.95, title, fontsize=9, color=t["ink"], family=MONO)
        for i in range(15):
            color = t["adapter"] if i in kept else t["grid"]
            ax.add_patch(Rectangle((i * 3.4, y), 2.9, 0.7, facecolor=color, edgecolor="none"))
        ax.text(54, y + 0.35, note, fontsize=8.5, color=t["ink2"], va="center")
    save(fig, out, "sampling", theme)


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "wiki/assets")
    out.mkdir(parents=True, exist_ok=True)
    for theme in THEMES:
        fastq_record(out, theme)
        pair_layouts(out, theme)
        sampling(out, theme)
        adapter_trimming(out, theme)
        quality_tail(out, theme)


if __name__ == "__main__":
    main()
