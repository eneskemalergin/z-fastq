#!/usr/bin/env python3
"""Generate bench/count/REPORT.md and figures from zebrac JSON."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

MPL_CONFIG_DIR = Path(__file__).resolve().parents[2] / "tools" / ".local" / "matplotlib"
MPL_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CONFIG_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import tabulate  # noqa: F401  -- pd.to_markdown()

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

DATASET_ORDER = ["Dense", "Variable", "Long"]
BASELINE = "z-fastq"
REFERENCE_TOOLS = frozenset({"seqfu"})
NO_RSS_RATIO = frozenset({"seqfu"})

# z-fastq, native inflate, rust adapters, then C/Nim CLIs.
PLAIN_TOOLS = [
    "z-fastq",
    "needletail",
    "helicase",
    "seqtk",
    "seqfu",
    "fqtools",
]
GZIP_TOOLS = [
    "z-fastq",
    "z-fastq-native",
    "needletail",
    "helicase",
    "seqtk",
    "seqfu",
    "fqtools",
]
SCALING_TOOLS = ["z-fastq", "needletail", "helicase", "seqtk", "fqtools"]
SCALING_GZIP_TOOLS = ["z-fastq", "z-fastq-native", "needletail", "helicase", "seqtk", "fqtools"]

COLORS = {
    "z-fastq": "#F7A41D",
    "z-fastq-native": "#FFB74D",
    "needletail": "#C45C26",
    "helicase": "#8B3A2A",
    "seqtk": "#6A1B9A",
    "seqfu": "#009485",
    "fqtools": "#555555",
}

DISPLAY = {
    "z-fastq": "z-fastq (ISA-L)",
    "z-fastq-native": "z-fastq (native)",
    "needletail": "Needletail",
    "helicase": "Helicase",
    "seqtk": "seqtk size",
    "seqfu": "SeqFu (descriptive)",
    "fqtools": "fqtools count",
}

FACET_WSPACE = 0.28
SIZE_MB_ORDER = [1, 5, 10, 25, 50, 100, 250, 500]
READ_COUNT_ORDER = [100000, 250000, 500000, 1000000]
MARKDOWNLINT_DISABLE = "<!-- markdownlint-disable MD024 MD032 MD033 MD036 MD041 MD049 -->"


class ReportCounters:
    def __init__(self) -> None:
        self.table = 1
        self.figure = 1

    def next_table(self) -> int:
        n = self.table
        self.table += 1
        return n

    def next_figure(self) -> int:
        n = self.figure
        self.figure += 1
        return n


def join(parts: list[str]) -> str:
    return "\n".join(parts)


def align_pipe_table(text: str) -> str:
    lines = [line.rstrip() for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        return text
    rows: list[list[str]] = []
    for line in lines:
        body = line.strip()
        if body.startswith("|"):
            body = body[1:]
        if body.endswith("|"):
            body = body[:-1]
        rows.append([cell.strip() for cell in body.split("|")])
    width = len(rows[0])
    if width == 0 or any(len(row) != width for row in rows):
        return text
    col_w = [max(len(row[i]) for row in rows) for i in range(width)]
    out: list[str] = []
    for ri, row in enumerate(rows):
        cells: list[str] = []
        for i, cell in enumerate(row):
            if ri == 1 and cell and set(cell) <= set("-:"):
                left = cell.startswith(":")
                right = cell.endswith(":")
                if left and right:
                    piece = ":" + ("-" * max(1, col_w[i] - 2)) + ":"
                elif left:
                    piece = ":" + ("-" * max(1, col_w[i] - 1))
                elif right:
                    piece = ("-" * max(1, col_w[i] - 1)) + ":"
                else:
                    piece = "-" * max(3, col_w[i])
                cells.append(piece)
            else:
                cells.append(cell.ljust(col_w[i]))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def to_markdown_aligned(df: pd.DataFrame, *, index: bool = True) -> str:
    return align_pipe_table(df.to_markdown(index=index))


def display_tool(tool: str) -> str:
    return DISPLAY.get(tool, tool)


def format_ratio(ratio: float | None) -> str:
    if ratio is None:
        return "n/a"
    if abs(ratio - 1.0) < 1e-9:
        return "1x"
    if ratio >= 100:
        return f"{ratio:.0f}x"
    if ratio >= 10:
        return f"{ratio:.1f}x"
    if ratio >= 1:
        return f"{ratio:.2f}x"
    return f"{ratio:.3f}x"


def _cell_float(row, key: str, default: float = 0.0) -> float:
    try:
        val = row[key]
    except (KeyError, IndexError):
        return default
    if pd.isna(val):
        return default
    return float(val)


def fmt_wall(row) -> str:
    mean = _cell_float(row, "mean")
    std = _cell_float(row, "stddev")
    return f"{mean:.3f}±{std:.3f} s"


def fmt_rss(row) -> str:
    mean = _cell_float(row, "peak_rss_mb")
    std = _cell_float(row, "peak_rss_stddev_mb")
    return f"{mean:.2f}±{std:.2f} MB"


def fmt_faults(row) -> str:
    mean = _cell_float(row, "minor_faults")
    std = _cell_float(row, "minor_faults_stddev")
    return f"{mean:.0f}±{std:.0f}"


def fmt_decoded_mibs(row) -> str:
    try:
        val = row["throughput_mibs"]
    except (KeyError, IndexError):
        return "n/a"
    if val is None or pd.isna(val):
        return "n/a"
    return f"{float(val):.1f} MiB/s"


def fmt_compressed_mibs(row) -> str:
    try:
        val = row["compressed_throughput_mibs"]
    except (KeyError, IndexError):
        val = None
    if val is not None and not pd.isna(val):
        return f"{float(val):.1f} MiB/s"
    try:
        ib = row["input_bytes"]
        mean = row["mean"]
    except (KeyError, IndexError):
        return "n/a"
    if ib is None or pd.isna(ib) or mean is None or pd.isna(mean) or float(mean) <= 0:
        return "n/a"
    return f"{(float(ib) / (1024.0 * 1024.0)) / float(mean):.1f} MiB/s"


def short_peer_version(name: str, raw: str) -> str:
    text = (raw or "").strip()
    if not text:
        return ""
    if name == "seqfu" and not text.lower().startswith("seqfu"):
        return f"seqfu {text.split()[0]}"
    if name in {"needletail", "helicase"}:
        for part in text.split(";"):
            part = part.strip()
            if part.startswith("engine="):
                return part.split("=", 1)[1].strip()
        return text.split(";")[0].strip()
    return text


def dataset_sort_key(name):
    try:
        return DATASET_ORDER.index(name), str(name)
    except ValueError:
        return len(DATASET_ORDER), str(name)


def tools_in_run(df: pd.DataFrame | None, order: list[str]) -> list[str]:
    if df is None or df.empty:
        return []
    present = set(df["tool"].unique())
    return [t for t in order if t in present]


def peer_tools(tools: list[str], baseline: str = BASELINE) -> list[str]:
    return [t for t in tools if t != baseline]


def _save(fig, path: Path) -> Path:
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def load_latest_manifest(results_dir: Path) -> dict:
    latest = results_dir / "LATEST"
    if not latest.is_file():
        raise SystemExit(f"error: missing {latest}")
    ts = latest.read_text(encoding="utf-8").strip()
    path = results_dir / f"run_{ts}.json"
    if not path.is_file():
        raise SystemExit(f"error: missing {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_metadata(results_dir: Path, manifest: dict) -> pd.DataFrame | None:
    name = manifest.get("metadata")
    if not name:
        return None
    path = results_dir / name
    if not path.is_file():
        return None
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            rows.append(json.loads(line))
    return pd.DataFrame(rows) if rows else None


def load_zebrac_json(path: Path, metadata_df: pd.DataFrame | None) -> pd.DataFrame:
    data = json.loads(path.read_text())
    meta_rows = []
    if metadata_df is not None and not metadata_df.empty:
        meta_rows = metadata_df[metadata_df["raw_json"] == str(path)].to_dict("records")
        if not meta_rows:
            meta_rows = metadata_df[
                metadata_df["raw_json"].astype(str).str.endswith(path.name)
            ].to_dict("records")
    meta_by_command = {row.get("command"): row for row in meta_rows}
    rows = []
    for idx, result in enumerate(data.get("results", [])):
        command = result.get("command", "")
        meta = meta_by_command.get(command, {})
        if not meta and idx < len(meta_rows):
            meta = meta_rows[idx]
        wall = result.get("wall_time", {})
        peak = result.get("peak_rss", {})
        minor = result.get("minor_faults", {})
        mean_s = float(wall.get("mean", 0) or 0) / 1_000_000_000.0
        std_s = float(wall.get("std_dev", 0) or 0) / 1_000_000_000.0
        rss = float(peak.get("mean", 0) or 0) / (1024.0 * 1024.0)
        rss_std = float(peak.get("std_dev", 0) or 0) / (1024.0 * 1024.0)
        faults = float(minor.get("mean", 0) or 0)
        faults_std = float(minor.get("std_dev", 0) or 0)
        workload = meta.get("workload", path.stem.split("__", 1)[0])
        section = meta.get("section", "")
        tool = meta.get("tool")
        if not tool:
            stem = path.stem
            tool = stem.split("__", 1)[1] if "__" in stem else "unknown"
        input_bytes = meta.get("input_bytes")
        decoded_bytes = meta.get("decoded_bytes")
        decoded_mib = None
        throughput = None
        compressed_throughput = None
        if decoded_bytes is not None:
            decoded_mib = float(decoded_bytes) / (1024.0 * 1024.0)
            if mean_s > 0:
                throughput = decoded_mib / mean_s
        if input_bytes is not None and mean_s > 0:
            compressed_throughput = (float(input_bytes) / (1024.0 * 1024.0)) / mean_s
        dataset = workload if workload in DATASET_ORDER else None
        rows.append(
            {
                "tool": tool,
                "section": section,
                "workload": workload,
                "dataset": dataset,
                "mean": mean_s,
                "stddev": std_s,
                "peak_rss_mb": rss,
                "peak_rss_stddev_mb": rss_std,
                "minor_faults": faults,
                "minor_faults_stddev": faults_std,
                "input_bytes": input_bytes,
                "decoded_bytes": decoded_bytes,
                "decoded_mib": decoded_mib,
                "throughput_mibs": throughput,
                "compressed_throughput_mibs": compressed_throughput,
                "command": command,
            }
        )
    return pd.DataFrame(rows)


def load_section(results_dir: Path, manifest: dict, key: str) -> pd.DataFrame | None:
    rel = (manifest.get("sections") or {}).get(key)
    if not rel:
        return None
    section_dir = results_dir / rel
    if not section_dir.is_dir():
        return None
    metadata = load_metadata(results_dir, manifest)
    frames = [load_zebrac_json(path, metadata) for path in sorted(section_dir.glob("*.json"))]
    frames = [frame for frame in frames if not frame.empty]
    return pd.concat(frames, ignore_index=True) if frames else None


def enrich_size(df: pd.DataFrame | None) -> pd.DataFrame | None:
    if df is None or df.empty:
        return df
    out = df.copy()
    out["size_mb"] = out["workload"].astype(str).str.replace("mb", "", regex=False).astype(float).astype(int)
    return out


def enrich_reads(df: pd.DataFrame | None) -> pd.DataFrame | None:
    if df is None or df.empty:
        return df
    out = df.copy()
    out["read_count"] = out["workload"].astype(int)
    return out


def build_ratio_comparisons(
    work: pd.DataFrame,
    value_col: str,
    *,
    baseline: str,
    peers: list[str],
    group_col: str = "dataset",
    group_sort=None,
) -> pd.DataFrame:
    if group_sort is None:
        group_sort = dataset_sort_key if group_col == "dataset" else lambda x: x
    groups = sorted(work[group_col].unique(), key=group_sort)
    rows: list[dict] = []
    for group in groups:
        base = work[(work[group_col] == group) & (work["tool"] == baseline)]
        if base.empty:
            continue
        base_v = float(base[value_col].values[0])
        for tool in peers:
            if tool == baseline:
                continue
            hit = work[(work[group_col] == group) & (work["tool"] == tool)]
            if hit.empty:
                continue
            peer_v = float(hit[value_col].values[0])
            ratio = peer_v / base_v if base_v > 0 else None
            rows.append(
                {
                    "dataset": group,
                    "tool": tool,
                    "competitor": display_tool(tool),
                    "zfasta_v": base_v,
                    "comp_v": peer_v,
                    "ratio": ratio,
                }
            )
    return pd.DataFrame(rows)


def md_pivot(df: pd.DataFrame, tools: list[str], index_col: str, fmt, *, index_order=None) -> str:
    work = df[df["tool"].isin(tools)].copy()
    if work.empty:
        return "_No data._"
    work["cell"] = work.apply(fmt, axis=1)
    pivot = work.pivot(index=index_col, columns="tool", values="cell")
    cols = [c for c in tools if c in pivot.columns]
    pivot = pivot[cols]
    if index_order is not None:
        pivot = pivot.reindex([x for x in index_order if x in pivot.index])
    elif index_col == "dataset":
        pivot = pivot.reindex([d for d in DATASET_ORDER if d in pivot.index])
    pivot = pivot.rename(columns={c: display_tool(c) for c in pivot.columns})
    if index_col == "dataset":
        pivot.index.name = "Dataset"
    return to_markdown_aligned(pivot)


def md_ratio_table(
    comparisons: pd.DataFrame,
    *,
    zf_label: str,
    peer_label: str,
    ratio_label: str,
    fmt_zf,
    fmt_comp,
    group_label: str = "Dataset",
) -> str:
    if comparisons.empty:
        return "_No comparisons._"
    rows = []
    for row in comparisons.itertuples(index=False):
        rows.append(
            {
                group_label: row.dataset,
                "z-fastq vs": row.competitor,
                zf_label: fmt_zf(row),
                peer_label: fmt_comp(row),
                ratio_label: format_ratio(row.ratio),
            }
        )
    return to_markdown_aligned(pd.DataFrame(rows), index=False)


def _std_col(value_col: str) -> str | None:
    if value_col == "mean":
        return "stddev"
    if value_col == "peak_rss_mb":
        return "peak_rss_stddev_mb"
    if value_col == "minor_faults":
        return "minor_faults_stddev"
    return None


def _bar_patches(tools: list[str]) -> list:
    patches = []
    for tool in tools:
        color = COLORS.get(tool, "#888888")
        kw: dict = {
            "facecolor": color,
            "edgecolor": color if tool in REFERENCE_TOOLS else "none",
            "label": display_tool(tool),
            "alpha": 0.75 if tool in REFERENCE_TOOLS else 0.88,
        }
        if tool in REFERENCE_TOOLS:
            kw["hatch"] = "///"
        patches.append(mpatches.Patch(**kw))
    return patches


def _annotate_ratios(
    ax,
    dataset,
    tools,
    bar_tops,
    comparisons,
    baseline,
    width,
    *,
    skip_near_one: bool = False,
) -> None:
    base_color = COLORS.get(baseline, "#F7A41D")
    base_key = (dataset, baseline)
    if base_key in bar_tops:
        bx, by = bar_tops[base_key]
        xs = [bar_tops[(dataset, t)][0] for t in tools if (dataset, t) in bar_tops]
        if xs:
            ax.hlines(
                by,
                min(xs) - width * 0.65,
                max(xs) + width * 0.65,
                colors=base_color,
                linestyles=(0, (4, 3)),
                linewidth=1.0,
                alpha=0.45,
                zorder=1,
            )
        ax.annotate(
            "1x",
            (bx, by),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            rotation=90,
            fontsize=7,
            fontweight="bold",
            color=base_color,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor=base_color, linewidth=0.9),
            zorder=5,
        )
    if comparisons.empty:
        return
    for row in comparisons.itertuples(index=False):
        key = (dataset, row.tool)
        if key not in bar_tops:
            continue
        if skip_near_one and row.ratio is not None and 0.92 <= float(row.ratio) <= 1.08:
            continue
        xpos, mean_t = bar_tops[key]
        color = COLORS.get(row.tool, "#666666")
        ax.annotate(
            format_ratio(row.ratio),
            (xpos, mean_t),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            rotation=90,
            fontsize=7,
            fontweight="bold",
            color="#111111",
            bbox=dict(
                boxstyle="round,pad=0.18",
                facecolor="white",
                edgecolor=color,
                linewidth=0.9,
                alpha=0.96,
            ),
            zorder=5,
        )


def _draw_metric_facet(ax, work, tools, datasets, *, value_col, ylabel, value_floor, log_y, baseline) -> None:
    std_col = _std_col(value_col)
    width = min(0.095, 0.80 / max(1, len(tools)))
    x = list(range(len(datasets)))
    bar_tops: dict[tuple[str, str], tuple[float, float]] = {}
    for ti, tool in enumerate(tools):
        color = COLORS.get(tool, "#888888")
        is_ref = tool in REFERENCE_TOOLS
        for di, ds in enumerate(datasets):
            row = work[(work["dataset"] == ds) & (work["tool"] == tool)]
            if row.empty:
                continue
            val = max(float(row[value_col].iloc[0]), value_floor)
            std = 0.0
            if std_col and std_col in row.columns and pd.notna(row[std_col].iloc[0]):
                std = max(float(row[std_col].iloc[0]), 0.0)
            xpos = di + (ti - len(tools) / 2 + 0.5) * width
            bar_kw: dict = {"width": width, "color": color, "alpha": 0.75 if is_ref else 0.88, "zorder": 2}
            if is_ref:
                bar_kw.update({"hatch": "///", "edgecolor": color, "linewidth": 0.6})
            if std > 0:
                bar_kw.update(
                    {
                        "yerr": std,
                        "capsize": 2,
                        "error_kw": {"elinewidth": 0.8, "capthick": 0.8, "ecolor": "#333333"},
                    }
                )
            ax.bar(xpos, val, **bar_kw)
            bar_tops[(ds, tool)] = (xpos, val)
    peers = [t for t in peer_tools(tools, baseline) if t not in NO_RSS_RATIO or value_col != "peak_rss_mb"]
    for ds in datasets:
        ds_work = work[work["dataset"] == ds]
        present = [t for t in tools if (ds, t) in bar_tops]
        if not present:
            continue
        comparisons = build_ratio_comparisons(
            ds_work,
            value_col,
            baseline=baseline,
            peers=[t for t in peers if t in present],
        )
        _annotate_ratios(
            ax,
            ds,
            present,
            bar_tops,
            comparisons,
            baseline,
            width,
            skip_near_one=(value_col == "peak_rss_mb"),
        )
    if log_y:
        ax.set_yscale("log")
        lo, hi = ax.get_ylim()
        if hi > lo > 0:
            ax.set_ylim(lo, hi * 2.0)
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.grid(axis="y", alpha=0.28, which="both")
    ax.set_axisbelow(True)


def fig_metric_facets(work: pd.DataFrame, out: Path, tools: list[str], title: str, note: str) -> Path:
    filtered = work[work["tool"].isin(tools)].copy()
    legend_tools = [t for t in tools if t in filtered["tool"].unique()]
    datasets = [d for d in DATASET_ORDER if d in filtered["dataset"].unique()]
    facets = [
        ("mean", "Wall Time (s)", True, 1e-6),
        ("peak_rss_mb", "Peak RSS (MB)", False, 0.1),
        ("minor_faults", "Minor Page Faults", True, 1.0),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.8), squeeze=False)
    for ax, (value_col, ylabel, log_y, floor) in zip(axes[0], facets):
        _draw_metric_facet(
            ax,
            filtered,
            legend_tools,
            datasets,
            value_col=value_col,
            ylabel=ylabel,
            value_floor=floor,
            log_y=log_y,
            baseline=BASELINE,
        )
    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.98)
    fig.text(0.5, 0.93, note, ha="center", va="top", fontsize=9, color="#444444", style="italic")
    patches = _bar_patches(legend_tools)
    fig.legend(
        handles=patches,
        labels=[p.get_label() for p in patches],
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.06),
        ncol=min(len(legend_tools), 7),
        frameon=False,
        columnspacing=1.0,
        handletextpad=0.4,
    )
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.18, top=0.86, wspace=FACET_WSPACE)
    return _save(fig, out)


def _scaling_style(tool: str) -> dict:
    is_ref = tool in REFERENCE_TOOLS
    is_native = tool == "z-fastq-native"
    return {
        "color": COLORS.get(tool, "#888888"),
        "marker": "D" if is_native else "o",
        "linestyle": "--" if is_ref else "-",
        "linewidth": 2.4 if tool == BASELINE else (1.9 if is_native else 1.8),
        "markersize": 5.5 if is_native else 6,
        "alpha": 0.92 if is_native else (0.75 if is_ref else 0.95),
        "zorder": 3 if tool == BASELINE else (2.6 if is_native else 2),
    }


SCALING_METRICS = (
    ("mean", "Wall time (s)", True, 1e-6),
    ("peak_rss_mb", "Peak RSS (MB)", True, 0.1),
    ("minor_faults", "Minor page faults", True, 1.0),
)


def fig_scaling_abs(df: pd.DataFrame, out: Path, *, param_col: str, xlabel: str, title: str, tools: list[str]) -> Path:
    work = df[df["tool"].isin(tools)].copy()
    tools = [t for t in tools if t in work["tool"].unique()]
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.4), squeeze=False)
    for ax, (col, ylab, log_y, floor) in zip(axes[0], SCALING_METRICS):
        for tool in tools:
            tdf = work[work["tool"] == tool].sort_values(param_col)
            if tdf.empty:
                continue
            st = _scaling_style(tool)
            ax.plot(
                tdf[param_col].astype(float),
                tdf[col].astype(float).clip(lower=floor),
                marker=st["marker"],
                linestyle=st["linestyle"],
                color=st["color"],
                linewidth=st["linewidth"],
                markersize=st["markersize"],
                markeredgecolor="white",
                markeredgewidth=0.6,
                alpha=st["alpha"],
                zorder=st["zorder"],
            )
        ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")
        ax.grid(alpha=0.28, which="both")
        ax.set_axisbelow(True)
        ax.set_ylabel(ylab, fontsize=10)
        ax.set_title(ylab, fontsize=11, fontweight="bold", pad=8)
        ax.set_xlabel(xlabel, fontsize=10)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.18, top=0.84, wspace=FACET_WSPACE)
    pos = axes[0, 0].get_position()
    fig.text(pos.x0, 0.965, title, ha="left", va="top", fontsize=12, fontweight="bold")
    fig.text(
        pos.x0,
        0.922,
        "Log-log absolute values. Circle = ISA-L; diamond = native inflate.",
        ha="left",
        va="top",
        fontsize=8.5,
        color="#444444",
        style="italic",
    )
    handles = []
    for tool in tools:
        st = _scaling_style(tool)
        handles.append(
            mlines.Line2D(
                [],
                [],
                color=st["color"],
                marker=st["marker"],
                linestyle=st["linestyle"],
                linewidth=st["linewidth"],
                markersize=st["markersize"],
                markeredgecolor="white",
                markeredgewidth=0.6,
                alpha=st["alpha"],
                label=display_tool(tool),
            )
        )
    fig.legend(
        handles=handles,
        fontsize=8.5,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.08),
        ncol=min(len(tools), 6),
        frameon=False,
    )
    return _save(fig, out)


def fig_scaling_slope(df: pd.DataFrame, out: Path, *, param_col: str, title: str, tools: list[str]) -> Path:
    work = df[df["tool"].isin(tools)].copy()
    tools = [t for t in tools if t in work["tool"].unique() and t != BASELINE]
    xs = sorted(work[param_col].unique())
    if len(xs) < 2:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "Need at least two scale points", ha="center")
        return _save(fig, out)
    x0, x1 = xs[0], xs[-1]
    fig, axes = plt.subplots(1, 3, figsize=(15, 6.2), squeeze=False)
    for ax, (col, ylab, _log_y, floor) in zip(axes[0], SCALING_METRICS):
        b0 = work[(work["tool"] == BASELINE) & (work[param_col] == x0)]
        b1 = work[(work["tool"] == BASELINE) & (work[param_col] == x1)]
        if b0.empty or b1.empty:
            continue
        base0 = float(b0[col].iloc[0])
        base1 = float(b1[col].iloc[0])
        for tool in tools:
            r0 = work[(work["tool"] == tool) & (work[param_col] == x0)]
            r1 = work[(work["tool"] == tool) & (work[param_col] == x1)]
            if r0.empty or r1.empty:
                continue
            v0 = max(float(r0[col].iloc[0]), floor) / max(base0, floor)
            v1 = max(float(r1[col].iloc[0]), floor) / max(base1, floor)
            st = _scaling_style(tool)
            ax.plot(
                [0, 1],
                [v0, v1],
                color=st["color"],
                linewidth=2.0,
                marker=st["marker"],
                markersize=7,
                markeredgecolor="white",
                markeredgewidth=0.7,
                alpha=0.92,
            )
            ax.text(1.04, v1, display_tool(tool), fontsize=7.5, va="center", color=st["color"])
        ax.axhline(1.0, color=COLORS[BASELINE], linestyle=(0, (4, 3)), linewidth=1.2, alpha=0.7)
        ax.set_yscale("log")
        ax.set_xlim(-0.08, 1.55)
        ax.set_xticks([0, 1])
        if param_col == "size_mb":
            ax.set_xticklabels([f"{int(x0)} MB", f"{int(x1)} MB"])
        else:
            ax.set_xticklabels([f"{int(x0):,}", f"{int(x1):,}"])
        ax.set_ylabel("x vs z-fastq ISA-L", fontsize=10)
        ax.set_title(ylab, fontsize=11, fontweight="bold", pad=8)
        ax.grid(axis="y", alpha=0.28, which="both")
        ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.07, right=0.90, bottom=0.12, top=0.84, wspace=0.32)
    pos = axes[0, 0].get_position()
    fig.text(pos.x0, 0.965, title, ha="left", va="top", fontsize=12, fontweight="bold")
    fig.text(
        pos.x0,
        0.922,
        "Left = smallest fixture; right = largest. Gold dashed = 1x. Rising = the peer gets relatively slower or heavier.",
        ha="left",
        va="top",
        fontsize=8.5,
        color="#444444",
        style="italic",
    )
    return _save(fig, out)


def md_metric_block(
    *,
    heading: str,
    intro: str,
    work: pd.DataFrame,
    tools: list[str],
    value_col: str,
    fmt,
    nums: ReportCounters,
    ratio_summary: str,
    zf_label: str,
    peer_label: str,
    ratio_label: str,
    fmt_zf,
    fmt_comp,
    skip_rss_tools: bool = False,
) -> str:
    t_main = nums.next_table()
    pivot = md_pivot(work, tools, "dataset", fmt)
    peers = peer_tools(tools)
    if skip_rss_tools:
        peers = [t for t in peers if t not in NO_RSS_RATIO]
    comparisons = build_ratio_comparisons(work, value_col, baseline=BASELINE, peers=peers)
    t_ratio = nums.next_table()
    ratio = md_ratio_table(
        comparisons,
        zf_label=zf_label,
        peer_label=peer_label,
        ratio_label=ratio_label,
        fmt_zf=fmt_zf,
        fmt_comp=fmt_comp,
    )
    return join(
        [
            f"## {heading}",
            "",
            intro,
            "",
            f"**Table {t_main}:** Mean ± stddev by dataset and tool.",
            "",
            pivot,
            "",
            f"<details><summary><strong>Table {t_ratio}:</strong> {ratio_summary}</summary>",
            "",
            ratio,
            "",
            "</details>",
            "",
            '<div style="margin: 1.5em 0"></div>',
            "",
        ]
    )


def md_figure_block(nums: ReportCounters, rel: str, caption: str, reading: list[str]) -> str:
    fnum = nums.next_figure()
    bullets = "\n".join(f"- {item}" for item in reading)
    return join(
        [
            f"**Figure {fnum}:** {caption}",
            "",
            f"![Figure {fnum}]({rel})",
            "",
            f"**Reading Figure {fnum}**",
            "",
            bullets,
            "",
            '<div style="margin: 1.5em 0"></div>',
            "",
        ]
    )


def md_overview(manifest: dict) -> str:
    real_set = manifest.get("real_set", "full")
    dense_note = (
        "Dense in this run is MiniSeq R1 SRR10325788 (bring-up set; not the large SRR1810900 file)."
        if real_set == "small"
        else "Dense is the large fixed-length Illumina file SRR1810900 (24.5 million 50 bp reads)."
    )
    verify = manifest.get("verify_pass")
    if manifest.get("verify_skipped"):
        verify_line = "The count-agreement check was skipped (`--skip-tests`). Do not treat this as a published run."
    elif verify:
        verify_line = (
            "seqtk `size`, Needletail, Helicase, fqtools, and both z-fastq binaries agreed on the record count "
            "for every timed file before any duration measurement. Details are in Count agreement below. "
            "I am not repeating those integers here."
        )
    else:
        verify_line = "The count-agreement check was not recorded."
    return join(
        [
            "This report times `z-fastq count` on matched plain and gzip FASTQ. The three real inputs "
            "are a dense fixed-length Illumina run, a variable-length Illumina run, and an ONT long-read run. "
            + dense_note,
            "",
            "**What is timed**",
            "",
            "- One process, one file. Zebrac starts that argv; it does not start a shell.",
            "- Plain files use the ISA-L product binary only. Gzip files also time the native inflate binary (`-Disa-l=false`).",
            "- seqtk `size`, Needletail, Helicase, and fqtools are complete count peers. "
            "SeqFu may launch helper processes; its wall time is shown and its RSS is not used for ratios.",
            "",
            verify_line,
        ]
    )


def md_capability() -> str:
    matrix = pd.DataFrame(
        [
            {
                "Tool": "z-fastq (ISA-L)",
                "Command": "`z-fastq count`",
                "Same count job": "yes",
                "Timed as": "product default",
            },
            {
                "Tool": "z-fastq (native)",
                "Command": "`z-fastq-native count`",
                "Same count job": "yes, gzip only",
                "Timed as": "second inflate backend",
            },
            {
                "Tool": "Needletail",
                "Command": "`needletail-adapter count`",
                "Same count job": "yes, same one-integer line",
                "Timed as": "equal-work wrapper",
            },
            {
                "Tool": "Helicase",
                "Command": "`helicase-adapter count`",
                "Same count job": "yes, same one-integer line",
                "Timed as": "equal-work wrapper",
            },
            {
                "Tool": "seqtk",
                "Command": "`seqtk size`",
                "Same count job": "reads field; also prints bases",
                "Timed as": "well-known reference, complete peer",
            },
            {
                "Tool": "SeqFu",
                "Command": "`seqfu count`",
                "Same count job": "numeric count on one file",
                "Timed as": "hatched; RSS not comparable",
            },
            {
                "Tool": "fqtools",
                "Command": "`fqtools count`",
                "Same count job": "yes on four-line files",
                "Timed as": "complete peer",
            },
        ]
    )
    return join(
        [
            "`z-fastq count` prints one decimal record count. That is the compared fact.",
            "",
            to_markdown_aligned(matrix, index=False),
            "",
            "seqtk stdout is `N` then a tab then bases. The agreement check parses `N`. It does not rewrite seqtk to look like z-fastq. SeqFu is not part of that check.",
        ]
    )


def md_correctness(manifest: dict) -> str:
    log = manifest.get("verify_log", "")
    status = manifest.get("verify_pass") or "not recorded"
    skipped = bool(manifest.get("verify_skipped"))
    if skipped:
        body = "This run skipped the count-agreement check."
    else:
        body = (
            f"Status: **{status}**. The check required byte-identical count lines from z-fastq, "
            "z-fastq native (gzip files), Needletail, Helicase, and the independent four-line counter, plus a matching "
            "seqtk reads field and fqtools integer. The first mismatch would have stopped the run."
        )
    return join(
        [
            "Valid four-line FASTQ only. seqtk is not an error-status reference. Malformed fixtures are a "
            "separate z-fastq reject check.",
            "",
            body,
            "",
            f"Log: `{log}`." if log else "No verify log name was recorded.",
        ]
    )


def md_provenance(manifest: dict) -> str:
    tools = manifest.get("tools") or {}
    tool_lines = [
        f"- **{name}:** {short_peer_version(name, tools[name])}"
        for name in ("needletail", "helicase", "seqtk", "seqfu", "fqtools")
        if tools.get(name)
    ]
    isa_bytes = int(manifest.get("z_fastq_bytes") or 0)
    nat_bytes = int(manifest.get("z_fastq_native_bytes") or 0)
    isa_line = f"- **z-fastq ISA-L:** {manifest.get('z_fastq', '')}"
    nat_line = f"- **z-fastq native:** {manifest.get('z_fastq_native', '')}"
    if isa_bytes:
        isa_line += f" ({isa_bytes:,} bytes)"
    if nat_bytes:
        nat_line += f" ({nat_bytes:,} bytes)"
    real_set = manifest.get("real_set", "full")
    dense_line = (
        "Dense in this run is MiniSeq R1 (bring-up). Variable is ERR164407. Long is ERR5404926."
        if real_set == "small"
        else "Dense is SRR1810900 (fixed 50 bp Illumina). Variable is ERR164407 (variable Illumina). Long is ERR5404926 (ONT)."
    )
    sections = manifest.get("sections") or {}
    has_scale = any(key.startswith("scale_") for key in sections)
    scale_note = (
        " Synthetic size files hold 100k uniform reads. Synthetic read-count files hold 150 bp reads. "
        "Gzip copies of the read-count files use compresslevel 6."
        if has_scale
        else " This run skipped the synthetic size and read-count sweeps."
    )
    verify = manifest.get("verify_pass")
    log = manifest.get("verify_log") or ""
    if manifest.get("verify_skipped"):
        check_line = "- **Count check:** skipped (`--skip-tests`)"
    elif verify and log:
        check_line = f"- **Count check:** {verify} (`{log}`)"
    elif verify:
        check_line = f"- **Count check:** {verify}"
    else:
        check_line = "- **Count check:** not recorded"
    parts = [
        f"- **Timestamp:** `{manifest.get('timestamp', '')}`",
        f"- **Runner:** zebrac, warm page cache, runs={manifest.get('runs')}, warmup={manifest.get('warmup')}, duration_ms={manifest.get('duration_ms')}",
        f"- **zebrac:** {manifest.get('zebrac', '')}",
        isa_line,
        nat_line,
        check_line,
        "",
        "**Peer versions**",
        "",
    ]
    parts += tool_lines or ["- none recorded"]
    parts += ["", dense_line + scale_note]
    return join(parts)


def md_perf_section(
    df: pd.DataFrame,
    *,
    title: str,
    intro: str,
    tools: list[str],
    fig_name: str,
    fig_title: str,
    fig_note: str,
    nums: ReportCounters,
    figures_dir: Path,
    include_throughput: bool = False,
) -> str:
    present = tools_in_run(df, tools)
    fig_metric_facets(df, figures_dir / fig_name, present, fig_title, fig_note)
    blocks = [
        f"## {title}",
        "",
        intro,
        "",
        md_metric_block(
            heading="Wall time",
            intro="Zebrac mean wall time for one `count` / peer process with stdout discarded.",
            work=df,
            tools=present,
            value_col="mean",
            fmt=fmt_wall,
            nums=nums,
            ratio_summary="Time x = lane wall / z-fastq ISA-L. Same ratios as bar labels.",
            zf_label="z-fastq",
            peer_label="Peer",
            ratio_label="Time x",
            fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
            fmt_comp=lambda r: f"{r.comp_v:.4f}s",
        ),
    ]
    if include_throughput:
        t_dec = nums.next_table()
        t_comp = nums.next_table()
        blocks.extend(
            [
                "## Decoded throughput",
                "",
                "Decoded FASTQ MiB/s is uncompressed sibling size divided by mean wall. That is the gzip unit that matches inflate work. Compressed MiB/s (on-disk gzip size / wall) is under the fold.",
                "",
                f"**Table {t_dec}:** Decoded FASTQ MiB/s.",
                "",
                md_pivot(df, present, "dataset", fmt_decoded_mibs),
                "",
                f"<details><summary><strong>Table {t_comp}:</strong> Compressed on-disk MiB/s (supporting column).</summary>",
                "",
                md_pivot(df, present, "dataset", fmt_compressed_mibs),
                "",
                "</details>",
                "",
                '<div style="margin: 1.5em 0"></div>',
                "",
            ]
        )
    blocks.extend(
        [
            md_metric_block(
                heading="Memory",
                intro="Peak RSS of the timed process. Zebrac starts the tool directly and does not start a shell. Count streams, so RSS does not grow with file size. The scanner does not heap-allocate; remaining RSS is the process image plus a 256 KiB read buffer. SeqFu is omitted from RSS ratios because child RSS is not measured.",
                work=df,
                tools=present,
                value_col="peak_rss_mb",
                fmt=fmt_rss,
                nums=nums,
                ratio_summary="RSS x = lane peak RSS / z-fastq ISA-L.",
                zf_label="z-fastq",
                peer_label="Peer",
                ratio_label="RSS x",
                fmt_zf=lambda r: f"{r.zfasta_v:.2f} MB",
                fmt_comp=lambda r: f"{r.comp_v:.2f} MB",
                skip_rss_tools=True,
            ),
            md_metric_block(
                heading="Page faults",
                intro="Minor page faults for the same commands.",
                work=df,
                tools=present,
                value_col="minor_faults",
                fmt=fmt_faults,
                nums=nums,
                ratio_summary="Faults x = lane minor faults / z-fastq ISA-L.",
                zf_label="z-fastq",
                peer_label="Peer",
                ratio_label="Faults x",
                fmt_zf=lambda r: f"{r.zfasta_v:.0f}",
                fmt_comp=lambda r: f"{r.comp_v:.0f}",
            ),
            md_figure_block(
                nums,
                f"results/figures/{fig_name}",
                "Wall, peak RSS, and minor page faults for the real datasets.",
                [
                    "Facets: wall time (log y) | peak RSS (linear) | minor page faults (log y).",
                    "X-axis: Dense, Variable, Long.",
                    (
                        "Gold bars: z-fastq ISA-L (1x). Amber bars: native inflate."
                        if include_throughput
                        else "Gold bars: z-fastq ISA-L (1x). One z-fastq lane on plain FASTQ."
                    ),
                    "Hatched bars: SeqFu (descriptive; helper processes). They are not dominance comparisons.",
                    "Error bars show standard deviation.",
                ],
            ),
        ]
    )
    return join(blocks)


def _scaling_pivot(df: pd.DataFrame, tools: list[str], param_col: str, param_order: list, xlabel: str, fmt) -> str:
    work = df[df["tool"].isin(tools)].copy()
    if work.empty:
        return "_No data._"
    work["cell"] = work.apply(fmt, axis=1)
    pivot = work.pivot(index=param_col, columns="tool", values="cell")
    pivot = pivot.reindex([p for p in param_order if p in pivot.index])
    pivot = pivot[[c for c in tools if c in pivot.columns]]
    pivot = pivot.rename(columns={c: display_tool(c) for c in pivot.columns})
    if param_col == "size_mb":
        pivot.index = [f"{int(v)} MB" for v in pivot.index]
    else:
        pivot.index = [f"{int(v):,}" for v in pivot.index]
    pivot.index.name = xlabel
    return to_markdown_aligned(pivot)


def _scaling_ratio(
    df: pd.DataFrame,
    tools: list[str],
    param_col: str,
    value_col: str,
    *,
    fmt_zf,
    fmt_comp,
    ratio_label: str,
    xlabel: str,
) -> str:
    comparisons = build_ratio_comparisons(
        df,
        value_col,
        baseline=BASELINE,
        peers=peer_tools(tools),
        group_col=param_col,
        group_sort=lambda x: float(x),
    )
    if comparisons.empty:
        return "_No comparisons._"
    comparisons = comparisons.copy()
    comparisons["dataset"] = comparisons["dataset"].map(
        lambda v: f"{int(v)} MB" if param_col == "size_mb" else f"{int(v):,}"
    )
    return md_ratio_table(
        comparisons,
        zf_label="z-fastq",
        peer_label="Peer",
        ratio_label=ratio_label,
        fmt_zf=fmt_zf,
        fmt_comp=fmt_comp,
        group_label=xlabel,
    )


def md_scaling(
    df: pd.DataFrame,
    *,
    param_col: str,
    param_order: list,
    xlabel: str,
    title: str,
    fig_stem: str,
    fig_title: str,
    tools: list[str],
    nums: ReportCounters,
    figures_dir: Path,
) -> str:
    present = tools_in_run(df, tools)
    fig_scaling_abs(
        df,
        figures_dir / f"{fig_stem}.png",
        param_col=param_col,
        xlabel=xlabel,
        title=fig_title,
        tools=present,
    )
    fig_scaling_slope(
        df,
        figures_dir / f"{fig_stem}_slope.png",
        param_col=param_col,
        title=f"{fig_title}: x at smallest vs largest",
        tools=present,
    )
    t_wall = nums.next_table()
    t_wall_ratio = nums.next_table()
    t_rss = nums.next_table()
    t_rss_ratio = nums.next_table()
    t_pf = nums.next_table()
    t_pf_ratio = nums.next_table()
    f_abs = nums.next_figure()
    f_slope = nums.next_figure()
    return join(
        [
            f"## {title}",
            "",
            "Complete count peers only. SeqFu is excluded so scale lines stay equivalent-work.",
            "",
            f"**Table {t_wall}:** Mean ± stddev wall time.",
            "",
            _scaling_pivot(df, present, param_col, param_order, xlabel, fmt_wall),
            "",
            f"<details><summary><strong>Table {t_wall_ratio}:</strong> Time x = peer / z-fastq ISA-L at each point.</summary>",
            "",
            _scaling_ratio(
                df,
                present,
                param_col,
                "mean",
                fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
                fmt_comp=lambda r: f"{r.comp_v:.4f}s",
                ratio_label="Time x",
                xlabel=xlabel,
            ),
            "",
            "</details>",
            "",
            f"<details><summary><strong>Table {t_rss}:</strong> Peak RSS (MB) by {xlabel}.</summary>",
            "",
            _scaling_pivot(df, present, param_col, param_order, xlabel, fmt_rss),
            "",
            "</details>",
            "",
            f"<details><summary><strong>Table {t_rss_ratio}:</strong> RSS x = peer / z-fastq ISA-L.</summary>",
            "",
            _scaling_ratio(
                df,
                present,
                param_col,
                "peak_rss_mb",
                fmt_zf=lambda r: f"{r.zfasta_v:.2f} MB",
                fmt_comp=lambda r: f"{r.comp_v:.2f} MB",
                ratio_label="RSS x",
                xlabel=xlabel,
            ),
            "",
            "</details>",
            "",
            f"<details><summary><strong>Table {t_pf}:</strong> Minor page faults by {xlabel}.</summary>",
            "",
            _scaling_pivot(df, present, param_col, param_order, xlabel, fmt_faults),
            "",
            "</details>",
            "",
            f"<details><summary><strong>Table {t_pf_ratio}:</strong> Faults x = peer / z-fastq ISA-L.</summary>",
            "",
            _scaling_ratio(
                df,
                present,
                param_col,
                "minor_faults",
                fmt_zf=lambda r: f"{r.zfasta_v:.0f}",
                fmt_comp=lambda r: f"{r.comp_v:.0f}",
                ratio_label="Faults x",
                xlabel=xlabel,
            ),
            "",
            "</details>",
            "",
            '<div style="margin: 1.5em 0"></div>',
            "",
            f"**Figure {f_abs}:** Absolute wall / RSS / page faults vs {xlabel} (log-log).",
            "",
            f"![Figure {f_abs}](results/figures/{fig_stem}.png)",
            "",
            f"**Reading Figure {f_abs}**",
            "",
            "- Facets: wall time | peak RSS | minor page faults.",
            "- Gold circle = ISA-L. Amber diamond = native inflate (gzip sweep only).",
            "",
            '<div style="margin: 1.5em 0"></div>',
            "",
            f"**Figure {f_slope}:** x vs z-fastq ISA-L at the smallest vs largest fixture.",
            "",
            f"![Figure {f_slope}](results/figures/{fig_stem}_slope.png)",
            "",
            f"**Reading Figure {f_slope}**",
            "",
            "- Left = smallest fixture; right = largest. Gold dashed = 1x.",
            "- Rising = the peer gets relatively worse with scale.",
            "",
            '<div style="margin: 1.5em 0"></div>',
        ]
    )


def md_run_banner(manifest: dict) -> str | None:
    notes: list[str] = []
    if manifest.get("real_set") == "small":
        notes.append("Dense is MiniSeq R1, not SRR1810900")
    try:
        duration = int(manifest.get("duration_ms") or 0)
        runs = int(manifest.get("runs") or 0)
        warmup = int(manifest.get("warmup") or 0)
    except (TypeError, ValueError):
        duration, runs, warmup = 0, 0, 0
    if (duration and duration < 5000) or (runs and runs < 25) or (warmup and warmup < 5):
        notes.append(
            f"zebrac duration_ms={duration}, runs={runs}, warmup={warmup} "
            "(published default is 5000 ms, 25 runs, 5 warmups)"
        )
    sections = manifest.get("sections") or {}
    if manifest.get("skip_scale") or not any(key.startswith("scale_") for key in sections):
        notes.append("synthetic sweeps skipped")
    if not notes:
        return None
    return "_This is a bring-up measurement: " + "; ".join(notes) + ". Read it as a check of the suite, not a published ranking._"


def squeeze_blank_lines(text: str) -> str:
    while "\n\n\n" in text:
        text = text.replace("\n\n\n", "\n\n")
    return text


def generate(results_dir: Path, allow_incomplete: bool) -> None:
    manifest = load_latest_manifest(results_dir)
    if not allow_incomplete and manifest.get("verify_skipped"):
        raise SystemExit("error: publishable report requires a passing count-agreement check")
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    nums = ReportCounters()
    banner = md_run_banner(manifest)
    lines = [
        MARKDOWNLINT_DISABLE,
        "",
        "# z-fastq Count Benchmark Report",
        "",
        "_Generated by `bench/count/generate_report.py` from zebrac results._",
        "",
    ]
    if banner:
        lines += [banner, ""]
    lines += [
        "## Overview",
        "",
        md_overview(manifest),
        "",
        "## What each tool does",
        "",
        md_capability(),
        "",
        "## Count agreement",
        "",
        md_correctness(manifest),
        "",
        "## Run provenance",
        "",
        md_provenance(manifest),
        "",
    ]

    plain = load_section(results_dir, manifest, "perf_plain")
    gzip_df = load_section(results_dir, manifest, "perf_gzip")
    size_df = enrich_size(load_section(results_dir, manifest, "scale_size"))
    reads_df = enrich_reads(load_section(results_dir, manifest, "scale_reads"))
    scale_gz = enrich_reads(load_section(results_dir, manifest, "scale_gzip"))

    required = []
    if not allow_incomplete:
        if plain is None:
            required.append("perf_plain")
        if gzip_df is None:
            required.append("perf_gzip")
    if required:
        raise SystemExit("error: missing sections: " + ", ".join(required) + " (pass --allow-incomplete for a draft)")

    sample_n = manifest.get("runs", "?")
    if plain is not None and not plain.empty:
        lines.append(
            md_perf_section(
                plain,
                title="Performance: plain FASTQ",
                intro=(
                    "Uncompressed siblings of the same records. One z-fastq lane (the ISA-L product binary). "
                    "Ratios use z-fastq ISA-L as 1x."
                ),
                tools=PLAIN_TOOLS,
                fig_name="perf_plain.png",
                fig_title="Plain FASTQ count: wall, RSS, page faults",
                fig_note=f"Error bars = zebrac stddev (n={sample_n}). Hatched = SeqFu.",
                nums=nums,
                figures_dir=figures_dir,
            )
        )
        lines.append("")
    if gzip_df is not None and not gzip_df.empty:
        lines.append(
            md_perf_section(
                gzip_df,
                title="Performance: gzip FASTQ",
                intro=(
                    "The same records as gzip. Gold bars are z-fastq ISA-L (product default, 1x). "
                    "Amber bars are native Zig inflate (`-Disa-l=false`). "
                    "Decoded FASTQ MiB/s (uncompressed sibling size / wall) is the gzip unit that matches inflate work."
                ),
                tools=GZIP_TOOLS,
                fig_name="perf_gzip.png",
                fig_title="Gzip FASTQ count: wall, RSS, page faults",
                fig_note=f"Error bars = zebrac stddev (n={sample_n}). Amber = native inflate. Hatched = SeqFu.",
                nums=nums,
                figures_dir=figures_dir,
                include_throughput=True,
            )
        )
        lines.append("")
    if size_df is not None and not size_df.empty:
        lines.append(
            md_scaling(
                size_df,
                param_col="size_mb",
                param_order=SIZE_MB_ORDER,
                xlabel="File size (MB)",
                title="Scaling: file size",
                fig_stem="scaling_size",
                fig_title="Count vs file size (100k uniform reads)",
                tools=SCALING_TOOLS,
                nums=nums,
                figures_dir=figures_dir,
            )
        )
        lines.append("")
    if reads_df is not None and not reads_df.empty:
        lines.append(
            md_scaling(
                reads_df,
                param_col="read_count",
                param_order=READ_COUNT_ORDER,
                xlabel="Read count",
                title="Scaling: read count",
                fig_stem="scaling_reads",
                fig_title="Count vs read count (150 bp dense)",
                tools=SCALING_TOOLS,
                nums=nums,
                figures_dir=figures_dir,
            )
        )
        lines.append("")
    if scale_gz is not None and not scale_gz.empty:
        lines.append(
            md_scaling(
                scale_gz,
                param_col="read_count",
                param_order=READ_COUNT_ORDER,
                xlabel="Read count",
                title="Scaling: gzip read count",
                fig_stem="scaling_gzip",
                fig_title="Gzip count vs read count (150 bp, compresslevel 6)",
                tools=SCALING_GZIP_TOOLS,
                nums=nums,
                figures_dir=figures_dir,
            )
        )
        lines.append("")

    report_path = SCRIPT_DIR / "REPORT.md"
    text = squeeze_blank_lines("\n".join(lines).rstrip() + "\n")
    report_path.write_text(text, encoding="utf-8")
    print(f"wrote {report_path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results_dir", nargs="?", default=str(RESULTS_DIR))
    ap.add_argument(
        "--allow-incomplete",
        default="false",
        help="true/false: allow missing perf sections or a skipped count-agreement check",
    )
    args = ap.parse_args()
    allow = str(args.allow_incomplete).lower() in {"1", "true", "yes"}
    generate(Path(args.results_dir), allow)


if __name__ == "__main__":
    main()
