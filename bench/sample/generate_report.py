#!/usr/bin/env python3
"""Generate bench/sample/REPORT.md and figures from zebrac JSON."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

MPL_CONFIG_DIR = Path(__file__).resolve().parents[2] / "tools" / ".local" / "matplotlib"
MPL_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CONFIG_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import pandas as pd
import tabulate  # noqa: F401 -- pd.to_markdown()

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"

BASELINE = "z-fastq"
PLAIN_TOOLS = ["z-fastq", "seqtk", "seqkit", "rasusa", "fq", "fqkit", "irma", "bbtools"]
GZIP_TOOLS = ["z-fastq", "z-fastq-native", "seqtk", "seqkit", "rasusa", "fq", "fqkit", "irma", "bbtools"]
DESCRIPTIVE_TOOLS = frozenset({"seqkit", "rasusa", "fq", "fqkit", "irma", "bbtools"})
NO_OCCUPANCY = frozenset({"bbtools", "seqkit"})
ROLE_MATRIX_ORDER = [
    "se_fraction",
    "se_count",
    "paired_fraction",
    "paired_count",
    "interleaved_fraction",
    "interleaved_count",
]
ROLE_MATRIX_LABELS = {
    "se_fraction": "SE fraction",
    "se_count": "SE count",
    "paired_fraction": "paired fraction",
    "paired_count": "paired count",
    "interleaved_fraction": "interleaved fraction",
    "interleaved_count": "interleaved count",
}
FAMILY_ROLES = {
    "se_fraction": ("se_fraction",),
    "se_count": ("se_count",),
    "pairs_fraction": ("paired_fraction", "interleaved_fraction"),
    "pairs_count": ("paired_count", "interleaved_count"),
}
FIXTURE_ORDER = [
    "screened_se_fraction",
    "screened_se_count",
    "slash_pair_fraction",
    "slash_pair_count",
    "interleaved_fraction",
    "interleaved_count",
]
FIXTURE_LABELS = {
    "screened_se_fraction": "screened SE fraction 0.5",
    "screened_se_count": "screened SE count 3",
    "slash_pair_fraction": "slash /1 /2 pair fraction 0.5",
    "slash_pair_count": "slash /1 /2 pair count 3",
    "interleaved_fraction": "interleaved pair fraction 0.5",
    "interleaved_count": "interleaved pair count 3",
}

COLORS = {
    "z-fastq": "#F7A41D",
    "z-fastq-native": "#FFB74D",
    "seqtk": "#6A1B9A",
    "seqkit": "#1565C0",
    "rasusa": "#2E7D32",
    "fq": "#0277BD",
    "fqkit": "#6D4C41",
    "irma": "#00838F",
    "bbtools": "#C62828",
}
DISPLAY = {
    "z-fastq": "z-fastq (ISA-L)",
    "z-fastq-native": "z-fastq (native)",
    "seqtk": "seqtk sample",
    "seqkit": "SeqKit sample",
    "rasusa": "Rasusa reads",
    "fq": "fq subsample",
    "fqkit": "fqkit subfq",
    "irma": "IRMA Core sampler",
    "bbtools": "BBTools reformat",
}

WORKLOAD_ORDER: list[str] = []
WORKLOAD_LABELS: dict[str, str] = {}
WORKLOAD_ROLES: dict[str, str] = {}
WORKLOAD_MARKERS = ("o", "s", "^", "D", "P", "X", "v", "<", ">", "h")
WORKLOAD_LINESTYLES = ("-", "--", "-.", ":", (0, (8, 3, 2, 3, 2, 3)), (0, (5, 2, 1, 2)))
WORKLOAD_HATCHES = ("", "//", "xx", "..", "\\\\", "++", "oo", "**", "--", "||")
MARKDOWNLINT_DISABLE = "<!-- markdownlint-disable MD024 MD032 MD033 MD036 MD041 MD049 -->"


class ReportCounters:
    def __init__(self) -> None:
        self.table = 1
        self.figure = 1

    def next_table(self) -> int:
        value = self.table
        self.table += 1
        return value

    def next_figure(self) -> int:
        value = self.figure
        self.figure += 1
        return value


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
    col_widths = [max(len(row[index]) for row in rows) for index in range(width)]
    output: list[str] = []
    for row_index, row in enumerate(rows):
        cells: list[str] = []
        for index, cell in enumerate(row):
            if row_index == 1 and cell and set(cell) <= set("-:"):
                left = cell.startswith(":")
                right = cell.endswith(":")
                if left and right:
                    cell = ":" + "-" * max(1, col_widths[index] - 2) + ":"
                elif left:
                    cell = ":" + "-" * max(1, col_widths[index] - 1)
                elif right:
                    cell = "-" * max(1, col_widths[index] - 1) + ":"
                else:
                    cell = "-" * max(3, col_widths[index])
            cells.append(cell.ljust(col_widths[index]))
        output.append("| " + " | ".join(cells) + " |")
    return "\n".join(output)


def to_markdown_aligned(frame: pd.DataFrame, *, index: bool = True) -> str:
    return align_pipe_table(frame.to_markdown(index=index))


def display_tool(tool: str) -> str:
    return DISPLAY.get(tool, tool)


def display_workload(workload: str) -> str:
    return WORKLOAD_LABELS.get(workload, workload)


def family_workloads(roles: tuple[str, ...]) -> list[str]:
    wanted = set(roles)
    return [key for key in WORKLOAD_ORDER if WORKLOAD_ROLES.get(key) in wanted]


def occupancy_tools(tools: list[str]) -> list[str]:
    return [tool for tool in tools if tool not in NO_OCCUPANCY]


def drop_disallowed_lanes(frame: pd.DataFrame) -> pd.DataFrame:
    """IRMA on one interleaved file is SE records, not pair units."""
    if frame is None or frame.empty:
        return frame
    roles = frame["workload"].map(lambda key: WORKLOAD_ROLES.get(str(key), ""))
    drop = (frame["tool"] == "irma") & roles.str.startswith("interleaved")
    return frame.loc[~drop].copy()


def overview_shapes(roles: tuple[str, ...]) -> str:
    prefixes = (
        "SE fraction: ",
        "SE count: ",
        "paired fraction: ",
        "paired count: ",
        "interleaved fraction: ",
        "interleaved count: ",
    )
    parts: list[str] = []
    for key in family_workloads(roles):
        label = display_workload(key)
        for prefix in prefixes:
            if label.startswith(prefix):
                label = label[len(prefix) :]
                break
        parts.append(label)
    return ", ".join(parts) if parts else "none"


def workload_style(workload: str) -> tuple[str, str, str]:
    try:
        index = WORKLOAD_ORDER.index(workload)
    except ValueError:
        index = 0
    return (
        WORKLOAD_MARKERS[index % len(WORKLOAD_MARKERS)],
        WORKLOAD_LINESTYLES[index % len(WORKLOAD_LINESTYLES)],
        WORKLOAD_HATCHES[index % len(WORKLOAD_HATCHES)],
    )


def format_ratio(value: float | None) -> str:
    if value is None:
        return "n/a"
    if abs(value - 1.0) < 1e-9:
        return "1x"
    if value >= 100:
        return f"{value:.0f}x"
    if value >= 10:
        return f"{value:.1f}x"
    if value >= 1:
        return f"{value:.2f}x"
    return f"{value:.3f}x"


def _cell_float(row, key: str, default: float = 0.0) -> float:
    try:
        value = row[key]
    except (KeyError, IndexError):
        return default
    if pd.isna(value):
        return default
    return float(value)


def fmt_wall(row) -> str:
    return f"{_cell_float(row, 'mean'):.3f}±{_cell_float(row, 'stddev'):.3f} s"


def fmt_rss(row) -> str:
    return f"{_cell_float(row, 'peak_rss_mb'):.2f}±{_cell_float(row, 'peak_rss_stddev_mb'):.2f} MB"


def fmt_faults(row) -> str:
    return f"{_cell_float(row, 'minor_faults'):.0f}±{_cell_float(row, 'minor_faults_stddev'):.0f}"


def fmt_decoded_mibs(row) -> str:
    value = row.get("throughput_mibs")
    if value is None or pd.isna(value):
        return "n/a"
    return f"{float(value):.1f} MiB/s"


def fmt_compressed_mibs(row) -> str:
    value = row.get("compressed_throughput_mibs")
    if value is None or pd.isna(value):
        return "n/a"
    return f"{float(value):.1f} MiB/s"


def _std_col(value_col: str) -> str | None:
    return {
        "mean": "stddev",
        "peak_rss_mb": "peak_rss_stddev_mb",
        "minor_faults": "minor_faults_stddev",
    }.get(value_col)


def _save(fig, path: Path, *, dpi: int = 150) -> Path:
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def load_latest_manifest(results_dir: Path) -> dict:
    latest = results_dir / "LATEST"
    if not latest.is_file():
        raise SystemExit(f"error: missing {latest}")
    timestamp = latest.read_text(encoding="utf-8").strip()
    path = results_dir / f"run_{timestamp}.json"
    if not path.is_file():
        raise SystemExit(f"error: missing {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def configure_workloads(manifest: dict) -> None:
    global WORKLOAD_ORDER, WORKLOAD_LABELS, WORKLOAD_ROLES
    rows = manifest.get("workloads") or []
    WORKLOAD_ORDER = []
    WORKLOAD_LABELS = {}
    WORKLOAD_ROLES = {}
    for row in rows:
        key = str(row["key"])
        role = str(row.get("role", "workload"))
        ids = [str(value) for value in row.get("ids", [])]
        if key in WORKLOAD_LABELS:
            continue
        if role.startswith("se_"):
            kind = "fraction" if "fraction" in role else "count"
            label = f"SE {kind}: " + ", ".join(ids)
        elif role.startswith("paired_"):
            kind = "fraction" if "fraction" in role else "count"
            label = f"paired {kind}: " + " + ".join(ids)
        elif role.startswith("interleaved_"):
            kind = "fraction" if "fraction" in role else "count"
            label = f"interleaved {kind}: " + ", ".join(ids)
        else:
            label = f"{role}: " + " + ".join(ids)
        WORKLOAD_ORDER.append(key)
        WORKLOAD_LABELS[key] = label
        WORKLOAD_ROLES[key] = role


def load_metadata(results_dir: Path, manifest: dict) -> pd.DataFrame | None:
    name = manifest.get("metadata")
    if not name:
        return None
    path = results_dir / name
    if not path.is_file():
        return None
    rows = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    return pd.DataFrame(rows) if rows else None


def load_zebrac_json(path: Path, metadata: pd.DataFrame | None) -> pd.DataFrame:
    data = json.loads(path.read_text(encoding="utf-8"))
    meta_rows: list[dict] = []
    if metadata is not None and not metadata.empty:
        matches = metadata[metadata["raw_json"] == str(path)]
        if matches.empty:
            matches = metadata[metadata["raw_json"].astype(str).str.endswith(path.name)]
        meta_rows = matches.to_dict("records")
    by_command = {row.get("command"): row for row in meta_rows}
    rows: list[dict] = []
    for index, result in enumerate(data.get("results", [])):
        command = result.get("command", "")
        meta = by_command.get(command, {})
        if not meta and index < len(meta_rows):
            meta = meta_rows[index]
        wall = result.get("wall_time", {})
        peak = result.get("peak_rss", {})
        faults = result.get("minor_faults", {})
        mean = float(wall.get("mean", 0) or 0) / 1_000_000_000.0
        stddev = float(wall.get("std_dev", 0) or 0) / 1_000_000_000.0
        rss = float(peak.get("mean", 0) or 0) / (1024.0 * 1024.0)
        rss_stddev = float(peak.get("std_dev", 0) or 0) / (1024.0 * 1024.0)
        input_bytes = meta.get("input_bytes")
        decoded_bytes = meta.get("decoded_bytes")
        decoded_mib = None if decoded_bytes is None else float(decoded_bytes) / (1024.0 * 1024.0)
        compressed_mib = None if input_bytes is None else float(input_bytes) / (1024.0 * 1024.0)
        rows.append(
            {
                "tool": meta.get("tool") or path.stem.split("__", 1)[-1],
                "section": meta.get("section", ""),
                "workload": meta.get("workload", path.stem.split("__", 1)[0]),
                "mean": mean,
                "stddev": stddev,
                "peak_rss_mb": rss,
                "peak_rss_stddev_mb": rss_stddev,
                "minor_faults": float(faults.get("mean", 0) or 0),
                "minor_faults_stddev": float(faults.get("std_dev", 0) or 0),
                "input_bytes": input_bytes,
                "decoded_bytes": decoded_bytes,
                "decoded_mib": decoded_mib,
                "throughput_mibs": decoded_mib / mean if decoded_mib is not None and mean > 0 else None,
                "compressed_throughput_mibs": compressed_mib / mean if compressed_mib is not None and mean > 0 else None,
                "command": command,
            }
        )
    return pd.DataFrame(rows)


def load_section(results_dir: Path, manifest: dict, key: str) -> pd.DataFrame | None:
    relative = (manifest.get("sections") or {}).get(key)
    if not relative:
        return None
    directory = results_dir / relative
    if not directory.is_dir():
        return None
    metadata = load_metadata(results_dir, manifest)
    frames = [
        load_zebrac_json(path, metadata)
        for path in sorted(directory.glob("*.json"))
    ]
    frames = [frame for frame in frames if not frame.empty]
    return pd.concat(frames, ignore_index=True) if frames else None


def tools_in_run(frame: pd.DataFrame | None, order: list[str]) -> list[str]:
    if frame is None or frame.empty:
        return []
    present = set(frame["tool"].unique())
    return [tool for tool in order if tool in present]


def build_ratio_comparisons(
    frame: pd.DataFrame,
    value_col: str,
    *,
    peers: list[str],
) -> pd.DataFrame:
    rows: list[dict] = []
    for workload in WORKLOAD_ORDER:
        base = frame[(frame["workload"] == workload) & (frame["tool"] == BASELINE)]
        if base.empty:
            continue
        base_value = float(base[value_col].iloc[0])
        for tool in peers:
            hit = frame[(frame["workload"] == workload) & (frame["tool"] == tool)]
            if hit.empty:
                continue
            peer_value = float(hit[value_col].iloc[0])
            rows.append(
                {
                    "workload": workload,
                    "tool": tool,
                    "competitor": display_tool(tool),
                    "reference": base_value,
                    "peer": peer_value,
                    "ratio": peer_value / base_value if base_value > 0 else None,
                }
            )
    return pd.DataFrame(rows)


def md_pivot(frame: pd.DataFrame, tools: list[str], fmt) -> str:
    work = frame[frame["tool"].isin(tools)].copy()
    if work.empty:
        return "_No data._"
    work["cell"] = work.apply(fmt, axis=1)
    pivot = work.pivot(index="workload", columns="tool", values="cell")
    pivot = pivot[[tool for tool in tools if tool in pivot.columns]]
    pivot = pivot.reindex([key for key in WORKLOAD_ORDER if key in pivot.index])
    pivot = pivot.fillna("")
    pivot.index = [display_workload(key) for key in pivot.index]
    pivot.index.name = "Workload"
    pivot = pivot.rename(columns={tool: display_tool(tool) for tool in pivot.columns})
    return to_markdown_aligned(pivot)


def md_ratio_table(frame: pd.DataFrame, tools: list[str], fmt_reference, fmt_peer) -> str:
    comparisons = build_ratio_comparisons(frame, "mean", peers=[tool for tool in tools if tool != BASELINE])
    if comparisons.empty:
        return "_No comparable peer lanes._"
    rows = []
    for row in comparisons.itertuples(index=False):
        rows.append(
            {
                "Workload": display_workload(row.workload),
                "z-fastq vs": row.competitor,
                "z-fastq": fmt_reference(row),
                "Peer": fmt_peer(row),
                "Time relative to z-fastq": format_ratio(row.ratio),
            }
        )
    return to_markdown_aligned(pd.DataFrame(rows), index=False)


def _bar_patches(tools: list[str]) -> list:
    patches = []
    for tool in tools:
        descriptive = tool in DESCRIPTIVE_TOOLS
        patches.append(
            mpatches.Patch(
                facecolor=COLORS.get(tool, "#888888"),
                edgecolor=COLORS.get(tool, "#888888") if descriptive else "none",
                hatch="///" if descriptive else None,
                alpha=0.78 if descriptive else 0.88,
                label=display_tool(tool),
            )
        )
    return patches


def _draw_metric_facet(ax, frame, tools, workloads, *, value_col, ylabel, value_floor, log_y) -> None:
    std_col = {
        "mean": "stddev",
        "peak_rss_mb": "peak_rss_stddev_mb",
        "minor_faults": "minor_faults_stddev",
    }[value_col]
    width = min(0.12, 0.80 / max(1, len(tools)))
    x = list(range(len(workloads)))
    tops: dict[tuple[str, str], tuple[float, float]] = {}
    for tool_index, tool in enumerate(tools):
        for workload_index, workload in enumerate(workloads):
            row = frame[(frame["workload"] == workload) & (frame["tool"] == tool)]
            if row.empty:
                continue
            value = max(float(row[value_col].iloc[0]), value_floor)
            stddev = max(float(row[std_col].iloc[0]), 0.0)
            xpos = workload_index + (tool_index - len(tools) / 2 + 0.5) * width
            kwargs: dict = {
                "width": width,
                "color": COLORS.get(tool, "#888888"),
                "alpha": 0.78 if tool in DESCRIPTIVE_TOOLS else 0.88,
                "zorder": 2,
            }
            if tool in DESCRIPTIVE_TOOLS:
                kwargs.update({"hatch": "///", "edgecolor": COLORS.get(tool, "#888888"), "linewidth": 0.5})
            if stddev > 0:
                kwargs.update(
                    {
                        "yerr": stddev,
                        "capsize": 2,
                        "error_kw": {"elinewidth": 0.8, "capthick": 0.8, "ecolor": "#333333"},
                    }
                )
            ax.bar(xpos, value, **kwargs)
            tops[(workload, tool)] = (xpos, value)

    comparisons = build_ratio_comparisons(
        frame,
        value_col,
        peers=[tool for tool in tools if tool != BASELINE],
    )
    for row in comparisons.itertuples(index=False):
        if value_col != "mean":
            continue
        key = (row.workload, row.tool)
        if key not in tops:
            continue
        xpos, value = tops[key]
        ax.annotate(
            format_ratio(row.ratio),
            (xpos, value),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            rotation=90,
            fontsize=7,
            color="#111111",
            bbox=dict(
                boxstyle="round,pad=0.16",
                facecolor="white",
                edgecolor=COLORS.get(row.tool, "#666666"),
                linewidth=0.8,
                alpha=0.95,
            ),
            zorder=5,
        )

    if log_y:
        ax.set_yscale("log")
        low, high = ax.get_ylim()
        if high > low > 0:
            ax.set_ylim(low, high * 2.0)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [display_workload(workload) for workload in workloads],
        rotation=48,
        ha="right",
        fontsize=8,
    )
    ax.set_ylabel(ylabel, fontsize=10)
    ax.grid(axis="y", alpha=0.28, which="both")
    ax.set_axisbelow(True)


def fig_metric_facets(frame: pd.DataFrame, out: Path, tools: list[str], title: str, note: str) -> Path:
    filtered = frame[frame["tool"].isin(tools)].copy()
    workloads = [key for key in WORKLOAD_ORDER if key in set(filtered["workload"])]
    if filtered.empty or not workloads:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No performance rows", ha="center", va="center")
        ax.set_axis_off()
        return _save(fig, out)
    facets = [
        ("mean", "Wall Time (s)", True, 1e-6),
        ("peak_rss_mb", "Peak RSS (MB)", False, 0.1),
        ("minor_faults", "Minor Page Faults", True, 1.0),
    ]
    width = max(12.0 if len(workloads) <= 3 else 17.0, 2.4 * len(workloads) + 6.0)
    fig, axes = plt.subplots(1, 3, figsize=(width, 7.6), squeeze=False)
    for ax, (value_col, ylabel, log_y, floor) in zip(axes[0], facets):
        _draw_metric_facet(
            ax,
            filtered,
            tools,
            workloads,
            value_col=value_col,
            ylabel=ylabel,
            value_floor=floor,
            log_y=log_y,
        )
    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.98)
    fig.text(0.5, 0.935, note, ha="center", va="top", fontsize=9, color="#444444", style="italic")
    patches = _bar_patches(tools)
    fig.legend(
        handles=patches,
        labels=[patch.get_label() for patch in patches],
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.055),
        ncol=min(len(patches), 6),
        frameon=False,
        columnspacing=1.0,
        handletextpad=0.4,
    )
    fig.subplots_adjust(left=0.045, right=0.99, bottom=0.28, top=0.86, wspace=0.28)
    return _save(fig, out)


def occupancy_frame(frame: pd.DataFrame, tools: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    scored = occupancy_tools(tools)
    for workload in WORKLOAD_ORDER:
        base = frame[(frame["workload"] == workload) & (frame["tool"] == BASELINE)]
        if base.empty:
            continue
        z_wall = float(base["mean"].iloc[0])
        z_rss = float(base["peak_rss_mb"].iloc[0])
        if z_wall <= 0 or z_rss <= 0:
            continue
        z_occupancy = z_wall * z_rss
        for tool in scored:
            hit = frame[(frame["workload"] == workload) & (frame["tool"] == tool)]
            if hit.empty:
                continue
            wall = float(hit["mean"].iloc[0])
            rss = float(hit["peak_rss_mb"].iloc[0])
            occupancy = wall * rss
            rows.append(
                {
                    "workload": workload,
                    "tool": tool,
                    "wall": wall,
                    "rss": rss,
                    "occupancy": occupancy,
                    "wall_x": wall / z_wall,
                    "rss_x": rss / z_rss,
                    "cost_x": occupancy / z_occupancy,
                    "efficiency": z_occupancy / occupancy,
                }
            )
    return pd.DataFrame(rows)


def _plot_isocost_xy(
    ax,
    xs: list[float],
    ys: list[float],
    *,
    color: str,
    linestyle,
    linewidth: float,
    alpha: float,
    zorder: float,
) -> None:
    if len(xs) < 2:
        return
    ax.plot(
        xs,
        ys,
        color=color,
        linestyle=linestyle,
        linewidth=linewidth,
        alpha=alpha,
        zorder=zorder,
        solid_capstyle="round",
        dash_capstyle="round",
        solid_joinstyle="round",
        clip_on=True,
    )


def fig_occupancy_scatter(
    frame: pd.DataFrame,
    out: Path,
    title: str,
    tools: list[str],
) -> Path:
    workloads = [key for key in WORKLOAD_ORDER if key in set(frame["workload"])]
    if frame.empty or not workloads:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No occupancy rows", ha="center", va="center")
        ax.set_axis_off()
        return _save(fig, out)

    walls = [float(value) for value in frame["wall"]]
    rsses = [float(value) for value in frame["rss"]]
    t_min = max(min(walls) / 1.22, 1e-4)
    t_max = max(walls) * 1.32
    z_rsses = [float(value) for value in frame.loc[frame["tool"] == BASELINE, "rss"]]
    if not z_rsses:
        raise SystemExit("error: occupancy scatter needs z-fastq ISA-L")
    y_max = max(max(rsses) * 1.16, max(z_rsses) * 4.05)
    label_stroke = [pe.withStroke(linewidth=3.2, foreground="white", alpha=0.92)]

    fig, ax = plt.subplots(figsize=(10.2, 7.6))
    steps = 160
    isocosts = (
        (1.0, 2.25, 0.96, 1.55),
        (2.0, 1.35, 0.72, 1.25),
        (4.0, 1.12, 0.55, 1.15),
        (8.0, 0.95, 0.42, 1.05),
    )
    for workload in workloads:
        subset = frame[frame["workload"] == workload]
        base = subset[subset["tool"] == BASELINE]
        if base.empty:
            continue
        z_wall = float(base["wall"].iloc[0])
        occupancy = z_wall * float(base["rss"].iloc[0])
        marker, linestyle, _ = workload_style(workload)
        t0 = max(float(subset["wall"].min()) / 1.16, 1e-4)
        t1 = float(subset["wall"].max()) * 1.18
        span = math.log(t1 / t0)
        grid = [t0 * math.exp(span * index / (steps - 1)) for index in range(steps)]
        for multiplier, linewidth, shade, zorder in isocosts:
            color = COLORS[BASELINE] if multiplier == 1 else "#5a5a5a"
            xs: list[float] = []
            ys: list[float] = []
            for wall in grid:
                rss = occupancy * multiplier / wall
                if not math.isfinite(rss) or rss <= 0 or rss > y_max * 1.01:
                    _plot_isocost_xy(
                        ax,
                        xs,
                        ys,
                        color=color,
                        linestyle=linestyle,
                        linewidth=linewidth,
                        alpha=shade,
                        zorder=zorder,
                    )
                    xs = []
                    ys = []
                    continue
                xs.append(wall)
                ys.append(rss)
            _plot_isocost_xy(
                ax,
                xs,
                ys,
                color=color,
                linestyle=linestyle,
                linewidth=linewidth,
                alpha=shade,
                zorder=zorder,
            )
            label_wall = min(max(z_wall * 1.18, t0 * 1.05), t1)
            label_rss = occupancy * multiplier / label_wall
            if 0.12 * y_max < label_rss < 0.96 * y_max:
                ax.text(
                    label_wall,
                    label_rss,
                    f"{multiplier:.0f}x",
                    fontsize=7.5,
                    color="#3f3f3f" if multiplier > 1 else "#7a4a10",
                    ha="left",
                    va="bottom",
                    zorder=2.4,
                    path_effects=label_stroke,
                )

    for tool in tools:
        color = COLORS.get(tool, "#888888")
        is_z = tool == BASELINE
        is_native = tool == "z-fastq-native"
        for workload in workloads:
            hit = frame[(frame["tool"] == tool) & (frame["workload"] == workload)]
            if hit.empty:
                continue
            marker, _, _ = workload_style(workload)
            ax.scatter(
                float(hit["wall"].iloc[0]),
                float(hit["rss"].iloc[0]),
                s=118 if is_z else (92 if is_native else 86),
                marker=marker,
                facecolor=color,
                edgecolor="#111111" if is_z else "white",
                linewidths=1.25 if is_z else 0.85,
                zorder=5 if is_z else (4.4 if is_native else 4),
                label="_nolegend_",
            )

    ax.set_xscale("log")
    ax.set_xlim(t_min, t_max)
    ax.set_ylim(0.0, y_max)
    ax.set_xlabel("Mean wall time (s), log")
    ax.set_ylabel("Peak RSS (MB)")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.grid(True, which="major", alpha=0.32, linewidth=0.7)
    ax.grid(True, which="minor", alpha=0.14, linewidth=0.45)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("#4a4a4a")
        spine.set_linewidth(0.85)
    ax.tick_params(colors="#333333", length=4, width=0.8)

    tool_handles = [
        mlines.Line2D(
            [],
            [],
            color=COLORS.get(tool, "#888888"),
            marker="o",
            linestyle="None",
            markersize=8.5,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=display_tool(tool),
        )
        for tool in tools
    ]
    workload_handles = [
        mlines.Line2D(
            [],
            [],
            color="#333333",
            marker=workload_style(workload)[0],
            linestyle=workload_style(workload)[1],
            linewidth=1.7,
            markersize=8.5,
            markerfacecolor="#333333",
            markeredgecolor="white",
            markeredgewidth=0.6,
            label=display_workload(workload),
        )
        for workload in workloads
    ]
    tool_legend = fig.legend(
        handles=tool_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.095),
        ncol=min(len(tool_handles), 4),
        fontsize=8,
        frameon=False,
        columnspacing=1.15,
        handletextpad=0.4,
    )
    fig.add_artist(tool_legend)
    fig.legend(
        handles=workload_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.035),
        ncol=min(len(workload_handles), 3),
        fontsize=7.5,
        frameon=False,
        columnspacing=1.0,
        handletextpad=0.4,
    )
    fig.text(
        0.5,
        0.968,
        "Isocosts: RSS = k × (z-fastq wall × RSS) / wall. Dash matches workload shape. Gold 1x, gray 2x/4x/8x. Lower-left is better.",
        ha="center",
        va="top",
        fontsize=8,
        color="#444444",
        style="italic",
    )
    fig.subplots_adjust(left=0.11, right=0.985, bottom=0.31, top=0.885)
    return _save(fig, out, dpi=180)


def fig_efficiency_bars(
    frame: pd.DataFrame,
    out: Path,
    title: str,
    tools: list[str],
) -> Path:
    workloads = [key for key in WORKLOAD_ORDER if key in set(frame["workload"])]
    if frame.empty or not workloads:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No occupancy peers", ha="center", va="center")
        ax.set_axis_off()
        return _save(fig, out)

    width = min(0.18, 0.72 / max(1, len(workloads)))
    x = list(range(len(tools)))
    fig, ax = plt.subplots(figsize=(13.0, 5.9))
    ax.axhline(1.0, color=COLORS[BASELINE], linestyle=(0, (4, 3)), linewidth=1.2, alpha=0.7, zorder=1)
    for workload_index, workload in enumerate(workloads):
        hatch = workload_style(workload)[2]
        positions = [
            index + (workload_index - len(workloads) / 2 + 0.5) * width
            for index in x
        ]
        for tool, position in zip(tools, positions):
            hit = frame[(frame["tool"] == tool) & (frame["workload"] == workload)]
            if hit.empty:
                continue
            value = float(hit["efficiency"].iloc[0])
            ax.bar(
                position,
                value,
                width=width,
                color=COLORS.get(tool, "#888888"),
                alpha=0.88,
                hatch=hatch,
                edgecolor="#333333",
                linewidth=0.4,
                zorder=2,
            )
            if value > 0:
                ax.annotate(
                    f"{value:.2f}x",
                    (position, value),
                    xytext=(0, 4),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=6.5,
                    rotation=90,
                )
    ax.set_xticks(x)
    ax.set_xticklabels([display_tool(tool) for tool in tools], fontsize=9)
    ax.set_ylabel("Efficiency vs z-fastq ISA-L (occupancy)", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.28)
    ax.set_axisbelow(True)
    workload_handles = [
        mpatches.Patch(
            facecolor="#888888",
            hatch=workload_style(workload)[2],
            edgecolor="#333333",
            label=display_workload(workload),
        )
        for workload in workloads
    ]
    fig.legend(
        handles=workload_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=min(len(workload_handles), 3),
        fontsize=7.5,
        frameon=False,
        columnspacing=1.0,
    )
    fig.text(
        0.5,
        0.955,
        "Higher is better. 1.00 = the same wall × RSS as z-fastq ISA-L on that workload. Color is the tool; hatch is the workload.",
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.27, top=0.86)
    return _save(fig, out)


def md_occupancy_tables(
    frame: pd.DataFrame,
    tools: list[str],
    nums: ReportCounters,
) -> str:
    scored = occupancy_tools(tools)
    if frame.empty:
        return "_No occupancy rows._"
    work = frame[frame["tool"].isin(scored)].copy()
    if work.empty:
        return "_No occupancy rows._"
    occupancy = work.pivot(index="workload", columns="tool", values="occupancy")
    occupancy = occupancy.reindex([key for key in WORKLOAD_ORDER if key in occupancy.index])
    occupancy = occupancy[[tool for tool in scored if tool in occupancy.columns]]
    occupancy = occupancy.map(lambda value: f"{value:.4f}" if pd.notna(value) else "")
    occupancy = occupancy.rename(columns={tool: display_tool(tool) for tool in occupancy.columns})
    occupancy.index = [display_workload(key) for key in occupancy.index]
    occupancy.index.name = "Workload"

    efficiency = work.pivot(index="workload", columns="tool", values="efficiency")
    efficiency = efficiency.reindex([key for key in WORKLOAD_ORDER if key in efficiency.index])
    efficiency = efficiency[[tool for tool in scored if tool in efficiency.columns]]

    def format_efficiency(value) -> str:
        if pd.isna(value):
            return ""
        return "1x" if abs(float(value) - 1.0) < 1e-9 else format_ratio(float(value))

    efficiency = efficiency.map(format_efficiency)
    efficiency = efficiency.rename(columns={tool: display_tool(tool) for tool in efficiency.columns})
    efficiency.index = [display_workload(key) for key in efficiency.index]
    efficiency.index.name = "Workload"
    occupancy_number = nums.next_table()
    efficiency_number = nums.next_table()
    return join(
        [
            f"**Table {occupancy_number}:** Occupancy (MB·s) = mean wall × peak RSS.",
            "",
            to_markdown_aligned(occupancy),
            "",
            f"**Table {efficiency_number}:** Occupancy efficiency vs z-fastq ISA-L. Higher is better.",
            "",
            to_markdown_aligned(efficiency),
        ]
    )


def md_figure_block(nums: ReportCounters, path: str, description: str, bullets: list[str]) -> str:
    number = nums.next_figure()
    return join(
        [
            f"**Figure {number}:** {description}",
            "",
            f"![Figure {number}]({path})",
            "",
            f"**Reading Figure {number}**",
            "",
            "\n".join(f"- {bullet}" for bullet in bullets),
            "",
            '<div style="margin: 1.5em 0"></div>',
            "",
        ]
    )


def md_metric_block(
    frame: pd.DataFrame,
    tools: list[str],
    *,
    heading: str,
    intro: str,
    value_col: str,
    formatter,
    nums: ReportCounters,
    include_time_ratios: bool = False,
) -> str:
    table_number = nums.next_table()
    blocks = [
        f"### {heading}",
        "",
        intro,
        "",
        f"**Table {table_number}:** Mean ± standard deviation by workload and tool.",
        "",
        md_pivot(frame, tools, formatter),
        "",
    ]
    if include_time_ratios:
        ratio_number = nums.next_table()
        blocks.extend(
            [
                f"<details><summary><strong>Table {ratio_number}:</strong> Time relative to z-fastq ISA-L.</summary>",
                "",
                md_ratio_table(
                    frame,
                    tools,
                    lambda row: f"{row.reference:.4f}s",
                    lambda row: f"{row.peer:.4f}s",
                ),
                "",
                "</details>",
                "",
            ]
        )
    return join(blocks)


def md_throughput(frame: pd.DataFrame, tools: list[str], nums: ReportCounters) -> str:
    decoded_number = nums.next_table()
    compressed_number = nums.next_table()
    return join(
        [
            "### Decoded throughput",
            "",
            "Decoded FASTQ MiB/s is the uncompressed input size divided by wall time. That is the gzip unit that matches inflate work; compressed on-disk throughput is supporting context.",
            "",
            f"**Table {decoded_number}:** Decoded FASTQ MiB/s.",
            "",
            md_pivot(frame, tools, fmt_decoded_mibs),
            "",
            f"<details><summary><strong>Table {compressed_number}:</strong> Compressed on-disk MiB/s.</summary>",
            "",
            md_pivot(frame, tools, fmt_compressed_mibs),
            "",
            "</details>",
            "",
        ]
    )


def lane_role(workload: str) -> str:
    return WORKLOAD_ROLES.get(workload, workload.split("__", 1)[0])


def summarize_lane_reason(reason: str) -> str:
    if "accepted" in reason:
        return "timed"
    if "probed, not timed" in reason:
        return "probed, not timed"
    if "not supported" in reason:
        return "unsupported"
    if reason.startswith("exit "):
        return "rejected"
    if "unavailable" in reason:
        return "unavailable"
    return reason


def md_role_matrix(manifest: dict) -> str:
    reasons = manifest.get("lane_reasons") or {}
    by_tool: dict[str, dict[str, dict[str, str]]] = {}
    for key, reason in reasons.items():
        parts = key.split("|", 2)
        if len(parts) != 3:
            continue
        section, workload, tool = parts
        role = lane_role(workload)
        if tool == "irma" and role.startswith("interleaved"):
            continue
        compression = "gzip" if "gzip" in section else "plain"
        by_tool.setdefault(tool, {}).setdefault(role, {})[compression] = summarize_lane_reason(str(reason))

    def cell(tool: str, role: str) -> str:
        if tool == "z-fastq":
            return "timed"
        if tool == "z-fastq-native":
            return "timed (gzip)"
        if tool == "irma" and role.startswith("interleaved"):
            return "unsupported"
        statuses = by_tool.get(tool, {}).get(role, {})
        plain = statuses.get("plain")
        gzip = statuses.get("gzip")
        values = {status for status in (plain, gzip) if status}
        if not values:
            return "unsupported" if tool == "irma" and role.startswith("interleaved") else ""
        if values == {"timed"}:
            return "timed"
        if plain == "timed" and gzip == "probed, not timed":
            return "timed (plain)"
        if gzip == "probed, not timed" and plain in (None, "unsupported", "unavailable"):
            return "probed, not timed"
        if "timed" in values and "probed, not timed" in values:
            return "timed (plain)"
        if "timed" in values:
            return "timed"
        if "probed, not timed" in values:
            return "probed, not timed"
        if "rejected" in values:
            return "rejected"
        if "unavailable" in values:
            return "unavailable"
        if "unsupported" in values:
            return "unsupported"
        return next(iter(values))

    columns = ["Tool"] + [ROLE_MATRIX_LABELS[role] for role in ROLE_MATRIX_ORDER]
    rows = []
    for tool in ["z-fastq", "z-fastq-native", *PLAIN_TOOLS[1:]]:
        row = {"Tool": display_tool(tool)}
        for role in ROLE_MATRIX_ORDER:
            row[ROLE_MATRIX_LABELS[role]] = cell(tool, role)
        rows.append(row)
    return to_markdown_aligned(pd.DataFrame(rows, columns=columns), index=False)


def md_lane_fold(manifest: dict) -> str:
    reasons = manifest.get("lane_reasons") or {}
    rows = []
    for key, reason in sorted(reasons.items()):
        section, workload, tool = key.split("|", 2)
        result = summarize_lane_reason(str(reason))
        if tool == "irma" and lane_role(workload).startswith("interleaved"):
            result = "unsupported"
        rows.append(
            {
                "Section": section.removeprefix("perf_"),
                "Workload": display_workload(workload),
                "Tool": display_tool(tool),
                "Result": result,
            }
        )
    if not rows:
        return "_No peer probes recorded._"
    return to_markdown_aligned(pd.DataFrame(rows), index=False)


def md_fixture_table(manifest: dict) -> str:
    rows_in = manifest.get("peer_fixtures") or []
    if not rows_in:
        return "_No peer fixture probes were recorded in this run. Re-run `bench/sample/run.sh` to fill this matrix._"
    by_key = {(str(row.get("fixture")), str(row.get("tool"))): str(row.get("result", "")) for row in rows_in}
    tools = list(PLAIN_TOOLS[1:])
    table_rows = []
    for fixture in FIXTURE_ORDER:
        row = {"Fixture": FIXTURE_LABELS.get(fixture, fixture)}
        for tool in tools:
            row[display_tool(tool)] = by_key.get((fixture, tool), "")
        table_rows.append(row)
    return to_markdown_aligned(pd.DataFrame(table_rows), index=False)


def md_capability(manifest: dict) -> str:
    peers = pd.DataFrame(
        [
            {
                "Tool": "z-fastq (ISA-L)",
                "Command": "`z-fastq sample --fraction P|--count K --seed 11`",
                "Same sample job": "yes; source-ordered subset, pair as unit",
                "Timed as": "contract reference",
            },
            {
                "Tool": "z-fastq (native)",
                "Command": "`z-fastq-native sample ...`",
                "Same sample job": "yes, gzip only",
                "Timed as": "second inflate backend",
            },
            {
                "Tool": "seqtk",
                "Command": "`seqtk sample -s 11` and `seqtk sample -2 -s 11` for exact count",
                "Same sample job": "SE only, on screened LF/nonempty/bare-plus inputs",
                "Timed as": "exact SE reference",
            },
            {
                "Tool": "SeqKit",
                "Command": "`seqkit sample -j 1 -p P -s 11` / `seqkit sample2 -j 1 -2 -n K -s 11`",
                "Same sample job": "no; different RNG; SE only",
                "Timed as": "hatched; descriptive",
            },
            {
                "Tool": "Rasusa",
                "Command": "`rasusa reads -p P|-n K -s 11` (SE stdout; two `-o /dev/null` when paired)",
                "Same sample job": "no; `-p` is one-pass probability; paired writes split files",
                "Timed as": "hatched; descriptive",
            },
            {
                "Tool": "fq",
                "Command": "`fq subsample -p P|-n K -s 11 --r1-dst /dev/null` (and `--r2-dst` when paired)",
                "Same sample job": "no; split outputs, different RNG; no interleaved input",
                "Timed as": "hatched; paired behavior reference",
            },
            {
                "Tool": "fqkit",
                "Command": "`fqkit subfq -@ 1 -2 -n K -s 11`",
                "Same sample job": "no; exact SE only, no fraction or pairs",
                "Timed as": "hatched; descriptive",
            },
            {
                "Tool": "IRMA Core",
                "Command": "`irma-core sampler -s 11 -p 10|-t K`",
                "Same sample job": "no; integer percent; gzip and plain select different records; one interleaved file is SE records",
                "Timed as": "hatched; plain SE and paired only",
            },
            {
                "Tool": "BBTools",
                "Command": "`reformat.sh samplerate=P|samplereadstarget=K sampleseed=11 out=/dev/null`",
                "Same sample job": "no; JVM; can rewrite qualities",
                "Timed as": "hatched; RSS not occupancy",
            },
        ]
    )
    omitted = pd.DataFrame(
        [
            {
                "Tool": "IRMA Core on interleaved",
                "Why it is not here": "one-file sampler treats interleaved FASTQ as SE records, not pair units",
            },
            {
                "Tool": "Fasten `fasten_sample`",
                "Why it is not here": "stdin-only, no seed, no gzip; zebrac cannot feed it without a shell",
            },
            {
                "Tool": "SeqFu",
                "Why it is not here": "no sample command",
            },
            {
                "Tool": "seqtk on pairs",
                "Why it is not here": "record-independent, not a pair unit",
            },
            {
                "Tool": "Needletail / Helicase adapters",
                "Why it is not here": "count and stats wrappers; no sample command",
            },
            {
                "Tool": "fastp",
                "Why it is not here": "QC and transformation, not a sampler",
            },
        ]
    )
    return join(
        [
            "A valid z-fastq sample writes a source-ordered FASTQ subset to stdout and exits 0. `--fraction` is one pass and approximate; `--count` is two passes over a regular file and emits `min(K, N)` units. `--paired` and `--interleaved` sample the pair, never a lone mate, and emit interleaved stdout. seqtk is the screened SE exact-reference peer. Other samplers are descriptive: they are timed only when they exit 0 on a positive workload they claim to support as the same sampling unit.",
            "",
            to_markdown_aligned(peers, index=False),
            "",
            "**Mode coverage**",
            "",
            md_role_matrix(manifest),
            "",
            "**Not timed**",
            "",
            to_markdown_aligned(omitted, index=False),
            "",
            "<details><summary><strong>Peer lane decisions</strong> (plain and gzip, every workload)</summary>",
            "",
            md_lane_fold(manifest),
            "",
            "</details>",
        ]
    )


def md_overview(manifest: dict) -> str:
    set_name = manifest.get("real_set", "publication")
    fraction = manifest.get("sample_fraction", "0.1")
    count = manifest.get("sample_count", 1000)
    seed = manifest.get("sample_seed", 11)
    se_frac = overview_shapes(FAMILY_ROLES["se_fraction"])
    se_count = overview_shapes(FAMILY_ROLES["se_count"])
    pair_frac = overview_shapes(FAMILY_ROLES["pairs_fraction"])
    pair_count = overview_shapes(FAMILY_ROLES["pairs_count"])
    return join(
        [
            "This report times `z-fastq sample` on the shared FASTQ corpus. Coverage is kept in one suite; the showcase is four families so SE fraction, SE exact-count, pair-layout fraction, and pair-layout count are not one race.",
            "",
            f"Set: **{set_name}**. Parameters: **fraction {fraction}**, **count {count}**, **seed {seed}**.",
            "",
            f"- **SE fraction:** {se_frac}. One-pass Bernoulli. Occupancy lives here.",
            f"- **SE count:** {se_count}. Two-pass reservoir, `min(K, N)` records. Occupancy lives here separately; it is not mixed with fraction.",
            f"- **Pair layout, fraction:** {pair_frac}. `--paired` vs `--interleaved` on the matching mates. SE is half the records and is not on this plot.",
            f"- **Pair layout, count:** {pair_count}. Exact `min(K, M)` complete pairs. Tables and 3-facet only.",
            "",
            "**What is timed**",
            "",
            "- One process, one sampling workload. Zebrac starts the argv directly; no shell or pipeline is measured. Stdout is discarded by zebrac; file-writing peers sink to `/dev/null`.",
            "- Plain and gzip inputs are separate sections. Gzip includes native z-fastq only as the second inflate backend.",
            "- seqtk is timed as the screened SE exact reference. Other samplers are descriptive peers and are hatched.",
            "- IRMA Core gzip is probed and not timed: gzip and plain select different records. Fasten is not timed: stdin-only, no seed. IRMA interleaved is not timed: one file is SE records, not pair units.",
        ]
    )


def md_correctness(manifest: dict) -> str:
    status = manifest.get("verify_pass") or "not recorded"
    log = manifest.get("verify_log") or ""
    if manifest.get("verify_skipped"):
        body = "The sample contract preflight was skipped. This run is not publishable."
    else:
        body = (
            f"Status: **{status}**. Both z-fastq binaries passed fraction 0/1 and count 0/1 fixtures, stdin fraction, "
            "screened seqtk byte identity (fraction 0.5 and count 3, seed 11), paired vs interleaved identical "
            "pair indexes on fixtures and retained real outputs, even pair-complete output, stdin exact-count rejection, and `--json` rejection. "
            "Exact-count real files were required to emit `min(K, N)` records (or `2 × min(K, M)` pairs). "
            "ISA-L vs native gzip outputs were byte-compared when retained. "
            "seqtk byte identity on real SE files is required except Dense (named plus). Positive peer lanes "
            "were probed next; incompatible peer behavior was recorded rather than treated as a z-fastq failure."
        )
    return join(
        [
            body,
            "",
            f"Log: `{log}`." if log else "No verification log was recorded.",
            "",
            "**Peer fixture probes**",
            "",
            "These rows are descriptive. A fail does not stop the run. `pass` is exit 0, and for stdout writers a matching record count on count fixtures or even pair output on pair-fraction fixtures; `fail` is a nonzero exit or a cardinality mismatch on a fixture z-fastq accepts; `unsupported` means that tool is not invoked for that mode. File-writing peers (`fq`, Rasusa paired, BBTools) are exit 0 only.",
            "",
            md_fixture_table(manifest),
        ]
    )


def md_provenance(manifest: dict) -> str:
    tools = manifest.get("tools") or {}
    lines = [
        f"- **{name}:** {tools[name]}"
        for name in ("seqtk", "seqkit", "rasusa", "fq", "fqkit", "irma", "bbtools")
        if tools.get(name)
    ]
    isa = f"- **z-fastq ISA-L:** {manifest.get('z_fastq', '')}"
    native = f"- **z-fastq native:** {manifest.get('z_fastq_native', '')}"
    if manifest.get("z_fastq_bytes"):
        isa += f" ({int(manifest['z_fastq_bytes']):,} bytes)"
    if manifest.get("z_fastq_native_bytes"):
        native += f" ({int(manifest['z_fastq_native_bytes']):,} bytes)"
    return join(
        [
            f"- **Timestamp:** `{manifest.get('timestamp', '')}`",
            f"- **Runner:** zebrac, warm page cache, runs={manifest.get('runs')}, warmup={manifest.get('warmup')}, duration_ms={manifest.get('duration_ms')}",
            f"- **Parameters:** fraction={manifest.get('sample_fraction')}, count={manifest.get('sample_count')}, seed={manifest.get('sample_seed')}",
            f"- **zebrac:** {manifest.get('zebrac', '')}",
            isa,
            native,
            "",
            "**Peer versions**",
            "",
            "\n".join(lines) if lines else "- none recorded",
        ]
    )


def md_perf_section(
    frame: pd.DataFrame,
    *,
    title: str,
    intro: str,
    tools: list[str],
    fig_name: str,
    fig_title: str,
    fig_note: str,
    nums: ReportCounters,
    figures_dir: Path,
    roles: tuple[str, ...],
    include_throughput: bool = False,
    include_occupancy: bool = False,
    figure_caption: str | None = None,
    figure_bullets: list[str] | None = None,
) -> str:
    global WORKLOAD_ORDER
    keys = family_workloads(roles)
    if not keys:
        return ""
    saved_order = WORKLOAD_ORDER
    WORKLOAD_ORDER = keys
    try:
        work = drop_disallowed_lanes(frame[frame["workload"].isin(set(keys))].copy())
        present = tools_in_run(work, tools)
        if work.empty or not present:
            return ""
        figure = figures_dir / fig_name
        fig_metric_facets(work, figure, present, fig_title, fig_note)
        bullets = figure_bullets or [
            "Facets: wall time (log y) | peak RSS (linear) | minor page faults (log y).",
            "X-axis: workload role and catalog ids; paired inputs are one workload.",
            "Gold bars are z-fastq ISA-L. Hatched bars are descriptive samplers with different RNGs or output layouts.",
        ]
        blocks = [
            f"## {title}",
            "",
            intro,
            "",
            md_metric_block(
                work,
                present,
                heading="Wall time",
                intro="Zebrac mean wall time for one sample process. Vertical annotations are wall time relative to z-fastq ISA-L, including hatched descriptive lanes. Hatched ratios are not a ranking.",
                value_col="mean",
                formatter=fmt_wall,
                nums=nums,
                include_time_ratios=True,
            ),
            md_metric_block(
                work,
                present,
                heading="Peak memory",
                intro="Peak RSS of the timed process. These are raw process measurements; sampler implementations do not share one memory contract.",
                value_col="peak_rss_mb",
                formatter=fmt_rss,
                nums=nums,
            ),
            md_metric_block(
                work,
                present,
                heading="Minor page faults",
                intro="Minor page faults for the same direct commands. These are descriptive process measurements.",
                value_col="minor_faults",
                formatter=fmt_faults,
                nums=nums,
            ),
        ]
        if include_throughput:
            blocks.append(md_throughput(work, present, nums))
        blocks.append(
            md_figure_block(
                nums,
                f"results/figures/{fig_name}",
                figure_caption or "Wall time, peak RSS, and minor page faults for this family.",
                bullets,
            )
        )
        if include_occupancy:
            occupancy = occupancy_frame(work, present)
            scored = occupancy_tools(present)
            if not occupancy.empty and scored:
                stem = fig_name.removesuffix(".png")
                occupancy_name = f"{stem}_occupancy.png"
                efficiency_name = f"{stem}_efficiency.png"
                kind = "Gzip" if include_throughput else "Plain"
                fig_occupancy_scatter(
                    occupancy,
                    figures_dir / occupancy_name,
                    f"{kind} SE sample occupancy: peak RSS vs wall",
                    scored,
                )
                fig_efficiency_bars(
                    occupancy,
                    figures_dir / efficiency_name,
                    f"{kind} SE sample occupancy efficiency vs z-fastq",
                    scored,
                )
                workload_ids = [key for key in WORKLOAD_ORDER if key in set(occupancy["workload"])]
                blocks.extend(
                    [
                        "### Occupancy (50/50 wall and RSS)",
                        "",
                        "Occupancy is mean wall time × peak RSS, and only SE families use it (fraction and count separately). These figures use accepted positive lanes; blank cells mean unsupported or rejected. BBTools stays on the wall/RSS facets and is omitted here because its RSS is a JVM image, not a single sampler process. SeqKit stays on the wall/RSS facets and is omitted here because its RSS is the Go image. Occupancy is descriptive resource accounting, not a semantic-equivalence ranking.",
                        "",
                        md_occupancy_tables(occupancy, scored, nums),
                        "",
                        md_figure_block(
                            nums,
                            f"results/figures/{occupancy_name}",
                            "Peak RSS vs mean wall. Per-workload occupancy isocosts through z-fastq ISA-L.",
                            [
                                "X: mean wall (s, log). Y: peak RSS (MB, linear).",
                                "Color = tool. Shape and line dash = workload; the legend names every SE shape.",
                                "Curves are 1x/2x/4x/8x memory-seconds of that workload's z-fastq ISA-L lane. Gold 1x, gray 2x/4x/8x.",
                                "Lower-left is better. Missing points are not treated as zero.",
                            ],
                        ),
                        md_figure_block(
                            nums,
                            f"results/figures/{efficiency_name}",
                            "Occupancy efficiency vs z-fastq ISA-L. Higher is better. 1.00 matches z-fastq memory-seconds.",
                            [
                                f"X-axis is the tool. Grouped bars are {' / '.join(display_workload(key) for key in workload_ids)}; hatches identify workloads.",
                                "Bar color is the tool. Gold dashed line is 1x.",
                                "This is the inverse of occupancy cost on the scatter and is descriptive for accepted positive SE sample lanes.",
                            ],
                        ),
                    ]
                )
        return join(blocks)
    finally:
        WORKLOAD_ORDER = saved_order


def md_run_banner(manifest: dict) -> str | None:
    try:
        duration = int(manifest.get("duration_ms") or 0)
        runs = int(manifest.get("runs") or 0)
        warmup = int(manifest.get("warmup") or 0)
    except (TypeError, ValueError):
        duration, runs, warmup = 0, 0, 0
    if (duration and duration < 5000) or (runs and runs < 25) or (warmup and warmup < 5):
        return (
            "_This is a bring-up measurement: "
            f"zebrac duration_ms={duration}, runs={runs}, warmup={warmup} "
            "(published default is 5000 ms, 25 runs, 5 warmups). Read it as a suite check, not a published ranking._"
        )
    return None


def squeeze_blank_lines(text: str) -> str:
    while "\n\n\n" in text:
        text = text.replace("\n\n\n", "\n\n")
    return text


def generate(results_dir: Path, allow_incomplete: bool) -> None:
    manifest = load_latest_manifest(results_dir)
    configure_workloads(manifest)
    if not allow_incomplete and manifest.get("verify_skipped"):
        raise SystemExit("error: publishable report requires a passing sample contract preflight")
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    nums = ReportCounters()
    banner = md_run_banner(manifest)
    lines = [
        MARKDOWNLINT_DISABLE,
        "",
        "# z-fastq Sample Benchmark Report",
        "",
        "_Generated by `bench/sample/generate_report.py` from zebrac results._",
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
        md_capability(manifest),
        "",
        "## Sample agreement",
        "",
        md_correctness(manifest),
        "",
        "## Run provenance",
        "",
        md_provenance(manifest),
        "",
    ]

    plain = load_section(results_dir, manifest, "perf_plain")
    gzip = load_section(results_dir, manifest, "perf_gzip")
    required = []
    if not allow_incomplete:
        if plain is None or plain.empty:
            required.append("perf_plain")
        if gzip is None or gzip.empty:
            required.append("perf_gzip")
    if required:
        raise SystemExit("error: missing sections: " + ", ".join(required) + " (pass --allow-incomplete for a draft)")

    sample_count = manifest.get("runs", "?")
    shared_kw = {"nums": nums, "figures_dir": figures_dir}
    hifi_note = ""
    if any("HiFi" in display_workload(key) for key in family_workloads(FAMILY_ROLES["se_count"])):
        hifi_note = " HiFi has 1,114 records, so K=1000 keeps almost every record."
    se_frac_bullets_plain = [
        "Facets: wall time (log y) | peak RSS (linear) | minor page faults (log y).",
        "X-axis: SE catalog shapes. One-pass `--fraction`. Occupancy follows this family.",
        "Gold bars are z-fastq ISA-L. Hatched bars are descriptive samplers. seqtk is the screened exact SE reference and is not hatched.",
    ]
    se_frac_bullets_gzip = [
        "Facets: wall time (log y) | peak RSS (linear) | minor page faults (log y).",
        "X-axis: SE catalog shapes. One-pass `--fraction`. Occupancy follows this family.",
        "Gold bars are z-fastq ISA-L. Amber bars are native inflate. Hatched bars are descriptive samplers. seqtk is the screened exact SE reference and is not hatched.",
    ]
    se_count_bullets_plain = [
        "Facets: wall time (log y) | peak RSS (linear) | minor page faults (log y).",
        "X-axis: SE catalog shapes. Two-pass `--count`. Occupancy follows this family and is not mixed with fraction.",
        "Gold bars are z-fastq ISA-L. Hatched bars are descriptive samplers. seqtk `sample -2` is the screened exact SE reference.",
    ]
    se_count_bullets_gzip = [
        "Facets: wall time (log y) | peak RSS (linear) | minor page faults (log y).",
        "X-axis: SE catalog shapes. Two-pass `--count`. Occupancy follows this family and is not mixed with fraction.",
        "Gold bars are z-fastq ISA-L. Amber bars are native inflate. Hatched bars are descriptive samplers. seqtk `sample -2` is the screened exact SE reference.",
    ]
    pair_frac_bullets_plain = [
        "Facets: wall time (log y) | peak RSS (linear) | minor page faults (log y).",
        "X-axis: `--paired` then `--interleaved` on the matching mates. SE is not on this plot.",
        "Gold bars are z-fastq ISA-L. Hatched bars are descriptive samplers. Blank cells are unsupported for that layout.",
    ]
    pair_frac_bullets_gzip = [
        "Facets: wall time (log y) | peak RSS (linear) | minor page faults (log y).",
        "X-axis: `--paired` then `--interleaved` on the matching mates. SE is not on this plot.",
        "Gold bars are z-fastq ISA-L. Amber bars are native inflate. Hatched bars are descriptive samplers. Blank cells are unsupported for that layout.",
    ]
    pair_count_bullets_plain = [
        "Facets: wall time (log y) | peak RSS (linear) | minor page faults (log y).",
        "X-axis: `--paired --count` then `--interleaved --count` on the matching mates.",
        "Gold bars are z-fastq ISA-L. Hatched bars are descriptive samplers. Blank cells are unsupported for that layout.",
    ]
    pair_count_bullets_gzip = [
        "Facets: wall time (log y) | peak RSS (linear) | minor page faults (log y).",
        "X-axis: `--paired --count` then `--interleaved --count` on the matching mates.",
        "Gold bars are z-fastq ISA-L. Amber bars are native inflate. Hatched bars are descriptive samplers. Blank cells are unsupported for that layout.",
    ]
    pair_blank = " Blank cells mean that tool is not invoked for that layout."

    if plain is not None and not plain.empty:
        lines.append(
            md_perf_section(
                plain,
                title="Performance: SE fraction",
                intro="Single-end `--fraction` on catalog SE shapes. Small set: DenseSmall, Variable, HiFi. Publication set: Dense, Variable, Long. Occupancy lives only here among fraction workloads. seqtk is the screened exact reference; other samplers are descriptive.",
                tools=PLAIN_TOOLS,
                fig_name="perf_plain_se_fraction.png",
                fig_title="Plain SE fraction sample: wall, RSS, page faults",
                fig_note=f"Error bars = zebrac standard deviation (n={sample_count}). Hatched = descriptive samplers.",
                roles=FAMILY_ROLES["se_fraction"],
                include_occupancy=True,
                figure_caption="Wall time, peak RSS, and minor page faults for SE fraction sampling.",
                figure_bullets=se_frac_bullets_plain,
                **shared_kw,
            )
        )
        lines.append("")
        lines.append(
            md_perf_section(
                plain,
                title="Performance: SE count",
                intro=(
                    "Single-end `--count` on the same SE shapes. Two-pass reservoir; output is "
                    f"`min(K, N)` records in source order.{hifi_note} Occupancy is separate from fraction because the jobs are not the same."
                ),
                tools=PLAIN_TOOLS,
                fig_name="perf_plain_se_count.png",
                fig_title="Plain SE exact-count sample: wall, RSS, page faults",
                fig_note=f"Error bars = zebrac standard deviation (n={sample_count}). Hatched = descriptive samplers.",
                roles=FAMILY_ROLES["se_count"],
                include_occupancy=True,
                figure_caption="Wall time, peak RSS, and minor page faults for SE exact-count sampling.",
                figure_bullets=se_count_bullets_plain,
                **shared_kw,
            )
        )
        lines.append("")
    if gzip is not None and not gzip.empty:
        lines.append(
            md_perf_section(
                gzip,
                title="Performance: SE fraction (gzip)",
                intro="The same SE `--fraction` shapes as the plain section. Gold is z-fastq ISA-L; amber is native Zig inflate. IRMA Core is not timed on gzip because gzip and plain select different records. Decoded FASTQ MiB/s uses the uncompressed sibling size divided by wall time.",
                tools=GZIP_TOOLS,
                fig_name="perf_gzip_se_fraction.png",
                fig_title="Gzip SE fraction sample: wall, RSS, page faults",
                fig_note=f"Error bars = zebrac standard deviation (n={sample_count}). Amber = native inflate. Hatched = descriptive samplers.",
                roles=FAMILY_ROLES["se_fraction"],
                include_throughput=True,
                include_occupancy=True,
                figure_caption="Wall time, peak RSS, and minor page faults for gzip SE fraction sampling.",
                figure_bullets=se_frac_bullets_gzip,
                **shared_kw,
            )
        )
        lines.append("")
        lines.append(
            md_perf_section(
                gzip,
                title="Performance: SE count (gzip)",
                intro=(
                    "The same SE `--count` shapes as the plain section. Gold is z-fastq ISA-L; amber is native Zig inflate."
                    f"{hifi_note}"
                ),
                tools=GZIP_TOOLS,
                fig_name="perf_gzip_se_count.png",
                fig_title="Gzip SE exact-count sample: wall, RSS, page faults",
                fig_note=f"Error bars = zebrac standard deviation (n={sample_count}). Amber = native inflate. Hatched = descriptive samplers.",
                roles=FAMILY_ROLES["se_count"],
                include_throughput=True,
                include_occupancy=True,
                figure_caption="Wall time, peak RSS, and minor page faults for gzip SE exact-count sampling.",
                figure_bullets=se_count_bullets_gzip,
                **shared_kw,
            )
        )
        lines.append("")
    if plain is not None and not plain.empty:
        lines.append(
            md_perf_section(
                plain,
                title="Performance: pair layout, fraction",
                intro=(
                    "`--paired --fraction` vs `--interleaved --fraction` on the matching mates (same pair count). "
                    "SE is half the records of this pair and is not plotted here as pair overhead. "
                    f"Tables and a 3-facet only; no occupancy.{pair_blank}"
                ),
                tools=PLAIN_TOOLS,
                fig_name="perf_plain_pairs_fraction.png",
                fig_title="Plain pair-layout fraction sample: wall, RSS, page faults",
                fig_note=f"Error bars = zebrac standard deviation (n={sample_count}). Hatched = descriptive samplers.",
                roles=FAMILY_ROLES["pairs_fraction"],
                figure_caption="Wall time, peak RSS, and minor page faults for paired vs interleaved fraction sampling.",
                figure_bullets=pair_frac_bullets_plain,
                **shared_kw,
            )
        )
        lines.append("")
        lines.append(
            md_perf_section(
                plain,
                title="Performance: pair layout, count",
                intro=(
                    "`--paired --count` vs `--interleaved --count` on the matching mates. "
                    "Output is `min(K, M)` complete pairs as interleaved records. "
                    f"Tables and a 3-facet only; no occupancy.{pair_blank}"
                ),
                tools=PLAIN_TOOLS,
                fig_name="perf_plain_pairs_count.png",
                fig_title="Plain pair-layout exact-count sample: wall, RSS, page faults",
                fig_note=f"Error bars = zebrac standard deviation (n={sample_count}). Hatched = descriptive samplers.",
                roles=FAMILY_ROLES["pairs_count"],
                figure_caption="Wall time, peak RSS, and minor page faults for paired vs interleaved exact-count sampling.",
                figure_bullets=pair_count_bullets_plain,
                **shared_kw,
            )
        )
        lines.append("")
    if gzip is not None and not gzip.empty:
        lines.append(
            md_perf_section(
                gzip,
                title="Performance: pair layout, fraction (gzip)",
                intro=(
                    "The same paired vs interleaved `--fraction` mates as the plain section. "
                    "Gold is z-fastq ISA-L; amber is native Zig inflate. IRMA Core is not timed on gzip."
                    f"{pair_blank}"
                ),
                tools=GZIP_TOOLS,
                fig_name="perf_gzip_pairs_fraction.png",
                fig_title="Gzip pair-layout fraction sample: wall, RSS, page faults",
                fig_note=f"Error bars = zebrac standard deviation (n={sample_count}). Amber = native inflate. Hatched = descriptive samplers.",
                roles=FAMILY_ROLES["pairs_fraction"],
                include_throughput=True,
                figure_caption="Wall time, peak RSS, and minor page faults for gzip paired vs interleaved fraction sampling.",
                figure_bullets=pair_frac_bullets_gzip,
                **shared_kw,
            )
        )
        lines.append("")
        lines.append(
            md_perf_section(
                gzip,
                title="Performance: pair layout, count (gzip)",
                intro=(
                    "The same paired vs interleaved `--count` mates as the plain section. "
                    "Gold is z-fastq ISA-L; amber is native Zig inflate."
                    f"{pair_blank}"
                ),
                tools=GZIP_TOOLS,
                fig_name="perf_gzip_pairs_count.png",
                fig_title="Gzip pair-layout exact-count sample: wall, RSS, page faults",
                fig_note=f"Error bars = zebrac standard deviation (n={sample_count}). Amber = native inflate. Hatched = descriptive samplers.",
                roles=FAMILY_ROLES["pairs_count"],
                include_throughput=True,
                figure_caption="Wall time, peak RSS, and minor page faults for gzip paired vs interleaved exact-count sampling.",
                figure_bullets=pair_count_bullets_gzip,
                **shared_kw,
            )
        )
        lines.append("")

    output = squeeze_blank_lines("\n".join(lines).rstrip() + "\n")
    report_path = SCRIPT_DIR / "REPORT.md"
    report_path.write_text(output, encoding="utf-8")
    print(f"wrote {report_path}")

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", nargs="?", default=str(RESULTS_DIR))
    parser.add_argument(
        "--allow-incomplete",
        default="false",
        help="true/false: allow missing perf sections or a skipped sample contract preflight",
    )
    args = parser.parse_args()
    allow = str(args.allow_incomplete).lower() in {"1", "true", "yes"}
    generate(Path(args.results_dir), allow)


if __name__ == "__main__":
    main()
