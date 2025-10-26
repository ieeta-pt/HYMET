#!/usr/bin/env python3
"""Generate combined contig/read figures for CAMI suite runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from matplotlib.ticker import FuncFormatter

PREFERRED_TOOL_ORDER = [
    "hymet",
    "hymet_reads",
    "kraken2",
    "centrifuge",
    "ganon2",
    "metaphlan4",
    "sourmash_gather",
    "camitax",
    "basta",
    "tama",
    "phabox",
    "phyloflash",
    "viwrap",
    "squeezemeta",
    "megapath_nano",
]
MODE_COLORS = {
    "contigs": "#2563eb",
    "reads": "#ec6f20",
}

RANK_ORDER = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
LINE_STYLES = ["solid", "dashdot", (0, (3, 1)), "dotted", (0, (1, 1))]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", required=True, help="Path to suite tables/ directory")
    parser.add_argument("--outdir", required=True, help="Directory for combined figures")
    parser.add_argument("--run", default="", help="Run directory (metadata label)")
    parser.add_argument("--suite", default="", help="Suite name for figure titles")
    return parser.parse_args()


def load_rows(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def ensure_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: F401

    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except Exception:
        plt.style.use("ggplot")

    plt.rcParams.update(
        {
            "axes.facecolor": "#fbfdff",
            "figure.facecolor": "white",
            "axes.edgecolor": "#d0d5dd",
            "grid.alpha": 0.3,
            "grid.linestyle": (0, (4, 4)),
            "axes.titleweight": "bold",
            "axes.labelweight": "semibold",
            "font.size": 11,
            "axes.titlesize": 15,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        }
    )
    return matplotlib.pyplot


def beautify_axis(ax, grid_axis="y"):
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#d2d7e3")
    ax.grid(axis=grid_axis, color="#dbe2f0", linewidth=0.7, alpha=0.8, linestyle="--")
    ax.tick_params(colors="#4a5568")


def order_tools(tools: Iterable[str]) -> List[str]:
    seen = []
    remaining = [tool for tool in tools if tool]
    for tool in PREFERRED_TOOL_ORDER:
        if tool in remaining:
            seen.append(tool)
            remaining.remove(tool)
    seen.extend(sorted(remaining))
    return seen


def load_tool_modes(run_path: Path | None) -> Dict[str, set]:
    mapping: Dict[str, set] = defaultdict(set)
    if not run_path:
        return mapping
    meta_path = run_path / "metadata.json"
    if not meta_path.is_file():
        return mapping
    try:
        metadata = json.loads(meta_path.read_text())
    except Exception:
        return mapping

    def register(field: str, mode_label: str):
        raw = metadata.get(field, "")
        if isinstance(raw, str):
            tokens = raw.split(",")
        elif isinstance(raw, list):
            tokens = raw
        else:
            return
        for token in tokens:
            tool = (token or "").strip()
            if tool:
                mapping[tool].add(mode_label)

    register("contig_tools", "contigs")
    register("read_tools", "reads")
    return mapping


def order_modes(modes: Iterable[str]) -> List[str]:
    ordered = []
    remaining = [mode for mode in modes if mode]
    for preferred in ("contigs", "reads"):
        if preferred in remaining:
            ordered.append(preferred)
            remaining.remove(preferred)
    ordered.extend(sorted(remaining))
    return ordered


def get_tool_colors(tools: Iterable[str]) -> Dict[str, str]:
    palette = [
        "#2563eb",
        "#f97316",
        "#059669",
        "#a855f7",
        "#ef4444",
        "#0ea5e9",
        "#f59e0b",
        "#14b8a6",
        "#6d28d9",
        "#d946ef",
        "#22c55e",
        "#f43f5e",
    ]
    ordered = order_tools(tools)
    return {tool: palette[i % len(palette)] for i, tool in enumerate(ordered)}


def collect_runtime(runtime_rows: List[Dict[str, str]]) -> Tuple[Dict[Tuple[str, str], float], List[str], List[str]]:
    totals = defaultdict(float)
    modes = set()
    tools = set()
    tool_totals = defaultdict(float)
    for row in runtime_rows:
        stage = (row.get("stage") or "").strip()
        if stage and stage not in {"run", "eval"}:
            continue
        if stage == "eval":
            continue
        tool = (row.get("tool") or "").strip() or "unknown"
        mode = (row.get("mode") or "").strip() or "unknown"
        wall_seconds = safe_float(row.get("wall_seconds"))
        hours = wall_seconds / 3600.0
        totals[(tool, mode)] += hours
        tool_totals[tool] += hours
        tools.add(tool)
        modes.add(mode)
    ordered_tools = sorted(order_tools(tools), key=lambda t: tool_totals.get(t, 0.0), reverse=True)
    return totals, ordered_tools, order_modes(modes)


def collect_rank_f1(summary_rows: List[Dict[str, str]], allowed_modes: Dict[str, set] | None = None) -> Tuple[Dict[str, Dict[str, Dict[str, float]]], Dict[str, List[str]]]:
    values = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for row in summary_rows:
        rank = (row.get("rank") or "").strip().lower()
        if rank not in RANK_ORDER:
            continue
        mode = (row.get("mode") or "").strip() or "unknown"
        tool = (row.get("tool") or "").strip() or "unknown"
        if allowed_modes:
            modes_expected = allowed_modes.get(tool)
            if modes_expected and mode not in modes_expected:
                continue
        score = safe_float(row.get("F1_%"))
        # Include zero scores as valid points; they matter for visibility
        values[mode][tool][rank].append(score)
    averages: Dict[str, Dict[str, Dict[str, float]]] = {}
    tool_lists: Dict[str, List[str]] = {}
    for mode, tool_map in values.items():
        averages[mode] = {}
        ordered_tools = order_tools(tool_map.keys())
        tool_lists[mode] = ordered_tools
        for tool in ordered_tools:
            rank_map = tool_map[tool]
            averages[mode][tool] = {rank: sum(vals) / len(vals) for rank, vals in rank_map.items()}
    return averages, tool_lists


def collect_contig_accuracy(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, float]]:
    per_tool = defaultdict(lambda: defaultdict(list))
    for row in rows:
        tool = (row.get("tool") or "").strip()
        rank = (row.get("rank") or "").strip().lower()
        if not tool or rank not in RANK_ORDER:
            continue
        acc = safe_float(row.get("accuracy_percent"))
        per_tool[tool][rank].append(acc)
    averages = {
        tool: {rank: sum(vals) / len(vals) for rank, vals in rank_map.items() if vals}
        for tool, rank_map in per_tool.items()
    }
    return {tool: averages[tool] for tool in order_tools(averages.keys())}


def collect_abundance_errors(summary_rows: List[Dict[str, str]], allowed_modes: Dict[str, set] | None = None) -> Dict[Tuple[str, str], Dict[str, float]]:
    metrics = defaultdict(lambda: defaultdict(list))
    for row in summary_rows:
        tool = (row.get("tool") or "").strip()
        mode = (row.get("mode") or "").strip() or "unknown"
        if allowed_modes:
            modes_expected = allowed_modes.get(tool)
            if modes_expected and mode not in modes_expected:
                continue
        if not tool:
            continue
        key = (tool, mode)
        metrics[key]["L1"].append(safe_float(row.get("L1_total_variation_pctpts")))
        metrics[key]["BrayCurtis"].append(safe_float(row.get("BrayCurtis_pct")))
    averages = {}
    for key, m in metrics.items():
        averages[key] = {
            metric: (sum(vals) / len(vals) if vals else 0.0)
            for metric, vals in m.items()
        }
    return averages


def resolve_color(mode: str, idx: int) -> str:
    palette = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
    ]
    if mode in MODE_COLORS:
        return MODE_COLORS[mode]
    return palette[idx % len(palette)]


def format_hours(value: float) -> str:
    if value >= 1:
        return f"{value:.1f} h"
    if value <= 0:
        return "0 m"
    minutes = value * 60
    return f"{minutes:.0f} m" if minutes >= 1 else f"{minutes * 60:.0f} s"


def format_minutes(value: float) -> str:
    if value >= 60:
        hours = value / 60.0
        return f"{hours:.1f} h" if hours < 10 else f"{hours:.0f} h"
    if value >= 1:
        return f"{value:.0f} m"
    if value <= 0:
        return "0 s"
    return f"{value*60:.0f} s"


def annotate_horizontal_bars(ax, bars, formatter):
    for bar in bars:
        width = bar.get_width()
        if width <= 0:
            continue
        label = formatter(width)
        ax.annotate(
            label,
            xy=(width, bar.get_y() + bar.get_height() / 2),
            xytext=(8, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=9,
            color="#1f2933",
            fontweight="semibold",
        )


def plot_runtime_bars(data: Dict[Tuple[str, str], float], tools: List[str], modes: List[str], out_path: Path):
    if not data or not tools or not modes:
        return False
    plt = ensure_matplotlib()
    # Flatten into separate rows per (tool, mode) so each bar is centered on its own label
    labels: List[str] = []
    y_for_mode: Dict[str, List[int]] = {mode: [] for mode in modes}
    vals_for_mode: Dict[str, List[float]] = {mode: [] for mode in modes}
    y = 0
    present_any = False
    for tool in tools:
        # For consistent grouping, keep contigs first
        for mode in modes:
            hours = data.get((tool, mode))
            if hours is None or hours <= 0:
                continue
            minutes = hours * 60.0
            labels.append(tool)
            y_for_mode[mode].append(y)
            vals_for_mode[mode].append(minutes)
            y += 1
            present_any = True
    if not present_any:
        return False

    fig_height = max(3.2, len(labels) * 0.45)
    fig, ax = plt.subplots(figsize=(9.6, fig_height))

    containers = []
    for idx, mode in enumerate(modes):
        if not y_for_mode[mode]:
            continue
        bars = ax.barh(
            y_for_mode[mode],
            vals_for_mode[mode],
            height=0.36,
            color=resolve_color(mode, idx),
            edgecolor="#e5e7eb",
            linewidth=0.9,
            label=mode.capitalize(),
        )
        for bar in bars:
            ax.scatter(bar.get_width(), bar.get_y() + bar.get_height() / 2.0, s=18, color="#ffffff", zorder=3)
            ax.scatter(bar.get_width(), bar.get_y() + bar.get_height() / 2.0, s=26, color=resolve_color(mode, idx), zorder=3, edgecolor="#ffffff", linewidth=0.6)
        containers.append(bars)

    # Axis scale/ticks (symlog so 0 is represented at the origin)
    all_minutes = [v for vs in vals_for_mode.values() for v in vs]
    max_pos = max(all_minutes) if all_minutes else 1.0
    ax.set_xscale("symlog", linthresh=0.1, linscale=1.0, base=10)
    left = 0.0
    right = max_pos * 1.15
    ax.set_xlim(left, right)

    def _select_ticks(lo: float, hi: float):
        base = [0.0, 0.1, 0.5, 1, 2, 5, 10, 20, 30, 60, 120, 240, 480, 960]
        ticks = [t for t in base if lo <= t <= hi]
        if len(ticks) > 8:
            step = (len(ticks) + 7) // 8
            ticks = ticks[::step]
        if ticks and ticks[0] != 0.0 and lo == 0.0:
            ticks = [0.0] + ticks
        return ticks

    ticks = _select_ticks(left, right)
    if ticks:
        ax.set_xticks(ticks)
        def _fmt(x, pos):
            if abs(x) < 1e-12:
                return "0 s"
            return format_minutes(x)
        ax.get_xaxis().set_major_formatter(FuncFormatter(_fmt))
    ax.minorticks_off()

    ax.set_yticks(list(range(len(labels))))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Wall time (minutes)")
    beautify_axis(ax, grid_axis="x")
    ax.tick_params(axis="y", labelsize=10, pad=6)
    for bars in containers:
        annotate_horizontal_bars(ax, bars, lambda m: format_minutes(m))

    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(0.98, 1.02), ncol=len(modes), handlelength=2.5)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return True


def plot_contig_accuracy_lines(data: Dict[str, Dict[str, float]], out_path: Path):
    if not data:
        return False
    plt = ensure_matplotlib()
    ranks = RANK_ORDER
    x_positions = list(range(len(ranks)))
    tools = list(data.keys())
    colors = get_tool_colors(tools)

    fig, ax = plt.subplots(figsize=(10.5, 5.0))
    for idx, tool in enumerate(tools):
        rank_map = data[tool]
        y_vals = [rank_map.get(rank, float("nan")) for rank in ranks]
        ax.plot(
            x_positions,
            y_vals,
            linewidth=2.0,
            color=colors.get(tool, "#4f46e5"),
            label=tool,
        )
        ax.scatter(
            x_positions,
            y_vals,
            s=24,
            color=colors.get(tool, "#4f46e5"),
            edgecolor="white",
            linewidth=0.6,
        )
    ax.set_xticks(x_positions)
    ax.set_xticklabels([rank.capitalize() for rank in ranks], rotation=20, ha="right")
    ax.set_ylabel("Contig accuracy (%)")
    beautify_axis(ax)
    ax.set_ylim(0, 105)
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(0.98, 1.02), ncol=2)
    ax.text(0.01, 1.02, "Higher = better; % of contigs correctly classified per rank", transform=ax.transAxes, fontsize=9, color="#475467")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return True


def plot_abundance_bars(metrics: Dict[Tuple[str, str], Dict[str, float]], out_path: Path):
    if not metrics:
        return False
    plt = ensure_matplotlib()
    entries = list(metrics.items())
    entries.sort(key=lambda kv: order_tools([kv[0][0]])[0])
    labels = [f"{tool} ({mode})" for (tool, mode), _ in entries]
    l1_vals = [kv[1].get("L1", 0.0) for kv in entries]
    bc_vals = [kv[1].get("BrayCurtis", 0.0) for kv in entries]
    y_positions = list(range(len(labels)))

    fig, ax = plt.subplots(figsize=(9.5, max(3.0, len(labels) * 0.45)))
    ax.barh([y - 0.2 for y in y_positions], l1_vals, height=0.35, color="#2563eb", alpha=0.9, label="L1 deviation")
    ax.barh([y + 0.2 for y in y_positions], bc_vals, height=0.35, color="#f97316", alpha=0.9, label="Bray–Curtis")
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Abundance error (percentage points)")
    beautify_axis(ax, grid_axis="x")
    ax.tick_params(axis="y", labelsize=10)
    ax.legend(frameon=False, loc="upper right")
    ax.text(0.01, 1.02, "Lower = better; abundance deviation from truth", transform=ax.transAxes, fontsize=9, color="#475467")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return True


def plot_rank_lines(rank_values: Dict[str, Dict[str, Dict[str, float]]], tool_lists: Dict[str, List[str]], out_path: Path):
    if not rank_values:
        return False
    plt = ensure_matplotlib()
    ranks = RANK_ORDER
    x_positions = list(range(len(ranks)))
    fig, ax = plt.subplots(figsize=(10.5, 5.0))

    label_offsets = {
        "contigs": (6, -2),
        "reads": (6, 2),
    }

    label_register: Dict[int, List[float]] = defaultdict(list)

    def nudge_for(x_idx: int, y_value: float) -> float:
        min_gap = 6.0
        offsets = [0, 6, -6, 12, -12, 18, -18, 24, -24, 30, -30]
        for delta in offsets:
            candidate = y_value + delta
            if all(abs(candidate - existing) >= min_gap for existing in label_register[x_idx]):
                label_register[x_idx].append(candidate)
                return delta
        label_register[x_idx].append(y_value + offsets[-1])
        return offsets[-1]

    ordered_modes = order_modes(rank_values.keys())
    best_mean = -1.0
    best_point = None
    best_color = None
    best_mode = None
    best_tool_name = None
    for idx, mode in enumerate(ordered_modes):
        tools = tool_lists.get(mode, [])
        if not tools:
            continue
        color = resolve_color(mode, idx)
        for tool_idx, tool in enumerate(tools):
            rank_map = rank_values[mode].get(tool, {})
            if not rank_map:
                continue
            y_vals = [rank_map.get(rank, float("nan")) for rank in ranks]
            style = LINE_STYLES[tool_idx % len(LINE_STYLES)]
            ax.plot(
                x_positions,
                y_vals,
                color=color,
                linewidth=2.0,
                linestyle=style,
                alpha=0.9,
            )
            ax.scatter(
                x_positions,
                y_vals,
                s=26,
                color=color,
                edgecolor="white",
                linewidth=0.6,
                alpha=0.95,
            )
            valid_vals = [val for val in y_vals if not math.isnan(val)]
            tool_mean = None
            if valid_vals:
                tool_mean = sum(valid_vals) / len(valid_vals)
                if tool_mean > best_mean:
                    for pos in reversed(range(len(y_vals))):
                        val = y_vals[pos]
                        if math.isnan(val):
                            continue
                        best_mean = tool_mean
                        best_point = (x_positions[pos], val)
                        best_color = color
                        best_mode = mode
                        best_tool_name = tool
                        break
            # Label each tool at its last finite point with mean F1
            for pos in reversed(range(len(y_vals))):
                val = y_vals[pos]
                if math.isnan(val):
                    continue
                dx, dy = label_offsets.get(mode, (6, -6))
                dy += nudge_for(pos, val)
                label_mean = tool_mean if tool_mean is not None else None
                label_text = f"{tool} {label_mean:.1f}%" if label_mean is not None else tool
                ax.annotate(
                    label_text,
                    xy=(x_positions[pos], val),
                    xytext=(dx, dy),
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    fontsize=9,
                    color=color,
                    fontweight="bold" if (best_tool_name == tool and best_mode == mode) else "normal",
                )
                break

    ax.set_xticks(x_positions)
    ax.set_xticklabels([rank.capitalize() for rank in ranks], rotation=20, ha="right")
    ax.set_ylabel("Mean F1 (%)")
    beautify_axis(ax)
    ax.set_ylim(0, 105)
    legend = ax.legend(
        [arg for arg in []],
        [],
    )
    legend_handles = [
        plt.Line2D([0], [0], color=resolve_color(mode, idx), linewidth=2.5)
        for idx, mode in enumerate(ordered_modes)
    ]
    if best_point and best_color:
        ax.scatter(
            best_point[0],
            best_point[1],
            marker="*",
            s=160,
            color=best_color,
            edgecolor="white",
            linewidth=0.9,
            zorder=5,
        )

    ax.legend(
        legend_handles,
        [mode.capitalize() for mode in ordered_modes],
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.98, 1.02),
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return True


def main():
    args = parse_args()
    tables_dir = Path(args.tables).resolve()
    combined_dir = tables_dir / "combined"
    outdir = Path(args.outdir).resolve()
    run_path = Path(args.run).resolve() if args.run else None
    tool_modes = load_tool_modes(run_path)

    runtime_rows = load_rows(tables_dir / "runtime_memory.tsv")
    summary_rows = load_rows(combined_dir / "summary_per_tool_per_sample.tsv")
    contig_rows = load_rows(combined_dir / "contig_accuracy_per_tool.tsv")

    generated_any = False

    runtime_data, runtime_tools, runtime_modes = collect_runtime(runtime_rows)
    runtime_path = outdir / "fig_runtime_by_tool_combined.png"
    if plot_runtime_bars(runtime_data, runtime_tools, runtime_modes, runtime_path):
        generated_any = True

    rank_values, tool_lists = collect_rank_f1(summary_rows, tool_modes)
    rank_path = outdir / "fig_f1_by_rank_combined.png"
    if plot_rank_lines(rank_values, tool_lists, rank_path):
        generated_any = True

    contig_data = collect_contig_accuracy(contig_rows)
    contig_path = outdir / "fig_contig_accuracy_by_rank_combined.png"
    if plot_contig_accuracy_lines(contig_data, contig_path):
        generated_any = True

    abundance_metrics = collect_abundance_errors(summary_rows, tool_modes)
    abundance_path = outdir / "fig_abundance_error_combined.png"
    if plot_abundance_bars(abundance_metrics, abundance_path):
        generated_any = True

    if not generated_any:
        print("[combined] No figures generated (missing data?)")


if __name__ == "__main__":
    main()
