#!/usr/bin/env python3
"""Create benchmark plots from aggregated CAMI metrics."""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict
from pathlib import Path

RANK_ORDER = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
SUPPORTED_RANKS = {
    "sourmash_gather": {"superkingdom", "species"},
}
PREFERRED_TOOL_ORDER = [
    "hymet",
    "kraken2",
    "centrifuge",
    "ganon2",
    "metaphlan4",
    "sourmash_gather",
]
PALETTE = [
    "#264653",
    "#2a9d8f",
    "#e9c46a",
    "#f4a261",
    "#e76f51",
    "#9d4edd",
    "#457b9d",
    "#f72585",
]


def load_table(path: Path, key_cols):
    if not path.is_file():
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        data = {}
        for row in reader:
            key = tuple(row[col] for col in key_cols)
            data[key] = row
    return data


def load_rows(path: Path):
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_runtime_rows(path: Path):
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
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: F401
        try:
            plt.style.use("seaborn-v0_8-whitegrid")
        except Exception:
            plt.style.use("ggplot")
        plt.rcParams.update(
            {
                "axes.facecolor": "#f4f6fb",
                "figure.facecolor": "white",
                "axes.edgecolor": "#d0d5dd",
                "grid.alpha": 0.3,
                "grid.linestyle": "--",
                "axes.titleweight": "semibold",
                "axes.labelweight": "semibold",
                "font.size": 11,
                "axes.titlesize": 15,
                "axes.labelsize": 12,
            }
        )
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("matplotlib is required to generate plots. Install it or skip plotting.") from exc


def order_ranks(ranks):
    ranks = list(ranks)
    ordered = [rank for rank in RANK_ORDER if rank in ranks]
    extras = [rank for rank in ranks if rank not in RANK_ORDER]
    return ordered + sorted(extras)


def get_tool_colors(tools):
    import matplotlib.pyplot as plt

    if not tools:
        return {}
    colors = PALETTE
    if len(tools) > len(colors):
        cmap = plt.get_cmap("tab20")
        colors = [cmap(i % cmap.N) for i in range(len(tools))]
    return {tool: colors[i % len(colors)] for i, tool in enumerate(tools)}


def order_tools(tools):
    ordered = []
    seen = set()
    for tool in PREFERRED_TOOL_ORDER:
        if tool in tools and tool not in seen:
            ordered.append(tool)
            seen.add(tool)
    for tool in tools:
        if tool not in seen:
            ordered.append(tool)
            seen.add(tool)
    return ordered


def mean(values):
    values = [v for v in values if v is not None]
    return sum(values) / len(values) if values else 0.0


def clean_axis(ax):
    ax.set_axisbelow(True)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.6, color="#d9dfe8")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#cccccc")
    ax.tick_params(axis="both", labelsize=10)


def bar_positions(x, idx, total, width):
    group_width = width * total
    return [xi - group_width / 2 + (idx + 0.5) * width for xi in x]


def add_bar_labels(ax, containers):
    total_bars = sum(len(container) for container in containers)
    if total_bars > 24:
        return
    ymin, ymax = ax.get_ylim()
    span = max(ymax - ymin, 1e-6)
    for container in containers:
        for bar in container:
            height = bar.get_height()
            if height is None or math.isnan(height):
                continue
            ax.annotate(
                f"{height:.1f}",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, -6 if height > ymin + span * 0.2 else 6),
                textcoords="offset points",
                ha="center",
                va="top" if height > ymin + span * 0.2 else "bottom",
                fontsize=9,
                color="#344054",
                fontweight="medium",
                clip_on=True,
            )


def place_legend(ax, pos="upper center", bbox=(0.5, 1.02), ncols=3):
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None
    legend = ax.legend(
        handles,
        labels,
        loc=pos,
        bbox_to_anchor=bbox,
        frameon=False,
        fontsize=10,
        ncol=ncols,
        columnspacing=0.8,
        handlelength=1.4,
        borderaxespad=0.6,
    )
    return legend


def place_grid_legend(fig, axes, tools, ncols=3, loc="upper center", bbox=(0.5, 1.05)):
    handles, labels = axes[0].get_legend_handles_labels()
    if not handles:
        return None
    legend = fig.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(0.88, 0.5),
        frameon=True,
        fontsize=10,
        ncol=1,
        columnspacing=0.8,
        handlelength=1.4,
    )
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.65)
    fig.subplots_adjust(right=0.84)
    return legend


def place_stack_legend(ax, tools, ncols=4):
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None
    legend = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.88, 1.0),
        frameon=True,
        fontsize=10,
        ncol=1,
        columnspacing=0.8,
        handlelength=1.4,
    )
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.65)
    ax.figure.subplots_adjust(right=0.84)
    return legend


def place_bar_legend(ax, handles=None, labels=None, loc="upper left", bbox=(0.02, 0.98)):
    if handles is None or labels is None:
        handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None
    legend = ax.legend(
        handles,
        labels,
        loc=loc,
        bbox_to_anchor=bbox,
        frameon=True,
        fontsize=10,
        ncol=1,
        columnspacing=0.8,
        handlelength=1.4,
    )
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.65)
    return legend



def plot_f1_by_rank(summary_rows, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    data = defaultdict(lambda: defaultdict(list))
    for row in summary_rows:
        tool = row["tool"]
        rank = row["rank"]
        allowed = SUPPORTED_RANKS.get(tool)
        if allowed and rank not in allowed:
            continue
        data[rank][tool].append(safe_float(row.get("F1_%")))

    ranks = order_ranks(data.keys())
    tools = order_tools([tool for tool in tool_colors if any(tool in rank_data for rank_data in data.values())])

    fig, ax = plt.subplots(figsize=(11, 5.6))
    width = 0.8 / max(1, len(tools))
    x = list(range(len(ranks)))
    containers = []
    for idx, tool in enumerate(tools):
        offsets = bar_positions(x, idx, len(tools), width)
        heights = []
        for rank in ranks:
            rank_values = data.get(rank, {}).get(tool)
            if rank_values:
                heights.append(mean(rank_values))
            else:
                heights.append(float('nan'))
        container = ax.bar(
            offsets,
            heights,
            width=width,
            label=tool,
            color=tool_colors.get(tool),
            edgecolor="white",
            linewidth=0.6,
        )
        containers.append(container)

    ax.set_xticks([xi for xi in x])
    ax.set_xticklabels(ranks, rotation=20)
    ax.set_ylabel("F1 (%)")
    ax.set_ylim(0, 100)
    clean_axis(ax)
    handles = [c[0] for c in containers]
    legend = place_bar_legend(ax, handles, tools, loc="upper right", bbox=(0.98, 0.98))
    if legend:
        legend.get_frame().set_linewidth(0.0)
        legend.borderpad = 0.35
    ax.set_xlim(min(x) - 0.6, max(x) + 0.6)
    add_bar_labels(ax, containers)
    fig.tight_layout(pad=0.4)
    fig.savefig(out_path)
    plt.close(fig)


def plot_l1_bray(summary_rows, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    tool_metrics = defaultdict(lambda: {"L1": defaultdict(list), "Bray": defaultdict(list)})
    for row in summary_rows:
        tool = row["tool"]
        rank = row["rank"]
        tool_metrics[tool]["L1"][rank].append(safe_float(row.get("L1_total_variation_pctpts")))
        tool_metrics[tool]["Bray"][rank].append(safe_float(row.get("BrayCurtis_pct")))

    tools = order_tools([tool for tool in tool_colors if tool in tool_metrics])
    ranks = order_ranks({rank for tm in tool_metrics.values() for rank in tm["L1"].keys()})

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, metric_key, label in zip(axes, ["L1", "Bray"], ["L1 total variation (pct-pts)", "Bray-Curtis (%)"]):
        for tool in tools:
            means = [mean(tool_metrics[tool][metric_key].get(rank, [])) for rank in ranks]
            ax.plot(
                ranks,
                means,
                marker="o",
                linewidth=2,
                markersize=5,
                label=tool,
                color=tool_colors.get(tool),
            )
        ax.set_title(label.title(), fontsize=13, pad=10)
        ax.set_xticks(range(len(ranks)))
        ax.set_xticklabels(ranks, rotation=25)
        clean_axis(ax)
        ax.set_ylim(bottom=0)
        ax.spines["left"].set_alpha(0.6)
        ax.spines["bottom"].set_alpha(0.6)
        ax.grid(axis="y", alpha=0.25, linestyle="--", linewidth=0.7)

    axes[0].set_ylabel("Mean value")
    legend = axes[0].legend(
        loc="upper left",
        bbox_to_anchor=(0.5, -0.25),
        frameon=True,
        fontsize=10,
        ncol=min(len(tools), 4),
        columnspacing=1.0,
        handlelength=1.4,
    )
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.65)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.26)
    fig.savefig(out_path)
    plt.close(fig)


def plot_accuracy(contig_rows, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    data = defaultdict(lambda: defaultdict(list))
    for row in contig_rows:
        if safe_float(row.get("n")) <= 0:
            continue
        tool = row["tool"]
        rank = row["rank"]
        allowed = SUPPORTED_RANKS.get(tool)
        if allowed and rank not in allowed:
            continue
        data[rank][tool].append(safe_float(row.get("accuracy_percent")))

    if not data:
        return

    ranks = order_ranks(data.keys())
    tools = order_tools([tool for tool in tool_colors if any(tool in rank_data for rank_data in data.values())])

    fig, ax = plt.subplots(figsize=(11, 5.6))
    width = 0.8 / max(1, len(tools))
    x = list(range(len(ranks)))
    containers = []
    for idx, tool in enumerate(tools):
        offsets = bar_positions(x, idx, len(tools), width)
        heights = []
        for rank in ranks:
            rank_values = data.get(rank, {}).get(tool)
            if rank_values:
                heights.append(mean(rank_values))
            else:
                heights.append(float("nan"))
        container = ax.bar(
            offsets,
            heights,
            width=width,
            label=tool,
            color=tool_colors.get(tool),
            edgecolor="white",
            linewidth=0.6,
        )
        containers.append(container)

    ax.set_xticks([xi for xi in x])
    ax.set_xticklabels(ranks, rotation=20)
    ax.set_ylabel("Contig Accuracy (%)")
    ax.set_ylim(0, 100)
    clean_axis(ax)
    handles = [c[0] for c in containers]
    place_bar_legend(ax, handles, tools, loc="upper right", bbox=(0.98, 0.98))
    ax.set_xlim(min(x) - 0.6, max(x) + 0.6)
    add_bar_labels(ax, containers)
    fig.tight_layout(pad=0.4)
    fig.savefig(out_path)
    plt.close(fig)


def plot_per_sample_stack(summary_rows, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    species_rows = [row for row in summary_rows if row.get("rank") == "species"]
    if not species_rows:
        species_rows = summary_rows

    samples = sorted({row["sample"] for row in species_rows})
    tools = order_tools([tool for tool in tool_colors if any(row["tool"] == tool for row in species_rows)])

    data = {sample: defaultdict(list) for sample in samples}
    for row in species_rows:
        data[row["sample"]][row["tool"]].append(safe_float(row.get("F1_%")))

    fig, ax = plt.subplots(figsize=(12, 5.6))
    bottoms = [0.0] * len(samples)
    x = list(range(len(samples)))
    for tool in tools:
        heights = []
        for sample in samples:
            values = data[sample].get(tool, [])
            heights.append(mean(values) if values else 0.0)
        ax.bar(x, heights, bottom=bottoms, label=tool, color=tool_colors.get(tool))
        bottoms = [b + h for b, h in zip(bottoms, heights)]

    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=25)
    ax.set_ylabel("Average F1 (%)")
    max_total = max(bottoms) if bottoms else 0.0
    ax.set_ylim(0, max_total * 1.05 if max_total else 100)
    ax.margins(y=0.04)
    legend = place_stack_legend(ax, tools, ncols=min(len(tools), 4))
    if legend:
        legend.set_bbox_to_anchor((0.98, 0.98))
        legend._loc = 1  # upper right
    clean_axis(ax)
    fig.tight_layout(pad=0.4)
    fig.savefig(out_path)
    plt.close(fig)


def summarise_runtime(runtime_rows):
    data = defaultdict(lambda: {"cpu_sec": [], "max_gb": []})
    for row in runtime_rows:
        if row.get("stage") != "run":
            continue
        tool = row.get("tool")
        if not tool:
            continue
        try:
            user = float(row.get("user_seconds", 0.0))
        except ValueError:
            user = 0.0
        try:
            sys_time = float(row.get("sys_seconds", 0.0))
        except ValueError:
            sys_time = 0.0
        try:
            rss = float(row.get("max_rss_gb", 0.0))
        except ValueError:
            rss = 0.0
        data[tool]["cpu_sec"].append(user + sys_time)
        data[tool]["max_gb"].append(rss)
    summary = {}
    for tool, metrics in data.items():
        cpu_vals = [v for v in metrics["cpu_sec"] if v >= 0]
        mem_vals = [v for v in metrics["max_gb"] if v >= 0]
        summary[tool] = {
            "cpu_min": mean(cpu_vals) / 60.0 if cpu_vals else 0.0,
            "max_gb": max(mem_vals) if mem_vals else 0.0,
        }
    return summary


def plot_runtime(summary, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    if not summary:
        return

    tools = order_tools(list(summary.keys()))
    values = [summary[t]["cpu_min"] for t in tools]
    x = list(range(len(tools)))

    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    bars = ax.bar(x, values, color=[tool_colors.get(t) for t in tools], edgecolor="white", linewidth=0.6)
    ax.set_ylabel("CPU time (minutes)")
    ax.set_xticks(x)
    ax.set_xticklabels(tools, rotation=15)
    clean_axis(ax)
    ax.set_ylim(0, max(values) * 1.15 if values else 1.0)
    ax.legend_.remove() if ax.get_legend() else None
    fig.tight_layout(pad=0.4)
    fig.savefig(out_path)
    plt.close(fig)


def plot_memory(summary, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    if not summary:
        return

    tools = order_tools(list(summary.keys()))
    values = [summary[t]["max_gb"] for t in tools]
    x = list(range(len(tools)))

    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    bars = ax.bar(x, values, color=[tool_colors.get(t) for t in tools], edgecolor="white", linewidth=0.6)
    ax.set_ylabel("Peak RSS (GB)")
    ax.set_xticks(x)
    ax.set_xticklabels(tools, rotation=15)
    clean_axis(ax)
    ax.set_ylim(0, max(values) * 1.15 if values else 1.0)
    ax.legend_.remove() if ax.get_legend() else None
    fig.tight_layout(pad=0.4)
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate plots from CAMI aggregate metrics.")
    ap.add_argument("--bench-root", default=str(Path(__file__).resolve().parent.parent), help="Bench directory root.")
    ap.add_argument("--outdir", default="out", help="Relative output directory for figures.")
    args = ap.parse_args()

    bench_root = Path(args.bench_root)
    out_root = bench_root / args.outdir
    out_root.mkdir(parents=True, exist_ok=True)

    summary_path = out_root / "summary_per_tool_per_sample.tsv"
    contig_path = out_root / "contig_accuracy_per_tool.tsv"

    summary_rows = load_rows(summary_path)
    contig_rows = load_rows(contig_path)
    runtime_rows = load_runtime_rows(out_root / "runtime_memory.tsv")
    if not summary_rows:
        print("[plot] No summary data available; skipping figure generation.")
        return

    ensure_matplotlib()

    tools = sorted({row["tool"] for row in summary_rows})
    if contig_rows:
        tools = sorted(set(tools).union({row["tool"] for row in contig_rows}))
    tool_colors = get_tool_colors(tools)

    plot_f1_by_rank(summary_rows, out_root / "fig_f1_by_rank.png", tool_colors)
    plot_l1_bray(summary_rows, out_root / "fig_l1_braycurtis.png", tool_colors)
    if contig_rows:
        plot_accuracy(contig_rows, out_root / "fig_accuracy_by_rank.png", tool_colors)
    plot_per_sample_stack(summary_rows, out_root / "fig_per_sample_f1_stack.png", tool_colors)

    runtime_summary = summarise_runtime(runtime_rows)
    if runtime_summary:
        plot_runtime(runtime_summary, out_root / "fig_cpu_time_by_tool.png", tool_colors)
        plot_memory(runtime_summary, out_root / "fig_peak_memory_by_tool.png", tool_colors)


if __name__ == "__main__":
    main()
