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
    "camitax",
    "basta",
    "phabox",
    "phyloflash",
    "viwrap",
    "squeezemeta",
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

# Whether an emoji-capable font was detected and enabled
HAS_EMOJI_FONT = False


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
        # Optional: allow higher-resolution PNGs via env var
        dpi_env = os.environ.get("HYMET_PLOTS_DPI")
        if dpi_env:
            try:
                dpi_val = int(float(dpi_env))
                if dpi_val > 0:
                    plt.rcParams.update({
                        "savefig.dpi": dpi_val,
                        "figure.dpi": dpi_val,
                    })
            except Exception:
                pass

        # Try to include an emoji-capable font for labels if available
        try:
            from matplotlib import font_manager as _fm

            available = {f.name for f in _fm.fontManager.ttflist}
            candidates = [
                "Noto Color Emoji",
                "Segoe UI Emoji",
                "Apple Color Emoji",
                "Twitter Color Emoji",
                "Symbola",
                "EmojiOne Color",
            ]
            chosen = next((name for name in candidates if name in available), None)
            if chosen:
                global HAS_EMOJI_FONT
                HAS_EMOJI_FONT = True
                existing = list(plt.rcParams.get("font.sans-serif", []))
                new_list = ["DejaVu Sans", chosen] + [n for n in existing if n not in ("DejaVu Sans", chosen)]
                plt.rcParams["font.sans-serif"] = new_list
                plt.rcParams.setdefault("font.family", ["sans-serif"])  # ensure sans fallback
        except Exception:
            # If font discovery fails, keep defaults; crown may fall back to mono
            pass
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

def plot_f1_by_rank_lines(summary_rows, out_path: Path, tool_colors):
    """Plot mean F1 by taxonomic rank with one line per tool.

    Enhancements:
    - Label each line at the last rank with the tool name
    - Highlight the top-3 tools (by overall mean F1) with thicker lines
    """
    import matplotlib.pyplot as plt
    import math as _math
    from matplotlib import patheffects as _pe

    # Collect F1 values per (rank, tool) across samples
    data = defaultdict(lambda: defaultdict(list))
    for row in summary_rows:
        tool = row["tool"]
        rank = row["rank"]
        allowed = SUPPORTED_RANKS.get(tool)
        if allowed and rank not in allowed:
            continue
        data[rank][tool].append(safe_float(row.get("F1_%")))

    if not data:
        return

    ranks = order_ranks(data.keys())
    tools = order_tools([t for t in tool_colors if any(t in rank_data for rank_data in data.values())])
    x = list(range(len(ranks)))

    # Prepare y-series per tool and compute overall means for ranking
    yseries = {}
    for tool in tools:
        vals = []
        for rank in ranks:
            vlist = data.get(rank, {}).get(tool)
            vals.append(mean(vlist) if vlist else float("nan"))
        yseries[tool] = vals

    def _nanmean(arr):
        arr2 = [v for v in arr if not _math.isnan(v)]
        return sum(arr2) / len(arr2) if arr2 else 0.0

    tool_means = {t: _nanmean(y) for t, y in yseries.items()}
    ranked_tools = sorted(tool_means.items(), key=lambda kv: kv[1], reverse=True)
    top_tools = [t for t, _ in ranked_tools[:3]]
    best_tool = ranked_tools[0][0] if ranked_tools else None

    fig, ax = plt.subplots(figsize=(11.5, 6.0))

    # Draw a smooth, elegant line with markers for each tool
    for tool in tools:
        y_vals = yseries[tool]
        lw = 3.2 if tool in top_tools else 1.9
        alpha = 1.0 if tool in top_tools else 0.9
        ax.plot(
            x,
            y_vals,
            label=tool,
            color=tool_colors.get(tool),
            linewidth=lw,
            marker="o",
            markersize=6.5 if tool in top_tools else 5.8,
            markerfacecolor=tool_colors.get(tool),
            markeredgecolor="white",
            markeredgewidth=1.0 if tool in top_tools else 0.9,
            alpha=alpha,
        )

    # Label each line at the last available point with non-overlapping labels
    # 1) Gather endpoints
    end_info = []  # (tool, end_x, end_y)
    for tool in tools:
        y_vals = yseries[tool]
        valid_idx = [i for i, v in enumerate(y_vals) if not _math.isnan(v)]
        if not valid_idx:
            continue
        end_x = valid_idx[-1]
        end_y = y_vals[end_x]
        end_info.append((tool, end_x, end_y))

    if end_info:
        # 2) Compute separated target y-positions in data coords
        end_info.sort(key=lambda t: t[2])  # sort by desired y (ascending)
        ymin, ymax = 0.0, 100.0
        margin = 2.0
        avail = max((ymax - ymin) - 2 * margin, 1.0)
        nlab = len(end_info)
        base_sep = 2.4
        min_sep = min(base_sep, avail / max(nlab - 1, 1)) if nlab > 1 else 0

        y_adj = []
        for i, (_, _, y) in enumerate(end_info):
            if i == 0:
                y_adj.append(max(y, ymin + margin))
            else:
                y_adj.append(max(y, y_adj[i - 1] + min_sep))
        # Pull back into bounds if we exceeded the top
        overflow = y_adj[-1] - (ymax - margin)
        if overflow > 0:
            y_adj = [yi - overflow for yi in y_adj]
            # Ensure bottom bound after shift
            under = (ymin + margin) - min(y_adj)
            if under > 0:
                y_adj = [yi + under for yi in y_adj]

        # 3) Place labels to the right with connectors
        label_x = (len(ranks) - 1) + 0.65
        # Replace crown with a widely supported black star
        star_symbol = "★"
        label_annotations = {}
        for (tool, ex, ey), ly in zip(end_info, y_adj):
            # Best tool gets a star after the name
            label_text = f"{tool} {star_symbol}" if (best_tool and tool == best_tool) else tool
            txt = ax.annotate(
                label_text,
                xy=(ex, ey),
                xytext=(label_x, ly),
                textcoords="data",
                ha="left",
                va="center",
                fontsize=10.5 if tool in top_tools else 9.8,
                color=tool_colors.get(tool),
                fontweight="semibold" if tool in top_tools else "normal",
                bbox=dict(facecolor="white", edgecolor="none", boxstyle="round,pad=0.2", alpha=0.75),
                arrowprops=dict(
                    arrowstyle="-",
                    color=tool_colors.get(tool),
                    lw=1.2 if tool in top_tools else 0.8,
                    alpha=0.6,
                    shrinkA=0,
                    shrinkB=0,
                    relpos=(0, 0.5),
                ),
            )
            txt.set_path_effects([_pe.withStroke(linewidth=2.5, foreground="white", alpha=0.9)])
            label_annotations[tool] = txt

    # Axes styling and room for labels on the right
    ax.set_xticks(x)
    ax.set_xticklabels(ranks, rotation=20)
    ax.set_ylabel("Mean F1 (%)")
    ax.set_ylim(0, 100)
    ax.set_xlim(min(x) - 0.2, (len(ranks) - 1) + 1.25)
    clean_axis(ax)
    # Labels make the legend redundant here; hide if present
    if ax.get_legend():
        ax.get_legend().remove()

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
    if not tools or not ranks:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), sharey=True, gridspec_kw={"wspace": 0.08})

    metric_series = {
        metric_key: {
            tool: [mean(tool_metrics[tool][metric_key].get(rank, [])) for rank in ranks]
            for tool in tools
        }
        for metric_key in ("L1", "Bray")
    }
    global_max = max(
        (
            max(values)
            for series in metric_series.values()
            for values in series.values()
            if values
        ),
        default=0.0,
    )
    top_limit = global_max * 1.05 if global_max > 0 else 1.0

    for ax, metric_key, label in zip(axes, ["L1", "Bray"], ["L1 total variation (pct-pts)", "Bray-Curtis (%)"]):
        for tool in tools:
            means = metric_series[metric_key][tool]
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
        ax.set_ylim(0, top_limit)
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


def plot_l1_bray_lines(summary_rows, out_path: Path, tool_colors):
    """Two-panel line plot for L1 total variation and Bray-Curtis with end labels.

    Style mirrors other line plots:
    - One line per tool across ranks
    - Top-3 tools highlighted (per metric)
    - Best tool marked with a star at end label (per metric)
    - Non-overlapping end labels with connectors
    """
    import matplotlib.pyplot as plt
    import math as _math
    from matplotlib import patheffects as _pe

    # Aggregate metric values per tool per rank
    tool_metrics = defaultdict(lambda: {"L1": defaultdict(list), "Bray": defaultdict(list)})
    for row in summary_rows:
        tool = row["tool"]
        rank = row["rank"]
        tool_metrics[tool]["L1"][rank].append(safe_float(row.get("L1_total_variation_pctpts")))
        tool_metrics[tool]["Bray"][rank].append(safe_float(row.get("BrayCurtis_pct")))

    tools = order_tools([t for t in tool_colors if t in tool_metrics])
    ranks = order_ranks({r for tm in tool_metrics.values() for r in tm["L1"].keys() | tm["Bray"].keys()})
    if not tools or not ranks:
        return

    # Build series per metric
    def _mean_seq(metric_key):
        return {
            tool: [mean(tool_metrics[tool][metric_key].get(rank, [])) for rank in ranks]
            for tool in tools
        }

    series = {k: _mean_seq(k) for k in ("L1", "Bray")}

    def _nanmean(arr):
        arr2 = [v for v in arr if not _math.isnan(v)]
        return sum(arr2) / len(arr2) if arr2 else 0.0

    # For error/distance measures, lower is better. Rank accordingly per metric.
    higher_is_better = {"L1": False, "Bray": False}
    rankings = {
        k: sorted(
            ((t, _nanmean(vals)) for t, vals in ser.items()),
            key=lambda kv: kv[1],
            reverse=higher_is_better.get(k, True),
        )
        for k, ser in series.items()
    }

    x = list(range(len(ranks)))
    # Increase wspace to keep end labels from overlapping across panels
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.6, 6.0),
        sharey=True,
        gridspec_kw={"wspace": 0.22},
    )

    # Establish a shared y-limit that comfortably accommodates end labels
    def _finite_max(seq):
        vals = [v for v in seq if not _math.isnan(v)]
        return max(vals) if vals else 0.0
    global_max = 0.0
    for mkey, ser in series.items():
        for _tool, vals in ser.items():
            global_max = max(global_max, _finite_max(vals))
    top_limit = (max(100.0, global_max) * 1.05) if global_max > 0 else 100.0

    labels = {"L1": "L1 total variation (pct-pts)", "Bray": "Bray-Curtis (%)"}
    star_symbol = "★"

    for ax, key in zip(axes, ("L1", "Bray")):
        ser = series[key]
        ranked = rankings[key]
        top_tools = [t for t, _ in ranked[:3]]
        best_tool = ranked[0][0] if ranked else None

        for tool in tools:
            y_vals = ser[tool]
            lw = 3.2 if tool in top_tools else 1.9
            alpha = 1.0 if tool in top_tools else 0.9
            ax.plot(
                x,
                y_vals,
                label=tool,
                color=tool_colors.get(tool),
                linewidth=lw,
                marker="o",
                markersize=6.5 if tool in top_tools else 5.8,
                markerfacecolor=tool_colors.get(tool),
                markeredgecolor="white",
                markeredgewidth=1.0 if tool in top_tools else 0.9,
                alpha=alpha,
            )

        # Use a shared, bottom legend to avoid text overlap in dense plots

        ax.set_title(labels[key], fontsize=13, pad=10)
        ax.set_xticks(x)
        ax.set_xticklabels(ranks, rotation=25)
        ax.set_xlim(min(x) - 0.2, (len(ranks) - 1) + 1.25)
        ax.set_ylim(0, top_limit)
        ax.grid(axis="y", alpha=0.25, linestyle="--", linewidth=0.7)
        ax.spines["left"].set_alpha(0.6)
        ax.spines["bottom"].set_alpha(0.6)

    axes[0].set_ylabel("Mean value")
    for ax in axes:
        clean_axis(ax)

    # Build a clean, de-duplicated legend from the first panel
    best_set = set()
    if rankings.get("L1"):
        best_set.add(rankings["L1"][0][0])
    if rankings.get("Bray"):
        best_set.add(rankings["Bray"][0][0])

    handles, labels = axes[0].get_legend_handles_labels()
    uniq = {}
    new_h, new_l = [], []
    for h, l in zip(handles, labels):
        if l in uniq:
            continue
        uniq[l] = True
        new_h.append(h)
        new_l.append(f"{l} {star_symbol}" if l in best_set else l)

    fig.subplots_adjust(bottom=0.26)
    fig_legend = fig.legend(
        new_h,
        new_l,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),
        frameon=True,
        fontsize=9,
        ncol=min(len(new_l), 6),
        columnspacing=0.9,
        handlelength=1.4,
    )
    if fig_legend:
        fig_legend.get_frame().set_facecolor("white")
        fig_legend.get_frame().set_alpha(0.75)
        fig_legend.get_frame().set_linewidth(0.0)

    fig.savefig(out_path, bbox_inches="tight")
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


def plot_accuracy_by_rank_lines(contig_rows, out_path: Path, tool_colors):
    """Line plot: mean contig accuracy (%) by rank with one line per tool.

    Mirrors the style of plot_f1_by_rank_lines:
    - Non-overlapping end labels with connectors
    - Top-3 tools highlighted with thicker lines
    - Best tool marked with a star after its name
    """
    import matplotlib.pyplot as plt
    import math as _math
    from matplotlib import patheffects as _pe

    # Collect accuracy per (rank, tool); only include rows with positive n
    data = defaultdict(lambda: defaultdict(list))
    for row in contig_rows:
        try:
            n_val = float(row.get("n", 0))
        except Exception:
            n_val = 0.0
        if n_val <= 0:
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
    tools = order_tools([t for t in tool_colors if any(t in rank_data for rank_data in data.values())])
    x = list(range(len(ranks)))

    # Prepare series and compute overall tool means
    def _nanmean(arr):
        arr2 = [v for v in arr if not _math.isnan(v)]
        return sum(arr2) / len(arr2) if arr2 else 0.0

    yseries = {}
    for tool in tools:
        vals = []
        for rank in ranks:
            vlist = data.get(rank, {}).get(tool)
            vals.append(mean(vlist) if vlist else float("nan"))
        yseries[tool] = vals

    tool_means = {t: _nanmean(y) for t, y in yseries.items()}
    ranked_tools = sorted(tool_means.items(), key=lambda kv: kv[1], reverse=True)
    top_tools = [t for t, _ in ranked_tools[:3]]
    best_tool = ranked_tools[0][0] if ranked_tools else None

    fig, ax = plt.subplots(figsize=(11.5, 6.0))

    # Draw lines
    for tool in tools:
        y_vals = yseries[tool]
        lw = 3.2 if tool in top_tools else 1.9
        alpha = 1.0 if tool in top_tools else 0.9
        ax.plot(
            x,
            y_vals,
            label=tool,
            color=tool_colors.get(tool),
            linewidth=lw,
            marker="o",
            markersize=6.5 if tool in top_tools else 5.8,
            markerfacecolor=tool_colors.get(tool),
            markeredgecolor="white",
            markeredgewidth=1.0 if tool in top_tools else 0.9,
            alpha=alpha,
        )

    # End labels with connectors; compute endpoints
    end_info = []  # (tool, end_x, end_y)
    for tool in tools:
        y_vals = yseries[tool]
        valid_idx = [i for i, v in enumerate(y_vals) if not _math.isnan(v)]
        if not valid_idx:
            continue
        end_x = valid_idx[-1]
        end_y = y_vals[end_x]
        end_info.append((tool, end_x, end_y))

    if end_info:
        end_info.sort(key=lambda t: t[2])
        ymin, ymax = 0.0, 100.0
        margin = 2.0
        avail = max((ymax - ymin) - 2 * margin, 1.0)
        nlab = len(end_info)
        base_sep = 2.4
        min_sep = min(base_sep, avail / max(nlab - 1, 1)) if nlab > 1 else 0

        y_adj = []
        for i, (_, _, y) in enumerate(end_info):
            if i == 0:
                y_adj.append(max(y, ymin + margin))
            else:
                y_adj.append(max(y, y_adj[i - 1] + min_sep))
        overflow = y_adj[-1] - (ymax - margin)
        if overflow > 0:
            y_adj = [yi - overflow for yi in y_adj]
            under = (ymin + margin) - min(y_adj)
            if under > 0:
                y_adj = [yi + under for yi in y_adj]

        label_x = (len(ranks) - 1) + 0.65
        star_symbol = "★"
        for (tool, ex, ey), ly in zip(end_info, y_adj):
            label_text = f"{tool} {star_symbol}" if (best_tool and tool == best_tool) else tool
            txt = ax.annotate(
                label_text,
                xy=(ex, ey),
                xytext=(label_x, ly),
                textcoords="data",
                ha="left",
                va="center",
                fontsize=10.5 if tool in top_tools else 9.8,
                color=tool_colors.get(tool),
                fontweight="semibold" if tool in top_tools else "normal",
                bbox=dict(facecolor="white", edgecolor="none", boxstyle="round,pad=0.2", alpha=0.75),
                arrowprops=dict(
                    arrowstyle="-",
                    color=tool_colors.get(tool),
                    lw=1.2 if tool in top_tools else 0.8,
                    alpha=0.6,
                    shrinkA=0,
                    shrinkB=0,
                    relpos=(0, 0.5),
                ),
            )
            txt.set_path_effects([_pe.withStroke(linewidth=2.5, foreground="white", alpha=0.9)])

    # Axes and limits
    ax.set_xticks(x)
    ax.set_xticklabels(ranks, rotation=20)
    ax.set_ylabel("Contig Accuracy (%)")
    ax.set_ylim(0, 100)
    ax.set_xlim(min(x) - 0.2, (len(ranks) - 1) + 1.25)
    clean_axis(ax)
    if ax.get_legend():
        ax.get_legend().remove()

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


def summarise_runtime(runtime_rows, bench_root: Path | None = None):
    data = defaultdict(lambda: {"cpu_sec": [], "wall_sec": [], "max_gb": []})
    samples_by_tool = defaultdict(set)
    for row in runtime_rows:
        if row.get("stage") != "run":
            continue
        tool = row.get("tool")
        if not tool:
            continue
        sample = row.get("sample") or row.get("dataset")
        if sample:
            samples_by_tool[tool].add(sample)
        try:
            user = float(row.get("user_seconds", 0.0))
        except ValueError:
            user = 0.0
        try:
            sys_time = float(row.get("sys_seconds", 0.0))
        except ValueError:
            sys_time = 0.0
        try:
            wall = float(row.get("wall_seconds", 0.0))
        except ValueError:
            wall = 0.0
        try:
            rss = float(row.get("max_rss_gb", 0.0))
        except ValueError:
            rss = 0.0
        data[tool]["cpu_sec"].append(user + sys_time)
        data[tool]["wall_sec"].append(wall)
        data[tool]["max_gb"].append(rss)
    # Special handling: CAMITAX Nextflow runs can be under-reported by /usr/bin/time
    # Try to parse wall time from Nextflow logs when available, falling back to measured values.
    if bench_root is not None and samples_by_tool.get("camitax"):
        cami_out = bench_root / "out"
        parsed_wall = []
        import re
        from datetime import datetime
        ts_re = re.compile(r"^(?P<ts>[A-Za-z]{3}-\d{2} \d{2}:\d{2}:\d{2}(?:\.\d+)?) ")

        def parse_ts(s: str):
            m = ts_re.match(s)
            if not m:
                return None
            txt = m.group("ts")
            for fmt in ("%b-%d %H:%M:%S.%f", "%b-%d %H:%M:%S"):
                try:
                    return datetime.strptime(txt, fmt)
                except Exception:
                    pass
            return None

        for sample in sorted(samples_by_tool["camitax"]):
            log_path = cami_out / sample / "camitax" / "run" / ".nextflow.log"
            if not log_path.is_file():
                continue
            try:
                with log_path.open("r", encoding="utf-8", errors="ignore") as fh:
                    first_ts = None
                    last_ts = None
                    for line in fh:
                        if "Session start" in line or "Session start" in line:
                            ts = parse_ts(line)
                            if ts and first_ts is None:
                                first_ts = ts
                        if "Execution complete" in line or "Goodbye" in line:
                            ts = parse_ts(line)
                            if ts:
                                last_ts = ts
                    if first_ts and last_ts and last_ts >= first_ts:
                        parsed_wall.append((last_ts - first_ts).total_seconds())
            except Exception:
                pass
        if parsed_wall:
            # Replace measured wall times for camitax with parsed ones (in seconds)
            data["camitax"]["wall_sec"] = parsed_wall

    summary = {}
    for tool, metrics in data.items():
        cpu_vals = [v for v in metrics["cpu_sec"] if v >= 0]
        wall_vals = [v for v in metrics["wall_sec"] if v >= 0]
        mem_vals = [v for v in metrics["max_gb"] if v >= 0]
        summary[tool] = {
            "cpu_min": mean(cpu_vals) / 60.0 if cpu_vals else 0.0,
            "wall_min": mean(wall_vals) / 60.0 if wall_vals else 0.0,
            "max_gb": max(mem_vals) if mem_vals else 0.0,
        }
    return summary


def plot_runtime(summary, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt
    import math as _math

    def _fmt_minutes(m):
        """Human-friendly time from minutes.

        Rules:
        - < 90s: show seconds ("Xs")
        - < 10m: show one decimal minute ("X.Ym")
        - < 60m: whole minutes ("Xm")
        - < 24h: hours with optional minutes ("Xh Ym")
        - >= 24h: days with optional hours ("Xd Yh")
        """
        try:
            m = float(m)
        except Exception:
            return "NA"
        if m <= 0:
            return "NA"
        s = m * 60.0
        if s < 90:
            # Avoid showing "0s" for tiny non-zero values
            return f"{max(1, int(round(s)))}s"
        if m < 10:
            return f"{m:.1f}m"
        if m < 60:
            return f"{int(round(m))}m"
        h = m / 60.0
        if h < 24:
            whole = int(h)
            rem_m = int(round((h - whole) * 60))
            return f"{whole}h" + (f" {rem_m}m" if rem_m and whole < 8 else "")
        d = h / 24.0
        whole = int(d)
        rem_h = int(round((d - whole) * 24))
        return f"{whole}d" + (f" {rem_h}h" if rem_h and whole < 5 else "")

    # Use all known tools (color key) so the figure stays consistent
    tools_all = order_tools(list(tool_colors.keys() or summary.keys()))
    if not tools_all:
        return

    # Build list and sort by ascending wall-clock time (faster tools first)
    items = []  # (tool, value or None)
    for t in tools_all:
        v = summary.get(t, {}).get("wall_min")
        if v is None or v <= 0:
            items.append((t, None))
        else:
            items.append((t, float(v)))

    # Separate available and missing, sort available ascending (lower is better)
    avail = [(t, v) for t, v in items if v is not None]
    missing = [(t, v) for t, v in items if v is None]
    avail.sort(key=lambda kv: kv[1])
    order = avail + missing
    tools = [t for t, _ in order]
    values = [v for _, v in order]

    # Prepare figure: horizontal lollipop on log scale
    n = len(tools)
    fig_h = max(4.5, min(0.5 * n + 1.6, 9.5))
    fig, ax = plt.subplots(figsize=(10.5, fig_h))

    y = list(range(n))
    y.reverse()  # top = rank 1 (fastest)
    tools = list(reversed(tools))
    values = list(reversed(values))

    xmin_candidates = [v for v in values if v is not None and v > 0]
    if not xmin_candidates:
        return
    xmin = min(xmin_candidates)
    xmax = max(xmin_candidates)

    # Draw stems and markers
    for yi, (t, v) in enumerate(zip(tools, values)):
        color = tool_colors.get(t)
        if v is None:
            # Place a subtle NA label to the left of the y-label row
            ax.text(
                xmax * 0.7 if xmax > 0 else 1.0,
                y[yi],
                "NA",
                va="center",
                ha="right",
                fontsize=9,
                color="#9ca3af",
            )
            continue
        # Stem from xmin reference to value for visual scale perception
        ax.hlines(y[yi], xmin, v, color="#e0e7ef", linewidth=2.0)
        ax.plot(
            v,
            y[yi],
            marker="o",
            markersize=9,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=1.2,
        )
        # Value label with humanized units
        ax.text(
            v,
            y[yi],
            f"  {_fmt_minutes(v)}",
            va="center",
            ha="left",
            fontsize=10,
            color="#374151",
            fontweight="medium",
        )

    # Aesthetics and scales
    ax.set_yticks(y)
    ax.set_yticklabels(tools)
    ax.set_xlabel("Wall time (minutes, log scale)")
    ax.set_xscale("log")
    # Pad a bit on both ends for labels
    left = xmin / 1.6 if xmin > 0 else 0.1
    right = xmax * 1.35 if xmax > 0 else 10
    ax.set_xlim(left, right)
    clean_axis(ax)
    # Emphasize that lower is better
    ax.set_title("Wall-clock Time by Tool (lower is better)", pad=10)
    # Light grid on x only for readability
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    # Remove legend; color already encodes tool
    if ax.get_legend():
        ax.get_legend().remove()

    fig.tight_layout(pad=0.6)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_memory(summary, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    def _fmt_gb(g):
        try:
            g = float(g)
        except Exception:
            return "NA"
        if g <= 0:
            return "NA"
        if g < 1:
            return f"{g*1024:.0f} MB"
        if g < 10:
            return f"{g:.1f} GB"
        return f"{g:.0f} GB"

    tools_all = order_tools(list(tool_colors.keys() or summary.keys()))
    if not tools_all:
        return

    items = []
    for t in tools_all:
        v = summary.get(t, {}).get("max_gb")
        if v is None or v <= 0:
            items.append((t, None))
        else:
            items.append((t, float(v)))

    avail = [(t, v) for t, v in items if v is not None]
    missing = [(t, v) for t, v in items if v is None]
    # Sort ascending (lower memory usage first)
    avail.sort(key=lambda kv: kv[1])
    order = avail + missing
    tools = [t for t, _ in order]
    values = [v for _, v in order]

    n = len(tools)
    fig_h = max(4.5, min(0.5 * n + 1.6, 9.5))
    fig, ax = plt.subplots(figsize=(10.5, fig_h))

    y = list(range(n))
    y.reverse()
    tools = list(reversed(tools))
    values = list(reversed(values))

    xmin_candidates = [v for v in values if v is not None and v > 0]
    if not xmin_candidates:
        return
    xmin = min(xmin_candidates)
    xmax = max(xmin_candidates)

    for yi, (t, v) in enumerate(zip(tools, values)):
        color = tool_colors.get(t)
        if v is None:
            ax.text(
                xmax * 0.7 if xmax > 0 else 1.0,
                y[yi],
                "NA",
                va="center",
                ha="right",
                fontsize=9,
                color="#9ca3af",
            )
            continue
        ax.hlines(y[yi], xmin, v, color="#e0e7ef", linewidth=2.0)
        ax.plot(
            v,
            y[yi],
            marker="o",
            markersize=9,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=1.2,
        )
        ax.text(
            v,
            y[yi],
            f"  {_fmt_gb(v)}",
            va="center",
            ha="left",
            fontsize=10,
            color="#374151",
            fontweight="medium",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(tools)
    ax.set_xlabel("Peak RSS (GB, log scale)")
    ax.set_xscale("log")
    left = xmin / 1.6 if xmin > 0 else 0.1
    right = xmax * 1.35 if xmax > 0 else 10
    ax.set_xlim(left, right)
    clean_axis(ax)
    ax.set_title("Peak Memory by Tool (lower is better)", pad=10)
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    if ax.get_legend():
        ax.get_legend().remove()

    fig.tight_layout(pad=0.6)
    fig.savefig(out_path, bbox_inches="tight")
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

    runtime_summary = summarise_runtime(runtime_rows, bench_root)
    tools = sorted({row["tool"] for row in summary_rows})
    if contig_rows:
        tools = sorted(set(tools).union({row["tool"] for row in contig_rows}))
    if runtime_summary:
        tools = sorted(set(tools).union(runtime_summary.keys()))
    tool_colors = get_tool_colors(tools)

    plot_f1_by_rank(summary_rows, out_root / "fig_f1_by_rank.png", tool_colors)
    # Alternative, easier-to-compare line chart version
    plot_f1_by_rank_lines(summary_rows, out_root / "fig_f1_by_rank_lines.png", tool_colors)
    plot_l1_bray(summary_rows, out_root / "fig_l1_braycurtis.png", tool_colors)
    plot_l1_bray_lines(summary_rows, out_root / "fig_l1_braycurtis_lines.png", tool_colors)
    if contig_rows:
        plot_accuracy(contig_rows, out_root / "fig_accuracy_by_rank.png", tool_colors)
        # Line version with end labels and highlights
        plot_accuracy_by_rank_lines(contig_rows, out_root / "fig_accuracy_by_rank_lines.png", tool_colors)
    plot_per_sample_stack(summary_rows, out_root / "fig_per_sample_f1_stack.png", tool_colors)

    if runtime_summary:
        # Keep legacy filename for continuity
        plot_runtime(runtime_summary, out_root / "fig_cpu_time_by_tool.png", tool_colors)
        # Also emit an explicit wall-time figure for clarity
        plot_runtime(runtime_summary, out_root / "fig_wall_time_by_tool.png", tool_colors)
        plot_memory(runtime_summary, out_root / "fig_peak_memory_by_tool.png", tool_colors)


if __name__ == "__main__":
    main()
