#!/usr/bin/env python3
"""Create benchmark plots from aggregated CAMI metrics."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from collections import Counter, defaultdict
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
    "#0f4c81",  # deep sapphire
    "#4d908e",  # muted teal
    "#c8b6ff",  # lavender bloom
    "#ff6b6b",  # coral rose
    "#ffd166",  # golden hour
    "#5f0f40",  # mulberry
    "#277da1",  # blue lagoon
    "#90be6d",  # fresh fern
    "#f28482",  # salmon blush
    "#118ab2",  # refined cyan
    "#ef709d",  # modern magenta
    "#7b6ef6",  # periwinkle pop
]

# Whether an emoji-capable font was detected and enabled
HAS_EMOJI_FONT = False

TOOL_NAME_OVERRIDES = {
    "basta": "BASTA",
    "camitax": "CAMITAX",
    "centrifuge": "Centrifuge",
    "ganon2": "Ganon 2",
    "hymet": "HYMET",
    "kraken2": "Kraken 2",
    "megapath_nano": "MegaPath-Nano",
    "metaphlan4": "MetaPhlAn 4",
    "phabox": "PHABOX",
    "phyloflash": "phyloFlash",
    "snakemags": "SnakeMAGs",
    "sourmash_gather": "sourmash gather",
    "squeezemeta": "SqueezeMeta",
    "tama": "TAMA",
    "viwrap": "ViWrap",
}


def soften_color(color, mix: float = 0.18):
    """Blend a color toward white for softer fills."""
    from matplotlib.colors import to_rgb

    if not color:
        color = "#4b5563"
    try:
        r, g, b = to_rgb(color)
    except ValueError:
        r, g, b = to_rgb("#4b5563")
    mix = max(0.0, min(1.0, mix))
    return tuple(r + (1.0 - r) * mix for r in (r, g, b))


def format_tool_label(tool: str) -> str:
    return TOOL_NAME_OVERRIDES.get(tool, tool.replace("_", " ").title())


def format_seconds(value, include_long: bool = False) -> str:
    try:
        seconds = float(value)
    except Exception:
        return "NA"
    if seconds <= 0:
        return "NA"
    base = f"{seconds:.0f} s" if seconds >= 1 else f"{seconds:.2f} s"
    if not include_long or seconds < 60:
        return base
    minutes = seconds / 60.0
    if minutes < 60:
        return f"{base} ({minutes:.2f} min)"
    hours = minutes / 60.0
    return f"{base} ({hours:.2f} h)"


def format_gib(value) -> str:
    try:
        gib = float(value)
    except Exception:
        return "NA"
    if gib <= 0:
        return "NA"
    if gib < 1:
        return f"{gib:.2f} GiB"
    if gib < 10:
        return f"{gib:.2f} GiB"
    return f"{gib:.1f} GiB"


def env_flag(name: str, default: bool = False) -> bool:
    val = os.environ.get(name)
    if val is None:
        return default
    val = val.strip().lower()
    if not val:
        return default
    return val in {"1", "true", "yes", "on"}


def percentile(values, q):
    if not values:
        return 0.0
    arr = sorted(values)
    idx = (len(arr) - 1) * q
    lower = int(idx)
    upper = min(lower + 1, len(arr) - 1)
    weight = idx - lower
    return arr[lower] * (1 - weight) + arr[upper] * weight


def describe_series(values):
    clean = [v for v in values if v and v > 0]
    if not clean:
        return {}
    clean.sort()
    mid = clean[len(clean) // 2] if len(clean) % 2 else (clean[len(clean) // 2 - 1] + clean[len(clean) // 2]) / 2
    return {
        "count": len(clean),
        "min": clean[0],
        "max": clean[-1],
        "median": mid,
        "p10": percentile(clean, 0.10),
        "p90": percentile(clean, 0.90),
    }


def parse_threads(command: str | None):
    if not command:
        return None
    match = re.search(r"--threads(?:\s+|=)(\d+)", command)
    if match:
        try:
            return int(match.group(1))
        except Exception:
            return None
    return None


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
    threads_by_tool = defaultdict(list)
    all_samples = set()
    for row in runtime_rows:
        if row.get("stage") != "run":
            continue
        tool = row.get("tool")
        if not tool:
            continue
        sample = row.get("sample") or row.get("dataset")
        if sample:
            samples_by_tool[tool].add(sample)
            all_samples.add(sample)
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
        threads = parse_threads(row.get("command"))
        if threads:
            threads_by_tool[tool].append(threads)
    # Special handling: CAMITAX Nextflow runs can be under-reported by /usr/bin/time
    # Try to parse wall time from Nextflow logs when available, falling back to measured values.
    if bench_root is not None and samples_by_tool.get("camitax"):
        cami_out = bench_root / "out"
        parsed_wall = []
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

    thread_meta = {}
    for tool, values in threads_by_tool.items():
        if values:
            thread_meta[tool] = Counter(values).most_common(1)[0][0]

    summary = {}
    for tool, metrics in data.items():
        cpu_vals = [v for v in metrics["cpu_sec"] if v > 0]
        wall_vals = [v for v in metrics["wall_sec"] if v > 0]
        mem_vals = [v for v in metrics["max_gb"] if v > 0]
        cpu_stats = describe_series(cpu_vals)
        wall_stats = describe_series(wall_vals)
        mem_stats = describe_series(mem_vals)
        summary[tool] = {
            "cpu_min": mean(cpu_vals) / 60.0 if cpu_vals else 0.0,
            "wall_min": mean(wall_vals) / 60.0 if wall_vals else 0.0,
            "max_gb": max(mem_vals) if mem_vals else 0.0,
            "cpu_sec_stats": cpu_stats,
            "wall_sec_stats": wall_stats,
            "mem_gb_stats": mem_stats,
            "sample_count": len(samples_by_tool.get(tool, [])),
            "threads": thread_meta.get(tool),
        }

    meta = {
        "thread_counts": sorted({v for v in thread_meta.values() if v}),
        "total_samples": sorted(all_samples),
    }
    return summary, meta



def plot_runtime(
    summary,
    out_path: Path,
    tool_colors,
    *,
    metric_key: str = 'wall_min',
    stats_key: str = 'wall_sec_stats',
    title: str = 'Runtime budget by tool',
    xlabel: str = 'Runtime per CAMI sample (seconds, log10 scale)',
    subtitle: str = 'Lower is better - aggregated across CAMI samples',
    runtime_meta: dict | None = None,
    show_notes: bool = False,
):
    import matplotlib.pyplot as plt
    from matplotlib import patheffects as _pe
    from matplotlib.ticker import FuncFormatter, LogLocator

    runtime_meta = runtime_meta or {}

    def _tick_label(value, _pos=None):
        label = format_seconds(value)
        return '' if label == 'NA' else label

    tools_all = order_tools(list(tool_colors.keys() or summary.keys()))
    if not tools_all:
        return

    rows = []
    missing_labels = []
    for tool in tools_all:
        metrics = summary.get(tool)
        label = format_tool_label(tool)
        base_entry = {
            'tool': tool,
            'label': label,
            'value': None,
            'stats': {},
            'color': tool_colors.get(tool, '#4b5563'),
            'sample_count': 0,
        }
        if not metrics:
            missing_labels.append(label)
            rows.append(base_entry)
            continue
        base_entry['sample_count'] = metrics.get('sample_count', 0)
        value = metrics.get(metric_key)
        stats = metrics.get(stats_key) or {}
        if value is None or value <= 0:
            missing_labels.append(label)
            rows.append(base_entry)
            continue
        if metric_key.endswith('_min'):
            value = float(value) * 60.0
        else:
            value = float(value)
        base_entry['value'] = value
        base_entry['stats'] = stats
        rows.append(base_entry)

    if not rows:
        return

    rows.sort(key=lambda r: r['value'] if r['value'] is not None else float('inf'))

    n = len(rows)
    fig_h = max(4.8, min(0.58 * n + 1.6, 10.5))
    fig, ax = plt.subplots(figsize=(11.8, fig_h))
    fig.patch.set_facecolor('#ffffff')
    ax.set_facecolor('#fdfdff')

    values = [r['value'] for r in rows if r['value'] is not None]
    if values:
        xmin = min(values)
        xmax = max(values)
        left = max(xmin / 1.8 if xmin > 0 else 0.5, 0.3)
        right = xmax * 1.8 if xmax > 0 else 10.0
    else:
        left, right = 0.5, 10.0

    y_positions = list(range(n))
    ax.set_yticks(y_positions)
    ax.set_yticklabels([])
    ax.tick_params(axis='y', length=0)
    y_label_transform = ax.get_yaxis_transform()

    for idx, row in enumerate(rows):
        y = idx
        base_color = row['color']
        value = row['value']
        if value is not None:
            face = soften_color(base_color, 0.25)
            ax.barh(
                y,
                value,
                height=0.52,
                color=face,
                edgecolor=base_color,
                linewidth=1.0,
                zorder=3,
            )
        ax.text(
            -0.03,
            y,
            row['label'],
            transform=y_label_transform,
            ha='right',
            va='center',
            fontsize=11,
            color='#0f172a',
            zorder=6,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.9, pad=0.2),
        )
        if value is None:
            ax.text(
                left * 1.05,
                y,
                'not reported',
                va='center',
                ha='left',
                fontsize=10,
                color='#9ca3af',
            )
            continue
        stats = row['stats']
        if stats.get('p10') and stats.get('p90') and stats['p90'] > stats['p10']:
            ax.hlines(
                y,
                stats['p10'],
                stats['p90'],
                color=base_color,
                linewidth=3.2,
                alpha=0.75,
                zorder=4,
            )
        min_label_x = left * 1.12
        text_x = max(value * 1.15, min_label_x)
        text_x = min(text_x, right / 1.03)
        ax.text(
            text_x,
            y,
            format_seconds(value, include_long=True),
            va='center',
            ha='left',
            fontsize=10.2,
            color='#0f172a',
            zorder=6,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=0.15),
        ).set_path_effects([_pe.withStroke(linewidth=1.6, foreground='white', alpha=0.9)])

    ax.set_xscale('log')
    ax.set_xlim(left, right)
    ax.grid(axis='x', linestyle='--', linewidth=0.6, alpha=0.4, color='#d5dceb')
    ax.tick_params(axis='both', labelsize=10)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_color('#d3d8e3')
    ax.spines['bottom'].set_color('#d3d8e3')
    ax.xaxis.set_major_locator(LogLocator(base=10.0))
    ax.xaxis.set_major_formatter(FuncFormatter(_tick_label))
    ax.set_xlabel(xlabel)
    if show_notes:
        notes = []
        total_samples = runtime_meta.get('total_samples') or []
        if total_samples:
            notes.append(f"Benchmark covers {len(total_samples)} CAMI samples.")
        sample_counts = sorted({r['sample_count'] for r in rows if r['sample_count']})
        if sample_counts:
            if len(sample_counts) == 1:
                notes.append(f"Per-tool sample count: {sample_counts[0]}.")
            else:
                notes.append(f"Per-tool sample count range: {sample_counts[0]}–{sample_counts[-1]}.")
        thread_counts = runtime_meta.get('thread_counts') or []
        if thread_counts:
            notes.append(f"Bench scripts invoked with --threads={', '.join(str(t) for t in thread_counts)}.")
        else:
            notes.append('Thread counts were not recorded in runtime metadata.')
        notes.append('Hardware model and RAM details were not captured in RUN_0 metadata.')
        notes.append('Bar colors follow the same tool palette across all runtime/memory plots.')
        notes.append('Bars show mean runtime per CAMI sample; whiskers span the 10th–90th percentile across samples.')
        if missing_labels:
            notes.append('Not reported: ' + ', '.join(missing_labels))

        fig.text(
            0.01,
            0.02,
            "\n".join(notes),
            ha='left',
            fontsize=9,
            color='#4b5563',
        )

    fig.subplots_adjust(left=0.36, right=0.96, top=0.85, bottom=0.2)
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)


def plot_memory(summary, out_path: Path, tool_colors, runtime_meta: dict | None = None, show_notes: bool = False):
    import matplotlib.pyplot as plt
    from matplotlib import patheffects as _pe
    from matplotlib.ticker import FuncFormatter, LogLocator

    runtime_meta = runtime_meta or {}

    def _tick_label(value, _pos=None):
        label = format_gib(value)
        return '' if label == 'NA' else label

    tools_all = order_tools(list(tool_colors.keys() or summary.keys()))
    if not tools_all:
        return

    rows = []
    missing_labels = []
    for tool in tools_all:
        metrics = summary.get(tool)
        label = format_tool_label(tool)
        base_entry = {
            'tool': tool,
            'label': label,
            'value': None,
            'stats': {},
            'color': tool_colors.get(tool, '#4b5563'),
        }
        if not metrics:
            missing_labels.append(label)
            rows.append(base_entry)
            continue
        value = metrics.get('max_gb')
        stats = metrics.get('mem_gb_stats') or {}
        if value is None or value <= 0:
            missing_labels.append(label)
            rows.append(base_entry)
            continue
        base_entry['value'] = float(value)
        base_entry['stats'] = stats
        rows.append(base_entry)

    if not rows:
        return

    rows.sort(key=lambda r: r['value'] if r['value'] is not None else float('inf'))
    n = len(rows)
    fig_h = max(4.8, min(0.58 * n + 1.6, 10.5))
    fig, ax = plt.subplots(figsize=(11.5, fig_h))
    fig.patch.set_facecolor('#ffffff')
    ax.set_facecolor('#fdfdff')

    values = [r['value'] for r in rows if r['value'] is not None]
    if values:
        xmin = min(values)
        xmax = max(values)
        left = max(xmin / 1.9 if xmin > 0 else 0.05, 0.02)
        right = xmax * 1.85 if xmax > 0 else 10.0
    else:
        left, right = 0.05, 10.0

    y_positions = list(range(n))
    ax.set_yticks(y_positions)
    ax.set_yticklabels([])
    ax.tick_params(axis='y', length=0)
    y_label_transform = ax.get_yaxis_transform()

    for idx, row in enumerate(rows):
        y = idx
        base_color = row['color']
        value = row['value']
        if value is not None:
            face = soften_color(base_color, 0.28)
            ax.barh(
                y,
                value,
                height=0.52,
                color=face,
                edgecolor=base_color,
                linewidth=1.0,
                zorder=3,
            )
        ax.text(
            -0.03,
            y,
            row['label'],
            transform=y_label_transform,
            ha='right',
            va='center',
            fontsize=11,
            color='#0f172a',
            zorder=6,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.9, pad=0.2),
        )
        if value is None:
            ax.text(
                left * 1.05,
                y,
                'not reported',
                va='center',
                ha='left',
                fontsize=10,
                color='#9ca3af',
            )
            continue
        stats = row['stats']
        if stats.get('p10') and stats.get('p90') and stats['p90'] > stats['p10']:
            ax.hlines(
                y,
                stats['p10'],
                stats['p90'],
                color=base_color,
                linewidth=3.0,
                alpha=0.75,
                zorder=4,
            )
        min_label_x = left * 1.12
        text_x = max(value * 1.15, min_label_x)
        text_x = min(text_x, right / 1.03)
        ax.text(
            text_x,
            y,
            format_gib(value),
            va='center',
            ha='left',
            fontsize=10.2,
            color='#0f172a',
            zorder=6,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=0.15),
        ).set_path_effects([_pe.withStroke(linewidth=1.6, foreground='white', alpha=0.9)])

    ax.set_xscale('log')
    ax.set_xlim(left, right)
    ax.grid(axis='x', linestyle='--', linewidth=0.6, alpha=0.4, color='#d5dceb')
    ax.tick_params(axis='both', labelsize=10)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_color('#d3d8e3')
    ax.spines['bottom'].set_color('#d3d8e3')
    ax.xaxis.set_major_locator(LogLocator(base=10.0))
    ax.xaxis.set_major_formatter(FuncFormatter(_tick_label))
    ax.set_xlabel('Peak resident memory (GiB, log10 scale)')
    # Titles removed per styling guidance; rely on axis labels and footnotes.

    if show_notes:
        notes = [
            'Bars show the maximum RSS observed per CAMI sample; whiskers mark the 10th–90th percentile across samples.'
        ]
        notes.append('Bar colors follow the same tool palette across all runtime/memory plots.')
        if runtime_meta.get('total_samples'):
            notes.append(f"Benchmark covers {len(runtime_meta['total_samples'])} CAMI samples.")
        thread_counts = runtime_meta.get('thread_counts') or []
        if thread_counts:
            notes.append(f"Benchmark threads per tool: {', '.join(str(t) for t in thread_counts)}.")
        else:
            notes.append('Thread counts were not recorded in runtime metadata.')
        notes.append('Hardware model and RAM details were not captured in RUN_0 metadata.')
        if missing_labels:
            notes.append('Not reported: ' + ', '.join(missing_labels))

        fig.text(0.01, 0.02, "\n".join(notes), ha='left', fontsize=9, color='#4b5563')

    fig.subplots_adjust(left=0.36, right=0.96, top=0.85, bottom=0.2)
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)

def plot_runtime_memory_panels(summary, out_path: Path, tool_colors, runtime_meta: dict | None = None, show_notes: bool = False):
    import matplotlib.pyplot as plt
    from matplotlib import patheffects as _pe
    from matplotlib.ticker import FuncFormatter, LogLocator

    runtime_meta = runtime_meta or {}
    tools_all = order_tools(list(tool_colors.keys() or summary.keys()))
    if not tools_all:
        return

    rows = []
    missing_cpu = []
    missing_mem = []
    for tool in tools_all:
        metrics = summary.get(tool)
        label = format_tool_label(tool)
        if not metrics:
            missing_cpu.append(label)
            missing_mem.append(label)
            rows.append(
                {
                    'tool': tool,
                    'label': label,
                    'cpu': None,
                    'cpu_stats': {},
                    'mem': None,
                    'mem_stats': {},
                    'color': tool_colors.get(tool, '#4b5563'),
                }
            )
            continue
        cpu_val = metrics.get('cpu_min')
        if cpu_val and cpu_val > 0:
            cpu_val = float(cpu_val) * 60.0
        else:
            cpu_val = None
            missing_cpu.append(label)
        mem_val = metrics.get('max_gb')
        if mem_val and mem_val > 0:
            mem_val = float(mem_val)
        else:
            mem_val = None
            missing_mem.append(label)
        rows.append(
            {
                'tool': tool,
                'label': label,
                'cpu': cpu_val,
                'cpu_stats': metrics.get('cpu_sec_stats') or {},
                'mem': mem_val,
                'mem_stats': metrics.get('mem_gb_stats') or {},
                'color': tool_colors.get(tool, '#4b5563'),
            }
        )

    if not rows:
        return

    rows.sort(key=lambda r: (float('inf') if r['cpu'] is None else r['cpu']))
    n = len(rows)
    fig_h = max(4.6, min(0.45 * n + 2.4, 9.5))
    fig, (ax_cpu, ax_mem) = plt.subplots(1, 2, sharey=True, figsize=(13.5, fig_h))
    fig.patch.set_facecolor('#ffffff')
    for ax in (ax_cpu, ax_mem):
        ax.set_facecolor('#fdfdff')

    y_positions = list(range(n))
    ax_cpu.set_yticks(y_positions)
    ax_cpu.set_yticklabels([])
    ax_cpu.tick_params(axis='y', length=0)
    ax_mem.set_yticks(y_positions)
    ax_mem.set_yticklabels([])
    y_label_transform = ax_cpu.get_yaxis_transform()
    for idx, row in enumerate(rows):
        ax_cpu.text(
            -0.06,
            idx,
            row['label'],
            transform=y_label_transform,
            ha='right',
            va='center',
            fontsize=10.5,
            color='#0f172a',
            zorder=6,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.9, pad=0.2),
        )

    def _draw_axis(ax, values_key, stats_key, formatter, xlabel, left_pad, right_pad):
        values = [r[values_key] for r in rows if r[values_key]]
        if values:
            xmin = min(values)
            xmax = max(values)
            left = max(xmin / left_pad if xmin > 0 else 0.05, 0.02)
            right = xmax * right_pad if xmax > 0 else 10.0
        else:
            left, right = 0.5, 1.0
        ax.set_xscale('log')
        ax.set_xlim(left, right)
        ax.grid(axis='x', linestyle='--', linewidth=0.55, alpha=0.35, color='#d5dceb')
        ax.tick_params(axis='x', labelsize=9)
        ax.xaxis.set_major_locator(LogLocator(base=10.0))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda val, pos=None: formatter(val)))
        for idx, row in enumerate(rows):
            y = idx
            value = row[values_key]
            if not value:
                ax.text(left * 1.05, y, 'not reported', va='center', ha='left', fontsize=9, color='#9ca3af')
                continue
            base_color = row['color']
            face = soften_color(base_color, 0.25 if values_key == 'cpu' else 0.3)
            ax.barh(y, value, height=0.45, color=face, edgecolor=base_color, linewidth=0.9, zorder=3)
            stats = row[stats_key]
            if stats.get('p10') and stats.get('p90') and stats['p90'] > stats['p10']:
                ax.hlines(y, stats['p10'], stats['p90'], color=base_color, linewidth=2.4, alpha=0.75, zorder=4)
            min_label_x = left * 1.12
            text_x = max(value * 1.12, min_label_x)
            text_x = min(text_x, right / 1.03)
            align = 'left'
            ax.text(
                text_x,
                y,
                formatter(value, include_long=True) if values_key == 'cpu' else format_gib(value),
                va='center',
                ha=align,
                fontsize=9.5,
                color='#111827',
                zorder=6,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=0.1),
            ).set_path_effects([_pe.withStroke(linewidth=1.4, foreground='white', alpha=0.9)])
        ax.set_xlabel(xlabel)
        ax.spines['left'].set_color('#d3d8e3')
        ax.spines['bottom'].set_color('#d3d8e3')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    _draw_axis(
        ax_cpu,
        'cpu',
        'cpu_stats',
        lambda val, include_long=False: format_seconds(val) if not include_long else format_seconds(val, include_long=True),
        'CPU time per sample (seconds, log10 scale)',
        left_pad=1.8,
        right_pad=1.7,
    )
    _draw_axis(
        ax_mem,
        'mem',
        'mem_stats',
        lambda val, include_long=False: format_gib(val),
        'Peak memory per sample (GiB, log10 scale)',
        left_pad=1.9,
        right_pad=1.75,
    )

    ax_mem.tick_params(axis='y', left=False, labelleft=False, right=False)

    if show_notes:
        notes = []
        if runtime_meta.get('total_samples'):
            notes.append(f"Data aggregated across {len(runtime_meta['total_samples'])} CAMI samples.")
        thread_counts = runtime_meta.get('thread_counts') or []
        if thread_counts:
            notes.append(f"Bench scripts ran with --threads={', '.join(str(t) for t in thread_counts)}.")
        else:
            notes.append('Thread counts were not recorded in runtime metadata.')
        notes.append('Hardware model and RAM details were not captured in RUN_0 metadata.')
        if missing_cpu:
            notes.append('CPU not reported: ' + ', '.join(sorted(set(missing_cpu))))
        if missing_mem:
            notes.append('Memory not reported: ' + ', '.join(sorted(set(missing_mem))))

        fig.text(0.01, 0.02, "\n".join(notes), ha='left', fontsize=9, color='#4b5563')
    fig.subplots_adjust(left=0.36, right=0.98, top=0.9, bottom=0.2, wspace=0.18)
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)

def main() -> None:
    ap = argparse.ArgumentParser(description="Generate plots from CAMI aggregate metrics.")
    ap.add_argument("--bench-root", default=str(Path(__file__).resolve().parent.parent), help="Bench directory root.")
    ap.add_argument("--outdir", default="out", help="Relative or absolute output directory for figures.")
    ap.add_argument(
        "--tables",
        default=None,
        help="Directory containing summary_per_tool_per_sample.tsv, runtime_memory.tsv, etc. "
        "If not provided, defaults to the figure output directory (previous behaviour).",
    )
    args = ap.parse_args()

    bench_root = Path(args.bench_root)
    out_root = Path(args.outdir)
    if not out_root.is_absolute():
        out_root = bench_root / out_root
    out_root.mkdir(parents=True, exist_ok=True)

    tables_root = Path(args.tables) if args.tables else out_root
    if not tables_root.is_absolute():
        tables_root = bench_root / tables_root
    if not tables_root.exists():
        raise FileNotFoundError(f"Tables directory not found: {tables_root}")

    summary_path = tables_root / "summary_per_tool_per_sample.tsv"
    contig_path = tables_root / "contig_accuracy_per_tool.tsv"

    summary_rows = load_rows(summary_path)
    contig_rows = load_rows(contig_path)
    runtime_rows = load_runtime_rows(tables_root / "runtime_memory.tsv")
    if not summary_rows:
        print("[plot] No summary data available; skipping figure generation.")
        return

    ensure_matplotlib()

    runtime_summary, runtime_meta = summarise_runtime(runtime_rows, bench_root)
    show_notes = env_flag("HYMET_PLOTS_SHOW_NOTES", False)
    tools = sorted({row["tool"] for row in summary_rows})
    if contig_rows:
        tools = sorted(set(tools).union({row["tool"] for row in contig_rows}))
    if runtime_summary:
        tools = sorted(set(tools).union(runtime_summary.keys()))
    tool_colors = get_tool_colors(order_tools(tools))

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
        plot_runtime(
            runtime_summary,
            out_root / "fig_cpu_time_by_tool.png",
            tool_colors,
            metric_key="cpu_min",
            stats_key="cpu_sec_stats",
            title="CPU time by tool",
            xlabel="Runtime per CAMI sample (seconds, log10 scale)",
            subtitle="User + system CPU seconds aggregated across CAMI samples - lower is better",
            runtime_meta=runtime_meta,
            show_notes=show_notes,
        )
        # Also emit an explicit wall-time figure for clarity
        plot_runtime(
            runtime_summary,
            out_root / "fig_wall_time_by_tool.png",
            tool_colors,
            metric_key="wall_min",
            stats_key="wall_sec_stats",
            title="Wall-clock time by tool",
            xlabel="Wall time per CAMI sample (seconds, log10 scale)",
            subtitle="Lower is better - aggregated across CAMI samples",
            runtime_meta=runtime_meta,
            show_notes=show_notes,
        )
        plot_memory(
            runtime_summary,
            out_root / "fig_peak_memory_by_tool.png",
            tool_colors,
            runtime_meta=runtime_meta,
            show_notes=show_notes,
        )
        plot_runtime_memory_panels(
            runtime_summary,
            out_root / "fig_runtime_cpu_mem.png",
            tool_colors,
            runtime_meta=runtime_meta,
            show_notes=show_notes,
        )


if __name__ == "__main__":
    main()
