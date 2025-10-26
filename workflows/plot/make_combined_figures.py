#!/usr/bin/env python3
"""Generate combined contig/read figures for CAMI suite runs."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

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
    "contigs": "#1f77b4",
    "reads": "#f97316",
}


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
            "axes.facecolor": "#f7f9fb",
            "figure.facecolor": "white",
            "axes.edgecolor": "#d0d5dd",
            "grid.alpha": 0.25,
            "grid.linestyle": "--",
            "axes.titleweight": "semibold",
            "axes.labelweight": "semibold",
            "font.size": 11,
        }
    )
    return matplotlib.pyplot


def order_tools(tools: Iterable[str]) -> List[str]:
    seen = []
    remaining = [tool for tool in tools if tool]
    for tool in PREFERRED_TOOL_ORDER:
        if tool in remaining:
            seen.append(tool)
            remaining.remove(tool)
    seen.extend(sorted(remaining))
    return seen


def collect_runtime(runtime_rows: List[Dict[str, str]]) -> Tuple[Dict[Tuple[str, str], float], List[str], List[str]]:
    totals = defaultdict(float)
    modes = set()
    tools = set()
    for row in runtime_rows:
        stage = (row.get("stage") or "").strip()
        if stage and stage not in {"run", "eval"}:
            continue
        if stage == "eval":
            continue
        tool = (row.get("tool") or "").strip() or "unknown"
        mode = (row.get("mode") or "").strip() or "unknown"
        wall_seconds = safe_float(row.get("wall_seconds"))
        totals[(tool, mode)] += wall_seconds / 3600.0  # hours
        tools.add(tool)
        modes.add(mode)
    return totals, order_tools(tools), sorted(modes)


def collect_genus_f1(summary_rows: List[Dict[str, str]]) -> Tuple[Dict[Tuple[str, str], float], List[str], List[str]]:
    values = defaultdict(list)
    modes = set()
    tools = set()
    for row in summary_rows:
        rank = (row.get("rank") or "").strip().lower()
        if rank and rank != "genus":
            continue
        tool = (row.get("tool") or "").strip() or "unknown"
        mode = (row.get("mode") or "").strip() or "unknown"
        score = safe_float(row.get("F1_%"))
        if score == 0.0:
            continue
        values[(tool, mode)].append(score / 100.0)
        tools.add(tool)
        modes.add(mode)
    averages = {key: (sum(vals) / len(vals)) for key, vals in values.items()}
    return averages, order_tools(tools), sorted(modes)


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


def plot_grouped_bars(data: Dict[Tuple[str, str], float], tools: List[str], modes: List[str], ylabel: str, title: str, out_path: Path):
    if not data or not tools or not modes:
        return False
    plt = ensure_matplotlib()
    fig, ax = plt.subplots(figsize=(max(6, len(tools) * 0.7), 5.0))
    width = min(0.25, 0.8 / max(len(modes), 1))
    offsets = [width * (i - (len(modes) - 1) / 2.0) for i in range(len(modes))]
    x_positions = range(len(tools))

    for idx, mode in enumerate(modes):
        heights = [data.get((tool, mode), 0.0) for tool in tools]
        ax.bar(
            [pos + offsets[idx] for pos in x_positions],
            heights,
            width=width,
            color=resolve_color(mode, idx),
            label=mode.capitalize(),
        )

    ax.set_xticks(list(x_positions))
    ax.set_xticklabels(tools, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.set_ylim(bottom=0)
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

    runtime_rows = load_rows(tables_dir / "runtime_memory.tsv")
    summary_rows = load_rows(combined_dir / "summary_per_tool_per_sample.tsv")

    generated_any = False
    title_suffix = args.suite or Path(args.run).name or "CAMI suite"

    runtime_data, runtime_tools, runtime_modes = collect_runtime(runtime_rows)
    runtime_path = outdir / "fig_runtime_by_tool_combined.png"
    if plot_grouped_bars(
        runtime_data,
        runtime_tools,
        runtime_modes,
        ylabel="Wall time (hours)",
        title=f"Runtime by tool and mode – {title_suffix}",
        out_path=runtime_path,
    ):
        generated_any = True

    genus_data, genus_tools, genus_modes = collect_genus_f1(summary_rows)
    genus_path = outdir / "fig_genus_f1_combined.png"
    if plot_grouped_bars(
        genus_data,
        genus_tools,
        genus_modes,
        ylabel="Mean genus F1",
        title=f"Genus-level F1 – {title_suffix}",
        out_path=genus_path,
    ):
        generated_any = True

    if not generated_any:
        print("[combined] No figures generated (missing data?)")


if __name__ == "__main__":
    main()
