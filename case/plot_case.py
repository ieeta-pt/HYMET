#!/usr/bin/env python3
"""Generate publication-quality summary figures for the HYMET case study outputs."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


CASES_ROOT = Path(__file__).resolve().parent
DEFAULT_OUT_ROOT = CASES_ROOT / "out"
DEFAULT_FIG_ROOT = DEFAULT_OUT_ROOT / "figures"
DEFAULT_MAX_TAXA = 12  # show top-N taxa per sample


@dataclass
class SampleResult:
    name: str
    runtime_wall: float = 0.0
    runtime_cpu: float = 0.0
    runtime_rss: float = 0.0
    top_taxa: List[Dict[str, str]] = field(default_factory=list)

    @property
    def display_name(self) -> str:
        return self.name.replace("_", " ")


def load_runtime(sample: str, table: Path) -> Dict[str, float]:
    if not table.exists():
        return {}
    with table.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["sample"] == sample and row["tool"] == "hymet" and row["stage"] == "run":
                return {
                    "wall": float(row["wall_seconds"]),
                    "cpu": float(row["user_seconds"]) + float(row["sys_seconds"]),
                    "rss": float(row["max_rss_gb"]),
                }
    return {}


def load_top_taxa(sample_dir: Path) -> List[Dict[str, str]]:
    tsv = sample_dir / "top_taxa.tsv"
    if not tsv.exists():
        return []
    taxa: List[Dict[str, str]] = []
    with tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                pct = float(row.get("Percentage", "0"))
            except ValueError:
                continue
            label = row.get("TaxPathSN") or row.get("TaxPath") or row.get("TaxID") or "Unknown"
            taxa.append({
                "rank": row.get("Rank", ""),
                "label": label,
                "pct": pct,
            })
    taxa.sort(key=lambda item: item["pct"], reverse=True)
    return taxa


def collect_case_results(out_root: Path, max_taxa: int) -> List[SampleResult]:
    results: List[SampleResult] = []
    runtime_table = out_root / "runtime_memory.tsv"
    if not runtime_table.exists():
        return results
    for sample_dir in sorted(out_root.iterdir()):
        if not sample_dir.is_dir():
            continue
        sample = SampleResult(name=sample_dir.name)
        runtime = load_runtime(sample_dir.name, runtime_table)
        if runtime:
            sample.runtime_wall = runtime["wall"]
            sample.runtime_cpu = runtime["cpu"]
            sample.runtime_rss = runtime["rss"]
        taxa = load_top_taxa(sample_dir)
        sample.top_taxa = taxa[:max_taxa]
        results.append(sample)
    return results


def format_duration(seconds: float) -> str:
    if seconds >= 3600:
        hours = seconds / 3600.0
        return f"{hours:.1f} h"
    if seconds >= 60:
        mins = seconds / 60.0
        return f"{mins:.1f} min"
    return f"{seconds:.0f} s"


def plot_runtime(results: List[SampleResult], fig_root: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    palette = sns.color_palette("deep", len(results))
    names = [r.display_name for r in results]

    wall_minutes = np.array([r.runtime_wall / 60.0 for r in results])
    mem_gb = np.array([r.runtime_rss for r in results])

    def annotate(ax, bars, values, formatter):
        threshold = values.max() * 0.25 if values.size else 0
        for bar, val in zip(bars, values):
            y = bar.get_y() + bar.get_height() / 2
            text = formatter(val)
            if val > threshold:
                ax.text(val - values.max() * 0.02, y, text, va="center", ha="right",
                        fontsize=10, weight="semibold", color="white")
            else:
                ax.text(val + values.max() * 0.05, y, text, va="center", ha="left",
                        fontsize=10, color="#1c2333")

    bars = axes[0].barh(names, wall_minutes, color=palette, edgecolor="white")
    axes[0].set_xlabel("Wall-clock Time (minutes)")
    axes[0].set_title("HYMET runtime per sample")
    annotate(axes[0], bars, wall_minutes, lambda v: format_duration(v * 60.0))

    bars_mem = axes[1].barh(names, mem_gb, color=palette, edgecolor="white")
    axes[1].set_xlabel("Peak Resident Set Size (GB)")
    axes[1].set_title("Memory footprint")
    annotate(axes[1], bars_mem, mem_gb, lambda v: f"{v:.2f} GB")

    sns.despine(fig)
    fig.suptitle("HYMET Case Study – Performance Summary", fontsize=15, weight="semibold")
    fig_root.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_root / "fig_case_runtime.png", dpi=300)
    plt.close(fig)


def taxon_label(entry: Dict[str, str]) -> str:
    label = entry["label"]
    if "|" in label:
        label = label.split("|")[-1]
    label = label.replace("_", " ")
    return textwrap.shorten(label, width=32, placeholder="…")


def plot_top_taxa_panels(results: List[SampleResult], fig_root: Path) -> None:
    if not results:
        return
    n = len(results)
    fig, axes = plt.subplots(n, 1, figsize=(11, 3.2 * n), constrained_layout=True)
    if n == 1:
        axes = [axes]
    palette = sns.color_palette("Spectral", MAX_TAXA)

    for ax, sample in zip(axes, results):
        taxa = sample.top_taxa
        if not taxa:
            ax.axis("off")
            continue
        labels = [taxon_label(t) for t in taxa][::-1]
        values = [t["pct"] for t in taxa][::-1]
        bars = ax.barh(labels, values, color=palette[:len(values)], edgecolor="white")
        ax.set_xlabel("Relative abundance (%)")
        ax.set_xlim(0, max(values) * 1.05)
        ax.set_title(f"{sample.display_name} – Top {len(values)} taxa", loc="left", weight="semibold")
        for bar, val in zip(bars, values):
            xpos = bar.get_width()
            y = bar.get_y() + bar.get_height() / 2
            if xpos > max(values) * 0.6:
                ax.text(xpos - max(values) * 0.02, y, f"{val:.2f}", va="center", ha="right",
                        fontsize=9, color="white", weight="semibold")
            else:
                ax.text(xpos + max(values) * 0.02, y, f"{val:.2f}", va="center", ha="left",
                        fontsize=9, color="#1c2333")
        sns.despine(ax=ax, left=True)

    fig_root.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_root / "fig_case_top_taxa_panels.png", dpi=300)
    plt.close(fig)


def plot_taxa_heatmap(results: List[SampleResult], fig_root: Path) -> None:
    taxa_names = sorted({taxon_label(t) for r in results for t in r.top_taxa})
    if not taxa_names:
        return
    data = np.zeros((len(taxa_names), len(results)))
    for j, sample in enumerate(results):
        for t in sample.top_taxa:
            i = taxa_names.index(taxon_label(t))
            data[i, j] = t["pct"]

    fig, ax = plt.subplots(figsize=(10, 0.45 * len(taxa_names) + 1.5), constrained_layout=True)
    sns.heatmap(data,
                annot=True,
                fmt=".1f",
                cmap="rocket_r",
                annot_kws={"fontsize": 9, "color": "black"},
                linewidths=0.4,
                linecolor="white",
                cbar_kws={"label": "Relative abundance (%)"},
                ax=ax)
    ax.set_xticklabels([r.display_name for r in results], rotation=30, ha="right")
    ax.set_yticklabels(taxa_names, rotation=0)
    ax.set_title("Top taxa overlap across case-study samples", loc="left", weight="semibold")
    fig_root.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_root / "fig_case_top_taxa_heatmap.png", dpi=300)
    plt.close(fig)


def save_metadata(results: List[SampleResult], fig_root: Path) -> None:
    payload = {
        "samples": [
            {
                "name": r.name,
                "runtime": {
                    "wall_seconds": r.runtime_wall,
                    "cpu_seconds": r.runtime_cpu,
                    "rss_gb": r.runtime_rss,
                },
                "top_taxa": r.top_taxa,
            }
            for r in results
        ]
    }
    fig_root.mkdir(parents=True, exist_ok=True)
    (fig_root / "case_figures_metadata.json").write_text(json.dumps(payload, indent=2))


def generate_figures(args: argparse.Namespace) -> None:
    sns.set_theme(context="talk", style="whitegrid", font="DejaVu Sans")
    out_root = Path(args.case_root).resolve()
    fig_root = Path(args.figures_dir).resolve()
    results = collect_case_results(out_root, args.max_taxa)
    if not results:
        raise SystemExit(f"No case-study outputs found under {out_root}.")
    plot_runtime(results, fig_root)
    plot_top_taxa_panels(results, fig_root)
    plot_taxa_heatmap(results, fig_root)
    save_metadata(results, fig_root)
    print(f"Generated figures in {fig_root}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", default=str(DEFAULT_OUT_ROOT), help="Directory containing case outputs (default: case/out)")
    parser.add_argument("--figures-dir", help="Directory to store generated figures (default: <case-root>/figures)")
    parser.add_argument("--max-taxa", type=int, default=DEFAULT_MAX_TAXA, help="Top-N taxa to display per sample (default: 12)")
    args = parser.parse_args()
    if args.figures_dir is None:
        args.figures_dir = str(Path(args.case_root).resolve() / "figures")
    return args


def main() -> None:
    args = parse_args()
    generate_figures(args)


if __name__ == "__main__":
    main()
