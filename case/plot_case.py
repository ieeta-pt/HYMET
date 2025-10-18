#!/usr/bin/env python3
"""Generate publication-quality summary figures for the HYMET case study outputs."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


CASES_ROOT = Path(__file__).resolve().parent
OUT_ROOT = CASES_ROOT / "out"
FIG_ROOT = OUT_ROOT / "figures"
MAX_TAXA = 12  # show top-N taxa per sample


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
    return taxa[:MAX_TAXA]


def collect_case_results() -> List[SampleResult]:
    results: List[SampleResult] = []
    runtime_table = OUT_ROOT / "runtime_memory.tsv"
    for sample_dir in sorted(OUT_ROOT.iterdir()):
        if not sample_dir.is_dir() or sample_dir.name == "figures":
            continue
        sample = SampleResult(name=sample_dir.name)
        runtime = load_runtime(sample_dir.name, runtime_table)
        if runtime:
            sample.runtime_wall = runtime["wall"]
            sample.runtime_cpu = runtime["cpu"]
            sample.runtime_rss = runtime["rss"]
        sample.top_taxa = load_top_taxa(sample_dir)
        results.append(sample)
    return results


def ensure_output_dir() -> None:
    FIG_ROOT.mkdir(parents=True, exist_ok=True)


def format_duration(seconds: float) -> str:
    if seconds >= 3600:
        hours = seconds / 3600.0
        return f"{hours:.1f} h"
    if seconds >= 60:
        mins = seconds / 60.0
        return f"{mins:.1f} min"
    return f"{seconds:.0f} s"


def plot_runtime(results: List[SampleResult]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), gridspec_kw={"wspace": 0.25})
    palette = sns.color_palette("deep", len(results))
    names = [r.display_name for r in results]

    wall_minutes = np.array([r.runtime_wall / 60.0 for r in results])
    cpu_minutes = np.array([r.runtime_cpu / 60.0 for r in results])
    mem_gb = np.array([r.runtime_rss for r in results])

    axes[0].barh(names, wall_minutes, color=palette)
    axes[0].set_xlabel("Wall-clock Time (minutes)")
    axes[0].set_title("HYMET runtime per sample")
    for bar, val in zip(axes[0].patches, wall_minutes):
        axes[0].text(val + wall_minutes.max() * 0.02, bar.get_y() + bar.get_height() / 2,
                     format_duration(val * 60.0), va="center", fontsize=10)

    axes[1].barh(names, mem_gb, color=palette)
    axes[1].set_xlabel("Peak Resident Set Size (GB)")
    axes[1].set_title("Memory footprint")
    for bar, val in zip(axes[1].patches, mem_gb):
        axes[1].text(val + mem_gb.max() * 0.03, bar.get_y() + bar.get_height() / 2,
                     f"{val:.2f} GB", va="center", fontsize=10)

    sns.despine(fig)
    fig.suptitle("HYMET Case Study – Performance Summary", fontsize=14, weight="semibold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_ROOT / "fig_case_runtime.png", dpi=300)
    plt.close(fig)


def taxon_label(entry: Dict[str, str]) -> str:
    label = entry["label"]
    if "|" in label:
        label = label.split("|")[-1]
    return label.replace("_", " ")


def plot_top_taxa_panels(results: List[SampleResult]) -> None:
    if not results:
        return
    n = len(results)
    fig, axes = plt.subplots(n, 1, figsize=(11, 3.5 * n))
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
        bars = ax.barh(labels, values, color=palette[:len(values)])
        ax.set_xlabel("Relative abundance (%)")
        ax.set_xlim(0, max(values) * 1.1)
        ax.set_title(f"{sample.display_name}: Top {len(values)} taxa", loc="left", weight="semibold")
        for bar, val in zip(bars, values):
            ax.text(val + max(values) * 0.01,
                    bar.get_y() + bar.get_height() / 2,
                    f"{val:.2f}", va="center", fontsize=9)
        sns.despine(ax=ax, left=True)

    fig.tight_layout()
    fig.savefig(FIG_ROOT / "fig_case_top_taxa_panels.png", dpi=300)
    plt.close(fig)


def plot_taxa_heatmap(results: List[SampleResult]) -> None:
    taxa_names = sorted({taxon_label(t) for r in results for t in r.top_taxa})
    if not taxa_names:
        return
    data = np.zeros((len(taxa_names), len(results)))
    for j, sample in enumerate(results):
        for t in sample.top_taxa:
            i = taxa_names.index(taxon_label(t))
            data[i, j] = t["pct"]

    fig, ax = plt.subplots(figsize=(10, 0.35 * len(taxa_names) + 2))
    sns.heatmap(data,
                annot=True,
                fmt=".1f",
                cmap="mako",
                linewidths=0.5,
                linecolor="white",
                cbar_kws={"label": "Relative abundance (%)"},
                ax=ax)
    ax.set_xticklabels([r.display_name for r in results], rotation=45, ha="right")
    ax.set_yticklabels(taxa_names, rotation=0)
    ax.set_title("Top taxa overlap across case-study samples", loc="left", weight="semibold")
    fig.tight_layout()
    fig.savefig(FIG_ROOT / "fig_case_top_taxa_heatmap.png", dpi=300)
    plt.close(fig)


def save_metadata(results: List[SampleResult]) -> None:
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
    (FIG_ROOT / "case_figures_metadata.json").write_text(json.dumps(payload, indent=2))


def generate_figures(_: argparse.Namespace) -> None:
    sns.set_theme(context="talk", style="whitegrid", font="DejaVu Sans")
    ensure_output_dir()
    results = collect_case_results()
    if not results:
        raise SystemExit("No case-study outputs found under case/out/.")
    plot_runtime(results)
    plot_top_taxa_panels(results)
    plot_taxa_heatmap(results)
    save_metadata(results)
    print(f"Generated figures in {FIG_ROOT}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    generate_figures(args)


if __name__ == "__main__":
    main()
