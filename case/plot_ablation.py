#!/usr/bin/env python3
"""Generate plots for database ablation experiments."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import pandas as pd


def load_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "level_fraction" in df.columns:
        df = df.sort_values("level_fraction")
    # Break out incremental fallback contributions for clearer plots
    df["species_pct"] = df["assigned_species_pct"]
    df["genus_only_pct"] = (df["assigned_genus_pct"] - df["assigned_species_pct"]).clip(lower=0)
    df["family_only_pct"] = (df["assigned_family_pct"] - df["assigned_genus_pct"]).clip(lower=0)
    df["higher_pct"] = df["assigned_higher_pct"].clip(lower=0)
    return df


def _active_categories(df: pd.DataFrame) -> List[tuple[str, pd.Series, str]]:
    categories = [
        ("Species/strain", df["species_pct"], "#1f77b4"),
        ("Genus fallback", df["genus_only_pct"], "#9467bd"),
        ("Family fallback", df["family_only_pct"], "#ff7f0e"),
        ("Higher (root)", df["higher_pct"], "#2ca02c"),
    ]
    if df["genus_only_pct"].max() <= 1e-6:
        categories = [c for c in categories if c[0] != "Genus fallback"]
    if df["family_only_pct"].max() <= 1e-6:
        categories = [c for c in categories if c[0] != "Family fallback"]
    return categories


def plot_rank_fallback(df: pd.DataFrame, out_path: Path) -> None:
    plt.figure(figsize=(8, 4.5))
    x = df["level_fraction"]
    cats = _active_categories(df)
    labels = [label for label, _, _ in cats]
    data = [series for _, series, _ in cats]
    colors = [color for _, _, color in cats]
    plt.stackplot(x, data, labels=labels, colors=colors, alpha=0.85)
    plt.xlabel("Fraction of dominant taxa removed")
    plt.ylabel("Assignments retained (%)")
    plt.title("Rank fallback under database ablation")
    plt.legend(loc="upper right")
    plt.grid(True, linestyle="--", alpha=0.35)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_rank_stack(df: pd.DataFrame, out_path: Path) -> None:
    plt.figure(figsize=(8, 4.5))
    levels = df["level_label"].astype(str)
    species = df["species_pct"]
    bottom = species.copy()
    plt.bar(levels, species, label="Species/strain", color="#1f77b4")
    if df["genus_only_pct"].max() > 1e-6:
        genus = df["genus_only_pct"]
        plt.bar(levels, genus, bottom=bottom, label="Genus fallback", color="#9467bd")
        bottom = bottom + genus
    if df["family_only_pct"].max() > 1e-6:
        family = df["family_only_pct"]
        plt.bar(levels, family, bottom=bottom, label="Family fallback", color="#ff7f0e")
        bottom = bottom + family
    higher = df["higher_pct"]
    plt.bar(levels, higher, bottom=bottom, label="Higher (root)", color="#2ca02c")
    plt.xlabel("Ablation level (%)")
    plt.ylabel("Assignments (%)")
    plt.title("Assignment distribution by rank")
    plt.legend()
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_eval_metrics(df: pd.DataFrame, out_path: Path) -> None:
    if df.empty:
        return
    ranks = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
    plt.figure(figsize=(9, 5))
    for rank in ranks:
        subset = df[df["rank"] == rank]
        if subset.empty:
            continue
        plt.plot(subset["level_fraction"], subset["F1"].astype(float), marker="o", label=rank.title())
    if plt.gca().has_data():
        plt.xlabel("Fraction of dominant taxa removed")
        plt.ylabel("F1 score (%)")
        plt.title("F1 by rank under database ablation")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.4)
        plt.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=200)
        plt.close()
    else:
        plt.close()


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot case-study ablation results.")
    ap.add_argument("--summary", required=True, help="case/ablation_summary.tsv")
    ap.add_argument("--eval", help="case/ablation_eval_summary.tsv (optional)")
    ap.add_argument("--outdir", required=True, help="Directory for output figures")
    args = ap.parse_args()

    summary_path = Path(args.summary)
    eval_path = Path(args.eval) if args.eval else None
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_summary(summary_path)
    plot_rank_fallback(df, outdir / "fig_ablation_rank_fallback.png")
    plot_rank_stack(df, outdir / "fig_ablation_rank_stack.png")

    if eval_path and eval_path.is_file():
        eval_df = pd.read_csv(eval_path, sep="\t")
        plot_eval_metrics(eval_df, outdir / "fig_ablation_f1_by_rank.png")


if __name__ == "__main__":
    main()
