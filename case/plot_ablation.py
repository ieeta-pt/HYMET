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
    return df


def plot_rank_fallback(df: pd.DataFrame, out_path: Path) -> None:
    plt.figure(figsize=(8, 4.5))
    x = df["level_fraction"]
    plt.plot(x, df["assigned_species_pct"], marker="o", label="Species")
    plt.plot(x, df["assigned_genus_pct"], marker="o", label="≤ Genus")
    plt.plot(x, df["assigned_family_pct"], marker="o", label="≤ Family")
    plt.plot(x, df["assigned_higher_pct"], marker="o", label="Higher ranks")
    plt.xlabel("Fraction of dominant taxa removed")
    plt.ylabel("Assignments retained (%)")
    plt.title("Rank fallback under database ablation")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.4)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_rank_stack(df: pd.DataFrame, out_path: Path) -> None:
    plt.figure(figsize=(8, 4.5))
    levels = df["level_label"].astype(str)
    species = df["assigned_species_pct"]
    genus = df["assigned_genus_pct"] - df["assigned_species_pct"]
    family = df["assigned_family_pct"] - df["assigned_genus_pct"]
    higher = df["assigned_higher_pct"]
    bottom = species
    plt.bar(levels, species, label="Species")
    plt.bar(levels, genus, bottom=species, label="Genus (fallback)")
    plt.bar(levels, family, bottom=species + genus, label="Family (fallback)")
    plt.bar(levels, higher, bottom=species + genus + family, label="Higher ranks")
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
