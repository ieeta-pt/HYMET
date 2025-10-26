#!/usr/bin/env python3
"""Aggregate evaluation metrics across CAMI benchmark runs."""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

PROFILE_SUMMARY = "profile_summary.tsv"
CONTIG_SUMMARY = "contigs_per_rank.tsv"


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def collect_eval(sample: str, tool: str, eval_dir: Path) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    prof_path = eval_dir / PROFILE_SUMMARY
    contig_path = eval_dir / CONTIG_SUMMARY
    prof_rows: List[Dict[str, str]] = []
    contig_rows: List[Dict[str, str]] = []
    if prof_path.is_file() and prof_path.stat().st_size > 0:
        prof_rows = read_tsv(prof_path)
        for row in prof_rows:
            row.update({"sample": sample, "tool": tool})
    if contig_path.is_file() and contig_path.stat().st_size > 0:
        raw_rows = read_tsv(contig_path)
        filtered_rows: List[Dict[str, str]] = []
        for row in raw_rows:
            n_val = (row.get("n") or "").strip()
            try:
                n_float = float(n_val)
            except ValueError:
                continue
            if n_float <= 0:
                continue
            row.update({"sample": sample, "tool": tool})
            filtered_rows.append(row)
        contig_rows = filtered_rows
    return prof_rows, contig_rows


def average_metrics(rows: List[Dict[str, str]], keys: List[str]) -> Dict[str, float]:
    agg = {k: 0.0 for k in keys}
    count = {k: 0 for k in keys}
    for row in rows:
        for key in keys:
            val = row.get(key)
            if val is None or val == "":
                continue
            try:
                agg[key] += float(val)
                count[key] += 1
            except ValueError:
                continue
    return {k: (agg[k] / count[k] if count[k] else 0.0) for k in keys}


def write_summary(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def main() -> None:
    ap = argparse.ArgumentParser(description="Aggregate CAMI evaluation metrics across samples and tools.")
    ap.add_argument("--bench-root", default=str(Path(__file__).resolve().parent), help="Bench directory (default: script parent).")
    ap.add_argument("--samples-root", help="Directory that contains per-sample outputs (default: <bench-root>/out).")
    ap.add_argument("--outdir", default="out", help="Relative or absolute output directory for aggregated tables.")
    args = ap.parse_args()

    bench_root = Path(args.bench_root)
    samples_root = Path(args.samples_root) if args.samples_root else bench_root / "out"
    out_path = Path(args.outdir)
    out_root = out_path if out_path.is_absolute() else bench_root / out_path
    out_root.mkdir(parents=True, exist_ok=True)

    per_sample_rows: List[Dict[str, str]] = []
    contig_rows: List[Dict[str, str]] = []

    if not samples_root.is_dir():
        print(f"[aggregate] No benchmark outputs under {samples_root}; skipping aggregation.")
        return

    sample_dirs = [p for p in samples_root.iterdir() if p.is_dir()]
    for sdir in sorted(sample_dirs, key=lambda p: p.name):
        sample = sdir.name
        for tool_dir in sorted([p for p in sdir.iterdir() if p.is_dir()]):
            tool = tool_dir.name
            eval_dir = tool_dir / "eval"
            if not eval_dir.is_dir():
                continue
            prof_rows, cont_rows = collect_eval(sample, tool, eval_dir)
            per_sample_rows.extend(prof_rows)
            contig_rows.extend(cont_rows)

    if per_sample_rows:
        summary_fields = [
            "sample",
            "tool",
            "rank",
            "L1_total_variation_pctpts",
            "BrayCurtis_pct",
            "Precision_%",
            "Recall_%",
            "F1_%",
            "TP",
            "FP",
            "FN",
        ]
        write_summary(out_root / "summary_per_tool_per_sample.tsv", per_sample_rows, summary_fields)

        metrics_by_tool_rank: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
        for row in per_sample_rows:
            metrics_by_tool_rank[(row["tool"], row["rank"])].append(row)

        leaderboard_rows: List[Dict[str, str]] = []
        metric_keys = [
            "L1_total_variation_pctpts",
            "BrayCurtis_pct",
            "Precision_%",
            "Recall_%",
            "F1_%",
        ]
        for (tool, rank), rows in sorted(metrics_by_tool_rank.items()):
            avg = average_metrics(rows, metric_keys)
            leaderboard_rows.append(
                {
                    "tool": tool,
                    "rank": rank,
                    "samples": str(len(rows)),
                    "mean_L1_total_variation_pctpts": f"{avg['L1_total_variation_pctpts']:.4f}",
                    "mean_BrayCurtis_pct": f"{avg['BrayCurtis_pct']:.4f}",
                    "mean_Precision_%": f"{avg['Precision_%']:.2f}",
                    "mean_Recall_%": f"{avg['Recall_%']:.2f}",
                    "mean_F1_%": f"{avg['F1_%']:.2f}",
                }
            )

        leaderboard_fields = [
            "tool",
            "rank",
            "samples",
            "mean_L1_total_variation_pctpts",
            "mean_BrayCurtis_pct",
            "mean_Precision_%",
            "mean_Recall_%",
            "mean_F1_%",
        ]
        write_summary(out_root / "leaderboard_by_rank.tsv", leaderboard_rows, leaderboard_fields)

    if contig_rows:
        contig_fields = ["sample", "tool", "rank", "n", "correct", "accuracy_percent"]
        write_summary(out_root / "contig_accuracy_per_tool.tsv", contig_rows, contig_fields)


if __name__ == "__main__":
    main()
