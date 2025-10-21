#!/usr/bin/env python3
"""Compute per-group precision/recall/F1 for the HYMET mutation experiment."""
from __future__ import annotations

import argparse
import json
import csv
from pathlib import Path
from typing import Dict, List

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_ALIAS = {
    "domain": "superkingdom",
    "kingdom": "superkingdom",
    "sk": "superkingdom",
    "k": "superkingdom",
    "superkingdom": "superkingdom",
    "phylum": "phylum",
    "p": "phylum",
    "class": "class",
    "c": "class",
    "order": "order",
    "o": "order",
    "family": "family",
    "f": "family",
    "genus": "genus",
    "g": "genus",
    "species": "species",
    "s": "species",
}

GROUPS = [
    "Viruses",
    "Archaea",
    "Bacteria",
    "Fungi",
    "Plants",
    "Protozoa",
    "Invertebrates",
    "Vertebrate Mammals",
    "Other Vertebrates",
]


def parse_lineage(lineage: str) -> Dict[str, str]:
    result = {rank: "" for rank in RANKS}
    if not lineage or lineage == "Unknown":
        return result
    for part in lineage.split(";"):
        part = part.strip()
        if not part or ":" not in part:
            continue
        rank_raw, name = part.split(":", 1)
        rank = RANK_ALIAS.get(rank_raw.strip().lower())
        if rank:
            result[rank] = name.strip()
    return result


def load_predictions(path: Path) -> Dict[str, Dict[str, str]]:
    predictions: Dict[str, Dict[str, str]] = {}
    if not path.exists():
        return predictions
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            contig = row["Query"].strip()
            lineage = row["Lineage"].strip()
            rank_map = parse_lineage(lineage)
            predictions[contig] = rank_map
    return predictions


def evaluate(
    metadata: Dict[str, dict],
    predictions: Dict[str, Dict[str, str]],
) -> Dict[str, Dict[str, Dict[str, int]]]:
    stats = {
        group: {rank: {"TP": 0, "FP": 0, "FN": 0} for rank in RANKS}
        for group in GROUPS
    }
    for contig, info in metadata.items():
        group = info["group"]
        truth_names: List[str] = info["names"]
        pred_map = predictions.get(contig, {})
        for idx, rank in enumerate(RANKS):
            truth = truth_names[idx]
            if not truth:
                continue
            pred = pred_map.get(rank, "")
            counters = stats[group][rank]
            if pred == truth:
                counters["TP"] += 1
            elif pred:
                counters["FP"] += 1
                counters["FN"] += 1
            else:
                counters["FN"] += 1
    return stats


def compute_metrics(tp: int, fp: int, fn: int):
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    return precision, recall, f1


def format_float(value: float) -> str:
    return f"{value:.2f}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute HYMET mutation metrics")
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--outputs", required=True, type=Path, help="Directory holding hymet_outputs/rate_*/classified_sequences.tsv")
    parser.add_argument("--rates", nargs="*", type=str, required=True, help="List of rates identifiers (e.g. 0.00 0.02 0.05 ...)")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    metadata = json.loads(args.metadata.read_text())
    args.outdir.mkdir(parents=True, exist_ok=True)

    summary_rows = []

    for rate in args.rates:
        rate_dir = args.outputs / f"rate_{rate}"
        classified_path = rate_dir / "classified_sequences.tsv"
        predictions = load_predictions(classified_path)
        stats = evaluate(metadata, predictions)
        for group in GROUPS:
            for rank in RANKS:
                tp = stats[group][rank]["TP"]
                fp = stats[group][rank]["FP"]
                fn = stats[group][rank]["FN"]
                precision, recall, f1 = compute_metrics(tp, fp, fn)
                summary_rows.append(
                    {
                        "rate": rate,
                        "group": group,
                        "rank": rank,
                        "precision": precision,
                        "recall": recall,
                        "f1": f1,
                        "tp": tp,
                        "fp": fp,
                        "fn": fn,
                    }
                )

    output_csv = args.outdir / "mutation_metrics.tsv"
    with output_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["rate", "group", "rank", "precision", "recall", "f1", "tp", "fp", "fn"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    print(f"[metrics] wrote {output_csv}")


if __name__ == "__main__":
    main()
