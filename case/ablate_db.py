#!/usr/bin/env python3
"""
Create ablated FASTA references by removing fractions of sequences that belong
to specified TaxIDs. Fractions are applied per-taxid to allow progressive
database incompleteness experiments.
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple


def parse_levels(text: str) -> List[float]:
    parts = [p.strip() for p in text.split(",") if p.strip()]
    levels: List[float] = []
    for p in parts:
        try:
            val = float(p)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"Invalid level '{p}'") from exc
        if not (0.0 <= val <= 1.0):
            raise argparse.ArgumentTypeError(f"Ablation levels must be between 0 and 1 inclusive (got {val}).")
        levels.append(val)
    if not levels:
        raise argparse.ArgumentTypeError("At least one ablation level is required.")
    return sorted(set(levels))


def load_seqmap(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            seq, tax = row[0].strip(), row[1].strip()
            if seq and tax:
                mapping[seq] = tax
    return mapping


def group_sequences_by_taxa(mapping: Dict[str, str], targets: Set[str]) -> Dict[str, List[str]]:
    grouped: Dict[str, List[str]] = {tax: [] for tax in targets}
    for seq, tax in mapping.items():
        if tax in grouped:
            grouped[tax].append(seq)
    return grouped


def determine_removals(grouped: Dict[str, List[str]], level: float, rng: random.Random) -> Set[str]:
    to_remove: Set[str] = set()
    for taxid, seqs in grouped.items():
        if not seqs:
            continue
        count = int(round(level * len(seqs), 0))
        if count <= 0:
            continue
        choices = rng.sample(seqs, min(count, len(seqs)))
        to_remove.update(choices)
    return to_remove


def write_ablated_fasta(
    fasta_path: Path,
    out_path: Path,
    removal_set: Set[str],
) -> Tuple[int, int]:
    total = 0
    removed = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with fasta_path.open() as fin, out_path.open("w") as fout:
        keep = True
        current_id = None
        for line in fin:
            if line.startswith(">"):
                total += 1
                current_id = line[1:].strip().split()[0]
                keep = current_id not in removal_set
                if not keep:
                    removed += 1
                else:
                    fout.write(line)
            else:
                if keep:
                    fout.write(line)
    return total, removed


def main() -> None:
    ap = argparse.ArgumentParser(description="Create ablated FASTA files removing fractions of selected taxa.")
    ap.add_argument("--fasta", required=True, help="Input FASTA used for baseline databases.")
    ap.add_argument("--seqmap", required=True, help="TSV mapping sequence IDs to TaxIDs (seqid<TAB>taxid).")
    ap.add_argument("--taxa", required=True, help="Comma-separated list of TaxIDs to ablate.")
    ap.add_argument("--levels", default="0,0.25,0.5,0.75,1.0", help="Comma-separated fractions (0-1). Default: 0,0.25,0.5,0.75,1.0")
    ap.add_argument("--out-dir", required=True, help="Directory where ablated FASTA files will be written.")
    ap.add_argument("--prefix", default="combined_subset", help="Prefix for generated FASTA files.")
    ap.add_argument("--seed", type=int, default=1337, help="Random seed for reproducible sampling.")
    args = ap.parse_args()

    fasta_path = Path(args.fasta).resolve()
    seqmap_path = Path(args.seqmap).resolve()
    out_dir = Path(args.out_dir).resolve()

    if not fasta_path.is_file():
        raise SystemExit(f"Input FASTA not found: {fasta_path}")
    if not seqmap_path.is_file():
        raise SystemExit(f"Sequence map missing: {seqmap_path}")

    levels = parse_levels(args.levels)
    targets = {tax.strip() for tax in args.taxa.split(",") if tax.strip()}
    if not targets:
        raise SystemExit("No target TaxIDs provided.")

    mapping = load_seqmap(seqmap_path)
    grouped = group_sequences_by_taxa(mapping, targets)
    rng = random.Random(args.seed)

    summary_path = out_dir / "ablation_summary.tsv"
    out_dir.mkdir(parents=True, exist_ok=True)
    if not summary_path.exists():
        with summary_path.open("w") as summary:
            summary.write("level_fraction\tlevel_label\ttarget_taxid\ttotal_sequences\tdropped_sequences\n")

    for level in levels:
        label = f"{int(level*100):03d}"
        removal_set = determine_removals(grouped, level, rng)
        out_path = out_dir / f"{args.prefix}.ablate{label}.fasta"
        total, removed = write_ablated_fasta(fasta_path, out_path, removal_set)

        with summary_path.open("a") as summary:
            for taxid, seqs in grouped.items():
                count = int(round(level * len(seqs), 0))
                summary.write(f"{level}\t{label}\t{taxid}\t{len(seqs)}\t{min(count, len(seqs))}\n")

        print(f"[ablate] level={level:.2f} ({label}) → wrote {out_path.name} (removed {removed}/{total} sequences)")


if __name__ == "__main__":
    main()
