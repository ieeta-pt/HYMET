#!/usr/bin/env python3
"""
Rewrite superkingdom rows in CAMI profile TSVs so they match the identifiers used
in a reference (truth) profile. Converters that rely on GTDB taxonomy often emit
`Bacillati`/`Pseudomonadati` (taxids 1783272/3379134) while the CAMI ground truth
still references the canonical NCBI `Bacteria` (taxid 2). By rebuilding the
superkingdom rows from lower ranks, this helper aligns predicted profiles with
the truth (or falls back to canonical NCBI IDs) so evaluation metrics at the
highest rank are meaningful again.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple

RANK_PRIORITY = [
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

CANONICAL_SUPERKINGDOMS = {
    "2",      # Bacteria
    "2157",   # Archaea
    "2759",   # Eukaryota
    "10239",  # Viruses
    "12884",  # Viroids
}


def load_taxonomy(path: Path) -> Dict[str, Tuple[str, str, str]]:
    """Return mapping TaxID -> (parent, rank, name)."""
    taxonomy: Dict[str, Tuple[str, str, str]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            taxonomy[row["TaxID"]] = (row["ParentTaxID"], row["Rank"].lower(), row["Name"])
    return taxonomy


def canonical_superkingdom(taxid: str, taxonomy: Dict[str, Tuple[str, str, str]]) -> str:
    """Ascend the taxonomy tree to the canonical superkingdom/domain identifier."""
    current = taxid
    visited = set()
    while current and current not in visited:
        visited.add(current)
        if current in CANONICAL_SUPERKINGDOMS:
            return current
        parent, rank, _ = taxonomy.get(current, ("", "", ""))
        if not parent or parent == current:
            break
        current = parent
    return taxid


def align_to_targets(taxid: str, targets: set[str], taxonomy: Dict[str, Tuple[str, str, str]]) -> str | None:
    """Ascend until a taxid present in the target set is found."""
    current = taxid
    visited = set()
    while current and current not in visited:
        visited.add(current)
        if current in targets:
            return current
        parent, _, _ = taxonomy.get(current, ("", "", ""))
        if not parent or parent == current:
            break
        current = parent
    return None


def load_truth_superkingdoms(path: Path) -> set[str]:
    targets: set[str] = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith(("@", "#")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1].lower() == "superkingdom":
                targets.add(parts[0])
    return targets


def rewrite_profile(profile: Path, taxonomy: Dict[str, Tuple[str, str, str]], targets: set[str]) -> None:
    lines: List[str] = []
    with profile.open() as handle:
        lines = handle.readlines()

    if not lines:
        return

    header_prefixes = ("@", "#")
    body: List[List[str]] = []
    header_lines: List[str] = []
    for line in lines:
        if line.startswith(header_prefixes):
            header_lines.append(line)
        else:
            body.append(line.rstrip("\n").split("\t"))

    if not body:
        return

    remainder: List[List[str]] = []
    aggregates: Dict[str, float] = {}

    original_super: List[List[str]] = [row for row in body if len(row) >= 2 and row[1].lower() == "superkingdom"]

    available_ranks = {row[1].lower() for row in body if len(row) >= 2 and row[1].lower() != "superkingdom"}
    aggregate_rank = next((rank for rank in RANK_PRIORITY if rank in available_ranks), None)

    for row in body:
        if len(row) < 5:
            continue
        taxid, rank = row[0], row[1].lower()
        if rank == "superkingdom":
            continue  # rebuild later
        remainder.append(row)
        if aggregate_rank and rank != aggregate_rank:
            continue
        try:
            perc = float(row[4])
        except ValueError:
            perc = 0.0
        target = align_to_targets(taxid, targets, taxonomy)
        if not target:
            target = canonical_superkingdom(taxid, taxonomy)
        aggregates[target] = aggregates.get(target, 0.0) + perc

    # Ensure targets exist in aggregates even if zero to avoid dropping strata
    for target in targets:
        aggregates.setdefault(target, 0.0)

    super_rows: List[List[str]] = []
    if aggregates:
        for taxid, perc in sorted(aggregates.items()):
            if perc <= 0:
                continue
            parent, _, name = taxonomy.get(taxid, ("", "", ""))
            width = 7
            path = [taxid] + ["NA"] * (width - 1)
            names = [name or "NA"] + ["NA"] * (width - 1)
            super_rows.append([
                taxid,
                "superkingdom",
                "|".join(path),
                "|".join(names),
                f"{perc:.6f}",
            ])
        existing_keys = {row[0] for row in super_rows}
        for row in original_super:
            if row[0] not in existing_keys:
                super_rows.append(row)
    else:
        super_rows = original_super

    with profile.open("w") as handle:
        handle.writelines(header_lines)
        for row in super_rows + remainder:
            handle.write("\t".join(row) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Canonicalize superkingdom taxids in CAMI profile TSV files.")
    parser.add_argument("--profile", required=True, help="Path to the CAMI profile TSV to rewrite in-place.")
    parser.add_argument(
        "--truth-profile",
        required=True,
        help="Truth profile used for evaluation; determines which superkingdom IDs to align with.",
    )
    parser.add_argument(
        "--taxonomy",
        default=None,
        help="Path to taxonomy_hierarchy.tsv (defaults to HYMET/data/taxonomy_hierarchy.tsv relative to this script).",
    )
    args = parser.parse_args()

    profile_path = Path(args.profile)
    if not profile_path.exists():
        raise SystemExit(f"profile not found: {profile_path}")

    truth_path = Path(args.truth_profile)
    if not truth_path.exists():
        raise SystemExit(f"truth profile not found: {truth_path}")

    taxonomy_path = Path(args.taxonomy) if args.taxonomy else (Path(__file__).resolve().parents[2] / "data" / "taxonomy_hierarchy.tsv")
    if not taxonomy_path.exists():
        raise SystemExit(f"taxonomy hierarchy not found: {taxonomy_path}")

    taxonomy = load_taxonomy(taxonomy_path)
    targets = load_truth_superkingdoms(truth_path)
    if not targets and not CANONICAL_SUPERKINGDOMS:
        raise SystemExit("no superkingdoms found in truth profile")
    rewrite_profile(profile_path, taxonomy, targets or CANONICAL_SUPERKINGDOMS)


if __name__ == "__main__":
    main()
