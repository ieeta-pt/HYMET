#!/usr/bin/env python3
"""Convert phyloFlash NTU abundance outputs into CAMI profile format."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

if __package__ is None or __package__ == "":  # pragma: no cover
    import sys

    sys.path.append(os.path.dirname(__file__))
    from common import (  # type: ignore
        RANKS,
        taxonkit_name2taxid,
        taxonkit_taxpath,
        rollup_to_ancestors,
        write_cami_profile,
    )
else:  # pragma: no cover
    from .common import RANKS, taxonkit_name2taxid, taxonkit_taxpath, rollup_to_ancestors, write_cami_profile

RankPath = Dict[str, str]
RANK_ORDER = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


def _clean_name(raw: str) -> str:
    name = raw.strip()
    if not name:
        return ""
    if name.startswith("(") and name.endswith(")"):
        name = name[1:-1].strip()
    if not name or name.lower() in {"unclassified", "unassigned", "unknown", "na"}:
        return ""
    return name


def _read_phyloflash(path: Path) -> List[Tuple[RankPath, float]]:
    entries: List[Tuple[RankPath, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row:
                continue
            lineage_raw = row[0].strip()
            if not lineage_raw:
                continue
            try:
                count = float(row[1]) if len(row) > 1 else 0.0
            except ValueError:
                continue
            parts = [p.strip() for p in lineage_raw.split(";")]
            while len(parts) < len(RANK_ORDER):
                parts.append("")
            lineage: RankPath = {}
            for rank, value in zip(RANK_ORDER, parts[: len(RANK_ORDER)]):
                cleaned = _clean_name(value)
                if cleaned:
                    lineage[rank] = cleaned
            if lineage and count > 0:
                entries.append((lineage, count))
    return entries


def _choose_rank(lineage: RankPath) -> Optional[Tuple[str, str]]:
    for rank in reversed(RANK_ORDER):
        if rank in lineage:
            return lineage[rank], rank
    return None


def build_profile(
    entries: List[Tuple[RankPath, float]],
    sample_id: str,
    tool: str,
    taxdb: str,
    profile_path: Path,
    classified_path: Path,
) -> None:
    name_rank_counts: Dict[Tuple[str, str], float] = {}
    names: set[str] = set()
    for lineage, count in entries:
        chosen = _choose_rank(lineage)
        if not chosen:
            continue
        name, rank = chosen
        names.add(name)
        key = (name, rank)
        name_rank_counts[key] = name_rank_counts.get(key, 0.0) + count

    if not name_rank_counts:
        write_cami_profile([], str(profile_path), sample_id, tool)
        classified_path.unlink(missing_ok=True)
        return

    name_to_taxid = taxonkit_name2taxid(names, taxdb)
    taxid_counts: Dict[Tuple[str, str], float] = {}
    for (name, rank), count in name_rank_counts.items():
        hit = name_to_taxid.get(name)
        if not hit:
            continue
        taxid, real_rank = hit
        taxid_counts[(taxid, real_rank or rank)] = taxid_counts.get((taxid, real_rank or rank), 0.0) + count

    if not taxid_counts:
        write_cami_profile([], str(profile_path), sample_id, tool)
        classified_path.unlink(missing_ok=True)
        return

    taxids = [taxid for taxid, _ in taxid_counts.keys()]
    taxpaths = taxonkit_taxpath(taxids, taxdb)
    total = sum(taxid_counts.values())
    cami_rows = []
    for (taxid, rank), count in taxid_counts.items():
        ids_raw, names_raw = taxpaths.get(taxid, ("|".join(["NA"] * len(RANKS)), "|".join(["NA"] * len(RANKS))))
        id_vec = ids_raw.split("|")
        name_vec = names_raw.split("|")
        if len(id_vec) < len(RANKS):
            id_vec += ["NA"] * (len(RANKS) - len(id_vec))
        if len(name_vec) < len(RANKS):
            name_vec += ["NA"] * (len(RANKS) - len(name_vec))
        cami_rows.append(
            {
                "taxid": taxid,
                "rank": rank or "species",
                "taxpath": id_vec[: len(RANKS)],
                "taxpathsn": name_vec[: len(RANKS)],
                "percentage": 100.0 * count / total if total else 0.0,
            }
        )

    cami_rows = rollup_to_ancestors(cami_rows)
    write_cami_profile(cami_rows, str(profile_path), sample_id, tool, normalise=True)
    with classified_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Query", "TaxID"])


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert phyloFlash NTU abundance outputs to CAMI format.")
    parser.add_argument("--input", required=True, help="Path to <lib>.phyloFlash.NTUfull_abundance.csv.")
    parser.add_argument("--out", required=True, help="Output CAMI profile TSV.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--tool", default="phyloflash", help="Tool identifier for metadata.")
    parser.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="TaxonKit database directory.")
    parser.add_argument("--classified-out", default="", help="Optional path for classified_sequences.tsv.")
    args = parser.parse_args()

    entries = _read_phyloflash(Path(args.input))
    classified_path = Path(args.classified_out) if args.classified_out else Path(args.out).with_name("classified_sequences.tsv")
    build_profile(entries, args.sample_id, args.tool, args.taxdb, Path(args.out), classified_path)


if __name__ == "__main__":
    main()
