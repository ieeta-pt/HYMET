#!/usr/bin/env python3
"""Convert SnakeMAGs GTDB-Tk summaries into CAMI-compatible profiles."""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Tuple

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
    import sys

    sys.path.append(os.path.dirname(__file__))
    from common import (  # type: ignore
        RANKS,
        rollup_to_ancestors,
        taxonkit_name2taxid,
        taxonkit_taxpath,
        write_cami_profile,
    )
else:  # pragma: no cover
    from .common import RANKS, rollup_to_ancestors, taxonkit_name2taxid, taxonkit_taxpath, write_cami_profile


RANK_ALIAS = {
    "d": "superkingdom",
    "k": "superkingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}


def load_mapping(path: Path) -> Dict[str, Tuple[str, float]]:
    mapping: Dict[str, Tuple[str, float]] = {}
    if not path.is_file():
        return mapping
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for idx, row in enumerate(reader):
            if idx == 0 and row and row[0] == "mag_id":
                continue
            if len(row) < 3:
                continue
            mag_id = row[0].strip()
            contig_id = row[1].strip()
            try:
                length = float(row[2])
            except ValueError:
                length = 0.0
            if mag_id and contig_id:
                mapping[mag_id] = (contig_id, max(length, 0.0))
    return mapping


def parse_lineage(classification: str) -> Dict[str, str]:
    ranks: Dict[str, str] = {}
    if not classification:
        return ranks
    parts = classification.split(";")
    for part in parts:
        if "__" not in part:
            continue
        prefix, name = part.split("__", 1)
        prefix = prefix.strip().lower()
        name = name.strip()
        if not name or name in {"unclassified", "uncultured"}:
            continue
        clean_name = name.replace("_", " ")
        rank = RANK_ALIAS.get(prefix)
        if rank:
            ranks[rank] = clean_name
    return ranks


def load_summary(path: Path) -> List[Tuple[str, Dict[str, str]]]:
    if not path.is_file():
        return []
    rows: List[Tuple[str, Dict[str, str]]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            genome = (row.get("user_genome") or row.get("user_genome_name") or row.get("genome") or "").strip()
            if not genome:
                continue
            classification = (row.get("classification") or row.get("classification_taxonomy") or "").strip()
            ranks = parse_lineage(classification)
            rows.append((genome, ranks))
    return rows


def select_taxon(
    ranks: Dict[str, str],
    name_to_taxid: Dict[str, Tuple[str, str]],
) -> Tuple[str, str]:
    for rank in ("species", "genus", "family", "order", "class", "phylum", "superkingdom"):
        name = ranks.get(rank)
        if not name:
            continue
        hit = name_to_taxid.get(name)
        if hit:
            taxid, actual_rank = hit
            mapped_rank = actual_rank if actual_rank in RANKS else rank
            if taxid and taxid != "0":
                return taxid, mapped_rank
    return "", ""


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert SnakeMAGs GTDB summaries to CAMI format.")
    ap.add_argument("--summary", required=True, help="Path to gtdbtk.summary.tsv.")
    ap.add_argument("--mapping", required=True, help="MAG-to-contig mapping with lengths.")
    ap.add_argument("--out", required=True, help="Output CAMI profile TSV.")
    ap.add_argument("--sample-id", required=True, help="Sample identifier.")
    ap.add_argument("--tool", default="snakemags", help="Tool identifier for CAMI metadata.")
    ap.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="TaxonKit database directory.")
    ap.add_argument("--classified-out", help="Optional path for classified_sequences.tsv.")
    args = ap.parse_args()

    mapping = load_mapping(Path(args.mapping))
    summary_rows = load_summary(Path(args.summary))

    if not mapping or not summary_rows:
        write_cami_profile([], args.out, args.sample_id, args.tool, normalise=False)
        classified_path = Path(args.classified_out) if args.classified_out else None
        if classified_path and classified_path.exists():
            classified_path.unlink()
        return

    names: List[str] = []
    for _, rank_map in summary_rows:
        for name in rank_map.values():
            if name and name not in names:
                names.append(name)

    name_to_taxid = taxonkit_name2taxid(names, args.taxdb)

    weight_by_taxon: DefaultDict[Tuple[str, str], float] = defaultdict(float)
    assignments: Dict[str, str] = {}
    selected_taxids: List[str] = []

    for genome, rank_map in summary_rows:
        mag_id = genome
        if mag_id.endswith(".fa") or mag_id.endswith(".fna"):
            mag_id = Path(mag_id).stem
        contig_info = mapping.get(mag_id)
        if not contig_info:
            continue
        contig_id, length = contig_info
        if length <= 0:
            length = 1.0
        taxid, rank = select_taxon(rank_map, name_to_taxid)
        if not taxid:
            continue
        weight_by_taxon[(taxid, rank or "species")] += length
        assignments[contig_id] = taxid
        if taxid not in selected_taxids:
            selected_taxids.append(taxid)

    if not weight_by_taxon:
        write_cami_profile([], args.out, args.sample_id, args.tool, normalise=False)
        classified_path = Path(args.classified_out) if args.classified_out else None
        if classified_path and classified_path.exists():
            classified_path.unlink()
        return

    tax_paths = taxonkit_taxpath(selected_taxids, args.taxdb)
    default_path = "|".join(["NA"] * len(RANKS))
    cami_rows: List[Dict[str, object]] = []

    for (taxid, rank), weight in weight_by_taxon.items():
        ids_str, names_str = tax_paths.get(taxid, (default_path, default_path))
        taxpath = ids_str.split("|")
        names = names_str.split("|")
        if len(taxpath) < len(RANKS):
            taxpath += ["NA"] * (len(RANKS) - len(taxpath))
        if len(names) < len(RANKS):
            names += ["NA"] * (len(RANKS) - len(names))
        cami_rows.append(
            {
                "taxid": taxid,
                "rank": rank,
                "taxpath": taxpath[: len(RANKS)],
                "taxpathsn": names[: len(RANKS)],
                "percentage": weight,
            }
        )

    cami_rows = rollup_to_ancestors(cami_rows)
    write_cami_profile(cami_rows, args.out, args.sample_id, args.tool, normalise=True)

    if args.classified_out:
        classified_path = Path(args.classified_out)
        if assignments:
            classified_path.parent.mkdir(parents=True, exist_ok=True)
            with classified_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["Query", "TaxID"])
                for contig, taxid in sorted(assignments.items()):
                    writer.writerow([contig, taxid])
        elif classified_path.exists():
            classified_path.unlink()


if __name__ == "__main__":
    main()
