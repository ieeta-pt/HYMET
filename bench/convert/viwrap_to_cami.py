#!/usr/bin/env python3
"""Convert geNomad (ViWrap) virus summary into CAMI profile format."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

if __package__ is None or __package__ == "":  # pragma: no cover
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

RANK_ORDER = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


def _parse_taxonomy(taxonomy: str) -> Dict[str, str]:
    lineage = [part.strip() for part in taxonomy.split(";")]
    rank_map: Dict[str, str] = {}
    for rank, value in zip(RANK_ORDER, lineage):
        if not value or value in {"NA", "Unclassified"}:
            continue
        rank_map[rank] = value
    return rank_map


def load_virus_summary(path: Path, score_cutoff: float) -> List[Tuple[str, Dict[str, str]]]:
    entries: List[Tuple[str, Dict[str, str]]] = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            seq_name = (row.get("seq_name") or "").strip()
            if not seq_name:
                continue
            try:
                score = float(row.get("virus_score") or "0")
            except ValueError:
                score = 0.0
            if score < score_cutoff:
                continue
            taxonomy = (row.get("taxonomy") or "").strip()
            if not taxonomy:
                continue
            rank_map = _parse_taxonomy(taxonomy)
            if rank_map:
                entries.append((seq_name, rank_map))
    return entries


def build_profile(
    entries: List[Tuple[str, Dict[str, str]]],
    sample_id: str,
    tool: str,
    taxdb: str,
    profile_path: Path,
    classified_path: Path,
) -> None:
    if not entries:
        write_cami_profile([], str(profile_path), sample_id, tool)
        classified_path.unlink(missing_ok=True)
        return

    names = {rank_map[next(reversed(rank_map))] for _, rank_map in entries if rank_map}
    name_to_taxid = taxonkit_name2taxid(names, taxdb)

    counts: Dict[Tuple[str, str], int] = {}
    classified_rows: List[Tuple[str, str]] = []
    for seq_name, rank_map in entries:
        for rank in reversed(RANK_ORDER):
            name = rank_map.get(rank)
            if not name:
                continue
            hit = name_to_taxid.get(name)
            if not hit:
                continue
            taxid, actual_rank = hit
            counts[(taxid, actual_rank or rank)] = counts.get((taxid, actual_rank or rank), 0) + 1
            classified_rows.append((seq_name, taxid))
            break

    if not counts:
        write_cami_profile([], str(profile_path), sample_id, tool)
        classified_path.unlink(missing_ok=True)
        return

    taxids = [taxid for taxid, _ in counts.keys()]
    taxpaths = taxonkit_taxpath(taxids, taxdb)

    total = sum(counts.values())
    cami_rows: List[Dict[str, object]] = []
    for (taxid, rank), count in counts.items():
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
        for seq_name, taxid in classified_rows:
            writer.writerow([seq_name, taxid])


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert geNomad virus summary to CAMI format.")
    parser.add_argument("--input", required=True, help="Path to <sample>_summary/<sample>_virus_summary.tsv")
    parser.add_argument("--out", required=True, help="Output CAMI profile TSV.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--tool", default="viwrap", help="Tool identifier.")
    parser.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="TaxonKit database directory.")
    parser.add_argument("--score-cutoff", type=float, default=0.5, help="Minimum virus_score to retain (default 0.5).")
    parser.add_argument("--classified-out", default="", help="Optional path for classified_sequences.tsv.")
    args = parser.parse_args()

    summary_entries = load_virus_summary(Path(args.input), args.score_cutoff)
    classified_path = Path(args.classified_out) if args.classified_out else Path(args.out).with_name("classified_sequences.tsv")
    build_profile(summary_entries, args.sample_id, args.tool, args.taxdb, Path(args.out), classified_path)


if __name__ == "__main__":
    main()
