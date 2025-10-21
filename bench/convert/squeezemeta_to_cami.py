#!/usr/bin/env python3
"""Convert SqueezeMeta contig taxonomy outputs into CAMI profile format."""

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
        taxonkit_name2taxid,
        taxonkit_taxpath,
        write_cami_profile,
    )
else:  # pragma: no cover
    from .common import RANKS, taxonkit_name2taxid, taxonkit_taxpath, write_cami_profile

RANK_KEYS = {
    "superkingdom": {"superkingdom", "domain", "kingdom"},
    "phylum": {"phylum"},
    "class": {"class"},
    "order": {"order"},
    "family": {"family"},
    "genus": {"genus"},
    "species": {"species"},
}


def _detect_columns(header: List[str]) -> Dict[str, int]:
    indices: Dict[str, int] = {}
    lowered = [h.strip().lower() for h in header]
    for rank, candidates in RANK_KEYS.items():
        for candidate in candidates:
            if candidate in lowered:
                indices[rank] = lowered.index(candidate)
                break
    return indices


def load_taxonomy(path: Path) -> Dict[str, Dict[str, str]]:
    taxonomy: Dict[str, Dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if not header:
            return taxonomy
        indices = _detect_columns(header)
        try:
            contig_idx = header.index("contig")
        except ValueError:
            contig_idx = 0
        for row in reader:
            if len(row) <= contig_idx:
                continue
            contig_id = row[contig_idx].strip()
            if not contig_id:
                continue
            rank_map: Dict[str, str] = {}
            for rank, idx in indices.items():
                if len(row) <= idx:
                    continue
                value = row[idx].strip()
                if not value or value in {"NA", "Unclassified"}:
                    continue
                rank_map[rank] = value
            if rank_map:
                taxonomy[contig_id] = rank_map
    return taxonomy


def build_profile(
    taxonomy: Dict[str, Dict[str, str]],
    sample_id: str,
    tool: str,
    taxdb: str,
    profile_path: Path,
    classified_path: Path,
) -> None:
    if not taxonomy:
        write_cami_profile([], str(profile_path), sample_id, tool)
        classified_path.unlink(missing_ok=True)
        return

    names = {name for rank_map in taxonomy.values() for name in rank_map.values() if name}
    name_to_taxid = taxonkit_name2taxid(names, taxdb)

    counts: Dict[Tuple[str, str], int] = {}
    classified_rows: List[Tuple[str, str]] = []
    for contig, rank_map in taxonomy.items():
        for rank in ("species", "genus", "family", "order", "class", "phylum", "superkingdom"):
            name = rank_map.get(rank)
            if not name:
                continue
            hit = name_to_taxid.get(name)
            if not hit:
                continue
            taxid, actual_rank = hit
            counts[(taxid, actual_rank or rank)] = counts.get((taxid, actual_rank or rank), 0) + 1
            classified_rows.append((contig, taxid))
            break

    if not counts:
        write_cami_profile([], str(profile_path), sample_id, tool)
        classified_path.unlink(missing_ok=True)
        return

    taxids = [taxid for taxid, _ in counts.keys()]
    taxpaths = taxonkit_taxpath(taxids, taxdb)

    total = sum(counts.values())
    cami_rows = []
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

    write_cami_profile(cami_rows, str(profile_path), sample_id, tool, normalise=False)

    with classified_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Query", "TaxID"])
        for contig, taxid in classified_rows:
            writer.writerow([contig, taxid])


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert SqueezeMeta contig taxonomy to CAMI format.")
    parser.add_argument("--input", required=True, help="Path to contig_taxonomy.summary (TSV).")
    parser.add_argument("--out", required=True, help="Output CAMI profile TSV.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--tool", default="squeezemeta", help="Tool identifier for metadata.")
    parser.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="TaxonKit database directory.")
    parser.add_argument("--classified-out", default="", help="Optional path for classified_sequences.tsv.")
    args = parser.parse_args()

    taxonomy = load_taxonomy(Path(args.input))
    classified_path = Path(args.classified_out) if args.classified_out else Path(args.out).with_name("classified_sequences.tsv")
    build_profile(taxonomy, args.sample_id, args.tool, args.taxdb, Path(args.out), classified_path)


if __name__ == "__main__":
    main()
