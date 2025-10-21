#!/usr/bin/env python3
"""Convert PhaBOX taxonomy predictions into CAMI profile format."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
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


def _read_id_map(path: Optional[Path]) -> Dict[str, str]:
    if not path:
        return {}
    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        header_lower = [h.lower() for h in header]
        try:
            new_idx = header_lower.index("new_id")
        except ValueError:
            new_idx = 0
        try:
            orig_idx = header_lower.index("original_id")
        except ValueError:
            if len(header_lower) > 1:
                orig_idx = 1
            else:
                orig_idx = 0
        for row in reader:
            if len(row) <= max(new_idx, orig_idx):
                continue
            new_id = row[new_idx].strip()
            orig_id = row[orig_idx].strip()
            if new_id and orig_id:
                mapping[new_id] = orig_id
    return mapping


def _parse_lineage(raw: str) -> RankPath:
    lineage: RankPath = {}
    for part in raw.split(";"):
        part = part.strip()
        if not part or ":" not in part:
            continue
        level, name = part.split(":", 1)
        level = level.strip().lower()
        name = name.strip()
        if not name or name.lower() in {"unclassified", "unassigned", "na"}:
            continue
        # normalize synonyms
        if level == "kingdom":
            level = "superkingdom"
        lineage[level] = name
    return lineage


def _detect_columns(header: List[str]) -> Tuple[int, int, Optional[int]]:
    lowered = [h.lower() for h in header]
    id_idx = 0
    lineage_idx = 1
    status_idx: Optional[int] = None
    candidates_id = ["contig_id", "contigid", "query", "contig", "accession"]
    for candidate in candidates_id:
        if candidate in lowered:
            id_idx = lowered.index(candidate)
            break
    lineage_candidates = ["lineage", "taxonomy", "taxon", "best_match", "tax_info"]
    for candidate in lineage_candidates:
        if candidate in lowered:
            lineage_idx = lowered.index(candidate)
            break
    status_candidates = ["status", "prediction", "note"]
    for candidate in status_candidates:
        if candidate in lowered:
            status_idx = lowered.index(candidate)
            break
    return id_idx, lineage_idx, status_idx


def load_phabox(path: Path) -> Dict[str, RankPath]:
    predictions: Dict[str, RankPath] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if not header:
            return predictions
        id_idx, lineage_idx, status_idx = _detect_columns(header)
        for row in reader:
            if len(row) <= max(id_idx, lineage_idx):
                continue
            contig_id = row[id_idx].strip()
            lineage_raw = row[lineage_idx].strip()
            if not contig_id or not lineage_raw:
                continue
            if status_idx is not None and status_idx < len(row):
                status = row[status_idx].strip().lower()
                if status in {"filtered", "unpredicted", "unassigned"}:
                    continue
            lineage = _parse_lineage(lineage_raw)
            if not lineage:
                continue
            predictions[contig_id] = lineage
    return predictions


def _choose_rank(lineage: RankPath) -> Optional[Tuple[str, str]]:
    for rank in reversed(RANK_ORDER):  # prefer the most specific available
        if rank in lineage:
            return lineage[rank], rank
    return None


def build_profile(
    predictions: Dict[str, RankPath],
    id_map: Dict[str, str],
    sample_id: str,
    tool_name: str,
    taxdb: str,
    out_path: Path,
    classified_out: Path,
) -> None:
    contig_to_name_rank: Dict[str, Tuple[str, str]] = {}
    names = set()
    for contig_id, lineage in predictions.items():
        chosen = _choose_rank(lineage)
        if not chosen:
            continue
        name, rank = chosen
        names.add(name)
        real_id = id_map.get(contig_id, contig_id)
        contig_to_name_rank[real_id] = (name, rank)

    if not contig_to_name_rank:
        write_cami_profile([], str(out_path), sample_id, tool_name)
        classified_out.unlink(missing_ok=True)
        return

    name_to_tax = taxonkit_name2taxid(names, taxdb)
    counts: Counter[Tuple[str, str]] = Counter()
    classified_rows: List[Tuple[str, str]] = []
    for contig_id, (name, rank) in contig_to_name_rank.items():
        hit = name_to_tax.get(name)
        if not hit:
            continue
        taxid, actual_rank = hit
        counts[(taxid, actual_rank or rank)] += 1
        classified_rows.append((contig_id, taxid))

    if not counts:
        write_cami_profile([], str(out_path), sample_id, tool_name)
        classified_out.unlink(missing_ok=True)
        return

    taxid_to_paths = taxonkit_taxpath([taxid for taxid, _ in counts.keys()], taxdb)
    total = sum(counts.values())
    cami_rows = []
    for (taxid, rank), count in counts.items():
        ids_raw, names_raw = taxid_to_paths.get(taxid, ("|".join(["NA"] * len(RANKS)), "|".join(["NA"] * len(RANKS))))
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
    write_cami_profile(cami_rows, str(out_path), sample_id, tool_name, normalise=True)

    with classified_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Query", "TaxID"])
        for contig_id, taxid in classified_rows:
            writer.writerow([contig_id, taxid])


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert PhaBOX phagcn predictions to CAMI format.")
    parser.add_argument("--input", required=True, help="phagcn_prediction.tsv (or similar) produced by PhaBOX.")
    parser.add_argument("--out", required=True, help="Output CAMI profile TSV.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--tool", default="phabox", help="Tool identifier for metadata.")
    parser.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="TaxonKit database directory.")
    parser.add_argument("--id-map", default="", help="Optional TSV mapping new IDs to original contig IDs.")
    parser.add_argument("--classified-out", default="", help="Optional path for classified_sequences.tsv.")
    args = parser.parse_args()

    predictions = load_phabox(Path(args.input))
    id_map = _read_id_map(Path(args.id_map)) if args.id_map else {}
    classified_path = Path(args.classified_out) if args.classified_out else Path(args.out).with_name("classified_sequences.tsv")
    build_profile(predictions, id_map, args.sample_id, args.tool, args.taxdb, Path(args.out), classified_path)


if __name__ == "__main__":
    main()
