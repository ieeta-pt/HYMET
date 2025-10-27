#!/usr/bin/env python3
"""Convert MegaPath-Nano microbe statistics into CAMI profile format."""

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
        taxonkit_taxpath,
        write_cami_profile,
    )
else:  # pragma: no cover
    from .common import RANKS, rollup_to_ancestors, taxonkit_taxpath, write_cami_profile


def _parse_float(value: str | None) -> float:
    if value is None:
        return 0.0
    value = value.strip()
    if not value:
        return 0.0
    try:
        return float(value)
    except ValueError:
        try:
            return float(value.replace(",", ""))
        except ValueError:
            return 0.0


def _normalise_taxid(value: str | None) -> str:
    if value is None:
        return ""
    token = value.strip()
    if not token or token == "0":
        return ""
    if token.startswith("taxid|"):
        token = token.split("|", 1)[-1]
    return token


def _select_taxon(row: Dict[str, str]) -> Tuple[str, str, str]:
    candidates = (
        ("species", row.get("species_tax_id"), row.get("species_tax_name")),
        ("species", row.get("tax_id"), row.get("tax_name")),
        ("genus", row.get("genus_tax_id"), row.get("genus_tax_name")),
    )
    for rank, taxid, name in candidates:
        taxid_norm = _normalise_taxid(taxid)
        if taxid_norm:
            return taxid_norm, (name or "").strip(), rank
    return "", "", ""


def _load_mapping(path: Path | None) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if path is None or not path.is_file():
        return mapping
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for idx, row in enumerate(reader):
            if idx == 0 and row and row[0] in {"read_id", "sequence_id"}:
                continue
            if len(row) < 2:
                continue
            read_id = row[0].strip()
            contig = row[1].strip()
            if read_id and contig:
                mapping[read_id] = contig
    return mapping


def _read_microbe_sequence(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return rows
        for row in reader:
            rows.append({k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()})
    return rows


def _read_microbe_assembly(path: Path | None) -> Dict[str, Dict[str, str]]:
    if path is None or not path.is_file():
        return {}
    assemblies: Dict[str, Dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            assembly_id = (row.get("assembly_id") or "").strip()
            if assembly_id and assembly_id not in assemblies:
                assemblies[assembly_id] = {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
    return assemblies


def _aggregate_contig_calls(
    seq_rows: Iterable[Dict[str, str]],
    id_map: Dict[str, str],
) -> Dict[str, str]:
    votes: DefaultDict[str, DefaultDict[str, float]] = defaultdict(lambda: defaultdict(float))
    for row in seq_rows:
        read_id = (row.get("sequence_id") or "").strip()
        if not read_id:
            continue
        contig = id_map.get(read_id)
        if not contig:
            contig = read_id.split("|", 1)[0]
        taxid, _, _ = _select_taxon(row)
        if not taxid:
            continue
        weight = _parse_float(row.get("sequence_total_aligned_bp"))
        if weight <= 0.0:
            weight = _parse_float(row.get("sequence_length"))
        if weight <= 0.0:
            weight = 1.0
        votes[contig][taxid] += weight

    assignments: Dict[str, str] = {}
    for contig, weight_map in votes.items():
        if not weight_map:
            continue
        best_taxid, _ = max(weight_map.items(), key=lambda item: (item[1], item[0]))
        assignments[contig] = best_taxid
    return assignments


def _build_profile_rows(
    assemblies: Dict[str, Dict[str, str]],
    seq_rows: Iterable[Dict[str, str]],
) -> List[Tuple[str, str, float]]:
    records: Dict[str, Tuple[str, str, float]] = {}
    for assembly in assemblies.values():
        taxid, _, rank = _select_taxon(assembly)
        if not taxid:
            continue
        weight = _parse_float(assembly.get("adjusted_total_aligned_bp"))
        if weight <= 0.0:
            weight = _parse_float(assembly.get("assembly_adjusted_total_aligned_bp"))
        if weight <= 0.0:
            weight = _parse_float(assembly.get("pre_total_aligned_bp"))
        if weight <= 0.0:
            weight = _parse_float(assembly.get("assembly_length"))
        if weight <= 0.0:
            weight = 1.0
        key = taxid
        existing = records.get(key)
        if existing:
            records[key] = (taxid, rank, existing[2] + weight)
        else:
            records[key] = (taxid, rank, weight)

    if not records and seq_rows:
        fallback: DefaultDict[str, float] = defaultdict(float)
        for row in seq_rows:
            taxid, _, rank = _select_taxon(row)
            if not taxid:
                continue
            weight = _parse_float(row.get("sequence_total_aligned_bp"))
            if weight <= 0.0:
                weight = _parse_float(row.get("sequence_length"))
            if weight <= 0.0:
                weight = 1.0
            fallback[(taxid, rank)] += weight
        for (taxid, rank), weight in fallback.items():
            records[taxid] = (taxid, rank, weight)

    return list(records.values())


def _prepare_cami_rows(
    entries: List[Tuple[str, str, float]],
    taxdb: str,
) -> List[Dict[str, object]]:
    if not entries:
        return []
    unique_taxids = [taxid for taxid, _rank, _weight in entries if taxid and taxid != "NA"]
    tax_paths = taxonkit_taxpath(unique_taxids, taxdb)
    default_path = "|".join(["NA"] * len(RANKS))
    cami_rows: List[Dict[str, object]] = []
    for taxid, rank, weight in entries:
        if weight <= 0.0 or not taxid:
            continue
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
                "rank": rank or "species",
                "taxpath": taxpath[: len(RANKS)],
                "taxpathsn": names[: len(RANKS)],
                "percentage": weight,
            }
        )
    return cami_rows


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert MegaPath-Nano outputs to CAMI profile format.")
    ap.add_argument("--input", required=True, help="Path to *.microbe_stat_by_sequence_id_assembly_info file.")
    ap.add_argument("--assembly-stat", help="Optional path to *.microbe_stat summary table.")
    ap.add_argument("--out", required=True, help="Destination CAMI profile TSV.")
    ap.add_argument("--sample-id", required=True, help="Sample identifier.")
    ap.add_argument("--tool", default="megapath_nano", help="Tool identifier.")
    ap.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="TaxonKit database directory.")
    ap.add_argument("--classified-out", help="Optional path for classified_sequences.tsv")
    ap.add_argument("--id-map", help="Optional read→contig mapping file.")
    args = ap.parse_args()

    seq_rows = _read_microbe_sequence(Path(args.input))
    assemblies = _read_microbe_assembly(Path(args.assembly_stat)) if args.assembly_stat else {}

    entries = _build_profile_rows(assemblies, seq_rows)
    cami_rows = _prepare_cami_rows(entries, args.taxdb)
    cami_rows = rollup_to_ancestors(cami_rows)
    write_cami_profile(cami_rows, args.out, args.sample_id, args.tool, normalise=True)

    classified_path = Path(args.classified_out) if args.classified_out else None
    if classified_path:
        mapping = _load_mapping(Path(args.id_map)) if args.id_map else {}
        assignments = _aggregate_contig_calls(seq_rows, mapping)
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
