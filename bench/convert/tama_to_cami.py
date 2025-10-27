#!/usr/bin/env python3
"""Convert TAMA abundance outputs into CAMI-compatible taxonomic profiles."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
    sys.path.append(os.path.dirname(__file__))
    from common import (  # type: ignore
        RANKS,
        default_taxpath,
        rollup_to_ancestors,
        taxonkit_taxpath,
        write_cami_profile,
    )
else:  # pragma: no cover
    from .common import RANKS, default_taxpath, rollup_to_ancestors, taxonkit_taxpath, write_cami_profile


def parse_abundance(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return rows
        header_map = {name.strip().lower(): idx for idx, name in enumerate(header)}
        name_idx = header_map.get("scientific name")
        taxid_idx = header_map.get("taxonomy id") or header_map.get("taxon id")
        abundance_idx = header_map.get("abundance")
        if abundance_idx is None:
            # Some versions emit "relative abundance"
            abundance_idx = header_map.get("relative abundance")
        if abundance_idx is None or taxid_idx is None:
            return rows
        for line in reader:
            if not line:
                continue
            taxid = (line[taxid_idx].strip() if taxid_idx < len(line) else "").strip()
            name = (line[name_idx].strip() if name_idx is not None and name_idx < len(line) else "").strip()
            abundance_str = (line[abundance_idx].strip() if abundance_idx < len(line) else "").strip()
            if not abundance_str:
                continue
            try:
                abundance = float(abundance_str)
            except ValueError:
                continue
            rows.append(
                {
                    "taxid": taxid or "NA",
                    "name": name or "NA",
                    "abundance": abundance,
                }
            )
    return rows


def build_profile(
    records: Iterable[Dict[str, str]],
    sample_id: str,
    tool_name: str,
    out_path: Path,
    taxdb: str,
    rank: str,
) -> None:
    entries = list(records)
    if not entries:
        write_cami_profile([], str(out_path), sample_id, tool_name)
        return

    taxids = [entry["taxid"] for entry in entries if entry["taxid"] not in {"", "NA", "NaN"}]
    taxpaths: Dict[str, Tuple[str, str]] = taxonkit_taxpath(taxids, taxdb) if taxids else {}

    cami_rows: List[Dict[str, object]] = []

    for entry in entries:
        taxid = entry["taxid"] if entry["taxid"] not in {"", "NaN"} else "NA"
        abundance = float(entry.get("abundance", 0.0))
        ids_default, names_default = default_taxpath()
        taxpath_ids_list = list(ids_default)
        taxpath_names_list = list(names_default)
        if taxid != "NA" and taxid in taxpaths:
            taxpath_ids, taxpath_names = taxpaths[taxid]
            taxpath_ids_list = taxpath_ids.split("|")
            taxpath_names_list = taxpath_names.split("|")
        target_len = len(ids_default)
        if len(taxpath_ids_list) < target_len:
            taxpath_ids_list.extend(["NA"] * (target_len - len(taxpath_ids_list)))
        if len(taxpath_names_list) < target_len:
            taxpath_names_list.extend(["NA"] * (target_len - len(taxpath_names_list)))
        cami_rows.append(
            {
                "taxid": taxid,
                "rank": rank,
                "taxpath": taxpath_ids_list,
                "taxpathsn": taxpath_names_list,
                "percentage": 100.0 * abundance,
            }
        )
    cami_rows = rollup_to_ancestors(cami_rows)
    write_cami_profile(cami_rows, str(out_path), sample_id, tool_name, normalise=True)


def convert_classifications(read_path: Optional[Path], out_path: Path) -> None:
    if not read_path or not read_path.exists():
        return
    with read_path.open("r", encoding="utf-8", errors="ignore") as inp, out_path.open(
        "w", encoding="utf-8", newline=""
    ) as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["Query", "TaxID"])
        for line in inp:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            query = parts[0].strip()
            taxon_field = parts[1].strip()
            if not query or not taxon_field:
                continue
            taxon_candidates = [tok.strip() for tok in taxon_field.replace(";", ",").split(",") if tok.strip()]
            if not taxon_candidates:
                continue
            taxid = taxon_candidates[0]
            if taxid.lower() == "na":
                continue
            writer.writerow([query, taxid])


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert TAMA outputs into CAMI profile format.")
    parser.add_argument("--profile", required=True, help="Path to TAMA abundance_profile*.out file.")
    parser.add_argument("--out", required=True, help="Destination CAMI profile TSV.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--tool", default="tama", help="Tool identifier for CAMI metadata.")
    parser.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="Optional TaxonKit DB directory.")
    parser.add_argument("--rank", default="species", help="Taxonomic rank TAMA analysed (default: species).")
    parser.add_argument("--read-classi", default="", help="Optional path to read_classi*.out file.")
    parser.add_argument("--classified-out", default="", help="Optional path to classified_sequences.tsv output.")
    args = parser.parse_args()

    profile_path = Path(args.profile)
    if not profile_path.is_file():
        raise FileNotFoundError(f"TAMA abundance profile not found: {profile_path}")

    rows = parse_abundance(profile_path)
    out_path = Path(args.out)
    build_profile(rows, args.sample_id, args.tool, out_path, args.taxdb, args.rank)

    if args.classified_out:
        read_path = Path(args.read_classi) if args.read_classi else None
        convert_classifications(read_path, Path(args.classified_out))


if __name__ == "__main__":
    main()
