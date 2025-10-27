#!/usr/bin/env python3
"""Convert CAMITAX classification summary into CAMI profile format."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict, List

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
    sys.path.append(os.path.dirname(__file__))
    from common import RANKS, rollup_to_ancestors, taxonkit_taxpath, write_cami_profile  # type: ignore
else:  # pragma: no cover
    from .common import RANKS, rollup_to_ancestors, taxonkit_taxpath, write_cami_profile


def _load_camitax(path: str, sample_id: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(path, "r") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for raw in reader:
            if not raw:
                continue
            cells = [cell.strip() for cell in raw if cell is not None]
            if not cells:
                continue
            if cells[0].lower().startswith("genome"):
                # Header row
                continue
            genome_field = cells[0]
            if sample_id and sample_id not in genome_field:
                # Skip rows from other samples when using shared report files.
                continue
            taxid = cells[1] if len(cells) > 1 else ""
            if not taxid or not taxid.isdigit():
                continue
            rank = cells[3].lower() if len(cells) > 3 and cells[3] else "species"
            genome_name = genome_field
            if sample_id and genome_name.startswith(sample_id):
                genome_name = genome_name[len(sample_id) :]
                if genome_name.startswith(sample_id):
                    genome_name = genome_name[len(sample_id) :]
            rows.append(
                {
                    "genome": genome_name or sample_id,
                    "taxid": taxid,
                    "rank": rank or "species",
                }
            )
    return rows


def _pad_path(raw: str) -> List[str]:
    values = raw.split("|") if raw else []
    if len(values) < len(RANKS):
        values.extend(["NA"] * (len(RANKS) - len(values)))
    return values[: len(RANKS)]


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert CAMITAX camitax.tsv to CAMI profile format.")
    ap.add_argument("--input", required=True, help="camitax.tsv file produced by CAMITAX.")
    ap.add_argument("--out", required=True, help="Output CAMI TSV path.")
    ap.add_argument("--sample-id", required=True, help="Sample identifier.")
    ap.add_argument("--tool", default="camitax", help="Tool identifier for metadata.")
    ap.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="Optional TaxonKit data directory.")
    args = ap.parse_args()

    records = _load_camitax(args.input, args.sample_id)
    if not records:
        write_cami_profile([], args.out, args.sample_id, args.tool)
        return

    taxids = [entry["taxid"] for entry in records if entry["taxid"].isdigit()]
    taxid_to_paths = taxonkit_taxpath(taxids, args.taxdb)
    weight = 100.0 / len(records)

    cami_rows: List[Dict[str, object]] = []
    for entry in records:
        taxid = entry["taxid"]
        if not taxid.isdigit():
            continue
        taxpath_ids = ["NA"] * len(RANKS)
        taxpath_names = ["NA"] * len(RANKS)
        if taxid in taxid_to_paths:
            ids_raw, names_raw = taxid_to_paths[taxid]
            taxpath_ids = _pad_path(ids_raw)
            taxpath_names = _pad_path(names_raw)

        cami_rows.append(
            {
                "taxid": taxid,
                "rank": entry["rank"],
                "taxpath": taxpath_ids,
                "taxpathsn": taxpath_names,
                "percentage": weight,
            }
        )

    cami_rows = rollup_to_ancestors(cami_rows)
    write_cami_profile(cami_rows, args.out, args.sample_id, args.tool, normalise=True)


if __name__ == "__main__":
    main()
