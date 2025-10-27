#!/usr/bin/env python3
"""Convert Centrifuge kreport-style output into CAMI format."""

from __future__ import annotations

import argparse
import os
import sys

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
    sys.path.append(os.path.dirname(__file__))
    from kreport import parse_kreport  # type: ignore
    from common import rollup_to_ancestors, taxonkit_taxpath, RANKS, write_cami_profile  # type: ignore
else:  # pragma: no cover
    from .kreport import parse_kreport
    from .common import rollup_to_ancestors, taxonkit_taxpath, RANKS, write_cami_profile


def _parse_standard_report(path: str, taxdb: str) -> list[dict[str, object]]:
    rows = []
    taxids = []
    abundances = {}
    with open(path, "r", newline="") as handle:
        header = handle.readline()
        if not header:
            return []
        parts = header.strip().split("\t")
        if "taxID" not in parts or "abundance" not in parts:
            return []
        idx_taxid = parts.index("taxID")
        idx_abund = parts.index("abundance")
        handle.seek(0)
        import csv

        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) <= max(idx_taxid, idx_abund):
                continue
            taxid = row[idx_taxid].strip()
            if not taxid or taxid in {"NA", "0"}:
                continue
            try:
                abundance = float(row[idx_abund]) * 100.0
            except ValueError:
                continue
            if abundance <= 0:
                continue
            abundances[taxid] = abundances.get(taxid, 0.0) + abundance
            taxids.append(taxid)
    if not abundances:
        return []
    taxpaths = taxonkit_taxpath(list(abundances.keys()), taxdb)
    for taxid, percentage in abundances.items():
        ids_str, names_str = taxpaths.get(taxid, ("|".join(["NA"] * len(RANKS)), "|".join(["NA"] * len(RANKS))))
        id_vec = ids_str.split("|")
        name_vec = names_str.split("|")
        if len(id_vec) < len(RANKS):
            id_vec += ["NA"] * (len(RANKS) - len(id_vec))
        if len(name_vec) < len(RANKS):
            name_vec += ["NA"] * (len(RANKS) - len(name_vec))
        rank = "species"
        for idx in range(len(RANKS) - 1, -1, -1):
            if idx < len(id_vec) and id_vec[idx] not in {"", "NA"}:
                rank = RANKS[idx]
                break
        rows.append(
            {
                "taxid": taxid,
                "rank": rank,
                "taxpath": id_vec[: len(RANKS)],
                "taxpathsn": name_vec[: len(RANKS)],
                "percentage": percentage,
            }
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert Centrifuge reports to CAMI format.")
    ap.add_argument("--report", required=True, help="Path to centrifuge-kreport output.")
    ap.add_argument("--out", required=True, help="Output CAMI TSV path.")
    ap.add_argument("--sample-id", required=True, help="Sample identifier.")
    ap.add_argument("--tool", default="centrifuge", help="Tool identifier.")
    ap.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="TaxonKit database directory for abundance reports.")
    args = ap.parse_args()

    rows = parse_kreport(args.report)
    if not rows:
        rows = _parse_standard_report(args.report, args.taxdb)
    rows = rollup_to_ancestors(rows)
    write_cami_profile(rows, args.out, sample_id=args.sample_id, tool_name=args.tool, normalise=True)


if __name__ == "__main__":
    main()
