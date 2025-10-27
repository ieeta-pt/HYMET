#!/usr/bin/env python3
"""Convert sourmash gather CSV output into CAMI profile format."""

from __future__ import annotations

import argparse
import csv
import os
import sys
import re
from collections import defaultdict
from typing import Dict, Iterable, List

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
    sys.path.append(os.path.dirname(__file__))
    from common import RANKS, normalise_rows, rollup_to_ancestors, taxonkit_taxpath, write_cami_profile  # type: ignore
else:  # pragma: no cover
    from .common import RANKS, normalise_rows, rollup_to_ancestors, taxonkit_taxpath, write_cami_profile


def load_seqid_map(path: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if not path or not os.path.exists(path):
        raise FileNotFoundError(f"seqid2taxid map not found: {path}")
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            seqid, taxid = parts[0].strip(), parts[1].strip()
            if not seqid or not taxid:
                continue
            mapping.setdefault(seqid, taxid)
            if "." in seqid:
                mapping.setdefault(seqid.split(".", 1)[0], taxid)
    return mapping


def lookup_taxid(name: str, seqmap: Dict[str, str]) -> str | None:
    candidates: List[str] = []
    if not name:
        return None
    cleaned = name.strip()
    if not cleaned:
        return None
    # direct match and whitespace/pipe delimiters
    tokens = re.split(r"[\s\|,;]+", cleaned)
    candidates.extend(tokens)
    # also consider removing trailing descriptions
    candidates.append(cleaned.split()[0])
    for cand in candidates:
        cand = cand.strip()
        if not cand:
            continue
        if cand in seqmap:
            return seqmap[cand]
        if "." in cand:
            base = cand.split(".", 1)[0]
            if base in seqmap:
                return seqmap[base]
    return None


def gather_rows(gather_csv: str, seqmap: Dict[str, str]) -> Dict[str, float]:
    totals: Dict[str, float] = defaultdict(float)
    with open(gather_csv, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        frac_keys = [
            "f_unique_to_query",
            "fraction_unique_to_query",
            "unique_fraction",
        ]
        name_keys = ["name", "match_name", "filename"]
        for row in reader:
            if not row:
                continue
            frac = None
            for key in frac_keys:
                if key in row and row[key]:
                    try:
                        frac = float(row[key])
                        break
                    except ValueError:
                        continue
            if frac is None or frac <= 0.0:
                continue
            name_val = None
            for key in name_keys:
                if key in row and row[key]:
                    name_val = row[key]
                    break
            taxid = lookup_taxid(name_val or "", seqmap)
            if not taxid:
                continue
            totals[taxid] += frac * 100.0
    return totals


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert sourmash gather output to CAMI profile format.")
    ap.add_argument("--gather", required=True, help="sourmash gather CSV output")
    ap.add_argument("--seqmap", required=True, help="seqid2taxid map (from make_seqid_map.py)")
    ap.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="TaxonKit database directory")
    ap.add_argument("--out", required=True, help="Output CAMI profile path")
    ap.add_argument("--sample-id", required=True, help="Sample identifier")
    ap.add_argument("--tool", default="sourmash_gather", help="Tool name for CAMI header")
    args = ap.parse_args()

    seqmap = load_seqid_map(args.seqmap)
    totals = gather_rows(args.gather, seqmap)
    if not totals:
        write_cami_profile([], args.out, args.sample_id, args.tool, normalise=False)
        return

    tax_paths = taxonkit_taxpath(totals.keys(), args.taxdb)

    rows = []
    for taxid, percentage in totals.items():
        ids_str, names_str = tax_paths.get(taxid, ("|".join(["NA"] * len(RANKS)), "|".join(["NA"] * len(RANKS))))
        tax_ids = ids_str.split("|")
        tax_names = names_str.split("|")
        rank = "species"
        for idx in range(len(RANKS) - 1, -1, -1):
            if idx < len(tax_ids) and tax_ids[idx] not in {"", "NA"}:
                rank = RANKS[idx]
                break
        rows.append(
            {
                "taxid": taxid,
                "rank": rank,
                "taxpath": tax_ids,
                "taxpathsn": tax_names,
                "percentage": percentage,
            }
        )

    rows = normalise_rows(rows)
    rows = rollup_to_ancestors(rows)
    write_cami_profile(rows, args.out, args.sample_id, args.tool, normalise=True)


if __name__ == "__main__":
    main()
