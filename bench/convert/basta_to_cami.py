#!/usr/bin/env python3
"""Convert BASTA taxonomy assignments into CAMI profile format."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
    sys.path.append(os.path.dirname(__file__))
    from common import (  # type: ignore
        RANKS,
        taxonkit_name2taxid,
        taxonkit_taxpath,
        write_cami_profile,
    )
else:  # pragma: no cover
    from .common import RANKS, taxonkit_name2taxid, taxonkit_taxpath, write_cami_profile


Assignment = Tuple[str, List[str], str]
CANONICAL_TAXA: Dict[str, Tuple[str, str]] = {
    "bacteria": ("2", "superkingdom"),
    "archaea": ("2157", "superkingdom"),
    "eukaryota": ("2759", "superkingdom"),
    "viruses": ("10239", "superkingdom"),
}


def parse_basta(path: Path) -> List[Assignment]:
    assignments: List[Assignment] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            query = row[0].strip()
            if not query:
                continue
            lineage_raw = row[1].strip() if len(row) > 1 else ""
            if not lineage_raw:
                continue
            lineage = [token.strip() for token in lineage_raw.split(";") if token.strip()]
            if not lineage:
                continue
            provided_taxid = ""
            for cell in row[2:]:
                value = cell.strip()
                if value.isdigit():
                    provided_taxid = value
                    break
            assignments.append((query, lineage, provided_taxid))
    return assignments


def build_profile(
    assignments_iter: Iterable[Assignment],
    sample_id: str,
    tool_name: str,
    taxdb: str,
    out_path: Path,
    classified_out: Path | None,
) -> None:
    assignments = list(assignments_iter)
    label_counts: Counter[Tuple[str, str]] = Counter()
    taxid_to_label: Dict[str, str] = {}
    classified_rows: List[Tuple[str, str]] = []

    lineage_names = {
        token.replace("_", " ")
        for _, lineage, _ in assignments
        for token in lineage
    }
    name_to_taxid = taxonkit_name2taxid(sorted(lineage_names), taxdb) if lineage_names else {}

    for query, lineage, provided_taxid in assignments:
        taxid = ""
        rank = ""

        chosen_lineage = lineage[-1]
        if provided_taxid:
            taxid = provided_taxid
            index = min(len(lineage), len(RANKS)) - 1
            rank = RANKS[index] if index >= 0 else "species"
        else:
            for name in reversed(lineage):
                normalized = name.replace("_", " ")
                canonical = CANONICAL_TAXA.get(normalized.lower())
                if canonical:
                    taxid, rank = canonical
                    break
                info = name_to_taxid.get(normalized)
                if info:
                    taxid, rank = info
                    break

        if not taxid:
            continue

        label_counts[(taxid, rank)] += 1
        taxid_to_label[taxid] = chosen_lineage
        classified_rows.append((query, taxid))

    if not label_counts:
        write_cami_profile([], str(out_path), sample_id, tool_name)
        if classified_out:
            classified_out.unlink(missing_ok=True)
        return

    taxid_paths = taxonkit_taxpath(taxid_to_label.keys(), taxdb)
    total = sum(label_counts.values())
    cami_rows = []
    for (taxid, rank), count in sorted(label_counts.items(), key=lambda kv: kv[0]):
        ids_raw, names_raw = taxid_paths.get(taxid, ("|".join(["NA"] * len(RANKS)), "|".join(["NA"] * len(RANKS))))
        ids_vec = ids_raw.split("|")
        names_vec = names_raw.split("|")
        if len(ids_vec) < len(RANKS):
            ids_vec += ["NA"] * (len(RANKS) - len(ids_vec))
        if len(names_vec) < len(RANKS):
            names_vec += ["NA"] * (len(RANKS) - len(names_vec))
        cami_rows.append(
            {
                "taxid": taxid,
                "rank": rank or "species",
                "taxpath": ids_vec[: len(RANKS)],
                "taxpathsn": names_vec[: len(RANKS)],
                "percentage": 100.0 * count / total if total else 0.0,
            }
        )

    write_cami_profile(cami_rows, str(out_path), sample_id, tool_name, normalise=False)
    if classified_out:
        with classified_out.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Query", "TaxID"])
            for query, taxid in classified_rows:
                writer.writerow([query, taxid])


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert BASTA outputs into CAMI format.")
    parser.add_argument("--input", required=True, help="BASTA taxonomy file (tab-separated).")
    parser.add_argument("--out", required=True, help="Path to write CAMI profile TSV.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--tool", default="basta", help="Tool identifier for metadata.")
    parser.add_argument("--taxdb", default=os.environ.get("TAXONKIT_DB", ""), help="TaxonKit database directory.")
    parser.add_argument("--classified-out", default="", help="Optional classified_sequences.tsv output path.")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        raise FileNotFoundError(f"BASTA taxonomy file not found: {input_path}")

    assignments = parse_basta(input_path)
    out_path = Path(args.out)
    classified_path = Path(args.classified_out) if args.classified_out else None
    build_profile(assignments, args.sample_id, args.tool, args.taxdb, out_path, classified_path)


if __name__ == "__main__":
    main()
