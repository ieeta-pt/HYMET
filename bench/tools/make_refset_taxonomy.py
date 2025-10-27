#!/usr/bin/env python3
"""
Build a lightweight taxonomy lookup table for a reference FASTA.

Given a FASTA file used by the CAMI benchmark harness, gather the canonical
names from the sequence headers, resolve them to TaxIDs via taxonkit, and emit
an identifiers table suitable for HYMET (GCF / TaxID / Identifiers).
"""

from __future__ import annotations

import argparse
import subprocess
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="Reference FASTA used by the benchmark (e.g. bench/refsets/combined_subset.fasta).",
    )
    parser.add_argument(
        "--taxonkit-db",
        required=True,
        type=Path,
        help="Directory containing taxonkit database files (typically HYMET/taxonomy_files).",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Destination TSV path (e.g. data/detailed_taxonomy.tsv).",
    )
    parser.add_argument(
        "--taxid-prefix",
        default="FAKE",
        help="Prefix for synthetic accession IDs in the output table (default: FAKE).",
    )
    return parser.parse_args()


def canonical_name(header: str) -> str:
    """
    Extract a moderately stable name from a FASTA header.
    Mimic the behaviour previously embedded in the benchmark instructions.
    """
    tokens = header.split()
    if not tokens:
        return ""
    seq_id = tokens[0]
    desc_tokens = tokens[1:]
    if desc_tokens:
        clean = desc_tokens[0].split(",", 1)[0]
        words = [tok for tok in clean.split() if tok]
        if len(words) >= 2:
            return " ".join(words[:2])
        if words:
            return words[0]
    return seq_id


def gather_names(fasta_path: Path) -> OrderedDict[str, List[str]]:
    mapping: OrderedDict[str, List[str]] = OrderedDict()
    with fasta_path.open("r") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            rest = line[1:].strip()
            if not rest:
                continue
            seq_id = rest.split()[0]
            name = canonical_name(rest) or seq_id
            mapping.setdefault(name, []).append(seq_id)
    return mapping


def resolve_taxids(names: Iterable[str], taxonkit_db: Path) -> Dict[str, str]:
    """
    Resolve names to TaxIDs using taxonkit name2taxid --show-rank.
    Returns a dict mapping input name -> taxid (string).
    """
    names = [name for name in names if name]
    if not names:
        return {}
    cmd = [
        "taxonkit",
        "name2taxid",
        "--data-dir",
        str(taxonkit_db),
        "--show-rank",
    ]
    proc = subprocess.run(
        cmd,
        input="\n".join(names) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    mapping: Dict[str, str] = {}
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        name, taxid = parts[0], parts[1]
        mapping[name] = taxid or "0"
    return mapping


def main() -> None:
    args = parse_args()
    if not args.fasta.exists():
        raise SystemExit(f"FASTA not found: {args.fasta}")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    names_to_ids = gather_names(args.fasta)
    if not names_to_ids:
        raise SystemExit(f"No FASTA headers found in {args.fasta}")
    name2tax = resolve_taxids(names_to_ids.keys(), args.taxonkit_db)

    with args.output.open("w") as out:
        out.write("GCF\tTaxID\tIdentifiers\n")
        for name, seq_ids in names_to_ids.items():
            taxid = name2tax.get(name, "0")
            accession = f"{args.taxid_prefix}_{taxid or '0'}"
            out.write(f"{accession}\t{taxid}\t{';'.join(seq_ids)}\n")
    print(f"Wrote {args.output} (entries: {len(names_to_ids)})")


if __name__ == "__main__":
    main()
