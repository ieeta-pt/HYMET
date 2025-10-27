#!/usr/bin/env python3
"""
Populate a HYMET cache directory with a curated FASTA and taxonomy map.

This script copies a curated reference FASTA into an existing cache entry,
rebuilds the identifier -> TaxID table from a seqid→taxid mapping, and removes
any stale minimap2 indices so they will be regenerated on the next run.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache-dir",
        required=True,
        type=Path,
        help="Target cache directory (e.g. data/downloaded_genomes/cache_case/<sha1>).",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="Curated reference FASTA (plain or gzipped).",
    )
    parser.add_argument(
        "--seqmap",
        required=True,
        type=Path,
        help="Tab-separated file mapping sequence IDs to TaxIDs (two columns).",
    )
    parser.add_argument(
        "--taxid-prefix",
        default="Curated",
        help="Prefix for synthetic accession IDs in detailed_taxonomy.tsv (default: Curated).",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def copy_fasta(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open_text(src) as fin, dest.open("w") as fout:
        shutil.copyfileobj(fin, fout)


def load_seqmap(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            seq, taxid = row[0].strip(), row[1].strip()
            if seq and taxid:
                mapping[seq] = taxid
    if not mapping:
        raise ValueError(f"No sequence→taxid entries found in {path}")
    return mapping


def build_taxonomy(fasta: Path, seq2tax: dict[str, str], prefix: str) -> list[tuple[str, str, list[str]]]:
    buckets = defaultdict(list)
    with open_text(fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            seq_id = line[1:].split()[0]
            taxid = seq2tax.get(seq_id)
            if taxid:
                buckets[taxid].append(seq_id)
    if not buckets:
        raise ValueError("No FASTA headers matched seqmap entries; check inputs.")
    records = []
    def sort_key(value: str):
        try:
            return (0, int(value))
        except ValueError:
            return (1, value)

    for taxid in sorted(buckets.keys(), key=sort_key):
        accession = f"{prefix}_{taxid}"
        records.append((accession, taxid, buckets[taxid]))
    return records


def write_taxonomy(records: list[tuple[str, str, list[str]]], dest: Path) -> None:
    with dest.open("w") as out:
        out.write("GCF\tTaxID\tIdentifiers\n")
        for accession, taxid, identifiers in records:
            out.write(f"{accession}\t{taxid}\t{';'.join(identifiers)}\n")


def remove_indices(cache_dir: Path) -> None:
    for idx in cache_dir.glob("*.mmi"):
        idx.unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    cache_dir = args.cache_dir
    fasta_dest = cache_dir / "combined_genomes.fasta"
    taxonomy_dest = cache_dir / "detailed_taxonomy.tsv"

    cache_dir.mkdir(parents=True, exist_ok=True)
    seq2tax = load_seqmap(args.seqmap)
    copy_fasta(args.fasta, fasta_dest)
    taxonomy_records = build_taxonomy(fasta_dest, seq2tax, args.taxid_prefix)
    write_taxonomy(taxonomy_records, taxonomy_dest)
    remove_indices(cache_dir)

    print(f"Updated cache at {cache_dir}")
    print(f"- combined_genomes.fasta (source: {args.fasta})")
    print(f"- detailed_taxonomy.tsv (entries: {len(taxonomy_records)})")
    print("- reference.mmi removed (will be rebuilt on next run)")


if __name__ == "__main__":
    main()
