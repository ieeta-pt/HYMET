#!/usr/bin/env python3
"""Create a seqid -> TaxID mapping using HYMET's detailed taxonomy file."""

from __future__ import annotations

import argparse
import csv
import re
from typing import Dict


def load_id_map(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            key = row[0].strip()
            tax = row[1].strip()
            if key and tax and key not in out:
                out[key] = tax
                if "." in key:
                    out.setdefault(key.split(".", 1)[0], tax)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate seqid2taxid map for Centrifuge/Ganon.")
    ap.add_argument("--fasta", required=True, help="Reference FASTA file.")
    ap.add_argument("--taxonomy-map", required=True, help="Detailed taxonomy map TSV (from build_id_map.py).")
    ap.add_argument("--out", required=True, help="Output path for seqid2taxid.map.")
    args = ap.parse_args()

    id_map = load_id_map(args.taxonomy_map)

    def lookup(token: str) -> str | None:
        if not token:
            return None
        if token in id_map:
            return id_map[token]
        if "." in token:
            base = token.split(".", 1)[0]
            if base in id_map:
                return id_map[base]
        return None

    missing = 0
    total = 0
    with open(args.fasta, "r") as fin, open(args.out, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                continue
            total += 1
            header = line[1:].strip()
            seq_id = header.split()[0]
            tax = lookup(seq_id)
            if not tax:
                tokens = re.split(r"[\\s\\|,;]+", header)
                for tok in tokens:
                    tax = lookup(tok)
                    if tax:
                        break
            if tax:
                fout.write(f"{seq_id}\t{tax}\n")
            else:
                missing += 1

    print(f"[INFO] seqid2taxid: mapped {total - missing} / {total} sequences (missing {missing}).")


if __name__ == "__main__":
    main()
