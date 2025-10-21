#!/usr/bin/env python3
"""Split a mutated combined FASTA into per-group FASTA files."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict

GROUPS = [
    "Viruses",
    "Archaea",
    "Bacteria",
    "Fungi",
    "Plants",
    "Protozoa",
    "Invertebrates",
    "Vertebrate Mammals",
    "Other Vertebrates",
]


def stream_fasta(path: Path):
    with path.open() as fh:
        header = None
        seq_lines = []
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line.strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            yield header, "".join(seq_lines)


def contig_id(header: str) -> str:
    return header[1:].split()[0]


def split_fasta(fasta: Path, metadata: Dict[str, dict], outdir: Path) -> None:
    handles = {}
    try:
        for group in GROUPS:
            path = outdir / f"{group.replace(' ', '_').lower()}.fna"
            path.parent.mkdir(parents=True, exist_ok=True)
            handles[group] = path.open("w")
        missing = 0
        for header, sequence in stream_fasta(fasta):
            cid = contig_id(header)
            info = metadata.get(cid) or metadata.get(cid.split(".", 1)[0])
            if not info:
                missing += 1
                continue
            group = info.get("group", "Protozoa")
            handle = handles.get(group)
            if handle is None:
                missing += 1
                continue
            handle.write(f"{header}\n")
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i : i + 80] + "\n")
        if missing:
            print(f"[split] warning: {missing} contigs missing metadata were skipped.")
    finally:
        for handle in handles.values():
            handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Split mutated FASTA by HYMET group")
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metadata = json.loads(args.metadata.read_text())
    args.outdir.mkdir(parents=True, exist_ok=True)
    split_fasta(args.fasta, metadata, args.outdir)


if __name__ == "__main__":
    main()
