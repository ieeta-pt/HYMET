#!/usr/bin/env python3
"""Convert contig FASTA files into synthetic single-end reads for TAMA.

The TAMA workflow expects FASTQ reads produced by Kraken/Centrifuge/CLARK.
This helper slices each contig into fixed-size windows and emits high-quality
single-end reads so that assembled contigs can be analysed without needing the
original raw data.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Tuple


def iter_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    name = None
    seq_parts = []
    with path.open("r") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts).upper()
                name = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if name is not None:
            yield name, "".join(seq_parts).upper()


def write_reads(
    fasta: Path,
    fastq: Path,
    chunk_size: int = 250,
    min_chunk: int = 100,
) -> int:
    total = 0
    with fastq.open("w") as out:
        for contig_id, seq in iter_fasta(fasta):
            if not seq:
                continue
            seq = seq.replace("N", "A")
            if len(seq) <= chunk_size:
                total += 1
                qual = "I" * len(seq)
                out.write(f"@{contig_id}|0|{len(seq)}\n{seq}\n+\n{qual}\n")
                continue
            start = 0
            while start < len(seq):
                end = min(start + chunk_size, len(seq))
                chunk = seq[start:end]
                if len(chunk) < min_chunk and start != 0:
                    break
                total += 1
                qual = "I" * len(chunk)
                out.write(f"@{contig_id}|{start}|{end}\n{chunk}\n+\n{qual}\n")
                start += chunk_size
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", required=True, help="Input contigs FASTA")
    parser.add_argument("--out", required=True, help="Output FASTQ path")
    parser.add_argument("--chunk-size", type=int, default=250, help="Window size (bp)")
    parser.add_argument("--min-chunk", type=int, default=100, help="Minimum tail chunk size to emit")
    args = parser.parse_args()

    contigs = Path(args.contigs).resolve()
    out_fastq = Path(args.out).resolve()
    out_fastq.parent.mkdir(parents=True, exist_ok=True)

    total = write_reads(contigs, out_fastq, chunk_size=args.chunk_size, min_chunk=args.min_chunk)
    print(f"[contigs_to_reads] wrote {total} reads → {out_fastq}")


if __name__ == "__main__":
    main()
