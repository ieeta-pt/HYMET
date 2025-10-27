#!/usr/bin/env python3
"""Extract a size-limited subset from a FASTA file.

To support benchmarking on machines with constrained disk or RAM,
this tool streams a multi-FASTA input and writes the first N sequences
and/or bases to an output path. The subset can be fed to the DB build
scripts via the REF_FASTA environment variable.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def iter_fasta(path: Path):
    name = None
    seq_lines = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_lines)
                name = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if name is not None:
            yield name, "".join(seq_lines)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Create a sequence-limited subset from a FASTA file."
    )
    ap.add_argument("--input", required=True, help="Input FASTA file.")
    ap.add_argument("--output", required=True, help="Output FASTA path.")
    ap.add_argument(
        "--max-seqs",
        type=int,
        default=1000,
        help="Maximum number of sequences to emit (default: 1000).",
    )
    ap.add_argument(
        "--max-bases",
        type=int,
        default=500_000_000,
        help="Maximum number of bases to emit across all sequences (default: 500M).",
    )
    args = ap.parse_args()

    inp = Path(args.input)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    emitted_seqs = 0
    emitted_bases = 0

    with out.open("w", encoding="utf-8") as handle:
        for header, seq in iter_fasta(inp):
            if emitted_seqs >= args.max_seqs or emitted_bases >= args.max_bases:
                break
            remaining_bases = args.max_bases - emitted_bases
            subseq = seq if len(seq) <= remaining_bases else seq[:remaining_bases]
            handle.write(f"{header}\n")
            for i in range(0, len(subseq), 80):
                handle.write(subseq[i : i + 80] + "\n")
            emitted_seqs += 1
            emitted_bases += len(subseq)
            if len(subseq) < len(seq):
                break

    print(
        f"[subset_fasta] wrote {emitted_seqs} sequences / {emitted_bases} bases to {out}"
    )


if __name__ == "__main__":
    main()
