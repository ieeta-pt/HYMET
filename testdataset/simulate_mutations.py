#!/usr/bin/env python3
"""Generate mutated FASTA sequences for HYMET benchmarking."""

from __future__ import annotations

import argparse
import random
from pathlib import Path
from typing import Iterable, List

DNA_BASES = ("A", "C", "G", "T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Introduce substitutions and short indels into FASTA records."
    )
    parser.add_argument("--fasta", required=True, type=Path, help="Input FASTA.")
    parser.add_argument("--output", required=True, type=Path, help="Output FASTA.")
    parser.add_argument(
        "--sub-rate",
        type=float,
        default=0.1,
        help="Per-base substitution probability (default: 0.1).",
    )
    parser.add_argument(
        "--indel-rate",
        type=float,
        default=0.0,
        help="Per-base chance of an insertion or deletion event (default: 0.0).",
    )
    parser.add_argument(
        "--max-indel-length",
        type=int,
        default=3,
        help="Maximum length (bp) for insertions/deletions (default: 3).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed to ensure reproducibility (default: 42).",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> Iterable[tuple[str, str]]:
    header = None
    seq_chunks: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line.upper())
        if header is not None:
            yield header, "".join(seq_chunks)


def mutate_sequence(
    sequence: str,
    rng: random.Random,
    sub_rate: float,
    indel_rate: float,
    max_indel: int,
) -> str:
    mutated: List[str] = []
    i = 0
    while i < len(sequence):
        base = sequence[i]
        if base not in DNA_BASES:
            mutated.append(base)
            i += 1
            continue

        # deletion event
        if indel_rate > 0 and rng.random() < indel_rate / 2:
            del_len = rng.randint(1, max_indel)
            i += del_len
            continue

        # insertion event
        if indel_rate > 0 and rng.random() < indel_rate / 2:
            ins_len = rng.randint(1, max_indel)
            insert = "".join(rng.choice(DNA_BASES) for _ in range(ins_len))
            mutated.append(insert)

        # substitution
        if rng.random() < sub_rate:
            choices = [b for b in DNA_BASES if b != base]
            mutated.append(rng.choice(choices))
        else:
            mutated.append(base)

        i += 1
    return "".join(mutated)


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as out:
        for header, seq in read_fasta(args.fasta):
            mutated = mutate_sequence(
                seq,
                rng,
                sub_rate=args.sub_rate,
                indel_rate=args.indel_rate,
                max_indel=max(1, args.max_indel_length),
            )
            out.write(f"{header}\n")
            for i in range(0, len(mutated), 80):
                out.write(mutated[i : i + 80] + "\n")


if __name__ == "__main__":
    main()
