#!/usr/bin/env python3
"""
limit_candidates.py

Post-process the Mash candidate list to (optionally) deduplicate at the species
level and cap the number of candidates retained for downstream reference
building. When species deduplication is enabled the script keeps the
highest-scoring assembly per species (based on Mash screen scores) and trims
the list to the requested maximum.

This helper is designed to be idempotent — it never mutates the input file in
place. The final candidate list is written to the path supplied via
`--output`.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import os
import pathlib
import sys
import tempfile
import urllib.request
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

SUMMARY_URLS = {
    "assembly_summary_refseq.txt": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt",
    "assembly_summary_genbank.txt": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt",
}

DEFAULT_MAX_CANDIDATES = 5000
REFRESH_DAYS = 14  # refresh assembly summaries every two weeks


@dataclass
class Candidate:
    name: str
    score: float
    species_key: str
    species_label: str
    order: int


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Limit Mash candidate genomes with optional species-level deduplication."
    )
    parser.add_argument("--selected", required=True, help="Input file with one candidate per line.")
    parser.add_argument("--output", required=True, help="Path to write the filtered candidate list.")
    parser.add_argument(
        "--score-file",
        action="append",
        default=[],
        dest="score_files",
        help="Tab-delimited Mash screen file (column 1 score, column 5 accession). "
        "May be supplied multiple times.",
    )
    parser.add_argument(
        "--assembly-dir",
        default=None,
        help="Directory holding assembly_summary_*.txt (downloaded if missing). "
        "Defaults to HYMET/data/downloaded_genomes/assembly_summaries relative to this script.",
    )
    parser.add_argument(
        "--max",
        type=int,
        default=DEFAULT_MAX_CANDIDATES,
        help=f"Maximum number of candidates to retain (default: {DEFAULT_MAX_CANDIDATES}).",
    )
    parser.add_argument(
        "--dedupe",
        action="store_true",
        help="Enable species-level deduplication before applying the candidate cap.",
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Optional path to append a log line summarising how many candidates were kept.",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Skip downloading assembly summaries if they are missing. Species deduplication "
        "will fall back to per-accession unique keys.",
    )
    return parser.parse_args(argv)


def read_candidates(path: pathlib.Path) -> List[str]:
    with path.open("r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def load_scores(files: Iterable[pathlib.Path]) -> Dict[str, float]:
    scores: Dict[str, float] = {}
    for file in files:
        if not file.exists():
            continue
        try:
            with file.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    candidate = parts[4].strip()
                    if not candidate:
                        continue
                    try:
                        score = float(parts[0])
                    except ValueError:
                        continue
                    # Keep the best score encountered for the candidate.
                    if candidate not in scores or score > scores[candidate]:
                        scores[candidate] = score
        except OSError:
            continue
    return scores


def needs_refresh(path: pathlib.Path) -> bool:
    if not path.exists():
        return True
    try:
        mtime = path.stat().st_mtime
    except OSError:
        return True
    age_days = ( _dt.datetime.now().timestamp() - mtime ) / 86400.0
    return age_days > REFRESH_DAYS


def ensure_assembly_summary(name: str, target_dir: pathlib.Path, allow_download: bool) -> Optional[pathlib.Path]:
    target_dir.mkdir(parents=True, exist_ok=True)
    path = target_dir / name
    if path.exists() and not needs_refresh(path):
        return path
    if not allow_download and not path.exists():
        return None
    url = SUMMARY_URLS[name]
    tmp_name: Optional[str] = None
    try:
        with tempfile.NamedTemporaryFile("wb", delete=False) as tmp:
            tmp_name = tmp.name
            urllib.request.urlretrieve(url, tmp_name)
            tmp.flush()
            os.replace(tmp_name, path)
    except Exception:
        # On failure keep existing file (if any) and proceed.
        try:
            if tmp_name:
                os.unlink(tmp_name)  # best effort cleanup
        except Exception:
            pass
        if not path.exists():
            return None
    return path if path.exists() else None


def load_species_map(directory: pathlib.Path, allow_download: bool) -> Dict[str, Tuple[str, str]]:
    mapping: Dict[str, Tuple[str, str]] = {}
    for name in SUMMARY_URLS:
        path = ensure_assembly_summary(name, directory, allow_download)
        if not path:
            continue
        try:
            with path.open("r", encoding="utf-8", errors="ignore") as handle:
                reader = csv.reader(handle, delimiter="\t")
                for row in reader:
                    if not row or row[0].startswith("#"):
                        continue
                    if len(row) < 8:
                        continue
                    accession = row[0].strip()
                    # Prefer species_taxid (column index 6), fall back to taxid (index 5).
                    species_taxid = (row[6] or row[5]).strip() if len(row) > 6 else row[5].strip()
                    organism_name = row[7].strip() if len(row) > 7 else ""
                    if accession:
                        mapping[accession] = (species_taxid or accession, organism_name or accession)
        except OSError:
            continue
    return mapping


def accession_from_filename(candidate: str) -> str:
    pieces = candidate.split("_", 2)
    if len(pieces) >= 2:
        return f"{pieces[0]}_{pieces[1]}"
    return candidate


def build_candidate_objects(
    names: Sequence[str],
    scores: Dict[str, float],
    species_map: Dict[str, Tuple[str, str]],
    dedupe: bool,
) -> List[Candidate]:
    candidates: List[Candidate] = []
    for idx, name in enumerate(names):
        score = scores.get(name, float("-inf"))
        accession = accession_from_filename(name)
        species_key, species_label = species_map.get(accession, (accession, accession))
        if not dedupe:
            species_key = name  # treat each candidate as unique
        candidates.append(
            Candidate(
                name=name,
                score=score,
                species_key=species_key,
                species_label=species_label,
                order=idx,
            )
        )
    # Sort by score descending, tie-break on original order to preserve determinism.
    candidates.sort(key=lambda c: (-c.score, c.order))
    return candidates


def choose_candidates(candidates: Sequence[Candidate], limit: int) -> Tuple[List[Candidate], Dict[str, int]]:
    chosen: List[Candidate] = []
    seen_species: Dict[str, int] = {}
    for cand in candidates:
        if cand.species_key in seen_species:
            continue
        seen_species[cand.species_key] = len(chosen)
        chosen.append(cand)
        if limit > 0 and len(chosen) >= limit:
            break
    return chosen, seen_species


def write_candidates(path: pathlib.Path, candidates: Sequence[Candidate]) -> None:
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        for cand in candidates:
            handle.write(f"{cand.name}\n")
    os.replace(tmp_path, path)


def append_log(log_path: pathlib.Path, message: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip("\n") + "\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    selected_path = pathlib.Path(args.selected)
    output_path = pathlib.Path(args.output)

    if args.max <= 0:
        raise SystemExit("The --max value must be greater than zero.")

    candidates_raw = read_candidates(selected_path)
    if not candidates_raw:
        raise SystemExit(f"No candidates found in {selected_path}")

    score_files = [pathlib.Path(p) for p in args.score_files]
    score_map = load_scores(score_files)

    assembly_dir = (
        pathlib.Path(args.assembly_dir)
        if args.assembly_dir
        else (pathlib.Path(__file__).resolve().parent.parent / "data" / "downloaded_genomes" / "assembly_summaries")
    )
    species_map = load_species_map(assembly_dir, allow_download=not args.no_download) if args.dedupe else {}

    candidate_objs = build_candidate_objects(candidates_raw, score_map, species_map, args.dedupe)
    chosen, seen_species = choose_candidates(candidate_objs, args.max)

    write_candidates(output_path, chosen)

    kept = len(chosen)
    deduped = len(seen_species)
    original = len(candidates_raw)
    summary = (
        f"[limit_candidates] kept {kept} / {original} candidates "
        f"({deduped if args.dedupe else kept} unique keys) "
        f"{'(species dedupe)' if args.dedupe else ''}"
    )
    print(summary)

    if args.log:
        append_log(pathlib.Path(args.log), summary)

    return 0


if __name__ == "__main__":
    sys.exit(main())
