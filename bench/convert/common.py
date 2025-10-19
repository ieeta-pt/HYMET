#!/usr/bin/env python3
"""Helpers for converting third-party tool outputs into CAMI profile format."""

from __future__ import annotations

import csv
import os
import pathlib
import subprocess
from typing import Iterable, List, Dict, Tuple

RANKS: List[str] = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_CODES: Dict[str, str] = {
    "U": None,
    "R": None,
    "D": "superkingdom",
    "K": "superkingdom",
    "P": "phylum",
    "C": "class",
    "O": "order",
    "F": "family",
    "G": "genus",
    "S": "species",
}


def ensure_parent(path: str) -> None:
    pathlib.Path(path).expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)


def default_taxpath() -> Tuple[List[str], List[str]]:
    ids = ["NA"] * len(RANKS)
    names = ["NA"] * len(RANKS)
    return ids, names


def normalise_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    total = sum(float(row.get("percentage", 0.0)) for row in rows)
    if total <= 0:
        return rows
    for row in rows:
        row["percentage"] = 100.0 * float(row.get("percentage", 0.0)) / total
    return rows


def _format_path(value) -> str:
    if value is None:
        return "|".join(["NA"] * len(RANKS))
    if isinstance(value, str):
        if value:
            return value
        return "|".join(["NA"] * len(RANKS))
    return "|".join(str(v) for v in value)


def write_cami_profile(
    rows: Iterable[Dict[str, str]],
    out_path: str,
    sample_id: str,
    tool_name: str,
    normalise: bool = False,
) -> None:
    rows = list(rows)
    if normalise:
        rows = normalise_rows(rows)
    ensure_parent(out_path)
    with open(out_path, "w", newline="") as handle:
        handle.write(f"@SampleID:\t{sample_id}\n")
        handle.write("@Version:\t0.9.1\n")
        handle.write("@Ranks:\t" + "|".join(RANKS) + "\n")
        handle.write(f"@ToolID:\t{tool_name}\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["@@TAXID", "RANK", "TAXPATH", "TAXPATHSN", "PERCENTAGE"])
        for row in rows:
            taxid = str(row.get("taxid", "NA")).strip() or "NA"
            rank = str(row.get("rank", "unknown")).lower()
            taxpath = _format_path(row.get("taxpath"))
            taxpathsn = _format_path(row.get("taxpathsn"))
            perc = float(row.get("percentage", 0.0))
            writer.writerow([
                taxid,
                rank,
                taxpath,
                taxpathsn,
                f"{perc:.6f}",
            ])


def _run_taxonkit(args, stdin, taxdb: str) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    if taxdb:
        env["TAXONKIT_DB"] = taxdb
    try:
        proc = subprocess.run(
            args,
            input=stdin,
            text=True,
            capture_output=True,
            check=True,
            env=env,
        )
    except FileNotFoundError as exc:
        raise RuntimeError("taxonkit executable not found. Please install taxonkit and ensure TAXONKIT_DB points to the NCBI taxonomy directory.") from exc
    return proc


def _normalise_path_str(path: str) -> str:
    parts = (path or "").split("|") if path else []
    if len(parts) < len(RANKS):
        parts.extend(["NA"] * (len(RANKS) - len(parts)))
    return "|".join(parts[:len(RANKS)])


def taxonkit_taxpath(taxids: Iterable[str], taxdb: str) -> Dict[str, Tuple[str, str]]:
    tids = [t for t in dict.fromkeys(taxids) if t and t != "NA"]
    if not tids:
        return {}
    proc = _run_taxonkit(
        [
            "taxonkit",
            "reformat",
            "-I",
            "1",
            "-f",
            "{d}|{p}|{c}|{o}|{f}|{g}|{s}",
            "-t",
            "-T",
        ]
        + (["--data-dir", taxdb] if taxdb else []),
        "\n".join(tids) + "\n",
        taxdb,
    )
    mapping: Dict[str, Tuple[str, str]] = {}
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            names = _normalise_path_str(parts[1])
            ids = _normalise_path_str(parts[2])
            mapping[parts[0]] = (ids, names)
    return mapping


def taxonkit_name2taxid(names: Iterable[str], taxdb: str) -> Dict[str, Tuple[str, str]]:
    unique = [n for n in dict.fromkeys(names) if n and n.lower() != "unclassified"]
    if not unique:
        return {}
    proc = _run_taxonkit(
        ["taxonkit", "name2taxid", "--show-rank"] + (["--data-dir", taxdb] if taxdb else []),
        "\n".join(unique) + "\n",
        taxdb,
    )
    mapping: Dict[str, Tuple[str, str]] = {}
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3 and parts[1].isdigit():
            mapping[parts[0]] = (parts[1], parts[2])
    return mapping
