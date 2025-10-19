#!/usr/bin/env python3
"""Convert HYMET classified_sequences.tsv into a CAMI taxonomic profile."""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_ALIAS = {
    "domain": "superkingdom",
    "kingdom": "superkingdom",
    "sk": "superkingdom",
    "k": "superkingdom",
    "phylum": "phylum",
    "p": "phylum",
    "class": "class",
    "c": "class",
    "order": "order",
    "o": "order",
    "family": "family",
    "f": "family",
    "genus": "genus",
    "g": "genus",
    "species": "species",
    "s": "species",
}


def _run(cmd: List[str], input_text: str | None = None) -> str:
    proc = subprocess.run(
        cmd,
        input=input_text,
        text=True,
        capture_output=True,
        check=True,
    )
    return proc.stdout


def parse_lineage(lineage: str) -> Dict[str, str]:
    out = {rank: "" for rank in RANKS}
    if not lineage:
        return out
    for part in lineage.split(";"):
        part = part.strip()
        if not part or ":" not in part:
            continue
        rk, name = part.split(":", 1)
        rk = RANK_ALIAS.get(rk.strip().lower(), rk.strip().lower())
        if rk in out:
            out[rk] = name.strip()
    return out

def load_records(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    taxid_records: List[str] = []
    lineage_records: List[Dict[str, str]] = []
    with path.open(encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        has_taxid = 'TaxID' in reader.fieldnames if reader.fieldnames else False
        for row in reader:
            taxid = (row.get("TaxID") or "").strip() if has_taxid else ""
            if taxid and taxid.lower() != "unknown":
                taxid_records.append(taxid)
                continue
            parsed = parse_lineage(row.get("Lineage", ""))
            if any(parsed.values()):
                lineage_records.append(parsed)
    return taxid_records, lineage_records


def batch_name2taxid(names: Iterable[str], taxdb: str) -> Dict[str, str]:
    names = sorted({n for n in names if n})
    if not names:
        return {}
    cmd = ["taxonkit", "name2taxid", "--data-dir", taxdb, "--show-rank"]
    input_text = "\n".join(names) + "\n"
    output = _run(cmd, input_text)
    mapping: Dict[str, str] = {}
    for line in output.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            mapping[parts[0]] = parts[1]
    return mapping


def batch_taxpath(taxids: Iterable[str], taxdb: str) -> Dict[str, Dict[str, object]]:
    taxids = sorted({t for t in taxids if t})
    if not taxids:
        return {}
    cmd = [
        "taxonkit",
        "reformat",
        "--data-dir",
        taxdb,
        "-I",
        "1",
        "-f",
        "{d}|{p}|{c}|{o}|{f}|{g}|{s}",
        "-t",
        "-T",
    ]
    input_text = "\n".join(taxids) + "\n"
    output = _run(cmd, input_text)
    mapping: Dict[str, Dict[str, object]] = {}
    for line in output.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            names = parts[1]
            ids = parts[2]
            def _pad(values: List[str]) -> List[str]:
                padded = [v if v and v.upper() != "NA" else "" for v in values]
                if len(padded) < len(RANKS):
                    padded.extend([""] * (len(RANKS) - len(padded)))
                return padded[: len(RANKS)]
            names_list = _pad(names.split("|"))
            ids_list = _pad(ids.split("|"))
            mapping[parts[0]] = {
                "names": names,
                "ids": ids,
                "names_list": names_list,
                "ids_list": ids_list,
            }
    return mapping


def accumulate(
    taxid_records: List[str],
    lineage_records: List[Dict[str, str]],
    taxdb: str,
) -> Tuple[
    Dict[str, Dict[str, float]],
    Dict[str, float],
    Dict[str, Dict[str, dict]],
    Dict[str, str],
]:
    counts = {rank: defaultdict(float) for rank in RANKS}
    totals = {rank: 0.0 for rank in RANKS}
    meta = {rank: {} for rank in RANKS}

    all_names = set()
    for parsed in lineage_records:
        for name in parsed.values():
            if name:
                all_names.add(name)
    name2tid = batch_name2taxid(all_names, taxdb)

    taxids_needed: set[str] = set()
    pending: List[Tuple[str, int, float]] = []

    for parsed in lineage_records:
        for idx, rank in enumerate(RANKS):
            name = parsed.get(rank)
            if not name:
                continue
            tid = name2tid.get(name)
            if not tid:
                continue
            taxids_needed.add(tid)
            pending.append((tid, idx, 1.0))

    for tid in taxid_records:
        if not tid:
            continue
        taxids_needed.add(tid)
        for idx in range(len(RANKS)):
            pending.append((tid, idx, 1.0))

    taxid2path = batch_taxpath(taxids_needed, taxdb)

    def add_count(tid: str, idx: int, weight: float) -> None:
        info = taxid2path.get(tid)
        if not info:
            return
        ids_list = info.get("ids_list", [])
        names_list = info.get("names_list", [])
        ancestor_tid = tid
        if idx < len(ids_list) and ids_list[idx]:
            ancestor_tid = ids_list[idx]
        rank = RANKS[idx]
        counts[rank][ancestor_tid] += weight
        totals[rank] += weight
        if ancestor_tid not in meta[rank]:
            taxpath_ids = "|".join(ids_list[: idx + 1]) if ids_list else ancestor_tid
            taxpath_names = "|".join(names_list[: idx + 1]) if names_list else ""
            meta[rank][ancestor_tid] = {
                "taxpath_ids": taxpath_ids or ancestor_tid,
                "taxpath_names": taxpath_names or (names_list[idx] if idx < len(names_list) else ""),
            }

    for tid, idx, weight in pending:
        if 0 <= idx < len(RANKS):
            add_count(tid, idx, weight)

    return counts, totals, meta, name2tid

def emit_cami(counts: Dict[str, Dict[str, float]], totals: Dict[str, float], meta: Dict[str, Dict[str, dict]]) -> None:
    print("#CAMI Submission for Taxonomic Profiling")
    print("@Version:0.9.1 @Ranks:superkingdom|phylum|class|order|family|genus|species @SampleID:sample_0")
    print("@@TAXID RANK TAXPATH TAXPATHSN PERCENTAGE")
    for rank in RANKS:
        total = totals.get(rank, 0.0)
        if total <= 0:
            continue
        rank_counts = counts.get(rank, {})
        rank_meta = meta.get(rank, {})
        for tid, count in sorted(rank_counts.items(), key=lambda kv: kv[1], reverse=True):
            info = rank_meta.get(tid)
            if not info:
                continue
            ids = info.get("taxpath_ids", "")
            names = info.get("taxpath_names", "")
            pct = 100.0 * count / total
            print(f"{tid}	{rank}	{ids}	{names}	{pct:.6f}")




def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: hymet2cami.py <classified_sequences.tsv>", file=sys.stderr)
        sys.exit(1)

    input_path = Path(sys.argv[1]).resolve()
    if not input_path.is_file():
        print(f"Missing classified_sequences TSV: {input_path}", file=sys.stderr)
        sys.exit(1)

    taxdb = os.environ.get("TAXONKIT_DB", str(Path(__file__).resolve().parents[1] / "taxonomy_files"))
    print(f"[hymet2cami] using taxonomy DB at {taxdb}", file=sys.stderr)

    taxid_records, lineage_records = load_records(input_path)
    print(
        f"[hymet2cami] parsed {len(lineage_records)} lineage rows + {len(taxid_records)} direct taxid rows",
        file=sys.stderr,
    )

    all_names = set()
    for parsed in lineage_records:
        for name in parsed.values():
            if name:
                all_names.add(name)

    print(f"[hymet2cami] converting {len(all_names)} unique taxon names", file=sys.stderr)

    counts, totals, meta, name2tid = accumulate(taxid_records, lineage_records, taxdb)
    print(f"[hymet2cami] mapped {len(name2tid)} names to taxids", file=sys.stderr)
    total_taxids = sum(len(v) for v in meta.values())
    print(f"[hymet2cami] converting {total_taxids} taxids to paths", file=sys.stderr)

    emit_cami(counts, totals, meta)
    print("[hymet2cami] done", file=sys.stderr)



if __name__ == "__main__":
    main()
