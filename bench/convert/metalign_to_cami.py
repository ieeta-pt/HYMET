#!/usr/bin/env python3
"""Normalise Metalign CAMI-like TSV to canonical CAMI profile.

Metalign typically emits a CAMI-format profile already. This helper:
- Accepts Metalign's output TSV
- If it already contains CAMI metadata/header, it rewrites SampleID/ToolID
- Otherwise, tries to interpret a minimal table with the expected columns

Outputs a canonical CAMI TSV compatible with the HYMET evaluation helpers.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict, Iterable, List

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
    sys.path.append(os.path.dirname(__file__))
    from common import RANKS, write_cami_profile  # type: ignore
else:  # pragma: no cover
    from .common import RANKS, write_cami_profile


def looks_like_cami(lines: List[str]) -> bool:
    if not lines:
        return False
    meta = [l for l in lines if l.startswith("@")]
    header = [l for l in lines if l.startswith("@@")]
    return any("@SampleID:" in l for l in meta) and any("@@TAXID" in l for l in header)


def rewrite_cami(in_lines: List[str], sample_id: str, tool: str) -> List[str]:
    out: List[str] = []
    meta_done = set()
    for line in in_lines:
        if line.startswith("@SampleID:"):
            out.append(f"@SampleID:\t{sample_id}\n")
            meta_done.add("SampleID")
        elif line.startswith("@ToolID:"):
            out.append(f"@ToolID:\t{tool}\n")
            meta_done.add("ToolID")
        elif line.startswith("@Ranks:"):
            out.append("@Ranks:\t" + "|".join(RANKS) + "\n")
        else:
            out.append(line)
    # Ensure required meta lines exist
    existing = "".join(in_lines)
    if "@SampleID:" not in existing:
        out.insert(0, f"@SampleID:\t{sample_id}\n")
    if "@Version:" not in existing:
        out.insert(1, "@Version:\t0.9.1\n")
    if "@Ranks:" not in existing:
        out.insert(2, "@Ranks:\t" + "|".join(RANKS) + "\n")
    if "@ToolID:" not in existing:
        out.insert(3, f"@ToolID:\t{tool}\n")
    return out


def maybe_parse_minimal_table(path: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    # Accept lowercase/uppercase headers and minor variations
    field_map = {
        "taxid": {"taxid", "@@taxid", "id"},
        "rank": {"rank", "@@rank"},
        "taxpath": {"taxpath", "@@taxpath"},
        "taxpathsn": {"taxpathsn", "@@taxpathsn"},
        "percentage": {"percentage", "@@percentage", "abundance", "%"},
    }
    with open(path, "r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return rows
        header_norm = [h.strip().lstrip("@").lower() for h in header]
        idx: Dict[str, int] = {}
        for key, aliases in field_map.items():
            for i, name in enumerate(header_norm):
                if name in aliases:
                    idx.setdefault(key, i)
        if not {"taxid", "rank", "taxpath", "taxpathsn", "percentage"}.issubset(idx.keys()):
            return []
        for parts in reader:
            if not parts or all(not p.strip() for p in parts):
                continue
            try:
                perc = float(parts[idx["percentage"]])
            except Exception:
                continue
            rows.append(
                {
                    "taxid": parts[idx["taxid"]].strip() or "NA",
                    "rank": parts[idx["rank"]].strip().lower() or "species",
                    "taxpath": parts[idx["taxpath"]].strip(),
                    "taxpathsn": parts[idx["taxpathsn"]].strip(),
                    "percentage": perc,
                }
            )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="Normalise Metalign TSV to CAMI profile.")
    ap.add_argument("--input", required=True, help="Metalign output TSV (abundance/profile)")
    ap.add_argument("--out", required=True, help="Output CAMI TSV path")
    ap.add_argument("--sample-id", required=True, help="Sample identifier")
    ap.add_argument("--tool", default="metalign", help="Tool name for CAMI header")
    args = ap.parse_args()

    with open(args.input, "r", encoding="utf-8", errors="ignore") as fin:
        lines = fin.readlines()

    if looks_like_cami(lines):
        out_lines = rewrite_cami(lines, args.sample_id, args.tool)
        with open(args.out, "w", encoding="utf-8") as fout:
            fout.writelines(out_lines)
        return

    # Fall back to minimal table reader
    rows = maybe_parse_minimal_table(args.input)
    write_cami_profile(rows, args.out, sample_id=args.sample_id, tool_name=args.tool, normalise=True)


if __name__ == "__main__":
    main()

