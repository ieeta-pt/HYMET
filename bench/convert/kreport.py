#!/usr/bin/env python3
"""Utilities to convert Kraken-style kreport tables to CAMI format."""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, Iterable, List, Tuple

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
    sys.path.append(os.path.dirname(__file__))
    from common import RANKS, RANK_CODES, default_taxpath, write_cami_profile  # type: ignore
else:  # pragma: no cover
    from .common import RANKS, RANK_CODES, default_taxpath, write_cami_profile


def parse_kreport(report_path: str) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    rank_state: Dict[str, Dict[str, str]] = {
        rank: {"taxid": "NA", "name": "NA"} for rank in RANKS
    }
    stack: List[Tuple[str, str, str]] = []

    with open(report_path, "r") as handle:
        for raw_line in handle:
            if not raw_line.strip():
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                perc = float(parts[0])
            except ValueError:
                continue
            rank_code = parts[3].strip()
            if not rank_code:
                continue
            rank_code = rank_code[0].upper()
            taxid = parts[4].strip()
            name_field = parts[5]
            leading = len(name_field) - len(name_field.lstrip())
            depth = max(0, leading // 2)
            name = name_field.strip() or taxid or "unknown"

            while len(stack) > depth:
                code, _, _ = stack.pop()
                mapped = RANK_CODES.get(code)
                if mapped in rank_state:
                    rank_state[mapped] = {"taxid": "NA", "name": "NA"}

            stack.append((rank_code, taxid, name))
            mapped_rank = RANK_CODES.get(rank_code)

            if mapped_rank:
                rank_state[mapped_rank] = {"taxid": taxid, "name": name}
                for lower in RANKS[RANKS.index(mapped_rank) + 1 :]:
                    rank_state[lower] = {"taxid": "NA", "name": "NA"}

            if mapped_rank and perc > 0.0 and taxid not in {"0", "", "NA"}:
                taxpath = [rank_state[r]["taxid"] for r in RANKS]
                names = [rank_state[r]["name"] for r in RANKS]
                rows.append(
                    {
                        "taxid": taxid,
                        "rank": mapped_rank,
                        "taxpath": taxpath,
                        "taxpathsn": names,
                        "percentage": perc,
                    }
                )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert Kraken-style reports to CAMI format.")
    ap.add_argument("--report", required=True, help="Input kreport file.")
    ap.add_argument("--out", required=True, help="Output CAMI TSV path.")
    ap.add_argument("--sample-id", required=True, help="Sample identifier.")
    ap.add_argument("--tool", default="unknown", help="Tool identifier.")
    args = ap.parse_args()

    rows = parse_kreport(args.report)
    write_cami_profile(rows, args.out, sample_id=args.sample_id, tool_name=args.tool, normalise=True)


if __name__ == "__main__":
    main()
