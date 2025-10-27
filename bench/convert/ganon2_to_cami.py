#!/usr/bin/env python3
"""Convert ganon2 kreport output into CAMI format."""

from __future__ import annotations

import argparse
import os
import sys

if __package__ is None or __package__ == "":  # pragma: no cover - CLI fallback
    sys.path.append(os.path.dirname(__file__))
    from kreport import parse_kreport  # type: ignore
    from common import write_cami_profile  # type: ignore
else:  # pragma: no cover
    from .kreport import parse_kreport
    from .common import write_cami_profile


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert ganon2 reports to CAMI format.")
    ap.add_argument("--report", required=True, help="Path to ganon report (kreport format).")
    ap.add_argument("--out", required=True, help="Output CAMI TSV path.")
    ap.add_argument("--sample-id", required=True, help="Sample identifier.")
    ap.add_argument("--tool", default="ganon2", help="Tool identifier.")
    args = ap.parse_args()

    rows = parse_kreport(args.report)
    write_cami_profile(rows, args.out, sample_id=args.sample_id, tool_name=args.tool, normalise=True)


if __name__ == "__main__":
    main()
