#!/usr/bin/env python3
"""
prune_cache.py

Utility to prune HYMET's downloaded genome cache based on age and/or total
size. Designed for the benchmark harness which defaults to
data/downloaded_genomes/cache_bench, but works with any cache directory that
contains hashed run subdirectories.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import shutil
import sys
import time
from dataclasses import dataclass
from typing import Iterable, List


@dataclass
class CacheEntry:
    path: pathlib.Path
    size_bytes: int
    mtime: float

    @property
    def age_days(self) -> float:
        return max(0.0, (time.time() - self.mtime) / 86400.0)

    def human_size(self) -> str:
        units = ["B", "KiB", "MiB", "GiB", "TiB"]
        size = float(self.size_bytes)
        for unit in units:
            if size < 1024.0:
                return f"{size:.1f} {unit}"
            size /= 1024.0
        return f"{size:.1f} PiB"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prune HYMET downloaded genome cache directories by age or total size."
    )
    parser.add_argument(
        "cache_root",
        nargs="?",
        default=None,
        help="Cache root to prune (default: HYMET/data/downloaded_genomes/cache_bench relative to this script).",
    )
    parser.add_argument(
        "--max-age-days",
        type=float,
        default=None,
        help="Remove cache directories whose newest file is older than this many days.",
    )
    parser.add_argument(
        "--max-size-gb",
        type=float,
        default=None,
        help="Ensure cache size stays below this threshold (GiB). Oldest entries removed first.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show which directories would be removed without deleting them.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="List cache entries before pruning.",
    )
    return parser.parse_args()


def compute_directory_size(path: pathlib.Path) -> int:
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            try:
                total += (pathlib.Path(root) / name).stat().st_size
            except OSError:
                continue
    return total


def scan_cache(root: pathlib.Path) -> List[CacheEntry]:
    entries: List[CacheEntry] = []
    if not root.exists():
        return entries
    for entry in root.iterdir():
        if not entry.is_dir():
            continue
        try:
            stat = entry.stat()
        except OSError:
            continue
        size_bytes = compute_directory_size(entry)
        entries.append(CacheEntry(path=entry, size_bytes=size_bytes, mtime=stat.st_mtime))
    return entries


def remove_entry(entry: CacheEntry, dry_run: bool) -> None:
    if dry_run:
        print(f"[dry-run] would remove {entry.path.name} ({entry.human_size()}, {entry.age_days:.1f} days old)")
        return
    print(f"Removing {entry.path.name} ({entry.human_size()}, {entry.age_days:.1f} days old)")
    shutil.rmtree(entry.path, ignore_errors=True)


def prune_by_age(entries: Iterable[CacheEntry], max_age_days: float, dry_run: bool) -> List[CacheEntry]:
    remaining: List[CacheEntry] = []
    for entry in entries:
        if entry.age_days > max_age_days:
            remove_entry(entry, dry_run)
        else:
            remaining.append(entry)
    return remaining


def prune_by_size(entries: List[CacheEntry], max_size_gb: float, dry_run: bool) -> None:
    if max_size_gb is None or max_size_gb <= 0:
        return
    limit_bytes = max_size_gb * (1024 ** 3)
    # Sort by most recent modification (newest last) so we remove oldest first.
    entries.sort(key=lambda e: e.mtime)  # oldest first
    total = sum(e.size_bytes for e in entries)
    if total <= limit_bytes:
        print(f"Cache size {total / (1024 ** 3):.2f} GiB within limit {max_size_gb:.2f} GiB")
        return
    print(f"Cache size {total / (1024 ** 3):.2f} GiB exceeds limit {max_size_gb:.2f} GiB → pruning oldest entries")
    for entry in list(entries):
        if total <= limit_bytes:
            break
        remove_entry(entry, dry_run)
        total -= entry.size_bytes


def main() -> int:
    args = parse_args()
    repo_root = pathlib.Path(__file__).resolve().parents[2]
    default_cache_root = repo_root / "data" / "downloaded_genomes" / "cache_bench"
    cache_root = pathlib.Path(args.cache_root) if args.cache_root else default_cache_root

    if not cache_root.exists():
        print(f"Cache root {cache_root} does not exist; nothing to prune.")
        return 0

    entries = scan_cache(cache_root)
    if not entries:
        print(f"No cache directories found under {cache_root}")
        return 0

    if args.verbose:
        for entry in sorted(entries, key=lambda e: e.mtime):
            print(
                f"{entry.path.name}\t{entry.human_size()}\t{entry.age_days:.1f} days"
            )

    if args.max_age_days is not None and args.max_age_days >= 0:
        entries = prune_by_age(entries, args.max_age_days, args.dry_run)

    prune_by_size(entries, args.max_size_gb, args.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())

