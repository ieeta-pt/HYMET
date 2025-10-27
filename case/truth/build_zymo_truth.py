#!/usr/bin/env python3
"""
Derive Zymo mock community truth labels by mapping contigs against a multi-strain
reference panel and assigning per-contig TaxIDs using a simple best-hit/LCA rule.

Outputs:
  - Per-contig truth table (contig_id, taxid, rank, match_bases, identity_pct)
  - CAMI-style truth profile (length-weighted abundance per rank)

Usage:
  build_zymo_truth.py --contigs /path/to/contigs.fna \
                      --seqmap HYMET/case/truth/zymo_refs/seqid2taxid.tsv \
                      --paf HYMET/case/truth/zymo_mc/zymo_mc_vs_refs.paf \
                      --out-contigs HYMET/case/truth/zymo_mc/truth_contigs.tsv \
                      --out-profile HYMET/case/truth/zymo_mc/truth_profile.cami.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


def load_seqmap(path: Path) -> Dict[str, int]:
    mapping: Dict[str, int] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            seqid, taxid = line.strip().split("\t")[:2]
            mapping[seqid] = int(taxid)
    return mapping


def load_nodes(path: Path) -> Tuple[Dict[int, int], Dict[int, str]]:
    parent: Dict[int, int] = {}
    rank: Dict[int, str] = {}
    with path.open() as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3:
                continue
            tid = int(parts[0])
            parent[tid] = int(parts[1])
            rank[tid] = parts[2]
    return parent, rank


def climb_to_rank(taxid: int, target_rank: str, parent: Dict[int, int], rank: Dict[int, str]) -> Optional[int]:
    seen = set()
    current = taxid
    while current not in seen:
        seen.add(current)
        if rank.get(current) == target_rank:
            return current
        nxt = parent.get(current)
        if nxt is None or nxt == current:
            return None
        current = nxt
    return None


def load_contig_lengths(path: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with path.open() as fh:
        current = None
        acc = 0
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = acc
                current = line[1:].split()[0]
                acc = 0
            else:
                acc += len(line)
        if current is not None:
            lengths[current] = acc
    return lengths


def parse_paf(
    path: Path,
    seq2tax: Dict[str, int],
    min_match: int,
    min_identity: float,
    min_coverage: float,
) -> Dict[str, List[Tuple[int, str, int, float, float]]]:
    hits: Dict[str, List[Tuple[int, str, int, float, float]]] = defaultdict(list)
    with path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            query = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2])
            qend = int(parts[3])
            target = parts[5]
            match = int(parts[9])
            block = int(parts[10])
            mapq = int(parts[11])
            if match < min_match or block <= 0:
                continue
            dv = None
            for tag in parts[12:]:
                if tag.startswith("dv:f:"):
                    dv = float(tag.split(":", 2)[2])
                    break
            identity = 1.0 - dv if dv is not None else match / block
            if identity < min_identity:
                continue
            cov = (qend - qstart) / qlen if qlen > 0 else 0.0
            if cov < min_coverage:
                continue
            taxid = seq2tax.get(target)
            if taxid is None:
                continue
            hits[query].append((taxid, target, match, identity, cov))
    return hits


def assign_taxids(
    hits: Dict[str, List[Tuple[int, str, int, float, float]]],
    parent: Dict[int, int],
    rank: Dict[int, str],
    tolerance: float,
) -> Dict[str, Tuple[int, str, int, float, float]]:
    assignments: Dict[str, Tuple[int, str, int, float, float]] = {}
    for contig, rows in hits.items():
        if not rows:
            continue
        best_match = max(row[2] for row in rows)
        threshold = best_match * (1.0 - tolerance)
        kept = [row for row in rows if row[2] >= threshold]
        species_taxids = {row[0] for row in kept}
        chosen_rank = "species"
        chosen_taxid: Optional[int] = None
        if len(species_taxids) == 1:
            chosen_taxid = next(iter(species_taxids))
            chosen_rank = "species"
        else:
            genus_taxids = set()
            for taxid in species_taxids:
                genus = climb_to_rank(taxid, "genus", parent, rank)
                if genus:
                    genus_taxids.add(genus)
            if len(genus_taxids) == 1:
                chosen_taxid = next(iter(genus_taxids))
                chosen_rank = "genus"
        if chosen_taxid:
            # Representative match/identity from best hit
            primary = max(kept, key=lambda r: r[2])
            assignments[contig] = (chosen_taxid, chosen_rank, primary[2], primary[3] * 100.0, primary[4] * 100.0)
    return assignments


def taxonkit_paths(taxids: Iterable[int], taxdb: Path) -> Dict[int, Tuple[str, str]]:
    taxids = sorted({tid for tid in taxids if tid})
    if not taxids:
        return {}
    cmd = [
        "taxonkit",
        "reformat",
        "--data-dir",
        str(taxdb),
        "-I",
        "1",
        "-f",
        "{d}|{p}|{c}|{o}|{f}|{g}|{s}",
        "-t",
    ]
    proc = subprocess.run(cmd, input="\n".join(map(str, taxids)) + "\n", text=True, capture_output=True, check=True)
    mapping: Dict[int, Tuple[str, str]] = {}
    for line in proc.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            tid = int(parts[0])
            mapping[tid] = (parts[1], parts[2])
    return mapping


def build_profile(
    assignments: Dict[str, Tuple[int, str, int, float, float]],
    lengths: Dict[str, int],
    taxdb: Path,
) -> List[Tuple[str, str, str, str, float]]:
    totals = Counter()  # rank -> total length
    accum: Dict[str, Counter] = {rank: Counter() for rank in RANKS}
    needed_taxids: set[int] = set()
    for contig, (taxid, _, _, _, _) in assignments.items():
        length = lengths.get(contig, 1)
        if length <= 0:
            continue
        needed_taxids.add(taxid)
    tax_paths = taxonkit_paths(needed_taxids, taxdb)

    # Expand to include ancestor taxids from paths
    ancestor_taxids: set[int] = set()
    for names_ids in tax_paths.values():
        _, ids = names_ids
        for tid in ids.split("|"):
            if tid and tid != "NA":
                ancestor_taxids.add(int(tid))
    tax_paths.update(taxonkit_paths(ancestor_taxids - set(tax_paths), taxdb))

    for contig, (taxid, _, _, _, _) in assignments.items():
        length = lengths.get(contig, 1)
        names_ids = tax_paths.get(taxid)
        if not names_ids:
            continue
        _, ids = names_ids
        id_list = ids.split("|")
        for idx, rank in enumerate(RANKS):
            if idx >= len(id_list):
                continue
            tid = id_list[idx]
            if not tid or tid == "NA":
                continue
            tid_int = int(tid)
            accum[rank][tid_int] += length
            totals[rank] += length

    profile_rows: List[Tuple[str, str, str, str, float]] = []
    for rank in RANKS:
        total = totals.get(rank, 0)
        if total <= 0:
            continue
        for tid, length in accum[rank].items():
            names_ids = tax_paths.get(tid)
            if not names_ids:
                continue
            names, ids = names_ids
            pct = 100.0 * length / total if total else 0.0
            profile_rows.append((str(tid), rank, ids, names, pct))
    return profile_rows


def format_lineage(names_ids: Tuple[str, str]) -> str:
    names, _ = names_ids
    parts = names.split("|")
    lineage_parts = []
    for rank, name in zip(RANKS, parts):
        if name and name != "NA":
            lineage_parts.append(f"{rank}:{name}")
    return ";".join(lineage_parts)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build Zymo mock truth data from contig mappings.")
    ap.add_argument("--contigs", required=True, help="FASTA with assembly contigs (query).")
    ap.add_argument("--seqmap", required=True, help="TSV mapping reference seqid to species TaxID.")
    ap.add_argument("--paf", required=True, help="PAF from mapping contigs to reference panel.")
    ap.add_argument("--out-contigs", required=True, help="Output TSV with contig truth assignments.")
    ap.add_argument("--out-profile", required=True, help="Output CAMI-style profile TSV.")
    ap.add_argument("--taxonomy-dir", default=str(Path(__file__).resolve().parents[2] / "taxonomy_files"), help="NCBI taxonomy directory for taxonkit.")
    ap.add_argument("--min-match", type=int, default=1000, help="Minimum matching bases per alignment.")
    ap.add_argument("--min-identity", type=float, default=0.95, help="Minimum identity per alignment (0-1).")
    ap.add_argument("--tolerance", type=float, default=0.01, help="Fractional tolerance for considering alternative hits.")
    ap.add_argument("--min-coverage", type=float, default=0.5, help="Minimum fraction of contig covered by an alignment.")
    args = ap.parse_args()

    seqmap = load_seqmap(Path(args.seqmap))
    parent, rank = load_nodes(Path(args.taxonomy_dir) / "nodes.dmp")
    lengths = load_contig_lengths(Path(args.contigs))
    hits = parse_paf(Path(args.paf), seqmap, args.min_match, args.min_identity, args.min_coverage)
    assignments = assign_taxids(hits, parent, rank, args.tolerance)

    assigned_species = sum(1 for _, (_, r, _, _, _) in assignments.items() if r == "species")
    assigned_genus = sum(1 for _, (_, r, _, _, _) in assignments.items() if r == "genus")
    print(f"[truth] Assigned {len(assignments)} contigs ({assigned_species} species-level, {assigned_genus} genus-level)")

    taxdb = Path(args.taxonomy_dir)
    tax_paths = taxonkit_paths({tid for tid, *_ in assignments.values()}, taxdb)

    # Write per-contig truth table
    with open(args.out_contigs, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["contig_id", "taxid", "rank", "match_bases", "identity_percent", "coverage_percent"])
        for contig, (taxid, assigned_rank, match, ident, cov) in sorted(assignments.items()):
            writer.writerow([contig, taxid, assigned_rank, match, f"{ident:.2f}", f"{cov:.2f}"])

    # Build profile (length-weighted)
    profile_rows = build_profile(assignments, lengths, taxdb)
    with open(args.out_profile, "w") as out:
        out.write("#CAMI Submission for Taxonomic Profiling\n")
        out.write("@Version:0.9.1 @Ranks:superkingdom|phylum|class|order|family|genus|species @SampleID:zymo_mc_truth\n")
        out.write("@@TAXID RANK TAXPATH TAXPATHSN PERCENTAGE\n")
        ALT_SUPERKINGDOM = {"Bacteria": 3379134}
        for tid, rank_name, ids, names, pct in sorted(profile_rows, key=lambda x: (RANKS.index(x[1]), -x[4])):
            write_tid = str(tid)
            if rank_name == "superkingdom":
                first_name = names.split("|")[0] if names else ""
                alt_tid = ALT_SUPERKINGDOM.get(first_name)
                if alt_tid:
                    write_tid = str(alt_tid)
            out.write(f"{write_tid}\t{rank_name}\t{ids}\t{names}\t{pct:.6f}\n")


if __name__ == "__main__":
    main()
