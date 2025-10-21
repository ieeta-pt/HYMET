#!/usr/bin/env python3
"""
Prepare per-group truth data and metadata for the HYMET mutation experiment.

Outputs:
  * Per-group FASTA files (original sequences, 0% mutation).
  * Truth profile + contig mapping TSVs per group (CAMI format).
  * contig_metadata.json with group, taxonomy path, ancestry and length.
"""
from __future__ import annotations

import argparse
import csv
import re
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple
csv.field_size_limit(1024 * 1024 * 1024)

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_ALIAS = {
    "domain": "superkingdom",
    "kingdom": "superkingdom",
    "sk": "superkingdom",
    "k": "superkingdom",
    "superkingdom": "superkingdom",
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
    "subspecies": "species",
}

GROUPS = [
    "Viruses",
    "Archaea",
    "Bacteria",
    "Fungi",
    "Plants",
    "Protozoa",
    "Invertebrates",
    "Vertebrate Mammals",
    "Other Vertebrates",
]


def normalise_rank(rank: str) -> Optional[str]:
    if not rank:
        return None
    return RANK_ALIAS.get(rank.strip().lower())


def load_taxonomy_hierarchy(path: Path) -> Tuple[Dict[str, str], Dict[str, Optional[str]], Dict[str, str]]:
    tax_name: Dict[str, str] = {}
    tax_parent: Dict[str, Optional[str]] = {}
    tax_rank: Dict[str, str] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tid = row["TaxID"].strip()
            tax_name[tid] = row["Name"].strip()
            parent = row["ParentTaxID"].strip()
            tax_parent[tid] = parent if parent and parent != tid else None
            rank = normalise_rank(row["Rank"])
            if rank:
                tax_rank[tid] = rank
    return tax_name, tax_parent, tax_rank


def build_lineage(
    taxid: str,
    tax_name: Dict[str, str],
    tax_parent: Dict[str, Optional[str]],
    tax_rank: Dict[str, str],
) -> Tuple[List[Optional[str]], List[str]]:
    taxids: List[Optional[str]] = [None] * len(RANKS)
    names: List[str] = [""] * len(RANKS)
    current = taxid
    visited = set()
    while current and current not in visited:
        visited.add(current)
        rank = tax_rank.get(current)
        if rank in RANKS:
            idx = RANKS.index(rank)
            if taxids[idx] is None:
                taxids[idx] = current
                names[idx] = tax_name.get(current, "")
        current = tax_parent.get(current)
    return taxids, names


def compute_ancestry(taxid: str, tax_parent: Dict[str, Optional[str]]) -> List[str]:
    ancestry = []
    current = taxid
    visited = set()
    while current and current not in visited:
        visited.add(current)
        ancestry.append(current)
        current = tax_parent.get(current)
    return ancestry


def assign_group(ancestry: List[str]) -> str:
    anc_set = set(ancestry)
    if "10239" in anc_set:
        return "Viruses"
    if "2" in anc_set:
        return "Bacteria"
    if "2157" in anc_set:
        return "Archaea"
    if "2759" in anc_set:
        if "4751" in anc_set:
            return "Fungi"
        if "33090" in anc_set:
            return "Plants"
        if "33208" in anc_set:
            if "7711" in anc_set:
                if "40674" in anc_set:
                    return "Vertebrate Mammals"
                return "Other Vertebrates"
            return "Invertebrates"
        return "Protozoa"
    return "Protozoa"


def load_id_to_taxid(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            taxid = row.get("TaxID", "").strip()
            if not taxid:
                continue
            identifiers = row.get("Identifiers", "")
            if not identifiers:
                continue
            for token in re.split(r"[;,|\s]+", identifiers):
                token = token.strip()
                if token:
                    mapping[token] = taxid
    return mapping


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_truth_contigs(path: Path, entries: List[Tuple[str, str, int]]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "contig_id",
            "taxid",
            "rank",
            "match_bases",
            "identity_percent",
            "coverage_percent",
        ])
        for contig_id, taxid, length in entries:
            writer.writerow([contig_id, taxid, "species", length, "100.0", "100.0"])


def write_truth_profile(path: Path, rows: List[Tuple[str, str, List[str], List[str], float]]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="") as fh:
        fh.write("#CAMI Submission for Taxonomic Profiling\n")
        fh.write("@Version:0.9.1 @Ranks:" + "|".join(RANKS) + " @SampleID:mutation_truth\n")
        fh.write("@@TAXID RANK TAXPATH TAXPATHSN PERCENTAGE\n")
        writer = csv.writer(fh, delimiter="\t")
        for taxid, rank, taxpath_ids, taxpath_names, pct in rows:
            writer.writerow([taxid, rank, "|".join(taxpath_ids), "|".join(taxpath_names), f"{pct:.6f}"])


def stream_fasta(path: Path):
    with path.open() as fh:
        header = None
        seq_lines = []
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line.strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            yield header, "".join(seq_lines)


def contig_id_from_header(header: str) -> str:
    return header[1:].split()[0]


def prepare_dataset(
    fasta_path: Path,
    taxonomy_path: Path,
    hierarchy_path: Path,
    output_dir: Path,
) -> None:
    id_to_taxid = load_id_to_taxid(taxonomy_path)
    tax_name, tax_parent, tax_rank = load_taxonomy_hierarchy(hierarchy_path)

    fasta_handles = {}
    for group in GROUPS:
        path = output_dir / "original" / f"{group.replace(' ', '_').lower()}.fna"
        ensure_dir(path.parent)
        fasta_handles[group] = path.open("w")

    contig_metadata = {}
    group_entries: Dict[str, List[Tuple[str, str, int]]] = defaultdict(list)
    group_profile: Dict[str, Dict[str, Dict[str, Dict[str, object]]]] = {
        group: {rank: {} for rank in RANKS} for group in GROUPS
    }
    total_length: Dict[str, int] = defaultdict(int)

    missing_taxid = 0
    missing_lineage = 0

    for header, sequence in stream_fasta(fasta_path):
        contig = contig_id_from_header(header)
        length = len(sequence)

        taxid = id_to_taxid.get(contig)
        if not taxid and "." in contig:
            taxid = id_to_taxid.get(contig.split(".", 1)[0])
        if not taxid:
            missing_taxid += 1
            continue

        taxids, names = build_lineage(taxid, tax_name, tax_parent, tax_rank)
        if all(t is None for t in taxids):
            missing_lineage += 1
            continue

        ancestry = compute_ancestry(taxid, tax_parent)
        group = assign_group(ancestry)
        if group not in GROUPS:
            group = "Protozoa"

        out = fasta_handles[group]
        out.write(f"{header}\n")
        for i in range(0, len(sequence), 80):
            out.write(sequence[i : i + 80] + "\n")

        contig_metadata[contig] = {
            "group": group,
            "length": length,
            "taxids": taxids,
            "names": names,
            "species_taxid": taxid,
            "ancestry": ancestry,
        }

        group_entries[group].append((contig, taxid, length))
        total_length[group] += length

        for idx, rank in enumerate(RANKS):
            rank_taxid = taxids[idx]
            if not rank_taxid:
                continue
            entry = group_profile[group][rank].setdefault(
                rank_taxid,
                {
                    "taxpath_ids": [tid or "" for tid in taxids],
                    "taxpath_names": names[:],
                    "length": 0,
                },
            )
            entry["length"] += length

    for handle in fasta_handles.values():
        handle.close()

    if missing_taxid:
        print(f"[prepare] warning: {missing_taxid} contigs lacked a taxid mapping and were skipped.")
    if missing_lineage:
        print(f"[prepare] warning: {missing_lineage} contigs lacked lineage information and were skipped.")

    truth_root = output_dir / "truth"
    for group in GROUPS:
        entries = group_entries[group]
        if not entries:
            continue
        group_dir = truth_root / group.replace(" ", "_").lower()
        write_truth_contigs(group_dir / "truth_contigs.tsv", entries)

        total_len = total_length[group]
        profile_rows: List[Tuple[str, str, List[str], List[str], float]] = []
        for rank in RANKS:
            for taxid, info in group_profile[group][rank].items():
                pct = (info["length"] / total_len) * 100 if total_len else 0.0
                profile_rows.append(
                    (
                        taxid,
                        rank,
                        [tid or "" for tid in info["taxpath_ids"]],
                        info["taxpath_names"],
                        pct,
                    )
                )
        write_truth_profile(group_dir / "truth_profile.cami.tsv", profile_rows)

    metadata_path = output_dir / "contig_metadata.json"
    ensure_dir(metadata_path.parent)
    with metadata_path.open("w") as fh:
        json.dump(contig_metadata, fh)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare HYMET mutation truth artefacts")
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--taxonomy", required=True, type=Path)
    parser.add_argument("--hierarchy", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    prepare_dataset(args.fasta, args.taxonomy, args.hierarchy, args.outdir)


if __name__ == "__main__":
    main()
