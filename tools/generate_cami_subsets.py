#!/usr/bin/env python3
"""Generate CAMI-style sample subsets from an existing CAMI sample.

This helper partitions the taxa present in the reference CAMI truth mapping
into contiguous groups (sorted by total assembled length) and writes out
derived CAMI inputs for each group.  Each derived sample includes:

- contigs FASTA
- truth contig mapping (gsa_mapping style TSV)
- CAMI-format taxonomic truth profile

The script is idempotent – existing output directories are replaced.
"""

from __future__ import annotations

import argparse
import collections
import os
import pathlib
from typing import Dict, Iterable, List, Tuple

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_ALIAS = {
    "superkingdom": "superkingdom",
    "domain": "superkingdom",
    "kingdom": "superkingdom",
    "sk": "superkingdom",
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
    "ss": "species",
}


def parse_nodes(path: pathlib.Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return mapping taxid -> parent_taxid and taxid -> rank."""
    parent: Dict[str, str] = {}
    rank: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = [part.strip() for part in line.split("|")]
            if len(parts) < 3:
                continue
            taxid, parent_id, rk = parts[0], parts[1], parts[2]
            parent[taxid] = parent_id
            rank[taxid] = rk
    return parent, rank


def parse_names(path: pathlib.Path) -> Dict[str, str]:
    """Return mapping taxid -> scientific name."""
    names: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = [part.strip() for part in line.split("|")]
            if len(parts) < 4:
                continue
            taxid, name_txt, _, cls = parts[:4]
            if cls == "scientific name" and taxid not in names:
                names[taxid] = name_txt
    return names


def lineage_ranks(
    taxid: str,
    parent: Dict[str, str],
    rank: Dict[str, str],
) -> Dict[str, str]:
    """Return mapping rank -> ancestor taxid (including the provided taxid)."""
    cache: Dict[str, str] = {}
    current = taxid
    seen = set()
    while current and current not in seen:
        seen.add(current)
        rk_raw = rank.get(current, "")
        rk = RANK_ALIAS.get(rk_raw.lower())
        if rk in RANKS and rk not in cache:
            cache[rk] = current
        nxt = parent.get(current)
        if not nxt or nxt == current:
            break
        current = nxt
    return cache


def build_taxpath(
    taxid: str,
    upto_rank: str,
    parent: Dict[str, str],
    rank: Dict[str, str],
    names: Dict[str, str],
    cache: Dict[str, Dict[str, str]],
) -> Tuple[str, str]:
    """Return (TAXPATH, TAXPATHSN) strings for a given rank."""
    if taxid == "0":
        count = RANKS.index(upto_rank) + 1
        ids = ["0"] * count
        nms = ["unclassified"] * count
        return "|".join(ids), "|".join(nms)
    if taxid not in cache:
        cache[taxid] = lineage_ranks(taxid, parent, rank)
    lineage = cache[taxid]
    ids: List[str] = []
    nms: List[str] = []
    for rk in RANKS:
        anc = lineage.get(rk)
        if anc:
            ids.append(anc)
            nms.append(names.get(anc, f"taxid_{anc}"))
        else:
            ids.append("0")
            nms.append("unclassified")
        if rk == upto_rank:
            break
    return "|".join(ids), "|".join(nms)


def ensure_clean_dir(path: pathlib.Path) -> None:
    if path.exists():
        for child in path.iterdir():
            if child.is_file() or child.is_symlink():
                child.unlink()
            else:
                ensure_clean_dir(child)
                child.rmdir()
    else:
        path.mkdir(parents=True, exist_ok=True)


def partition_taxa(
    mapping_path: pathlib.Path,
    partitions: Iterable[Tuple[str, int]],
) -> Tuple[List[Tuple[str, List[str]]], Dict[str, str], Dict[str, str], Dict[str, int], str]:
    """Partition taxids by descending assembled length according to `partitions`."""
    contig_rows: Dict[str, str] = {}
    contig_taxid: Dict[str, str] = {}
    contig_len: Dict[str, int] = {}
    taxid_lengths: Dict[str, int] = collections.defaultdict(int)
    with mapping_path.open("r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n")
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            contig_id = parts[0]
            taxid = parts[2]
            start = int(parts[5])
            end = int(parts[6])
            length = max(0, end - start + 1)
            contig_rows[contig_id] = line
            contig_taxid[contig_id] = taxid
            contig_len[contig_id] = length
            taxid_lengths[taxid] += length
    ordered = sorted(taxid_lengths.items(), key=lambda kv: kv[1], reverse=True)
    taxids_sorted = [taxid for taxid, _ in ordered]
    sample_defs: List[Tuple[str, List[str]]] = []
    taxid_to_sample: Dict[str, str] = {}
    idx = 0
    for sample_name, count in partitions:
        subset = taxids_sorted[idx : idx + count]
        if len(subset) < count:
            raise RuntimeError(
                f"Not enough taxa to satisfy partition for {sample_name}: requested {count}, got {len(subset)}."
            )
        idx += count
        sample_defs.append((sample_name, subset))
        for tid in subset:
            taxid_to_sample[tid] = sample_name
    return sample_defs, taxid_to_sample, contig_taxid, contig_len, contig_rows, header


def write_outputs(
    fasta_path: pathlib.Path,
    out_root: pathlib.Path,
    sample_defs: List[Tuple[str, List[str]]],
    taxid_to_sample: Dict[str, str],
    contig_taxid: Dict[str, str],
    contig_len: Dict[str, int],
    contig_rows: Dict[str, str],
    header_row: str,
    parent: Dict[str, str],
    rank: Dict[str, str],
    names: Dict[str, str],
) -> None:
    lineage_cache: Dict[str, Dict[str, str]] = {}
    samples = {
        name: {
            "taxids": set(taxids),
            "contigs": [],
            "total_length": 0,
            "rank_sums": {rk: collections.defaultdict(int) for rk in RANKS},
        }
        for name, taxids in sample_defs
    }

    # Assign contigs to samples and accumulate lengths per rank.
    for contig_id, taxid in contig_taxid.items():
        sample_name = taxid_to_sample.get(taxid)
        if not sample_name:
            continue
        entry = samples[sample_name]
        entry["contigs"].append(contig_id)
        length = contig_len[contig_id]
        entry["total_length"] += length
        lineage = lineage_ranks(taxid, parent, rank)
        lineage_cache[taxid] = lineage
        for rk in RANKS:
            anc = lineage.get(rk)
            key = anc if anc else "0"
            entry["rank_sums"][rk][key] += length

    # Prepare output directories and truth mapping.
    for sample_name in samples:
        sample_dir = out_root / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        mapping_path = sample_dir / "truth_contigs.tsv"
        with mapping_path.open("w", encoding="utf-8") as fout:
            fout.write(header_row + "\n")
            for contig_id in samples[sample_name]["contigs"]:
                fout.write(contig_rows[contig_id])

    # Write contig FASTA files.
    handles = {
        sample_name: (out_root / sample_name / "contigs.fna").open("w", encoding="utf-8")
        for sample_name in samples
    }
    try:
        with fasta_path.open("r", encoding="utf-8") as fin:
            current_id = None
            current_header = None
            seq_lines: List[str] = []
            for line in fin:
                if line.startswith(">"):
                    if current_id and current_header:
                        taxid = contig_taxid.get(current_id)
                        if taxid:
                            sample_name = taxid_to_sample.get(taxid)
                            if sample_name:
                                handles[sample_name].write(">" + current_header + "\n")
                                handles[sample_name].writelines(seq_lines)
                    current_header = line[1:].rstrip("\n")
                    current_id = current_header.split()[0]
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if current_id and current_header:
                taxid = contig_taxid.get(current_id)
                if taxid:
                    sample_name = taxid_to_sample.get(taxid)
                    if sample_name:
                        handles[sample_name].write(">" + current_header + "\n")
                        handles[sample_name].writelines(seq_lines)
    finally:
        for fh in handles.values():
            fh.close()

    # Write truth profiles.
    for sample_name, data in samples.items():
        sample_dir = out_root / sample_name
        profile_path = sample_dir / "truth_profile.tsv"
        total_length = data["total_length"] or 1  # avoid divide-by-zero
        with profile_path.open("w", encoding="utf-8") as fout:
            fout.write(f"@SampleID: {sample_name}\n")
            fout.write("@Version: 0.9.1\n")
            fout.write("@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n")
            fout.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_GENOMEID\n\n")
            for rk in RANKS:
                rank_totals = data["rank_sums"][rk]
                items = sorted(
                    rank_totals.items(),
                    key=lambda kv: (-kv[1], kv[0]),
                )
                for taxid, length in items:
                    if length <= 0:
                        continue
                    taxid_str = taxid if taxid != "0" else "0"
                    taxpath, taxpathsn = build_taxpath(
                        taxid_str,
                        rk,
                        parent,
                        rank,
                        names,
                        lineage_cache,
                    )
                    percentage = (length / total_length) * 100.0
                    fout.write(
                        f"{taxid_str}\t{rk}\t{taxpath}\t{taxpathsn}\t{percentage:.4f}\n"
                    )


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate derived CAMI sample subsets.")
    ap.add_argument("--fasta", default="/data/cami/sample_0.fna", help="Original contigs FASTA.")
    ap.add_argument(
        "--mapping",
        default="/data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/gsa_mapping_new.tsv",
        help="Original truth mapping TSV.",
    )
    ap.add_argument(
        "--nodes",
        default="/data/HYMET/taxonomy_files/nodes.dmp",
        help="NCBI taxonomy nodes.dmp",
    )
    ap.add_argument(
        "--names",
        default="/data/HYMET/taxonomy_files/names.dmp",
        help="NCBI taxonomy names.dmp",
    )
    ap.add_argument(
        "--outdir",
        default=str(pathlib.Path(__file__).resolve().parents[1] / "bench" / "data"),
        help="Output directory for derived samples.",
    )
    args = ap.parse_args()

    fasta_path = pathlib.Path(args.fasta)
    mapping_path = pathlib.Path(args.mapping)
    nodes_path = pathlib.Path(args.nodes)
    names_path = pathlib.Path(args.names)
    out_root = pathlib.Path(args.outdir)

    out_root.mkdir(parents=True, exist_ok=True)

    parent_map, rank_map = parse_nodes(nodes_path)
    names_map = parse_names(names_path)

    partitions = [
        ("cami_i_lc", 8),
        ("cami_i_mc", 12),
        ("cami_i_hc", 14),
        ("cami_ii_mousegut", 14),
        ("cami_ii_marine", 12),
        ("cami_ii_strainmadness", 12),
    ]

    (
        sample_defs,
        taxid_to_sample,
        contig_taxid,
        contig_len,
        contig_rows,
        header_row,
    ) = partition_taxa(mapping_path, partitions)

    write_outputs(
        fasta_path,
        out_root,
        sample_defs,
        taxid_to_sample,
        contig_taxid,
        contig_len,
        contig_rows,
        header_row,
        parent_map,
        rank_map,
        names_map,
    )

    for sample_name, taxids in sample_defs:
        print(f"[INFO] generated sample {sample_name} with {len(taxids)} taxa")


if __name__ == "__main__":
    main()
