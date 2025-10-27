# !/usr/bin/env python3
import os
import csv
from collections import defaultdict
import argparse
from multiprocessing import Pool
import logging
import sys

csv.field_size_limit(sys.maxsize)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_taxonomy_file(taxonomy_file):
    taxonomy = {}
    with open(taxonomy_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            taxid = row["TaxID"]
            for identifier in row["Identifiers"].split(";"):
                cleaned_id = identifier.strip()
                if cleaned_id:
                    taxonomy[cleaned_id] = taxid
    logging.info(f"Loaded {len(taxonomy)} taxonomy mappings")
    return taxonomy

def load_taxonomy_hierarchy_file(taxonomy_hierarchy_file):
    hierarchy = {}
    with open(taxonomy_hierarchy_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            taxid = row["TaxID"]
            hierarchy[taxid] = row["Lineage"].strip()
    logging.info(f"Loaded {len(hierarchy)} taxonomy hierarchies")
    return hierarchy

def parse_paf_file(paf_file):
    query_map = defaultdict(list)
    ref_counts = defaultdict(int)
    
    with open(paf_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
                
            query_id = parts[0]
            query_len = int(parts[1])
            ref_id = parts[5]
            align_len = int(parts[10])
            
            coverage = align_len / query_len if query_len > 0 else 0
            is_exact = (query_id == ref_id) and (coverage >= 0.99)
            
            query_map[query_id].append((ref_id, coverage, is_exact))
            ref_counts[ref_id] += 1

    logging.info(f"Processed {len(query_map)} queries from PAF file")
    return query_map, ref_counts

def determine_taxonomic_level(lineage):
    rank_order = [
        'superkingdom', 'phylum', 'class', 'order',
        'family', 'genus', 'species', 'strain'
    ]
    
    current_level = None
    for part in lineage.split(';'):
        part = part.strip()
        if ':' in part:
            rank, name = part.split(':', 1)
            rank = rank.strip().lower()
            if rank in rank_order:
                if current_level is None:
                    current_level = rank
                else:
                    current_index = rank_order.index(current_level)
                    new_index = rank_order.index(rank)
                    if new_index > current_index:
                        current_level = rank
    return current_level if current_level is not None else 'root'

def calculate_weighted_lineage(refs, ref_abundance, taxonomy):
    taxid_weights = defaultdict(float)
    total_weight = 0.0
    
    for ref_id, coverage, _ in refs:
        if ref_id not in taxonomy:
            continue
            
        taxid = taxonomy[ref_id]
        weight = coverage * ref_abundance.get(ref_id, 1)
        taxid_weights[taxid] += weight
        total_weight += weight
    
    return taxid_weights, total_weight

def determine_lca(taxid_weights, total_weight, taxonomy_hierarchy):
    if total_weight == 0:
        return "Unknown", "root", 0.0, None

    lineages = []
    for taxid, weight in taxid_weights.items():
        if taxid in taxonomy_hierarchy:
            lineage = taxonomy_hierarchy[taxid].split(";")
            lineages.append((lineage, weight/total_weight))

    if not lineages:
        return "Unknown", "root", 0.0, None

    consensus = {}
    confidence = 1.0
    rank_order = [
        'superkingdom', 'phylum', 'class', 'order',
        'family', 'genus', 'species', 'strain'
    ]
    best_taxid = None
    best_weight = -1.0

    for rank in rank_order:
        level_counts = defaultdict(float)
        for lineage, weight in lineages:
            for part in lineage:
                if part.startswith(f"{rank}:"):
                    level_counts[part] += weight
                    break
        
        if level_counts:
            best_match, conf = max(level_counts.items(), key=lambda x: x[1])
            consensus[rank] = best_match
            confidence *= conf
        else:
            break  # Stop at first missing rank

    # Determine representative taxid (highest weight among contributing taxids)
    for taxid, weight in taxid_weights.items():
        if taxid in taxonomy_hierarchy and weight > best_weight:
            best_taxid = taxid
            best_weight = weight

    lineage_parts = [consensus.get(rank) for rank in rank_order if consensus.get(rank)]
    if not lineage_parts:
        return "Unknown", "root", 0.0, best_taxid

    full_lineage = ";".join(lineage_parts)
    level = determine_taxonomic_level(full_lineage)
    return full_lineage, level, min(confidence, 1.0), best_taxid

def process_query(args):
    query, refs, ref_abundance, taxonomy, taxonomy_hierarchy = args
    
    # Check for exact matches first
    exact_matches = [ref for ref, _, is_exact in refs if is_exact and ref in taxonomy]
    if exact_matches:
        taxid = taxonomy[exact_matches[0]]
        if taxid in taxonomy_hierarchy:
            lineage = taxonomy_hierarchy[taxid]
            level = determine_taxonomic_level(lineage)
            return (query, lineage, level, taxid, 1.0)

    # Calculate LCA for non-exact matches
    taxid_weights, total_weight = calculate_weighted_lineage(refs, ref_abundance, taxonomy)
    lineage, level, confidence, representative_taxid = determine_lca(taxid_weights, total_weight, taxonomy_hierarchy)
    taxid_out = representative_taxid if representative_taxid is not None else "Unknown"
    
    return (query, lineage, level, taxid_out, confidence)

def main_process(paf_file, taxonomy_file, hierarchy_file, output_file, processes=4):
    taxonomy = load_taxonomy_file(taxonomy_file)
    taxonomy_hierarchy = load_taxonomy_hierarchy_file(hierarchy_file)
    query_map, ref_abundance = parse_paf_file(paf_file)

    tasks = [(query, refs, ref_abundance, taxonomy, taxonomy_hierarchy) 
             for query, refs in query_map.items()]

    results = []
    with Pool(processes) as pool:
        results = pool.map(process_query, tasks)

    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Query', 'Lineage', 'Taxonomic Level', 'TaxID', 'Confidence'])
        
        classified = 0
        for query, lineage, level, taxid, confidence in results:
            if lineage != 'Unknown' and taxid != "Unknown":
                classified += 1
            taxid_str = str(taxid) if taxid is not None else "Unknown"
            writer.writerow([query, lineage, level, taxid_str, f"{confidence:.4f}"])

    logging.info(f"Classification complete. Results saved to {output_file}")
    logging.info(f"Classified: {classified}/{len(results)} ({classified/len(results):.1%})")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Advanced LCA/Best Match Taxonomic Classifier")
    parser.add_argument("--paf", required=True, help="Input PAF file")
    parser.add_argument("--taxonomy", required=True, help="Taxonomy mapping file")
    parser.add_argument("--hierarchy", required=True, help="Taxonomy hierarchy file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--processes", type=int, default=4, help="Number of parallel processes")
    
    args = parser.parse_args()
    
    main_process(
        args.paf,
        args.taxonomy,
        args.hierarchy,
        args.output,
        args.processes
    )
