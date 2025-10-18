#!/usr/bin/env python3
import os
import re
import csv
import gzip
import argparse
import logging
import sys
from collections import defaultdict, Counter
from multiprocessing import Pool

# --- config / logging ---
csv.field_size_limit(1024 * 1024 * 1024)  # tolerate giant Identifier fields
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
RANK_ALIAS = {
    'domain': 'superkingdom', 'kingdom': 'superkingdom', 'sk': 'superkingdom', 'k': 'superkingdom',
    'phylum': 'phylum', 'p': 'phylum',
    'class': 'class', 'c': 'class',
    'order': 'order', 'o': 'order',
    'family': 'family', 'f': 'family',
    'genus': 'genus', 'g': 'genus',
    'species': 'species', 's': 'species',
    'subspecies': 'strain', 'ss': 'strain', 'strain': 'strain'
}
GCFA_RE = re.compile(r'GC[AF]_\d+(?:\.\d+)?(?:_PRJ[A-Z]+\d+)?')
ACC_RE = re.compile(r'(NC_\d+\.\d+|NZ_[A-Z]{2}\d+\.\d+|NZ_[A-Z]{5}\d+\.\d+|CP\d+\.\d+|CM\d+\.\d+|[A-Z]{2}_\d+\.\d+)')

# --- globals for worker processes ---
_TAX = None                 # dict: identifier -> taxid
_HIER = None                # dict: taxid -> tuple/list of names by RANKS
_REF_ABUND = None           # dict: ref_id -> count

def _init_worker(tax_map, hier_map, ref_abund):
    global _TAX, _HIER, _REF_ABUND
    _TAX = tax_map
    _HIER = hier_map
    _REF_ABUND = ref_abund

# -------- taxonomy loaders --------

def _add_token(m, tok, taxid):
    """Add token and versionless variants into map if not present."""
    if not tok:
        return
    tok = tok.strip()
    if not tok:
        return
    m.setdefault(tok, taxid)
    # versionless (e.g., NC_014743.1 -> NC_014743 ; GCF_000001.1 -> GCF_000001)
    if '.' in tok:
        m.setdefault(tok.split('.', 1)[0], taxid)

def _split_identifiers(s):
    """Split identifiers field on common separators; keep meaningful tokens."""
    if not s:
        return []
    # split on ; | , whitespace
    parts = re.split(r'[;|,\s]+', s)
    return [p for p in (x.strip() for x in parts) if p]

def load_taxonomy_file(taxonomy_file):
    """
    Build identifier -> TaxID map from detailed_taxonomy.tsv.
    Uses columns:
      - TaxID (required)
      - Identifiers (tokens split; also versionless)
      - Any field containing a GCF/GCA pattern (added as keys)
    """
    m = {}
    with open(taxonomy_file, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'TaxID' not in reader.fieldnames:
            raise RuntimeError("TaxID column not found in taxonomy file")
        for row in reader:
            taxid = (row.get('TaxID') or '').strip()
            if not taxid:
                continue

            # 1) capture any GCF/GCA-like accessions present in any column
            for v in row.values():
                if not v:
                    continue
                for acc in GCFA_RE.findall(v):
                    _add_token(m, acc, taxid)

            # 2) Identifiers column (broad token list)
            ids = row.get('Identifiers') or ''
            for tok in _split_identifiers(ids):
                _add_token(m, tok, taxid)

            # 3) Also capture contig-style accessions embedded in Identifiers or other fields
            #    (we already split 'ids', but patterns may be embedded)
            for v in (ids,) + tuple(row.get(k) or '' for k in row.keys()):
                if not v:
                    continue
                for mm in ACC_RE.findall(v):
                    _add_token(m, mm, taxid)

    logging.info(f"Loaded {len(m):,} taxonomy mappings")
    return m

def _parse_lineage_to_names(lineage_raw):
    """
    Convert various lineage encodings to a normalized list of names per RANKS.
    Accepted forms:
      - 'rank:name; rank:name; ...'
      - 'k__Bacteria; p__Firmicutes; ...'
      - 'Bacteria; Firmicutes; Bacilli; ...'
      - 'name1|name2|...' (pipe-separated)
    """
    names_by_rank = [''] * len(RANKS)
    if not lineage_raw:
        return names_by_rank

    s = lineage_raw.strip()

    # Case A: rank:name pairs
    if ':' in s:
        parts = re.split(r'[;|]+', s)
        for part in parts:
            part = part.strip()
            if not part or ':' not in part:
                continue
            rk, nm = part.split(':', 1)
            rk = RANK_ALIAS.get(rk.strip().lower(), None)
            nm = nm.strip()
            if not rk or not nm:
                continue
            idx = RANKS.index(rk)
            names_by_rank[idx] = nm
        return names_by_rank

    # Case B: k__/p__ style labels
    if '__' in s:
        parts = re.split(r'[;|]+', s)
        for part in parts:
            part = part.strip()
            if not part or '__' not in part:
                continue
            rk_tag, nm = part.split('__', 1)
            rk = RANK_ALIAS.get(rk_tag.strip().lower(), None)
            nm = nm.strip()
            if not rk or not nm:
                continue
            idx = RANKS.index(rk)
            names_by_rank[idx] = nm
        return names_by_rank

    # Case C: plain names, assume ordered from superkingdom downward
    parts = re.split(r'[;|]+', s)
    seq = [p.strip() for p in parts if p.strip() and p.strip().upper() != 'NA']
    for i, nm in enumerate(seq[:len(RANKS)]):
        names_by_rank[i] = nm
    return names_by_rank

def load_taxonomy_hierarchy_file(taxonomy_hierarchy_file):
    """
    Map TaxID -> names-by-rank (list aligned to RANKS).
    """
    hierarchy = {}
    with open(taxonomy_hierarchy_file, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'TaxID' not in reader.fieldnames or 'Lineage' not in reader.fieldnames:
            raise RuntimeError("Hierarchy file must have TaxID and Lineage columns")
        for row in reader:
            tid = (row.get('TaxID') or '').strip()
            lin = (row.get('Lineage') or '').strip()
            if not tid:
                continue
            hierarchy[tid] = _parse_lineage_to_names(lin)
    logging.info(f"Loaded {len(hierarchy):,} taxonomy hierarchies")
    return hierarchy

# -------- PAF parsing --------

def _opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

def parse_paf_file(paf_file):
    """
    Return:
      query_map: dict q -> list of (tname, coverage)
      ref_counts: dict tname -> total alignments count (for weighting)
    """
    query_map = defaultdict(list)
    ref_counts = defaultdict(int)
    with _opener(paf_file) as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 11:
                continue
            qname = parts[0]
            try:
                qlen = int(parts[1])
                aln_block = int(parts[10])  # PAF col 11: alignment block length
            except Exception:
                qlen = 0
                aln_block = 0
            tname = parts[5]
            cov = (aln_block / qlen) if qlen > 0 else 0.0
            query_map[qname].append((tname, cov))
            ref_counts[tname] += 1
    logging.info(f"Processed {len(query_map):,} queries from PAF file")
    return query_map, ref_counts

# -------- classification core --------

def _generate_lookup_candidates(tname):
    """
    Produce a ranked list of candidate keys to look up in _TAX.
    Includes:
      - original tname
      - versionless
      - first token before '|' and its versionless
      - any embedded GCFA/ACC patterns
    """
    cands = []
    def add(x):
        if x and x not in cands:
            cands.append(x)
        if x and '.' in x:
            xv = x.split('.', 1)[0]
            if xv not in cands:
                cands.append(xv)

    add(tname)
    # split by whitespace and pipe
    head = re.split(r'[|\s]+', tname)[0]
    add(head)

    # embedded accession patterns
    for g in GCFA_RE.findall(tname):
        add(g)
    for a in ACC_RE.findall(tname):
        add(a)

    return cands

def _lookup_taxid(tname):
    """Try multiple normalized forms against the taxonomy map."""
    for cand in _generate_lookup_candidates(tname):
        tid = _TAX.get(cand)
        if tid:
            return tid
    return None

def _weighted_lca(taxid_weights):
    """
    taxid_weights: dict taxid -> weight
    Use _HIER (taxid -> names-by-rank) to compute a weighted consensus lineage.
    Returns: (lineage_str, level, confidence, representative_taxid)
    """
    total_w = sum(taxid_weights.values())
    if total_w <= 0:
        return "Unknown", "root", 0.0, None

    chosen = []
    conf_product = 1.0
    best_taxid = None
    best_weight = -1.0

    for r_idx, rank in enumerate(RANKS):
        # gather weights per name at this rank
        name_w = defaultdict(float)
        denom = 0.0
        for tid, w in taxid_weights.items():
            names = _HIER.get(tid)
            if not names:
                continue
            nm = names[r_idx] if r_idx < len(names) else ''
            if nm:
                name_w[nm] += w
                denom += w
        if denom <= 0 or not name_w:
            break
        best_name, best_w = max(name_w.items(), key=lambda kv: kv[1])
        conf_i = best_w / denom
        chosen.append(best_name)
        conf_product *= conf_i

    if not chosen:
        return "Unknown", "root", 0.0, None

    lineage_str = "; ".join(f"{RANKS[i]}:{name}" for i, name in enumerate(chosen))
    level = RANKS[len(chosen) - 1]
    for tid, weight in taxid_weights.items():
        if weight > best_weight:
            best_taxid = tid
            best_weight = weight
    return lineage_str, level, min(conf_product, 1.0), best_taxid

def _process_one(task):
    """
    task: (query, refs) where refs = [(tname, cov), ...]
    Uses globals _TAX, _HIER, _REF_ABUND.
    """
    q, refs = task
    tw = defaultdict(float)
    any_hit = False
    for tname, cov in refs:
        tid = _lookup_taxid(tname)
        if not tid:
            continue
        any_hit = True
        w = cov * _REF_ABUND.get(tname, 1)
        tw[tid] += w
    if not any_hit:
        return (q, "Unknown", "root", "Unknown", 0.0)
    lineage, level, conf, rep_tid = _weighted_lca(tw)
    tid_out = rep_tid if rep_tid is not None else "Unknown"
    return (q, lineage, level, tid_out, conf)

# -------- main driver --------

def main_process(paf_file, taxonomy_file, hierarchy_file, output_file, processes=4):
    tax_map = load_taxonomy_file(taxonomy_file)
    hier_map = load_taxonomy_hierarchy_file(hierarchy_file)
    query_map, ref_abund = parse_paf_file(paf_file)

    tasks = list(query_map.items())
    results = []
    with Pool(processes, initializer=_init_worker, initargs=(tax_map, hier_map, ref_abund)) as pool:
        for res in pool.imap_unordered(_process_one, tasks, chunksize=200):
            results.append(res)

    classified = sum(1 for _, lin, tid, _, _ in results if lin != "Unknown" and tid != "Unknown")
    resdict = {q: (lin, lvl, tid, conf) for q, lin, lvl, tid, conf in results}
    with open(output_file, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Query', 'Lineage', 'Taxonomic Level', 'TaxID', 'Confidence'])
        for q in query_map.keys():
            lin, lvl, tid, conf = resdict.get(q, ("Unknown", "root", "Unknown", 0.0))
            w.writerow([q, lin, lvl, tid, f"{conf:.4f}"])

    total = len(results)
    logging.info(f"Classification complete. Results saved to {output_file}")
    logging.info(f"Classified: {classified}/{total} ({(classified/total if total else 0):.1%})")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="HYMET taxonomy classifier (robust identifiers + weighted LCA)")
    p.add_argument("--paf", required=True, help="Input PAF file (.paf or .paf.gz)")
    p.add_argument("--taxonomy", required=True, help="detailed_taxonomy.tsv")
    p.add_argument("--hierarchy", required=True, help="taxonomy_hierarchy.tsv")
    p.add_argument("--output", required=True, help="Output TSV")
    p.add_argument("--processes", type=int, default=4, help="Parallel workers")
    args = p.parse_args()

    main_process(args.paf, args.taxonomy, args.hierarchy, args.output, args.processes)
