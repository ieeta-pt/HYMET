#!/usr/bin/env python3
import os, sys, csv, argparse, math, collections, subprocess, pathlib, re, hashlib
csv.field_size_limit(1024*1024*1024)

# -------- Defaults bound to your layout --------
DEF_PRED_PROFILE = "/data/hymet_out/sample_0/hymet.sample_0.cami.tsv"
DEF_TRUTH_PROFILE= "/data/cami/sample_0/taxonomic_profile_0.txt"
DEF_PRED_CONTIGS = "/data/hymet_out/sample_0/work/classified_sequences.tsv"  # HYMET per-contig predictions
DEF_TRUTH_CONTIGS = "/data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/gsa_mapping_new.tsv"
DEF_TRUTH_CONTIGS_FALLBACK = "/data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/gsa_mapping.tsv"
DEF_PRED_FASTA = "/data/cami/sample_0.fna"  # the FASTA HYMET was run on
DEF_TRUTH_FASTA= "/data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/anonymous_gsa.fasta"
DEF_TAXDB = "/data/HYMET/taxonomy_files"
DEF_OUTDIR = "/data/hymet_out/sample_0/eval"
RANKS = ["superkingdom","phylum","class","order","family","genus","species"]
RANKC = ["k","p","c","o","f","g","s"]

# ----------------- utils -----------------
def ensure_dir(p):
    pathlib.Path(p).mkdir(parents=True, exist_ok=True)

def is_num(s): 
    return bool(re.fullmatch(r"[0-9]+", s or ""))

def chunked(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i+n]

def fasta_lengths(paths):
    lens = {}
    for path in paths:
        if not path or not os.path.isfile(path):
            continue
        name=None; L=0
        with open(path) as f:
            for ln in f:
                if ln.startswith(">"):
                    if name is not None: lens.setdefault(name, L)
                    name = ln[1:].strip().split()[0]; L = 0
                else:
                    L += len(ln.strip())
        if name is not None: lens.setdefault(name, L)
    return lens

def fasta_hashes(path):
    """Return dict: contig_id -> md5(sequence)."""
    hmap = {}
    if not path or not os.path.isfile(path): 
        return hmap
    name=None; chunks=[]
    md = None
    with open(path) as f:
        for ln in f:
            if ln.startswith(">"):
                if name is not None:
                    md = hashlib.md5()
                    for c in chunks: md.update(c)
                    hmap[name] = md.hexdigest()
                name = ln[1:].strip().split()[0]
                chunks = []
            else:
                s = ln.strip().encode("utf-8")
                if s: chunks.append(s)
        if name is not None:
            md = hashlib.md5()
            for c in chunks: md.update(c)
            hmap[name] = md.hexdigest()
    return hmap

def taxonkit_name2taxid(names, taxdb):
    out = {}
    if not names: return out
    env = os.environ.copy(); env["TAXONKIT_DB"] = taxdb
    for ch in chunked(list(names), 50000):
        proc = subprocess.run(
            ["taxonkit","name2taxid","--data-dir",taxdb,"--show-rank"],
            input="\n".join(ch)+"\n", text=True, capture_output=True, check=True, env=env
        )
        for line in proc.stdout.strip().splitlines():
            ps = line.split("\t")
            if len(ps) >= 2 and is_num(ps[1]):
                out[ps[0]] = ps[1]
    return out

def taxonkit_taxpath(taxids, taxdb):
    out = {}
    if not taxids: return out
    env = os.environ.copy(); env["TAXONKIT_DB"] = taxdb
    for ch in chunked(list(taxids), 50000):
        proc = subprocess.run(
            ["taxonkit","reformat","--data-dir",taxdb,"-I","1","-f","{d}|{p}|{c}|{o}|{f}|{g}|{s}","-t","-T"],
            input="\n".join(ch)+"\n", text=True, capture_output=True, check=True, env=env
        )
        for line in proc.stdout.strip().splitlines():
            ps = line.split("\t")
            if len(ps) >= 3:
                tid, names, ids = ps[0], ps[1], ps[2]
                def _norm(path):
                    parts = (path or "").split("|") if path else []
                    if len(parts) < len(RANKS):
                        parts.extend(["NA"] * (len(RANKS) - len(parts)))
                    return "|".join(parts[:len(RANKS)])
                out[tid] = (_norm(ids), _norm(names))
    return out

# -------------- profiles (CAMI) ---------------
def _parse_cami_like(lines, taxdb):
    """Robust CAMI profile parser. Returns dict[rank]->Counter(taxid->percent)."""
    prof = {r: collections.Counter() for r in RANKS}
    # try standard CAMI: TAXID RANK TAXPATH TAXPATHSN PERCENTAGE
    for ln in lines:
        if not ln.strip() or ln[0] in "#@": 
            continue
        ps = ln.rstrip("\n").split("\t")
        if len(ps) >= 5 and is_num(ps[0]):   # looks like standard
            tid, rank, perc = ps[0], ps[1].strip().lower(), ps[4]
            if rank in prof:
                try: prof[rank][tid] += float(perc)
                except: pass
            continue
        # alt formats: try to locate headers
        break
    if any(prof[r] for r in RANKS):
        return prof  # parsed ok

    # try with header row and different column names
    rdr = csv.reader([ln for ln in lines if ln.strip() and ln[0] not in "#@"], delimiter="\t")
    try:
        hdr = next(rdr)
    except StopIteration:
        return prof
    h = [c.strip().lower() for c in hdr]
    def idx(*names):
        for n in names:
            if n in h: return h.index(n)
        return -1
    i_taxid = idx("taxid","taxon_id","ncbi_taxid","ncbi_tax_id")
    i_rank  = idx("rank")
    i_perc  = idx("percentage","abundance","rel_abundance","fraction_total_reads")
    i_taxpath = idx("taxpath")
    i_taxpathsn = idx("taxpathsn","taxpath_sn","taxpath_names","lineage")

    rows = list(rdr)
    if i_taxid >= 0 and i_rank >=0 and i_perc >= 0:
        for ps in rows:
            try:
                tid = ps[i_taxid].strip()
                rk  = ps[i_rank].strip().lower()
                val = float(ps[i_perc]) * (100.0 if "abundance" in h[i_perc] or "fraction" in h[i_perc] else 1.0)
                if rk in prof and is_num(tid):
                    prof[rk][tid] += val
            except: 
                continue
        return prof

    # derive TAXID from TAXPATHSN (names) or TAXPATH (ids)
    if i_rank >= 0 and (i_taxpath >= 0 or i_taxpathsn >= 0):
        rk_to_idx = dict(zip(RANKS, range(len(RANKS))))
        name2tid = {}
        for ps in rows:
            try:
                rk = ps[i_rank].strip().lower()
                if rk not in rk_to_idx: continue
                if i_perc < 0: continue
                val = float(ps[i_perc]) * (100.0 if "abundance" in h[i_perc] or "fraction" in h[i_perc] else 1.0)
                if i_taxpath >= 0:
                    path = ps[i_taxpath].strip().split("|")
                    if path and all((x.isdigit() or x=="NA") for x in path):
                        idx_r = rk_to_idx[rk]
                        if idx_r < len(path) and path[idx_r]!="NA":
                            prof[rk][path[idx_r]] += val
                        continue
                if i_taxpathsn >= 0:
                    pathn = [p.strip() for p in ps[i_taxpathsn].split("|")]
                    idx_r = rk_to_idx[rk]
                    if idx_r < len(pathn):
                        nm = pathn[idx_r]
                        if nm: name2tid[nm] = None
            except:
                continue
        if name2tid:
            mapped = taxonkit_name2taxid(list(name2tid.keys()), taxdb)
            for k,v in mapped.items(): name2tid[k]=v
            for ps in rows:
                try:
                    rk = ps[i_rank].strip().lower()
                    if rk not in rk_to_idx: continue
                    val = float(ps[i_perc]) * (100.0 if "abundance" in h[i_perc] or "fraction" in h[i_perc] else 1.0)
                    pathn = [p.strip() for p in ps[i_taxpathsn].split("|")] if i_taxpathsn>=0 else []
                    idx_r = rk_to_idx[rk]
                    if idx_r < len(pathn):
                        nm = pathn[idx_r]
                        tid = name2tid.get(nm)
                        if tid and is_num(tid):
                            prof[rk][tid] += val
                except:
                    continue
    return prof

def load_profile_any(path, taxdb):
    if not os.path.isfile(path):
        return {r: collections.Counter() for r in RANKS}
    with open(path) as f:
        lines = f.readlines()
    return _parse_cami_like(lines, taxdb)

# -------------- contig mapping ----------------
def load_pred_contigs(pred_file):
    """Return dict contig -> last-name-from-Lineage (string)."""
    out = {}
    if not os.path.isfile(pred_file):
        return out
    with open(pred_file) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            q = row.get("Query") or row.get("qname") or row.get("q")
            lin = (row.get("Lineage","") or "").strip()
            if not q or not lin:
                continue
            last = lin.split(";")[-1].strip()
            if ":" in last: last = last.split(":",1)[1].strip()
            if last: out[q] = last
    return out

def load_gt_contigs(gt_file, taxdb):
    """Return dict contig -> TaxID (str). Robust to column names."""
    out = {}
    if not os.path.isfile(gt_file):
        return out
    with open(gt_file) as f:
        rdr = csv.reader(f, delimiter="\t")
        hdr = next(rdr)
        h = [c.strip().lower() for c in hdr]
        # contig/sequence column
        cand_contig = ["contig","contig_name","sequence","sequence_name","name","id"]
        ci = next((h.index(c) for c in cand_contig if c in h), 0)
        # prefer explicit taxid columns
        cand_taxid = ["taxid","taxon_id","ncbi_tax_id","ncbi_taxid","species_taxid"]
        ti = next((h.index(c) for c in cand_taxid if c in h), -1)
        # taxonomy names/path fallbacks
        ti_taxonomy = next((h.index(c) for c in ["taxonomy","taxpathsn","lineage","taxpath_names","clade_name","organism_name"] if c in h), -1)
        ti_taxpath  = next((h.index(c) for c in ["taxpath","taxid_path","ncbi_taxpath"] if c in h), -1)

        rows = list(rdr)
        if ti >= 0:
            for ps in rows:
                if len(ps) <= max(ci,ti): continue
                cont, val = ps[ci], ps[ti].strip()
                if is_num(val): out[cont] = val
            return out

        if ti_taxpath >= 0:
            for ps in rows:
                if len(ps) <= max(ci,ti_taxpath): continue
                cont, path = ps[ci], ps[ti_taxpath].strip()
                ids = [x for x in path.split("|") if x and x!="NA"]
                if ids and is_num(ids[-1]): out[cont] = ids[-1]
            if out: return out

        if ti_taxonomy >= 0:
            names = set(); vals = {}
            for ps in rows:
                if len(ps) <= max(ci,ti_taxonomy): continue
                cont, tax = ps[ci], ps[ti_taxonomy].strip()
                if not tax: continue
                parts = re.split(r"\|+|;+", tax)
                nm = parts[-1].strip()
                vals[cont] = nm; names.add(nm)
            mapnm = taxonkit_name2taxid(names, taxdb)
            for cont, nm in vals.items():
                tid = mapnm.get(nm)
                if tid and is_num(tid): out[cont] = tid
            return out

        for ps in rows:
            if len(ps) < 2: continue
            if is_num(ps[1]): out[ps[0]] = ps[1]
    return out

# ---------- build profiles from contigs ----------
def profiles_from_contig_maps(contig2tid, lengths, taxdb):
    """Return dict[rank]->Counter(taxid->percent) weighted by contig length."""
    prof = {r: collections.Counter() for r in RANKS}
    if not contig2tid:
        return prof
    tids = set(contig2tid.values())
    paths = taxonkit_taxpath(tids, taxdb)
    acc = collections.Counter()
    for cont, tid in contig2tid.items():
        w = lengths.get(cont, 1)
        names_ids = paths.get(tid)
        if not names_ids: 
            continue
        ids = names_ids[0].split("|")
        for i, code in enumerate(RANKC):
            if i < len(ids) and ids[i] != "NA":
                prof[RANKS[i]][ids[i]] += w
                acc[RANKS[i]] += w
    for r in RANKS:
        s = acc[r]
        if s > 0:
            for k in list(prof[r].keys()):
                prof[r][k] = 100.0 * prof[r][k] / s
    return prof

# ---------------- metrics ----------------
def l1_and_braycurtis(a: dict, b: dict):
    keys = set(a) | set(b)
    if not keys: return 0.0, 0.0
    sum_abs = sum(abs(a.get(k,0.0)-b.get(k,0.0)) for k in keys)
    l1 = 0.5 * sum_abs
    sump = sum(a.get(k,0.0) for k in keys)
    sumt = sum(b.get(k,0.0) for k in keys)
    shared = sum(min(a.get(k,0.0), b.get(k,0.0)) for k in keys)
    bc = 1.0 - (2.0*shared / (sump + sumt if (sump+sumt)>0 else 1.0))
    return l1, bc*100.0

def prf_presence(a: dict, b: dict, thr=0.1):
    A = {k for k,v in a.items() if v >= thr}
    B = {k for k,v in b.items() if v >= thr}
    tp = len(A & B); fp = len(A - B); fn = len(B - A)
    prec = tp/(tp+fp) if (tp+fp)>0 else 0.0
    rec  = tp/(tp+fn) if (tp+fn)>0 else 0.0
    f1   = 2*prec*rec/(prec+rec) if (prec+rec)>0 else 0.0
    return prec*100.0, rec*100.0, f1*100.0, tp, fp, fn

# -------------- contig-level eval --------------
def eval_contigs(pred_file, gt_file, taxdb, outdir, pred_fasta=None, gt_fasta=None):
    # predicted: contig -> last-name
    pred_name = load_pred_contigs(pred_file)
    gt_map = load_gt_contigs(gt_file, taxdb)

    # map predicted names -> TaxID
    mapped = taxonkit_name2taxid(set(pred_name.values()), taxdb)
    pred_tid = {q: mapped.get(nm) for q, nm in pred_name.items() if mapped.get(nm)}

    # direct name intersection
    pairs = [(q, pred_tid[q], gt_map[q]) for q in pred_tid.keys() if q in gt_map]
    if not pairs and pred_fasta and gt_fasta and os.path.isfile(pred_fasta) and os.path.isfile(gt_fasta):
        # build cross-map via MD5 sequence hashes
        pred_hash = fasta_hashes(pred_fasta)
        gt_hash = fasta_hashes(gt_fasta)
        # invert GT hashes → contig name
        inv_gt = collections.defaultdict(list)
        for gname, h in gt_hash.items():
            inv_gt[h].append(gname)
        cross = {}
        for q in pred_tid.keys():
            h = pred_hash.get(q)
            if not h: continue
            names = inv_gt.get(h, [])
            if not names: 
                continue
            # if multiple contigs share the same hash, pick the first deterministically
            cross[q] = sorted(names)[0]
        # remap with cross-name
        pairs = [(q, pred_tid[q], gt_map.get(cross[q])) for q in cross.keys() if gt_map.get(cross[q])]

    usable = len(pairs)
    exact = sum(1 for _,pt,gtid in pairs if pt == gtid)

    # per-rank accuracy
    tids = {pt for _,pt,_ in pairs} | {gtid for *_,gtid in pairs}
    tpaths = taxonkit_taxpath(tids, taxdb)

    per_rank = {}
    for i, r in enumerate(RANKS):
        tot = 0; ok = 0
        for _, pt, gtid in pairs:
            pids = tpaths.get(pt, ("",""))[0]
            gids = tpaths.get(gtid,("",""))[0]
            if not pids or not gids: 
                continue
            pvec = pids.split("|"); gvec = gids.split("|")
            if i >= len(pvec) or i >= len(gvec): 
                continue
            pid = pvec[i]; gid = gvec[i]
            if pid=="NA" or gid=="NA": continue
            tot += 1
            if pid == gid: ok += 1
        acc = 100.0*ok/tot if tot else 0.0
        per_rank[r] = {"n": tot, "acc": acc, "correct": ok}

    # write files
    with open(os.path.join(outdir,"contigs_exact.tsv"),"w",newline="") as w:
        wr=csv.writer(w, delimiter="\t")
        wr.writerow(["metric","value"])
        wr.writerow(["usable_pairs", usable])
        wr.writerow(["exact_taxid_matches", exact])
        wr.writerow(["exact_taxid_accuracy_percent", 100.0*exact/usable if usable else 0.0])

    with open(os.path.join(outdir,"contigs_per_rank.tsv"),"w",newline="") as w:
        wr=csv.writer(w, delimiter="\t")
        wr.writerow(["rank","n","correct","accuracy_percent"])
        for r in RANKS:
            m=per_rank.get(r,{"n":0,"correct":0,"acc":0.0})
            wr.writerow([r, m["n"], m["correct"], f"{m['acc']:.4f}"])
    return {"usable_pairs": usable, "exact": exact, "per_rank": per_rank}

# ------------------- main -------------------
def main():
    ap = argparse.ArgumentParser(description="Evaluate HYMET vs CAMI ground truth (robust, name+sequence reconciliation).")
    ap.add_argument("--pred-profile", default=DEF_PRED_PROFILE, help="Predicted CAMI profile TSV")
    ap.add_argument("--truth-profile", default=DEF_TRUTH_PROFILE, help="Ground-truth CAMI profile TSV")
    ap.add_argument("--pred-contigs", default=DEF_PRED_CONTIGS, help="Predicted per-contig file (classified_sequences.tsv)")
    ap.add_argument("--truth-contigs", default="", help="Ground-truth contig mapping TSV (gsa_mapping*.tsv)")
    ap.add_argument("--pred-fasta", default=DEF_PRED_FASTA, help="FASTA used by HYMET (for MD5 reconciliation)")
    ap.add_argument("--truth-fasta", default=DEF_TRUTH_FASTA, help="Ground-truth contigs FASTA (for MD5 reconciliation)")
    ap.add_argument("--taxdb", default=DEF_TAXDB, help="TaxonKit DB dir")
    ap.add_argument("--outdir", default=DEF_OUTDIR, help="Output directory for evaluation")
    ap.add_argument("--presence-thresh", type=float, default=0.1, help="Presence threshold (%%) for P/R/F1 on profiles")
    args = ap.parse_args()

    ensure_dir(args.outdir)
    gt_contigs_path = args.truth_contigs or (DEF_TRUTH_CONTIGS if os.path.isfile(DEF_TRUTH_CONTIGS) else DEF_TRUTH_CONTIGS_FALLBACK)

    # --------- Try to load profiles directly ----------
    pred_prof = load_profile_any(args.pred_profile, args.taxdb)
    truth_prof = load_profile_any(args.truth_profile, args.taxdb)

    # --------- Fallback: rebuild profiles from contigs ----------
    need_pred_fb  = all(not pred_prof[r] for r in RANKS)
    need_truth_fb = all(not truth_prof[r] for r in RANKS)

    lens = {}
    if need_pred_fb or need_truth_fb:
        lens = fasta_lengths([args.pred_fasta, args.truth_fasta])

    if need_pred_fb:
        pred_name = load_pred_contigs(args.pred_contigs)
        mapped = taxonkit_name2taxid(set(pred_name.values()), args.taxdb)
        pred_tid = {q: mapped.get(nm) for q, nm in pred_name.items() if mapped.get(nm)}
        pred_prof = profiles_from_contig_maps(pred_tid, lens, args.taxdb)

    if need_truth_fb:
        gt_map = load_gt_contigs(gt_contigs_path, args.taxdb)
        truth_prof = profiles_from_contig_maps(gt_map, lens, args.taxdb)

    # --------- Profile metrics ----------
    def write_rank_diffs(outdir, rank, a: dict, b: dict, top=200):
        p = pathlib.Path(outdir, f"profile_diffs_{rank}.tsv")
        rows = []
        keys = set(a)|set(b)
        for k in keys:
            rows.append((k, a.get(k,0.0), b.get(k,0.0), a.get(k,0.0)-b.get(k,0.0)))
        rows.sort(key=lambda x: abs(x[3]), reverse=True)
        with p.open("w", newline="") as w:
            wr = csv.writer(w, delimiter="\t")
            wr.writerow(["taxid","pred_percent","truth_percent","delta_pred_minus_truth"])
            for r in rows[:top]:
                wr.writerow(r)

    summ = []
    summ.append("# Profile-level metrics (per rank)")
    with open(os.path.join(args.outdir,"profile_summary.tsv"),"w",newline="") as w:
        wr = csv.writer(w, delimiter="\t")
        wr.writerow(["rank","L1_total_variation_pctpts","BrayCurtis_pct","Precision_%","Recall_%","F1_%","TP","FP","FN"])
        for r in RANKS:
            l1, bc = l1_and_braycurtis(pred_prof[r], truth_prof[r])
            pr, rc, f1, tp, fp, fn = prf_presence(pred_prof[r], truth_prof[r], args.presence_thresh)
            wr.writerow([r, f"{l1:.4f}", f"{bc:.4f}", f"{pr:.2f}", f"{rc:.2f}", f"{f1:.2f}", tp, fp, fn])
            summ.append(f"{r:14s}  L1={l1:.3f}  BC={bc:.3f}%  P/R/F1={pr:.1f}/{rc:.1f}/{f1:.1f}% (TP={tp}, FP={fp}, FN={fn})")
            write_rank_diffs(args.outdir, r, pred_prof[r], truth_prof[r], top=200)

    # --------- Contig-level metrics (with MD5 reconciliation) ----------
    summ.append("\n# Contig-level accuracy")
    cres = eval_contigs(args.pred_contigs, gt_contigs_path, args.taxdb, args.outdir,
                        pred_fasta=args.pred_fasta, gt_fasta=args.truth_fasta)
    usable = cres.get("usable_pairs", 0); exact = cres.get("exact", 0)
    summ.append(f"Exact TaxID: {exact}/{usable} ({(100.0*exact/usable if usable else 0.0):.2f}%)")
    for r in RANKS:
        m = cres.get("per_rank", {}).get(r, {"n":0,"acc":0.0,"correct":0})
        summ.append(f"{r:14s}  n={m['n']:<8d}  acc={m['acc']:.2f}%")

    with open(os.path.join(args.outdir,"summary.txt"),"w") as w:
        w.write("\n".join(summ)+"\n")
    print("\n".join(summ))

    with open(os.path.join(args.outdir,"_debug_info.txt"),"w") as w:
        w.write(f"pred_profile_path: {args.pred_profile}\n")
        w.write(f"truth_profile_path: {args.truth_profile}\n")
        w.write(f"pred_contigs_path: {args.pred_contigs}\n")
        w.write(f"truth_contigs_path: {gt_contigs_path}\n")
        w.write(f"pred_fasta: {args.pred_fasta}\n")
        w.write(f"truth_fasta: {args.truth_fasta}\n")
        w.write(f"taxdb: {args.taxdb}\n")
        w.write(f"presence_thresh: {args.presence_thresh}\n")

    print(f"\n[WROTE] {os.path.join(args.outdir,'summary.txt')}")
    print(f"[WROTE] {os.path.join(args.outdir,'profile_summary.tsv')}")
    print(f"[WROTE] {os.path.join(args.outdir,'contigs_exact.tsv')}")
    print(f"[WROTE] {os.path.join(args.outdir,'contigs_per_rank.tsv')}")
    print(f"[WROTE] per-rank diffs: {os.path.join(args.outdir,'profile_diffs_<rank>.tsv')}")
    print(f"[WROTE] debug: {os.path.join(args.outdir,'_debug_info.txt')}")

if __name__ == "__main__":
    main()
