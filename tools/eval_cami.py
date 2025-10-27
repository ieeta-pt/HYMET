#!/usr/bin/env python3
import os, sys, csv, argparse, collections, subprocess, pathlib, re, hashlib, shutil, gzip
csv.field_size_limit(1024*1024*1024)

# -------- Defaults for your layout --------
DEF_PRED_PROFILE = "/data/hymet_out/sample_0/hymet.sample_0.cami.tsv"
DEF_TRUTH_PROFILE= "/data/cami/sample_0/taxonomic_profile_0.txt"
DEF_PRED_CONTIGS = "/data/hymet_out/sample_0/work/classified_sequences.tsv"
DEF_TRUTH_CONTIGS_NEW = "/data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/gsa_mapping_new.tsv"
DEF_TRUTH_CONTIGS_OLD = "/data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/gsa_mapping.tsv"
DEF_PRED_FASTA = "/data/cami/sample_0.fna"
DEF_TRUTH_FASTA= "/data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/anonymous_gsa.fasta"
DEF_TAXDB = "/data/HYMET/taxonomy_files"
DEF_TAXMAP = "/data/HYMET/data/detailed_taxonomy.tsv"
DEF_PAF = "/data/hymet_out/sample_0/work/resultados.paf"
DEF_OUTDIR = "/data/hymet_out/sample_0/eval"
RANKS = ["superkingdom","phylum","class","order","family","genus","species"]
RANKC = ["k","p","c","o","f","g","s"]
LINEAGE_ALIAS = {
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
    "ss": "species",
    "strain": "species"
}

# ----------------- utils -----------------
def ensure_dir(p): pathlib.Path(p).mkdir(parents=True, exist_ok=True)
def is_num(s):
    s = (s or "").strip()
    if not s:
        return False
    if s.isdigit():
        return True
    if re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", s):
        return True
    return False

def normalize_taxid(val: str) -> str:
    if not val:
        return ""
    match = re.search(r"[0-9]+", val)
    return match.group(0) if match else ""
def open_any(path):
    if path.endswith(".gz"): return gzip.open(path, "rt")
    return open(path, "r")
def chunked(seq, n):
    for i in range(0, len(seq), n): yield seq[i:i+n]


def parse_lineage_string(lineage_raw: str) -> dict:
    names_by_rank = {}
    if not lineage_raw:
        return names_by_rank
    parts = [seg.strip() for seg in lineage_raw.split(";") if seg.strip()]
    for part in parts:
        if ":" not in part:
            continue
        rk, nm = part.split(":", 1)
        rk = LINEAGE_ALIAS.get(rk.strip().lower())
        nm = nm.strip()
        if rk and nm:
            names_by_rank[rk] = nm
    return names_by_rank

def fasta_lengths(paths):
    lens = {}
    for path in paths:
        if not path or not os.path.isfile(path): continue
        name=None; L=0
        with open(path) as f:
            for ln in f:
                if ln.startswith(">"):
                    if name is not None: lens.setdefault(name, L)
                    name = ln[1:].strip().split()[0]; L = 0
                else: L += len(ln.strip())
        if name is not None: lens.setdefault(name, L)
    return lens

def fasta_hashes(path):
    hmap = {}
    if not path or not os.path.isfile(path): return hmap
    name=None; md=None
    with open(path) as f:
        for ln in f:
            if ln.startswith(">"):
                if name is not None: hmap[name]=md.hexdigest()
                name = ln[1:].strip().split()[0]; md = hashlib.md5()
            else:
                s=ln.strip().encode()
                if s: md.update(s)
        if name is not None: hmap[name]=md.hexdigest()
    return hmap

def taxonkit_name2taxid(names, taxdb):
    out={}
    if not names: return out
    env=os.environ.copy(); env["TAXONKIT_DB"]=taxdb
    for ch in chunked(list(names), 50000):
        p=subprocess.run(["taxonkit","name2taxid","--data-dir",taxdb,"--show-rank"],
                         input="\n".join(ch)+"\n", text=True, capture_output=True, check=True, env=env)
        for line in p.stdout.strip().splitlines():
            ps=line.split("\t")
            if len(ps)>=2 and is_num(ps[1]): out[ps[0]]=ps[1]
    return out

def taxonkit_taxpath(taxids, taxdb):
    out={}
    if not taxids: return out
    env=os.environ.copy(); env["TAXONKIT_DB"]=taxdb
    for ch in chunked(list(taxids), 50000):
        p=subprocess.run(["taxonkit","reformat","--data-dir",taxdb,"-I","1","-f","{d}|{p}|{c}|{o}|{f}|{g}|{s}","-t","-T"],
                         input="\n".join(ch)+"\n", text=True, capture_output=True, check=True, env=env)
        for line in p.stdout.strip().splitlines():
            ps=line.split("\t")
            if len(ps)>=3:
                def _norm(path):
                    parts=(path or "").split("|") if path else []
                    if len(parts)<len(RANKS): parts.extend(["NA"]*(len(RANKS)-len(parts)))
                    return "|".join(parts[:len(RANKS)])
                names=_norm(ps[1])
                ids=_norm(ps[2])
                out[ps[0]]=(ids,names)
    return out

# -------------- taxonomy id map ---------------
GCFA_RE = re.compile(r'GC[AF]_\d+(?:\.\d+)?(?:_PRJ[A-Z]+\d+)?')
ACC_RE  = re.compile(r'(NC_\d+\.\d+|NZ_[A-Z]{2}\d+\.\d+|NZ_[A-Z]{5}\d+\.\d+|CP\d+\.\d+|CM\d+\.\d+|[A-Z]{2}_\d+\.\d+)')

def _add_tok(m, tok, taxid):
    if not tok: return
    tok=tok.strip()
    if not tok: return
    m.setdefault(tok, taxid)
    if "." in tok: m.setdefault(tok.split(".",1)[0], taxid)

def load_id_map(taxmap_path):
    id2tax={}
    if not os.path.isfile(taxmap_path): return id2tax
    with open(taxmap_path, newline='') as f:
        r=csv.DictReader(f, delimiter="\t")
        for row in r:
            tax=normalize_taxid(row.get('TaxID') or '')
            if not tax: continue
            for key in ('GCF','GCA'):
                v=(row.get(key) or '').strip()
                if v:
                    _add_tok(id2tax, v, tax)
            ids=(row.get('Identifiers') or '')
            for tok in re.split(r'[;|,\s]+', ids):
                _add_tok(id2tax, tok, tax)
            # embedded patterns anywhere
            for v in row.values():
                if not v: continue
                for g in GCFA_RE.findall(v): _add_tok(id2tax, g, tax)
                for a in ACC_RE.findall(v):  _add_tok(id2tax, a, tax)
    return id2tax

# -------------- profiles (CAMI) ---------------
def _parse_cami_like(lines, taxdb):
    prof = {r: collections.Counter() for r in RANKS}
    ok=False
    for ln in lines:
        if not ln.strip() or ln[0] in "#@": continue
        ps=ln.rstrip("\n").split("\t")
        if len(ps)>=5 and is_num(ps[0]):
            rk=ps[1].strip().lower(); 
            if rk in prof:
                try: prof[rk][ps[0]]+=float(ps[4]); ok=True
                except: pass
            continue
        break
    if ok: return prof
    rdr=csv.reader([ln for ln in lines if ln.strip() and ln[0] not in "#@"], delimiter="\t")
    try: hdr=next(rdr)
    except StopIteration: return prof
    h=[c.strip().lower() for c in hdr]
    def idx(*names):
        for n in names:
            if n in h: return h.index(n)
        return -1
    i_taxid = idx("taxid","taxon_id","ncbi_taxid","ncbi_tax_id")
    i_rank  = idx("rank")
    i_perc  = idx("percentage","abundance","rel_abundance","fraction_total_reads")
    i_taxpath=idx("taxpath")
    i_taxpathsn=idx("taxpathsn","taxpath_sn","taxpath_names","lineage")
    rows=list(rdr)
    if i_taxid>=0 and i_rank>=0 and i_perc>=0:
        mul = 100.0 if "abundance" in h[i_perc] or "fraction" in h[i_perc] else 1.0
        for ps in rows:
            try:
                tid=normalize_taxid(ps[i_taxid]); rk=ps[i_rank].strip().lower(); val=float(ps[i_perc])*mul
                if rk in prof and tid: prof[rk][tid]+=val
            except: pass
        return prof
    if i_rank>=0 and (i_taxpath>=0 or i_taxpathsn>=0) and i_perc>=0:
        rk_to_idx=dict(zip(RANKS, range(len(RANKS))))
        mul = 100.0 if "abundance" in h[i_perc] or "fraction" in h[i_perc] else 1.0
        if i_taxpath>=0:
            for ps in rows:
                try:
                    rk=ps[i_rank].strip().lower()
                    ids=[x for x in ps[i_taxpath].strip().split("|") if x and x!="NA"]
                    idx_r=rk_to_idx.get(rk, -1)
                    if idx_r>=0 and idx_r<len(ids): prof[rk][ids[idx_r]]+=float(ps[i_perc])*mul
                except: pass
            return prof
        names=set(); keep=[]
        for ps in rows:
            try:
                rk=ps[i_rank].strip().lower()
                pathn=[p.strip() for p in ps[i_taxpathsn].split("|")]
                idx_r=rk_to_idx.get(rk, -1)
                if idx_r>=0 and idx_r<len(pathn) and pathn[idx_r]: names.add(pathn[idx_r])
                keep.append(ps)
            except: pass
        m=taxonkit_name2taxid(names, taxdb)
        for ps in keep:
            try:
                rk=ps[i_rank].strip().lower(); pathn=[p.strip() for p in ps[i_taxpathsn].split("|")]
                idx_r=rk_to_idx.get(rk, -1)
                if idx_r>=0 and idx_r<len(pathn):
                    tid=m.get(pathn[idx_r])
                    if tid: prof[rk][tid]+=float(ps[i_perc])*mul
            except: pass
    return prof

def load_profile_any(path, taxdb):
    if not os.path.isfile(path):
        return {r: collections.Counter() for r in RANKS}
    with open(path) as f: lines=f.readlines()
    return _parse_cami_like(lines, taxdb)

# -------------- contig truth mapping --------------
def load_gt_contigs(gt_file, taxdb):
    out={}
    if not os.path.isfile(gt_file): return out
    with open(gt_file) as fh:
        first_line = fh.readline()
    if "\t" in first_line:
        sep = "\t"
        header = [c.strip() for c in first_line.strip().split("\t")]
        with open_any(gt_file) as f:
            rdr = csv.reader(f, delimiter="\t")
            hdr = next(rdr)
            rows = [row for row in rdr]
    elif "," in first_line:
        sep = ","
        header = [c.strip() for c in first_line.strip().split(",")]
        with open_any(gt_file) as f:
            rdr = csv.reader(f, delimiter=",")
            hdr = next(rdr)
            rows = [row for row in rdr]
    else:
        sep = None
        header = [c.strip() for c in first_line.strip().split()]
        rows = []
        with open_any(gt_file) as f:
            next(f)
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rows.append(line.split())
        hdr = header
    h=[c.strip().lstrip('#').lower() for c in hdr]
    contig_keys = [k for k in h if any(x in k for x in ("contig","sequence","scaffold"))]
    taxid_keys  = [k for k in h if ("tax" in k and "path" not in k)] + [k for k in h if k in ("ncbi_taxid","ncbi_tax_id","taxid","tax_id","species_taxid","genome_taxid")]
    ci = h.index(contig_keys[0]) if contig_keys else 0
    ti = h.index(taxid_keys[0])  if taxid_keys  else -1
    if ti>=0:
        for ps in rows:
            if len(ps)<=max(ci,ti): continue
            raw=(ps[ti] or "").strip()
            if not raw: continue
            if not is_num(raw):
                raw = normalize_taxid(raw)
            if raw:
                out[ps[ci]]=normalize_taxid(raw)
    else:
        try:
            tpi = h.index("taxpath")
            for ps in rows:
                ids=[x for x in ps[tpi].split("|") if x and x!="NA"]
                if ids:
                    tid = normalize_taxid(ids[-1])
                    if tid:
                        out[ps[ci]] = tid
        except ValueError:
            for ps in rows:
                for x in ps[1:]:
                    if is_num(x):
                        out[ps[0]]=normalize_taxid(x)
                        break
    return out

# ---------- build profiles from contigs ----------
def profiles_from_contig_maps(contig2tid, lengths, taxdb):
    prof = {r: collections.Counter() for r in RANKS}
    if not contig2tid: return prof
    normalized = {}
    for cont, tid in contig2tid.items():
        ntid = normalize_taxid(tid)
        if ntid:
            normalized[cont] = ntid
    if not normalized:
        return prof
    tids=set(normalized.values()); paths=taxonkit_taxpath(tids, taxdb)
    acc=collections.Counter()
    for cont,tid in normalized.items():
        w=lengths.get(cont,1); names_ids=paths.get(tid)
        if not names_ids: continue
        ids=names_ids[0].split("|")
        for i,_ in enumerate(RANKC):
            if i<len(ids) and ids[i]!="NA":
                prof[RANKS[i]][ids[i]]+=w; acc[RANKS[i]]+=w
    for r in RANKS:
        s=acc[r]
        if s>0:
            for k in list(prof[r].keys()): prof[r][k]=100.0*prof[r][k]/s
    return prof

# ---------------- minimap2 fallback ----------------
def have_minimap2(): return shutil.which("minimap2") is not None
def run_minimap2(pred_fa, truth_fa, paf_out, threads=8):
    cmd=["minimap2","-x","asm10","--secondary=no","-t",str(threads), truth_fa, pred_fa]
    with open(paf_out,"w") as w: subprocess.run(cmd, check=True, stdout=w)

def besthit_map_from_paf(paf_path, min_cov=0.95, min_id=0.95):
    out={}; best={}
    with open(paf_path) as f:
        for ln in f:
            if not ln.strip() or ln[0]=='#': continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<12: continue
            q,qlen,qstart,qend=p[0],int(p[1]),int(p[2]),int(p[3])
            t=p[5]; nmatch,alen=int(p[9]),int(p[10])
            cov=(qend-qstart)/qlen if qlen>0 else 0.0
            iden=nmatch/alen if alen>0 else 0.0
            if cov<min_cov or iden<min_id: continue
            score=nmatch
            cur=best.get(q)
            if cur is None or score>cur[0]: best[q]=(score,t,cov,iden)
    for q,(_,t,_,_) in best.items(): out[q]=t
    return out

# -------------- PAF (HYMET) first-hit --------------
def paf_firsthit_q2t(paf_path):
    q2t={}
    if not os.path.isfile(paf_path): return q2t
    with open(paf_path) as f:
        for ln in f:
            if not ln.strip() or ln[0]=='#': continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<6: continue
            q,t=p[0],p[5]
            if q not in q2t: q2t[q]=t
    return q2t

# ---------------- metrics ----------------
def l1_and_braycurtis(a: dict, b: dict):
    keys=set(a)|set(b)
    if not keys: return 0.0,0.0
    sum_abs=sum(abs(a.get(k,0.0)-b.get(k,0.0)) for k in keys)
    l1=0.5*sum_abs
    sump=sum(a.get(k,0.0) for k in keys); sumt=sum(b.get(k,0.0) for k in keys)
    shared=sum(min(a.get(k,0.0),b.get(k,0.0)) for k in keys)
    bc = 1.0 - (2.0*shared / (sump+sumt if (sump+sumt)>0 else 1.0))
    return l1, bc*100.0

def prf_presence(a: dict, b: dict, thr=0.1):
    A={k for k,v in a.items() if v>=thr}; B={k for k,v in b.items() if v>=thr}
    tp=len(A&B); fp=len(A-B); fn=len(B-A)
    prec=tp/(tp+fp) if (tp+fp)>0 else 0.0
    rec =tp/(tp+fn) if (tp+fn)>0 else 0.0
    f1 =2*prec*rec/(prec+rec) if (prec+rec)>0 else 0.0
    return prec*100.0, rec*100.0, f1*100.0, tp, fp, fn

# -------------- read HYMET predictions into TaxIDs --------------
def preds_taxid_from_classified(classified_tsv, taxdb, idmap, paf_path):
    """Map each contig to the most specific resolvable TaxID."""
    cont2tid: dict[str, str] = {}
    lineage_records: dict[str, dict] = {}
    fallback_info: dict[str, dict] = {}
    all_names: set[str] = set()

    if os.path.isfile(classified_tsv):
        with open(classified_tsv, encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f, delimiter="	")
            raw_fields = reader.fieldnames or []
            headers = [(h or "").strip().lower() for h in raw_fields]
            idx_query = headers.index("query") if "query" in headers else None
            idx_taxid = headers.index("taxid") if "taxid" in headers else None
            idx_target = next((headers.index(c) for c in ("target", "tname") if c in headers), None)
            idx_lineage = headers.index("lineage") if "lineage" in headers else None
            key_query = raw_fields[idx_query] if idx_query is not None else None
            key_taxid = raw_fields[idx_taxid] if idx_taxid is not None else None
            key_target = raw_fields[idx_target] if idx_target is not None else None
            key_lineage = raw_fields[idx_lineage] if idx_lineage is not None else None

            for row in reader:
                q = row.get(key_query) if key_query else (row.get("Query") or row.get("qname") or row.get("q"))
                if not q:
                    continue

                lineage = parse_lineage_string(row.get(key_lineage, "") if key_lineage else row.get("Lineage", ""))
                if lineage:
                    lineage_records[q] = lineage
                    for nm in lineage.values():
                        if nm:
                            all_names.add(nm)

                fallback_info[q] = {
                    "taxid": row.get(key_taxid) if key_taxid else row.get("TaxID"),
                    "target": row.get(key_target) if key_target else (row.get("Target") or row.get("tname"))
                }

    name_map = taxonkit_name2taxid(all_names, taxdb) if all_names else {}

    for q, lineage in lineage_records.items():
        for rank in reversed(RANKS):
            nm = lineage.get(rank)
            if not nm:
                continue
            entry = name_map.get(nm)
            if not entry:
                continue
            tid = entry[0] if isinstance(entry, tuple) else entry
            ntid = normalize_taxid(tid)
            if ntid:
                cont2tid[q] = ntid
                break

    for q, info in fallback_info.items():
        if q in cont2tid:
            continue
        raw_tid = info.get("taxid") or ""
        ntid = normalize_taxid(raw_tid)
        if ntid:
            cont2tid[q] = ntid

    for q, info in fallback_info.items():
        if q in cont2tid:
            continue
        target = (info.get("target") or "").strip()
        if not target:
            continue
        base = target.split("|", 1)[0]
        candidates = [target, base]
        if "." in base:
            candidates.append(base.split(".", 1)[0])
        for cand in candidates:
            if cand in idmap:
                ntid = normalize_taxid(idmap[cand])
                if ntid:
                    cont2tid[q] = ntid
                    break

    if paf_path and os.path.isfile(paf_path):
        paf_hits = paf_firsthit_q2t(paf_path)
        for q, target in paf_hits.items():
            if q in cont2tid:
                continue
            base = target.split("|", 1)[0]
            candidates = [target, base]
            if "." in base:
                candidates.append(base.split(".", 1)[0])
            for cand in candidates:
                if cand in idmap:
                    ntid = normalize_taxid(idmap[cand])
                    if ntid:
                        cont2tid[q] = ntid
                        break

    return cont2tid

# -------------- contig-level eval --------------
def eval_contigs(pred_file, gt_files, taxdb, outdir, pred_fasta=None, gt_fasta=None, threads=8,
                 taxmap_path=DEF_TAXMAP, paf_path=DEF_PAF):
    idmap=load_id_map(taxmap_path)
    pred_tid = preds_taxid_from_classified(pred_file, taxdb, idmap, paf_path)

    gt_map={}
    for gtf in gt_files:
        m=load_gt_contigs(gtf, taxdb) if gtf else {}
        gt_map.update(m)

    print(f"[DEBUG] loaded pred contigs with TaxID: {len(pred_tid)}", file=sys.stderr)
    print(f"[DEBUG] loaded truth contigs with TaxID: {len(gt_map)}", file=sys.stderr)

    pairs=[]
    n_direct=0
    for q,tid in pred_tid.items():
        if q in gt_map:
            pairs.append((q, tid, gt_map[q])); n_direct+=1

    if not pairs and pred_fasta and gt_fasta and os.path.isfile(pred_fasta) and os.path.isfile(gt_fasta):
        pred_hash=fasta_hashes(pred_fasta); gt_hash=fasta_hashes(gt_fasta)
        inv_gt=collections.defaultdict(list)
        for gname,h in gt_hash.items(): inv_gt[h].append(gname)
        n_md5=0
        for q in list(pred_tid.keys()):
            h=pred_hash.get(q)
            if not h: continue
            for t in inv_gt.get(h, []):
                gtid=gt_map.get(t)
                if gtid:
                    pairs.append((q, pred_tid[q], gtid)); n_md5+=1
        print(f"[DEBUG] MD5‑paired contigs: {n_md5}", file=sys.stderr)

    if not pairs and have_minimap2() and pred_fasta and gt_fasta and os.path.isfile(pred_fasta) and os.path.isfile(gt_fasta):
        paf_tmp=os.path.join(outdir,"pred_vs_truth.paf")
        run_minimap2(pred_fasta, gt_fasta, paf_tmp, threads=threads)
        q2t=besthit_map_from_paf(paf_tmp, min_cov=0.95, min_id=0.95)
        n_map=0
        for q,t in q2t.items():
            pt=pred_tid.get(q); gtid=gt_map.get(t)
            if pt and gtid:
                pairs.append((q, pt, gtid)); n_map+=1
        print(f"[DEBUG] minimap‑paired contigs: {n_map}", file=sys.stderr)

    usable=len(pairs); exact=sum(1 for _,pt,gtid in pairs if pt==gtid)

    tids={pt for _,pt,_ in pairs} | {gtid for *_,gtid in pairs}
    tpaths=taxonkit_taxpath(tids, taxdb)

    per_rank={}
    for i,r in enumerate(RANKS):
        tot=0; ok=0
        for _,pt,gtid in pairs:
            pids=tpaths.get(pt,("",""))[0]; gids=tpaths.get(gtid,("",""))[0]
            if not pids or not gids: continue
            pvec=pids.split("|"); gvec=gids.split("|")
            if i>=len(pvec) or i>=len(gvec): continue
            pid, gid = pvec[i], gvec[i]
            if pid=="NA" or gid=="NA": continue
            tot+=1
            if pid==gid: ok+=1
        per_rank[r]={"n":tot,"acc":(100.0*ok/tot if tot else 0.0),"correct":ok}

    exact_path=os.path.join(outdir,"contigs_exact.tsv")
    perrank_path=os.path.join(outdir,"contigs_per_rank.tsv")

    if usable > 0:
        with open(exact_path,"w",newline="") as w:
            wr=csv.writer(w, delimiter="\t"); wr.writerow(["metric","value"])
            wr.writerow(["usable_pairs", usable]); wr.writerow(["exact_taxid_matches", exact])
            wr.writerow(["exact_taxid_accuracy_percent", 100.0*exact/usable if usable else 0.0])

        with open(perrank_path,"w",newline="") as w:
            wr=csv.writer(w, delimiter="\t"); wr.writerow(["rank","n","correct","accuracy_percent"])
            for r in RANKS:
                m=per_rank.get(r,{"n":0,"correct":0,"acc":0.0})
                wr.writerow([r, m["n"], m["correct"], f"{m['acc']:.4f}"])
    else:
        for p in (exact_path, perrank_path):
            if os.path.exists(p):
                os.remove(p)

    return {"usable_pairs": usable, "exact": exact, "per_rank": per_rank, "pred_n": len(pred_tid), "gt_n": len(gt_map)}

# ------------------- main -------------------
def main():
    ap=argparse.ArgumentParser(description="Evaluate HYMET vs CAMI ground truth (post-processing only; classifier unchanged).")
    ap.add_argument("--pred-profile", default=DEF_PRED_PROFILE)
    ap.add_argument("--truth-profile", default=DEF_TRUTH_PROFILE)
    ap.add_argument("--pred-contigs", default=DEF_PRED_CONTIGS)
    ap.add_argument("--truth-contigs", default="")
    ap.add_argument("--pred-fasta", default=DEF_PRED_FASTA)
    ap.add_argument("--truth-fasta", default=DEF_TRUTH_FASTA)
    ap.add_argument("--taxdb", default=DEF_TAXDB)
    ap.add_argument("--taxmap", default=DEF_TAXMAP)
    ap.add_argument("--paf", default=DEF_PAF)
    ap.add_argument("--outdir", default=DEF_OUTDIR)
    ap.add_argument("--presence-thresh", type=float, default=0.1)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("OMP_NUM_THREADS","8")))
    args=ap.parse_args()

    ensure_dir(args.outdir)
    gt_files=[args.truth_contigs] if args.truth_contigs else [DEF_TRUTH_CONTIGS_NEW, DEF_TRUTH_CONTIGS_OLD]

    pred_prof = load_profile_any(args.pred_profile, args.taxdb)
    truth_prof= load_profile_any(args.truth_profile, args.taxdb)
    need_pred_fb  = all(not pred_prof[r]  for r in RANKS)
    need_truth_fb = all(not truth_prof[r] for r in RANKS)

    lens={}
    if need_pred_fb or need_truth_fb:
        lens=fasta_lengths([args.pred_fasta, args.truth_fasta])

    if need_pred_fb:
        print("[INFO] Rebuilding predicted profile from per‑contig classifications.", file=sys.stderr)
        cont2tid=preds_taxid_from_classified(args.pred_contigs, args.taxdb, load_id_map(args.taxmap), args.paf)
        pred_prof=profiles_from_contig_maps(cont2tid, lens, args.taxdb)

    if need_truth_fb:
        print("[INFO] Rebuilding truth profile from contig mapping.", file=sys.stderr)
        gt_map={}
        for g in gt_files:
            gt_map.update(load_gt_contigs(g, args.taxdb))
        truth_prof=profiles_from_contig_maps(gt_map, lens, args.taxdb)

    def l1_bc_row(rank):
        a,b=pred_prof[rank], truth_prof[rank]
        keys=set(a)|set(b)
        sum_abs=sum(abs(a.get(k,0.0)-b.get(k,0.0)) for k in keys)
        l1=0.5*sum_abs
        sump=sum(a.get(k,0.0) for k in keys); sumt=sum(b.get(k,0.0) for k in keys)
        shared=sum(min(a.get(k,0.0),b.get(k,0.0)) for k in keys)
        bc = 1.0 - (2.0*shared / (sump+sumt if (sump+sumt)>0 else 1.0))
        pr,rc,f1,tp,fp,fn=prf_presence(a,b,args.presence_thresh)
        return l1, bc*100.0, pr, rc, f1, tp, fp, fn

    with open(os.path.join(args.outdir,"profile_summary.tsv"),"w",newline="") as w:
        wr=csv.writer(w, delimiter="\t")
        wr.writerow(["rank","L1_total_variation_pctpts","BrayCurtis_pct","Precision_%","Recall_%","F1_%","TP","FP","FN"])
        for r in RANKS:
            l1,bc,pr,rc,f1,tp,fp,fn=l1_bc_row(r)
            wr.writerow([r, f"{l1:.4f}", f"{bc:.4f}", f"{pr:.2f}", f"{rc:.2f}", f"{f1:.2f}", tp, fp, fn])

    print("# Profile-level metrics (per rank)")
    for r in RANKS:
        l1,bc,pr,rc,f1,tp,fp,fn=l1_bc_row(r)
        print(f"{r:14s}  L1={l1:.3f}  BC={bc:.3f}%  P/R/F1={pr:.1f}/{rc:.1f}/{f1:.1f}% (TP={tp}, FP={fp}, FN={fn})")

    print("\n# Contig-level accuracy")
    cres = eval_contigs(args.pred_contigs, gt_files, args.taxdb, args.outdir,
                        pred_fasta=args.pred_fasta, gt_fasta=args.truth_fasta, threads=args.threads,
                        taxmap_path=args.taxmap, paf_path=args.paf)
    usable=cres.get("usable_pairs",0); exact=cres.get("exact",0)
    print(f"Exact TaxID: {exact}/{usable} ({(100.0*exact/usable if usable else 0.0):.2f}%)")
    for r in RANKS:
        m=cres.get("per_rank",{}).get(r,{"n":0,"acc":0.0,"correct":0})
        print(f"{r:14s}  n={m['n']:<8d}  acc={m['acc']:.2f}%")

    with open(os.path.join(args.outdir,"_debug_info.txt"),"w") as w:
        w.write(f"pred_profile_path: {args.pred_profile}\n")
        w.write(f"truth_profile_path: {args.truth_profile}\n")
        w.write(f"pred_contigs_path: {args.pred_contigs}\n")
        w.write("truth_contigs_paths:\n  " + "\n  ".join([g for g in gt_files if g]) + "\n")
        w.write(f"pred_fasta: {args.pred_fasta}\n")
        w.write(f"truth_fasta: {args.truth_fasta}\n")
        w.write(f"taxdb: {args.taxdb}\n")
        w.write(f"taxmap: {args.taxmap}\n")
        w.write(f"paf: {args.paf}\n")

    print(f"\n[WROTE] {os.path.join(args.outdir,'profile_summary.tsv')}")
    print(f"[WROTE] {os.path.join(args.outdir,'contigs_exact.tsv')}")
    print(f"[WROTE] {os.path.join(args.outdir,'contigs_per_rank.tsv')}")
    print(f"[WROTE] debug: {os.path.join(args.outdir,'_debug_info.txt')}")

if __name__ == "__main__":
    main()
