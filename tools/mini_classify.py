#!/usr/bin/env python3
import sys, csv, gzip
if len(sys.argv)!=4:
    sys.exit("usage: mini_classify.py input.paf id_to_taxid.tsv out.tsv")
paf,mapf,outp=sys.argv[1],sys.argv[2],sys.argv[3]

idmap={}
with open(mapf, encoding='utf-8', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if not row: continue
        idmap.setdefault(row[0], row[1])

def opener(p): 
    return gzip.open(p,'rt',encoding='utf-8',errors='ignore') if p.endswith('.gz') else open(p,'r',encoding='utf-8',errors='ignore')

seen=set(); n=0; tot=0
with opener(paf) as f, open(outp,'w',encoding='utf-8',newline='') as w:
    wr=csv.writer(w, delimiter='\t')
    wr.writerow(["qname","tname","taxid"])
    for ln in f:
        if not ln or ln[0]=='#': continue
        p=ln.rstrip('\n').split('\t')
        if len(p)<6: continue
        q, t = p[0], p[5]
        tot += 1
        if q in seen: 
            continue
        tax = idmap.get(t) or idmap.get(t.split('.',1)[0])
        if tax:
            wr.writerow([q,t,tax]); n+=1; seen.add(q)
print(f"[mini] Classified {n}/{tot} alignments (first-hit per query) -> {outp}")
