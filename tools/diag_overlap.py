#!/usr/bin/env python3
import sys, gzip, csv

if len(sys.argv)!=3:
    sys.exit("usage: diag_overlap.py input.paf id_to_taxid.tsv")
paf,mapf=sys.argv[1],sys.argv[2]

idmap={}
with open(mapf, encoding='utf-8', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if not row: continue
        idmap[row[0]]=row[1]
keys=set(idmap.keys())

def opener(p):
    return gzip.open(p,'rt',encoding='utf-8',errors='ignore') if p.endswith('.gz') else open(p,'r',encoding='utf-8',errors='ignore')

targets=set()
with opener(paf) as f:
    for ln in f:
        if not ln or ln[0]=='#': continue
        p=ln.rstrip('\n').split('\t')
        if len(p)>=6:
            targets.add(p[5])

direct=sum(1 for t in targets if t in keys)
versionless=sum(1 for t in targets if (t.split('.',1)[0] in keys) and (t not in keys))

print(f"[diag] id_map_keys={len(keys):,} unique_paf_targets={len(targets):,} direct_overlap={direct} versionless_overlap={versionless}")
