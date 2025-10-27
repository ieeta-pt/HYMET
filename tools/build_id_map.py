#!/usr/bin/env python3
import sys, csv
from collections import OrderedDict

if len(sys.argv)!=3:
    sys.exit("usage: build_id_map.py detailed_taxonomy.tsv out_map.tsv")

# lift csv field size limit to handle huge "Identifiers" column
try:
    csv.field_size_limit(sys.maxsize)
except Exception:
    csv.field_size_limit(10**9)

inp, outp = sys.argv[1], sys.argv[2]
id2tax = OrderedDict()

def emit(k, tax):
    if not k: return
    id2tax.setdefault(k, tax)
    # versionless (e.g., NC_XXXXXX from NC_XXXXXX.1)
    if '.' in k:
        id2tax.setdefault(k.split('.', 1)[0], tax)

with open(inp, 'r', encoding='utf-8', errors='ignore', newline='') as f:
    first = f.readline()
    if not first:
        sys.exit("empty taxonomy file")
    hdr = first.rstrip('\n').split('\t')
    # Expect: GCF, TaxID, Identifiers
    try:
        i_gcf = hdr.index('GCF'); i_tax = hdr.index('TaxID'); i_ids = hdr.index('Identifiers')
    except ValueError:
        i_gcf, i_tax, i_ids = 0, 1, 2

    for line in f:
        if not line.strip(): continue
        row = line.rstrip('\n').split('\t')
        # guard
        if len(row) <= max(i_gcf, i_tax): continue
        gcf = row[i_gcf].strip()
        tax = row[i_tax].strip()
        if gcf:
            emit(gcf, tax)
        # identifiers list may be missing
        ids = row[i_ids].strip() if len(row) > i_ids else ''
        if ids:
            for tok in ids.split(';'):
                emit(tok.strip(), tax)

with open(outp, 'w', encoding='utf-8', newline='') as w:
    wr = csv.writer(w, delimiter='\t')
    for k, v in id2tax.items():
        wr.writerow([k, v])

print(f"wrote {len(id2tax):,} ids -> TaxID to {outp}")
