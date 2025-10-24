#!/usr/bin/env bash
# Execute HYMET on real-world case-study samples defined in case/manifest.tsv.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${CASE_ROOT}/manifest.tsv"
OUT_ROOT="${CASE_ROOT}/out"
THREADS="${THREADS:-8}"
TOP_N="${TOP_N:-10}"
SANITY_METAPHLAN=0
METAPHLAN_CMD="${METAPHLAN_CMD:-metaphlan}"
METAPHLAN_OPTS="${METAPHLAN_OPTS:-}"
MEASURE="${CASE_ROOT}/lib/measure.sh"
RUNTIME_TSV="${OUT_ROOT}/runtime_memory.tsv"

usage(){
  cat <<'USAGE'
Usage: run_case.sh [--manifest TSV] [--out DIR] [--threads N]
                   [--top-n K] [--sanity-metaphlan]
                   [--metaphlan-cmd PATH] [--metaphlan-opts "..."]
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2;;
    --out) OUT_ROOT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --top-n) TOP_N="$2"; shift 2;;
    --sanity-metaphlan) SANITY_METAPHLAN=1; shift;;
    --metaphlan-cmd) METAPHLAN_CMD="$2"; shift 2;;
    --metaphlan-opts) METAPHLAN_OPTS="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

MANIFEST="$(resolve_path "${MANIFEST}")"
OUT_ROOT="$(resolve_path "${OUT_ROOT}")"
ensure_dir "${OUT_ROOT}"

[[ -s "${MANIFEST}" ]] || die "Manifest not found: ${MANIFEST}"

append_summary_header(){
  local path="$1"
  if [[ ! -s "${path}" ]]; then
    ensure_dir "$(dirname "${path}")"
    cat <<'EOF' >"${path}"
sample	metric	value
EOF
  fi
}

TOP_SUMMARY="${OUT_ROOT}/top_taxa_summary.tsv"
append_summary_header "${TOP_SUMMARY}"

while IFS= read -r line || [[ -n "${line}" ]]; do
  [[ -z "${line}" || "${line}" == \#* || "${line}" == sample_id* ]] && continue
  IFS=$'\x1f' read -r sample_id contigs truth_contigs truth_profile expected_taxa citation <<<"$(manifest_split_line "${line}")"
  log "[case] parsed row: sample='${sample_id}' contigs='${contigs}'"

  if [[ -z "${sample_id}" ]]; then
    log "Skipping manifest line with empty sample_id."
    continue
  fi

  contigs_abs="$(resolve_path "${contigs}")"
  if [[ ! -s "${contigs_abs}" ]]; then
    log "WARNING: sample ${sample_id} missing contigs at ${contigs_abs}; skipping."
    continue
  fi

  sample_out="${OUT_ROOT}/${sample_id}"
  hymet_out="${sample_out}/hymet"
  ensure_dir "${hymet_out}"

  log "[case] Running HYMET for ${sample_id}"
  THREADS="${THREADS}" "${MEASURE}" \
    --tool hymet \
    --sample "${sample_id}" \
    --stage run \
    --out "${RUNTIME_TSV}" \
    -- "${BENCH_ROOT}/run_hymet.sh" \
         --sample "${sample_id}" \
         --contigs "${contigs_abs}" \
         --out "${hymet_out}" \
         --threads "${THREADS}"

  profile="${hymet_out}/profile.cami.tsv"
  classified="${hymet_out}/classified_sequences.tsv"
  [[ -s "${profile}" ]] || die "HYMET profile missing for ${sample_id}: ${profile}"

  top_table="${sample_out}/top_taxa.tsv"
  python3 - "${profile}" "${top_table}" "${TOP_N}" <<'PY'
import csv, sys
profile_path, out_path, top_n = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = []
with open(profile_path) as fh:
    for line in fh:
        line = line.strip()
        if not line or line[0] in "#@":
            continue
        taxid, rank, taxpath, taxpathsn, pct = line.split("\t")
        try:
            pct_val = float(pct)
        except ValueError:
            continue
        rows.append((pct_val, rank, taxid, taxpathsn, taxpath))
rows.sort(reverse=True)
with open(out_path, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Rank", "TaxID", "TaxPathSN", "TaxPath", "Percentage"])
    for pct, rank, taxid, taxpathsn, taxpath in rows[:top_n]:
        writer.writerow([rank, taxid, taxpathsn, taxpath, f"{pct:.6f}"])
PY

  python3 - "${TOP_SUMMARY}" "${sample_id}" "${top_table}" <<'PY'
import csv, sys, pathlib
summary_path, sample_id, table_path = sys.argv[1:4]
rows = []
with open(table_path) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        rows.append(row)
with open(summary_path, "a", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    for entry in rows:
        writer.writerow([sample_id, f"top_{entry['Rank']}", f"{entry['TaxPathSN']} ({entry['Percentage']})"])
PY

  git_commit="$(git -C "${HYMET_ROOT}" rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
  git_dirty="$(git -C "${HYMET_ROOT}" status --porcelain 2>/dev/null || echo '')"
  if [[ -n "${git_dirty}" && "${git_commit}" != "unknown" ]]; then git_commit="${git_commit}-dirty"; fi

  cat > "${sample_out}/metadata.json" <<EOF
{
  "sample_id": "${sample_id}",
  "hymet_commit": "${git_commit}",
  "threads": ${THREADS},
  "contigs": "${contigs_abs}",
  "truth_contigs": "${truth_contigs}",
  "truth_profile": "${truth_profile}",
  "expected_taxa": "${expected_taxa}",
  "citation": "${citation}",
  "hymet_profile": "${profile}",
  "hymet_classified": "${classified}"
}
EOF

  normalize_metadata_json "${sample_out}/metadata.json" "${sample_out}"

  if [[ "${SANITY_METAPHLAN}" -eq 1 ]]; then
    metaphlan_out="${sample_out}/metaphlan"
    ensure_dir "${metaphlan_out}"
    mp_profile="${metaphlan_out}/profile.tsv"
    log "[case] Running MetaPhlAn sanity check for ${sample_id}"
    "${MEASURE}" \
      --tool metaphlan4 \
      --sample "${sample_id}" \
      --stage run \
      --out "${RUNTIME_TSV}" \
      -- "${METAPHLAN_CMD}" "${contigs_abs}" \
            --input_type fasta \
            --nproc "${THREADS}" \
            ${METAPHLAN_OPTS} \
            -o "${mp_profile}"

    if [[ -s "${mp_profile}" ]]; then
      comparison="${metaphlan_out}/comparison.tsv"
      metrics="${metaphlan_out}/metrics.tsv"
      python3 - "${profile}" "${mp_profile}" "${comparison}" "${metrics}" "${sample_id}" <<'PY'
import csv, sys, math
from collections import OrderedDict

hymet_path, mp_path, comp_path, metrics_path, sample_id = sys.argv[1:6]
EPS = 1e-8

def load_profile(path):
    prof = OrderedDict()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "#@":
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            taxpathsn = parts[3]
            pct = float(parts[4])
            prof[taxpathsn] = pct
    return prof

def build_distribution(hymet, meta):
    taxa = sorted(set(hymet) | set(meta))
    hp = []
    mp = []
    for tax in taxa:
        hp.append(max(hymet.get(tax, 0.0), 0.0))
        mp.append(max(meta.get(tax, 0.0), 0.0))
    total_h = sum(hp)
    total_m = sum(mp)
    if total_h <= 0:
        hp = [1.0 / len(hp) for _ in hp]
    else:
        hp = [v / total_h for v in hp]
    if total_m <= 0:
        mp = [1.0 / len(mp) for _ in mp]
    else:
        mp = [v / total_m for v in mp]
    return taxa, hp, mp

def symmetric_kl(p, q):
    kl_pq = 0.0
    kl_qp = 0.0
    for pi, qi in zip(p, q):
        pi = max(pi, EPS)
        qi = max(qi, EPS)
        kl_pq += pi * math.log(pi / qi)
        kl_qp += qi * math.log(qi / pi)
    return 0.5 * (kl_pq + kl_qp)

def ranks(values):
    n = len(values)
    order = sorted(enumerate(values), key=lambda x: x[1], reverse=True)
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and order[j + 1][1] == order[i][1]:
            j += 1
        rank_value = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            idx = order[k][0]
            ranks[idx] = rank_value
        i = j + 1
    return ranks

def spearman(p, q):
    if len(p) < 2:
        return float("nan")
    rp = ranks(p)
    rq = ranks(q)
    mean_p = sum(rp) / len(rp)
    mean_q = sum(rq) / len(rq)
    num = sum((a - mean_p) * (b - mean_q) for a, b in zip(rp, rq))
    den_p = math.sqrt(sum((a - mean_p) ** 2 for a in rp))
    den_q = math.sqrt(sum((b - mean_q) ** 2 for b in rq))
    if den_p == 0 or den_q == 0:
        return float("nan")
    return num / (den_p * den_q)

hymet = load_profile(hymet_path)
meta = load_profile(mp_path)
taxa, hp, mp = build_distribution(hymet, meta)

with open(comp_path, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["TaxPathSN", "HYMET_Percent", "MetaPhlAn_Percent", "Absolute_Difference"])
    for tax, h, m in zip(taxa, hp, mp):
        writer.writerow([tax, f"{h*100:.6f}", f"{m*100:.6f}", f"{abs(h-m)*100:.6f}"])

sym_kl = symmetric_kl(hp, mp)
spearman_corr = spearman(hp, mp)
with open(metrics_path, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Sample", "Symmetric_KL_Divergence", "Spearman_Rank"])
    writer.writerow([sample_id, f"{sym_kl:.6f}", f"{spearman_corr:.6f}"])
PY
      if [[ -s "${metrics}" ]]; then
        summary_metrics="${OUT_ROOT}/metaphlan_metrics.tsv"
        if [[ ! -s "${summary_metrics}" ]]; then
          ensure_dir "$(dirname "${summary_metrics}")"
          echo -e "sample\tSymmetric_KL_Divergence\tSpearman_Rank" > "${summary_metrics}"
        fi
        tail -n +2 "${metrics}" >> "${summary_metrics}"
      fi
    else
      log "WARNING: MetaPhlAn profile missing for ${sample_id}; comparison skipped."
    fi
  fi
done < "${MANIFEST}"

log "[case] Completed case-study run. Outputs in ${OUT_ROOT}"
