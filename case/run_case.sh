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
SCENARIO="${CASE_SCENARIO:-cases}"
SUITE_NAME="${CASE_SUITE:-canonical}"
SUITE_PATH=""
PUBLISH_RESULTS=1
CUSTOM_OUT=0
RUN_STAMP=""
RUN_DIR=""
TABLES_DIR=""
FIGURES_DIR=""
declare -a CASE_SAMPLES=()

usage(){
  cat <<'USAGE'
Usage: run_case.sh [--manifest TSV] [--out DIR] [--threads N]
                   [--top-n K] [--sanity-metaphlan]
                   [--metaphlan-cmd PATH] [--metaphlan-opts "..."]
                   [--scenario NAME] [--suite NAME]
                   [--suite-path REL/PATH] [--no-publish]
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
    --scenario) SCENARIO="$2"; shift 2;;
    --suite) SUITE_NAME="$2"; shift 2;;
    --suite-path) SUITE_PATH="$2"; shift 2;;
    --no-publish) PUBLISH_RESULTS=0; shift;;
    -h|--help) usage;;
    *) usage;;
  esac
done

MANIFEST="$(resolve_path "${MANIFEST}")"

if [[ "${OUT_ROOT}" != "${CASE_ROOT}/out" ]]; then
  CUSTOM_OUT=1
  PUBLISH_RESULTS=0
fi

if [[ ${CUSTOM_OUT} -eq 1 || "${PUBLISH_RESULTS}" -eq 0 ]]; then
  OUT_ROOT="$(resolve_path "${OUT_ROOT}")"
  ensure_dir "${OUT_ROOT}"
else
  RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
  if [[ -n "${SUITE_PATH}" ]]; then
    if [[ "${SUITE_PATH}" = /* ]]; then
      run_base="${SUITE_PATH}"
    else
      run_base="${HYMET_ROOT}/results/${SUITE_PATH}"
    fi
  else
    run_base="${HYMET_ROOT}/results/${SCENARIO}/${SUITE_NAME}"
  fi
  RUN_DIR="${run_base}/run_${RUN_STAMP}"
  OUT_ROOT="$(resolve_path "${RUN_DIR}/raw")"
  TABLES_DIR="$(resolve_path "${RUN_DIR}/tables")"
  FIGURES_DIR="$(resolve_path "${RUN_DIR}/figures")"
  ensure_dir "${OUT_ROOT}"
  ensure_dir "${TABLES_DIR}"
  ensure_dir "${FIGURES_DIR}"
fi

RUNTIME_TSV="${OUT_ROOT}/runtime_memory.tsv"

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
SUMMARY_METAPHLAN="${OUT_ROOT}/metaphlan_metrics.tsv"
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

  CASE_SAMPLES+=("${sample_id}")
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
        if [[ ! -s "${SUMMARY_METAPHLAN}" ]]; then
          ensure_dir "$(dirname "${SUMMARY_METAPHLAN}")"
          echo -e "sample\tSymmetric_KL_Divergence\tSpearman_Rank" > "${SUMMARY_METAPHLAN}"
        fi
        tail -n +2 "${metrics}" >> "${SUMMARY_METAPHLAN}"
      fi
    else
      log "WARNING: MetaPhlAn profile missing for ${sample_id}; comparison skipped."
    fi
  fi
done < "${MANIFEST}"

if [[ "${PUBLISH_RESULTS}" -eq 1 && -n "${RUN_DIR}" ]]; then
  copy_if_exists(){
    local src="$1"
    local dst="$2"
    if [[ -f "${src}" ]]; then
      install -m 0644 -D "${src}" "${dst}"
    fi
  }

  copy_if_exists "${RUNTIME_TSV}" "${TABLES_DIR}/runtime_memory.tsv"
  copy_if_exists "${TOP_SUMMARY}" "${TABLES_DIR}/top_taxa_summary.tsv"
  copy_if_exists "${SUMMARY_METAPHLAN}" "${TABLES_DIR}/metaphlan_metrics.tsv"

  python3 "${CASE_ROOT}/plot_case.py" \
    --case-root "${OUT_ROOT}" \
    --figures-dir "${FIGURES_DIR}" \
    --max-taxa "${TOP_N}"

  CASE_SAMPLE_LIST=""
  if [[ ${#CASE_SAMPLES[@]} -gt 0 ]]; then
    CASE_SAMPLE_LIST="$(printf "%s\n" "${CASE_SAMPLES[@]}")"
  fi
  export CASE_SAMPLE_LIST
  export RUN_DIR SCENARIO SUITE_NAME RUN_STAMP MANIFEST THREADS TOP_N SANITY_METAPHLAN HYMET_ROOT
  python3 - <<'PY'
import json, os, pathlib, subprocess
run_dir = pathlib.Path(os.environ["RUN_DIR"])
manifest = pathlib.Path(os.environ["MANIFEST"]).resolve()
samples = [s for s in os.environ.get("CASE_SAMPLE_LIST", "").splitlines() if s]
meta = {
    "scenario": os.environ["SCENARIO"],
    "suite": os.environ["SUITE_NAME"],
    "run_id": run_dir.name,
    "timestamp": os.environ.get("RUN_STAMP") or "manual",
    "manifest": str(manifest),
    "threads": int(os.environ.get("THREADS") or 0),
    "top_n": int(os.environ.get("TOP_N") or 0),
    "samples": samples,
    "sanity_metaphlan": bool(int(os.environ.get("SANITY_METAPHLAN") or "0")),
    "source": "case/run_case.sh",
    "git_commit": subprocess.run(
        ["git", "-C", os.environ["HYMET_ROOT"], "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip() or "unknown",
}
(run_dir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
PY
  log "[case] Published case-study artefacts to ${RUN_DIR}"
else
  log "[case] Completed case-study run. Outputs in ${OUT_ROOT}"
fi
