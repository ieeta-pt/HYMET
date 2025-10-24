# HYMET Reproducibility Playbook

This playbook expands the minimal quick-start instructions in the README into a detailed, GigaScience-aligned protocol. Follow it to recreate every experiment reported in the manuscript, regenerate tables/figures, and capture provenance so that independent reviewers can validate the submission.

> **Scope.** Commands are written for GNU/Linux (x86\_64) with Bash ≥ 4. They were tested on Ubuntu 22.04 LTS, but any recent distribution with the listed prerequisites should work.

---

## 1. Requirements

### 1.1 Hardware
- **CPU:** ≥8 physical cores recommended (16 for full CAMI reruns).
- **RAM:** 32 GB minimum; 64 GB recommended when running MetaPhlAn4 alongside HYMET.
- **Disk:** ~220 GB free space for references, intermediates, and outputs:
  - Mash sketches: 2.6 GB
  - Tool databases (CAMI harness): ≈120 GB
  - Benchmark outputs & logs: 15–25 GB
  - Case-study runs & ablation: ≈8 GB

### 1.2 Software
- `git`, `curl`, `wget`, `tar`, `unzip`
- `micromamba` (preferred) or `conda`
- Python ≥ 3.9
- Optional: `aria2c` (for faster downloads), Docker/Singularity (for container confirmation)

### 1.3 Network
- Access to NCBI/GTDB FTP mirrors, Zenodo (doi:10.5281/zenodo.17428354), MGnify, and CAMI archives.
- When working behind a proxy, configure `http_proxy`/`https_proxy` before running the scripts.

---

## 2. Provenance checklist (fill as you go)

| Item | Command / Location | Notes |
|------|-------------------|-------|
| Repository commit | `git rev-parse HEAD` | Record in lab notebook |
| HYMET environment file | `environment.lock.yml` | Keep copy with submission |
| Sketch checksums | `data/sketch_checksums.tsv` | Compare with Zenodo |
| Bench aggregate hash | `sha256sum bench/out/summary_per_tool_per_sample.tsv` | Match Zenodo |
| Case runtime log hash | `sha256sum case/out/runtime_memory.tsv` | Match Zenodo |

Use this table as a cover sheet for the reproducibility package.

---

## 3. Repository checkout

```bash
git clone https://github.com/ieeta-pt/HYMET.git
cd HYMET
```

Optional: pin to the exact revision referenced in `bench/out/*/metadata.json`:

```bash
git checkout <commit-hash>
```

Record the resulting hash.

---

## 4. Environment creation

### 4.1 Micromamba (preferred)

```bash
micromamba env create -f environment.lock.yml -n hymet-env
micromamba activate hymet-env
```

> If the lockfile cannot be solved (e.g., on unsupported architectures), fall back to `environment.yml` and regenerate the lock:
> ```bash
> micromamba env create -f environment.yml -n hymet-env
> micromamba env export -n hymet-env > environment.lock.yml
> ```

### 4.2 Verify toolchain

```bash
python -V
for tool in mash minimap2 taxonkit pytest; do
  command -v "$tool" >/dev/null || echo "Missing $tool"
done
```

Capture the Python version and any warnings in your notes.

---

## 5. Reference assets

### 5.1 Fetch Mash sketches (Zenodo)

```bash
tools/fetch_sketches.sh                # downloads sketch{1,2,3}.msh
tools/verify_sketches.sh               # validates SHA256 against Zenodo metadata
```

Expected files:
- `data/sketch{1,2,3}.msh`
- `data/sketch_checksums.tsv`

Store the console output confirming the checksum match.

### 5.2 NCBI taxonomy dump

```bash
mkdir -p taxonomy_files
curl -L ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip -o /tmp/taxdmp.zip
unzip -j /tmp/taxdmp.zip names.dmp nodes.dmp merged.dmp delnodes.dmp -d taxonomy_files
```

Check file modification dates to ensure you captured the latest release.

---

## 6. Dataset acquisition

### 6.1 CAMI subsets (lightweight panels)

```bash
cd bench
./fetch_cami.sh --manifest cami_manifest.tsv --dest data --threads 8
cd ..
```

Outputs appear under `bench/data/<sample>/`. Confirm that each folder contains:
- `contigs.fna`
- `truth/` tables (`contigs.tsv`, `profile.tsv`)

### 6.2 Case-study contamination panels

```bash
cd case
./fetch_case_data.sh --dest data case/manifest.tsv
cd ..
```

Verify `case/data/gut_case/*.fna` and `case/data/zymo_mc/*.fna` exist, plus truth sets under `case/truth/`.

---

## 7. CLI smoke tests

Run dry-run invocations to ensure the CLI resolves paths and dependencies.

```bash
bin/hymet version
bin/hymet run   --contigs bench/data/cami_i_lc/contigs.fna --out out/smoke --threads 1 --dry-run
bin/hymet bench --manifest bench/cami_manifest.tsv --tools hymet --max-samples 1 --dry-run
bin/hymet case  --manifest case/manifest_zymo.tsv --dry-run
```

Expected: commands print the planned shell invocation and exit with status 0.

---

## 8. Full CAMI benchmark reproduction

### 8.1 Resource preparation

```bash
export THREADS=16
export CACHE_ROOT="$(pwd)/bench/data/downloaded_genomes/cache_bench"
mkdir -p "$CACHE_ROOT"
```

> Each HYMET run hashes selected references and stores them in `$CACHE_ROOT/<sha1>/`. Keep this directory for re-runs.

### 8.2 Execute benchmark

The command below reproduces the tool panel reported in the manuscript. Adjust `--tools` to include or exclude specific competitors.

```bash
bin/hymet bench \
  --manifest bench/cami_manifest.tsv \
  --tools hymet,kraken2,centrifuge,ganon2,metaphlan4,sourmash_gather,megapath_nano,tama,basta,phabox,phyloflash,viwrap,squeezemeta \
  --resume
```

Key outputs:
- `bench/out/<sample>/<tool>/` – raw PAF files, CAMI profiles, logs, metadata.
- `bench/out/summary_per_tool_per_sample.tsv` – per-rank precision/recall/F1.
- `bench/out/runtime_memory.tsv` – `/usr/bin/time -v` metrics.
- `bench/out/fig_*.png` – aggregated plots (copied to `results/bench/`).

### 8.3 Partial reruns

To refresh a single tool on two samples:

```bash
THREADS=8 bin/hymet bench --tools metaphlan4 --max-samples 2 --resume
```

The harness skips tools/samples with existing success markers unless `FORCE_DOWNLOAD=1` or `--no-build` is specified.

---

## 9. Case studies & ablation

### 9.1 Standard case run

```bash
export CACHE_ROOT="$(pwd)/case/data/downloaded_genomes/cache_case"
export THREADS=8
bin/hymet case --manifest case/manifest.tsv --sanity-metaphlan
```

Outputs:
- `case/out/<sample>/hymet/` – CAMI profile, classified contigs, metadata.
- `case/out/<sample>/metaphlan/` – optional sanity comparison.
- `case/out/top_taxa_summary.tsv` – combined view.
- `case/out/runtime_memory.tsv` – resource usage.

### 9.2 Zymo curated-panel ablation

Steps:

1. Preload the curated reference into the cache (optional but recommended to limit downloads):
   ```bash
   python case/tools/preload_cache_from_fasta.py \
     --cache-dir "${CACHE_ROOT}/<sha1-from-log>" \
     --fasta case/truth/zymo_refs/zymo_refs.fna.gz \
     --seqmap case/truth/zymo_refs/seqid2taxid.tsv \
     --taxid-prefix ZymoTax
   ```
   Remove `reference.mmi` in that directory to force a fresh index.

2. Run ablation:
   ```bash
   THREADS=4 bin/hymet ablation \
     --sample zymo_mc \
     --taxa 1423,562,28901,1639,1280,1351,287,1613,4932,5207 \
     --levels 0,0.25,0.5,0.75,1.0 \
     --fasta case/truth/zymo_refs/zymo_refs.fna.gz \
     --seqmap case/truth/zymo_refs/seqid2taxid.tsv
   ```

Outputs:
- `case/ablation/<sample>/level_*/hymet/` – per-level predictions + evaluation.
- `case/ablation_summary.tsv`, `case/ablation_eval_summary.tsv`.
- `case/ablation/figures/*.png` – fallback curves and F1 summaries.

---

## 10. Regenerate tables and figures

You can rebuild all plots purely from the TSV aggregates.

```bash
# CAMI aggregates and figures
python bench/aggregate_metrics.py --bench-root bench --outdir out
python bench/plot/make_figures.py --bench-root bench --outdir out

# Case-study figures
python case/plot_case.py
python case/plot_ablation.py --summary case/ablation_summary.tsv --eval case/ablation_eval_summary.tsv --outdir case/ablation/figures
```

If the scripts complain about missing dependencies (e.g., matplotlib), ensure the environment is activated or install the required packages. When only a subset of artefacts is available, run the relevant commands selectively.

**Quick refresh:** `bin/hymet artifacts` orchestrates the aggregation and plotting steps above (skipping stages when inputs are missing) to rebuild supplementary artefacts in one command.

---

## 11. Validation and QA

### 11.1 Check hashes against Zenodo

The Zenodo record includes SHA256 sums for the canonical outputs. Compare your reruns:

```bash
sha256sum bench/out/summary_per_tool_per_sample.tsv
sha256sum bench/out/leaderboard_by_rank.tsv
sha256sum bench/out/runtime_memory.tsv
sha256sum case/out/runtime_memory.tsv
sha256sum results/bench/fig_f1_by_rank_lines.png
sha256sum results/case/fig_case_top_taxa_panels.png
```

Record the hashes and confirm they match the deposited versions (see `zenodo_checksums.txt` in the archive).

### 11.2 Inspect metadata

Each HYMET run writes `metadata.json` noting commit hash, sketch checksums, tool versions, and cache key.

```bash
jq '.' bench/out/cami_i_lc/hymet/metadata.json
```

Ensure `hymet_commit` matches the checked-out revision and that the sketch SHA256 values agree with Zenodo.

### 11.3 Unit and smoke tests

```bash
pytest -q
```

All tests should pass; a failure typically indicates missing dependencies or path misconfiguration.

---

## 12. Packaging for submission

1. Create archives of key outputs:
   ```bash
   tar -cf bench_outputs.tar.gz bench/out results/bench
   tar -cf case_outputs.tar.gz case/out case/ablation results/case
   ```
2. Export the environment:
   ```bash
   micromamba env export -n hymet-env > repro_env.yml
   ```
3. Assemble a manifest:
   - `manifest.txt` summarising tarball contents
   - Hash list (`sha256sum bench_outputs.tar.gz case_outputs.tar.gz > hashes.txt`)
4. Store logs (`bench/run_all_cami.log`, any terminal transcripts) alongside the tarballs.

Upload the bundles to Zenodo (or GigaDB) and reference them in the manuscript’s Data Availability statement.

---

## 13. Troubleshooting appendix

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `curl: (7) Failed to connect` when fetching sketches | Proxy/firewall | Set `http_proxy`/`https_proxy` env vars or mirror the Zenodo files manually. |
| `micromamba` cannot resolve `environment.lock.yml` | Platform mismatch | Install from `environment.yml`, then re-lock. |
| `metaphlan` database download stalls | Large DB (>30 GB) | Pre-download via `metaphlan --install` with `--nproc 4`, or skip MetaPhlAn by removing it from `--tools`. |
| `Permission denied` writing to cache | CACHE_ROOT on read-only volume | Point `CACHE_ROOT` to a writable directory. |
| CAMITAX/Nextflow failures | Missing Nextflow/java | Install Nextflow, ensure alternative tools still execute; remove `camitax` from `--tools` if not required. |
| Excessive disk usage | Large `.mmi` indices in cache | Periodically prune with `python bench/tools/prune_cache.py --max-age-days 30 --max-size-gb 150`. |

---

## 14. Minimal completion checklist

☐ Repository cloned, commit hash recorded  
☐ Environment created from `environment.lock.yml`, `pytest -q` passes  
☐ Sketches downloaded, `tools/verify_sketches.sh` reports success  
☐ CAMI benchmark rerun (or verified against existing outputs)  
☐ Case-study and ablation pipelines rerun (or verified)  
☐ Figures regenerated via plotting scripts  
☐ SHA256 sums documented and matched to archival record  
☐ Tarballs + environment export created for submission  

Sign off each item and include the checklist in the reproducibility dossier submitted to reviewers or the journal.
