# HYMET CAMI Benchmark Harness

This directory documents the automation used to benchmark HYMET and selected baseline profilers on CAMI-style datasets. It includes database preparation, tool execution, evaluation against CAMI truth, aggregation, and figure generation. All scripts target a standard Linux workstation or container with sufficient CPU, RAM, and disk.

## 1. Directory Map

```
bench/
├── aggregate_metrics.py      # Merge evaluation outputs → summary TSVs
├── cami_manifest.tsv         # Sample manifest (local paths + optional URLs)
├── convert/                  # Tool output → CAMI format converters
├── data/                     # Lightweight CAMI subsets (generated locally)
├── db/                       # Database builders + cached indices
├── environment.yml           # Conda/micromamba environment spec
├── environment.lock.yml      # Frozen environment capture (optional)
├── fetch_cami.sh             # Optional downloader driven by manifest URLs
├── lib/                      # Shared helpers (common.sh, measure.sh, run_eval.sh, ...)
├── out/                      # Per-run outputs (one folder per sample/tool)
├── plot/                     # Figure generator using matplotlib
├── refsets/                  # Shared reference FASTA subsets
└── run_*.sh                  # Tool-specific wrappers (HYMET, Kraken2, Centrifuge, ganon2, sourmash gather, MetaPhlAn4)
```

Key support scripts:

| Path | Purpose |
|------|---------|
| `lib/common.sh` | Logging, path resolution, root discovery; sourced by every runner. |
| `lib/measure.sh` | Wraps commands with `/usr/bin/time -v`, appends to `out/runtime_memory.tsv`. |
| `lib/run_eval.sh` | Invokes `HYMET/tools/eval_cami.py` and removes empty contig reports. |
| `convert/*.py` | Convert raw outputs into CAMI-compliant profiles. |
| `aggregate_metrics.py` | Builds `summary_per_tool_per_sample.tsv`, `leaderboard_by_rank.tsv`, `contig_accuracy_per_tool.tsv`. |

## 2. Prerequisites

| Requirement | Notes |
|-------------|-------|
| OS | Linux/x86_64 with Bash ≥ 4 and Python ≥ 3.9. |
| Package manager | `micromamba` (recommended) or `conda` for environment creation. |
| CPU/RAM | ~16 threads and ≥ 32 GB RAM. MetaPhlAn4 fits < 20 GB with `METAPHLAN_THREADS=4` + `--split_reads`. |
| Disk | Allocate ~160 GB (MetaPhlAn DB ≈34 GB, HYMET data ≈45 GB, tool outputs ≈2 GB). Remove HYMET minimap2 indices to save ~52 GB per sample. |
| Taxonomy dump | `HYMET/taxonomy_files/` must hold NCBI `names.dmp`, `nodes.dmp`, etc. |
| MetaPhlAn DB | Install with `metaphlan --install -x mpa_vJun23_CHOCOPhlAnSGB_202307 --db_dir bench/db/metaphlan`. |

## 3. Input Data

### 3.1 Manifest (`cami_manifest.tsv`)

Columns:
- `sample_id`
- `contigs_fa`
- `truth_contigs_tsv`
- `truth_profile_tsv`
- Optional `*_url` columns used by `fetch_cami.sh`.

Relative paths are resolved against `bench/`; absolute paths (e.g., `/data/cami/...`) are allowed.

### 3.2 Lightweight CAMI subsets

`bench/data/` stores derived samples (e.g., `cami_i_lc`). Regenerate with:

```bash
python ../tools/generate_cami_subsets.py \
  --fasta /data/cami/sample_0.fna \
  --mapping /data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/gsa_mapping_new.tsv \
  --outdir $(pwd)/data
```

### 3.3 Shared reference FASTA

Database builders expect a shared FASTA (`refsets/combined_subset.fasta`). Build or override via `REF_FASTA=/path/to/fasta`.

## 4. Environment Setup

```bash
cd HYMET/bench
micromamba env create -f environment.yml -n hymet-benchmark
micromamba activate hymet-benchmark
```

Capture an exact environment snapshot after modifications:

```bash
micromamba env export -n hymet-benchmark > environment.lock.yml
```

## 5. Database Builders

All builders live in `bench/db/` and write `.build.stamp` files to avoid redundant work.

| Script | Purpose | Key env variables |
|--------|---------|-------------------|
| `build_kraken2.sh` | Build Kraken2 DB | `REF_FASTA`, `THREADS`, `KRAKEN2_DB_DIR` |
| `build_centrifuge.sh` | Build Centrifuge index | `REF_FASTA`, `THREADS`, `CENTRIFUGE_BMAX`, `CENTRIFUGE_DCV` |
| `build_ganon2.sh` | Build ganon2 HIBF | `REF_FASTA`, `THREADS`, `GANON_FILTER_SIZE` |
| `build_sourmash.sh` | Build sourmash sketch + SBT | `REF_FASTA`, `SOURMASH_KSIZE`, `SOURMASH_SCALED` |

Generate a smaller shared FASTA when RAM is limited:

```bash
python lib/subset_fasta.py \
  --input ../data/downloaded_genomes/combined_genomes.fasta \
  --output refsets/combined_subset.fasta \
  --max-seqs 1000 --max-bases 500000000
export REF_FASTA=$(pwd)/refsets/combined_subset.fasta
```

## 6. Running the Benchmark

### 6.1 One-button driver

```bash
THREADS=16 METAPHLAN_THREADS=4 METAPHLAN_OPTS="--split_reads" \
  ./run_all_cami.sh --tools hymet,kraken2,centrifuge,ganon2,sourmash_gather,metaphlan4
```

Useful options:
- `--tools` – comma-separated list or `all`.
- `--no-build` – reuse existing databases.
- `--threads N` – override default thread count.
- `--max-samples N` – process first N manifest rows.
- `--resume` – keep existing `out/runtime_memory.tsv`.

Important environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `THREADS` | 8 | Fallback thread count. |
| `KEEP_HYMET_WORK` | 0 | Delete HYMET minimap2 indices unless set to 1. |
| `KEEP_GANON_RAW`, etc. | 0 | Preserve heavy raw files when set to 1. |
| `GANON_REL_CUTOFF` | 0 | Disable ganon2 relative cutoff (increase to tighten). |
| `GANON_REL_FILTER` | 1 | Keep all hits relative to best match. |
| `SOURMASH_TOP_HITS` | 500 | Limit gather hits. |
| `METAPHLAN_DB_DIR` | `bench/db/metaphlan` | Location of MetaPhlAn indices. |
| `METAPHLAN_INDEX` | `mpa_vJun23_CHOCOPhlAnSGB_202307` | MetaPhlAn index name. |
| `METAPHLAN_THREADS` | `THREADS` | Bowtie2 threads. |
| `METAPHLAN_OPTS` | (empty) | Additional MetaPhlAn options. |

### 6.2 Individual tool runs

Call the runner directly:

```bash
./run_ganon2.sh --sample cami_i_lc --contigs data/cami_i_lc/contigs.fna --threads 16
```

Then regenerate aggregates:

```bash
python aggregate_metrics.py --bench-root . --outdir out
python plot/make_figures.py  --bench-root . --outdir out  # optional
```

### 6.3 Evaluation-only reruns

Useful after manual tweaks to converter scripts:

```bash
./lib/run_eval.sh --sample cami_sample_0 --tool ganon2 \
  --pred-profile out/cami_sample_0/ganon2/profile.cami.tsv \
  --truth-profile /data/cami/sample_0/taxonomic_profile_0.txt \
  --pred-contigs out/cami_sample_0/ganon2/classified_sequences.tsv \
  --truth-contigs /data/cami/sample_0/.../gsa_mapping_new.tsv \
  --pred-fasta /data/cami/sample_0.fna --threads 16
python aggregate_metrics.py --bench-root . --outdir out
```

## 7. Output Layout

`out/<sample>/<tool>/` typically contains:

| File/Dir | Description |
|----------|-------------|
| `profile.cami.tsv` | CAMI profile predicted by the tool. |
| `classified_sequences.tsv` | Per-contig assignments (when available). |
| `metadata.json` | Provenance info (sample, tool, key paths). |
| `resultados.paf` | HYMET PAF alignment (if produced). |
| `eval/profile_summary.tsv` | Rank-wise abundance metrics. |
| `eval/contigs_exact.tsv` + `eval/contigs_per_rank.tsv` | Contig metrics (omitted if no usable pairs). |
| `eval/_debug_info.txt` | Diagnostic input summary. |

`aggregate_metrics.py` writes:
- `summary_per_tool_per_sample.tsv`
- `leaderboard_by_rank.tsv`
- `contig_accuracy_per_tool.tsv` (rows with `n <= 0` filtered)
- `runtime_memory.tsv`

`plot/make_figures.py` produces:
- `fig_accuracy_by_rank.png`
- `fig_f1_by_rank.png`
- `fig_l1_braycurtis.png`
- `fig_per_sample_f1_stack.png`

## 8. Resource Tips & Troubleshooting

- **Disk usage**: HYMET minimap2 index (~52 GB) is removed by default; set `KEEP_HYMET_WORK=1` to keep it. The MetaPhlAn `.tar` download (12 GB) can be deleted after extraction. Clear `out/` between runs to reclaim space.
- **ganon2 coverage**: Relaxed defaults (`GANON_REL_CUTOFF=0`, `GANON_REL_FILTER=1`) keep long-contig matches. Tighten thresholds if precision is a concern.
- **MetaPhlAn memory**: Use `METAPHLAN_THREADS=4` and `--split_reads` to stay under 20 GB.
- **No contig output**: MetaPhlAn4 and sourmash gather do not emit contig assignments. The evaluator now removes empty reports, and aggregates omit those rows.
- **Logs**: `runtime_memory.tsv` captures command lines and resources. `_debug_info.txt` in each `eval/` folder lists the evaluation inputs.

## 9. Suggested Workflow

1. Activate environment.
2. Validate `cami_manifest.tsv` paths.
3. Build databases (`db/build_*.sh`).
4. Ensure MetaPhlAn DB exists under `bench/db/metaphlan`.
5. Run benchmark (see §6.1).
6. Regenerate aggregates/plots if needed.
7. Review outputs in `bench/out/`.

## 10. Extending the Harness

- **New tool**: add `run_<tool>.sh`, converter in `convert/`, register in `run_all_cami.sh`, document here.
- **Custom datasets**: update `cami_manifest.tsv` and regenerate subsets with `tools/generate_cami_subsets.py`.
- **Plots**: extend `plot/make_figures.py` or add new scripts and call them after aggregation.
- **Environment**: update `environment.yml`, regenerate `environment.lock.yml`, mention changes in README.

## 11. Case Study Integration

The case-study toolkit under `case/` reuses the HYMET runner and database filters described here. After completing a CAMI benchmark, run `case/run_case.sh` to produce real-data summaries and `case/run_ablation.sh` to quantify the impact of removing dominant taxa from the shared FASTA.

## 12. Checklist Before Sharing Results

- [ ] Environment created from `environment.yml` (or `environment.lock.yml`).
- [ ] MetaPhlAn database installed and referenced via `METAPHLAN_DB_DIR`.
- [ ] Manifest paths verified (`ls` each entry).
- [ ] Database builders completed (check `.build.stamp`).
- [ ] Benchmark run (`run_all_cami.sh` exit status 0).
- [ ] Aggregates regenerated (`aggregate_metrics.py`).
- [ ] Disk usage monitored (`df -h`).
- [ ] `bench/out/` archived with TSVs, figures, `runtime_memory.tsv`.

With these steps, the harness yields reproducible CAMI benchmark results suitable for regression tracking or publication.
