# HYMET CAMI Benchmark Harness

This directory documents the automation used to benchmark HYMET and selected baseline profilers on CAMI-style datasets. It includes database preparation, tool execution, evaluation against CAMI truth, aggregation, and figure generation. All scripts target a standard Linux workstation or container with sufficient CPU, RAM, and disk.

## 1. Directory Map

```
bench/
├── aggregate_metrics.py      # Merge evaluation outputs → summary TSVs
├── bin/                      # Helper entry points (e.g., Nextflow launcher for CAMITAX)
├── cami_manifest.tsv         # Sample manifest (local paths + optional URLs)
├── config/                   # Parameter templates for external profilers (TAMA, etc.)
├── convert/                  # Tool output → CAMI format converters
├── data/                     # Lightweight CAMI subsets (generated locally)
├── db/                       # Database builders + cached indices
├── docs/                     # Supplemental notes and walkthroughs
├── environment.yml           # Conda/micromamba environment spec
├── environment.lock.yml      # Frozen environment capture (optional)
├── fetch_cami.sh             # Optional downloader driven by manifest URLs
├── lib/                      # Shared helpers (common.sh, measure.sh, run_eval.sh, ...)
├── nextflow/                 # Local Nextflow assets and cached work directories
├── out/                      # Per-run outputs (one folder per sample/tool)
├── plot/                     # Figure generator using matplotlib
├── refsets/                  # Shared reference FASTA subsets
├── results_summary.md        # Rolling benchmark status and tuning notes
├── run_all_cami.sh           # Batch driver that orchestrates every configured tool
├── run_*.sh                  # Tool-specific wrappers (HYMET, Kraken2, Centrifuge, ganon2, sourmash gather, MetaPhlAn4, CAMITAX, BASTA, PhaBOX, phyloFlash, ViWrap/geNomad, SqueezeMeta, MegaPath-Nano, SnakeMAGs, TAMA, etc.)
├── tmp_downloads/            # Staging area for large third-party downloads
└── tools/                    # One-off utilities (cache pruning, subset generation, taxonomy helpers)
```

Key support scripts:

| Path | Purpose |
|------|---------|
| `lib/common.sh` | Logging, path resolution, root discovery; sourced by every runner. |
| `lib/measure.sh` | Wraps commands with `/usr/bin/time -v`, appends to `out/runtime_memory.tsv`. |
| `lib/run_eval.sh` | Invokes `HYMET/tools/eval_cami.py` and removes empty contig reports. |
| `convert/*.py` | Convert raw outputs into CAMI-compliant profiles. |
| `aggregate_metrics.py` | Builds `summary_per_tool_per_sample.tsv`, `leaderboard_by_rank.tsv`, `contig_accuracy_per_tool.tsv`. |
| `plot/make_figures.py` | Generates accuracy/F1/abundance/resource figures from aggregated metrics. |

## 2. Prerequisites

| Requirement | Notes |
|-------------|-------|
| OS | Linux/x86_64 with Bash ≥ 4 and Python ≥ 3.9. |
| Package manager | `micromamba` (recommended) or `conda` for environment creation. |
| CPU/RAM | ~16 threads and ≥ 32 GB RAM. MetaPhlAn4 fits < 20 GB with `METAPHLAN_THREADS=4` + `--split_reads`. |
| Disk | Allocate ~160 GB (MetaPhlAn DB ≈34 GB, HYMET data ≈45 GB, tool outputs ≈2 GB). Remove HYMET minimap2 indices to save ~52 GB per sample. |
| Taxonomy dump | `HYMET/taxonomy_files/` must hold NCBI `names.dmp`, `nodes.dmp`, etc. |
| Nextflow | Installable via `curl -s https://get.nextflow.io \| bash` (place `nextflow` on PATH or keep in `bench/`). Required for CAMITAX integration. |
| BASTA + BLAST | Install BLAST+ (`blastx`) and the BASTA CLI; provide a compatible BLAST database for LCA assignment. |
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

### 4.1 Reproducibility

- Prefer running with the locked environment (`environment.lock.yml`) for exact versions.
- Database locations are passed via environment variables; avoid hard‑coding paths in scripts.
- The harness writes all per‑run artefacts under `bench/out/` and regenerates aggregate tables and figures deterministically from those artefacts.
- Recovery utilities that attempted to “rebuild” or “restore” runtime logs from previous commits have been removed to keep results provenance clear. If a runtime table is missing, re‑run the affected tool or stage instead of reconstructing from history.

## 5. Database Builders & External Tool Setup

All builders live in `bench/db/` and write `.build.stamp` files to avoid redundant work.

| Script | Purpose | Key env variables |
|--------|---------|-------------------|
| `build_kraken2.sh` | Build Kraken2 DB | `REF_FASTA`, `THREADS`, `KRAKEN2_DB_DIR` |
| `build_centrifuge.sh` | Build Centrifuge index | `REF_FASTA`, `THREADS`, `CENTRIFUGE_BMAX`, `CENTRIFUGE_DCV` |
| `build_ganon2.sh` | Build ganon2 HIBF | `REF_FASTA`, `THREADS`, `GANON_FILTER_SIZE` |
| `build_sourmash.sh` | Build sourmash sketch + SBT | `REF_FASTA`, `SOURMASH_KSIZE`, `SOURMASH_SCALED` |

### 5.1 CAMITAX Reference & Runtime

CAMITAX relies on a large reference bundle and several third-party binaries outside the default benchmark environment.

1. **Create the runtime environment** (uses micromamba for reproducibility):
   ```bash
   micromamba create -y -n camitax -c conda-forge -c bioconda \
     python=3.10 kaiju prodigal checkm-genome hmmer pplacer bioconductor-dada2 tbb=2020.2
   ```
   The explicit `tbb` pin avoids a missing `tbb::task` symbol when loading `dada2`.

2. **Download the CAMITAX reference database** (≈27 GB):
   ```bash
   ./nextflow run CAMI-challenge/CAMITAX/init.nf --db /path/to/camitax_db -work-dir /path/to/camitax_work
   ```
   If `nextflow pull CAMI-challenge/CAMITAX` is blocked, clone the repository manually and run `nextflow run /path/to/CAMITAX/init.nf …`.

3. **Wrapper script**: `bin/nextflow_camitax.sh` executes Nextflow inside the `camitax` environment. Point `NEXTFLOW_CMD` to this script when invoking the benchmark:
   ```bash
   export NEXTFLOW_CMD="$(pwd)/bin/nextflow_camitax.sh"
   export CAMITAX_DB=/path/to/camitax_db
   export CAMITAX_PIPELINE=/path/to/CAMITAX/main.nf   # optional, defaults to remote ID
   export CAMITAX_EXTRA_OPTS="-c $(pwd)/lib/camitax_local.config -without-docker -resume"
   ```
   `lib/camitax_local.config` disables container usage and reduces default CPU/memory footprints for local execution.

4. **Run CAMITAX via the harness** (other profilers remain untouched):
   ```bash
   THREADS=8 KEEP_CAMITAX_WORK=1 ./run_all_cami.sh --tools camitax --no-build --resume
  ```
  The `KEEP_CAMITAX_WORK` flag preserves per-sample Nextflow work directories inside `out/<sample>/camitax/run/` for easier debugging. Omit it to reclaim disk space.

### 5.2 BASTA Workflow

BASTA consumes translated alignments and assigns taxonomy via an LCA strategy. The wrapper prefers `diamond blastx` whenever `diamond` is available, but you can force classic BLAST+ by exporting `BASTA_USE_DIAMOND=0`.

1. **Install prerequisites**:
   - BLAST+ (`blastx`) available on `PATH` for the fallback path.
   - [`diamond`](https://github.com/bbuchfink/diamond) (recommended); set `BASTA_DIAMOND_DB` if the DIAMOND database does not share the BLAST prefix.
   - BASTA CLI (`basta` command).
   - A BLAST/DIAMOND-formatted protein database (e.g. UniProt).
   - Optional: `TAXONKIT_DB` pointing to `HYMET/taxonomy_files` to accelerate lineage lookups.

2. **Execute on a single sample**:
   ```bash
   export BASTA_BLAST_DB=/path/to/blast/db/prefix
   export BASTA_DIAMOND_DB=/path/to/diamond/db.dmnd   # optional if DIAMOND uses a different basename
   export DIAMOND_EXTRA_OPTS="--fast"                 # forwarded to diamond blastx
   export BASTA_TAXON_MODE=uni
   THREADS=8 ./run_basta.sh --sample cami_sample_0 --contigs /data/cami/sample_0.fna
   ```
   Additional environment toggles:
   - `BASTA_USE_DIAMOND=1` – force DIAMOND even when both databases exist.
   - `BLAST_MAX_TARGETS=50` – adjust translated hit fan-out.
   - `BLASTX_EXTRA_OPTS="--max_hsps 5"`
   - `BASTA_EXTRA_OPTS="-d /root/.basta/taxonomy"` – custom taxonomy directory.
   - `BLASTX_CMD` / `BASTA_CMD` – override executable paths.

3. **Batch mode** (preserves existing profiler outputs):
   ```bash
   THREADS=8 ./run_all_cami.sh --tools basta --no-build --resume
   ```

4. **Outputs**:
   - `out/<sample>/basta/profile.cami.tsv`
   - `out/<sample>/basta/classified_sequences.tsv`
   - Intermediate files cached under `out/<sample>/basta/run/`

### 5.3 PhaBOX Workflow

PhaBOX performs phage lifestyle/host predictions. The harness relabels contigs, calls the CLI, and converts `phagcn_prediction.tsv` into CAMI-ready outputs.

1. **Install prerequisites**:
   - Ensure `phabox2` (or equivalent) is on `PATH`, or set `PHABOX_CMD`.
   - Download the PhaBOX database and set `PHABOX_DB_DIR`.
   - Optional: export `PHABOX_WORKDIR` if the CLI must run from a specific directory.
   - `run_phabox.sh` auto-installs [`prodigal-gv`](https://github.com/apcamargo/prodigal-gv) and [`taxonkit`](https://bioinf.shenwei.me/taxonkit/) via `micromamba` (Bioconda) when missing; install them manually if `micromamba` is unavailable.

2. **Execute on a single sample**:
   ```bash
   export PHABOX_DB_DIR=/path/to/phabox_db_v2.0.0
   export PHABOX_CMD=phabox2                     # e.g., "python /opt/phabox/bin/phabox2"
   THREADS=8 ./run_phabox.sh --sample cami_sample_0 --contigs /data/cami/sample_0.fna
   ```
   Additional toggles:
   - `PHABOX_TASK=phagcn` (default task)
   - `PHABOX_EXTRA_OPTS="--minlen 5000"` (forwarded to the CLI)

3. **Batch mode**:
   ```bash
   THREADS=8 ./run_all_cami.sh --tools phabox --no-build --resume
   ```

4. **Outputs**:
   - `out/<sample>/phabox/profile.cami.tsv`
   - `out/<sample>/phabox/classified_sequences.tsv`
   - Intermediate artefacts under `out/<sample>/phabox/run/` (converted FASTA, ID map, raw prediction TSV)

Generate a smaller shared FASTA when RAM is limited:

```bash
python lib/subset_fasta.py \
  --input ../data/downloaded_genomes/combined_genomes.fasta \
  --output refsets/combined_subset.fasta \
  --max-seqs 1000 --max-bases 500000000
export REF_FASTA=$(pwd)/refsets/combined_subset.fasta
```

### 5.4 ViWrap (geNomad) Workflow

The ViWrap integration invokes `genomad end-to-end` on assembled contigs and converts the virus summary into CAMI outputs.

1. **Install prerequisites**:
   - Create a dedicated environment (default `/opt/envs/genomad`):
     ```bash
     micromamba create -y -p /opt/envs/genomad -c conda-forge -c bioconda genomad
     ```
   - Download the geNomad database to shared storage (e.g. `/data/ref/viwrap/genomad_db`):
     ```bash
     micromamba run -p /opt/envs/genomad genomad download-database /data/ref/viwrap/genomad_db
     ```
2. **Execute on a single sample**:
   ```bash
   export VIWRAP_DB_DIR=/data/ref/viwrap/genomad_db/genomad_db
   THREADS=8 ./run_viwrap.sh --sample cami_sample_0 --contigs /data/cami/sample_0.fna
   ```
   Optional toggles:
   - `VIWRAP_ENV_PREFIX=/custom/genomad/env`
   - `VIWRAP_SCORE_CUTOFF=0.7` (minimum `virus_score` retained)
   - `VIWRAP_EXTRA_OPTS="--lenient-taxonomy"` (forwarded to `genomad end-to-end`)
3. **Batch mode**:
   ```bash
   THREADS=8 ./run_all_cami.sh --tools viwrap --no-build --resume
   ```
4. **Outputs**:
   - `out/<sample>/viwrap/profile.cami.tsv`
   - `out/<sample>/viwrap/classified_sequences.tsv`
   - Raw geNomad artefacts under `out/<sample>/viwrap/run/`.

### 5.5 SqueezeMeta Workflow

SqueezeMeta sequential mode can operate directly on pre-assembled contigs. The wrapper fabricates lightweight reads to satisfy the pipeline, runs `SqueezeMeta.pl -extassembly`, and converts the resulting contig taxonomy table.

1. **Install prerequisites**:
   - Create the environment (default `/opt/envs/squeezemeta`):
     ```bash
     micromamba create -y -p /opt/envs/squeezemeta -c conda-forge -c bioconda -c anaconda -c fpusan squeezemeta=1.6 --no-channel-priority
     ```
   - Populate the SqueezeMeta database somewhere with plenty of space (recommend `/data/ref/squeezemeta`). Use the upstream `download_databases.pl` helper and point `SQUEEZEMETA_DB_DIR` at the extracted directory.
2. **Execute on a single sample**:
   ```bash
   export SQUEEZEMETA_DB_DIR=/data/ref/squeezemeta
   THREADS=8 ./run_squeezemeta.sh --sample cami_sample_0 --contigs /data/cami/sample_0.fna
   ```
   Optional toggles:
   - `SQUEEZEMETA_ENV_PREFIX=/custom/squeezemeta/env`
   - `SQUEEZEMETA_SYNTH_FRAG_LEN=200` to adjust synthetic fragment size.
   - `SQUEEZEMETA_EXTRA_OPTS="--nobins --nodiamond"` to tailor the workflow.
3. **Batch mode**:
   ```bash
   THREADS=8 ./run_all_cami.sh --tools squeezemeta --no-build --resume
   ```
4. **Outputs**:
   - `out/<sample>/squeezemeta/profile.cami.tsv`
   - `out/<sample>/squeezemeta/classified_sequences.tsv`
   - Full SqueezeMeta project retained under `out/<sample>/squeezemeta/run/<sample>/`.

### 5.6 PhyloFlash Workflow

PhyloFlash detects and profiles SSU rRNA from assembled contigs by synthesizing pseudo-reads and running `phyloFlash.pl`.

1. **Install prerequisites**:
   - Install phyloFlash and its dependencies into a dedicated environment (default path `/opt/envs/phyloflash`). Example:
     ```bash
     micromamba create -y -p /opt/envs/phyloflash -c conda-forge -c bioconda phyloflash barrnap bedtools seqtk
     ```
   - Download the SILVA-derived database (e.g. `138.1`) and unpack it somewhere with plenty of space, such as `/data/ref/phyloflash/138.1`.
2. **Execute on a single sample**:
   ```bash
   export PHYLOFLASH_DB_DIR=/data/ref/phyloflash/138.1
   THREADS=8 ./run_phyloflash.sh --sample cami_sample_0 --contigs /data/cami/sample_0.fna
   ```
   Optional environment toggles:
   - `PHYLOFLASH_ENV_PREFIX=/custom/env/prefix`
   - `PHYLOFLASH_EXTRA_OPTS="--zip"` (passed straight to `phyloFlash.pl`)
   - `PHYLOFLASH_FRAGMENT_LEN=250` and `PHYLOFLASH_MIN_FRAGMENT_LEN=60` to adjust pseudo-read tiling.
3. **Batch mode**:
   ```bash
   THREADS=8 ./run_all_cami.sh --tools phyloflash --no-build --resume
   ```
4. **Outputs**:
   - `out/<sample>/phyloflash/profile.cami.tsv`
   - `out/<sample>/phyloflash/classified_sequences.tsv` (header only when no per-read assignments are available)
   - Intermediate artefacts (GFF, rRNA FASTA, pseudo-reads, raw phyloFlash output) under `out/<sample>/phyloflash/run/`.

### 5.7 TAMA Workflow

The harness includes a reproducible TAMA setup driven by static parameter files.

1. **Install / locate TAMA** and its dependent databases. Note the directory containing `TAMA.pl` (set as `TAMA_ROOT` below).
2. **Generate synthetic reads** from the benchmark contigs so TAMA receives FASTQ input:
   ```bash
   mkdir -p data/tama_reads
   while IFS=$'\t' read -r sample contigs _; do
     [[ $sample == sample_id || $sample == \#* ]] && continue
     python tools/contigs_to_reads.py \
       --contigs "$(realpath "$contigs")" \
       --out "data/tama_reads/${sample}.fastq"
   done < cami_manifest.tsv
   ```
   This slices contigs into 250 bp windows (tail chunks ≥100 bp) and writes high-quality single-end reads for each CAMI sample.
3. **Source the helper environment script** to register shared locations:
   ```bash
   source config/tama_env.sh          # sets BENCH_ROOT and TAMA_PARAM_DIR
   export TAMA_ROOT=/abs/path/to/TAMA # update to your installation
   ```
   `config/tama_params/` already contains one parameter file per CAMI sample, each now pointing to `data/tama_reads/<sample>.fastq`. The helper script prepends the micromamba `perl` toolchain (which supplies `perl-sort-key`) so that abundance estimation succeeds without extra tweaks. The default params enable only the Centrifuge and Kraken backends because the CLARK database requires >150 GB RAM to load on this host.
4. **Run a single sample** (keeps other tools untouched):
   ```bash
   THREADS=8 TAMA_PARAM_FILE="" ./run_tama.sh \
     --sample cami_sample_0 \
     --contigs /data/cami/sample_0.fna
   ```
   Leaving `TAMA_PARAM_FILE` empty allows the wrapper to pick the matching file from `config/tama_params/<sample>.params.txt`. Set `KEEP_TAMA_WORK=1` to preserve the large intermediate directory under `out/<sample>/tama/run/`.
5. **Batch mode**:
   ```bash
   THREADS=8 ./run_all_cami.sh --tools tama --no-build --resume
   ```
   Combine with other profilers as required: `--tools hymet,kraken2,tama`.
6. **Outputs**:
   - `out/<sample>/tama/profile.cami.tsv`
   - `out/<sample>/tama/classified_sequences.tsv` (if `read_classi*.out` is present)
   - The parameter file used for the run is copied to `out/<sample>/tama/params.txt` for traceability.

### 5.8 MegaPath-Nano Workflow

MegaPath-Nano targets long-read (ONT) data. The CAMI harness fabricates long-read surrogates from contigs, runs the taxonomy-only pipeline, and converts the resulting microbe statistics.

1. **Install prerequisites**:
   - Clone the upstream repository and compile its bundled dependencies (see `Syst_Review/Configuration/megapathnano.sh` for a scripted walkthrough):
     ```bash
     git clone https://github.com/HKU-BAL/MegaPath-Nano /opt/tools/MegaPath-Nano
     ```
   - Populate the MegaPath-Nano databases (`./install_db.sh` within the repo) and ensure `minimap2`, `porechop`, and other bundled binaries are operational inside the tree.
   - Optional: create a dedicated Conda environment with the packages referenced by `megapath_nano.py`.
2. **Single-sample execution**:
   ```bash
   export MPN_ROOT=/opt/tools/MegaPath-Nano
   # Optional overrides:
   # export MPN_CMD="micromamba run -n mpn python3"
   # export MPN_CHUNK_SIZE=12000   # synthetic long-read length
   THREADS=8 ./run_megapath_nano.sh \
     --sample cami_sample_0 \
     --contigs /data/cami/sample_0.fna
   ```
   The wrapper forces `--taxon_module_only`, disables heavy AMR/noise outputs, caps aligner threads (≤64), and prunes temporary files unless `MPN_KEEP_WORK=1`.
3. **Batch mode**:
   ```bash
   THREADS=8 ./run_all_cami.sh --tools megapath_nano --no-build --resume
   ```
4. **Outputs**:
   - `out/<sample>/megapath_nano/profile.cami.tsv`
   - `out/<sample>/megapath_nano/classified_sequences.tsv`
   - Minimal run artefacts (ID map, synthetic FASTQ, raw stats) under `out/<sample>/megapath_nano/run/` when `MPN_KEEP_WORK=1`.

### 5.9 SnakeMAGs Classification Path

SnakeMAGs is a Snakemake workflow that normally consumes raw Illumina reads. For CAMI contig benchmarks we expose a lightweight classification path that treats long contigs as provisional MAGs and runs GTDB-Tk via Snakemake.

1. **Install prerequisites**:
   - Provide a working `snakemake` binary (`micromamba install -c conda-forge -c bioconda snakemake`).
   - Download and unpack the GTDB-Tk data bundle; export `SNAKEMAGS_GTDB` to that directory (e.g., `/data/ref/gtdb/release214`).
2. **Single-sample execution**:
   ```bash
   export SNAKEMAGS_GTDB=/data/ref/gtdb/release214
   # Optional overrides:
   # export SNAKEMAGS_MIN_CONTIG=8000          # retain only long contigs
   # export SNAKEMAGS_SNKMK_CMD="snakemake -p" # print shell commands
   THREADS=8 ./run_snakemags.sh \
     --sample cami_sample_0 \
     --contigs /data/cami/sample_0.fna
   ```
   Contigs shorter than `SNAKEMAGS_MIN_CONTIG` are skipped to reduce false positives and runtime. Intermediate MAG FASTA files and the Snakemake workdir live under `out/<sample>/snakemags/run/`.
3. **Batch mode**:
   ```bash
   THREADS=8 ./run_all_cami.sh --tools snakemags --no-build --resume
   ```
4. **Outputs**:
   - `out/<sample>/snakemags/profile.cami.tsv`
   - `out/<sample>/snakemags/classified_sequences.tsv`
  - GTDB-Tk summaries retained under `out/<sample>/snakemags/run/` unless `SNAKEMAGS_KEEP_WORK=0` (default cleans up).

### 5.10 Metalign Workflow

Metalign is a read-based profiler. The harness synthesizes single-end reads from contigs, runs Metalign, and normalizes its CAMI-style output.

1. Install prerequisites and database:
   - Preferred: install via Bioconda. The builder will install Metalign into `/opt/conda` if missing and fetch the reference data into `bench/db/metalign/data`:
     ```bash
     ./db/build_metalign.sh
     ```
   - Override locations via environment:
     - `METALIGN_ENV_PREFIX=/custom/conda/prefix`
     - `METALIGN_DB_DIR=/data/ref/metalign/data` (to avoid re-downloading)

2. Execute on a single sample:
   ```bash
   THREADS=8 ./run_metalign.sh --sample cami_sample_0 --contigs /data/cami/sample_0.fna
   ```
   Optional toggles:
   - `METALIGN_PRESET=--precise` (or `--sensitive`)
   - `METALIGN_OPTS="--length_normalize --rank_renormalize"`
   - `METALIGN_CHUNK_SIZE=250 METALIGN_MIN_CHUNK=100` to adjust synthetic reads.

3. Batch mode:
   ```bash
   THREADS=8 ./run_all_cami.sh --tools metalign --no-build --resume
   ```

4. Outputs:
   - `out/<sample>/metalign/profile.cami.tsv`
   - `out/<sample>/metalign/metalign_abundances.tsv`
   - Synthetic reads under `out/<sample>/metalign/run/`

## 6. Running the Benchmark

### 6.1 HYMET CLI shortcuts

The `bin/hymet` wrapper provides friendlier entry points around the harness while inheriting the same environment variables as the shell scripts.

```bash
# Single-sample HYMET run
bin/hymet run \
  --contigs /path/to/contigs.fna \
  --out /path/to/output \
  --threads 16

# Full CAMI manifest across multiple profilers
bin/hymet bench \
  --manifest bench/cami_manifest.tsv \
  --tools hymet,kraken2,centrifuge
```

Run `bin/hymet bench --help` for advanced options (`--samples`, `--resume`, etc.).

### 6.2 One-button driver

```bash
THREADS=16 METAPHLAN_THREADS=4 METAPHLAN_OPTS="--split_reads" \
  ./run_all_cami.sh --tools hymet,kraken2,centrifuge,ganon2,sourmash_gather,metaphlan4,camitax
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

### 6.3 Individual tool runs

Call the runner directly:

```bash
./run_ganon2.sh --sample cami_i_lc --contigs data/cami_i_lc/contigs.fna --threads 16
```

Then regenerate aggregates and sync the `results/bench/` snapshots:

```bash
./publish_results.sh
```

### 6.4 Evaluation-only reruns

Useful after manual tweaks to converter scripts:

```bash
./lib/run_eval.sh --sample cami_sample_0 --tool ganon2 \
  --pred-profile out/cami_sample_0/ganon2/profile.cami.tsv \
  --truth-profile /data/cami/sample_0/taxonomic_profile_0.txt \
  --pred-contigs out/cami_sample_0/ganon2/classified_sequences.tsv \
  --truth-contigs /data/cami/sample_0/.../gsa_mapping_new.tsv \
  --pred-fasta /data/cami/sample_0.fna --threads 16
./publish_results.sh
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
- `fig_cpu_time_by_tool.png`
- `fig_wall_time_by_tool.png`
- `fig_peak_memory_by_tool.png`
- `fig_contig_accuracy_heatmap.png`

`run_all_cami.sh` now invokes `publish_results.sh`, so aggregate tables and figures are refreshed under both `bench/out/` and `results/bench/`.

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
5. Run benchmark (see §6.2).
6. Regenerate aggregates/plots if needed.
7. Review outputs in `bench/out/`.

## 10. Extending the Harness

- **New tool**: add `run_<tool>.sh`, converter in `convert/`, register in `run_all_cami.sh`, document here.
- **Custom datasets**: update `cami_manifest.tsv` and regenerate subsets with `tools/generate_cami_subsets.py`.
- **Plots**: extend `plot/make_figures.py` or add new scripts and call them after aggregation.
- **Environment**: update `environment.yml`, regenerate `environment.lock.yml`, mention changes in README.

## 11. Case Study Integration

The case-study toolkit under `case/` reuses the HYMET runner and database filters described here. After completing a CAMI benchmark, run `case/run_case.sh` to produce real-data summaries and `case/run_ablation.sh` to quantify the impact of removing dominant taxa from the shared FASTA.
