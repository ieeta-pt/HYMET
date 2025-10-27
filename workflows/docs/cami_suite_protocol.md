# HYMET CAMI Suite Protocol

This protocol describes how to use `workflows/run_cami_suite.sh` to generate reproducible CAMI benchmark suites. Each execution produces `results/<scenario>/<suite>/run_<timestamp>/` containing:

- `raw/` – full bench outputs grouped by mode (e.g. `raw/contigs/<sample>/<tool>/`).
- `tables/<mode>/` – `summary_per_tool_per_sample.tsv`, `leaderboard_by_rank.tsv`, `runtime_memory.tsv`, `contig_accuracy_per_tool.tsv`, `manifest.snapshot.tsv`.
- `figures/<mode>/` – `fig_*.png` regenerated from the raw data.
- `metadata.json` – manifest path, git commit, commands, tool roster, environment.

The canonical manuscript run resides in `results/cami/canonical/run_<timestamp>/`. Additional suites (contig-only, contig vs reads, reviewer-specific panels) are stored alongside it without ever touching `bench/out/`.

---

## 1. Prerequisites

| Requirement | Recommendation | Notes |
|-------------|----------------|-------|
| CPU | ≥16 threads | HYMET/Minimap2 scale with more cores. |
| RAM | 64 GB | Kraken2, Centrifuge, and SqueezeMeta may use >32 GB. |
| Disk | ≥250 GB free | Reference caches + multiple suite runs. |
| OS | Linux (Ubuntu 22.04 tested) | macOS works with Conda/Mamba. |
| Software | `git`, `curl`, `wget`, `tar`, `sha256sum`, `python ≥3.9`, `micromamba` | Install before starting. |

Clone the repository, create the `hymet-env`, and download Mash sketches + manifests as described in `docs/reproducibility.md`.

---

## 2. Running a suite

### 2.1 Default contig+read panel

```bash
cd HYMET
THREADS=16 CACHE_ROOT=data/downloaded_genomes/cache_bench \
workflows/run_cami_suite.sh \
  --scenario cami \
  --suite contig_full
```

This runs the default contig tools (`hymet,kraken2,centrifuge,ganon2,viwrap,tama,squeezemeta,megapath_nano`) plus the `hymet_reads` pseudo-read variant. Outputs appear in `results/cami/contig_full/run_<timestamp>/`.

### 2.2 Custom panels

```bash
workflows/run_cami_suite.sh \
  --scenario cami \
  --suite reviewer_set \
  --contig-tools hymet,kraken2,ganon2 \
  --read-tools hymet_reads \
  --threads 24
```

Use `--manifest` to point at alternate CAMI subsets, `--modes` to drop read mode, and `--bench-extra` to forward additional flags to `bin/hymet bench` (e.g. `--resume`).

### 2.3 Dry-run mode

Add `--dry-run` to record metadata and planned commands without executing the benchmark. This is useful when previewing tool lists or verifying manifests.

---

## 3. Outputs & inspection

For a given run `results/<scenario>/<suite>/run_<timestamp>/`:

| Path | Contents |
|------|----------|
| `raw/contigs/<sample>/<tool>/` | Per-tool CAMI outputs (classified_sequences.tsv, profile.cami.tsv, eval/…). |
| `raw/reads/...` | Pseudo-read outputs (if `hymet_reads` enabled). |
| `tables/<mode>/summary_per_tool_per_sample.tsv` | Precision/recall/F1 for all tools in that mode. |
| `tables/<mode>/leaderboard_by_rank.tsv` | Rank-level averages. |
| `tables/<mode>/runtime_memory.tsv` | `/usr/bin/time` measurements per tool/stage. |
| `figures/<mode>/fig_*.png` | Accuracy/L1/peak-memory plots regenerated from the raw data. |
| `metadata.json` | Timestamp, git commit, manifest path, tool roster, commands, environment variables. |

To regenerate figures for a specific mode, point the plotting helper at that mode’s `tables/` directory:

```bash
python bench/plot/make_figures.py \
  --bench-root bench \
  --tables results/cami/<suite>/run_<timestamp>/tables/contigs \
  --outdir results/cami/<suite>/run_<timestamp>/figures/contigs
```

---

## 4. Provenance checklist

| Item | Command | Location |
|------|---------|----------|
| Git revision | `git rev-parse HEAD` | `metadata.json` (`git_commit`). |
| Environment export | `micromamba env export -n hymet-env > hymet-env.yml` | Archive alongside results. |
| Manifest snapshot | Already copied to `tables/<mode>/manifest.snapshot.tsv`. |
| Aggregate checksum | `sha256sum results/<scenario>/<suite>/run_<timestamp>/tables/<mode>/summary_per_tool_per_sample.tsv` | Record per run. |
| Command log | Stored in `metadata.json[
| Command log | Stored in `metadata.json["commands"]` | Auto |

---

## 5. Packaging & sharing

```bash
SUITE_DIR=results/cami/contig_full/run_<timestamp>
mkdir -p share
cp "$SUITE_DIR"/metadata.json share/
cp "$SUITE_DIR"/tables/contigs/summary_per_tool_per_sample.tsv share/
cp "$SUITE_DIR"/tables/contigs/leaderboard_by_rank.tsv share/
cp "$SUITE_DIR"/tables/contigs/runtime_memory.tsv share/
cp "$SUITE_DIR"/figures/contigs/fig_* share/
micromamba env export -n hymet-env > share/hymet-env.yml
tar -C share -czf cami_suite_results.tar.gz .
```

Publish the tarball alongside the manifest snapshot stored in `tables/<mode>/manifest.snapshot.tsv` and, if necessary, the cached `raw/` tree.

---

## 6. Troubleshooting

| Symptom | Resolution |
|---------|------------|
| `manifest not found` | Pass `--manifest /path/to/cami_manifest.tsv` or run from the repo root. |
| Suite directory already exists | Use a new `--suite`/`--scenario` combination or remove the previous `run_<timestamp>/` folder. |
| Missing figures | Re-run `python bench/plot/make_figures.py --bench-root bench --tables results/.../tables/<mode> --outdir results/.../figures/<mode>`. |
| Missing contig accuracy | Only tools that emit per-contig assignments populate `contig_accuracy_per_tool.tsv` (HYMET, Kraken2, Centrifuge, ganon2, ViWrap, TAMA, SqueezeMeta, MegaPath-Nano). |

---

With this workflow, every CAMI rerun is permanently archived under `results/<scenario>/<suite>/run_<timestamp>/`, making it easy for reviewers to trace individual tools, inspect raw files, and regenerate figures without ever overwriting the canonical manuscript outputs.
