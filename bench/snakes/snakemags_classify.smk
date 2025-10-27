"""Minimal SnakeMAGs-style classification workflow for CAMI benchmarking."""

import os
from glob import glob
from pathlib import Path

SAMPLE = config["sample"]
MAG_DIR = config["mag_dir"]
WORK_DIR = config["work_dir"]
OUT_DIR = config["out_dir"]
GTDB = config["gtdb"]
THREADS = int(config.get("threads", 4))
PPLACER_THREADS = int(config.get("pplacer_threads", min(THREADS, 4)))

MAG_FASTAS = sorted(glob(os.path.join(MAG_DIR, "*.fa")))


rule all:
    input:
        os.path.join(OUT_DIR, "gtdbtk.summary.tsv")


rule run_gtdbtk:
    input:
        MAG_FASTAS
    output:
        touch(os.path.join(WORK_DIR, "gtdbtk.done"))
    threads: THREADS
    run:
        os.makedirs(Path(output[0]).parent, exist_ok=True)
        out_dir = Path(WORK_DIR) / "gtdbtk"
        out_dir.mkdir(parents=True, exist_ok=True)
        env = os.environ.copy()
        env["GTDBTK_DATA_PATH"] = GTDB
        cmd = [
            "gtdbtk",
            "classify_wf",
            # use mash if DB provides it; otherwise skip ANI screen which older
            # GTDB releases omit. We default to the safe flag here; GTDB-Tk
            # accepts it even when mash is available.
            "--skip_ani_screen",
            "--genome_dir",
            MAG_DIR,
            "--out_dir",
            str(out_dir),
            "--extension",
            "fa",
            "--cpus",
            str(threads),
            "--pplacer_cpus",
            str(PPLACER_THREADS),
        ]
        shell(" ".join(cmd), env=env)
        Path(output[0]).touch()


rule collect_summary:
    input:
        os.path.join(WORK_DIR, "gtdbtk.done")
    output:
        summary=os.path.join(OUT_DIR, "gtdbtk.summary.tsv")
    run:
        out_path = Path(output.summary)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        gtdb_dir = Path(WORK_DIR) / "gtdbtk"
        summaries = sorted(gtdb_dir.glob("gtdbtk.*.summary.tsv"))
        if not summaries:
            out_path.write_text("")
            return
        with out_path.open("w", encoding="utf-8") as dest:
            header_written = False
            for summary in summaries:
                with summary.open("r", encoding="utf-8", errors="ignore") as src:
                    for line_idx, line in enumerate(src):
                        if line_idx == 0:
                            if not header_written:
                                dest.write(line)
                                header_written = True
                        else:
                            dest.write(line)
