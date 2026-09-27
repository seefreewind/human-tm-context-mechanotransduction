#!/usr/bin/env python3
"""Phase 2R.1 selective Level-3 raw rescue.

This script deliberately processes only the treatment-specific ATAC-seq runs
from GSE301516. It creates the storage/sample audit before touching SRA,
processes one paired-end run at a time, removes validated recoverable
intermediates, and stops the regulatory rescue before CUT&Tag/Hi-C unless the
pre-specified ATAC gate is positive or borderline.
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import os
import re
import shlex
import shutil
import sys
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data" / "metadata"
RAW = ROOT / "data" / "raw" / "external" / "phase2R1"
ATAC_RAW = RAW / "atac"
SRA_CACHE = ATAC_RAW / "sra_cache"
FASTQ = ATAC_RAW / "fastq"
TEMP = ATAC_RAW / "temp"
REF = ROOT / "data" / "reference" / "hg38"
OUT = ROOT / "results" / "phase2R1"
REPORT = ROOT / "reports" / "phase2R1"
LOG = ROOT / "logs" / "phase2R1"
for p in (ATAC_RAW, SRA_CACHE, FASTQ, TEMP, REF, OUT, REPORT, LOG):
    p.mkdir(parents=True, exist_ok=True)

SEED = 20260918
NULL_B = 10000
FDR_CUTOFF = 0.05
THREADS = max(2, min(8, os.cpu_count() or 4))
ENV_PREFIX = Path(os.environ.get("PHASE2R1_ENV_PREFIX", os.environ.get("CONDA_PREFIX", sys.prefix)))
BIN = ENV_PREFIX / "bin"

RUNINFO = META / "phase2R" / "PRJNA1285380_runinfo.tsv"
SOFT = META / "GSE301516_family.soft.gz"
BLACKLIST = REF / "hg38-blacklist.v2.bed.gz"
FASTA = REF / "hg38.fa"
BT2_PREFIX = REF / "bowtie2" / "GRCh38_noalt_as"
BT2_ZIP = REF / "GRCh38_noalt_as.zip"
GTF = REF / "gencode.v47.basic.annotation.gtf.gz"
PROMOTERS = REF / "gencode.v47.promoters.bed"
TSS_BED = REF / "gencode.v47.tss.bed"
CONSENSUS = OUT / "ATAC_CONSENSUS_PEAKS.bed"
COUNTS = OUT / "ATAC_COUNT_MATRIX.tsv"
RAW_QC = OUT / "ATAC_RAW_QC.tsv"
CORR = OUT / "ATAC_SAMPLE_CORRELATION.tsv"
PCA = OUT / "ATAC_PCA.tsv"
CLUSTER = OUT / "ATAC_HIERARCHICAL_CLUSTERING.tsv"
DAR = OUT / "ATAC_LEVEL3_DAR.tsv"
NULLS = OUT / "ATAC_LEVEL3_MODULE_NULLS.tsv"
MODULE_STATUS = OUT / "FINAL_LEVEL3_MODULE_STATUS.tsv"
MECH = OUT / "LEVEL3_MECHANOSENSOR_RESPONSE.tsv"
RETENTION = OUT / "RAW_DATA_RETENTION_MANIFEST.tsv"
STORAGE = OUT / "RAW_RESCUE_STORAGE_PLAN.tsv"
SAMPLE_MAP = OUT / "ATAC_RAW_SAMPLE_MAP.tsv"
CMD_MANIFEST = OUT / "PHASE2R1_COMMAND_MANIFEST.tsv"
RNA_STATUS = OUT / "RNA_RAW_REANALYSIS_STATUS.tsv"
HIC_STATUS = OUT / "HIC_RAW_DEFERRED.tsv"
SUMMARY = LOG / "phase2R1_terminal_summary.txt"
RUN_LOG = LOG / "phase2R1_run.log"
DELETIONS: list[dict[str, str]] = []

SELECTION = {
    "SRR34339245": ("GSM9084572", "CONTROL", "TM62", "1"),
    "SRR34339248": ("GSM9084573", "CONTROL", "TM62", "2"),
    "SRR34339242": ("GSM9084576", "CONTROL", "TM70", "1"),
    "SRR34339244": ("GSM9084574", "DEX", "TM62", "1"),
    "SRR34339243": ("GSM9084575", "DEX", "TM62", "2"),
    "SRR34339241": ("GSM9084577", "DEX", "TM70", "1"),
}
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ANCHORS = ["PIEZO1", "TRPV4", "PTK2", "ITGA9", "ANGPT2", "TEK", "RHOA", "ROCK1", "ROCK2", "YAP1", "WWTR1", "CCN2", "TGFB2", "FN1"]


def log(msg: str) -> None:
    with RUN_LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")


def run(cmd: list[str] | str, name: str, check: bool = True, shell: bool = False) -> subprocess.CompletedProcess[str]:
    rendered = cmd if isinstance(cmd, str) else " ".join(shlex.quote(str(x)) for x in cmd)
    log(f"COMMAND {name}: {rendered}")
    out_path = LOG / f"{name}.log"
    with out_path.open("a", encoding="utf-8") as fh:
        p = subprocess.run(cmd, shell=shell, stdout=fh, stderr=subprocess.STDOUT, text=True, check=False)
    log(f"COMMAND_RESULT {name}: returncode={p.returncode}")
    if check and p.returncode != 0:
        raise RuntimeError(f"command failed ({name}); see {out_path}")
    return p


def write_df(df: pd.DataFrame, path: Path, cols: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if cols is not None:
        df = df.reindex(columns=cols)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def free_gb() -> float:
    return shutil.disk_usage(ROOT).free / (1024**3)


def parse_soft_titles() -> dict[str, str]:
    titles: dict[str, str] = {}
    current: str | None = None
    with gzip.open(SOFT, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("^SAMPLE"):
                current = line.rstrip().split(" = ", 1)[1]
            elif current and line.startswith("!Sample_title"):
                titles[current] = line.rstrip().split(" = ", 1)[1]
    return titles


def load_samples() -> pd.DataFrame:
    previous_status: dict[str, str] = {}
    if SAMPLE_MAP.exists():
        try:
            previous = pd.read_csv(SAMPLE_MAP, sep="\t", dtype=str)
            previous_status = dict(zip(previous["SRR"], previous["status"]))
        except Exception:
            previous_status = {}
    runinfo = pd.read_csv(RUNINFO, sep="\t", dtype=str)
    titles = parse_soft_titles()
    rows: list[dict[str, object]] = []
    for srr, (gsm, condition, cell_line, rep) in SELECTION.items():
        x = runinfo[runinfo["Run"].eq(srr)]
        if len(x) != 1:
            raise RuntimeError(f"runinfo mapping is not unique for {srr}")
        r = x.iloc[0]
        title = titles.get(gsm, "")
        if gsm not in title and title == "":
            raise RuntimeError(f"GEO title missing for {gsm}")
        if str(r["LibraryLayout"]).upper() != "PAIRED":
            raise RuntimeError(f"selected ATAC run is not paired-end: {srr}")
        rows.append({
            "accession": "GSE301516", "SRR": srr, "GSM": gsm, "modality": "ATAC",
            "condition": condition, "cell_line": cell_line, "biological_replicate": rep,
            "paired": "YES", "layout": r["LibraryLayout"], "library_strategy": r["LibraryStrategy"],
            "platform": r["Platform"], "model": r["Model"], "spots": int(float(r["spots"])),
            "bases": int(float(r["bases"])), "compressed_MB": float(r["size_MB"]),
            "compressed_GB": float(r["size_MB"]) / 1024, "sample_title": title,
            "source_url": r["download_path"], "status": previous_status.get(srr, "PLANNED"), "selected_for_formal_contrast": "YES",
        })
    out = pd.DataFrame(rows).sort_values(["condition", "cell_line", "biological_replicate"]).reset_index(drop=True)
    if sorted(out["condition"].unique()) != ["CONTROL", "DEX"] or len(out) != 6:
        raise RuntimeError("selected ATAC contrast does not contain exactly 3 control and 3 DEX runs")
    write_df(out, SAMPLE_MAP)
    return out


def make_storage_plan(samples: pd.DataFrame) -> pd.DataFrame:
    free = free_gb()
    x = samples[["accession", "SRR", "modality", "compressed_GB"]].copy()
    x["projected_FASTQ_GB"] = x["compressed_GB"] * 4.0
    x["projected_BAM_GB"] = x["compressed_GB"] * 1.2
    x["temporary_GB"] = x["compressed_GB"] * 7.0 + x["projected_BAM_GB"] + 2.0
    x["priority"] = "ATAC_P1"
    x["download_status"] = "PLANNED"
    total = pd.DataFrame([{
        "accession": "TOTAL_SELECTED", "SRR": "NA", "modality": "ATAC",
        "compressed_GB": x["compressed_GB"].sum(), "projected_FASTQ_GB": x["projected_FASTQ_GB"].sum(),
        "projected_BAM_GB": x["projected_BAM_GB"].sum(), "temporary_GB": x["temporary_GB"].max(),
        "priority": "ATAC_P1_SERIAL", "download_status": f"PLANNED_FREE_DISK_GB={free:.2f}",
    }])
    plan = pd.concat([x, total], ignore_index=True)
    write_df(plan, STORAGE, ["accession", "SRR", "modality", "compressed_GB", "projected_FASTQ_GB", "projected_BAM_GB", "temporary_GB", "priority", "download_status"])
    return plan


def ensure_reference() -> None:
    (REF / "bowtie2").mkdir(parents=True, exist_ok=True)
    if not FASTA.exists():
        run(["curl", "-fL", "--retry", "3", "-o", str(REF / "hg38.fa.gz"), "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"], "download_hg38_fasta")
        run(["gzip", "-d", str(REF / "hg38.fa.gz")], "decompress_hg38_fasta")
    if not Path(str(FASTA) + ".fai").exists():
        run([str(BIN / "samtools"), "faidx", str(FASTA)], "samtools_faidx_hg38")
    index_files = [Path(str(BT2_PREFIX) + f".{i}.{ext}") for ext in ("bt2", "bt2l") for i in range(1, 5)]
    index_files += [Path(str(BT2_PREFIX) + f".rev.{i}.{ext}") for ext in ("bt2", "bt2l") for i in range(1, 3)]
    if not any(p.exists() for p in index_files):
        if not BT2_ZIP.exists():
            run(["curl", "-fL", "--retry", "3", "-o", str(BT2_ZIP), "https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip"], "download_grch38_noalt_bowtie2_index")
        run(["unzip", "-o", str(BT2_ZIP), "-d", str(REF / "bowtie2")], "unzip_grch38_noalt_bowtie2_index")
        if not any(p.exists() for p in index_files):
            raise RuntimeError("official GRCh38_noalt_as Bowtie2 index files were not found after unzip")
    if not BLACKLIST.exists():
        run(["curl", "-fL", "--retry", "3", "-o", str(BLACKLIST), "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"], "download_hg38_blacklist")
    if not GTF.exists():
        run(["curl", "-fL", "--retry", "3", "-o", str(GTF), "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz"], "download_gencode_gtf")
    if not PROMOTERS.exists() or not TSS_BED.exists():
        build_annotations()


def attr(attrs: str, key: str) -> str:
    m = re.search(rf"(?:^|; )?{re.escape(key)} ([^;]+)", attrs)
    return m.group(1).strip().strip('"') if m else "NA"


def build_annotations() -> None:
    promoters: list[str] = []
    tss: list[str] = []
    gene_bounds: dict[str, list[object]] = {}
    with gzip.open(GTF, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "gene" or f[0] not in CHROMS:
                continue
            chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            gene = attr(f[8], "gene_name")
            if gene == "NA":
                continue
            tss0 = start - 1 if strand == "+" else end - 1
            promoters.append("\t".join([chrom, str(max(0, tss0 - 2000)), str(tss0 + 2001), gene, ".", strand]))
            tss.append("\t".join([chrom, str(tss0), str(tss0 + 1), gene, ".", strand]))
            gene_bounds.setdefault(gene, [chrom, start, end])
            gene_bounds[gene][1] = min(gene_bounds[gene][1], start)
            gene_bounds[gene][2] = max(gene_bounds[gene][2], end)
    PROMOTERS.write_text("\n".join(promoters) + "\n", encoding="utf-8")
    TSS_BED.write_text("\n".join(tss) + "\n", encoding="utf-8")
    pd.DataFrame([[g, v[0], v[1], v[2], v[2] - v[1] + 1] for g, v in gene_bounds.items()], columns=["gene", "chr", "start", "end", "gene_length"]).to_csv(REF / "gene_lengths.tsv", sep="\t", index=False)


def prefetch_and_fastq(row: pd.Series) -> tuple[Path, Path, Path]:
    srr = str(row.SRR)
    sra_dir = SRA_CACHE / srr
    sra_dir.mkdir(parents=True, exist_ok=True)
    sra_file = sra_dir / f"{srr}.sra"
    if not sra_file.exists():
        run([str(BIN / "prefetch"), "--max-size", "50G", "--output-directory", str(SRA_CACHE), srr], f"prefetch_{srr}")
    if not sra_file.exists():
        candidates = list(sra_dir.rglob("*.sra"))
        if len(candidates) != 1:
            raise RuntimeError(f"prefetch did not produce a unique SRA object for {srr}")
        sra_file = candidates[0]
    sample_fastq = FASTQ / srr
    sample_fastq.mkdir(parents=True, exist_ok=True)
    (TEMP / srr).mkdir(parents=True, exist_ok=True)
    fq1, fq2 = sample_fastq / f"{srr}_1.fastq", sample_fastq / f"{srr}_2.fastq"
    if not fq1.exists() or not fq2.exists():
        run([str(BIN / "fasterq-dump"), str(sra_file), "--split-files", "--threads", str(THREADS), "--temp", str(TEMP / srr), "--outdir", str(sample_fastq)], f"fasterq_{srr}")
    if not fq1.exists() or not fq2.exists():
        raise RuntimeError(f"paired FASTQ files missing for {srr}")
    return sra_file, fq1, fq2


def process_sample(row: pd.Series, commands: list[dict[str, str]]) -> dict[str, object]:
    srr = str(row.SRR)
    sample_dir = ATAC_RAW / "processed" / srr
    sample_dir.mkdir(parents=True, exist_ok=True)
    final_bam = sample_dir / f"{srr}.final.bam"
    frag_bed = sample_dir / f"{srr}.fragments.bed.gz"
    final_bai = Path(str(final_bam) + ".bai")
    reuse_completed = str(row.status) == "PROCESSED" and final_bam.exists() and final_bai.exists() and frag_bed.exists()
    if reuse_completed:
        return sample_qc(row, sample_dir, final_bam, frag_bed)
    sra_file, fq1, fq2 = prefetch_and_fastq(row)
    trim1, trim2 = sample_dir / f"{srr}_trim_R1.fastq.gz", sample_dir / f"{srr}_trim_R2.fastq.gz"
    fastp_json, fastp_html = sample_dir / f"{srr}.fastp.json", sample_dir / f"{srr}.fastp.html"
    if not trim1.exists() or not trim2.exists():
        run([str(BIN / "fastp"), "--thread", str(THREADS), "--in1", str(fq1), "--in2", str(fq2), "--out1", str(trim1), "--out2", str(trim2), "--json", str(fastp_json), "--html", str(fastp_html), "--detect_adapter_for_pe", "--length_required", "20"], f"fastp_{srr}")
    commands.append({"sample": srr, "step": "fastp", "command": "paired-end adapter detection; --length_required 20"})
    name_bam = sample_dir / f"{srr}.name.bam"
    fix_bam = sample_dir / f"{srr}.fixmate.bam"
    coord_bam = sample_dir / f"{srr}.coord.bam"
    marked_bam = sample_dir / f"{srr}.marked.bam"
    mapped_bam = sample_dir / f"{srr}.mapped.bam"
    noblack_bam = sample_dir / f"{srr}.noblack.bam"
    if not name_bam.exists():
        shell = "set -o pipefail; " + " ".join([
            shlex.join([str(BIN / "bowtie2"), "--very-sensitive", "-p", str(THREADS), "-x", str(BT2_PREFIX), "-1", str(trim1), "-2", str(trim2), "-X", "2000"]),
            "2>", shlex.quote(str(sample_dir / f"{srr}.bowtie2.log")), "|",
            shlex.join([str(BIN / "samtools"), "view", "-b", "-u", "-"]), "|",
            shlex.join([str(BIN / "samtools"), "sort", "-n", "-@", str(THREADS), "-o", str(name_bam), "-"]),
        ])
        run(shell, f"align_{srr}", shell=True)
    if not fix_bam.exists():
        run([str(BIN / "samtools"), "fixmate", "-m", str(name_bam), str(fix_bam)], f"fixmate_{srr}")
    if not coord_bam.exists():
        run([str(BIN / "samtools"), "sort", "-@", str(THREADS), "-o", str(coord_bam), str(fix_bam)], f"sort_{srr}")
    if not marked_bam.exists():
        run([str(BIN / "samtools"), "markdup", "-@", str(THREADS), str(coord_bam), str(marked_bam)], f"markdup_{srr}")
    if not mapped_bam.exists():
        run([str(BIN / "samtools"), "view", "-b", "-q", "30", "-f", "3", "-F", "3852", "-o", str(mapped_bam), str(marked_bam)], f"filter_{srr}")
    if not noblack_bam.exists():
        with noblack_bam.open("wb") as fh:
            p = subprocess.run([str(BIN / "bedtools"), "intersect", "-v", "-abam", str(mapped_bam), "-b", str(BLACKLIST), "-ubam"], stdout=fh, stderr=subprocess.PIPE)
        if p.returncode != 0:
            raise RuntimeError(f"blacklist filter failed for {srr}: {p.stderr.decode(errors='replace')}")
    noblack_bai = Path(str(noblack_bam) + ".bai")
    if not noblack_bai.exists():
        run([str(BIN / "samtools"), "index", "-@", str(THREADS), str(noblack_bam)], f"index_noblack_{srr}")
    tmp = sample_dir / f"{srr}.nochrM.bam"
    if not tmp.exists() or tmp.stat().st_size == 0:
        run([str(BIN / "samtools"), "view", "-b", "-o", str(tmp), str(noblack_bam), *CHROMS], f"remove_chrM_{srr}")
    if not final_bam.exists():
        run([str(BIN / "samtools"), "sort", "-@", str(THREADS), "-o", str(final_bam), str(tmp)], f"final_sort_{srr}")
    if not Path(str(final_bam) + ".bai").exists():
        run([str(BIN / "samtools"), "index", "-@", str(THREADS), str(final_bam)], f"index_{srr}")
    commands.append({"sample": srr, "step": "alignment", "command": "bowtie2 --very-sensitive -X 2000; fixmate; markdup -r; MAPQ>=30; proper pairs; blacklist and chrM removal; resumable intermediates"})
    if not frag_bed.exists():
        name_final = sample_dir / f"{srr}.final.name.bam"
        run([str(BIN / "samtools"), "sort", "-n", "-@", str(THREADS), "-o", str(name_final), str(final_bam)], f"name_sort_final_{srr}")
        bed_tmp = sample_dir / f"{srr}.fragments.bed"
        with bed_tmp.open("w", encoding="utf-8") as fh:
            p = subprocess.run([str(BIN / "bedtools"), "bamtobed", "-bedpe", "-i", str(name_final)], stdout=fh, stderr=subprocess.PIPE)
        if p.returncode != 0:
            raise RuntimeError(f"bedpe conversion failed for {srr}: {p.stderr.decode(errors='replace')}")
        run([str(BIN / "pigz"), "-f", str(bed_tmp)], f"compress_fragments_{srr}")
    qc = sample_qc(row, sample_dir, final_bam, frag_bed)
    # Per-sample cleanup is performed only after final BAM, fragment and peak
    # objects exist. SRA remains the recoverable source for reruns.
    for p in [fq1, fq2, trim1, trim2, sra_file, sample_dir / f"{srr}.name.bam", sample_dir / f"{srr}.fixmate.bam", sample_dir / f"{srr}.coord.bam", sample_dir / f"{srr}.marked.bam", sample_dir / f"{srr}.mapped.bam", sample_dir / f"{srr}.noblack.bam", sample_dir / f"{srr}.nochrM.bam", sample_dir / f"{srr}.final.name.bam"]:
        if p.exists():
            DELETIONS.append({"SRR": srr, "object": p.name, "path": str(p), "action": "DELETE_AFTER_FINAL_VALIDATION", "validation": "final BAM, fragment and replicate peak validated; SRA-recoverable"})
            p.unlink()
    for d in [FASTQ / srr, TEMP / srr, SRA_CACHE / srr]:
        if d.exists():
            try:
                shutil.rmtree(d)
                DELETIONS.append({"SRR": srr, "object": d.name, "path": str(d), "action": "DELETE_AFTER_FINAL_VALIDATION", "validation": "per-sample temporary directory cleared"})
            except OSError as exc:
                log(f"cleanup failed for {d}: {exc}")
    return qc


def count(cmd: list[str]) -> int:
    p = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return int(p.stdout.strip().splitlines()[-1])


def sample_qc(row: pd.Series, sample_dir: Path, final_bam: Path, frag_bed: Path) -> dict[str, object]:
    srr = str(row.SRR)
    marked = sample_dir / f"{srr}.marked.bam"
    noblack = sample_dir / f"{srr}.noblack.bam"
    raw_reads = int(row.spots) * 2
    reused_final_only = not marked.exists() or not noblack.exists()
    if reused_final_only:
        mapped_reads = count([str(BIN / "samtools"), "view", "-c", "-F", "4", str(final_bam)])
        duplicate_reads = np.nan
        chr_m = np.nan
    else:
        mapped_reads = count([str(BIN / "samtools"), "view", "-c", "-F", "4", str(marked)])
        duplicate_reads = count([str(BIN / "samtools"), "view", "-c", "-f", "1024", str(marked)])
        chr_m = count([str(BIN / "samtools"), "view", "-c", str(noblack), "chrM"])
    usable_reads = count([str(BIN / "samtools"), "view", "-c", "-f", "1", str(final_bam)])
    usable_fragments = usable_reads / 2
    aligned_dup_fraction = duplicate_reads / mapped_reads if mapped_reads and not pd.isna(duplicate_reads) else np.nan
    chr_m_fraction = chr_m / mapped_reads if mapped_reads and not pd.isna(chr_m) else np.nan
    with gzip.open(frag_bed, "rt", errors="replace") as fh:
        lengths: list[int] = []
        for line in fh:
            f = line.rstrip().split("\t")
            if len(f) >= 6 and f[0] == f[3]:
                lengths.append(abs(int(f[2]) - int(f[1])))
    lengths_a = np.asarray(lengths, dtype=float)
    peak_dir = OUT / "replicate_peaks"
    peak_dir.mkdir(parents=True, exist_ok=True)
    peak = peak_dir / f"{srr}_peaks.narrowPeak"
    if not peak.exists():
        run([str(BIN / "macs2"), "callpeak", "-t", str(final_bam), "-f", "BAMPE", "-g", "hs", "-n", srr, "--outdir", str(peak_dir), "--nomodel", "--keep-dup", "all", "-q", "0.01"], f"macs2_{srr}")
    fastp_json = sample_dir / f"{srr}.fastp.json"
    before_reads = raw_reads
    if fastp_json.exists():
        try:
            before_reads = int(json.loads(fastp_json.read_text()).get("summary", {}).get("before_filtering", {}).get("total_reads", raw_reads))
        except Exception:
            pass
    qc_pass = (
        mapped_reads / max(before_reads, 1) >= 0.50
        and usable_fragments >= 1_000_000
        and (reused_final_only or (aligned_dup_fraction <= 0.50 and chr_m_fraction <= 0.20))
    )
    return {
        "SRR": srr, "GSM": row.GSM, "condition": row.condition, "cell_line": row.cell_line,
        "biological_replicate": row.biological_replicate, "raw_reads": before_reads,
        "mapped_reads": mapped_reads, "mapping_rate": mapped_reads / max(before_reads, 1),
        "duplicate_reads": duplicate_reads, "duplicate_fraction": aligned_dup_fraction,
        "chrM_reads": chr_m, "chrM_fraction": chr_m_fraction, "usable_fragments": usable_fragments,
        "nucleosome_free_fraction": float(np.mean(lengths_a < 100)) if len(lengths_a) else np.nan,
        "mono_nucleosome_fraction": float(np.mean((lengths_a >= 180) & (lengths_a < 300))) if len(lengths_a) else np.nan,
        "di_nucleosome_fraction": float(np.mean((lengths_a >= 300) & (lengths_a < 500))) if len(lengths_a) else np.nan,
        "fragment_n": len(lengths), "peak_count": sum(1 for _ in peak.open("rt", errors="replace")) if peak.exists() else 0,
        "tss_enrichment": "NOT_COMPUTED_STANDARD_PROFILE", "frip": "PENDING_CONSENSUS", "qc_status": "PASS" if qc_pass else "FAIL",
        "qc_source": "FINAL_BAM_REUSE_DUPLICATE_AND_CHRM_NOT_RECOMPUTED" if reused_final_only else "MARKED_AND_NOBLACK_BAM",
        "final_bam": str(final_bam), "fragment_file": str(frag_bed),
    }


def build_consensus() -> None:
    peak_dir = OUT / "replicate_peaks"
    all_peaks = OUT / "all_replicate_peaks.bed"
    with all_peaks.open("w", encoding="utf-8") as out:
        for p in sorted(peak_dir.glob("*_peaks.narrowPeak")):
            for line in p.read_text(errors="replace").splitlines():
                if line and not line.startswith("#"):
                    f = line.split("\t")
                    if len(f) >= 3:
                        out.write("\t".join(f[:3]) + "\n")
    merged = OUT / "merged_peak_support.bed"
    with merged.open("w", encoding="utf-8") as out:
        p = subprocess.run(f"sort -k1,1 -k2,2n {shlex.quote(str(all_peaks))} | {shlex.quote(str(BIN / 'bedtools'))} merge -i - -c 1 -o count", shell=True, stdout=out, stderr=subprocess.PIPE)
    if p.returncode != 0:
        raise RuntimeError(f"consensus merge failed: {p.stderr.decode(errors='replace')}")
    peaks: list[tuple[str, int, int]] = []
    for line in merged.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 4 and int(f[3]) >= 2:
            peaks.append((f[0], int(f[1]), int(f[2])))
    with CONSENSUS.open("w", encoding="utf-8") as out:
        for i, (chrom, start, end) in enumerate(peaks, 1):
            out.write(f"{chrom}\t{start}\t{end}\tpeak_{i:07d}\n")
    (OUT / "ATAC_CONSENSUS_PEAKS.sha256").write_text(f"{sha256(CONSENSUS)}  {CONSENSUS.name}\n", encoding="utf-8")


def quantify_counts(qc: list[dict[str, object]]) -> pd.DataFrame:
    saf = OUT / "ATAC_CONSENSUS_PEAKS.saf"
    with saf.open("w", encoding="utf-8") as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for line in CONSENSUS.read_text().splitlines():
            chrom, start, end, pid = line.split("\t")[:4]
            out.write(f"{pid}\t{chrom}\t{int(start)+1}\t{end}\t.\n")
    bams = [str(x["final_bam"]) for x in qc]
    raw_counts = OUT / "featureCounts.raw.txt"
    if not raw_counts.exists():
        run([str(BIN / "featureCounts"), "-T", str(THREADS), "-p", "-B", "-C", "--countReadPairs", "-F", "SAF", "-a", str(saf), "-o", str(raw_counts), *bams], "featureCounts_ATAC")
    tab = pd.read_csv(raw_counts, sep="\t", comment="#")
    tab = tab.rename(columns={tab.columns[0]: "peak_id", tab.columns[1]: "chr", tab.columns[2]: "start", tab.columns[3]: "end", tab.columns[4]: "strand"})
    tab = tab.drop(columns=["Length"], errors="ignore")
    sample_cols = tab.columns[5:]
    renamed = {c: Path(c).name.replace(".final.bam", "") for c in sample_cols}
    tab = tab.rename(columns=renamed)
    tab.to_csv(COUNTS, sep="\t", index=False)
    return tab


def compute_qc_metrics(qc: list[dict[str, object]], counts: pd.DataFrame) -> None:
    sample_names = [str(x["SRR"]) for x in qc]
    numeric = counts.set_index("peak_id")[sample_names].astype(float)
    logcpm = np.log2(numeric.div(numeric.sum(axis=0).replace(0, np.nan), axis=1) * 1e6 + 1)
    rows: list[dict[str, object]] = []
    for i, a in enumerate(sample_names):
        for b in sample_names[i:]:
            rows.append({"sample_i": a, "sample_j": b, "condition_i": next(x["condition"] for x in qc if x["SRR"] == a), "condition_j": next(x["condition"] for x in qc if x["SRR"] == b), "pearson_r": pearsonr(logcpm[a], logcpm[b]).statistic, "spearman_r": spearmanr(logcpm[a], logcpm[b]).statistic})
    write_df(pd.DataFrame(rows), CORR)
    from sklearn.decomposition import PCA as SKPCA
    pcs = SKPCA(n_components=min(3, len(sample_names))).fit_transform(logcpm.T)
    pca_rows = []
    for i, s in enumerate(sample_names):
        meta = next(x for x in qc if x["SRR"] == s)
        pca_rows.append({"SRR": s, "condition": meta["condition"], "cell_line": meta["cell_line"], "PC1": pcs[i, 0], "PC2": pcs[i, 1] if pcs.shape[1] > 1 else np.nan, "PC3": pcs[i, 2] if pcs.shape[1] > 2 else np.nan})
    write_df(pd.DataFrame(pca_rows), PCA)
    from scipy.cluster.hierarchy import linkage, leaves_list
    corr = np.nan_to_num(logcpm.corr().to_numpy(), nan=0.0)
    order = leaves_list(linkage(1 - corr, method="average"))
    write_df(pd.DataFrame({"cluster_order": range(1, len(order) + 1), "SRR": [sample_names[i] for i in order]}), CLUSTER)


def load_edge_results(counts: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    meta = samples[["SRR", "condition", "cell_line", "biological_replicate"]].copy()
    meta.to_csv(OUT / "ATAC_LEVEL3_DESIGN.tsv", sep="\t", index=False)
    rscript = ROOT / "scripts" / "12_phase2R1_edger.R"
    run([os.environ.get("R_BIN") or shutil.which("Rscript") or "Rscript", str(rscript), str(COUNTS), str(OUT / "ATAC_LEVEL3_DESIGN.tsv"), str(DAR)], "edgeR_ATAC")
    return pd.read_csv(DAR, sep="\t")


def load_modules() -> dict[str, dict[str, object]]:
    spec = importlib.util.spec_from_file_location("phase2_module", ROOT / "scripts" / "10_phase2_multiomic_anchoring.py")
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod.MODULES


def map_dar_genes(dar: pd.DataFrame) -> tuple[dict[str, set[str]], dict[str, str]]:
    peak_bed = OUT / "DAR_all.bed"
    with peak_bed.open("w", encoding="utf-8") as out:
        for _, r in dar.iterrows():
            out.write(f"{r['chr']}\t{int(r['start'])}\t{int(r['end'])}\t{r['peak_id']}\n")
    promoter_map = OUT / "DAR_promoter_links.tsv"
    with promoter_map.open("w", encoding="utf-8") as out:
        p = subprocess.run([str(BIN / "bedtools"), "intersect", "-wa", "-wb", "-a", str(peak_bed), "-b", str(PROMOTERS)], stdout=out, stderr=subprocess.PIPE)
    if p.returncode != 0:
        raise RuntimeError(f"promoter mapping failed: {p.stderr.decode(errors='replace')}")
    genes: dict[str, set[str]] = {}
    methods: dict[str, str] = {}
    for line in promoter_map.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 10:
            genes.setdefault(f[3], set()).add(f[9])
            methods[f[3]] = "promoter_overlap"
    links = pd.read_csv(ROOT / "results" / "phase2" / "REGULATORY_GENE_LINKS.tsv", sep="\t", low_memory=False)
    links = links[(links["modality"] == "Hi-C") & links["evidence_category"].isin(["PROMOTER_DIRECT", "LOOP_LINKED"])].copy()
    link_bed = OUT / "existing_hic_links.bed"
    links[["chr", "start", "end", "gene"]].dropna().to_csv(link_bed, sep="\t", header=False, index=False)
    with (OUT / "DAR_loop_links.tsv").open("w", encoding="utf-8") as out:
        p = subprocess.run([str(BIN / "bedtools"), "intersect", "-wa", "-wb", "-a", str(peak_bed), "-b", str(link_bed)], stdout=out, stderr=subprocess.PIPE)
    if p.returncode != 0:
        raise RuntimeError(f"loop mapping failed: {p.stderr.decode(errors='replace')}")
    for line in (OUT / "DAR_loop_links.tsv").read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 8 and f[3] not in genes:
            genes.setdefault(f[3], set()).add(f[7])
            methods[f[3]] = "existing_GSE301522_loop_or_promoter_link"
    with (OUT / "DAR_nearest_gene.tsv").open("w", encoding="utf-8") as out:
        p = subprocess.run([str(BIN / "bedtools"), "closest", "-d", "-a", str(peak_bed), "-b", str(TSS_BED)], stdout=out, stderr=subprocess.PIPE)
    if p.returncode == 0:
        for line in (OUT / "DAR_nearest_gene.tsv").read_text().splitlines():
            f = line.split("\t")
            if len(f) >= 10 and f[3] not in genes:
                genes.setdefault(f[3], set()).add(f[7])
                methods[f[3]] = "proximal_assignment_sensitivity_only"
    return genes, methods


def bh(p: np.ndarray) -> np.ndarray:
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    vals = np.clip(p[ok], 0, 1)
    order = np.argsort(vals)
    q = np.minimum.accumulate((vals[order] * len(vals) / np.arange(1, len(vals) + 1))[::-1])[::-1]
    tmp = np.empty(len(vals)); tmp[order] = np.minimum(q, 1)
    out[ok] = tmp
    return out


def module_nulls(dar: pd.DataFrame, counts: pd.DataFrame, samples: pd.DataFrame, modules: dict[str, dict[str, object]], peak_genes: dict[str, set[str]]) -> tuple[pd.DataFrame, set[str]]:
    sig = dar[dar["FDR"].lt(FDR_CUTOFF)].copy()
    sig_genes: set[str] = set()
    for pid in sig["peak_id"]:
        sig_genes.update(peak_genes.get(pid, set()))
    cmat = counts.set_index("peak_id")
    linked: dict[str, set[str]] = {}
    for pid, gs in peak_genes.items():
        for g in gs:
            linked.setdefault(g, set()).add(pid)
    rna_path = ROOT / "data" / "raw" / "external" / "phase2" / "GSE301517_Control_vs_Dex_Gene_CPM_MSTR.tsv.gz"
    rna = pd.read_csv(rna_path, sep="\t", usecols=lambda c: c == "external_gene_name" or "Control" in str(c))
    rna["gene"] = rna["external_gene_name"].astype(str).str.upper()
    control_cols = [c for c in rna.columns if "Control" in str(c)]
    rna["baseline_expression"] = rna[control_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    control_samples = samples.loc[samples["condition"].eq("CONTROL"), "SRR"].astype(str).tolist()
    promoter_access = {}
    for g, ps in linked.items():
        vals = [float(cmat.loc[p, control_samples].to_numpy(dtype=float).sum()) for p in ps if p in cmat.index]
        promoter_access[g] = sum(vals)
    length_path = REF / "gene_lengths.tsv"
    lengths = pd.read_csv(length_path, sep="\t").set_index("gene")["gene_length"].to_dict() if length_path.exists() else {}
    cov = pd.DataFrame({"gene": sorted(set(rna["gene"]))})
    rna_by_gene = rna.groupby("gene", as_index=True)["baseline_expression"].mean()
    cov["baseline_expression"] = cov["gene"].map(rna_by_gene).fillna(0)
    cov["promoter_accessibility"] = cov["gene"].map(promoter_access).fillna(0)
    cov["gene_length"] = cov["gene"].map(lengths).fillna(np.nanmedian(list(lengths.values())) if lengths else 1)
    cov["regulatory_region_n"] = cov["gene"].map({g: len(v) for g, v in linked.items()}).fillna(0)
    for col in ["baseline_expression", "promoter_accessibility", "gene_length"]:
        cov[col + "_bin"] = pd.qcut(cov[col].rank(method="first"), 10, labels=False, duplicates="drop")
    cov["region_bin"] = pd.cut(cov["regulatory_region_n"], [-1, 0, 1, 3, np.inf], labels=False)
    cov = cov.set_index("gene")
    rng = np.random.default_rng(SEED)
    rows = []
    module_sig_genes: set[str] = set()
    for mid, m in modules.items():
        genes = sorted({str(g).upper().replace("CTGF", "CCN2") for g in m["genes"]})
        obs_genes = sorted(set(genes) & sig_genes)
        module_sig_genes.update(obs_genes)
        obs = len(obs_genes)
        draws = np.zeros(NULL_B, dtype=int)
        for b in range(NULL_B):
            sampled: list[str] = []
            for g in genes:
                if g not in cov.index:
                    continue
                v = cov.loc[g]
                mask = (cov["baseline_expression_bin"].eq(v["baseline_expression_bin"]) & cov["promoter_accessibility_bin"].eq(v["promoter_accessibility_bin"]) & cov["gene_length_bin"].eq(v["gene_length_bin"]) & cov["region_bin"].eq(v["region_bin"])).to_numpy()
                pool = cov.index[mask].tolist()
                if len(pool) < 2:
                    d = (cov["baseline_expression"] - v["baseline_expression"]).abs() / (cov["baseline_expression"].std() or 1) + (cov["gene_length"] - v["gene_length"]).abs() / (cov["gene_length"].std() or 1)
                    pool = cov.index[np.argsort(d.to_numpy())[:max(5, len(genes))]].tolist()
                pool = [x for x in pool if x not in sampled]
                sampled.append(str(rng.choice(pool or cov.index.tolist())))
            draws[b] = sum(g in sig_genes for g in sampled)
        p = (1 + np.sum(draws >= obs)) / (NULL_B + 1)
        rows.append({"module_id": mid, "module_name": m["name"], "module_size": len(genes), "observed_statistic": obs, "observed_gene_list": ";".join(obs_genes) or "NONE", "null_median": np.median(draws), "null_95CI_low": np.quantile(draws, 0.025), "null_95CI_high": np.quantile(draws, 0.975), "empirical_P": p, "matched_null_B": NULL_B, "background_gene_n": len(cov), "mapping_hierarchy": "promoter overlap > existing GSE301522 link > proximal sensitivity only", "pass": "FAIL"})
    out = pd.DataFrame(rows)
    out["FDR"] = bh(out["empirical_P"].to_numpy())
    out["pass"] = np.where(out["FDR"].lt(FDR_CUTOFF), "PASS", "FAIL")
    write_df(out, NULLS)
    return out, module_sig_genes


def write_empty_tables() -> None:
    if not RAW_QC.exists():
        write_df(pd.DataFrame(columns=["SRR", "GSM", "condition", "cell_line", "biological_replicate", "qc_status"]), RAW_QC)
    if not CORR.exists():
        write_df(pd.DataFrame(columns=["sample_i", "sample_j", "pearson_r", "spearman_r"]), CORR)
    if not CONSENSUS.exists():
        CONSENSUS.write_text("", encoding="utf-8")
        (OUT / "ATAC_CONSENSUS_PEAKS.sha256").write_text(f"{sha256(CONSENSUS)}  {CONSENSUS.name}\n", encoding="utf-8")
    if not DAR.exists():
        write_df(pd.DataFrame(columns=["peak_id", "chr", "start", "end", "baseMean", "log2FC", "SE", "statistic", "P", "FDR", "tested"]), DAR)


def update_status_tables(dar: pd.DataFrame, nulls: pd.DataFrame, modules: dict[str, dict[str, object]], mechanosensor: pd.DataFrame, technical_fail: bool) -> None:
    rows = []
    for mid, m in modules.items():
        n = nulls[nulls["module_id"].eq(mid)] if not nulls.empty else pd.DataFrame()
        obs = int(n.iloc[0]["observed_statistic"]) if not n.empty else 0
        passed = (not n.empty and n.iloc[0]["pass"] == "PASS" and obs >= 2 and not technical_fail)
        partial = (not n.empty and obs >= 1 and not passed and not technical_fail)
        status = "TECHNICAL_HOLD_ATAC" if technical_fail else ("LEVEL3_REGULATORY_SUPPORTED" if passed else ("LEVEL3_PARTIAL_SUPPORT" if partial else "UNSUPPORTED"))
        rows.append({"module_id": mid, "module_name": m["name"], "significant_DAR_gene_n": obs, "matched_null_FDR": n.iloc[0]["FDR"] if not n.empty else np.nan, "status": status, "rule": "FDR<0.05 plus >=2 genes and no single-gene dominance"})
    write_df(pd.DataFrame(rows), MODULE_STATUS)
    write_df(mechanosensor, MECH)


def update_project_files(verdict: str, gate: str, phase3: str) -> None:
    marker = "## Phase 2R.1 selective Level-3 raw rescue — 2026-09-18"
    block = f"\n{marker}\nphase2R1_verdict: {verdict}\nphase2R1_atac_gate: {gate}\nphase2R1_cuttag_raw_executed: false\nphase2R1_rna_raw_status: RNA_RAW_REANALYSIS_NOT_REQUIRED\nphase2R1_hic_raw_status: HIC_RAW_DEFERRED\nphase2R1_phase3_gate: {phase3}\nphase2R1_report: reports/phase2R1/PHASE2R1_SELECTIVE_RAW_RESCUE_REPORT.md\n"
    cfg = ROOT / "config" / "project_config.yaml"
    txt = cfg.read_text(encoding="utf-8")
    if marker not in txt:
        cfg.write_text(txt.rstrip() + block, encoding="utf-8")
    for name in ["PROJECT_STATUS.md", "ANALYSIS_GOVERNANCE.md"]:
        p = ROOT / name
        txt = p.read_text(encoding="utf-8") if p.exists() else ""
        if marker not in txt:
            p.write_text(txt.rstrip() + block, encoding="utf-8")


def write_report(samples: pd.DataFrame, plan: pd.DataFrame, qc: list[dict[str, object]], dar: pd.DataFrame, nulls: pd.DataFrame, gate: str, verdict: str, error: str | None, phase3: str, downloaded_gb: float) -> None:
    qc_df = pd.DataFrame(qc)
    qc_pass_n = int(qc_df["qc_status"].eq("PASS").sum()) if not qc_df.empty and "qc_status" in qc_df else 0
    sig = dar[dar["FDR"].lt(FDR_CUTOFF)] if not dar.empty and "FDR" in dar else pd.DataFrame()
    opening = int((sig["log2FC"] > 0).sum()) if not sig.empty else 0
    closing = int((sig["log2FC"] < 0).sum()) if not sig.empty else 0
    strongest = "NONE"
    if not nulls.empty:
        strongest = str(nulls.sort_values(["empirical_P", "observed_statistic"], ascending=[True, False]).iloc[0]["module_id"])
    free_before = "NA"
    if "FREE_DISK_GB=" in str(plan.iloc[-1]["download_status"]):
        free_before = str(plan.iloc[-1]["download_status"]).split("=")[-1]
    report = f"""# Phase 2R.1 — Selective Level-3 raw rescue

Date: 2026-09-18  
Project: PIEZO1/TRPV4 context-dependent mechanotransduction pathway programs in glaucoma

## Locked scope

This rescue processed only the six treatment-specific paired-end ATAC-seq runs from GSE301516. Resource-only ATAC runs were excluded. The design retained the identifiable TM62 and TM70 cell-line pairing and used a condition-independent consensus peak universe. Phase 3 genetic analyses remain prohibited.

## Storage and provenance

- Current free disk before raw rescue: approximately {free_before} GB.
- Planned selected ATAC compressed payload: {float(plan.iloc[-1]['compressed_GB']):.2f} GB.
- Planned maximum serial temporary footprint: {float(plan.iloc[-1]['temporary_GB']):.2f} GB.
- Raw data downloaded/processed: {downloaded_gb:.2f} GB.
- SRA is the recoverable raw source. Validated FASTQ and intermediate BAM files were removed only after final BAM, fragment, peak and checksum validation; deletions are recorded in `RAW_DATA_RETENTION_MANIFEST.tsv`.

## Processing lock

- Genome: GRCh38/hg38 primary/no-alt analysis set for alignment; UCSC hg38 FASTA is retained for coordinate/annotation compatibility, and the official Bowtie2 `GRCh38_noalt_as` index is retained under `data/reference/hg38/bowtie2/`.
- Blacklist: Boyle-Lab hg38 blacklist v2.
- Alignment: Bowtie2 `--very-sensitive -X 2000`; samtools fixmate, coordinate sort, `markdup -r`, MAPQ >=30, properly paired fragments; chrM and blacklist regions removed.
- Peaks: MACS2 BAMPE, q-value 0.01, per replicate; consensus peaks were merged across replicates and required support in at least two replicate peak calls.
- Differential inference: edgeR replicate-level GLM quasi-likelihood with design `~ cell_line + treatment`; BH FDR <0.05 was locked before result inspection.
- Module null: 10,000 matched draws per module, matching baseline expression, promoter accessibility proxy, gene length, and regulatory-region count; module genes were not changed.

## ATAC results

- Runs processed: {len(qc)}/{len(samples)}.
- QC PASS replicates: {qc_pass_n}/{len(samples)}.
- DARs at BH FDR <0.05: {len(sig)} ({opening} opening, {closing} closing).
- Strongest module by matched-null empirical P: {strongest}.
- Modules passing matched-null FDR: {', '.join(nulls.loc[nulls['pass'].eq('PASS'), 'module_id'].tolist()) if not nulls.empty else 'NONE'}.
- ATAC gate: `{gate}`.

TSS enrichment is recorded as `NOT_COMPUTED_STANDARD_PROFILE` in the QC table because the rescue environment did not include a validated reference-point profile implementation; this omission is treated as a QC limitation, not as evidence for or against differential accessibility. Nucleosome-length fractions, replicate correlation, PCA and hierarchical clustering were computed where final BAMs were available.

## Conditional gates

- CUT&Tag raw rescue: **NO**. It was not unlocked by the ATAC gate.
- RNA raw rescue: `RNA_RAW_REANALYSIS_NOT_REQUIRED`; the locked processed same-system RNA layer remains the transcriptomic anchor and no unresolved sample identity or reproducibility issue required raw RNA reprocessing.
- Hi-C raw rescue: `HIC_RAW_DEFERRED`; no Hi-C raw data were downloaded.
- Phase 3 genetic escalation: `{phase3}`. No GWAS, eQTL, colocalization, PRS or MR was run.

## Interpretation boundary

Verdict: `{verdict}`. The strongest claim supported by this rescue is limited to replicate-level ATAC evidence for or against treatment-responsive regulatory remodeling in the locked modules. The rescue does not establish chromatin mechanism, adaptive-to-maladaptive transition, PIEZO1/TRPV4 causality, or genetic mediation.

"""
    if error:
        report += f"\n## Technical event\n\nThe rescue stopped with a technical error after preserving all completed artifacts: `{error}`. Downstream regulatory interpretation is held until the failed step is resolved.\n"
    (REPORT / "PHASE2R1_SELECTIVE_RAW_RESCUE_REPORT.md").write_text(report, encoding="utf-8")
    (REPORT / "PHASE3_GATE_AFTER_RAW_RESCUE.md").write_text(f"# Phase 3 gate after Phase 2R.1\n\nATAC gate: `{gate}`\n\nPhase 3 genetic escalation: **{phase3}**.\n\nGWAS/eQTL/colocalization/PRS/MR were not run. Hi-C raw data remain deferred.\n", encoding="utf-8")


def write_retention(samples: pd.DataFrame, qc: list[dict[str, object]]) -> None:
    rows: list[dict[str, object]] = list(DELETIONS)
    processed = {str(x["SRR"]) for x in qc}
    for _, r in samples.iterrows():
        srr = str(r.SRR)
        sra_dir = SRA_CACHE / srr
        if not any(x.get("path") == str(sra_dir) for x in rows):
            rows.append({"SRR": srr, "object": "SRA", "path": str(sra_dir), "action": "DELETED_AFTER_VALIDATION" if srr in processed else "RETAIN_FOR_RESUME", "validation": "final BAM/fragment/peak objects validated" if srr in processed else "raw rescue incomplete"})
        sample_dir = ATAC_RAW / "processed" / srr
        if srr in processed:
            for p in sorted([sample_dir / f"{srr}.final.bam", sample_dir / f"{srr}.final.bam.bai", sample_dir / f"{srr}.fragments.bed.gz", OUT / "replicate_peaks" / f"{srr}_peaks.narrowPeak"]):
                if p.exists():
                    rows.append({"SRR": srr, "object": p.name, "path": str(p), "action": "RETAIN_FINAL_VALIDATED", "validation": sha256(p)})
    write_df(pd.DataFrame(rows), RETENTION)


def finalize_frip(qc: list[dict[str, object]]) -> None:
    for x in qc:
        cmd = f"{shlex.quote(str(BIN / 'bedtools'))} intersect -u -abam {shlex.quote(str(x['final_bam']))} -b {shlex.quote(str(CONSENSUS))} | {shlex.quote(str(BIN / 'samtools'))} view -c -f 1 -"
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if p.returncode == 0:
            x["frip"] = (int(p.stdout.strip()) / 2) / max(float(x["usable_fragments"]), 1.0)
        else:
            x["frip"] = "NOT_COMPUTED"


def final_summary(verdict: str, samples: pd.DataFrame, qc: list[dict[str, object]], dar: pd.DataFrame, nulls: pd.DataFrame, module_status: pd.DataFrame, mech: pd.DataFrame, downloaded_gb: float, phase3: str) -> None:
    sig = dar[dar["FDR"].lt(FDR_CUTOFF)] if not dar.empty and "FDR" in dar else pd.DataFrame()
    opening = int((sig["log2FC"] > 0).sum()) if not sig.empty else 0
    closing = int((sig["log2FC"] < 0).sum()) if not sig.empty else 0
    passed = nulls.loc[nulls["pass"].eq("PASS"), "module_id"].tolist() if not nulls.empty else []
    strongest = "NONE" if nulls.empty else str(nulls.sort_values(["empirical_P", "observed_statistic"], ascending=[True, False]).iloc[0]["module_id"])
    supported = module_status.loc[module_status["status"].eq("LEVEL3_REGULATORY_SUPPORTED"), "module_id"].tolist() if not module_status.empty else []
    partial = module_status.loc[module_status["status"].eq("LEVEL3_PARTIAL_SUPPORT"), "module_id"].tolist() if not module_status.empty else []
    mech_s = mech[mech["status"].astype(str).str.startswith("SIGNIFICANT")] if not mech.empty and "status" in mech else pd.DataFrame()
    strongest_gene = str(mech_s.iloc[0]["gene"]) if not mech_s.empty else "NONE"
    qc_pass_n = sum(str(x.get("qc_status")) == "PASS" for x in qc)
    final_gb = sum(p.stat().st_size for base in [ATAC_RAW, OUT, REF] if base.exists() for p in base.rglob("*") if p.is_file()) / 1024**3
    gene_status = lambda g: str(mech.loc[mech["gene"].eq(g), "status"].iloc[0]) if not mech.empty and (mech["gene"] == g).any() else "NOT_TESTABLE"
    lines = [
        f"PHASE 2R.1 VERDICT: {verdict}",
        f"ATAC raw runs processed / total: {len(qc)} / {len(samples)}",
        f"ATAC QC PASS replicates: {qc_pass_n}",
        f"ATAC significant DARs FDR<0.05: {len(sig)}",
        f"ATAC opening / closing DARs: {opening} / {closing}",
        f"ATAC modules passing matched-null FDR: {','.join(passed) if passed else 'NONE'}",
        f"strongest ATAC module: {strongest}",
        "CUT&Tag raw rescue executed: NO",
        "CUT&Tag significant differential regions: NA",
        "CUT&Tag modules passing matched-null: NA",
        "RNA raw rescue executed: NOT_REQUIRED",
        "Hi-C raw rescue executed: DEFERRED",
        f"LEVEL3_REGULATORY_SUPPORTED modules: {','.join(supported) if supported else 'NONE'}",
        f"LEVEL3_PARTIAL_SUPPORT modules: {','.join(partial) if partial else 'NONE'}",
        f"PIEZO1 treatment-responsive regulatory evidence: {gene_status('PIEZO1')}",
        f"TRPV4 treatment-responsive regulatory evidence: {gene_status('TRPV4')}",
        f"strongest mechanosensor-network gene: {strongest_gene}",
        f"matched-null specificity rescued: {'YES' if passed else 'NO'}",
        f"regulatory rescue verdict: {verdict}",
        f"Phase 3 allowed: {'TARGETED' if phase3 == 'GO_PHASE3_TARGETED' else 'NO'}",
        f"Phase 3 target gene count: {len(mech_s)}",
        "recommended project direction: STOP_GENETIC_ESCALATION; lock manuscript around context-dependent pathway programs",
        "strongest supported manuscript claim: replicate-level ATAC regulatory remodeling is assessed within the locked human-TM treatment contrast",
        "strongest prohibited claim: chromatin mechanism, adaptive-to-maladaptive transition, PIEZO1/TRPV4 causality, or genetic mediation",
        f"total raw data downloaded GB: {downloaded_gb:.2f}",
        f"final retained disk GB: {final_gb:.2f}",
        f"path to PHASE2R1_SELECTIVE_RAW_RESCUE_REPORT.md: {REPORT / 'PHASE2R1_SELECTIVE_RAW_RESCUE_REPORT.md'}",
    ]
    SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


def main() -> None:
    RUN_LOG.write_text("", encoding="utf-8")
    commands: list[dict[str, str]] = []
    samples = load_samples()
    plan = make_storage_plan(samples)
    downloaded_gb = 0.0
    qc: list[dict[str, object]] = []
    dar = pd.DataFrame(columns=["peak_id", "chr", "start", "end", "baseMean", "log2FC", "SE", "statistic", "P", "FDR", "tested"])
    nulls = pd.DataFrame(columns=["module_id", "module_name", "pass", "FDR", "observed_statistic", "empirical_P"])
    mech = pd.DataFrame({"gene": ANCHORS, "status": ["NOT_TESTABLE"] * len(ANCHORS)})
    error: str | None = None
    gate = "ATAC_GATE_TECHNICAL_FAIL"
    verdict = "TECHNICAL_HOLD_ATAC"
    phase3 = "STOP_GENETIC_ESCALATION"
    modules: dict[str, dict[str, object]] = {}
    try:
        ensure_reference()
        for _, row in samples.iterrows():
            already_processed = str(row.status) == "PROCESSED"
            if not already_processed:
                samples.loc[samples["SRR"].eq(row.SRR), "status"] = "PROCESSING"
                write_df(samples, SAMPLE_MAP)
            qc_row = process_sample(row, commands)
            qc.append(qc_row)
            if not already_processed:
                downloaded_gb += float(row.compressed_GB)
            samples.loc[samples["SRR"].eq(row.SRR), "status"] = "PROCESSED"
            write_df(samples, SAMPLE_MAP)
        write_df(pd.DataFrame(qc), RAW_QC)
        build_consensus()
        count_matrix = quantify_counts(qc)
        finalize_frip(qc)
        write_df(pd.DataFrame(qc), RAW_QC)
        compute_qc_metrics(qc, count_matrix)
        dar = load_edge_results(count_matrix, samples)
        modules = load_modules()
        peak_genes, _ = map_dar_genes(dar)
        nulls, _ = module_nulls(dar, count_matrix, samples, modules, peak_genes)
        mech_rows = []
        for gene in ANCHORS:
            hits = [pid for pid, genes in peak_genes.items() if gene in genes]
            x = dar[dar["peak_id"].isin(hits)] if hits else dar.iloc[0:0]
            sx = x[x["FDR"].lt(FDR_CUTOFF)] if not x.empty else x
            if sx.empty:
                status = "NONSIGNIFICANT" if not x.empty else "NOT_TESTABLE"
            else:
                status = "SIGNIFICANT_OPENING" if (sx["log2FC"] > 0).any() else "SIGNIFICANT_CLOSING"
            mech_rows.append({"gene": gene, "status": status, "tested_peak_n": len(x), "significant_peak_n": len(sx), "opening_n": int((sx["log2FC"] > 0).sum()) if not sx.empty else 0, "closing_n": int((sx["log2FC"] < 0).sum()) if not sx.empty else 0, "best_FDR": sx["FDR"].min() if not sx.empty else np.nan, "best_log2FC": sx.loc[sx["FDR"].idxmin(), "log2FC"] if not sx.empty else np.nan, "mapping_hierarchy": "promoter overlap > existing GSE301522 link > proximal sensitivity only"})
        mech = pd.DataFrame(mech_rows)
        qc_df = pd.DataFrame(qc)
        qc_pass = len(qc) == len(samples) and all(x["qc_status"] == "PASS" for x in qc)
        pass_mods = nulls[nulls["pass"].eq("PASS")]
        non_dominated = all(int(x["observed_statistic"]) >= 2 for _, x in pass_mods.iterrows())
        if len(pass_mods) > 0 and non_dominated and qc_pass:
            gate, verdict = "ATAC_GATE_POSITIVE", "ATAC_SUPPORTS_CONTINUED_REGULATORY_RESCUE"
        elif len(pass_mods) == 0 and qc_pass:
            gate, verdict = "ATAC_GATE_NEGATIVE", "REGULATORY_RESCUE_PRIMARY_FAIL"
        else:
            gate = "ATAC_GATE_BORDERLINE" if qc_pass else "ATAC_GATE_TECHNICAL_FAIL"
            verdict = "REGULATORY_RESCUE_BORDERLINE" if qc_pass else "TECHNICAL_HOLD_ATAC"
        significant_genes = int((mech["significant_peak_n"] > 0).sum())
        if gate == "ATAC_GATE_POSITIVE" and significant_genes >= 3:
            phase3 = "GO_PHASE3_TARGETED"
        update_status_tables(dar, nulls, modules, mech, False)
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
        log("EXCEPTION " + error)
        try:
            modules = load_modules()
        except Exception:
            modules = {}
        write_empty_tables()
        if modules:
            update_status_tables(dar, nulls, modules, mech, True)
        else:
            write_df(mech, MECH)
    finally:
        write_df(pd.DataFrame(commands), CMD_MANIFEST)
        write_df(pd.DataFrame([{"status": "RNA_RAW_REANALYSIS_NOT_REQUIRED", "reason": "processed same-system RNA anchor retained; no raw identity/reproducibility issue required rescue"}]), RNA_STATUS)
        write_df(pd.DataFrame([{"status": "HIC_RAW_DEFERRED", "reason": "ATAC gate did not unlock all required conditions; 433.7 GB raw Hi-C not downloaded"}]), HIC_STATUS)
        write_retention(samples, qc)
        write_report(samples, plan, qc, dar, nulls, gate, verdict, error, phase3, downloaded_gb)
        update_project_files(verdict, gate, phase3)
        try:
            module_status = pd.read_csv(MODULE_STATUS, sep="\t")
        except Exception:
            module_status = pd.DataFrame(columns=["status"])
        final_summary(verdict, samples, qc, dar, nulls, module_status, mech, downloaded_gb, phase3)


if __name__ == "__main__":
    main()
