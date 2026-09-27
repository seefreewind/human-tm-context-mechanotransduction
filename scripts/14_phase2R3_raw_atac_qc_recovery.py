#!/usr/bin/env python3
"""Phase 2R.3 raw ATAC QC recovery and project closure.

This script is deliberately limited to technical QC recovery.  It audits
retained raw/pre-filter objects first, reuses them when available, and only
then reacquires missing public SRA objects one sample at a time.  Temporary
FASTQ/alignment objects are removed after metric validation; metrics, commands,
logs and checksums are retained.  No peaks, DARs, modules or new biological
discoveries are generated here.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "results" / "phase2R1"
P2 = ROOT / "results" / "phase2R3"
REPORT = ROOT / "reports" / "phase2R3"
LOG = ROOT / "logs" / "phase2R3"
RAW = ROOT / "data" / "raw" / "external" / "phase2R1" / "atac"
PROCESSED = RAW / "processed"
SRA_CACHE = RAW / "sra_cache"
FASTQ = RAW / "fastq"
TEMP = RAW / "temp"
REF = ROOT / "data" / "reference" / "hg38"
SAMPLE_MAP = P1 / "ATAC_RAW_SAMPLE_MAP.tsv"
OLD_QC = P2 / "ATAC_RAW_QC_COMPLETED.tsv"
FASTQ_INPUT_MAP = P2 / "FASTQ_INPUT_MAP.tsv"
SOURCE_AUDIT = P2 / "ATAC_QC_RECOVERY_SOURCE.current.tsv"
LEGACY_PARTIAL = P2 / ".ATAC_FINAL_QC.partial.tsv"
# The legacy checkpoint is read-only provenance.  New recovery state is written
# separately so no prior Phase 2R.3 intermediate is overwritten.
PARTIAL = P2 / "ATAC_FINAL_QC.phase2r3a.partial.tsv"
DAR = P1 / "ATAC_LEVEL3_DAR.tsv"
NULLS = ROOT / "results" / "phase2R2" / "ATAC_LEVEL3_MODULE_NULLS_FINAL.tsv"

ENV_PREFIX = Path(os.environ.get("PHASE2R1_ENV_PREFIX", os.environ.get("CONDA_PREFIX", sys.prefix)))
BIN = ENV_PREFIX / "bin"
THREADS = max(2, min(8, os.cpu_count() or 4))
MAPPING_THRESHOLD = 0.50
FRAGMENT_THRESHOLD = 1_000_000
MAX_DUPLICATE_FRACTION = 0.50
MAX_CHRM_FRACTION = 0.20
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
RUN_ID = time.strftime("%Y%m%d_%H%M%S")
RUN_LOG = LOG / "phase2r3_run.log"
COMMANDS = P2 / "ATAC_QC_RECOVERY_COMMANDS.current.tsv"
REGENERATED_VALIDATION = P2 / "REGENERATED_BAM_VALIDATION.tsv"

for d in (P2, REPORT, LOG, SRA_CACHE, FASTQ, TEMP):
    d.mkdir(parents=True, exist_ok=True)


def log(message: str) -> None:
    with RUN_LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}\n")


def run(cmd: list[str] | str, name: str, check: bool = True, shell: bool = False) -> subprocess.CompletedProcess[str]:
    rendered = cmd if isinstance(cmd, str) else " ".join(shlex.quote(str(x)) for x in cmd)
    log(f"COMMAND {name}: {rendered}")
    out = LOG / f"{name}.log"
    with out.open("a", encoding="utf-8") as fh:
        p = subprocess.run(cmd, shell=shell, stdout=fh, stderr=subprocess.STDOUT, text=True, check=False)
    log(f"COMMAND_RESULT {name}: returncode={p.returncode}")
    if check and p.returncode != 0:
        raise RuntimeError(f"command failed: {name}; see {out}")
    return p


def write_df(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def load_fastq_input_map(srr: str) -> dict[str, Path] | None:
    """Return an explicit, read-only authoritative FASTQ mapping if registered."""
    if not FASTQ_INPUT_MAP.exists():
        return None
    df = pd.read_csv(FASTQ_INPUT_MAP, sep="\t", dtype=str)
    x = df.loc[df.SRR.astype(str).eq(str(srr))].copy()
    if x.empty or not {"read", "mapped_input"}.issubset(x.columns):
        return None
    out: dict[str, Path] = {}
    for _, row in x.iterrows():
        read = str(row["read"])
        path = Path(str(row["mapped_input"]))
        if read in {"R1", "R2"}:
            out[read] = path
    if set(out) != {"R1", "R2"}:
        raise RuntimeError(f"incomplete FASTQ input map for {srr}")
    for read, path in out.items():
        if not path.is_file():
            raise RuntimeError(f"mapped FASTQ missing for {srr} {read}: {path}")
    return out


def current_path(value: object) -> str:
    """Normalize legacy absolute paths only in newly generated artifacts."""
    text = str(value)
    marker = "/piezo1_trpv4_glaucoma/"
    if marker in text:
        suffix = text.split(marker, 1)[-1]
        return str(ROOT / suffix)
    return text


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def size_bytes(path: Path) -> int:
    try:
        return int(path.stat().st_size)
    except OSError:
        return 0


def free_gib() -> float:
    return shutil.disk_usage(ROOT).free / (1024 ** 3)


def locate_one(directory: Path, patterns: list[str]) -> Path | None:
    if not directory.exists():
        return None
    candidates: list[Path] = []
    for pattern in patterns:
        candidates.extend(p for p in directory.rglob(pattern) if p.is_file() and not p.name.startswith("._"))
    return sorted(set(candidates), key=lambda p: p.stat().st_size, reverse=True)[0] if candidates else None


def audit_sources(samples: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in samples.iterrows():
        srr = str(row.SRR)
        mapped = load_fastq_input_map(srr)
        sra = locate_one(SRA_CACHE / srr, [f"{srr}.sra", f"{srr}.sralite", "*.sra", "*.sralite"])
        fq1 = mapped["R1"] if mapped else locate_one(FASTQ / srr, [f"{srr}_1.fastq", f"{srr}_1.fastq.gz"])
        fq2 = mapped["R2"] if mapped else locate_one(FASTQ / srr, [f"{srr}_2.fastq", f"{srr}_2.fastq.gz"])
        pre_bam = locate_one(PROCESSED / srr, [f"{srr}.marked.bam", f"{srr}.mapped.bam", f"{srr}.noblack.bam", f"{srr}.nochrM.bam", "*.pre_filter.bam"])
        logs = sorted(str(p) for p in (PROCESSED / srr).glob("*log*") if p.is_file() and not p.name.startswith("._"))
        final_bam = PROCESSED / srr / f"{srr}.final.bam"
        sources = []
        if sra:
            sources.append("SRA")
        if fq1 and fq2:
            sources.append("AUDITED_FASTQ_MAP" if mapped else "FASTQ")
        if pre_bam:
            sources.append("PRE_FILTER_BAM")
        if logs:
            sources.append("ALIGNMENT_LOG")
        if final_bam.exists() and not (sra or (fq1 and fq2) or pre_bam):
            sources.append("FINAL_BAM_ONLY")
        local_path = mapped["R1"] if mapped else (sra or fq1 or pre_bam or final_bam)
        rows.append({
            "SRR": srr, "GSM": row.GSM, "condition": row.condition, "cell_line": row.cell_line,
            "biological_replicate": row.biological_replicate,
            "local_SRA": "YES" if sra else "NO", "local_FASTQ_paired": "YES" if fq1 and fq2 else "NO",
            "local_pre_filter_BAM": "YES" if pre_bam else "NO", "alignment_log": "YES" if logs else "NO",
            "final_BAM": "YES" if final_bam.exists() else "NO",
            "local_source_status": ";".join(sources) if sources else "NOT_RETAINED",
            "recovery_source": "AUDITED_RUN3_FASTQ_REBUILD" if mapped else ("LOCAL_REUSE" if (sra or (fq1 and fq2) or pre_bam) else "FINAL_BAM_ONLY_INCOMPLETE"),
            "local_path": current_path(local_path),
            "local_size_bytes": size_bytes(local_path) if local_path else 0,
            "initial_audit_run": RUN_ID,
        })
    out = pd.DataFrame(rows)
    write_df(out, SOURCE_AUDIT)
    return out


def find_sra(srr: str) -> Path | None:
    return locate_one(SRA_CACHE / srr, [f"{srr}.sra", f"{srr}.sralite", "*.sra", "*.sralite"])


def recover_fastq_from_temp_chunks(srr: str, temp_dir: Path, fq1: Path, fq2: Path) -> bool:
    """Recover completed fasterq-dump chunks after an interrupted assembly."""
    chunk_dirs = sorted(
        p for p in temp_dir.glob("fasterq.tmp.*")
        if p.is_dir() and not p.name.startswith("._")
    )
    if len(chunk_dirs) != 1:
        return False
    chunk_dir = chunk_dirs[0]
    pattern = re.compile(r"\.(\d+)\.(1|2)$")
    by_read: dict[str, list[tuple[int, Path]]] = {"1": [], "2": []}
    for path in chunk_dir.iterdir():
        if not path.is_file() or path.name.startswith("._"):
            continue
        match = pattern.search(path.name)
        if match:
            by_read[match.group(2)].append((int(match.group(1)), path))
    if not by_read["1"] or {i for i, _ in by_read["1"]} != {i for i, _ in by_read["2"]}:
        return False
    for read in ("1", "2"):
        target = fq1 if read == "1" else fq2
        temporary = target.with_suffix(target.suffix + ".recovered.tmp")
        with temporary.open("wb") as out:
            for _, chunk in sorted(by_read[read]):
                with chunk.open("rb") as src:
                    shutil.copyfileobj(src, out, length=1024 * 1024)
        os.replace(temporary, target)
    log(f"RECOVERED_FASTQ_FROM_TEMP {srr}: chunks={len(by_read['1'])} fq1_bytes={fq1.stat().st_size} fq2_bytes={fq2.stat().st_size}")
    return True


def obtain_fastq(row: pd.Series) -> tuple[Path | None, Path, Path, str]:
    """Use only retained/local inputs; public reacquisition and fasterq-dump are closed."""
    srr = str(row.SRR)
    mapped = load_fastq_input_map(srr)
    if mapped:
        sra = find_sra(srr)
        log(f"AUTHORITATIVE_FASTQ_MAP {srr}: R1={mapped['R1']} R2={mapped['R2']} sra_retained={bool(sra)}")
        return sra, mapped["R1"], mapped["R2"], "AUDITED_RUN3_FASTQ_REBUILD"
    sample_fastq = FASTQ / srr
    fq1 = sample_fastq / f"{srr}_1.fastq"
    fq2 = sample_fastq / f"{srr}_2.fastq"
    if fq1.is_file() and fq2.is_file():
        return find_sra(srr), fq1, fq2, "LOCAL_FASTQ_REUSE"
    raise RuntimeError(f"NO_RETAINED_PAIRED_FASTQ_NO_PUBLIC_REACQUISITION_ALLOWED:{srr}")


def count(cmd: list[str], name: str) -> int:
    rendered = " ".join(shlex.quote(str(x)) for x in cmd)
    log(f"COMMAND {name}: {rendered}")
    q = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (LOG / f"{name}.log").write_text((q.stdout or "") + (q.stderr or ""), encoding="utf-8")
    log(f"COMMAND_RESULT {name}: returncode={q.returncode}")
    if q.returncode != 0:
        raise RuntimeError(f"command failed ({name}); see {LOG / f'{name}.log'}")
    try:
        return int(q.stdout.strip().splitlines()[-1])
    except (IndexError, ValueError) as exc:
        raise RuntimeError(f"could not parse scalar output for {name}: {q.stdout!r}") from exc


def validate_regenerated_bam(srr: str, metrics: dict[str, object]) -> pd.DataFrame:
    """Validate the regenerated name-sorted BAM before any final QC merge."""
    bam = Path(str(metrics["regenerated_name_bam"]))
    sample_dir = bam.parent
    fastp_json = sample_dir / f"{srr}.phase2r3.regenerated.fastp.json"
    mapped = load_fastq_input_map(srr)
    rows: list[dict[str, object]] = []

    def add(check: str, observed: object, expected: object, status: str, evidence: str) -> None:
        rows.append({"SRR": srr, "object": str(bam), "check": check, "observed": observed, "expected": expected, "status": status, "evidence": evidence})

    quick = subprocess.run([str(BIN / "samtools"), "quickcheck", "-v", str(bam)], capture_output=True, text=True, check=False)
    add("samtools_quickcheck", "PASS" if quick.returncode == 0 else (quick.stderr or quick.stdout).strip(), "returncode=0", "PASS" if quick.returncode == 0 else "FAIL", "samtools quickcheck -v")

    header = subprocess.run([str(BIN / "samtools"), "view", "-H", str(bam)], capture_output=True, text=True, check=False)
    hd = next((line for line in header.stdout.splitlines() if line.startswith("@HD")), "")
    # `samtools view -H` returns literal tab-delimited header fields.
    sort_order = next((x.split(":", 1)[1] for x in hd.split("\t") if x.startswith("SO:")), "MISSING")
    add("header_inspection", "PASS" if header.returncode == 0 else "FAIL", "header readable", "PASS" if header.returncode == 0 else "FAIL", hd or "samtools view -H failed")
    add("sort_order", sort_order, "queryname", "PASS" if sort_order == "queryname" else "FAIL", "@HD SO field")

    bam_total = count([str(BIN / "samtools"), "view", "-c", str(bam)], f"validate_total_{srr}")
    bam_paired = count([str(BIN / "samtools"), "view", "-c", "-f", "1", str(bam)], f"validate_paired_{srr}")
    bam_proper = count([str(BIN / "samtools"), "view", "-c", "-f", "2", "-F", "4", str(bam)], f"validate_proper_{srr}")
    add("paired_read_structure", f"paired_records={bam_paired};total_records={bam_total};properly_paired_mapped={bam_proper}", "all BAM records paired; paired count even", "PASS" if bam_total > 0 and bam_paired == bam_total and bam_paired % 2 == 0 else "FAIL", "samtools flag counts")

    audited_r1 = audited_r2 = np.nan
    if mapped:
        m = pd.read_csv(FASTQ_INPUT_MAP, sep="\t", dtype=str)
        for read in ("R1", "R2"):
            x = m.loc[(m.SRR.astype(str) == srr) & (m.read == read), "record_count"]
            if not x.empty:
                if read == "R1": audited_r1 = int(x.iloc[0])
                else: audited_r2 = int(x.iloc[0])
    fastp_before = fastp_after = np.nan
    if fastp_json.exists():
        d = json.loads(fastp_json.read_text())
        summary = d.get("summary", {})
        fastp_before = summary.get("before_filtering", {}).get("total_reads", np.nan)
        fastp_after = summary.get("after_filtering", {}).get("total_reads", np.nan)
    audited_status = audited_r1 == audited_r2 == 142540485 and fastp_before == 285080970
    add("audited_FASTQ_to_input_reconciliation", f"R1={audited_r1};R2={audited_r2};fastp_before={fastp_before}", "R1=R2=142540485;fastp_before=285080970", "PASS" if audited_status else "FAIL", "FASTQ_INPUT_MAP.tsv plus regenerated fastp JSON")
    bam_vs_fastp = bool(pd.notna(fastp_after) and bam_total == int(fastp_after))
    add("BAM_read_count_reconciliation", f"bam_total={bam_total};fastp_after={fastp_after}", "BAM total equals fastp after_filtering total_reads", "PASS" if bam_vs_fastp else "FAIL", "samtools view -c versus fastp JSON")
    out = pd.DataFrame(rows)
    write_df(out, REGENERATED_VALIDATION)
    if not out.status.eq("PASS").all():
        raise RuntimeError("regenerated BAM validation failed; see results/phase2R3/REGENERATED_BAM_VALIDATION.tsv")
    return out


def fastp_and_align(row: pd.Series, sra: Path | None, fq1: Path, fq2: Path) -> dict[str, object]:
    srr = str(row.SRR)
    sample_dir = PROCESSED / srr
    sample_dir.mkdir(parents=True, exist_ok=True)
    # All fallback outputs carry an explicit regenerated label and are kept
    # alongside, never in place of, historical/final objects.
    trim1 = sample_dir / f"{srr}.phase2r3.regenerated_trim_R1.fastq.gz"
    trim2 = sample_dir / f"{srr}.phase2r3.regenerated_trim_R2.fastq.gz"
    fastp_json = sample_dir / f"{srr}.phase2r3.regenerated.fastp.json"
    fastp_html = sample_dir / f"{srr}.phase2r3.regenerated.fastp.html"
    if not trim1.exists() or not trim2.exists() or not fastp_json.exists():
        run([str(BIN / "fastp"), "--thread", str(THREADS), "--in1", str(fq1), "--in2", str(fq2), "--out1", str(trim1), "--out2", str(trim2), "--json", str(fastp_json), "--html", str(fastp_html), "--detect_adapter_for_pe", "--length_required", "20"], f"fastp_{srr}")
    name_bam = sample_dir / f"{srr}.phase2r3.regenerated.name.bam"
    fix_bam = sample_dir / f"{srr}.phase2r3.regenerated.fixmate.bam"
    coord_bam = sample_dir / f"{srr}.phase2r3.regenerated.coord.bam"
    marked_bam = sample_dir / f"{srr}.phase2r3.regenerated.marked.bam"
    if not name_bam.exists():
        cmd = "set -o pipefail; " + " ".join([
            shlex.join([str(BIN / "bowtie2"), "--very-sensitive", "-p", str(THREADS), "-x", str(REF / "bowtie2" / "GRCh38_noalt_as"), "-1", str(trim1), "-2", str(trim2), "-X", "2000"]),
            "2>", shlex.quote(str(sample_dir / f"{srr}.phase2r3.regenerated.bowtie2.log")), "|",
            shlex.join([str(BIN / "samtools"), "view", "-b", "-u", "-"]), "|",
            shlex.join([str(BIN / "samtools"), "sort", "-n", "-@", str(THREADS), "-o", str(name_bam), "-"]),
        ])
        run(cmd, f"align_{srr}", shell=True)
    if not fix_bam.exists():
        run([str(BIN / "samtools"), "fixmate", "-m", str(name_bam), str(fix_bam)], f"fixmate_{srr}")
    if not coord_bam.exists():
        run([str(BIN / "samtools"), "sort", "-@", str(THREADS), "-o", str(coord_bam), str(fix_bam)], f"sort_{srr}")
    if not marked_bam.exists():
        run([str(BIN / "samtools"), "markdup", "-@", str(THREADS), str(coord_bam), str(marked_bam)], f"markdup_{srr}")
    raw_reads = int(json.loads(fastp_json.read_text()).get("summary", {}).get("before_filtering", {}).get("total_reads", int(row.spots) * 2))
    mapped = count([str(BIN / "samtools"), "view", "-c", "-F", "4", str(marked_bam)], f"mapped_{srr}")
    # Region-restricted samtools view requires an index.  The index is a new
    # regenerated auxiliary object; no historical/final BAM is modified.
    marked_bai = Path(str(marked_bam) + ".bai")
    if not marked_bai.exists():
        run([str(BIN / "samtools"), "index", "-@", str(THREADS), str(marked_bam)], f"index_regenerated_marked_{srr}")
    chrm = count([str(BIN / "samtools"), "view", "-c", "-F", "4", str(marked_bam), "chrM"], f"chrM_{srr}")
    duplicates = count([str(BIN / "samtools"), "view", "-c", "-f", "1024", str(marked_bam)], f"duplicates_{srr}")
    properly_paired = count([str(BIN / "samtools"), "view", "-c", "-f", "3", "-F", "4", str(marked_bam)], f"properly_paired_{srr}")
    usable_nuclear = count([str(BIN / "samtools"), "view", "-c", "-q", "30", "-f", "3", "-F", "3852", str(marked_bam), *CHROMS], f"usable_nuclear_{srr}")
    mapping_rate = mapped / max(raw_reads, 1)
    duplicate_fraction = duplicates / max(mapped, 1)
    chrm_fraction = chrm / max(mapped, 1)
    return {
        "SRR": srr, "GSM": row.GSM, "condition": row.condition, "cell_line": row.cell_line,
        "biological_replicate": row.biological_replicate, "raw_reads": raw_reads,
        "mapped_reads_raw_alignment": mapped, "mapping_rate_raw_alignment": mapping_rate,
        "chrM_aligned_reads": chrm, "chrM_fraction": chrm_fraction,
        "duplicate_reads": duplicates, "duplicate_fraction": duplicate_fraction,
        "properly_paired_reads": properly_paired, "usable_nuclear_reads": usable_nuclear,
        "usable_nuclear_fragments": usable_nuclear / 2.0,
        "tss_enrichment": "FROM_PHASE2R2_FINAL_BAM", "tss_enrichment_status": "RECOVERED_FINAL_BAM_PROFILE",
        "frip": "FROM_PHASE2R2_FINAL_BAM", "source_used": "AUDITED_FASTQ_OR_LOCAL_FASTQ_PRE_FILTER_ALIGNMENT",
        "raw_alignment_bam": str(marked_bam), "regenerated_name_bam": str(name_bam), "status": "RECOVERED",
    }


def cleanup_sample(srr: str, keep_sra: bool = True) -> list[dict[str, object]]:
    """Retain all recovery inputs/intermediates for audit and safe resume."""
    sample_dir = PROCESSED / srr
    retained = sorted(
        p for p in sample_dir.glob(f"{srr}.phase2r3.regenerated*")
        if p.is_file() and not p.name.startswith("._")
    )
    rows = [{
        "SRR": srr, "object": p.name, "path": str(p),
        "action": "RETAINED_FOR_PHASE2R3_AUDIT_AND_SAFE_RESUME",
        "size_bytes": size_bytes(p), "sha256_before_delete": sha256(p),
    } for p in retained]
    sra = find_sra(srr)
    if sra:
        rows.append({"SRR": srr, "object": sra.name, "path": str(sra), "action": "RETAINED_AUTHORITATIVE_SRA_SOURCE", "size_bytes": size_bytes(sra), "sha256_before_delete": sha256(sra)})
    return rows


def final_bam_only_metrics(row: pd.Series, old_qc: pd.DataFrame) -> dict[str, object]:
    """Carry forward final-BAM metrics while explicitly retaining missing raw QC."""
    srr = str(row.SRR)
    old = old_qc.loc[old_qc.SRR.astype(str).eq(srr)] if not old_qc.empty else pd.DataFrame()
    old_row = old.iloc[0].to_dict() if not old.empty else {}
    final_bam = Path(current_path(old_row.get("final_bam", PROCESSED / srr / f"{srr}.final.bam")))
    raw_reads = old_row.get("raw_reads", float(row.spots) * 2)
    mapped = old_row.get("mapped_reads_final_bam", np.nan)
    mapping = old_row.get("mapping_rate_final_bam_vs_raw_reads", np.nan)
    usable_fragments = old_row.get("usable_fragments", np.nan)
    return {
        "SRR": srr, "GSM": row.GSM, "condition": row.condition, "cell_line": row.cell_line,
        "biological_replicate": row.biological_replicate, "raw_reads": raw_reads,
        "mapped_reads_raw_alignment": mapped, "mapping_rate_raw_alignment": mapping,
        "chrM_aligned_reads": np.nan, "chrM_fraction": np.nan,
        "duplicate_reads": np.nan, "duplicate_fraction": np.nan,
        "properly_paired_reads": np.nan,
        "usable_nuclear_reads": float(usable_fragments) * 2 if pd.notna(usable_fragments) else np.nan,
        "usable_nuclear_fragments": usable_fragments,
        "tss_enrichment": old_row.get("tss_enrichment", np.nan),
        "tss_enrichment_status": old_row.get("tss_enrichment_status", "UNAVAILABLE"),
        "frip": old_row.get("frip", np.nan),
        "source_used": "EXISTING_FINAL_BAM_ONLY_NO_PRE_FILTER_RECOVERY",
        "raw_alignment_bam": str(final_bam),
        "status": "INCOMPLETE_FINAL_BAM_ONLY",
        "recovery_source": "FINAL_BAM_ONLY_INCOMPLETE",
    }


def load_completed_qc() -> pd.DataFrame:
    if OLD_QC.exists():
        return pd.read_csv(OLD_QC, sep="\t", dtype={"SRR": str})
    return pd.DataFrame()


def merge_final_qc(samples: pd.DataFrame, recovered: pd.DataFrame, old_qc: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, s in samples.iterrows():
        srr = str(s.SRR)
        r = recovered.loc[recovered.SRR.eq(srr)]
        old = old_qc.loc[old_qc.SRR.astype(str).eq(srr)] if not old_qc.empty else pd.DataFrame()
        reasons: list[str] = []
        if r.empty:
            reasons.append("raw_or_pre_filter_source_not_recovered")
            raw = {}
        else:
            raw = r.iloc[0].to_dict()
        old_row = old.iloc[0].to_dict() if not old.empty else {}
        mapping = float(raw.get("mapping_rate_raw_alignment", np.nan))
        fragments = float(raw.get("usable_nuclear_fragments", np.nan))
        dup = float(raw.get("duplicate_fraction", np.nan))
        chrm = float(raw.get("chrM_fraction", np.nan))
        for val, label in [(mapping, "mapping_rate"), (fragments, "usable_nuclear_fragments"), (dup, "duplicate_fraction"), (chrm, "chrM_fraction")]:
            if not np.isfinite(val):
                reasons.append(f"{label}_missing")
        if np.isfinite(mapping) and mapping < MAPPING_THRESHOLD:
            reasons.append(f"mapping_rate<{MAPPING_THRESHOLD}")
        if np.isfinite(fragments) and fragments < FRAGMENT_THRESHOLD:
            reasons.append(f"usable_nuclear_fragments<{FRAGMENT_THRESHOLD}")
        if np.isfinite(dup) and dup > MAX_DUPLICATE_FRACTION:
            reasons.append(f"duplicate_fraction>{MAX_DUPLICATE_FRACTION}")
        if np.isfinite(chrm) and chrm > MAX_CHRM_FRACTION:
            reasons.append(f"chrM_fraction>{MAX_CHRM_FRACTION}")
        if str(old_row.get("tss_enrichment_status", "")).startswith("INCOMPLETE"):
            reasons.append("TSS_enrichment_unavailable")
        if str(old_row.get("frip", "NA")) in {"NA", "NOT_COMPUTED", ""}:
            reasons.append("FRiP_unavailable")
        if not r.empty and not reasons:
            status = "PASS"
        elif r.empty or any(x.endswith("_missing") or "unavailable" in x for x in reasons):
            status = "INCOMPLETE"
        else:
            status = "FAIL"
        rows.append({
            "SRR": srr, "GSM": s.GSM, "condition": s.condition, "cell_line": s.cell_line,
            "biological_replicate": s.biological_replicate, "raw_reads": raw.get("raw_reads", np.nan),
            "mapped_reads_raw_alignment": raw.get("mapped_reads_raw_alignment", np.nan),
            "mapping_rate_raw_alignment": mapping, "chrM_aligned_reads": raw.get("chrM_aligned_reads", np.nan),
            "chrM_fraction": chrm, "duplicate_reads": raw.get("duplicate_reads", np.nan),
            "duplicate_fraction": dup, "properly_paired_reads": raw.get("properly_paired_reads", np.nan),
            "usable_nuclear_reads": raw.get("usable_nuclear_reads", np.nan),
            "usable_nuclear_fragments": fragments, "tss_enrichment": old_row.get("tss_enrichment", np.nan),
            "tss_enrichment_status": old_row.get("tss_enrichment_status", "UNAVAILABLE"),
            "frip": old_row.get("frip", np.nan),
            "duplicate_fraction_status": "RECOVERED_REGENERATED_MARKED_BAM" if np.isfinite(dup) else "INCOMPLETE_PREDEDUP_OBJECT_NOT_RETAINED",
            "chrM_fraction_status": "RECOVERED_REGENERATED_PRE_FILTER_ALIGNMENT" if np.isfinite(chrm) else "INCOMPLETE_PRE_FILTER_OBJECT_NOT_RETAINED",
            "raw_metric_source": raw.get("source_used", "NO_RECOVERED_RAW_SOURCE"),
            "QC_status": status,
            "reason": ";".join(reasons) if reasons else "all_locked_QC_rules_met",
            "locked_thresholds": "mapping>=0.50;usable_nuclear_fragments>=1000000;duplicate_fraction<=0.50;chrM_fraction<=0.20",
            "final_bam": current_path(old_row.get("final_bam", str(PROCESSED / srr / f"{srr}.final.bam"))),
        })
    out = pd.DataFrame(rows)
    write_df(out, P2 / "ATAC_FINAL_QC.tsv")
    return out


def qc_balance(qc: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for factor in ["condition", "cell_line", "biological_replicate"]:
        for level in sorted(qc[factor].astype(str).unique()):
            x = qc[qc[factor].astype(str).eq(level)]
            rows.append({"record_type": "summary", "factor": factor, "level": level, "sample_n": len(x), "pass_n": int(x.QC_status.eq("PASS").sum()), "fail_n": int(x.QC_status.eq("FAIL").sum()), "incomplete_n": int(x.QC_status.eq("INCOMPLETE").sum()), "pass_rate": float(x.QC_status.eq("PASS").mean()), "audit_note": "technical QC balance only; no treatment interpretation and no sample removal"})
    rows.append({"record_type": "overall", "factor": "none", "level": "all", "sample_n": len(qc), "pass_n": int(qc.QC_status.eq("PASS").sum()), "fail_n": int(qc.QC_status.eq("FAIL").sum()), "incomplete_n": int(qc.QC_status.eq("INCOMPLETE").sum()), "pass_rate": float(qc.QC_status.eq("PASS").mean()) if len(qc) else np.nan, "audit_note": "technical QC balance only"})
    for factor in ["condition", "cell_line"]:
        grouped = qc.groupby(factor)["QC_status"].apply(lambda x: float(x.eq("PASS").mean()))
        incomplete_rates = qc.groupby(factor)["QC_status"].apply(lambda x: float(x.eq("INCOMPLETE").mean()))
        imbalance = bool(len(grouped) > 1 and (grouped.max() != grouped.min() or incomplete_rates.max() != incomplete_rates.min()))
        rows.append({"record_type": "imbalance", "factor": factor, "level": "all", "sample_n": len(qc), "pass_n": int(qc.QC_status.eq("PASS").sum()), "fail_n": int(qc.QC_status.eq("FAIL").sum()), "incomplete_n": int(qc.QC_status.eq("INCOMPLETE").sum()), "pass_rate": float(qc.QC_status.eq("PASS").mean()) if len(qc) else np.nan, "qc_imbalance": "YES" if imbalance else "NO", "audit_note": "descriptive technical missingness/pass balance only; no treatment or cell-line interpretation"})
    out = pd.DataFrame(rows)
    write_df(out, P2 / "ATAC_FINAL_QC_BALANCE.tsv")
    return out


def dar_validity(samples: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    dar = pd.read_csv(DAR, sep="\t")
    nulls = pd.read_csv(NULLS, sep="\t") if NULLS.exists() else pd.DataFrame()
    pass_samples = samples[samples.SRR.astype(str).isin(qc.loc[qc.QC_status.eq("PASS"), "SRR"].astype(str))].copy()
    cond_n = pass_samples.condition.value_counts().to_dict() if not pass_samples.empty else {}
    design = pd.DataFrame({"intercept": 1.0, "cell_line_TM70": (pass_samples.cell_line == "TM70").astype(int), "treatment_DEX": (pass_samples.condition == "DEX").astype(int)}) if not pass_samples.empty else pd.DataFrame()
    rank = int(np.linalg.matrix_rank(design.to_numpy(dtype=float))) if not design.empty else 0
    full_rank = rank == design.shape[1] if not design.empty else False
    module_pass = int((nulls.FDR < 0.05).sum()) if not nulls.empty and "FDR" in nulls else 0
    if len(qc) == len(samples) and qc.QC_status.eq("PASS").all():
        result = "RETAIN_EXISTING_FULL_SET_DAR"
        reason = "all six samples satisfy locked QC; existing edgeR full-set DAR remains primary and no new discovery was run"
    elif full_rank and cond_n.get("CONTROL", 0) > 0 and cond_n.get("DEX", 0) > 0:
        result = "QC_PASS_ONLY_SENSITIVITY_ALLOWED_BUT_NOT_REQUIRED_FOR_PRIMARY"
        reason = "QC-pass-only design is identifiable; the locked full-set result remains primary and no replacement analysis was selected"
    else:
        result = "QC_PASS_ONLY_NOT_IDENTIFIABLE"
        reason = "QC-pass-only subset lacks a full-rank cell_line+treatment design; no sensitivity differential analysis was run"
    rows = [
        {"analysis": "locked_full_set", "sample_n": len(samples), "control_n": int((samples.condition == "CONTROL").sum()), "dex_n": int((samples.condition == "DEX").sum()), "significant_DAR_n": int((dar.FDR < 0.05).sum()), "opening_n": int(((dar.FDR < 0.05) & (dar.log2FC > 0)).sum()), "closing_n": int(((dar.FDR < 0.05) & (dar.log2FC < 0)).sum()), "design_rank": "3/3", "status": "RETAINED_PRIMARY", "module_pass_n": module_pass, "reason": "verified existing locked edgeR result under ~ cell_line + treatment"},
        {"analysis": "QC_pass_only_sensitivity", "sample_n": len(pass_samples), "control_n": int(cond_n.get("CONTROL", 0)), "dex_n": int(cond_n.get("DEX", 0)), "significant_DAR_n": np.nan, "opening_n": np.nan, "closing_n": np.nan, "design_rank": f"{rank}/{design.shape[1] if not design.empty else 0}", "status": result, "module_pass_n": module_pass, "reason": reason},
    ]
    out = pd.DataFrame(rows)
    write_df(out, P2 / "ATAC_FINAL_DAR_VALIDITY.tsv")
    return out


def retention(rows: list[dict[str, object]], source_audit: pd.DataFrame, recovered: pd.DataFrame) -> pd.DataFrame:
    all_rows = list(rows)
    for _, row in source_audit.iterrows():
        srr = str(row.SRR)
        rr = recovered.loc[recovered.SRR.eq(srr)]
        all_rows.append({"SRR": srr, "object": "retained_source_or_final_bam", "path": current_path(row.local_path), "action": "RETAINED_OR_FINAL_BAM_ONLY_NO_REACQUISITION", "size_bytes": row.local_size_bytes, "sha256_before_delete": "NOT_DELETED;SOURCE_OR_FINAL_OBJECT_RETAINED"})
    out = pd.DataFrame(all_rows)
    write_df(out, P2 / "ATAC_QC_RECOVERY_RETENTION.tsv")
    return out


def update_project_files(verdict: str, classification: str, source_status: str) -> None:
    marker = "## Phase 2R.3A input recovery and safe BAM regeneration — 2026-09-26"
    status = f"{marker}\nphase2R3_historical_name_bam: HISTORICAL_BAM_UNRECOVERABLE\nphase2R3_audited_run3_fastq_fallback: AUDITED_RUN3_FASTQ_REBUILD\nphase2R3_verdict: {verdict}\nphase2R3_atac_classification: {classification}\nphase2R3_raw_source_status: {source_status}\nphase2R3_module_specific_regulatory_support: NO_MODULE_SPECIFIC_REGULATORY_SUPPORT\nphase2R3_regulatory_gate: STOP_REGULATORY_ESCALATION\nphase2R3_genetic_gate: STOP_GENETIC_ESCALATION\nphase2R3_cuttag_raw_executed: false\nphase2R3_hic_raw_executed: false\nphase2R3_phase3_gate: CLOSED\nphase2R3_report: reports/phase2R3/PHASE2R3_FINAL_ATAC_QC_REPORT.md\n"
    for path in [ROOT / "PROJECT_STATUS.md", ROOT / "ANALYSIS_GOVERNANCE.md"]:
        text = path.read_text(encoding="utf-8")
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + status, encoding="utf-8")
    cfg = ROOT / "config" / "project_config.yaml"
    text = cfg.read_text(encoding="utf-8")
    cfg_marker = f"# {marker}"
    if cfg_marker not in text:
        # Keep the YAML syntactically valid while preserving the same audit
        # addendum recorded in the Markdown status files.
        cfg_status = "\n".join(f"# {line}" for line in status.splitlines())
        cfg.write_text(text.rstrip() + "\n\n" + cfg_status + "\n", encoding="utf-8")


def manuscript_hierarchy(classification: str) -> None:
    text = f"""# Final manuscript evidence hierarchy — Phase 2R.3

ATAC final classification: **{classification}**

## Locked evidence hierarchy

1. External acute mechanical-stretch transcriptomics.
2. Pathological TGFβ2 / DEX transcriptomic contexts.
3. GSE313736 large human trabecular-meshwork transportability.
4. Donor-paired descriptive concordance.
5. Replicate-level exploratory ATAC accessibility changes, retained only at the confidence level supported by final QC.
6. Absence of module-specific matched-null regulatory enrichment: 0/9 modules pass.
7. PIEZO1/TRPV4 as mechanosensory interpretation anchors.

## Allowed interpretation

The ATAC layer can support exploratory accessibility differences only to the extent indicated by the final QC table. The matched-null analysis provides **NO_MODULE_SPECIFIC_REGULATORY_SUPPORT**.

## Prohibited interpretation

Do not claim established chromatin rewiring, established mechanotransduction regulatory remodeling, an adaptive-to-maladaptive transition, PIEZO1/TRPV4 causality, genetic mediation, or a `PRS × DEX` effect.
"""
    (REPORT / "FINAL_MANUSCRIPT_EVIDENCE_HIERARCHY.md").write_text(text, encoding="utf-8")


def final_report(samples: pd.DataFrame, source: pd.DataFrame, qc: pd.DataFrame, balance: pd.DataFrame, validity: pd.DataFrame, classification: str, verdict: str, source_status: str) -> None:
    nulls = pd.read_csv(NULLS, sep="\t") if NULLS.exists() else pd.DataFrame()
    module_pass = int((nulls.FDR < 0.05).sum()) if not nulls.empty else 0
    dar = pd.read_csv(DAR, sep="\t")
    report = f"""# Phase 2R.3 — final ATAC QC recovery and project closure

Date: {time.strftime('%Y-%m-%d')}
Project: PIEZO1/TRPV4 context-dependent mechanotransduction pathway programs in glaucoma

## Final phase verdict

**{verdict}**

ATAC classification: **{classification}**

## Raw-source audit and recovery

- Samples audited: {len(samples)}/6.
- Initial local source status: {source_status}.
- Historical `SRR34339248.phase2r3.name.bam` was not found or verifiable: `HISTORICAL_BAM_UNRECOVERABLE`.
- The retained, audited run3 paired FASTQ was registered by direct source-path mapping and used without copying, overwriting, SRA conversion or public re-download.
- Samples without retained pre-filter/raw inputs remain explicitly `FINAL_BAM_ONLY_INCOMPLETE`; their final-BAM TSS/FRiP/mapping values were not substituted for missing duplicate/chrM recovery metrics.
- Regenerated inputs and intermediates were retained for audit and safe resume; no old final BAM, SRA cache, run1/run2/run3 object or historical log was deleted.

Source audit: `results/phase2R3/ATAC_QC_RECOVERY_SOURCE.current.tsv`
Retention audit: `results/phase2R3/ATAC_QC_RECOVERY_RETENTION.tsv`

## Final technical QC

- PASS: {int(qc.QC_status.eq('PASS').sum())}; FAIL: {int(qc.QC_status.eq('FAIL').sum())}; INCOMPLETE: {int(qc.QC_status.eq('INCOMPLETE').sum())}.
- Locked thresholds: mapping rate ≥0.50; usable nuclear fragments ≥1,000,000; duplicate fraction ≤0.50; chrM fraction ≤0.20.
- Duplicate fraction recovered for {int(qc.duplicate_fraction.notna().sum())}/6; chrM fraction recovered for {int(qc.chrM_fraction.notna().sum())}/6; TSS recovered for {int(qc.tss_enrichment.notna().sum())}/6; FRiP recovered for {int(qc.frip.notna().sum())}/6.
- TSS and FRiP values are carried from the completed final-BAM technical audit; raw/pre-filter values are not substituted for them.

Final QC: `results/phase2R3/ATAC_FINAL_QC.tsv`
Balance audit: `results/phase2R3/ATAC_FINAL_QC_BALANCE.tsv`

## DAR validity

- Existing locked full-set DARs retained: {int((dar.FDR < 0.05).sum())} (opening {int(((dar.FDR < 0.05) & (dar.log2FC > 0)).sum())}; closing {int(((dar.FDR < 0.05) & (dar.log2FC < 0)).sum())}).
- No differential-accessibility rerun was used to replace the primary result.
- QC-pass sensitivity was considered only after design-rank checking; a non-identifiable subset was not analyzed.

DAR validity: `results/phase2R3/ATAC_FINAL_DAR_VALIDITY.tsv`

## Module and regulatory conclusion

- Frozen matched-null modules passing FDR<0.05: {module_pass}/9.
- Final conclusion: **NO_MODULE_SPECIFIC_REGULATORY_SUPPORT**.
- Regulatory gate: **STOP_REGULATORY_ESCALATION**.
- Genetic gate: **STOP_GENETIC_ESCALATION**.
- Raw CUT&Tag, raw Hi-C, GWAS, eQTL, colocalization, PRS and MR remain closed/prohibited.

## Manuscript route

Use the context-dependent mechanotransduction pathway-program framing. Given the final QC classification, retain ATAC only as supplementary exploratory/audit evidence and do not make it a primary mechanistic pillar. Explicitly state that module-specific regulatory enrichment is absent.

Evidence hierarchy: `reports/phase2R3/FINAL_MANUSCRIPT_EVIDENCE_HIERARCHY.md`
"""
    (REPORT / "PHASE2R3_FINAL_ATAC_QC_REPORT.md").write_text(report, encoding="utf-8")
    (REPORT / "FINAL_PROJECT_GATE.md").write_text(f"# Final project gate — Phase 2R.3\n\nFINAL_PROJECT_GATE: **{verdict}**\n\nATAC classification: **{classification}**\nHistorical name BAM: **HISTORICAL_BAM_UNRECOVERABLE**\nRun3 FASTQ fallback: **AUDITED_RUN3_FASTQ_REBUILD**\nQC complete table: **{len(qc)}/6 samples**\nMatched-null modules passing FDR<0.05: **{module_pass}/9**\nNO_MODULE_SPECIFIC_REGULATORY_SUPPORT\nSTOP_REGULATORY_ESCALATION\nSTOP_GENETIC_ESCALATION\nCUT&Tag raw rescue: NO\nRaw Hi-C rescue: NO\nGWAS/eQTL/colocalization/PRS/MR: CLOSED\nRecommended next phase: manuscript phase; no regulatory or genetic escalation\n", encoding="utf-8")


def main() -> None:
    log(f"PHASE2R3_START run_id={RUN_ID} free_gib={free_gib():.2f}")
    samples = pd.read_csv(SAMPLE_MAP, sep="\t", dtype=str)
    source = audit_sources(samples)
    source_status = ";".join(sorted(source.local_source_status.unique()))
    old_qc = load_completed_qc()
    recovered_rows: list[dict[str, object]] = []
    retention_rows: list[dict[str, object]] = []
    checkpoint = PARTIAL if PARTIAL.exists() else LEGACY_PARTIAL
    if checkpoint.exists():
        try:
            prior = pd.read_csv(checkpoint, sep="\t", dtype={"SRR": str})
            recovered_rows.extend(prior.loc[prior["status"].isin(["RECOVERED", "INCOMPLETE_FINAL_BAM_ONLY"])].to_dict("records"))
            log(f"RESUME checkpoint={checkpoint.name} recovered_checkpoint_samples={len(recovered_rows)}")
        except Exception as exc:
            log(f"CHECKPOINT_READ_FAILED {checkpoint}: {exc}")
    done = {str(r["SRR"]) for r in recovered_rows}
    for _, row in samples.iterrows():
        srr = str(row.SRR)
        if srr in done:
            continue
        log(f"SAMPLE_START {srr}")
        sra = fq1 = fq2 = None
        try:
            sra, fq1, fq2, source_kind = obtain_fastq(row)
            metrics = fastp_and_align(row, sra, fq1, fq2)
            metrics["recovery_source"] = source_kind
            if source_kind == "AUDITED_RUN3_FASTQ_REBUILD":
                validation = validate_regenerated_bam(srr, metrics)
                metrics["regenerated_bam_validation"] = "PASS" if validation.status.eq("PASS").all() else "FAIL"
            recovered_rows = [r for r in recovered_rows if str(r.get("SRR")) != srr] + [metrics]
            write_df(pd.DataFrame(recovered_rows), PARTIAL)
            retention_rows.extend(cleanup_sample(srr))
            write_df(pd.DataFrame(retention_rows), P2 / "ATAC_QC_RECOVERY_RETENTION.partial.tsv")
            log(f"SAMPLE_DONE {srr} mapping={metrics['mapping_rate_raw_alignment']:.6f} dup={metrics['duplicate_fraction']:.6f} chrM={metrics['chrM_fraction']:.6f}")
        except RuntimeError as exc:
            if str(exc).startswith("NO_RETAINED_PAIRED_FASTQ_NO_PUBLIC_REACQUISITION_ALLOWED:"):
                metrics = final_bam_only_metrics(row, old_qc)
                recovered_rows = [r for r in recovered_rows if str(r.get("SRR")) != srr] + [metrics]
                write_df(pd.DataFrame(recovered_rows), PARTIAL)
                log(f"SAMPLE_INCOMPLETE_FINAL_BAM_ONLY {srr}: no retained raw/pre-filter input; no public reacquisition")
                continue
            log(f"SAMPLE_FAILED {srr}: {type(exc).__name__}: {exc}")
            raise
        except Exception as exc:
            log(f"SAMPLE_FAILED {srr}: {type(exc).__name__}: {exc}")
            # Keep all inputs/intermediates for safe resume; do not delete a failed sample.
            raise
    recovered = pd.DataFrame(recovered_rows)
    qc = merge_final_qc(samples, recovered, old_qc)
    balance = qc_balance(qc)
    validity = dar_validity(samples, qc)
    retention(retention_rows, source, recovered)
    if (qc.QC_status == "PASS").all():
        classification = "FORMALLY_USABLE_EXPLORATORY_ATAC"
    elif int(qc.QC_status.eq("PASS").sum()) >= 4 and int(qc.QC_status.eq("INCOMPLETE").sum()) == 0:
        classification = "LIMITED_EXPLORATORY_ATAC"
    else:
        classification = "TECHNICALLY_UNRELIABLE_ATAC"
    verdict = classification
    update_project_files(verdict, classification, source_status)
    manuscript_hierarchy(classification)
    final_report(samples, source, qc, balance, validity, classification, verdict, source_status)
    # Keep the completed checkpoint as an auditable recovery object; it is not
    # the final QC table and is never used as a final result.
    write_df(pd.DataFrame(recovered_rows), P2 / "ATAC_QC_RECOVERY_METRICS.tsv")
    print("PHASE 2R.3 VERDICT")
    print(verdict)
    print(f"ATAC CLASSIFICATION: {classification}")
    print(f"QC PASS / FAIL / INCOMPLETE: {int(qc.QC_status.eq('PASS').sum())} / {int(qc.QC_status.eq('FAIL').sum())} / {int(qc.QC_status.eq('INCOMPLETE').sum())}")
    dar_count = int(pd.read_csv(DAR, sep="\t").FDR.lt(0.05).sum())
    print(f"DARs RETAINED: {dar_count}")
    print("MODULE-SPECIFIC REGULATORY SUPPORT: NO_MODULE_SPECIFIC_REGULATORY_SUPPORT")
    print("TARGETED CUT&Tag RESCUE: NO")
    print("RAW Hi-C RESCUE: NO")
    print("GENETIC ESCALATION: STOP_GENETIC_ESCALATION")
    print("PHASE 3: CLOSED")
    print(f"REPORT: {REPORT / 'PHASE2R3_FINAL_ATAC_QC_REPORT.md'}")


if __name__ == "__main__":
    main()
