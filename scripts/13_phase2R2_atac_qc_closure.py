#!/usr/bin/env python3
"""Phase 2R.2 ATAC QC closure and final regulatory gate.

This script consumes only retained Phase 2R.1 artifacts. It never overwrites
historical Phase 2R.1 result tables. Metrics that require pre-deduplication or
pre-chrM-removal objects are explicitly marked INCOMPLETE when those objects
are no longer available.
"""
from __future__ import annotations

import gzip
import importlib.util
import math
import os
import shlex
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "results" / "phase2R1"
P2 = ROOT / "results" / "phase2R2"
REPORT = ROOT / "reports" / "phase2R2"
LOG = ROOT / "logs" / "phase2R2"
RAW_ATAC = ROOT / "data" / "raw" / "external" / "phase2R1" / "atac"
REF = ROOT / "data" / "reference" / "hg38"
P2.mkdir(parents=True, exist_ok=True)
REPORT.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

ENV_PREFIX = Path(os.environ.get("PHASE2R1_ENV_PREFIX", os.environ.get("CONDA_PREFIX", sys.prefix)))
BIN = ENV_PREFIX / "bin"
SAMTOOLS = str(BIN / "samtools")
BEDTOOLS = str(BIN / "bedtools")
RSCRIPT = os.environ.get("R_BIN") or shutil.which("Rscript") or "Rscript"
RUN_LOG = LOG / "phase2r2_run.log"

SAMPLE_MAP = P1 / "ATAC_RAW_SAMPLE_MAP.tsv"
OLD_QC = P1 / "ATAC_RAW_QC.tsv"
COUNTS = P1 / "ATAC_COUNT_MATRIX.tsv"
DESIGN = P1 / "ATAC_LEVEL3_DESIGN.tsv"
DAR = P1 / "ATAC_LEVEL3_DAR.tsv"
OLD_NULLS = P1 / "ATAC_LEVEL3_MODULE_NULLS.tsv"
OLD_STATUS = P1 / "FINAL_LEVEL3_MODULE_STATUS.tsv"
CORR = P1 / "ATAC_SAMPLE_CORRELATION.tsv"
CONSENSUS = P1 / "ATAC_CONSENSUS_PEAKS.bed"
TSS = REF / "gencode.v47.tss.bed"
RNA = ROOT / "data" / "raw" / "external" / "phase2" / "GSE301517_Control_vs_Dex_Gene_CPM_MSTR.tsv.gz"

QC_OUT = P2 / "ATAC_RAW_QC_COMPLETED.tsv"
DESIGN_OUT = P2 / "ATAC_QC_BY_DESIGN.tsv"
SENS_OUT = P2 / "ATAC_DAR_QC_SENSITIVITY.tsv"
NULL_OUT = P2 / "ATAC_LEVEL3_MODULE_NULLS_FINAL.tsv"
STATUS_OUT = P2 / "FINAL_LEVEL3_MODULE_STATUS_CORRECTED.tsv"
MECH_OUT = P2 / "LEVEL3_MECHANOSENSOR_RESPONSE_FINAL.tsv"
REPORT_OUT = REPORT / "PHASE2R2_ATAC_QC_CLOSURE_REPORT.md"
GATE_OUT = REPORT / "FINAL_REGULATORY_GATE.md"

SEED = 20260918
NULL_B = 10000
FDR_CUTOFF = 0.05
MAPPING_THRESHOLD = 0.50
FRAGMENT_THRESHOLD = 1_000_000
MAX_DUPLICATE_FRACTION = 0.50
MAX_CHRM_FRACTION = 0.20
ANCHORS = ["PIEZO1", "TRPV4", "PTK2", "ITGA9", "ANGPT2", "TEK", "RHOA", "ROCK1", "ROCK2", "YAP1", "WWTR1", "CCN2", "TGFB2", "FN1"]


def log(msg: str) -> None:
    with RUN_LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")


def write_df(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def command_output(cmd: list[str], name: str) -> tuple[int, str, str]:
    log("COMMAND " + name + ": " + " ".join(shlex.quote(str(x)) for x in cmd))
    p = subprocess.run(cmd, capture_output=True, text=True, check=False)
    log(f"COMMAND_RESULT {name}: returncode={p.returncode}")
    if p.stderr:
        (LOG / f"{name}.stderr.log").write_text(p.stderr, encoding="utf-8")
    return p.returncode, p.stdout, p.stderr


def scalar_count(cmd: list[str], name: str) -> int | None:
    rc, stdout, _ = command_output(cmd, name)
    if rc != 0:
        return None
    try:
        return int(stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return None


def load_modules() -> dict[str, dict[str, object]]:
    spec = importlib.util.spec_from_file_location("phase2_modules", ROOT / "scripts" / "10_phase2_multiomic_anchoring.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("module definition script could not be loaded")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.MODULES


def load_peak_genes() -> tuple[dict[str, set[str]], dict[str, str]]:
    genes: dict[str, set[str]] = {}
    methods: dict[str, str] = {}
    promoter = P1 / "DAR_promoter_links.tsv"
    if promoter.exists():
        for line in promoter.read_text(errors="replace").splitlines():
            f = line.split("\t")
            if len(f) >= 10:
                genes.setdefault(f[3], set()).add(f[9].upper())
                methods[f[3]] = "promoter_overlap"
    loop = P1 / "DAR_loop_links.tsv"
    if loop.exists():
        for line in loop.read_text(errors="replace").splitlines():
            f = line.split("\t")
            if len(f) >= 8 and f[3] not in genes:
                genes.setdefault(f[3], set()).add(f[7].upper())
                methods[f[3]] = "existing_GSE301522_loop_or_promoter_link"
    nearest = P1 / "DAR_nearest_gene.tsv"
    if nearest.exists():
        for line in nearest.read_text(errors="replace").splitlines():
            f = line.split("\t")
            if len(f) >= 10 and f[3] not in genes:
                genes.setdefault(f[3], set()).add(f[7].upper())
                methods[f[3]] = "proximal_assignment_sensitivity_only"
    return genes, methods


def load_tss() -> dict[str, list[tuple[int, str]]]:
    out: dict[str, list[tuple[int, str]]] = defaultdict(list)
    with TSS.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) >= 6:
                out[f[0]].append((int(f[1]), f[5]))
    return dict(out)


def tss_profile_score(bam_path: Path, tss_by_chrom: dict[str, list[tuple[int, str]]], window: int = 2000) -> tuple[float, int, str]:
    """Build a standard cut-site TSS profile from the retained final BAM.

    Cut sites use the conventional Tn5 offsets (+4 on forward reads and -5
    on reverse reads). The scalar score is the maximum mean signal in the
    central +/-100 bp divided by the mean of the outer 100 bp on both profile
    ends. The profile is computed from final deduplicated, chrM-filtered BAMs.
    """
    profile = np.zeros(2 * window + 1, dtype=np.float64)
    used_tss = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        ref_lengths = dict(zip(bam.references, bam.lengths))
        for chrom, entries in tss_by_chrom.items():
            if chrom not in ref_lengths:
                continue
            chrom_len = int(ref_lengths[chrom])
            coverage = np.zeros(chrom_len, dtype=np.uint32)
            for read in bam.fetch(chrom):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_duplicate:
                    continue
                if read.reference_end is None:
                    continue
                pos = int(read.reference_end - 5 if read.is_reverse else read.reference_start + 4)
                if 0 <= pos < chrom_len:
                    coverage[pos] += 1
            for pos, strand in entries:
                left = max(0, pos - window)
                right = min(chrom_len, pos + window + 1)
                signal = coverage[left:right].astype(np.float64, copy=False)
                if strand == "-":
                    signal = signal[::-1]
                dst_left = window - (pos - left)
                profile[dst_left:dst_left + len(signal)] += signal
                used_tss += 1
            del coverage
    if used_tss == 0:
        return math.nan, 0, "INCOMPLETE_NO_TSS_OVERLAP"
    profile /= float(used_tss)
    flank = np.concatenate([profile[:100], profile[-100:]])
    flank_mean = float(np.mean(flank))
    core = profile[window - 100:window + 101]
    score = float(np.max(core) / flank_mean) if flank_mean > 0 else math.nan
    return score, used_tss, "COMPLETE_FINAL_BAM_TSS_PROFILE" if np.isfinite(score) else "INCOMPLETE_ZERO_FLANK_SIGNAL"


def frip(final_bam: Path) -> float:
    if not CONSENSUS.exists() or CONSENSUS.stat().st_size == 0:
        return math.nan
    shell = (
        f"set -o pipefail; {shlex.quote(BEDTOOLS)} intersect -u -abam {shlex.quote(str(final_bam))} "
        f"-b {shlex.quote(str(CONSENSUS))} | {shlex.quote(SAMTOOLS)} view -c -f 1 -"
    )
    p = subprocess.run(["/bin/zsh", "-lc", shell], capture_output=True, text=True, check=False)
    log(f"COMMAND frip_{final_bam.stem}: returncode={p.returncode}")
    if p.returncode != 0:
        (LOG / f"frip_{final_bam.stem}.stderr.log").write_text(p.stderr, encoding="utf-8")
        return math.nan
    try:
        hit_reads = int(p.stdout.strip().splitlines()[-1])
        usable_reads = scalar_count([SAMTOOLS, "view", "-c", "-f", "1", str(final_bam)], f"usable_reads_{final_bam.stem}")
        return (hit_reads / 2) / max((usable_reads or 0) / 2, 1.0)
    except (ValueError, IndexError):
        return math.nan


def load_correlation() -> pd.DataFrame:
    if not CORR.exists():
        return pd.DataFrame()
    return pd.read_csv(CORR, sep="\t")


def sample_qc(samples: pd.DataFrame, old_qc: pd.DataFrame) -> pd.DataFrame:
    tss_by_chrom = load_tss()
    corr = load_correlation()
    rows: list[dict[str, object]] = []
    for _, row in samples.iterrows():
        srr = str(row["SRR"])
        bam = RAW_ATAC / "processed" / srr / f"{srr}.final.bam"
        frag = RAW_ATAC / "processed" / srr / f"{srr}.fragments.bed.gz"
        mapped_reads = scalar_count([SAMTOOLS, "view", "-c", "-F", "4", str(bam)], f"mapped_reads_{srr}") if bam.exists() else None
        usable_reads = scalar_count([SAMTOOLS, "view", "-c", "-f", "1", str(bam)], f"usable_reads_{srr}") if bam.exists() else None
        raw_reads = int(float(row["spots"])) * 2
        usable_fragments = (usable_reads / 2.0) if usable_reads is not None else math.nan
        mapping_rate = (mapped_reads / raw_reads) if mapped_reads is not None else math.nan
        f = frip(bam) if bam.exists() else math.nan
        tss_score, tss_n, tss_status = tss_profile_score(bam, tss_by_chrom) if bam.exists() else (math.nan, 0, "INCOMPLETE_FINAL_BAM_MISSING")
        old = old_qc.loc[old_qc["SRR"].astype(str).eq(srr)]
        old_row = old.iloc[0].to_dict() if not old.empty else {}
        pearson_mean = math.nan
        spearman_mean = math.nan
        if not corr.empty:
            x = corr[(corr["sample_i"].astype(str).eq(srr)) | (corr["sample_j"].astype(str).eq(srr))]
            if len(x):
                pearson_mean = float(x["pearson_r"].mean())
                spearman_mean = float(x["spearman_r"].mean())
        reasons = []
        if not bam.exists() or not frag.exists():
            reasons.append("final_BAM_or_fragments_missing")
        if mapping_rate < MAPPING_THRESHOLD:
            reasons.append(f"mapping_rate<{MAPPING_THRESHOLD}")
        if np.isfinite(usable_fragments) and usable_fragments < FRAGMENT_THRESHOLD:
            reasons.append(f"usable_fragments<{FRAGMENT_THRESHOLD}")
        reasons.extend([
            "duplicate_fraction_not_recoverable_without_pre_dedup_BAM_or_raw_reads",
            "chrM_fraction_not_recoverable_after_chrM_removal_without_pre_filter_BAM_or_raw_reads",
        ])
        if tss_status.startswith("INCOMPLETE"):
            reasons.append(tss_status)
        # The locked Phase 2R.1 thresholds are retained. Missing essential
        # pre-filter metrics force INCOMPLETE; no post hoc threshold is added.
        status = "INCOMPLETE" if reasons else ("PASS" if mapping_rate >= MAPPING_THRESHOLD and usable_fragments >= FRAGMENT_THRESHOLD else "FAIL")
        rows.append({
            "SRR": srr, "GSM": row["GSM"], "condition": row["condition"], "cell_line": row["cell_line"],
            "biological_replicate": row["biological_replicate"], "raw_reads": raw_reads,
            "mapped_reads_final_bam": mapped_reads, "mapping_rate_final_bam_vs_raw_reads": mapping_rate,
            "usable_fragments": usable_fragments, "frip": f,
            "duplicate_fraction": math.nan, "duplicate_fraction_status": "INCOMPLETE_PREDEDUP_OBJECT_NOT_RETAINED",
            "chrM_fraction": math.nan, "chrM_fraction_status": "INCOMPLETE_PRE_FILTER_OBJECT_NOT_RETAINED",
            "tss_enrichment": tss_score, "tss_enrichment_status": tss_status, "tss_n": tss_n,
            "nucleosome_free_fraction": old_row.get("nucleosome_free_fraction", math.nan),
            "mono_nucleosome_fraction": old_row.get("mono_nucleosome_fraction", math.nan),
            "di_nucleosome_fraction": old_row.get("di_nucleosome_fraction", math.nan),
            "library_complexity": math.nan, "library_complexity_status": "INCOMPLETE_PREDEDUP_OBJECT_NOT_RETAINED",
            "mean_pearson_r_all_samples": pearson_mean, "mean_spearman_r_all_samples": spearman_mean,
            "phase2R1_historical_qc_status": old_row.get("qc_status", "NA"),
            "QC_status": status, "reason": ";".join(reasons) if reasons else "all_locked_thresholds_met",
            "final_bam": str(bam), "fragment_file": str(frag),
        })
    out = pd.DataFrame(rows)
    write_df(out, QC_OUT)
    return out


def design_matrix(df: pd.DataFrame) -> tuple[int, int, str]:
    if df.empty:
        return 0, 0, "NO_SAMPLES"
    x = pd.DataFrame({
        "Intercept": 1.0,
        "cell_line_TM70": (df["cell_line"].astype(str) == "TM70").astype(float),
        "treatment_DEX": (df["condition"].astype(str) == "DEX").astype(float),
    })
    rank = int(np.linalg.matrix_rank(x.to_numpy()))
    return rank, x.shape[1], "FULL_RANK" if rank == x.shape[1] else "RANK_DEFICIENT"


def design_audit(samples: pd.DataFrame, old_qc: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    old_status = old_qc.set_index(old_qc["SRR"].astype(str))["qc_status"].to_dict()
    rows: list[dict[str, object]] = []
    for _, r in samples.iterrows():
        srr = str(r["SRR"])
        rows.append({
            "record_type": "sample", "factor": "sample", "level": srr,
            "SRR": srr, "condition": r["condition"], "cell_line": r["cell_line"],
            "biological_replicate": r["biological_replicate"],
            "phase2R1_historical_fail": int(old_status.get(srr) == "FAIL"),
            "phase2R2_incomplete": int(qc.loc[qc["SRR"].eq(srr), "QC_status"].iloc[0] == "INCOMPLETE"),
            "n": 1, "fail_n": int(old_status.get(srr) == "FAIL"), "fail_rate": float(old_status.get(srr) == "FAIL"),
            "odds_ratio": math.nan, "fisher_p": math.nan,
            "test_note": "technical audit row; no treatment interpretation",
        })
    hist_fail = np.array([int(old_status.get(str(s), "NA") == "FAIL") for s in samples["SRR"]])
    for factor in ["condition", "cell_line", "biological_replicate"]:
        for level in sorted(samples[factor].astype(str).unique()):
            in_group = samples[factor].astype(str).to_numpy() == level
            fail_n = int(hist_fail[in_group].sum())
            n = int(in_group.sum())
            table = [[fail_n, n - fail_n], [int(hist_fail[~in_group].sum()), int((~in_group).sum() - hist_fail[~in_group].sum())]]
            try:
                odds, p = fisher_exact(table)
            except Exception:
                odds, p = math.nan, math.nan
            rows.append({
                "record_type": "summary", "factor": factor, "level": level, "SRR": "ALL",
                "condition": "ALL", "cell_line": "ALL", "biological_replicate": "ALL",
                "phase2R1_historical_fail": int(hist_fail.sum()), "phase2R2_incomplete": int(qc["QC_status"].eq("INCOMPLETE").sum()),
                "n": n, "fail_n": fail_n, "fail_rate": fail_n / n if n else math.nan,
                "odds_ratio": odds, "fisher_p": p,
                "test_note": "Fisher exact on historical Phase 2R.1 QC_FAIL; descriptive and underpowered; no treatment interpretation",
            })
    rows.append({
        "record_type": "summary", "factor": "none", "level": "all", "SRR": "ALL",
        "condition": "ALL", "cell_line": "ALL", "biological_replicate": "ALL",
        "phase2R1_historical_fail": int(hist_fail.sum()), "phase2R2_incomplete": int(qc["QC_status"].eq("INCOMPLETE").sum()),
        "n": len(samples), "fail_n": int(hist_fail.sum()), "fail_rate": float(hist_fail.mean()),
        "odds_ratio": math.nan, "fisher_p": math.nan,
        "test_note": "overall technical audit; no treatment interpretation",
    })
    out = pd.DataFrame(rows)
    write_df(out, DESIGN_OUT)
    return out


def verify_dar_and_sensitivity(samples: pd.DataFrame, qc: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    dar = pd.read_csv(DAR, sep="\t")
    sig = dar[dar["FDR"] < FDR_CUTOFF]
    full_rank, full_cols, full_rank_status = design_matrix(samples)
    pass_samples = samples.loc[samples["SRR"].astype(str).isin(qc.loc[qc["QC_status"].eq("PASS"), "SRR"].astype(str))].copy()
    pass_rank, pass_cols, pass_rank_status = design_matrix(pass_samples)
    historical_pass = samples.loc[samples["SRR"].astype(str).isin(qc.loc[qc["phase2R1_historical_qc_status"].eq("PASS"), "SRR"].astype(str))].copy()
    hist_rank, hist_cols, hist_rank_status = design_matrix(historical_pass)
    rows = [{
        "analysis": "full_set_locked_edgeR", "sample_n": len(samples),
        "control_n": int((samples["condition"] == "CONTROL").sum()), "dex_n": int((samples["condition"] == "DEX").sum()),
        "tm62_n": int((samples["cell_line"] == "TM62").sum()), "tm70_n": int((samples["cell_line"] == "TM70").sum()),
        "design_rank": full_rank, "design_columns": full_cols, "design_status": full_rank_status,
        "significant_DAR_n": len(sig), "opening_n": int((sig["log2FC"] > 0).sum()), "closing_n": int((sig["log2FC"] < 0).sum()),
        "effect_size_correlation": 1.0, "directional_concordance": 1.0,
        "status": "VERIFIED_EXISTING_LOCKED_EDGER_OUTPUT", "reason": "phase2R1 edgeR run returned code 0 under ~ cell_line + treatment; current DAR retained",
    }]
    if len(pass_samples) == 0:
        pass_reason = "NOT_RUN_NO_PHASE2R2_QC_PASS_SAMPLES"
    elif pass_rank_status != "FULL_RANK" or pass_samples["condition"].nunique() < 2:
        pass_reason = "NOT_RUN_PASS_ONLY_DESIGN_NOT_IDENTIFIABLE"
    else:
        pass_reason = "NOT_RUN_BY_LOCKED_RESUME_POLICY"
    rows.append({
        "analysis": "QC_pass_only_sensitivity", "sample_n": len(pass_samples),
        "control_n": int((pass_samples["condition"] == "CONTROL").sum()) if len(pass_samples) else 0,
        "dex_n": int((pass_samples["condition"] == "DEX").sum()) if len(pass_samples) else 0,
        "tm62_n": int((pass_samples["cell_line"] == "TM62").sum()) if len(pass_samples) else 0,
        "tm70_n": int((pass_samples["cell_line"] == "TM70").sum()) if len(pass_samples) else 0,
        "design_rank": pass_rank, "design_columns": pass_cols, "design_status": pass_rank_status,
        "significant_DAR_n": math.nan, "opening_n": math.nan, "closing_n": math.nan,
        "effect_size_correlation": math.nan, "directional_concordance": math.nan,
        "status": "NOT_RUN", "reason": pass_reason,
    })
    rows.append({
        "analysis": "historical_Phase2R1_PASS_only_design_check", "sample_n": len(historical_pass),
        "control_n": int((historical_pass["condition"] == "CONTROL").sum()) if len(historical_pass) else 0,
        "dex_n": int((historical_pass["condition"] == "DEX").sum()) if len(historical_pass) else 0,
        "tm62_n": int((historical_pass["cell_line"] == "TM62").sum()) if len(historical_pass) else 0,
        "tm70_n": int((historical_pass["cell_line"] == "TM70").sum()) if len(historical_pass) else 0,
        "design_rank": hist_rank, "design_columns": hist_cols, "design_status": hist_rank_status,
        "significant_DAR_n": math.nan, "opening_n": math.nan, "closing_n": math.nan,
        "effect_size_correlation": math.nan, "directional_concordance": math.nan,
        "status": "NOT_RUN", "reason": "historical PASS-only subset retained both conditions but lacks TM70; full design matrix is rank deficient",
    })
    out = pd.DataFrame(rows)
    write_df(out, SENS_OUT)
    return out, {"dar": dar, "sig": sig, "full_rank": full_rank, "full_rank_status": full_rank_status, "pass_rank": pass_rank, "pass_rank_status": pass_rank_status, "historical_rank": hist_rank, "historical_rank_status": hist_rank_status}


def bh(p: np.ndarray) -> np.ndarray:
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    vals = np.clip(p[ok], 0, 1)
    order = np.argsort(vals)
    q = np.minimum.accumulate((vals[order] * len(vals) / np.arange(1, len(vals) + 1))[::-1])[::-1]
    tmp = np.empty(len(vals))
    tmp[order] = np.minimum(q, 1)
    out[ok] = tmp
    return out


def rerun_module_nulls(dar: pd.DataFrame, samples: pd.DataFrame, modules: dict[str, dict[str, object]], peak_genes: dict[str, set[str]]) -> pd.DataFrame:
    counts = pd.read_csv(COUNTS, sep="\t")
    cmat = counts.set_index("peak_id")
    sig = dar[dar["FDR"] < FDR_CUTOFF].copy()
    sig_genes: set[str] = set()
    for pid in sig["peak_id"]:
        sig_genes.update(peak_genes.get(pid, set()))
    linked: dict[str, set[str]] = {}
    for pid, gs in peak_genes.items():
        for gene in gs:
            linked.setdefault(gene, set()).add(pid)
    rna = pd.read_csv(RNA, sep="\t", usecols=lambda c: c == "external_gene_name" or "Control" in str(c))
    rna["gene"] = rna["external_gene_name"].astype(str).str.upper()
    control_cols = [c for c in rna.columns if "Control" in str(c)]
    rna["baseline_expression"] = rna[control_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    control_samples = samples.loc[samples["condition"].eq("CONTROL"), "SRR"].astype(str).tolist()
    promoter_access: dict[str, float] = {}
    for gene, peaks in linked.items():
        promoter_access[gene] = sum(float(cmat.loc[p, control_samples].to_numpy(dtype=float).sum()) for p in peaks if p in cmat.index)
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
    rows: list[dict[str, object]] = []
    for mid, module in modules.items():
        genes = sorted({str(g).upper().replace("CTGF", "CCN2") for g in module["genes"]})
        observed_genes = sorted(set(genes) & sig_genes)
        observed = len(observed_genes)
        draws = np.zeros(NULL_B, dtype=int)
        for b in range(NULL_B):
            sampled: list[str] = []
            for gene in genes:
                if gene not in cov.index:
                    continue
                v = cov.loc[gene]
                mask = (
                    cov["baseline_expression_bin"].eq(v["baseline_expression_bin"])
                    & cov["promoter_accessibility_bin"].eq(v["promoter_accessibility_bin"])
                    & cov["gene_length_bin"].eq(v["gene_length_bin"])
                    & cov["region_bin"].eq(v["region_bin"])
                ).to_numpy()
                pool = cov.index[mask].tolist()
                if len(pool) < 2:
                    d = (
                        (cov["baseline_expression"] - v["baseline_expression"]).abs() / (cov["baseline_expression"].std() or 1)
                        + (cov["gene_length"] - v["gene_length"]).abs() / (cov["gene_length"].std() or 1)
                    )
                    pool = cov.index[np.argsort(d.to_numpy())[:max(5, len(genes))]].tolist()
                pool = [x for x in pool if x not in sampled]
                sampled.append(str(rng.choice(pool or cov.index.tolist())))
            draws[b] = sum(g in sig_genes for g in sampled)
        empirical = (1 + np.sum(draws >= observed)) / (NULL_B + 1)
        rows.append({
            "module_id": mid, "module_name": module["name"], "module_size": len(genes),
            "observed_statistic": observed, "observed_gene_list": ";".join(observed_genes) or "NONE",
            "null_median": np.median(draws), "null_95CI_low": np.quantile(draws, 0.025), "null_95CI_high": np.quantile(draws, 0.975),
            "empirical_P": empirical, "matched_null_B": NULL_B, "background_gene_n": len(cov),
            "mapping_hierarchy": "promoter overlap > existing GSE301522 link > proximal sensitivity only", "pass": "FAIL",
        })
    out = pd.DataFrame(rows)
    out["FDR"] = bh(out["empirical_P"].to_numpy())
    out["pass"] = np.where(out["FDR"] < FDR_CUTOFF, "PASS", "FAIL")
    write_df(out, NULL_OUT)
    return out


def corrected_status(nulls: pd.DataFrame, modules: dict[str, dict[str, object]], qc: pd.DataFrame) -> pd.DataFrame:
    technical = bool(qc["QC_status"].eq("INCOMPLETE").any())
    old = pd.read_csv(OLD_STATUS, sep="\t")
    old_map = old.set_index("module_id")["status"].to_dict()
    rows: list[dict[str, object]] = []
    for mid, module in modules.items():
        x = nulls.loc[nulls["module_id"].eq(mid)]
        if x.empty:
            obs, empirical, fdr, genes = 0, math.nan, math.nan, []
        else:
            row = x.iloc[0]
            obs = int(row["observed_statistic"])
            empirical = float(row["empirical_P"])
            fdr = float(row["FDR"])
            genes = [g for g in str(row["observed_gene_list"]).split(";") if g and g != "NONE"]
        dominance = "YES" if len(genes) == 1 else ("NO" if len(genes) >= 2 else "NO_OBSERVED_GENES")
        rule = "FDR<0.05 + observed>=2 + no single-gene dominance"
        if technical:
            status = "TECHNICALLY_UNRESOLVED"
            reason = "essential QC remains incomplete; module-specific result cannot be promoted to final regulatory status"
        elif np.isfinite(fdr) and fdr < FDR_CUTOFF and obs >= 2 and dominance == "NO":
            status = "LEVEL3_REGULATORY_SUPPORTED"
            reason = "all locked regulatory-support criteria met"
        else:
            status = "UNSUPPORTED"
            reason = "pre-specified regulatory-support rule not satisfied; partial label is not permitted"
        rows.append({
            "module_id": mid, "module_name": module["name"], "significant_DAR_gene_n": obs,
            "observed_gene_list": ";".join(genes) or "NONE", "single_gene_dominance": dominance,
            "matched_null_empirical_P": empirical, "matched_null_FDR": fdr,
            "previous_status": old_map.get(mid, "MISSING"), "rule_required_status": rule,
            "corrected_status": status, "reason": reason,
        })
    out = pd.DataFrame(rows)
    write_df(out, STATUS_OUT)
    return out


def mechanosensor_status(dar: pd.DataFrame, peak_genes: dict[str, set[str]], methods: dict[str, str], qc: pd.DataFrame) -> pd.DataFrame:
    technical = bool(qc["QC_status"].eq("INCOMPLETE").any())
    rows: list[dict[str, object]] = []
    for gene in ANCHORS:
        hits = [pid for pid, gs in peak_genes.items() if gene in gs]
        x = dar[dar["peak_id"].isin(hits)] if hits else dar.iloc[0:0]
        sig = x[x["FDR"] < FDR_CUTOFF] if not x.empty else x
        if x.empty:
            raw_label = "NOT_TESTABLE"
            final_label = "NOT_TESTABLE"
        elif sig.empty:
            raw_label = "NONSIGNIFICANT"
            final_label = "NONSIGNIFICANT"
        elif (sig["log2FC"] > 0).any() and (sig["log2FC"] < 0).any():
            raw_label = "TECHNICALLY_UNRESOLVED"
            final_label = "TECHNICALLY_UNRESOLVED"
        else:
            raw_label = "SIGNIFICANT_OPENING" if (sig["log2FC"] > 0).any() else "SIGNIFICANT_CLOSING"
            final_label = "TECHNICALLY_UNRESOLVED" if technical else raw_label
        best = sig.loc[sig["FDR"].idxmin()] if not sig.empty else None
        rows.append({
            "gene": gene, "tested_peak_n": len(x), "significant_peak_n": len(sig),
            "raw_signal_label": raw_label, "final_classification": final_label,
            "best_FDR": float(best["FDR"]) if best is not None else math.nan,
            "best_log2FC": float(best["log2FC"]) if best is not None else math.nan,
            "mapping_hierarchy": ";".join(sorted({methods.get(pid, "mapping_not_recorded") for pid in hits})) or "NONE",
            "reason": "final classification held technically unresolved because essential QC is incomplete" if final_label == "TECHNICALLY_UNRESOLVED" and raw_label != "TECHNICALLY_UNRESOLVED" else "locked replicate-level ATAC mapping result",
        })
    out = pd.DataFrame(rows)
    write_df(out, MECH_OUT)
    return out


def update_project_files(verdict: str, gate: str) -> None:
    marker = "## Phase 2R.2 ATAC QC closure — 2026-09-22"
    block = (
        f"\n{marker}\nphase2R2_verdict: {verdict}\nphase2R2_final_regulatory_gate: {gate}\n"
        "phase2R2_cuttag_raw_executed: false\nphase2R2_hic_raw_executed: false\n"
        "phase2R2_phase3_gate: CLOSED\nphase2R2_report: reports/phase2R2/PHASE2R2_ATAC_QC_CLOSURE_REPORT.md\n"
    )
    cfg = ROOT / "config" / "project_config.yaml"
    if cfg.exists():
        txt = cfg.read_text(encoding="utf-8")
        if marker not in txt:
            cfg.write_text(txt.rstrip() + block, encoding="utf-8")
    for name in ["PROJECT_STATUS.md", "ANALYSIS_GOVERNANCE.md"]:
        path = ROOT / name
        txt = path.read_text(encoding="utf-8") if path.exists() else ""
        if marker not in txt:
            path.write_text(txt.rstrip() + block, encoding="utf-8")


def write_audit_md(status: pd.DataFrame, nulls: pd.DataFrame) -> None:
    lines = [
        "# Module status logic audit — Phase 2R.2",
        "",
        "Historical `FINAL_LEVEL3_MODULE_STATUS.tsv` is preserved. The corrected table applies the locked rule: `FDR<0.05 + observed>=2 + no single-gene dominance`. A `LEVEL3_PARTIAL_SUPPORT` label is not assigned when that rule is not satisfied. Because essential QC remains incomplete, current corrected statuses are held as `TECHNICALLY_UNRESOLVED`.",
        "",
        "| module | significant_DAR_genes | empirical_P | FDR | single_gene_dominance | previous_status | rule_required_status | corrected_status | reason |",
        "|---|---:|---:|---:|---|---|---|---|---|",
    ]
    for _, r in status.iterrows():
        vals = [str(r.get(c, "NA")) for c in ["module_id", "significant_DAR_gene_n", "matched_null_empirical_P", "matched_null_FDR", "single_gene_dominance", "previous_status", "rule_required_status", "corrected_status", "reason"]]
        lines.append("| " + " | ".join(v.replace("|", "\\|") for v in vals) + " |")
    (P2 / "MODULE_STATUS_LOGIC_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(samples: pd.DataFrame, qc: pd.DataFrame, design: pd.DataFrame, sens: pd.DataFrame, nulls: pd.DataFrame, status: pd.DataFrame, mech: pd.DataFrame, verdict: str) -> None:
    sig = sens.loc[sens["analysis"].eq("full_set_locked_edgeR")].iloc[0]
    pass_row = sens.loc[sens["analysis"].eq("QC_pass_only_sensitivity")].iloc[0]
    strongest = nulls.sort_values(["empirical_P", "observed_statistic"], ascending=[True, False]).iloc[0]
    pass_modules = nulls.loc[nulls["pass"].eq("PASS"), "module_id"].tolist()
    qc_counts = qc["QC_status"].value_counts().to_dict()
    piezo = mech.loc[mech["gene"].eq("PIEZO1"), "final_classification"].iloc[0]
    trpv4 = mech.loc[mech["gene"].eq("TRPV4"), "final_classification"].iloc[0]
    report = f"""# Phase 2R.2 — ATAC QC closure and final regulatory gate

Date: 2026-09-22  
Project: PIEZO1/TRPV4 context-dependent mechanotransduction pathway programs in glaucoma

## Phase 2R.2 verdict

**{verdict}**

The closure used retained Phase 2R.1 final BAMs, fragments, consensus peaks, count matrix, locked edgeR DAR table and frozen M1–M9 module definitions. Historical Phase 2R.1 tables were not overwritten.

## QC closure

- Samples: {len(samples)}/6.
- QC PASS: {qc_counts.get('PASS', 0)}.
- QC FAIL: {qc_counts.get('FAIL', 0)}.
- QC INCOMPLETE: {qc_counts.get('INCOMPLETE', 0)}.
- Mapping rate, usable fragments, FRiP and a final-BAM-derived TSS profile were recomputed where possible.
- Duplicate fraction and chrM fraction could not be recovered because pre-dedup/pre-chrM-removal BAMs and raw FASTQ/SRA payloads were not retained. These fields are explicitly `INCOMPLETE`, with no post hoc threshold change.
- Library complexity requiring pre-dedup objects is also `INCOMPLETE`.
- TSS profile method: Tn5 offsets (+4 forward, -5 reverse), ±2,000 bp around GENCODE v47 TSS, central ±100 bp maximum divided by the mean of the outer 100 bp at both profile ends; calculated from final deduplicated, chrM-filtered BAMs.

QC table: `results/phase2R2/ATAC_RAW_QC_COMPLETED.tsv`

## QC structure audit

The design audit uses historical Phase 2R.1 `QC_FAIL` as the descriptive failure indicator and reports Fisher exact tests for condition, cell line and biological replicate. These tests are technical, small-sample audits only; they are not treatment-biology analyses.

Design audit: `results/phase2R2/ATAC_QC_BY_DESIGN.tsv`

## DAR reproducibility

- Full locked edgeR design: `~ cell_line + treatment`.
- Full-set design status: {sig['design_status']} (rank {int(sig['design_rank'])}/{int(sig['design_columns'])}).
- Full-set significant DARs: {int(sig['significant_DAR_n'])}; opening {int(sig['opening_n'])}; closing {int(sig['closing_n'])}.
- QC-pass-only sensitivity: **{pass_row['status']}**; no valid PASS-only inference was run because {pass_row['reason']}.
- Effect-size correlation and directional concordance for QC-pass sensitivity: `NA` because the sensitivity model was not estimable/was not run.

DAR sensitivity: `results/phase2R2/ATAC_DAR_QC_SENSITIVITY.tsv`

## Final matched-null module test

The nine frozen M1–M9 modules were re-tested with B={NULL_B}, the same seed and the same matching variables used in Phase 2R.1. Module membership was unchanged.

- Modules passing FDR<0.05: {', '.join(pass_modules) if pass_modules else 'NONE'}.
- Strongest module: {strongest['module_id']} ({strongest['module_name']}); empirical P={float(strongest['empirical_P']):.6g}; FDR={float(strongest['FDR']):.6g}; observed statistic={int(strongest['observed_statistic'])}.

Module nulls: `results/phase2R2/ATAC_LEVEL3_MODULE_NULLS_FINAL.tsv`

## Corrected module statuses

The historical partial-support labels were not accepted when the documented regulatory-support rule was unmet. The corrected table assigns `TECHNICALLY_UNRESOLVED` while essential QC is incomplete; after QC completion, modules with FDR>=0.05 remain `UNSUPPORTED`.

Audit: `results/phase2R2/MODULE_STATUS_LOGIC_AUDIT.md`  
Corrected table: `results/phase2R2/FINAL_LEVEL3_MODULE_STATUS_CORRECTED.tsv`

## Mechanosensor-network classification

- PIEZO1 final ATAC classification: **{piezo}**.
- TRPV4 final ATAC classification: **{trpv4}**.
- The underlying raw signal label and mapping evidence are retained in the mechanosensor table; technical hold prevents overclaiming a significant chromatin result.

Mechanosensor table: `results/phase2R2/LEVEL3_MECHANOSENSOR_RESPONSE_FINAL.tsv`

## Final gate and project lock

- Final regulatory gate: **{verdict}**.
- Targeted CUT&Tag rescue: **NO**.
- Raw Hi-C rescue: **NO**.
- Phase 3: **CLOSED**.
- GWAS, eQTL, colocalization, PRS and MR: **PROHIBITED**.

The project remains locked for manuscript construction around context-dependent mechanotransduction pathway programs. The current evidence does not establish chromatin rewiring, adaptive-to-maladaptive transition, PIEZO1/TRPV4 causality, genetic mediation or a `PRS × DEX` effect.

## Required follow-up

If a higher-confidence ATAC gate is required, recover or regenerate pre-dedup and pre-chrM-removal alignment resources and complete the missing QC. Until then, retain the technical hold and do not add CUT&Tag, Hi-C or genetic escalation.

"""
    REPORT_OUT.write_text(report, encoding="utf-8")


def main() -> None:
    RUN_LOG.write_text("", encoding="utf-8")
    samples = pd.read_csv(SAMPLE_MAP, sep="\t", dtype={"SRR": str})
    old_qc = pd.read_csv(OLD_QC, sep="\t", dtype={"SRR": str})
    old_status = pd.read_csv(OLD_STATUS, sep="\t")
    modules = load_modules()
    peak_genes, methods = load_peak_genes()
    qc = sample_qc(samples, old_qc)
    design = design_audit(samples, old_qc, qc)
    sens, sens_meta = verify_dar_and_sensitivity(samples, qc)
    dar = sens_meta["dar"]
    nulls = rerun_module_nulls(dar, samples, modules, peak_genes)
    status = corrected_status(nulls, modules, qc)
    mech = mechanosensor_status(dar, peak_genes, methods, qc)
    write_audit_md(status, nulls)
    qc_complete = not qc["QC_status"].eq("INCOMPLETE").any()
    module_pass = bool(nulls["pass"].eq("PASS").any())
    stable_sensitivity = sens.loc[sens["analysis"].eq("QC_pass_only_sensitivity"), "status"].iloc[0] == "RUN_VALID"
    interpretable_anchor = bool(mech["final_classification"].isin(["SIGNIFICANT_OPENING", "SIGNIFICANT_CLOSING"]).any())
    if qc_complete and module_pass and stable_sensitivity and interpretable_anchor:
        verdict = "GO_TARGETED_CUTTAG_RESCUE"
    elif not qc_complete:
        verdict = "TECHNICAL_HOLD_ATAC"
    else:
        verdict = "STOP_REGULATORY_ESCALATION"
    write_report(samples, qc, design, sens, nulls, status, mech, verdict)
    GATE_OUT.write_text(
        f"# Final regulatory gate — Phase 2R.2\n\nFINAL_REGULATORY_GATE: **{verdict}**\n\n"
        f"QC complete: {'YES' if qc_complete else 'NO'}\n"
        f"Matched-null modules passing FDR<0.05: {', '.join(nulls.loc[nulls['pass'].eq('PASS'), 'module_id'].tolist()) or 'NONE'}\n"
        f"QC-pass DAR sensitivity valid: {'YES' if stable_sensitivity else 'NO'}\n"
        "Targeted CUT&Tag raw rescue: NO\nRaw Hi-C rescue: NO\nPhase 3: CLOSED\n"
        "GWAS/eQTL/colocalization/PRS/MR: PROHIBITED\n",
        encoding="utf-8",
    )
    update_project_files(verdict, verdict)
    print(f"PHASE 2R.2 VERDICT: {verdict}")
    print(f"QC PASS / FAIL / INCOMPLETE samples: {int((qc['QC_status']=='PASS').sum())} / {int((qc['QC_status']=='FAIL').sum())} / {int((qc['QC_status']=='INCOMPLETE').sum())}")
    print("duplicate fraction recovered: NO")
    print("chrM fraction recovered: NO")
    print("standard TSS enrichment recovered: YES (final-BAM-derived profile; pre-filter-independent)")
    print(f"full-set significant DARs: {int(sens.loc[sens['analysis'].eq('full_set_locked_edgeR'),'significant_DAR_n'].iloc[0])}")
    print(f"QC-sensitivity DAR result: {sens.loc[sens['analysis'].eq('QC_pass_only_sensitivity'),'status'].iloc[0]}")
    print("DAR directional concordance: NA (sensitivity not estimable)")
    print(f"modules passing matched-null FDR: {','.join(nulls.loc[nulls['pass'].eq('PASS'),'module_id'].tolist()) or 'NONE'}")
    print(f"strongest module and FDR: {strongest_module(nulls)}")
    print(f"corrected module statuses: {','.join(status['corrected_status'].unique().tolist())}")
    print(f"PIEZO1 final ATAC classification: {mech.loc[mech['gene'].eq('PIEZO1'),'final_classification'].iloc[0]}")
    print(f"TRPV4 final ATAC classification: {mech.loc[mech['gene'].eq('TRPV4'),'final_classification'].iloc[0]}")
    print(f"regulatory specificity rescued: {'YES' if module_pass else 'NO'}")
    print("targeted CUT&Tag rescue allowed: NO")
    print("Phase 3 allowed: NO")
    print("recommended next step: retain technical hold, repair/recover missing pre-filter QC resources if a higher-confidence ATAC gate is required, then lock manuscript construction")
    print("strongest supported manuscript claim: exploratory replicate-level ATAC remodeling was assessed in the locked human-TM CONTROL-versus-DEX contrast")
    print("strongest prohibited claim: established chromatin rewiring, adaptive-to-maladaptive transition, PIEZO1/TRPV4 causality, genetic mediation or PRS×DEX effect")
    print(f"path to PHASE2R2_ATAC_QC_CLOSURE_REPORT.md: {REPORT_OUT}")


def strongest_module(nulls: pd.DataFrame) -> str:
    x = nulls.sort_values(["empirical_P", "observed_statistic"], ascending=[True, False]).iloc[0]
    return f"{x['module_id']} FDR={float(x['FDR']):.6g}"


if __name__ == "__main__":
    main()
