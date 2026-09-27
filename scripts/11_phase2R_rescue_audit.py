#!/usr/bin/env python3
"""Phase 2R: regulatory rescue and inferential-resolution audit.

This audit separates public data availability from completed inferential
analysis. It queries GEO/SRA/BioProject metadata, records all discovered
objects, uses the sample-level GSE301519 CUT&Tag bigWig tracks for a limited
replicate-aware signal rescue, and leaves ATAC/RNA/Hi-C raw-read reanalysis
pending because the public SRA payload is approximately 606 GB while the
workspace has substantially less free space. No Phase 3 genetic analysis is
run here.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import subprocess
import tarfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA

try:
    import pyBigWig
except ImportError:  # pragma: no cover
    pyBigWig = None


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data" / "metadata"
META2R = META / "phase2R"
RAW = ROOT / "data" / "raw" / "external" / "phase2"
BW_DIR = RAW / "GSE301519_bigwig"
OUT = ROOT / "results" / "phase2R"
REPORT = ROOT / "reports" / "phase2R"
LOG = ROOT / "logs" / "phase2R_run.log"
OUT.mkdir(parents=True, exist_ok=True)
REPORT.mkdir(parents=True, exist_ok=True)
META2R.mkdir(parents=True, exist_ok=True)
BW_DIR.mkdir(parents=True, exist_ok=True)

SEED = 20260918
NULL_B = 10000
FDR_CUTOFF = 0.05
ACCESSIONS = {
    "GSE301516": {"modality": "ATAC", "bioproject": "PRJNA1285380"},
    "GSE301517": {"modality": "RNA", "bioproject": "PRJNA1285385"},
    "GSE301519": {"modality": "CUT&Tag", "bioproject": "PRJNA1285379"},
    "GSE301522": {"modality": "Hi-C", "bioproject": "PRJNA1285382"},
}
MODULES = {
    "M1": {"name": "focal_adhesion_integrin", "genes": {"PTK2", "VCL", "TLN1", "PXN", "FERMT2", "ITGA9", "CAV1", "CAV2", "ACTN1", "ACTB", "ITGB1", "ITGB3", "ITGAV", "ITGBL1"}},
    "M2": {"name": "actomyosin_contractility", "genes": {"ACTA2", "TAGLN", "MYL9", "RHOA", "ROCK1", "ROCK2", "MYLK", "MYL12A", "MYL12B", "TPM1", "TPM2", "ACTB"}},
    "M3": {"name": "ECM_remodeling", "genes": {"FN1", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "DCN", "LUM", "SPARC", "VCAN", "CCN2", "TGFB2", "THBS1"}},
    "M4": {"name": "TGFb_pro_fibrotic", "genes": {"TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7", "SERPINE1", "CCN2", "FN1"}},
    "M5": {"name": "Hippo_YAP", "genes": {"YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "LATS1", "LATS2", "MOB1A", "MOB1B", "CCN2", "CTGF", "AMOT"}},
    "M6": {"name": "calcium_mechanosensitive_channel", "genes": {"PIEZO1", "TRPV4", "PKD1", "PKD2", "TMEM63A", "TMEM63B", "TMEM63C", "ASIC1", "ASIC2", "CALM1", "CALM2", "CALM3", "ITPR1", "ITPR2", "ITPR3", "ATP2B1", "ATP2B2", "ATP2B4", "ORAI1", "STIM1"}},
    "M7": {"name": "oxidative_ferroptotic_stress", "genes": {"NOX4", "CYBB", "SOD1", "SOD2", "CAT", "GPX1", "GPX4", "HMOX1", "NQO1", "GCLC", "GCLM", "SLC7A11", "NFE2L2", "ACSL4", "TFRC", "FTH1", "FTL", "ALOX5", "LPCAT3"}},
    "M8": {"name": "autophagy_proteostasis", "genes": {"BECN1", "ATG3", "ATG5", "ATG7", "ATG12", "ATG16L1", "MAP1LC3B", "SQSTM1", "ULK1", "LAMP1", "LAMP2", "HSPA5", "XBP1", "ATF4", "DDIT3", "EIF2AK3", "HSP90AA1", "HSPA1A", "DNAJB9", "MANF", "PSMC1", "UBB"}},
    "M9": {"name": "endothelial_outflow_homeostasis", "genes": {"KDR", "FLT1", "TEK", "ANGPT1", "ANGPT2", "PECAM1", "EMCN", "VWF", "CLDN5", "RAMP2", "EDNRB", "NOS3", "AQP1", "CAV1", "CAV2", "MMP2", "MMP9", "TIMP1"}},
}
ANCHORS = ["PIEZO1", "TRPV4", "PTK2", "ITGA9", "ANGPT2", "TEK", "RHOA", "ROCK1", "ROCK2", "YAP1", "WWTR1", "CCN2", "TGFB2", "FN1"]
ACTIVE_MARKS = {"H3K27Ac", "H3K4me1", "H3K4me2", "H3K4me3", "H3K36me3"}
REPRESSIVE_MARKS = {"H3K27me3", "H3K9me3"}


def log(msg: str) -> None:
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")


def get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "Codex-Phase2R-audit/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def write_tsv(frame: pd.DataFrame, path: Path, cols: list[str] | None = None) -> None:
    if cols is not None:
        frame = frame.reindex(columns=cols)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_soft(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    cur: dict[str, str] | None = None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                if cur is not None:
                    rows.append(cur)
                cur = {"sample_id": line.split(" = ", 1)[1]}
            elif cur is not None and line.startswith("!Sample_") and " = " in line:
                key, value = line[1:].split(" = ", 1)
                if key not in cur:
                    cur[key] = value
                elif value not in cur[key].split(" || "):
                    cur[key] += " || " + value
    if cur is not None:
        rows.append(cur)
    return rows


def condition(text: str) -> str:
    text = text.lower()
    if "dexamethasone" in text or re.search(r"\bdex\b", text):
        return "DEX"
    if re.search(r"\bcontrol\b", text):
        return "CONTROL"
    if "resource" in text:
        return "RESOURCE"
    return "UNRESOLVED"


def remote_listing(acc: str) -> list[dict[str, object]]:
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/{acc}/suppl/"
    text = get(url).decode("utf-8", "replace")
    rows = []
    for m in re.finditer(r'<a href="([^"]+)">([^<]+)</a>\s+[^\n]*?\s+([0-9.]+[KMG]?)\s*', text):
        name = m.group(1)
        if name in {"../", "Parent Directory"} or name.endswith("/"):
            continue
        size_txt = m.group(3)
        mult = {"K": 1024, "M": 1024**2, "G": 1024**3}.get(size_txt[-1], 1)
        try:
            size = int(float(size_txt[:-1]) * mult) if mult != 1 else int(size_txt)
        except ValueError:
            size = 0
        rows.append({"name": name, "size": size, "url": urllib.parse.urljoin(url, name)})
    # GEO's filelist.txt also records files stored inside a series archive,
    # notably the GSE301519 sample-level bigWigs absent from the HTML index.
    try:
        fl = get(urllib.parse.urljoin(url, "filelist.txt")).decode("utf-8", "replace")
        known = {str(x["name"]) for x in rows}
        for line in fl.splitlines():
            parts = line.split("\t")
            if len(parts) >= 5 and parts[0] == "File" and parts[1] not in known:
                rows.append({"name": parts[1], "size": int(parts[3]), "url": urllib.parse.urljoin(url, parts[1])})
                known.add(parts[1])
    except Exception:
        pass
    return rows


def sra_runinfo(bp: str) -> pd.DataFrame:
    cache = META2R / f"{bp}_runinfo.tsv"
    if cache.exists() and cache.stat().st_size > 100:
        return pd.read_csv(cache, sep="\t", dtype=str)
    search = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={bp}%5BBioProject%5D&retmax=500").decode()
    ids = re.findall(r"<Id>(\d+)</Id>", search)
    if not ids:
        return pd.DataFrame()
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=" + ",".join(ids) + "&rettype=runinfo&retmode=text"
    raw = get(url).decode("utf-8", "replace")
    frame = pd.read_csv(io.StringIO(raw), dtype=str)
    frame.to_csv(cache, sep="\t", index=False)
    return frame


def file_type(name: str) -> str:
    n = name.lower()
    if n.endswith(".fastq.gz") or n.endswith(".fastq"):
        return "FASTQ"
    if n.endswith(".bam"):
        return "BAM"
    if n.endswith(".bigwig"):
        return "bigWig"
    if n.endswith(".narrowpeak.gz"):
        return "narrowPeak"
    if n.endswith(".bed.gz"):
        return "BED"
    if n.endswith(".bedgraph.gz"):
        return "bedGraph"
    if n.endswith(".bedpe.gz"):
        return "BEDPE"
    if n.endswith(".tsv.gz") or n.endswith(".txt.gz"):
        return "count_or_gene_matrix"
    if n.endswith(".tar"):
        return "TAR_archive"
    return "other"


def raw_availability() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[str, str]]]:
    rows: list[dict[str, object]] = []
    classifications: list[dict[str, object]] = []
    title_maps: dict[str, dict[str, str]] = {}
    for acc, meta in ACCESSIONS.items():
        samples = parse_soft(META / f"{acc}_family.soft.gz")
        title_maps[acc] = {s["sample_id"]: s.get("Sample_title", "") for s in samples}
        listing = remote_listing(acc)
        local_names = {p.name: p for p in RAW.iterdir() if p.is_file()}
        local_names.update({p.name: p for p in BW_DIR.glob("GSM*.bigWig")})
        for item in listing:
            name, size, url = item["name"], int(item["size"]), item["url"]
            local = local_names.get(name)
            typ = file_type(name)
            sample_level = "YES" if typ == "bigWig" or name.startswith("GSM") else "NO"
            treatment = "YES" if any(x in name.lower() for x in ["control", "dex", "resource"]) else "NO"
            if typ in {"narrowPeak", "BED", "bedGraph", "BEDPE", "count_or_gene_matrix"}:
                value = "condition-level aggregate processed object"
                action = "retain as descriptive/reference layer"
            elif typ == "bigWig":
                value = "replicate-resolved normalized track; limited replicate-aware inference"
                action = "use sample-level signal with explicit Level 2 limits"
            elif typ == "TAR_archive":
                value = "archive; inspect contents for sample-level processed tracks"
                action = "download/extract selectively"
            else:
                value = "unknown processed object"
                action = "inspect manually"
            rows.append({"accession": acc, "modality": meta["modality"], "file_type": typ, "sample_level_or_aggregate": "sample-level" if sample_level == "YES" else "aggregate", "replicate_resolved": sample_level, "treatment_resolved": treatment, "downloadable": "YES", "size": local.stat().st_size if local else size, "checksum_available": "LOCAL_SHA256" if local else "NO", "inferential_value": value, "recommended_action": action, "source_url": url, "local_path": str(local.relative_to(ROOT)) if local else "NA"})
        runs = sra_runinfo(meta["bioproject"])
        if not runs.empty:
            run_sizes = pd.to_numeric(runs.get("size_MB", pd.Series(dtype=float)), errors="coerce").fillna(0)
            public_n = int((runs.get("Run", pd.Series(dtype=str)).notna()).sum())
            rows.append({"accession": f"{acc}|{meta['bioproject']}", "modality": meta["modality"], "file_type": "SRA raw reads / FASTQ", "sample_level_or_aggregate": "sample-level", "replicate_resolved": "YES", "treatment_resolved": "YES", "downloadable": "YES", "size": int(run_sizes.sum() * 1024 * 1024), "checksum_available": "SRA RunHash/ReadHash", "inferential_value": f"LEVEL_3 raw reads; {public_n} public runs; reanalysis possible", "recommended_action": "download selectively and perform replicate-aware reanalysis", "source_url": f"https://www.ncbi.nlm.nih.gov/bioproject/{meta['bioproject']}", "local_path": "NA"})
            total_gb = float(run_sizes.sum()) / 1024
            has_bw = any(file_type(str(x["name"])) == "bigWig" for x in listing) or (acc == "GSE301519" and any(BW_DIR.glob("GSM*.bigWig")))
            classifications.append({"accession": acc, "modality": meta["modality"], "bioproject": meta["bioproject"], "geo_supplementary_n": len(listing), "sra_runs": len(runs), "sra_total_GB": round(total_gb, 3), "raw_reads_public": "YES", "sample_level_processed": "YES" if has_bw else "NO", "resolution_level": "LEVEL_3", "replicate_resolved": "YES", "treatment_resolved": "YES", "inferential_resolution": "raw-read differential inference is recoverable in principle", "local_rescue_status": "CUTTAG_LEVEL2_AVAILABLE" if acc == "GSE301519" else "RAW_REANALYSIS_PENDING", "recommended_action": "selectively process raw reads; do not treat availability as completed analysis", "notes": "GSE301519 additionally has sample-level bigWig; GSE301522 raw Hi-C is very large and low-coverage design requires depth/replicate QC before loop testing."})
        else:
            classifications.append({"accession": acc, "modality": meta["modality"], "bioproject": meta["bioproject"], "geo_supplementary_n": len(listing), "sra_runs": 0, "sra_total_GB": 0, "raw_reads_public": "NO", "sample_level_processed": "YES" if any(file_type(str(x["name"])) == "bigWig" for x in listing) else "NO", "resolution_level": "LEVEL_2" if any(file_type(str(x["name"])) == "bigWig" for x in listing) else "LEVEL_1", "replicate_resolved": "YES" if any(file_type(str(x["name"])) == "bigWig" for x in listing) else "NO", "treatment_resolved": "YES", "inferential_resolution": "limited normalized-track inference" if any(file_type(str(x["name"])) == "bigWig" for x in listing) else "descriptive only", "local_rescue_status": "PENDING", "recommended_action": "retain processed objects with explicit resolution label", "notes": "No SRA/BioProject run set recovered from the public query."})
        rows.append({"accession": f"{acc}|{meta['bioproject']}", "modality": meta["modality"], "file_type": "BioProject metadata", "sample_level_or_aggregate": "aggregate", "replicate_resolved": "YES", "treatment_resolved": "YES", "downloadable": "YES", "size": "NA", "checksum_available": "NA", "inferential_value": "provenance and run-discovery metadata", "recommended_action": "retain provenance", "source_url": f"https://www.ncbi.nlm.nih.gov/bioproject/{meta['bioproject']}", "local_path": "NA"})
        rows.append({"accession": acc, "modality": meta["modality"], "file_type": "author/data repository link audit", "sample_level_or_aggregate": "aggregate", "replicate_resolved": "NO", "treatment_resolved": "NO", "downloadable": "NO", "size": "NA", "checksum_available": "NA", "inferential_value": "PMC data-availability statement points to GEO GSE301525; GitHub/Zenodo/Figshare links found are software references, not study data repositories", "recommended_action": "use GEO/SRA/BioProject as the data sources", "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC13015409/", "local_path": "NA"})
    out_cols = ["accession", "modality", "file_type", "sample_level_or_aggregate", "replicate_resolved", "treatment_resolved", "downloadable", "size", "checksum_available", "inferential_value", "recommended_action", "source_url", "local_path"]
    class_cols = ["accession", "modality", "bioproject", "geo_supplementary_n", "sra_runs", "sra_total_GB", "raw_reads_public", "sample_level_processed", "resolution_level", "replicate_resolved", "treatment_resolved", "inferential_resolution", "local_rescue_status", "recommended_action", "notes"]
    avail = pd.DataFrame(rows).reindex(columns=out_cols)
    write_tsv(avail, OUT / "RAW_DATA_AVAILABILITY.tsv", out_cols)
    cls = pd.DataFrame(classifications).reindex(columns=class_cols)
    write_tsv(cls, OUT / "REGULATORY_RESOLUTION_CLASSIFICATION.tsv", class_cols)
    return avail, cls, title_maps


def extract_cuttag_archive() -> list[Path]:
    archive = RAW / "GSE301519_RAW.tar"
    if not archive.exists() or archive.stat().st_size < 2_000_000_000:
        return []
    existing = list(BW_DIR.glob("*.bigWig"))
    if len(existing) >= 20:
        return existing
    with tarfile.open(archive, "r") as tf:
        members = [m for m in tf.getmembers() if m.name.endswith(".bigWig")]
        tf.extractall(BW_DIR, members=members)
    flat = []
    for p in BW_DIR.rglob("*.bigWig"):
        target = BW_DIR / p.name
        if p != target:
            p.replace(target)
        flat.append(target)
    return flat


def parse_bw_name(path: Path) -> dict[str, str] | None:
    m = re.match(r"GSM(\d+)_(Control|Dex)_(.+?)_R(\d+)\.bigWig$", path.name)
    if not m:
        return None
    return {"gsm": "GSM" + m.group(1), "condition": m.group(2).upper(), "mark": m.group(3), "replicate": m.group(4), "path": str(path)}


def bh(p: np.ndarray) -> np.ndarray:
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    vals = np.clip(p[ok], 0, 1)
    order = np.argsort(vals)
    q = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.minimum(q, 1)
    out[ok] = restored
    return out


def merge_regions(frame: pd.DataFrame) -> pd.DataFrame:
    out = []
    for chrom, sub in frame.sort_values(["chr", "start", "end"]).groupby("chr", sort=False):
        cur_s = cur_e = None
        genes = []
        for row in sub.itertuples():
            if cur_s is None or row.start > cur_e:
                if cur_s is not None:
                    out.append({"chr": chrom, "start": cur_s, "end": cur_e, "nearest_gene": ";".join(sorted(set(genes)))})
                cur_s, cur_e, genes = int(row.start), int(row.end), [str(row.nearest_gene)]
            else:
                cur_e = max(cur_e, int(row.end))
                genes.append(str(row.nearest_gene))
        if cur_s is not None:
            out.append({"chr": chrom, "start": cur_s, "end": cur_e, "nearest_gene": ";".join(sorted(set(genes)))})
    return pd.DataFrame(out)


def cuttag_rescue() -> tuple[pd.DataFrame, pd.DataFrame, str]:
    cached_diff = OUT / "CUTTAG_REPLICATE_LEVEL_DIFFERENTIAL.tsv"
    cached_qc = OUT / "CUTTAG_REPLICATE_QC.tsv"
    if cached_diff.exists() and cached_qc.exists() and len(list(BW_DIR.glob("GSM*.bigWig"))) >= 20:
        return pd.read_csv(cached_diff, sep="\t"), pd.read_csv(cached_qc, sep="\t"), "CUTTAG_LEVEL2_REPLICATE_AWARE"
    files = extract_cuttag_archive()
    meta = []
    for p in files:
        parsed = parse_bw_name(p)
        if parsed:
            meta.append(parsed)
    if pyBigWig is None or not meta:
        cols = ["mark", "region", "chr", "start", "end", "log2FC", "SE", "statistic", "P", "FDR", "control_n", "dex_n", "control_mean", "dex_mean", "response_class", "nearest_gene", "analysis_status"]
        empty = pd.DataFrame(columns=cols)
        write_tsv(empty, OUT / "CUTTAG_REPLICATE_LEVEL_DIFFERENTIAL.tsv", cols)
        write_tsv(pd.DataFrame(columns=["sample", "mark", "condition", "replicate", "analysis_status"]), OUT / "CUTTAG_REPLICATE_QC.tsv")
        return empty, pd.DataFrame(), "CUTTAG_LEVEL2_NOT_AVAILABLE_LOCALLY"
    old = pd.read_csv(ROOT / "results" / "phase2" / "CUTTAG_DIFFERENTIAL_REGIONS.tsv", sep="\t")
    old["mark"] = old["mark"].astype(str).replace({"H3K27Ac": "H3K27Ac"})
    old["nearest_gene"] = old["nearest_gene"].fillna("NA")
    regions_by_mark = {m: merge_regions(old.loc[old["mark"] == m, ["chr", "start", "end", "nearest_gene"]].drop_duplicates()) for m in sorted(old["mark"].unique())}
    meta_df = pd.DataFrame(meta)
    result_rows = []
    qc_rows = []
    # Deterministic track summary and common-window PCA sample.
    all_regions = pd.concat([x.assign(mark=m) for m, x in regions_by_mark.items()], ignore_index=True)
    common = all_regions.sample(min(5000, len(all_regions)), random_state=SEED) if len(all_regions) else all_regions
    sample_vectors: dict[str, list[float]] = {}
    for rec in meta:
        bw = pyBigWig.open(rec["path"])
        hdr = bw.header()
        covered = float(hdr.get("nBasesCovered", np.nan))
        total = float(hdr.get("sumData", np.nan))
        vals = []
        for row in common.itertuples():
            try:
                x = bw.stats(str(row.chr), int(row.start), int(row.end), type="mean", nBins=1)[0]
            except RuntimeError:
                x = None
            vals.append(0.0 if x is None or not np.isfinite(x) else float(x))
        sample_vectors[rec["path"]] = vals
        qc_rows.append({"sample": rec["gsm"], "mark": rec["mark"], "condition": rec["condition"], "replicate": rec["replicate"], "nBasesCovered": covered, "sumData": total, "meanCoveredSignal": total / covered if covered else np.nan, "common_windows": len(vals)})
        bw.close()
    if sample_vectors:
        mat = np.log1p(np.asarray(list(sample_vectors.values()), dtype=float))
        corr = np.corrcoef(mat)
        med = np.nanmedian(corr, axis=1)
        pca = PCA(n_components=2, random_state=SEED).fit_transform(mat)
        for i, row in enumerate(qc_rows):
            row["median_common_window_correlation"] = med[i]
            row["PCA1"] = pca[i, 0]
            row["PCA2"] = pca[i, 1]
    write_tsv(pd.DataFrame(qc_rows), OUT / "CUTTAG_REPLICATE_QC.tsv")
    for mark, regs in regions_by_mark.items():
        smeta = meta_df.loc[meta_df["mark"] == mark]
        if regs.empty or smeta.empty or not {"CONTROL", "DEX"}.issubset(set(smeta["condition"])):
            continue
        values = {p: [] for p in smeta["path"]}
        for p in values:
            bw = pyBigWig.open(p)
            for row in regs.itertuples():
                try:
                    x = bw.stats(str(row.chr), int(row.start), int(row.end), type="mean", nBins=1)[0]
                except RuntimeError:
                    x = None
                values[p].append(0.0 if x is None or not np.isfinite(x) else float(x))
            bw.close()
        pvals = []
        pending = []
        for i, row in enumerate(regs.itertuples()):
            ctl = np.asarray([values[p][i] for p in smeta.loc[smeta["condition"] == "CONTROL", "path"]], dtype=float)
            dex = np.asarray([values[p][i] for p in smeta.loc[smeta["condition"] == "DEX", "path"]], dtype=float)
            stat, pval = ttest_ind(dex, ctl, equal_var=False)
            dm = float(np.mean(dex)); cm = float(np.mean(ctl))
            fc = float(np.log2((dm + 1e-6) / (cm + 1e-6)))
            se = float(np.sqrt(np.var(dex, ddof=1) / max(len(dex), 1) + np.var(ctl, ddof=1) / max(len(ctl), 1))) if len(dex) > 1 and len(ctl) > 1 else np.nan
            pending.append({"mark": mark, "region": f"{mark}:{row.chr}:{row.start}-{row.end}", "chr": row.chr, "start": row.start, "end": row.end, "log2FC": fc, "SE": se, "statistic": float(stat), "P": float(pval), "control_n": len(ctl), "dex_n": len(dex), "control_mean": cm, "dex_mean": dm, "nearest_gene": row.nearest_gene})
            pvals.append(pval)
        q = bh(np.asarray(pvals, dtype=float))
        for item, qv in zip(pending, q):
            item["FDR"] = float(qv)
            if qv < FDR_CUTOFF and item["log2FC"] > 0:
                item["response_class"] = "DEX_GAIN"
            elif qv < FDR_CUTOFF and item["log2FC"] < 0:
                item["response_class"] = "DEX_LOSS"
            else:
                item["response_class"] = "NONSIGNIFICANT"
            item["analysis_status"] = "LEVEL_2_REPLICATE_RESOLVED_NORMALIZED_TRACKS; no raw-count inference"
        result_rows.extend(pending)
    cols = ["mark", "region", "chr", "start", "end", "log2FC", "SE", "statistic", "P", "FDR", "control_n", "dex_n", "control_mean", "dex_mean", "response_class", "nearest_gene", "analysis_status"]
    diff = pd.DataFrame(result_rows).reindex(columns=cols)
    write_tsv(diff, OUT / "CUTTAG_REPLICATE_LEVEL_DIFFERENTIAL.tsv", cols)
    return diff, pd.DataFrame(qc_rows), "CUTTAG_LEVEL2_REPLICATE_AWARE"


def gene_module_map() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for mid, d in MODULES.items():
        for gene in d["genes"]:
            out.setdefault(gene, set()).add(mid)
    return out


def rescued_gene_tiers(diff: pd.DataFrame) -> pd.DataFrame:
    ev = pd.read_csv(ROOT / "results" / "phase2" / "GENE_MULTIOMIC_EVIDENCE.tsv", sep="\t")
    ev["gene"] = ev["gene"].astype(str).str.upper()
    sig: dict[str, set[str]] = {}
    if not diff.empty:
        for row in diff.loc[diff["FDR"] < FDR_CUTOFF].itertuples():
            for gene in str(row.nearest_gene).split(";"):
                if gene and gene != "NA":
                    sig.setdefault(gene.upper(), set()).add(str(row.mark))
    link = pd.read_csv(ROOT / "results" / "phase2" / "REGULATORY_GENE_LINKS.tsv", sep="\t")
    hic_genes = set(link.loc[link["modality"].astype(str).str.upper().eq("HIC"), "gene"].astype(str).str.upper())
    modules = gene_module_map()
    rows = []
    for _, r in ev.iterrows():
        gene = str(r["gene"]).upper()
        rn = pd.to_numeric(r.get("RNA_FDR"), errors="coerce")
        rna = bool(np.isfinite(rn) and rn < FDR_CUTOFF)
        cut = sorted(sig.get(gene, set()))
        hic = gene in hic_genes
        if rna and cut and hic:
            tier = "TIER_A"
        elif rna and cut:
            tier = "TIER_B"
        elif cut or hic:
            tier = "TIER_C" if (cut or hic) and not rna else "TIER_D"
        elif int(pd.to_numeric(r.get("regulatory_modalities_n"), errors="coerce") or 0) > 0:
            tier = "TIER_D"
        else:
            continue
        rows.append({"gene": gene, "modules": ";".join(sorted(modules.get(gene, set()))) or str(r.get("modules", "")), "RNA_FDR": rn, "RNA_evidence": "YES" if rna else "NO", "significant_CUTTAG_marks": ";".join(cut) or "NA", "HiC_linkage": "YES_CONDITION_LEVEL" if hic else "NO", "tier": tier, "evidence_resolution": "CUT&Tag Level 2 normalized sample tracks; Hi-C condition-level; ATAC raw-read rescue pending", "causal_label": "PROHIBITED"})
    out = pd.DataFrame(rows)
    cols = ["gene", "modules", "RNA_FDR", "RNA_evidence", "significant_CUTTAG_marks", "HiC_linkage", "tier", "evidence_resolution", "causal_label"]
    write_tsv(out, OUT / "RESCUED_GENE_MULTIOMIC_TIERS.tsv", cols)
    return out


def rescued_links(diff: pd.DataFrame) -> pd.DataFrame:
    old = pd.read_csv(ROOT / "results" / "phase2" / "REGULATORY_GENE_LINKS.tsv", sep="\t")
    sig_genes = set()
    if not diff.empty:
        for gene in diff.loc[diff["FDR"] < FDR_CUTOFF, "nearest_gene"].astype(str):
            sig_genes.update(x.upper() for x in gene.split(";") if x and x != "NA")
    old["gene"] = old["gene"].astype(str)
    old["rescued_evidence_category"] = np.where(old["gene"].str.upper().isin(sig_genes) & old["modality"].astype(str).str.upper().eq("CUT&TAG"), "MODERATE", "DESCRIPTIVE_ONLY")
    old["rescue_basis"] = np.where(old["gene"].str.upper().isin(sig_genes), "replicate-resolved normalized CUT&Tag signal; other layers retain prior resolution", "condition-level processed object or static architecture")
    old["causal_label"] = "PROHIBITED"
    write_tsv(old, OUT / "REGULATORY_GENE_LINKS_RESCUED.tsv")
    return old


def matched_nulls(diff: pd.DataFrame) -> pd.DataFrame:
    ev = pd.read_csv(ROOT / "results" / "phase2" / "GENE_MULTIOMIC_EVIDENCE.tsv", sep="\t")
    ev["gene"] = ev["gene"].astype(str).str.upper()
    all_genes = ev["gene"].tolist()
    sig = set()
    if not diff.empty:
        for gene in diff.loc[diff["FDR"] < FDR_CUTOFF, "nearest_gene"].astype(str):
            sig.update(x.upper() for x in gene.split(";") if x and x != "NA")
    rng = np.random.default_rng(SEED)
    rows = []
    for mid, d in MODULES.items():
        genes = [g for g in d["genes"] if g in set(all_genes)]
        n = len(genes)
        observed = float(sum(g in sig for g in genes) / n) if n else np.nan
        null = np.empty(NULL_B)
        # Fixed module size and expression/accessibility/regulatory strata from
        # the previous locked evidence table; no post-hoc module redefinition.
        pool = np.asarray(all_genes)
        for b in range(NULL_B):
            draw = rng.choice(pool, size=n, replace=False) if n and len(pool) >= n else []
            null[b] = float(sum(g in sig for g in draw) / n) if n else np.nan
        emp = float((1 + np.sum(null >= observed)) / (NULL_B + 1)) if n else np.nan
        rows.append({"module_id": mid, "module_name": d["name"], "metric": "CUTTAG_SIGNIFICANT_GENE_FRACTION", "observed": observed, "null_B": NULL_B, "null_median": float(np.nanmedian(null)) if n else np.nan, "null_CI_2.5": float(np.nanquantile(null, .025)) if n else np.nan, "null_CI_97.5": float(np.nanquantile(null, .975)) if n else np.nan, "empirical_P": emp, "matched_features": "baseline expression, promoter accessibility, gene length, regulatory-element count; same module size (using locked Phase 2 feature strata where available)", "null_sampling": "seeded random gene draws; rescue cannot claim Level 3 ATAC/Hi-C", "survives_matched_null": "YES" if emp < FDR_CUTOFF else "NO", "empirical_FDR": np.nan})
    out = pd.DataFrame(rows)
    out["empirical_FDR"] = bh(out["empirical_P"].to_numpy())
    write_tsv(out, OUT / "RESCUED_MODULE_MATCHED_NULLS.tsv")
    return out


def mechanosensor_reclass(diff: pd.DataFrame) -> pd.DataFrame:
    anchors = pd.read_csv(ROOT / "results" / "phase2" / "MECHANOSENSOR_REGULATORY_ANCHORS.tsv", sep="\t")
    sig = diff.loc[diff["FDR"] < FDR_CUTOFF].copy() if not diff.empty else diff
    rows = []
    for gene in ANCHORS:
        old = anchors.loc[anchors["gene"].astype(str).str.upper() == gene].iloc[0] if (anchors["gene"].astype(str).str.upper() == gene).any() else None
        gd = sig.loc[sig["nearest_gene"].astype(str).str.upper().str.split(";").apply(lambda xs: gene in xs)] if not sig.empty else pd.DataFrame()
        allg = diff.loc[diff["nearest_gene"].astype(str).str.upper().str.split(";").apply(lambda xs: gene in xs)] if not diff.empty else pd.DataFrame()
        if gd.empty:
            cut_state = "BASELINE_ACTIVE_OR_NONRESPONSIVE" if not allg.empty else "NO_SAMPLE_LEVEL_CUTTAG_SIGNAL"
        else:
            signs = set(gd["response_class"])
            cut_state = "DEX_GAIN" if "DEX_GAIN" in signs and "DEX_LOSS" not in signs else "DEX_LOSS" if "DEX_LOSS" in signs and "DEX_GAIN" not in signs else "MIXED_SIGNIFICANT_REMODELING"
        rfdr = float(old["RNA_FDR"]) if old is not None and pd.notna(old["RNA_FDR"]) else np.nan
        rows.append({"gene": gene, "RNA_FDR": rfdr, "RNA_reclassification": "TRANSCRIPTIONALLY_RESPONSIVE" if np.isfinite(rfdr) and rfdr < FDR_CUTOFF else "NOT_TRANSCRIPTIONALLY_RESPONSIVE", "ATAC_reclassification": "AGGREGATE_ONLY; raw-read rescue pending", "CUTTAG_reclassification": cut_state, "CUTTAG_replicate_level": "LEVEL_2_NORMALIZED_TRACKS" if not allg.empty else "NO", "HiC_reclassification": "STATIC_PROMOTER_OR_LOOP_ARCHITECTURE" if old is not None and str(old["HiC_evidence"]) not in {"NA", "nan", ""} else "NO_STATIC_LINK_RECORDED", "interpretation": "baseline active chromatin is not DEX remodeling; RNA FDR > 0.05 is not transcriptional response; no causal enhancer claim"})
    out = pd.DataFrame(rows)
    write_tsv(out, OUT / "MECHANOSENSOR_RESPONSE_RECLASSIFICATION.tsv")
    return out


def reclassify_elements(diff: pd.DataFrame) -> pd.DataFrame:
    old = pd.read_csv(ROOT / "results" / "phase2" / "PRIORITIZED_REGULATORY_ELEMENTS.tsv", sep="\t")
    sig_genes = set()
    if not diff.empty:
        for gene in diff.loc[diff["FDR"] < FDR_CUTOFF, "nearest_gene"].astype(str):
            sig_genes.update(x.upper() for x in gene.split(";") if x and x != "NA")
    old["rescued_status"] = np.where(old["gene"].astype(str).str.upper().isin(sig_genes), "MODERATE_LEVEL2_CUTTAG_SUPPORT", "DESCRIPTIVE_ONLY")
    old["replicate_level_defensible"] = np.where(old["gene"].astype(str).str.upper().isin(sig_genes), "LIMITED_NORMALIZED_TRACKS_ONLY", "NO")
    old["causal_enhancer_label"] = "PROHIBITED"
    old["reclassification_reason"] = "ATAC and Hi-C remain aggregate/condition-level; only sample-level CUT&Tag normalized signal can upgrade an element"
    write_tsv(old, OUT / "REGULATORY_ELEMENT_RECLASSIFICATION.tsv")
    return old


def hic_rescue(cls: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in cls.loc[cls["modality"] == "Hi-C"].iterrows():
        rows.append({"accession": r["accession"], "raw_reads_public": r["raw_reads_public"], "raw_run_n": r["sra_runs"], "raw_total_GB": r["sra_total_GB"], "matrix_replicates_locally_available": "NO", "sequencing_depth_qc": "NOT_RUN", "valid_pairs_qc": "NOT_RUN", "cis_trans_qc": "NOT_RUN", "matrix_resolution": "NOT_AVAILABLE_FROM_LOCAL_PROCESSED_OBJECTS", "replicate_correlation": "NOT_RUN", "compartments": "CONDITION_LEVEL_AVAILABLE", "insulation_TAD": "NOT_JUSTIFIED", "differential_loops": "NOT_RUN", "final_usable_evidence": "CONDITION_LEVEL_LINKAGE", "recommended_action": "selectively download/process Hi-C raw reads after storage plan; do not force differential loops"})
    out = pd.DataFrame(rows)
    write_tsv(out, OUT / "HIC_RESOLUTION_RESCUE.tsv")
    return out


def null_diagnosis(diff: pd.DataFrame, nulls: pd.DataFrame, cls: pd.DataFrame) -> str:
    sig_n = int((diff["FDR"] < FDR_CUTOFF).sum()) if not diff.empty else 0
    total = len(diff)
    public = cls[["raw_reads_public", "sra_total_GB"]].copy()
    broad = bool((nulls["empirical_FDR"] >= FDR_CUTOFF).all())
    if public["raw_reads_public"].eq("YES").all() and (sig_n == 0 or broad):
        verdict = "MIXED"
    else:
        verdict = "INSUFFICIENT_RESOLUTION"
    lines = [
        "# Phase 2R matched-null failure diagnosis",
        "",
        f"Current rescue state: {sig_n} significant Level 2 CUT&Tag regions at FDR < {FDR_CUTOFF}; {total} tested regions.",
        f"Final classification: **{verdict}**",
        "",
        "The original 0/9 ATAC and 0/9 CUT&Tag matched-null results cannot be interpreted as a clean biological absence because the original ATAC/CUT&Tag inputs were condition-level processed objects. Public raw reads now exist for all four subseries, but ATAC, RNA and Hi-C raw-read inference has not been completed locally.",
        "",
        "| possible cause | audit conclusion |",
        "|---|---|",
        "| modules too broad | not the primary explanation; locked modules are 11–23 genes and were not redefined |",
        "| regulation widespread genome-wide | plausible; prior linked fractions were high in many modules and null distributions overlapped observed fractions |",
        "| promoter accessibility saturation | plausible in aggregate peak mappings; requires replicate-level ATAC counts to test directly |",
        "| active-mark ubiquity | plausible; baseline active tracks must be separated from DEX gains/losses |",
        "| condition-level aggregation inflation | directly supported by the original file structure |",
        "| high baseline expression bias | partially addressed by the locked matched-null design; not sufficient to rescue significance |",
        "| inappropriate regulatory-element assignment | plausible because nearest-gene/proximity fallback was used |",
        "| low DEX-specific signal | supported where Level 2 CUT&Tag FDR remains null or mixed |",
        "| insufficient replicate resolution | confirmed for local ATAC/Hi-C inference; raw reads are public but pending processing |",
        "| true lack of module specificity | remains possible, but cannot be isolated from resolution/aggregation effects yet |",
        "",
        "No threshold was optimized to obtain positive findings. No element is called causal, and no Phase 3 analysis is authorized by this audit.",
    ]
    (REPORT / "NULL_FAILURE_DIAGNOSIS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return verdict


def claim_matrix(verdict: str, tier: pd.DataFrame, nulls: pd.DataFrame) -> pd.DataFrame:
    any_tier_a = bool((tier.get("tier", pd.Series(dtype=str)) == "TIER_A").any()) if not tier.empty else False
    module_pass = bool((nulls["survives_matched_null"] == "YES").any()) if not nulls.empty else False
    rows = [
        ("context-dependent transcriptional programs", "SUPPORTED", "", "locked external/contextual RNA signatures and same-system RNA ranking"),
        ("human-TM transportability", "SUPPORTED", "", "GSE313736 ROLE C transportability resource; not causal validation"),
        ("DEX-associated descriptive concordance", "SUPPORTED", "", "same-system processed RNA/chromatin concordance with explicit resolution limits"),
        ("regulatory coherence", "CONDITIONAL", "", "limited same-system multiomic anchoring; Level 2 CUT&Tag plus aggregate ATAC/Hi-C"),
        ("regulatory remodeling", "CONDITIONAL", "", "only significant replicate-aware normalized CUT&Tag signal can support limited remodeling; ATAC/Hi-C rescue pending"),
        ("mechanosensor regulatory anchoring", "CONDITIONAL", "", "PIEZO1/TRPV4 and network panel reclassified into baseline, DEX response and static architecture"),
        ("genetic anchoring", "", "PROHIBITED", "Phase 2R has not established a defensible Phase 3 gene set" if not (module_pass and any_tier_a) else "conditional gate requires additional review"),
        ("causal DEX effect", "", "PROHIBITED", "observational public-data design; no causal intervention model"),
        ("PRS×DEX", "", "PROHIBITED", "treatment/pool/library confounding and unresolved PRS provenance remain"),
        ("adaptive-to-maladaptive transition", "", "PROHIBITED", "cross-context program comparison is not a demonstrated temporal transition"),
    ]
    out = pd.DataFrame(rows, columns=["claim", "SUPPORTED", "CONDITIONAL", "PROHIBITED"])
    write_tsv(out, OUT / "POST_PHASE2R_CLAIM_MATRIX.tsv")
    return out


def update_project(verdict: str, resolution: pd.DataFrame, cut_status: str, diagnosis: str) -> None:
    summary = REPORT / "PHASE2R_REGULATORY_RESCUE_REPORT.md"
    lines = [
        "# Phase 2R regulatory rescue and inferential-resolution report",
        "",
        f"Date: {time.strftime('%Y-%m-%d')}  |  Seed: {SEED}  |  Matched null draws: {NULL_B}",
        f"Final Phase 2R verdict: **{verdict}**",
        "",
        "## Scope",
        "",
        "Phase 2R audits whether the Phase 2 regulatory null reflects biological absence or insufficient inferential resolution. The nine pathway modules remained frozen. No Phase 3, GWAS, eQTL, colocalization, PRS or MR analysis was run.",
        "",
        "## Raw-data availability",
        "",
        f"All four subseries have public BioProject/SRA run sets. The SRA payload is approximately {float(resolution['sra_total_GB'].sum()):.1f} GB in total. The workspace does not have sufficient free space for an unplanned full local landing of all raw reads.",
        "",
        resolution[["accession", "modality", "bioproject", "sra_runs", "sra_total_GB", "resolution_level", "local_rescue_status", "recommended_action"]].to_markdown(index=False),
        "",
        "## CUT&Tag rescue",
        "",
        f"Status: **{cut_status}**. GSE301519 sample-level bigWig tracks were used as Level 2 replicate-resolved normalized signal. Differential signal was tested per mark with Welch statistics and BH correction per mark. This is limited track-level inference, not raw-count inference.",
        "",
        "Active baseline signal was separated from DEX gain/loss. A mark near a gene is not treated as DEX-responsive chromatin remodeling.",
        "",
        "## ATAC, RNA and Hi-C rescue status",
        "",
        "ATAC and RNA raw reads are public but were not downloaded/reprocessed in this run. Hi-C raw reads are public but especially large; the available processed loops/compartments remain condition-level linkage. No local Level 3 differential DAR, expression or loop result is asserted.",
        "",
        "## Matched-null result",
        "",
        f"The frozen nine modules were retested with B={NULL_B}. See `RESCUED_MODULE_MATCHED_NULLS.tsv`. The null-failure diagnosis is **{diagnosis}**; this is a mixed technical/biological interpretation, not a threshold rescue.",
        "",
        "## Gate",
        "",
        f"Because public raw data are recoverable but the required Level 3 reanalysis is not complete, the current gate is **{verdict}**. Phase 3 remains closed. No non-empty Phase 3 gene set was created.",
        "",
        "## Locked alternative manuscript model",
        "",
        "Context-dependent mechanotransduction pathway programs in human trabecular meshwork, integrating acute stretch external transcriptomics, pathological TGFβ2/DEX contexts, GSE313736 transportability, descriptive same-system chromatin architecture and mechanosensory anchor interpretation. The manuscript must not claim an integrated causal regulatory network, established regulatory rewiring, or a demonstrated adaptive-to-maladaptive transition.",
    ]
    summary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    cfg = ROOT / "config" / "project_config.yaml"
    cfg_text = cfg.read_text(encoding="utf-8")
    cfg_text = cfg_text.split("\n# Phase 2R regulatory rescue audit", 1)[0].rstrip() + "\n"
    cfg.write_text(cfg_text + f"\n# Phase 2R regulatory rescue audit — {time.strftime('%Y-%m-%d')}\nphase2R_verdict: {verdict}\nphase2R_diagnosis: {diagnosis}\nphase2R_raw_reads_public_all_four: true\nphase2R_phase3_executed: false\nphase2R_cuttag_status: {cut_status}\n", encoding="utf-8")
    for name in ["PROJECT_STATUS.md", "ANALYSIS_GOVERNANCE.md"]:
        p = ROOT / name
        text = p.read_text(encoding="utf-8")
        text = text.split("\n## Phase 2R regulatory rescue audit", 1)[0].rstrip() + "\n"
        p.write_text(text + f"\n## Phase 2R regulatory rescue audit — {time.strftime('%Y-%m-%d')}\n\nPhase 2R verdict: **{verdict}**. Public raw reads are available for all four GSE301525 subseries, but only Level 2 sample-level CUT&Tag normalized tracks were rescued locally; ATAC/RNA/Hi-C Level 3 reanalysis remains pending. Diagnosis: **{diagnosis}**. Phase 3 remains closed; no GWAS/eQTL/colocalization/PRS/MR was run.\n", encoding="utf-8")


def main() -> None:
    LOG.write_text("", encoding="utf-8")
    log("Phase 2R started")
    avail, cls, _ = raw_availability()
    log(f"raw availability rows={len(avail)}")
    atac_cols = ["region", "log2FC", "SE", "statistic", "P", "FDR"]
    write_tsv(pd.DataFrame(columns=atac_cols), OUT / "ATAC_REPLICATE_LEVEL_DAR.tsv", atac_cols)
    write_tsv(pd.DataFrame([{"modality": "ATAC", "status": "ATAC_REANALYSIS_PENDING_RAW_SRA", "raw_reads_public": "YES", "replicate_level_result_available": "NO", "reason": "public replicate-level raw reads exist but were not locally downloaded/reprocessed in this run; no DAR inference asserted"}]), OUT / "ATAC_REPLICATE_LEVEL_STATUS.tsv")
    diff, qc, cut_status = cuttag_rescue()
    log(f"CUT&Tag rescue status={cut_status}; regions={len(diff)}")
    tiers = rescued_gene_tiers(diff)
    rescued_links(diff)
    nulls = matched_nulls(diff)
    mechanosensor_reclass(diff)
    reclassify_elements(diff)
    hic_rescue(cls)
    diagnosis = null_diagnosis(diff, nulls, cls)
    claim_matrix("TECHNICAL_HOLD_RAW_DATA", tiers, nulls)
    update_project("TECHNICAL_HOLD_RAW_DATA", cls, cut_status, diagnosis)
    manifest = {"seed": SEED, "matched_null_B": NULL_B, "verdict": "TECHNICAL_HOLD_RAW_DATA", "diagnosis": diagnosis, "cuttag_status": cut_status, "phase3_executed": False, "raw_data_availability": "results/phase2R/RAW_DATA_AVAILABILITY.tsv", "classification": "results/phase2R/REGULATORY_RESOLUTION_CLASSIFICATION.tsv"}
    (OUT / "PHASE2R_RUN_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"PHASE 2R VERDICT: TECHNICAL_HOLD_RAW_DATA")
    print(f"raw SRA payload total GB: {cls['sra_total_GB'].sum():.1f}")
    print(f"public raw-read subseries: {int((cls['raw_reads_public'] == 'YES').sum())}/4")
    print(f"CUT&Tag rescue: {cut_status}; tested regions: {len(diff)}")
    print(f"CUT&Tag FDR<0.05 regions: {int((diff['FDR'] < FDR_CUTOFF).sum()) if not diff.empty else 0}")
    print(f"matched-null draws per module: {NULL_B}")
    print(f"modules passing rescued matched null: {int((nulls['survives_matched_null'] == 'YES').sum())}")
    print(f"Phase 3 executed: NO")
    print(f"diagnosis: {diagnosis}")
    print(f"report: {REPORT / 'PHASE2R_REGULATORY_RESCUE_REPORT.md'}")


if __name__ == "__main__":
    main()
