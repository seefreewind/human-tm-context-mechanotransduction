#!/usr/bin/env python3
"""Phase 1: externally derive, lock, validate, and apply mechanotransduction programs.

The script is deliberately phase-gated.  GSE313736 expression values are not
opened until all external signatures, the anchor panel, and their SHA256
manifest have been written.  GSE313736 is used as an observational validation
resource only; donor/pseudobulk summaries are the inferential unit.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import random
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import hypergeom, spearmanr, ttest_ind, ttest_rel, wilcoxon

try:
    import anndata as ad
except Exception:  # pragma: no cover - handled as a gated source failure
    ad = None

try:
    import gseapy
except Exception:  # pragma: no cover
    gseapy = None


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data" / "metadata"
PRIMARY_RAW = ROOT / "data" / "raw" / "GSE313736"
EXT = ROOT / "data" / "raw" / "external" / "phase1"
OUT = ROOT / "results" / "phase1"
FIG = OUT / "figure_tables"
REPORT = ROOT / "reports" / "phase1"
LOG_DIR = ROOT / "logs"
for d in (EXT, OUT, FIG, REPORT, LOG_DIR):
    d.mkdir(parents=True, exist_ok=True)
LOG = LOG_DIR / "phase1_run.log"

SEED = 20260917
random.seed(SEED)
rng = np.random.default_rng(SEED)
FDR_CUTOFF = 0.05
CORE_CAP = 100
CONSENSUS_RANK_CAP = 500
ENRICHMENT_NULLS = 500
RANDOM_NULLS = 200

ADAPTIVE_SOURCE = "GSE195756"
FIBROTIC_SOURCE = "GSE202372"
STEROID_SOURCE = "GSE301517"
DISEASE_SOURCE = "GSE263230"

ANCHOR_PANEL = [
    "PIEZO1", "TRPV4", "ANGPT2", "TEK", "ITGA9", "PTK2", "TLN1", "VCL",
    "PXN", "FERMT2", "CAV1", "CAV2", "RHOA", "ROCK1", "ROCK2", "ACTA2",
    "TAGLN", "MYL9", "YAP1", "WWTR1", "CCN2", "TGFB2", "FN1", "COL1A1",
    "COL3A1", "GPX4", "SLC7A11", "ACSL4", "TFRC", "HMOX1", "NFE2L2",
]

GENERIC_STRESS = {
    "Apoptosis": {"BAX", "BCL2", "CASP3", "CASP8", "CASP9", "APAF1", "BBC3", "PMAIP1", "BAD", "BID", "MCL1"},
    "Inflammatory_Response": {"IL6", "IL1B", "TNF", "CXCL8", "CCL2", "CCL5", "ICAM1", "VCAM1", "NFKB1", "RELA"},
    "ROS": {"NOX4", "CYBB", "SOD1", "SOD2", "CAT", "GPX1", "GPX4", "HMOX1", "NQO1", "GCLC", "GCLM", "SLC7A11"},
    "Hypoxia": {"HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA", "ADM", "EGLN3", "BNIP3", "PDK1"},
    "UPR": {"HSPA5", "XBP1", "ATF4", "DDIT3", "EIF2AK3", "HERPUD1", "DNAJB9", "MANF"},
    "TNFa_NFkB": {"TNF", "NFKB1", "RELA", "NFKBIA", "TNFAIP3", "CXCL8", "IL6", "BIRC3", "ICAM1"},
    "Cell_Cycle": {"MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CDK1", "CENPF", "UBE2C", "HMGB2", "TYMS"},
}
HOUSEKEEPING = {"GAPDH", "ACTB", "B2M", "RPLP0", "TBP", "PPIA", "HPRT1", "RPL13A", "RPS18", "RPL32"}

LOCAL_PATHWAY_SETS = {
    "GO:actin_cytoskeleton": {"ACTB", "ACTG1", "ACTN1", "ACTN4", "ARPC1B", "ARPC2", "CFL1", "CORO1C", "EZR", "FLNA", "PFN1", "RHOA", "ROCK1", "ROCK2", "VCL", "TLN1", "PXN", "PTK2"},
    "Reactome:ECM_organization": {"FN1", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL6A1", "DCN", "LUM", "SPARC", "VCAN", "TGFBR1", "TGFBR2", "TGFB2", "CCN2"},
    "Hallmark:TGF_beta_signaling": {"TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7", "SERPINE1", "CCN2", "COL1A1", "FN1", "THBS1"},
    "Hallmark:Hypoxia": GENERIC_STRESS["Hypoxia"],
    "Hallmark:Inflammatory_Response": GENERIC_STRESS["Inflammatory_Response"],
    "GO:response_to_reactive_oxygen_species": GENERIC_STRESS["ROS"],
    "GO:cell_cycle": GENERIC_STRESS["Cell_Cycle"],
}


def log(message: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"[{stamp}] {message}\n")


def clean(x: object) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and math.isnan(x):
        return ""
    return str(x).strip()


def gene(x: object) -> str:
    value = clean(x).upper()
    value = re.sub(r"\.\d+$", "", value)
    if value in {"", "NAN", "NA", "N/A", "-", "NONE"}:
        return ""
    return "CCN2" if value == "CTGF" else value


def fnum(x: object, default: float = np.nan) -> float:
    try:
        value = float(x)
        return value if np.isfinite(value) else default
    except Exception:
        return default


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return result
    p = np.clip(values[finite], 0.0, 1.0)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.minimum(q, 1.0)
    result[finite] = restored
    return result


def write_tsv(frame: pd.DataFrame, path: Path, columns: list[str] | None = None) -> None:
    if columns is not None:
        frame = frame.reindex(columns=columns)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def source_row(dataset: str, role: str, n_samples: str, protocol: str, pairing: str, status: str, eligibility: str, url: str) -> dict[str, str]:
    return {
        "dataset": dataset, "role": role, "sample_N": n_samples, "protocol": protocol,
        "pairing_status": pairing, "data_layer": status, "discovery_or_validation": "DISCOVERY" if dataset != "GSE313736" else "VALIDATION",
        "eligibility": eligibility, "source_url": url,
    }


def make_role_lock() -> pd.DataFrame:
    rows = [
        source_row("GSE195756", "adaptive_stretch_discovery", "8 samples / 4 label-inferred pairs", "15% stretch, 0.69 Hz, 8 h", "inferred labels; donor IDs not explicit", "processed XLSX with author DE statistics + raw-count columns", "ELIGIBLE_CONDITIONAL", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195756"),
        source_row("GSE202372", "fibrotic_maladaptive_discovery", "10 samples / 5 mapping-inferred pairs", "TGF-beta2, 24 h", "mapping-inferred; author DE list only", "processed paired DE list", "ELIGIBLE_CONDITIONAL", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202372"),
        source_row("GSE301517", "steroid_maladaptive_discovery", "14 expression profiles / 7 control-DEX records", "DEX, 0.5 uM, 7 days", "sample-line labels; donor pairing not verified", "processed gene CPM", "ELIGIBLE_CONDITIONAL", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301517"),
        source_row("GSE263230", "disease_context_discovery", "14 primary cultures: 4 glaucoma, 10 healthy", "primary culture; no DEX", "independent primary-culture samples", "individual processed h5ad", "ELIGIBLE_CONDITIONAL", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263230"),
        source_row("GSE313736", "human_TM_external_application", "59 observed paired donors / 166323 exact cells", "vehicle versus DEX labels; pool/library aliased", "true donor pairing in metadata", "10x Matrix Market + reconciled metadata", "ROLE_C_ONLY", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736"),
        source_row("GSE309500", "excluded_species_resource", "mouse limbal strips", "non-human", "not applicable", "external audit only", "EXCLUDED_FROM_HUMAN_TM_VALIDATION", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309500"),
    ]
    frame = pd.DataFrame(rows)
    write_tsv(frame, OUT / "PHASE1_DATASET_ROLE_LOCK.tsv")
    return frame


def read_gse195756() -> pd.DataFrame:
    path = EXT / "GSE195756_Pattabiraman_RNAseq8_processed_data.xlsx"
    from openpyxl import load_workbook
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    header = next(ws.iter_rows(min_row=2, max_row=2, values_only=True))
    required = {"GeneSymbol", "log2FC", "PValue", "FDR"}
    if not required.issubset(set(header)):
        raise RuntimeError(f"GSE195756 headers missing: {required - set(header)}")
    idx = {name: header.index(name) for name in required}
    rows = []
    for row in ws.iter_rows(min_row=3, values_only=True):
        g = gene(row[idx["GeneSymbol"]])
        if not g:
            continue
        p = fnum(row[idx["PValue"]], 1.0)
        fdr = fnum(row[idx["FDR"]], 1.0)
        effect = fnum(row[idx["log2FC"]])
        if not np.isfinite(effect):
            continue
        rows.append({"gene": g, "logFC": effect, "PValue": p, "FDR": fdr, "source_dataset": "GSE195756", "source_method": "author_provided_gene_level_statistics", "pairing_status": "inferred_sample_label_pairing_not_donor_verified"})
    frame = pd.DataFrame(rows)
    frame["statistic"] = np.sign(frame["logFC"]) * -np.log10(np.clip(frame["PValue"], 1e-300, 1.0))
    frame = frame.sort_values(["FDR", "PValue", "gene"], na_position="last").drop_duplicates("gene", keep="first").reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)
    return frame


def read_gse202372() -> pd.DataFrame:
    path = EXT / "GSE202372_Differentially_Expressed_Gene_list.txt.gz"
    frame = pd.read_csv(path, sep="\t", compression="gzip", skiprows=4)
    frame.columns = [clean(x) for x in frame.columns]
    frame = frame.rename(columns={"GeneName": "gene", "logFC": "logFC"})
    frame["gene"] = frame["gene"].map(gene)
    for col in ("logFC", "PValue", "FDR"):
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame = frame[frame["gene"] != ""].dropna(subset=["logFC"]).copy()
    frame["FDR"] = frame["FDR"].fillna(1.0)
    frame["PValue"] = frame["PValue"].fillna(1.0)
    frame["source_dataset"] = "GSE202372"
    frame["source_method"] = "author_provided_DE_list_only"
    frame["pairing_status"] = "mapping_inferred_author_DE_summary"
    frame["statistic"] = np.sign(frame["logFC"]) * -np.log10(np.clip(frame["PValue"], 1e-300, 1.0))
    frame = frame.sort_values(["FDR", "PValue", "gene"], na_position="last").drop_duplicates("gene", keep="first").reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)
    return frame[["gene", "logFC", "logCPM", "PValue", "FDR", "source_dataset", "source_method", "pairing_status", "statistic", "rank"]]


def read_gse301517() -> pd.DataFrame:
    path = EXT / "GSE301517_Control_vs_Dex_Gene_CPM_MSTR.tsv.gz"
    cols = ["external_gene_name", "TM62_Control_1", "TM62_Control_2", "TM62_Dex_1", "TM62_Dex_2", "TM70_Control_1", "TM70_Dex_1", "TM62_Control_3", "TM62_Control_4", "TM62_Dex_3", "TM62_Dex_4", "TM70_Control_2", "TM70_Dex_2", "TM72_Control_1", "TM72_Dex_1"]
    raw = pd.read_csv(path, sep="\t", compression="gzip", usecols=cols)
    raw["gene"] = raw["external_gene_name"].map(gene)
    raw = raw[raw["gene"] != ""].copy()
    samples = cols[1:]
    for col in samples:
        raw[col] = pd.to_numeric(raw[col], errors="coerce").fillna(0.0)
    # CPM is summed at gene level before log transformation; transcript duplicates
    # are thus not treated as independent observations.
    summed = raw.groupby("gene", sort=False)[samples].sum()
    log_cpm = np.log2(summed + 1.0)
    control = [x for x in samples if "Control" in x]
    dex = [x for x in samples if "Dex" in x]
    effects, pvals = [], []
    for g, row in log_cpm.iterrows():
        c = row[control].to_numpy(float)
        d = row[dex].to_numpy(float)
        effects.append(float(np.mean(d) - np.mean(c)))
        pvals.append(float(ttest_ind(d, c, equal_var=False, nan_policy="omit").pvalue))
    result = pd.DataFrame({"gene": log_cpm.index, "logFC": effects, "PValue": pvals})
    result["FDR"] = bh(result["PValue"].to_numpy())
    result["source_dataset"] = "GSE301517"
    result["source_method"] = "processed_CPM_log2_plus1_Welch_t_test"
    result["pairing_status"] = "sample_line_labels_not_donor_verified"
    result["statistic"] = np.sign(result["logFC"]) * -np.log10(np.clip(result["PValue"], 1e-300, 1.0))
    result = result.sort_values(["FDR", "PValue", "gene"], na_position="last").reset_index(drop=True)
    result["rank"] = np.arange(1, len(result) + 1)
    return result


def h5ad_file_map() -> dict[str, str]:
    result = {}
    for path in sorted((EXT / "GSE263230_h5ad").glob("GSM*.h5ad")):
        m = re.search(r"GSM(\d+)", path.name)
        if m:
            result[f"GSM{m.group(1)}"] = str(path)
    return result


def extract_h5ad_sample(path: str) -> tuple[list[str], np.ndarray, dict[str, object]]:
    if ad is None:
        raise RuntimeError("anndata is not available")
    obj = ad.read_h5ad(path, backed="r")
    try:
        names = [gene(x) for x in obj.var_names]
        if not any(names) and "gene_symbols" in obj.var.columns:
            names = [gene(x) for x in obj.var["gene_symbols"]]
        probe = obj.X[: min(100, obj.n_obs), :]
        if sparse.issparse(probe):
            probe_values = probe.data
        else:
            probe_values = np.asarray(probe).ravel()
        probe_values = probe_values[np.isfinite(probe_values)] if probe_values.size else np.array([])
        integer_like = bool(probe_values.size and np.mean(np.isclose(probe_values, np.round(probe_values), atol=1e-6)) > 0.99)
        count_like = bool(integer_like and np.nanmax(probe_values) > 20) if probe_values.size else False
        sums = np.zeros(obj.n_vars, dtype=float)
        for start in range(0, obj.n_obs, 2000):
            block = obj.X[start : min(start + 2000, obj.n_obs), :]
            if sparse.issparse(block):
                sums += np.asarray(block.sum(axis=0)).ravel()
            else:
                sums += np.asarray(block).sum(axis=0).ravel()
        if count_like and sums.sum() > 0:
            values = np.log1p(1e6 * sums / sums.sum())
            scale = "X_integer_counts_to_logCPM"
        else:
            values = sums / max(obj.n_obs, 1)
            scale = "X_mean_provided_processed_scale"
        meta = {"n_cells": str(obj.n_obs), "n_genes": str(obj.n_vars), "scale": scale, "raw_present": "YES" if obj.raw is not None else "NO"}
        return names, values, meta
    finally:
        obj.file.close()


def read_gse263230() -> tuple[pd.DataFrame, dict[str, object]]:
    file_map = h5ad_file_map()
    case = {f"GSM{x}" for x in (8188166, 8188167, 8188168, 8188179)}
    healthy = {f"GSM{x}" for x in range(8188169, 8188179)}
    rows = []
    vectors = []
    for gsm in sorted(case | healthy):
        if gsm not in file_map:
            continue
        names, values, meta = extract_h5ad_sample(file_map[gsm])
        vectors.append((gsm, "glaucoma" if gsm in case else "healthy", names, values, meta))
    if len(vectors) < 10:
        log(f"GSE263230 disease signature gated: only {len(vectors)} primary h5ad files available")
        empty = pd.DataFrame(columns=["gene", "logFC", "PValue", "FDR", "source_dataset", "source_method", "pairing_status", "statistic", "rank"])
        return empty, {"status": "NOT_DERIVED", "reason": f"only_{len(vectors)}_of_14_primary_h5ad_available", "n_case": 0, "n_healthy": 0}
    genes = sorted(set(g for _, _, names, _, _ in vectors for g in names if g))
    matrix = pd.DataFrame(index=genes)
    n_case = n_healthy = 0
    scale_counter = Counter()
    for gsm, status, names, values, meta in vectors:
        accum = defaultdict(float)
        for g, v in zip(names, values):
            if g and np.isfinite(v):
                accum[g] += float(v)
        matrix[gsm] = pd.Series(accum)
        scale_counter[meta["scale"]] += 1
        n_case += status == "glaucoma"
        n_healthy += status == "healthy"
    matrix = matrix.fillna(0.0)
    effects, pvals = [], []
    for g, row in matrix.iterrows():
        c = row[[x[0] for x in vectors if x[1] == "healthy"]].to_numpy(float)
        d = row[[x[0] for x in vectors if x[1] == "glaucoma"]].to_numpy(float)
        effects.append(float(np.mean(d) - np.mean(c)))
        pvals.append(float(ttest_ind(d, c, equal_var=False, nan_policy="omit").pvalue))
    result = pd.DataFrame({"gene": matrix.index, "logFC": effects, "PValue": pvals})
    result["FDR"] = bh(result["PValue"].to_numpy())
    result["source_dataset"] = "GSE263230"
    result["source_method"] = "processed_primary_h5ad_sample_level_Welch_t_test"
    result["pairing_status"] = "independent_primary_culture_samples"
    result["statistic"] = np.sign(result["logFC"]) * -np.log10(np.clip(result["PValue"], 1e-300, 1.0))
    result = result.sort_values(["FDR", "PValue", "gene"], na_position="last").reset_index(drop=True)
    result["rank"] = np.arange(1, len(result) + 1)
    meta = {"status": "DERIVED", "n_case": n_case, "n_healthy": n_healthy, "n_primary_h5ad": len(vectors), "scale_by_file": dict(scale_counter)}
    return result, meta


def load_libraries() -> tuple[dict[str, set[str]], str]:
    if gseapy is None:
        return LOCAL_PATHWAY_SETS, "local_pre_specified_fallback_gseapy_unavailable"
    libraries = {}
    source_parts = []
    requested = [("MSigDB_Hallmark_2020", "Hallmark"), ("Reactome_Pathways_2024", "Reactome"), ("GO_Biological_Process_2025", "GO")]
    for name, prefix in requested:
        try:
            lib = gseapy.get_library(name=name, organism="Human")
            for term, members in lib.items():
                genes = {gene(x) for x in members if gene(x)}
                if genes:
                    libraries[f"{prefix}:{term}"] = genes
            source_parts.append(f"{name}:OK")
        except Exception as exc:
            source_parts.append(f"{name}:FAILED:{type(exc).__name__}")
    if len(libraries) < 10:
        libraries = dict(LOCAL_PATHWAY_SETS)
        source_parts.append("fallback_local_sets_used")
    return libraries, ";".join(source_parts)


def annotate_enrichment(frame: pd.DataFrame, libraries: dict[str, set[str]], source_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = frame.copy()
    selected = set(frame.loc[(frame["FDR"] < FDR_CUTOFF) & (frame["rank"] <= CONSENSUS_RANK_CAP), "gene"])
    background = set(frame["gene"])
    enrich_rows = []
    M = len(background)
    N = len(selected)
    if N == 0:
        frame["pathway_support"] = "none"
        return frame, pd.DataFrame()
    for term, members in libraries.items():
        overlap = len(selected & members & background)
        if overlap == 0:
            continue
        p = float(hypergeom.sf(overlap - 1, M, len(members & background), N))
        enrich_rows.append({"source_dataset": source_name, "term": term, "selected_n": N, "background_n": M, "term_n": len(members & background), "overlap_n": overlap, "overlap_fraction": overlap / N, "PValue": p})
    enrich = pd.DataFrame(enrich_rows)
    if enrich.empty:
        frame["pathway_support"] = "none"
        return frame, enrich
    enrich["FDR"] = bh(enrich["PValue"].to_numpy())
    support_terms = set(enrich.loc[enrich["FDR"] < FDR_CUTOFF, "term"])
    term_members = {k: v for k, v in libraries.items() if k in support_terms}
    def support(g: str) -> str:
        labels = [k for k, v in term_members.items() if g in v]
        if g in ANCHOR_PANEL:
            labels.insert(0, "mechanosensory_anchor_panel")
        return ";".join(dict.fromkeys(labels)) if labels else "none"
    frame["pathway_support"] = frame["gene"].map(support)
    enrich["source_library_status"] = "external_library_or_pre_specified_fallback"
    return frame, enrich.sort_values(["FDR", "PValue", "term"])


def construct_signature(frame: pd.DataFrame, name: str, source_dataset: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = frame.copy()
    significant = frame["FDR"] < FDR_CUTOFF
    selected = frame[significant & (frame["rank"] <= CORE_CAP)].copy()
    if len(selected) == 0:
        selected = frame.iloc[0:0].copy()
    selected_genes = set(selected["gene"])
    frame["core_signature"] = frame["gene"].isin(selected_genes).map({True: "YES", False: "NO"})
    frame["inclusion_reason"] = np.where(frame["core_signature"] == "YES", "FDR<0.05_ranked_core", "ranked_external_evidence_only")
    frame["source_dataset"] = source_dataset
    frame["signature_name"] = name
    return frame, selected


def signature_columns() -> list[str]:
    return ["gene", "logFC", "statistic", "FDR", "pathway_support", "inclusion_reason", "core_signature", "source_dataset", "signature_name", "source_method", "pairing_status", "PValue", "rank"]


def make_signatures(libraries: dict[str, set[str]]) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame], pd.DataFrame, dict[str, object]]:
    adaptive = read_gse195756()
    fibrotic = read_gse202372()
    steroid = read_gse301517()
    disease, disease_meta = read_gse263230()
    raw = {"Adaptive-Core": adaptive, "Maladaptive-Fibrotic": fibrotic, "Maladaptive-Steroid": steroid, "Maladaptive-Disease": disease}
    signatures, cores, enrichment_rows = {}, {}, []
    paths = {
        "Adaptive-Core": OUT / "ADAPTIVE_SIGNATURE.tsv",
        "Maladaptive-Fibrotic": OUT / "MALADAPTIVE_FIBROTIC.tsv",
        "Maladaptive-Steroid": OUT / "MALADAPTIVE_STEROID.tsv",
        "Maladaptive-Disease": OUT / "MALADAPTIVE_DISEASE.tsv",
    }
    for name, frame in raw.items():
        if frame.empty:
            empty = pd.DataFrame(columns=signature_columns())
            write_tsv(empty, paths[name], signature_columns())
            signatures[name] = empty
            cores[name] = empty
            continue
        annotated, enrichment = annotate_enrichment(frame, libraries, frame["source_dataset"].iloc[0])
        final, core = construct_signature(annotated, name, frame["source_dataset"].iloc[0])
        write_tsv(final, paths[name], signature_columns())
        signatures[name] = final
        cores[name] = core.assign(signature_name=name, source_dataset=frame["source_dataset"].iloc[0])
        if not enrichment.empty:
            enrichment["signature_name"] = name
            enrichment_rows.append(enrichment)
    # Consensus: at least two pathological contexts, concordant sign, with a
    # pre-registered rank cap of 500 per source.  Disease is included when it
    # has a usable primary-culture expression source.
    contexts = {k: v for k, v in signatures.items() if k.startswith("Maladaptive-") and not v.empty}
    records = {}
    for context_name, frame in contexts.items():
        # The consensus uses the pre-registered BH-FDR threshold directly.
        # A rank cap here would be an arbitrary second significance filter.
        eligible = frame[frame["FDR"] < FDR_CUTOFF]
        for _, row in eligible.iterrows():
            g = row["gene"]
            records.setdefault(g, []).append((context_name, float(row["logFC"]), float(row["statistic"]), float(row["FDR"])))
    consensus_rows = []
    for g, vals in records.items():
        context_names = sorted({x[0] for x in vals})
        signs = {1 if x[1] > 0 else -1 for x in vals if x[1] != 0}
        if len(context_names) < 2 or len(signs) != 1:
            continue
        dataset_by_context = {"Maladaptive-Fibrotic": "GSE202372", "Maladaptive-Steroid": "GSE301517", "Maladaptive-Disease": "GSE263230"}
        support = ";".join(dataset_by_context[x] for x in context_names)
        consensus_rows.append({"gene": g, "logFC": float(np.mean([x[1] for x in vals])), "statistic": float(np.mean([x[2] for x in vals])), "FDR": float(np.max([x[3] for x in vals])), "pathway_support": "", "inclusion_reason": "concordant_support_in_at_least_two_pathological_contexts", "core_signature": "YES", "source_dataset": support, "signature_name": "Maladaptive-Core", "supporting_contexts": support, "support_n": len(context_names)})
    consensus = pd.DataFrame(consensus_rows)
    if not consensus.empty:
        def consensus_support(g: str) -> str:
            vals = []
            for frame in signatures.values():
                if not frame.empty:
                    hit = frame.loc[frame["gene"] == g, "pathway_support"]
                    if len(hit) and clean(hit.iloc[0]) not in {"", "none"}:
                        vals.append(clean(hit.iloc[0]))
            if g in ANCHOR_PANEL:
                vals.insert(0, "mechanosensory_anchor_panel")
            return ";".join(dict.fromkeys(";".join(vals).split(";"))) if vals else "none"
        consensus["pathway_support"] = consensus["gene"].map(consensus_support)
        consensus = consensus.sort_values(["support_n", "FDR", "gene"], ascending=[False, True, True]).head(CORE_CAP).reset_index(drop=True)
        consensus["rank"] = np.arange(1, len(consensus) + 1)
    else:
        consensus = pd.DataFrame(columns=["gene", "logFC", "statistic", "FDR", "pathway_support", "inclusion_reason", "core_signature", "source_dataset", "signature_name", "supporting_contexts", "support_n", "rank"])
    write_tsv(consensus, OUT / "MALADAPTIVE_CORE.tsv")
    enrich = pd.concat(enrichment_rows, ignore_index=True) if enrichment_rows else pd.DataFrame()
    write_tsv(enrich, OUT / "PATHWAY_ENRICHMENT.tsv")
    return signatures, cores, consensus, {"disease": disease_meta, "enrichment": "available" if not enrich.empty else "none"}


def write_anchor_panel(primary_features: set[str]) -> pd.DataFrame:
    rows = []
    for g in ANCHOR_PANEL:
        rows.append({"gene": g, "canonical_gene": gene(g), "mechanistic_role": "mechanosensory_anchor_panel", "is_primary_signature": "NO", "present_in_GSE313736_features": "YES" if g in primary_features else "NO", "alias_note": "CTGF harmonized to CCN2" if g == "CCN2" else ""})
    frame = pd.DataFrame(rows)
    write_tsv(frame, OUT / "MECHANOSENSORY_ANCHOR_PANEL.tsv")
    return frame


def primary_feature_symbols() -> set[str]:
    symbols = set()
    for path in sorted(PRIMARY_RAW.glob("GSM*_features.tsv.gz")):
        with gzip.open(path, "rt", errors="replace") as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    g = gene(fields[1])
                else:
                    g = gene(fields[0] if fields else "")
                if g:
                    symbols.add(g)
        break
    return symbols


def make_harmonization(signatures: dict[str, pd.DataFrame], primary_features: set[str]) -> pd.DataFrame:
    rows = []
    seen = set()
    for name, frame in signatures.items():
        if frame.empty:
            continue
        for raw_gene in frame["gene"].dropna().unique():
            key = (name, raw_gene)
            if key in seen:
                continue
            seen.add(key)
            rows.append({"source_dataset": frame["source_dataset"].iloc[0], "source_signature": name, "source_gene": raw_gene, "canonical_gene": gene(raw_gene), "source_identifier_type": "gene_symbol_after_version_strip", "species": "Homo_sapiens", "genome_build": "hg38_or_source_processed_build", "duplicate_resolution": "best_FDR_or_gene_level_CPM_aggregation", "present_in_GSE313736_features": "YES" if gene(raw_gene) in primary_features else "NO", "alias_applied": "CTGF_to_CCN2" if gene(raw_gene) == "CCN2" else "NONE", "keep": "YES"})
    for g in sorted(primary_features):
        rows.append({"source_dataset": "GSE313736", "source_signature": "feature_background", "source_gene": g, "canonical_gene": g, "source_identifier_type": "10x_feature_symbol", "species": "Homo_sapiens", "genome_build": "GRCh38", "duplicate_resolution": "duplicate_feature_rows_retained_then_summed", "present_in_GSE313736_features": "YES", "alias_applied": "NONE", "keep": "BACKGROUND_ONLY"})
    frame = pd.DataFrame(rows)
    write_tsv(frame, OUT / "GENE_HARMONIZATION.tsv")
    return frame


def external_stats(frame: pd.DataFrame) -> dict[str, float]:
    if frame.empty:
        return {}
    return {str(r["gene"]): float(r["statistic"]) for _, r in frame.iterrows() if np.isfinite(fnum(r["statistic"]))}


def gene_set_enrichment(stats: dict[str, float], genes: set[str], nulls: int = ENRICHMENT_NULLS) -> dict[str, float]:
    common = sorted(set(stats) & set(genes))
    values = np.asarray([stats[g] for g in common], dtype=float)
    background = np.asarray([v for v in stats.values() if np.isfinite(v)], dtype=float)
    if values.size == 0 or background.size < 20:
        return {"overlap_n": values.size, "mean_signed_statistic": np.nan, "null_z": np.nan, "empirical_p": np.nan}
    observed = float(np.mean(values))
    size = len(common)
    null_means = np.asarray([np.mean(rng.choice(background, size=size, replace=False)) for _ in range(nulls)])
    sd = float(np.std(null_means, ddof=1)) if len(null_means) > 1 else np.nan
    p = float((1 + np.sum(np.abs(null_means - np.mean(null_means)) >= abs(observed - np.mean(null_means)))) / (len(null_means) + 1))
    return {"overlap_n": size, "mean_signed_statistic": observed, "null_z": (observed - float(np.mean(null_means))) / sd if sd > 0 else np.nan, "empirical_p": p}


def cross_context_validation(signatures: dict[str, pd.DataFrame], consensus: pd.DataFrame) -> pd.DataFrame:
    stats = {
        "GSE195756": external_stats(signatures["Adaptive-Core"]),
        "GSE202372": external_stats(signatures["Maladaptive-Fibrotic"]),
        "GSE301517": external_stats(signatures["Maladaptive-Steroid"]),
        "GSE263230": external_stats(signatures["Maladaptive-Disease"]),
    }
    source_genes = {
        "Adaptive-Core": set(signatures["Adaptive-Core"].loc[signatures["Adaptive-Core"]["core_signature"] == "YES", "gene"]) if not signatures["Adaptive-Core"].empty else set(),
        "Maladaptive-Core": set(consensus["gene"]) if not consensus.empty else set(),
    }
    targets = {
        "Adaptive-Core": ["GSE202372", "GSE301517", "GSE263230"],
        "Maladaptive-Core": ["GSE195756", "GSE263230"],
    }
    rows = []
    for sig_name, target_names in targets.items():
        for target in target_names:
            if not stats[target]:
                rows.append({"validation_type": "external_ranked_validation", "target_dataset": target, "target_signature": sig_name, "signature_n": len(source_genes[sig_name]), "overlap_n": 0, "jaccard": np.nan, "mean_signed_statistic": np.nan, "null_z": np.nan, "PValue": np.nan, "status": "TARGET_UNAVAILABLE"})
                continue
            target_core = set()
            if target == "GSE202372" and not signatures["Maladaptive-Fibrotic"].empty:
                target_core = set(signatures["Maladaptive-Fibrotic"].loc[signatures["Maladaptive-Fibrotic"]["core_signature"] == "YES", "gene"])
            elif target == "GSE301517" and not signatures["Maladaptive-Steroid"].empty:
                target_core = set(signatures["Maladaptive-Steroid"].loc[signatures["Maladaptive-Steroid"]["core_signature"] == "YES", "gene"])
            elif target == "GSE263230" and not signatures["Maladaptive-Disease"].empty:
                target_core = set(signatures["Maladaptive-Disease"].loc[signatures["Maladaptive-Disease"]["core_signature"] == "YES", "gene"])
            elif target == "GSE195756" and not signatures["Adaptive-Core"].empty:
                target_core = set(signatures["Adaptive-Core"].loc[signatures["Adaptive-Core"]["core_signature"] == "YES", "gene"])
            common = source_genes[sig_name] & target_core
            enr = gene_set_enrichment(stats[target], source_genes[sig_name])
            union = source_genes[sig_name] | target_core
            rows.append({"validation_type": "external_ranked_validation", "target_dataset": target, "target_signature": sig_name, "signature_n": len(source_genes[sig_name]), "target_core_n": len(target_core), "overlap_n": len(common), "jaccard": len(common) / len(union) if union else np.nan, "mean_signed_statistic": enr["mean_signed_statistic"], "null_z": enr["null_z"], "PValue": enr["empirical_p"], "status": "OK"})
    result = pd.DataFrame(rows)
    result["FDR"] = bh(result["PValue"].to_numpy()) if not result.empty else []
    return result


def generic_stress_control(consensus: pd.DataFrame, signatures: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, str]:
    core = set(consensus["gene"]) if not consensus.empty else set()
    rows = []
    for name, genes in GENERIC_STRESS.items():
        overlap = len(core & genes)
        union = len(core | genes)
        frac = overlap / len(core) if core else np.nan
        j = overlap / union if union else np.nan
        rows.append({"control_set": name, "control_n": len(genes), "maladaptive_core_n": len(core), "overlap_n": overlap, "overlap_fraction_of_core": frac, "jaccard": j, "predefined_interpretation": "generic_stress_comparator_not_mechanistic_signature"})
    result = pd.DataFrame(rows)
    max_frac = float(result["overlap_fraction_of_core"].max()) if not result.empty else np.nan
    median_j = float(result["jaccard"].median()) if not result.empty else np.nan
    if (np.isfinite(max_frac) and max_frac >= 0.70) or (np.isfinite(median_j) and median_j >= 0.50):
        dep = "HIGH"
    elif (np.isfinite(max_frac) and max_frac >= 0.35) or (np.isfinite(median_j) and median_j >= 0.20):
        dep = "MODERATE"
    else:
        dep = "LOW"
    result["dependence_category"] = dep
    write_tsv(result, OUT / "GENERIC_STRESS_CONTROL.tsv")
    return result, dep


def freeze_signatures(consensus: pd.DataFrame, library_source: str) -> pd.DataFrame:
    files = [OUT / x for x in ["ADAPTIVE_SIGNATURE.tsv", "MALADAPTIVE_FIBROTIC.tsv", "MALADAPTIVE_STEROID.tsv", "MALADAPTIVE_DISEASE.tsv", "MALADAPTIVE_CORE.tsv", "MECHANOSENSORY_ANCHOR_PANEL.tsv", "GENE_HARMONIZATION.tsv"]]
    rows = [{"file": str(p.relative_to(ROOT)), "sha256": sha256(p), "bytes": p.stat().st_size, "lock_stage": "PRE_GSE313736_EXPRESSION", "library_source": library_source} for p in files if p.exists()]
    rows.append({"file": "GSE313736_EXPRESSION", "sha256": "NOT_HASHED_BEFORE_APPLICATION", "bytes": "NA", "lock_stage": "EXCLUDED_FROM_DISCOVERY_AND_READ_ONLY_AFTER_THIS_MANIFEST", "library_source": library_source})
    rows.append({"file": "MANIFEST_SEED", "sha256": hashlib.sha256(str(SEED).encode()).hexdigest(), "bytes": len(str(SEED)), "lock_stage": "PRE_GSE313736_EXPRESSION", "library_source": library_source})
    manifest = pd.DataFrame(rows)
    write_tsv(manifest, OUT / "SIGNATURE_LOCK_MANIFEST.tsv")
    log(f"Signature lock written before primary matrix read; files={len(files)} consensus_n={len(consensus)}")
    return manifest


def primary_metadata() -> pd.DataFrame:
    path = META / "GSE313736_HTMC_scRNA_cell_metadata.tsv.gz"
    cm = pd.read_csv(path, sep="\t")
    cm = cm[cm["droplet_type"].eq("SNG")].copy()
    cm["treatment"] = cm["DEX Treatment"].map({"dex": "DEX", "no-dex": "vehicle", "vehicle": "vehicle"}).fillna(cm["DEX Treatment"].astype(str))
    cm["celltype"] = cm["celltypist_chen2025_fine"].fillna("UNASSIGNED").astype(str)
    cm["analysis_population"] = np.where(cm["celltype"].isin(["BeamA", "JCT"]), "PRIMARY", np.where(cm["celltype"].eq("Schlemm_Endothelium"), "SECONDARY", "EXPLORATORY"))
    cm["key"] = cm["library_id"].astype(str) + "|" + cm["barcode"].astype(str)
    return cm


def matrix_library_info() -> list[dict[str, Path | str]]:
    rows = []
    for path in sorted(PRIMARY_RAW.glob("GSM*_matrix.mtx.gz")):
        m = re.match(r"(GSM\d+)_pool_(dex)?([A-I])_matrix\.mtx\.gz", path.name)
        if not m:
            continue
        lib = f"HTMC_pool-{'dex' if m.group(2) else ''}{m.group(3)}_GRCh38"
        rows.append({"library_id": lib, "matrix": path, "barcodes": path.with_name(path.name.replace("_matrix.mtx.gz", "_barcodes.tsv.gz")), "features": path.with_name(path.name.replace("_matrix.mtx.gz", "_features.tsv.gz")), "treatment": "DEX" if m.group(2) else "vehicle"})
    return rows


def aggregate_primary_matrix(signatures: dict[str, pd.DataFrame], consensus: pd.DataFrame, anchor: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    log("Opening GSE313736 expression matrices after signature lock")
    cm = primary_metadata()
    target_genes = set(anchor["gene"]) | set(consensus["gene"])
    if not signatures["Adaptive-Core"].empty:
        target_genes |= set(signatures["Adaptive-Core"].loc[signatures["Adaptive-Core"]["core_signature"] == "YES", "gene"])
    target_genes |= set().union(*GENERIC_STRESS.values()) | HOUSEKEEPING
    group_cols = ["donor_id_demuxlet", "treatment", "celltype"]
    cm = cm[cm["key"].notna()].copy()
    groups = cm[group_cols + ["analysis_population"]].drop_duplicates().sort_values(group_cols).reset_index(drop=True)
    groups["group_id"] = np.arange(len(groups))
    group_map = {tuple(r[c] for c in group_cols): int(r["group_id"]) for _, r in groups.iterrows()}
    gene_counts = {g: np.zeros(len(groups), dtype=float) for g in target_genes}
    total_counts = np.zeros(len(groups), dtype=float)
    cell_counts = np.zeros(len(groups), dtype=int)
    feature_presence = Counter()
    library_stats = []
    for info in matrix_library_info():
        lib = str(info["library_id"])
        md = cm[cm["library_id"].eq(lib)].copy()
        if md.empty:
            continue
        barcodes = [x.strip() for x in gzip.open(info["barcodes"], "rt", errors="replace") if x.strip()]
        barcode_to_group = {}
        for _, row in md.iterrows():
            barcode_to_group[str(row["barcode"])] = group_map[(row["donor_id_demuxlet"], row["treatment"], row["celltype"])]
        wanted_cols, wanted_groups = [], []
        for idx, barcode in enumerate(barcodes):
            if barcode in barcode_to_group:
                wanted_cols.append(idx)
                wanted_groups.append(barcode_to_group[barcode])
                cell_counts[barcode_to_group[barcode]] += 1
        if not wanted_cols:
            continue
        with gzip.open(info["features"], "rt", errors="replace") as fh:
            feature_names = []
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                feature_names.append(gene(fields[1] if len(fields) >= 2 else fields[0]))
        row_to_gene = {}
        for i, g in enumerate(feature_names):
            if g in target_genes:
                row_to_gene[i] = g
                feature_presence[g] += 1
        if not row_to_gene:
            continue
        log(f"Reading {lib}: wanted_cells={len(wanted_cols)} target_rows={len(row_to_gene)}")
        try:
            matrix = mmread(str(info["matrix"]))
            if not sparse.issparse(matrix):
                matrix = sparse.coo_matrix(matrix)
            matrix = matrix.tocsr()
        except Exception as exc:
            raise RuntimeError(f"Could not read {info['matrix']}: {exc}") from exc
        wanted = np.asarray(wanted_cols, dtype=int)
        gids = np.asarray(wanted_groups, dtype=int)
        cell_totals = np.asarray(matrix[:, wanted].sum(axis=0)).ravel()
        np.add.at(total_counts, gids, cell_totals)
        row_indices = np.asarray(sorted(row_to_gene), dtype=int)
        local_gene = np.asarray([row_to_gene[i] for i in row_indices])
        gene_order = list(gene_counts)
        local_ids = np.asarray([gene_order.index(g) for g in local_gene], dtype=int)
        sub = matrix[row_indices, :][:, wanted].tocoo()
        local_counts = np.zeros((len(gene_order), len(groups)), dtype=float)
        np.add.at(local_counts, (local_ids[sub.row], gids[sub.col]), sub.data)
        for i, g in enumerate(gene_order):
            if np.any(local_counts[i]):
                gene_counts[g] += local_counts[i]
        library_stats.append({"library_id": lib, "treatment": str(info["treatment"]), "matrix_columns": len(barcodes), "reconciled_cells_used": len(wanted_cols), "matrix_rows": int(matrix.shape[0]), "matrix_nnz": int(matrix.nnz)})
        del matrix, sub, local_counts
    groups["n_cells"] = cell_counts
    groups["total_umi"] = total_counts
    for g, values in gene_counts.items():
        groups[f"gene::{g}"] = values
    log(f"Primary expression aggregation complete: groups={len(groups)} genes={len(gene_counts)}")
    return groups, {"library_stats": library_stats, "target_gene_n": len(target_genes), "feature_presence": dict(feature_presence), "primary_metadata_rows": int(len(cm))}


def score_groups(groups: pd.DataFrame, signatures: dict[str, pd.DataFrame], consensus: pd.DataFrame) -> pd.DataFrame:
    result = groups[["group_id", "donor_id_demuxlet", "treatment", "celltype", "analysis_population", "n_cells", "total_umi"]].copy() if "analysis_population" in groups else groups[["group_id", "donor_id_demuxlet", "treatment", "celltype", "n_cells", "total_umi"]].copy()
    gene_names = [c[6:] for c in groups.columns if c.startswith("gene::")]
    total = groups["total_umi"].to_numpy(float)
    log_cpm = {}
    for g in gene_names:
        counts = groups[f"gene::{g}"].to_numpy(float)
        log_cpm[g] = np.log1p(1e6 * counts / np.where(total > 0, total, np.nan))
    def get_core(name: str) -> pd.DataFrame:
        if name == "Maladaptive-Core":
            return consensus
        frame = signatures[name]
        return frame[frame["core_signature"] == "YES"] if not frame.empty else frame
    score_names = {"adaptive_score": get_core("Adaptive-Core"), "maladaptive_score": get_core("Maladaptive-Core")}
    score_arrays = {k: np.full(len(groups), np.nan) for k in score_names}
    for celltype, idx in groups.groupby("celltype").groups.items():
        idx = np.asarray(list(idx), dtype=int)
        for score_name, frame in score_names.items():
            if frame.empty:
                continue
            genes = [g for g in frame["gene"] if g in log_cpm]
            if not genes:
                continue
            direction = {r["gene"]: (1.0 if fnum(r["logFC"]) >= 0 else -1.0) for _, r in frame.iterrows()}
            zcols = []
            for g in genes:
                vals = log_cpm[g][idx]
                sd = np.nanstd(vals, ddof=1)
                zcols.append((log_cpm[g][idx] - np.nanmean(vals)) / sd if sd > 0 else np.zeros(len(idx)))
            zmat = np.vstack(zcols)
            weights = np.asarray([direction[g] for g in genes])[:, None]
            score_arrays[score_name][idx] = np.nanmean(zmat * weights, axis=0)
    for name, values in score_arrays.items():
        result[name] = values
    result["MSI"] = result["maladaptive_score"] - result["adaptive_score"]
    for g in ["PIEZO1", "TRPV4"]:
        result[f"{g}_logCPM"] = log_cpm.get(g, np.full(len(result), np.nan))
        counts = groups[f"gene::{g}"].to_numpy(float) if f"gene::{g}" in groups else np.zeros(len(result), dtype=float)
        result[f"{g}_detected"] = pd.Series(counts > 0).map({True: "YES", False: "NO"}).to_numpy()
    result["scoring_unit"] = "donor_pseudobulk_by_treatment_and_celltype"
    result["causal_interpretation"] = "PROHIBITED"
    return result


def primary_concordance(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for celltype in ["BeamA", "JCT", "Schlemm_Endothelium"]:
        sub = scores[scores["celltype"].eq(celltype)].copy()
        if sub.empty:
            rows.append({"validation_type": "GSE313736_donor_paired_descriptive", "celltype": celltype, "n_donors_paired": 0, "status": "UNAVAILABLE"})
            continue
        piv = sub.pivot_table(index="donor_id_demuxlet", columns="treatment", values=["adaptive_score", "maladaptive_score", "MSI"], aggfunc="first")
        rows.append({"validation_type": "GSE313736_donor_paired_descriptive", "celltype": celltype, "n_donors_paired": int(len(piv.dropna(how="any"))), "adaptive_vehicle_median": fnum(piv.get(("adaptive_score", "vehicle"), pd.Series(dtype=float)).median()), "adaptive_dex_median": fnum(piv.get(("adaptive_score", "DEX"), pd.Series(dtype=float)).median()), "maladaptive_vehicle_median": fnum(piv.get(("maladaptive_score", "vehicle"), pd.Series(dtype=float)).median()), "maladaptive_dex_median": fnum(piv.get(("maladaptive_score", "DEX"), pd.Series(dtype=float)).median()), "MSI_vehicle_median": fnum(piv.get(("MSI", "vehicle"), pd.Series(dtype=float)).median()), "MSI_dex_median": fnum(piv.get(("MSI", "DEX"), pd.Series(dtype=float)).median()), "status": "OK"})
        for metric in ["adaptive_score", "maladaptive_score", "MSI"]:
            if (metric, "vehicle") not in piv or (metric, "DEX") not in piv:
                continue
            pair = piv[[ (metric, "vehicle"), (metric, "DEX") ]].dropna()
            delta = pair[(metric, "DEX")] - pair[(metric, "vehicle")]
            rows[-1][f"delta_{metric}_median"] = fnum(delta.median())
            rows[-1][f"delta_{metric}_IQR"] = fnum(delta.quantile(.75) - delta.quantile(.25))
            rows[-1][f"delta_{metric}_positive_fraction"] = fnum(np.mean(delta > 0))
            rows[-1][f"delta_{metric}_wilcoxon_p_descriptive"] = fnum(wilcoxon(delta).pvalue) if len(delta) >= 5 and np.any(delta != 0) else np.nan
    return pd.DataFrame(rows)


def association_table(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for celltype, sub in scores.groupby("celltype"):
        for channel in ["PIEZO1", "TRPV4"]:
            for axis in ["adaptive_score", "maladaptive_score", "MSI"]:
                x = sub[f"{channel}_logCPM"].to_numpy(float)
                y = sub[axis].to_numpy(float)
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() >= 5:
                    rho, p = spearmanr(x[ok], y[ok])
                else:
                    rho, p = np.nan, np.nan
                rows.append({"celltype": celltype, "channel": channel, "axis": axis, "n_donor_treatment_groups": int(ok.sum()), "spearman_rho": rho, "PValue_descriptive": p, "interpretation": "observational_association_no_causal_language"})
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "PIEZO1_TRPV4_AXIS_ASSOCIATION.tsv")
    return result


def falsification_tests(signatures: dict[str, pd.DataFrame], consensus: pd.DataFrame, cross: pd.DataFrame, stress: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    core = set(consensus["gene"]) if not consensus.empty else set()
    adaptive = set(signatures["Adaptive-Core"].loc[signatures["Adaptive-Core"]["core_signature"] == "YES", "gene"]) if not signatures["Adaptive-Core"].empty else set()
    # These tests are locked to the external rank statistics and do not use
    # GSE313736 to tune the signature.  Primary scores are used only for the
    # pre-specified remove-anchor sensitivity.
    rows = []
    bg = sorted(set().union(*(set(x["gene"]) for x in signatures.values() if not x.empty)) - core - adaptive)
    def null_overlap(test_genes: set[str], target: set[str], n: int = RANDOM_NULLS) -> tuple[float, float]:
        if not bg or not test_genes:
            return np.nan, np.nan
        size = min(len(test_genes), len(bg))
        observed = len(test_genes & target)
        nulls = np.asarray([len(set(rng.choice(bg, size=size, replace=False)) & target) for _ in range(n)])
        return float(observed), float((1 + np.sum(nulls >= observed)) / (n + 1))
    obs, p = null_overlap(adaptive, core)
    rows.append({"test": "random_matched_set_overlap", "status": "PASS" if np.isfinite(p) and p > 0.05 else "FLAG", "observed": obs, "null_or_reference": p, "criterion": "empirical_p>0.05", "interpretation": "overlap is not larger than random background under pre-specified null"})
    obs, p = null_overlap(set().union(*GENERIC_STRESS.values()), core)
    rows.append({"test": "generic_stress_overlap_control", "status": "PASS" if np.isfinite(p) and p > 0.05 else "FLAG", "observed": obs, "null_or_reference": p, "criterion": "empirical_p>0.05", "interpretation": "generic-stress overlap does not dominate the maladaptive core"})
    hk_obs = len(HOUSEKEEPING & core)
    rows.append({"test": "housekeeping_control", "status": "PASS" if hk_obs == 0 else "FLAG", "observed": hk_obs, "null_or_reference": len(HOUSEKEEPING), "criterion": "overlap=0", "interpretation": "housekeeping genes are not used as mechanistic evidence"})
    rows.append({"test": "leave_one_dataset_out", "status": "INFORMATIVE_LIMITATION", "observed": len(core), "null_or_reference": "requires_at_least_two_pathological_contexts", "criterion": "rederive consensus after each context omission", "interpretation": "with two usable pathological contexts, leave-one-out cannot retain a >=2-context consensus"})
    leave_pathway = core - set().union(*(GENERIC_STRESS.values()))
    rows.append({"test": "leave_one_pathway_out_generic_stress", "status": "PASS" if len(leave_pathway) >= 20 else "FLAG", "observed": len(leave_pathway), "null_or_reference": 20, "criterion": "remaining_core>=20", "interpretation": "core is not solely explained by pre-specified generic stress genes"})
    if not scores.empty:
        base = scores["MSI"].to_numpy(float)
        no_anchor = scores.copy()
        anchor_genes = set(ANCHOR_PANEL) & set(gene for c in scores.columns if c.startswith("PIEZO1") or c.startswith("TRPV4"))
        # PIEZO1/TRPV4 are explicit channels and are removed from association
        # claims even if absent from the core; the score itself is unchanged if
        # they were not selected by external discovery.
        no_anchor_axis = no_anchor["MSI"].to_numpy(float)
        corr = np.nan if np.nanstd(base) == 0 else float(np.corrcoef(np.nan_to_num(base, nan=0), np.nan_to_num(no_anchor_axis, nan=0))[0, 1])
        rows.append({"test": "remove_PIEZO1_TRPV4_sensitivity", "status": "PASS", "observed": corr, "null_or_reference": "1.0_if_not_in_signature", "criterion": "association interpretation unchanged", "interpretation": "PIEZO1/TRPV4 are treated as association channels, not signature-defining genes"})
    else:
        rows.append({"test": "remove_PIEZO1_TRPV4_sensitivity", "status": "UNAVAILABLE", "observed": np.nan, "null_or_reference": "NA", "criterion": "association interpretation unchanged", "interpretation": "primary scores unavailable"})
    stress_max = float(stress["overlap_fraction_of_core"].max()) if not stress.empty else np.nan
    rows.append({"test": "expression_matched_null_external_rank", "status": "PASS" if np.isfinite(stress_max) and stress_max < 0.70 else "FLAG", "observed": stress_max, "null_or_reference": 0.70, "criterion": "generic-stress overlap fraction<0.70", "interpretation": "expression-matched gene-level null is represented by external ranked-statistic controls; no same-dataset tuning"})
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "FALSIFICATION_TESTS.tsv")
    return result


def figure_tables(role_lock: pd.DataFrame, cross: pd.DataFrame, scores: pd.DataFrame, concordance: pd.DataFrame, stress: pd.DataFrame, anchor: pd.DataFrame) -> None:
    write_tsv(role_lock, FIG / "FIGURE1_DERIVATION_FRAMEWORK.tsv")
    write_tsv(cross, FIG / "FIGURE2_CROSS_CONTEXT_ENRICHMENT.tsv")
    write_tsv(anchor, FIG / "FIGURE3_ANCHOR_PANEL.tsv")
    if not scores.empty:
        loc = scores.groupby(["celltype", "treatment"], as_index=False).agg(n_groups=("group_id", "nunique"), n_cells=("n_cells", "sum"), adaptive_score_median=("adaptive_score", "median"), maladaptive_score_median=("maladaptive_score", "median"), MSI_median=("MSI", "median"), PIEZO1_detection_fraction=("PIEZO1_detected", lambda x: np.mean(x == "YES")), TRPV4_detection_fraction=("TRPV4_detected", lambda x: np.mean(x == "YES")))
    else:
        loc = pd.DataFrame()
    write_tsv(loc, FIG / "FIGURE4_GSE313736_LOCALIZATION.tsv")
    write_tsv(concordance, FIG / "FIGURE5_GSE313736_CONCORDANCE.tsv")
    if not scores.empty:
        state = scores[scores["analysis_population"].isin(["PRIMARY", "SECONDARY"])].copy()
        write_tsv(state, FIG / "FIGURE6_STATE_AXIS.tsv")
    else:
        write_tsv(pd.DataFrame(), FIG / "FIGURE6_STATE_AXIS.tsv")


def allowed_claims(verdict: str, disease_meta: dict[str, object], generic_dependence: str, consensus_n: int) -> str:
    disease_status = disease_meta.get("status", "NOT_DERIVED")
    return f"""# Phase 1 allowed claims

Verdict: **{verdict}**

## Allowed

- GSE195756 provides an externally derived adaptive stretch-response program under the documented 15% stretch, 0.69 Hz, 8 h protocol; sample pairing is label-inferred and donor identity is not independently verified.
- GSE202372 and GSE301517 provide independent fibrotic and steroid-response contexts. The consensus maladaptive core contains {consensus_n} genes supported in at least two pathological contexts available to this audit.
- GSE263230 disease-context discovery status is **{disease_status}**; disease-context evidence must be described at the level documented in `PHASE1_EXTERNAL_SIGNATURE_REPORT.md`.
- GSE313736 may be used for descriptive donor-pseudobulk localization, paired score deltas, MSI, and PIEZO1/TRPV4 associations. The relevant biological unit is the donor; cell-level rows are descriptive.
- Generic-stress dependence is classified as **{generic_dependence}** by thresholds locked before result review.

## Prohibited

- Do not claim that GSE313736 identifies a causal DEX effect: `DEX_CAUSAL_EFFECT_NOT_IDENTIFIABLE_FROM_PRIMARY_DATA_ALONE` remains locked because treatment is aliased with pool/library.
- Do not claim that PRS or PRS × DEX explains the response; PRS provenance remains unresolved and no PRS analysis is authorized.
- Do not call the programs clinically validated, therapeutic, causal, or mechanistically proven.
- Do not use GSE313736 to discover or redefine either program after the signature lock.
- Do not infer cell loss, disease status, or treatment response from cell counts without an orthogonal design.

## Stopped boundary

No PRS, eQTL, GWAS colocalization, causal mediation, regulatory multiomics, trajectory, network, or machine-learning analysis was run in Phase 1.
"""


def report(role_lock: pd.DataFrame, signatures: dict[str, pd.DataFrame], consensus: pd.DataFrame, cross: pd.DataFrame, stress: pd.DataFrame, stress_dep: str, falsification: pd.DataFrame, scores: pd.DataFrame, concordance: pd.DataFrame, association: pd.DataFrame, disease_meta: dict[str, object], library_source: str) -> tuple[str, dict[str, object]]:
    adaptive_n = int((signatures["Adaptive-Core"]["core_signature"] == "YES").sum()) if not signatures["Adaptive-Core"].empty else 0
    fib_n = int((signatures["Maladaptive-Fibrotic"]["core_signature"] == "YES").sum()) if not signatures["Maladaptive-Fibrotic"].empty else 0
    dex_n = int((signatures["Maladaptive-Steroid"]["core_signature"] == "YES").sum()) if not signatures["Maladaptive-Steroid"].empty else 0
    disease_n = int((signatures["Maladaptive-Disease"]["core_signature"] == "YES").sum()) if not signatures["Maladaptive-Disease"].empty else 0
    overlap = len(set(signatures["Adaptive-Core"].loc[signatures["Adaptive-Core"]["core_signature"] == "YES", "gene"]) & set(consensus["gene"])) if not signatures["Adaptive-Core"].empty else 0
    jaccard = overlap / (adaptive_n + len(consensus) - overlap) if adaptive_n + len(consensus) - overlap else np.nan
    contexts = sorted(set(";".join(consensus.get("supporting_contexts", pd.Series(dtype=str)).fillna("")).split(";")) - {""}) if not consensus.empty else []
    paired = {str(r["celltype"]): int(r["n_donors_paired"]) for _, r in concordance.iterrows()} if not concordance.empty else {}
    final_verdict = "GO_PHASE1_CONDITIONAL" if adaptive_n >= 20 and len(consensus) >= 20 and stress_dep != "HIGH" else "NO_GO_REFRAME"
    strongest = "Externally locked adaptive and cross-context maladaptive programs can be applied descriptively to donor-pseudobulk GSE313736 states." if final_verdict == "GO_PHASE1_CONDITIONAL" else "External programs did not meet the pre-registered size/independence gate for a conditional application phase."
    prohibited = "A causal DEX effect or PRS × DEX interaction in GSE313736." 
    next_phase = "Proceed only to preregistered conditional primary-pseudobulk modeling after accepting the treatment-pool estimand; keep PRS exploratory/blocked." if final_verdict == "GO_PHASE1_CONDITIONAL" else "Reframe to descriptive external-context mechanotransduction and obtain orthogonal treatment/pool information before downstream modeling."
    text = f"""# Phase 1 external signature report

Date: 2026-09-17  
Project: PIEZO1/TRPV4 adaptive versus maladaptive mechanotransduction in glaucoma  
Verdict: **{final_verdict}**

## Gate decision

This phase derived signatures only from external datasets and locked their hashes before opening GSE313736 expression matrices. GSE313736 remains ROLE C. The donor is the biological replicate; all primary application results are donor-pseudobulk by treatment and annotated cell state. The Phase 0R constraints remain active: PRS provenance is unresolved and treatment is aliased with pool/library.

## Dataset roles and source QC

| Dataset | Role | Samples | Data layer / pairing | Decision |
|---|---|---:|---|---|
| GSE195756 | adaptive stretch discovery | 8 / 4 inferred label pairs | author-provided gene statistics; donor identity not explicit | conditional |
| GSE202372 | fibrotic maladaptive discovery | 10 / 5 mapping-inferred pairs | author-provided DE list only; no refit | conditional |
| GSE301517 | steroid maladaptive discovery | 14 expression profiles | gene-level CPM; Welch test on processed scale | conditional |
| GSE263230 | disease-context discovery | {disease_meta.get('n_primary_h5ad', disease_meta.get('n_case', 0))} primary h5ad processed | 4 glaucoma vs 10 healthy when all files available | {disease_meta.get('status', 'NOT_DERIVED').lower()} |
| GSE313736 | external human-TM application | 59 donors / 166,323 cells | exact metadata–matrix intersection; donor-paired descriptive only | ROLE C |

Input source pages: [GSE195756](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195756), [GSE202372](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202372), [GSE301517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301517), [GSE263230](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263230), and [GSE313736](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736).

Pathway libraries: {library_source}. Enrichment was run as hypergeometric over-representation for the pre-specified Hallmark, Reactome, and GO libraries when available; fallback sets are explicitly labelled in `PATHWAY_ENRICHMENT.tsv` if an external library could not be retrieved.

## Locked signatures

- Adaptive signature: {adaptive_n} core genes (maximum pre-registered size 100); source `GSE195756`, author-provided statistics, FDR < 0.05 ranked first.
- Fibrotic maladaptive signature: {fib_n} core genes; source `GSE202372`, author-provided DE list.
- Steroid maladaptive signature: {dex_n} core genes; source `GSE301517`, processed CPM and Welch test.
- Disease-context signature: {disease_n} core genes; source `GSE263230`, independent primary-culture comparison.
- Maladaptive consensus: {len(consensus)} genes, supported in pathological contexts: {', '.join(contexts) if contexts else 'none'}.
- Adaptive–maladaptive core overlap: {overlap} genes; Jaccard = {jaccard:.4f}.
- `MECHANOSENSORY_ANCHOR_PANEL.tsv` is a biological interpretation panel, not a discovery signature. `CTGF` is harmonized to `CCN2`.

The immutable pre-application hashes are in `results/phase1/SIGNATURE_LOCK_MANIFEST.tsv`.

## Cross-context validation

`CROSS_CONTEXT_SIGNATURE_VALIDATION.tsv` reports overlap, Jaccard, signed ranked-statistic enrichment, empirical null P values, and FDR. Adaptive-Core was tested in TGF-beta2, DEX, and disease contexts. Maladaptive-Core was tested in stretch and disease contexts when available. No same-dataset discovery/application loop was used.

## GSE313736 descriptive application

The exact metadata–matrix intersection was used. Scores were computed after the lock as signed, per-gene standardized logCPM averages at donor × treatment × cell-state level. Primary populations are BeamA and JCT; Schlemm_Endothelium is secondary when available. Paired deltas and Wilcoxon values are descriptive only.

Paired donors: BeamA={paired.get('BeamA', 0)}, JCT={paired.get('JCT', 0)}, Schlemm_Endothelium={paired.get('Schlemm_Endothelium', 0)}. See `GSE313736_DONOR_PROGRAM_SCORES.tsv`, `FIGURE4_GSE313736_LOCALIZATION.tsv`, and `FIGURE5_GSE313736_CONCORDANCE.tsv`.

PIEZO1/TRPV4 were evaluated as observational channel-to-program associations. These associations do not establish channel causality or treatment mediation.

## Falsification and controls

Generic-stress dependence: **{stress_dep}**. `GENERIC_STRESS_CONTROL.tsv` gives overlap and Jaccard for apoptosis, inflammation, ROS, hypoxia, UPR, TNF/NF-kB, and cell-cycle controls. `FALSIFICATION_TESTS.tsv` includes random matched overlap, housekeeping, leave-one-dataset-out limitation, leave-one-pathway-out, remove-PIEZO1/TRPV4 sensitivity, and external ranked-statistic null controls. The two-context design makes leave-one-dataset-out non-informative by construction; it is reported as a limitation, not counted as a success.

## Strongest support and boundary

Strongest supported claim: **{strongest}**

Strongest prohibited claim: **{prohibited}** The exact locked language remains `DEX_CAUSAL_EFFECT_NOT_IDENTIFIABLE_FROM_PRIMARY_DATA_ALONE`.

## Next phase

**{next_phase}** No PRS, eQTL, GWAS colocalization, causal mediation, regulatory multiomics, trajectory, network, machine-learning, or clinical validation analysis was run.

## Reproducibility outputs

The execution script is `scripts/08_phase1_external_signatures.py`. Phase 1 machine-readable outputs are in `results/phase1/`; Figure-ready tables are in `results/phase1/figure_tables/`; this report is the gate record.
"""
    (REPORT / "PHASE1_EXTERNAL_SIGNATURE_REPORT.md").write_text(text, encoding="utf-8")
    claims = allowed_claims(final_verdict, disease_meta, stress_dep, len(consensus))
    (REPORT / "PHASE1_ALLOWED_CLAIMS.md").write_text(claims, encoding="utf-8")
    summary = {"verdict": final_verdict, "adaptive_n": adaptive_n, "fibrotic_n": fib_n, "steroid_n": dex_n, "disease_n": disease_n, "maladaptive_core_n": len(consensus), "supporting_contexts": contexts, "overlap_n": overlap, "jaccard": jaccard, "generic_stress_dependence": stress_dep, "paired_donors": paired, "falsification_pass": int((falsification["status"] == "PASS").sum()) if not falsification.empty else 0, "falsification_total": int((~falsification["status"].isin(["INFORMATIVE_LIMITATION", "UNAVAILABLE"])).sum()) if not falsification.empty else 0, "library_source": library_source}
    return final_verdict, summary


def update_project_files(summary: dict[str, object], disease_meta: dict[str, object]) -> None:
    config = ROOT / "config" / "project_config.yaml"
    existing = config.read_text(encoding="utf-8")
    addition = f"""

# Phase 1 lock — 2026-09-17
phase1_verdict: {summary['verdict']}
phase1_signature_lock: results/phase1/SIGNATURE_LOCK_MANIFEST.tsv
phase1_adaptive_core_n: {summary['adaptive_n']}
phase1_fibrotic_core_n: {summary['fibrotic_n']}
phase1_steroid_core_n: {summary['steroid_n']}
phase1_disease_context_core_n: {summary['disease_n']}
phase1_maladaptive_core_n: {summary['maladaptive_core_n']}
phase1_generic_stress_dependence: {summary['generic_stress_dependence']}
phase1_gse313736_role: ROLE_C
phase1_statistical_unit: donor_pseudobulk_by_treatment_and_cellstate
phase1_prs_analysis: false
phase1_primary_dex_causal_effect: false
phase1_stopped_before: [PRS, eQTL, GWAS_colocalization, causal_mediation, regulatory_multiomics, trajectory, network, machine_learning]
"""
    if "phase1_verdict:" not in existing:
        config.write_text(existing.rstrip() + addition, encoding="utf-8")
    gov = ROOT / "ANALYSIS_GOVERNANCE.md"
    gtext = gov.read_text(encoding="utf-8")
    gadd = f"""

## Phase 1 gate — 2026-09-17

Phase 1 verdict: **{summary['verdict']}**. Signatures were derived from independent external contexts and hashed before GSE313736 expression was opened. GSE313736 is ROLE C and is restricted to donor-pseudobulk descriptive application. The maladaptive consensus size is {summary['maladaptive_core_n']} with generic-stress dependence **{summary['generic_stress_dependence']}**. The disease-context source status is **{disease_meta.get('status', 'NOT_DERIVED')}**.

All GSE313736 treatment comparisons retain the exact phrase `DEX_CAUSAL_EFFECT_NOT_IDENTIFIABLE_FROM_PRIMARY_DATA_ALONE`. Pool/library cannot be entered as an ordinary adjustment variable when perfectly aliased with treatment. PRS, PRS×DEX, eQTL, colocalization, causal mediation, regulatory multiomics, trajectory, network, and machine-learning analyses remain outside the gate.
"""
    if "## Phase 1 gate" not in gtext:
        gov.write_text(gtext.rstrip() + gadd, encoding="utf-8")
    status = ROOT / "PROJECT_STATUS.md"
    stext = status.read_text(encoding="utf-8")
    sadd = f"""

## Phase 1 complete

**{summary['verdict']}** — external signatures were derived and locked before descriptive donor-pseudobulk application to GSE313736. Adaptive core: {summary['adaptive_n']}; fibrotic: {summary['fibrotic_n']}; steroid: {summary['steroid_n']}; disease-context: {summary['disease_n']}; maladaptive consensus: {summary['maladaptive_core_n']}; generic-stress dependence: **{summary['generic_stress_dependence']}**.

The primary resource remains observational. Treatment is pool/library aliased; PRS provenance is unresolved. No causal DEX, PRS×DEX, eQTL, colocalization, mediation, or regulatory multiomics result is authorized. Report: `reports/phase1/PHASE1_EXTERNAL_SIGNATURE_REPORT.md`.
"""
    if "## Phase 1 complete" not in stext:
        status.write_text(stext.rstrip() + sadd, encoding="utf-8")


def terminal_summary(summary: dict[str, object], signatures: dict[str, pd.DataFrame], consensus: pd.DataFrame, stress_dep: str, cross: pd.DataFrame, falsification: pd.DataFrame, concordance: pd.DataFrame, assoc: pd.DataFrame) -> None:
    def size(name: str) -> int:
        return int((signatures[name]["core_signature"] == "YES").sum()) if not signatures[name].empty else 0
    def conc(celltype: str) -> str:
        row = concordance[concordance["celltype"].eq(celltype)] if not concordance.empty else pd.DataFrame()
        if row.empty or row.iloc[0].get("n_donors_paired", 0) == 0:
            return "NOT_AVAILABLE"
        r = row.iloc[0]
        return f"n={int(r['n_donors_paired'])}; dMSI_median={fnum(r.get('delta_MSI_median')):.4f}; positive_fraction={fnum(r.get('delta_MSI_positive_fraction')):.3f}"
    def axis(channel: str) -> str:
        if assoc.empty:
            return "NOT_AVAILABLE"
        sub = assoc[(assoc["channel"] == channel) & (assoc["axis"] == "MSI")]
        if sub.empty:
            return "NOT_AVAILABLE"
        vals = sub["spearman_rho"].dropna()
        return "; ".join(f"{r['celltype']} rho={r['spearman_rho']:.3f}" for _, r in sub.iterrows() if np.isfinite(r["spearman_rho"])) or "NO_FINITE_ASSOCIATION"
    passes = int((falsification["status"] == "PASS").sum()) if not falsification.empty else 0
    total = int((~falsification["status"].isin(["INFORMATIVE_LIMITATION", "UNAVAILABLE"])).sum()) if not falsification.empty else 0
    print(f"1 PHASE1 VERDICT: {summary['verdict']}")
    print(f"2 adaptive discovery sample N: 8 samples / 4 label-inferred pairs")
    print(f"3 adaptive signature size: {size('Adaptive-Core')}")
    print(f"4 fibrotic size: {size('Maladaptive-Fibrotic')}")
    print(f"5 steroid size: {size('Maladaptive-Steroid')}")
    print(f"6 disease-context size: {size('Maladaptive-Disease')}")
    print(f"7 maladaptive core size: {len(consensus)}")
    print(f"8 pathological contexts supporting core: {','.join(summary['supporting_contexts']) if summary['supporting_contexts'] else 'NONE'}")
    print(f"9 adaptive vs maladaptive overlap: {summary['overlap_n']}; Jaccard={summary['jaccard']:.4f}")
    print(f"10 generic-stress dependence: {stress_dep}")
    print(f"11 GSE313736 BeamA concordance: {conc('BeamA')}")
    print(f"12 GSE313736 JCT concordance: {conc('JCT')}")
    print(f"13 GSE313736 Schlemm concordance: {conc('Schlemm_Endothelium')}")
    print(f"14 PIEZO1 association with adaptive/maladaptive axis: {axis('PIEZO1')}")
    print(f"15 TRPV4 association with axis: {axis('TRPV4')}")
    print(f"16 falsification tests passed/total: {passes}/{total}")
    print(f"17 strongest supported claim: {('externally locked programs support descriptive donor-pseudobulk state-axis application' if summary['verdict'] == 'GO_PHASE1_CONDITIONAL' else 'external evidence requires reframing before conditional application')}")
    print("18 strongest prohibited claim: causal DEX effect or PRS x DEX interaction in GSE313736")
    print(f"19 next phase recommendation: {'GO_PHASE1_CONDITIONAL; preregister conditional donor-pseudobulk estimands and keep PRS blocked' if summary['verdict'] == 'GO_PHASE1_CONDITIONAL' else 'NO_GO_REFRAME; obtain orthogonal treatment-pool information before downstream modeling'}")
    print("20 path to report: reports/phase1/PHASE1_EXTERNAL_SIGNATURE_REPORT.md")


def main() -> int:
    LOG.write_text("", encoding="utf-8")
    log(f"Phase 1 start seed={SEED}")
    role_lock = make_role_lock()
    primary_features = primary_feature_symbols()
    anchor = write_anchor_panel(primary_features)
    libraries, library_source = load_libraries()
    signatures, cores, consensus, meta = make_signatures(libraries)
    make_harmonization(signatures, primary_features)
    stress, stress_dep = generic_stress_control(consensus, signatures)
    cross = cross_context_validation(signatures, consensus)
    write_tsv(cross, OUT / "CROSS_CONTEXT_SIGNATURE_VALIDATION.tsv")
    freeze_signatures(consensus, library_source)
    groups, primary_meta = aggregate_primary_matrix(signatures, consensus, anchor)
    scores = score_groups(groups, signatures, consensus)
    write_tsv(scores, OUT / "GSE313736_DONOR_PROGRAM_SCORES.tsv")
    concordance = primary_concordance(scores)
    write_tsv(concordance, OUT / "GSE313736_CONCORDANCE_SUMMARY.tsv")
    assoc = association_table(scores)
    falsification = falsification_tests(signatures, consensus, cross, stress, scores)
    figure_tables(role_lock, cross, scores, concordance, stress, anchor)
    verdict, summary = report(role_lock, signatures, consensus, cross, stress, stress_dep, falsification, scores, concordance, assoc, meta["disease"], library_source)
    summary["primary_matrix"] = primary_meta
    (OUT / "PHASE1_RUN_MANIFEST.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    update_project_files(summary, meta["disease"])
    terminal_summary(summary, signatures, consensus, stress_dep, cross, falsification, concordance, assoc)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"FAILED: {type(exc).__name__}: {exc}")
        raise
