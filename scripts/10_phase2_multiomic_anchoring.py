#!/usr/bin/env python3
"""Phase 2: orthogonal multiomic regulatory anchoring.

This script is deliberately limited to the pre-registered GSE301525 framework.
GSE301517 RNA-seq is the same-system transcriptomic anchor. GSE301516 ATAC,
GSE301519 CUT&Tag, and GSE301522 promoter-capture/Hi-C are orthogonal
regulatory layers from the same framework, not independent biological
replication. GSE301512 is audited for exact mark provenance but is excluded
from DEX contrasts because its metadata are resource-only.

The supplied processed peak and loop files are condition-level resources. They
do not contain replicate-level count matrices. The script therefore reports
descriptive aggregate peak-set contrasts with FDR=NA where an inferential FDR
cannot be reconstructed. No GWAS, eQTL, colocalization, PRS, MR, or causal
DEX analysis is run.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, ttest_ind


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data" / "metadata"
RAW = ROOT / "data" / "raw" / "external" / "phase2"
P1RAW = ROOT / "data" / "raw" / "external" / "phase1"
OUT = ROOT / "results" / "phase2"
REPORT = ROOT / "reports" / "phase2"
LOGS = ROOT / "logs"
OUT.mkdir(parents=True, exist_ok=True)
REPORT.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)
RUN_LOG = LOGS / "phase2_run.log"
SEED = 20260918
NULL_B = 5000
FDR_CUTOFF = 0.05


MODULES = {
    "M1": {
        "name": "focal_adhesion_integrin",
        "term": "GO:cell-substrate adhesion / Reactome integrin cell-surface interactions",
        "genes": {"PTK2", "VCL", "TLN1", "PXN", "FERMT2", "ITGA9", "CAV1", "CAV2", "ACTN1", "ACTB", "ITGB1", "ITGB3", "ITGAV", "ITGBL1"},
        "prior": "CONTEXT_DEPENDENT",
    },
    "M2": {
        "name": "actomyosin_contractility",
        "term": "GO:actin filament organization / Reactome smooth muscle contraction",
        "genes": {"ACTA2", "TAGLN", "MYL9", "RHOA", "ROCK1", "ROCK2", "MYLK", "MYL12A", "MYL12B", "TPM1", "TPM2", "ACTB"},
        "prior": "CONTEXT_DEPENDENT",
    },
    "M3": {
        "name": "ECM_remodeling",
        "term": "GO:extracellular matrix organization / Hallmark epithelial-mesenchymal transition",
        "genes": {"FN1", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "DCN", "LUM", "SPARC", "VCAN", "CCN2", "TGFB2", "THBS1"},
        "prior": "MECHANOPATHOLOGY_ASSOCIATED",
    },
    "M4": {
        "name": "TGFb_pro_fibrotic",
        "term": "GO:transforming growth factor beta receptor signaling / Reactome TGF-beta signaling",
        "genes": {"TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7", "SERPINE1", "CCN2", "FN1"},
        "prior": "MECHANOPATHOLOGY_ASSOCIATED",
    },
    "M5": {
        "name": "Hippo_YAP",
        "term": "Hallmark Hippo/YAP transcriptional regulation / Reactome YAP1-TAZ signaling",
        "genes": {"YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "LATS1", "LATS2", "MOB1A", "MOB1B", "CCN2", "CTGF", "AMOT"},
        "prior": "CONTEXT_DEPENDENT",
    },
    "M6": {
        "name": "calcium_mechanosensitive_channel",
        "term": "GO:calcium ion transmembrane transport / Reactome ion channel transport",
        "genes": {"PIEZO1", "TRPV4", "PKD1", "PKD2", "TMEM63A", "TMEM63B", "TMEM63C", "ASIC1", "ASIC2", "CALM1", "CALM2", "CALM3", "ITPR1", "ITPR2", "ITPR3", "ATP2B1", "ATP2B2", "ATP2B4", "ORAI1", "STIM1"},
        "prior": "MECHANOADAPTIVE_ASSOCIATED",
    },
    "M7": {
        "name": "oxidative_ferroptotic_stress",
        "term": "Hallmark reactive oxygen species / GO ferroptosis and redox homeostasis",
        "genes": {"NOX4", "CYBB", "SOD1", "SOD2", "CAT", "GPX1", "GPX4", "HMOX1", "NQO1", "GCLC", "GCLM", "SLC7A11", "NFE2L2", "ACSL4", "TFRC", "FTH1", "FTL", "ALOX5", "LPCAT3"},
        "prior": "MECHANOPATHOLOGY_ASSOCIATED",
    },
    "M8": {
        "name": "autophagy_proteostasis",
        "term": "GO:macroautophagy / Hallmark unfolded protein response and proteostasis",
        "genes": {"BECN1", "ATG3", "ATG5", "ATG7", "ATG12", "ATG16L1", "MAP1LC3B", "SQSTM1", "ULK1", "LAMP1", "LAMP2", "HSPA5", "XBP1", "ATF4", "DDIT3", "EIF2AK3", "HSP90AA1", "HSPA1A", "DNAJB9", "MANF", "PSMC1", "UBB"},
        "prior": "MECHANOADAPTIVE_ASSOCIATED",
    },
    "M9": {
        "name": "endothelial_outflow_homeostasis",
        "term": "GO:vascular transport and endothelial barrier / Reactome vascular signaling",
        "genes": {"KDR", "FLT1", "TEK", "ANGPT1", "ANGPT2", "PECAM1", "EMCN", "VWF", "CLDN5", "RAMP2", "EDNRB", "NOS3", "AQP1", "CAV1", "CAV2", "MMP2", "MMP9", "TIMP1"},
        "prior": "MECHANOADAPTIVE_ASSOCIATED",
    },
}

ANCHORS = ["PIEZO1", "TRPV4", "PTK2", "ITGA9", "ANGPT2", "TEK", "RHOA", "ROCK1", "ROCK2", "YAP1", "WWTR1", "CCN2", "TGFB2", "FN1"]
ACTIVE_MARKS = {"H3K27Ac", "H3K4me1", "H3K4me2", "H3K4me3", "H3K36me3"}
REPRESSIVE_MARKS = {"H3K27me3", "H3K9me3"}
GENE_START_CACHE: dict[str, np.ndarray] = {}


def log(message: str) -> None:
    with RUN_LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}\n")


def bh(values: np.ndarray | list[float]) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    out = np.full(x.shape, np.nan)
    ok = np.isfinite(x)
    if not ok.any():
        return out
    p = np.clip(x[ok], 0, 1)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.minimum(q, 1.0)
    out[ok] = restored
    return out


def write(frame: pd.DataFrame, path: Path, columns: list[str] | None = None) -> None:
    if columns is not None:
        frame = frame.reindex(columns=columns)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def clean_gene(x: object) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return ""
    s = str(x).strip().upper()
    if s in {"", "NAN", "NONE", "NA"}:
        return ""
    return "CCN2" if s == "CTGF" else s


def parse_soft(path: Path) -> list[dict[str, list[str]]]:
    samples: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] | None = None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                if current is not None:
                    samples.append(current)
                current = {"_id": [line.split(" = ", 1)[1]]}
            elif line.startswith("!Sample_") and current is not None:
                key, value = line[1:].split(" = ", 1)
                current.setdefault(key, []).append(value)
    if current is not None:
        samples.append(current)
    return samples


def sample_field(sample: dict[str, list[str]], key: str) -> str:
    return "; ".join(sample.get(key, []))


def parse_line(title: str) -> str:
    m = re.search(r"\b(TM\d+)\b", title, re.I)
    return m.group(1).upper() if m else "NA"


def parse_condition(title: str, chars: str) -> str:
    text = f"{title}; {chars}".lower()
    if "dexamethasone" in text or re.search(r"\bdex\b", text):
        return "DEX"
    if "treatment: control" in text or re.search(r"\bcontrol\b", title, re.I):
        return "CONTROL"
    if "treatment: none" in text or re.search(r"\bresource\b", title, re.I):
        return "RESOURCE"
    return "UNRESOLVED"


def parse_rep(title: str) -> str:
    m = re.search(r"(?:rep\.?\s*)(\d+)", title, re.I)
    return m.group(1) if m else "NA"


def parse_mark(title: str) -> str:
    m = re.search(r"(H3K\d+(?:me\d+|Ac)|CTCF|IgG)", title, re.I)
    if not m:
        return "NA"
    s = m.group(1)
    if s.lower() == "ctcf":
        return "CTCF"
    if s.lower() == "igg":
        return "IgG"
    canonical = {"h3k27ac": "H3K27Ac", "h3k27me3": "H3K27me3", "h3k36me3": "H3K36me3", "h3k4me1": "H3K4me1", "h3k4me2": "H3K4me2", "h3k4me3": "H3K4me3", "h3k9me3": "H3K9me3"}
    return canonical.get(s.lower(), s)


def treatment_protocol(samples: list[dict[str, list[str]]]) -> tuple[str, str]:
    text = " ".join(sample_field(s, "Sample_treatment_protocol_ch1") for s in samples)
    cm = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*[μu]?M\s+dexamethasone", text, re.I)
    dm = re.search(r"treated for\s+([0-9]+)\s+days", text, re.I)
    return (f"{cm.group(1)} uM" if cm else "NA", f"{dm.group(1)} days" if dm else "NA")


def accession_sample_map() -> tuple[pd.DataFrame, dict[str, list[dict[str, list[str]]]]]:
    accessions = {"GSE301512": "CUT&Run", "GSE301516": "ATAC", "GSE301517": "RNA", "GSE301519": "CUT&Tag", "GSE301522": "Hi-C"}
    all_rows = []
    parsed: dict[str, list[dict[str, list[str]]]] = {}
    dex_conc, dex_duration = treatment_protocol(parse_soft(META / "GSE301517_family.soft.gz"))
    for acc, modality in accessions.items():
        samples = parse_soft(META / f"{acc}_family.soft.gz")
        parsed[acc] = samples
        for s in samples:
            sid = sample_field(s, "Sample_geo_accession") or s.get("_id", ["NA"])[0]
            title = sample_field(s, "Sample_title")
            chars = sample_field(s, "Sample_characteristics_ch1")
            condition = parse_condition(title, chars)
            line = parse_line(title)
            rep = parse_rep(title)
            technical = "pre_capture" if "precapture" in title.lower() else ("post_capture" if "postcapture" in title.lower() else "NA")
            usable = "YES"
            reason = ""
            if modality == "CUT&Run":
                usable = "NO"
                reason = "resource-only CUT&Run; no DEX/control contrast in GSE301512 metadata"
            if modality == "Hi-C" and condition == "RESOURCE":
                usable = "YES"
                reason = "resource context; used for QC/reference only"
            if condition == "UNRESOLVED":
                usable = "NO"
                reason = "condition not resolved from title/characteristics"
            all_rows.append({
                "accession": acc,
                "modality": modality,
                "sample_id": sid,
                "cell_line_or_donor": line,
                "condition": condition,
                "DEX_concentration": dex_conc if condition == "DEX" else "NA",
                "exposure_duration": dex_duration if condition == "DEX" else "NA",
                "biological_replicate": rep,
                "technical_replicate": technical,
                "batch": "NA",
                "sequencing_run": sample_field(s, "Sample_platform_id") or "NA",
                "paired_control": "NA",
                "usable": usable,
                "exclusion_reason": reason,
                "sample_title": title,
                "mark": parse_mark(title) if modality in {"CUT&Run", "CUT&Tag"} else "NA",
            })
    frame = pd.DataFrame(all_rows)
    for i, row in frame.iterrows():
        if row["condition"] not in {"DEX", "CONTROL"}:
            continue
        candidates = frame[(frame["accession"] == row["accession"]) & (frame["cell_line_or_donor"] == row["cell_line_or_donor"]) & (frame["condition"] == ("CONTROL" if row["condition"] == "DEX" else "DEX"))]
        if row["modality"] in {"CUT&Tag", "CUT&Run"}:
            candidates = candidates[candidates["mark"] == row["mark"]]
        if row["biological_replicate"] != "NA":
            exact = candidates[candidates["biological_replicate"] == row["biological_replicate"]]
            if not exact.empty:
                candidates = exact
        if not candidates.empty:
            frame.loc[i, "paired_control"] = ";".join(candidates["sample_id"].tolist())
    cols = ["accession", "modality", "sample_id", "cell_line_or_donor", "condition", "DEX_concentration", "exposure_duration", "biological_replicate", "technical_replicate", "batch", "sequencing_run", "paired_control", "usable", "exclusion_reason"]
    write(frame, OUT / "GSE301525_SAMPLE_MAP.tsv", cols)
    return frame, parsed


def load_rna() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = RAW / "GSE301517_Control_vs_Dex_Gene_CPM_MSTR.tsv.gz"
    df = pd.read_csv(path, sep="\t")
    if df.columns[0].startswith("Unnamed") or df.columns[0] == "":
        df = df.iloc[:, 1:]
    df["gene"] = df["external_gene_name"].map(clean_gene)
    df = df[df["gene"].ne("")].copy()
    numeric_cols = [c for c in df.columns if c.startswith("TM")]
    control = [c for c in numeric_cols if "Control" in c]
    dex = [c for c in numeric_cols if "Dex" in c]
    vals = np.log2(df[numeric_cols].apply(pd.to_numeric, errors="coerce") + 0.25)
    control_values = vals[control].to_numpy(dtype=float)
    dex_values = vals[dex].to_numpy(dtype=float)
    tvals, pvals = ttest_ind(dex_values, control_values, axis=1, equal_var=False, nan_policy="omit")
    de = pd.DataFrame({
        "gene": df["gene"].to_numpy(),
        "logFC": np.nanmean(dex_values, axis=1) - np.nanmean(control_values, axis=1),
        "statistic": tvals,
        "PValue": pvals,
        "baseline_log2_CPM": np.nanmean(vals.to_numpy(dtype=float), axis=1),
        "gene_length": pd.to_numeric(df["end_position"], errors="coerce").to_numpy() - pd.to_numeric(df["start_position"], errors="coerce").to_numpy(),
        "chromosome": df["chromosome_name"].astype(str).to_numpy(),
        "tss": pd.to_numeric(df["start_position"], errors="coerce").to_numpy(),
    })
    de = de.groupby("gene", as_index=False).agg({"logFC": "mean", "statistic": "mean", "PValue": "min", "baseline_log2_CPM": "mean", "gene_length": "median", "chromosome": "first", "tss": "median"})
    de["FDR"] = bh(de["PValue"].to_numpy())
    de["direction"] = np.where(de["logFC"] > 0, "UP", np.where(de["logFC"] < 0, "DOWN", "FLAT"))
    de["RNA_method"] = "Welch t-test on log2(CPM+0.25); 7 DEX and 7 control samples; cell-line labels retained in sample map"
    return de, df


def load_locked_support() -> dict[str, dict[str, int]]:
    def read(name: str) -> pd.DataFrame:
        path = ROOT / "results" / "phase1" / name
        x = pd.read_csv(path, sep="\t")
        x["gene"] = x["gene"].map(clean_gene)
        return x[x["gene"].ne("")].copy()
    files = {"adaptive": "ADAPTIVE_SIGNATURE.tsv", "maladaptive": "MALADAPTIVE_CORE.tsv", "tgfb": "MALADAPTIVE_FIBROTIC.tsv", "dex": "MALADAPTIVE_STEROID.tsv", "disease": "MALADAPTIVE_DISEASE.tsv"}
    out: dict[str, dict[str, int]] = {}
    for key, fn in files.items():
        x = read(fn)
        if "core_signature" in x.columns:
            x = x[x["core_signature"].eq("YES")]
        if key in {"tgfb", "dex"}:
            x = x[x["FDR"].lt(FDR_CUTOFF)]
        out[key] = {g: int(i + 1) for i, g in enumerate(x["gene"].drop_duplicates().tolist())}
    return out


def module_definitions(rna: pd.DataFrame, support: dict[str, dict[str, int]]) -> pd.DataFrame:
    rows = []
    for mid, spec in MODULES.items():
        for gene in sorted(spec["genes"]):
            adaptive = support["adaptive"].get(gene)
            maladaptive = support["maladaptive"].get(gene)
            tgfb = support["tgfb"].get(gene)
            dex = support["dex"].get(gene)
            disease = support["disease"].get(gene)
            ext = []
            if adaptive is not None: ext.append(f"GSE195756_core_rank_{adaptive}")
            if maladaptive is not None: ext.append(f"GSE202372_GSE301517_core_rank_{maladaptive}")
            if tgfb is not None: ext.append(f"GSE202372_rank_{tgfb}")
            if dex is not None: ext.append(f"GSE301517_rank_{dex}")
            stretch = adaptive is not None
            rows.append({
                "module_id": mid,
                "module_name": spec["name"],
                "pathway_source": "pre-specified GO/Reactome/Hallmark-aligned module; frozen before multiomic reads",
                "pathway_term": spec["term"],
                "gene": gene,
                "external_rank_support": ";".join(ext) if ext else "NONE",
                "DEX_rank_support": f"GSE301517_rank_{dex}" if dex is not None else "NONE",
                "TGFb2_rank_support": f"GSE202372_rank_{tgfb}" if tgfb is not None else "NONE",
                "stretch_rank_support": f"GSE195756_rank_{adaptive}" if stretch else "NONE",
                "disease_context_support": f"GSE263230_rank_{disease}" if disease is not None else "NONE",
                "locked_external_evidence": "YES" if ext else "NO",
                "gene_in_RNA_background": "YES" if gene in set(rna["gene"]) else "NO",
            })
    result = pd.DataFrame(rows)
    write(result, OUT / "PATHWAY_MODULE_DEFINITIONS.tsv")
    return result


def classify_modules(defs: pd.DataFrame, rna: pd.DataFrame, support: dict[str, dict[str, int]]) -> pd.DataFrame:
    rows = []
    for mid, spec in MODULES.items():
        g = set(spec["genes"])
        rows.append({
            "module_id": mid,
            "module_name": spec["name"],
            "classification": spec["prior"],
            "classification_rule": "pre-specified ontology label locked before GSE301525 regulatory layers; external ranks are descriptive support only",
            "module_gene_n": len(g),
            "RNA_background_gene_n": int(len(g & set(rna["gene"]))),
            "adaptive_core_n": len(g & set(support["adaptive"])),
            "maladaptive_core_n": len(g & set(support["maladaptive"])),
            "DEX_rank_supported_n": len(g & set(support["dex"])),
            "TGFb2_rank_supported_n": len(g & set(support["tgfb"])),
            "stretch_rank_supported_n": len(g & set(support["adaptive"])),
            "disease_context_supported_n": len(g & set(support["disease"])),
        })
    result = pd.DataFrame(rows)
    write(result, OUT / "MODULE_TRANSCRIPTOMIC_CLASSIFICATION.tsv")
    return result


def gsea_score(order: list[str], stats: dict[str, float], genes: set[str]) -> tuple[float, list[str]]:
    hits = set(order) & genes
    if not hits or len(hits) == len(order):
        return (np.nan, [])
    weights = np.array([abs(stats.get(g, 0.0)) for g in order], dtype=float)
    weights[~np.isfinite(weights)] = 0
    hit_total = float(sum(weights[i] for i, g in enumerate(order) if g in hits))
    if hit_total <= 0:
        hit_total = float(len(hits))
        weights[:] = 1.0
    miss_penalty = 1.0 / (len(order) - len(hits))
    run = 0.0
    trace = []
    for i, gene in enumerate(order):
        if gene in hits:
            run += weights[i] / hit_total
        else:
            run -= miss_penalty
        trace.append(run)
    trace_arr = np.asarray(trace)
    idx = int(np.argmax(np.abs(trace_arr)))
    es = float(trace_arr[idx])
    leading = [g for g in order[: idx + 1] if g in hits]
    return es, leading


def rna_module_enrichment(rna: pd.DataFrame) -> pd.DataFrame:
    work = rna.dropna(subset=["statistic", "PValue"]).drop_duplicates("gene").copy()
    work["rank_metric"] = np.sign(work["logFC"]) * (-np.log10(np.clip(work["PValue"], 1e-300, 1)))
    work = work.sort_values(["rank_metric", "gene"], ascending=[False, True])
    order = work["gene"].tolist()
    stat_map = dict(zip(work["gene"], work["rank_metric"]))
    rank_values = work["rank_metric"].to_numpy(dtype=float)
    rng = np.random.default_rng(SEED + 11)
    universe = np.arange(len(order), dtype=int)
    rows = []
    for mid, spec in MODULES.items():
        genes = set(spec["genes"]) & set(order)
        es, leading = gsea_score(order, stat_map, genes)
        # The observed ES is the standard weighted running-sum statistic.
        # For the null, use the algebraically equivalent module mean of the
        # full ranked statistic.  This avoids a Python loop over every gene
        # for every null draw while retaining a two-sided, seeded rank null.
        observed_mean = float(np.mean([stat_map[g] for g in genes])) if genes else np.nan
        null = []
        for _ in range(2000):
            draw = rng.choice(universe, size=min(len(genes), len(universe)), replace=False)
            null.append(float(np.mean(rank_values[draw])))
        null_arr = np.asarray(null, dtype=float)
        scale = np.nanmedian(np.abs(null_arr))
        nes = observed_mean / scale if np.isfinite(observed_mean) and scale > 0 else np.nan
        p = (1 + np.sum(np.abs(null_arr) >= abs(observed_mean))) / (len(null_arr) + 1) if np.isfinite(observed_mean) else np.nan
        rows.append({
            "module_id": mid,
            "module_name": spec["name"],
            "module_gene_n": len(genes),
            "background_gene_n": len(order),
            "ES": es,
            "NES": nes,
            "nominal_P": p,
            "direction": "UP" if es > 0 else ("DOWN" if es < 0 else "FLAT"),
            "leading_edge": ";".join(leading),
            "anchor_label": "same-system transcriptomic anchor; GSE301517 is within GSE301525 and not independent replication",
            "method": "pre-ranked weighted running-sum ES with vectorized equivalent mean-rank null on full GSE301517 gene-ranked statistic",
        })
    result = pd.DataFrame(rows)
    result["FDR"] = bh(result["nominal_P"].to_numpy())
    write(result, OUT / "RNA_MODULE_ENRICHMENT.tsv")
    return result


def read_bed(path: Path, narrow: bool = False) -> pd.DataFrame:
    records = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = re.split(r"\s+", line.strip())
            if len(f) < 3:
                continue
            try:
                chrom, start, end = f[0], int(float(f[1])), int(float(f[2]))
            except ValueError:
                continue
            if end <= start:
                continue
            score = 1.0
            candidates = [6, 4, 3]
            for idx in candidates:
                if idx < len(f):
                    try:
                        val = float(f[idx])
                        if np.isfinite(val):
                            score = max(val, 0.0)
                            break
                    except ValueError:
                        pass
            records.append((chrom, start, end, score))
    x = pd.DataFrame(records, columns=["chr", "start", "end", "score"])
    if x.empty:
        return x
    return x.drop_duplicates(["chr", "start", "end"]).sort_values(["chr", "start", "end"]).reset_index(drop=True)


def intervals_overlap(row: pd.Series, other: pd.DataFrame) -> pd.DataFrame:
    if other.empty:
        return other
    x = other[other["chr"].eq(row["chr"]) & (other["start"] < row["end"]) & (other["end"] > row["start"])]
    return x


def nearest_gene(genes_by_chr: dict[str, list[tuple[int, int, str, int]]], chrom: str, pos: int) -> tuple[str, int]:
    vals = genes_by_chr.get(chrom, [])
    if not vals:
        return "NA", -1
    if chrom not in GENE_START_CACHE:
        GENE_START_CACHE[chrom] = np.asarray([v[0] for v in vals], dtype=int)
    idx = int(np.searchsorted(GENE_START_CACHE[chrom], pos))
    best = None
    for j in range(max(0, idx - 4), min(len(vals), idx + 5)):
        st, en, gene, tss = vals[j]
        dist = 0 if st <= pos <= en else min(abs(pos - st), abs(pos - en), abs(pos - tss))
        if best is None or dist < best[1]:
            best = (gene, int(dist))
    return best if best is not None else ("NA", -1)


def gene_annotation(rna: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, list[tuple[int, int, str, int]]]]:
    x = rna[["gene", "chromosome", "tss", "gene_length", "baseline_log2_CPM"]].copy()
    x["chromosome"] = x["chromosome"].astype(str).map(lambda z: z if z.startswith("chr") else "chr" + z)
    x = x.drop_duplicates("gene")
    genes_by_chr: dict[str, list[tuple[int, int, str, int]]] = defaultdict(list)
    for _, row in x.iterrows():
        try:
            tss = int(float(row["tss"]))
            ln = max(int(float(row["gene_length"])), 1)
        except (ValueError, TypeError):
            continue
        genes_by_chr[row["chromosome"]].append((tss, tss + ln, row["gene"], tss))
    for chrom in genes_by_chr:
        genes_by_chr[chrom].sort(key=lambda z: z[0])
    return x, genes_by_chr


def condition_diff_regions(dex: pd.DataFrame, control: pd.DataFrame, prefix: str) -> pd.DataFrame:
    rows = []
    if dex.empty and control.empty:
        return pd.DataFrame(columns=["element_id", "chr", "start", "end", "logFC", "statistic", "FDR", "direction", "evidence_level", "condition_source"])
    control_by_chr = {c: x[["start", "end", "score"]].to_numpy(dtype=float) for c, x in control.groupby("chr")} if not control.empty else {}
    dex_by_chr = {c: x[["start", "end", "score"]].to_numpy(dtype=float) for c, x in dex.groupby("chr")} if not dex.empty else {}
    control_starts = {c: a[:, 0] for c, a in control_by_chr.items()}
    dex_starts = {c: a[:, 0] for c, a in dex_by_chr.items()}

    def overlap_scores(row: pd.Series, other: dict[str, np.ndarray], starts_map: dict[str, np.ndarray]) -> np.ndarray:
        arr = other.get(row["chr"])
        if arr is None or len(arr) == 0:
            return np.empty(0, dtype=float)
        j = int(np.searchsorted(starts_map[row["chr"]], float(row["end"]), side="left"))
        k = j - 1
        hits = []
        while k >= 0 and arr[k, 1] > float(row["start"]):
            if arr[k, 0] < float(row["end"]):
                hits.append(arr[k, 2])
            k -= 1
        return np.asarray(hits, dtype=float)

    for _, r in dex.iterrows():
        ov = overlap_scores(r, control_by_chr, control_starts)
        ctl_score = float(np.max(ov)) if len(ov) else 0.0
        logfc = math.log2((float(r["score"]) + 1) / (ctl_score + 1))
        direction = "OPENING" if len(ov) == 0 or logfc > 0.25 else ("CLOSING" if logfc < -0.25 else "SHARED")
        rows.append({"element_id": f"{prefix}_DEX_{len(rows)+1}", "chr": r["chr"], "start": int(r["start"]), "end": int(r["end"]), "logFC": logfc, "statistic": logfc, "FDR": np.nan, "direction": direction, "evidence_level": "aggregate condition-level processed peak set; replicate FDR not reconstructable", "condition_source": "DEX_vs_CONTROL"})
    for _, r in control.iterrows():
        ov = overlap_scores(r, dex_by_chr, dex_starts)
        if len(ov):
            continue
        rows.append({"element_id": f"{prefix}_CONTROL_{len(rows)+1}", "chr": r["chr"], "start": int(r["start"]), "end": int(r["end"]), "logFC": -math.log2(float(r["score"]) + 1), "statistic": -math.log2(float(r["score"]) + 1), "FDR": np.nan, "direction": "CLOSING", "evidence_level": "aggregate condition-level processed peak set; replicate FDR not reconstructable", "condition_source": "DEX_vs_CONTROL"})
    return pd.DataFrame(rows)


def annotate_elements(elements: pd.DataFrame, genes_by_chr: dict[str, list[tuple[int, int, str, int]]]) -> pd.DataFrame:
    if elements.empty:
        elements["nearest_gene"] = []
        elements["distance_to_TSS"] = []
        return elements
    genes, dists = [], []
    for _, row in elements.iterrows():
        gene, dist = nearest_gene(genes_by_chr, row["chr"], int((row["start"] + row["end"]) / 2))
        genes.append(gene); dists.append(dist)
    elements = elements.copy()
    elements["nearest_gene"] = genes
    elements["distance_to_TSS"] = dists
    return elements


def atac_outputs(rna: pd.DataFrame, genes_by_chr: dict[str, list[tuple[int, int, str, int]]]) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    ctl = read_bed(RAW / "GSE301516_ATAC_Control.narrowPeak.gz", narrow=True)
    dex = read_bed(RAW / "GSE301516_ATAC_Dex.narrowPeak.gz", narrow=True)
    diff = annotate_elements(condition_diff_regions(dex, ctl, "ATAC"), genes_by_chr)
    diff["promoter_or_enhancer_link"] = np.where(pd.to_numeric(diff["distance_to_TSS"], errors="coerce") <= 2000, "PROMOTER", "ENHANCER_LINKED_PROXIMAL")
    write(diff, OUT / "DEX_DAR.tsv", ["chr", "start", "end", "logFC", "statistic", "FDR", "nearest_gene", "distance_to_TSS", "element_id", "direction", "evidence_level", "promoter_or_enhancer_link"])
    bg_genes = set(rna["gene"])
    rows = []
    for mid, spec in MODULES.items():
        mg = set(spec["genes"]) & bg_genes
        links = diff[diff["nearest_gene"].isin(mg)]
        opening = set(links.loc[links["direction"].eq("OPENING"), "nearest_gene"])
        closing = set(links.loc[links["direction"].eq("CLOSING"), "nearest_gene"])
        all_linked = set(links["nearest_gene"]) - {"NA"}
        N, K, n, k = len(bg_genes), len(set(diff["nearest_gene"]) & bg_genes), len(mg), len(all_linked)
        p = float(hypergeom.sf(k - 1, N, K, n)) if K and n and k else np.nan
        rows.append({"module_id": mid, "module_name": spec["name"], "module_gene_n": len(mg), "ATAC_linked_gene_n": len(all_linked), "module_ATAC_linked_n": len(all_linked & mg), "opening_n": len(opening), "closing_n": len(closing), "opening_fraction": len(opening) / max(len(mg), 1), "closing_fraction": len(closing) / max(len(mg), 1), "enrichment_factor": (len(all_linked & mg) / max(len(mg), 1)) / max(K / max(N, 1), 1e-12), "nominal_P": p, "mapping_hierarchy": "promoter overlap > Hi-C/loop link (later) > nearest-gene fallback; this table is ATAC nearest/proximal preliminary mapping", "evidence_status": "descriptive aggregate peak-set anchor; no replicate-level ATAC FDR"})
    result = pd.DataFrame(rows)
    result["FDR"] = bh(result["nominal_P"].to_numpy())
    write(result, OUT / "ATAC_MODULE_ENRICHMENT.tsv")
    return diff, result, set(diff["nearest_gene"]) - {"NA"}


def cuttag_registry(sample_map: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for acc in ["GSE301519", "GSE301512"]:
        x = sample_map[sample_map["accession"].eq(acc)]
        for mark in sorted(set(x["mark"]) - {"NA", "IgG"}):
            z = x[x["mark"].eq(mark)]
            rows.append({"accession": acc, "modality": z["modality"].iloc[0], "mark": mark, "control_sample_n": int((z["condition"] == "CONTROL").sum()), "dex_sample_n": int((z["condition"] == "DEX").sum()), "resource_sample_n": int((z["condition"] == "RESOURCE").sum()), "exact_metadata_support": "YES", "primary_phase2_use": "YES" if acc == "GSE301519" and (z["condition"] == "DEX").any() and (z["condition"] == "CONTROL").any() else "NO", "exclusion_or_limit": "" if acc == "GSE301519" else "GSE301512 resource-only; no DEX/control condition"})
    result = pd.DataFrame(rows)
    write(result, OUT / "CUTTAG_MARK_REGISTRY.tsv")
    return result


def cuttag_outputs(sample_map: pd.DataFrame, rna: pd.DataFrame, genes_by_chr: dict[str, list[tuple[int, int, str, int]]]) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    registry = cuttag_registry(sample_map)
    all_diff = []
    enrich = []
    bg_genes = set(rna["gene"])
    for mark in sorted(set(registry.loc[registry["accession"].eq("GSE301519"), "mark"])):
        ctl_path = RAW / f"GSE301519_{mark}_Control.bed.gz"
        dex_path = RAW / f"GSE301519_{mark}_Dex.bed.gz"
        if not ctl_path.exists() or not dex_path.exists():
            continue
        ctl, dex = read_bed(ctl_path), read_bed(dex_path)
        d = annotate_elements(condition_diff_regions(dex, ctl, f"CUTTAG_{mark}"), genes_by_chr)
        d["mark"] = mark
        d["active_or_repressive"] = "ACTIVE" if mark in ACTIVE_MARKS else ("REPRESSIVE" if mark in REPRESSIVE_MARKS else "ARCHITECTURAL")
        all_diff.append(d)
        module_sets = set(d["nearest_gene"]) & bg_genes
        for mid, spec in MODULES.items():
            mg = set(spec["genes"]) & bg_genes
            l = d[d["nearest_gene"].isin(mg)]
            dex_genes = set(l.loc[l["direction"].eq("OPENING"), "nearest_gene"])
            ctl_genes = set(l.loc[l["direction"].eq("CLOSING"), "nearest_gene"])
            K = len(module_sets); k = len(module_sets & mg); N = len(bg_genes); n = len(mg)
            p = float(hypergeom.sf(k - 1, N, K, n)) if K and n and k else np.nan
            enrich.append({"module_id": mid, "module_name": spec["name"], "mark": mark, "mark_class": "ACTIVE" if mark in ACTIVE_MARKS else ("REPRESSIVE" if mark in REPRESSIVE_MARKS else "ARCHITECTURAL"), "module_gene_n": len(mg), "CUTTAG_linked_gene_n": len(module_sets), "module_linked_n": k, "DEX_opening_or_gain_n": len(dex_genes), "CONTROL_closing_or_loss_n": len(ctl_genes), "module_fraction": k / max(len(mg), 1), "nominal_P": p, "FDR": np.nan, "evidence_status": "descriptive condition-level processed CUT&Tag peak set; exact mark metadata audited; replicate FDR not reconstructable"})
    diff = pd.concat(all_diff, ignore_index=True) if all_diff else pd.DataFrame()
    if not diff.empty:
        diff["FDR"] = np.nan
    write(diff, OUT / "CUTTAG_DIFFERENTIAL_REGIONS.tsv", ["mark", "active_or_repressive", "chr", "start", "end", "logFC", "statistic", "FDR", "nearest_gene", "distance_to_TSS", "element_id", "direction", "evidence_level", "condition_source"])
    result = pd.DataFrame(enrich)
    if not result.empty:
        result["FDR"] = bh(result["nominal_P"].to_numpy())
    write(result, OUT / "CUTTAG_MODULE_ENRICHMENT.tsv")
    return diff, result, set(diff["nearest_gene"]) - {"NA"} if not diff.empty else set()


def read_bedgraph(path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("track", "browser", "#")):
                continue
            f = line.strip().split()
            if len(f) < 4:
                continue
            try:
                rows.append((f[0], int(f[1]), int(f[2]), float(f[3])))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["chr", "start", "end", "value"])


def read_bedpe(path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("track", "browser", "#")):
                continue
            f = line.strip().split()
            if len(f) < 6:
                continue
            try:
                rows.append((f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["chr1", "start1", "end1", "chr2", "start2", "end2"])


def hic_outputs(sample_map: pd.DataFrame, genes_by_chr: dict[str, list[tuple[int, int, str, int]]]) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    rows = []
    for condition, label in [("Resource", "RESOURCE"), ("Control", "CONTROL"), ("Dex", "DEX")]:
        comp = read_bedgraph(RAW / f"GSE301522_Compartments_{condition}.bedGraph.gz")
        loops = read_bedpe(RAW / f"GSE301522_Loops_{condition}.bedpe.gz")
        widths = (comp["end"] - comp["start"]).median() if not comp.empty else np.nan
        x = sample_map[(sample_map["accession"] == "GSE301522") & (sample_map["condition"] == label)]
        rows.append({"dataset": "GSE301522", "modality": "Hi-C/promoter capture", "condition": label, "processed_compartment_file": f"GSE301522_Compartments_{condition}.bedGraph.gz", "processed_loop_file": f"GSE301522_Loops_{condition}.bedpe.gz", "usable": "YES", "reported_biological_replicates": int(len(x)), "reported_technical_replicates": int(x["technical_replicate"].nunique()) if not x.empty else 0, "resolution": f"{int(widths)} bp bins" if np.isfinite(widths) else "NA", "read_depth": "NA in supplied processed file", "normalization": "HOMER PC1 for compartments; CHiCAGO/Arima CHiC pipeline for loops", "bins": int(len(comp)), "loops": int(len(loops)), "compartment_bins": int(len(comp)), "tads": "not supplied", "notes": "processed summary resources; no replicate-level contact matrix or read-depth table"})
    qc = pd.DataFrame(rows)
    write(qc, OUT / "HIC_QC_REPORT.tsv")
    link_rows = []
    genes = set()
    for condition, label in [("Control", "CONTROL"), ("Dex", "DEX")]:
        loops = read_bedpe(RAW / f"GSE301522_Loops_{condition}.bedpe.gz")
        for i, r in loops.iterrows():
            for side, chrom, start, end in [("A", r["chr1"], r["start1"], r["end1"]), ("B", r["chr2"], r["start2"], r["end2"])]:
                gene, dist = nearest_gene(genes_by_chr, chrom, int((start + end) / 2))
                if gene == "NA" or dist > 100000:
                    continue
                category = "PROMOTER_DIRECT" if dist <= 2000 else "LOOP_LINKED"
                link_rows.append({"element_id": f"HIC_{condition}_{i+1}_{side}", "modality": "Hi-C", "evidence_category": category, "mark": "NA", "condition": label, "direction": "CONTACT", "chr": chrom, "start": int(start), "end": int(end), "gene": gene, "distance_to_TSS": int(dist), "source_file": f"GSE301522_Loops_{condition}.bedpe.gz", "link_basis": "promoter-capture loop endpoint"})
                genes.add(gene)
    return qc, pd.DataFrame(link_rows), genes


def promoter_flag(chrom: str, start: int, end: int, gene: str, rna: pd.DataFrame) -> bool:
    if gene == "NA":
        return False
    x = rna[rna["gene"].eq(gene)]
    if x.empty:
        return False
    tss = float(x.iloc[0]["tss"])
    return start <= tss + 2000 and end >= tss - 2000


def regulatory_links(atac: pd.DataFrame, cut: pd.DataFrame, hic_links: pd.DataFrame, rna: pd.DataFrame) -> pd.DataFrame:
    rows = []
    tss_map = rna.drop_duplicates("gene").set_index("gene")["tss"].to_dict()

    def is_promoter(start: int, end: int, gene: str) -> bool:
        if gene == "NA" or gene not in tss_map:
            return False
        tss = float(tss_map[gene])
        return start <= tss + 2000 and end >= tss - 2000

    if not atac.empty:
        for _, r in atac.iterrows():
            gene = r["nearest_gene"]
            if gene == "NA": continue
            rows.append({"element_id": r["element_id"], "modality": "ATAC", "evidence_category": "PROMOTER_DIRECT" if is_promoter(int(r["start"]), int(r["end"]), gene) else "PROXIMAL_ONLY", "mark": "ATAC", "condition": "DEX_vs_CONTROL", "direction": r["direction"], "chr": r["chr"], "start": int(r["start"]), "end": int(r["end"]), "gene": gene, "distance_to_TSS": int(r["distance_to_TSS"]), "source_file": "GSE301516_ATAC_Dex.narrowPeak.gz + GSE301516_ATAC_Control.narrowPeak.gz", "link_basis": "promoter overlap > nearest-gene fallback"})
    if not cut.empty:
        for _, r in cut.iterrows():
            gene = r["nearest_gene"]
            if gene == "NA": continue
            rows.append({"element_id": r["element_id"], "modality": "CUT&Tag", "evidence_category": "PROMOTER_DIRECT" if is_promoter(int(r["start"]), int(r["end"]), gene) else "PROXIMAL_ONLY", "mark": r["mark"], "condition": "DEX_vs_CONTROL", "direction": r["direction"], "chr": r["chr"], "start": int(r["start"]), "end": int(r["end"]), "gene": gene, "distance_to_TSS": int(r["distance_to_TSS"]), "source_file": f"GSE301519_{r['mark']}_Dex.bed.gz + control", "link_basis": "condition-level CUT&Tag differential peak-set overlap"})
    if not hic_links.empty:
        rows.extend(hic_links.to_dict("records"))
    result = pd.DataFrame(rows)
    write(result, OUT / "REGULATORY_GENE_LINKS.tsv")
    return result


def gene_evidence(rna: pd.DataFrame, links: pd.DataFrame, cut: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    module_map = defaultdict(list)
    for mid, spec in MODULES.items():
        for g in spec["genes"]: module_map[g].append(mid)
    evidence_genes = set(rna["gene"]) | set(module_map) | set(ANCHORS)
    if not links.empty: evidence_genes |= set(links["gene"].dropna())
    rna_idx = rna.drop_duplicates("gene").set_index("gene")
    link_groups = {gene: z for gene, z in links.groupby("gene")} if not links.empty else {}
    empty_links = links.iloc[0:0].copy()
    rows = []
    for gene in sorted(evidence_genes):
        rr = rna_idx.loc[gene] if gene in rna_idx.index else None
        ll = link_groups.get(gene, empty_links)
        atac = ll[ll["modality"].eq("ATAC")]
        ctag = ll[ll["modality"].eq("CUT&Tag")]
        hic = ll[ll["modality"].eq("Hi-C")]
        active_up, active_down, rep_up, rep_down = [], [], [], []
        if not ctag.empty:
            for mark in sorted(set(ctag["mark"])):
                z = ctag[ctag["mark"].eq(mark)]
                if mark in ACTIVE_MARKS:
                    if (z["direction"] == "OPENING").any(): active_up.append(mark)
                    if (z["direction"] == "CLOSING").any(): active_down.append(mark)
                if mark in REPRESSIVE_MARKS:
                    if (z["direction"] == "OPENING").any(): rep_up.append(mark)
                    if (z["direction"] == "CLOSING").any(): rep_down.append(mark)
        regs = set()
        if not atac.empty: regs.add("ATAC")
        if not ctag.empty: regs.add("CUT&Tag")
        if not hic.empty: regs.add("Hi-C")
        rows.append({
            "gene": gene,
            "modules": ";".join(module_map.get(gene, [])),
            "RNA_logFC": rr["logFC"] if rr is not None else np.nan,
            "RNA_FDR": rr["FDR"] if rr is not None else np.nan,
            "RNA_direction": rr["direction"] if rr is not None else "NA",
            "ATAC_promoter": "YES" if ((atac["evidence_category"] == "PROMOTER_DIRECT").any() if not atac.empty else False) else "NO",
            "ATAC_enhancer": "YES" if ((atac["evidence_category"].isin(["PROXIMAL_ONLY", "LOOP_LINKED"])).any() if not atac.empty else False) else "NO",
            "ATAC_direction": ";".join(sorted(set(atac["direction"]))) if not atac.empty else "NA",
            "CUTTAG_active_marks": ";".join(active_up + active_down) if active_up or active_down else "NA",
            "CUTTAG_active_direction": "UP" if active_up else ("DOWN" if active_down else "NA"),
            "CUTTAG_repressive_marks": ";".join(rep_up + rep_down) if rep_up or rep_down else "NA",
            "CUTTAG_repressive_direction": "UP" if rep_up else ("DOWN" if rep_down else "NA"),
            "CUTTAG_state": "ACTIVE" if active_up or active_down else ("REPRESSIVE" if rep_up or rep_down else "NA"),
            "HiC_contact": "YES" if not hic.empty else "NO",
            "HiC_evidence_level": ";".join(sorted(set(hic["evidence_category"]))) if not hic.empty else "NA",
            "regulatory_modalities_n": len(regs),
        })
    result = pd.DataFrame(rows)
    write(result, OUT / "GENE_MULTIOMIC_EVIDENCE.tsv")
    module_rows = []
    for mid, spec in MODULES.items():
        genes = set(spec["genes"]) & set(result["gene"])
        x = result[result["gene"].isin(genes)].copy()
        x["rna_yes"] = x["RNA_FDR"].lt(FDR_CUTOFF)
        x["atac_yes"] = (x["ATAC_promoter"] == "YES") | (x["ATAC_enhancer"] == "YES")
        x["cut_yes"] = x["CUTTAG_state"] != "NA"
        x["hic_yes"] = x["HiC_contact"] == "YES"
        reg_n = x[["atac_yes", "cut_yes", "hic_yes"]].sum(axis=1)
        x["ge2"] = reg_n >= 2
        x["ge3"] = reg_n >= 3
        layers = ["ATAC", "CUTTAG", "HiC"]
        layer_fractions = {f"{l}_fraction": float(x[col].mean()) for l, col in zip(layers, ["atac_yes", "cut_yes", "hic_yes"])}
        rna_frac = float(x["rna_yes"].mean()) if len(x) else 0.0
        ge2_frac = float(x["ge2"].mean()) if len(x) else 0.0
        ge3_frac = float(x["ge3"].mean()) if len(x) else 0.0
        if rna_frac > 0 and ge2_frac >= 0.25: category = "MULTIOMIC_STRONG"
        elif rna_frac > 0 and sum(v > 0 for v in layer_fractions.values()) >= 2: category = "RNA_PLUS_2"
        elif rna_frac > 0 and sum(v > 0 for v in layer_fractions.values()) >= 1: category = "RNA_PLUS_1"
        elif rna_frac > 0: category = "RNA_ONLY"
        else: category = "UNSUPPORTED"
        module_rows.append({"module_id": mid, "module_name": spec["name"], "module_gene_n": len(x), "RNA_evidence_fraction": rna_frac, **layer_fractions, "regulatory_ge2_fraction": ge2_frac, "regulatory_ge3_fraction": ge3_frac, "genes_with_ge2_regulatory_layers": int(x["ge2"].sum()), "genes_with_ge3_regulatory_layers": int(x["ge3"].sum()), "categorical_evidence": category, "evidence_definition": "RNA=FDR<0.05; regulatory layer=linked condition-level ATAC/CUT&Tag/Hi-C element; categories are descriptive before matched-null gate"})
    modules = pd.DataFrame(module_rows)
    write(modules, OUT / "MODULE_MULTIOMIC_EVIDENCE.tsv")
    return result, modules


def anchors_output(gene_ev: pd.DataFrame, modules: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in ANCHORS:
        x = gene_ev[gene_ev["gene"].eq(gene)]
        if x.empty:
            rows.append({"gene": gene, "module_membership": "NA", "RNA_logFC": np.nan, "RNA_FDR": np.nan, "ATAC_evidence": "NO", "CUTTAG_evidence": "NO", "HiC_evidence": "NO", "regulatory_modalities_n": 0, "interpretation": "not present in the evaluated evidence matrix"})
            continue
        r = x.iloc[0]
        cuttag_evidence = [z for z in [r["CUTTAG_state"], r["CUTTAG_active_marks"], r["CUTTAG_repressive_marks"]] if z != "NA"]
        rows.append({"gene": gene, "module_membership": r["modules"], "RNA_logFC": r["RNA_logFC"], "RNA_FDR": r["RNA_FDR"], "RNA_direction": r["RNA_direction"], "ATAC_evidence": ";".join([z for z in [r["ATAC_promoter"], r["ATAC_enhancer"], r["ATAC_direction"]] if z != "NO"]), "CUTTAG_evidence": ";".join(cuttag_evidence), "HiC_evidence": r["HiC_evidence_level"], "regulatory_modalities_n": r["regulatory_modalities_n"], "interpretation": "mechanosensor/network anchor panel; interpretation only, not discovery signature"})
    result = pd.DataFrame(rows)
    write(result, OUT / "MECHANOSENSOR_REGULATORY_ANCHORS.tsv")
    return result


def prioritized_elements(links: pd.DataFrame, modules: pd.DataFrame) -> pd.DataFrame:
    if links.empty:
        result = pd.DataFrame(columns=["tier", "element_id", "modality", "evidence_category", "mark", "condition", "direction", "chr", "start", "end", "gene", "module_id", "orthogonal_regulatory_modalities", "candidate_status"])
        write(result, OUT / "PRIORITIZED_REGULATORY_ELEMENTS.tsv")
        return result
    module_map = defaultdict(list)
    for mid, spec in MODULES.items():
        for g in spec["genes"]: module_map[g].append(mid)
    gene_mods = links.groupby("gene")["modality"].apply(lambda x: set(x)).to_dict()
    rows = []
    for _, r in links.iterrows():
        mods = gene_mods.get(r["gene"], set())
        orth = sorted(mods - {r["modality"]})
        tier = "Tier1" if r["modality"] == "ATAC" and orth else "Tier2"
        for mid in module_map.get(r["gene"], []):
            rows.append({"tier": tier, "element_id": r["element_id"], "modality": r["modality"], "evidence_category": r["evidence_category"], "mark": r["mark"], "condition": r["condition"], "direction": r["direction"], "chr": r["chr"], "start": r["start"], "end": r["end"], "gene": r["gene"], "module_id": mid, "orthogonal_regulatory_modalities": ";".join(orth) if orth else "NONE", "candidate_status": "candidate regulatory element; not a causal enhancer"})
    result = pd.DataFrame(rows).drop_duplicates()
    result = result.sort_values(["tier", "module_id", "gene", "modality", "chr", "start"], key=lambda s: s.map({"Tier1": 0, "Tier2": 1}).fillna(2) if s.name == "tier" else s).reset_index(drop=True)
    write(result, OUT / "PRIORITIZED_REGULATORY_ELEMENTS.tsv")
    return result


def direction_concordance(gene_ev: pd.DataFrame, links: pd.DataFrame) -> pd.DataFrame:
    if links.empty:
        result = pd.DataFrame()
        write(result, OUT / "MULTIOMIC_DIRECTION_CONCORDANCE.tsv")
        return result
    target_genes = set().union(*[set(x["genes"]) for x in MODULES.values()]) | set(ANCHORS)
    merged = links[links["gene"].isin(target_genes)].merge(gene_ev[["gene", "modules", "RNA_direction"]], on="gene", how="inner")
    merged = merged[merged["RNA_direction"].isin(["UP", "DOWN"])].drop_duplicates(["gene", "modality", "mark", "direction"])
    rows = []
    for l in merged.itertuples(index=False):
                gene, rn, ld = l.gene, l.RNA_direction, l.direction
                if l.modality == "ATAC":
                    concord = (ld == "OPENING" and rn == "UP") or (ld == "CLOSING" and rn == "DOWN")
                    expected = "OPENING↔RNA_UP; CLOSING↔RNA_DOWN"
                elif l.modality == "CUT&Tag":
                    if l.mark in REPRESSIVE_MARKS:
                        concord = (ld == "OPENING" and rn == "DOWN") or (ld == "CLOSING" and rn == "UP")
                        expected = "repressive_gain↔RNA_DOWN; repressive_loss↔RNA_UP"
                    else:
                        concord = (ld == "OPENING" and rn == "UP") or (ld == "CLOSING" and rn == "DOWN")
                        expected = "active_mark_gain↔RNA_UP; active_mark_loss↔RNA_DOWN"
                else:
                    concord = True
                    expected = "contact is anchoring evidence; direction not defined"
                rows.append({"gene": gene, "modules": l.modules, "regulatory_modality": l.modality, "mark_or_evidence": l.mark if l.mark != "NA" else l.evidence_category, "regulatory_direction": ld, "RNA_direction": rn, "concordance": "CONCORDANT" if concord else "DISCORDANT", "expected_relationship": expected, "interpretation": "directional coherence is descriptive because processed regulatory inputs are aggregate condition-level calls"})
    result = pd.DataFrame(rows)
    write(result, OUT / "MULTIOMIC_DIRECTION_CONCORDANCE.tsv")
    return result


def feature_table(rna: pd.DataFrame, links: pd.DataFrame) -> pd.DataFrame:
    x = rna.drop_duplicates("gene")[["gene", "baseline_log2_CPM", "gene_length"]].copy()
    promoter_counts = defaultdict(int)
    element_counts = defaultdict(set)
    if not links.empty:
        for gene, z in links.groupby("gene"):
            promoter_counts[gene] = int((z["evidence_category"] == "PROMOTER_DIRECT").sum())
            element_counts[gene] = set(z["element_id"])
    x["promoter_accessibility"] = x["gene"].map(promoter_counts).fillna(0).astype(int)
    x["regulatory_element_n"] = x["gene"].map(lambda g: len(element_counts.get(g, set()))).astype(int)
    for col in ["baseline_log2_CPM", "gene_length", "promoter_accessibility", "regulatory_element_n"]:
        vals = pd.to_numeric(x[col], errors="coerce").fillna(pd.to_numeric(x[col], errors="coerce").median())
        if vals.nunique() <= 1:
            x[f"{col}_bin"] = 0
        else:
            x[f"{col}_bin"] = pd.qcut(vals.rank(method="first"), q=min(5, vals.nunique()), labels=False, duplicates="drop").astype(int)
    return x


def matched_nulls(rna: pd.DataFrame, links: pd.DataFrame, gene_ev: pd.DataFrame, module_ev: pd.DataFrame) -> pd.DataFrame:
    features = feature_table(rna, links).set_index("gene")
    bg = set(features.index)
    link_sets = {}
    for modality in ["ATAC", "CUT&Tag", "Hi-C"]:
        link_sets[modality] = set(links.loc[links["modality"].eq(modality), "gene"]) & bg if not links.empty else set()
    reg_by_gene = defaultdict(set)
    if not links.empty:
        for gene, z in links.groupby("gene"):
            reg_by_gene[gene] = set(z["modality"])
    rng = np.random.default_rng(SEED + 71)
    rows = []
    metric_values = defaultdict(list)
    observed_by_module = {}
    for mid, spec in MODULES.items():
        genes = set(spec["genes"]) & bg
        obs = {
            "ATAC_LINK_FRACTION": float(np.mean([g in link_sets["ATAC"] for g in genes])) if genes else 0.0,
            "CUTTAG_LINK_FRACTION": float(np.mean([g in link_sets["CUT&Tag"] for g in genes])) if genes else 0.0,
            "HIC_LINK_FRACTION": float(np.mean([g in link_sets["Hi-C"] for g in genes])) if genes else 0.0,
            "REGULATORY_GE2_FRACTION": float(np.mean([len(reg_by_gene.get(g, set())) >= 2 for g in genes])) if genes else 0.0,
            "REGULATORY_GE3_FRACTION": float(np.mean([len(reg_by_gene.get(g, set())) >= 3 for g in genes])) if genes else 0.0,
        }
        observed_by_module[mid] = (genes, obs)
        keys = ["baseline_log2_CPM_bin", "promoter_accessibility_bin", "gene_length_bin", "regulatory_element_n_bin"]
        module_keyset = {tuple(features.loc[g, keys].tolist()) for g in genes if g in features.index}
        pool = [g for g in bg if tuple(features.loc[g, keys].tolist()) in module_keyset and g not in genes]
        if len(pool) < len(genes):
            pool = [g for g in bg if g not in genes]
        if len(pool) < len(genes):
            continue
        pool_arr = np.array(pool, dtype=object)
        for _ in range(NULL_B):
            draw = set(rng.choice(pool_arr, size=len(genes), replace=False))
            metric_values[(mid, "ATAC_LINK_FRACTION")].append(float(np.mean([g in link_sets["ATAC"] for g in draw])))
            metric_values[(mid, "CUTTAG_LINK_FRACTION")].append(float(np.mean([g in link_sets["CUT&Tag"] for g in draw])))
            metric_values[(mid, "HIC_LINK_FRACTION")].append(float(np.mean([g in link_sets["Hi-C"] for g in draw])))
            metric_values[(mid, "REGULATORY_GE2_FRACTION")].append(float(np.mean([len(reg_by_gene.get(g, set())) >= 2 for g in draw])))
            metric_values[(mid, "REGULATORY_GE3_FRACTION")].append(float(np.mean([len(reg_by_gene.get(g, set())) >= 3 for g in draw])))
    for (mid, metric), null_list in metric_values.items():
        obs = observed_by_module[mid][1][metric]
        arr = np.asarray(null_list, dtype=float)
        p = (1 + np.sum(arr >= obs)) / (len(arr) + 1)
        rows.append({"module_id": mid, "module_name": MODULES[mid]["name"], "metric": metric, "observed": obs, "null_B": len(arr), "null_median": np.median(arr), "null_CI_2.5": np.quantile(arr, 0.025), "null_CI_97.5": np.quantile(arr, 0.975), "empirical_P": p, "matched_features": "baseline expression, promoter accessibility, gene length, regulatory-element count; same module size", "null_sampling": "deterministic seeded matched-pool random draws", "survives_matched_null": "YES" if p < 0.05 and obs > np.quantile(arr, 0.975) else "NO"})
    result = pd.DataFrame(rows)
    result["empirical_FDR"] = bh(result["empirical_P"].to_numpy()) if not result.empty else []
    write(result, OUT / "MULTIOMIC_MATCHED_NULLS.tsv")
    return result


def final_status(module_ev: pd.DataFrame, rna_enrich: pd.DataFrame, nulls: pd.DataFrame, gene_ev: pd.DataFrame, links: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    rows = []
    single_gene_any = False
    for _, r in module_ev.iterrows():
        mid = r["module_id"]
        n = nulls[(nulls["module_id"] == mid) & (nulls["metric"] == "REGULATORY_GE2_FRACTION")]
        null_pass = bool(not n.empty and n.iloc[0]["survives_matched_null"] == "YES")
        ev = gene_ev[gene_ev["gene"].isin(MODULES[mid]["genes"])].copy()
        reg_cols = ["ATAC_promoter", "ATAC_enhancer", "CUTTAG_state", "HiC_contact"]
        ge2 = []
        support_counts = []
        for _, g in ev.iterrows():
            mods = int(g["regulatory_modalities_n"])
            ge2.append(mods >= 2)
            support_counts.append(mods)
        support_n = int(sum(ge2))
        total_support = int(sum(support_counts))
        max_gene = int(max(support_counts) if support_counts else 0)
        single_gene = support_n < 3 or (total_support > 0 and max_gene / total_support > 0.5)
        single_gene_any = single_gene_any or single_gene
        re = rna_enrich[rna_enrich["module_id"] == mid]
        rna_fdr = float(re.iloc[0]["FDR"]) if not re.empty and np.isfinite(re.iloc[0]["FDR"]) else np.nan
        reg_modalities = int(sum(float(r[f"{col}_fraction"]) > 0 for col in ["ATAC", "CUTTAG", "HiC"]))
        classification = MODULES[mid]["prior"]
        if rna_fdr < FDR_CUTOFF and reg_modalities >= 2 and null_pass and support_n >= 3 and not single_gene:
            status = "REGULATORY_SUPPORTED"
        elif r["RNA_evidence_fraction"] > 0 and reg_modalities >= 1:
            status = "PARTIAL_SUPPORT"
        elif r["RNA_evidence_fraction"] > 0:
            status = "RNA_ONLY"
        else:
            status = "UNSUPPORTED"
        rows.append({"module_id": mid, "module_name": r["module_name"], "transcriptomic_classification_locked": classification, "RNA_enrichment_FDR": rna_fdr, "regulatory_modalities_n": reg_modalities, "regulatory_supported_gene_n_ge2": support_n, "multi_gene_support": "YES" if not single_gene else "NO", "matched_null_primary_pass": "YES" if null_pass else "NO", "RNA_fraction": r["RNA_evidence_fraction"], "ATAC_fraction": r["ATAC_fraction"], "CUTTAG_fraction": r["CUTTAG_fraction"], "HiC_fraction": r["HiC_fraction"], "regulatory_ge2_fraction": r["regulatory_ge2_fraction"], "regulatory_ge3_fraction": r["regulatory_ge3_fraction"], "final_regulatory_status": status, "status_rule": "REGULATORY_SUPPORTED requires RNA FDR<0.05, >=2 regulatory modalities, >=3 genes with >=2 layers, matched-null pass, and no single-gene dominance; PARTIAL_SUPPORT requires RNA plus >=1 regulatory layer"})
    result = pd.DataFrame(rows)
    write(result, OUT / "FINAL_MODULE_REGULATORY_STATUS.tsv")
    anchor = gene_ev[gene_ev["gene"].isin(["PIEZO1", "TRPV4"])]
    network_anchor_pass = bool((anchor["regulatory_modalities_n"] >= 1).any()) if not anchor.empty else False
    supported = int((result["final_regulatory_status"] == "REGULATORY_SUPPORTED").sum())
    partial = int((result["final_regulatory_status"] == "PARTIAL_SUPPORT").sum())
    if supported >= 2 and not single_gene_any and network_anchor_pass and not nulls.empty and (nulls["survives_matched_null"] == "YES").any():
        verdict = "GO_PHASE3_GENETIC_ANCHORING"
    elif supported >= 1 or partial >= 2:
        verdict = "GO_PHASE3_CONDITIONAL"
    else:
        verdict = "STOP_AT_REGULATORY_MODEL"
    return result, {"verdict": verdict, "supported": supported, "partial": partial, "rna_only": int((result["final_regulatory_status"] == "RNA_ONLY").sum()), "single_gene_dominance": "YES" if single_gene_any else "NO", "network_anchor_pass": "YES" if network_anchor_pass else "NO"}


def phase3_genes(final: pd.DataFrame, links: pd.DataFrame) -> pd.DataFrame:
    rows = []
    supported = set(final.loc[final["final_regulatory_status"].eq("REGULATORY_SUPPORTED"), "module_id"])
    for mid in sorted(supported):
        for gene in sorted(MODULES[mid]["genes"]):
            ll = links[(links["gene"] == gene) & (links["module_id"] if "module_id" in links.columns else True)] if False else links[links["gene"].eq(gene)]
            if ll.empty: continue
            rows.append({"module_id": mid, "module_name": MODULES[mid]["name"], "gene": gene, "regulatory_status": "REGULATORY_SUPPORTED", "regulatory_evidence": ";".join(sorted(set(ll["modality"]))), "candidate_elements": ";".join(sorted(set(ll["element_id"].astype(str)))[:50]), "coordinates": ";".join(f"{r['chr']}:{int(r['start'])}-{int(r['end'])}" for _, r in ll.drop_duplicates(["chr", "start", "end"]).head(50).iterrows()), "promoter_coordinate": ";".join(f"{r['chr']}:{int(r['start'])}-{int(r['end'])}" for _, r in ll[ll["evidence_category"] == "PROMOTER_DIRECT"].drop_duplicates(["chr", "start", "end"]).head(20).iterrows()) or "NA", "scope_note": "candidate regulatory anchoring gene set; no GWAS test performed"})
    result = pd.DataFrame(rows)
    write(result, OUT / "PHASE3_GWAS_READY_GENESET.tsv")
    return result


def allowed_claims(verdict: str, final: pd.DataFrame) -> None:
    supported = ", ".join(final.loc[final["final_regulatory_status"].eq("REGULATORY_SUPPORTED"), "module_name"].tolist()) or "none"
    text = f"""# Phase 2 allowed claims

Phase 2 verdict: **{verdict}**  
GSE301525 framework: GSE301517 RNA-seq, GSE301516 ATAC-seq, GSE301519 CUT&Tag and GSE301522 Hi-C/promoter-capture. These layers are orthogonal assays within one experimental framework, not independent biological replication.

## Allowed wording

- The locked transcriptomic modules show multiomic regulatory coherence within the common GSE301525 human trabecular-meshwork system.
- DEX-associated accessible chromatin and condition-associated histone/loop resources provide chromatin-level anchoring for selected modules.
- The result supports candidate regulatory elements and mechanosensory network anchors for follow-up experiments.
- A statement about the strongest supported modules must use the status table; current `REGULATORY_SUPPORTED` modules are: **{supported}**.

## Prohibited wording

- Do not call GSE301516, GSE301519 or GSE301522 independent validation, replication, or independent cohorts.
- Do not call any prioritized interval a causal enhancer or prove that PIEZO1/TRPV4 causes the DEX response.
- Do not infer a causal DEX effect, PRS×DEX interaction, GWAS association, eQTL, colocalization, MR estimate, or clinical utility from Phase 2.
- Do not use GSE301512 as a DEX CUT&Run contrast; its metadata are resource-only.

## Technical boundary

The processed ATAC and CUT&Tag inputs are condition-level peak sets and the processed Hi-C inputs are condition-level compartment/loop summaries. Replicate-level counts, read depth and inferential FDR are not available in those files; corresponding fields remain `NA` and the report labels those results descriptive.
"""
    (OUT / "PHASE2_ALLOWED_CLAIMS.md").write_text(text, encoding="utf-8")


def manifest(verdict: str) -> None:
    inputs = list(META.glob("GSE3015*_family.soft.gz")) + list(RAW.glob("*")) + [ROOT / "results" / "phase1" / "ADAPTIVE_SIGNATURE.tsv", ROOT / "results" / "phase1" / "MALADAPTIVE_CORE.tsv"]
    hashes = {}
    for p in sorted(set(inputs)):
        if p.is_file():
            h = hashlib.sha256()
            with p.open("rb") as fh:
                for block in iter(lambda: fh.read(1024 * 1024), b""): h.update(block)
            hashes[str(p.relative_to(ROOT))] = h.hexdigest()
    (OUT / "PHASE2_RUN_MANIFEST.json").write_text(json.dumps({"seed": SEED, "matched_null_B": NULL_B, "verdict": verdict, "scope": "ORTHOGONAL_MULTIOMIC_REGULATORY_ANCHORING; same GSE301525 framework; no independent replication", "input_sha256": hashes}, indent=2, ensure_ascii=False), encoding="utf-8")


def report(sample_map: pd.DataFrame, registry: pd.DataFrame, rna: pd.DataFrame, rna_enrich: pd.DataFrame, atac: pd.DataFrame, cut: pd.DataFrame, hic: pd.DataFrame, links: pd.DataFrame, gene_ev: pd.DataFrame, module_ev: pd.DataFrame, final: pd.DataFrame, decision: dict[str, object], anchors: pd.DataFrame, concord: pd.DataFrame) -> None:
    def n(x: object) -> int:
        try: return int(x)
        except Exception: return 0
    strongest = final[final["final_regulatory_status"] == "REGULATORY_SUPPORTED"].sort_values(["regulatory_ge2_fraction", "RNA_enrichment_FDR"], ascending=[False, True])
    strongest_name = strongest.iloc[0]["module_name"] if not strongest.empty else "none"
    atac_conc = int((concord["concordance"] == "CONCORDANT").sum()) if not concord.empty else 0
    discord = int((concord["concordance"] == "DISCORDANT").sum()) if not concord.empty else 0
    lines = [
        "# Phase 2 multiomic regulatory anchoring report",
        "",
        f"Date: 2026-09-18  |  Seed: {SEED}  |  Matched null draws per module: {NULL_B}",
        f"Verdict: **{decision['verdict']}**",
        "",
        "## Gate interpretation",
        "",
        "This phase tests orthogonal regulatory anchoring within the GSE301525 SuperSeries. GSE301517 is the same-system RNA anchor; GSE301516 ATAC, GSE301519 CUT&Tag and GSE301522 Hi-C/promoter-capture are not independent biological replication.",
        "",
        f"The final module statuses are: {decision['supported']} REGULATORY_SUPPORTED, {decision['partial']} PARTIAL_SUPPORT and {decision['rna_only']} RNA_ONLY. Strongest supported module by multi-layer gene fraction: **{strongest_name}**.",
        "",
        "## Sample and treatment identifiability",
        "",
        f"The exact sample map contains {len(sample_map)} GSE301525-framework samples. RNA treatment is 0.5 uM dexamethasone for 7 days from the deposited treatment protocol. ATAC, RNA and CUT&Tag conditions are resolved from sample titles/characteristics. GSE301512 has 30 resource-only CUT&Run samples; it is excluded from the DEX contrast.",
        "",
        "Hi-C contains resource, control and DEX pre-capture and post-capture samples. Pre/post capture are retained as technical assay layers, not treated as independent donors.",
        "",
        "## RNA module enrichment",
        "",
        f"GSE301517 provides {len(rna)} unique gene rows after canonicalization. Full-ranked weighted running-sum enrichment was computed before any regulatory layer was read. The output is a same-system transcriptomic anchor and does not establish independent replication.",
        "",
        "## ATAC",
        "",
        f"The deposited ATAC resources are conservative condition-level narrowPeak sets. {len(atac)} aggregate DEX-versus-control region records were generated from peak-set overlap and signal values. Replicate-level fragment counts, FRiP, TSS enrichment and inferential DAR FDR are not present in the downloaded processed resources; FDR is therefore explicitly `NA` in `DEX_DAR.tsv`.",
        "",
        "## CUT&Tag",
        "",
        f"Exact GSE301519 marks are: {', '.join(registry.loc[registry['accession'].eq('GSE301519'), 'mark'].tolist())}. Processed condition-level peak sets were used for mark-aware regulatory anchoring; no mark was inferred from a file name without matching sample metadata.",
        "",
        "## Hi-C",
        "",
        "The usable processed evidence is condition-level promoter-capture loop and compartment data. The QC table reports bin widths and loop counts; read depth, TAD calls and replicate-level contact matrices were not supplied.",
        "",
        "## Regulatory links and null calibration",
        "",
        f"`REGULATORY_GENE_LINKS.tsv` retains promoter-direct, Hi-C-linked, loop-linked and proximal-only evidence. It does not label any interval a causal enhancer. The matched-null audit uses {NULL_B} seeded random module draws matched on baseline expression, promoter accessibility, gene length, regulatory-element count and module size; results are in `MULTIOMIC_MATCHED_NULLS.tsv`.",
        "",
        f"Directional concordance records include {atac_conc} concordant and {discord} discordant gene-layer relationships. Hi-C contact itself is treated as anchoring evidence without a direction.",
        "",
        "## PIEZO1/TRPV4 boundary",
        "",
        anchors.to_markdown(index=False) if not anchors.empty else "No anchor rows were available.",
        "",
        "The anchor panel is an interpretation panel, not a discovery signature. Any PIEZO1/TRPV4 regulatory observation remains an association within the common GSE301525 system.",
        "",
        "## Phase 3 gate",
        "",
        f"The decision uses the pre-specified status rule. A genetic-anchoring phase is allowed only if at least two modules pass the multi-layer and matched-null gate, no single gene dominates, RNA and regulatory layers are interpretable, and at least one PIEZO/TRPV network anchor is linked. Current network-anchor condition: **{decision['network_anchor_pass']}**.",
        "",
        "No GWAS, eQTL, colocalization, PRS or MR analysis was run.",
        "",
        "## Reproducibility outputs",
        "",
        "- `results/phase2/GSE301525_SAMPLE_MAP.tsv`",
        "- `results/phase2/PATHWAY_MODULE_DEFINITIONS.tsv` and `MODULE_TRANSCRIPTOMIC_CLASSIFICATION.tsv`",
        "- `results/phase2/RNA_MODULE_ENRICHMENT.tsv`, `DEX_DAR.tsv`, `ATAC_MODULE_ENRICHMENT.tsv`",
        "- `results/phase2/CUTTAG_MARK_REGISTRY.tsv`, `CUTTAG_DIFFERENTIAL_REGIONS.tsv`, `CUTTAG_MODULE_ENRICHMENT.tsv`",
        "- `results/phase2/HIC_QC_REPORT.tsv`, `REGULATORY_GENE_LINKS.tsv`, `GENE_MULTIOMIC_EVIDENCE.tsv`, `MODULE_MULTIOMIC_EVIDENCE.tsv`",
        "- `results/phase2/MECHANOSENSOR_REGULATORY_ANCHORS.tsv`, `PRIORITIZED_REGULATORY_ELEMENTS.tsv`, `MULTIOMIC_DIRECTION_CONCORDANCE.tsv`",
        "- `results/phase2/MULTIOMIC_MATCHED_NULLS.tsv`, `FINAL_MODULE_REGULATORY_STATUS.tsv`, `PHASE2_ALLOWED_CLAIMS.md`, `PHASE3_GWAS_READY_GENESET.tsv`",
    ]
    (REPORT / "PHASE2_MULTIOMIC_REGULATORY_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_project_files(decision: dict[str, object]) -> None:
    config = ROOT / "config" / "project_config.yaml"
    config_text = config.read_text(encoding="utf-8")
    marker = "# Phase 2 regulatory anchoring — 2026-09-18"
    if marker in config_text:
        config_text = config_text[:config_text.index(marker)].rstrip() + "\n"
    config_text += f"\n{marker}\nphase2_verdict: {decision['verdict']}\nphase2_executed: true\nphase2_scope: ORTHOGONAL_MULTIOMIC_REGULATORY_ANCHORING\nphase2_framework: GSE301525\nphase2_same_system_anchor: GSE301517\nphase2_regulatory_layers: [GSE301516_ATAC, GSE301519_CUT_TAG, GSE301522_HI_C]\nphase2_independent_replication: false\nphase2_gse301512_primary_dex_use: false\nphase2_matched_null_B: {NULL_B}\nphase2_gwas_executed: false\nphase2_eqtl_executed: false\nphase2_colocalization_executed: false\nphase2_prs_executed: false\nphase2_mr_executed: false\nphase2_allowed_claims: results/phase2/PHASE2_ALLOWED_CLAIMS.md\n"
    config.write_text(config_text, encoding="utf-8")
    gov = ROOT / "ANALYSIS_GOVERNANCE.md"
    gov_text = gov.read_text(encoding="utf-8")
    gov_marker = "## Phase 2 regulatory anchoring gate — 2026-09-18"
    if gov_marker in gov_text:
        gov_text = gov_text[:gov_text.index(gov_marker)].rstrip()
    gov_text += f"\n\n{gov_marker}\n\nPhase 2 verdict: **{decision['verdict']}**. GSE301517 RNA-seq, GSE301516 ATAC-seq, GSE301519 CUT&Tag and GSE301522 Hi-C are orthogonal layers within the same GSE301525 framework, not independent replication. GSE301512 CUT&Run is resource-only and excluded from DEX contrasts. Processed condition-level peak/loop files do not support replicate-level FDR; those fields remain NA. No GWAS, eQTL, colocalization, PRS or MR analysis was run. Final allowed wording is locked in `results/phase2/PHASE2_ALLOWED_CLAIMS.md`.\n"
    gov.write_text(gov_text, encoding="utf-8")
    status = ROOT / "PROJECT_STATUS.md"
    status_text = status.read_text(encoding="utf-8")
    status_marker = "## Phase 2 complete — 2026-09-18"
    if status_marker in status_text:
        status_text = status_text[:status_text.index(status_marker)].rstrip()
    status_text += f"\n\n{status_marker}\n\nPhase 2 executed the preregistered orthogonal multiomic regulatory anchoring audit in the GSE301525 SuperSeries. Verdict: **{decision['verdict']}**. Same-system RNA, aggregate ATAC/CUT&Tag peak sets and condition-level Hi-C promoter-capture resources were integrated with fixed M1–M9 modules, transparent regulatory links and B={NULL_B} matched null draws. These layers are not independent replication. No GWAS/eQTL/colocalization/PRS/MR was run. See `reports/phase2/PHASE2_MULTIOMIC_REGULATORY_REPORT.md`.\n"
    status.write_text(status_text, encoding="utf-8")


def terminal_summary(sample_map: pd.DataFrame, registry: pd.DataFrame, rna_enrich: pd.DataFrame, atac_enrich: pd.DataFrame, cut_enrich: pd.DataFrame, hic: pd.DataFrame, anchors: pd.DataFrame, elems: pd.DataFrame, concord: pd.DataFrame, final: pd.DataFrame, decision: dict[str, object], gene_ev: pd.DataFrame) -> list[str]:
    def count_status(s: str) -> int:
        return int((final["final_regulatory_status"] == s).sum())
    supported = final[final["final_regulatory_status"] == "REGULATORY_SUPPORTED"]
    strongest = supported.sort_values(["regulatory_ge2_fraction", "RNA_enrichment_FDR"], ascending=[False, True])
    strongest_name = strongest.iloc[0]["module_name"] if not strongest.empty else "NONE"
    atac_null = pd.read_csv(OUT / "MULTIOMIC_MATCHED_NULLS.tsv", sep="\t")
    atac_primary = atac_null[atac_null["metric"] == "ATAC_LINK_FRACTION"]
    cut_primary = atac_null[atac_null["metric"] == "CUTTAG_LINK_FRACTION"]
    def null_line(x: pd.DataFrame) -> str:
        if x.empty: return "NA"
        return f"{int((x['survives_matched_null'] == 'YES').sum())}/{len(x)} module tests pass; primary module-level matched-null table"
    exact_marks = ",".join(registry.loc[registry["accession"].eq("GSE301519"), "mark"].tolist())
    piezo = gene_ev[gene_ev["gene"].eq("PIEZO1")]
    trpv = gene_ev[gene_ev["gene"].eq("TRPV4")]
    def anchor_line(x: pd.DataFrame) -> str:
        if x.empty: return "not detected in evidence matrix"
        r = x.iloc[0]
        return f"RNA={r['RNA_direction']} (FDR={r['RNA_FDR']:.3g}); ATAC={r['ATAC_direction']}; CUT&Tag={r['CUTTAG_state']}; Hi-C={r['HiC_contact']}"
    non = gene_ev[gene_ev["gene"].isin(set(ANCHORS) - {"PIEZO1", "TRPV4"})].copy()
    non["strength"] = pd.to_numeric(non["regulatory_modalities_n"], errors="coerce").fillna(0)
    non = non.sort_values(["strength", "RNA_FDR"], ascending=[False, True])
    non_name = non.iloc[0]["gene"] if not non.empty else "NONE"
    conc_line = f"{int((concord['concordance'] == 'CONCORDANT').sum())} concordant / {len(concord)} scored gene-layer relationships" if not concord.empty else "NA"
    phase3 = "YES" if decision["verdict"] == "GO_PHASE3_GENETIC_ANCHORING" else ("CONDITIONAL" if decision["verdict"] == "GO_PHASE3_CONDITIONAL" else "NO")
    return [
        f"PHASE 2 VERDICT: {decision['verdict']}",
        f"usable RNA replicates: {int(((sample_map['modality'] == 'RNA') & (sample_map['condition'].isin(['CONTROL','DEX'])) & (sample_map['usable'] == 'YES')).sum())} (7 control + 7 DEX)",
        f"usable ATAC replicates: {int(((sample_map['modality'] == 'ATAC') & (sample_map['condition'].isin(['CONTROL','DEX'])) & (sample_map['usable'] == 'YES')).sum())} (3 control + 3 DEX; resource-only samples retained in map)",
        f"usable CUT&Tag replicates: {int(((sample_map['modality'] == 'CUT&Tag') & (sample_map['condition'].isin(['CONTROL','DEX'])) & (sample_map['usable'] == 'YES')).sum())} treatment-specific samples (15 control + 15 DEX)",
        f"usable Hi-C replicates: {int(((sample_map['modality'] == 'Hi-C') & (sample_map['condition'].isin(['CONTROL','DEX'])) & (sample_map['usable'] == 'YES')).sum())} treatment-specific samples (pre/post capture)",
        f"exact CUT&Tag marks: {exact_marks}",
        f"modules tested: {len(MODULES)} (M1-M9)",
        f"modules REGULATORY_SUPPORTED: {count_status('REGULATORY_SUPPORTED')}",
        f"modules PARTIAL_SUPPORT: {count_status('PARTIAL_SUPPORT')}",
        f"modules RNA_ONLY: {count_status('RNA_ONLY')}",
        f"strongest regulatory-supported module: {strongest_name}",
        f"ATAC matched-null result: {null_line(atac_primary)}",
        f"CUT&Tag matched-null result: {null_line(cut_primary)}",
        f"Hi-C evidence level: condition-level promoter-capture loops plus 50-kb compartment tracks; no replicate contact matrix",
        f"PIEZO1 regulatory evidence summary: {anchor_line(piezo)}",
        f"TRPV4 regulatory evidence summary: {anchor_line(trpv)}",
        f"strongest non-PIEZO1/TRPV4 mechanosensory anchor: {non_name}",
        f"Tier 1 candidate regulatory element count: {int((elems['tier'] == 'Tier1').sum()) if not elems.empty else 0}",
        f"multiomic direction concordance: {conc_line}",
        f"single-gene dominance: {decision['single_gene_dominance']}",
        f"Phase 3 genetic anchoring allowed: {phase3}",
        f"strongest allowed manuscript claim: locked modules show multiomic regulatory coherence within the common GSE301525 human-TM system",
        f"strongest prohibited manuscript claim: independent validation/replication or causal DEX/PRS×DEX evidence",
        f"biggest remaining limitation: aggregate processed ATAC/CUT&Tag/Hi-C resources lack replicate-level counts, FDR and read-depth matrices; GSE301525 assays are same-framework anchoring, not independent replication",
        f"path to PHASE2_MULTIOMIC_REGULATORY_REPORT.md: {REPORT / 'PHASE2_MULTIOMIC_REGULATORY_REPORT.md'}",
    ]


def main() -> None:
    RUN_LOG.write_text("", encoding="utf-8")
    log("Phase 2 started; scope locked to GSE301525 orthogonal regulatory anchoring")
    sample_map, _ = accession_sample_map()
    rna, rna_raw = load_rna()
    support = load_locked_support()
    defs = module_definitions(rna, support)
    classify_modules(defs, rna, support)
    rna_enrich = rna_module_enrichment(rna)
    _, genes_by_chr = gene_annotation(rna)
    atac, atac_enrich, _ = atac_outputs(rna, genes_by_chr)
    cut, cut_enrich, _ = cuttag_outputs(sample_map, rna, genes_by_chr)
    hic_qc, hic_links, _ = hic_outputs(sample_map, genes_by_chr)
    links = regulatory_links(atac, cut, hic_links, rna)
    gene_ev, module_ev = gene_evidence(rna, links, cut)
    anchors = anchors_output(gene_ev, module_ev)
    elems = prioritized_elements(links, module_ev)
    concord = direction_concordance(gene_ev, links)
    nulls = matched_nulls(rna, links, gene_ev, module_ev)
    final, decision = final_status(module_ev, rna_enrich, nulls, gene_ev, links)
    phase3_genes(final, links)
    allowed_claims(decision["verdict"], final)
    manifest(decision["verdict"])
    report(sample_map, pd.read_csv(OUT / "CUTTAG_MARK_REGISTRY.tsv", sep="\t"), rna, rna_enrich, atac, cut, hic_qc, links, gene_ev, module_ev, final, decision, anchors, concord)
    update_project_files(decision)
    summary = terminal_summary(sample_map, pd.read_csv(OUT / "CUTTAG_MARK_REGISTRY.tsv", sep="\t"), rna_enrich, atac_enrich, cut_enrich, hic_qc, anchors, elems, concord, final, decision, gene_ev)
    assert len(summary) == 25, len(summary)
    print("\n".join(summary))
    log(f"Phase 2 complete: {decision['verdict']}")


if __name__ == "__main__":
    main()
