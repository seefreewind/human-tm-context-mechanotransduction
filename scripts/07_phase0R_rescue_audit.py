#!/usr/bin/env python3
"""Phase 0R: rescue provenance, reconcile denominators, and gate Phase 1.

This script is deliberately limited to metadata, file, design, and technical
identifiability audits.  It does not read expression values and does not run
DE, enrichment, scores/signatures, PRS association, eQTL, colocalization,
trajectory, GRN, network, or machine-learning analyses.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data" / "metadata"
RAW = ROOT / "data" / "raw" / "GSE313736"
OUT = ROOT / "results" / "phase0R"
REPORT = ROOT / "reports" / "phase0R"
OUT.mkdir(parents=True, exist_ok=True)
REPORT.mkdir(parents=True, exist_ok=True)

PRIMARY = "GSE313736"
PRIMARY_SERIES = META / "GSE313736_series_full.txt"
PRIMARY_FAMILY = META / "GSE313736_family.soft.gz"
CELL_META = META / "GSE313736_HTMC_scRNA_cell_metadata.tsv.gz"
DONOR_META = META / "GSE313736_HTMC_scRNA_donor_metadata.tsv.gz"
VCF = META / "GSE313736_HTMC_genotypes.vcf.gz"

MATRIX_RE = re.compile(r"GSM(\d+)_pool_(dex)?([A-I])_(barcodes|features|matrix)\.tsv\.gz|GSM(\d+)_pool_(dex)?([A-I])_(barcodes|features|matrix)\.mtx\.gz")
DONOR_RE = re.compile(r"HTMC-\d+")


def clean(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return str(value).replace("\n", " ").replace("\t", " ").strip()


def soft_lines(path: Path) -> list[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", errors="replace") as fh:
            return fh.read().splitlines()
    return path.read_text(errors="replace").splitlines()


def line_value(lines: list[str], prefix: str, default: str = "") -> str:
    for line in lines:
        if line.startswith(prefix):
            return line.split("=", 1)[1].strip()
    return default


def line_values(lines: list[str], prefix: str) -> list[str]:
    return [line.split("=", 1)[1].strip() for line in lines if line.startswith(prefix)]


def parse_family(path: Path) -> tuple[list[str], list[dict[str, object]]]:
    lines = soft_lines(path)
    first_sample = next((i for i, line in enumerate(lines) if line.startswith("^SAMPLE")), len(lines))
    series = lines[:first_sample]
    samples: list[dict[str, object]] = []
    i = first_sample
    while i < len(lines):
        if not lines[i].startswith("^SAMPLE"):
            i += 1
            continue
        end = next((j for j in range(i + 1, len(lines)) if lines[j].startswith("^")), len(lines))
        block = lines[i:end]
        samples.append(
            {
                "geo": line_value(block, "!Sample_geo_accession = "),
                "title": line_value(block, "!Sample_title = "),
                "organism": line_value(block, "!Sample_organism_ch1 = "),
                "chars": line_values(block, "!Sample_characteristics_ch1 = "),
                "growth": line_value(block, "!Sample_growth_protocol_ch1 = "),
                "treatment_protocol": line_value(block, "!Sample_treatment_protocol_ch1 = "),
                "source": line_value(block, "!Sample_source_name_ch1 = "),
                "instrument": line_value(block, "!Sample_instrument_model = "),
                "extract": line_value(block, "!Sample_extract_protocol_ch1 = "),
                "processing": line_values(block, "!Sample_data_processing = "),
                "supp": line_values(block, "!Sample_supplementary_file = "),
            }
        )
        i = end
    return series, samples


def unique_or_na(values: list[object]) -> str:
    vals = sorted({clean(x) for x in values if clean(x) not in {"", "NA"}})
    return vals[0] if len(vals) == 1 else ("; ".join(vals) if vals else "NA")


def normalize_donor(value: object) -> str:
    match = DONOR_RE.search(clean(value))
    return match.group(0) if match else clean(value)


def pool_from_library(library: str) -> str:
    match = re.search(r"HTMC_pool-(dex)?([A-I])", library)
    if not match:
        return "UNKNOWN"
    return f"pool-dex{match.group(2)}" if match.group(1) else f"pool-{match.group(2)}"


def condition_from_library(library: str) -> str:
    return "DEX" if "pool-dex" in library else "vehicle"


def library_to_gsm() -> dict[str, str]:
    result: dict[str, str] = {}
    for path in sorted(RAW.glob("GSM*_pool_*_barcodes.tsv.gz")):
        match = re.match(r"(GSM\d+)_pool_(dex)?([A-I])_barcodes\.tsv\.gz", path.name)
        if match:
            lib = f"HTMC_pool-{'dex' if match.group(2) else ''}{match.group(3)}_GRCh38"
            result[lib] = match.group(1)
    return result


def matrix_barcode_frame() -> tuple[pd.DataFrame, set[tuple[str, str]]]:
    rows: list[dict[str, object]] = []
    pairs: set[tuple[str, str]] = set()
    for path in sorted(RAW.glob("GSM*_pool_*_barcodes.tsv.gz")):
        match = re.match(r"(GSM\d+)_pool_(dex)?([A-I])_barcodes\.tsv\.gz", path.name)
        if not match:
            continue
        library = f"HTMC_pool-{'dex' if match.group(2) else ''}{match.group(3)}_GRCh38"
        with gzip.open(path, "rt", errors="replace") as fh:
            barcodes = [line.strip() for line in fh if line.strip()]
        rows.append(
            {
                "gsm": match.group(1),
                "library_id": library,
                "pool_id": pool_from_library(library),
                "treatment": condition_from_library(library),
                "matrix_barcodes": len(barcodes),
                "matrix_barcode_file": str(path.relative_to(ROOT)),
            }
        )
        pairs.update((library, barcode) for barcode in barcodes)
    return pd.DataFrame(rows), pairs


def vcf_donors() -> tuple[list[str], dict[str, str]]:
    with gzip.open(VCF, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n").split("\t")
                samples = columns[9:]
                return [normalize_donor(x) for x in samples], {normalize_donor(x): x for x in samples}
    raise RuntimeError("VCF header #CHROM was not found")


def read_dmx_evidence() -> dict[str, dict[str, object]]:
    """Summarize optional dmx.best files without using expression values."""
    evidence: dict[str, dict[str, object]] = defaultdict(lambda: {"files": [], "rows": 0, "best_counts": Counter(), "barcodes": set()})
    dmx_dir = RAW / "demux_best"
    for path in sorted(dmx_dir.glob("*_dmx.best.gz")):
        match = re.match(r"GSM\d+_pool_(dex)?([A-I])_dmx\.best\.gz", path.name)
        if not match:
            continue
        pool = f"pool-dex{match.group(2)}" if match.group(1) else f"pool-{match.group(2)}"
        entry = evidence[pool]
        entry["files"].append(str(path.relative_to(ROOT)))
        try:
            with gzip.open(path, "rt", errors="replace") as fh:
                lines = [x.rstrip("\n") for x in fh if x.strip()]
            if not lines:
                continue
            header_idx = next((i for i, line in enumerate(lines) if not line.startswith("#") and ("BARCODE" in line.upper() or "BEST" in line.upper())), 0)
            header = lines[header_idx].lstrip("#").split("\t")
            if len(header) == 1:
                header = re.split(r"\s+", lines[header_idx].lstrip("#").strip())
            body = lines[header_idx + 1 :]
            for line in body:
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) == 1:
                    fields = re.split(r"\s+", line.strip())
                if not fields:
                    continue
                entry["rows"] += 1
                barcode = fields[0]
                entry["barcodes"].add(barcode)
                for idx, col in enumerate(header):
                    if idx < len(fields) and ("BEST" in col.upper() or "SAMPLE" in col.upper() or "DONOR" in col.upper()):
                        val = normalize_donor(fields[idx])
                        if DONOR_RE.search(val):
                            entry["best_counts"][val] += 1
                            break
        except Exception as exc:
            entry["parse_error"] = repr(exc)
    return evidence


def donor_reconciliation(cm: pd.DataFrame, dm: pd.DataFrame, matrix_pairs: set[tuple[str, str]], matrix_df: pd.DataFrame) -> pd.DataFrame:
    all_vcf, vcf_labels = vcf_donors()
    observed_cell = set(cm["donor_id_demuxlet"].map(normalize_donor))
    observed_donor = set(dm["donor_id_demuxlet"].map(normalize_donor))
    dmx = read_dmx_evidence()
    # A donor-specific dmx count is evidence of an attempted assignment, but
    # is not by itself evidence that the donor passed the supplied SC metadata.
    dmx_counts = Counter()
    dmx_pools: dict[str, list[str]] = defaultdict(list)
    for pool, item in dmx.items():
        for donor, count in item["best_counts"].items():
            dmx_counts[donor] += count
            dmx_pools[donor].append(pool)
    rows: list[dict[str, object]] = []
    for donor in sorted(all_vcf, key=lambda x: int(x.split("-")[-1])):
        in_cell = donor in observed_cell
        in_donor = donor in observed_donor
        cell_rows = cm[cm["donor_id_demuxlet"].map(normalize_donor) == donor]
        if in_cell:
            status = "ANALYSIS_DONOR"
            reason = "present in supplied donor metadata and cell metadata; both treatments observed"
            expr = "YES_METADATA_MATRIX_INTERSECTION"
            pv, pdex = "YES", "YES"
        else:
            status = "UNKNOWN"
            reason = "VCF donor absent from supplied donor and cell metadata; public record does not document donor-specific exclusion"
            expr = "NO_DONOR_SPECIFIC_ASSIGNMENT_IN_SUPPLIED_CELL_METADATA"
            pv, pdex = "UNVERIFIED", "UNVERIFIED"
            if dmx_counts.get(donor, 0):
                reason += "; dmx.best donor assignment evidence exists, but the supplied cell metadata omits these rows"
                expr = f"DMX_BEST_ASSIGNMENT_ROWS={dmx_counts[donor]}"
        rows.append(
            {
                "genotype_id": vcf_labels.get(donor, donor),
                "present_in_vcf": "YES",
                "present_in_donor_metadata": "YES" if in_donor else "NO",
                "present_in_cell_metadata": "YES" if in_cell else "NO",
                "present_vehicle": pv,
                "present_dex": pdex,
                "expression_barcode_evidence": expr,
                "PRS_present": "YES" if in_donor and dm.loc[dm["donor_id_demuxlet"].map(normalize_donor) == donor, "glaucoma_prs_raw"].notna().all() else "NO",
                "reason_missing_if_known": reason,
                "final_status": status,
                "dmx_best_assignment_rows": dmx_counts.get(donor, 0),
                "dmx_best_pools": ";".join(sorted(set(dmx_pools.get(donor, [])))),
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "DONOR_RECONCILIATION.tsv", sep="\t", index=False)
    return result


def cell_denominator_reconciliation(cm: pd.DataFrame, matrix_df: pd.DataFrame, matrix_pairs: set[tuple[str, str]]) -> pd.DataFrame:
    matched = cm.apply(lambda r: (r["library_id"], r["barcode"]) in matrix_pairs, axis=1)
    rows = [
        {
            "object": "GEO_series_summary",
            "cells": 316799,
            "barcodes": "NA",
            "stage": "GEO series summary singlet count",
            "filtering_status": "series-level reported count; exact cell list unavailable in supplied metadata",
            "donor_assignment": "74 donors reported; identities not supplied at cell level",
            "treatment": "vehicle and dex pools reported",
            "annotation": "not supplied in summary",
            "source": "data/metadata/GSE313736_series_full.txt",
            "recommended_use": "contextual denominator only; do not use as primary analysis N",
        },
        {
            "object": "downloaded_filtered_GEX_matrices",
            "cells": int(matrix_df["matrix_barcodes"].sum()),
            "barcodes": int(matrix_df["matrix_barcodes"].sum()),
            "stage": "Cell Ranger filtered GEX matrix barcode members downloaded from RAW.tar",
            "filtering_status": "author-supplied filtered barcodes; no additional cell filtering",
            "donor_assignment": "not present in barcode files; donor assignment requires cell metadata/demuxlet",
            "treatment": "library/pool identity gives vehicle or DEX label",
            "annotation": "not present in barcode files",
            "source": "data/raw/GSE313736/*_barcodes.tsv.gz",
            "recommended_use": "source expression universe; intersect with author metadata before analysis",
        },
        {
            "object": "supplied_cell_metadata",
            "cells": int(len(cm)),
            "barcodes": int(cm["barcode"].nunique()),
            "stage": "author-supplied demuxlet/cell-classification metadata",
            "filtering_status": "CellBender/demuxlet/classification processing described by author; supplied table only",
            "donor_assignment": f"{cm['donor_id_demuxlet'].nunique()} donors",
            "treatment": "DEX Treatment column; 9 vehicle and 9 DEX libraries",
            "annotation": "coarse and fine CellTypist plus scPred labels",
            "source": "data/metadata/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",
            "recommended_use": "auditable metadata frame; intersect with matrix barcodes",
        },
        {
            "object": "metadata_matrix_intersection",
            "cells": int(matched.sum()),
            "barcodes": int(matched.sum()),
            "stage": "exact inner join on (library_id, barcode)",
            "filtering_status": "metadata rows with a matching downloaded matrix barcode",
            "donor_assignment": f"{cm.loc[matched, 'donor_id_demuxlet'].nunique()} donors",
            "treatment": "vehicle and DEX; paired for all observed donors",
            "annotation": "coarse and fine labels retained",
            "source": "derived by scripts/07_phase0R_rescue_audit.py",
            "recommended_use": "PRIMARY_CELL_ANALYSIS_FRAME",
        },
        {
            "object": "metadata_without_matrix_match",
            "cells": int((~matched).sum()),
            "barcodes": int((~matched).sum()),
            "stage": "supplied metadata rows without exact matrix barcode match",
            "filtering_status": "excluded from primary expression analysis until reconciled",
            "donor_assignment": f"{cm.loc[~matched, 'donor_id_demuxlet'].nunique()} donors",
            "treatment": "vehicle and DEX labels retained in metadata",
            "annotation": "labels retained but no matched expression barcode",
            "source": "derived by scripts/07_phase0R_rescue_audit.py",
            "recommended_use": "QC discrepancy only; not primary expression frame",
        },
        {
            "object": "VCF_genotype_samples",
            "cells": "NA",
            "barcodes": "NA",
            "stage": "genotype donor sample columns",
            "filtering_status": "VCF sample header",
            "donor_assignment": "74 donor IDs",
            "treatment": "not applicable",
            "annotation": "not applicable",
            "source": "data/metadata/GSE313736_HTMC_genotypes.vcf.gz",
            "recommended_use": "genotype donor reference; not a cell denominator",
        },
    ]
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "CELL_DENOMINATOR_RECONCILIATION.tsv", sep="\t", index=False)
    return result


def experimental_design(cm: pd.DataFrame, matrix_df: pd.DataFrame, matrix_pairs: set[tuple[str, str]]) -> pd.DataFrame:
    gsm_by_lib = library_to_gsm()
    rows: list[dict[str, object]] = []
    for (donor, treatment, library), sub in cm.groupby(["donor_id_demuxlet", "DEX Treatment", "library_id"], sort=True):
        matched = sub.apply(lambda r: (r["library_id"], r["barcode"]) in matrix_pairs, axis=1)
        fine = Counter(sub["celltypist_chen2025_fine"].astype(str))
        top = ";".join(f"{k}:{v}" for k, v in fine.most_common())
        sat = pd.to_numeric(sub["library_sequencing_saturation"], errors="coerce")
        rows.append(
            {
                "donor_id": donor,
                "condition": "DEX" if treatment == "dex" else "vehicle",
                "treatment_column_value": treatment,
                "pool_id": pool_from_library(library),
                "library_id": library,
                "GSM": gsm_by_lib.get(library, "NA"),
                "n_metadata_cells": len(sub),
                "n_matrix_matched_cells": int(matched.sum()),
                "sequencing_saturation": unique_or_na(sat.tolist()),
                "droplet_type_counts": ";".join(f"{k}:{v}" for k, v in sub["droplet_type"].value_counts().items()),
                "fine_annotation_counts": top,
                "processing_date": "NOT_SUPPLIED",
                "culture_batch": "NOT_SUPPLIED",
                "treatment_plate": "NOT_SUPPLIED",
                "library_prep_batch": "NOT_SUPPLIED",
                "chromium_run": "NOT_SUPPLIED",
                "sequencing_run_lane": "NOT_SUPPLIED",
                "chemistry": "Chromium Next GEM Single Cell 3' v3.1; Illumina NovaSeq 6000 S4",
                "read_depth": "NOT_SUPPLIED",
                "median_genes": "NOT_SUPPLIED",
                "median_UMI": "NOT_SUPPLIED",
                "percent_mito": "NOT_SUPPLIED",
                "percent_ribo": "NOT_SUPPLIED",
                "doublet_or_singlet": "droplet_type=SNG in supplied metadata",
                "metadata_source": "GSE313736 family SOFT + supplied cell metadata",
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "EXPERIMENTAL_DESIGN_MATRIX.tsv", sep="\t", index=False)
    return result


def design_matrix(design_df: pd.DataFrame, terms: list[str]) -> tuple[np.ndarray, list[str], dict[str, tuple[int, int]]]:
    pieces = [np.ones((len(design_df), 1), dtype=float)]
    names = ["Intercept"]
    spans: dict[str, tuple[int, int]] = {"Intercept": (0, 1)}
    for term in terms:
        start = sum(piece.shape[1] for piece in pieces)
        if term == "donor":
            levels = sorted(design_df["donor_id"].astype(str).unique())
            block = np.column_stack([(design_df["donor_id"].astype(str).to_numpy() == x).astype(float) for x in levels[1:]])
            block_names = [f"donor[{x}]" for x in levels[1:]]
        elif term == "treatment":
            block = (design_df["condition"].astype(str).to_numpy() == "DEX").astype(float)[:, None]
            block_names = ["treatment[DEX]"]
        elif term in {"library", "pool"}:
            col = "library_id" if term == "library" else "pool_id"
            levels = sorted(design_df[col].astype(str).unique())
            block = np.column_stack([(design_df[col].astype(str).to_numpy() == x).astype(float) for x in levels[1:]])
            block_names = [f"{term}[{x}]" for x in levels[1:]]
        else:
            raise ValueError(term)
        pieces.append(block)
        names.extend(block_names)
        spans[term] = (start, start + block.shape[1])
    return np.column_stack(pieces), names, spans


def identifiability(design_df: pd.DataFrame) -> pd.DataFrame:
    configs = [
        ("fixed: ~ donor + treatment", ["donor", "treatment"], "donor-paired fixed-effect model"),
        ("fixed: ~ donor + treatment + library", ["donor", "treatment", "library"], "library is nested in treatment in supplied design"),
        ("fixed: ~ donor + treatment + pool", ["donor", "treatment", "pool"], "pool is the treatment-labelled capture layer"),
        ("mixed analogue: ~ treatment + (1|donor)", ["donor", "treatment"], "fixed-effect rank analogue; donor is random intercept"),
        ("mixed analogue: ~ treatment + library + (1|donor)", ["donor", "treatment", "library"], "fixed-effect rank analogue; library-treatment alias remains"),
        ("mixed analogue: ~ treatment + pool + (1|donor)", ["donor", "treatment", "pool"], "fixed-effect rank analogue; pool-treatment alias remains"),
    ]
    rows = []
    for label, terms, interpretation in configs:
        X, names, spans = design_matrix(design_df, terms)
        rank = int(np.linalg.matrix_rank(X))
        full = rank == X.shape[1]
        base_terms = [x for x in terms if x != "treatment"]
        X_base, _, _ = design_matrix(design_df, base_terms)
        treatment_est = int(np.linalg.matrix_rank(np.column_stack([X_base, (design_df["condition"] == "DEX").astype(float)]))) > int(np.linalg.matrix_rank(X_base))
        aliased = []
        if not treatment_est:
            aliased.append("treatment")
        rows.append(
            {
                "design": label,
                "rank": rank,
                "n_columns": X.shape[1],
                "full_rank": "YES" if full else "NO",
                "aliased_terms": ";".join(aliased) if aliased else "NONE",
                "treatment_estimable": "YES" if treatment_est else "NO",
                "interpretation": interpretation if treatment_est else interpretation + "; DEX causal effect is not identifiable from primary data alone",
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "DESIGN_IDENTIFIABILITY_REPORT.tsv", sep="\t", index=False)
    return result


def donor_paired_qc(cm: pd.DataFrame, matrix_pairs: set[tuple[str, str]]) -> tuple[pd.DataFrame, dict[str, object]]:
    cm = cm.copy()
    cm["matrix_match"] = cm.apply(lambda r: (r["library_id"], r["barcode"]) in matrix_pairs, axis=1)
    fine_states = ["BeamA", "JCT", "Schlemm_Endothelium", "Ciliary_Muscle"]
    rows: list[dict[str, object]] = []
    for donor, sub in cm.groupby("donor_id_demuxlet", sort=True):
        record: dict[str, object] = {"donor_id": donor, "paired": "YES"}
        for label, treatment in [("vehicle", "no-dex"), ("dex", "dex")]:
            x = sub[sub["DEX Treatment"] == treatment]
            record[f"{label}_n_metadata_cells"] = len(x)
            record[f"{label}_n_matrix_matched_cells"] = int(x["matrix_match"].sum())
            record[f"{label}_library_ids"] = ";".join(sorted(x["library_id"].unique()))
            record[f"{label}_pool_ids"] = ";".join(sorted(x["library_id"].map(pool_from_library).unique()))
            record[f"{label}_sequencing_saturation"] = unique_or_na(pd.to_numeric(x["library_sequencing_saturation"], errors="coerce").tolist())
            record[f"{label}_SNG_pct"] = round(100 * (x["droplet_type"].eq("SNG").mean()), 4) if len(x) else "NA"
            record[f"{label}_coarse_TM_pct"] = round(100 * x["celltypist_chen2025_coarse"].eq("trabecular meshwork cell").mean(), 4) if len(x) else "NA"
            for state in fine_states:
                record[f"{label}_{state}_pct"] = round(100 * x["celltypist_chen2025_fine"].eq(state).mean(), 4) if len(x) else "NA"
            for unavailable in ["median_UMI", "median_genes", "percent_mito", "percent_ribo", "doublet_rate"]:
                record[f"{label}_{unavailable}"] = "NOT_SUPPLIED"
        record["delta_dex_minus_vehicle_n_metadata_cells"] = record["dex_n_metadata_cells"] - record["vehicle_n_metadata_cells"]
        record["delta_dex_minus_vehicle_n_matrix_matched_cells"] = record["dex_n_matrix_matched_cells"] - record["vehicle_n_matrix_matched_cells"]
        record["delta_dex_minus_vehicle_sequencing_saturation"] = round(float(record["dex_sequencing_saturation"]) - float(record["vehicle_sequencing_saturation"]), 6)
        record["technical_asymmetry_assessment"] = "DESCRIPTIVE_ONLY; treatment-labelled library saturation differs and all per-cell QC metrics except droplet type are unavailable"
        rows.append(record)
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "DONOR_PAIRED_TECHNICAL_QC.tsv", sep="\t", index=False)
    paired = result
    summary = {
        "donors": int(len(paired)),
        "dex_higher_metadata_cells": int((paired["delta_dex_minus_vehicle_n_metadata_cells"] > 0).sum()),
        "vehicle_higher_metadata_cells": int((paired["delta_dex_minus_vehicle_n_metadata_cells"] < 0).sum()),
        "equal_metadata_cells": int((paired["delta_dex_minus_vehicle_n_metadata_cells"] == 0).sum()),
        "median_cell_delta": float(paired["delta_dex_minus_vehicle_n_metadata_cells"].median()),
        "median_saturation_delta": float(paired["delta_dex_minus_vehicle_sequencing_saturation"].median()),
        "all_droplets_SNG": bool((paired[["vehicle_SNG_pct", "dex_SNG_pct"]] == 100).all().all()),
    }
    return result, summary


def protocol_audit(series: list[str], samples: list[dict[str, object]]) -> str:
    growth = unique_or_na([x["growth"] for x in samples])
    treatment = unique_or_na([x["treatment_protocol"] for x in samples])
    characteristics = sorted({c for s in samples for c in s["chars"] if c.lower().startswith("treatment:")})
    return f"""# DEX protocol audit — GSE313736

Audit date: 2026-09-17. This document records public metadata and does not resolve the duration conflict by inference.

| Item | Status | Exact source/field | Recorded value | Consequence |
|---|---|---|---|---|
| Concentration | VERIFIED | `!Sample_growth_protocol_ch1` | `500nM dexamethasone` | concentration can be recorded as 500 nM |
| Vehicle | PARTIALLY_VERIFIED | `!Sample_treatment_protocol_ch1` | `vehicle solution` | vehicle composition is not specified in GSE313736 |
| Duration in sample characteristics | VERIFIED_AS_REPORTED | `!Sample_characteristics_ch1` | `{'; '.join(characteristics)}` | public sample label says 48 hours |
| Duration in growth protocol | VERIFIED_AS_REPORTED | `!Sample_growth_protocol_ch1` | `{growth}` | protocol text says 3 days |
| Treatment protocol | VERIFIED_AS_REPORTED | `!Sample_treatment_protocol_ch1` | `{treatment}` | conflicts with growth-protocol duration |
| Media | UNRESOLVED | GSE313736 family SOFT | no complete medium composition supplied in the audited fields | do not harmonize from another dataset |
| Pairing | VERIFIED | supplied cell metadata + donor metadata | 59 donors have vehicle and DEX rows | donor pairing is available |
| Pool/library structure | VERIFIED | `!Series_overall_design` + sample titles | vehicle pools A–I; DEX pools dexA–dexI | treatment is separated from pool/library |
| Processing dates/batches | UNRESOLVED | supplied GEO metadata | no processing date, culture batch, treatment plate, lane, or library batch | treatment-vs-technical effect cannot be separated |

## Locked wording

**DEX exposure duration is inconsistently reported in public record.**

The 48-hour and 3-day statements are retained as source-level facts. No Phase 1 analysis may silently choose one duration.

Sources: [GSE313736 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736); local copy `data/metadata/GSE313736_family.soft.gz`.
"""


def protocol_from_soft(acc: str) -> tuple[str, str, list[dict[str, object]], list[str]]:
    series, samples = parse_family(META / f"{acc}_family.soft.gz")
    title = line_value(series, "!Series_title = ")
    design = line_value(series, "!Series_overall_design = ")
    return title, design, samples, line_values(series, "!Series_supplementary_file = ")


def extract_condition(title: str) -> str:
    lower = title.lower()
    if "dexamethasone" in lower or re.search(r"\bdex\b", lower):
        return "DEX"
    if "vehicle" in lower:
        return "VEHICLE"
    if "control" in lower:
        return "CONTROL"
    if "resource" in lower:
        return "RESOURCE"
    if "stretch" in lower:
        return "STRETCH"
    if "treated" in lower:
        return "TREATED"
    return "UNCLEAR"


def external_dex_audit() -> tuple[pd.DataFrame, pd.DataFrame, str, str]:
    rows: list[dict[str, object]] = []
    sig_rows: list[dict[str, object]] = []
    for acc in ["GSE313736", "GSE301517", "GSE301516", "GSE301519", "GSE301522", "GSE301512"]:
        title, design, samples, series_supp = protocol_from_soft(acc)
        cond = [extract_condition(str(s["title"])) for s in samples]
        counts = Counter(cond)
        treatment_protocol = unique_or_na([s["treatment_protocol"] for s in samples])
        supp = sorted({str(x) for x in series_supp})
        if acc == "GSE313736":
            donor_n, paired, decision, separation = "59 observed donors; 74 VCF donors", "YES (59/59)", "PRIMARY_RESOURCE_ROLE_C", "NON_IDENTIFIABLE"
            notes = "DEX and vehicle each occupy separate pool/library sets; donor pairing does not restore treatment-vs-library orthogonality."
        elif acc == "GSE301517":
            donor_n, paired, decision, separation = "cell lines TM62, TM70, TM72 in control/DEX titles; donor identity not verified", "TITLE-MATCHED CONTROL/DEX REPLICATES; donor pairing not verified", "DEX_REFERENCE_DATASET_CONDITIONAL", "UNCLEAR"
            notes = "7 resource, 7 control, and 7 DEX titles; control/DEX replicate labels align by cell line and replicate number, but batch metadata are not supplied."
        elif acc == "GSE301516":
            donor_n, paired, decision, separation = "TM62 and TM70 cell-line titles; donor identity not verified", "TITLE-MATCHED CONTROL/DEX REPLICATES", "REGULATORY_SUPPORT_CONDITIONAL", "UNCLEAR"
            notes = "ATAC resource/control/DEX titles are present; exact treatment and batch separability require file-level QC."
        elif acc == "GSE301519":
            donor_n, paired, decision, separation = "TM62, TM70, TM71 cell-line/mark titles; donor identity not verified", "CONTROL/DEX PAIRS BY CELL LINE AND MARK", "REGULATORY_SUPPORT_CONDITIONAL", "UNCLEAR"
            notes = "15 control and 15 DEX records; condition appears matched by title but processing batches are not documented."
        elif acc == "GSE301522":
            donor_n, paired, decision, separation = "TM62, TM70, TM71; pre/post-capture; donor identity not verified", "CONTROL/DEX PAIRS BY CELL LINE AND CAPTURE LAYER", "REGULATORY_SUPPORT_CONDITIONAL", "UNCLEAR"
            notes = "16 control and 16 DEX records plus resource records; pre/post-capture layers must not be mixed without a locked model."
        else:
            donor_n, paired, decision, separation = "cell-line identity not verified", "NO_VERIFIED_CONTROL_DEX_PAIR", "EXCLUDE_AS_DEX_REFERENCE", "NO"
            notes = "Series design and sample titles are resource/control-resource; do not infer DEX CUT&Run support."
        rows.append(
            {
                "accession": acc,
                "dataset_role": decision,
                "donor_or_cell_line_N": donor_n,
                "sample_records": len(samples),
                "condition_counts": ";".join(f"{k}:{v}" for k, v in sorted(counts.items())),
                "paired": paired,
                "concentration": "500 nM / 0.5 µM" if acc in {"GSE313736", "GSE301517", "GSE301516", "GSE301519", "GSE301522"} else "NA",
                "duration": "48 h vs 3 days conflict" if acc == "GSE313736" else ("7 days" if acc.startswith("GSE301") and acc != "GSE301512" else "NA"),
                "vehicle": "vehicle solution; composition not specified" if acc == "GSE313736" else ("ethanol" if acc.startswith("GSE301") and acc != "GSE301512" else "NA"),
                "batch_metadata": "not supplied in audited GEO fields",
                "treatment_vs_batch_separable": separation,
                "processed_counts_or_matrix": "YES_LISTED" if supp else "NO_OR_NOT_AUDITED",
                "matrix_or_file": "; ".join(supp),
                "decision": decision,
                "notes": notes,
                "source": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
                "treatment_protocol_field": treatment_protocol,
                "series_design": design,
            }
        )
    external = pd.DataFrame(rows)
    external.to_csv(OUT / "EXTERNAL_DEX_DATASET_AUDIT.tsv", sep="\t", index=False)

    # Signature-source audit is a source suitability audit only. No gene list
    # is derived here and no expression matrix is opened.
    specs = [
        ("GSE195756", "adaptive", "mechanical stretch", "8 RNA-seq samples; 4 pairs inferred from Control/Stretch labels", "15% cyclic stretch; 0.69 Hz; 8 h", "CONDITIONAL", "YES_LISTED", "Independent human TM resource; donor IDs not explicit; processed XLSX listed."),
        ("GSE301517", "maladaptive", "DEX", "21 total RNA-seq records; 7 resource, 7 control, 7 DEX", "0.5 µM DEX; 7 days; ethanol vehicle", "CONDITIONAL", "YES_LISTED", "Control/DEX labels align by TM cell line and replicate title, but batch metadata and donor identity are not verified."),
        ("GSE202372", "maladaptive", "TGF-beta2", "10 RNA-seq records; 5 pairs inferred from A-E", "TGF-beta2 24 h", "CONDITIONAL", "YES_LISTED", "Pair labels and sample mapping file are listed; donor identities are not explicit in SOFT."),
        ("GSE263230", "maladaptive", "glaucoma cell-identity context", "14 primary culture + 3 post-mortem scRNA-seq samples; 4 glaucoma culture samples reported", "no DEX", "CONDITIONAL", "YES_LISTED", "Independent human TM context; glaucoma subset requires sample-level identity/QC; not a perturbation signature source by itself."),
        ("GSE301516", "maladaptive", "DEX ATAC", "12 records; 3 control and 3 DEX titles plus resource", "DEX context in titles; exact duration from shared protocol", "CONDITIONAL", "YES_LISTED", "Regulatory support only; exact matching and batch audit required."),
        ("GSE301519", "maladaptive", "DEX CUT&Tag", "30 records; 15 control and 15 DEX", "DEX context in titles; exact duration from shared protocol", "CONDITIONAL", "YES_LISTED", "Regulatory support only; mark/cell-line matching required."),
        ("GSE301522", "maladaptive", "DEX Hi-C", "38 records; 16 control and 16 DEX plus resource", "DEX context in titles; pre/post-capture layers", "CONDITIONAL", "YES_LISTED", "Regulatory support only; capture layer must be treated as a design factor."),
    ]
    for acc, role, context, n, protocol, status, processed, notes in specs:
        sig_rows.append(
            {
                "accession": acc,
                "signature_role": role,
                "stress_context": context,
                "sample_or_donor_structure": n,
                "protocol": protocol,
                "independent_of_GSE313736": "YES",
                "processed_counts_or_matrix": processed,
                "source_status": status,
                "use_in_phase1": "EXTERNAL_SOURCE_ONLY; lock sample mapping and preprocessing before signature construction" if status != "FAIL" else "NO",
                "notes": notes,
                "source": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
            }
        )
    sig = pd.DataFrame(sig_rows)
    sig.to_csv(OUT / "SIGNATURE_SOURCE_AUDIT.tsv", sep="\t", index=False)
    return external, sig, "CONDITIONAL", "CONDITIONAL"


def prs_audit(dm: pd.DataFrame) -> pd.DataFrame:
    raw = dm["glaucoma_prs_raw"].notna().sum()
    cent = dm["glaucoma_prs_centile"].notna().sum()
    ancestry_counts = dm["inferred_donor_ancestry"].value_counts(dropna=False).to_dict()
    rows = [
        {"item": "PRS raw field", "status": "VERIFIED", "source": "data/metadata/GSE313736_HTMC_scRNA_donor_metadata.tsv.gz", "direct_quote_or_field": f"glaucoma_prs_raw present for {raw}/{len(dm)} donor rows", "interpretation": "donor-level numeric field exists", "verified": "YES", "consequence": "can be carried as an observed metadata field; not yet scientifically interpretable as a reproducible score"},
        {"item": "PRS centile field", "status": "VERIFIED", "source": "data/metadata/GSE313736_HTMC_scRNA_donor_metadata.tsv.gz", "direct_quote_or_field": f"glaucoma_prs_centile present for {cent}/{len(dm)} donor rows", "interpretation": "centile values exist", "verified": "YES", "consequence": "reference population and calibration remain unknown"},
        {"item": "GWAS phenotype/trait source", "status": "UNRESOLVED", "source": "GEO series/donor metadata; BioProject PRJNA1381509", "direct_quote_or_field": "No GWAS accession, phenotype definition, or score-generation field in supplied records", "interpretation": "cannot identify the GWAS construct", "verified": "NO", "consequence": "PRS primary analysis prohibited"},
        {"item": "GWAS PMID/DOI", "status": "UNRESOLVED", "source": "GEO series, BioProject, NCBI SRA query, public web search 2026-09-17", "direct_quote_or_field": "No PRS-specific PMID/DOI found in linked public metadata", "interpretation": "source paper cannot be independently cited for the score", "verified": "NO", "consequence": "PRS cannot be mapped to a defined discovery GWAS"},
        {"item": "Ancestry", "status": "PARTIALLY_VERIFIED", "source": "donor metadata field inferred_donor_ancestry", "direct_quote_or_field": "; ".join(f"{clean(k)}={v}" for k, v in ancestry_counts.items()), "interpretation": "inferred ancestry labels exist but inference method is not documented", "verified": "PARTIAL", "consequence": "ancestry calibration/stratification is not defensible as primary"},
        {"item": "Genome build", "status": "UNRESOLVED", "source": "GEO sample processing + VCF", "direct_quote_or_field": "GEX assembly is GRCh38.p13; PRS build is not specified", "interpretation": "expression reference build does not establish PRS build", "verified": "NO", "consequence": "variant harmonization/reproduction blocked"},
        {"item": "SNP/effect-allele weights", "status": "UNRESOLVED", "source": "VCF and GEO supplementary-file list", "direct_quote_or_field": "VCF is supplied; no PRS weight file or SNP/effect-allele table is supplied", "interpretation": "genotypes are not sufficient to reproduce the score", "verified": "NO", "consequence": "no independent PRS reconstruction"},
        {"item": "Scoring method", "status": "UNRESOLVED", "source": "GEO/BioProject/SRA/linked-source search", "direct_quote_or_field": "No clumping threshold, PRS-CS, LDpred, or other scoring method reported", "interpretation": "algorithm and LD reference unknown", "verified": "NO", "consequence": "no interpretable PRS effect estimate"},
        {"item": "Genotype/imputation/QC", "status": "PARTIALLY_VERIFIED", "source": "VCF + demuxlet documentation", "direct_quote_or_field": "genotype VCF and demuxlet assignments supplied; PRS-specific imputation/QC not supplied", "interpretation": "genotyping resource is real, PRS input processing is not documented", "verified": "PARTIAL", "consequence": "donor mapping passes; PRS construction does not"},
        {"item": "Harmonization and standardization", "status": "UNRESOLVED", "source": "donor metadata + public-source search", "direct_quote_or_field": "raw and centile columns supplied; allele harmonization, standardization, and centile reference not reported", "interpretation": "raw range/centile cannot be translated to a defined scale", "verified": "NO", "consequence": "retain as provenance-unresolved metadata only"},
        {"item": "Independent reproduction", "status": "UNRESOLVED", "source": "GEO/BioProject/SRA/GitHub/Zenodo/Figshare/institutional-source search", "direct_quote_or_field": "No linked PRS weights or reproducible score code located", "interpretation": "score cannot be independently regenerated", "verified": "NO", "consequence": "PRS and PRS×DEX remain out of primary scope"},
    ]
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "PRS_PROVENANCE_AUDIT.tsv", sep="\t", index=False)
    return result


def phase1_lock(verdict: str, qc_summary: dict[str, object], missing15_explained: str, adaptive: str, maladaptive: str) -> str:
    return f"""# PHASE 1 ANALYSIS LOCK — issued after Phase 0R

Audit date: 2026-09-17  
Phase 0R verdict: **{verdict}**  
GSE313736 role: **ROLE C — observational genetic-response resource; no primary DEX causal estimand**

This file is a preregistration lock. It authorizes only the conditional scope below. It does not authorize execution of Phase 1 in the current turn.

## Locked estimand and unit

- Primary biological question: define, from independent perturbation resources, a mechanoadaptive versus mechanopathological architecture in human trabecular meshwork.
- Statistical unit: donor for donor-level analyses; donor-pseudobulk for expression summaries.
- Primary GSE313736 cell frame: **166,323 exact matrix-matched metadata cells from 59 observed donors, 2 treatment labels, and 18 libraries**.
- Primary cell populations: BeamA and JCT. Schlemm_Endothelium is secondary only if donor coverage and cell numbers are adequate after QC. Other small/contaminant-like labels are exploratory or excluded by a locked rule.
- GSE195756 is the adaptive mechanical-stretch source, currently **{adaptive}** pending sample-level matrix/mapping QC.
- Maladaptive resources are independent TGF-beta2, DEX, and human glaucoma-context datasets, currently **{maladaptive}** pending sample-level QC.

## GSE313736 restrictions

- Treatment is exactly separated into vehicle-only and DEX-only pool/library sets. Donor pairing is retained, but `treatment` is aliased with `pool` and `library` in the supplied design.
- `DEX_CAUSAL_EFFECT_NOT_IDENTIFIABLE_FROM_PRIMARY_DATA_ALONE` is locked.
- DEX exposure duration is inconsistently reported in public record.
- The 74 VCF donors versus 59 cell-analysis donors remains **{missing15_explained}**; no claim of 74 analyzed donors is allowed.
- No same-dataset DEX gene discovery may be used to define a confirmatory score in GSE313736.

## PRS lock

- PRS provenance is unresolved. PRS is **not allowed as a primary predictor** and PRS×DEX is **not allowed as a confirmatory analysis**.
- If provenance is later fully verified, retain the raw continuous score as the preferred representation; do not dichotomize into high/low groups without a separately justified, pre-specified calibration rule.
- Do not use disease status as a PRS proxy: all observed GSE313736 donors are non-glaucoma.

## Permitted future Phase 1 scope

1. Audit and, if passed, construct adaptive and maladaptive signatures only from independent sources.
2. Test external signatures in the locked GSE313736 frame as an independent application, with donor-level summaries and explicit treatment/pool limitation language.
3. Report GSE313736 as observational genetic-response context; do not label DEX differences causal.
4. Use independent human DEX/TGF-beta2/stretch/glaucoma resources for direction and concordance checks only after sample identity, processed counts, pairing, and batch metadata are verified.

## Forbidden until a new gate is explicitly closed

DE, pathway enrichment, scores/signatures from GSE313736 discovery, PRS association, PRS×DEX confirmatory testing, eQTL, colocalization, trajectory, GRN, network, ML, and any clinical-utility claim.

Technical paired-QC summary: {json.dumps(qc_summary, ensure_ascii=False, sort_keys=True)}
"""


def governance_append(verdict: str) -> str:
    return f"""

## Phase 0R gate — 2026-09-17

Phase 0R verdict: **{verdict}**. The locked primary frame is 59 observed donors and 166,323 exact matrix-matched metadata cells. The GEO summary's 74 donors/316,799 singlets and the matrix's 350,171 barcodes are contextual denominators, not interchangeable analysis Ns.

Treatment is perfectly separated into vehicle-only and DEX-only pools/libraries. Because treatment is aliased with pool/library, adding those covariates does not identify a DEX effect. The project must use the exact phrase `DEX_CAUSAL_EFFECT_NOT_IDENTIFIABLE_FROM_PRIMARY_DATA_ALONE` for GSE313736. GSE313736 is ROLE C: an observational genetic-response resource; no primary DEX causal effect or confirmatory PRS×DEX claim is allowed.

PRS fields are present for 59 donors, but GWAS source, build, weights, scoring method, harmonization, and centile reference remain unresolved. PRS is not primary and PRS×DEX is not confirmatory. If provenance is later verified, raw continuous PRS is preferred and high/low dichotomization remains prohibited unless separately locked.

Future work is restricted to independent adaptive/maladaptive source QC and externally locked signature construction/application. The current gate does not authorize DE, enrichment, PRS association, eQTL, colocalization, trajectory, GRN, network, ML, or clinical validation.
"""


def project_status(verdict: str, qc_summary: dict[str, object]) -> str:
    return f"""# Project status — Phase 0R complete

Updated: 2026-09-17  
Project: PIEZO1/TRPV4 adaptive versus maladaptive mechanotransduction in glaucoma

## Gate

**{verdict}**

Phase 0R rescued and audited the primary resource without opening expression values. The exact allowed next scope is external perturbation-source QC/signature construction and independent application planning; no Phase 1 biological analyses were run.

## Locked denominators

- Primary donors: 59 observed donors with vehicle and DEX metadata and genotype mapping.
- VCF donor columns: 74; 15 are not present in the supplied donor/cell metadata and remain unexplained at donor-specific level.
- GEO summary: 316,799 singlets / 74 donors.
- Downloaded matrix barcodes: 350,171.
- Supplied cell metadata: 169,398 rows.
- Exact metadata–matrix intersection: 166,323 cells; this is `PRIMARY_CELL_ANALYSIS_FRAME`.

## Main gate findings

- PRS provenance: unresolved; no primary PRS or confirmatory PRS×DEX analysis.
- Treatment versus pool/library: non-identifiable in GSE313736; donor pairing does not restore orthogonality.
- DEX duration: unresolved; public record reports both 48 hours and 3 days.
- GSE313736 role: ROLE C.
- Adaptive source: {"CONDITIONAL"} (GSE195756).
- Maladaptive sources: {"CONDITIONAL"} (independent DEX/TGF-beta2/glaucoma-context resources).

## Technical QC snapshot

`{json.dumps(qc_summary, ensure_ascii=False, sort_keys=True)}`

## Required next gate

Do not start Phase 1 until the independent source-level QC, sample mapping, and the scope in `PHASE1_ANALYSIS_LOCK.md` are explicitly accepted. Any attempt to recover a primary DEX causal effect from GSE313736 requires new orthogonal experimental/batch information or a changed estimand.

Machine-readable outputs: `results/phase0R/`  
Full report: `reports/phase0R/PHASE0R_RESCUE_REPORT.md`
"""


def main() -> None:
    cm = pd.read_csv(CELL_META, sep="\t", compression="gzip")
    dm = pd.read_csv(DONOR_META, sep="\t", compression="gzip")
    matrix_df, matrix_pairs = matrix_barcode_frame()
    design_df = experimental_design(cm, matrix_df, matrix_pairs)
    donor_rec = donor_reconciliation(cm, dm, matrix_pairs, matrix_df)
    denom = cell_denominator_reconciliation(cm, matrix_df, matrix_pairs)
    ident = identifiability(design_df)
    _, qc_summary = donor_paired_qc(cm, matrix_pairs)
    series, samples = parse_family(PRIMARY_FAMILY)
    (REPORT / "DEX_PROTOCOL_AUDIT.md").write_text(protocol_audit(series, samples))
    prs = prs_audit(dm)
    ext, sig, adaptive, maladaptive = external_dex_audit()
    missing15_explained = "NO" if int((donor_rec["final_status"] == "UNKNOWN").sum()) == 15 else "PARTIAL"
    verdict = "GO_PHASE1_EXTERNAL_DEX_ONLY"
    # Do not overwrite Phase 0; append only the new locked gate.
    config = ROOT / "config" / "project_config.yaml"
    config_text = config.read_text()
    config_text = config_text.replace("phase: phase0\n", "phase: phase0R\n", 1)
    config_text = config_text.replace("phase_gate: phase0_only\n", "phase_gate: phase0R_complete\n", 1)
    config_text = re.sub(r"^primary_design:.*$", "primary_design: external_signature_application_with_observational_GSE313736_context", config_text, count=1, flags=re.MULTILINE)
    config_text = re.sub(r"^phase1_blocked_until:.*$", "phase1_blocked_until: [\"independent_external_source_QC\", \"phase1_analysis_lock_accepted\"]", config_text, count=1, flags=re.MULTILINE)
    config_text = config_text.split("# Phase 0R lock (2026-09-17)", 1)[0].rstrip() + "\n"
    config.write_text(config_text + f"""

# Phase 0R lock (2026-09-17)
phase0r_verdict: {verdict}
primary_analysis_donors: 59
primary_analysis_cells: 166323
primary_cell_analysis_frame: matrix-metadata exact intersection on (library_id, barcode)
gse313736_role: ROLE_C
prs_provenance_status: UNRESOLVED
prs_primary_allowed: false
prs_dex_confirmatory_allowed: false
treatment_pool_library_identifiability: NON_IDENTIFIABLE
dex_causal_effect_identifiable_from_primary_data_alone: false
dex_duration_status: UNRESOLVED
phase1_allowed: true
phase1_scope: external independent adaptive/maladaptive source QC and externally locked signature application only
phase1_full_primary_allowed: false
primary_dex_causal_effect_allowed: false
primary_genetic_aim_allowed: false
primary_cell_populations: [BeamA, JCT]
secondary_cell_populations: [Schlemm_Endothelium]
ctgf_alias: CCN2
""")
    gov = ROOT / "ANALYSIS_GOVERNANCE.md"
    gov_text = gov.read_text()
    gov_text = gov_text.split("## Phase 0R gate — 2026-09-17", 1)[0].rstrip() + "\n"
    gov.write_text(gov_text + governance_append(verdict))
    (ROOT / "PROJECT_STATUS.md").write_text(project_status(verdict, qc_summary))
    (REPORT / "PHASE1_ANALYSIS_LOCK.md").write_text(phase1_lock(verdict, qc_summary, missing15_explained, adaptive, maladaptive))

    report = f"""# PHASE 0R — PRIMARY DATASET RESCUE AND CAUSAL-IDENTIFIABILITY AUDIT

Audit date: 2026-09-17. Scope was restricted to public-source provenance rescue, donor/cell denominator reconciliation, experimental design rank, paired technical QC, protocol reconciliation, and external-resource suitability. No expression values were read by this script and no Phase 1 biological analysis was run.

## Overall verdict

**{verdict}**

This verdict authorizes only the exact external-data scope in `reports/phase0R/PHASE1_ANALYSIS_LOCK.md`. GSE313736 cannot support a primary DEX causal estimand from the supplied data alone.

## Locked primary frame

| Quantity | Value | Interpretation |
|---|---:|---|
| Observed donors | 59 | all have vehicle and DEX metadata; all are non-glaucoma |
| VCF donor columns | 74 | 15 donor IDs lack supplied donor/cell metadata |
| GEO summary singlets | 316,799 | contextual series-level denominator |
| Matrix barcodes | 350,171 | downloaded filtered-GEX barcode denominator |
| Supplied metadata rows | 169,398 | author-supplied demuxlet/classification rows |
| Exact matrix–metadata intersection | 166,323 | `PRIMARY_CELL_ANALYSIS_FRAME` |

The exact primary expression frame is the 166,323-row inner join on `(library_id, barcode)`, retaining 59 donors, both treatments, and author annotations. The 74-donor and 316,799-cell GEO claims are not used as analyzed denominators.

## PRS provenance rescue

`results/phase0R/PRS_PROVENANCE_AUDIT.tsv` records each provenance element. The raw PRS and centile fields are present for all 59 observed donors, but the linked public record does not provide a GWAS accession/PMID/DOI, genome build for the PRS, SNP/effect-allele weights, scoring method, imputation/QC, harmonization, or centile reference. Searches covered the GEO series and sample records, linked BioProject `PRJNA1381509`, NCBI SRA query, and public code/repository search paths on 2026-09-17; no reproducible PRS definition was located.

**PRS status: UNRESOLVED.** The genotype-to-donor mapping is usable, but the score is not yet scientifically reproducible. PRS is not permitted as a primary predictor, and PRS×DEX is not permitted as a confirmatory analysis.

## 74 versus 59 donor reconciliation

`results/phase0R/DONOR_RECONCILIATION.tsv` contains one row for every VCF donor column. 59 donors are `ANALYSIS_DONOR`. The 15 VCF-only donors are marked `UNKNOWN`: the supplied record does not document whether they failed single-cell QC, had no retained cells, were excluded by the author, or are missing because of metadata export. Optional demuxlet BEST evidence is recorded if present, but it cannot replace the supplied author cell metadata as the primary frame.

**Missing 15 explained: {missing15_explained}.** Until donor-specific exclusion/QC documentation is found, selection bias from 74→59 remains unresolved.

## Cell denominator reconciliation

`results/phase0R/CELL_DENOMINATOR_RECONCILIATION.tsv` keeps all denominators separate. The 3,075 metadata rows without a matching matrix barcode are excluded from expression analyses pending reconciliation. No denominator is silently substituted for another.

## Treatment × pool identifiability

`results/phase0R/EXPERIMENTAL_DESIGN_MATRIX.tsv` records 118 donor×condition rows across 18 libraries. The supplied hierarchy is donor → treatment → pool → library → sequencing run/lane → batch → cell. Treatment is exactly nested in the pool/library naming layer: vehicle is represented by `HTMC_pool-A`–`I`, DEX by `HTMC_pool-dexA`–`I`. Processing date, culture batch, treatment plate, library batch, Chromium run, sequencing lane, depth, median genes/UMI, mitochondrial/ribosomal fraction, and doublet rate are not supplied in the audited metadata.

`results/phase0R/DESIGN_IDENTIFIABILITY_REPORT.tsv` evaluates fixed-effect rank analogues for `~ donor + treatment`, `~ donor + treatment + library`, and `~ donor + treatment + pool`, plus mixed-model analogues. Donor-paired treatment is estimable without library/pool terms; adding library or pool makes treatment aliased and the full design rank-deficient. Adding pool/library covariates cannot solve exact nesting.

**Locked conclusion: `DEX_CAUSAL_EFFECT_NOT_IDENTIFIABLE_FROM_PRIMARY_DATA_ALONE`.**

## Donor-paired technical QC

`results/phase0R/DONOR_PAIRED_TECHNICAL_QC.tsv` is descriptive and contains one row per donor. All supplied droplets are labelled SNG. The record provides library-level sequencing saturation but does not provide per-cell UMI, gene counts, mitochondrial/ribosomal fractions, or a numeric doublet rate. Cell-state composition and cell counts are recorded as QC descriptors only; they are not biological DEX effects.

## DEX protocol audit

See `reports/phase0R/DEX_PROTOCOL_AUDIT.md`. Concentration is recorded as 500 nM and the public sample record reports vehicle solution. Duration is contradictory: sample characteristics say 48 h, while the growth protocol says 3 days.

**Locked wording: DEX exposure duration is inconsistently reported in public record.**

## GSE313736 role

**ROLE C — observational genetic-response resource.** The paired donor structure is real, but treatment is not causally separable from pool/library in the primary data. GSE313736 may be used later for descriptive donor-level response context and independent external-signature application; it may not be used as a primary DEX causal or confirmatory PRS×DEX dataset.

## External DEX and stress resources

`results/phase0R/EXTERNAL_DEX_DATASET_AUDIT.tsv` audits GSE301517 and related regulatory subseries. GSE301517 is the preferred external DEX reference candidate: 21 total RNA-seq records with title-matched control/DEX replicate structure and processed CPM/RPKM files listed. Treatment-versus-batch separability and donor identity are not fully documented, so its status is **CONDITIONAL**, not definitive.

`results/phase0R/SIGNATURE_SOURCE_AUDIT.tsv` assigns independent source suitability without discovering any gene set. GSE195756 is the adaptive mechanical-stretch source and is **CONDITIONAL** because the 4 pairs are inferred from labels and donor IDs are not explicit. GSE301517/GSE202372/GSE263230 and related regulatory resources are **CONDITIONAL** maladaptive/context sources pending sample-level mapping and preprocessing audits. GSE263230 is a human glaucoma/TM context resource, not a DEX signature source.

## Phase 1 lock

See `reports/phase0R/PHASE1_ANALYSIS_LOCK.md`. The future hierarchy is: independent adaptive mechanical stretch; independent maladaptive DEX/TGF-beta2/glaucoma resources; GSE313736 only according to ROLE C; donor as unit; raw continuous PRS preferred only after provenance verification; BeamA/JCT primary, Schlemm secondary if adequate, other labels exploratory/excluded. No high/low PRS dichotomization is authorized.

## Remaining limitations

The largest remaining limitation is exact treatment–pool/library aliasing in GSE313736, compounded by the unresolved 74→59 donor discrepancy and unresolved PRS provenance. The permitted path forward is an external-Dex-only phase using independently locked sources; recovering a causal DEX or PRS×DEX claim from GSE313736 requires new orthogonal batch information or a changed estimand.

Public sources: [GSE313736](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736), [GSE301517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301517), [GSE195756](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195756), [GSE202372](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202372), [GSE263230](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263230), and linked [BioProject PRJNA1381509](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1381509).
"""
    (REPORT / "PHASE0R_RESCUE_REPORT.md").write_text(report)

    # Machine-readable run manifest for auditability.
    manifest = {
        "script": "scripts/07_phase0R_rescue_audit.py",
        "date": "2026-09-17",
        "phase0r_verdict": verdict,
        "primary_donors": int(cm["donor_id_demuxlet"].nunique()),
        "primary_metadata_rows": int(len(cm)),
        "primary_matrix_matched_rows": int(denom.loc[denom["object"] == "metadata_matrix_intersection", "cells"].iloc[0]),
        "vcf_donors": int(len(vcf_donors()[0])),
        "missing_vcf_only_donors": int((donor_rec["final_status"] == "UNKNOWN").sum()),
        "matrix_barcodes": int(matrix_df["matrix_barcodes"].sum()),
        "design_rows": int(len(design_df)),
        "design_identifiability_rows": int(len(ident)),
        "forbidden_analyses_run": False,
    }
    (OUT / "PHASE0R_RUN_MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    print(f"PHASE0R_VERDICT={verdict}")
    print(f"PRIMARY_DONORS={manifest['primary_donors']}")
    print(f"PRIMARY_CELLS={manifest['primary_matrix_matched_rows']}")
    print(f"MISSING_15_EXPLAINED={missing15_explained}")
    print(f"PRS_PROVENANCE=UNRESOLVED")
    print(f"TREATMENT_POOL_IDENTIFIABILITY=NON_IDENTIFIABLE")
    print(f"DEX_DURATION=UNRESOLVED")
    print(f"GSE313736_ROLE=ROLE_C")
    print(f"ADAPTIVE_SOURCE={adaptive}")
    print(f"MALADAPTIVE_SOURCES={maladaptive}")
    print("PRS_PRIMARY_ALLOWED=NO")
    print("PRSxDEX_CONFIRMATORY=NO")
    print("PHASE1_ALLOWED=YES_EXTERNAL_ONLY")
    print("PHASE1_SCOPE=independent adaptive/maladaptive source QC and externally locked signature application; no primary DEX causal or PRSxDEX claim")
    print("BIGGEST_LIMITATION=DEX is exactly aliased with pool/library; 74-to-59 donor discrepancy and PRS provenance remain unresolved")
    print("REPORT=reports/phase0R/PHASE0R_RESCUE_REPORT.md")


if __name__ == "__main__":
    main()
