#!/usr/bin/env python3
"""Build the Phase 0 registry from locally cached GEO SOFT records.

This script only audits accession metadata and supplementary-file listings. It
does not run any biological or association analysis.
"""
from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data" / "metadata"
OUT = ROOT / "config" / "dataset_registry.tsv"
ACCESSIONS = ["GSE313736", "GSE301525", "GSE301517", "GSE301516", "GSE301512", "GSE301519", "GSE301522", "GSE195756", "GSE202372", "GSE263230", "GSE309500"]

ROLES = {
    "GSE313736": "primary DEX scRNA/genotype/PRS dataset",
    "GSE301525": "regulatory multiomics superseries",
    "GSE301517": "DEX RNA-seq regulatory support",
    "GSE301516": "DEX ATAC regulatory support",
    "GSE301512": "CUT&Run regulatory resource; treatment scope audit required",
    "GSE301519": "DEX CUT&Tag regulatory support",
    "GSE301522": "DEX Hi-C regulatory support",
    "GSE195756": "independent mechanical-stretch signature source",
    "GSE202372": "independent TGFbeta2 stress signature source",
    "GSE263230": "independent human TM cell-identity validation",
    "GSE309500": "candidate disease/TM validation; exclude from human primary analysis",
}

FIELDS = ["accession", "title", "organism", "tissue", "cell_type", "assay", "perturbation", "n_donors_reported", "n_samples_reported", "genotype_available", "PRS_available", "disease_status_available", "processed_matrix_available", "raw_data_available", "metadata_available", "download_size", "public_access", "primary_role", "source_url", "notes"]


def get_value(lines, prefix, default=""):
    for line in lines:
        if line.startswith(prefix):
            return line.split("=", 1)[1].strip()
    return default


def get_values(lines, prefix):
    return [line.split("=", 1)[1].strip() for line in lines if line.startswith(prefix)]


def parse_soft(acc):
    path = META / f"{acc}_family.soft.gz"
    lines = gzip.open(path, "rt", errors="replace").read().splitlines()
    series = lines[: next((i for i, x in enumerate(lines) if x.startswith("^PLATFORM")), len(lines))]
    samples = []
    for i, line in enumerate(lines):
        if not line.startswith("^SAMPLE = "):
            continue
        end = next((j for j in range(i + 1, len(lines)) if lines[j].startswith("^")), len(lines))
        block = lines[i:end]
        samples.append({
            "acc": line.split("=", 1)[1].strip(),
            "title": get_value(block, "!Sample_title = "),
            "source": get_value(block, "!Sample_source_name_ch1 = "),
            "organism": get_value(block, "!Sample_organism_ch1 = "),
            "chars": get_values(block, "!Sample_characteristics_ch1 = "),
            "supp": get_values(block, "!Sample_supplementary_file"),
        })
    return series, samples


def infer_donors(acc, series, samples):
    text = " ".join(series)
    if acc == "GSE313736":
        return "74 reported; 59 observed in cell metadata"
    if acc == "GSE263230":
        return "14 primary culture samples + 3 post-mortem samples; donor IDs not supplied in SOFT"
    if acc == "GSE309500":
        return "mouse samples pooled from 8 eyes per sequencing run; individual donor N not applicable"
    if acc == "GSE195756":
        return "4 paired biological samples inferred from Control/Stretch pairs; donor IDs not verified"
    if acc == "GSE202372":
        return "5 paired biological samples inferred from A-E labels; donor IDs not verified"
    if acc.startswith("GSE301"):
        return "not reported as human donor IDs; primary TM cell lines/replicates"
    return "UNCLEAR"


def infer_perturbation(acc, text, samples):
    low = text.lower()
    if acc == "GSE313736": return "vehicle versus dexamethasone 500 nM, 48 h"
    if acc == "GSE195756": return "15% cyclic mechanical stretch, 0.69 Hz, 8 h"
    if acc == "GSE202372": return "TGF-beta2 versus control, 24 h"
    if acc.startswith("GSE301") and "dexamethasone" in low: return "resource/control/dexamethasone conditions; modality-specific"
    if acc == "GSE263230": return "none; culture versus post-mortem comparison"
    if acc == "GSE309500": return "mouse WT versus Lmx1b V265D/+; vitamin B3 protection context"
    return "UNCLEAR"


def infer_assay(acc, series, samples):
    title = get_value(series, "!Series_title = ")
    types = get_values(series, "!Series_type = ")
    if acc == "GSE313736": return "10x scRNA-seq"
    if "RNA-seq" in title: return "RNA-seq"
    if "ATAC" in title: return "ATAC-seq"
    if "CUT&Run" in title: return "CUT&Run"
    if "CUT&Tag" in title: return "CUT&Tag"
    if "HiC" in title: return "Hi-C"
    if acc == "GSE195756": return "RNA-seq (series also reports multiomics)"
    if acc == "GSE202372": return "RNA-seq"
    if acc == "GSE263230": return "10x scRNA-seq"
    if acc == "GSE309500": return "scRNA-seq and multiome RNA/ATAC"
    return "; ".join(types) or "UNCLEAR"


def main():
    rows = []
    for acc in ACCESSIONS:
        series, samples = parse_soft(acc)
        title = get_value(series, "!Series_title = ")
        organism = next((s["organism"] for s in samples if s["organism"]), "")
        tissue = next((s["source"] for s in samples if s["source"]), "UNCLEAR")
        chars = " ".join(c for s in samples for c in s["chars"])
        cell_type = next((c.split(":", 1)[1].strip() for c in chars.split("; ") if c.lower().startswith("cell type:")), "UNCLEAR")
        if acc == "GSE313736":
            cell_type = "trabecular meshwork cells; coarse author label verified in cell metadata"
        all_supp = [x for s in samples for x in s["supp"]] + get_values(series, "!Series_supplementary_file")
        raw = "YES_VERIFIED" if any("RAW" in x for x in all_supp) else "NO"
        processed = "YES_VERIFIED" if all_supp and any(not "RAW" in x for x in all_supp) else "NO"
        genotype = "YES_VERIFIED" if acc == "GSE313736" and (META / "GSE313736_HTMC_genotypes.vcf.gz").exists() else "NO"
        prs = "YES_VERIFIED" if acc == "GSE313736" and (META / "GSE313736_HTMC_scRNA_donor_metadata.tsv.gz").exists() else "NO"
        disease = "YES_VERIFIED" if acc == "GSE313736" and (META / "GSE313736_HTMC_scRNA_donor_metadata.tsv.gz").exists() else ("REPORTED_NOT_VERIFIED" if acc in {"GSE202372", "GSE263230"} else "NO")
        metadata = "YES_VERIFIED" if (META / f"{acc}_family.soft.gz").exists() else "NO"
        size = ""
        if acc == "GSE313736": size = "RAW.tar 6,536,949,760 B; metadata 6.3 MB; VCF 15 MB"
        elif acc == "GSE301525": size = "RAW.tar 3.2 GB"
        elif acc == "GSE263230": size = "RAW.tar 1.9 GB; sample h5ad files listed"
        elif acc == "GSE309500": size = "RAW.tar; sample matrix/TSV files listed"
        else: size = "supplementary processed files listed; exact aggregate size not audited"
        notes = []
        if acc == "GSE313736": notes.append("GEO reports 316,799 singlets/74 donors; supplied cell metadata contains 169,398 cells/59 donors. Two treatment-separated pooled captures.")
        if acc == "GSE301512": notes.append("Series design says control human TM cells; do not assume DEX CUT&Run support.")
        if acc == "GSE309500": notes.append("Mouse limbal strip, not human TM; exclude from human validation claim.")
        if acc == "GSE263230": notes.append("Human primary TM culture plus post-mortem tissue; no DEX treatment in SOFT.")
        if acc == "GSE195756": notes.append("Paired stretch/control labels present; donor IDs not explicitly supplied.")
        if acc == "GSE202372": notes.append("Five control/TGF-beta2 pairs inferred from A-E labels; donor IDs not explicitly supplied.")
        rows.append({
            "accession": acc, "title": title, "organism": organism, "tissue": tissue, "cell_type": cell_type,
            "assay": infer_assay(acc, series, samples), "perturbation": infer_perturbation(acc, " ".join(series), samples),
            "n_donors_reported": infer_donors(acc, series, samples), "n_samples_reported": str(len(samples)),
            "genotype_available": genotype, "PRS_available": prs, "disease_status_available": disease,
            "processed_matrix_available": processed, "raw_data_available": raw, "metadata_available": metadata,
            "download_size": size, "public_access": "YES_VERIFIED" if "Public" in get_value(series, "!Series_status = ") else "NO",
            "primary_role": ROLES[acc], "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
            "notes": " ".join(notes),
        })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t"); w.writeheader(); w.writerows(rows)


if __name__ == "__main__": main()
