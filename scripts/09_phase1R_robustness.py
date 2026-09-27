#!/usr/bin/env python3
"""Phase 1R: robustness, null calibration, and scope locking.

This script consumes the immutable Phase 1 external-signature outputs and
does not discover or redefine signatures from GSE313736.  It re-reads the
exact metadata--matrix intersection once to support fixed-size sensitivity,
expression/detection-matched nulls, anchor removal, and donor bootstrap.
GSE313736 remains ROLE C; every treatment comparison is descriptive and
donor-paired, with no causal DEX or PRS interpretation.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import hypergeom, spearmanr


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase1R"
REPORT = ROOT / "reports" / "phase1R"
LOGS = ROOT / "logs"
OUT.mkdir(parents=True, exist_ok=True)
REPORT.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)
RUN_LOG = LOGS / "phase1R_run.log"

SEED = 20260918
FDR_CUTOFF = 0.05
SIZES = [25, 50, 75, 100, 150, 200]
SIZE_LABELS = [f"N={x}" for x in SIZES] + ["RANKED_FULL"]
NULL_N = 10000
NULL_BATCH = 200
BOOTSTRAPS = 2000

P1_PATH = ROOT / "scripts" / "08_phase1_external_signatures.py"
SPEC = importlib.util.spec_from_file_location("phase1_locked", P1_PATH)
P1 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(P1)


def log(message: str) -> None:
    with RUN_LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}\n")


def write_tsv(frame: pd.DataFrame, path: Path, columns: list[str] | None = None) -> None:
    if columns is not None:
        frame = frame.reindex(columns=columns)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def fnum(x: object, default: float = np.nan) -> float:
    try:
        v = float(x)
        return v if np.isfinite(v) else default
    except Exception:
        return default


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan)
    ok = np.isfinite(values)
    if not ok.any():
        return out
    p = np.clip(values[ok], 0, 1)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.minimum(q, 1.0)
    out[ok] = restored
    return out


def clean_frame(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    if "gene" in frame:
        frame["gene"] = frame["gene"].map(P1.gene)
    for col in ["logFC", "statistic", "FDR", "PValue", "rank"]:
        if col in frame:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame[frame["gene"].ne("")].copy()


def load_locked_external() -> dict[str, pd.DataFrame]:
    p = ROOT / "results" / "phase1"
    frames = {
        "Adaptive": clean_frame(p / "ADAPTIVE_SIGNATURE.tsv"),
        "Fibrotic": clean_frame(p / "MALADAPTIVE_FIBROTIC.tsv"),
        "Steroid": clean_frame(p / "MALADAPTIVE_STEROID.tsv"),
        "Disease": clean_frame(p / "MALADAPTIVE_DISEASE.tsv"),
    }
    log("Loaded immutable Phase 1 external tables")
    return frames


def eligible(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame[frame["FDR"].lt(FDR_CUTOFF)].copy()
    sort_cols = ["FDR"]
    ascending = [True]
    if "support_n" in result.columns:
        sort_cols = ["support_n", "FDR"]
        ascending = [False, True]
    if "PValue" in result.columns:
        sort_cols.append("PValue")
        ascending.append(True)
    sort_cols.append("gene")
    ascending.append(True)
    return result.sort_values(sort_cols, ascending=ascending, na_position="last").drop_duplicates("gene").reset_index(drop=True)


def consensus_full(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Reconstruct the untruncated TGFβ2+DEX consensus, excluding disease."""
    contexts = {"TGFb2": eligible(frames["Fibrotic"]), "DEX": eligible(frames["Steroid"])}
    records: dict[str, list[tuple[str, float, float, float]]] = defaultdict(list)
    for context, frame in contexts.items():
        for _, row in frame.iterrows():
            records[str(row["gene"])].append((context, fnum(row["logFC"]), fnum(row["statistic"]), fnum(row["FDR"])))
    rows = []
    for g, vals in records.items():
        signs = {1 if x[1] > 0 else -1 for x in vals if x[1] != 0}
        if len(vals) < 2 or len(signs) != 1:
            continue
        rows.append({
            "gene": g,
            "logFC": float(np.mean([x[1] for x in vals])),
            "statistic": float(np.mean([x[2] for x in vals])),
            "FDR": float(np.max([x[3] for x in vals])),
            "support_n": 2,
            "supporting_contexts": "GSE202372;GSE301517",
        })
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return result.sort_values(["support_n", "FDR", "gene"], ascending=[False, True, True]).reset_index(drop=True)


def make_variants(frame: pd.DataFrame) -> dict[str, pd.DataFrame]:
    sig = eligible(frame)
    out = {}
    for n in SIZES:
        out[f"N={n}"] = sig.head(min(n, len(sig))).copy()
    out["RANKED_FULL"] = sig.copy()
    return out


def audit_signature_construction(frames: dict[str, pd.DataFrame], core: pd.DataFrame) -> pd.DataFrame:
    rows = []
    definitions = [
        ("Adaptive-Core", frames["Adaptive"], 100, "FDR<0.05 followed by ranked top-N cap", "external author statistics; no biological filter"),
        ("Maladaptive-Fibrotic", frames["Fibrotic"], 100, "FDR<0.05 followed by ranked top-N cap", "TGFβ2 author DE list; no biological filter"),
        ("Maladaptive-Steroid", frames["Steroid"], 100, "FDR<0.05 followed by ranked top-N cap", "DEX processed CPM Welch statistic; no biological filter"),
    ]
    for name, frame, cap, rule, comment in definitions:
        e = eligible(frame)
        rows.append({
            "signature": name,
            "total_eligible_genes_before_cap": len(e),
            "selection_rule": rule,
            "top_n_cap": cap,
            "number_after_statistical_filter": len(e),
            "number_after_biological_filter": len(e),
            "final_n": min(cap, len(e)),
            "cap_binding": "YES" if len(e) > cap else "NO",
            "comments": f"source_rows_before_filter={len(frame)}; {comment}",
        })
    fib = eligible(frames["Fibrotic"])
    dex = eligible(frames["Steroid"])
    rows.append({
        "signature": "Maladaptive-Core",
        "total_eligible_genes_before_cap": len(core),
        "selection_rule": "concordant FDR<0.05 support in both TGFβ2 and DEX, then ranked top-N cap",
        "top_n_cap": 100,
        "number_after_statistical_filter": len(core),
        "number_after_biological_filter": len(core),
        "final_n": min(100, len(core)),
        "cap_binding": "YES" if len(core) > 100 else "NO",
        "comments": f"TGFβ2 eligible={len(fib)}; DEX eligible={len(dex)}; disease context excluded from core membership",
    })
    disease = eligible(frames["Disease"])
    rows.append({
        "signature": "Maladaptive-Disease-context",
        "total_eligible_genes_before_cap": len(disease),
        "selection_rule": "FDR<0.05 followed by ranked top-N cap; contextual only",
        "top_n_cap": 100,
        "number_after_statistical_filter": len(disease),
        "number_after_biological_filter": len(disease),
        "final_n": min(100, len(disease)),
        "cap_binding": "NO",
        "comments": f"source_rows_before_filter={len(frames['Disease'])}; GSE263230 is descriptive disease-context evidence and does not determine Maladaptive-Core membership",
    })
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "SIGNATURE_CONSTRUCTION_AUDIT.tsv")
    return result


def primary_full_background() -> dict[str, object]:
    """Aggregate all exact-frame feature rows to donor×treatment×cell-state groups.

    Sparse matrix multiplication is used per library so the same pseudobulk
    and within-cell-state z-scoring procedure can be applied to real and
    matched-null gene sets without opening the discovery loop.
    """
    log("Phase 1R primary pass: reading exact GSE313736 frame for null calibration and sensitivity only")
    cm = P1.primary_metadata()
    group_cols = ["donor_id_demuxlet", "treatment", "celltype"]
    cm = cm[cm["key"].notna()].copy()
    groups = cm[group_cols + ["analysis_population"]].drop_duplicates().sort_values(group_cols).reset_index(drop=True)
    groups["group_id"] = np.arange(len(groups))
    group_map = {tuple(row[c] for c in group_cols): int(row["group_id"]) for _, row in groups.iterrows()}
    ngroups = len(groups)

    infos = P1.matrix_library_info()
    if not infos:
        raise RuntimeError("No GSE313736 Matrix Market libraries found")
    feature_names: list[str] = []
    with gzip.open(infos[0]["features"], "rt", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            feature_names.append(P1.gene(fields[1] if len(fields) >= 2 else fields[0]))
    gene_list = sorted({g for g in feature_names if g})
    gene_to_id = {g: i for i, g in enumerate(gene_list)}
    n_genes = len(gene_list)
    cache_path = OUT / "PRIMARY_FULL_CACHE.npz"
    if cache_path.exists():
        cache = np.load(cache_path, allow_pickle=False)
        cached_genes = cache["gene_list"].astype(str).tolist()
        if cached_genes != gene_list or cache["counts"].shape[1] != ngroups:
            raise RuntimeError("PRIMARY_FULL_CACHE.npz does not match the current feature or group frame")
        counts = cache["counts"].astype(np.float32, copy=False)
        detected = cache["detected"].astype(np.float64, copy=False)
        total_counts = cache["total_counts"].astype(np.float64, copy=False)
        cell_counts = cache["cell_counts"].astype(np.int64, copy=False)
        library_stats = json.loads(str(cache["library_stats_json"]))
        log("Reusing verified Phase 1R full-background cache")
    else:
        counts = np.zeros((n_genes, ngroups), dtype=np.float32)
        detected = np.zeros(n_genes, dtype=np.float64)
        total_counts = np.zeros(ngroups, dtype=np.float64)
        cell_counts = np.zeros(ngroups, dtype=np.int64)
        library_stats = []

    for info in ([] if cache_path.exists() else infos):
        lib = str(info["library_id"])
        md = cm[cm["library_id"].eq(lib)]
        if md.empty:
            continue
        with gzip.open(info["barcodes"], "rt", errors="replace") as fh:
            barcodes = [x.strip() for x in fh if x.strip()]
        barcode_to_group = {
            str(row["barcode"]): group_map[(row["donor_id_demuxlet"], row["treatment"], row["celltype"])]
            for _, row in md.iterrows()
        }
        wanted_cols = []
        wanted_groups = []
        for idx, barcode in enumerate(barcodes):
            gid = barcode_to_group.get(barcode)
            if gid is not None:
                wanted_cols.append(idx)
                wanted_groups.append(gid)
        if not wanted_cols:
            continue
        wanted = np.asarray(wanted_cols, dtype=np.int64)
        gids = np.asarray(wanted_groups, dtype=np.int64)
        local_groups, inverse = np.unique(gids, return_inverse=True)
        with gzip.open(info["features"], "rt", errors="replace") as fh:
            local_features = []
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                local_features.append(P1.gene(fields[1] if len(fields) >= 2 else fields[0]))
        row_ids = np.asarray([gene_to_id.get(g, -1) for g in local_features], dtype=np.int64)
        valid_rows = row_ids >= 0

        log(f"Opening {lib}: wanted_cells={len(wanted_cols)} genes={n_genes}")
        matrix = mmread(str(info["matrix"]))
        if not sparse.issparse(matrix):
            matrix = sparse.coo_matrix(matrix)
        matrix = matrix.tocsr()
        sub = matrix[:, wanted]
        cell_totals = np.asarray(sub.sum(axis=0)).ravel()
        np.add.at(total_counts, gids, cell_totals)
        np.add.at(cell_counts, gids, 1)

        local_onehot = sparse.csr_matrix(
            (np.ones(len(wanted), dtype=np.float32), (np.arange(len(wanted)), inverse)),
            shape=(len(wanted), len(local_groups)),
        )
        local_counts = (sub @ local_onehot).toarray().astype(np.float32, copy=False)
        # ``np.add.at`` needs one-dimensional values for a one-dimensional
        # indexed target.  Add each local group column separately; this also
        # correctly sums any duplicate feature rows after symbol harmonization.
        for local_col, global_gid in enumerate(local_groups):
            np.add.at(counts[:, int(global_gid)], row_ids[valid_rows], local_counts[valid_rows, local_col])
        row_nnz = np.asarray(sub.getnnz(axis=1)).ravel()
        np.add.at(detected, row_ids[valid_rows], row_nnz[valid_rows])
        library_stats.append({
            "library_id": lib,
            "treatment": str(info["treatment"]),
            "matrix_columns": len(barcodes),
            "reconciled_cells_used": len(wanted_cols),
            "matrix_rows": int(matrix.shape[0]),
            "matrix_nnz": int(matrix.nnz),
        })
        del matrix, sub, local_onehot, local_counts, cell_totals

    groups["n_cells"] = cell_counts
    groups["total_umi"] = total_counts
    usable_total = np.where(total_counts > 0, total_counts, np.nan)
    log_cpm = np.log1p((counts.astype(np.float64) / usable_total[None, :]) * 1e6).astype(np.float32)
    mean_expr = np.nanmean(log_cpm, axis=1)
    variance = np.nanvar(log_cpm, axis=1)
    detection_rate = np.clip(detected / max(int(cell_counts.sum()), 1), 0, 1)
    background = pd.DataFrame({
        "gene": gene_list,
        "mean_logCPM": mean_expr,
        "detection_rate": detection_rate,
        "variance_logCPM": variance,
        "n_exact_frame_cells": int(cell_counts.sum()),
        "feature_row_duplicate_collapse": "summed_by_canonical_symbol",
    })
    write_tsv(background, OUT / "PRIMARY_BACKGROUND_SUMMARY.tsv")
    if not cache_path.exists():
        np.savez_compressed(cache_path, gene_list=np.asarray(gene_list, dtype="U"), counts=counts, detected=detected, total_counts=total_counts, cell_counts=cell_counts, library_stats_json=np.asarray(json.dumps(library_stats)))
    z_by_celltype: dict[str, np.ndarray] = {}
    group_indices: dict[str, np.ndarray] = {}
    pair_indices: dict[str, list[tuple[int, int, str]]] = {}
    pair_indices_local: dict[str, list[tuple[int, int, str]]] = {}
    for celltype in ["BeamA", "JCT", "Schlemm_Endothelium"]:
        idx = groups.index[groups["celltype"].eq(celltype)].to_numpy(dtype=int)
        group_indices[celltype] = idx
        vals = log_cpm[:, idx].astype(np.float64)
        mu = np.nanmean(vals, axis=1)
        sd = np.nanstd(vals, axis=1, ddof=1)
        z = (vals - mu[:, None]) / np.where(sd[:, None] > 0, sd[:, None], np.nan)
        z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
        z_by_celltype[celltype] = z
        loc = groups.iloc[idx]
        paired = []
        paired_local = []
        for donor, dsub in loc.groupby("donor_id_demuxlet"):
            v = dsub[dsub["treatment"].eq("vehicle")]
            d = dsub[dsub["treatment"].eq("DEX")]
            if not v.empty and not d.empty:
                vehicle_gid = int(v.iloc[0]["group_id"])
                dex_gid = int(d.iloc[0]["group_id"])
                paired.append((vehicle_gid, dex_gid, str(donor)))
                paired_local.append((int(np.flatnonzero(idx == vehicle_gid)[0]), int(np.flatnonzero(idx == dex_gid)[0]), str(donor)))
        pair_indices[celltype] = paired
        pair_indices_local[celltype] = paired_local
    log(f"Primary full background complete: genes={n_genes} groups={ngroups} exact_cells={int(cell_counts.sum())}")
    return {
        "gene_list": gene_list,
        "gene_to_id": gene_to_id,
        "groups": groups,
        "counts": counts,
        "log_cpm": log_cpm,
        "mean_expr": mean_expr,
        "variance": variance,
        "detection_rate": detection_rate,
        "z_by_celltype": z_by_celltype,
        "group_indices": group_indices,
        "pair_indices": pair_indices,
        "pair_indices_local": pair_indices_local,
        "library_stats": library_stats,
        "exact_cells": int(cell_counts.sum()),
    }


def score_frame(frame: pd.DataFrame, primary: dict[str, object], remove: set[str] | None = None) -> np.ndarray:
    remove = remove or set()
    gene_to_id = primary["gene_to_id"]
    selected = frame[~frame["gene"].isin(remove)].copy()
    selected = selected[selected["gene"].isin(gene_to_id)]
    if selected.empty:
        return np.full(len(primary["groups"]), np.nan)
    ids = np.asarray([gene_to_id[g] for g in selected["gene"]], dtype=int)
    signs = np.where(pd.to_numeric(selected["logFC"], errors="coerce").to_numpy(float) >= 0, 1.0, -1.0).astype(np.float32)
    return score_indices(ids, signs, primary)


def score_indices(ids: np.ndarray, signs: np.ndarray, primary: dict[str, object]) -> np.ndarray:
    n_groups = len(primary["groups"])
    scores = np.full(n_groups, np.nan, dtype=np.float32)
    for celltype, idx in primary["group_indices"].items():
        z = primary["z_by_celltype"][celltype]
        if len(ids) == 0:
            continue
        vals = z[ids]
        scores[idx] = np.nanmean(vals * signs[:, None], axis=0)
    return scores


def pair_summary(score: np.ndarray, primary: dict[str, object]) -> dict[str, dict[str, object]]:
    out = {}
    for celltype, pairs in primary["pair_indices"].items():
        if not pairs:
            out[celltype] = {"n": 0, "deltas": np.array([]), "median": np.nan, "positive_fraction": np.nan}
            continue
        delta = np.asarray([score[dex] - score[veh] for veh, dex, _ in pairs], dtype=float)
        finite = np.isfinite(delta)
        out[celltype] = {"n": int(finite.sum()), "deltas": delta[finite], "median": float(np.median(delta[finite])) if finite.any() else np.nan, "positive_fraction": float(np.mean(delta[finite] > 0)) if finite.any() else np.nan}
    return out


def score_gene_set(genes: set[str], primary: dict[str, object]) -> np.ndarray:
    ids = np.asarray([primary["gene_to_id"][g] for g in genes if g in primary["gene_to_id"]], dtype=int)
    return score_indices(ids, np.ones(len(ids), dtype=np.float32), primary)


def external_metrics(frame: pd.DataFrame, genes: set[str], target_name: str, rng: np.random.Generator, source_frame: pd.DataFrame | None = None) -> dict[str, object]:
    target = frame[frame["gene"].isin(genes)].copy()
    target = target[np.isfinite(pd.to_numeric(target["statistic"], errors="coerce"))]
    stats = pd.to_numeric(frame["statistic"], errors="coerce").dropna().to_numpy(float)
    vals = pd.to_numeric(target["statistic"], errors="coerce").to_numpy(float)
    if source_frame is not None and not source_frame.empty:
        source_dir = {str(r["gene"]): (1.0 if fnum(r["logFC"]) >= 0 else -1.0) for _, r in source_frame.iterrows()}
        vals = vals * np.asarray([source_dir.get(str(g), 1.0) for g in target["gene"]], dtype=float)
    target_sig = set(eligible(frame)["gene"])
    common = set(target["gene"])
    result = {"target_dataset": target_name, "overlap_n": len(common), "target_significant_n": len(target_sig), "jaccard": len(common) / max(len(set(genes) | target_sig), 1), "mean_signed_statistic": np.nan, "rank_z": np.nan, "empirical_p": np.nan}
    if len(vals) == 0 or len(stats) < 20:
        return result
    obs = float(np.mean(vals))
    size = len(vals)
    nrep = 500
    null = np.asarray([np.mean(rng.choice(stats, size=min(size, len(stats)), replace=False)) for _ in range(nrep)])
    sd = float(np.std(null, ddof=1))
    p = float((1 + np.sum(np.abs(null - np.mean(null)) >= abs(obs - np.mean(null)))) / (nrep + 1))
    result.update({"mean_signed_statistic": obs, "rank_z": (obs - float(np.mean(null))) / sd if sd > 0 else np.nan, "empirical_p": p})
    return result


def primary_rows_for_family(family: str, label: str, frame: pd.DataFrame, primary: dict[str, object], score: np.ndarray, paired_msi: dict[str, float] | None = None) -> list[dict[str, object]]:
    summary = pair_summary(score, primary)
    rows = []
    for celltype in ["BeamA", "JCT", "Schlemm_Endothelium"]:
        s = summary[celltype]
        rows.append({
            "signature_family": family,
            "size": label,
            "requested_n": int(label.split("=")[1]) if label.startswith("N=") else "RANKED_FULL",
            "realized_n": int(len(frame)),
            "layer": "GSE313736_donor_pseudobulk",
            "target_or_celltype": celltype,
            "n_donors": s["n"],
            "median_delta_score": s["median"],
            "positive_fraction": s["positive_fraction"],
            "median_delta_MSI": paired_msi.get(celltype, np.nan) if paired_msi else np.nan,
            "status": "OK" if s["n"] else "UNAVAILABLE",
            "notes": "same signed within-cell-state z-score procedure as Phase 1; donor is unit; causal DEX interpretation prohibited",
        })
    return rows


def build_size_sensitivity(frames: dict[str, pd.DataFrame], core_all: pd.DataFrame, primary: dict[str, object], rng: np.random.Generator) -> tuple[pd.DataFrame, dict[str, dict[str, pd.DataFrame]], dict[str, dict[str, np.ndarray]]]:
    family_frames = {
        "Adaptive": frames["Adaptive"],
        "Maladaptive-Fibrotic": frames["Fibrotic"],
        "Maladaptive-Steroid": frames["Steroid"],
        "Maladaptive-Disease": frames["Disease"],
        "Maladaptive-Core": core_all,
    }
    variants = {name: make_variants(frame) for name, frame in family_frames.items()}
    scores: dict[str, dict[str, np.ndarray]] = defaultdict(dict)
    rows: list[dict[str, object]] = []
    adaptive_scores = {}
    core_scores = {}
    target_map = {
        "Adaptive": ["Maladaptive-Fibrotic", "Maladaptive-Steroid", "Maladaptive-Disease"],
        "Maladaptive-Fibrotic": ["Maladaptive-Steroid", "Maladaptive-Disease"],
        "Maladaptive-Steroid": ["Maladaptive-Fibrotic", "Maladaptive-Disease"],
        "Maladaptive-Disease": ["Maladaptive-Fibrotic", "Maladaptive-Steroid"],
        "Maladaptive-Core": ["Adaptive"],
    }
    rng_external = np.random.default_rng(SEED + 11)
    for family, family_variants in variants.items():
        for label, variant in family_variants.items():
            score = score_frame(variant, primary)
            scores[family][label] = score
            if family == "Adaptive":
                adaptive_scores[label] = score
            if family == "Maladaptive-Core":
                core_scores[label] = score
            paired_msi = None
            if label in adaptive_scores and label in core_scores:
                a = pair_summary(adaptive_scores[label], primary)
                m = pair_summary(core_scores[label], primary)
                paired_msi = {ct: float(np.nanmedian(m[ct]["deltas"] - a[ct]["deltas"])) if a[ct]["n"] and m[ct]["n"] else np.nan for ct in primary["pair_indices"]}
            rows.extend(primary_rows_for_family(family, label, variant, primary, score, paired_msi))
            for target in target_map[family]:
                target_frame = family_frames[target]
                e = external_metrics(target_frame, set(variant["gene"]), target, rng_external, source_frame=variant)
                rows.append({
                    "signature_family": family,
                    "size": label,
                    "requested_n": int(label.split("=")[1]) if label.startswith("N=") else "RANKED_FULL",
                    "realized_n": int(len(variant)),
                    "layer": "external_ranked_validation",
                    "target_or_celltype": target,
                    "n_donors": "NA",
                    "median_delta_score": "NA",
                    "positive_fraction": "NA",
                    "median_delta_MSI": "NA",
                    "external_overlap_n": e["overlap_n"],
                    "external_target_significant_n": e["target_significant_n"],
                    "external_jaccard": e["jaccard"],
                    "external_mean_signed_statistic": e["mean_signed_statistic"],
                    "external_rank_z": e["rank_z"],
                    "external_empirical_p": e["empirical_p"],
                    "status": "OK",
                    "notes": "external projection; no GSE313736 tuning",
                })
            stress = set().union(*P1.GENERIC_STRESS.values())
            rows.append({
                "signature_family": family,
                "size": label,
                "requested_n": int(label.split("=")[1]) if label.startswith("N=") else "RANKED_FULL",
                "realized_n": int(len(variant)),
                "layer": "generic_stress_control",
                "target_or_celltype": "all_predefined_generic_stress_sets",
                "n_donors": "NA",
                "median_delta_score": "NA",
                "positive_fraction": "NA",
                "median_delta_MSI": "NA",
                "generic_stress_overlap_n": len(set(variant["gene"]) & stress),
                "generic_stress_overlap_fraction": len(set(variant["gene"]) & stress) / max(len(variant), 1),
                "status": "OK",
                "notes": "apoptosis/inflammation/ROS/hypoxia/UPR/TNF-NFkB/cell-cycle control union",
            })
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "SIGNATURE_SIZE_SENSITIVITY.tsv")
    return result, variants, scores


def null_strata(primary: dict[str, object], excluded: set[str]) -> tuple[dict[tuple[int, int, int], np.ndarray], np.ndarray, dict[str, tuple[int, int, int]]]:
    bg = pd.read_csv(OUT / "PRIMARY_BACKGROUND_SUMMARY.tsv", sep="\t")
    genes = bg["gene"].tolist()
    values = bg[["mean_logCPM", "detection_rate", "variance_logCPM"]].to_numpy(float)
    bins = []
    for col in range(values.shape[1]):
        v = values[:, col]
        finite = np.isfinite(v)
        qs = np.quantile(v[finite], np.linspace(0, 1, 9)) if finite.any() else np.arange(9)
        qs = np.unique(qs)
        bins.append(np.searchsorted(qs, v, side="right").clip(0, 8))
    strata: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    gene_to_stratum = {}
    for i, g in enumerate(genes):
        key = (int(bins[0][i]), int(bins[1][i]), int(bins[2][i]))
        gene_to_stratum[g] = key
        if g not in excluded and g in primary["gene_to_id"]:
            strata[key].append(primary["gene_to_id"][g])
    pools = {k: np.asarray(v, dtype=int) for k, v in strata.items() if v}
    candidates = np.asarray(sorted(set().union(*[set(x) for x in pools.values()])), dtype=int)
    return pools, candidates, gene_to_stratum


def nearest_pool(key: tuple[int, int, int], pools: dict[tuple[int, int, int], np.ndarray]) -> np.ndarray:
    if key in pools:
        return pools[key]
    choices = sorted(pools.items(), key=lambda item: sum(abs(a - b) for a, b in zip(key, item[0])))
    return choices[0][1]


def random_matched_ids(selected: pd.DataFrame, primary: dict[str, object], pools: dict[tuple[int, int, int], np.ndarray], gene_strata: dict[str, tuple[int, int, int]], candidates: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    genes = selected["gene"].tolist()
    signs = np.where(pd.to_numeric(selected["logFC"], errors="coerce").to_numpy(float) >= 0, 1.0, -1.0).astype(np.float32)
    out = []
    used = set()
    for g in genes:
        pool = nearest_pool(gene_strata.get(g, (4, 4, 4)), pools)
        pick = None
        for _ in range(12):
            candidate = int(rng.choice(pool))
            if candidate not in used:
                pick = candidate
                break
        if pick is None:
            remaining = np.asarray([x for x in candidates if int(x) not in used], dtype=int)
            pick = int(rng.choice(remaining if len(remaining) else candidates))
        out.append(pick)
        used.add(pick)
    return np.asarray(out, dtype=int), signs


def batch_null_deltas(ids_batch: np.ndarray, signs_batch: np.ndarray, primary: dict[str, object]) -> dict[str, np.ndarray]:
    result = {}
    for celltype, pairs in primary["pair_indices_local"].items():
        z = primary["z_by_celltype"][celltype]
        vals = z[ids_batch]
        score = np.nanmean(vals * signs_batch[:, :, None], axis=1)
        if not pairs:
            result[celltype] = np.full(len(ids_batch), np.nan)
            continue
        delta = np.asarray([score[:, dex] - score[:, veh] for veh, dex, _ in pairs]).T
        result[celltype] = np.nanmedian(delta, axis=1)
    return result


def matched_null_tests(frames: dict[str, pd.DataFrame], variants: dict[str, dict[str, pd.DataFrame]], primary: dict[str, object], rng: np.random.Generator) -> tuple[pd.DataFrame, dict[str, object]]:
    excluded = set().union(*[set(eligible(x)["gene"]) for x in frames.values()]) | set(P1.HOUSEKEEPING) | set().union(*P1.GENERIC_STRESS.values())
    pools, candidates, gene_strata = null_strata(primary, excluded)
    log(f"Matched-null background: candidates={len(candidates)} strata={len(pools)} sets_per_family_size={NULL_N}")
    tested = {"Adaptive": variants["Adaptive"], "Maladaptive-Core": variants["Maladaptive-Core"]}
    rows = []
    diagnostics = []
    for label in SIZE_LABELS:
        real_a = tested["Adaptive"][label]
        real_m = tested["Maladaptive-Core"][label]
        observed_scores = {"Adaptive": score_frame(real_a, primary), "Maladaptive-Core": score_frame(real_m, primary)}
        obs_delta = {fam: pair_summary(sc, primary) for fam, sc in observed_scores.items()}
        null_a = {ct: np.empty(NULL_N, dtype=float) for ct in primary["pair_indices"]}
        null_m = {ct: np.empty(NULL_N, dtype=float) for ct in primary["pair_indices"]}
        null_msi = {ct: np.empty(NULL_N, dtype=float) for ct in primary["pair_indices"]}
        for start in range(0, NULL_N, NULL_BATCH):
            stop = min(start + NULL_BATCH, NULL_N)
            a_ids, a_sign = [], []
            m_ids, m_sign = [], []
            for _ in range(start, stop):
                ai, ass = random_matched_ids(real_a, primary, pools, gene_strata, candidates, rng)
                mi, mss = random_matched_ids(real_m, primary, pools, gene_strata, candidates, rng)
                a_ids.append(ai)
                a_sign.append(ass)
                m_ids.append(mi)
                m_sign.append(mss)
            a_ids = np.vstack(a_ids)
            m_ids = np.vstack(m_ids)
            a_sign = np.vstack(a_sign)
            m_sign = np.vstack(m_sign)
            ad = batch_null_deltas(a_ids, a_sign, primary)
            md = batch_null_deltas(m_ids, m_sign, primary)
            for ct in primary["pair_indices"]:
                null_a[ct][start:stop] = ad[ct]
                null_m[ct][start:stop] = md[ct]
                null_msi[ct][start:stop] = md[ct] - ad[ct]
        for family, nulls, obs_family in [("Adaptive", null_a, obs_delta["Adaptive"]), ("Maladaptive-Core", null_m, obs_delta["Maladaptive-Core"])]:
            for celltype in ["BeamA", "JCT", "Schlemm_Endothelium"]:
                observed = fnum(obs_family[celltype]["median"])
                vals = nulls[celltype][np.isfinite(nulls[celltype])]
                p = float((1 + np.sum(vals >= observed)) / (len(vals) + 1)) if len(vals) else np.nan
                rows.append({
                    "signature_family": family,
                    "size": label,
                    "realized_n": len(tested[family][label]),
                    "metric": "median_delta_score",
                    "celltype": celltype,
                    "null_sets": NULL_N,
                    "observed_statistic": observed,
                    "null_median": float(np.median(vals)) if len(vals) else np.nan,
                    "null_ci95_low": float(np.quantile(vals, .025)) if len(vals) else np.nan,
                    "null_ci95_high": float(np.quantile(vals, .975)) if len(vals) else np.nan,
                    "empirical_p_direction_consistent": p,
                    "observed_percentile": float(np.mean(vals <= observed)) if len(vals) else np.nan,
                    "status": "PASS" if np.isfinite(p) and p < 0.05 else "FLAG",
                    "matching": "mean_logCPM+detection_rate+variance bins; gene length and GC unavailable",
                })
                if family == "Maladaptive-Core":
                    observed_msi = fnum(obs_delta["Maladaptive-Core"][celltype]["median"] - obs_delta["Adaptive"][celltype]["median"])
                    vals = null_msi[celltype][np.isfinite(null_msi[celltype])]
                    p = float((1 + np.sum(vals >= observed_msi)) / (len(vals) + 1)) if len(vals) else np.nan
                    rows.append({
                        "signature_family": "Adaptive_vs_Maladaptive",
                        "size": label,
                        "realized_n": f"{len(real_a)}/{len(real_m)}",
                        "metric": "median_delta_MSI",
                        "celltype": celltype,
                        "null_sets": NULL_N,
                        "observed_statistic": observed_msi,
                        "null_median": float(np.median(vals)) if len(vals) else np.nan,
                        "null_ci95_low": float(np.quantile(vals, .025)) if len(vals) else np.nan,
                        "null_ci95_high": float(np.quantile(vals, .975)) if len(vals) else np.nan,
                        "empirical_p_direction_consistent": p,
                        "observed_percentile": float(np.mean(vals <= observed_msi)) if len(vals) else np.nan,
                        "status": "PASS" if np.isfinite(p) and p < 0.05 else "FLAG",
                        "matching": "paired independent matched nulls; same signed pseudobulk scoring and MSI contrast",
                    })
        diagnostics.append({"size": label, "adaptive_mean_expr": float(np.mean([primary["mean_expr"][primary["gene_to_id"][g]] for g in real_a["gene"] if g in primary["gene_to_id"]])), "maladaptive_mean_expr": float(np.mean([primary["mean_expr"][primary["gene_to_id"][g]] for g in real_m["gene"] if g in primary["gene_to_id"]]))})
    result = pd.DataFrame(rows)
    result["FDR_across_phase1R_null_tests"] = bh(pd.to_numeric(result["empirical_p_direction_consistent"], errors="coerce").to_numpy())
    write_tsv(result, OUT / "MATCHED_NULL_EMPIRICAL_TESTS.tsv")
    return result, {"pools": len(pools), "candidates": len(candidates), "diagnostics": diagnostics}


def random_null_diagnosis(frames: dict[str, pd.DataFrame], core_all: pd.DataFrame, nulls: pd.DataFrame, size_result: pd.DataFrame) -> str:
    phase1 = pd.read_csv(ROOT / "results" / "phase1" / "FALSIFICATION_TESTS.tsv", sep="\t")
    original = phase1[phase1["test"].eq("random_matched_set_overlap")]
    original_p = fnum(original.iloc[0]["null_or_reference"]) if not original.empty else np.nan
    primary_rows = nulls[nulls["metric"].eq("median_delta_score") & nulls["size"].eq("N=100")]
    pass_fraction = float(np.mean(primary_rows["status"].eq("PASS"))) if not primary_rows.empty else 0.0
    msis = nulls[nulls["metric"].eq("median_delta_MSI") & nulls["size"].eq("N=100")]
    msi_pass = float(np.mean(msis["status"].eq("PASS"))) if not msis.empty else 0.0
    stress = set().union(*P1.GENERIC_STRESS.values())
    core_genes = set(core_all["gene"])
    modules = {
        "highly_expressed_ECM": {"COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "DCN", "LUM", "SPARC", "VCAN", "CCN2"},
        "generic_DEx_response": stress,
    }
    module_lines = "\n".join(f"- {name}: overlap_with_maladaptive_core={len(core_genes & genes)}" for name, genes in modules.items())
    if original_p < 0.05 and pass_fraction >= 0.5 and msi_pass >= 0.5:
        conclusion = "RESOLVED_NULL_MISSPECIFICATION"
    elif original_p < 0.05 and pass_fraction >= 0.5:
        conclusion = "RESOLVED_BUT_EFFECT_ATTENUATED"
    elif pass_fraction == 0 and msi_pass == 0:
        conclusion = "FAIL_SIGNATURE_SPECIFICITY"
    else:
        conclusion = "UNRESOLVED_NULL_OVERLAP"
    text = f"""# Phase 1R random-null flag diagnosis

## Locked question

The original Phase 1 flag was not silently accepted or discarded. It was
reconstructed and compared with a new null that matches the primary
GSE313736 baseline as closely as feasible.

## What the original test actually measured

Phase 1 used a small uniform external-gene overlap null for the observed
Adaptive-Core versus Maladaptive-Core intersection. It did not match primary
TM expression, detection, variance, gene length, or GC content, and it did
not apply the full donor-pseudobulk signed-score procedure to each null set.
The original empirical reference was {original_p:.6g}; the observed overlap
was therefore a sparse-set/background calibration result, not a specificity
test for the GSE313736 score.

## Hypothesis audit

- Expression mismatch: Phase 1 did not match mean logCPM; Phase 1R matched mean logCPM bins.
- Detection mismatch: Phase 1 did not match TM detection; Phase 1R matched detection-rate bins from the exact frame.
- Signature-size artifact: Phase 1R tested N=25, 50, 75, 100, 150, 200 and ranked-full versions.
- Highly expressed ECM genes: tested as a pre-specified overlap control, not as a post-hoc pathway.
- Generic DEX-responsive genes: compared with the frozen apoptosis, inflammation, ROS, hypoxia, UPR, TNF/NF-kB and cell-cycle controls.
- Correlated gene modules: evaluated with fixed pathway-module removal and cross-context rank concordance.
- Donor-level global transcriptional shift: donor-paired null Δscore and ΔMSI were computed with the same within-cell-state z-scoring procedure.
- Normalization/scoring bias: real and null sets used identical pseudobulk totals, within-cell-state standardization, signed averaging, and donor pairing.

## Observed diagnostics

- Original Phase 1 overlap reference: {original_p:.6g}.
- At N=100, matched-null primary score rows passing the direction-consistent empirical test: {pass_fraction:.3f}.
- At N=100, matched-null ΔMSI rows passing: {msi_pass:.3f}.
{module_lines}

Gene length and GC content were unavailable in the supplied project inputs and
are explicitly not claimed as matched. The null therefore controls expression,
detection, variance, and gene-set size, with baseline TM detectability captured
by the exact GSE313736 frame.

## Required conclusion

**{conclusion}**

This label is a scope decision. It does not convert the observational
GSE313736 treatment contrast into a causal DEX estimate.
"""
    (REPORT / "RANDOM_NULL_FLAG_DIAGNOSIS.md").write_text(text, encoding="utf-8")
    return conclusion


def context_support(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    all_genes = sorted(set().union(*[set(eligible(x)["gene"]) for x in frames.values()]))
    lookups = {name: frame.set_index("gene") for name, frame in frames.items()}
    rows = []
    for g in all_genes:
        vals = {}
        for name in ["Fibrotic", "Steroid", "Disease"]:
            row = lookups[name].loc[g] if g in lookups[name].index else None
            vals[name] = {"logFC": fnum(row["logFC"]) if row is not None else np.nan, "FDR": fnum(row["FDR"], 1.0) if row is not None else 1.0}
        tf = vals["Fibrotic"]["FDR"] < FDR_CUTOFF
        dx = vals["Steroid"]["FDR"] < FDR_CUTOFF
        ds = vals["Disease"]["FDR"] < FDR_CUTOFF
        td_conc = tf and dx and np.sign(vals["Fibrotic"]["logFC"]) == np.sign(vals["Steroid"]["logFC"]) and vals["Fibrotic"]["logFC"] != 0 and vals["Steroid"]["logFC"] != 0
        disease_conc = ds and ((tf and np.sign(vals["Disease"]["logFC"]) == np.sign(vals["Fibrotic"]["logFC"])) or (dx and np.sign(vals["Disease"]["logFC"]) == np.sign(vals["Steroid"]["logFC"])))
        if td_conc and disease_conc:
            label = "BOTH_PLUS_DISEASE_CONCORDANT"
        elif td_conc:
            label = "BOTH_TGFb2_DEX"
        elif tf:
            label = "TGFb2_ONLY"
        elif dx:
            label = "DEX_ONLY"
        elif ds:
            label = "DISEASE_CONTEXT_SUPPORT_ONLY"
        else:
            label = "NONE"
        rows.append({
            "gene": g,
            "TGFb2_status": "SIGNIFICANT" if tf else "NOT_SIGNIFICANT",
            "DEX_status": "SIGNIFICANT" if dx else "NOT_SIGNIFICANT",
            "disease_status": "SIGNIFICANT" if ds else "NOT_SIGNIFICANT",
            "context_support_label": label,
            "TGFb2_logFC": vals["Fibrotic"]["logFC"],
            "DEX_logFC": vals["Steroid"]["logFC"],
            "disease_logFC": vals["Disease"]["logFC"],
            "TGFb2_FDR": vals["Fibrotic"]["FDR"],
            "DEX_FDR": vals["Steroid"]["FDR"],
            "disease_FDR": vals["Disease"]["FDR"],
            "direction_concordance_TGFb2_DEX": "YES" if td_conc else "NO",
            "eligible_for_maladaptive_core": "YES" if td_conc else "NO",
        })
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "MALADAPTIVE_CONTEXT_SUPPORT.tsv")
    return result


def primary_signature_summary(frame: pd.DataFrame, primary: dict[str, object]) -> dict[str, dict[str, object]]:
    return pair_summary(score_frame(frame, primary), primary)


def loco_validation(frames: dict[str, pd.DataFrame], core_all: pd.DataFrame, primary: dict[str, object]) -> tuple[pd.DataFrame, dict[str, bool]]:
    tgf = eligible(frames["Fibrotic"]).head(100)
    dex = eligible(frames["Steroid"]).head(100)
    models = [
        ("M1_TGFb2_derived_tested_in_DEX", tgf, "Steroid", "GSE301517", "TGFb2-only signature in DEX"),
        ("M2_DEX_derived_tested_in_TGFb2", dex, "Fibrotic", "GSE202372", "DEX-only signature in TGFb2"),
        ("M3_consensus_excluding_disease_tested_in_GSE263230", core_all.head(100), "Disease", "GSE263230", "consensus excluding disease; descriptive disease context"),
        ("M4_TGFb2_only_applied_to_GSE313736", tgf, None, "GSE313736", "TGFb2-only primary application"),
        ("M5_DEX_only_applied_to_GSE313736", dex, None, "GSE313736", "DEX-only primary application"),
    ]
    rows = []
    gates = {}
    for name, frame, target_family, target_name, note in models:
        e = external_metrics(frames[target_family], set(frame["gene"]), target_name, np.random.default_rng(SEED + len(rows) + 80), source_frame=frame) if target_family else {}
        ps = primary_signature_summary(frame, primary) if target_family is None else {}
        if ps:
            ok = all(ps[ct]["n"] >= 50 and ps[ct]["median"] > 0 and ps[ct]["positive_fraction"] >= 0.5 for ct in ["BeamA", "JCT", "Schlemm_Endothelium"])
            for ct in ["BeamA", "JCT", "Schlemm_Endothelium"]:
                rows.append({"model": name, "signature_n": len(frame), "target": target_name, "celltype": ct, "median_delta_score": ps[ct]["median"], "positive_fraction": ps[ct]["positive_fraction"], "n_donors": ps[ct]["n"], "external_mean_signed_statistic": np.nan, "external_rank_z": np.nan, "external_empirical_p": np.nan, "status": "PASS" if ok else "FLAG", "notes": note})
        else:
            ok = np.isfinite(fnum(e.get("mean_signed_statistic"))) and fnum(e.get("rank_z"), -np.inf) > 0
            rows.append({"model": name, "signature_n": len(frame), "target": target_name, "celltype": "EXTERNAL", "median_delta_score": np.nan, "positive_fraction": np.nan, "n_donors": "NA", "external_mean_signed_statistic": e.get("mean_signed_statistic", np.nan), "external_rank_z": e.get("rank_z", np.nan), "external_empirical_p": e.get("empirical_p", np.nan), "status": "PASS" if ok else "FLAG", "notes": note})
        gates[name] = bool(ok)
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "LEAVE_ONE_CONTEXT_OUT_VALIDATION.tsv")
    return result, gates


def rank_concordance(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    names = ["Adaptive", "Fibrotic", "Steroid", "Disease"]
    rows = []
    for i, a_name in enumerate(names):
        a = frames[a_name].set_index("gene")[["logFC", "statistic", "rank"]].rename(columns=lambda x: f"a_{x}")
        for b_name in names[i + 1:]:
            b = frames[b_name].set_index("gene")[["logFC", "statistic", "rank"]].rename(columns=lambda x: f"b_{x}")
            m = a.join(b, how="inner").replace([np.inf, -np.inf], np.nan).dropna()
            if m.empty:
                continue
            rho_stat = spearmanr(m["a_statistic"], m["b_statistic"]).statistic
            rho_fc = spearmanr(m["a_logFC"], m["b_logFC"]).statistic
            signs = np.sign(m["a_logFC"]) * np.sign(m["b_logFC"])
            nonzero = signs.ne(0)
            directional = float(np.mean(signs[nonzero] > 0)) if nonzero.any() else np.nan
            overlap_rows = []
            M = len(m)
            best_p = 1.0
            best_k = np.nan
            best_overlap = 0
            for k in [25, 50, 100, 200, 500]:
                ka = set(m.sort_values("a_rank").head(min(k, M)).index)
                kb = set(m.sort_values("b_rank").head(min(k, M)).index)
                ov = len(ka & kb)
                p = float(hypergeom.sf(ov - 1, M, len(ka), len(kb)))
                if p < best_p:
                    best_p, best_k, best_overlap = p, k, ov
            rows.append({"dataset_a": a_name, "dataset_b": b_name, "n_shared_genes": len(m), "spearman_statistic_rho": rho_stat, "spearman_logFC_rho": rho_fc, "directional_concordance": directional, "best_rank_overlap_k": best_k, "best_rank_overlap_n": best_overlap, "best_rank_overlap_hypergeom_p": best_p, "interpretation": "rank/effect concordance; not top-gene overlap alone"})
    result = pd.DataFrame(rows)
    result["best_rank_overlap_FDR"] = bh(result["best_rank_overlap_hypergeom_p"].to_numpy())
    write_tsv(result, OUT / "CROSS_DATASET_RANK_CONCORDANCE.tsv")
    return result


def overlap_table(frames: dict[str, pd.DataFrame], core_all: pd.DataFrame) -> pd.DataFrame:
    av = make_variants(frames["Adaptive"])
    mv = make_variants(core_all)
    fib = frames["Fibrotic"].set_index("gene")
    dex = frames["Steroid"].set_index("gene")
    rows = []
    for label in SIZE_LABELS:
        ag = set(av[label]["gene"])
        mg = set(mv[label]["gene"])
        inter = sorted(ag & mg)
        cats = []
        for g in inter:
            asgn = np.sign(fnum(av[label].set_index("gene").loc[g]["logFC"]))
            ctx = [np.sign(fnum(fib.loc[g]["logFC"])) if g in fib.index else np.nan, np.sign(fnum(dex.loc[g]["logFC"])) if g in dex.index else np.nan]
            ctx = [x for x in ctx if np.isfinite(x) and x != 0]
            if len(set(ctx)) > 1:
                cat = "CONTEXT_DEPENDENT"
            elif ctx and asgn == ctx[0]:
                cat = "SAME_DIRECTION"
            elif ctx:
                cat = "OPPOSITE_DIRECTION"
            else:
                cat = "CONTEXT_DEPENDENT"
            cats.append(cat)
        rows.append({"size": label, "adaptive_n": len(ag), "maladaptive_n": len(mg), "intersection_n": len(inter), "jaccard": len(inter) / max(len(ag | mg), 1), "overlap_coefficient": len(inter) / max(min(len(ag), len(mg)), 1), "same_direction_n": cats.count("SAME_DIRECTION"), "opposite_direction_n": cats.count("OPPOSITE_DIRECTION"), "context_dependent_n": cats.count("CONTEXT_DEPENDENT"), "overlap_genes": ";".join(inter), "direction_classification": ";".join(f"{g}:{c}" for g, c in zip(inter, cats))})
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "ADAPTIVE_MALADAPTIVE_OVERLAP.tsv")
    return result


PATHWAY_MODULES = {
    "focal_adhesion": {"PTK2", "VCL", "TLN1", "PXN", "FERMT2", "ITGA9", "CAV1", "CAV2", "ACTN1", "ACTB"},
    "integrin_signaling": {"ITGA9", "ITGB1", "ITGB3", "ITGAV", "ITGBL1", "PTK2", "TLN1", "FERMT2", "PXN", "VCL"},
    "actomyosin_contractility": {"ACTA2", "TAGLN", "MYL9", "RHOA", "ROCK1", "ROCK2", "MYLK", "MYL12A", "MYL12B", "TPM1", "TPM2", "ACTB"},
    "ECM_organization": {"FN1", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "DCN", "LUM", "SPARC", "VCAN", "CCN2", "TGFB2", "THBS1"},
    "TGF_beta_signaling": {"TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7", "SERPINE1", "CCN2", "FN1"},
    "Hippo_YAP": {"YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "LATS1", "LATS2", "MOB1A", "MOB1B", "CCN2", "CTGF", "AMOT"},
    "calcium_signaling": {"PIEZO1", "TRPV4", "CALM1", "CALM2", "CALM3", "CAMK2A", "CAMK2B", "ITPR1", "ITPR2", "ITPR3", "ATP2B1", "ATP2B2", "ATP2B4", "ORAI1", "STIM1"},
    "oxidative_stress": {"NOX4", "CYBB", "SOD1", "SOD2", "CAT", "GPX1", "GPX4", "HMOX1", "NQO1", "GCLC", "GCLM", "SLC7A11", "NFE2L2"},
    "ferroptosis": {"GPX4", "SLC7A11", "ACSL4", "TFRC", "HMOX1", "NFE2L2", "FTH1", "FTL", "ALOX5", "LPCAT3"},
    "autophagy": {"BECN1", "ATG3", "ATG5", "ATG7", "ATG12", "ATG16L1", "MAP1LC3B", "SQSTM1", "ULK1", "LAMP1", "LAMP2"},
    "proteostasis": {"HSPA5", "XBP1", "ATF4", "DDIT3", "EIF2AK3", "HSP90AA1", "HSPA1A", "DNAJB9", "MANF", "PSMC1", "UBB"},
    "mechanosensitive_channels": {"PIEZO1", "TRPV4", "PKD1", "PKD2", "TMEM63A", "TMEM63B", "TMEM63C", "ASIC1", "ASIC2"},
}


def pathway_robustness(frames: dict[str, pd.DataFrame], core_all: pd.DataFrame, primary: dict[str, object], variants: dict[str, dict[str, pd.DataFrame]]) -> pd.DataFrame:
    families = {"Adaptive": variants["Adaptive"], "Maladaptive-Core": variants["Maladaptive-Core"]}
    rows = []
    primary_bg = set(primary["gene_list"])
    for module, module_genes in PATHWAY_MODULES.items():
        for family, family_variants in families.items():
            for label, frame in family_variants.items():
                sig_genes = set(frame["gene"]) & primary_bg
                term = module_genes & primary_bg
                overlap = len(sig_genes & term)
                p = float(hypergeom.sf(overlap - 1, len(primary_bg), len(term), len(sig_genes))) if sig_genes and term else np.nan
                score = score_gene_set(term, primary)
                ps = pair_summary(score, primary)
                for ct in ["BeamA", "JCT", "Schlemm_Endothelium"]:
                    rows.append({"module": module, "signature_family": family, "size": label, "source_library": "pre_specified_GO_Reactome_Hallmark_aligned_module; no GSE313736 tuning", "module_n_in_primary": len(term), "signature_n_in_primary": len(sig_genes), "overlap_n": overlap, "odds_ratio": (overlap * max(len(primary_bg) - len(term) - len(sig_genes) + overlap, 1)) / max((len(term) - overlap) * (len(sig_genes) - overlap), 1), "PValue": p, "primary_celltype": ct, "primary_median_delta": ps[ct]["median"], "primary_positive_fraction": ps[ct]["positive_fraction"], "status": "OK"})
    result = pd.DataFrame(rows)
    result["FDR"] = bh(result["PValue"].to_numpy())
    write_tsv(result, OUT / "PATHWAY_ROBUSTNESS.tsv")
    return result


def anchor_removal(frames: dict[str, pd.DataFrame], core_all: pd.DataFrame, primary: dict[str, object], scores: dict[str, dict[str, np.ndarray]]) -> tuple[pd.DataFrame, bool]:
    adaptive100 = make_variants(frames["Adaptive"])["N=100"]
    core100 = make_variants(core_all)["N=100"]
    tests = [("remove_PIEZO1_TRPV4", {"PIEZO1", "TRPV4"}), ("remove_anchor_panel", set(P1.ANCHOR_PANEL))]
    rows = []
    all_pass = True
    for removal, remove in tests:
        for family, frame, original in [("Adaptive", adaptive100, scores["Adaptive"]["N=100"]), ("Maladaptive-Core", core100, scores["Maladaptive-Core"]["N=100"])]:
            new = score_frame(frame, primary, remove)
            old_s = pair_summary(original, primary)
            new_s = pair_summary(new, primary)
            for ct in ["BeamA", "JCT", "Schlemm_Endothelium"]:
                oi, ni = old_s[ct], new_s[ct]
                keep = bool(ni["median"] > 0 and ni["positive_fraction"] >= 0.5)
                all_pass = all_pass and keep
                idx = primary["group_indices"][ct]
                a = original[idx]
                b = new[idx]
                corr = float(np.corrcoef(a, b)[0, 1]) if np.nanstd(a) > 0 and np.nanstd(b) > 0 else np.nan
                rows.append({"removal": removal, "removed_n_in_signature": len(set(frame["gene"]) & remove), "signature_family": family, "celltype": ct, "original_median_delta": oi["median"], "removed_median_delta": ni["median"], "original_positive_fraction": oi["positive_fraction"], "removed_positive_fraction": ni["positive_fraction"], "group_score_correlation": corr, "direction_preserved": "YES" if keep else "NO", "status": "PASS" if keep else "FLAG"})
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "ANCHOR_REMOVAL_SENSITIVITY.tsv")
    return result, all_pass


def single_pathway_dominance(frames: dict[str, pd.DataFrame], core_all: pd.DataFrame, primary: dict[str, object], scores: dict[str, dict[str, np.ndarray]]) -> tuple[pd.DataFrame, str]:
    core100 = make_variants(core_all)["N=100"]
    base = primary_signature_summary(core100, primary)
    removal_modules = {
        "ECM_collagen": PATHWAY_MODULES["ECM_organization"],
        "TGF_beta": PATHWAY_MODULES["TGF_beta_signaling"],
        "actomyosin": PATHWAY_MODULES["actomyosin_contractility"],
        "oxidative_stress": PATHWAY_MODULES["oxidative_stress"],
        "Hippo_YAP": PATHWAY_MODULES["Hippo_YAP"],
    }
    rows = []
    collapse = []
    module_overlap_total = 0
    for name, module in removal_modules.items():
        remaining = core100[~core100["gene"].isin(module)].copy()
        module_overlap = len(set(core100["gene"]) & module)
        module_overlap_total += module_overlap
        sc = score_frame(remaining, primary)
        ps = pair_summary(sc, primary)
        ext_tgf = external_metrics(frames["Fibrotic"], set(remaining["gene"]), "GSE202372", np.random.default_rng(SEED + len(rows) + 120), source_frame=remaining)
        ext_dex = external_metrics(frames["Steroid"], set(remaining["gene"]), "GSE301517", np.random.default_rng(SEED + len(rows) + 121), source_frame=remaining)
        preserved = all(ps[ct]["median"] > 0 and ps[ct]["positive_fraction"] >= 0.5 for ct in ["BeamA", "JCT", "Schlemm_Endothelium"])
        external_preserved = fnum(ext_tgf["mean_signed_statistic"], -np.inf) > 0 and fnum(ext_dex["mean_signed_statistic"], -np.inf) > 0
        collapsed = not preserved and not external_preserved
        collapse.append(collapsed)
        for ct in ["BeamA", "JCT", "Schlemm_Endothelium"]:
            rows.append({"removed_module": name, "removed_n": module_overlap, "remaining_n": len(remaining), "celltype": ct, "median_delta_score": ps[ct]["median"], "positive_fraction": ps[ct]["positive_fraction"], "external_TGFb2_mean_signed_statistic": ext_tgf["mean_signed_statistic"], "external_DEX_mean_signed_statistic": ext_dex["mean_signed_statistic"], "direction_preserved": "YES" if preserved else "NO", "status": "FLAG_SINGLE_MODULE_COLLAPSE" if collapsed else "PASS_MULTI_MODULE", "notes": "zero overlap with this pre-specified module; leave-one-module-out does not establish module independence" if module_overlap == 0 else "pre-specified leave-one-module-out"})
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "SINGLE_PATHWAY_DOMINANCE.tsv")
    verdict = "SINGLE_PATHWAY_DOMINATED" if any(collapse) else "MULTI_MODULE_ROBUST"
    log(f"Single-pathway test: verdict={verdict}; core-module overlap across tested modules={module_overlap_total}")
    return result, verdict


def donor_bootstrap(scores: dict[str, dict[str, np.ndarray]], primary: dict[str, object]) -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    actual_a = scores["Adaptive"]["N=100"]
    actual_m = scores["Maladaptive-Core"]["N=100"]
    metric_scores = {"Adaptive": actual_a, "Maladaptive-Core": actual_m, "MSI": actual_m - actual_a}
    rng = np.random.default_rng(SEED + 2000)
    rows = []
    summary = {}
    for ct, pairs in primary["pair_indices"].items():
        n = len(pairs)
        summary[ct] = {}
        for metric, score in metric_scores.items():
            delta = np.asarray([score[dex] - score[veh] for veh, dex, _ in pairs], dtype=float)
            if n == 0:
                continue
            boot = np.asarray([np.median(delta[rng.integers(0, n, size=n)]) for _ in range(BOOTSTRAPS)])
            record = {"celltype": ct, "metric": metric, "n_donors": n, "B": BOOTSTRAPS, "observed_median_delta": float(np.median(delta)), "bootstrap_median": float(np.median(boot)), "ci95_low": float(np.quantile(boot, .025)), "ci95_high": float(np.quantile(boot, .975)), "directional_concordance": float(np.mean(delta > 0)), "interpretation": "donor stability of descriptive concordance; not causal DEX inference"}
            rows.append(record)
            summary[ct][metric] = record
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "DONOR_BOOTSTRAP_STABILITY.tsv")
    return result, summary


def celltype_heterogeneity(scores: dict[str, dict[str, np.ndarray]], primary: dict[str, object]) -> pd.DataFrame:
    a = pair_summary(scores["Adaptive"]["N=100"], primary)
    m = pair_summary(scores["Maladaptive-Core"]["N=100"], primary)
    rows = []
    for ct in ["BeamA", "JCT", "Schlemm_Endothelium"]:
        shared = all(x[ct]["median"] > 0 and x[ct]["positive_fraction"] >= 0.5 for x in [a, m])
        if shared:
            label = "SHARED_ACROSS_THREE"
        else:
            vals = {"BeamA": a[ct]["median"] + m[ct]["median"], "JCT": a[ct]["median"] + m[ct]["median"], "Schlemm_Endothelium": a[ct]["median"] + m[ct]["median"]}
            label = "TM_SUBTYPE_SPECIFIC_OR_MIXED"
        rows.append({"celltype": ct, "role": "PRIMARY" if ct in {"BeamA", "JCT"} else "SECONDARY", "adaptive_median_delta": a[ct]["median"], "adaptive_positive_fraction": a[ct]["positive_fraction"], "maladaptive_median_delta": m[ct]["median"], "maladaptive_positive_fraction": m[ct]["positive_fraction"], "MSI_median_delta": m[ct]["median"] - a[ct]["median"], "classification": label, "interpretation": "cell states remain separate; no pooled significance claim"})
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "CELLTYPE_PROGRAM_HETEROGENEITY.tsv")
    return result


def final_decision(size_result: pd.DataFrame, nulls: pd.DataFrame, loco: pd.DataFrame, loco_gates: dict[str, bool], anchor_pass: bool, dominance: str, rank: pd.DataFrame) -> dict[str, object]:
    required = size_result[(size_result["layer"].eq("GSE313736_donor_pseudobulk")) & (size_result["signature_family"].isin(["Adaptive", "Maladaptive-Core"])) & (size_result["size"].isin(["N=50", "N=75", "N=100", "N=150"]))]
    size_ok = bool(not required.empty and np.all(pd.to_numeric(required["median_delta_score"], errors="coerce") > 0) and np.all(pd.to_numeric(required["positive_fraction"], errors="coerce") >= 0.5))
    null_primary = nulls[(nulls["metric"].eq("median_delta_score")) & (nulls["size"].isin(["N=50", "N=75", "N=100", "N=150"]))]
    null_fraction = float(np.mean(null_primary["status"].eq("PASS"))) if not null_primary.empty else 0.0
    null_gate = bool(null_fraction >= 0.80 and np.all(nulls[(nulls["metric"].eq("median_delta_score")) & (nulls["size"].eq("N=100")) & (nulls["celltype"].isin(["BeamA", "JCT"]))]["status"].eq("PASS")))
    loco_tgf = bool(loco_gates.get("M4_TGFb2_only_applied_to_GSE313736", False))
    loco_dex = bool(loco_gates.get("M5_DEX_only_applied_to_GSE313736", False))
    loco_pass = bool(all(loco_gates.get(x, False) for x in ["M1_TGFb2_derived_tested_in_DEX", "M2_DEX_derived_tested_in_TGFb2", "M4_TGFb2_only_applied_to_GSE313736", "M5_DEX_only_applied_to_GSE313736"]))
    rank_rows = rank[(rank["dataset_a"].isin(["Adaptive"]) | rank["dataset_b"].isin(["Adaptive"]))]
    rank_sep = bool(not rank_rows.empty and np.nanmedian(pd.to_numeric(rank_rows["directional_concordance"], errors="coerce")) < 0.75)
    if size_ok and null_gate and loco_pass and anchor_pass and dominance != "SINGLE_PATHWAY_DOMINATED":
        representation = "CORE_GENE_SET"
        verdict = "GO_PHASE2_REGULATORY"
        terminology = "adaptive–maladaptive state axis"
    elif size_ok and dominance != "SINGLE_PATHWAY_DOMINATED" and (null_fraction >= 0.50 or loco_tgf or loco_dex):
        representation = "RANKED_PROGRAM"
        verdict = "GO_PHASE2_CONDITIONAL"
        terminology = "mechanoadaptive-associated program / mechanopathology-associated program"
    elif dominance != "SINGLE_PATHWAY_DOMINATED":
        representation = "PATHWAY_MODULE"
        verdict = "GO_PHASE2_CONDITIONAL" if (loco_tgf or loco_dex) else "REFRAME_CONTEXT_DEPENDENT"
        terminology = "context-dependent mechanotransduction pathway modules"
    else:
        representation = "PATHWAY_MODULE"
        verdict = "REFRAME_CONTEXT_DEPENDENT"
        terminology = "context-dependent mechanotransduction pathway modules"
    scope = "ORTHOGONAL_MULTIOMIC_REGULATORY_ANCHORING (GSE301516/301519/301522 anchored to GSE301517; not independent biological replication)" if verdict != "REFRAME_CONTEXT_DEPENDENT" else "NO_REGULATORY_PHASE2_UNTIL_REFRAMED"
    return {"verdict": verdict, "top_n_dependence": "YES", "size_robust": "PASS" if size_ok else "FLAG", "matched_null_specificity": "PASS" if null_gate else ("FLAG" if null_fraction >= 0.50 else "FAIL"), "matched_null_pass_fraction": null_fraction, "loco": "PASS" if loco_pass else "FLAG", "loco_tgf_primary": "PASS" if loco_tgf else "FAIL", "loco_dex_primary": "PASS" if loco_dex else "FAIL", "rank_separation": "PASS" if rank_sep else "FLAG", "anchor_removal": "PASS" if anchor_pass else "FAIL", "single_pathway_dominated": "YES" if dominance == "SINGLE_PATHWAY_DOMINATED" else "NO", "final_representation": representation, "allowed_terminology": terminology, "phase2_allowed": "YES" if verdict != "REFRAME_CONTEXT_DEPENDENT" else "NO", "phase2_scope": scope, "biggest_limitation": "GSE313736 treatment remains perfectly aliased with pool/library; PRS provenance and 74-versus-59 donor discrepancy remain unresolved; only TGFβ2 and DEX are discovery pathological contexts", "size_required_rows": len(required)}


def final_signature_model(decision: dict[str, object], frames: dict[str, pd.DataFrame], core_all: pd.DataFrame) -> pd.DataFrame:
    """Write the single representation selected before any Phase 2 work."""
    rows = [
        {"component": "adaptive", "representation": decision["final_representation"], "source": "GSE195756", "locked_gene_set": "results/phase1/ADAPTIVE_SIGNATURE.tsv", "gene_n_reference": 100, "status": "PATHWAY_OR_PROGRAM_SCOPE_LOCKED"},
        {"component": "maladaptive", "representation": decision["final_representation"], "source": "GSE202372+GSE301517; GSE263230 contextual only", "locked_gene_set": "results/phase1/MALADAPTIVE_CORE.tsv", "gene_n_reference": 100, "status": "PATHWAY_OR_PROGRAM_SCOPE_LOCKED"},
        {"component": "adaptive_module_panel", "representation": decision["final_representation"], "source": "pre-specified GO/Reactome/Hallmark-aligned modules", "locked_gene_set": "results/phase1R/PATHWAY_ROBUSTNESS.tsv", "gene_n_reference": "NA", "status": "DESCRIPTIVE_REGULATORY_ANCHOR_ONLY"},
        {"component": "maladaptive_module_panel", "representation": decision["final_representation"], "source": "pre-specified GO/Reactome/Hallmark-aligned modules", "locked_gene_set": "results/phase1R/PATHWAY_ROBUSTNESS.tsv", "gene_n_reference": "NA", "status": "DESCRIPTIVE_REGULATORY_ANCHOR_ONLY"},
        {"component": "terminology", "representation": decision["allowed_terminology"], "source": "Phase 1R scope decision", "locked_gene_set": "reports/phase1R/PHASE1R_SIGNATURE_ROBUSTNESS_REPORT.md", "gene_n_reference": "NA", "status": "LOCKED"},
    ]
    result = pd.DataFrame(rows)
    write_tsv(result, OUT / "FINAL_SIGNATURE_MODEL.tsv")
    return result


def phase2_preregistration(decision: dict[str, object]) -> None:
    text = f"""# Phase 2 preregistration proposal — regulatory anchoring only

Phase 1R verdict: **{decision['verdict']}**  
Final representation: **{decision['final_representation']}**  
Phase 2 allowed in this proposal: **{decision['phase2_allowed']}**

This file is a preregistration boundary, not an execution record. Phase 2
regulatory multiomics was not run in Phase 1R.

## Exact scope

{decision['phase2_scope']}

The proposed regulatory resources are GSE301516 ATAC-seq, GSE301519 CUT&Tag,
GSE301522 promoter-capture/Hi-C, and GSE301517 RNA-seq from the same GSE301525
experimental framework. These assays may provide orthogonal regulatory
anchoring of the externally locked adaptive and mechanopathology-associated
programs. They must not be described as independent biological replication.

## Locked questions

1. Are DEX-associated accessible chromatin regions linked to the frozen
   maladaptive program and distinguishable from adaptive stretch-associated
   loci?
2. Do CUT&Tag marks and promoter contacts support candidate enhancers or
   regulatory contacts near PIEZO1/TRPV4-associated modules and the broader
   pathway representation?
3. Do RNA, accessibility, histone-mark and contact evidence converge on the
   pre-specified modules without redefining the signatures from GSE313736?
4. Can any glaucoma GWAS mapping be reported only as regulatory anchoring,
   with genome build, variant source, allele harmonization and locus rules
   fixed before inspection?

## Prohibited at Phase 2 entry

- no causal DEX claim from GSE313736;
- no ordinary covariate adjustment for the perfectly aliased treatment/pool;
- no PRS or PRS×DEX analysis until provenance, build, weights and donor mapping
  are closed;
- no same-super-series assay described as independent replication;
- no signature redefinition from human-TM primary expression;
- no clinical, therapeutic or mechanistically proven wording.

## Required analysis-unit rule

The donor remains the biological unit where donor-level information exists.
Regulatory peaks, marks and contacts are evidence objects anchored to the
pre-specified signature/module, not independent donor observations. All
multiomic maps require assay-specific QC, genome-build checks, duplicate and
blacklist handling, and a predeclared multiple-testing plan.
"""
    (REPORT / "PHASE2_PREREGISTRATION.md").write_text(text, encoding="utf-8")


def write_report(decision: dict[str, object], audit: pd.DataFrame, context: pd.DataFrame, overlap: pd.DataFrame, pathway: pd.DataFrame, bootstrap: pd.DataFrame, heterogeneity: pd.DataFrame, null_diag: str, primary: dict[str, object], rank: pd.DataFrame) -> None:
    cap_rows = audit[audit["signature"].isin(["Adaptive-Core", "Maladaptive-Fibrotic", "Maladaptive-Steroid", "Maladaptive-Core"])]
    core_context = context[context["eligible_for_maladaptive_core"].eq("YES")]
    null100 = pd.read_csv(OUT / "MATCHED_NULL_EMPIRICAL_TESTS.tsv", sep="\t")
    null100 = null100[null100["size"].eq("N=100")]
    lines = []
    for ct in ["BeamA", "JCT", "Schlemm_Endothelium"]:
        r = bootstrap[(bootstrap["celltype"].eq(ct)) & (bootstrap["metric"].eq("MSI"))]
        if not r.empty:
            x = r.iloc[0]
            lines.append(f"- {ct}: observed median ΔMSI={x['observed_median_delta']:.4f}; bootstrap median={x['bootstrap_median']:.4f}; 95% CI [{x['ci95_low']:.4f}, {x['ci95_high']:.4f}]; directional concordance={x['directional_concordance']:.3f}.")
    report = f"""# Phase 1R signature robustness, null calibration and scope lock

Date: 2026-09-18  
Project: PIEZO1/TRPV4 adaptive versus maladaptive mechanotransduction in glaucoma  
Verdict: **{decision['verdict']}**

## Executive gate

Phase 1R used the immutable Phase 1 external tables and preserved the Phase 0R
constraints. GSE313736 was reopened only for a single exact-frame aggregation
supporting predeclared size sensitivity, matched nulls, anchor removal and
donor bootstrap. It remains ROLE C. The perfectly aliased treatment/pool
design still prohibits a causal DEX estimate and PRS provenance remains
unresolved.

| Gate | Result |
|---|---|
| 100-gene cap dependence | **{decision['top_n_dependence']}** |
| Fixed-size robustness | **{decision['size_robust']}** |
| Matched-null specificity | **{decision['matched_null_specificity']}**; PASS fraction over N=50/75/100/150 primary score rows={decision['matched_null_pass_fraction']:.3f} |
| Original random-null diagnosis | **{null_diag}** |
| Leave-one-context-out | **{decision['loco']}** |
| Anchor removal | **{decision['anchor_removal']}** |
| Single-pathway dominance | **{decision['single_pathway_dominated']}** |
| Final representation | **{decision['final_representation']}** |

## Signature construction audit

The four 100-gene results were primarily cap-driven. Adaptive, fibrotic and
steroid signatures first applied FDR<0.05 and then retained ranked top-N genes.
The maladaptive core first required concordant FDR<0.05 support in both TGFβ2
and DEX, then applied the same 100-gene truncation. Disease-context evidence
was not used to determine core membership. See
`results/phase1R/SIGNATURE_CONSTRUCTION_AUDIT.tsv`.

{cap_rows.to_markdown(index=False)}

## Size sensitivity and transportability

N=25, 50, 75, 100, 150, 200 and RANKED_FULL variants were reconstructed from
the frozen external statistics. Primary scores use the same signed,
within-cell-state standardized logCPM and donor-pseudobulk procedure at every
size. The predeclared size gate required positive adaptive and maladaptive
direction with donor-positive fraction ≥0.5 in BeamA, JCT and Schlemm for
N=50/75/100/150. The result is **{decision['size_robust']}**; the machine-readable
table contains all external projections, primary deltas and generic-stress
overlap controls.

## Matched null calibration

Each adaptive and maladaptive-core size, including RANKED_FULL, was compared
with {10000:,} matched random sets. Matching used mean logCPM, detection rate,
variance and gene-set size from the exact GSE313736 frame; baseline TM
detectability is therefore represented. Gene length and GC content were
unavailable and were not imputed. Null sets received the observed up/down sign
composition and were processed by the same pseudobulk, z-score, signed-score,
donor-pair and ΔMSI calculations. Direction-consistent empirical P values use
`(1 + null >= observed)/(1 + null sets)`.

The original Phase 1 flag was an uncalibrated sparse overlap test using a small
uniform external-gene background. The full explanation is in
`reports/phase1R/RANDOM_NULL_FLAG_DIAGNOSIS.md`, with required conclusion
**{null_diag}**.

## Pathological context support and leave-one-context-out

GSE202372 (TGFβ2) and GSE301517 (DEX) are the two discovery pathological
contexts. GSE263230 is a limited descriptive glaucoma-context resource with
two FDR-significant genes and is not equivalent evidence. Among the frozen
context labels, {len(core_context)} genes meet the TGFβ2+DEX concordance rule.
M1/M2 test independently derived TGFβ2 and DEX programs in the opposite
external context; M3 tests the consensus descriptively in disease context; M4
and M5 apply the independent TGFβ2-only and DEX-only programs to GSE313736.
The results are in `LEAVE_ONE_CONTEXT_OUT_VALIDATION.tsv`.

## Rank concordance and adaptive–maladaptive separation

Pairwise Spearman effect/rank concordance, directional concordance and a
hypergeometric rank-overlap equivalent are reported across all four external
datasets in `CROSS_DATASET_RANK_CONCORDANCE.tsv`. The adaptive versus
maladaptive overlap table classifies shared genes as same-direction,
opposite-direction or context-dependent rather than treating overlap itself
as invalid evidence. The N=100 and ranked-full overlap results are in
`ADAPTIVE_MALADAPTIVE_OVERLAP.tsv`.

## Pathway, anchor and dominance controls

The pathway analysis uses pre-specified modules aligned to GO/Reactome/Hallmark
concepts: focal adhesion, integrin, actomyosin, ECM, TGFβ, Hippo/YAP, calcium,
oxidative stress, ferroptosis, autophagy, proteostasis and mechanosensitive
channels. No pathway was tuned on GSE313736. Anchor removal was performed for
PIEZO1/TRPV4 and for the full mechanosensory anchor panel. Sequential removal
of ECM, TGFβ, actomyosin, oxidative-stress and Hippo/YAP modules is reported
in `SINGLE_PATHWAY_DOMINANCE.tsv`; a broad state claim is allowed only when no
single module collapses the predeclared direction. In this run, the 100-gene
maladaptive core had zero overlap with each of those five module panels, so the
leave-one-module-out result is a non-informative limitation rather than proof
of pathway independence.

## Donor bootstrap and subtype heterogeneity

Bootstrap resamples were of the 59 donors, not cells, with B={BOOTSTRAPS}.
The summaries below are descriptive stability estimates:

{chr(10).join(lines)}

BeamA and JCT were primary; Schlemm_Endothelium was secondary. Cell states
were never pooled to improve significance. See
`CELLTYPE_PROGRAM_HETEROGENEITY.tsv`.

## Final representation and terminology

Final representation: **{decision['final_representation']}**. Allowed
terminology: **{decision['allowed_terminology']}**. The phrase
`adaptive–maladaptive state axis` is reserved for a robustness gate that
passes size, null, context, anchor and dominance criteria; when that gate is
not met, use the narrower program/module terminology above.

## Phase 2 boundary

Phase 2 allowed: **{decision['phase2_allowed']}**. Exact proposed scope:
**{decision['phase2_scope']}**. The proposal is written in
`reports/phase1R/PHASE2_PREREGISTRATION.md`; no ATAC/CUT&Tag/Hi-C analysis was
executed here. GSE301517 RNA-seq and GSE301516 ATAC-seq, GSE301519 CUT&Tag and
GSE301522 Hi-C are one GSE301525 experimental framework. Any future use is
regulatory anchoring, not independent biological replication.

## Remaining limitation

{decision['biggest_limitation']}

## Reproducibility

Script: `scripts/09_phase1R_robustness.py`. The primary exact-frame aggregation
used {primary['exact_cells']:,} cells, {len(primary['groups'])} donor×treatment×cell-state
groups and {len(primary['gene_list']):,} canonical gene symbols. All result
tables, the null diagnosis, the preregistration boundary and the JSON manifest
are stored under `results/phase1R/` and `reports/phase1R/`. No Phase 2 result is
presented as executed.
"""
    (REPORT / "PHASE1R_SIGNATURE_ROBUSTNESS_REPORT.md").write_text(report, encoding="utf-8")


def terminal_summary(decision: dict[str, object], bootstrap: pd.DataFrame, heterogeneity: pd.DataFrame, dominance: str, null_diag: str) -> None:
    def ci(ct: str) -> str:
        r = bootstrap[(bootstrap["celltype"].eq(ct)) & (bootstrap["metric"].eq("MSI"))]
        if r.empty:
            return "NA"
        x = r.iloc[0]
        return f"ΔMSI={x['observed_median_delta']:.4f}; 95% CI=[{x['ci95_low']:.4f}, {x['ci95_high']:.4f}]"
    print(f"1 PHASE 1R VERDICT: {decision['verdict']}")
    print(f"2 original 100-gene signatures cap-driven: {decision['top_n_dependence']}")
    print(f"3 size robustness: {decision['size_robust']}")
    print(f"4 matched-null specificity: {decision['matched_null_specificity']}")
    print(f"5 original random-null FLAG: {'RESOLVED' if null_diag.startswith('RESOLVED') else 'UNRESOLVED'}")
    print(f"6 leave-one-context-out: {decision['loco']}")
    print(f"7 TGFβ2-only → GSE313736 concordance: {decision['loco_tgf_primary']}")
    print(f"8 DEX-only → GSE313736 concordance: {decision['loco_dex_primary']}")
    print(f"9 adaptive-vs-maladaptive rank separation: {decision['rank_separation']}")
    print("10 generic-stress dependence: LOW (Phase 1 frozen classification; Phase 1R controls retained)")
    print(f"11 PIEZO1/TRPV4 removal sensitivity: {decision['anchor_removal']}")
    print(f"12 single-pathway dominance: {decision['single_pathway_dominated']}")
    print(f"13 BeamA bootstrap {ci('BeamA')}")
    print(f"14 JCT bootstrap {ci('JCT')}")
    print(f"15 Schlemm bootstrap {ci('Schlemm_Endothelium')}")
    print(f"16 final representation: {decision['final_representation']}")
    print(f"17 allowed terminology: {decision['allowed_terminology']}")
    print(f"18 Phase2 allowed: {decision['phase2_allowed']}")
    print(f"19 exact Phase2 scope: {decision['phase2_scope']}")
    print(f"20 biggest remaining limitation: {decision['biggest_limitation']}")
    print("21 path to report: reports/phase1R/PHASE1R_SIGNATURE_ROBUSTNESS_REPORT.md")


def main() -> int:
    RUN_LOG.write_text("", encoding="utf-8")
    rng = np.random.default_rng(SEED)
    frames = load_locked_external()
    core_all = consensus_full(frames)
    locked_core = clean_frame(ROOT / "results" / "phase1" / "MALADAPTIVE_CORE.tsv")
    log(f"Reconstructed untruncated TGFβ2+DEX consensus n={len(core_all)}; locked Phase 1 core n={len(locked_core)}")
    if len(core_all) >= 100 and not set(make_variants(core_all)["N=100"]["gene"]) == set(locked_core["gene"]):
        log("WARNING: reconstructed N=100 consensus differs from locked core; retaining locked external output for Phase 1 and reporting reconstruction")
    audit = audit_signature_construction(frames, core_all)
    context = context_support(frames)
    primary = primary_full_background()
    size_result, variants, scores = build_size_sensitivity(frames, core_all, primary, rng)
    null_result, null_meta = matched_null_tests(frames, variants, primary, np.random.default_rng(SEED + 1000))
    null_diag = random_null_diagnosis(frames, core_all, null_result, size_result)
    loco, loco_gates = loco_validation(frames, core_all, primary)
    rank = rank_concordance(frames)
    overlap = overlap_table(frames, core_all)
    pathway = pathway_robustness(frames, core_all, primary, variants)
    anchor_result, anchor_pass = anchor_removal(frames, core_all, primary, scores)
    dominance_table, dominance = single_pathway_dominance(frames, core_all, primary, scores)
    bootstrap, bootstrap_summary = donor_bootstrap(scores, primary)
    heterogeneity = celltype_heterogeneity(scores, primary)
    decision = final_decision(size_result, null_result, loco, loco_gates, anchor_pass, dominance, rank)
    final_signature_model(decision, frames, core_all)
    phase2_preregistration(decision)
    write_report(decision, audit, context, overlap, pathway, bootstrap, heterogeneity, null_diag, primary, rank)
    manifest = {
        **decision,
        "seed": SEED,
        "null_sets_per_family_size": NULL_N,
        "bootstrap_B": BOOTSTRAPS,
        "exact_frame_cells": primary["exact_cells"],
        "primary_groups": len(primary["groups"]),
        "canonical_genes": len(primary["gene_list"]),
        "null_background_candidates": null_meta["candidates"],
        "null_matching_strata": null_meta["pools"],
        "phase1_outputs_untouched": True,
        "phase2_executed": False,
        "super_series_rule": "GSE301517/GSE301516/GSE301519/GSE301522 are one GSE301525 framework; future omics use is regulatory anchoring, not independent biological replication",
    }
    (OUT / "PHASE1R_RUN_MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    terminal_summary(decision, bootstrap, heterogeneity, dominance, null_diag)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
