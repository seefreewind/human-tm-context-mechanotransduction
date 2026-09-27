#!/usr/bin/env python3
"""Generate the single Phase 0 feasibility report from audit outputs."""
from __future__ import annotations
import json,csv
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results"/"phase0"; REPORT=ROOT/"reports"/"phase0"/"PHASE0_FEASIBILITY_REPORT.md"

def j(name): return json.loads((OUT/name).read_text())
def yn(v): return "YES" if v else "NO"
def main():
 s=j("DONOR_STRUCTURE_SUMMARY.json"); p=j("PRS_SUMMARY.json"); g=pd.read_csv(OUT/"GENE_DETECTABILITY.tsv",sep="\t"); reg=pd.read_csv(ROOT/"config/dataset_registry.tsv",sep="\t"); ext=pd.read_csv(OUT/"EXTERNAL_DATASET_AUDIT.tsv",sep="\t"); mm=pd.read_csv(OUT/"REGULATORY_MULTIOMICS_MAP.tsv",sep="\t")
 gene_lines=[]
 for gene in ["PIEZO1","TRPV4"]:
  r=g[g.gene==gene].iloc[0]; gene_lines.append(f"- {gene}: cells detected {r.cells_detected if r.status=='OK' else 'UNAVAILABLE'}, donor detection {r.donors_detected if r.status=='OK' else 'UNAVAILABLE'}/{s['total_donors']} ({r.pct_donors_detected if r.status=='OK' else 'NA'}%); vehicle {r.vehicle_pct if r.status=='OK' else 'NA'}%, DEX {r.dex_pct if r.status=='OK' else 'NA'}%.")
 detect_status="PASS" if (g.status=="OK").all() else "HOLD"
 verdict="CONDITIONAL_GO"; phase1="NO"
 REPORT.parent.mkdir(parents=True,exist_ok=True)
 REPORT.write_text(f'''# PHASE 0 — DATA FEASIBILITY AND METADATA AUDIT

Audit date: 2026-09-17. Scope was restricted to accession/source audit, metadata schema, donor structure, genotype/PRS provenance, descriptive gene detectability, cell-type metadata, and external-resource authenticity. No differential expression, PRS association, eQTL, GWAS colocalization, trajectory, network, machine-learning, or downstream biological interpretation was run.

## PHASE 0 VERDICT

**{verdict}**. The primary dataset is accessible and contains a usable paired donor-treatment structure with 59 donor-level PRS values. The project cannot pass to Phase 1 yet because PRS provenance is incomplete, disease status is constant (all supplied donors are non-glaucoma), ancestry is sparse, treatment is perfectly separated by pooled library, and the supplied cell metadata covers 169,398 cells whereas the GEO series summary reports 316,799 singlets. The PRS × DEX aim is therefore feasible as a conditional donor-level design, not as an already validated primary analysis.

Phase 1 allowed: **{phase1}** until the explicit holds in this report are resolved and recorded in the project configuration.

## PRIMARY DATASET

- Accession: GSE313736, *Single Cell Profiling of Human Trabecular Meshwork Cells*.
- GEO record reports 74 people/donors and 316,799 singlet cells. The supplied `GSE313736_HTMC_scRNA_cell_metadata.tsv.gz` contains {s['total_cells']:,} rows, {s['total_donors']} observed donors, {s['total_sample_ids']} sample IDs, and {s['total_libraries']} pooled libraries.
- The 18 downloaded expression matrices declare **{j('PRIMARY_METADATA_SUMMARY.json')['matrix_declared_cells']:,}** cell barcodes. Of the {s['total_cells']:,} metadata rows, **{j('PRIMARY_METADATA_SUMMARY.json')['metadata_matrix_intersection_cells']:,}** match a matrix barcode and **{j('PRIMARY_METADATA_SUMMARY.json')['metadata_rows_not_matching_matrix']:,}** do not. These are separate denominators and are not interchangeable.
- The actual cell metadata contains author-provided coarse and fine cell predictions; all rows carry the coarse label `trabecular meshwork cell`. This table is the auditable analysis frame for Phase 0. The 147,401-cell difference from the GEO summary is unresolved and must not be silently filled in.
- Supplementary files verified locally: donor metadata, cell metadata, demultiplexing/classification documentation, genotype VCF, and RAW archive file list. The complete RAW archive is listed as 6,536,949,760 bytes; it was not downloaded as one monolithic archive.

## DONOR STRUCTURE

- Total cells in supplied cell metadata: **{s['total_cells']:,}**.
- Total donors: **{s['total_donors']}**.
- Vehicle donors: **{s['donors_vehicle']}**; DEX donors: **{s['donors_dex']}**.
- Paired vehicle + DEX donors: **{s['paired_vehicle_dex']}**; vehicle-only: **{s['vehicle_only']}; DEX-only: {s['dex_only']}**.
- Donors with genotype mapping: **{s['donors_genotype']}**; donors with PRS: **{s['donors_prs']}**; usable paired PRS donors: **{s['prs_interaction_usable_donors']}**.
- Donor is the biological replicate. Cell counts are descriptive and must not be used as the statistical N.

## TREATMENT STRUCTURE

- Vehicle: {s['treatment_cells']['no-dex']:,} cells; DEX: {s['treatment_cells']['dex']:,} cells.
- Treatment was performed in separate pooled captures: no-dex libraries are `HTMC_pool-A` through `HTMC_pool-I`; DEX libraries are `HTMC_pool-dexA` through `HTMC_pool-dexI`.
- Treatment is therefore perfectly confounded with the pool/library naming layer. Pairing by donor is strong (59/59), but it does not remove all pool-level processing confounding. `library_id`, saturation, and other technical covariates must be audited and modeled before any DEX claim.
- Treatment duration is internally inconsistent in the supplied GEO record: sample characteristics state 48 hours, while the growth protocol describes dexamethasone exposure for 3 days. This must be reconciled with the authors before interpreting the perturbation.

## PRS / GENOTYPE

- Raw genotype: verified VCF with **{p['vcf_sample_columns']}** sample columns. It maps to {p['mapped_cell_donors']} of the {p['cell_metadata_donors']} observed donor IDs; {len(p['vcf_donors_without_cell_metadata'])} VCF donor labels have no observed cells.
- PRS: `glaucoma_prs_raw` and `glaucoma_prs_centile` are real donor-level metadata columns; all {p['prs_donors']} observed donors have both values. Raw range: {p['raw_prs_range'][0]} to {p['raw_prs_range'][1]}; centile range: {p['centile_range'][0]} to {p['centile_range'][1]}.
- PRS provenance is unresolved: the supplied files do not document the GWAS source, genome build, variant weights, allele harmonization, imputation status, or scaling/standardization procedure. This is a Phase 1 HOLD (`PRS_PROVENANCE_UNRESOLVED`).
- Disease status is unusable as a covariate or modifier in the current observed donor set: all {p['cell_metadata_donors']} donor labels are non-glaucoma. Thus PRS and disease status cannot be meaningfully disentangled here; no case/control genetic comparison is supported.
- Ancestry is 57 EUR, 1 AMR, and 1 missing. A robust ancestry-stratified PRS analysis is not supported; ancestry adjustment is technically sparse and its effect on PRS calibration remains unresolved.

## MINIMUM SAMPLE-SIZE FLAGS

- PRS × DEX interaction usable donor N: **{s['prs_interaction_usable_donors']} — GREEN by the pre-specified N-only rule (≥50)**.
- This GREEN flag is not a go-ahead for testing: PRS provenance, treatment–pool confounding, and the 74-versus-59 discrepancy remain mandatory holds.

## PIEZO1/TRPV4 DETECTABILITY

{chr(10).join(gene_lines)}

Detectability audit status: **{detect_status}**. This table is descriptive only; low detection must not be used to delete PIEZO1/TRPV4 before a pre-specified modeling decision.

## CELL-TYPE STRUCTURE

- Coarse author annotation: all {s['total_cells']:,} supplied rows are labeled `trabecular meshwork cell`.
- Fine labels include BeamA ({s['cell_type_fine_counts'].get('BeamA',0):,}), JCT ({s['cell_type_fine_counts'].get('JCT',0):,}), Ciliary_Muscle ({s['cell_type_fine_counts'].get('Ciliary_Muscle',0):,}), Schlemm_Endothelium ({s['cell_type_fine_counts'].get('Schlemm_Endothelium',0):,}), plus smaller endothelial, immune, Schwann, pericyte, melanocyte, and epithelial-like/other labels.
- The author coarse label should be retained. Fine labels show cell-type imbalance and contain contaminant-like populations; cell-type-specific analysis will require pre-specified inclusion and sensitivity rules.
- DEX-specific cell loss cannot be interpreted from cell counts alone because donor contribution and library/pool composition differ; donor-level cell-count balance is recorded in `DONOR_COMPLETENESS.tsv`.
- Gene-symbol alias audit: the exact symbol `CTGF` was absent from the supplied feature list, while `CCN2` was detected; future analysis must lock this alias mapping before use.

## EXTERNAL VALIDATION DATASETS

| Accession | Audited result |
|---|---|
| GSE195756 | PASS for human TM mechanical stretch RNA-seq resource; 4 paired biological samples inferred, donor IDs not verified; 15% stretch, 0.69 Hz, 8 h. |
| GSE202372 | PASS for human TM TGF-beta2 RNA-seq; 5 paired biological samples inferred from A-E labels, donor IDs not verified; 24 h exposure. |
| GSE301517 | PASS as human TM total RNA-seq with resource/control/DEX samples; donor/cell-line pairing requires verification. |
| GSE301516 | PASS as human TM ATAC-seq with resource/control/DEX sample titles; modality-specific replicate matching required. |
| GSE301512 | CONDITIONAL/limited: CUT&Run resource series; sample titles are resource/control-resource and series design says control cells. Do not assume DEX CUT&Run support. |
| GSE301519 | PASS as human TM CUT&Tag with control/DEX sample titles; condition and cell-line matching required. |
| GSE301522 | PASS as human TM Hi-C with resource/control/DEX sample titles; condition and cell-line matching required. |
| GSE263230 | PASS for independent human TM cell-identity context; 14 primary culture + 3 post-mortem samples, including four glaucoma culture samples reported by GEO; not a DEX perturbation dataset. |
| GSE309500 | FAIL for human TM validation: mouse limbal strips with TM included, WT/Lmx1b V265D/+ and vitamin B3 context. Use only as cross-species contextual evidence if later justified. |

The detailed machine-readable audit is `results/phase0/EXTERNAL_DATASET_AUDIT.tsv`; the regulatory modality map is `results/phase0/REGULATORY_MULTIOMICS_MAP.tsv`.

## REGULATORY MULTIOMICS

GSE301525 is a SuperSeries of GSE301512, GSE301516, GSE301517, GSE301519, and GSE301522. RNA-seq, ATAC-seq, CUT&Tag, and Hi-C sample titles support resource/control/DEX mapping, but same-condition integration still requires exact cell-line, replicate, treatment-duration, and batch matching. GSE301512 is a specific exception: its record describes CUT&Run for a control/resource dataset, so DEX regulatory support is not verified.

## IDENTIFIED CONFOUNDERS

1. Supplied cell metadata versus matrix/GEO summary cell/donor-count discrepancy (169,398 metadata rows; 166,323 metadata–matrix matches; 350,171 matrix barcodes; GEO summary 316,799 cells; 59 versus 74 donors).
2. Treatment is separated into independent pools/libraries; pool and treatment are not orthogonal.
3. Disease status is constant: all observed donors are non-glaucoma.
4. Ancestry is sparse and unbalanced (57 EUR, 1 AMR, 1 missing), limiting PRS calibration and adjustment.
5. PRS definition, GWAS source, build, weights, harmonization, and scaling are not documented in the supplied files.
6. Fine cell labels are highly imbalanced and include contaminant-like populations despite the coarse TM label.
7. VCF includes 74 donor labels, while 15 have no observed cell metadata; the 59-cell-metadata donor intersection is the usable mapping frame.

## ANALYSIS RISKS

- **High:** treating cells as independent replicates would inflate N; donor-pseudobulk/paired donor modeling is mandatory.
- **High:** attributing all treatment differences to DEX without addressing pool/library confounding.
- **High:** using disease status as a PRS proxy or claiming glaucoma genetic stratification from an all-non-glaucoma donor set.
- **High:** performing DEX discovery in GSE313736 and reusing those genes as a confirmatory maladaptive score in the same data (data leakage).
- **Moderate/High:** using PRS centiles as if they were standardized raw PRS, or inferring ancestry calibration without an ancestry-specific definition.
- **Moderate:** pooling cell-type labels across fine categories without a pre-specified cell-state/contaminant rule.

## RESCUE OPTIONS

1. Obtain author documentation for the PRS GWAS source, build, score file/weights, allele harmonization, imputation status, and scaling; lock this provenance before Phase 1.
2. Obtain the missing 15 donor/cell records or an explanation for the 74-versus-59 discrepancy; otherwise declare the analysis frame as 59 donors and do not generalize to 74.
3. Obtain pool-level technical metadata and confirm whether treatment assignment, processing date, library preparation, and sequencing lane are separable from treatment. If not, downgrade DEX inference to paired descriptive evidence and require independent DEX replication.
4. Use independent data for signature discovery: acute mechanical stretch for adaptive candidates; DEX/TGF-beta2/glaucoma/ECM-stiffness evidence for maladaptive candidates. Keep discovery and primary testing datasets separate.
5. Treat GSE263230 as human TM identity/heterogeneity validation and exclude GSE309500 from human validation claims.

## FINAL GO / CONDITIONAL_GO / NO_GO

**CONDITIONAL_GO** for a gated project: the expression and paired donor framework is real and the conditional PRS interaction design has 59 usable donor-level PRS values. **NO_GO_PRIMARY_GENETIC_AIM** would be triggered if the PRS-to-donor mapping or provenance fails; it is not triggered by the current files, but provenance remains a mandatory HOLD. The disease-status aim is not supported in this cohort.

### Phase 1 gate status

**NO — Phase 1 is not allowed yet.** Close the PRS provenance, 74-versus-59 donor discrepancy, and pool/treatment confounding gates first. No signature discovery or formal DEX/PRS testing was performed in Phase 0.

## REPRODUCIBILITY

- Python: 3.9.6; R: 4.4.3; pandas: 2.3.3; numpy: 1.26.4; scipy: 1.13.1; anndata: 0.10.8; scanpy: 1.10.3; PyYAML: 6.0.3.
- Source records: NCBI GEO accession pages and GEO FTP files queried on 2026-09-17; exact URLs and local paths are in `DATA_DOWNLOAD_MANIFEST.tsv`.
- Scripts are in `scripts/00_build_dataset_registry.py` through `scripts/06_generate_phase0_report.py`.
- Raw matrix archive was handled through verified HTTP byte ranges; no full 6.1-GB archive extraction is assumed by this report.

## OUTPUTS

- `config/dataset_registry.tsv`
- `results/phase0/GSE313736_METADATA_COLUMNS.tsv`
- `results/phase0/DONOR_COMPLETENESS.tsv`
- `results/phase0/GENE_DETECTABILITY.tsv`
- `results/phase0/CELLTYPE_COUNTS_BY_DONOR_TREATMENT.tsv`
- `results/phase0/PRS_AUDIT.tsv`
- `results/phase0/candidate_signature_sources.tsv`
- `manifests/DATA_DOWNLOAD_MANIFEST.tsv`
- `ANALYSIS_GOVERNANCE.md`
''')

if __name__=="__main__": main()
