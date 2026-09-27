# Context-dependent mechanotransduction programs in human trabecular meshwork

This repository contains the analysis code and selected derived tables supporting a public-data study of context-dependent mechanotransduction pathway programs in human trabecular meshwork (TM). It is a curated research release, not a dump of the full working directory.

## Evidence snapshot

- Acute cyclic stretch, TGFβ2, dexamethasone-labelled, and limited glaucoma-context datasets showed context-specific source responses. Robustness analyses did not support a universal fixed adaptive–maladaptive gene signature; the locked representation is pathway modules.
- GSE313736 is used for descriptive donor-level application. Treatment labels are perfectly aliased with pool/library, so donor-paired score concordance is not a causal treatment effect.
- PIEZO1 and TRPV4 are treated as mechanosensory anchors, not established causal drivers.
- The final ATAC assessment is `TECHNICALLY_UNRELIABLE_ATAC`: all 6 samples are incomplete under the locked QC requirements and 0/9 frozen modules pass matched-null FDR control. ATAC is therefore audit/supplementary evidence, not module-specific regulatory support.
- The GSE301525 RNA/ATAC/CUT&Tag/Hi-C assays belong to one experimental framework and are not independent biological replication.
- PRS, GWAS, eQTL, colocalization, MR, raw CUT&Tag and raw Hi-C escalation are outside the final project gate.

The machine-readable gate and detailed final QC report are in [`reports/phase2R3/`](reports/phase2R3/). The selected tables in [`data/derived/`](data/derived/) are project-derived outputs; they are not independent validation data.

## Repository contents

- `scripts/`: phase-ordered analysis, audit, and acquisition utilities. `04_audit_prs.py` is a provenance audit only; it does not perform a PRS association analysis.
- `config/`: project scope/gate configuration and accession registry.
- `data/source_maps/`: compact sample/source mapping.
- `data/derived/`: selected, small result tables used to inspect the locked conclusions. The local BAM-path column was omitted from the public ATAC QC table.
- `provenance/`: dataset, script, seed, numerical, and software-version records.
- `reports/phase2R3/`: final ATAC quality, project-gate, and manuscript-evidence reports.

## Source data and reproducibility

The source studies are available through the GEO accessions listed in [`config/dataset_registry.tsv`](config/dataset_registry.tsv), including GSE195756, GSE202372, GSE263230, GSE301525 and its assay-specific subseries, and GSE313736. Retrieve source files from GEO/SRA and the relevant reference files from their original providers; check each source's terms and checksums before analysis. Source FASTQ/SRA files, BAMs, large cell/gene matrices, reference-genome files, temporary alignments, and run logs are intentionally not redistributed here.

Scripts are numbered by project phase. In broad order, scripts 00–07 audit dataset/sample structure and project gates; 08–09 derive and assess external signatures; 10 performs same-framework multiomic anchoring; 11–13 contain historical rescue/closure stages; and 14 implements the final ATAC recovery/gate logic. The Phase 2R.1 R model is `scripts/12_phase2R1_edger.R`.

The released copies of scripts 12–14 use portable environment discovery (`PHASE2R1_ENV_PREFIX`, `CONDA_PREFIX`, or the active Python prefix) instead of the original workstation-specific executable paths. These portability edits do not alter the analysis thresholds or statistical logic.

This release is **not yet a clean-environment, end-to-end reproduction package**. The contemporaneous downstream Python/R/edgeR environment was not fully captured; see [`provenance/EER_SOFTWARE_VERSION_LOCK.tsv`](provenance/EER_SOFTWARE_VERSION_LOCK.tsv). Do not substitute current package versions and describe them as historical. The selected outputs are archived evidence; their inclusion does not claim that a clean rerun has reproduced them.

## Integrity and licensing

`SHA256SUMS.txt` records checksums for the files in this release (excluding the checksum file itself). No explicit reuse license has yet been selected by the authors. The repository is publicly viewable, but public visibility alone does not grant permission to reuse or redistribute its contents. A license and Zenodo DOI will be added only after the authors resolve the applicable rights and complete Zenodo linking/release steps.

## Citation

See [`CITATION.cff`](CITATION.cff). No DOI has been assigned yet. Cite the original GEO studies for source data; do not cite this repository as an independent biological validation cohort.
