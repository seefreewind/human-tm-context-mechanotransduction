# Phase 2R.3 — final ATAC QC recovery and project closure

Date: 2026-09-27
Project: PIEZO1/TRPV4 context-dependent mechanotransduction pathway programs in glaucoma

## Final phase verdict

**TECHNICALLY_UNRELIABLE_ATAC**

ATAC classification: **TECHNICALLY_UNRELIABLE_ATAC**

## Raw-source audit and recovery

- Samples audited: 6/6.
- Initial local source status: ALIGNMENT_LOG;FINAL_BAM_ONLY;SRA;ALIGNMENT_LOG;SRA;AUDITED_FASTQ_MAP;ALIGNMENT_LOG.
- Historical `SRR34339248.phase2r3.name.bam` was not found or verifiable: `HISTORICAL_BAM_UNRECOVERABLE`.
- The retained, audited run3 paired FASTQ was registered by direct source-path mapping and used without copying, overwriting, SRA conversion or public re-download.
- Samples without retained pre-filter/raw inputs remain explicitly `FINAL_BAM_ONLY_INCOMPLETE`; their final-BAM TSS/FRiP/mapping values were not substituted for missing duplicate/chrM recovery metrics.
- Regenerated inputs and intermediates were retained for audit and safe resume; no old final BAM, SRA cache, run1/run2/run3 object or historical log was deleted.

Source audit: `results/phase2R3/ATAC_QC_RECOVERY_SOURCE.current.tsv`
Retention audit: `results/phase2R3/ATAC_QC_RECOVERY_RETENTION.tsv`

## Final technical QC

- PASS: 0; FAIL: 0; INCOMPLETE: 6.
- Locked thresholds: mapping rate ≥0.50; usable nuclear fragments ≥1,000,000; duplicate fraction ≤0.50; chrM fraction ≤0.20.
- Duplicate fraction recovered for 2/6; chrM fraction recovered for 2/6; TSS recovered for 0/6; FRiP recovered for 0/6.
- TSS and FRiP values are carried from the completed final-BAM technical audit; raw/pre-filter values are not substituted for them.

Final QC: `results/phase2R3/ATAC_FINAL_QC.tsv`
Balance audit: `results/phase2R3/ATAC_FINAL_QC_BALANCE.tsv`

## DAR validity

- Existing locked full-set DARs retained: 18342 (opening 6839; closing 11503).
- No differential-accessibility rerun was used to replace the primary result.
- QC-pass sensitivity was considered only after design-rank checking; a non-identifiable subset was not analyzed.

DAR validity: `results/phase2R3/ATAC_FINAL_DAR_VALIDITY.tsv`

## Module and regulatory conclusion

- Frozen matched-null modules passing FDR<0.05: 0/9.
- Final conclusion: **NO_MODULE_SPECIFIC_REGULATORY_SUPPORT**.
- Regulatory gate: **STOP_REGULATORY_ESCALATION**.
- Genetic gate: **STOP_GENETIC_ESCALATION**.
- Raw CUT&Tag, raw Hi-C, GWAS, eQTL, colocalization, PRS and MR remain closed/prohibited.

## Manuscript route

Use the context-dependent mechanotransduction pathway-program framing. Given the final QC classification, retain ATAC only as supplementary exploratory/audit evidence and do not make it a primary mechanistic pillar. Explicitly state that module-specific regulatory enrichment is absent.

Evidence hierarchy: `reports/phase2R3/FINAL_MANUSCRIPT_EVIDENCE_HIERARCHY.md`
