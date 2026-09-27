# Phase 2 allowed claims

Phase 2 verdict: **GO_PHASE3_CONDITIONAL**  
GSE301525 framework: GSE301517 RNA-seq, GSE301516 ATAC-seq, GSE301519 CUT&Tag and GSE301522 Hi-C/promoter-capture. These layers are orthogonal assays within one experimental framework, not independent biological replication.

## Allowed wording

- The locked transcriptomic modules show multiomic regulatory coherence within the common GSE301525 human trabecular-meshwork system.
- DEX-associated accessible chromatin and condition-associated histone/loop resources provide chromatin-level anchoring for selected modules.
- The result supports candidate regulatory elements and mechanosensory network anchors for follow-up experiments.
- A statement about the strongest supported modules must use the status table; current `REGULATORY_SUPPORTED` modules are: **none**.

## Prohibited wording

- Do not call GSE301516, GSE301519 or GSE301522 independent validation, replication, or independent cohorts.
- Do not call any prioritized interval a causal enhancer or prove that PIEZO1/TRPV4 causes the DEX response.
- Do not infer a causal DEX effect, PRS×DEX interaction, GWAS association, eQTL, colocalization, MR estimate, or clinical utility from Phase 2.
- Do not use GSE301512 as a DEX CUT&Run contrast; its metadata are resource-only.

## Technical boundary

The processed ATAC and CUT&Tag inputs are condition-level peak sets and the processed Hi-C inputs are condition-level compartment/loop summaries. Replicate-level counts, read depth and inferential FDR are not available in those files; corresponding fields remain `NA` and the report labels those results descriptive.
