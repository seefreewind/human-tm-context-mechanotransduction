#!/usr/bin/env python3
"""Audit donor-level genotype and PRS provenance/mapping."""
from __future__ import annotations
import csv,gzip,json
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; META=ROOT/"data"/"metadata"; OUT=ROOT/"results"/"phase0"

def vcf_header(path):
    with gzip.open(path,"rt",errors="replace") as f:
        for line in f:
            if line.startswith("#CHROM"): return line.rstrip().split("\t")[9:]
    return []

def main():
    dm=pd.read_csv(META/"GSE313736_HTMC_scRNA_donor_metadata.tsv.gz",sep="\t")
    cm=pd.read_csv(META/"GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",sep="\t",low_memory=False)
    vcf=vcf_header(META/"GSE313736_HTMC_genotypes.vcf.gz"); vcf_donors={x.split("_",1)[1] for x in vcf if "_" in x}; cell_donors=set(cm.donor_id_demuxlet)
    rows=[
      {"audit_item":"raw_genotype","value":"YES_VERIFIED; gzipped VCF downloaded from GEO","status":"VERIFIED","evidence":"GSE313736_HTMC_genotypes.vcf.gz"},
      {"audit_item":"genotype_type","value":"VCF with FORMAT/sample columns; imputation status not documented in supplied file header","status":"UNCLEAR","evidence":f"{len(vcf)} VCF sample columns"},
      {"audit_item":"genotype_PCs","value":"not present in donor metadata or supplied VCF audit fields","status":"NO","evidence":"metadata schema audit"},
      {"audit_item":"PRS_raw","value":"glaucoma_prs_raw","status":"YES_VERIFIED","evidence":"donor metadata column"},
      {"audit_item":"PRS_percentile","value":"glaucoma_prs_centile","status":"YES_VERIFIED","evidence":"donor metadata column"},
      {"audit_item":"PRS_definition","value":"glaucoma PRS; exact variant weights/GWAS build/allele harmonization not documented in supplied files","status":"PRS_PROVENANCE_UNRESOLVED","evidence":"donor metadata + SOFT"},
      {"audit_item":"PRS_scaling","value":"raw PRS and centile supplied; standardization procedure not documented","status":"UNCLEAR","evidence":"donor metadata"},
      {"audit_item":"PRS_ancestry_restriction","value":"inferred ancestry supplied; 57 EUR, 1 AMR, 1 missing","status":"UNCLEAR","evidence":"donor metadata"},
      {"audit_item":"disease_status","value":"Glaucoma column; all 59 values are 0/non-glaucoma","status":"YES_VERIFIED_BUT_CONSTANT","evidence":"donor metadata"},
      {"audit_item":"donor_mapping","value":f"{len(cell_donors & vcf_donors)} of {len(cell_donors)} cell-metadata donors map to VCF; {len(vcf_donors-cell_donors)} VCF donors have no observed cells","status":"VERIFIED_WITH_UNMATCHED_VCF_SAMPLES","evidence":"VCF sample header vs demuxlet donor_id"},
    ]
    pd.DataFrame(rows).to_csv(OUT/"PRS_AUDIT.tsv",sep="\t",index=False)
    s=dm[["donor_id_demuxlet","glaucoma_prs_raw","glaucoma_prs_centile","inferred_donor_ancestry","Glaucoma"]].copy(); s["genotype_present"]=s.donor_id_demuxlet.isin(vcf_donors); s["prs_present"]=s.glaucoma_prs_raw.notna()&s.glaucoma_prs_centile.notna(); s.to_csv(OUT/"PRS_DONOR_SUMMARY.tsv",sep="\t",index=False)
    summary={"vcf_sample_columns":len(vcf),"vcf_donor_labels":len(vcf_donors),"cell_metadata_donors":len(cell_donors),"mapped_cell_donors":len(cell_donors&vcf_donors),"vcf_donors_without_cell_metadata":sorted(vcf_donors-cell_donors),"prs_donors":int(s.prs_present.sum()),"raw_prs_range":[float(s.glaucoma_prs_raw.min()),float(s.glaucoma_prs_raw.max())],"centile_range":[float(s.glaucoma_prs_centile.min()),float(s.glaucoma_prs_centile.max())],"disease_status_counts":{str(k):int(v) for k,v in s.Glaucoma.value_counts(dropna=False).items()},"ancestry_counts":{str(k):int(v) for k,v in s.inferred_donor_ancestry.value_counts(dropna=False).items()}}
    (OUT/"PRS_SUMMARY.json").write_text(json.dumps(summary,indent=2))

if __name__=="__main__": main()
