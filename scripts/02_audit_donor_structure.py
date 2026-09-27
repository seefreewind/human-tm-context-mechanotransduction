#!/usr/bin/env python3
"""Compute donor/treatment structure and completeness for GSE313736."""
from __future__ import annotations
import csv, gzip, json, re
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; META=ROOT/"data"/"metadata"; OUT=ROOT/"results"/"phase0"

def main():
    cm=pd.read_csv(META/"GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",sep="\t",compression="gzip",low_memory=False)
    dm=pd.read_csv(META/"GSE313736_HTMC_scRNA_donor_metadata.tsv.gz",sep="\t",compression="gzip",low_memory=False)
    t=cm.groupby(["donor_id_demuxlet","DEX Treatment"]).size().unstack(fill_value=0)
    for c in ["no-dex","dex"]:
        if c not in t.columns:t[c]=0
    t["vehicle_present"]=t["no-dex"]>0; t["dex_present"]=t["dex"]>0
    d=dm.set_index("donor_id_demuxlet").copy(); d=d.join(t[["vehicle_present","dex_present","no-dex","dex"]],how="outer")
    d=d.rename(columns={"no-dex":"n_vehicle_cells","dex":"n_dex_cells"})
    d["vehicle_present"]=d["vehicle_present"].fillna(False).astype(bool); d["dex_present"]=d["dex_present"].fillna(False).astype(bool)
    d["genotype_present"]=d.index.isin(vcf_donors(META/"GSE313736_HTMC_genotypes.vcf.gz"))
    d["prs_present"]=d["glaucoma_prs_raw"].notna() & d["glaucoma_prs_centile"].notna()
    d["disease_status"]=d["Glaucoma"].map({0:"non-glaucoma",1:"glaucoma"}).fillna("missing/unknown")
    d["ancestry"]=d["inferred_donor_ancestry"].fillna("missing")
    d["sex"]=d["assigned_donor_sex"].map({1:"code_1",2:"code_2"}).fillna("missing")
    d["age"]=d["assigned_donor_age"]
    d["batch"]=cm.groupby("donor_id_demuxlet")["library_id"].agg(lambda x:";".join(sorted(set(x.astype(str))))).reindex(d.index)
    d["comments"]=""
    d.loc[~(d.vehicle_present&d.dex_present),"comments"]="unpaired donor in supplied cell metadata"
    d.loc[d.inferred_donor_ancestry.isna(),"comments"] += "; ancestry missing"
    d.loc[d.Glaucoma.eq(0),"comments"] += "; all supplied donor disease labels are non-glaucoma"
    cols=["vehicle_present","dex_present","n_vehicle_cells","n_dex_cells","genotype_present","prs_present","disease_status","ancestry","sex","age","batch","comments"]
    out=OUT/"DONOR_COMPLETENESS.tsv"; out.parent.mkdir(parents=True,exist_ok=True)
    dd=d.reset_index().rename(columns={"donor_id_demuxlet":"donor_id","index":"donor_id"}); dd.to_csv(out,sep="\t",index=False,columns=["donor_id"]+cols)
    ctab=cm.groupby(["donor_id_demuxlet","DEX Treatment","celltypist_chen2025_fine"]).size().reset_index(name="cells")
    ctab.to_csv(OUT/"CELLTYPE_COUNTS_BY_DONOR_TREATMENT.tsv",sep="\t",index=False)
    summary={
      "total_cells":int(len(cm)),"total_sample_ids":int(cm.sampleid.nunique()),"total_libraries":int(cm.library_id.nunique()),"total_donors":int(cm.donor_id_demuxlet.nunique()),
      "donors_vehicle":int(t.vehicle_present.sum()),"donors_dex":int(t.dex_present.sum()),"paired_vehicle_dex":int((t.vehicle_present&t.dex_present).sum()),"vehicle_only":int((t.vehicle_present&~t.dex_present).sum()),"dex_only":int((t.dex_present&~t.vehicle_present).sum()),
      "donors_genotype":int(d.genotype_present.sum()),"donors_prs":int(d.prs_present.sum()),"donors_disease_status":int(d.Glaucoma.notna().sum()),"donors_age":int(d.assigned_donor_age.notna().sum()),"donors_sex":int(d.assigned_donor_sex.notna().sum()),"donors_ancestry":int(d.inferred_donor_ancestry.notna().sum()),
      "prs_interaction_usable_donors":int((d.prs_present & d.vehicle_present & d.dex_present).sum()),"prs_raw_range":[float(d.glaucoma_prs_raw.min()),float(d.glaucoma_prs_raw.max())],"prs_centile_range":[float(d.glaucoma_prs_centile.min()),float(d.glaucoma_prs_centile.max())],
      "disease_status_counts":{str(k):int(v) for k,v in d.Glaucoma.value_counts(dropna=False).items()},"ancestry_counts":{str(k):int(v) for k,v in d.inferred_donor_ancestry.value_counts(dropna=False).items()},"sex_code_counts":{str(k):int(v) for k,v in d.assigned_donor_sex.value_counts(dropna=False).items()},
      "treatment_cells":{str(k):int(v) for k,v in cm["DEX Treatment"].value_counts().items()},"cell_type_coarse_counts":{str(k):int(v) for k,v in cm["celltypist_chen2025_coarse"].value_counts(dropna=False).items()},"cell_type_fine_counts":{str(k):int(v) for k,v in cm["celltypist_chen2025_fine"].value_counts(dropna=False).items()},
      "treatment_library_counts":{str(k):int(v) for k,v in cm.groupby("DEX Treatment").library_id.nunique().items()},"pool_donor_counts":cm.groupby(["DEX Treatment","library_id"]).donor_id_demuxlet.nunique().unstack(fill_value=0).to_dict(),
    }
    (OUT/"DONOR_STRUCTURE_SUMMARY.json").write_text(json.dumps(summary,indent=2,default=str))

def vcf_donors(path):
    with gzip.open(path,"rt",errors="replace") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return {x.split("_",1)[1] for x in line.rstrip().split("\t")[9:] if "_" in x}
    return set()

if __name__=="__main__": main()
