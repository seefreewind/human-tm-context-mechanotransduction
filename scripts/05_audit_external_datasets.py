#!/usr/bin/env python3
"""Create a resource-level external dataset and regulatory modality audit."""
from __future__ import annotations
import csv,gzip,re
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; META=ROOT/"data"/"metadata"; OUT=ROOT/"results"/"phase0"
ACCS=["GSE195756","GSE202372","GSE301525","GSE301517","GSE301516","GSE301512","GSE301519","GSE301522","GSE263230","GSE309500"]

def value(block,prefix,default=""):
 for x in block:
  if x.startswith(prefix):return x.split("=",1)[1].strip()
 return default
def vals(block,prefix):return [x.split("=",1)[1].strip() for x in block if x.startswith(prefix)]
def parse(acc):
 lines=gzip.open(META/f"{acc}_family.soft.gz","rt",errors="replace").read().splitlines(); series=lines[:next((i for i,x in enumerate(lines) if x.startswith("^PLATFORM")),len(lines))]; ss=[]
 for i,x in enumerate(lines):
  if x.startswith("^SAMPLE = "):
   end=next((j for j in range(i+1,len(lines)) if lines[j].startswith("^")),len(lines)); b=lines[i:end]; ss.append({"acc":x.split("=",1)[1].strip(),"title":value(b,"!Sample_title = "),"chars":vals(b,"!Sample_characteristics_ch1 = "),"source":value(b,"!Sample_source_name_ch1 = "),"org":value(b,"!Sample_organism_ch1 = "),"supp":vals(b,"!Sample_supplementary_file")})
 return series,ss

def main():
 rows=[]; mod=[]
 for acc in ACCS:
  series,ss=parse(acc); title=value(series,"!Series_title = "); design=value(series,"!Series_overall_design = "); text=" ".join([title,design]+[x["title"]+" "+" ".join(x["chars"]) for x in ss]); low=text.lower()
  if acc=="GSE195756": assay="RNA-seq"; protocol="15% stretch; 0.69 Hz; 8 h"; controls="unstretched controls"; donors="4 paired biological samples inferred; donor IDs not verified"; match="PASS for human TM mechanical stretch resource, donor identity UNCLEAR"
  elif acc=="GSE202372": assay="RNA-seq"; protocol="TGF-beta2 exposure 24 h"; controls="vehicle/control"; donors="5 paired biological samples inferred from A-E; donor IDs not verified"; match="PASS for human TM TGF-beta2 perturbation; donor identity UNCLEAR"
  elif acc.startswith("GSE301"): assay=title.split("[")[-1].rstrip("]") if "[" in title else title; protocol="primary human TM cell lines; exact condition per sample metadata"; controls="resource/control/dex conditions where listed"; donors="cell lines/replicates; donor-level identity not verified"; match="PASS for modality audit; same-condition regulatory integration is condition/modality dependent"
  elif acc=="GSE263230": assay="10x scRNA-seq"; protocol="primary culture and post-mortem TM"; controls="culture/tissue comparison"; donors="14 primary culture samples + 3 post-mortem samples"; match="PASS for human TM cell identity; not a DEX validation set"
  else: assay="scRNA-seq + multiome"; protocol="mouse limbal strips; WT/Lmx1b V265D/+; vitamin B3 context"; controls="WT"; donors="pooled mouse eyes; individual donor N not applicable"; match="FAIL for human TM validation: mouse limbal strip, TM included but not isolated human TM"
  rows.append({"accession":acc,"title":title,"organism":next((x["org"] for x in ss if x["org"]),""),"sample_records":len(ss),"assay":assay,"donor_n_or_replicates":donors,"protocol_duration_amplitude_frequency":protocol,"control_design":controls,"processed_or_raw_files":"; ".join(vals(series,"!Series_supplementary_file")),"human_TM_match":match,"notes":design})
  if acc.startswith("GSE301") and acc!="GSE301525":
   conditions=[]
   for s in ss:
    tl=s["title"].lower(); cond="dex" if "dex" in tl else ("control" if "control" in tl else ("resource" if "resource" in tl else "unclear")); conditions.append(cond)
   mod.append({"accession":acc,"modality":assay,"sample_records":len(ss),"resource_n":conditions.count("resource"),"control_n":conditions.count("control"),"dex_n":conditions.count("dex"),"same_condition_regulatory_integration":"YES_CONDITIONALLY" if conditions.count("dex") and conditions.count("control") else "NO_OR_UNCLEAR","notes":"GSE301512 sample titles are resource-only/control-resource; do not infer DEX CUT&Run from parent series." if acc=="GSE301512" else "Sample titles support resource/control/dex mapping; donor/cell-line matching must be verified before integration."})
 pd.DataFrame(rows).to_csv(OUT/"EXTERNAL_DATASET_AUDIT.tsv",sep="\t",index=False); pd.DataFrame(mod).to_csv(OUT/"REGULATORY_MULTIOMICS_MAP.tsv",sep="\t",index=False)

if __name__=="__main__":main()
