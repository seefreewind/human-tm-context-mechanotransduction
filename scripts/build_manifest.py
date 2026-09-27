#!/usr/bin/env python3
"""Build the provenance/download manifest for all Phase 0 inputs."""
from __future__ import annotations
import csv,datetime,hashlib,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; META=ROOT/"data/metadata"; RAW=ROOT/"data/raw/GSE313736"; OUT=ROOT/"manifests/DATA_DOWNLOAD_MANIFEST.tsv"; DATE="2026-09-17"
def digest(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for x in iter(lambda:f.read(1024*1024),b""):h.update(x)
 return h.hexdigest()
def add(rows,p,url,acc,kind):
 if p.exists():rows.append({"source_url":url,"accession":acc,"filename":p.name,"filesize_bytes":p.stat().st_size,"md5":"","sha256":digest(p),"download_date":DATE,"local_path":str(p.relative_to(ROOT)),"file_kind":kind,"status":"DOWNLOADED"})
 else:rows.append({"source_url":url,"accession":acc,"filename":p.name,"filesize_bytes":"","md5":"","sha256":"","download_date":DATE,"local_path":str(p.relative_to(ROOT)),"file_kind":kind,"status":"NOT_DOWNLOADED"})
def main():
 rows=[]
 for p in sorted(META.glob("GSE*_family.soft.gz")): add(rows,p,f"https://ftp.ncbi.nlm.nih.gov/geo/series/{p.name[:6]}nnn/{p.name[:-len('_family.soft.gz')]}/soft/{p.name}",p.name[:-len('_family.soft.gz')],"GEO family SOFT")
 for p in [META/"GSE313736_HTMC_genotypes.vcf.gz",META/"GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",META/"GSE313736_HTMC_scRNA_demultiplexing_and_cell_classification.txt.gz",META/"GSE313736_HTMC_scRNA_donor_metadata.tsv.gz",META/"GSE313736_filelist.txt"]:
  name=p.name; url=f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313736/suppl/{name}"; add(rows,p,url,"GSE313736","primary supplementary metadata")
 for p in sorted(RAW.glob("*")):
  if p.name.startswith("._"): continue
  name=p.name; url=f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313736/suppl/GSE313736_RAW.tar#member={name}"; add(rows,p,url,"GSE313736","selected tar member")
 OUT.parent.mkdir(parents=True,exist_ok=True)
 with OUT.open("w",newline="") as f:
  fields=["source_url","accession","filename","filesize_bytes","md5","sha256","download_date","local_path","file_kind","status"]; w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(rows)
if __name__=="__main__":main()
