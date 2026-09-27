#!/usr/bin/env python3
"""Create a byte-range index for the uncompressed remote GEO TAR archive."""
from __future__ import annotations
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; META=ROOT/"data/metadata"; FILELIST=META/"GSE313736_filelist.txt"; OUT=META/"GSE313736_remote_tar_index.json"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313736/suppl/GSE313736_RAW.tar"
def main():
 rows=[]
 for line in FILELIST.read_text().splitlines():
  a=line.split("\t")
  if a and a[0]=="File": rows.append((a[1],int(a[3])))
 offset=0; members=[]
 for name,size in rows:
  members.append({"name":name,"header_offset":offset,"data_offset":offset+512,"size":size})
  offset += 512 + ((size+511)//512)*512
 OUT.write_text(json.dumps({"url":URL,"total_bytes":6536949760,"members":members},indent=2))
 print(f"members={len(members)} computed_archive_bytes={offset}")
if __name__=="__main__":main()
