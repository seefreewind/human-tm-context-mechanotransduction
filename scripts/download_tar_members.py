#!/usr/bin/env python3
"""Download selected members from the remote GSE313736 tar archive by HTTP range.

The full 6.1-GB archive is intentionally not downloaded. Selected Matrix
Market members, features, and barcodes are stored under data/raw/GSE313736.
"""
from __future__ import annotations
import concurrent.futures, hashlib, json, os, sys
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]; META=ROOT/"data/metadata"; RAW=ROOT/"data/raw/GSE313736"; LOG=ROOT/"logs/phase0_matrix_download.tsv"
INDEX=META/"GSE313736_remote_tar_index.json"
URL=json.loads(INDEX.read_text())["url"]

def one(m):
    name=m["name"]; out=RAW/name; part=out.with_suffix(out.suffix+".part")
    if out.exists() and out.stat().st_size==m["size"]: return name,"EXISTS",m["size"],sha(out)
    RAW.mkdir(parents=True,exist_ok=True)
    try:
        # The GEO server can truncate a large 206 response before the requested
        # end byte. Use bounded chunks and retain a partial file for retry.
        chunk_size=64*1024*1024
        start=part.stat().st_size if part.exists() else 0
        if start>m["size"]: part.unlink(); start=0
        with part.open("ab") as f:
            while start<m["size"]:
                end=min(m["size"]-1,start+chunk_size-1)
                h={"Range":f"bytes={m['data_offset']+start}-{m['data_offset']+end}"}
                with requests.get(URL,headers=h,stream=True,timeout=(30,600)) as r:
                    r.raise_for_status()
                    if r.status_code!=206: raise RuntimeError(f"expected 206 range response, got {r.status_code}")
                    got=0
                    for chunk in r.iter_content(1024*1024):
                        if chunk: f.write(chunk); got+=len(chunk)
                if got<=0: raise RuntimeError("empty range response")
                start += got
                if start>m["size"]: raise RuntimeError(f"range response exceeded expected size: {start} > {m['size']}")
        if part.stat().st_size!=m["size"]: raise RuntimeError(f"size {part.stat().st_size} != expected {m['size']}")
        part.replace(out); return name,"DOWNLOADED",m["size"],sha(out)
    except Exception as e:
        return name,"ERROR:"+repr(e),0,""

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for x in iter(lambda:f.read(1024*1024),b""): h.update(x)
    return h.hexdigest()

def main():
    members=[m for m in json.loads(INDEX.read_text())["members"] if m["name"].endswith(("_matrix.mtx.gz","_features.tsv.gz","_barcodes.tsv.gz"))]
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        for r in ex.map(one,members): results.append(r)
    LOG.parent.mkdir(parents=True,exist_ok=True)
    with LOG.open("w") as f:
        f.write("filename\tstatus\tbytes\tsha256\n")
        for r in results:f.write("\t".join(map(str,r))+"\n")
    bad=[r for r in results if r[1].startswith("ERROR")]
    print(f"selected_members={len(results)} completed={len(results)-len(bad)} errors={len(bad)}")
    if bad:
        for r in bad: print(r,file=sys.stderr)
        raise SystemExit(1)

if __name__=="__main__":main()
