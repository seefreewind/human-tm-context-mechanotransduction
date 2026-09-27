#!/usr/bin/env python3
"""Audit candidate-gene detection from Cell Ranger Matrix Market members.

This is a descriptive detectability audit only; no expression comparison or
biological interpretation is performed.
"""
from __future__ import annotations
import csv, gzip, json, re, tempfile
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; META=ROOT/"data"/"metadata"; RAW=ROOT/"data"/"raw"/"GSE313736"; OUT=ROOT/"results"/"phase0"
GENES=["PIEZO1","ANGPT2","TEK","ITGA9","PTK2","TLN1","VCL","PXN","FERMT2","CAV1","CAV2","TRPV4","RHOA","ROCK1","ROCK2","ACTA2","TAGLN","MYL9","YAP1","WWTR1","CTGF","CCN2","TGFB2","FN1","COL1A1","COL3A1","GPX4","SLC7A11","ACSL4","TFRC","FTH1","HMOX1","NFE2L2"]

def main():
    cm=pd.read_csv(META/"GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",sep="\t",compression="gzip",low_memory=False)
    by_barcode={(str(row.library_id),str(row.barcode)):(str(row.donor_id_demuxlet),str(row["DEX Treatment"])) for _,row in cm.iterrows()}
    all_results={g:{"cells_detected":0,"donors":set(),"vehicle":0,"dex":0} for g in GENES}
    members=json.loads((META/"GSE313736_remote_tar_index.json").read_text())["members"]
    matrix_members=[m for m in members if m["name"].endswith("_matrix.mtx.gz")]
    if not matrix_members or not all((RAW/(m["name"])).exists() for m in matrix_members):
        rows=[{"gene":g,"status":"MATRIX_MEMBERS_NOT_DOWNLOADED","cells_detected":"NA","pct_cells_detected":"NA","donors_detected":"NA","pct_donors_detected":"NA","vehicle_pct":"NA","dex_pct":"NA"} for g in GENES]
        pd.DataFrame(rows).to_csv(OUT/"GENE_DETECTABILITY.tsv",sep="\t",index=False); return
    import scipy.io
    for m in matrix_members:
        stem=m["name"].replace("_matrix.mtx.gz","")
        feat=RAW/f"{stem}_features.tsv.gz"; bar=RAW/f"{stem}_barcodes.tsv.gz"; mat=RAW/m["name"]
        feats=pd.read_csv(feat,sep="\t",header=None)
        symbols=feats.iloc[:,1].astype(str).tolist() if feats.shape[1]>1 else feats.iloc[:,0].astype(str).tolist()
        rows_by_gene={g:i for i,g in enumerate(symbols) if g in GENES}
        if not rows_by_gene: continue
        bars=[x.strip() for x in gzip.open(bar,"rt",errors="replace")]
        X=scipy.io.mmread(gzip.open(mat,"rb")).tocsr()
        for g,row in rows_by_gene.items():
            idx=X.getrow(row).indices
            mapped=[]
            # GEO member names encode the pool as pool_A..pool_I or pool_dexA..pool_dexI;
            # cell metadata uses HTMC_pool-A... / HTMC_pool-dexA... identifiers.
            m = re.search(r"_pool_(dex)?([A-I])$", stem)
            if not m: continue
            library=f"HTMC_pool-{'dex' if m.group(1) else ''}{m.group(2)}_GRCh38"
            for j in idx:
                key=(library, bars[j])
                if key in by_barcode: mapped.append(by_barcode[key])
            all_results[g]["cells_detected"] += len(mapped)
            all_results[g]["donors"].update(x[0] for x in mapped)
            all_results[g]["vehicle"] += sum(x[1]=="no-dex" for x in mapped)
            all_results[g]["dex"] += sum(x[1]=="dex" for x in mapped)
    # Percentages use only cell metadata rows whose (library, barcode) is
    # present in the corresponding expression matrix.
    matrix_keys=set()
    for m in matrix_members:
        stem=m["name"].replace("_matrix.mtx.gz",""); mm=re.search(r"_pool_(dex)?([A-I])$",stem)
        if not mm: continue
        library=f"HTMC_pool-{'dex' if mm.group(1) else ''}{mm.group(2)}_GRCh38"
        for b in gzip.open(RAW/f"{stem}_barcodes.tsv.gz","rt",errors="replace"):
            matrix_keys.add((library,b.strip()))
    analysis_keys=set(by_barcode) & matrix_keys
    ncell=len(analysis_keys); ndonor=cm.donor_id_demuxlet.nunique(); vehicle_n=sum(by_barcode[k][1]=="no-dex" for k in analysis_keys); dex_n=sum(by_barcode[k][1]=="dex" for k in analysis_keys)
    rows=[]
    for g in GENES:
        x=all_results[g]; rows.append({"gene":g,"status":"OK","cells_detected":x["cells_detected"],"pct_cells_detected":round(100*x["cells_detected"]/ncell,6),"donors_detected":len(x["donors"]),"pct_donors_detected":round(100*len(x["donors"])/ndonor,6),"vehicle_pct":round(100*x["vehicle"]/vehicle_n,6),"dex_pct":round(100*x["dex"]/dex_n,6),"percentage_denominator":"metadata-matrix barcode intersection"})
    pd.DataFrame(rows).to_csv(OUT/"GENE_DETECTABILITY.tsv",sep="\t",index=False)

if __name__=="__main__": main()
