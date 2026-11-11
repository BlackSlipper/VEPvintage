#!/usr/bin/env python
"""Labels as they stood in each January release 2016-2026 and in November 2026 (P/LP or B/LB, >=1 star, germline missense
SNVs), with the November 2026 class of each allele -> results/download_date_labels.tsv.gz. Also lists the hg38 ids of
alleles not in the evaluation history (label since changed) for fetching scores and allele frequencies
(results/extra_alleles_hg38.tsv, results/ids_extra.txt; then 03_fetch_dbnsfp.py -> results/dbnsfp_extra.tsv.gz and
03c_fetch_af.py -> results/af_extra.tsv.gz). Needs data/snap_parsed (01) and the current variant_summary (00)."""
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np
P=os.environ.get("VEP_ROOT",os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
V=pd.read_csv(f"{P}/results/variant_history.tsv.gz",sep="\t",usecols=["AlleleID"]); known=set(V.AlleleID)
STAR=lambda r: r.fillna("").str.contains("criteria provided|expert panel|practice guideline")&~r.fillna("").str.contains("no assertion criteria|conflicting")
snaps=[f"{y}-01" for y in range(2016,2027)]+["current"]; L=[]
for s in snaps:
    d=pd.read_csv(f"{P}/data/snap_parsed/{s}.tsv.gz",sep="\t",dtype=str,usecols=["AlleleID","Class","Review","Origin","Gene"])
    d=d[d.Class.isin(["PLP","BLB"])&d.Origin.fillna("").str.contains("germline|de novo|inherited|maternal|paternal|biparental|uniparental")&STAR(d.Review)].copy()
    d["AlleleID"]=d.AlleleID.astype(int); d=d.drop_duplicates("AlleleID"); d["release"]="2026-11" if s=="current" else s
    L.append(d[["release","AlleleID","Class","Gene"]])
Lb=pd.concat(L)
cur=pd.read_csv(f"{P}/data/snap_parsed/current.tsv.gz",sep="\t",dtype=str,usecols=["AlleleID","Class"]); cur["AlleleID"]=cur.AlleleID.astype(int)
Lb["class_2026_11"]=Lb.AlleleID.map(cur.drop_duplicates("AlleleID").set_index("AlleleID").Class).fillna("removed")
Lb.to_csv(f"{P}/results/download_date_labels.tsv.gz",sep="\t",index=False,compression="gzip")
need=set(Lb.AlleleID)-known; print("rows",len(Lb),"alleles needed:",len(need))
cols=["#AlleleID","Assembly","Chromosome","PositionVCF","ReferenceAlleleVCF","AlternateAlleleVCF"]; out=[]
for ch in pd.read_csv(f"{P}/data/clinvar_tab/variant_summary_current.txt.gz",sep="\t",dtype=str,usecols=cols,chunksize=500000):
    out.append(ch[(ch.Assembly=="GRCh38")&ch["#AlleleID"].astype(int).isin(need)])
X=pd.concat(out); X=X[X.ReferenceAlleleVCF.str.len().eq(1)&X.AlternateAlleleVCF.str.len().eq(1)&(X.ReferenceAlleleVCF!="na")]
X["hgvs38"]="chr"+X.Chromosome+":g."+X.PositionVCF+X.ReferenceAlleleVCF+">"+X.AlternateAlleleVCF
X=X.rename(columns={"#AlleleID":"AlleleID"})[["AlleleID","hgvs38"]].drop_duplicates("AlleleID")
print("with GRCh38 SNV coordinates:",len(X),"of",len(need))
X.to_csv(f"{P}/results/extra_alleles_hg38.tsv",sep="\t",index=False)
X.hgvs38.drop_duplicates().to_csv(f"{P}/results/ids_extra.txt",index=False,header=False)
