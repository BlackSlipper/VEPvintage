#!/usr/bin/env python
"""Merge parsed snapshots -> per-variant label history, vintage, and coordinates (current GRCh38)."""
import glob, gzip, os, re
import pandas as pd, numpy as np
P=__import__("os").environ.get("VEP_ROOT",__import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))
fs=sorted(glob.glob(f"{P}/data/snap_parsed/*.tsv.gz"))
snaps=sorted([os.path.basename(f)[:-7] for f in fs if "current" not in f])+["current"]
order={s:i for i,s in enumerate(snaps)}
CUR_DATE="2026-11"
def sdate(s): return CUR_DATE if s=="current" else s
cur0=pd.read_csv(f"{P}/data/snap_parsed/current.tsv.gz",sep="\t",dtype=str,usecols=["AlleleID","Class"])
keep=set(cur0.loc[cur0.Class.isin(["PLP","BLB","CONF"]),"AlleleID"].astype(int))
print("alleles kept (final PLP/BLB/CONF):",len(keep),flush=True)
L=[]
for f in fs:
    d=pd.read_csv(f,sep="\t",dtype=str,usecols=["Snapshot","AlleleID","Class","Review","Gene","Protein"])
    d["AlleleID"]=d.AlleleID.astype(int)
    L.append(d[d.AlleleID.isin(keep)])
H=pd.concat(L,ignore_index=True)
H["ord"]=H.Snapshot.map(order)
H=H.sort_values(["AlleleID","ord"])
print("snapshots:",len(snaps),"rows:",len(H),"unique alleles:",H.AlleleID.nunique())
cur=H[H.Snapshot=="current"].set_index("AlleleID")
# history string
hist=H.groupby("AlleleID").apply(lambda g: ";".join(f"{s}:{c}" for s,c in zip(g.Snapshot,g.Class)))
first_seen=H.groupby("AlleleID").Snapshot.first()
def vint(g):
    fin=g.Class.iloc[-1]
    if g.Snapshot.iloc[-1]!="current": return pd.Series({"final":np.nan})
    cls=g.Class.tolist(); sn=g.Snapshot.tolist()
    first=next(s for s,c in zip(sn,cls) if c==fin)
    # start of final uninterrupted run
    i=len(cls)-1
    while i>0 and cls[i-1]==fin: i-=1
    run=sn[i]
    prior=set(cls[:cls.index(fin)])
    return pd.Series({"final":fin,"vintage_first":sdate(first),"vintage_run":sdate(run),
        "prior_VUS":"VUS" in prior or "CONF" in prior,
        "prior_opposite":("BLB" in prior) if fin=="PLP" else (("PLP" in prior) if fin=="BLB" else False),
        "n_changes":int(sum(1 for a,b in zip(cls,cls[1:]) if a!=b)),
        "first_seen":sdate(sn[0])})
V=H.groupby("AlleleID").apply(vint)
V=V[V.final.notna()].copy()
V=V.join(cur[["Gene","Protein","Review"]].rename(columns={"Review":"Review_current"}))
# coordinates and origin from current variant_summary
cols=["AlleleID","Assembly","Chromosome","PositionVCF","ReferenceAlleleVCF","AlternateAlleleVCF","VariationID","OriginSimple","ReviewStatus","NumberSubmitters","SubmitterCategories","PhenotypeList","Type"]
it=pd.read_csv(f"{P}/data/clinvar_tab/variant_summary_current.txt.gz",sep="\t",dtype=str,usecols=lambda c: c.lstrip("#") in cols,chunksize=500000)
C=[]
for ch in it:
    ch.columns=[c.lstrip("#") for c in ch.columns]
    ch=ch[(ch.Assembly=="GRCh38")&(ch.Type=="single nucleotide variant")]
    C.append(ch)
C=pd.concat(C); C["AlleleID"]=C.AlleleID.astype(int); C=C.drop_duplicates("AlleleID").set_index("AlleleID")
V=V.join(C[["Chromosome","PositionVCF","ReferenceAlleleVCF","AlternateAlleleVCF","VariationID","OriginSimple","NumberSubmitters","SubmitterCategories","PhenotypeList"]],how="left")
ok=V.PositionVCF.notna()&V.ReferenceAlleleVCF.str.fullmatch("[ACGT]")&V.AlternateAlleleVCF.str.fullmatch("[ACGT]")
V["hgvs38"]=np.where(ok,"chr"+V.Chromosome.astype(str)+":g."+V.PositionVCF.astype(str)+V.ReferenceAlleleVCF.astype(str)+">"+V.AlternateAlleleVCF.astype(str),"")
stars={"practice guideline":4,"reviewed by expert panel":3,"criteria provided, multiple submitters, no conflicts":2,"criteria provided, single submitter":1}
V["stars"]=V.Review_current.map(stars).fillna(0).astype(int)
V["history"]=hist.reindex(V.index)
V.to_csv(f"{P}/results/variant_history.tsv.gz",sep="\t",compression="gzip")
print(V.final.value_counts().to_dict())
print("with coords:",(V.hgvs38!="").sum())
E=V[V.final.isin(["PLP","BLB"])&(V.stars>=1)&V.OriginSimple.fillna("").str.contains("germline")&(V.hgvs38!="")]
print("evaluation set:",len(E),E.final.value_counts().to_dict())
print(E.groupby(E.vintage_first.str[:4]).final.value_counts().unstack().fillna(0).astype(int))
ids=V.loc[(V.hgvs38!="")&V.final.isin(["PLP","BLB"]),"hgvs38"].drop_duplicates()
ids.to_csv(f"{P}/results/ids_hg38_plpblb.txt",index=False,header=False)
ids2=V.loc[(V.hgvs38!="")&V.final.isin(["CONF"]),"hgvs38"].drop_duplicates()
ids2.to_csv(f"{P}/results/ids_hg38_conf.txt",index=False,header=False)
print("ids for annotation:",len(ids),"conf:",len(ids2))
