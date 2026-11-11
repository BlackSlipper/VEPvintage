#!/usr/bin/env python
"""Likelihood ratios of ClinGen PP3/BP4 score intervals by label-vintage period."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
D=load()
import os
if os.environ.get("RARE"):
    A=pd.read_csv(f"{P}/results/variant_level_annotated.tsv.gz",sep="\t",usecols=["query","maxaf"]).drop_duplicates("query")
    D=D.merge(A,on="query",how="left"); D=D[D.maxaf.fillna(0)<1e-4]
    OUT="calibration_LR_by_period_rare.tsv"
    if os.environ.get("CLEAN"):
        Tt=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","clean"])
        D=D.merge(Tt,on="hgvs38",how="left"); D=D[D.clean.fillna(False)]
        OUT="calibration_LR_by_period_rare_clean.tsv"
else: OUT="calibration_LR_by_period.tsv"
T=pd.read_csv(f"{P}/meta/clingen_thresholds.tsv",sep="\t")
T["lo"]=T.lo.astype(float); T["hi"]=T.hi.astype(float)
def period(y):
    if y<=2019: return "<=2019"
    if y<=2022: return "2020-2022"
    return "2023-2026"
D["period"]=D.vyear.map(period)
D.loc[D.legacy,"period"]="<=2019"
ODDS={"Supporting":2.08,"Moderate":4.33,"Strong":18.7,"VeryStrong":350}
rows=[]
rng=np.random.default_rng(1)
for (tool,rc),tt in T.groupby(["tool","rawcol"]):
    if rc not in D.columns: continue
    for per,g in list(D.groupby("period"))+[("all",D)]:
        g=g[g[rc].notna()]
        npl=(g.y==1).sum(); nbl=(g.y==0).sum()
        for _,r in tt.iterrows():
            inb=(g[rc]>=r.lo)&(g[rc]<r.hi)
            a=(inb&(g.y==1)).sum(); b=(inb&(g.y==0)).sum()
            # LR for PP3 = P(bin|P)/P(bin|B); for BP4 report inverse (benign LR) so larger = stronger
            lr=((a+0.5)/(npl+0.5))/((b+0.5)/(nbl+0.5))
            if r.evidence.startswith("BP4"): lr=1/lr
            # bootstrap CI
            bs=[]
            ya=g.y.values; ib=inb.values
            for _ in range(200):
                i=rng.integers(0,len(ya),len(ya)); yy=ya[i]; bb=ib[i]
                A=(bb&(yy==1)).sum(); Bn=(bb&(yy==0)).sum(); n1=yy.sum(); n0=len(yy)-n1
                l=((A+0.5)/(n1+0.5))/((Bn+0.5)/(n0+0.5)); bs.append(1/l if r.evidence.startswith("BP4") else l)
            lo,hi=np.percentile(bs,[2.5,97.5])
            rows.append(dict(tool=tool,evidence=r.evidence,period=per,n_plp=int(npl),n_blb=int(nbl),in_bin_plp=int(a),in_bin_blb=int(b),LR=lr,lo=lo,hi=hi,nominal=ODDS[r.evidence.split("_")[1]],source=r.source))
R=pd.DataFrame(rows); R.to_csv(f"{P}/results/"+OUT,sep="\t",index=False)
pd.set_option("display.width",250)
print(R.pivot_table(index=["tool","evidence"],columns="period",values="LR").round(2).to_string())
