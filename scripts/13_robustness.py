#!/usr/bin/env python
"""Sensitivity analyses: vintage definition, review stars, Submitter A exclusion, VARITY in-sample vs LOO, PHACTboost training split."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
D=load(); T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","clean","rare","explicit_comp","comp_nonambry","ambry_only"])
D=D.merge(T,on="hgvs38")
D["per_first"]=np.where(D.legacy|(D.vyear<=2019),"<=2019",np.where(D.vyear<=2022,"2020-2022","2023-2026"))
D["per_run"]=np.where(D.vyear_run<=2019,"<=2019",np.where(D.vyear_run<=2022,"2020-2022","2023-2026"))
key=["REVEL","BayesDel_noAF","AlphaMissense","ESM1b","CADD","phyloP100way","GERP++","ClinPred","MetaRNN","VARITY_R","VARITY_R_LOO","PHACTboost","EVE","MutPred2"]
rows=[]
def add(name,G,per_col):
    for p in ["<=2019","2020-2022","2023-2026"]:
        g=G[G[per_col]==p]
        for t in key:
            m=g[t].notna(); rows.append(dict(analysis=name,period=p,tool=t,n=int(m.sum()),n_plp=int(g.y[m].sum()),auc=auc(g.y.values,g[t].values.astype(float))))
R=D[D.rare]
add("rare_vintage_first",R,"per_first")
add("rare_vintage_run",R,"per_run")
add("rare_stars>=2",R[R.stars>=2],"per_first")
add("rare_excl_SubmitterA_only",R[~R.ambry_only],"per_first")
add("rare_clean",R[R.clean],"per_first")
add("rare_clean_vintage_run",R[R.clean],"per_run")
O=pd.DataFrame(rows); O.to_csv(f"{P}/results/robustness.tsv",sep="\t",index=False)
pd.set_option("display.width",250)
print(O.pivot_table(index=["analysis","tool"],columns="period",values="auc").round(3).loc[(slice(None),["REVEL","ESM1b","phyloP100way","GERP++","AlphaMissense"]),:].to_string())
print(O.groupby(["analysis","period"]).n.max().unstack().to_string())
# VARITY in-sample vs LOO by vintage year (all variants)
D["vb"]=np.where(D.legacy,2015,D.vyear)
v=[]
for y,g in D.groupby("vb"):
    v.append(dict(vyear=y,n=len(g),VARITY_R=auc(g.y.values,g.VARITY_R.values.astype(float)),VARITY_R_LOO=auc(g.y.values,g.VARITY_R_LOO.values.astype(float)),PHACTboost=auc(g.y.values,g.PHACTboost.values.astype(float)),ClinPred=auc(g.y.values,g.ClinPred.values.astype(float)),ESM1b=auc(g.y.values,g.ESM1b.values.astype(float))))
V=pd.DataFrame(v); V["VARITY_insample_gain"]=V.VARITY_R-V.VARITY_R_LOO; print(V.round(4).to_string())
V.to_csv(f"{P}/results/leakage_by_vintage.tsv",sep="\t",index=False)
