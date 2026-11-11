#!/usr/bin/env python
"""Selection null on the same footing as Fig. 3c,d: 38 core predictors, variants scored by all of them (common subset)."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
D=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",dtype={"VariationID":str})
E=load(); tools=[t for t in TOOLS if t in E.columns]; cov=E[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]; del E
D=D[D.rare].dropna(subset=core)
Pc=D[(D.y==1)&D.code_citing]; Nc=D[(D.y==0)&D.code_citing]; Pall=D[D.y==1]; Nall=D[D.y==0]
rows=[]
for side,grp,other,obs_flag,sim_flag in [("PP3",Pc,Nall,Pc.has_PP3.fillna(False),Pc.REVEL_raw>=0.644),("BP4",Nc,Pall,Nc.has_BP4.fillna(False),Nc.REVEL_raw<=0.290)]:
    for kind,flag in [("observed",obs_flag),("simulated_REVEL",sim_flag)]:
        for t in core:
            a1=pd.concat([grp[flag.values],other]); a0=pd.concat([grp[~flag.values],other])
            rows.append(dict(side=side,kind=kind,tool=t,gap=auc(a1.y.values,a1[t].values.astype(float))-auc(a0.y.values,a0[t].values.astype(float)),frac_flag=float(flag.mean()),n_grp=len(grp)))
R=pd.DataFrame(rows); R.to_csv(f"{P}/results/selection_null_common.tsv",sep="\t",index=False)
W=R.pivot_table(index=["side","tool"],columns="kind",values="gap")
print(len(core), R.groupby(["side","kind"])[["frac_flag","n_grp"]].first())
for s in ["PP3","BP4"]:
    w=W.loc[s]; print(s,"median observed %.3f simulated %.3f r=%.2f; sim>obs for %d of %d"%(w.observed.median(),w.simulated_REVEL.median(),w.corr().iloc[0,1],(w.simulated_REVEL>w.observed).sum(),len(w)))
