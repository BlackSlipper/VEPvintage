#!/usr/bin/env python
"""Core analysis: AUROC of each predictor stratified by label vintage; vintage slopes; paired differences."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
import statsmodels.api as sm
D=load()
print("evaluation set:",len(D),D.final.value_counts().to_dict())
tools=[t for t in TOOLS if t in D.columns]
cov=D[tools].notna().mean().sort_values()
print("coverage:\n",cov.round(3).to_string())
D["vbin"]=D.vyear.clip(upper=2026).astype(str)
D.loc[D.legacy,"vbin"]="<=2015Q1"
order=["<=2015Q1"]+[str(y) for y in range(2015,2027)]
# --- AUROC per tool x vintage bin (all available variants)
rows=[]
for t in tools:
    for vb in order:
        g=D[D.vbin==vb]
        if len(g)==0: continue
        y=g.y.values; s=g[t].values.astype(float)
        a=auc(y,s)
        lo,hi=boot_auc(y,s,B=100) if not np.isnan(a) else (np.nan,np.nan)
        m=~np.isnan(s)
        rows.append(dict(tool=t,vbin=vb,n=int(m.sum()),n_plp=int(y[m].sum()),auc=a,lo=lo,hi=hi))
A=pd.DataFrame(rows); A.to_csv(f"{P}/results/auc_by_vintage_all.tsv",sep="\t",index=False)
# --- common-coverage subset: tools with >=90% coverage
core=[t for t in tools if cov[t]>=0.90]
C=D.dropna(subset=core)
print("core tools:",len(core),"common subset n=",len(C))
rows=[]
for t in core:
    for vb in order:
        g=C[C.vbin==vb]
        a=auc(g.y.values,g[t].values.astype(float))
        rows.append(dict(tool=t,vbin=vb,n=len(g),n_plp=int(g.y.sum()),auc=a))
AC=pd.DataFrame(rows); AC.to_csv(f"{P}/results/auc_by_vintage_core.tsv",sep="\t",index=False)
W=AC.pivot(index="tool",columns="vbin",values="auc")[order]
print(W.round(3).to_string())
# --- vintage slope: AUROC ~ year over 2016..2026 bins (excluding legacy), weighted by n
yrs=[str(y) for y in range(2016,2027)]
sl=[]
for t in core:
    s=AC[(AC.tool==t)&AC.vbin.isin(yrs)].dropna()
    X=sm.add_constant(s.vbin.astype(int)-2016); w=s.n
    f=sm.WLS(s.auc,X,weights=w).fit()
    sl.append(dict(tool=t,slope_per_year=f.params.iloc[1],se=f.bse.iloc[1],auc_2016_19=np.average(s[s.vbin.astype(int)<=2019].auc,weights=s[s.vbin.astype(int)<=2019].n),auc_2023_26=np.average(s[s.vbin.astype(int)>=2023].auc,weights=s[s.vbin.astype(int)>=2023].n)))
SL=pd.DataFrame(sl).merge(META[["tool","category","clinvar_in_training","clingen_calibration","pub_year"]],on="tool")
SL["delta"]=SL.auc_2023_26-SL.auc_2016_19
SL=SL.sort_values("slope_per_year",ascending=False); SL.to_csv(f"{P}/results/vintage_slopes_core.tsv",sep="\t",index=False)
print(SL.round(4).to_string())
