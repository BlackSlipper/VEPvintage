#!/usr/bin/env python
"""Drift-robust analyses: AF-stratified AUROC, within-vintage ranks, relative logit-AUROC slopes, category contrasts."""
import sys, os; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
import statsmodels.formula.api as smf
D=load()
AF=pd.concat([pd.read_csv(f,sep="\t") for f in sorted(glob.glob(f"{P}/results/af_shard*.tsv.gz"))]).drop_duplicates("query")
for c in AF.columns[1:]: AF[c]=pd.to_numeric(AF[c],errors="coerce")
AF["maxaf"]=AF[["gnomad_ex_popmax","gnomad_ge_popmax","gnomad_ex_af","gnomad_ge_af","exac_af"]].max(axis=1).fillna(0)
D=D.merge(AF[["query","maxaf"]],on="query",how="left"); D["maxaf"]=D.maxaf.fillna(0)
D["afbin"]=pd.cut(D.maxaf,[-1,0,1e-4,1e-3,1e-2,1.1],labels=["absent","<1e-4","1e-4-1e-3","1e-3-1e-2",">=1e-2"])
tools=[t for t in TOOLS if t in D.columns]
cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]
D["vbin"]=np.where(D.legacy,"<=2015Q1",D.vyear.astype(str))
order=["<=2015Q1"]+[str(y) for y in range(2015,2027)]
# composition by vintage
comp=D.groupby("vbin").apply(lambda g: pd.Series({"n":len(g),"frac_PLP":g.y.mean(),"n_genes":g.Gene.nunique(),
    "BLB_af_ge1e-3":(g[g.y==0].maxaf>=1e-3).mean(),"BLB_af_absent":(g[g.y==0].maxaf==0).mean(),"PLP_af_absent":(g[g.y==1].maxaf==0).mean(),
    "stars>=2":(g.stars>=2).mean(),"prior_VUS":g.prior_VUS.mean()})).reindex(order)
comp.to_csv(f"{P}/results/composition_by_vintage.tsv",sep="\t"); print(comp.round(3).to_string())
def auc_table(G,name):
    rows=[]
    for t in core:
        for vb in order:
            g=G[G.vbin==vb]; rows.append(dict(subset=name,tool=t,vbin=vb,n=len(g),n_plp=int(g.y.sum()),auc=auc(g.y.values,g[t].values.astype(float))))
    return pd.DataFrame(rows)
C=D.dropna(subset=core)
R=pd.concat([auc_table(C,"all_common"),auc_table(C[C.maxaf<1e-4],"rare_common"),auc_table(C[C.maxaf==0],"absent_common")])
R.to_csv(f"{P}/results/auc_by_vintage_subsets.tsv",sep="\t",index=False)
for sub in ["rare_common","absent_common"]:
    W=R[R.subset==sub].pivot(index="tool",columns="vbin",values="auc")[order]; print(sub); print(W.round(3).to_string())
# relative logit AUROC: tool minus vintage median across core tools
def logit(p): p=np.clip(p,1e-4,1-1e-4); return np.log(p/(1-p))
R["lauc"]=logit(R.auc)
R["rel"]=R.lauc-R.groupby(["subset","vbin"]).lauc.transform("median")
R["rank"]=R.groupby(["subset","vbin"]).auc.rank(ascending=False)
R.to_csv(f"{P}/results/auc_by_vintage_subsets.tsv",sep="\t",index=False)
yrs=[str(y) for y in range(2016,2027)]
out=[]
for sub in R.subset.unique():
    for t in core:
        s=R[(R.subset==sub)&(R.tool==t)&R.vbin.isin(yrs)].dropna(subset=["rel"])
        if len(s)<5: continue
        f=smf.wls("rel~I(vbin.astype(int)-2016)",data=s,weights=s.n).fit()
        out.append(dict(subset=sub,tool=t,rel_slope=f.params.iloc[1],se=f.bse.iloc[1],p=f.pvalues.iloc[1],
            rel_early=np.average(s[s.vbin.astype(int)<=2019].rel,weights=s[s.vbin.astype(int)<=2019].n) if (s.vbin.astype(int)<=2019).any() else np.nan,
            rel_late=np.average(s[s.vbin.astype(int)>=2023].rel,weights=s[s.vbin.astype(int)>=2023].n)))
O=pd.DataFrame(out).merge(META[["tool","category","clinvar_in_training","clingen_calibration","pub_year"]],on="tool")
O.to_csv(f"{P}/results/relative_slopes.tsv",sep="\t",index=False)
pd.set_option("display.width",250)
print(O[O.subset=="rare_common"].sort_values("rel_slope",ascending=False).round(4).to_string())
print(O.groupby(["subset","category"]).rel_slope.mean().round(4).unstack(0).to_string())
