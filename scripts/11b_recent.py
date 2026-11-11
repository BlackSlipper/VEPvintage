#!/usr/bin/env python
"""Post-training benchmark: rare variants labelled 2023-2026, all vs provenance-filtered (and conservation-filtered) labels."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
from scipy.stats import kendalltau
D=load(); T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","clean","clean_cons","rare","period"])
D=D.merge(T,on="hgvs38")
tools=[t for t in TOOLS if t in D.columns]
def gene_boot(G,t,B=200,seed=0):
    rng=np.random.default_rng(seed)
    y=G.y.values.astype(float); sc=G[t].values.astype(float)
    lev=pd.Series(sc).rank(method="dense").values.astype(int)-1; K=lev.max()+1
    gcode,gun=pd.factorize(G.Gene.values); ng=len(gun); out=[]
    for _ in range(B):
        wg=np.bincount(rng.integers(0,ng,ng),minlength=ng).astype(float); w=wg[gcode]
        Pk=np.bincount(lev,weights=w*y,minlength=K); Nk=np.bincount(lev,weights=w*(1-y),minlength=K)
        cN=np.cumsum(Nk)-Nk; tp,tn=Pk.sum(),Nk.sum()
        out.append(np.nan if tp==0 or tn==0 else (Pk*(cN+0.5*Nk)).sum()/(tp*tn))
    return np.nanpercentile(out,[2.5,97.5])
rows=[]
sets={"full_rare_2023_26":D[D.rare&(D.period=="2023-2026")],"clean_rare_2023_26":D[D.clean&D.rare&(D.period=="2023-2026")],
      "cleancons_rare_2023_26":D[D.clean_cons&D.rare&(D.period=="2023-2026")],"clean_rare_postVUS":D[D.clean&D.rare&D.prior_VUS]}
for name,G in sets.items():
    for t in tools:
        g=G[G[t].notna()]; lo,hi=gene_boot(g,t) if name in ("full_rare_2023_26","clean_rare_2023_26") else (np.nan,np.nan)
        rows.append(dict(set=name,tool=t,n=len(g),n_plp=int(g.y.sum()),cov=len(g)/max(len(G),1),auc=auc(g.y.values,g[t].values.astype(float)),lo=lo,hi=hi))
R=pd.DataFrame(rows); R.to_csv(f"{P}/results/benchmark_clean_rare_recent.tsv",sep="\t",index=False)
W=R.pivot(index="tool",columns="set",values="auc"); N=R.pivot(index="tool",columns="set",values="cov")
W["rank_full"]=W.full_rare_2023_26.rank(ascending=False); W["rank_clean"]=W.clean_rare_2023_26.rank(ascending=False); W["cov_clean"]=N.clean_rare_2023_26
W["drop"]=W.full_rare_2023_26-W.clean_rare_2023_26
pd.set_option("display.width",220)
print(W.sort_values("clean_rare_2023_26",ascending=False).round(3).to_string())
print(R.groupby("set")[["n","n_plp"]].max().to_string())
print("drop range",round(W["drop"].min(),3),round(W["drop"].max(),3))
print("tau full vs clean (rare 2023-26):",round(kendalltau(W.full_rare_2023_26,W.clean_rare_2023_26)[0],3))
