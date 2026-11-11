#!/usr/bin/env python
"""Provenance-aware benchmark tiers; rankings full vs clean; gene-cluster bootstrap CIs."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
from scipy.stats import kendalltau, spearmanr
D=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",dtype={"VariationID":str})
S=pd.read_csv(f"{P}/results/scv_textmining.tsv.gz",sep="\t",dtype={"VariationID":str})
cs=S.ClinSig.str.lower()
S["scv_class"]=np.where(cs.str.contains("pathogenic"),"PLP",np.where(cs.str.contains("benign"),"BLB","OTHER"))
S=S.merge(D[["VariationID","final"]],on="VariationID"); S=S[S.scv_class==S.final]
S["ambry"]=S.Submitter.eq("Ambry Genetics")
S["comp_nonambry"]=(~S.ambry)&(S.any_insilico==1)
S["cons_nonambry"]=(~S.ambry)&(S.conservation_text==1)
V=S.groupby("VariationID").agg(comp_nonambry=("comp_nonambry","max"),cons_nonambry=("cons_nonambry","max"),n_nonambry=("ambry",lambda s:(~s).sum()))
D=D.merge(V,left_on="VariationID",right_index=True,how="left")
D["comp_nonambry"]=D.comp_nonambry.fillna(False).astype(bool); D["cons_nonambry"]=D.cons_nonambry.fillna(False).astype(bool); D["n_nonambry"]=D.n_nonambry.fillna(0)
D["explicit_comp"]=np.where(D.y==1,D.has_PP3.fillna(False),D.has_BP4.fillna(False)).astype(bool)
D["ambry_only"]=D.n_nonambry==0
D["clean"]=(~D.explicit_comp)&(~D.comp_nonambry)&(~D.ambry_only)
D["clean_strict"]=D.clean&((D.stars>=2)|(D.n_conc_scv.fillna(0)>=2))
D["clean_cons"]=D.clean&~D.cons_nonambry
print(pd.crosstab(D.period,[D.final,D.clean]))
print("clean_cons n:",D.clean_cons.sum(),D[D.clean_cons].final.value_counts().to_dict())
print("clean n:",D.clean.sum(),D[D.clean].final.value_counts().to_dict(),"strict:",D.clean_strict.sum(),D[D.clean_strict].final.value_counts().to_dict())
D[["AlleleID","VariationID","hgvs38","Gene","Protein","final","vintage_first","period","stars","maxaf","rare","explicit_comp","comp_nonambry","cons_nonambry","ambry_only","clean","clean_strict","clean_cons","prior_VUS"]].to_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",index=False,compression="gzip")
tools=[t for t in TOOLS if t in D.columns]
def gene_boot(G,t,B=200,seed=0):
    rng=np.random.default_rng(seed)
    y=G.y.values.astype(float); sc=G[t].values.astype(float)
    lev=pd.Series(sc).rank(method="dense").values.astype(int)-1; K=lev.max()+1
    gcode,gun=pd.factorize(G.Gene.values); ng=len(gun)
    out=[]
    for _ in range(B):
        wg=np.bincount(rng.integers(0,ng,ng),minlength=ng).astype(float); w=wg[gcode]
        Pk=np.bincount(lev,weights=w*y,minlength=K); Nk=np.bincount(lev,weights=w*(1-y),minlength=K)
        cN=np.cumsum(Nk)-Nk
        tp,tn=Pk.sum(),Nk.sum()
        out.append(np.nan if tp==0 or tn==0 else (Pk*(cN+0.5*Nk)).sum()/(tp*tn))
    return np.nanpercentile(out,[2.5,97.5])
rows=[]
sets={"full":D,"full_rare":D[D.rare],"clean":D[D.clean],"clean_rare":D[D.clean&D.rare],"clean_strict":D[D.clean_strict],
      "nonclean":D[~D.clean],"clean_cons":D[D.clean_cons],"clean_cons_rare":D[D.clean_cons&D.rare],"clean_2023_26":D[D.clean&(D.period=="2023-2026")],"full_2023_26":D[D.period=="2023-2026"]}
for name,G in sets.items():
    for t in tools:
        g=G[G[t].notna()]
        a=auc(g.y.values,g[t].values.astype(float))
        lo,hi=gene_boot(g,t,B=200) if name in ("full","clean","clean_rare","full_rare") else (np.nan,np.nan)
        rows.append(dict(set=name,tool=t,n=len(g),n_plp=int(g.y.sum()),coverage=len(g)/len(G),auc=a,lo=lo,hi=hi))
R=pd.DataFrame(rows); R["rank"]=R.groupby("set").auc.rank(ascending=False)
R.to_csv(f"{P}/results/benchmark_full_vs_clean.tsv",sep="\t",index=False)
W=R.pivot(index="tool",columns="set",values="auc"); RK=R.pivot(index="tool",columns="set",values="rank")
W=W.join(META.set_index("tool")[["category"]])
W["delta_full_minus_clean"]=W.full-W.clean
W["rank_full"]=RK.full; W["rank_clean"]=RK.clean; W["rank_shift"]=RK.full-RK.clean
pd.set_option("display.width",250)
print(W.sort_values("clean",ascending=False)[["category","full","clean","full_rare","clean_rare","clean_strict","delta_full_minus_clean","rank_full","rank_clean","rank_shift"]].round(3).to_string())
for a,b in [("full","clean"),("full_rare","clean_rare"),("full_2023_26","clean_2023_26"),("full","clean_strict")]:
    x=W[[a,b]].dropna(); print(a,"vs",b,"kendall tau=%.3f spearman=%.3f"%(kendalltau(x[a],x[b])[0],spearmanr(x[a],x[b])[0]))
print(W.groupby("category").delta_full_minus_clean.mean().round(4))
