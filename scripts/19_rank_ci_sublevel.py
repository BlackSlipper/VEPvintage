#!/usr/bin/env python
"""(iv) gene-cluster bootstrap rank intervals (500 reps) for rare 2023-26: full, per-variant filtered, submitter-level;
consistency of rank shifts between filters; (vi) submitter-level calibration LRs by period."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
from scipy.stats import spearmanr
D=load()
T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","clean","rare","period"])
D=D.merge(T,on="hgvs38")
S=pd.read_csv(f"{P}/results/scv_textmining.tsv.gz",sep="\t",dtype={"VariationID":str})
C=pd.read_csv(f"{P}/results/scv_acmg_codes.tsv.gz",sep="\t",dtype={"VariationID":str})
st=S[S.has_text==1].groupby("Submitter").agg(n_text=("SCV","size"),frac_comp=("any_insilico","mean"))
cc=C.assign(comp=C.codes.str.contains(r"\bPP3\b|\bBP4\b")).groupby("Submitter").comp.mean().rename("fcc")
st=st.join(cc).fillna({"fcc":0}); noncomp=set(st[(st.n_text>=20)&(st.frac_comp<0.05)&(st.fcc<0.05)].index)-{"Ambry Genetics"}
cs=S.ClinSig.str.lower(); S["cl"]=np.where(cs.str.contains("pathogenic"),"PLP",np.where(cs.str.contains("benign"),"BLB","O"))
S=S.merge(D[["VariationID","final"]].drop_duplicates(),on="VariationID"); S=S[S.cl==S.final]
D["subclean"]=D.VariationID.isin(set(S[S.Submitter.isin(noncomp)].VariationID))
print("Submitter B qualifies:", "Labcorp Genetics (formerly Invitae), Labcorp" in noncomp)
tools=[t for t in TOOLS if t in D.columns]
def boot_ranks(G,B=500,seed=1):
    rng=np.random.default_rng(seed); gcode,gun=pd.factorize(G.Gene.values); ng=len(gun); y=G.y.values.astype(float)
    prep={}
    for t in tools:
        s=G[t].values.astype(float); m=~np.isnan(s)
        lev=pd.Series(s[m]).rank(method="dense").values.astype(int)-1
        prep[t]=(m,lev,lev.max()+1)
    R=np.zeros((B,len(tools)))
    for b in range(B):
        wg=np.bincount(rng.integers(0,ng,ng),minlength=ng).astype(float); w=wg[gcode]
        a=[]
        for t in tools:
            m,lev,K=prep[t]; ww=w[m]; yy=y[m]
            Pk=np.bincount(lev,weights=ww*yy,minlength=K); Nk=np.bincount(lev,weights=ww*(1-yy),minlength=K)
            cN=np.cumsum(Nk)-Nk; a.append((Pk*(cN+0.5*Nk)).sum()/(Pk.sum()*Nk.sum()))
        R[b]=pd.Series(a).rank(ascending=False).values
    return pd.DataFrame({"rank_lo":np.percentile(R,2.5,axis=0),"rank_med":np.median(R,axis=0),"rank_hi":np.percentile(R,97.5,axis=0)},index=tools)
sets={"full":D[D.rare&(D.period=="2023-2026")],"clean":D[D.rare&(D.period=="2023-2026")&D.clean],"subclean":D[D.rare&(D.period=="2023-2026")&D.subclean]}
res=[]
for k,G in sets.items():
    r=boot_ranks(G); r["set"]=k; r["auc"]=[auc(G.y.values,G[t].values.astype(float)) for t in tools]; r["rank"]=r.auc.rank(ascending=False); res.append(r)
Rk=pd.concat(res).reset_index().rename(columns={"index":"tool"}); Rk.to_csv(f"{P}/results/rank_bootstrap_rare_2023_26.tsv",sep="\t",index=False)
W=Rk.pivot(index="tool",columns="set",values="rank"); Lo=Rk.pivot(index="tool",columns="set",values="rank_lo"); Hi=Rk.pivot(index="tool",columns="set",values="rank_hi")
W["d_clean"]=W.full-W.clean; W["d_sub"]=W.full-W.subclean
# 'robust' shift: rank interval of filtered tier does not overlap the full-tier interval, in the same direction under both filters
W["robust_up"]=(Hi.clean<Lo.full)&(Hi.subclean<Lo.full); W["robust_down"]=(Lo.clean>Hi.full)&(Lo.subclean>Hi.full)
for c in ["full","clean","subclean"]: W[c+"_CI"]=Lo[c].astype(int).astype(str)+"-"+Hi[c].astype(int).astype(str)
pd.set_option("display.width",250)
print(W.sort_values("clean")[["full","full_CI","clean","clean_CI","subclean","subclean_CI","d_clean","d_sub","robust_up","robust_down"]].to_string())
rho,p=spearmanr(W.d_clean,W.d_sub); print("Spearman of rank shifts (clean vs subclean): rho=%.3f p=%.2g"%(rho,p))
print("same-direction shifts (|d|>=3 both):",((np.sign(W.d_clean)==np.sign(W.d_sub))&(W.d_clean.abs()>=3)&(W.d_sub.abs()>=3)).sum())
# (vi) submitter-level calibration
Th=pd.read_csv(f"{P}/meta/clingen_thresholds.tsv",sep="\t")
G=D[D.rare&D.subclean].copy(); G["per"]=np.where(G.period=="2023-2026","2023-2026","<=2022")
rows=[]
for (tool,rc),tt in Th.groupby(["tool","rawcol"]):
    for per,g in G.groupby("per"):
        g=g[g[rc].notna()]; n1=(g.y==1).sum(); n0=(g.y==0).sum()
        for _,r in tt[tt.evidence.str.startswith("PP3")].iterrows():
            inb=(g[rc]>=float(r.lo))&(g[rc]<float(r.hi)); a=(inb&(g.y==1)).sum(); b=(inb&(g.y==0)).sum()
            rows.append(dict(tool=tool,evidence=r.evidence,period=per,n_plp=int(n1),n_blb=int(n0),LR=((a+.5)/(n1+.5))/((b+.5)/(n0+.5))))
L=pd.DataFrame(rows); L.to_csv(f"{P}/results/calibration_submitter_level.tsv",sep="\t",index=False)
print(L[L.evidence=="PP3_Supporting"].pivot(index="tool",columns="period",values="LR").round(2).to_string()); print(L.groupby("period")[["n_plp","n_blb"]].max())
