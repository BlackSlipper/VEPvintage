#!/usr/bin/env python
"""Composition analyses for the vintage contrast on rare variants.
(A) Gene-identity null (after Harrizi et al. 2026): each variant is scored by the leave-one-out fraction of P/LP labels among
    the other labelled variants of its gene in the same label set (variants alone in their gene get the set's overall P/LP
    fraction); pooled AUROC of this null by period and tier, with gene-cluster bootstrap intervals, beside each core
    predictor's AUROC on the same variants.
(B) Composition-matched AUROC: on the gene x class cells present in both 2020-2022 and 2023-2026, the 2023-2026 labels are
    reweighted to the 2020-2022 gene x class composition (direct standardisation) and vice versa; paired gene-cluster
    bootstrap (B=1000, seed 2026) for the differences. Rare variants scored by all 38 core predictors."""
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
from multiprocessing import Pool
B=int(os.environ.get("VEP_BOOT","1000")); OUT=f"{P}/results/r2"; os.makedirs(OUT,exist_ok=True)
pd.set_option("display.width",250); pd.set_option("display.max_rows",300)
D=load()
T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","rare","clean","period"])
D=D.merge(T,on="hgvs38"); D["rare"]=D.rare.fillna(False).astype(bool); D["clean"]=D.clean.fillna(False).astype(bool)
tools=[t for t in TOOLS if t in D.columns]; cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]
R=D[D.rare].copy(); RC=R[R[core].notna().all(axis=1)].copy()
genes=pd.Index(sorted(D.Gene.dropna().unique())); ng=len(genes)
rng=np.random.default_rng(2026); W=np.stack([np.bincount(rng.integers(0,ng,ng),minlength=ng) for _ in range(B)]).astype(float)
ONE=np.ones(ng)
def ci(x): x=np.asarray(x,float); x=x[np.isfinite(x)]; return np.percentile(x,[2.5,97.5]) if len(x) else (np.nan,np.nan)
def pprep(G,t,wv=None):
    g=G[["Gene","y",t]].dropna(); lev=pd.Series(g[t].values).rank(method="dense").values.astype(int)-1
    w=np.ones(len(g)) if wv is None else np.asarray(wv,float)[G[t].notna().values]
    return lev,lev.max()+1,g.y.values.astype(float),genes.get_indexer(g.Gene.values),w
def pauc(pr,w):
    lev,K,y,gi,wv=pr; ww=w[gi]*wv; Pk=np.bincount(lev,weights=ww*y,minlength=K); Nk=np.bincount(lev,weights=ww*(1-y),minlength=K)
    cN=np.cumsum(Nk)-Nk; d=Pk.sum()*Nk.sum(); return (Pk*(cN+0.5*Nk)).sum()/d if d>0 else np.nan
PER=["<=2019","2020-2022","2023-2026"]
# ---------- (A) gene-identity null ----------
def add_null(G):
    G=G.copy(); gn=G.groupby("Gene").y.agg(["sum","count"]); s=G.Gene.map(gn["sum"]).values; n=G.Gene.map(gn["count"]).values
    prior=G.y.mean(); loo=np.where(n>1,(s-G.y.values)/np.maximum(n-1,1),prior); G["gene_null_loo"]=loo; G["gene_null_insample"]=s/n
    return G
rows=[]
for nm,G0 in [("all_common",RC),("all_own",R),("filtered_common",RC[RC.clean]),("filtered_own",R[R.clean])]:
    for per in PER:
        G=add_null(G0[G0.period==per])
        for t in ["gene_null_loo","gene_null_insample"]+core:
            pr=pprep(G,t); a=pauc(pr,ONE); bs=[pauc(pr,W[b]) for b in range(B)] if t.startswith("gene_null") else None
            r=dict(set=nm,period=per,tool=t,n=len(pr[2]),n_plp=int(pr[2].sum()),auc=a)
            if bs is not None: r["lo"],r["hi"]=ci(bs)
            if t.startswith("gene_null"): r["genes"]=G.Gene.nunique(); r["singleton_share"]=(G.groupby("Gene").y.transform("count")==1).mean()
            rows.append(r)
M=pd.DataFrame(rows); M.to_csv(f"{OUT}/M_gene_null.tsv",sep="\t",index=False)
print("gene-identity null AUROC by period"); print(M[M.tool.str.startswith("gene_null")].round(3).to_string())
nl=M[M.tool=="gene_null_loo"].set_index(["set","period"]).auc
for nm in ["all_common","filtered_common"]:
    sub=M[(M.set==nm)&~M.tool.str.startswith("gene_null")]
    piv=sub.pivot(index="tool",columns="period",values="auc"); 
    print(nm,"median predictor AUROC by period",piv.median().round(3).to_dict(),"null",{p:round(nl[(nm,p)],3) for p in PER})
    print(nm,"median margin over LOO null",{p:round((piv[p]-nl[(nm,p)]).median(),3) for p in PER})
# ---------- (B) composition matching ----------
A20=RC[RC.period=="2020-2022"].copy(); A23=RC[RC.period=="2023-2026"].copy()
def cells(G,by): return G.groupby(by).size()
res=[]
for mode,by in [("gene_class",["Gene","y"]),("gene",["Gene"])]:
    c20=cells(A20,by); c23=cells(A23,by); shared=c20.index.intersection(c23.index)
    k20=A20.set_index(by).index.isin(shared); k23=A23.set_index(by).index.isin(shared)
    S20=A20[k20].copy(); S23=A23[k23].copy()
    key20=pd.MultiIndex.from_frame(S20[by]) if len(by)>1 else pd.Index(S20[by[0]]); key23=pd.MultiIndex.from_frame(S23[by]) if len(by)>1 else pd.Index(S23[by[0]])
    w23to20=(c20.reindex(key23).values/c23.reindex(key23).values)   # 2023-2026 reweighted to the 2020-2022 composition
    w20to23=(c23.reindex(key20).values/c20.reindex(key20).values)   # 2020-2022 reweighted to the 2023-2026 composition
    print(mode,"shared cells",len(shared),"2020-22 kept",len(S20),"of",len(A20),"2023-26 kept",len(S23),"of",len(A23),
          "P/LP share 2020-22 S",round(S20.y.mean(),3),"2023-26 S",round(S23.y.mean(),3),"effective n 2023->2020",round(w23to20.sum()**2/(w23to20**2).sum()),"2020->2023",round(w20to23.sum()**2/(w20to23**2).sum()))
    def job(t):
        p20=pprep(A20,t); p23=pprep(A23,t); s20=pprep(S20,t); s23=pprep(S23,t); m23=pprep(S23,t,w23to20); m20=pprep(S20,t,w20to23)
        est={k:pauc(p,ONE) for k,p in [("full20",p20),("full23",p23),("S20",s20),("S23",s23),("S23_to20",m23),("S20_to23",m20)]}
        bs={k:np.array([pauc(p,W[b]) for b in range(B)]) for k,p in [("full20",p20),("full23",p23),("S20",s20),("S23",s23),("S23_to20",m23),("S20_to23",m20)]}
        r=dict(mode=mode,tool=t,**{k:v for k,v in est.items()})
        for nm,a,b_ in [("d_full","full23","full20"),("d_S","S23","S20"),("d_matched_to20","S23_to20","S20"),("d_matched_to23","S23","S20_to23")]:
            r[nm]=est[a]-est[b_]; r[nm+"_lo"],r[nm+"_hi"]=ci(bs[a]-bs[b_])
        return r
    with Pool(20) as pool: res+=pool.map(job,core,chunksize=2)
N=pd.DataFrame(res); N.to_csv(f"{OUT}/N_composition_matched.tsv",sep="\t",index=False)
for mode in ["gene_class","gene"]:
    sub=N[N["mode"]==mode]; print("==",mode)
    print(sub[["full20","full23","S20","S23","S23_to20","S20_to23","d_full","d_S","d_matched_to20","d_matched_to23"]].median().round(3).to_string())
    print("CI>0 counts:",{k:int((sub[k+"_lo"]>0).sum()) for k in ["d_full","d_S","d_matched_to20","d_matched_to23"]},"CI<0:",{k:int((sub[k+"_hi"]<0).sum()) for k in ["d_full","d_S","d_matched_to20","d_matched_to23"]})
    print(sub[sub.tool.isin(["REVEL","AlphaMissense","ESM1b","GERP++","phyloP100way","CADD","MetaRNN","PrimateAI"])][["tool","full20","full23","S20","S23","S23_to20","d_matched_to20","d_matched_to20_lo","d_matched_to20_hi","d_matched_to23"]].round(3).to_string())
print("DONE")
