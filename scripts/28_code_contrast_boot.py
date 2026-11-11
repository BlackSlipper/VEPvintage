#!/usr/bin/env python
"""Paired gene-cluster bootstrap (B=1000, seed 2026) for the evidence-code contrasts on rare variants scored by all
38 core predictors: (i) rare B/LB labels citing BP4 vs rare B/LB labels citing only other codes, each against all rare
P/LP variants; (ii) rare P/LP labels citing PP3 vs rare P/LP labels citing only other codes, each against all rare B/LB."""
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
from multiprocessing import Pool
B=int(os.environ.get("VEP_BOOT","1000")); OUT=f"{P}/results/r2"
D=load()
A=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",usecols=["hgvs38","rare","code_citing","has_PP3","has_BP4"]).drop_duplicates("hgvs38")
D=D.merge(A,on="hgvs38"); D["rare"]=D.rare.fillna(False).astype(bool)
for c in ["code_citing","has_PP3","has_BP4"]: D[c]=D[c].fillna(False).astype(bool)
tools=[t for t in TOOLS if t in D.columns]; cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]
R=D[D.rare]; RC=R[R[core].notna().all(axis=1)]
genes=pd.Index(sorted(D.Gene.dropna().unique())); ng=len(genes)
rng=np.random.default_rng(2026); W=np.stack([np.bincount(rng.integers(0,ng,ng),minlength=ng) for _ in range(B)]).astype(float)
def pprep(G,t):
    g=G[["Gene","y",t]].dropna(); lev=pd.Series(g[t].values).rank(method="dense").values.astype(int)-1
    return lev,lev.max()+1,g.y.values.astype(float),genes.get_indexer(g.Gene.values)
def pauc(pr,w):
    lev,K,y,gi=pr; ww=w[gi]; Pk=np.bincount(lev,weights=ww*y,minlength=K); Nk=np.bincount(lev,weights=ww*(1-y),minlength=K)
    cN=np.cumsum(Nk)-Nk; return (Pk*(cN+0.5*Nk)).sum()/(Pk.sum()*Nk.sum())
ONE=np.ones(ng)
def ci(x): x=np.asarray(x,float); x=x[np.isfinite(x)]; return np.percentile(x,[2.5,97.5])
plp=RC[RC.y==1]; blb=RC[RC.y==0]
SETS={"BP4":(pd.concat([plp,blb[blb.code_citing&blb.has_BP4]]),pd.concat([plp,blb[blb.code_citing&~blb.has_BP4]])),
      "PP3":(pd.concat([plp[plp.code_citing&plp.has_PP3],blb]),pd.concat([plp[plp.code_citing&~plp.has_PP3],blb]))}
for k,(a,b) in SETS.items(): print(k,"with:",len(a),int(a.y.sum()),"without:",len(b),int(b.y.sum()))
def job(k):
    side,t=k; a,b=SETS[side]; pa,pb=pprep(a,t),pprep(b,t)
    ba=np.array([pauc(pa,W[i]) for i in range(B)]); bb=np.array([pauc(pb,W[i]) for i in range(B)])
    lo,hi=ci(ba-bb); return dict(side=side,tool=t,auc_with=pauc(pa,ONE),auc_without=pauc(pb,ONE),gap=pauc(pa,ONE)-pauc(pb,ONE),gap_lo=lo,gap_hi=hi)
with Pool(20) as pool: res=pool.map(job,[(s,t) for s in SETS for t in core],chunksize=2)
J=pd.DataFrame(res); J.to_csv(f"{OUT}/L_code_contrast_bootstrap.tsv",sep="\t",index=False)
g=J.groupby("side"); print(pd.concat([g.gap.median(),g.gap.min().rename("gap_min"),g.gap.max().rename("gap_max"),g.apply(lambda d:(d.gap_lo>0).sum()).rename("n_CI_gt0")],axis=1).round(3).to_string())
print(J[J.tool.isin(["REVEL","MetaRNN","PrimateAI","GERP++"])].round(3).to_string()); print("DONE")
