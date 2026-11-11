#!/usr/bin/env python
"""Round-2 follow-ups: (a) MutPred2 and EVE ranked with the 38 core predictors on their own common subsets;
(b) gene composition of rare B/LB labels by period and source; (c) check of period-level AUROC changes on the rare common subset."""
import sys, os; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
from scipy.stats import rankdata
B=int(__import__("os").environ.get("VEP_BOOT","1000")); OUT=f"{P}/results/r2"
pd.set_option("display.width",250)
D=load()
T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","clean","rare","period"])
D=D.merge(T,on="hgvs38"); D["clean"]=D.clean.fillna(False).astype(bool)
A_=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",usecols=["hgvs38","conc_submitters"]).drop_duplicates("hgvs38")
D=D.merge(A_,on="hgvs38",how="left")
SA="Ambry Genetics"; SB="Labcorp Genetics (formerly Invitae), Labcorp"
tools=[t for t in TOOLS if t in D.columns]; cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]
R=D[D.rare].copy()
# subclean (same definition as scripts 17/19)
S=pd.read_csv(f"{P}/results/scv_textmining.tsv.gz",sep="\t",dtype={"VariationID":str},usecols=["VariationID","SCV","Submitter","ClinSig","has_text","any_insilico"])
C=pd.read_csv(f"{P}/results/scv_acmg_codes.tsv.gz",sep="\t",dtype={"VariationID":str})
st=S[S.has_text==1].groupby("Submitter").agg(n_text=("SCV","size"),frac_comp=("any_insilico","mean"))
cc=C.assign(comp=C.codes.str.contains(r"\bPP3\b|\bBP4\b")).groupby("Submitter").comp.mean().rename("fcc")
st=st.join(cc).fillna({"fcc":0}); noncomp=set(st[(st.n_text>=20)&(st.frac_comp<0.05)&(st.fcc<0.05)].index)-{SA}
cs=S.ClinSig.str.lower(); S["cl"]=np.where(cs.str.contains("pathogenic"),"PLP",np.where(cs.str.contains("benign"),"BLB","O"))
S=S.merge(D[["VariationID","final"]].drop_duplicates("VariationID"),on="VariationID"); S=S[S.cl==S.final]
R["subclean"]=R.VariationID.isin(set(S[S.Submitter.isin(noncomp)].VariationID))
genes=pd.Index(sorted(D.Gene.dropna().unique())); ng=len(genes)
rng=np.random.default_rng(2026); W=np.stack([np.bincount(rng.integers(0,ng,ng),minlength=ng) for _ in range(B)]).astype(float)
def pauc(G,t,w):
    g=G[["Gene","y",t]].dropna(); lev=pd.Series(g[t].values).rank(method="dense").values.astype(int)-1; K=lev.max()+1
    y=g.y.values.astype(float); ww=w[genes.get_indexer(g.Gene.values)]
    Pk=np.bincount(lev,weights=ww*y,minlength=K); Nk=np.bincount(lev,weights=ww*(1-y),minlength=K); cN=np.cumsum(Nk)-Nk
    return (Pk*(cN+0.5*Nk)).sum()/(Pk.sum()*Nk.sum())
# (a)
rows=[]
for extra in ["MutPred2","EVE"]:
    tl=core+[extra]; Cc=R[(R.period=="2023-2026")&R[tl].notna().all(axis=1)]
    for s,G in {"full":Cc,"clean":Cc[Cc.clean],"subclean":Cc[Cc.subclean]}.items():
        pt=np.array([pauc(G,t,np.ones(ng)) for t in tl]); M=np.array([[pauc(G,t,W[b]) for t in tl] for b in range(B)])
        rk=rankdata(-M,axis=1); i=len(tl)-1
        rows.append(dict(extra=extra,set=s,n=len(G),n_plp=int(G.y.sum()),auc=pt[i],rank=rankdata(-pt)[i],rank_lo=np.percentile(rk[:,i],2.5),rank_hi=np.percentile(rk[:,i],97.5),
                         REVEL_auc=pt[tl.index("REVEL")],REVEL_rank=rankdata(-pt)[tl.index("REVEL")],of=len(tl)))
X=pd.DataFrame(rows); X.to_csv(f"{OUT}/B2_extra_tools_ranks.tsv",sep="\t",index=False); print(X.round(3).to_string())
# (b) gene composition
nplp=R[R.y==1].groupby("Gene").size()
R["gene_nplp"]=R.Gene.map(nplp).fillna(0)
R["src"]=np.where(R.conc_submitters.fillna("")==SA,"A_only",np.where(R.conc_submitters.fillna("")==SB,"B_only","other"))
b=R[R.y==0]
out=b.groupby(["period","src"]).agg(n=("y","size"),no_plp_gene=("gene_nplp",lambda s:(s==0).mean()),lt5_plp_gene=("gene_nplp",lambda s:(s<5).mean()))
out2=b.groupby("period").agg(n=("y","size"),no_plp_gene=("gene_nplp",lambda s:(s==0).mean()),lt5_plp_gene=("gene_nplp",lambda s:(s<5).mean()))
out3=b[b.clean].groupby("period").agg(n=("y","size"),no_plp_gene=("gene_nplp",lambda s:(s==0).mean()))
print(out.round(3).to_string()); print(out2.round(3).to_string()); print("clean B/LB:\n",out3.round(3).to_string())
pd.concat([out.reset_index(),out2.reset_index().assign(src="all"),out3.reset_index().assign(src="clean_all")]).to_csv(f"{OUT}/B3_gene_composition.tsv",sep="\t",index=False)
# (c) period changes on the rare common subset
Cc=R[R[core].notna().all(axis=1)]
rows=[]
for per,g in Cc.groupby("period"):
    for t in core: rows.append(dict(period=per,tool=t,auc=auc(g.y.values,g[t].values.astype(float)),n=len(g)))
Pp=pd.DataFrame(rows).pivot(index="tool",columns="period",values="auc"); d=Pp["2023-2026"]-Pp["2020-2022"]
print("median by period:",Pp.median().round(3).to_dict(),"| 2023-26 minus 2020-22: n>0",int((d>0).sum()),"of",len(d),"range %.3f-%.3f median %.3f"%(d.min(),d.max(),d.median()))
print(Pp.loc[["GERP++","phyloP100way","REVEL","ESM1b","SiPhy","MetaRNN"]].round(3).to_string())
B_=pd.read_csv(f"{P}/results/benchmark_clean_rare_recent.tsv",sep="\t"); b2=B_[B_.tool.isin(core)&(B_.set=="full_rare_2023_26")]
print("rare 2023-26 own-variant AUROC range (core):",round(b2.auc.min(),3),b2.loc[b2.auc.idxmin(),"tool"],round(b2.auc.max(),3),b2.loc[b2.auc.idxmax(),"tool"])
w=B_.pivot(index="tool",columns="set",values="auc").loc[core]; dd=w.full_rare_2023_26-w.clean_rare_2023_26
print("full-clean drop over core: %.3f-%.3f, median %.3f"%(dd.min(),dd.max(),dd.median()), dd.idxmin(), dd.idxmax())
w2=B_.pivot(index="tool",columns="set",values="auc"); dd2=w2.full_rare_2023_26-w2.clean_rare_2023_26; print("all 44: %.3f-%.3f"%(dd2.min(),dd2.max()),dd2.idxmin(),dd2.idxmax())
print("DONE")
