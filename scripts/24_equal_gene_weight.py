#!/usr/bin/env python
"""Sensitivity: equal-gene-weight within-gene AUROC (mean of per-gene AUROCs over genes with >=3 P/LP and >=3 B/LB),
versus the pair-weighted within-gene AUROC of 22_round2.py. Rare variants; full / clean / subclean tiers by period;
full-minus-filtered differences (own genes and genes qualifying in both tiers) and the 2023-26 minus <=2019 jump.
Paired gene-cluster bootstrap with the same weight matrix as 22_round2.py (seed 2026)."""
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
B=int(os.environ.get("VEP_BOOT","1000")); OUT=f"{P}/results/r2"; os.makedirs(OUT,exist_ok=True); MINC=3
pd.set_option("display.width",250); pd.set_option("display.max_rows",300)
D=load()
T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","clean","rare","period"])
D=D.merge(T,on="hgvs38")
for c in ["clean","rare"]: D[c]=D[c].fillna(False).astype(bool)
SA="Ambry Genetics"
S=pd.read_csv(f"{P}/results/scv_textmining.tsv.gz",sep="\t",dtype={"VariationID":str},usecols=["VariationID","SCV","Submitter","ClinSig","has_text","any_insilico","year"])
C=pd.read_csv(f"{P}/results/scv_acmg_codes.tsv.gz",sep="\t",dtype={"VariationID":str})
st=S[S.has_text==1].groupby("Submitter").agg(n_text=("SCV","size"),frac_comp=("any_insilico","mean"))
cc=C.assign(comp=C.codes.str.contains(r"\bPP3\b|\bBP4\b")).groupby("Submitter").comp.mean().rename("fcc")
st=st.join(cc).fillna({"fcc":0}); noncomp=set(st[(st.n_text>=20)&(st.frac_comp<0.05)&(st.fcc<0.05)].index)-{SA}
cs=S.ClinSig.str.lower(); S["cl"]=np.where(cs.str.contains("pathogenic"),"PLP",np.where(cs.str.contains("benign"),"BLB","O"))
S=S.merge(D[["VariationID","final","vyear"]].drop_duplicates("VariationID"),on="VariationID"); S=S[S.cl==S.final]
D["subclean"]=D.VariationID.isin(set(S[S.Submitter.isin(noncomp)].VariationID))
tools=[t for t in TOOLS if t in D.columns]; cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]
R=D[D.rare].copy(); print("rare",len(R),"core",len(core),"noncomp",len(noncomp))
genes=pd.Index(sorted(D.Gene.dropna().unique())); ng=len(genes)
rng=np.random.default_rng(2026); W=np.stack([np.bincount(rng.integers(0,ng,ng),minlength=ng) for _ in range(B)]).astype(float)
def gstats(G,t):
    g=G[["Gene","y",t]].dropna(); gi=genes.get_indexer(g.Gene.values)
    r=g.groupby("Gene")[t].rank().values; y=g.y.values.astype(float)
    n1=np.bincount(gi,weights=y,minlength=ng); n=np.bincount(gi,minlength=ng).astype(float)
    rs=np.bincount(gi,weights=r*y,minlength=ng)
    return rs-n1*(n1+1)/2, n1*(n-n1), n1, n-n1
def ci(x): x=np.asarray(x,float); x=x[np.isfinite(x)]; return np.percentile(x,[2.5,97.5]) if len(x) else (np.nan,np.nan)
def eqw(U,pa,m):  # point and bootstrap equal-weight mean of per-gene AUROC over mask m
    a=np.where(m,U/np.where(pa>0,pa,1),0.0); return a[m].mean(), (W[:,m]@a[m])/W[:,m].sum(1)
per_list=["<=2019","2020-2022","2023-2026"]; ST={}; rows=[]
for per in per_list:
    G=R[R.period==per]
    for s,GG in [("full",G),("clean",G[G.clean]),("subclean",G[G.subclean])]:
        for t in core:
            U,pa,n1,n0=gstats(GG,t); m=(n1>=MINC)&(n0>=MINC); ST[(s,per,t)]=(U,pa,m)
            e,eb=eqw(U,pa,m); lo,hi=ci(eb)
            rows.append(dict(set=s,period=per,tool=t,n_genes_eq=int(m.sum()),within_pairw=U.sum()/pa.sum(),within_eqw=e,eqw_lo=lo,eqw_hi=hi))
G1=pd.DataFrame(rows); G1.to_csv(f"{OUT}/G_equal_gene_weight.tsv",sep="\t",index=False)
print("=== medians over core ==="); print(G1.groupby(["set","period"])[["n_genes_eq","within_pairw","within_eqw"]].median().round(3).to_string())
rows=[]
for per in per_list:
    for x in ["clean","subclean"]:
        for t in core:
            Uf,pf,mf=ST[("full",per,t)]; Uc,pc,mc=ST[(x,per,t)]
            ef,bf=eqw(Uf,pf,mf); ec,bc=eqw(Uc,pc,mc); m=mf&mc
            gf,gbf=eqw(Uf,pf,m); gc,gbc=eqw(Uc,pc,m)
            a,b=ci(bf-bc); c,d=ci(gbf-gbc)
            rows.append(dict(period=per,filter=x,tool=t,d_eqw=ef-ec,d_eqw_lo=a,d_eqw_hi=b,n_genes_both=int(m.sum()),d_eqw_common=gf-gc,d_eqwc_lo=c,d_eqwc_hi=d))
G2=pd.DataFrame(rows); G2.to_csv(f"{OUT}/G_equal_gene_weight_diff.tsv",sep="\t",index=False)
g=G2.groupby(["filter","period"])
print("=== full minus filtered (equal gene weight) ===")
print(pd.concat([g[["d_eqw","d_eqw_common","n_genes_both"]].median(),g.apply(lambda d:(d.d_eqw_lo>0).sum()).rename("n_CI_gt0"),
                 g.apply(lambda d:(d.d_eqwc_lo>0).sum()).rename("n_common_CI_gt0"),g.apply(lambda d:(d.d_eqwc_hi<0).sum()).rename("n_common_CI_lt0")],axis=1).round(3).to_string())
rows=[]
for x in ["full","clean","subclean"]:
    for t in core:
        Ua,pa_,ma=ST[(x,"2023-2026",t)]; Ub,pb,mb=ST[(x,"<=2019",t)]
        ea,ba=eqw(Ua,pa_,ma); eb_,bb=eqw(Ub,pb,mb); lo,hi=ci(ba-bb)
        rows.append(dict(set=x,tool=t,jump_eqw=ea-eb_,lo=lo,hi=hi))
G3=pd.DataFrame(rows); G3.to_csv(f"{OUT}/G_equal_gene_weight_jump.tsv",sep="\t",index=False)
g=G3.groupby("set"); print("=== jump 2023-26 minus <=2019 (equal gene weight) ===")
print(pd.concat([g.jump_eqw.median(),g.apply(lambda d:(d.lo>0).sum()).rename("n_CI_gt0"),g.apply(lambda d:(d.hi<0).sum()).rename("n_CI_lt0")],axis=1).round(3).to_string())
key=["REVEL","AlphaMissense","ESM1b","CADD","MetaRNN","BayesDel_noAF","VEST4"]
print(G2[G2.tool.isin(key)&(G2.period=="2023-2026")].round(3).to_string())
