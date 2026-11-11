#!/usr/bin/env python
"""(A) Variant inclusion flow for the evaluation set (counts at each step).
(B) Paired gene-cluster bootstrap (B=1000, seed 2026) for the change in AUROC between 2020-2022 and 2023-2026 labels,
    pooled and within-gene, on rare variants: (i) each predictor on the rare variants it scores, (ii) the rare variants
    scored by all 38 core predictors; for all labels and the provenance-filtered tier.
(C) Clopper-Pearson intervals for the text-mining audit (30/30 flagged, 1/40 unflagged)."""
import sys, os, gzip; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
from multiprocessing import Pool
from scipy.stats import beta
B=int(os.environ.get("VEP_BOOT","1000")); OUT=f"{P}/results/r2"; os.makedirs(OUT,exist_ok=True)
pd.set_option("display.width",250); pd.set_option("display.max_rows",300)
# (A) flow
rows=[]
snap=f"{P}/data/snap_parsed/current.tsv.gz"   # parsed November 2026 snapshot (script 01); only the first count needs it
if os.path.exists(snap):
    cur=pd.read_csv(snap,sep="\t",dtype=str,usecols=["AlleleID"])
    rows.append(("Missense SNVs with any ClinVar classification, November 2026",cur.AlleleID.nunique()))
else:
    rows.append(("Missense SNVs with any ClinVar classification, November 2026 (needs data/snap_parsed/current.tsv.gz from script 01)",np.nan))
V=pd.read_csv(f"{P}/results/variant_history.tsv.gz",sep="\t",dtype={"VariationID":str,"Chromosome":str})
pb=V[V.final.isin(["PLP","BLB"])]; rows.append(("  classified P/LP or B/LB",pb.AlleleID.nunique()))
V0=pb[pb.stars>=1]; rows.append(("  with at least one review star",V0.AlleleID.nunique()))
V2=V0[V0.hgvs38.notna()]; rows.append(("  with single-base GRCh38 reference and alternate alleles",len(V2)))
V1=V2[V2.OriginSimple.fillna("").str.contains("germline")]; rows.append(("  and germline origin",len(V1)))
D=load(); rows.append(("  with a dbNSFP 4.8a record (evaluation set)",len(D)))
T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","rare","clean","period"])
D=D.merge(T,on="hgvs38"); D["rare"]=D.rare.fillna(False).astype(bool); D["clean"]=D.clean.fillna(False).astype(bool)
tools=[t for t in TOOLS if t in D.columns]; cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]
R=D[D.rare]; rows.append(("    rare (max gnomAD v2.1.1/ExAC AF < 1e-4)",len(R)))
rows.append(("      P/LP",int((R.final=="PLP").sum()))); rows.append(("      B/LB",int((R.final=="BLB").sum())))
RC=R[R[core].notna().all(axis=1)]; rows.append(("    rare and scored by all 38 core predictors",len(RC)))
for per in ["<=2019","2020-2022","2023-2026"]:
    rows.append((f"      labelled {per}",int((RC.period==per).sum())))
rows.append(("    rare, provenance-filtered tier",int(R.clean.sum()))); rows.append(("    rare, 2023-2026, provenance-filtered",int((R.clean&(R.period=="2023-2026")).sum())))
F=pd.DataFrame(rows,columns=["step","n"]); F.to_csv(f"{OUT}/K_inclusion_flow.tsv",sep="\t",index=False); print(F.to_string())
# (B) period contrast
genes=pd.Index(sorted(D.Gene.dropna().unique())); ng=len(genes)
rng=np.random.default_rng(2026); W=np.stack([np.bincount(rng.integers(0,ng,ng),minlength=ng) for _ in range(B)]).astype(float)
def gstats(G,t):
    g=G[["Gene","y",t]].dropna(); gi=genes.get_indexer(g.Gene.values); r=g.groupby("Gene")[t].rank().values; y=g.y.values.astype(float)
    n1=np.bincount(gi,weights=y,minlength=ng); n=np.bincount(gi,minlength=ng).astype(float); rs=np.bincount(gi,weights=r*y,minlength=ng)
    return rs-n1*(n1+1)/2, n1*(n-n1)
def pprep(G,t):
    g=G[["Gene","y",t]].dropna(); lev=pd.Series(g[t].values).rank(method="dense").values.astype(int)-1
    return lev,lev.max()+1,g.y.values.astype(float),genes.get_indexer(g.Gene.values)
def pauc(pr,w):
    lev,K,y,gi=pr; ww=w[gi]; Pk=np.bincount(lev,weights=ww*y,minlength=K); Nk=np.bincount(lev,weights=ww*(1-y),minlength=K)
    cN=np.cumsum(Nk)-Nk; return (Pk*(cN+0.5*Nk)).sum()/(Pk.sum()*Nk.sum())
ONE=np.ones(ng)
def ci(x): x=np.asarray(x,float); x=x[np.isfinite(x)]; return np.percentile(x,[2.5,97.5])
SETS={}
for nm,G in [("all_own",R),("all_common",RC),("filtered_own",R[R.clean]),("filtered_common",RC[RC.clean])]:
    for per in ["2020-2022","2023-2026"]: SETS[(nm,per)]=G[G.period==per]
def job(k):
    nm,t=k; out={}
    for per in ["2020-2022","2023-2026"]:
        G=SETS[(nm,per)]; pr=pprep(G,t); U,pa=gstats(G,t)
        out[per]=(pauc(pr,ONE),np.array([pauc(pr,W[b]) for b in range(B)]),U.sum()/pa.sum(),(W@U)/(W@pa),len(pr[2]),int(pr[2].sum()))
    a,b_=out["2023-2026"],out["2020-2022"]
    lp,hp=ci(a[1]-b_[1]); lw,hw=ci(a[3]-b_[3]); lb,hb=ci((a[1]-b_[1])-(a[3]-b_[3]))
    return dict(set=nm,tool=t,n_2020_2022=b_[4],n_2023_2026=a[4],pooled_2020_2022=b_[0],pooled_2023_2026=a[0],d_pooled=a[0]-b_[0],dp_lo=lp,dp_hi=hp,
                within_2020_2022=b_[2],within_2023_2026=a[2],d_within=a[2]-b_[2],dw_lo=lw,dw_hi=hw,
                d_between=(a[0]-b_[0])-(a[2]-b_[2]),db_lo=lb,db_hi=hb,within_share=(a[2]-b_[2])/(a[0]-b_[0]) if abs(a[0]-b_[0])>1e-9 else np.nan)
KEYS=[(nm,t) for nm in ["all_own","all_common","filtered_own","filtered_common"] for t in core]
with Pool(20) as pool: res=pool.map(job,KEYS,chunksize=2)
J=pd.DataFrame(res); J.to_csv(f"{OUT}/J_period_contrast_2020_22_vs_2023_26.tsv",sep="\t",index=False)
g=J.groupby("set")
print(pd.concat([g[["d_pooled","d_within","d_between","within_share"]].median(),g.apply(lambda d:(d.db_lo>0).sum()).rename("n_between_CI_gt0"),g.d_pooled.min().rename("dp_min"),g.d_pooled.max().rename("dp_max"),g.apply(lambda d:(d.dp_lo>0).sum()).rename("n_pooled_CI_gt0"),
                 g.apply(lambda d:(d.dw_lo>0).sum()).rename("n_within_CI_gt0"),g.apply(lambda d:(d.dw_hi<0).sum()).rename("n_within_CI_lt0"),g.apply(lambda d:((d.dw_lo<0)&(d.dw_hi>0)).sum()).rename("n_within_CI_incl0")],axis=1).round(3).to_string())
print(J[J.tool.isin(["REVEL","AlphaMissense","ESM1b","GERP++","phyloP100way","CADD","MetaRNN"])&(J.set=="all_common")].round(3).to_string())
# (C) audit intervals
def cp(k,n): lo=beta.ppf(0.025,k,n-k+1) if k>0 else 0.0; hi=beta.ppf(0.975,k+1,n-k) if k<n else 1.0; return lo,hi
print("audit flagged 30/30 mention computational evidence: 95% CI",np.round(cp(30,30),3)); print("unflagged 1/40 mention:",np.round(cp(1,40),3))
print("DONE")
