#!/usr/bin/env python
"""Benchmarks as they would have been built on a given ClinVar download date.
For each January release 2016-2026 and November 2026: missense SNVs classified P/LP or B/LB (>=1 star, germline) in that
release, rare (max gnomAD v2.1.1/ExAC AF < 1e-4). Pooled and within-gene AUROC of the 38 core predictors; share of rare
B/LB labels in genes without a rare P/LP label; Kendall tau of rankings (variants scored by all core predictors) against
the November 2026 ranking; labels of each release that no longer hold in November 2026. Gene-cluster bootstrap (seed 2026,
genes resampled jointly across releases) for the change from January 2020 to November 2026."""
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa
from scipy.stats import kendalltau
B=int(os.environ.get("VEP_BOOT","1000")); OUT=f"{P}/results/r2"; os.makedirs(OUT,exist_ok=True)
pd.set_option("display.width",250); pd.set_option("display.max_rows",200)
D0=load(); tools=[t for t in TOOLS if t in D0.columns]; cov=D0[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]
print("core",len(core))
V=pd.read_csv(f"{P}/results/variant_history.tsv.gz",sep="\t",usecols=["AlleleID","hgvs38","Gene"])
X=pd.read_csv(f"{P}/results/extra_alleles_hg38.tsv",sep="\t") if os.path.exists(f"{P}/results/extra_alleles_hg38.tsv") else pd.DataFrame(columns=["AlleleID","hgvs38"])
A2H=pd.concat([V[["AlleleID","hgvs38"]],X[["AlleleID","hgvs38"]]]).drop_duplicates("AlleleID").set_index("AlleleID").hgvs38
GENE=V.drop_duplicates("AlleleID").set_index("AlleleID").Gene
fs=sorted(glob.glob(f"{P}/results/dbnsfp_shard*.tsv.gz"))+sorted(glob.glob(f"{P}/results/dbnsfp_extra*.tsv.gz"))
SC=pd.concat([pd.read_csv(f,sep="\t",dtype={"query":str,"uniprot":str,"aapos":str}) for f in fs]).drop_duplicates("query")
SC=SC[SC.found==1][["query"]+core]
AF=pd.concat([pd.read_csv(f,sep="\t") for f in sorted(glob.glob(f"{P}/results/af_shard*.tsv.gz"))+sorted(glob.glob(f"{P}/results/af_extra*.tsv.gz"))]).drop_duplicates("query")
for c in AF.columns[1:]: AF[c]=pd.to_numeric(AF[c],errors="coerce")
AF["maxaf"]=AF[["gnomad_ex_popmax","gnomad_ge_popmax","gnomad_ex_af","gnomad_ge_af","exac_af"]].max(axis=1).fillna(0)
AFm=AF.set_index("query").maxaf
LAB=pd.read_csv(f"{P}/results/download_date_labels.tsv.gz",sep="\t",dtype={"Gene":str})  # from 25b (or data/ via 00)
snaps=[f"{y}-01" for y in range(2016,2027)]+["current"]
SETS={}; rows_n=[]
for s in snaps:
    d=LAB[LAB.release==("2026-11" if s=="current" else s)].copy()
    n_all=len(d)
    d["query"]=d.AlleleID.map(A2H); d["Gene2"]=d.AlleleID.map(GENE).fillna(d.Gene)
    d=d.merge(SC,on="query",how="inner"); d["maxaf"]=d["query"].map(AFm)
    n_sc=len(d); n_af=int(d.maxaf.notna().sum())
    d=d[d.maxaf.fillna(1)<1e-4].copy(); d["y"]=(d.Class=="PLP").astype(int); d["Gene"]=d.Gene2
    d["held"]=d.class_2026_11==d.Class
    SETS[s]=d
    nplp=d[d.y==1].groupby("Gene").size(); b=d[d.y==0]
    rows_n.append(dict(release=("2026-11" if s=="current" else s),n_labels=n_all,n_scored=n_sc,frac_scored=n_sc/n_all,n_rare=len(d),n_rare_plp=int(d.y.sum()),
                       n_rare_blb=int((1-d.y).sum()),blb_no_plp_gene=(~b.Gene.isin(nplp.index)).mean(),frac_changed_by_2026_11=1-d.held.mean()))
N=pd.DataFrame(rows_n); print(N.round(3).to_string())
genes=pd.Index(sorted(set().union(*[set(SETS[s].Gene.dropna()) for s in snaps]))); ng=len(genes)
rng=np.random.default_rng(2026); W=np.stack([np.bincount(rng.integers(0,ng,ng),minlength=ng) for _ in range(B)]).astype(float)
def pprep(G,t):
    g=G[["Gene","y",t]].dropna(); lev=pd.Series(g[t].values).rank(method="dense").values.astype(int)-1
    return lev,lev.max()+1,g.y.values.astype(float),genes.get_indexer(g.Gene.values)
def pauc(pr,w):
    lev,K,y,gi=pr; ww=w[gi]; Pk=np.bincount(lev,weights=ww*y,minlength=K); Nk=np.bincount(lev,weights=ww*(1-y),minlength=K)
    cN=np.cumsum(Nk)-Nk; return (Pk*(cN+0.5*Nk)).sum()/(Pk.sum()*Nk.sum())
def gstats(G,t):
    g=G[["Gene","y",t]].dropna(); gi=genes.get_indexer(g.Gene.values); r=g.groupby("Gene")[t].rank().values; y=g.y.values.astype(float)
    n1=np.bincount(gi,weights=y,minlength=ng); n=np.bincount(gi,minlength=ng).astype(float); rs=np.bincount(gi,weights=r*y,minlength=ng)
    return rs-n1*(n1+1)/2, n1*(n-n1)
ONE=np.ones(ng); rows=[]; BP={}; BW={}
for s in snaps:
    G=SETS[s]; rel="2026-11" if s=="current" else s
    for t in core:
        pr=pprep(G,t); U,pa=gstats(G,t)
        rows.append(dict(release=rel,tool=t,n=len(pr[2]),pooled=pauc(pr,ONE),within=U.sum()/pa.sum(),genes_with_pairs=int((pa>0).sum())))
        if s in ("2020-01","current"):
            BP[(s,t)]=np.array([pauc(pr,W[b]) for b in range(B)]); BW[(s,t)]=(W@U)/(W@pa)
R=pd.DataFrame(rows); R.to_csv(f"{OUT}/H_download_date_auc.tsv",sep="\t",index=False)
def ci(x): x=np.asarray(x,float); x=x[np.isfinite(x)]; return np.percentile(x,[2.5,97.5])
rows=[]
for t in core:
    a=R[(R.release=="2020-01")&(R.tool==t)].iloc[0]; c=R[(R.release=="2026-11")&(R.tool==t)].iloc[0]
    lp,hp=ci(BP[("current",t)]-BP[("2020-01",t)]); lw,hw=ci(BW[("current",t)]-BW[("2020-01",t)])
    rows.append(dict(tool=t,pooled_2020_01=a.pooled,pooled_2026_11=c.pooled,d_pooled=c.pooled-a.pooled,dp_lo=lp,dp_hi=hp,
                     within_2020_01=a.within,within_2026_11=c.within,d_within=c.within-a.within,dw_lo=lw,dw_hi=hw))
C=pd.DataFrame(rows); C.to_csv(f"{OUT}/H_download_date_change.tsv",sep="\t",index=False)
# rankings on variants scored by all core predictors, vs November 2026
def rank_on_common(G):
    g=G[G[core].notna().all(axis=1)]; pr={t:pprep(g,t) for t in core}
    pooled=pd.Series({t:pauc(pr[t],ONE) for t in core}); w={}
    for t in core: U,pa=gstats(g,t); w[t]=U.sum()/pa.sum()
    return pooled,pd.Series(w),len(g)
refp,refw,_=rank_on_common(SETS["current"]); rows=[]
for s in snaps:
    p_,w_,n=rank_on_common(SETS[s])
    rows.append(dict(release="2026-11" if s=="current" else s,n_common=n,median_pooled_common=p_.median(),median_within_common=w_.median(),
                     tau_pooled_vs_2026_11=kendalltau(p_,refp)[0],tau_within_vs_2026_11=kendalltau(w_,refw)[0],
                     max_rank_shift_pooled=int((p_.rank(ascending=False)-refp.rank(ascending=False)).abs().max())))
K=pd.DataFrame(rows); K.to_csv(f"{OUT}/H_download_date_ranks.tsv",sep="\t",index=False)
N.to_csv(f"{OUT}/H_download_date_sets.tsv",sep="\t",index=False)
print("=== median over core by release ===")
print(R.groupby("release")[["n","pooled","within","genes_with_pairs"]].median().round(3).to_string())
print(K.round(3).to_string())
print("=== change Jan 2020 -> Nov 2026 ===")
print("median d_pooled %.3f (CI>0 for %d/%d); median d_within %.3f (CI>0 %d, CI<0 %d)"%(C.d_pooled.median(),(C.dp_lo>0).sum(),len(C),C.d_within.median(),(C.dw_lo>0).sum(),(C.dw_hi<0).sum()))
print(C[C.tool.isin(["REVEL","AlphaMissense","ESM1b","GERP++","CADD","MetaRNN"])].round(3).to_string())
# labels of Jan 2020 that no longer hold: pooled REVEL AUROC held vs changed
G=SETS["2020-01"]
for lab,g in [("held",G[G.held]),("changed",G[~G.held])]:
    print(lab,len(g),int(g.y.sum()),{t:round(auc(g.y.values,g[t].values.astype(float)),3) for t in ["REVEL","AlphaMissense","GERP++"]})
print("DONE")
