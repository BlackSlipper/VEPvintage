import pandas as pd, numpy as np, glob
P=__import__("os").environ.get("VEP_ROOT",__import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))
META=pd.read_csv(f"{P}/meta/predictors.tsv",sep="\t")
TOOLS=META.tool.tolist()
def load(eval_only=True, min_stars=1):
    V=pd.read_csv(f"{P}/results/variant_history.tsv.gz",sep="\t",dtype={"VariationID":str,"Chromosome":str})
    import os
    fs=sorted(glob.glob(os.environ.get("DBNSFP_GLOB",f"{P}/results/dbnsfp_shard*.tsv.gz")))
    S=pd.concat([pd.read_csv(f,sep="\t",dtype={"query":str,"uniprot":str,"aapos":str}) for f in fs]).drop_duplicates("query")
    D=V.merge(S,left_on="hgvs38",right_on="query",how="left")
    if eval_only:
        D=D[D.final.isin(["PLP","BLB"])&(D.stars>=min_stars)&D.OriginSimple.fillna("").str.contains("germline")&(D.found==1)]
    D["y"]=(D.final=="PLP").astype(int)
    D["vyear"]=D.vintage_first.str[:4].astype(int)
    D["vyear_run"]=D.vintage_run.str[:4].astype(int)
    D["vq"]=pd.PeriodIndex(D.vintage_first,freq="M").asfreq("Q").astype(str)
    D["legacy"]=(D.vintage_first=="2015-02")  # already labelled in first snapshot
    return D
def auc(y,s):
    """Mann-Whitney AUROC with ties; NaNs dropped."""
    m=~np.isnan(s); y=y[m]; s=s[m]
    n1=y.sum(); n0=len(y)-n1
    if n1<10 or n0<10: return np.nan
    r=pd.Series(s).rank().values
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
def boot_auc(y,s,B=200,seed=0):
    rng=np.random.default_rng(seed); n=len(y); out=[]
    for _ in range(B):
        i=rng.integers(0,n,n); out.append(auc(y[i],s[i]))
    return np.nanpercentile(out,[2.5,97.5])
