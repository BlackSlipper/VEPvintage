#!/usr/bin/env python
"""Sensitivity: 'rare' defined by gnomAD v4.1 joint AF (<1e-4; absent counted as rare)."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
G4=pd.read_csv(f"{P}/results/gnomad4_af_merged.tsv",sep="\t",header=None,names=["vid","ac","an","faf"])
G4["af4"]=pd.to_numeric(G4.ac,errors="coerce")/pd.to_numeric(G4.an,errors="coerce").replace(0,np.nan); G4["faf"]=pd.to_numeric(G4.faf,errors="coerce")
D=load(); T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","clean","rare","period"])
D=D.merge(T,on="hgvs38")
D["vid"]=D.hgvs38.str.replace("chr","",regex=False).str.replace(":g.","-",regex=False).str.replace(r"(\d+)([ACGT])>([ACGT])$",r"\1-\2-\3",regex=True)
D=D.merge(G4[["vid","af4","faf"]],on="vid",how="left")
D["af4"]=D.af4.fillna(0); D["faf"]=D.faf.fillna(0)
D["rare4"]=(D.af4<1e-4)&(D.faf<1e-4)
print("rare v2:",D.rare.sum(),"rare v4:",D.rare4.sum(),"both:",(D.rare&D.rare4).sum(),"v2-only:",(D.rare&~D.rare4).sum(),"v4-only:",(~D.rare&D.rare4).sum())
print(pd.crosstab([D.period,D.final],[D.rare,D.rare4]).to_string())
key=["REVEL","AlphaMissense","ESM1b","CADD","GERP++","phyloP100way","MetaRNN","MutPred2"]
tools=[t for t in TOOLS if t in D.columns]; cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.9]
C=D.dropna(subset=core)
rows=[]
def wg(G,t):
    num=den=0.0
    for g,d in G.groupby("Gene"):
        a=d[d.y==1][t].dropna().values; b=d[d.y==0][t].dropna().values
        if len(a)==0 or len(b)==0: continue
        s=np.concatenate([a,b]); r=pd.Series(s).rank().values; num+=r[:len(a)].sum()-len(a)*(len(a)+1)/2; den+=len(a)*len(b)
    return num/den if den else np.nan
for per in ["<=2019","2020-2022","2023-2026"]:
    g=C[C.rare4&(C.period==per)]
    for t in key:
        rows.append(dict(analysis="pooled_rare4_common",period=per,tool=t,n=len(g),auc=auc(g.y.values,g[t].values.astype(float))))
        rows.append(dict(analysis="withingene_rare4",period=per,tool=t,n=len(g),auc=wg(g,t)))
R4=D[D.rare4&(D.period=="2023-2026")]
for name,G in {"full_rare4_2023_26":R4,"clean_rare4_2023_26":R4[R4.clean]}.items():
    for t in key: rows.append(dict(analysis=name,period="2023-2026",tool=t,n=int(G[t].notna().sum()),auc=auc(G.y.values,G[t].values.astype(float))))
O=pd.DataFrame(rows); O.to_csv(f"{P}/results/sensitivity_gnomad4.tsv",sep="\t",index=False)
pd.set_option("display.width",220)
print(O.pivot_table(index=["analysis","tool"],columns="period",values="auc").round(3).to_string())
print(O.groupby(["analysis","period"]).n.max().to_string())
Th=pd.read_csv(f"{P}/meta/clingen_thresholds.tsv",sep="\t"); rows=[]
G=D[D.rare4].copy(); G["per"]=np.where(G.legacy|(G.vyear<=2019),"<=2019",np.where(G.vyear<=2022,"2020-2022","2023-2026"))
for (tool,rc),tt in Th.groupby(["tool","rawcol"]):
    for per,g in G.groupby("per"):
        g=g[g[rc].notna()]; n1=(g.y==1).sum(); n0=(g.y==0).sum()
        r=tt[tt.evidence=="PP3_Supporting"]
        if len(r)==0: continue
        r=r.iloc[0]; inb=(g[rc]>=float(r.lo))&(g[rc]<float(r.hi)); a=(inb&(g.y==1)).sum(); b=(inb&(g.y==0)).sum()
        rows.append(dict(tool=tool,period=per,LR=((a+.5)/(n1+.5))/((b+.5)/(n0+.5))))
L=pd.DataFrame(rows); L.to_csv(f"{P}/results/sensitivity_gnomad4_LR.tsv",sep="\t",index=False)
W=L.pivot(index="tool",columns="period",values="LR"); W["ratio"]=W["2023-2026"]/W["<=2019"]; print(W.round(2).to_string())
