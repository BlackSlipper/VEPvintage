#!/usr/bin/env python
"""Further analyses: (i) composition of filtered tier by submitter, exclusion of Submitter-B-only labels;
(ii) within-gene vintage jump and Submitter-A contrast; (iii) rationale availability by period;
(iv) bootstrap rank intervals; (v) BH-FDR for relative slopes; (vi) submitter-level calibration LRs; (vii) resolved-VUS by provenance."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
from statsmodels.stats.multitest import multipletests
D=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",dtype={"VariationID":str})
T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","clean"])
D=D.merge(T,on="hgvs38",how="left")
SA="Ambry Genetics"; SB="Labcorp Genetics (formerly Invitae), Labcorp"
out={}
# (i)
R=D[D.rare&(D.period=="2023-2026")]
c=R[R.clean.fillna(False)]
src=np.where(c.conc_submitters==SB,"B_only",np.where(c.conc_submitters.fillna("").str.contains(SA),"includes_A","other"))
print("(i) clean rare 2023-26 by source & class:\n",pd.crosstab(src,c.final))
key=["REVEL","AlphaMissense","ESM1b","CADD","GERP++","MetaRNN","MutPred2","phyloP100way"]
rows=[]
for name,G in {"clean":c,"clean_excl_Bonly":c[src!="B_only"],"clean_Bonly_BLB+allPLP":pd.concat([c[c.y==1],c[(c.y==0)&(src=="B_only")]])}.items():
    for t in key:
        g=G[G[t].notna()]; rows.append(dict(set=name,tool=t,n=len(g),n_plp=int(g.y.sum()),auc=auc(g.y.values,g[t].values.astype(float))))
X=pd.DataFrame(rows); print(X.pivot(index="tool",columns="set",values="auc").round(3).to_string()); print(X.groupby("set")[["n","n_plp"]].max())
X.to_csv(f"{P}/results/review_clean_by_source.tsv",sep="\t",index=False)
# (ii) within-gene AUROC: pooled over genes of P/LP-B/LB pairs within gene
def wg(G,t):
    num=den=0.0
    for g,d in G.groupby("Gene"):
        a=d[d.y==1][t].dropna().values; b=d[d.y==0][t].dropna().values
        if len(a)==0 or len(b)==0: continue
        s=np.concatenate([a,b]); r=pd.Series(s).rank().values
        num+=r[:len(a)].sum()-len(a)*(len(a)+1)/2; den+=len(a)*len(b)
    return num/den if den else np.nan
rows=[]
Rr=D[D.rare]
for per in ["<=2019","2020-2022","2023-2026"]:
    G=Rr[Rr.period==per]
    for t in key: rows.append(dict(analysis="within_gene_vintage",period=per,tool=t,auc=wg(G,t)))
    G2=G[G.y==1]
    Ga=pd.concat([G2,G[(G.y==0)&(G.conc_submitters==SA)]]); Go=pd.concat([G2,G[(G.y==0)&~G.conc_submitters.fillna("").str.contains("Ambry Genetics|Labcorp Genetics")]])
    for t in key:
        rows.append(dict(analysis="within_gene_BLB_A_only",period=per,tool=t,auc=wg(Ga,t)))
        rows.append(dict(analysis="within_gene_BLB_excl_AB",period=per,tool=t,auc=wg(Go,t)))
W=pd.DataFrame(rows); W.to_csv(f"{P}/results/review_within_gene.tsv",sep="\t",index=False)
print("(ii)\n",W.pivot_table(index=["analysis","tool"],columns="period",values="auc").round(3).to_string())
# (iii) rationale availability
D["any_text"]=D.cite_has_text.fillna(0)>0
print("(iii) share with concordant free-text rationale / code-citing, by period & class:\n",D.groupby(["period","final"])[["any_text","code_citing"]].mean().round(3).to_string())
# (v) BH-FDR on relative slopes
O=pd.read_csv(f"{P}/results/relative_slopes.tsv",sep="\t"); o=O[O.subset=="rare_common"].copy()
o["q_BH"]=multipletests(o.p,method="fdr_bh")[1]; o.to_csv(f"{P}/results/relative_slopes_rare_fdr.tsv",sep="\t",index=False)
print("(v)\n",o.sort_values("p")[["tool","rel_slope","p","q_BH"]].head(8).round(4).to_string())
# (vii) resolved VUS (T=2024-01) by provenance
def state_at(h,Tt):
    st="ABSENT"
    for tok in h.split(";"):
        s,cc=tok.split(":")
        if s=="current" or s>Tt: break
        st=cc
    return st
H=pd.read_csv(f"{P}/results/variant_history.tsv.gz",sep="\t",usecols=["hgvs38","history"]).drop_duplicates("hgvs38")
D2=D.merge(H,on="hgvs38",how="left",suffixes=("","_h"))
hc="history" if "history" in D2.columns and D2["history"].notna().any() else "history_h"
D2["stT"]=D2[hc].fillna("").map(lambda h: state_at(h,"2024-01") if h else "ABSENT")
V=D2[D2.stT.isin(["VUS","CONF"])]
rows=[]
for name,G in {"resolved_all":V,"resolved_clean":V[V.clean.fillna(False)],"resolved_nonclean":V[~V.clean.fillna(False)]}.items():
    for t in key: rows.append(dict(set=name,tool=t,n=int(G[t].notna().sum()),n_plp=int(G.y.sum()),auc=auc(G.y.values,G[t].values.astype(float))))
Y=pd.DataFrame(rows); Y.to_csv(f"{P}/results/review_resolved_by_provenance.tsv",sep="\t",index=False)
print("(vii)\n",Y.pivot(index="tool",columns="set",values="auc").round(3).to_string()); print(Y.groupby("set")[["n","n_plp"]].max())
