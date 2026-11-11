#!/usr/bin/env python
"""Gene-matched sensitivity: within-gene AUROC (pairs of P/LP and B/LB from the same gene), pooled across genes,
for BP4-cited vs other-code B/LB and PP3-cited vs other-code P/LP (rare variants)."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
D=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",dtype={"VariationID":str})
D=D[D.rare]
tools=["REVEL","BayesDel_noAF","AlphaMissense","ESM1b","CADD","phyloP100way","GERP++","PrimateAI","MetaRNN","VEST4","SIFT","PolyPhen2_HVAR"]
def within_gene_auc(pos,neg,t):
    num=0.0; den=0.0
    for g,pn in neg.groupby("Gene"):
        pp=pos[pos.Gene==g]
        if len(pp)==0: continue
        a=pp[t].dropna().values; b=pn[t].dropna().values
        if len(a)==0 or len(b)==0: continue
        s=np.concatenate([a,b]); r=pd.Series(s).rank().values
        u=r[:len(a)].sum()-len(a)*(len(a)+1)/2
        num+=u; den+=len(a)*len(b)
    return num/den if den>0 else np.nan, den
rows=[]
Pall=D[D.y==1]; Nall=D[D.y==0]
groups={"BLB_with_BP4":(Pall,Nall[Nall.code_citing&Nall.has_BP4.fillna(False)]),
        "BLB_codes_without_BP4":(Pall,Nall[Nall.code_citing&~Nall.has_BP4.fillna(False)]),
        "PLP_with_PP3":(Pall[Pall.code_citing&Pall.has_PP3.fillna(False)],Nall),
        "PLP_codes_without_PP3":(Pall[Pall.code_citing&~Pall.has_PP3.fillna(False)],Nall)}
# restrict BLB contrasts to genes present in both BLB groups; PLP contrasts to genes present in both PLP groups
gB=set(groups["BLB_with_BP4"][1].Gene)&set(groups["BLB_codes_without_BP4"][1].Gene)
gP=set(groups["PLP_with_PP3"][0].Gene)&set(groups["PLP_codes_without_PP3"][0].Gene)
for k,(pos,neg) in groups.items():
    gs=gB if k.startswith("BLB") else gP
    pos=pos[pos.Gene.isin(gs)]; neg=neg[neg.Gene.isin(gs)]
    for t in tools:
        a,npairs=within_gene_auc(pos,neg,t)
        rows.append(dict(contrast=k,tool=t,within_gene_auc=a,pairs=int(npairs),n_genes=len(gs),n_pos=len(pos),n_neg=len(neg)))
R=pd.DataFrame(rows); R.to_csv(f"{P}/results/within_gene_code_contrasts.tsv",sep="\t",index=False)
W=R.pivot(index="tool",columns="contrast",values="within_gene_auc")
W["BP4_infl"]=W.BLB_with_BP4-W.BLB_codes_without_BP4; W["PP3_infl"]=W.PLP_with_PP3-W.PLP_codes_without_PP3
pd.set_option("display.width",200); print(W.round(3).to_string()); print(R.groupby("contrast")[["n_genes","n_pos","n_neg"]].first())
