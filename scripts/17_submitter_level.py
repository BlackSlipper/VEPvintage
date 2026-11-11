#!/usr/bin/env python
"""Submitter-level provenance filter (avoids conditioning on per-variant evidence that depends on the score):
keep labels supported by >=1 concordant SCV from submitters whose free-text records rarely mention computational evidence."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
from scipy.stats import kendalltau
S=pd.read_csv(f"{P}/results/scv_textmining.tsv.gz",sep="\t",dtype={"VariationID":str})
C=pd.read_csv(f"{P}/results/scv_acmg_codes.tsv.gz",sep="\t",dtype={"VariationID":str})
st=S[S.has_text==1].groupby("Submitter").agg(n_text=("SCV","size"),frac_comp=("any_insilico","mean"))
cc=C.assign(comp=C.codes.str.contains(r"\bPP3\b|\bBP4\b")).groupby("Submitter").comp.mean().rename("frac_code_comp")
st=st.join(cc).fillna({"frac_code_comp":0})
noncomp=set(st[(st.n_text>=20)&(st.frac_comp<0.05)&(st.frac_code_comp<0.05)].index)-{"Ambry Genetics"}
print("non-computational submitters:",len(noncomp))
D=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",dtype={"VariationID":str})
cs=S.ClinSig.str.lower()
S["scv_class"]=np.where(cs.str.contains("pathogenic"),"PLP",np.where(cs.str.contains("benign"),"BLB","OTHER"))
S=S.merge(D[["VariationID","final"]],on="VariationID"); S=S[S.scv_class==S.final]
keep=set(S[S.Submitter.isin(noncomp)].VariationID)
D["subclean"]=D.VariationID.isin(keep)
print("submitter-level tier:",D.subclean.sum(),D[D.subclean].final.value_counts().to_dict())
print(D[D.subclean].groupby("period").final.value_counts().unstack().to_string())
top=S[S.Submitter.isin(noncomp)].Submitter.value_counts().head(10); print(top)
tools=[t for t in TOOLS if t in D.columns]
rows=[]
for name,G in {"full_rare_2023_26":D[D.rare&(D.period=="2023-2026")],"subclean_rare_2023_26":D[D.subclean&D.rare&(D.period=="2023-2026")],
               "subclean_rare_all":D[D.subclean&D.rare],"full_rare_all":D[D.rare]}.items():
    for t in tools:
        g=G[G[t].notna()]; rows.append(dict(set=name,tool=t,n=len(g),n_plp=int(g.y.sum()),auc=auc(g.y.values,g[t].values.astype(float))))
R=pd.DataFrame(rows); R.to_csv(f"{P}/results/submitter_level_tier.tsv",sep="\t",index=False)
W=R.pivot(index="tool",columns="set",values="auc")
B=pd.read_csv(f"{P}/results/benchmark_clean_rare_recent.tsv",sep="\t").pivot(index="tool",columns="set",values="auc")
W=W.join(B[["clean_rare_2023_26"]])
for c in W.columns: W["r_"+c]=W[c].rank(ascending=False)
pd.set_option("display.width",250)
print(W.sort_values("subclean_rare_2023_26",ascending=False).round(3).to_string())
print(R.groupby("set")[["n","n_plp"]].max().to_string())
print("tau subclean vs clean (rare 23-26):",round(kendalltau(W.subclean_rare_2023_26,W.clean_rare_2023_26)[0],3),
      "subclean vs full:",round(kendalltau(W.subclean_rare_2023_26,W.full_rare_2023_26)[0],3))
