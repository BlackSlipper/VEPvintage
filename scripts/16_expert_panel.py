#!/usr/bin/env python
"""Expert-panel (3-star) labels: vintage-resolved AUROC and PP3/BP4 usage."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
D=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",dtype={"VariationID":str})
E=D[D.stars>=3].copy()
print("expert-panel labels:",len(E),E.final.value_counts().to_dict(),"genes:",E.Gene.nunique())
print(E.groupby("period").final.value_counts().unstack().to_string())
print("code-citing share:",E.code_citing.mean().round(3),"PP3 among PLP code-citing:",E[(E.y==1)&E.code_citing].has_PP3.mean().round(3),"BP4 among BLB code-citing:",E[(E.y==0)&E.code_citing].has_BP4.mean().round(3))
tools=["REVEL","BayesDel_noAF","AlphaMissense","ESM1b","CADD","phyloP100way","GERP++","MetaRNN","VEST4","MutPred2"]
rows=[]
for per,g in E.groupby("period"):
    for t in tools:
        m=g[t].notna(); rows.append(dict(period=per,tool=t,n=int(m.sum()),n_plp=int(g.y[m].sum()),auc=auc(g.y.values,g[t].values.astype(float))))
for flag in [True,False]:
    g=E[(E.code_citing)&(np.where(E.y==1,E.has_PP3.fillna(False),E.has_BP4.fillna(False))==flag)]
    for t in tools:
        m=g[t].notna(); rows.append(dict(period=f"codes_comp={flag}",tool=t,n=int(m.sum()),n_plp=int(g.y[m].sum()),auc=auc(g.y.values,g[t].values.astype(float))))
R=pd.DataFrame(rows); R.to_csv(f"{P}/results/expert_panel_auc.tsv",sep="\t",index=False)
pd.set_option("display.width",200)
print(R.pivot(index="tool",columns="period",values="auc").round(3).to_string())
print(R.groupby("period")[["n","n_plp"]].max().to_string())
