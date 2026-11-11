#!/usr/bin/env python
"""Prospective (time-split) benchmark: variants resolved after reference snapshot T vs those already labelled at T."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
D=load()
tools=[t for t in TOOLS if t in D.columns]
def state_at(hist,T):
    st="ABSENT"
    for tok in hist.split(";"):
        s,c=tok.split(":")
        if s=="current" or s>T: break
        st=c
    return st
out=[]
for T in ["2020-01","2022-01","2023-01","2024-01","2025-01"]:
    st=D.history.map(lambda h: state_at(h,T))
    known=D[st==D.final]          # already carrying final label at T
    vus=D[st.isin(["VUS","CONF"])]  # unresolved at T, resolved later
    new=D[st=="ABSENT"]             # first submitted after T
    for name,G in [("known_at_T",known),("VUS_resolved_after_T",vus),("new_after_T",new)]:
        for t in tools:
            s=G[t].values.astype(float); y=G.y.values
            m=~np.isnan(s)
            a=auc(y,s)
            out.append(dict(T=T,set=name,tool=t,n=int(m.sum()),n_plp=int(y[m].sum()),auc=a))
O=pd.DataFrame(out)
O.to_csv(f"{P}/results/prospective_auc.tsv",sep="\t",index=False)
W=O[O["T"]=="2024-01"].pivot(index="tool",columns="set",values="auc")
W["drop_vus"]=W["VUS_resolved_after_T"]-W["known_at_T"]
W=W.join(META.set_index("tool")[["category"]]).sort_values("VUS_resolved_after_T",ascending=False)
print(O[(O.tool=="REVEL")].to_string())
print(W.round(3).to_string())
