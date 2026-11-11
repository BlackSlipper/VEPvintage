#!/usr/bin/env python
"""(a) Selection null: simulate PP3/BP4 assignment from published REVEL thresholds among code-citing labels and compare simulated vs observed AUROC gaps.
(b) Residual false-negative audit of the text-mining regex on 40 random unflagged non-template SCVs with free text."""
import sys, gzip, re, random; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
D=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",dtype={"VariationID":str})
D=D[D.rare]
tools=[t for t in TOOLS if t in D.columns]; cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.9]
Pc=D[(D.y==1)&D.code_citing]; Nc=D[(D.y==0)&D.code_citing]; Pall=D[D.y==1]; Nall=D[D.y==0]
rows=[]
for side,grp,other,obs_flag,sim_flag in [("PP3",Pc,Nall,Pc.has_PP3.fillna(False),Pc.REVEL_raw>=0.644),("BP4",Nc,Pall,Nc.has_BP4.fillna(False),Nc.REVEL_raw<=0.290)]:
    for kind,flag in [("observed",obs_flag),("simulated_REVEL",sim_flag)]:
        for t in core:
            a1=pd.concat([grp[flag.values],other]); a0=pd.concat([grp[~flag.values],other])
            rows.append(dict(side=side,kind=kind,tool=t,gap=auc(a1.y.values,a1[t].values.astype(float))-auc(a0.y.values,a0[t].values.astype(float)),frac_flag=float(flag.mean())))
R=pd.DataFrame(rows); W=R.pivot_table(index=["side","tool"],columns="kind",values="gap")
R.to_csv(f"{P}/results/selection_null.tsv",sep="\t",index=False)
print(R.groupby(["side","kind"]).frac_flag.first())
for s in ["PP3","BP4"]:
    w=W.loc[s]; print(s,"median observed gap %.3f, simulated %.3f; corr %.2f"%(w.observed.median(),w.simulated_REVEL.median(),w.corr().iloc[0,1]))
    print(w.loc[["REVEL","ESM1b","AlphaMissense","GERP++","PrimateAI","MetaRNN"]].round(3).to_string())
# (b) audit
S=pd.read_csv(f"{P}/results/scv_textmining.tsv.gz",sep="\t",usecols=["SCV","Submitter","any_insilico","has_text","conservation_text"])
neg=S[(S.Submitter!="Ambry Genetics")&(S.has_text==1)&(S.any_insilico==0)].sample(40,random_state=11)
want=set(neg.SCV); out={}
K=re.compile(r"predict|silico|comput|algorithm|model|pp3|bp4|polyphen|sift|revel|cadd|align|conserv|score",re.I)
with gzip.open(f"{P}/data/clinvar_tab/submission_summary_2026-11.txt.gz","rt",errors="replace") as fh:
    for l in fh:
        if l.startswith("#"): continue
        x=l.split("\t")
        if x[10] in want: out[x[10]]=x[3]
        if len(out)==len(want): break
for k,t in out.items():
    m=K.search(t); print("AUDIT|",(t[max(0,m.start()-90):m.start()+110] if m else "(no keyword) "+t[:80]).replace("\n"," "))
