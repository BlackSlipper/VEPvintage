#!/usr/bin/env python
"""(i) ACMG-code-level natural experiment: labels supported by PP3/BP4 vs code-citing labels without them.
(ii) Label-source stratification: one class from a given submitter vs all of the other class."""
import sys, gzip, re; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
D=pd.read_csv(f"{P}/results/variant_level_annotated.tsv.gz",sep="\t",dtype={"VariationID":str})
vids=set(D.VariationID.dropna())
CODE=re.compile(r"\b(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(?:[_ -]?(?:strong|moderate|supporting|very ?strong|str|mod|sup|vs|s|m|p))?\b",re.I)
RAW_SS=f"{P}/data/clinvar_tab/submission_summary_2026-11.txt.gz"
if __import__("os").path.exists(RAW_SS):
    rows=[]
    with gzip.open(RAW_SS,"rt",errors="replace") as fh:
        for l in fh:
            if l.startswith("#"): continue
            x=l.rstrip("\n").split("\t")
            if x[0] not in vids: continue
            codes=sorted(set(m.group(1).upper() for m in CODE.finditer(x[3])))
            if codes: rows.append((x[0],x[1],x[9],x[2],",".join(codes)))
    C=pd.DataFrame(rows,columns=["VariationID","ClinSig","Submitter","Date","codes"])
    cs=C.ClinSig.str.lower()
    C["scv_class"]=np.where(cs.str.contains("pathogenic"),"PLP",np.where(cs.str.contains("benign"),"BLB","OTHER"))
    C=C.merge(D[["VariationID","final"]],on="VariationID"); C=C[C.scv_class==C.final]
    C.to_csv(f"{P}/results/scv_acmg_codes.tsv.gz",sep="\t",index=False,compression="gzip")
else:  # release path: concordant code-citing SCVs shipped in data/ (copied by 00_use_release_data.py)
    C=pd.read_csv(f"{P}/results/scv_acmg_codes.tsv.gz",sep="\t",dtype={"VariationID":str})
    print("submission_summary_2026-11 not found; using results/scv_acmg_codes.tsv.gz")
print("code-citing concordant SCVs:",len(C),"variants:",C.VariationID.nunique())
print(C.groupby("final").codes.apply(lambda s: s.str.split(",").explode().value_counts().head(12)).to_string())
V=C.groupby("VariationID").agg(codes=("codes",lambda s: ",".join(sorted(set(",".join(s).split(","))))),code_subs=("Submitter",lambda s:"|".join(sorted(set(s)))))
V["has_PP3"]=V.codes.str.contains(r"\bPP3\b"); V["has_BP4"]=V.codes.str.contains(r"\bBP4\b")
D=D.merge(V,left_on="VariationID",right_index=True,how="left")
D["code_citing"]=D.codes.notna()
D.to_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",index=False,compression="gzip")
tools=[t for t in TOOLS if t in D.columns]; cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]
def logit(p): p=np.clip(p,1e-4,1-1e-4); return np.log(p/(1-p))
def panel(pos,neg,label):
    G=pd.concat([pos,neg]).dropna(subset=core)
    if G.y.sum()<30 or (1-G.y).sum()<30: return []
    a={t:auc(G.y.values,G[t].values.astype(float)) for t in core}; med=np.median([logit(v) for v in a.values()])
    return [dict(contrast=label,n_pos=int(G.y.sum()),n_neg=int((1-G.y).sum()),tool=t,auc=v,rel=logit(v)-med) for t,v in a.items()]
out=[]
for rare_only in [False,True]:
    B=D if not rare_only else D[D.rare]
    tag="rare" if rare_only else "all"
    P_all=B[B.y==1]; N_all=B[B.y==0]
    # benign side: BP4-supported vs code-citing without BP4
    out+=panel(P_all,N_all[N_all.code_citing&N_all.has_BP4],f"{tag}:BLB_with_BP4")
    out+=panel(P_all,N_all[N_all.code_citing&~N_all.has_BP4.fillna(False)],f"{tag}:BLB_codes_without_BP4")
    out+=panel(P_all[P_all.code_citing&P_all.has_PP3],N_all,f"{tag}:PLP_with_PP3")
    out+=panel(P_all[P_all.code_citing&~P_all.has_PP3.fillna(False)],N_all,f"{tag}:PLP_codes_without_PP3")
    # label source (benign side), 2023-2026
    for src in ["Ambry Genetics","Labcorp Genetics (formerly Invitae), Labcorp","CeGaT Center for Human Genetics Tuebingen","Blueprint Genetics","GeneDx","Illumina Laboratory Services, Illumina"]:
        for per in ["<=2019","2020-2022","2023-2026"]:
            Pp=P_all[P_all.period==per]; Nn=N_all[(N_all.period==per)&(N_all.conc_submitters.fillna("")==src)]
            out+=panel(Pp,Nn,f"{tag}:{per}:BLB_only_from:{src.split(chr(44))[0]}")
        Nn=N_all[(N_all.period=="2023-2026")&~N_all.conc_submitters.fillna("").str.contains("Ambry|Labcorp")]
    for per in ["<=2019","2020-2022","2023-2026"]:
        Pp=P_all[P_all.period==per]; Nn=N_all[(N_all.period==per)&~N_all.conc_submitters.fillna("").str.contains("Ambry Genetics|Labcorp Genetics")]
        out+=panel(Pp,Nn,f"{tag}:{per}:BLB_excluding_Ambry_Labcorp")
        out+=panel(Pp,N_all[N_all.period==per],f"{tag}:{per}:BLB_all")
R=pd.DataFrame(out); R.to_csv(f"{P}/results/evidence_source_contrasts.tsv",sep="\t",index=False)
pd.set_option("display.width",300); pd.set_option("display.max_columns",30)
key=["REVEL","BayesDel_noAF","AlphaMissense","ESM1b","CADD","phyloP100way","GERP++","VEST4","MetaRNN","ClinPred","PrimateAI","SIFT"]
W=R[R.tool.isin(key)].pivot_table(index="contrast",columns="tool",values="auc")
n=R.groupby("contrast")[["n_pos","n_neg"]].first()
print(W.join(n).round(3).to_string())
