#!/usr/bin/env python
"""Link SCV-level in-silico citations to variant labels; test whether tool-specific citation inflates tool-specific AUROC."""
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
D=load()
AF=pd.concat([pd.read_csv(f,sep="\t") for f in sorted(glob.glob(f"{P}/results/af_shard*.tsv.gz"))]).drop_duplicates("query")
for c in AF.columns[1:]: AF[c]=pd.to_numeric(AF[c],errors="coerce")
AF["maxaf"]=AF[["gnomad_ex_popmax","gnomad_ge_popmax","gnomad_ex_af","gnomad_ge_af","exac_af"]].max(axis=1).fillna(0)
D=D.merge(AF[["query","maxaf"]],on="query",how="left"); D["maxaf"]=D.maxaf.fillna(0)
S=pd.read_csv(f"{P}/results/scv_textmining.tsv.gz",sep="\t",dtype={"VariationID":str})
cs=S.ClinSig.str.lower()
S["scv_class"]=np.where(cs.str.contains("pathogenic")&~cs.str.contains("conflict"),"PLP",np.where(cs.str.contains("benign"),"BLB",np.where(cs.str.contains("uncertain"),"VUS","OTHER")))
S=S.merge(D[["VariationID","final"]],on="VariationID")
S["concordant"]=S.scv_class==S.final
C=S[S.concordant]
flags=["any_insilico","any_tool","REVEL","PolyPhen","SIFT","CADD","AlphaMissense","BayesDel","PP3","BP4","insilico_generic","has_text"]
V=C.groupby("VariationID")[flags].max().add_prefix("cite_")
V["n_conc_scv"]=C.groupby("VariationID").size()
V["conc_submitters"]=C.groupby("VariationID").Submitter.apply(lambda s:"|".join(sorted(set(s))))
D=D.merge(V,left_on="VariationID",right_index=True,how="left")
for f in flags: D["cite_"+f]=D["cite_"+f].fillna(0).astype(int)
D["period"]=np.where(D.legacy|(D.vyear<=2019),"<=2019",np.where(D.vyear<=2022,"2020-2022","2023-2026"))
D["rare"]=D.maxaf<1e-4
D.to_csv(f"{P}/results/variant_level_annotated.tsv.gz",sep="\t",index=False,compression="gzip")
pd.set_option("display.width",250)
print(D.groupby(["period","final"])[["cite_any_insilico","cite_REVEL","cite_PolyPhen","cite_SIFT","cite_BP4","cite_PP3","cite_has_text"]].mean().round(3).to_string())
print(D[D.rare].groupby(["period","final"])[["cite_any_insilico","cite_REVEL","cite_BP4","cite_PP3"]].mean().round(3).to_string())
# top submitters among concordant SCVs for rare BLB 2023-2026
x=D[(D.period=="2023-2026")&D.rare&(D.final=="BLB")].conc_submitters.fillna("").str.split("|").explode().value_counts().head(12)
print("rare BLB 2023-26 submitters:\n",x)
tools=[t for t in TOOLS if t in D.columns]; cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]
def logit(p): p=np.clip(p,1e-4,1-1e-4); return np.log(p/(1-p))
rows=[]
for per in ["<=2019","2020-2022","2023-2026"]:
    for rare in [True,False]:
        G0=D[(D.period==per)&(D.rare==rare)].dropna(subset=core)
        for flag in ["cite_any_insilico","cite_REVEL","cite_PolyPhen","cite_SIFT","cite_BP4","cite_PP3"]:
            for val in [0,1]:
                G=G0[G0[flag]==val]
                if G.y.sum()<20 or (1-G.y).sum()<20: continue
                aucs={t:auc(G.y.values,G[t].values.astype(float)) for t in core}
                med=np.median([logit(a) for a in aucs.values()])
                for t,a in aucs.items():
                    rows.append(dict(period=per,rare=rare,flag=flag,val=val,n=len(G),n_plp=int(G.y.sum()),tool=t,auc=a,rel=logit(a)-med))
R=pd.DataFrame(rows); R.to_csv(f"{P}/results/citation_stratified_auc.tsv",sep="\t",index=False)
for flag,tl in [("cite_REVEL",["REVEL","BayesDel_noAF","VEST4","ESM1b","phyloP100way","AlphaMissense"]),("cite_PolyPhen",["PolyPhen2_HVAR","PolyPhen2_HDIV","SIFT","REVEL","ESM1b"]),("cite_SIFT",["SIFT","SIFT4G","PolyPhen2_HVAR","REVEL","ESM1b"]),("cite_any_insilico",["REVEL","BayesDel_noAF","ESM1b","phyloP100way","GERP++","CADD"])]:
    r=R[(R.flag==flag)&R.tool.isin(tl)]
    print(flag); print(r.pivot_table(index=["period","rare","tool"],columns="val",values=["auc","rel"]).round(3).to_string())
    print(r.groupby(["period","rare","val"]).n.first().to_string())
