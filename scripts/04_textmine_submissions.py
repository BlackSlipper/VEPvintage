#!/usr/bin/env python
"""Text-mine ClinVar SCV interpretation text for in-silico evidence mentions (missense variants in our history table)."""
import gzip, re, pandas as pd, numpy as np
P=__import__("os").environ.get("VEP_ROOT",__import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))
V=pd.read_csv(f"{P}/results/variant_history.tsv.gz",sep="\t",usecols=["AlleleID","VariationID","final"],dtype={"VariationID":str})
vids=set(V.VariationID.dropna())
TOOLS={
 "REVEL":r"\brevel\b","CADD":r"\bcadd\b","PolyPhen":r"poly-?phen","SIFT":r"\bsift\b","MutationTaster":r"mutation\s?taster",
 "AlphaMissense":r"alpha\s?missense","BayesDel":r"bayes\s?del","PrimateAI":r"primate\s?ai","Align-GVGD":r"align[- ]?gvgd",
 "VEST":r"\bvest\d?\b","MutPred":r"\bmutpred","ESM":r"\besm-?1[bv]?\b","EVE":r"\beve\b(?! of)","MetaRNN":r"\bmetarnn\b","ClinPred":r"\bclinpred\b",
 "MetaSVM/LR":r"\bmeta(svm|lr)\b","FATHMM":r"\bfathmm","PROVEAN":r"\bprovean\b","GERP":r"\bgerp","phyloP":r"\bphylop\b","MutationAssessor":r"mutation\s?assessor",
 "VARITY":r"\bvarity\b","DANN":r"\bdann\b","M-CAP":r"\bm-?cap\b","SpliceAI":r"\bspliceai\b",
}
GEN={
 "PP3":r"\bpp3(?![0-9])","BP4":r"\bbp4(?![0-9])",
 "conservation_text":r"conservation analys|conserved across (species|vertebrates|mammals)|evolutionar(il)?y conserved|highly conserved|conserved (residue|position|amino acid)",
 "insilico_generic":r"in[- ]?silico|computational (prediction|tool|evidence|analys|model|algorithm)|prediction (tool|algorithm|program)s?|predictive (algorithm|tool|model)|algorithms? (developed|designed) to predict|advanced modeling of protein|evidence modeling of protein|modeling of protein sequence|computer[- ]based|bioinformatic(s)? (tool|prediction|analys)|(variant |missense )?prediction programs?|predictive programs?|predicted to be (deleterious|damaging|tolerated|benign|neutral) by|silico (tool|predict|analys|model)",
}
REG={k:re.compile(v,re.I) for k,v in {**TOOLS,**GEN}.items()}
rows=[]
with gzip.open(f"{P}/data/clinvar_tab/submission_summary_2026-11.txt.gz","rt",errors="replace") as fh:
    for line in fh:
        if line.startswith("##"): continue
        if line.startswith("#"):
            hdr=line[1:].rstrip("\n").split("\t"); ix={h:i for i,h in enumerate(hdr)}; continue
        x=line.rstrip("\n").split("\t")
        if x[0] not in vids: continue
        desc=x[ix["Description"]] if "Description" in ix else ""
        expl=x[ix["ExplanationOfInterpretation"]] if "ExplanationOfInterpretation" in ix and len(x)>ix["ExplanationOfInterpretation"] else ""
        txt=(desc+" "+expl).strip()
        r={"VariationID":x[0],"SCV":x[ix["SCV"]],"Submitter":x[ix["Submitter"]],"ClinSig":x[ix["ClinicalSignificance"]],
           "DateLastEvaluated":x[ix["DateLastEvaluated"]],"ReviewStatus":x[ix["ReviewStatus"]],"CollectionMethod":x[ix["CollectionMethod"]],
           "has_text":int(len(txt)>20 and txt not in ("-","not provided")),"text_len":len(txt)}
        for k,rg in REG.items(): r[k]=int(bool(rg.search(txt)))
        rows.append(r)
S=pd.DataFrame(rows)
tool_cols=[k for k in TOOLS if k!="SpliceAI"]
S["any_tool"]=S[tool_cols].max(axis=1)
S["any_insilico"]=S[["any_tool","insilico_generic","PP3","BP4"]].max(axis=1)
S["year"]=pd.to_datetime(S.DateLastEvaluated,errors="coerce").dt.year
S.to_csv(f"{P}/results/scv_textmining.tsv.gz",sep="\t",index=False,compression="gzip")
print("SCVs:",len(S),"with text:",S.has_text.sum())
print(S.groupby("year")[["has_text","any_insilico","any_tool","REVEL","PolyPhen","SIFT","CADD","AlphaMissense","BayesDel","PP3","BP4"]].mean().round(3).tail(16))
