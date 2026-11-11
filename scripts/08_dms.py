#!/usr/bin/env python
"""Map ProteinGym human DMS to ClinVar evaluation variants; DMS as label-independent reference."""
import sys, os, re, json, requests; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
from scipy.stats import spearmanr
DM=pd.read_csv(f"{P}/data/dms/DMS_substitutions.csv")
DM=DM[DM.taxon=="Human"]
mp_f=f"{P}/meta/uniprot_entry2gene.json"
if os.path.exists(mp_f): mp=json.load(open(mp_f))
else:
    mp={}
    for e in sorted(DM.UniProt_ID.unique()):
        r=requests.get(f"https://rest.uniprot.org/uniprotkb/search?query=id:{e}&fields=accession,gene_primary&format=tsv",timeout=60)
        L=r.text.strip().split("\n")
        if len(L)>1:
            acc,g=L[1].split("\t")[:2]; mp[e]={"acc":acc,"gene":g}
    json.dump(mp,open(mp_f,"w"),indent=1)
AA3={"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
D=load()
m=D.Protein.str.extract(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")
D["mut1"]=m[0].map(AA3)+m[1]+m[2].map(AA3)
tools=[t for t in TOOLS if t in D.columns]
rows=[]; merged=[]
for _,d in DM.iterrows():
    g=mp.get(d.UniProt_ID,{}).get("gene")
    if not g: continue
    f=f"{P}/data/dms/DMS_ProteinGym_substitutions/{d.DMS_filename}"
    if not os.path.exists(f): continue
    X=pd.read_csv(f)
    X=X[~X.mutant.str.contains(":")]
    seq=d.target_seq
    # check numbering consistency
    ok=X.mutant.map(lambda s: int(s[1:-1])<=len(seq) and seq[int(s[1:-1])-1]==s[0]).mean()
    G=D[D.Gene==g]
    if len(G)==0: continue
    J=G.merge(X[["mutant","DMS_score","DMS_score_bin"]],left_on="mut1",right_on="mutant")
    # wt check between ClinVar protein numbering and DMS target seq
    wt_ok=J.mut1.map(lambda s: int(s[1:-1])<=len(seq) and seq[int(s[1:-1])-1]==s[0]).mean() if len(J) else np.nan
    rows.append(dict(DMS_id=d.DMS_id,gene=g,year=d.year,assay=d.selection_assay,n_dms=len(X),seq_ok=ok,n_clinvar=len(J),n_plp=int(J.y.sum()) if len(J) else 0,wt_ok=wt_ok))
    if len(J): J["DMS_id"]=d.DMS_id; J["dms_year"]=d.year; merged.append(J)
S=pd.DataFrame(rows).sort_values("n_clinvar",ascending=False)
S.to_csv(f"{P}/results/dms_overlap_summary.tsv",sep="\t",index=False)
print(S.head(40).to_string())
M=pd.concat(merged)
M.to_csv(f"{P}/results/dms_clinvar_merged.tsv.gz",sep="\t",index=False,compression="gzip")
