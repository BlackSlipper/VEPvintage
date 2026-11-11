#!/usr/bin/env python
"""Parse ClinVar variant_summary snapshots -> per-snapshot missense SNV labels keyed by AlleleID."""
import sys, re, gzip, os
import pandas as pd
P=__import__("os").environ.get("VEP_ROOT",__import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))
f=sys.argv[1]
snap=re.search(r"(\d{4}-\d{2}|current)", os.path.basename(f)).group(1)
out=f"{P}/data/snap_parsed/{snap}.tsv.gz"
if os.path.exists(out): sys.exit(0)
AA=r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)"
pat=re.compile(r"\(p\."+AA+r"(\d+)"+AA+r"\)$")
def cls(s):
    s=str(s).lower()
    if "conflicting" in s: return "CONF"
    has_p = ("pathogenic" in s) and ("conflict" not in s)
    has_b = "benign" in s
    if "uncertain" in s: return "VUS"
    if has_p and not has_b:
        if "risk" in s or "low penetrance" in s: return "OTHER"
        return "PLP"
    if has_b and not has_p: return "BLB"
    return "OTHER"
rows=[]
with gzip.open(f,"rt",errors="replace") as fh:
    hdr=fh.readline().lstrip("#").rstrip("\n").split("\t")
    ix={h:i for i,h in enumerate(hdr)}
    need=["AlleleID","Type","Name","GeneSymbol","ClinicalSignificance","ReviewStatus","Assembly","NumberSubmitters","LastEvaluated","Origin"]
    for line in fh:
        x=line.rstrip("\n").split("\t")
        if len(x)<len(hdr): continue
        if x[ix["Type"]]!="single nucleotide variant": continue
        asm=x[ix["Assembly"]]
        if asm not in ("GRCh38","GRCh37"): continue
        m=pat.search(x[ix["Name"]])
        if not m: continue
        if m.group(1)==m.group(3): continue
        rows.append((int(x[ix["AlleleID"]]),asm,x[ix["GeneSymbol"]],m.group(1)+m.group(2)+m.group(3),
                     x[ix["ClinicalSignificance"]],cls(x[ix["ClinicalSignificance"]]),x[ix["ReviewStatus"]],
                     x[ix["NumberSubmitters"]],x[ix["LastEvaluated"]],x[ix["Origin"]] if "Origin" in ix else ""))
df=pd.DataFrame(rows,columns=["AlleleID","Assembly","Gene","Protein","ClinSig","Class","Review","NSub","LastEval","Origin"])
# one row per AlleleID: prefer GRCh38
df=df.sort_values(["AlleleID","Assembly"],ascending=[True,False]).drop_duplicates("AlleleID")
df.insert(0,"Snapshot",snap)
os.makedirs(os.path.dirname(out),exist_ok=True)
df.to_csv(out,sep="\t",index=False,compression="gzip")
print(snap,len(df),df.Class.value_counts().to_dict())
