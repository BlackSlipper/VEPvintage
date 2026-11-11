import pandas as pd
import os; P = __import__("os").environ.get("VEP_ROOT", __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))
A=pd.read_csv(f"{P}/results/r2/R_audit_sample.tsv",sep="\t")
lab={}
for i in list(range(20,40))+list(range(60,80))+list(range(100,120))+list(range(140,160)): lab[i]="applied in support of the classification"
for i in [24,28,104,107,118,119,141,143,149]: lab[i]="inconclusive, conflicting or explicitly not applied"
for i in [31,32,156]: lab[i]="mentioned as pointing against the classification"
for i in [100,101]: lab[i]="generic template listing possible criteria"
for i in list(range(0,20))+list(range(40,60))+list(range(80,100))+list(range(120,140)): lab[i]="no computational missense evidence"
for i in [123,134]: lab[i]="computational evidence present, missed by the patterns"
for i in [41,50]: lab[i]="in-house structural modelling, missed by the patterns (uncertain)"
for i in [17,135]: lab[i]="splicing predictor only (not counted)"
for i in [52,59]: lab[i]="laboratory model of experimental data (not a predictor)"
for i in [81,83]: lab[i]="automated frequency score (not a predictor)"
A["annotation"]=[lab[i] for i in range(len(A))]
A["annotated_by"]="draft annotation, to be checked by the authors"
A.to_csv(f"{P}/results/r2/R_audit_annotated.tsv",sep="\t",index=False)
A["flag"]=A.stratum.str.split("|").str[2]; A["per"]=A.stratum.str.split("|").str[0]; A["cls"]=A.stratum.str.split("|").str[1]
print(pd.crosstab([A.flag,A.annotation],[A.per,A.cls]).to_string())
print("submitters", A.Submitter.nunique())
