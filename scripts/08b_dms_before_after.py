#!/usr/bin/env python3
"""DMS assays against ClinVar labels assigned up to and after the assay's publication year (Supplementary Methods S9,
Supplementary Table S11, Fig. S1). Reads results/dms_clinvar_merged.tsv.gz from 08_dms.py and writes
results/dms_before_after_pub.tsv. A label is 'before' when the year of its first P/LP or B/LB class (vintage_first) is no
later than the publication year recorded for the assay in ProteinGym. An assay period is listed when it holds at least
8 labels of each class, and AUROCs are computed when it holds at least 10 of each (DMS scores are oriented so that lower
fitness means pathogenic)."""
import os, sys
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
P = os.environ.get("VEP_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{P}/results/dms_before_after_pub.tsv"
M = pd.read_csv(f"{P}/results/dms_clinvar_merged.tsv.gz", sep="\t")
M["yr"] = pd.to_numeric(M.vintage_first.astype(str).str[:4], errors="coerce")
M["set"] = np.where(M.yr <= M.dms_year, "before_pub", "after_pub")
M["DMS"] = -M.DMS_score
rows = []
for (a, s), d in M.groupby(["DMS_id", "set"]):
    n1 = int(d.y.sum()); n0 = len(d) - n1
    if min(n1, n0) < 8: continue
    r = dict(assay=a, gene=d.Gene.iloc[0] if "Gene" in d else a.split("_")[0], pub=int(d.dms_year.iloc[0]), set=s, n=len(d), n_plp=n1)
    for c, lab in [("DMS", "auc_DMS"), ("REVEL", "auc_REVEL"), ("AlphaMissense", "auc_AlphaMissense"), ("ESM1b", "auc_ESM1b")]:
        x = d[[c, "y"]].dropna()
        r[lab] = roc_auc_score(x.y, x[c]) if min(n1, n0) >= 10 and x.y.nunique() == 2 else np.nan
    rows.append(r)
R = pd.DataFrame(rows)
R["o"] = R.set.map({"before_pub": 0, "after_pub": 1}); R = R.sort_values(["assay", "o"]).drop(columns="o")
R.to_csv(OUT, sep="\t", index=False); print(R.to_string())
