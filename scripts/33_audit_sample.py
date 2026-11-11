#!/usr/bin/env python
"""Stratified audit sample of submission-record text (extends the audit of Supplementary Methods S4).
Strata: year last evaluated (<=2021, 2023-2026) x class (P/LP, B/LB) x text-mining flag (flagged as computational evidence,
i.e. any_insilico incl. PP3/BP4; not flagged but with free text). Concordant records of evaluation variants, Ambry excluded
(template text), at most two records per submitter and stratum, 20 records per stratum (seed 33).
Writes results/r2/R_audit_sample.tsv (record, stratum, text: a window around the first match for flagged records, the whole text otherwise) for annotation: for flagged records whether computational
evidence is applied in support of the classification, mentioned without being applied, negated (not met / conflicting), or a
false match; for unflagged records whether computational evidence is used but missed by the patterns."""
import gzip, re, os, pandas as pd, numpy as np
import os; P = __import__("os").environ.get("VEP_ROOT", __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))
src = open(f"{P}/scripts/04_textmine_submissions.py").read()
exec(src.split("REG=")[0].split("V=pd.read_csv")[0]); exec("TOOLS=" + src.split("TOOLS=")[1].split("REG=")[0])
ANY = re.compile("|".join(f"(?:{v})" for k, v in {**TOOLS, **GEN}.items() if k != "conservation_text"), re.I)
S = pd.read_csv(f"{P}/results/scv_textmining.tsv.gz", sep="\t", usecols=["VariationID", "SCV", "Submitter", "ClinSig", "has_text", "any_insilico", "year"])
T = pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz", sep="\t", usecols=["VariationID", "final"]); T["VariationID"] = T.VariationID.astype(str)
S["VariationID"] = S.VariationID.astype(str); S = S.merge(T, on="VariationID")
cs = S.ClinSig.fillna("").str.lower(); S["cls"] = np.select([cs.str.contains("conflict|uncertain"), cs.str.contains("pathogenic"), cs.str.contains("benign")], ["X", "PLP", "BLB"], "X")
S = S[(S.cls == S.final) & (S.Submitter != "Ambry Genetics") & (S.has_text == 1)]
S["per"] = np.where(S.year <= 2021, "<=2021", np.where(S.year >= 2023, "2023-2026", "2022"))
S = S[S.per != "2022"]
rng = np.random.default_rng(33); picks = []
for (per, cls, fl), g in S.groupby(["per", "cls", "any_insilico"]):
    g = g.sample(frac=1, random_state=int(rng.integers(1e9))); g = g.groupby("Submitter").head(2).head(20)
    picks.append(g.assign(stratum=f"{per}|{cls}|{'flagged' if fl else 'not_flagged'}"))
A = pd.concat(picks); want = set(A.SCV)
txt = {}
with gzip.open(f"{P}/data/clinvar_tab/submission_summary_2026-11.txt.gz", "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("##"): continue
        if line.startswith("#"): hdr = line[1:].rstrip("\n").split("\t"); ix = {h: i for i, h in enumerate(hdr)}; continue
        x = line.rstrip("\n").split("\t")
        if len(x) > ix["SCV"] and x[ix["SCV"]] in want:
            d = x[ix["Description"]] if "Description" in ix else ""
            e = x[ix["ExplanationOfInterpretation"]] if "ExplanationOfInterpretation" in ix and len(x) > ix["ExplanationOfInterpretation"] else ""
            txt[x[ix["SCV"]]] = (d + " " + e).strip()
def window(t, flagged):
    m = ANY.search(t) if flagged else None
    if m: a = max(0, m.start() - 260); b = min(len(t), m.end() + 260); return ("…" if a else "") + t[a:b] + ("…" if b < len(t) else "")
    return t   # unflagged records are kept whole so that a missed mention can be checked
A["text"] = [window(txt.get(s, ""), f == 1) for s, f in zip(A.SCV, A.any_insilico)]
A["n_submitters_in_stratum"] = A.groupby("stratum").Submitter.transform("nunique")
A[["stratum", "SCV", "VariationID", "Submitter", "ClinSig", "year", "text"]].to_csv(f"{P}/results/r2/R_audit_sample.tsv", sep="\t", index=False)
print(A.groupby("stratum").agg(n=("SCV", "size"), submitters=("Submitter", "nunique")).to_string()); print("missing text", sum(1 for s in A.SCV if s not in txt))
