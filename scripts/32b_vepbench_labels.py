#!/usr/bin/env python
"""Build the label table read by vepbench.py (data/vepbench_labels.tsv.gz) and the worked-example score file
(examples/phyloP100way_example.tsv.gz: raw phyloP100way, a UCSC conservation track, for the evaluation variants).
Recorded-evidence groups as in script 32: Ambry only; explicit computational evidence (class-matching PP3/BP4 cited, or a
non-Ambry concordant record mentions computational evidence or names a predictor); text without computational evidence
(a non-Ambry concordant record has free text); no free text."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
SUBA = "Ambry Genetics"; SUBB = "Labcorp Genetics (formerly Invitae), Labcorp"
OUTD = os.environ.get("VEPBENCH_OUT", P)
D = load()
T = pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz", sep="\t", usecols=["hgvs38", "rare", "clean", "period", "explicit_comp", "comp_nonambry", "ambry_only"])
A = pd.read_csv(f"{P}/results/variant_level_annotated.tsv.gz", sep="\t", usecols=["hgvs38", "conc_submitters"]).drop_duplicates("hgvs38")
D = D.merge(T, on="hgvs38").merge(A, on="hgvs38", how="left")
for c in ["rare", "clean", "explicit_comp", "comp_nonambry", "ambry_only"]: D[c] = D[c].fillna(False).astype(bool)
S = pd.read_csv(f"{P}/results/scv_textmining.tsv.gz", sep="\t", usecols=["VariationID", "Submitter", "ClinSig", "has_text"])
cs = S.ClinSig.fillna("").str.lower(); S["cls"] = np.select([cs.str.contains("conflict|uncertain"), cs.str.contains("pathogenic"), cs.str.contains("benign")], ["X", "PLP", "BLB"], "X")
S["VariationID"] = S.VariationID.astype(str); D["VariationID"] = D.VariationID.astype(str)
S = S.merge(D[["VariationID", "final"]], on="VariationID"); S = S[(S.cls == S.final) & (S.Submitter != SUBA)]
txt = S.groupby("VariationID").has_text.max(); D["nonambry_text"] = D.VariationID.map(txt).fillna(0).astype(int) > 0
G = ["Explicit computational evidence", "Text without computational evidence", "No free text", "Ambry only (template)"]
D["prov_group"] = np.select([D.ambry_only, D.explicit_comp | D.comp_nonambry, D.nonambry_text], [G[3], G[0], G[1]], G[2])
D["source"] = np.where(D.conc_submitters == SUBA, "Ambry only", np.where(D.conc_submitters == SUBB, "Labcorp only", "Other or several"))
cols = ["hgvs38", "VariationID", "Gene", "final", "y", "vintage_first", "period", "stars", "rare", "clean", "source", "prov_group"]
os.makedirs(f"{OUTD}/data", exist_ok=True); os.makedirs(f"{OUTD}/examples", exist_ok=True)
D[cols].to_csv(f"{OUTD}/data/vepbench_labels.tsv.gz", sep="\t", index=False)
D[["hgvs38", "phyloP100way_raw"]].rename(columns={"phyloP100way_raw": "score"}).dropna().to_csv(f"{OUTD}/examples/phyloP100way_example.tsv.gz", sep="\t", index=False)
print(len(D), D.prov_group.value_counts().to_dict(), D.phyloP100way_raw.notna().sum())
