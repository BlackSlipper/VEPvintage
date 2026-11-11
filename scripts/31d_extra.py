#!/usr/bin/env python
"""Extra checks for the discontinuity analysis: leave-three-out, difference without two pre-2022 submitters,
and how often PP3/BP4 mentions are negated in the submission text. Outputs results/r2/P_rd_extra.tsv"""
import sys, os, re, gzip; sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)),"31_threshold_rd.py")).read().split("# ------------------------------------------------------------------ analysis sets")[0])
X["per"] = np.where(X.year <= 2021, "<=2021", np.where(X.year == 2022, "2022", "2023-2026"))
rng = np.random.default_rng(SEED); out = []
d = X[X.cls.isin(["PLP", "VUS"]) & X.REVEL_raw.notna()]
x = d.REVEL_raw.values.astype(float); y = d.PP3.values.astype(float); g = d.gi.values; per = d.per.values; sub = d.Submitter.values
three = d.Submitter.str.contains("Mayo Clinic Laboratories|University of Leipzig|Variantyx", regex=True).values
two = d.Submitter.str.contains("Laboratory for Molecular Medicine|Ljubljana", regex=True).values
late = per == "2023-2026"; early = per == "<=2021"
for name, m in [("PP3 2023-2026 all", late), ("PP3 2023-2026 without the three step-like labs", late & ~three)]:
    ct = Cut(x[m], y[m], g[m], 0.644, H); j = ct.jump()[0]
    bs = np.array([ct.jump(gene_weights(rng))[0] for _ in range(B)])
    out.append(dict(analysis=name, est=j, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), nL=ct.nL, nR=ct.nR))
for name, m0 in [("PP3 diff late-early, all", early), ("PP3 diff late-early, early without two submitters", early & ~two)]:
    c0 = Cut(x[m0], y[m0], g[m0], 0.644, H); c1 = Cut(x[late], y[late], g[late], 0.644, H)
    bs = []
    for _ in range(B):
        gw = gene_weights(rng); bs.append(c1.jump(gw)[0] - c0.jump(gw)[0])
    out.append(dict(analysis=name, est=c1.jump()[0] - c0.jump()[0], lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), nL=c0.nL, nR=c0.nR))
# negation of PP3/BP4 mentions in the free text
cited = set(X[(X.PP3 == 1) | (X.BP4 == 1)].SCV)
neg = re.compile(r"(?:\b(?:not|no|nor|without|absent|unmet|n/a)\b[^.;]{0,25}\b(pp3|bp4)\b|\b(pp3|bp4)\b[^.;]{0,25}\b(?:not (?:met|applied|used|applicable)|was not|is not|does not|unmet|not)\b)", re.I)
rows = []
with gzip.open(f"{P}/data/clinvar_tab/submission_summary_2026-11.txt.gz", "rt", errors="replace") as fh:
    hdr = None
    for line in fh:
        if line.startswith("#VariationID"):
            hdr = line.lstrip("#").rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(hdr)}; continue
        if hdr is None or line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        scv = f[ix["SCV"]]
        if scv in cited:
            desc = f[ix["Description"]] if "Description" in ix else ""
            rows.append((scv, bool(neg.search(desc))))
Ng = pd.DataFrame(rows, columns=["SCV", "negated"])
Z = X[X.SCV.isin(cited)].merge(Ng, on="SCV", how="left")
for code in ["PP3", "BP4"]:
    z = Z[Z[code] == 1]
    for p in ["<=2021", "2022", "2023-2026"]:
        zz = z[z.per == p]
        out.append(dict(analysis=f"{code} mentions with negation, {p}", est=zz.negated.mean(), nL=len(zz)))
    zz = z[(z.per == "<=2021") & z.Submitter.str.contains("Laboratory for Molecular Medicine|Ljubljana", regex=True)]
    out.append(dict(analysis=f"{code} mentions with negation, <=2021, the two submitters", est=zz.negated.mean(), nL=len(zz)))
# PP3 step after dropping negated mentions
nx = set(Ng[Ng.negated].SCV)
y2 = np.where(d.SCV.isin(nx).values, 0.0, y)
for name, m in [("PP3 2023-2026, negated mentions not counted", late), ("PP3 <=2021, negated mentions not counted", early)]:
    ct = Cut(x[m], y2[m], g[m], 0.644, H); j = ct.jump()[0]
    bs = np.array([ct.jump(gene_weights(rng))[0] for _ in range(B)])
    out.append(dict(analysis=name, est=j, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), nL=ct.nL, nR=ct.nR))
O = pd.DataFrame(out); O.to_csv(f"{R2}/P_rd_extra.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(O.round(4).to_string())
ex = Z[(Z.PP3 == 1) & Z.negated & (Z.per == "<=2021")].SCV.head(3).tolist(); print("examples", ex)
