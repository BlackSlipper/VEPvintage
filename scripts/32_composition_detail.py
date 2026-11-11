#!/usr/bin/env python
"""Detail for the composition and provenance analyses (rare variants scored by all 38 core predictors).
(A) Recent (2023-2026) labels inside and outside the gene x class cells shared with 2020-2022: counts by class, source
    (Ambry only, Labcorp only, other) and whether the gene has any P/LP label; median core AUROC on the full recent set,
    the shared cells and the excluded cells.
(B) Provenance groups by what the concordant submission records say: explicit (a non-Ambry record mentions computational
    evidence or a predictor, or the class-matching PP3/BP4 code is cited), text without computational evidence, no free
    text at all (non-Ambry records only), and Ambry-only labels (template text). Pooled and within-gene AUROC per group and
    period, with gene-cluster bootstrap intervals (B=VEP_BOOT, default 500, seed 2026) for REVEL and the median predictor.
Outputs results/r2/Q_excluded_cells.tsv, Q_excluded_genes.tsv, Q_excluded_auc.tsv, Q_provenance_groups.tsv"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
B = int(os.environ.get("VEP_BOOT", "500")); OUT = f"{P}/results/r2"
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
SUBA = "Ambry Genetics"; SUBB = "Labcorp Genetics (formerly Invitae), Labcorp"
D = load()
T = pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz", sep="\t", usecols=["hgvs38", "rare", "clean", "period", "explicit_comp", "comp_nonambry", "ambry_only"])
A = pd.read_csv(f"{P}/results/variant_level_annotated.tsv.gz", sep="\t", usecols=["hgvs38", "conc_submitters"]).drop_duplicates("hgvs38")
D = D.merge(T, on="hgvs38").merge(A, on="hgvs38", how="left")
for c in ["rare", "clean", "explicit_comp", "comp_nonambry", "ambry_only"]: D[c] = D[c].fillna(False).astype(bool)
tools = [t for t in TOOLS if t in D.columns]; cov = D[tools].notna().mean(); core = [t for t in tools if cov[t] >= 0.90]
RC = D[D.rare & D[core].notna().all(axis=1)].copy()
RC["source"] = np.where(RC.conc_submitters == SUBA, "Ambry only", np.where(RC.conc_submitters == SUBB, "Labcorp only", "Other or several"))
genes = pd.Index(sorted(D.Gene.dropna().unique())); ng = len(genes)
rng = np.random.default_rng(2026); W = np.stack([np.bincount(rng.integers(0, ng, ng), minlength=ng) for _ in range(B)]).astype(float)

def prep(G, t):
    g = G[["Gene", "y", t]].dropna(); lev = pd.Series(g[t].values).rank(method="dense").values.astype(int) - 1
    return lev, lev.max() + 1, g.y.values.astype(float), genes.get_indexer(g.Gene.values)
def pauc(pr, w):
    lev, K, y, gi = pr; ww = w[gi]; Pk = np.bincount(lev, weights=ww * y, minlength=K); Nk = np.bincount(lev, weights=ww * (1 - y), minlength=K)
    cN = np.cumsum(Nk) - Nk; d = Pk.sum() * Nk.sum(); return (Pk * (cN + 0.5 * Nk)).sum() / d if d > 0 else np.nan
def within(G, t):
    """within-gene AUROC pooled over P/LP-B/LB pairs from the same gene."""
    g = G[["Gene", "y", t]].dropna(); num = 0.0; den = 0.0
    for gene, h in g.groupby("Gene"):
        yy = h.y.values; s = h[t].values.astype(float); p = s[yy == 1]; n = s[yy == 0]
        if len(p) == 0 or len(n) == 0: continue
        num += (p[:, None] > n[None, :]).sum() + 0.5 * (p[:, None] == n[None, :]).sum(); den += len(p) * len(n)
    return num / den if den > 0 else np.nan
def summarize(G, label, extra, boot=True):
    prs = {t: prep(G, t) for t in core}
    pt = np.array([pauc(prs[t], np.ones(ng)) for t in core]); med = np.nanmedian(pt); rev = pt[core.index("REVEL")]
    row = dict(label=label, **extra, n=len(G), n_plp=int(G.y.sum()), n_blb=int((1 - G.y).sum()), genes=G.Gene.nunique(), auc_REVEL=rev, auc_median=med)
    if boot and G.y.sum() >= 20 and (1 - G.y).sum() >= 20:
        bm = []; br = []
        for b in range(B):
            v = np.array([pauc(prs[t], W[b]) for t in core]); bm.append(np.nanmedian(v)); br.append(v[core.index("REVEL")])
        row.update(auc_median_lo=np.nanpercentile(bm, 2.5), auc_median_hi=np.nanpercentile(bm, 97.5), auc_REVEL_lo=np.nanpercentile(br, 2.5), auc_REVEL_hi=np.nanpercentile(br, 97.5))
    wv = [within(G, t) for t in ["REVEL", "AlphaMissense", "ESM1b", "GERP++"]]
    pairs = G.groupby("Gene").y.agg(lambda s: (s.sum() > 0) & ((1 - s).sum() > 0))
    row.update(within_REVEL=wv[0], within_AlphaMissense=wv[1], within_ESM1b=wv[2], within_GERP=wv[3], genes_with_pairs=int(pairs.sum()),
               n_in_pair_genes=int(G.Gene.isin(pairs[pairs].index).sum()))
    return row

# ---------- (A) shared vs excluded cells ----------
E = RC[RC.period == "2020-2022"]; L = RC[RC.period == "2023-2026"].copy()
cells = set(zip(E.Gene, E.y)); L["shared_cell"] = [(g, y) in cells for g, y in zip(L.Gene, L.y)]
L["gene_has_PLP_2020_22"] = L.Gene.isin(set(E.loc[E.y == 1, "Gene"])); L["gene_has_PLP_any"] = L.Gene.isin(set(D.loc[D.y == 1, "Gene"]))
L["gene_in_2020_22"] = L.Gene.isin(set(E.Gene))
tab = L.groupby(["shared_cell", "final", "source", "gene_in_2020_22", "gene_has_PLP_2020_22", "gene_has_PLP_any"]).size().rename("n").reset_index()
tab.to_csv(f"{OUT}/Q_excluded_cells.tsv", sep="\t", index=False)
print("shared share", L.shared_cell.mean().round(4), len(L), "shared n", int(L.shared_cell.sum()), "cells", len(set(zip(L.Gene, L.y)) & cells))
print(L.groupby(["shared_cell", "final"]).size().unstack().to_string())
print(pd.crosstab([L.shared_cell, L.final], L.source, normalize="index").round(3).to_string())
print(pd.crosstab([L.shared_cell, L.final], L.source).to_string())
print(pd.crosstab([L.shared_cell, L.final], [L.gene_in_2020_22, L.gene_has_PLP_any]).to_string())
X = L[~L.shared_cell]
for cls in ["PLP", "BLB"]:
    g = X[X.final == cls].Gene.value_counts(); print("excluded", cls, "genes", len(g), "top10 share", round(g.head(10).sum() / g.sum(), 3)); print(g.head(15).to_string())
X.groupby(["final", "Gene"]).size().rename("n").reset_index().sort_values("n", ascending=False).to_csv(f"{OUT}/Q_excluded_genes.tsv", sep="\t", index=False)
ar = [summarize(L, "recent_all", {}), summarize(L[L.shared_cell], "recent_shared_cells", {}), summarize(L[~L.shared_cell], "recent_excluded_cells", {}), summarize(E, "early_all", {})]
pd.DataFrame(ar).to_csv(f"{OUT}/Q_excluded_auc.tsv", sep="\t", index=False)
print(pd.DataFrame(ar).round(4).to_string())

# ---------- (B) provenance groups ----------
S = pd.read_csv(f"{P}/results/scv_textmining.tsv.gz", sep="\t", usecols=["VariationID", "Submitter", "ClinSig", "has_text"])
cs = S.ClinSig.fillna("").str.lower(); S["cls"] = np.select([cs.str.contains("conflict|uncertain"), cs.str.contains("pathogenic"), cs.str.contains("benign")], ["X", "PLP", "BLB"], "X")
S["VariationID"] = S.VariationID.astype(str)
V = RC[["VariationID", "final"]].copy(); V["VariationID"] = V.VariationID.astype(str)
S = S.merge(V, on="VariationID"); S = S[(S.cls == S.final) & (S.Submitter != SUBA)]
txt = S.groupby("VariationID").has_text.max()
RC["VariationID"] = RC.VariationID.astype(str); RC["nonambry_text"] = RC.VariationID.map(txt).fillna(0).astype(int) > 0
GROUPS = ["Explicit computational evidence", "Text without computational evidence", "No free text", "Ambry only (template)"]
RC["prov_group"] = np.select([RC.ambry_only, RC.explicit_comp | RC.comp_nonambry, RC.nonambry_text], [GROUPS[3], GROUPS[0], GROUPS[1]], GROUPS[2])
RC[["hgvs38", "Gene", "final", "period", "prov_group", "source"]].to_csv(f"{OUT}/Q_provenance_membership.tsv.gz", sep="\t", index=False)
rows = []
for per in ["2020-2022", "2023-2026"]:
    G0 = RC[RC.period == per]
    print(per); print(G0.groupby(["prov_group", "final"]).size().unstack().to_string())
    for grp in GROUPS:
        G = G0[G0.prov_group == grp]
        if G.y.sum() >= 10 and (1 - G.y).sum() >= 10: rows.append(summarize(G, grp, dict(period=per, comparison="within group")))
        else: rows.append(dict(label=grp, period=per, comparison="within group", n=len(G), n_plp=int(G.y.sum()), n_blb=int((1 - G.y).sum())))
    for grp in GROUPS:   # all P/LP of the period against the B/LB of one group
        G = pd.concat([G0[G0.y == 1], G0[(G0.y == 0) & (G0.prov_group == grp)]])
        if (1 - G.y).sum() >= 20: rows.append(summarize(G, grp, dict(period=per, comparison="all P/LP vs group B/LB")))
Q = pd.DataFrame(rows); Q.to_csv(f"{OUT}/Q_provenance_groups.tsv", sep="\t", index=False)
print(Q.round(4).to_string()); print("DONE")
