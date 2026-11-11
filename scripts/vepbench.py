#!/usr/bin/env python
"""vepbench: evaluate one missense variant effect score on the vintage- and provenance-annotated ClinVar benchmark.

    python vepbench.py --scores my_scores.tsv.gz --out my_run/
    python vepbench.py --scores examples/phyloP100way_example.tsv.gz --out example_run/      # worked example

Input: a tab-separated file with a variant id column (default `hgvs38`, GRCh38 HGVS g. notation as in
data/vepbench_labels.tsv.gz; `VariationID` also accepted) and a score column (default `score`; higher = more damaging,
use --lower-is-damaging otherwise). Variants absent from the score file are left out, so every number refers to the
variants the score covers.

Steps (Table 1 of the paper): (1) join the scores to the labels and keep rare variants (max population AF < 1e-4);
(2) pooled and within-gene AUROC by label period, for all labels and the provenance-filtered tier, with a gene-identity
null; (3) the change from 2020-2022 to 2023-2026 split into full sets, gene x class cells shared by both periods, and
shared cells reweighted to the 2020-2022 mix; (4) B/LB labels of each recorded-evidence group against all P/LP labels of
the period. Intervals: percentile intervals from a gene-cluster bootstrap (--boot, default 200; seed 2026), the same
resampled genes used for every quantity so that differences are paired.

Outputs in --out: auc_by_period.tsv, composition.tsv, evidence_groups.tsv, summary.txt."""
import argparse, os, sys
import numpy as np, pandas as pd

PERIODS = ["<=2019", "2020-2022", "2023-2026"]
GROUPS = ["Explicit computational evidence", "Text without computational evidence", "No free text", "Ambry only (template)"]

def args_():
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--scores", required=True); p.add_argument("--out", required=True)
    p.add_argument("--labels", default=os.path.join(here, "..", "data", "vepbench_labels.tsv.gz"))
    p.add_argument("--id-col", default="hgvs38"); p.add_argument("--score-col", default="score")
    p.add_argument("--lower-is-damaging", action="store_true"); p.add_argument("--all-frequencies", action="store_true", help="do not restrict to rare variants")
    p.add_argument("--boot", type=int, default=200); p.add_argument("--seed", type=int, default=2026)
    return p.parse_args()

class Bench:
    """AUROCs with gene-level weights, so that a bootstrap replicate is a vector of gene multiplicities."""
    def __init__(self, D, genes):
        self.genes = genes; self.D = D
    def prep(self, G, col="s"):
        g = G[["Gene", "y", col]].dropna()
        lev = pd.Series(g[col].values).rank(method="dense").values.astype(int) - 1
        gi = self.genes.get_indexer(g.Gene.values)
        # within-gene Mann-Whitney numerator and pair count per gene
        num = np.zeros(len(self.genes)); den = np.zeros(len(self.genes))
        r = g.groupby("Gene")[col].rank(method="average").values; y = g.y.values
        t = pd.DataFrame({"gi": gi, "y": y, "r": r})
        a = t.groupby("gi").agg(n1=("y", "sum"), n=("y", "size")); a["n0"] = a.n - a.n1
        rs = t[t.y == 1].groupby("gi").r.sum().reindex(a.index).fillna(0)
        ok = (a.n1 > 0) & (a.n0 > 0)
        num[a.index[ok]] = (rs - a.n1 * (a.n1 + 1) / 2)[ok].values; den[a.index[ok]] = (a.n1 * a.n0)[ok].values
        return dict(lev=lev, K=lev.max() + 1 if len(lev) else 0, y=y.astype(float), gi=gi, num=num, den=den, wv=np.ones(len(g)))
    @staticmethod
    def pooled(pr, w):
        if pr["K"] == 0: return np.nan
        ww = w[pr["gi"]] * pr["wv"]; P = np.bincount(pr["lev"], weights=ww * pr["y"], minlength=pr["K"]); N = np.bincount(pr["lev"], weights=ww * (1 - pr["y"]), minlength=pr["K"])
        d = P.sum() * N.sum(); return (P * (np.cumsum(N) - N + 0.5 * N)).sum() / d if d > 0 else np.nan
    @staticmethod
    def within(pr, w):
        d = (w * pr["den"]).sum(); return (w * pr["num"]).sum() / d if d > 0 else np.nan

def gene_null(G):
    """leave-one-out share of P/LP among the other labelled variants of the same gene; singletons get the overall share."""
    s = G.groupby("Gene").y.transform("sum"); n = G.groupby("Gene").y.transform("size")
    out = ((s - G.y) / (n - 1)).where(n > 1, G.y.mean()); return out

def ci(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return (np.percentile(x, 2.5), np.percentile(x, 97.5)) if len(x) else (np.nan, np.nan)

def main():
    a = args_(); os.makedirs(a.out, exist_ok=True)
    L = pd.read_csv(a.labels, sep="\t", dtype={"VariationID": str})
    S = pd.read_csv(a.scores, sep="\t", dtype={a.id_col: str})
    if a.id_col not in S.columns or a.score_col not in S.columns: sys.exit(f"score file needs columns {a.id_col} and {a.score_col}")
    key = a.id_col if a.id_col in L.columns else sys.exit(f"--id-col must be hgvs38 or VariationID")
    S = S[[a.id_col, a.score_col]].rename(columns={a.score_col: "s"}).dropna().drop_duplicates(a.id_col)
    S["s"] = pd.to_numeric(S.s, errors="coerce") * (-1 if a.lower_is_damaging else 1)
    D = L.merge(S, on=key, how="inner"); D = D[D.s.notna()]
    n_join = len(D)
    if not a.all_frequencies: D = D[D.rare.astype(bool)]
    D["null"] = np.nan
    for per in PERIODS:   # gene-identity null within each label set (period x tier)
        m = D.period == per; D.loc[m, "null"] = gene_null(D[m]).values
    genes = pd.Index(sorted(D.Gene.dropna().unique())); B = Bench(D, genes); ng = len(genes)
    rng = np.random.default_rng(a.seed); W = [np.ones(ng)] + [np.bincount(rng.integers(0, ng, ng), minlength=ng).astype(float) for _ in range(a.boot)]

    # (2) AUROC by period and tier
    rows = []
    for tier, M in [("all labels", np.ones(len(D), bool)), ("provenance-filtered", D.clean.astype(bool).values)]:
        for per in PERIODS:
            G = D[M & (D.period == per).values]
            if G.y.sum() < 10 or (1 - G.y).sum() < 10: continue
            Gn = G.copy(); Gn["null"] = gene_null(G).values if tier != "all labels" else G["null"]
            pr = B.prep(G); pn = B.prep(Gn, "null")
            vp = [B.pooled(pr, w) for w in W]; vw = [B.within(pr, w) for w in W]; vn = [B.pooled(pn, w) for w in W]
            rows.append(dict(tier=tier, period=per, n=len(G), n_plp=int(G.y.sum()), n_blb=int((1 - G.y).sum()), genes=G.Gene.nunique(),
                             genes_with_both=int((pr["den"] > 0).sum()), pooled_auc=vp[0], pooled_lo=ci(vp[1:])[0], pooled_hi=ci(vp[1:])[1],
                             within_auc=vw[0], within_lo=ci(vw[1:])[0], within_hi=ci(vw[1:])[1], gene_null_auc=vn[0]))
    A = pd.DataFrame(rows); A.to_csv(os.path.join(a.out, "auc_by_period.tsv"), sep="\t", index=False, float_format="%.4f")

    # (3) composition: full sets, shared gene x class cells, shared cells reweighted to the 2020-2022 mix
    E = D[D.period == "2020-2022"]; R = D[D.period == "2023-2026"]
    cells = set(zip(E.Gene, E.y)) & set(zip(R.Gene, R.y))
    Es = E[[c in cells for c in zip(E.Gene, E.y)]].copy(); Rs = R[[c in cells for c in zip(R.Gene, R.y)]].copy()
    pE, pR, pEs, pRs = B.prep(E), B.prep(R), B.prep(Es), B.prep(Rs)
    ce = Es.groupby(["Gene", "y"]).size(); cr = Rs.groupby(["Gene", "y"]).size().reindex(ce.index)
    cidx = pd.MultiIndex.from_frame(Rs[["Gene", "y"]]); g_of_cell = genes.get_indexer(ce.index.get_level_values(0))
    out = {k: [] for k in ["full", "within", "shared", "reweighted"]}
    for w in W:
        we = ce * w[g_of_cell]; wr = cr * w[g_of_cell]
        ratio = (we / we.sum()) / (wr / wr.sum()).replace(0, np.nan)
        pRs["wv"] = ratio.reindex(cidx).fillna(0).values
        rw = B.pooled(pRs, w); pRs["wv"] = np.ones(len(Rs)); sh_r = B.pooled(pRs, w); sh_e = B.pooled(pEs, w)
        out["full"].append(B.pooled(pR, w) - B.pooled(pE, w)); out["within"].append(B.within(pR, w) - B.within(pE, w))
        out["shared"].append(sh_r - sh_e); out["reweighted"].append(rw - sh_e)
    C = pd.DataFrame([dict(step=k, gain=v[0], lo=ci(v[1:])[0], hi=ci(v[1:])[1]) for k, v in out.items()])
    C["note"] = ["pooled AUROC, full sets", "within-gene AUROC, full sets", "pooled, shared gene x class cells only",
                 "pooled, shared cells with 2023-2026 reweighted to the 2020-2022 cell mix"]
    C.to_csv(os.path.join(a.out, "composition.tsv"), sep="\t", index=False, float_format="%.4f")
    share_shared = len(Rs) / len(R) if len(R) else np.nan

    # (4) recorded-evidence groups of B/LB labels, against all P/LP of the period
    rows = []
    for per in ["2020-2022", "2023-2026"]:
        G0 = D[D.period == per]
        for grp in GROUPS:
            G = pd.concat([G0[G0.y == 1], G0[(G0.y == 0) & (G0.prov_group == grp)]])
            if (1 - G.y).sum() < 20: continue
            pr = B.prep(G); v = [B.pooled(pr, w) for w in W]
            rows.append(dict(period=per, blb_group=grp, n_blb=int((1 - G.y).sum()), n_plp=int(G.y.sum()), pooled_auc=v[0], lo=ci(v[1:])[0], hi=ci(v[1:])[1]))
    P = pd.DataFrame(rows); P.to_csv(os.path.join(a.out, "evidence_groups.tsv"), sep="\t", index=False, float_format="%.4f")

    # summary
    f = lambda r, k: "%.3f (%.3f-%.3f)" % (r[k + "_auc"], r[k + "_lo"], r[k + "_hi"])
    lines = [f"vepbench: {os.path.basename(a.scores)}; {n_join:,} labelled variants with a score; {len(D):,} analysed "
             f"({'rare only' if not a.all_frequencies else 'all frequencies'}); {a.boot} gene-cluster bootstrap replicates", ""]
    lines.append("AUROC by label period (pooled | within gene | gene-identity null)")
    for _, r in A.iterrows():
        lines.append(f"  {r.tier:20s} {r.period:10s} n={r.n:>6,} (P/LP {r.n_plp:,})  {f(r,'pooled')} | {f(r,'within')} | {r.gene_null_auc:.3f}")
    lines += ["", f"Change 2020-2022 -> 2023-2026 (shared cells hold {share_shared:.0%} of the 2023-2026 labels)"]
    for _, r in C.iterrows(): lines.append(f"  {r.note:75s} {r.gain:+.3f} ({r.lo:+.3f} to {r.hi:+.3f})")
    lines += ["", "B/LB labels by recorded evidence, against all P/LP of the period (pooled AUROC)"]
    for _, r in P.iterrows(): lines.append(f"  {r.period:10s} {r.blb_group:38s} n={r.n_blb:>6,}  {r.pooled_auc:.3f} ({r.lo:.3f}-{r.hi:.3f})")
    lines += ["", "Read pooled AUROC as accuracy across genes and within-gene AUROC as accuracy among variants of the same gene;",
              "a gain that disappears in the shared or reweighted rows reflects which genes and classes were labelled."]
    txt = "\n".join(lines); open(os.path.join(a.out, "summary.txt"), "w").write(txt + "\n"); print(txt)

if __name__ == "__main__":
    main()
