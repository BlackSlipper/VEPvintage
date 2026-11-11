#!/usr/bin/env python
"""Pooled and within-gene AUROC of REVEL, ESM1b and GERP++ by label period on the rare variants scored by all 38 core
predictors (the set behind Fig. 2a-b and the text), with gene-cluster bootstrap intervals (B=VEP_BOOT, default 500, seed 2026).
Output results/r2/Q_fig2c_common.tsv (used by Fig. 2c)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "32_composition_detail.py")).read()
exec(src.split("# ---------- (A) shared vs excluded cells ----------")[0])
def within_w(pr_g, w):
    num, den = pr_g; d = (w * den).sum(); return (w * num).sum() / d if d > 0 else np.nan
def within_prep(G, t):
    g = G[["Gene", "y", t]].dropna(); num = np.zeros(ng); den = np.zeros(ng)
    for gene, h in g.groupby("Gene"):
        yy = h.y.values; s = h[t].values.astype(float); p = s[yy == 1]; n = s[yy == 0]
        if len(p) == 0 or len(n) == 0: continue
        i = genes.get_loc(gene); num[i] = (p[:, None] > n[None, :]).sum() + 0.5 * (p[:, None] == n[None, :]).sum(); den[i] = len(p) * len(n)
    return num, den
rows = []
for per in ["<=2019", "2020-2022", "2023-2026"]:
    G = RC[RC.period == per]
    for t in ["REVEL", "ESM1b", "GERP++"]:
        pr = prep(G, t); wp = within_prep(G, t)
        bp = [pauc(pr, W[b]) for b in range(B)]; bw = [within_w(wp, W[b]) for b in range(B)]
        rows.append(dict(period=per, tool=t, n=len(G), n_plp=int(G.y.sum()), pooled=pauc(pr, np.ones(ng)), pooled_lo=np.nanpercentile(bp, 2.5), pooled_hi=np.nanpercentile(bp, 97.5),
                         within=within_w(wp, np.ones(ng)), within_lo=np.nanpercentile(bw, 2.5), within_hi=np.nanpercentile(bw, 97.5)))
Q = pd.DataFrame(rows); Q.to_csv(f"{OUT}/Q_fig2c_common.tsv", sep="\t", index=False); print(Q.round(4).to_string())
