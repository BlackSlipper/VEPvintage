#!/usr/bin/env python
"""Recalibration of PP3/BP4 thresholds on labels up to 2022 and on 2023-2026 labels, closer to Pejaver et al. 2022.
Changes from script 30: (i) evidence strengths use the likelihood-ratio thresholds of Pejaver et al. (prior 4.41%):
2.406 supporting, 5.790 moderate, 33.53 strong (benign: reciprocals); (ii) the main setting keeps only genes with at least
one P/LP label in the evaluation set, as ClinGen kept genes with at least one pathogenic ClinVar variant; the unrestricted
setting is kept as a sensitivity analysis; (iii) shares are reported for P/LP and B/LB; (iv) an out-of-sample check splits
genes into two halves, calibrates on one half and applies the thresholds to the rare 2023-2026 variants of the other half.
Local lr+: 400 score quantiles, symmetric window with at least 100 rare variants, pooled with a 0.5 correction; Mantel-Haenszel
within genes; gene-cluster bootstrap (B=VEP_BOOT, default 1000, seed 2026); ClinGen rule = one-sided 95% lower bound.
Outputs results/r2/O2_threshold_shift.tsv, O2_threshold_shift_categories.tsv, O2_threshold_oos.tsv."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from multiprocessing import Pool
B = int(os.environ.get("VEP_BOOT", "1000")); OUT = f"{P}/results/r2"; WIN = 100; NG = 400
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
D = load()
T = pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz", sep="\t", usecols=["hgvs38", "rare", "clean", "period"])
D = D.merge(T, on="hgvs38"); D["rare"] = D.rare.fillna(False).astype(bool); D["clean"] = D.clean.fillna(False).astype(bool)
TH = pd.read_csv(f"{P}/meta/clingen_thresholds.tsv", sep="\t"); TH["lo"] = TH.lo.astype(float); TH["hi"] = TH.hi.astype(float)
PLP_GENES = set(D.loc[D.y == 1, "Gene"].dropna())          # genes with >=1 P/LP label anywhere in the evaluation set
R0 = D[D.rare].copy(); R0["recent"] = R0.period == "2023-2026"
genes = pd.Index(sorted(D.Gene.dropna().unique())); ng = len(genes)
rng = np.random.default_rng(2026); W = np.stack([np.bincount(rng.integers(0, ng, ng), minlength=ng) for _ in range(B)]).astype(float)
FOLD = np.random.default_rng(7).permutation(ng) % 2       # gene halves for the out-of-sample check
LRT = {"Supporting": 2.406, "Moderate": 5.790, "Strong": 33.53}
SETTINGS = {"genes_with_PLP": R0[R0.Gene.isin(PLP_GENES)], "all_genes": R0}

def windows(s):
    n = len(s); grid = np.quantile(s, np.linspace(0.005, 0.995, NG)); lo = np.empty(NG, int); hi = np.empty(NG, int)
    step = (s[-1] - s[0]) / 2000 if s[-1] > s[0] else 1e-6
    for k, g in enumerate(grid):
        h = 0.0
        while True:
            a = np.searchsorted(s, g - h, side="left"); b = np.searchsorted(s, g + h, side="right")
            if b - a >= WIN or (a == 0 and b == n): break
            h += step
        lo[k] = a; hi[k] = b
    return grid, lo, hi
def pooled_lr(lo, hi, y, gi, w):
    ww = w[gi]; cP = np.concatenate([[0], np.cumsum(ww * y)]); cN = np.concatenate([[0], np.cumsum(ww * (1 - y))])
    a = cP[hi] - cP[lo]; b = cN[hi] - cN[lo]; return ((a + 0.5) / (cP[-1] + 0.5)) / ((b + 0.5) / (cN[-1] + 0.5))
def mh_terms(lo, hi, y, gi):
    n1 = np.bincount(gi, weights=y, minlength=ng); n0 = np.bincount(gi, weights=1 - y, minlength=ng); N = np.maximum(n1 + n0, 1)
    A = np.zeros((ng, NG)); Bm = np.zeros((ng, NG))
    for k in range(NG):
        g = gi[lo[k]:hi[k]]; yy = y[lo[k]:hi[k]]
        a = np.bincount(g, weights=yy, minlength=ng); b = np.bincount(g, weights=1 - yy, minlength=ng)
        A[:, k] = a * n0 / N; Bm[:, k] = b * n1 / N
    return A, Bm
def thresholds(grid, lr):
    out = {}
    for ev, t in LRT.items():
        suf = np.cumprod((lr >= t)[::-1])[::-1].astype(bool); idx = np.where(suf)[0]; out["PP3_" + ev] = grid[idx[0]] if len(idx) else np.nan
        pre = np.cumprod((1 / lr) >= t).astype(bool); idx = np.where(pre)[0]; out["BP4_" + ev] = grid[idx[-1]] if len(idx) else np.nan
    return out
def levels(s, th):
    lv = np.zeros(len(s), int)
    for i, ev in enumerate(["Supporting", "Moderate", "Strong"], 1):
        t = th["PP3_" + ev]
        if np.isfinite(t): lv[s >= t] = i
    for i, ev in enumerate(["Supporting", "Moderate", "Strong"], 1):
        t = th["BP4_" + ev]
        if np.isfinite(t): lv[s <= t] = -i
    return lv
def prep(G, rc, sign):
    g = G[["Gene", "y", rc]].dropna()
    s = sign * g[rc].values.astype(float); o = np.argsort(s, kind="mergesort")
    return s[o], g.y.values.astype(float)[o], genes.get_indexer(g.Gene.values)[o], len(g)
def shares(s, y, th):
    lv = levels(s, th); out = {}
    for cls, m in [("PLP", y == 1), ("BLB", y == 0)]:
        out[f"{cls}_pp3_any"] = (lv[m] > 0).mean(); out[f"{cls}_pp3_modplus"] = (lv[m] >= 2).mean(); out[f"{cls}_pp3_strong"] = (lv[m] >= 3).mean()
        out[f"{cls}_bp4_any"] = (lv[m] < 0).mean(); out[f"{cls}_bp4_modplus"] = (lv[m] <= -2).mean(); out[f"n_{cls}"] = int(m.sum())
    return out

tools = TH[["tool", "rawcol", "direction"]].drop_duplicates().values.tolist()
def job(k):
    tool, rc, direction = k; sign = 1 if direction == "+" else -1; out = []; cat = []; oos = []
    for setting, R in SETTINGS.items():
        SETS = {"le2022_all": R[~R.recent], "2023_2026_all": R[R.recent], "le2022_filtered": R[~R.recent & R.clean], "2023_2026_filtered": R[R.recent & R.clean]}
        thr = {}
        for nm, G in SETS.items():
            s, y, gi, n = prep(G, rc, sign)
            if n < WIN * 2: continue
            grid, lo, hi = windows(s)
            lrp = pooled_lr(lo, hi, y, gi, np.ones(ng)); A, Bm = mh_terms(lo, hi, y, gi); lrm = (A.sum(0) + 0.5) / (Bm.sum(0) + 0.5)
            LRP = np.stack([pooled_lr(lo, hi, y, gi, W[b]) for b in range(B)]); LRM = (W @ A + 0.5) / (W @ Bm + 0.5)
            est = {"pooled": thresholds(grid, lrp), "within_gene": thresholds(grid, lrm),
                   "pooled_lower95": thresholds(grid, np.percentile(LRP, 5, axis=0)), "within_gene_lower95": thresholds(grid, np.percentile(LRM, 5, axis=0))}
            bsP = [thresholds(grid, LRP[b]) for b in range(B)]; bsM = [thresholds(grid, LRM[b]) for b in range(B)]
            thr[nm] = est
            for kind, bs in [("pooled", bsP), ("within_gene", bsM), ("pooled_lower95", bsP), ("within_gene_lower95", bsM)]:
                for ev, v in est[kind].items():
                    arr = np.array([d[ev] for d in bs], float); fin = arr[np.isfinite(arr)]
                    lo_, hi_ = (np.percentile(fin, [2.5, 97.5]) if len(fin) else (np.nan, np.nan))
                    out.append(dict(tool=tool, setting=setting, set=nm, kind=kind, n=n, n_plp=int(y.sum()), evidence=ev,
                                    threshold=sign * v if np.isfinite(v) else np.nan, lo=sign * min(lo_, hi_) if np.isfinite(lo_) else np.nan,
                                    hi=sign * max(lo_, hi_) if np.isfinite(hi_) else np.nan, reached_share=len(fin) / B))
        if "le2022_all" in thr and "2023_2026_all" in thr:
            s, y, gi, n = prep(R[R.recent], rc, sign)
            for kind in ["pooled", "pooled_lower95", "within_gene", "within_gene_lower95"]:
                for src in ["le2022_all", "2023_2026_all"]:
                    cat.append(dict(tool=tool, setting=setting, kind=kind, thresholds_from=src, applied_to="rare_2023_2026", **shares(s, y, thr[src][kind])))
        for f in (0, 1):
            cal_genes = set(genes[FOLD == f]); app_genes = set(genes[FOLD != f])
            for src, G in [("le2022_all", R[~R.recent]), ("2023_2026_all", R[R.recent])]:
                s, y, gi, n = prep(G[G.Gene.isin(cal_genes)], rc, sign)
                if n < WIN * 2: continue
                grid, lo, hi = windows(s); th = thresholds(grid, pooled_lr(lo, hi, y, gi, np.ones(ng)))
                s2, y2, _, n2 = prep(R[R.recent & R.Gene.isin(app_genes)], rc, sign)
                oos.append(dict(tool=tool, setting=setting, fold=f, thresholds_from=src, n_cal=n, **{k: v for k, v in th.items()}, **shares(s2, y2, th)))
    return out, cat, oos

with Pool(min(22, len(tools))) as pool: res = pool.map(job, tools)
rows, cat_rows, oos_rows = [], [], []
for a, b, c in res: rows += a; cat_rows += b; oos_rows += c
Q = pd.DataFrame(rows); Q.to_csv(f"{OUT}/O2_threshold_shift.tsv", sep="\t", index=False)
C = pd.DataFrame(cat_rows); C.to_csv(f"{OUT}/O2_threshold_shift_categories.tsv", sep="\t", index=False)
O = pd.DataFrame(oos_rows); O.to_csv(f"{OUT}/O2_threshold_oos.tsv", sep="\t", index=False)
for setting in SETTINGS:
    piv = Q[(Q.setting == setting) & (Q.kind == "pooled") & Q.evidence.str.startswith("PP3")].pivot_table(index=["tool", "evidence"], columns="set", values="threshold")
    print("==", setting, "pooled PP3 thresholds"); print(piv.round(3).to_string())
    c = C[(C.setting == setting) & (C.kind == "pooled")].pivot_table(index="tool", columns="thresholds_from", values=["PLP_pp3_modplus", "BLB_pp3_modplus", "PLP_pp3_strong"])
    print(c.round(3).to_string()); print("median", c.median().round(3).to_dict())
    o = O[O.setting == setting].groupby(["tool", "thresholds_from"])[["PLP_pp3_modplus", "BLB_pp3_modplus"]].mean().unstack()
    print("== OOS", setting); print(o.round(3).to_string()); print("median", o.median().round(3).to_dict())
print("DONE")
