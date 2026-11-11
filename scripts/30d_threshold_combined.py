#!/usr/bin/env python
"""Recalibration check that combines held-out genes with calibration sets matched on class counts (main setting of 30b:
rare labels in genes with >=1 P/LP label; pooled point estimates). In each of NREP replicates (seed 31) the genes are split
at random into a calibration half and an evaluation half. In the calibration half, thresholds are estimated from labels up to
2022 and from 2023-2026 labels, either as they are (design 'split') or after subsampling both periods to the same numbers of
P/LP and of B/LB labels (design 'split_matched'). All four threshold sets are applied to the evaluation half, separately to
its rare 2023-2026 labels and to its rare labels up to 2022. Reported per replicate: thresholds, the shares of P/LP and B/LB
variants reaching each strength (30b shares()), and the likelihood ratio achieved in each PP3/BP4 evidence interval
(supporting, moderate, strong; BP4 as the inverse ratio) with the counts behind it.
Outputs results/r2/O3_threshold_combined.tsv; summary printed (per-tool median over replicates, then median over tools)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "30b_threshold_shift_v2.py")).read()
exec(src.split("tools = TH[")[0].replace('B = int(os.environ.get("VEP_BOOT", "1000"))', 'B = 2'))
NREP = int(os.environ.get("NREP", "30"))
tools = TH[["tool", "rawcol", "direction"]].drop_duplicates().values.tolist()
R = SETTINGS["genes_with_PLP"]
rngs = np.random.default_rng(31); SPLITS = [rngs.permutation(ng) % 2 for _ in range(NREP)]; SEEDS = rngs.integers(1e9, size=(NREP, 2))
EV = ["Supporting", "Moderate", "Strong"]
def interval_lrs(s, y, th):
    out = {}; n1 = y.sum(); n0 = len(y) - n1
    cut = [th["PP3_" + e] for e in EV]
    for k, e in enumerate(EV):
        lo = cut[k]; hi = next((c for c in cut[k + 1:] if np.isfinite(c)), np.inf)
        m = (s >= lo) & (s < hi) if np.isfinite(lo) else np.zeros(len(s), bool)
        a = y[m].sum(); b = m.sum() - a
        out[f"PP3_{e}_LR"] = ((a + 0.5) / (n1 + 0.5)) / ((b + 0.5) / (n0 + 0.5)) if m.any() else np.nan; out[f"PP3_{e}_n"] = int(m.sum())
    cut = [th["BP4_" + e] for e in EV]
    for k, e in enumerate(EV):
        hi = cut[k]; lo = next((c for c in cut[k + 1:] if np.isfinite(c)), -np.inf)
        m = (s <= hi) & (s > lo) if np.isfinite(hi) else np.zeros(len(s), bool)
        a = y[m].sum(); b = m.sum() - a
        out[f"BP4_{e}_LR"] = ((b + 0.5) / (n0 + 0.5)) / ((a + 0.5) / (n1 + 0.5)) if m.any() else np.nan; out[f"BP4_{e}_n"] = int(m.sum())
    return out
def job(k):
    tool, rc, direction = k; sign = 1 if direction == "+" else -1; out = []
    G = R[R[rc].notna()]
    for r in range(NREP):
        cal = set(genes[SPLITS[r] == 0]); evl = set(genes[SPLITS[r] == 1])
        EA = G[~G.recent & G.Gene.isin(cal)]; LA = G[G.recent & G.Gene.isin(cal)]
        EB = G[~G.recent & G.Gene.isin(evl)]; LB = G[G.recent & G.Gene.isin(evl)]
        npl = min(int(EA.y.sum()), int(LA.y.sum())); nbl = min(int((EA.y == 0).sum()), int((LA.y == 0).sum()))
        def match(X, seed):
            return pd.concat([X[X.y == 1].sample(n=npl, random_state=int(seed)), X[X.y == 0].sample(n=nbl, random_state=int(seed) + 1)])
        sets = [("le2022", "split", EA), ("2023_2026", "split", LA), ("le2022", "split_matched", match(EA, SEEDS[r, 0])), ("2023_2026", "split_matched", match(LA, SEEDS[r, 1]))]
        evals = {"eval_2023_2026": prep(LB, rc, sign), "eval_le2022": prep(EB, rc, sign)}
        for src_, design, X in sets:
            s, y, gi, n = prep(X, rc, sign); grid, lo, hi = windows(s); th = thresholds(grid, pooled_lr(lo, hi, y, gi, np.ones(ng)))
            for ev, (s2, y2, _, n2) in evals.items():
                out.append(dict(tool=tool, rep=r, thresholds_from=src_, design=design, eval_set=ev, n_cal=n, n_cal_plp=int(y.sum()), n_eval=n2, n_eval_plp=int(y2.sum()),
                                **{kk: (sign * v if np.isfinite(v) else np.nan) for kk, v in th.items()}, **shares(s2, y2, th), **interval_lrs(s2, y2, th)))
    return out
with Pool(min(11, len(tools))) as pool: res = pool.map(job, tools)
Q = pd.DataFrame([r for x in res for r in x]); Q.to_csv(f"{OUT}/O3_threshold_combined.tsv", sep="\t", index=False)
cols = ["PLP_pp3_any", "PLP_pp3_modplus", "PLP_pp3_strong", "BLB_pp3_any", "BLB_pp3_modplus", "BLB_bp4_modplus", "PP3_Supporting_LR", "PP3_Moderate_LR", "PP3_Strong_LR", "BP4_Supporting_LR", "BP4_Moderate_LR"]
pt = Q.groupby(["design", "eval_set", "thresholds_from", "tool"])[cols].median()
print(pt.groupby(["design", "eval_set", "thresholds_from"]).median().round(3).to_string())
w = Q[Q.eval_set == "eval_2023_2026"].pivot_table(index=["design", "tool", "rep"], columns="thresholds_from", values="PLP_pp3_modplus")
print("share of replicates with higher recent P/LP moderate+ share:"); print((w["2023_2026"] > w["le2022"]).groupby(["design", "tool"]).mean().unstack(0).round(2).to_string())
print(Q[Q.tool == "REVEL"].groupby(["design", "eval_set", "thresholds_from"])[["PP3_Supporting", "PP3_Moderate", "PP3_Strong", "PLP_pp3_modplus", "PP3_Supporting_LR", "PP3_Moderate_LR"]].quantile([0.025, 0.5, 0.975]).unstack().round(3).to_string())
print("DONE")
