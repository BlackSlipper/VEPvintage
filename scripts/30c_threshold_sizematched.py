#!/usr/bin/env python
"""Size-matched check for the recalibration contrast (script 30b, main setting: genes with >=1 P/LP label, pooled point
estimates). The 2023-2026 set has far more B/LB labels than the set up to 2022, and the threshold rule (every window above the
threshold must reach the likelihood ratio) is sensitive to noise, so the sets are matched on class counts: in each of 20
replicates (seed 11) the B/LB labels of 2023-2026 are subsampled to the number up to 2022 and the P/LP labels up to 2022 to the
number in 2023-2026. Thresholds from each matched set are applied to all rare 2023-2026 variants.
Output results/r2/O2_threshold_sizematched.tsv (per tool and replicate) and a printed summary."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "30b_threshold_shift_v2.py")).read()
exec(src.split("tools = TH[")[0].replace('B = int(os.environ.get("VEP_BOOT", "1000"))', 'B = 2'))
NREP = 20
tools = TH[["tool", "rawcol", "direction"]].drop_duplicates().values.tolist()
R = SETTINGS["genes_with_PLP"]; E = R[~R.recent]; L = R[R.recent]
def job(k):
    tool, rc, direction = k; sign = 1 if direction == "+" else -1; out = []
    rng = np.random.default_rng(11)
    Ev = E[E[rc].notna()]; Lv = L[L[rc].notna()]
    nb = int((Ev.y == 0).sum()); npl = int((Lv.y == 1).sum())
    s2, y2, _, _ = prep(Lv, rc, sign)
    for r in range(NREP):
        Lm = pd.concat([Lv[Lv.y == 1], Lv[Lv.y == 0].sample(n=min(nb, int((Lv.y == 0).sum())), random_state=int(rng.integers(1e9)))])
        Em = pd.concat([Ev[Ev.y == 1].sample(n=min(npl, int((Ev.y == 1).sum())), random_state=int(rng.integers(1e9))), Ev[Ev.y == 0]])
        for nm, G in [("le2022_matched", Em), ("2023_2026_matched", Lm)]:
            s, y, gi, n = prep(G, rc, sign); grid, lo, hi = windows(s); th = thresholds(grid, pooled_lr(lo, hi, y, gi, np.ones(ng)))
            out.append(dict(tool=tool, rep=r, thresholds_from=nm, n=n, n_plp=int(y.sum()), **{kk: (sign * v if np.isfinite(v) else np.nan) for kk, v in th.items()}, **shares(s2, y2, th)))
    return out
with Pool(min(11, len(tools))) as pool: res = pool.map(job, tools)
Q = pd.DataFrame([r for x in res for r in x]); Q.to_csv(f"{OUT}/O2_threshold_sizematched.tsv", sep="\t", index=False)
cols = ["PLP_pp3_any", "PLP_pp3_modplus", "PLP_pp3_strong", "BLB_pp3_any", "BLB_pp3_modplus", "BLB_bp4_any", "BLB_bp4_modplus"]
per_tool = Q.groupby(["tool", "thresholds_from"])[cols + ["PP3_Supporting", "PP3_Moderate", "PP3_Strong"]].median()
print(per_tool.round(3).to_string())
print("median over tools"); print(per_tool.groupby("thresholds_from")[cols].median().round(3).to_string())
q = Q[Q.tool == "REVEL"].groupby("thresholds_from")[["PP3_Supporting", "PP3_Moderate", "PP3_Strong", "BP4_Supporting", "BP4_Moderate"]].quantile([0.05, 0.5, 0.95]).unstack()
print(q.round(3).to_string())
lo_ = Q.pivot_table(index=["tool", "rep"], columns="thresholds_from", values="PLP_pp3_modplus"); print("reps with recent > early, per tool:"); print((lo_["2023_2026_matched"] > lo_["le2022_matched"]).groupby("tool").mean().to_string())
print("DONE")
