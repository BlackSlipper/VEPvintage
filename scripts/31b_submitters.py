#!/usr/bin/env python
"""Submitter-specific discontinuities at REVEL 0.644 (PP3) and 0.290 (BP4); appended to 31_threshold_rd outputs."""
import sys, os; sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)),"31_threshold_rd.py")).read().split("# ------------------------------------------------------------------ analysis sets")[0])
REVEL_LIST = None
X["per"] = np.where(X.year <= 2021, "<=2021", np.where(X.year == 2022, "2022", "2023-2026"))
rng = np.random.default_rng(SEED)
out = []
for side, cls_, c in [("PP3", ["PLP", "VUS"], 0.644), ("BP4", ["BLB", "VUS"], -0.290)]:
    d0 = X[X.cls.isin(cls_) & X.REVEL_raw.notna()]
    x0 = d0.REVEL_raw.values.astype(float); x0 = x0 if side == "PP3" else -x0
    win = np.abs(x0 - c) < H
    for per in ["<=2021", "2023-2026"]:
        mp = (d0.per == per).values
        cited = d0[mp & win & (d0[side] == 1).values].Submitter.value_counts()
        tot_cited = int(cited.sum())
        for sub, ncit in cited.head(25).items():
            m = mp & (d0.Submitter == sub).values
            ct = Cut(x0[m], d0[side].values[m].astype(float), d0.gi.values[m], c, H)
            if ct.nL < 15 or ct.nR < 15:
                out.append(dict(side=side, period=per, submitter=sub, n_cited_window=int(ncit), share_cited_window=ncit / tot_cited, nL=ct.nL, nR=ct.nR)); continue
            j, aL, aR = ct.jump()
            bs = np.array([ct.jump(gene_weights(rng))[0] for _ in range(300)])
            out.append(dict(side=side, period=per, submitter=sub, n_cited_window=int(ncit), share_cited_window=ncit / tot_cited, nL=ct.nL, nR=ct.nR,
                            rate_below=aL, rate_above=aR, jump=j, lo=np.nanpercentile(bs, 2.5), hi=np.nanpercentile(bs, 97.5)))
        out.append(dict(side=side, period=per, submitter="(all submitters)", n_cited_window=tot_cited, share_cited_window=1.0))
O = pd.DataFrame(out); O.to_csv(f"{R2}/P_rd_submitters.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 45)
print(O.round(3).to_string())
# step-like submitters: jump >= 0.3 with lower bound > 0
for side in ["PP3", "BP4"]:
    for per in ["<=2021", "2023-2026"]:
        s = O[(O.side == side) & (O.period == per) & (O.submitter != "(all submitters)")]
        st = s[(s.jump >= 0.3) & (s.lo > 0)]
        print(side, per, "step-like:", len(st), "of", s.jump.notna().sum(), "share of window citations:", round(st.share_cited_window.sum(), 3))
B = pd.read_csv(f"{R2}/P_rd_bins.tsv", sep="\t")
b = B[(B.level == "scv") & (B.score == "REVEL") & (B.side == "PP3") & (B.x_lo >= 0.55) & (B.x_lo <= 0.75)]
print(b.pivot(index="x_lo", columns="period", values="rate").round(3).to_string())
