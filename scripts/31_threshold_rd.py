#!/usr/bin/env python
"""Threshold fingerprints in ClinVar.

Regression-discontinuity estimates at the ClinGen PP3/BP4 score thresholds:
(A) the probability that a submission (SCV) cites PP3 or BP4, as a function of the
    variant's REVEL (or AlphaMissense) score, by the year the SCV was last evaluated;
(B) the probability that a rare labelled variant is P/LP rather than B/LB, as a
    function of its REVEL score, by label vintage.
Local linear fits with a triangular kernel on each side of the cut point; intervals
from gene-cluster bootstrap; placebo distribution from a scan over cut points that
lie at least one bandwidth away from any published threshold.
Outputs results/r2/P_rd_*.tsv
"""
import sys, os; sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
R2 = f"{P}/results/r2"
H = float(os.environ.get("RD_H", 0.05)); B = int(os.environ.get("RD_B", 500)); SEED = 2026

# ------------------------------------------------------------------ data
S = pd.read_csv(f"{P}/results/scv_textmining.tsv.gz", sep="\t",
                usecols=["VariationID", "SCV", "Submitter", "ClinSig", "year", "PP3", "BP4", "REVEL", "AlphaMissense"])
D = load(eval_only=False)
V = D[["VariationID", "Gene", "REVEL_raw", "AlphaMissense_raw"]].drop_duplicates("VariationID")
V["VariationID"] = V.VariationID.astype(str); S["VariationID"] = S.VariationID.astype(str)
X = S.merge(V, on="VariationID", how="inner")
cs = X.ClinSig.str.lower()
X["cls"] = np.select([cs.str.contains("conflict"), cs.str.contains("uncertain"), cs.str.contains("pathogenic"), cs.str.contains("benign")],
                     ["CONF", "VUS", "PLP", "BLB"], "OTHER")
X = X[X.cls.isin(["PLP", "VUS", "BLB"]) & X.year.notna()].copy()
X["year"] = X.year.astype(int)
X["Gene"] = X.Gene.fillna("NA")

Tt = pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz", sep="\t", usecols=["hgvs38", "rare", "period"])
L = load()
L = L.merge(Tt, on="hgvs38", how="left")
L = L[L.rare.fillna(False).astype(bool)].copy()

genes = pd.Index(sorted(set(X.Gene) | set(L.Gene.fillna("NA"))))
X["gi"] = genes.get_indexer(X.Gene); L["gi"] = genes.get_indexer(L.Gene.fillna("NA"))
NG = len(genes)

# ------------------------------------------------------------------ estimator
def side_fit(u, y, w):
    S0 = w.sum(); S1 = (w * u).sum(); S2 = (w * u * u).sum(); T0 = (w * y).sum(); T1 = (w * u * y).sum()
    den = S0 * S2 - S1 * S1
    return (S2 * T0 - S1 * T1) / den if den > 0 else np.nan

class Cut:
    """Pre-selects the observations within one bandwidth of c; x >= c is the evidence side."""
    def __init__(self, x, y, g, c, h):
        m = (x >= c - h) & (x < c + h)
        self.u = x[m] - c; self.y = y[m].astype(float); self.g = g[m]
        self.k = 1 - np.abs(self.u) / h
        self.r = self.u >= 0
        self.nL = int((~self.r).sum()); self.nR = int(self.r.sum())
    def jump(self, gw=None):
        w = self.k if gw is None else self.k * gw[self.g]
        aR = side_fit(self.u[self.r], self.y[self.r], w[self.r])
        aL = side_fit(self.u[~self.r], self.y[~self.r], w[~self.r])
        return aR - aL, aL, aR

def gene_weights(rng):
    return np.bincount(rng.integers(0, NG, NG), minlength=NG).astype(float)

# ------------------------------------------------------------------ analysis sets
# orientation: for BP4 the score is negated so that the evidence side is always x >= c
REVEL_PP3 = {"PP3_Supporting": 0.644, "PP3_Moderate": 0.773, "PP3_Strong": 0.932}
REVEL_BP4 = {"BP4_Supporting": 0.290, "BP4_Moderate": 0.183, "BP4_Strong": 0.016}
REVEL_OLD = {"legacy_0.5": 0.5, "legacy_0.75": 0.75}
AM_PP3 = {"PP3_Supporting": 0.792, "PP3_Moderate": 0.906, "PP3_Strong": 0.990}
AM_BP4 = {"BP4_Supporting": 0.169, "BP4_Moderate": 0.099}

def scv_period(y):
    return "<=2021" if y <= 2021 else ("2022" if y == 2022 else "2023-2026")
def scv_period_am(y):
    return "<=2023" if y <= 2023 else ("2024" if y == 2024 else "2025-2026")
X["per"] = X.year.map(scv_period); X["per_am"] = X.year.map(scv_period_am)
L["per"] = np.where(L.period == "2023-2026", "2023-2026", "<=2022")

specs = []   # (level, score, side, sample df, x, y, periodcol)
for side, cls_ in [("PP3", ["PLP", "VUS"]), ("BP4", ["BLB", "VUS"])]:
    for score, col, pc in [("REVEL", "REVEL_raw", "per"), ("AlphaMissense", "AlphaMissense_raw", "per_am")]:
        d = X[X.cls.isin(cls_) & X[col].notna()]
        x = d[col].values.astype(float); x = x if side == "PP3" else -x
        specs.append(("scv", score, side, d, x, d[side].values.astype(float), pc))
d = L[L.REVEL_raw.notna()]
specs.append(("label", "REVEL", "PP3", d, d.REVEL_raw.values.astype(float), d.y.values.astype(float), "per"))
specs.append(("label", "REVEL", "BP4", d, -d.REVEL_raw.values.astype(float), 1 - d.y.values.astype(float), "per"))
dv = L[L.REVEL_raw.notna() & L.prior_VUS.fillna(False).astype(bool)]
specs.append(("label_priorVUS", "REVEL", "PP3", dv, dv.REVEL_raw.values.astype(float), dv.y.values.astype(float), "per"))
specs.append(("label_priorVUS", "REVEL", "BP4", dv, -dv.REVEL_raw.values.astype(float), 1 - dv.y.values.astype(float), "per"))

def cuts_for(score, side):
    if score == "REVEL":
        base = REVEL_PP3 if side == "PP3" else REVEL_BP4
        extra = REVEL_OLD if side == "PP3" else {}
    else:
        base = AM_PP3 if side == "PP3" else AM_BP4; extra = {}
    out = {k: (v if side == "PP3" else -v) for k, v in {**base, **extra}.items()}
    return out
def all_thresholds(score, side):
    return np.array(list(cuts_for(score, side).values()))

rows, scan, bins = [], [], []
rng = np.random.default_rng(SEED)
for level, score, side, d, x, y, pc in specs:
    g = d.gi.values; per = d[pc].values
    for p in sorted(set(per)):
        m = per == p
        xs, ys, gs = x[m], y[m], g[m]
        # binned means for plotting
        lo_, hi_ = (0, 1) if side == "PP3" else (-1, 0)
        e = np.round(np.arange(lo_, hi_ + 1e-9, 0.01), 3)
        bi = np.digitize(xs, e) - 1
        for k in range(len(e) - 1):
            mm = bi == k
            if mm.sum() >= 20:
                bins.append(dict(level=level, score=score, side=side, period=p, x_lo=e[k], x_hi=e[k + 1], n=int(mm.sum()), rate=ys[mm].mean()))
        # scan of cut points
        grid = np.round(np.arange(lo_ + H, hi_ - H + 1e-9, 0.005), 3)
        thr = all_thresholds(score, side)
        for c in grid:
            ct = Cut(xs, ys, gs, c, H)
            if ct.nL < 30 or ct.nR < 30: continue
            j, aL, aR = ct.jump()
            scan.append(dict(level=level, score=score, side=side, period=p, c=c, jump=j, nL=ct.nL, nR=ct.nR,
                             placebo=bool(np.min(np.abs(thr - c)) > H)))
        sc = pd.DataFrame([s for s in scan if s["level"] == level and s["score"] == score and s["side"] == side and s["period"] == p])
        null = sc[sc.placebo].jump.values if len(sc) else np.array([])
        # named cut points with bootstrap
        for name, c in cuts_for(score, side).items():
            for h in ([H, 0.03, 0.08] if name in ("PP3_Supporting", "BP4_Supporting") else [H]):
                ct = Cut(xs, ys, gs, c, h)
                if ct.nL < 30 or ct.nR < 30: continue
                j, aL, aR = ct.jump()
                bs = np.array([ct.jump(gene_weights(rng))[0] for _ in range(B)])
                pp = (np.sum(np.abs(null) >= abs(j)) + 1) / (len(null) + 1) if (h == H and len(null)) else np.nan
                rows.append(dict(level=level, score=score, side=side, period=p, cut=name, c=abs(c), h=h, nL=ct.nL, nR=ct.nR,
                                 rate_below=aL, rate_above=aR, jump=j, lo=np.nanpercentile(bs, 2.5), hi=np.nanpercentile(bs, 97.5),
                                 p_placebo=pp, n_placebo=len(null)))
R = pd.DataFrame(rows); R.to_csv(f"{R2}/P_rd_jumps.tsv", sep="\t", index=False)
pd.DataFrame(scan).to_csv(f"{R2}/P_rd_scan.tsv", sep="\t", index=False)
pd.DataFrame(bins).to_csv(f"{R2}/P_rd_bins.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ difference in discontinuities (late minus early), paired bootstrap
diff = []
pairs = {"per": ("<=2021", "2023-2026"), "per_am": ("<=2023", "2025-2026")}
for level, score, side, d, x, y, pc in specs:
    e0, e1 = pairs[pc] if level == "scv" else ("<=2022", "2023-2026")
    g = d.gi.values; per = d[pc].values
    m0, m1 = per == e0, per == e1
    for name, c in cuts_for(score, side).items():
        c0 = Cut(x[m0], y[m0], g[m0], c, H); c1 = Cut(x[m1], y[m1], g[m1], c, H)
        if min(c0.nL, c0.nR, c1.nL, c1.nR) < 30: continue
        j0 = c0.jump()[0]; j1 = c1.jump()[0]; bs = []
        for _ in range(B):
            gw = gene_weights(rng); bs.append(c1.jump(gw)[0] - c0.jump(gw)[0])
        bs = np.array(bs)
        diff.append(dict(level=level, score=score, side=side, cut=name, c=abs(c), early=e0, late=e1, jump_early=j0, jump_late=j1,
                         diff=j1 - j0, lo=np.nanpercentile(bs, 2.5), hi=np.nanpercentile(bs, 97.5)))
pd.DataFrame(diff).to_csv(f"{R2}/P_rd_diff.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ leave-one-submitter-out (SCV level, REVEL supporting cuts, late period)
lo_rows = []
for side, cls_, c in [("PP3", ["PLP", "VUS"], 0.644), ("BP4", ["BLB", "VUS"], -0.290)]:
    d = X[X.cls.isin(cls_) & X.REVEL_raw.notna() & (X.per == "2023-2026")]
    x = d.REVEL_raw.values.astype(float); x = x if side == "PP3" else -x
    top = d[d[side] == 1].Submitter.value_counts().head(10).index.tolist()
    for sub in ["(none)"] + top:
        m = (d.Submitter != sub).values
        ct = Cut(x[m], d[side].values[m].astype(float), d.gi.values[m], c, H)
        j = ct.jump()[0]
        lo_rows.append(dict(side=side, dropped=sub, n_dropped_cited=int(((d.Submitter == sub) & (d[side] == 1)).sum()), jump=j, nL=ct.nL, nR=ct.nR))
    # submitter-specific jumps for the top citers
    for sub in top:
        m = (d.Submitter == sub).values
        ct = Cut(x[m], d[side].values[m].astype(float), d.gi.values[m], c, H)
        if ct.nL >= 30 and ct.nR >= 30:
            lo_rows.append(dict(side=side, dropped="only:" + sub, n_dropped_cited=int((d[m][side] == 1).sum()), jump=ct.jump()[0], nL=ct.nL, nR=ct.nR))
pd.DataFrame(lo_rows).to_csv(f"{R2}/P_rd_loso.tsv", sep="\t", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(R[R.h == H][["level", "score", "side", "period", "cut", "nL", "nR", "rate_below", "rate_above", "jump", "lo", "hi", "p_placebo"]].round(4).to_string())
print(pd.DataFrame(diff).round(4).to_string())
print(pd.DataFrame(lo_rows).round(4).to_string())
