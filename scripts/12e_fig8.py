#!/usr/bin/env python
"""Figure 8: fingerprints of the ClinGen REVEL thresholds in submissions and labels (script 31 outputs)."""
import sys; sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)),"12_figures.py")).read().split("def fig1():")[0])
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)),"31_threshold_rd.py")).read().split("# ------------------------------------------------------------------ analysis sets")[0])
R2 = f"{P}/results/r2"
from matplotlib.lines import Line2D
Bn = pd.read_csv(f"{R2}/P_rd_bins.tsv", sep="\t"); Sc = pd.read_csv(f"{R2}/P_rd_scan.tsv", sep="\t")
Su = pd.read_csv(f"{R2}/P_rd_submitters.tsv", sep="\t")
X["per"] = np.where(X.year <= 2021, "<=2021", np.where(X.year == 2022, "2022", "2023-2026"))
L["per"] = np.where(L.period == "2023-2026", "2023-2026", "<=2022")

def fitline(x, y, c, h, right):
    m = (x >= c) & (x < c + h) if right else (x >= c - h) & (x < c)
    u = x[m] - c; w = 1 - np.abs(u) / h; yy = y[m]
    S0, S1, S2 = w.sum(), (w * u).sum(), (w * u * u).sum(); T0, T1 = (w * yy).sum(), (w * u * yy).sum()
    den = S0 * S2 - S1 * S1; a = (S2 * T0 - S1 * T1) / den; b = (S0 * T1 - S1 * T0) / den
    uu = np.linspace(0, h, 20) if right else np.linspace(-h, 0, 20)
    return c + uu, a + b * uu

PERS = [("<=2021", "Up to 2021", SKY), ("2023-2026", "2023–2026", ROSE)]
fig = plt.figure(figsize=(7.1, 5.0)); gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.42)

def rate_panel(ax, level, side, pers, xlim, cuts, main, ylab, title):
    flip = side == "BP4" and level == "scv"
    for p, nm, col in pers:
        b = Bn[(Bn.level == level) & (Bn.score == "REVEL") & (Bn.side == side) & (Bn.period == p)].copy()
        xm = (b.x_lo + b.x_hi) / 2; xm = -xm if flip else xm
        k = (xm >= xlim[0]) & (xm <= xlim[1])
        ax.scatter(xm[k], b.rate[k] * 100, s=7, color=col, edgecolor=EDGE, lw=0.3, zorder=2, label=nm)
        if level == "scv":
            d = X[X.cls.isin(["PLP", "VUS"] if side == "PP3" else ["BLB", "VUS"]) & X.REVEL_raw.notna() & (X.per == p)]
            x = d.REVEL_raw.values.astype(float); y = d[side].values.astype(float)
        else:
            d = L[L.REVEL_raw.notna() & (L.per == p)]; x = d.REVEL_raw.values.astype(float); y = d.y.values.astype(float)
        for right in (False, True):
            if side == "BP4" and level == "scv":
                xx, yy = fitline(-x, y, -main, H, right); xx = -xx
            else:
                xx, yy = fitline(x, y, main, H, right)
            ax.plot(xx, yy * 100, color=col, lw=1.4, zorder=3)
    for c in cuts: ax.axvline(c, color="#9e9e9e", lw=0.6, ls="--", zorder=1)
    ax.set_xlim(*xlim); ax.set_xlabel("REVEL score"); ax.set_ylabel(ylab); ax.set_title(title, loc="left")
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)

a = fig.add_subplot(gs[0, 0])
rate_panel(a, "scv", "PP3", PERS, (0.45, 0.85), [0.644, 0.773], 0.644, "Submission records\nciting PP3 (%)", "PP3 citation")
a.set_ylim(0, 14); a.legend(loc="upper left", fontsize=5.4, handletextpad=0.1, borderaxespad=0.3, title="Submission records\nevaluated", title_fontsize=5.4, frameon=True, framealpha=1, facecolor="white", edgecolor="#d9d9d9")
a.text(0.648, 0.3, "0.644", fontsize=5.4, color=INK2, va="bottom"); a.text(0.777, 0.3, "0.773", fontsize=5.4, color=INK2, va="bottom"); lab(a, "a")
b = fig.add_subplot(gs[0, 1])
rate_panel(b, "scv", "BP4", PERS, (0.10, 0.45), [0.183, 0.290], 0.290, "Submission records\nciting BP4 (%)", "BP4 citation")
b.set_ylim(0, 9); b.text(0.294, 8.7, "0.290", fontsize=5.4, color=INK2, va="top"); b.text(0.187, 8.7, "0.183", fontsize=5.4, color=INK2, va="top"); lab(b, "b")
c = fig.add_subplot(gs[0, 2])
rate_panel(c, "label", "PP3", [("<=2022", "Labels up to 2022", SKY), ("2023-2026", "Labels 2023–2026", ROSE)], (0.45, 0.85), [0.644, 0.773], 0.644,
           "Rare labels that are P/LP (%)", "Aggregate labels")
c.set_ylim(0, 105); c.legend(loc="lower right", fontsize=5.4, handletextpad=0.1, borderaxespad=0.3, frameon=True, framealpha=1, facecolor="white", edgecolor="#d9d9d9"); lab(c, "c")

def scan_panel(ax, side, xlim, cuts, title):
    for p, nm, col in PERS:
        s = Sc[(Sc.level == "scv") & (Sc.score == "REVEL") & (Sc.side == side) & (Sc.period == p)].copy()
        s["cc"] = s.c.abs(); s = s.sort_values("cc"); k = (s.cc >= xlim[0]) & (s.cc <= xlim[1])
        q = s[s.placebo].jump.quantile([0.025, 0.975]).values * 100
        ax.axhspan(q[0], q[1], color=col, alpha=0.16, lw=0, zorder=0)
        for qv in q: ax.axhline(qv, color=col, lw=0.9, ls=":", zorder=1)  # band limits in the period colour
        ax.plot(s.cc[k], s.jump[k] * 100, color=col, lw=1.1, zorder=2, label=nm)
    for cc in cuts: ax.axvline(cc, color="#9e9e9e", lw=0.6, ls="--", zorder=1)
    ax.axhline(0, color=INK2, lw=0.5)
    ax.set_xlim(*xlim); ax.set_xlabel("Cut point on REVEL"); ax.set_ylabel("Estimated step (percentage points)"); ax.set_title(title, loc="left")
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
d = fig.add_subplot(gs[1, 0]); scan_panel(d, "PP3", (0.40, 0.95), [0.644, 0.773, 0.932], "Scan of cut points, PP3"); lab(d, "d")
e = fig.add_subplot(gs[1, 1]); scan_panel(e, "BP4", (0.05, 0.50), [0.183, 0.290], "Scan of cut points, BP4"); lab(e, "e")

f = fig.add_subplot(gs[1, 2])
s = Su[(Su.period == "2023-2026") & Su.jump.notna()].copy()
order = s.sort_values(["side", "n_cited_window"], ascending=[False, False])
names = {}; letters = iter("CDEFGHIJKLMNOP")
for sub in order.submitter:
    if sub not in names: names[sub] = next(letters)
y = 0; yt = []; yl = []
PTC = "#5a5a5a"; CIC = "#9a9a9a"  # one colour for all submitters; PP3 circles, BP4 squares
for side, col, thr, mk in [("PP3", PTC, "0.644", "o"), ("BP4", PTC, "0.290", "s")]:
    g = order[order.side == side]
    for _, r in g.iterrows():
        f.hlines(y, max(r.lo * 100, -105), min(r.hi * 100, 105), color=CIC, lw=0.8, alpha=0.9)
        if r.lo * 100 < -105: f.plot(-105, y, marker="<", ms=3.2, color=CIC, mec=CIC, clip_on=False, zorder=2)
        if r.hi * 100 > 105: f.plot(105, y, marker=">", ms=3.2, color=CIC, mec=CIC, clip_on=False, zorder=2)
        f.scatter(r.jump * 100, y, s=8 + 1.2 * np.sqrt(r.n_cited_window), color=col, marker=mk, edgecolor=EDGE, lw=0.35, zorder=3)
        yt.append(y); yl.append(f"{names[r.submitter]} ({side})"); y += 1
    y += 0.6
f.axvline(0, color=INK2, lw=0.5); f.set_xlim(-110, 110); f.set_yticks(yt); f.set_yticklabels(yl, fontsize=5.4); f.set_ylim(y - 0.4, -0.8)
f.set_xlabel("Step at the threshold (percentage points)"); f.set_title("Single submitters, 2023–2026", loc="left")
f.xaxis.grid(True, color=GRID, lw=0.5); f.set_axisbelow(True); lab(f, "f")
ks = lambda n: np.sqrt(8 + 1.2 * np.sqrt(n))  # legend marker size (pt) for the scatter area used above
f.add_artist(f.legend(handles=[Line2D([], [], marker="o", ls="", color=PTC, mec=EDGE, mew=0.35, ms=ks(50), label="PP3 (0.644)"), Line2D([], [], marker="s", ls="", color=PTC, mec=EDGE, mew=0.35, ms=ks(50), label="BP4 (0.290)"),
                               Line2D([], [], marker="<", ls="", color=CIC, mec=CIC, ms=3.2, label="CI truncated")], loc="upper center", bbox_to_anchor=(0.45, -0.235), ncol=3, fontsize=5.2, handletextpad=0.1, columnspacing=0.7, borderaxespad=0))
f.legend(handles=[Line2D([], [], marker="o", ls="", color=PTC, mec=EDGE, mew=0.35, ms=ks(n), label=format(n, ",")) for n in (20, 100, 1000)], loc="upper center", bbox_to_anchor=(0.45, -0.345), ncol=3, fontsize=5.2,
         title="Citing records in the window", title_fontsize=5.2, handletextpad=0.2, columnspacing=1.0, borderaxespad=0)
pd.Series(names).rename("letter").to_csv(f"{R2}/P_rd_submitter_letters.tsv", sep="\t")
save(fig, "Fig8_threshold_fingerprints"); print("saved")
