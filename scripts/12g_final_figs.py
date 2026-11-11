#!/usr/bin/env python
"""Figs 2 and 7: Fig. 2 with panel c on the rare variants scored by all 38 core predictors (the set used in the
text; script 32c), and Fig. 7 with the class-matched and held-out-gene recalibration checks (scripts 30b, 30c, 30d).
usage: 12g_final_figs.py fig2n fig7n"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.dirname(os.path.abspath(__file__))
exec(open(f"{SRC}/12f_revision_figs.py").read().split("# ---------------------------------------------------------------- Fig. 1: the label set")[0])

REPL2 = [
    ('Rb=pd.read_csv(f"{P}/results/robustness.tsv",sep="\\t"); Wg=pd.read_csv(f"{P}/results/review_within_gene.tsv",sep="\\t")',
     'Qc=pd.read_csv(f"{P}/results/r2/Q_fig2c_common.tsv",sep="\\t")'),
    ('pooled=Rb[(Rb.analysis=="rare_vintage_first")&(Rb.tool==t)].set_index("period").auc.reindex(per)',
     'qq=Qc[Qc.tool==t].set_index("period").reindex(per); pooled=qq.pooled'),
    ('within=Wg[(Wg.analysis=="within_gene_vintage")&(Wg.tool==t)].set_index("period").auc.reindex(per)', 'within=qq.within'),
    ('c.plot(range(3),within.values,color=col,lw=1.5,marker="s",ms=3,ls="--",label="Within gene")',
     'c.plot(range(3),within.values,color=col,lw=1.5,marker="s",ms=3,ls="--",label="Within gene")\n'
     '        c.errorbar(np.arange(3),pooled.values,yerr=[pooled.values-qq.pooled_lo.values,qq.pooled_hi.values-pooled.values],fmt="none",ecolor=col,elinewidth=0.7,capsize=1.2)\n'
     '        c.errorbar(np.arange(3),within.values,yerr=[within.values-qq.within_lo.values,qq.within_hi.values-within.values],fmt="none",ecolor=col,elinewidth=0.7,capsize=1.2,alpha=0.6)'),
    ('lo=min(pooled.min(),within.min()); hi=max(pooled.max(),within.max())', 'lo=min(qq.pooled_lo.min(),qq.within_lo.min()); hi=max(qq.pooled_hi.max(),qq.within_hi.max())'),
    # presentation (round 7): display names, axis labels, coloured and separated end labels in b, pale pink band, note in d
    ('fig=plt.figure(figsize=(7.1,6.6)); G0=fig.add_gridspec(2,1,height_ratios=[2.2,1],hspace=0.62)', 'fig=plt.figure(figsize=(7.1,7.95)); G0=fig.add_gridspec(2,1,height_ratios=[2.744,1],hspace=0.7267)'),
    ('bot=G0[1].subgridspec(1,5,wspace=0.62)', 'bot=G0[1].subgridspec(1,6,width_ratios=[1,1,1,0.05,1,1],wspace=0.62)'),
    ('d=fig.add_subplot(bot[3:5])', 'd=fig.add_subplot(bot[4:6])'),
    ('bbox_to_anchor=(0.76,0.415)', 'bbox_to_anchor=(0.76,0.3913),title="Predictor category",title_fontsize=6'),
    ('a.set_yticklabels(W.index,fontsize=5.4)', 'a.set_yticklabels([NM2(t) for t in W.index],fontsize=5.4)'),
    ('cax=a.inset_axes([0.0,-0.2,0.45,0.022])', 'cax=a.inset_axes([0.0,-0.2004,0.45,0.01764])'),
    ('a.set_title("AUROC by label vintage, 38 core predictors",loc="left"); lab(a,"a")', 'a.set_title("AUROC by label vintage, 38 core predictors",loc="left"); lab(a,"a"); a.set_xlabel("Label vintage")'),
    ('ly=repel([W.loc[t].values[-1] for t in rep],0.0075)', 'ly=spread([W.loc[t].values[-1] for t in rep],0.0105)'),
    ('b.text(x[-1]+0.35,yy,t,fontsize=5.8,va="center",color="#0b0b0b")',
     'b.text(x[-1]+0.42,yy,NM2(t),fontsize=5.8,va="center",color=CC[META2.loc[t,"category"]],fontweight="bold"); b.plot([x[-1]+0.1,x[-1]+0.36],[y[-1],yy],color=CC[META2.loc[t,"category"]],lw=0.5)'),
    ('b.axvspan(8.5,12.5,color="#fbfbe3",zorder=0,lw=0)', 'b.axvspan(8.5,12.5,color="#fdecef",zorder=0,lw=0)'),
    ('b.set_title("Representative predictors",loc="left"); lab(b,"b")', 'b.set_title("Representative predictors",loc="left"); lab(b,"b"); b.set_xlabel("Label vintage")'),
    ('c.set_ylabel("AUROC (rare)")', 'c.set_ylabel("AUROC (rare variants)")'),
    ('d.text(i23+0.6,d.get_ylim()[0]+0.004,"rare B/LB %.1fk ', 'd.text(i23+0.1,(max(d.lines[0].get_ydata()[i23-2:i23+1])+min(d.lines[1].get_ydata()[i23-2:i23+1]))/2,"rare B/LB\\n%.1fk '),
    ('%(nb.iloc[i23]/1e3,nb.iloc[i24]/1e3),fontsize=5.4,color=INK2,va="bottom")', '%(nb.iloc[i23]/1e3,nb.iloc[i24]/1e3),fontsize=5.4,color=INK2,va="center",ha="right")'),
    ('d.set_ylabel("Median AUROC (rare)")', 'd.set_ylabel("Median AUROC\\n(rare variants)"); d.set_xlabel("ClinVar release")'),
]
def fig2n():
    src = open(f"{SRC}/12_figures.py").read()
    f2 = "def fig2():" + src.split("def fig2():")[1].split('if __name__=="__main__":')[0]
    for o, n in REPL2:
        assert f2.count(o) == 1, o[:60]; f2 = f2.replace(o, n)
    exec(f2, globals()); globals()["fig2"]()

def fig7n():
    C = pd.read_csv(f"{R2}/O2_threshold_shift_categories.tsv", sep="\t"); SM = pd.read_csv(f"{R2}/O2_threshold_sizematched.tsv", sep="\t")
    O3 = pd.read_csv(f"{R2}/O3_threshold_combined.tsv", sep="\t")
    fig = plt.figure(figsize=(7.1, 6.0)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.62, wspace=0.62)
    # a: local likelihood ratio of REVEL
    a = fig.add_subplot(gs[0, 0])
    D = load(); T = pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz", sep="\t", usecols=["hgvs38", "rare", "period"]); D = D.merge(T, on="hgvs38")
    D["rare"] = D.rare.fillna(False).astype(bool); plp_genes = set(D.loc[D.y == 1, "Gene"].dropna())
    Rr = D[D.rare & D.Gene.isin(plp_genes) & D.REVEL_raw.notna()]
    def curve(G, WIN=100, NG=400):
        o = np.argsort(G.REVEL_raw.values.astype(float), kind="mergesort"); s = G.REVEL_raw.values.astype(float)[o]; y = G.y.values.astype(float)[o]
        n = len(s); grid = np.quantile(s, np.linspace(0.005, 0.995, NG)); step = (s[-1] - s[0]) / 2000; out = []
        cP = np.concatenate([[0], np.cumsum(y)]); cN = np.concatenate([[0], np.cumsum(1 - y)])
        for g in grid:
            h = 0.0
            while True:
                lo_ = np.searchsorted(s, g - h, side="left"); hi_ = np.searchsorted(s, g + h, side="right")
                if hi_ - lo_ >= WIN or (lo_ == 0 and hi_ == n): break
                h += step
            out.append(((cP[hi_] - cP[lo_] + 0.5) / (cP[-1] + 0.5)) / ((cN[hi_] - cN[lo_] + 0.5) / (cN[-1] + 0.5)))
        return grid, np.array(out)
    for G, lb, col in [(Rr[Rr.period != "2023-2026"], "Labels ≤2022", SKY), (Rr[Rr.period == "2023-2026"], "Labels 2023–26", ROSE)]:
        gx, lr = curve(G); a.plot(gx, lr, color=col, lw=1.0, label=lb)
    for v, t in [(2.406, "supporting"), (5.790, "moderate"), (33.53, "strong")]:
        a.axhline(v, color=INK2, lw=0.5, ls=":"); a.axhline(1 / v, color=INK2, lw=0.5, ls=":"); a.text(0.01, v * 1.12, t, fontsize=5.2, color=INK2)
    for xv in [0.644, 0.773, 0.932]: a.axvline(xv, color="#a6a6a6", lw=0.5, ls="--")
    a.set_yscale("log"); a.set_ylim(0.01, 2000); a.set_xlim(0, 1); a.set_xlabel("REVEL score"); a.set_ylabel("Local likelihood ratio")
    hh, ll = a.get_legend_handles_labels(); a.legend(hh + [Line2D([], [], color="#a6a6a6", lw=0.6, ls="--")], ll + ["Published PP3\nthresholds"], loc="upper left", fontsize=5.2, frameon=True, framealpha=1, facecolor="white", edgecolor="#d9d9d9", borderaxespad=0.3, handlelength=1.6)
    a.set_title("REVEL, genes with P/LP labels", loc="left"); lab(a, "a")
    cc = C[(C.setting == "genes_with_PLP") & (C.kind == "pooled")].pivot(index="tool", columns="thresholds_from", values=["PLP_pp3_modplus", "BLB_pp3_any"])
    def dumb(ax, key, xlabel, title, letter, xmax):
        dd = cc[key].sort_values("2023_2026_all"); y = np.arange(len(dd))
        ax.plot(np.vstack([dd.le2022_all * 100, dd["2023_2026_all"] * 100]), np.vstack([y + 0.2, y - 0.2]), color="#cfcfcf", lw=0.8, zorder=1)  # r9: markers dodged, as in d and e
        ax.scatter(dd.le2022_all * 100, y + 0.2, s=14, color=SKY, zorder=2, label="Labels ≤2022", edgecolor=EDGE, lw=0.35)
        ax.scatter(dd["2023_2026_all"] * 100, y - 0.2, s=14, color=ROSE, zorder=3, label="Labels 2023–26", edgecolor=EDGE, lw=0.35)
        ax.set_yticks(y); ax.set_yticklabels([NM(t) for t in dd.index], fontsize=5.6); ax.set_ylim(-0.5, len(dd) - 0.5); ax.set_xlim(0, xmax); ax.set_xlabel(xlabel); grid_x(ax)
        if letter == "b": xt_ = ax.get_xticks(); ax.set_xlim(-0.04 * xmax, xmax); ax.set_xticks(xt_)  # r9: small margin below 0 so markers near 0% clear the spine; ticks as before
        ax.legend(loc="lower right", fontsize=5.0, handletextpad=0.2, borderaxespad=0.2, title="Thresholds from", title_fontsize=5.0); ax.set_title(title, loc="left"); lab(ax, letter)
    dumb(fig.add_subplot(gs[0, 1]), "PLP_pp3_modplus", "% of rare 2023–26 P/LP\nat PP3 moderate or stronger", "P/LP at PP3 ≥ moderate", "b", 100)
    dumb(fig.add_subplot(gs[0, 2]), "BLB_pp3_any", "% of rare 2023–26 B/LB\nreaching any PP3 level", "B/LB reaching any PP3", "c", max(5, cc["BLB_pp3_any"].max().max() * 100 * 1.15))
    # d, e: four designs, median over predictors (IQR across predictors)
    def per_tool(df, src_col, key):
        return df.groupby(["tool", src_col])[key].median().unstack()
    designs = []
    for key in ["PLP_pp3_modplus", "BLB_pp3_any"]:
        c0 = C[(C.setting == "genes_with_PLP") & (C.kind == "pooled")].pivot(index="tool", columns="thresholds_from", values=key).rename(columns={"le2022_all": "le", "2023_2026_all": "re"})
        c1 = per_tool(SM, "thresholds_from", key).rename(columns={"le2022_matched": "le", "2023_2026_matched": "re"})
        o3 = O3[O3.eval_set == "eval_2023_2026"]
        c2 = per_tool(o3[o3.design == "split"], "thresholds_from", key).rename(columns={"le2022": "le", "2023_2026": "re"})
        c3 = per_tool(o3[o3.design == "split_matched"], "thresholds_from", key).rename(columns={"le2022": "le", "2023_2026": "re"})
        designs.append([("All genes", c0), ("Matched\ncounts", c1), ("Held-out\ngenes", c2), ("Held-out genes,\nmatched counts", c3)])
    rows = []
    for (key, title, letter, col), dl in zip([("PLP_pp3_modplus", "P/LP at PP3 ≥ moderate", "d", gs[1, 0]), ("BLB_pp3_any", "B/LB reaching any PP3", "e", gs[1, 1])], designs):
        ax = fig.add_subplot(col); y = np.arange(len(dl))
        for i, (nm, tb) in enumerate(dl):
            for src_, colr, dy in [("le", SKY, -0.12), ("re", ROSE, 0.12)]:
                v = tb[src_] * 100; m = v.median(); q1, q3 = v.quantile([0.25, 0.75])
                ax.hlines(i + dy, q1, q3, color=colr, lw=1.0, zorder=2); ax.scatter(m, i + dy, s=18, color=colr, edgecolor=EDGE, lw=0.35, zorder=3)
                rows.append(dict(metric=key, design=nm.replace("\n", " "), thresholds_from=src_, median=m, q1=q1, q3=q3))
            ax.plot([tb["le"].median() * 100, tb["re"].median() * 100], [i - 0.12, i + 0.12], color="#bdbdbd", lw=0.6, zorder=1)
        ax.set_yticks(y); ax.set_yticklabels([nm for nm, _ in dl], fontsize=5.6); ax.set_ylim(len(dl) - 0.4, -0.6); grid_x(ax)
        ax.set_xlabel("Median over 11 predictors (%; bar, IQR)"); ax.set_title(title, loc="left"); lab(ax, letter)
        if letter == "e":
            ax.legend(handles=[Line2D([], [], marker="o", ls="", color=SKY, mec=EDGE, mew=0.35, ms=4, label="Labels ≤2022"), Line2D([], [], marker="o", ls="", color=ROSE, mec=EDGE, mew=0.35, ms=4, label="Labels 2023–26")],
                      loc="upper center", bbox_to_anchor=(-0.31, -0.235), ncol=2, fontsize=5.0, title="Thresholds from", title_fontsize=5.0, handletextpad=0.1, columnspacing=0.8, borderaxespad=0)  # r9: below d and e
    # f: achieved likelihood ratio on held-out genes (matched design), per predictor, against the target
    f = fig.add_subplot(gs[1, 2]); o = O3[(O3.design == "split_matched") & (O3.eval_set == "eval_2023_2026")]
    med = o.groupby(["tool", "thresholds_from"])[["PP3_Supporting_LR", "PP3_Moderate_LR", "PP3_Strong_LR"]].median()
    rng = np.random.default_rng(3)
    for k, (key, tgt) in enumerate([("PP3_Supporting_LR", 2.406), ("PP3_Moderate_LR", 5.790), ("PP3_Strong_LR", 33.53)]):
        f.hlines(tgt, k - 0.42, k + 0.42, color=INK2, lw=0.9, ls="--", zorder=1)
        for src_, colr, dx in [("le2022", SKY, -0.16), ("2023_2026", ROSE, 0.16)]:
            v = med.xs(src_, level="thresholds_from")[key].dropna()
            f.scatter(k + dx + rng.uniform(-0.06, 0.06, len(v)), v, s=9, color=colr, edgecolor=EDGE, lw=0.3, zorder=3)
            f.hlines(v.median(), k + dx - 0.1, k + dx + 0.1, color="#333333", lw=1.2, zorder=4)
    f.set_yscale("log"); f.set_xticks(range(3)); f.set_xticklabels(["Supporting", "Moderate", "Strong"], fontsize=6); f.set_xlim(-0.5, 2.5)
    f.set_ylabel("Likelihood ratio in the interval"); grid_y(f)
    f.legend(handles=[Line2D([], [], color=INK2, lw=0.9, ls="--", label="Target (ClinGen)"), Line2D([], [], color="#333333", lw=1.2, label="Median predictor")], loc="upper left", fontsize=5.0, borderaxespad=0.2)
    f.add_artist(f.get_legend())
    f.legend(handles=[Line2D([], [], marker="o", ls="", color=SKY, mec=EDGE, mew=0.3, ms=3.2, label="Labels ≤2022"), Line2D([], [], marker="o", ls="", color=ROSE, mec=EDGE, mew=0.3, ms=3.2, label="Labels 2023–26")],
             loc="upper center", bbox_to_anchor=(0.5, -0.235), ncol=2, fontsize=5.0, title="Thresholds from", title_fontsize=5.0, handletextpad=0.1, columnspacing=0.8, borderaxespad=0)
    f.set_title("Achieved on held-out genes", loc="left"); lab(f, "f")
    save(fig, "Fig7_recalibration")
    pd.DataFrame(rows).to_csv(f"{R2}/Q_fig7de_designs.tsv", sep="\t", index=False)

if __name__ == "__main__":
    for w in sys.argv[1:]: globals()[w](); print("done", w)
