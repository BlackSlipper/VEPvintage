#!/usr/bin/env python
"""Figures for the revised manuscript (order: label set, performance, composition, provenance x2, calibration x2, fingerprints).
usage: 12f_revision_figs.py fig1n fig3n fig4n fig5n fig6n fig7n
Fig. 2 is unchanged (12_figures.py fig2); Fig. 8 is 12e_fig7.py saved under its new name."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.dirname(os.path.abspath(__file__))
exec(open(f"{SRC}/12_figures.py").read().split("def fig1():")[0])
from matplotlib.lines import Line2D
R2 = f"{P}/results/r2"
PER = ["<=2019", "2020-2022", "2023-2026"]; PERL = ["≤2019", "2020–22", "2023–26"]
NM = lambda t: t.replace("100way", "").replace("PolyPhen2_", "PolyPhen2-").replace("_noAF", "")
NM2 = lambda t: t.replace("100way", "").replace("PolyPhen2_", "PolyPhen2-")  # display names that keep BayesDel_addAF/_noAF apart
def spread(ys, mind):
    """Push label positions apart symmetrically until neighbours are >= mind apart (order kept)."""
    y = np.array(ys, float).copy(); o = np.argsort(y)
    for _ in range(500):
        moved = False
        for i0, i1 in zip(o[:-1], o[1:]):
            gap = y[i1] - y[i0]
            if gap < mind - 1e-9: y[i0] -= (mind - gap) / 2; y[i1] += (mind - gap) / 2; moved = True
        if not moved: break
    return y
DARK = "#1a1a1a"; OPEN_GREY = "#8c8c8c"; PER_GREY = "#8c8c8c"; BLUEGREY = "#7f9db5"
def grid_y(ax): ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
def grid_x(ax): ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)

# ---------------------------------------------------------------- Fig. 1: the label set
def fig1n():
    T = pd.read_csv(f"{P}/results/snapshot_class_counts.tsv", sep="\t"); T["date"] = pd.to_datetime(T.Snapshot.replace("current", "2026-11"))
    Tr = pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz", sep="\t")
    A = pd.read_csv(f"{P}/results/variant_level_annotated.tsv.gz", sep="\t", usecols=["hgvs38", "conc_submitters"])
    Tr = Tr.merge(A.drop_duplicates("hgvs38"), on="hgvs38", how="left"); Tr["rare"] = Tr.rare.fillna(False).astype(bool)
    Tr["vbin"] = np.where(Tr.vintage_first == "2015-02", "<=2015Q1", Tr.vintage_first.str[:4])
    Tr["src"] = np.where(Tr.conc_submitters == SUBA, "A", np.where(Tr.conc_submitters == SUBB, "L", "O"))
    S = pd.read_csv(f"{P}/results/scv_textmining.tsv.gz", sep="\t")
    cs = S.ClinSig.str.lower(); S = S[cs.str.contains("pathogenic|benign") & ~cs.str.contains("conflict|uncertain")]
    fig, ax = plt.subplots(2, 2, figsize=(7.1, 5.2), gridspec_kw={"width_ratios": [1, 1.25], "hspace": 0.78, "wspace": 0.32}); ax = ax.ravel()
    a = ax[0]
    for c, col, nm in [("VUS", VUS_C, "VUS"), ("BLB", BLB_C, "B/LB"), ("PLP", PLP_C, "P/LP"), ("CONF", CONF_C, "Conflicting")]:
        a.plot(T.date, T[c].where(T[c] > 0), color=col, lw=1.5, label=nm)  # r9: non-positive counts masked, so each line starts at its first positive value
    ends = np.log10([T[c].iloc[-1] for c in ["VUS", "BLB", "PLP", "CONF"]]); ly = spread(ends, 0.2)
    for nm, yy, col in zip(["VUS", "B/LB", "P/LP", "Conflicting"], ly, [VUS_C, BLB_C, PLP_C, CONF_C]): a.text(T.date.iloc[-1] + pd.Timedelta(days=60), 10 ** yy, nm, va="center", fontsize=6, color=col, fontweight="bold")
    a.set_yscale("log"); a.set_ylabel("Missense SNVs in ClinVar"); a.set_xlim(T.date.iloc[0], T.date.iloc[-1] + pd.Timedelta(days=820))
    a.set_xticks(pd.to_datetime(["2015", "2018", "2021", "2024", "2027"])); a.set_xticklabels(["2015", "2018", "2021", "2024", "2027"])
    a.set_title("ClinVar missense classifications", loc="left"); lab(a, "a"); grid_y(a)
    b = ax[1]; x = np.arange(len(ORDER)); bot = np.zeros(len(ORDER))
    parts = [("P/LP", (Tr.final == "PLP"), FILL["rose"], None), ("B/LB, other sources", (Tr.final == "BLB") & (Tr.src == "O"), FILL["sky"], None),
             ("B/LB, Labcorp only", (Tr.final == "BLB") & (Tr.src == "L"), FILL["mint"], None), ("B/LB, Ambry only", (Tr.final == "BLB") & (Tr.src == "A"), FILL["butter"], "////")]
    for nm, m, col, h in parts:
        v = Tr[m].groupby("vbin").size().reindex(ORDER).fillna(0).values / 1e3
        b.bar(x, v, bottom=bot, color=col, width=0.8, label=nm, edgecolor=EDGE, lw=0.4, hatch=h); bot += v
    b.set_xticks(x); b.set_xticklabels(XL, rotation=60, ha="right"); b.set_ylabel("Variants (×1,000)"); b.set_xlabel("Label vintage (first year with current class)")
    b.set_ylim(0, 40); b.legend(loc="upper left", fontsize=5.8); b.set_title("Evaluation set (n = %s) by label vintage and source" % format(len(Tr), ","), loc="left"); lab(b, "b"); grid_y(b)
    c = ax[2]
    R = Tr[Tr.rare & (Tr.final == "BLB")].copy(); plp_genes = set(Tr[Tr.rare & (Tr.final == "PLP")].Gene.dropna())
    R["nogene"] = ~R.Gene.isin(plp_genes); R["star2"] = R.stars >= 2; R["A"] = R.src == "A"; R["L"] = R.src == "L"
    g = R.groupby("vbin")[["A", "L", "nogene", "star2"]].mean().reindex(ORDER); n = R.groupby("vbin").size().reindex(ORDER).fillna(0)
    for col_, cc, nm in [("A", MUSTARD, "From Ambry only"), ("L", MINT, "From Labcorp only"), ("nogene", ROSE, "In genes without a rare P/LP label"), ("star2", PERI, "Two or more review stars")]:
        c.plot(x, g[col_] * 100, color=cc, lw=1.5, marker="o", ms=2.5, label=nm)
    c.set_xticks(x); c.set_xticklabels(XL, rotation=60, ha="right"); c.set_ylabel("% of rare B/LB labels"); c.set_xlabel("Label vintage"); c.set_ylim(0, 100)
    for i, v in enumerate(n.values): c.text(i, 101, ("%.1fk" % (v / 1e3)) if v >= 1000 else "%d" % v, ha="center", va="bottom", fontsize=6, color=INK2, rotation=90)
    c.legend(loc="upper left", bbox_to_anchor=(0.0, 0.93), fontsize=5.6); c.set_title("Rare B/LB labels by vintage", loc="left", pad=22); lab(c, "c", y=1.21); grid_y(c)
    d = ax[3]
    S["year"] = pd.to_numeric(S.year, errors="coerce"); S = S[(S.year >= 2013) & (S.year <= 2026)]
    S["code"] = S[["PP3", "BP4"]].max(axis=1); S["tool_named"] = S["any_tool"]; S["generic_nonA"] = (S.any_insilico == 1) & (S.Submitter != SUBA)
    gg = S.groupby("year")[["generic_nonA", "tool_named", "code", "REVEL"]].mean()
    for col_, cc, nm in [("generic_nonA", SKY, "Any computational evidence*"), ("code", MINT, "Explicit PP3/BP4"), ("tool_named", PEACH, "Named predictor"), ("REVEL", LAV, "REVEL named")]:
        d.plot(gg.index, gg[col_] * 100, color=cc, lw=1.5, marker="o", ms=2.5, label=nm)
    d.set_ylabel("% of P/LP and B/LB\nsubmission records"); d.set_xlabel("Date last evaluated (year)"); d.set_ylim(0, 45); d.legend(loc="upper left", fontsize=5.8, borderaxespad=0.2)
    d.set_title("Computational evidence cited in submission text", loc="left", pad=22); lab(d, "d", y=1.21); grid_y(d)
    save(fig, "Fig1_label_set")
    g.assign(n=n).to_csv(f"{R2}/Q_fig1c_rare_blb_by_vintage.tsv", sep="\t")

# ---------------------------------------------------------------- Fig. 3: gene composition
def null_panel(a):
    M = pd.read_csv(f"{R2}/M_gene_null.tsv", sep="\t"); rng = np.random.default_rng(1)
    for g, (nm, ttl) in enumerate([("all_common", "All labels"), ("filtered_common", "Provenance-filtered")]):
        for i, per in enumerate(PER):
            x = g * 3.6 + i; sub = M[(M.set == nm) & (M.period == per)]
            pr = sub[~sub.tool.str.startswith("gene_null")].auc.values; nl = sub[sub.tool == "gene_null_loo"]
            a.scatter(x + rng.uniform(-0.18, 0.18, len(pr)), pr, s=6, color="#c4c4c4", zorder=1, lw=0)
            a.hlines(np.median(pr), x - 0.3, x + 0.3, color="#333333", lw=1.2, zorder=3)
            a.errorbar(x, nl.auc.iloc[0], yerr=[[nl.auc.iloc[0] - nl.lo.iloc[0]], [nl.hi.iloc[0] - nl.auc.iloc[0]]], fmt="D", ms=4, color=DARK, ecolor=DARK, elinewidth=0.8, capsize=1.5, zorder=4)
        a.text(g * 3.6 + 1, 1.03, ttl, ha="center", va="top", fontsize=6.2, color=INK2)
    a.set_xticks([0, 1, 2, 3.6, 4.6, 5.6]); a.set_xticklabels(["≤2019", "2020–\n2022", "2023–\n2026"] * 2, fontsize=6); a.set_ylabel("Pooled AUROC (rare variants)")
    a.set_ylim(0.6, 1.04); a.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0]); a.axvline(2.8, color=GRID, lw=0.8); grid_y(a)
    a.legend(handles=[Line2D([], [], marker="o", color="#c4c4c4", ls="", ms=3, label="Core predictor"), Line2D([], [], color="#333333", lw=1.2, label="Median predictor"),
                      Line2D([], [], marker="D", color=DARK, ls="", ms=4, label="Gene-identity null")], loc="upper center", bbox_to_anchor=(0.5, -0.19), ncol=3, fontsize=5.6, handletextpad=0.3, columnspacing=1.2)
    a.set_title("Predictors and the gene-identity null", loc="left"); lab(a, "a")

def fig3n():
    N = pd.read_csv(f"{R2}/N_composition_matched.tsv", sep="\t"); N = N[N["mode"] == "gene_class"].set_index("tool")
    Xc = pd.read_csv(f"{R2}/Q_excluded_cells.tsv", sep="\t"); Xa = pd.read_csv(f"{R2}/Q_excluded_auc.tsv", sep="\t").set_index("label")
    fig = plt.figure(figsize=(7.1, 5.5)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.55, wspace=0.42)
    null_panel(fig.add_subplot(gs[0, 0]))
    # b: three steps
    b = fig.add_subplot(gs[0, 1]); steps = [("d_full", "Full label sets"), ("d_S", "Shared gene × class\ncells only"), ("d_matched_to20", "Shared cells,\nreweighted")]
    xs = np.arange(3)
    for t, r in N.iterrows():
        col = CC.get(META2.loc[t, "category"], "#999") if t in META2.index else "#999"
        b.plot(xs, [r[k] for k, _ in steps], color=col, lw=0.6, alpha=0.55, zorder=1)
    med = [N[k].median() for k, _ in steps]
    b.plot(xs, med, color="#333333", lw=1.8, marker="o", ms=4, zorder=3)
    for i, (k, _) in enumerate(steps):
        npos = int((N[k + "_lo"] > 0).sum()); b.text(i, med[i] + 0.012, "%.3f" % med[i], ha="center", va="bottom", fontsize=6.2, color="#1a1a1a", fontweight="bold")
        b.text(i, -0.0472, "CI > 0: %d/%d" % (npos, len(N)), ha="center", va="bottom", fontsize=5.4, color=INK2)  # r9: just above the -0.05 gridline
    b.axhline(0, color=INK2, lw=0.6, ls=":"); b.set_xticks(xs); b.set_xticklabels([s for _, s in steps], fontsize=5.8); b.set_xlim(-0.35, 2.35); b.set_ylim(-0.065, 0.215)
    b.set_ylabel("AUROC gain, 2023–26 over 2020–22"); grid_y(b)
    b.legend(handles=[Line2D([], [], color="#999", lw=0.8, label="Single predictor"), Line2D([], [], color="#333333", lw=1.8, marker="o", ms=4, label="Median of 38")], loc="upper right", fontsize=5.6)
    b.set_title("Gain at each step of composition matching", loc="left"); lab(b, "b")
    # c: composition of retained and excluded recent labels
    c = fig.add_subplot(gs[1, 0])
    Xc["part"] = np.where(Xc.final == "PLP", "P/LP", np.where(Xc.source == "Ambry only", "B/LB, Ambry only", np.where(Xc.source == "Labcorp only", "B/LB, Labcorp only", "B/LB, other sources")))
    ct = Xc.groupby(["shared_cell", "part"]).n.sum().unstack().fillna(0)
    order = ["P/LP", "B/LB, other sources", "B/LB, Labcorp only", "B/LB, Ambry only"]; cols = [FILL["rose"], FILL["sky"], FILL["mint"], FILL["butter"]]; hat = [None, None, None, "////"]
    tot = ct.sum(axis=1); rows = [(True, "Shared cells"), (False, "Excluded cells")]
    for j, (key, nm) in enumerate(rows):
        left = 0
        for p, col, h in zip(order, cols, hat):
            v = ct.loc[key, p] / 1e3 if p in ct.columns else 0
            c.barh(j, v, left=left, color=col, edgecolor=EDGE, lw=0.4, hatch=h, height=0.62, label=p if j == 0 else None)
            if v > 5.5: c.text(left + v / 2, j, "%.0f%%" % (100 * v * 1e3 / tot[key]), ha="center", va="center", fontsize=5.6, bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none") if h else None)
            left += v
        c.text(left + 0.4, j, "n = %s (%.0f%%)" % (format(int(tot[key]), ","), 100 * tot[key] / tot.sum()), va="center", fontsize=5.8, color=INK2)
    nog = Xc[(~Xc.shared_cell) & (Xc.final == "BLB")]; share_nog = nog[~nog.gene_has_PLP_any].n.sum() / nog.n.sum()
    c.set_yticks([0, 1]); c.set_yticklabels([nm for _, nm in rows]); c.set_ylim(1.95, -1.25); c.set_xlim(0, tot.max() / 1e3 * 1.42)
    c.set_xlabel("Rare 2023–26 labels scored by all core predictors (×1,000)"); grid_x(c)
    c.legend(loc="upper left", ncol=2, fontsize=5.4, handlelength=1.4, columnspacing=0.8, borderaxespad=0.2)
    c.text(1.0, 1.42, "%.0f%% of excluded B/LB labels lie in genes\nwith no P/LP label in the evaluation set" % (100 * share_nog), fontsize=5.6, color=INK2, va="top")
    c.set_title("Which recent labels the shared cells leave out", loc="left"); lab(c, "c")
    # d: AUROC on the subsets
    d = fig.add_subplot(gs[1, 1]); sets = [("early_all", "2020–22, all"), ("recent_shared_cells", "2023–26,\nretained"), ("recent_excluded_cells", "2023–26,\nexcluded cells"), ("recent_all", "2023–26, all")]
    for k, (key, nm) in enumerate(sets):
        r = Xa.loc[key]
        for dx, val, lo, hi, col, mk, lbl in [(-0.2, r.auc_median, r.auc_median_lo, r.auc_median_hi, "#333333", "o", "Median core predictor, pooled"),
                                              (0.0, r.auc_REVEL, r.auc_REVEL_lo, r.auc_REVEL_hi, MUSTARD, "s", "REVEL, pooled")]:
            d.errorbar(k + dx, val, yerr=[[val - lo], [hi - val]], fmt=mk, ms=3.6, color=col, ecolor=col, elinewidth=0.8, capsize=1.5, label=lbl if k == 0 else None, zorder=3)
        d.scatter(k + 0.2, r.within_REVEL, marker="s", s=14, facecolor="white", edgecolor=MUSTARD, lw=1.0, zorder=4, label="REVEL, within gene" if k == 0 else None)
    d.set_xticks(range(4)); d.set_xticklabels([nm for _, nm in sets], fontsize=5.8); d.set_ylabel("AUROC (rare variants)"); d.set_xlim(-0.5, 3.5); grid_y(d)
    lo_, hi_ = d.get_ylim(); d.set_ylim(lo_ - (1.0 - lo_) * 0.3, 1.0); hh_, ll_ = d.get_legend_handles_labels(); oo_ = [ll_.index(s) for s in ["Median core predictor, pooled", "REVEL, pooled", "REVEL, within gene"]]
    d.legend([hh_[i] for i in oo_], [ll_[i] for i in oo_], loc="lower right", fontsize=5.5)
    d.set_title("Accuracy on retained and excluded labels", loc="left"); lab(d, "d")
    save(fig, "Fig3_gene_composition")

# ---------------------------------------------------------------- Fig. 4: sources, provenance groups and codes
exec(open(f"{SRC}/12b_figures.py").read().split("def dumbbell")[1].split("def fig3():")[0].join(["def dumbbell", ""]))
def fig4n():
    E = pd.read_csv(f"{P}/results/evidence_source_contrasts.tsv", sep="\t"); Q = pd.read_csv(f"{R2}/Q_provenance_groups.tsv", sep="\t")
    fig = plt.figure(figsize=(7.1, 7.25)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 2.47], hspace=0.34, wspace=0.6)  # r9: c-e taller, row pitch ~1.17 line heights
    for k, t in enumerate(["REVEL", "phyloP100way"]):
        a = fig.add_subplot(gs[0, k])
        for j, (key, nm, col) in enumerate([("BLB_all", "All B/LB", "#4d4d4d"), ("BLB_excluding_Ambry_Labcorp", "B/LB excl. Ambry, Labcorp", BLUEGREY), ("BLB_only_from:Ambry Genetics", "B/LB from Ambry only", MUSTARD)]):
            ys = []; xs = []
            for i, p in enumerate(PER):
                r = E[(E.contrast == f"rare:{p}:{key}") & (E.tool == t)]
                if len(r): xs.append(i + (j - 1) * 0.12); ys.append(r.auc.iloc[0])
            a.plot(xs, ys, marker="o", ms=3.5, lw=1.2, color=col, label=nm)
        a.set_xticks(range(3)); a.set_xticklabels(PERL); a.set_ylabel("AUROC (rare variants)"); a.set_title(NM(t), loc="left"); grid_y(a); a.set_xlim(-0.4, 2.4)
        if k == 0: lab(a, "a"); lo_, hi_ = a.get_ylim(); a.set_ylim(lo_, hi_ + (hi_ - lo_) * 0.38); a.legend(loc="upper left", fontsize=5.5)
    b = fig.add_subplot(gs[0, 2]); G = ["Explicit computational evidence", "Text without computational evidence", "No free text", "Ambry only (template)"]
    GL = ["Explicit", "Text, no\nmention", "No text", "Ambry\nonly"]
    q = Q[Q.comparison == "all P/LP vs group B/LB"]
    for j, (per, col) in enumerate([("2020-2022", SKY), ("2023-2026", ROSE)]):
        r = q[q.period == per].set_index("label").reindex(G); x = np.arange(4) + (j - 0.5) * 0.24
        ok = r.auc_median.notna()
        b.errorbar(x[ok], r.auc_median[ok], yerr=[(r.auc_median - r.auc_median_lo)[ok], (r.auc_median_hi - r.auc_median)[ok]], fmt="o", ms=3.4, color=col, elinewidth=0.8, capsize=1.2, label=PERL[j + 1])
    b.set_xticks(range(4)); b.set_xticklabels(GL, fontsize=5.4); b.set_ylabel("Median core AUROC"); grid_y(b); b.legend(loc="upper left", fontsize=5.5, ncol=2, title="Label period", title_fontsize=5.5, columnspacing=0.6, handletextpad=0.1); lo_, hi_ = b.get_ylim(); b.set_ylim(lo_, hi_ + (hi_ - lo_) * 0.12)
    b.set_title("B/LB labels by recorded evidence", loc="left"); lab(b, "b")
    c = fig.add_subplot(gs[1, 0]); d = fig.add_subplot(gs[1, 1])
    piv = E[E.contrast.str.startswith("rare:") & ~E.contrast.str.contains(":20|:<=")].pivot(index="tool", columns="contrast", values="auc"); piv.index = [NM2(t) for t in piv.index]
    for ax_, A_, B_, la_, lb_, ttl_, xl_, let_ in [(c, piv["rare:BLB_with_BP4"], piv["rare:BLB_codes_without_BP4"], "B/LB citing BP4", "B/LB citing other codes only", "Benign side: BP4", "AUROC (vs all rare P/LP)", "c"),
                                                   (d, piv["rare:PLP_with_PP3"], piv["rare:PLP_codes_without_PP3"], "P/LP citing PP3", "P/LP citing other codes only", "Pathogenic side: PP3", "AUROC (vs all rare B/LB)", "d")]:
        dumbbell(ax_, A_, B_, la_, lb_, DARK, "white", ttl_)  # code cited: filled black; other codes only: open grey
        ax_.collections[1].set_edgecolor(OPEN_GREY); ax_.collections[1].set_linewidth(0.8)
        ax_.set_title(ttl_, loc="left", pad=plt.rcParams["axes.titlepad"]); hh_, ll_ = ax_.get_legend_handles_labels()
        ax_.legend(hh_[::-1], ll_[::-1], loc="upper center", bbox_to_anchor=(0.4, -0.11), fontsize=5.4, handletextpad=0.2, borderaxespad=0)
        ax_.set_xlabel(xl_); lab(ax_, let_)
    e = fig.add_subplot(gs[1, 2])
    piv0 = E[E.contrast.str.startswith("rare:") & ~E.contrast.str.contains(":20|:<=")].pivot(index="tool", columns="contrast", values="auc")
    dl = pd.DataFrame({"BP4": piv0["rare:BLB_with_BP4"] - piv0["rare:BLB_codes_without_BP4"], "PP3": piv0["rare:PLP_with_PP3"] - piv0["rare:PLP_codes_without_PP3"]})
    dl["cat"] = META2.loc[dl.index, "category"]
    for cn in CATS:
        s = dl[dl.cat == cn]; e.scatter(s.BP4, s.PP3, s=14, color=CC[cn], edgecolor=EDGE, lw=0.35, label=cn)
    lpos = {"ESM1b": (0.068, 0.076, "right"), "REVEL": (0.0255, 0.062, "left"), "ClinPred": (0.028, 0.012, "left"), "AlphaMissense": (0.095, 0.012, "left"),
            "GERP++": (0.165, 0.024, "left"), "PrimateAI": (0.235, 0.083, "center"), "MetaLR": (0.124, 0.086, "left"), "MetaSVM": (0.113, 0.101, "left"), "FATHMM": (0.124, 0.125, "left")}
    for t, (tx, ty, ha_) in lpos.items():  # labels with leader lines; the three Supervised-curated points above 0.08 included
        if t in dl.index: e.annotate(NM2(t), (dl.loc[t, "BP4"], dl.loc[t, "PP3"]), xytext=(tx, ty), textcoords="data", fontsize=5.3, ha=ha_, va="center",
                                     arrowprops=dict(arrowstyle="-", lw=0.4, color=INK2, shrinkA=0.5, shrinkB=2.2))
    e.set_xlim(0.02, 0.27); e.set_ylim(0.003, 0.135)
    e.set_xlabel("AUROC gap, BP4-cited B/LB"); e.set_ylabel("AUROC gap, PP3-cited P/LP"); e.set_title("Gap by predictor", loc="left"); lab(e, "e")
    grid_x(e); grid_y(e); e.legend(loc="upper center", bbox_to_anchor=(0.45, -0.11), borderaxespad=0, ncol=2, fontsize=5.2, handletextpad=0.1, columnspacing=0.6)
    save(fig, "Fig4_sources_provenance")

# ---------------------------------------------------------------- Fig. 5 and 6: earlier Fig. 4 and 5 under new names
_b = open(f"{SRC}/12b_figures.py").read()
_f5 = _b.split("def fig4r2():")[1].split("def fig5r2():")[0]
for old, new in [('save(fig,"Fig4_rankings")', 'save(fig,"Fig5_filtering_rankings")'),
                 # wider gutters: b's rank labels at 6.4 pt and >= 1 cm more between b and c; taller figure for the larger labels
                 ('fig=plt.figure(figsize=(7.1,5.0)); gs=fig.add_gridspec(1,3,width_ratios=[1.0,0.95,0.9],wspace=1.15)',
                  # r9: 6.56 in wide (tight width <= 17.4 cm), names in a and c at 6.0 pt, c's legend wrapped under the panel
                  'fig=plt.figure(figsize=(6.56,5.8)); gs=fig.add_gridspec(1,5,width_ratios=[1.10,0.82,0.92,1.62,1.20],wspace=0,left=0.122,right=0.985)'),
                 ('b=fig.add_subplot(gs[1])', 'b=fig.add_subplot(gs[2])'), ('c=fig.add_subplot(gs[2])', 'c=fig.add_subplot(gs[4])'),
                 # a: pooled = filled black, within gene = open grey
                 ('color=ROSE,lw=0.7,alpha=0.8); a.scatter(d.d_pooled,y-0.15,s=9,color=ROSE,zorder=3,label="Pooled",edgecolor=EDGE,lw=0.3)',
                  'color="#333333",lw=0.7,alpha=0.8); a.scatter(d.d_pooled,y-0.15,s=9,color=DARK,zorder=3,label="Pooled",edgecolor=DARK,lw=0.3)'),
                 ('color=SKY,lw=0.7,alpha=0.8); a.scatter(d.d_within,y+0.15,s=9,color=SKY,zorder=3,label="Within gene",edgecolor=EDGE,lw=0.3)',
                  'color="#a6a6a6",lw=0.7,alpha=0.9); a.scatter(d.d_within,y+0.15,s=10,facecolor="white",zorder=3,label="Within gene",edgecolor="#808080",lw=0.8)'),
                 ('a.set_yticklabels([NM(t) for t in d.index],fontsize=5.0)', 'a.set_yticklabels([NM2(t) for t in d.index],fontsize=6.0)'),
                 ('NM(t),ha="right",va="center",fontsize=4.8); b.text(1.04,X.loc[t,"clean"],NM(t),ha="left",va="center",fontsize=4.8)',
                  'NM2(t),ha="right",va="center",fontsize=6.5); b.text(1.04,X.loc[t,"clean"],NM2(t),ha="left",va="center",fontsize=6.5)'),
                 ('p.index=[NM(t) for t in p.index]', 'p.index=[NM2(t) for t in p.index]'),
                 ('a.legend(loc="lower right",fontsize=5.6,handletextpad=0.2,frameon=False)', 'a.legend(loc="lower right",fontsize=6.0,handletextpad=0.2,frameon=False)'),
                 ('c.set_xlabel("AUROC"); lab(c,"c")', 'c.set_xlabel("AUROC"); lab(c,"c"); c.tick_params(axis="y",labelsize=6.0); hh_,ll_=c.get_legend_handles_labels(); c.legend(hh_,ll_,loc="upper center",bbox_to_anchor=(0.45,-0.09),fontsize=6.0,handletextpad=0.2)'),
                 ('"Labelled before 2024","VUS resolved after 2024-01"', '"Labelled before Jan 2024","VUS/conflicting resolved\\nafter Jan 2024"')]:
    assert _f5.count(old) == 1, old; _f5 = _f5.replace(old, new)
exec("def fig5n():" + _f5)
_f6 = _b.split("def fig5r2():")[1].split("def figS3():")[0]
for old, new in [('nom={"PP3_Supporting":2.08,"PP3_Moderate":4.33}', 'nom={"PP3_Supporting":2.406,"PP3_Moderate":5.790}'), ('save(fig,"Fig5_calibration")', 'save(fig,"Fig6_likelihood_ratios")'),
                 ('color=PERC[j]', 'color=[PER_GREY,SKY,ROSE][j]'), ('a.legend(title="Label vintage"', 'a.legend(title="Label period"'), ('rotation=50,ha="right",fontsize=5.2)', 'rotation=50,ha="right",fontsize=6.0)'),
                 ('a.set_title((ev.replace(chr(95),chr(32))', 'a.set_title((ev.replace(chr(95),chr(32)).replace("Supporting","supporting").replace("Moderate","moderate")')]:
    assert old in _f6, old; _f6 = _f6.replace(old, new)
exec("def fig6n():" + _f6)

# ---------------------------------------------------------------- Fig. 7: recalibration
def fig7n():
    C = pd.read_csv(f"{R2}/O2_threshold_shift_categories.tsv", sep="\t"); O = pd.read_csv(f"{R2}/O2_threshold_oos.tsv", sep="\t")
    TH = pd.read_csv(f"{P}/meta/clingen_thresholds.tsv", sep="\t")
    fig = plt.figure(figsize=(7.1, 5.6)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.55, wspace=0.42)
    a = fig.add_subplot(gs[0, 0])
    D = load(); T = pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz", sep="\t", usecols=["hgvs38", "rare", "period"]); D = D.merge(T, on="hgvs38")
    D["rare"] = D.rare.fillna(False).astype(bool); plp_genes = set(D.loc[D.y == 1, "Gene"].dropna())
    R = D[D.rare & D.Gene.isin(plp_genes) & D.REVEL_raw.notna()]
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
    for G, lb, col in [(R[R.period != "2023-2026"], "Labels up to 2022", SKY), (R[R.period == "2023-2026"], "Labels 2023–2026", ROSE)]:
        gx, lr = curve(G); a.plot(gx, lr, color=col, lw=1.2, label=lb)
    for v, t in [(2.406, "supporting"), (5.790, "moderate"), (33.53, "strong")]:
        a.axhline(v, color=INK2, lw=0.5, ls=":"); a.axhline(1 / v, color=INK2, lw=0.5, ls=":"); a.text(0.01, v * 1.12, t, fontsize=5.4, color=INK2)
    for xv in [0.644, 0.773, 0.932]: a.axvline(xv, color="#a6a6a6", lw=0.5, ls="--")
    a.set_yscale("log"); a.set_ylim(0.01, 400); a.set_xlim(0, 1); a.set_xlabel("REVEL score"); a.set_ylabel("Local likelihood ratio")
    hh, ll = a.get_legend_handles_labels(); a.legend(hh + [Line2D([], [], color="#a6a6a6", lw=0.6, ls="--")], ll + ["Published PP3 thresholds"], loc="lower right", fontsize=5.6)
    a.set_title("REVEL, rare variants in genes with P/LP labels", loc="left"); lab(a, "a")
    cc = C[(C.setting == "genes_with_PLP") & (C.kind == "pooled")].pivot(index="tool", columns="thresholds_from", values=["PLP_pp3_modplus", "BLB_pp3_any"])
    def dumb(ax, key, xlabel, title, letter, xmax):
        dd = cc[key].sort_values("2023_2026_all"); y = np.arange(len(dd))
        ax.hlines(y, dd.le2022_all * 100, dd["2023_2026_all"] * 100, color="#cfcfcf", lw=0.8, zorder=1)
        ax.scatter(dd.le2022_all * 100, y, s=16, color=SKY, zorder=2, label="Thresholds from labels ≤2022", edgecolor=EDGE, lw=0.35)
        ax.scatter(dd["2023_2026_all"] * 100, y, s=16, color=ROSE, zorder=3, label="Thresholds from labels 2023–26", edgecolor=EDGE, lw=0.35)
        ax.set_yticks(y); ax.set_yticklabels([NM(t) for t in dd.index], fontsize=6); ax.set_xlim(0, xmax); ax.set_xlabel(xlabel); grid_x(ax)
        ax.legend(loc="lower right", fontsize=5.2, handletextpad=0.2, borderaxespad=0.2); ax.set_title(title, loc="left"); lab(ax, letter)
    dumb(fig.add_subplot(gs[0, 1]), "PLP_pp3_modplus", "% of rare 2023–26 P/LP variants at PP3 moderate or stronger", "Pathogenic variants gaining evidence", "b", 100)
    dumb(fig.add_subplot(gs[1, 0]), "BLB_pp3_any", "% of rare 2023–26 B/LB variants reaching any PP3 level", "Benign variants given pathogenic evidence", "c", max(5, cc["BLB_pp3_any"].max().max() * 100 * 1.15))
    d = fig.add_subplot(gs[1, 1])
    cats = [("PLP_pp3_modplus", "P/LP: PP3\n≥ moderate"), ("PLP_pp3_strong", "P/LP: PP3\nstrong"), ("BLB_pp3_any", "B/LB: any\nPP3"), ("BLB_bp4_modplus", "B/LB: BP4\n≥ moderate")]
    ins = C[(C.setting == "genes_with_PLP") & (C.kind == "pooled")].groupby("thresholds_from")[[k for k, _ in cats]].median()
    oos = O[O.setting == "genes_with_PLP"].groupby(["tool", "thresholds_from"])[[k for k, _ in cats]].mean().groupby("thresholds_from").median()
    x = np.arange(len(cats)); w = 0.19
    for j, (src, col) in enumerate([("le2022_all", FILL["sky"]), ("2023_2026_all", FILL["rose"])]):
        d.bar(x + (j * 2 - 1.5) * w, ins.loc[src].values * 100, width=w, color=col, edgecolor=EDGE, lw=0.4, label=("Labels ≤2022" if j == 0 else "Labels 2023–26") + ", all genes")
        d.bar(x + (j * 2 - 0.5) * w, oos.loc[src].values * 100, width=w, color=col, edgecolor=EDGE, lw=0.4, hatch="////", label=("Labels ≤2022" if j == 0 else "Labels 2023–26") + ", other gene half")
    d.set_xticks(x); d.set_xticklabels([l for _, l in cats], fontsize=5.8); d.set_ylabel("Median over 11 predictors (%)"); grid_y(d)
    d.set_ylim(0, 100); d.legend(loc="upper right", fontsize=5.0, handlelength=1.4, title="Thresholds calibrated on", title_fontsize=5.2)
    d.set_title("Medians, all genes and held-out genes", loc="left"); lab(d, "d")
    save(fig, "Fig7_recalibration")
    pd.concat({"in_sample": ins, "held_out_genes": oos}).to_csv(f"{R2}/Q_fig7d_medians.tsv", sep="\t")

if __name__ == "__main__":
    for w in sys.argv[1:]: globals()[w](); print("done", w)
