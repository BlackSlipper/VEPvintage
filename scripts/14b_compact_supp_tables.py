#!/usr/bin/env python3
"""Build the compact Supplementary Tables S1-S9 of the paper from the full result workbook.

Usage: 14b_compact_supp_tables.py FULL.xlsx OUT.xlsx [MAPPED_FULL.xlsx]
FULL.xlsx is the workbook written by 14_supp_tables.py (sheets S1a ... S20, earlier numbering).
OUT.xlsx gets the nine tables of the paper; MAPPED_FULL.xlsx (optional) is a copy of the full workbook
with a first sheet that maps every full sheet to its compact table or marks it as kept only in the full set.
Values are copied unchanged; S7e is a per-predictor summary computed here from the replicate-level sheets.
"""
import sys
import numpy as np, pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

FULL, OUT = sys.argv[1], sys.argv[2]
MAPPED = sys.argv[3] if len(sys.argv) > 3 else None
x = pd.ExcelFile(FULL)
rd = lambda s: x.parse(s)

TABLES = [
 ("S1", "Evaluation set and predictors"),
 ("S2", "Accuracy by label vintage (Fig. 2)"),
 ("S3", "Sensitivity analyses of the vintage effect"),
 ("S4", "Gene composition and the pooled gain (Fig. 3)"),
 ("S5", "Label sources, recorded evidence and evidence codes (Fig. 4)"),
 ("S6", "Provenance tiers, rankings and time split (Fig. 5)"),
 ("S7", "Likelihood ratios and recalibration (Figs 6 and 7)"),
 ("S8", "Discontinuities at the published REVEL thresholds (Fig. 8)"),
 ("S9", "Audit of submission-record text"),
]

def merged_release():
    a, b = rd("S4k_release_sets"), rd("S4n_release_ranks")
    m = a.merge(b, on="release", how="outer", validate="one_to_one")
    assert len(m) == len(a) == len(b)
    return m

def per_predictor_designs():
    """Medians over replicates per predictor, threshold set and evaluation labels, for the four calibration designs."""
    keep = ["PLP_pp3_any", "PLP_pp3_modplus", "PLP_pp3_strong", "BLB_pp3_any", "BLB_pp3_modplus", "BLB_bp4_modplus"]
    thr = ["PP3_Supporting", "PP3_Moderate", "PP3_Strong", "BP4_Supporting", "BP4_Moderate"]
    lrs = ["PP3_Supporting_LR", "PP3_Moderate_LR", "PP3_Strong_LR", "BP4_Supporting_LR"]
    TF = {"le2022_all": "up to 2022", "2023_2026_all": "2023-2026", "le2022_matched": "up to 2022",
          "2023_2026_matched": "2023-2026", "le2022": "up to 2022", "2023_2026": "2023-2026"}
    rows = []
    # all genes, single calibration (in-sample), pooled point estimates, genes with a P/LP label
    i = rd("S16d_shares_in_sample")
    i = i[(i.setting == "genes_with_PLP") & (i.kind == "pooled")]
    for _, r in i.iterrows():
        d = dict(tool=r.tool, design="All genes", thresholds_from=TF[r.thresholds_from], evaluated_on="2023-2026", n_replicates=1)
        for k in keep: d[k + " (%)"] = 100 * r[k]
        rows.append(d)
    # class-count matching on all genes, 20 replicates
    s = rd("S16d_shares_size_matched")
    for (tool, tf), g in s.groupby(["tool", "thresholds_from"], sort=False):
        d = dict(tool=tool, design="Matched counts", thresholds_from=TF[tf], evaluated_on="2023-2026", n_replicates=g.rep.nunique())
        for k in keep: d[k + " (%)"] = 100 * g[k].median()
        for t in thr:
            d[t + " threshold, median"] = g[t].median(); d[t + " threshold, 5th pct"] = g[t].quantile(0.05); d[t + " threshold, 95th pct"] = g[t].quantile(0.95)
        rows.append(d)
    # held-out genes, with and without class-count matching, 30 splits
    h = rd("S16d_shares_held_out_30_splits")
    DES = {"split": "Held-out genes", "split_matched": "Held-out genes, matched counts"}
    EV = {"eval_2023_2026": "2023-2026", "eval_le2022": "up to 2022"}
    for (tool, des, ev, tf), g in h.groupby(["tool", "design", "eval_set", "thresholds_from"], sort=False):
        d = dict(tool=tool, design=DES[des], thresholds_from=TF[tf], evaluated_on=EV[ev], n_replicates=g.rep.nunique())
        for k in keep: d[k + " (%)"] = 100 * g[k].median()
        for t in thr:
            d[t + " threshold, median"] = g[t].median(); d[t + " threshold, 5th pct"] = g[t].quantile(0.05); d[t + " threshold, 95th pct"] = g[t].quantile(0.95)
        for l in lrs: d[l.replace("_LR", " achieved LR, median")] = g[l].median()
        rows.append(d)
    out = pd.DataFrame(rows)
    # replicates in which the 2023-2026 thresholds gave the higher P/LP share at PP3 moderate or above (rare 2023-2026 labels)
    def higher(df, keycols):
        a = df[df.thresholds_from.str.startswith("le2022")].set_index(keycols)["PLP_pp3_modplus"]
        b = df[df.thresholds_from.str.startswith("2023")].set_index(keycols)["PLP_pp3_modplus"]
        return (b > a.reindex(b.index)).groupby(level=0).sum()
    out["replicates with higher P/LP share at moderate or above under 2023-2026 thresholds"] = np.nan
    hm = higher(s, ["tool", "rep"])
    for tool, n in hm.items():
        out.loc[(out.tool == tool) & (out.design == "Matched counts") & (out.thresholds_from == "2023-2026"), out.columns[-1]] = n
    for des, lab in DES.items():
        hh = h[(h.design == des) & (h.eval_set == "eval_2023_2026")]
        for tool, n in higher(hh, ["tool", "rep"]).items():
            out.loc[(out.tool == tool) & (out.design == lab) & (out.evaluated_on == "2023-2026") & (out.thresholds_from == "2023-2026"), out.columns[-1]] = n
    return out

S = [  # (new sheet, source, description)
 ("S1a_inclusion_flow", "S1a_inclusion_flow", "Variants retained at each step of building the evaluation set and its rare subsets"),
 ("S1b_composition", "S1b_composition", "Composition of the evaluation set by label vintage"),
 ("S1c_predictors", "S2_predictors", "The 44 predictors: category, training labels, ClinVar in training, ClinGen calibration, coverage of the evaluation set, core status (coverage at least 90%) and reference"),
 ("S2a_period_contrast", "S4c2_period_contrast_2020s", "2023-2026 minus 2020-2022 AUROC: pooled, within-gene and between-gene components with paired gene-cluster 95% intervals, for each predictor on its own rare variants and on the rare variants scored by all core predictors"),
 ("S2b_contrast_le2019", "S4c_vintage_jump", "2023-2026 minus up-to-2019 AUROC, pooled and within-gene, with paired 95% and 90% intervals; equivalence margin +/-0.02"),
 ("S2c_fig2c_common_set", "S4c3_fig2c_common_set", "Values of Fig. 2c: pooled and within-gene AUROC of REVEL, ESM1b and GERP++ by label period on the 96,636 rare variants scored by all 38 core predictors, with gene-cluster 95% intervals (500 replicates)"),
 ("S2d_release_benchmarks", None, "Benchmarks built from each January release 2016-2026 and from November 2026: labels, scored and rare variants, share of rare B/LB labels in genes without a rare P/LP label, share of labels changed by November 2026, and agreement of the rankings with November 2026 (Kendall's tau, maximum rank shift)"),
 ("S2e_release_auc", "S4l_release_auc", "Pooled and within-gene AUROC of the 38 core predictors on the benchmark built from each release (rare variants; Fig. 2d)"),
 ("S3a_gnomad4_auc", "S6a_gnomad4_auc", "Rarity defined by gnomAD v4.1: pooled and within-gene AUROC by label period, and all versus provenance-filtered labels"),
 ("S3b_gnomad4_LR", "S6b_gnomad4_LR", "Rarity defined by gnomAD v4.1: pooled likelihood ratios of the PP3-supporting interval by label period"),
 ("S3c_sensitivity", "S7_sensitivity", "AUROC by label period under alternative definitions: vintage definition, review stars, exclusion of Ambry's labels and the provenance-filtered tier"),
 ("S4a_gene_composition", "S4g_gene_composition", "Share of rare B/LB labels in genes without (or with fewer than five) rare P/LP labels, by label period and source"),
 ("S4b_gene_identity_null", "S4o_gene_identity_null", "Gene-identity null (leave-one-out P/LP fraction of the variant's gene) and each core predictor: pooled AUROC by label period on all and provenance-filtered rare labels, with gene-cluster intervals for the null"),
 ("S4c_composition_matched", "S4p_composition_matched", "2023-2026 minus 2020-2022 AUROC on rare variants scored by all core predictors: full sets, shared gene-by-class (or gene) cells, and after reweighting one period to the other's cell frequencies, with paired gene-cluster intervals"),
 ("S4d_excluded_cells", "S18a_excluded_cells", "Rare 2023-2026 labels scored by all core predictors inside or outside the gene-by-class cells present in 2020-2022: counts by class, source and gene status"),
 ("S4e_retained_excluded_auc", "S18c_retained_excluded_auc", "AUROC of REVEL and of the median core predictor (pooled, 500 gene-cluster replicates) and within-gene AUROC of four predictors on 2020-2022 labels and on retained, excluded and all 2023-2026 labels (Fig. 3d)"),
 ("S5a_evidence_groups", "S19_evidence_groups", "Recorded-evidence groups by label period: AUROC of all P/LP labels against the B/LB labels of each group, and within each group; pooled with 500 gene-cluster replicates, and within genes (Fig. 4a, b)"),
 ("S5b_code_contrasts", "S8f_code_contrast_boot", "Evidence-code contrasts on rare variants scored by all core predictors: AUROC with and without the class-matching code, gap and paired gene-cluster 95% interval (Fig. 4c-e)"),
 ("S5c_within_gene_codes", "S8e_within_gene_codes", "Within-gene evidence-code contrasts"),
 ("S5d_selection_null", "S9_selection_null", "Simulated versus observed PP3/BP4 gaps when codes are assigned from REVEL thresholds, on variants scored by all 38 core predictors"),
 ("S6a_tiers_pooled_within", "S4a_pooled_within_tiers", "Pooled and within-gene AUROC on all labels and in the provenance-filtered, submitter-level and filtered-without-B-only tiers, by label period (1,000 gene-cluster replicates)"),
 ("S6b_filter_drop", "S4b_filter_decomposition", "All minus provenance-filtered AUROC by label period: pooled, within-gene and gene-matched within-gene differences with paired gene-cluster intervals (Fig. 5a)"),
 ("S6c_ranks", "S12a_ranks_core_common", "Ranks of the 38 core predictors on shared rare 2023-2026 variants, pooled and within genes, on all and filtered labels, with gene-cluster bootstrap rank intervals (Fig. 5b)"),
 ("S6d_submitter_level_tier", "S13_submitter_level_tier", "Submitter-level provenance tier: AUROC by label period"),
 ("S6e_time_split", "S8c_time_split", "Time-split benchmark at each reference date: AUROC on variants already labelled, on VUS or conflicting variants resolved later and on variants first submitted later (Fig. 5c)"),
 ("S7a_LR_pooled_MH", "S16a_LR_pooled_MH", "Pooled and Mantel-Haenszel within-gene likelihood ratios for all calibrated intervals by tier and label period, with 1,000 gene-cluster replicates (Fig. 6)"),
 ("S7b_LR_period_ratios", "S16b_LR_period_ratios", "Ratios of likelihood ratios between label periods, pooled and within-gene, with paired intervals"),
 ("S7c_thresholds", "S16c_threshold_shift", "PP3/BP4 thresholds recalibrated with the ClinGen likelihood ratios (2.406, 5.790, 33.53; reciprocals for BP4) on rare labels up to 2022 and from 2023-2026, all and provenance-filtered, by setting and estimate kind"),
 ("S7d_shares_by_design", "S16d_shares_design_summary", "Median over the 11 calibrated predictors, with the interquartile range across predictors, of the shares of rare 2023-2026 P/LP and B/LB variants reaching PP3 (and BP4) under each calibration design (Fig. 7d, e)"),
 ("S7e_per_predictor_designs", None, "Per predictor and calibration design (all genes; class-count matching, 20 replicates; held-out genes and held-out genes with matching, 30 random splits): median shares, thresholds with 5th and 95th percentiles, average likelihood ratios achieved in each interval on the held-out genes (Fig. 7f), and the number of replicates in which the 2023-2026 thresholds gave the higher P/LP share at moderate or above"),
 ("S8a_discontinuities", "S16e_threshold_discontinuity", "Regression-discontinuity steps at the ClinGen PP3/BP4 thresholds (REVEL; AlphaMissense) in the share of submission records citing PP3 or BP4, by year last evaluated, and in the share of rare labels that are P/LP, with intervals and placebo ranks"),
 ("S8b_discontinuity_difference", "S16f_discontinuity_difference", "Difference in discontinuities between periods at each threshold, with paired gene-cluster intervals"),
 ("S8c_discontinuity_submitters", "S16g_discontinuity_submitters", "Steps at REVEL 0.644 (PP3) and 0.290 (BP4) for the submitters citing the code most often near the threshold, by period; letter = label in Fig. 8f"),
 ("S8d_discontinuity_checks", "S16j_discontinuity_extra", "Further checks: PP3 step without the three submitters with threshold-like steps, difference between periods without the two submitters behind the opposite early step, and shares of negated code mentions"),
 ("S8e_discontinuity_loso", "S16h_discontinuity_loso", "Steps at REVEL 0.644 and 0.290 in 2023-2026 records after leaving out each of the ten submitters citing the code most often, and for each of them alone"),
 ("S9_text_audit", "S20_text_audit", "Stratified audit of 160 submission records: stratum (year last evaluated, class, flag), submitter, record text, draft annotation prepared with AI assistance and the authors' annotation"),
]

frames = {}
for new, src, desc in S:
    if new == "S2d_release_benchmarks": df = merged_release(); src = "S4k_release_sets + S4n_release_ranks"
    elif new == "S7e_per_predictor_designs": df = per_predictor_designs(); src = "summary of S16d_shares_in_sample, S16d_shares_size_matched, S16d_shares_held_out_30_splits"
    else: df = rd(src)
    frames[new] = (df, src, desc)

contents = pd.DataFrame([(n.split("_")[0], n, d, s) for n, (df, s, d) in frames.items()], columns=["table", "sheet", "description", "source in the full result workbook"])
with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    head = pd.DataFrame({"Supplementary Tables S1-S9": [f"Table {t}. {title}" for t, title in TABLES] + ["",
        "Further result tables (per-gene lists, replicate-level calibration results, the full discontinuity scan and further sensitivity analyses) are in Supplementary Data 1 (results/Supplementary_Tables_full.xlsx and results/*.tsv)."]})
    head.to_excel(w, sheet_name="Contents", index=False, startrow=0)
    contents.to_excel(w, sheet_name="Contents", index=False, startrow=len(head) + 2)
    for n, (df, s, d) in frames.items(): df.to_excel(w, sheet_name=n, index=False)

# light formatting
wb = load_workbook(OUT)
hdr_fill = PatternFill("solid", fgColor="EAF1F8")
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for c in row: c.font = Font(name="Arial", size=10, bold=c.font.bold)
    if ws.title != "Contents":
        for c in ws[1]: c.font = Font(name="Arial", size=10, bold=True); c.fill = hdr_fill
        ws.freeze_panes = "A2"
    for j, col in enumerate(ws.iter_cols(min_row=1, max_row=min(ws.max_row, 300)), 1):
        L = max((len(str(c.value)) for c in col if c.value is not None), default=8)
        ws.column_dimensions[get_column_letter(j)].width = min(max(10, L + 2), 60)
    if ws.title == "S9_text_audit":
        for row in ws.iter_rows(min_row=2):
            for c in row: c.alignment = Alignment(wrap_text=True, vertical="top")
ws = wb["Contents"]; ws.column_dimensions["A"].width = 14; ws.column_dimensions["B"].width = 30; ws.column_dimensions["C"].width = 110; ws.column_dimensions["D"].width = 45
for c in ws["A"]: c.alignment = Alignment(wrap_text=False)
for r in range(1, ws.max_row + 1):
    ws.cell(r, 3).alignment = Alignment(wrap_text=True, vertical="top")
ws["A1"].font = Font(name="Arial", size=11, bold=True)
wb.save(OUT)

if MAPPED:
    # copy of the full workbook with a mapping sheet in front
    used = {}
    for n, (df, s, d) in frames.items():
        for part in s.replace("summary of ", "").replace(" + ", ", ").split(", "): used[part] = n
    wbf = load_workbook(FULL)
    m = wbf.create_sheet("Mapping", 0)
    m.append(["sheet in this workbook (earlier numbering)", "in the paper's Supplementary Tables"])
    for ws in wbf.worksheets[1:]:
        if ws.title == "Contents": continue
        m.append([ws.title, used.get(ws.title, "full result set only (Supplementary Data 1)")])
    for c in m[1]: c.font = Font(name="Arial", size=10, bold=True)
    m.column_dimensions["A"].width = 38; m.column_dimensions["B"].width = 48
    wbf.save(MAPPED)
print("sheets:", len(frames), "rows:", sum(len(df) for df, s, d in frames.values()))
