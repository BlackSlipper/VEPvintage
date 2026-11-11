#!/usr/bin/env python
"""Assemble supplementary tables S1-S20 into one Excel workbook (with a contents sheet)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import os
R2=f"{P}/results/r2"
D=load()
cov=D[[t for t in TOOLS if t in D.columns]].notna().mean().rename("coverage_eval_set")
S1=META.set_index("tool").join(cov); S1["core_(>=90%)"]=S1.coverage_eval_set>=0.90
REF={"SIFT":"Ng PC, Henikoff S. Nucleic Acids Res 2003;31:3812-3814","SIFT4G":"Vaser R, et al. Nat Protoc 2016;11:1-9",
"PolyPhen2_HDIV":"Adzhubei IA, et al. Nat Methods 2010;7:248-249","PolyPhen2_HVAR":"Adzhubei IA, et al. Nat Methods 2010;7:248-249",
"LRT":"Chun S, Fay JC. Genome Res 2009;19:1553-1561","MutationTaster":"Schwarz JM, et al. Nat Methods 2014;11:361-362",
"MutationAssessor":"Reva B, et al. Nucleic Acids Res 2011;39:e118","FATHMM":"Shihab HA, et al. Hum Mutat 2013;34:57-65",
"PROVEAN":"Choi Y, et al. PLoS One 2012;7:e46688","VEST4":"Carter H, et al. BMC Genomics 2013;14(Suppl 3):S3",
"MetaSVM":"Dong C, et al. Hum Mol Genet 2015;24:2125-2137","MetaLR":"Dong C, et al. Hum Mol Genet 2015;24:2125-2137",
"MetaRNN":"Li C, et al. Genome Med 2022;14:115","M-CAP":"Jagadeesh KA, et al. Nat Genet 2016;48:1581-1586",
"REVEL":"Ioannidis NM, et al. Am J Hum Genet 2016;99:877-885","MutPred2":"Pejaver V, et al. Nat Commun 2020;11:5918",
"MVP":"Qi H, et al. Nat Commun 2021;12:510","MPC":"Samocha KE, et al. bioRxiv 2017;148353",
"PrimateAI":"Sundaram L, et al. Nat Genet 2018;50:1161-1170","DEOGEN2":"Raimondi D, et al. Nucleic Acids Res 2017;45:W201-W206",
"BayesDel_addAF":"Feng BJ. Hum Mutat 2017;38:243-251","BayesDel_noAF":"Feng BJ. Hum Mutat 2017;38:243-251",
"ClinPred":"Alirezaie N, et al. Am J Hum Genet 2018;103:474-483","LIST-S2":"Malhis N, et al. Nucleic Acids Res 2020;48:W154-W161",
"VARITY_R":"Wu Y, et al. Am J Hum Genet 2021;108:1891-1906","VARITY_ER":"Wu Y, et al. Am J Hum Genet 2021;108:1891-1906",
"VARITY_R_LOO":"Wu Y, et al. Am J Hum Genet 2021;108:1891-1906","VARITY_ER_LOO":"Wu Y, et al. Am J Hum Genet 2021;108:1891-1906",
"ESM1b":"Brandes N, et al. Nat Genet 2023;55:1512-1522","EVE":"Frazer J, et al. Nature 2021;599:91-95",
"AlphaMissense":"Cheng J, et al. Science 2023;381:eadg7492","CADD":"Rentzsch P, et al. Nucleic Acids Res 2019;47:D886-D894",
"DANN":"Quang D, et al. Bioinformatics 2015;31:761-763","FATHMM-MKL":"Shihab HA, et al. Bioinformatics 2015;31:1536-1543",
"FATHMM-XF":"Rogers MF, et al. Bioinformatics 2018;34:511-513","Eigen":"Ionita-Laza I, et al. Nat Genet 2016;48:214-220",
"Eigen-PC":"Ionita-Laza I, et al. Nat Genet 2016;48:214-220","GERP++":"Davydov EV, et al. PLoS Comput Biol 2010;6:e1001025",
"phyloP100way":"Pollard KS, et al. Genome Res 2010;20:110-121","phastCons100way":"Siepel A, et al. Genome Res 2005;15:1034-1050",
"SiPhy":"Garber M, et al. Bioinformatics 2009;25:i54-i62"}
S1["reference"]=[REF.get(t,"see dbNSFP 4.8a documentation (Liu et al. Genome Med 2020;12:103)") for t in S1.index]
V=pd.read_csv(f"{P}/results/variant_level_annotated.tsv.gz",sep="\t",usecols=["final","rare","period","conc_submitters","found"])
V=V[V.found==1]
sub=V.conc_submitters.fillna("(none concordant)").str.split("|").explode()
S2=pd.DataFrame({"submitter":sub}).join(V[["final","period","rare"]]).groupby(["submitter","final","period"]).size().unstack(["final","period"]).fillna(0).astype(int)
S2["total"]=S2.sum(axis=1); S2=S2.sort_values("total",ascending=False).head(60)
ps={"Ambry Genetics":"Ambry","Labcorp Genetics (formerly Invitae), Labcorp":"Labcorp"}
S2.insert(0,"name_in_text",[ps.get(i,"") for i in S2.index])
S2.columns=["_".join(map(str,c)).strip("_") if isinstance(c,tuple) else c for c in S2.columns]
rd=lambda f: pd.read_csv(f,sep="\t")
S5a=rd(f"{P}/results/calibration_LR_by_period_rare.tsv"); S5a.insert(0,"subset","rare")
S5b=rd(f"{P}/results/calibration_LR_by_period_rare_clean.tsv"); S5b.insert(0,"subset","rare_provenance_filtered")
S5c=rd(f"{P}/results/calibration_LR_by_period.tsv"); S5c.insert(0,"subset","all")
S12=rd(f"{R2}/B_ranks_core_common.tsv"); S12b=rd(f"{R2}/B2_extra_tools_ranks.tsv") if os.path.exists(f"{R2}/B2_extra_tools_ranks.tsv") else pd.DataFrame()
sheets=[
("S1_predictors",S1.reset_index(),"Predictor list: category, training labels, ClinGen calibration, coverage, core status (>=90% coverage), reference"),
("S2_submitters",S2.reset_index(),"Top 60 submitters by concordant labels, by class and period; short names used in the text"),
("S3a_inclusion_flow",rd(f"{R2}/K_inclusion_flow.tsv"),"Variants retained at each step of building the evaluation set and its rare subsets"),
("S3b_composition",rd(f"{P}/results/composition_by_vintage.tsv"),"Composition of the evaluation set by label vintage"),
("S4_sensitivity",rd(f"{P}/results/robustness.tsv").pivot_table(index=["analysis","tool"],columns="period",values="auc").reset_index(),"Sensitivity analyses: vintage definition, review stars, Ambry exclusion, filtered tiers"),
("S5_calibration_LR",pd.concat([S5a,S5b,S5c]),"Pooled likelihood ratios (variant-level bootstrap) for all calibrated intervals, periods and subsets"),
("S6_DMS",rd(f"{P}/results/dms_before_after_pub.tsv"),"DMS assays matched to ClinVar labels before and after publication"),
("S7a_full_vs_filtered",rd(f"{P}/results/benchmark_full_vs_clean.tsv"),"AUROC on all versus provenance-filtered labels (all vintages)"),
("S7b_rare_2023_26",rd(f"{P}/results/benchmark_clean_rare_recent.tsv"),"Rare variants labelled 2023-2026: all, filtered and conservation-filtered labels (200-replicate CIs; see S15a for 1,000-replicate CIs)"),
("S7c_time_split",rd(f"{P}/results/prospective_auc.tsv"),"Time-split benchmark by reference date"),
("S7d_provenance_contrasts",rd(f"{P}/results/evidence_source_contrasts.tsv"),"Evidence-code and label-source contrasts (variants scored by all core predictors)"),
("S7e_within_gene_codes",rd(f"{P}/results/within_gene_code_contrasts.tsv"),"Within-gene evidence-code contrasts"),
("S7f_code_contrast_boot",rd(f"{R2}/L_code_contrast_bootstrap.tsv"),"Evidence-code contrasts on rare variants scored by all core predictors: AUROC with and without the class-matching code, gap and paired gene-cluster 95% interval"),
("S8_expert_panel",rd(f"{P}/results/expert_panel_auc.tsv"),"Expert-panel (three-star) labels by period and by PP3/BP4 citation"),
("S9_submitter_level_tier",rd(f"{P}/results/submitter_level_tier.tsv"),"Submitter-level provenance tier"),
("S10a_within_gene",rd(f"{P}/results/review_within_gene.tsv"),"Within-gene AUROC by vintage and benign-label source"),
("S10b_filtered_by_source",rd(f"{P}/results/review_clean_by_source.tsv"),"Filtered tier by source, with and without Labcorp-only labels"),
("S11_selection_null",rd(f"{P}/results/selection_null_common.tsv"),"Simulated versus observed PP3/BP4 contrasts on variants scored by all 38 core predictors"),
("S12a_ranks_core_common",S12,"Ranks of the 38 core predictors on shared variants (rare, 2023-2026), pooled and within-gene, with 1,000-replicate gene-cluster bootstrap intervals"),
("S12b_ranks_extra_tools",S12b,"MutPred2 and EVE ranked with the 38 core predictors on the variants they score"),
("S13_calib_submitter_level",rd(f"{P}/results/calibration_submitter_level.tsv"),"Submitter-level tier: PP3 likelihood ratios"),
("S14a_gnomad4_auc",rd(f"{P}/results/sensitivity_gnomad4.tsv"),"gnomAD v4.1 rarity: pooled (variants scored by all core predictors) and within-gene AUROC by period; full vs filtered (all variants scored by each predictor)"),
("S14b_gnomad4_LR",rd(f"{P}/results/sensitivity_gnomad4_LR.tsv"),"gnomAD v4.1 rarity: pooled PP3-supporting likelihood ratios"),
("S15a_pooled_within_tiers",rd(f"{R2}/A_tiers_pooled_within.tsv"),"Pooled and within-gene AUROC for all, filtered, submitter-level and filtered-without-B-only tiers by period (1,000 gene-cluster replicates)"),
("S15b_filter_decomposition",rd(f"{R2}/A_full_minus_filtered.tsv"),"All minus filtered AUROC: pooled, within-gene and gene-matched within-gene differences with paired intervals"),
("S15c_vintage_jump",rd(f"{R2}/E_jump.tsv"),"2023-2026 minus <=2019 AUROC, pooled and within-gene, paired intervals (95% and 90%); equivalence at +/-0.02"),
("S15c2_period_contrast_2020s",rd(f"{R2}/J_period_contrast_2020_22_vs_2023_26.tsv"),"2023-2026 minus 2020-2022 AUROC: pooled, within-gene and between-gene components with paired gene-cluster intervals, for each predictor on its own rare variants and on the rare variants scored by all core predictors, all labels and the provenance-filtered tier"),
("S15c3_fig2c_common_set",rd(f"{R2}/Q_fig2c_common.tsv"),"Pooled and within-gene AUROC of REVEL, ESM1b and GERP++ by label period on the 96,636 rare variants scored by all 38 core predictors, with gene-cluster bootstrap 95% intervals (500 replicates, seed 2026; Fig. 2c; script 32c)"),
("S15h_equal_gene_weight",rd(f"{R2}/G_equal_gene_weight_diff.tsv"),"Sensitivity: all minus filtered within-gene AUROC with genes weighted equally (mean per-gene AUROC, genes with >=3 P/LP and >=3 B/LB), own genes and genes qualifying in both tiers"),
("S15i_equal_gene_weight_jump",rd(f"{R2}/G_equal_gene_weight_jump.tsv"),"Sensitivity: 2023-2026 minus <=2019 within-gene AUROC with genes weighted equally"),
("S15j_gene_drivers",rd(f"{R2}/I_gene_drivers_blb_2023_26.tsv"),"Rare B/LB labels of 2023-2026 in genes without a rare P/LP label: labels per gene, share supported only by Ambry or Labcorp, share kept by the per-variant filter"),
("S15k_release_sets",rd(f"{R2}/H_download_date_sets.tsv"),"Benchmarks built from each January release 2016-2026 and November 2026: labels, scored, rare P/LP and B/LB, share of rare B/LB in genes without rare P/LP, share of labels changed by November 2026"),
("S15l_release_auc",rd(f"{R2}/H_download_date_auc.tsv"),"Pooled and within-gene AUROC of the 38 core predictors on benchmarks built from each release (rare variants)"),
("S15m_release_change",rd(f"{R2}/H_download_date_change.tsv"),"Change from the January 2020 to the November 2026 benchmark, pooled and within-gene, with gene-cluster bootstrap intervals"),
("S15n_release_ranks",rd(f"{R2}/H_download_date_ranks.tsv"),"Rankings on variants scored by all core predictors in each release: Kendall tau against November 2026"),
("S15d_pre2020_reasons",rd(f"{R2}/D_removal_reasons.tsv"),"Reasons for removal by the filter and timing of flagging records, labels <=2019 and 2023-2026"),
("S15e_pre2020_subsets",rd(f"{R2}/D_subsets_le2019.tsv"),"Labels <=2019: AUROC after removing each provenance criterion, and as-of-vintage filter"),
("S15f_recent_subsets",rd(f"{R2}/D_subsets_2023-2026.tsv"),"Labels 2023-2026: AUROC after removing each provenance criterion, and as-of-vintage filter"),
("S15g_gene_composition",rd(f"{R2}/B3_gene_composition.tsv") if os.path.exists(f"{R2}/B3_gene_composition.tsv") else pd.DataFrame(),"Share of rare B/LB labels in genes without (or with fewer than 5) rare P/LP labels, by period and source"),
("S16a_LR_pooled_MH",rd(f"{R2}/C_LR_pooled_MH.tsv"),"Pooled and Mantel-Haenszel within-gene likelihood ratios for all calibrated intervals by tier and period (1,000 gene-cluster replicates)"),
("S16b_LR_period_ratios",rd(f"{R2}/C_LR_ratios.tsv"),"Ratios of likelihood ratios between periods, pooled and within-gene, with paired intervals"),
("S17_resolved_VUS",rd(f"{R2}/F_resolved_VUS.tsv"),"VUS/conflicting at January 2024 and resolved afterwards, by provenance tier, with gene-cluster CIs"),
("S15o_gene_identity_null",rd(f"{R2}/M_gene_null.tsv"),"Gene-identity null (leave-one-out P/LP fraction of the variant's gene) and each core predictor: pooled AUROC by period and tier on rare variants, with gene-cluster intervals for the null"),
("S15p_composition_matched",rd(f"{R2}/N_composition_matched.tsv"),"2023-2026 minus 2020-2022 AUROC on rare variants scored by all core predictors: full sets, shared gene x class (or gene) cells, and after reweighting one period to the other's composition, with paired gene-cluster intervals"),
("S16c_threshold_shift",rd(f"{R2}/O2_threshold_shift.tsv"),"PP3/BP4 thresholds recalibrated with the ClinGen likelihood ratios (2.406, 5.790, 33.53; BP4 reciprocals) on rare labels up to 2022 and from 2023-2026, all and provenance-filtered; setting genes_with_PLP (main) or all_genes; kind pooled or within_gene (Mantel-Haenszel), point estimate or 5th percentile of 1,000 gene-cluster replicates (_lower95); lo/hi = bootstrap range of the threshold; reached_share = share of replicates in which the strength was reached (script 30b)"),
("S16d_shares_in_sample",rd(f"{R2}/O2_threshold_shift_categories.tsv"),"Shares of rare 2023-2026 P/LP and B/LB variants reaching each PP3/BP4 strength under thresholds calibrated on labels up to 2022 (le2022_all) or 2023-2026 (2023_2026_all), by setting and estimate kind (script 30b)"),
("S16d_shares_held_out_genes",rd(f"{R2}/O2_threshold_oos.tsv"),"Held-out genes: thresholds (pooled point estimates) calibrated on one random half of the genes (fold) and applied to the rare 2023-2026 variants of the other half; shares by class and strength (script 30b, single split, seed 7; the 30-split analysis is in S16d_shares_held_out_30_splits)"),
("S16d_shares_size_matched",rd(f"{R2}/O2_threshold_sizematched.tsv"),"Calibration sets matched on class counts (20 replicates, seed 11): thresholds and shares of rare 2023-2026 variants per replicate (script 30c)"),
("S16d_shares_held_out_30_splits",rd(f"{R2}/O3_threshold_combined.tsv"),"Held-out genes with and without class-count matching, 30 random splits of the genes with at least one P/LP label into halves (seed 31; script 30d). Thresholds (pooled point estimates) were calibrated on the rare labels up to 2022 or from 2023-2026 of one half, as they are (design split) or matched on class counts (design split_matched), and applied to the rare labels of the other half from 2023-2026 (eval_2023_2026) or up to 2022 (eval_le2022). Columns give the thresholds, the shares of P/LP and B/LB variants reaching each strength, and the likelihood ratio achieved in each PP3/BP4 evidence interval (_LR; P/LP:B/LB for PP3, B/LB:P/LP for BP4) with the number of variants in the interval (_n)"),
("S16d_shares_design_summary",rd(f"{R2}/Q_fig7de_designs.tsv"),"Median over the 11 calibrated predictors, with the interquartile range across predictors, of the share of rare 2023-2026 P/LP variants at PP3 moderate or stronger and of B/LB variants reaching any PP3 level, by calibration design (Fig. 7d, e). For replicated designs each predictor was first summarised by its median over replicates (20 for matched counts, 30 for held-out genes), then the median over predictors was taken"),
("S16e_threshold_discontinuity",rd(f"{R2}/P_rd_jumps.tsv"),"Regression-discontinuity steps at ClinGen PP3/BP4 thresholds (REVEL; AlphaMissense) in the share of submission records citing PP3 or BP4, by year last evaluated, and in the share of rare labels that are P/LP, by label period (all labels; previously VUS): local linear fits with a triangular kernel (bandwidths 0.03, 0.05, 0.08), gene-cluster bootstrap intervals, and the rank of each step among placebo cut points more than one bandwidth from any published threshold"),
("S16f_discontinuity_difference",rd(f"{R2}/P_rd_diff.tsv"),"Difference in discontinuities (later minus earlier period) at each threshold, with paired gene-cluster bootstrap intervals"),
("S16g_discontinuity_submitters",rd(f"{R2}/P_rd_submitters.tsv").merge(pd.read_csv(f"{R2}/P_rd_submitter_letters.tsv",sep="\t").rename(columns={"Unnamed: 0":"submitter"}),on="submitter",how="left"),"Steps at REVEL 0.644 (PP3) and 0.290 (BP4) for the submitters citing the code most often within 0.05 of the threshold, by period; letter = label used in Fig. 8f"),
("S16h_discontinuity_loso",rd(f"{R2}/P_rd_loso.tsv"),"Steps at REVEL 0.644 (PP3) and 0.290 (BP4) in submission records evaluated 2023-2026 after leaving out each of the ten submitters citing the code most often, and for each of them alone"),
("S16i_discontinuity_scan",rd(f"{R2}/P_rd_scan.tsv"),"Step estimates at every cut point of the scan (grid 0.005, bandwidth 0.05), by level, score, side and period; placebo = more than one bandwidth from any published or legacy threshold"),
("S16j_discontinuity_extra",rd(f"{R2}/P_rd_extra.tsv"),"Further checks of the discontinuity analysis: PP3 step at 0.644 in 2023-2026 records without the three submitters with threshold-like steps; difference between periods without the two submitters behind the negative pre-2022 step; share of PP3/BP4 mentions that are negated in the text, by period; PP3 steps with negated mentions counted as non-citations"),
("S18a_excluded_cells",rd(f"{R2}/Q_excluded_cells.tsv"),"Rare 2023-2026 labels scored by all core predictors, inside (shared_cell True) or outside the gene x class cells present in 2020-2022: counts by class, source and gene status (script 32)"),
("S18b_excluded_genes",rd(f"{R2}/Q_excluded_genes.tsv"),"Genes of the 2023-2026 labels outside the shared cells, by class (script 32)"),
("S18c_retained_excluded_auc",rd(f"{R2}/Q_excluded_auc.tsv"),"AUROC of REVEL and the median core predictor (pooled, 500 gene-cluster replicates) and within-gene AUROC of four predictors on 2020-2022 labels and on retained, excluded and all 2023-2026 labels (script 32)"),
("S19_evidence_groups",rd(f"{R2}/Q_provenance_groups.tsv"),"Recorded-evidence groups by period: AUROC within each group, and of all P/LP labels against the B/LB labels of each group; pooled with 500 gene-cluster replicates, and within genes (script 32)"),
("S20_text_audit",rd(f"{R2}/R_audit_annotated.tsv").rename(columns={"annotation":"draft_annotation"}).drop(columns=["annotated_by"],errors="ignore").assign(author_annotation="",annotation_status="draft prepared with AI assistance; author_annotation to be entered by the authors"),"Stratified audit of 160 submission records (scripts 33, 33b): stratum (year last evaluated | class | flag), submitter, record text, draft annotation (prepared with AI assistance) and the authors' annotation. The sample has 20 records per stratum and is not representative of all records"),
]

import re as _re
MAP={1:2,2:3,3:1,4:7,10:5,15:4,14:6,7:8,11:9,8:10,6:11,12:12,9:13,17:14,5:15,16:16,13:17,18:18,19:19,20:20}
def _rn(x): return _re.sub(r"\bS(\d+)(?=[a-z_]|\b)",lambda m:"S"+str(MAP[int(m.group(1))]),x)
sheets=[(_rn(n),df,_rn(d)) for n,df,d in sheets]
sheets.sort(key=lambda t:(int(_re.match(r"S(\d+)",t[0]).group(1)),t[0]))
out=f"{P}/results/Supplementary_Tables.xlsx"
with pd.ExcelWriter(out) as w:
    pd.DataFrame({"sheet":[s[0] for s in sheets],"description":[s[2] for s in sheets]}).to_excel(w,"Contents",index=False)
    for nm,df,_ in sheets: df.to_excel(w,nm[:31],index=False)
print("written",out,len(sheets))
