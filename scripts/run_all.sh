#!/bin/bash
# Re-run the analyses from the released tables, in the order tested (see README).
# Needs results/dbnsfp_shard*.tsv.gz from scripts/03_fetch_dbnsfp.py. VEP_BOOT=100 gives a quick check.
# 08 (ProteinGym DMS) runs only if data/dms/ holds the ProteinGym v1.3 files; 01-04 and 20 need the raw ClinVar files.
set -e
cd "$(dirname "$0")"
PY=${PY:-python3}; L=../logs; mkdir -p $L
run(){ s=$1; shift; echo "$(date +%T) $s $*"; env "$@" $PY $s.py ${ARGS:-} > $L/$s${TAG:-}.log 2>&1; }
for s in 00_use_release_data 09_mechanism 10_evidence_codes 11_clean_benchmark 05_vintage_auc 05b_relative 06_prospective; do run $s; done
run 07_calibration_drift
TAG=_rare run 07_calibration_drift RARE=1
TAG=_rareclean run 07_calibration_drift RARE=1 CLEAN=1
[ -d ../data/dms/DMS_ProteinGym_substitutions ] && run 08_dms
[ -f ../results/dms_clinvar_merged.tsv.gz ] && run 08b_dms_before_after || true
for s in 11b_recent 13_robustness 15_within_gene 16_expert_panel 17_submitter_level 18_review_fixes 19_rank_ci_sublevel \
         20b_selection_null_common 21_gnomad4_sensitivity 22_round2 23_round2b 24_equal_gene_weight 26_gene_drivers \
         27_rigor 28_code_contrast_boot 29_composition 30b_threshold_shift_v2 30c_threshold_sizematched 30d_threshold_combined 31_threshold_rd 31b_submitters \
         32b_vepbench_labels 32_composition_detail 32c_fig2c_common; do run $s; done
if [ -s ../data/clinvar_tab/submission_summary_2026-11.txt.gz ]; then run 31d_extra; run 33_audit_sample; run 33b_audit_annotations; fi  # need the raw submission_summary file
[ -s ../results/dbnsfp_extra.tsv.gz ] && run 25_download_date  # needs 03_fetch_dbnsfp.py on results/ids_extra.txt
# figures in the order of the manuscript (12f: Figs 1 and 3-6; 12g: Figs 2 and 7, reading 12_figures and 12f;
# 12e_fig8: Fig. 8; 12b figS3 writes FigS2_relative_slopes; 12c figS2 writes FigS3_leakage) and the supplementary tables
ARGS="fig1n fig3n fig4n fig5n fig6n" run 12f_revision_figs
ARGS="fig2n fig7n" run 12g_final_figs
run 12e_fig8
ARGS="figS3" run 12b_figures
ARGS="figS1 figS2" run 12c_suppfigs
run 14_supp_tables
ARGS="../results/Supplementary_Tables.xlsx ../results/Supplementary_Tables_S1-S9.xlsx ../results/Supplementary_Tables_full.xlsx" run 14b_compact_supp_tables
# worked example of the evaluation tool; compare with ../examples/expected_output
$PY vepbench.py --scores ../examples/phyloP100way_example.tsv.gz --out ../examples/example_run > $L/vepbench_example.log 2>&1
echo "all done; logs in $L"
