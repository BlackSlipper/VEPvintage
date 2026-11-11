#!/usr/bin/env python
"""Place the released derived tables under results/ with the column names the analysis scripts expect.
Run once after unpacking the release. dbNSFP scores are not redistributed: create
results/dbnsfp_shard*.tsv.gz with scripts/03_fetch_dbnsfp.py first, then run this script, then 09, 10 and 11
(10 falls back to the released concordant code-citing SCV table when the raw ClinVar submission_summary is absent)."""
import os, shutil, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.environ.get("VEP_ROOT",os.path.dirname(HERE))
D=os.path.join(ROOT,"data"); R=os.path.join(ROOT,"results"); os.makedirs(R,exist_ok=True)
T=pd.read_csv(os.path.join(D,"benchmark_tiers.tsv.gz"),sep="\t",dtype={"VariationID":str})
T=T.rename(columns={"explicit_PP3_or_BP4":"explicit_comp","computational_text_nonA":"comp_nonambry",
                    "conservation_text_nonA":"cons_nonambry","submitterA_only":"ambry_only"})
T.to_csv(os.path.join(R,"benchmark_tiers.tsv.gz"),sep="\t",index=False,compression="gzip")
for f in ["variant_history.tsv.gz","scv_textmining.tsv.gz","scv_acmg_codes.tsv.gz","gnomad4_af_merged.tsv"]:
    src=os.path.join(D,f)
    if os.path.exists(src): shutil.copy(src,os.path.join(R,f))
src=os.path.join(D,"gnomad2_exac_af.tsv.gz")  # gnomAD v2.1.1 / ExAC AFs used for the rarity filter (read by 05b and 09)
if os.path.exists(src): shutil.copy(src,os.path.join(R,"af_shard_release.tsv.gz"))
for f in ["download_date_labels.tsv.gz","extra_alleles_hg38.tsv"]:  # benchmarks by ClinVar release (25)
    src=os.path.join(D,f)
    if os.path.exists(src): shutil.copy(src,os.path.join(R,f))
src=os.path.join(D,"gnomad2_exac_af_extra.tsv.gz")
if os.path.exists(src): shutil.copy(src,os.path.join(R,"af_extra.tsv.gz"))
src=os.path.join(D,"extra_alleles_hg38.tsv")  # ids to fetch with 03_fetch_dbnsfp.py -> results/dbnsfp_extra.tsv.gz
if os.path.exists(src): pd.read_csv(src,sep="\t").hgvs38.drop_duplicates().to_csv(os.path.join(R,"ids_extra.txt"),index=False,header=False)
print("results/ prepared:",sorted(os.listdir(R)))
