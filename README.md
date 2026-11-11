# Label vintage and gene composition in ClinVar benchmarks of missense variant effect predictors

Code for the manuscript "Label vintage and gene composition shape ClinVar benchmarks of missense variant effect
predictors" (Shin HD, Kim H; Briefings in Bioinformatics, Problem Solving Protocol, submitted 2026), with the label
table and worked example needed to run the evaluation tool `vepbench.py`. The derived tables, result tables, figures
and supplementary tables are in Supplementary Data 1 of the paper (release v1.3.1), which also contains this code.

## Quick start: evaluating your own score
`scripts/vepbench.py` runs steps 3-5 and 7 of the protocol (Table 1 of the paper) for any missense score:

    pip install -r requirements.txt
    python scripts/vepbench.py --scores my_scores.tsv.gz --out my_run/
    python scripts/vepbench.py --scores examples/phyloP100way_example.tsv.gz --out example_run/   # worked example, under 10 s

The score file is tab-separated with columns `hgvs38` (GRCh38 HGVS g. notation, e.g. `chr11:g.126277517A>G`; or
`--id-col VariationID`) and `score` (higher = more damaging; otherwise `--lower-is-damaging`). The tool joins the scores
to `data/vepbench_labels.tsv.gz`, keeps rare variants (max population AF < 1e-4; `--all-frequencies` to keep all) and
writes `auc_by_period.tsv` (pooled and within-gene AUROC and a gene-identity null by label period, all labels and the
provenance-filtered tier), `composition.tsv` (change 2020-2022 -> 2023-2026 on full sets, within genes, on shared gene x
class cells and after reweighting), `evidence_groups.tsv` (B/LB labels of each recorded-evidence group against all P/LP
labels of the period) and `summary.txt`. Intervals come from a gene-cluster bootstrap (`--boot`, default 200; seed 2026)
applied to every quantity, so differences are paired. `examples/expected_output/` holds the output of the worked example.

## Contents
- `scripts/vepbench.py`: the evaluation tool (steps 3-5 and 7 of the protocol).
- `scripts/`: the analysis pipeline, numbered in run order. `common.py` holds shared settings and functions;
  `run_all.sh` runs the analyses, figures and supplementary tables from the released tables in the tested order.
- `data/vepbench_labels.tsv.gz`: per-variant labels read by `vepbench.py` (class, vintage and period, stars, rare flag,
  provenance-filtered flag, source and recorded-evidence group; script 32b). `data/clinvar_snapshots_used.txt`: the
  archived ClinVar `variant_summary` files (read by `00_download.sh`).
- `meta/`: predictor categories and training-label information (`predictors.tsv`), ClinGen PP3/BP4 score intervals
  (`clingen_thresholds.tsv`) and the UniProt-to-gene map used with ProteinGym (`uniprot_entry2gene.json`).
- `examples/`: `phyloP100way_example.tsv.gz` (raw phyloP100way, a UCSC conservation score, for the evaluation variants)
  and `expected_output/` from the worked example.

## Re-running the analyses
Copy the `data/` folder of Supplementary Data 1 into this directory; it adds the label histories, benchmark tiers,
text-mining and evidence-code tables, allele frequencies and the label sets of earlier releases. dbNSFP scores are not
redistributed (see Licences), so they have to be fetched first:

    pip install -r requirements.txt
    mkdir -p results
    python -c "import pandas as pd; V=pd.read_csv('data/variant_history.tsv.gz',sep='\t',usecols=['hgvs38','final']); V[V.final.isin(['PLP','BLB'])].hgvs38.dropna().drop_duplicates().to_csv('results/ids_hg38_plpblb.txt',index=False,header=False)"
    for i in 0 1 2 3; do python scripts/03_fetch_dbnsfp.py results/ids_hg38_plpblb.txt results/dbnsfp_shard$i.tsv.gz $i 4 & done; wait
    python scripts/00_use_release_data.py     # copies data/ tables to results/ under the names the scripts use
    python scripts/03_fetch_dbnsfp.py results/ids_extra.txt results/dbnsfp_extra.tsv.gz   # 2,737 extra ids (script 25)
    python scripts/09_mechanism.py            # -> results/variant_level_annotated.tsv.gz
    python scripts/10_evidence_codes.py       # -> results/variant_level_annotated2.tsv.gz
    python scripts/11_clean_benchmark.py      # rebuilds results/benchmark_tiers.tsv.gz and the filtered benchmark
    python scripts/22_round2.py               # within-gene decomposition, ranks, Mantel-Haenszel LRs (results/r2/)

`scripts/run_all.sh` runs these and all remaining analyses in the tested order, then the figure scripts (12f, 12g,
12e_fig8, 12b, 12c), the supplementary tables (14, 14b) and the worked example of `vepbench.py`; scripts run by hand
should follow the same order (23, for example, reads the output of 11b). Script 07 reads the environment variables
`RARE=1` and `CLEAN=1` for the rare and filtered calibration tables. Script 08 needs the ProteinGym v1.3 substitution
files, and 31d, 33 and 33b the raw `submission_summary` file.

Environment variables: `VEP_ROOT` (project root; defaults to the directory above `scripts/`), `VEP_BOOT` (number of
bootstrap replicates; default 1000, the value used in the paper; 100 gives a quick check with identical point estimates),
`DBNSFP_GLOB` (location of the dbNSFP shards). Random seeds are fixed (gene-cluster bootstrap: 2026).

MyVariant.info serves the dbNSFP version current at the time of the query. The paper used dbNSFP 4.8a (retrieved
November 2026); a later version may change some scores slightly.

With 100 bootstrap replicates each script takes a few minutes or less on a 24-core server (script 18 about 11 min);
with the paper's 1,000 replicates, script 23 takes about 40 min. The dbNSFP download takes several hours. The analyses
were re-run on a second server from the released tables and a copy of the dbNSFP shards used for the paper: every
point estimate matched, and only bootstrap intervals differed with the number of replicates (the README of
Supplementary Data 1 lists the scripts checked).

## Full pipeline from the raw ClinVar archives
00 download ClinVar archives (00_download.sh) -> 01 parse snapshots (01_parse_snapshots.py, one archive file per call) -> 02 build label history -> 03 fetch dbNSFP 4.8a scores via
MyVariant.info (03c gnomAD v2.1.1/ExAC allele frequencies) -> 04 text-mine submission records -> 05/05b vintage AUROC and
relative slopes -> 06 time-split benchmark -> 07 calibration likelihood ratios -> 08/08b ProteinGym DMS (08b: labels before and after assay publication) -> 09/10 provenance and
ACMG-code contrasts -> 11/11b provenance-filtered benchmarks -> 13 sensitivity analyses -> 15-19 within-gene, expert-panel,
submitter-level tier, rank intervals -> 20/20b selection simulation and text-mining audit -> 21 gnomAD v4.1 rarity (03d
fetch) -> 22/23 within-gene decomposition, ranks on shared variants, Mantel-Haenszel LRs, pre-2020 diagnosis -> 24
equal-gene-weight sensitivity -> 25b/25 benchmarks rebuilt from each ClinVar release (25b needs the parsed snapshots from 01)
-> 26 genes behind the composition effect -> 27 inclusion flow, paired intervals for the 2020-2022 vs 2023-2026 contrast
and the text-mining audit -> 28 paired intervals for the evidence-code contrasts -> 29 gene-identity null and composition
matching -> 30b/30c/30d recalibration of PP3/BP4 thresholds on recent labels with the ClinGen likelihood ratios -> 31 regression discontinuities at the ClinGen REVEL and
AlphaMissense thresholds in submission records and labels (31b single submitters; 31d leave-three-out, negation and period-difference checks) -> 32/32b/32c recorded-evidence groups, the label table read by vepbench.py and the
common variant set of Fig. 2c -> 33/33b stratified audit sample of submission-record text -> 12f/12g/12e_fig8/12b/12c
figures -> 14/14b supplementary tables.
Script 00 downloads the 48 archived `variant_summary` files listed in `data/clinvar_snapshots_used.txt`, the November 2026
release and the ProteinGym v1.3 files (needed by 08). NCBI's tab-delimited archive does not keep past `submission_summary`
files, so a fresh download will not match the November 2026 file used here (`submission_summary_2026-11.txt.gz`); the tables
`scv_textmining.tsv.gz` and `scv_acmg_codes.tsv.gz` of Supplementary Data 1 contain what the paper used. Scripts 04 and 20 (text mining
and its audit) need the raw file; script 10 reads it when present and otherwise uses the code table of Supplementary Data 1.

## Licences and data access
Code and the label table in this repository: CC BY 4.0. The label table is derived from ClinVar and follows NCBI's
policies for ClinVar data. dbNSFP 4.8a (the academic release) scores are subject to the licence terms of dbNSFP and of
its component predictors and are not redistributed here; re-create them with `scripts/03_fetch_dbnsfp.py`. ProteinGym
v1.3 is available from its project site. Submitter names: the paper calls Ambry Genetics "Ambry" and Labcorp Genetics
(formerly Invitae) "Labcorp".

## Use of an AI assistant
The analysis code in this repository was written with an AI assistant (Claude, Anthropic) under the authors' direction;
Supplementary Methods S11 of the paper describes how it was used and checked.
