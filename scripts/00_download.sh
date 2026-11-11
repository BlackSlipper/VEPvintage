#!/bin/bash
P=${VEP_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}
mkdir -p $P/data/clinvar_tab $P/data/dms $P/data/snap_parsed $P/logs
cd $P/data/clinvar_tab
B=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive
cat $P/data/clinvar_snapshots_used.txt | xargs -P 4 -I{} sh -c "f=\$(basename {}); [ -s \$f ] || curl -s --retry 5 -o \$f $B/{}"
# the November 2026 release is the "current" snapshot of the paper
[ -s variant_summary_current.txt.gz ] || curl -s --retry 5 -o variant_summary_current.txt.gz $B/variant_summary_2026-11.txt.gz
# The tab_delimited archive does not keep past submission_summary files. The paper used the November 2026 file (submission_summary_2026-11.txt.gz);
# a later download differs, and data/scv_textmining.tsv.gz and data/scv_acmg_codes.tsv.gz hold what the paper used.
[ -s submission_summary_2026-11.txt.gz ] || { curl -s --retry 5 -o submission_summary_current.txt.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz; echo "NOTE: submission_summary_current.txt.gz is not the November 2026 file used in the paper"; }
cd $P/data/dms
curl -s --retry 5 -O https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip
curl -s --retry 5 -O https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_substitutions.csv
echo DONE > $P/logs/download.done
