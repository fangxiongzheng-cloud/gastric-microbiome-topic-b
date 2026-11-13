#!/bin/bash
set -uo pipefail
source /REPLACE_WITH_ENVIRONMENT_ROOT/activate.sh
B=/REPLACE_WITH_PRIVATE_WORKDIR; W=$B/data/external/third_cohort/PRJNA1061213; A=$W/art; mkdir -p $A
cd $W
python3 - <<'PY'
import pandas as pd
d=pd.read_csv('/REPLACE_WITH_PRIVATE_WORKDIR/data/external/third_cohort/PRJNA1061213_runs.tsv',sep='\t')
d=d[d.sample_alias.str.match(r'^(CT|NT)\d+')]
with open('/REPLACE_WITH_PRIVATE_WORKDIR/data/external/third_cohort/PRJNA1061213/manifest.tsv','w') as f:
    f.write('sample-id\tabsolute-filepath\n')
    for r in d.run_accession: f.write(f'{r}\t/REPLACE_WITH_PRIVATE_WORKDIR/data/external/third_cohort/PRJNA1061213/raw/{r}.fastq.gz\n')
print(len(d))
PY
echo "$(date) import"
[ -f $A/demux.qza ] || q2 qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path manifest.tsv --output-path $A/demux.qza --input-format SingleEndFastqManifestPhred33V2 || exit 1
echo "$(date) dada2"
[ -f $A/table.qza ] || q2 qiime dada2 denoise-single --i-demultiplexed-seqs $A/demux.qza --p-trunc-len 0 --p-max-ee 2.0 --p-n-threads 3 --o-table $A/table.qza --o-representative-sequences $A/rep.qza --o-denoising-stats $A/stats.qza --o-base-transition-stats $A/bts.qza || exit 1
for cl in silva:/REPLACE_WITH_ENVIRONMENT_ROOT/db/qiime2/silva-138-99-nb-classifier.qza gg2v4:/REPLACE_WITH_ENVIRONMENT_ROOT/db/qiime2/gg2-2024.09-backbone-v4-nb.qza; do
  n=${cl%%:*}; c=${cl#*:}; echo "$(date) classify $n"
  [ -f $A/tax_$n.qza ] || q2 qiime feature-classifier classify-sklearn --i-classifier $c --i-reads $A/rep.qza --p-n-jobs 1 --p-reads-per-batch 2000 --o-classification $A/tax_$n.qza || exit 1
  q2 qiime tools export --input-path $A/tax_$n.qza --output-path $A/exp_tax_$n
done
q2 qiime tools export --input-path $A/table.qza --output-path $A/exp_table
q2 biom convert -i $A/exp_table/feature-table.biom -o $A/exp_table/table.tsv --to-tsv
q2 qiime tools export --input-path $A/stats.qza --output-path $A/exp_stats
echo F2PIPEDONE $(date)
