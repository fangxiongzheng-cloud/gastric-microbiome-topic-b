#!/bin/bash
set -euo pipefail
source /REPLACE_WITH_ENVIRONMENT_ROOT/activate.sh
W=/REPLACE_WITH_PRIVATE_WORKDIR/data/external/korea16s; O=/REPLACE_WITH_PRIVATE_WORKDIR/results/korea16s
cd $W/art
q2 qiime feature-table merge $(for f in b*.table.qza; do echo --i-tables $f; done) --o-merged-table ../table.qza
q2 qiime feature-table merge-seqs $(for f in b*.rep.qza; do echo --i-data $f; done) --o-merged-data ../rep.qza
cd $W
for cl in silva:/REPLACE_WITH_ENVIRONMENT_ROOT/db/qiime2/silva-138-99-nb-classifier.qza gg2v4:/REPLACE_WITH_ENVIRONMENT_ROOT/db/qiime2/gg2-2024.09-backbone-v4-nb.qza; do
  n=${cl%%:*}; c=${cl#*:}
  [ -f tax_$n.qza ] || q2 qiime feature-classifier classify-sklearn --i-classifier $c --i-reads rep.qza --p-n-jobs 1 --p-reads-per-batch 2000 --o-classification tax_$n.qza
  : # Candidate: shared QIIME2 temporary cleanup omitted.
  q2 qiime tools export --input-path tax_$n.qza --output-path exp_tax_$n
done
q2 qiime tools export --input-path table.qza --output-path exp_table
q2 biom convert -i exp_table/feature-table.biom -o exp_table/table.tsv --to-tsv
for f in art/b*.stats.qza; do q2 qiime tools export --input-path $f --output-path exp_stats/$(basename $f .qza); done
echo FINALIZED
