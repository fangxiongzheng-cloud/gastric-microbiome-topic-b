#!/bin/bash
# PRJNA310127 China (CC#N/CC#T) — same QIIME2/DADA2 pipeline as Korea; primers Fadrosh 2014 319F/806R; trunc 260/190 fixed from quality/length profile before denoising
set -uo pipefail; export KRAKEN2_DB=/REPLACE_WITH_ENVIRONMENT_ROOT/db/k2_std8
E=/REPLACE_WITH_ENVIRONMENT_ROOT/envs/qiime2; q2(){ env PATH=$E/bin:/usr/bin:/bin CONDA_PREFIX=$E R_HOME=$E/lib/R "$@"; }
W=/REPLACE_WITH_PRIVATE_WORKDIR/data/external/yu310127; cd $W
[ -f demux.qza ] || q2 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.tsv --output-path demux.qza --input-format PairedEndFastqManifestPhred33V2 || exit 1
[ -f trim.qza ] || q2 qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza --p-front-f ACTCCTACGGGAGGCAGCAG --p-front-r GGACTACHVGGGTWTCTAAT --p-discard-untrimmed --p-cores 3 --o-trimmed-sequences trim.qza || exit 1
[ -f table.qza ] || q2 qiime dada2 denoise-paired --i-demultiplexed-seqs trim.qza --p-trunc-len-f 260 --p-trunc-len-r 190 --p-n-threads 3 --o-table table.qza --o-representative-sequences rep.qza --o-denoising-stats stats.qza --o-base-transition-stats bts.qza --verbose || exit 1
for cl in silva:/REPLACE_WITH_ENVIRONMENT_ROOT/db/qiime2/silva-138-99-nb-classifier.qza gg2v4:/REPLACE_WITH_ENVIRONMENT_ROOT/db/qiime2/gg2-2024.09-backbone-v4-nb.qza; do
  n=${cl%%:*}; c=${cl#*:}
  [ -f tax_$n.qza ] || q2 qiime feature-classifier classify-sklearn --i-classifier $c --i-reads rep.qza --p-n-jobs 1 --p-reads-per-batch 2000 --o-classification tax_$n.qza
  q2 qiime tools export --input-path tax_$n.qza --output-path exp_tax_$n
done
q2 qiime tools export --input-path table.qza --output-path exp_table
q2 biom convert -i exp_table/feature-table.biom -o exp_table/table.tsv --to-tsv
q2 qiime tools export --input-path stats.qza --output-path exp_stats
echo FINALIZED
