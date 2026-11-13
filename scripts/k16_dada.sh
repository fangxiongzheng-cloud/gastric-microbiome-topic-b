#!/bin/bash
# consumer: DADA2 on trims produced by k16_prod.sh, then finalize+analysis
set -uo pipefail
source /REPLACE_WITH_ENVIRONMENT_ROOT/activate.sh
cd /REPLACE_WITH_PRIVATE_WORKDIR/data/external/korea16s
for B in batches/b*; do
  b=$(basename $B)
  while pgrep -f "dada2 denoise-paired" >/dev/null; do sleep 30; done
  [ -f art/$b.table.qza ] && continue
  while [ ! -f art/$b.trim.qza ]; do pgrep -f k16_prod.sh >/dev/null || { echo "prod gone, no trim $b"; exit 1; }; sleep 30; done
  echo "$(date) dada2 $b"
  q2 qiime dada2 denoise-paired --i-demultiplexed-seqs art/$b.trim.qza --p-trunc-len-f 260 --p-trunc-len-r 220 --p-n-threads 4 --o-table art/$b.table.qza --o-representative-sequences art/$b.rep.qza --o-denoising-stats art/$b.stats.qza --o-base-transition-stats art/$b.bts.qza || { echo "dada2 fail $b"; exit 1; }
  rm art/$b.trim.qza; : # Candidate: shared QIIME2 temporary cleanup omitted.
done
cd /REPLACE_WITH_PRIVATE_WORKDIR
bash scripts/k16_finalize.sh || exit 1
python scripts/k16_genus.py && python scripts/k16_analysis.py && echo PIPEDONE
