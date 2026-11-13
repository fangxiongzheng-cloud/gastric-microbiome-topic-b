#!/bin/bash
# Batch download -> cutadapt -> DADA2 per batch; raw fastq deleted after artifacts exist
set -uo pipefail
source /REPLACE_WITH_ENVIRONMENT_ROOT/activate.sh
W=/REPLACE_WITH_PRIVATE_WORKDIR/data/external/korea16s
cd $W; mkdir -p batches art log
[ -d batches ] || split -l 24 -d runs16s.tsv batches/b
for B in batches/b*; do
  b=$(basename $B); [ "$b" \< "b02" ] && continue
  [ -f art/$b.table.qza ] && continue
  while [ $(ls art/*.trim.qza 2>/dev/null | wc -l) -ge 2 ]; do sleep 30; done
  avail=$(df --output=avail -BG / | tail -1 | tr -dc 0-9); echo "$(date) prod $b free=${avail}G"
  if [ "$avail" -lt 12 ]; then echo "LOW DISK, stop"; exit 1; fi
  if [ ! -f art/$b.trim.qza ]; then
  F=fq_$b; mkdir -p $F
  printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > $F/manifest.tsv
  while IFS=$'\t' read run s lib st ins rc ftp bytes; do
    u1=${ftp%%;*}; u2=${ftp##*;}
    echo "$run https://$u1 https://$u2"
    printf "%s\t%s\t%s\n" $run $W/$F/${run}_1.fastq.gz $W/$F/${run}_2.fastq.gz >> $F/manifest.tsv
  done < $B > $F/urls.txt
  cut -d' ' -f1 $F/urls.txt | xargs -P 8 -I{} bash -c 'cd '$W/$F' && for t in 1 2 3; do curl -sf -o {}.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/{}/{} && /REPLACE_WITH_LEGACY_ENVIRONMENT_ROOT/bin/micromamba run -p /REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun fasterq-dump --split-files -e 2 {}.sra -O . >/dev/null 2>&1 && break; done; rm -f {}.sra; mv {}.sra_1.fastq {}_1.fastq 2>/dev/null; mv {}.sra_2.fastq {}_2.fastq 2>/dev/null; gzip -1 {}_1.fastq {}_2.fastq'
  q2 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $F/manifest.tsv --output-path art/$b.demux.qza --input-format PairedEndFastqManifestPhred33V2 || { echo "import fail $b"; exit 1; }
  rm -rf $F
  q2 qiime cutadapt trim-paired --i-demultiplexed-sequences art/$b.demux.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACHVGGGTATCTAATCC --p-discard-untrimmed --p-cores 8 --o-trimmed-sequences art/$b.trimtmp.qza && rm art/$b.demux.qza
  mv art/$b.trimtmp.qza art/$b.trim.qza
  fi
done
echo PRODDONE
