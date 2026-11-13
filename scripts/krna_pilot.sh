#!/bin/bash
# Resumable Korean RNA-seq pilot (protocol v0.4 A1/A8). Usage: krna_pilot.sh runlist.tsv
set -uo pipefail
export KRAKEN2_DB=/REPLACE_WITH_ENVIRONMENT_ROOT/db/k2_std8; export PATH=/usr/bin:/bin
L=${1:-/REPLACE_WITH_PRIVATE_WORKDIR/data/external/rna_pilot40.tsv}
O=/REPLACE_WITH_PRIVATE_WORKDIR/results/korea_rna/per_run; W=/REPLACE_WITH_PRIVATE_SCRATCH/scratch_rna; mkdir -p $W $O
HP=/REPLACE_WITH_REFERENCE_ROOT/hp/hp12; VI=/REPLACE_WITH_REFERENCE_ROOT/hp/virul; IDX=/REPLACE_WITH_REFERENCE_ROOT/gencode/salmon_v46
S=/REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin; SAL=/REPLACE_WITH_ENVIRONMENT_ROOT/envs/salmon/bin/salmon
dl(){ r=$1; [ -e $W/$r.sra.ok ] && return; for t in 1 2 3; do /usr/bin/curl -sf -o $W/$r.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r && touch $W/$r.sra.ok && return; done; }
runs=$(tail -n +2 $L | cut -f1)
prev=""
for r in $runs; do
  [ -f $O/$r/DONE ] && continue
  avail=$(df --output=avail -BG / | tail -1 | tr -dc 0-9); [ "$avail" -lt 12 ] && { echo "$(date) LOW DISK $avail"; exit 1; }
  echo "$(date) START $r free=${avail}G"; mkdir -p $O/$r; cd $W
  dl $r
  # prefetch next run in background
  nxt=$(echo "$runs" | grep -A1 -x $r | tail -1); [ "$nxt" != "$r" ] && [ ! -f $O/$nxt/DONE ] && (dl $nxt &)
  $S/fastq-dump --split-spot -Z $W/$r.sra 2>/dev/null | $S/fastp --stdin --interleaved_in -o $W/$r.1.fq.gz -O $W/$r.2.fq.gz -z 1 -w 8 -j $O/$r/fastp.json -h /dev/null 2> $O/$r/fastp.log || { echo "fastp fail $r"; continue; }
  rm -f $W/$r.sra $W/$r.sra.ok
  $SAL quant -i $IDX -l A -1 $W/$r.1.fq.gz -2 $W/$r.2.fq.gz -p 8 --validateMappings -o $O/$r/salmon > $O/$r/salmon.log 2>&1 && gzip -f $O/$r/salmon/quant.sf
  $S/bowtie2 --very-sensitive -p 8 -x $HP -1 $W/$r.1.fq.gz -2 $W/$r.2.fq.gz --no-unal 2> $O/$r/bt2.log | $S/samtools view -b -o $O/$r/hp.bam -
  $S/samtools fastq -1 $W/$r.hp1.fq -2 $W/$r.hp2.fq -s /dev/null -0 /dev/null -F 0x900 $O/$r/hp.bam 2>/dev/null
  $S/kraken2 --db $KRAKEN2_DB --memory-mapping --paired $W/$r.hp1.fq $W/$r.hp2.fq --threads 4 --report $O/$r/k2.report --output $O/$r/k2.out.txt 2> $O/$r/k2.log; gzip -f $O/$r/k2.out.txt
  $S/bowtie2 --very-sensitive -p 4 -x $VI -1 $W/$r.hp1.fq -2 $W/$r.hp2.fq --no-unal 2> $O/$r/vir.log | $S/samtools view -b -o $O/$r/vir.bam -
  rm -f $W/$r.1.fq.gz $W/$r.2.fq.gz $W/$r.hp1.fq $W/$r.hp2.fq
  touch $O/$r/DONE; echo "$(date) DONE $r"
done
echo ALLDONE
