#!/usr/bin/env bash
# hp-hostrna: sequencing-depth-normalised H. pylori transcript-read load from host bulk RNA-seq.
# Packaged from the exact per-run commands of scripts/krna_pilot.sh used on 352 Korean runs (PRJNA1174735).
# Usage: hp_hostrna_quant.sh <config.sh>      (see config.example.sh)
# Resumable: a run with <outdir>/<run>/DONE is skipped.
set -uo pipefail
CFG=${1:?usage: hp_hostrna_quant.sh config.sh}; source "$CFG"
: "${RUNLIST:?}" "${OUTDIR:?}" "${WORKDIR:?}" "${SALMON_INDEX:?}" "${HP_BT2_INDEX:?}" "${KRAKEN2_DB:?}"
THREADS=${THREADS:-8}; K2_THREADS=${K2_THREADS:-4}; MIN_FREE_GB=${MIN_FREE_GB:-15}; VIR_BT2_INDEX=${VIR_BT2_INDEX:-}
mkdir -p "$OUTDIR" "$WORKDIR"
dl(){ local r=$1; for t in 1 2 3; do curl -sf --connect-timeout 20 --speed-time 30 --speed-limit 10000 --retry-all-errors \
      -o "$WORKDIR/$r.sra" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r" && return 0; done; return 1; }
for r in $(grep -v '^#' "$RUNLIST" | cut -f1 | grep -v '^run'); do
  o=$OUTDIR/$r; [ -f "$o/DONE" ] && continue; mkdir -p "$o"
  avail=$(df --output=avail -BG "$WORKDIR" | tail -1 | tr -dc 0-9); [ "$avail" -lt "$MIN_FREE_GB" ] && { echo "$(date) LOW DISK ${avail}G"; exit 1; }
  echo "$(date) START $r"; dl "$r" || { echo "download failed $r"; continue; }
  fastq-dump --split-spot -Z "$WORKDIR/$r.sra" 2>/dev/null | fastp --stdin --interleaved_in -o "$WORKDIR/$r.1.fq.gz" -O "$WORKDIR/$r.2.fq.gz" \
      -z 1 -w "$THREADS" -j "$o/fastp.json" -h /dev/null 2> "$o/fastp.log" || { echo "fastp fail $r"; continue; }
  rm -f "$WORKDIR/$r.sra"
  salmon quant -i "$SALMON_INDEX" -l A -1 "$WORKDIR/$r.1.fq.gz" -2 "$WORKDIR/$r.2.fq.gz" -p "$THREADS" --validateMappings -o "$o/salmon" > "$o/salmon.log" 2>&1 && gzip -f "$o/salmon/quant.sf"
  bowtie2 --very-sensitive -p "$THREADS" -x "$HP_BT2_INDEX" -1 "$WORKDIR/$r.1.fq.gz" -2 "$WORKDIR/$r.2.fq.gz" --no-unal 2> "$o/bt2.log" | samtools view -b -o "$o/hp.bam" -
  samtools fastq -1 "$WORKDIR/$r.hp1.fq" -2 "$WORKDIR/$r.hp2.fq" -s /dev/null -0 /dev/null -F 0x900 "$o/hp.bam" 2>/dev/null
  kraken2 --db "$KRAKEN2_DB" --memory-mapping --paired "$WORKDIR/$r.hp1.fq" "$WORKDIR/$r.hp2.fq" --threads "$K2_THREADS" \
      --report "$o/k2.report" --output "$o/k2.out.txt" 2> "$o/k2.log"; gzip -f "$o/k2.out.txt"
  if [ -n "$VIR_BT2_INDEX" ]; then bowtie2 --very-sensitive -p "$K2_THREADS" -x "$VIR_BT2_INDEX" -1 "$WORKDIR/$r.hp1.fq" -2 "$WORKDIR/$r.hp2.fq" --no-unal 2> "$o/vir.log" | samtools view -b -o "$o/vir.bam" -; fi
  rm -f "$WORKDIR/$r.1.fq.gz" "$WORKDIR/$r.2.fq.gz" "$WORKDIR/$r.hp1.fq" "$WORKDIR/$r.hp2.fq"
  touch "$o/DONE"; echo "$(date) DONE $r"
done
echo ALLDONE
