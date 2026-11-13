#!/bin/bash
# v0.6e wrapper of krna_pilot.sh for PRJEB9749.
# Download: AWS SRA Open Data Program (official same runs), not ENA FASTQ.
# Reason: ENA ftp.sra.ebi.ac.uk measured ~0.4-0.9 MB/s and stalled;
#         AWS ODP measured ~24 MB/s on this box (2026-11-10). All 30 runs present (~118 GB SRA).
# Ops only: same analysis commands/parameters/tools/refs as krna_pilot.sh
# (fastq-dump | fastp --stdin --interleaved_in, then Salmon/bowtie2/Kraken2).
# No prefetch, disk floor 15G, delete SRA/FASTQ after verified outputs.
# Usage: krna_prjeb9749.sh [runlist.tsv]
set -uo pipefail
export KRAKEN2_DB=/REPLACE_WITH_ENVIRONMENT_ROOT/db/k2_std8
export PATH=/usr/bin:/bin
ROOT=/REPLACE_WITH_PRIVATE_WORKDIR
L=${1:-$ROOT/data/external/rna_prjeb9749.tsv}
O=$ROOT/results/prjeb9749/per_run
W=/REPLACE_WITH_PRIVATE_SCRATCH/scratch_prjeb
mkdir -p "$W" "$O"
HP=/REPLACE_WITH_REFERENCE_ROOT/hp/hp12
VI=/REPLACE_WITH_REFERENCE_ROOT/hp/virul
IDX=/REPLACE_WITH_REFERENCE_ROOT/gencode/salmon_v46
S=/REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin
SAL=/REPLACE_WITH_ENVIRONMENT_ROOT/envs/salmon/bin/salmon

expect_size() {
  # $1 run; print content-length or empty
  /usr/bin/curl -sI -L --max-time 15 "https://sra-pub-run-odp.s3.amazonaws.com/sra/$1/$1" \
    | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r' | tail -1
}

dl() {
  local r=$1 dest="$W/$r.sra" exp="$2" t
  for t in 1 2 3; do
    /usr/bin/curl -fL --connect-timeout 20 --retry 2 --retry-delay 10 --retry-all-errors \
      --speed-time 30 --speed-limit 10000 -C - -o "$dest" \
      "https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r" || continue
    if [ -n "$exp" ] && [ -s "$dest" ]; then
      local got
      got=$(stat -c %s "$dest")
      [ "$got" = "$exp" ] && return 0
      echo "$(date) size mismatch $r got=$got expected=$exp, retry"
      rm -f "$dest"
    elif [ -s "$dest" ]; then
      return 0
    fi
    sleep 10
  done
  return 1
}

outputs_ok() {
  local r=$1
  [ -s "$O/$r/fastp.json" ] && [ -s "$O/$r/salmon/quant.sf.gz" ] && \
  [ -s "$O/$r/salmon/aux_info/meta_info.json" ] && [ -s "$O/$r/hp.bam" ] && \
  [ -s "$O/$r/k2.report" ] && [ -s "$O/$r/k2.out.txt.gz" ] && [ -s "$O/$r/vir.bam" ]
}

for r in $(tail -n +2 "$L" | cut -f1); do
  [ -f "$O/$r/DONE" ] && continue
  [ -f "$ROOT/results/prjeb9749/PAUSE" ] && { echo "$(date) PAUSE"; exit 0; }
  avail=$(df --output=avail -BG / | tail -1 | tr -dc 0-9)
  [ "$avail" -lt 15 ] && { echo "$(date) LOW DISK ${avail}G"; exit 1; }
  exp=$(expect_size "$r")
  echo "$(date) START $r free=${avail}G sra_bytes=${exp:-unknown}"
  mkdir -p "$O/$r"
  cd "$W" || exit 1
  rm -f "$W/$r.sra" "$W/$r.1.fq.gz" "$W/$r.2.fq.gz" "$W/$r.hp1.fq" "$W/$r.hp2.fq"
  if ! dl "$r" "$exp"; then
    echo "$(date) FAIL $r download"
    rm -f "$W/$r.sra"
    continue
  fi
  $S/fastq-dump --split-spot -Z "$W/$r.sra" 2> "$O/$r/fqdump.err" | $S/fastp --stdin --interleaved_in -o "$W/$r.1.fq.gz" -O "$W/$r.2.fq.gz" -z 1 -w 8 -j "$O/$r/fastp.json" -h /dev/null 2> "$O/$r/fastp.log" \
    || { echo "$(date) FAIL $r fastp"; rm -f "$W/$r.sra" "$W/$r.1.fq.gz" "$W/$r.2.fq.gz"; continue; }
  rm -f "$W/$r.sra"
  $SAL quant -i "$IDX" -l A -1 "$W/$r.1.fq.gz" -2 "$W/$r.2.fq.gz" -p 8 --validateMappings -o "$O/$r/salmon" > "$O/$r/salmon.log" 2>&1 && gzip -f "$O/$r/salmon/quant.sf"
  $S/bowtie2 --very-sensitive -p 8 -x "$HP" -1 "$W/$r.1.fq.gz" -2 "$W/$r.2.fq.gz" --no-unal 2> "$O/$r/bt2.log" | $S/samtools view -b -o "$O/$r/hp.bam" -
  $S/samtools fastq -1 "$W/$r.hp1.fq" -2 "$W/$r.hp2.fq" -s /dev/null -0 /dev/null -F 0x900 "$O/$r/hp.bam" 2>/dev/null
  $S/kraken2 --db "$KRAKEN2_DB" --memory-mapping --paired "$W/$r.hp1.fq" "$W/$r.hp2.fq" --threads 4 --report "$O/$r/k2.report" --output "$O/$r/k2.out.txt" 2> "$O/$r/k2.log"
  gzip -f "$O/$r/k2.out.txt"
  $S/bowtie2 --very-sensitive -p 4 -x "$VI" -1 "$W/$r.hp1.fq" -2 "$W/$r.hp2.fq" --no-unal 2> "$O/$r/vir.log" | $S/samtools view -b -o "$O/$r/vir.bam" -
  rm -f "$W/$r.1.fq.gz" "$W/$r.2.fq.gz" "$W/$r.hp1.fq" "$W/$r.hp2.fq"
  if outputs_ok "$r"; then
    touch "$O/$r/DONE"
    echo "$(date) DONE $r"
  else
    echo "$(date) FAIL $r incomplete outputs"
  fi
done
echo ALLDONE
