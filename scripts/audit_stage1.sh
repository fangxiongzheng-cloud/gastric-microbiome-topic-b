#!/bin/bash
# v0.6 §3 measurement audit, stage 1 (no re-download; reads only from existing per_run hp.bam).
# Deviation (implementation, user-approved t807u 2026-11-10): HP reads re-aligned to single genome 26695 (NC_000915.1) for MAPQ/breadth.
# Literal v0.6 MAPQ>=30 on multi-genome hp.bam is computed in stage 2 from the original hp.bam.
set -uo pipefail
B=/REPLACE_WITH_PRIVATE_WORKDIR; P=$B/results/korea_rna/per_run; A=$B/results/audit
S=/REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin; SAL=/REPLACE_WITH_ENVIRONMENT_ROOT/envs/salmon/bin/salmon; R=/REPLACE_WITH_REFERENCE_ROOT/hp/single26695/hp26695
[ -e $R.1.bt2 ] || $S/bowtie2-build -q --threads 8 /REPLACE_WITH_REFERENCE_ROOT/hp/NC_000915.1.fa $R
: > $A/all/hp_R1.fq; : > $A/all/hp_R2.fq
for d in $P/SRR*; do r=$(basename $d); o=$A/per_run/$r; mkdir -p $o
  [ -e $o/OK ] && { $S/samtools fastq -F 0x900 -1 >(cat >> $A/all/hp_R1.fq) -2 >(cat >> $A/all/hp_R2.fq) -0 /dev/null -s /dev/null $d/hp.bam 2>/dev/null; continue; }
  $S/samtools fastq -F 0x900 -1 $o/r1.fq -2 $o/r2.fq -0 /dev/null -s /dev/null $d/hp.bam 2>/dev/null
  cat $o/r1.fq >> $A/all/hp_R1.fq; cat $o/r2.fq >> $A/all/hp_R2.fq
  $S/bowtie2 --very-sensitive -p 8 -x $R -1 $o/r1.fq -2 $o/r2.fq --no-unal 2> $o/bt2_26695.log | $S/samtools sort -@2 -o $o/hp26695.bam - 2>/dev/null
  $S/samtools index $o/hp26695.bam; rm -f $o/r1.fq $o/r2.fq; touch $o/OK
done
echo "STAGE1 realign done $(date)"
$SAL quant -i /REPLACE_WITH_REFERENCE_ROOT/gencode/salmon_v46 -l A -1 $A/all/hp_R1.fq -2 $A/all/hp_R2.fq -p 8 --validateMappings --writeMappings=$A/all/human_map.sam -o $A/all/salmon_hpreads > $A/log/salmon_hpreads.log 2>&1
cut -f1 $A/all/human_map.sam | grep -v '^@' | sort -u > $A/all/human_hit_readnames.txt; rm -f $A/all/human_map.sam
echo "STAGE1 ALLDONE $(date)"
