#!/bin/bash
# Step10 E8 spike-in simulation runner (v0.6b E8). Detached; waits for E7 DONE. See freeze_log E8 entry.
set -uo pipefail
B=/REPLACE_WITH_PRIVATE_WORKDIR; O=$B/results/step10/E8; W=/REPLACE_WITH_PRIVATE_SCRATCH/scratch_e8; mkdir -p $W $O/bg $O/sim $O/log
export KRAKEN2_DB=/REPLACE_WITH_ENVIRONMENT_ROOT/db/k2_std8; S=/REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin; SAL=/REPLACE_WITH_ENVIRONMENT_ROOT/envs/salmon/bin/salmon
HP=/REPLACE_WITH_REFERENCE_ROOT/hp/hp12; IDX=/REPLACE_WITH_REFERENCE_ROOT/gencode/salmon_v46; R=$O/refs; TH=6
log(){ echo "$(date '+%F %T') $*" | tee -a $O/log/run.log; }
while [ ! -f $B/results/step10/E7/DONE ]; do
  if grep -q "RUNNER END" $B/results/step10/E7/log/status.log 2>/dev/null; then log "E7 runner ended without DONE -> E8 starts anyway"; break; fi
  sleep 300; done
log "E8 START"
quant(){ # $1=r1 $2=r2 $3=outdir ; bowtie2->samtools fastq->kraken2
  mkdir -p $3; $S/bowtie2 --very-sensitive -p $TH -x $HP -1 $1 -2 $2 --no-unal 2> $3/bt2.log | $S/samtools view -b -o $3/hp.bam -
  $S/samtools fastq -1 $3/h1.fq -2 $3/h2.fq -s /dev/null -0 /dev/null -F 0x900 $3/hp.bam 2>/dev/null
  $S/kraken2 --db $KRAKEN2_DB --memory-mapping --paired $3/h1.fq $3/h2.fq --threads 4 --report $3/k2.report --output $3/k2.out.txt 2> $3/k2.log
  gzip -f $3/k2.out.txt; rm -f $3/h1.fq $3/h2.fq; }
# ---------- simulations (small; independent of backgrounds)
python3 - <<PY
import os
seq=''.join(l.strip() for l in open('/REPLACE_WITH_REFERENCE_ROOT/hp/NC_000915.1.fa') if not l.startswith('>'))
reg=[(445337,448576),(1207581,1209081),(1473565,1476805),(1511135,1512634)]
with open('$R/rRNA_regions_26695.fa','w') as f:
    for a,b in reg: f.write(f'>rRNA_{a}_{b}\n'+seq[a-1-150:b+150]+'\n')
PY
declare -A REF=([T1_26695_uniform]=/REPLACE_WITH_REFERENCE_ROOT/hp/NC_000915.1.fa [T3_F30_offref]=$R/NC_017365.1.fa [T4a_acinonychis]=$R/NC_008229.1.fa [T4b_hepaticus]=$R/NC_004917.1.fa [T2_26695_rRNA90]=/REPLACE_WITH_REFERENCE_ROOT/hp/NC_000915.1.fa)
LEV="10 30 100 300 1000 3000 10000 30000"
sim1(){ # target ref n seed out1 out2 tag
  /REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin/wgsim -N $3 -1 151 -2 151 -d 300 -s 50 -e 0.002 -r 0.001 -R 0 -S $4 $2 $5.tmp1 $6.tmp2 > /dev/null 2>&1
  awk -v t=$7 'NR%4==1{sub(/^@/,"@"t"|")}1' $5.tmp1 >> $5; awk -v t=$7 'NR%4==1{sub(/^@/,"@"t"|")}1' $6.tmp2 >> $6; rm -f $5.tmp1 $6.tmp2; }
for T in T1_26695_uniform T2_26695_rRNA90 T3_F30_offref T4a_acinonychis T4b_hepaticus; do for SD in 20261030 20261031 20261032; do
  d=$O/sim/${T}_$SD; [ -f $d/DONE ] && continue; mkdir -p $d; rm -f $W/s1.fq $W/s2.fq
  for n in $LEV; do
    if [ $T = T2_26695_rRNA90 ]; then nr=$(python3 -c "print(round(0.9*$n))"); ng=$((n-nr))
      sim1 $T $R/rRNA_regions_26695.fa $nr $SD $W/s1.fq $W/s2.fq "L${n}_S${SD}_r"; [ $ng -gt 0 ] && sim1 $T ${REF[$T]} $ng $SD $W/s1.fq $W/s2.fq "L${n}_S${SD}_g"
    else sim1 $T ${REF[$T]} $n $SD $W/s1.fq $W/s2.fq "L${n}_S${SD}_g"; fi
  done
  grep -c "^@L" $W/s1.fq > $d/n_reads.txt; quant $W/s1.fq $W/s2.fq $d
  $S/samtools view -F 0x904 $d/hp.bam | cut -f1,5 > $d/bt2_primary_name_mapq.tsv; touch $d/DONE; log "sim $T $SD done"
done; done
# ---------- backgrounds
for r in $(tail -n +2 $O/backgrounds.csv | cut -d, -f1); do
  d=$O/bg/$r; [ -f $d/DONE ] && continue; mkdir -p $d
  av=$(df --output=avail -BG / | tail -1 | tr -dc 0-9); [ $av -lt 12 ] && { log "LOW DISK $av G -> stop"; exit 1; }
  log "bg $r START free=${av}G"
  for t in 1 2 3 4 5; do /usr/bin/curl -sf --connect-timeout 20 --speed-time 60 --speed-limit 10000 --retry 3 --retry-all-errors -o $W/$r.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r && break; log "download retry $t"; done
  $S/fastq-dump --split-spot -Z $W/$r.sra 2>/dev/null | $S/fastp --stdin --interleaved_in -o $W/$r.1.fq.gz -O $W/$r.2.fq.gz -z 1 -w $TH -j $d/fastp.json -h /dev/null 2> $d/fastp.log || { log "fastp fail $r"; continue; }
  rm -f $W/$r.sra
  zcat $W/$r.1.fq.gz | head -n 40000000 > $W/$r.b1.fq; zcat $W/$r.2.fq.gz | head -n 40000000 > $W/$r.b2.fq; rm -f $W/$r.1.fq.gz $W/$r.2.fq.gz
  echo $(( $(wc -l < $W/$r.b1.fq) / 4 )) > $d/n_pairs.txt
  $SAL quant -i $IDX -l A -1 $W/$r.b1.fq -2 $W/$r.b2.fq -p $TH --validateMappings -o $d/salmon > $d/salmon.log 2>&1 && gzip -f $d/salmon/quant.sf
  quant $W/$r.b1.fq $W/$r.b2.fq $d
  if [ $r = SRR31048547 ] && [ ! -f $O/equivalence_check.txt ]; then
    rm -f $W/e1.fq $W/e2.fq; sim1 x /REPLACE_WITH_REFERENCE_ROOT/hp/NC_000915.1.fa 1000 20261030 $W/e1.fq $W/e2.fq "EQ"
    cp $W/e1.fq $W/eq_sim1.fq; cp $W/e2.fq $W/eq_sim2.fq; cat $W/$r.b1.fq >> $W/e1.fq; cat $W/$r.b2.fq >> $W/e2.fq
    quant $W/e1.fq $W/e2.fq $O/eq_merged; quant $W/eq_sim1.fq $W/eq_sim2.fq $O/eq_simonly
    c(){ awk -F'\t' '$5=="210"{print $2}' $1/k2.report | head -1; }
    m=$(c $O/eq_merged); a=$(c $d); s=$(c $O/eq_simonly); a=${a:-0}; s=${s:-0}; m=${m:-0}
    echo "merged=$m bg=$a sim=$s sum=$((a+s))" > $O/equivalence_check.txt; log "EQUIVALENCE $(cat $O/equivalence_check.txt)"
    rm -f $W/e1.fq $W/e2.fq $W/eq_sim1.fq $W/eq_sim2.fq
    if [ "$m" != "$((a+s))" ]; then touch $O/BLOCKED_equivalence; log "EQUIVALENCE FAILED -> stop"; exit 2; fi
  fi
  rm -f $W/$r.b1.fq $W/$r.b2.fq; touch $d/DONE; log "bg $r DONE"
done
log "E8 processing ALLDONE"; touch $O/PROCESSING_DONE
source /REPLACE_WITH_ENVIRONMENT_ROOT/activate.sh; python $B/scripts/s10_e8_analyze.py > $O/log/analyze.log 2>&1 && touch $O/DONE; log "E8 analyze exit $?"
