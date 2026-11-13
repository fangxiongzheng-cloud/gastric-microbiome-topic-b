#!/usr/bin/env python3
"""Summarise per-run outputs of hp_hostrna_quant.sh into a table.
Metric: sequencing-depth-normalised H. pylori transcript-read load (RPM_host) =
  Kraken2-confirmed fragments in the taxid-210 clade / Salmon-mapped host fragments x 1e6.
Usage: summarize_hp_load.py <outdir> <out.csv> [--samtools PATH]"""
import sys,os,json,glob,subprocess,math,argparse
ap=argparse.ArgumentParser(); ap.add_argument('outdir'); ap.add_argument('out'); ap.add_argument('--samtools',default='samtools'); a=ap.parse_args()
def k2clade(f,taxid='210'):
    for l in open(f):
        x=l.rstrip('\n').split('\t')
        if len(x)>4 and x[4]==taxid: return int(x[1])
    return 0
rows=[]
for d in sorted(glob.glob(os.path.join(a.outdir,'*',''))):
    if not os.path.exists(d+'DONE'): continue
    r=os.path.basename(d.rstrip('/'))
    fp=json.load(open(d+'fastp.json'))['summary']['after_filtering']['total_reads']//2
    mi=json.load(open(d+'salmon/aux_info/meta_info.json'))
    bt=int(subprocess.run(f"{a.samtools} view -F 0x904 {d}hp.bam | cut -f1 | sort -u | wc -l",shell=True,capture_output=True,text=True).stdout)
    k=k2clade(d+'k2.report'); rpm=k/mi['num_mapped']*1e6
    rows.append([r,fp,mi['num_mapped'],round(mi['percent_mapped'],3),bt,k,rpm,math.log1p(rpm),int(k>=10)])
with open(a.out,'w') as f:
    f.write('run,pairs_postQC,salmon_mapped,salmon_pct,hp_bt2_fragments,hp_k2_confirmed,hp_rpm_host,hp_log1p_rpm,detected_ge10\n')
    for x in rows: f.write(','.join(map(str,x))+'\n')
print(f'{len(rows)} runs written to {a.out}')
