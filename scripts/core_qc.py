import pandas as pd,os,json,gzip
B='/REPLACE_WITH_PRIVATE_WORKDIR/';O=B+'results/korea_rna/per_run/'
L=pd.read_csv(B+'data/external/rna_full352.tsv',sep='\t')
rows=[]
for r in L.run_accession:
    d=O+r+'/'; g=lambda f:os.path.exists(d+f) and os.path.getsize(d+f)>0
    x=dict(run=r,DONE=os.path.exists(d+'DONE'),fastp=g('fastp.json'),salmon_quant=g('salmon/quant.sf.gz') or g('salmon/quant.sf'),
      meta=g('salmon/aux_info/meta_info.json'),hp_bam=g('hp.bam'),k2_report=os.path.exists(d+'k2.report'),k2_out=g('k2.out.txt.gz') or os.path.exists(d+'k2.out.txt'),vir_bam=g('vir.bam'),bt2log=g('bt2.log'))
    x['k2_report_nonempty']=g('k2.report')
    if x['meta']:
        mi=json.load(open(d+'salmon/aux_info/meta_info.json')); x['salmon_mapped']=mi['num_mapped'];x['salmon_pct']=mi['percent_mapped']
    if x['fastp']:
        x['pairs_postQC']=json.load(open(d+'fastp.json'))['summary']['after_filtering']['total_reads']//2
    if g('bt2.log'):
        t=open(d+'bt2.log').read().split('\n'); x['bt2_overall']=[l for l in t if 'overall alignment' in l][0].split('%')[0].strip() if any('overall alignment' in l for l in t) else None
    rows.append(x)
D=pd.DataFrame(rows);D=D.merge(L[['run_accession','library_name','read_count']],left_on='run',right_on='run_accession').drop(columns='run_accession')
D['ratio_postQC_vs_ena']=D.pairs_postQC/D.read_count
D.to_csv(B+'results/korea_rna/QC_final_per_run.csv',index=False)
print(len(D));print(D[['DONE','fastp','salmon_quant','meta','hp_bam','k2_report','k2_report_nonempty','k2_out','vir_bam','bt2log']].sum())
print(D.salmon_pct.describe(),D.ratio_postQC_vs_ena.describe())
print(D[~(D.DONE&D.fastp&D.salmon_quant&D.meta&D.hp_bam&D.k2_report&D.vir_bam)].run.tolist())
print(D[D.k2_report&~D.k2_report_nonempty].run.tolist())
print(D.library_name.str[-1].value_counts(), D.library_name.duplicated().sum())
