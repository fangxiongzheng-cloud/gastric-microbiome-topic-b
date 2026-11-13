# v0.6e E12/E13 (EXPLORATORY) on PRJEB9749. Implements plans/addendum_v0.6e_correa.md sections 3-4 exactly.
# QC rule = QC_final (files exist and non-empty). Metric = log1p(Kraken2 taxid210 clade / Salmon num_mapped x1e6).
# JT: J=sum_{i<j} U_ij (earlier>later, ties 0.5); one-sided P=(1+#J*>=Jobs)/(1+10000); label permutation w/o replacement,
# fresh numpy.random.default_rng(20261030) per test. No pairwise tests.
import pandas as pd,numpy as np,json,os,gzip,subprocess,hashlib
from scipy.stats import spearmanr
B='/REPLACE_WITH_PRIVATE_WORKDIR/'; O=B+'results/prjeb9749/'; R=O+'per_run/'; OUT=O+'E12_E13/'; os.makedirs(OUT,exist_ok=True)
SEED=20261030; NPERM=10000; NBOOT=2000
sha=lambda f:hashlib.sha256(open(f,'rb').read()).hexdigest()
FROZEN={'plans/addendum_v0.6e_correa.md':'d37aac254bf8237ad71c970eaed4dc1b5f645323db421b4cd1fdfe34889959f1','data/external/rna_prjeb9749.tsv':'2434ab442fb33647a505fd2ce4677594803888dcf1bcfc917528e2a2eedadd7c'}
for f,h in FROZEN.items(): assert sha(B+f)==h, f+' sha mismatch'
L=pd.read_csv(B+'data/external/rna_prjeb9749.tsv',sep='\t')
def k2clade(f,taxid='210'):
    for l in open(f):
        x=l.rstrip('\n').split('\t')
        if len(x)>4 and x[4].strip()==taxid: return int(x[1])
    return 0
t2g={}
with gzip.open('/REPLACE_WITH_REFERENCE_ROOT/gencode/gencode.v46.transcripts.fa.gz','rt') as f:
    for l in f:
        if l[0]=='>': x=l[1:].split('|'); t2g[x[0]]=x[5]
rows=[];ex={}
for _,r in L.iterrows():
    d=R+r.run_accession+'/'; g=lambda f:os.path.exists(d+f) and os.path.getsize(d+f)>0
    x=dict(run=r.run_accession,source_name=r.source_name,group5=r.group5,group5_order=int(r.group5_order),factor4=r.factor4,factor4_order=int(r.factor4_order),
      DONE=os.path.exists(d+'DONE'),fastp=g('fastp.json'),salmon_quant=g('salmon/quant.sf.gz') or g('salmon/quant.sf'),meta=g('salmon/aux_info/meta_info.json'),
      hp_bam=g('hp.bam'),k2_report=g('k2.report'),k2_out=g('k2.out.txt.gz'),vir_bam=g('vir.bam'))
    x['QC_pass']=all(x[k] for k in ['DONE','fastp','salmon_quant','meta','hp_bam','k2_report','k2_out','vir_bam'])
    if x['fastp']: x['pairs_postQC']=json.load(open(d+'fastp.json'))['summary']['after_filtering']['total_reads']//2
    if x['meta']:
        mi=json.load(open(d+'salmon/aux_info/meta_info.json')); x['salmon_mapped']=mi['num_mapped']; x['salmon_pct']=mi['percent_mapped']
        x['libtype']=json.load(open(d+'salmon/lib_format_counts.json'))['expected_format']
    if g('bt2.log'):
        t=[l for l in open(d+'bt2.log') if 'overall alignment' in l]; x['bt2_overall_pct']=float(t[0].split('%')[0]) if t else None
    if x['k2_report']: x['hp_k2_clade']=k2clade(d+'k2.report')
    if x['hp_bam']: x['hp_bt2_fragments']=int(subprocess.run(f"/REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin/samtools view -F 0x904 {d}hp.bam | cut -f1 | sort -u | wc -l",shell=True,capture_output=True,text=True).stdout)
    x['ena_read_count']=int(r.ena_read_count)
    rows.append(x)
    if x['QC_pass']:
        q=pd.read_csv(d+'salmon/quant.sf.gz',sep='\t',usecols=['Name','TPM']); q['g']=q.Name.map(t2g); ex[r.run_accession]=q.groupby('g').TPM.sum()
D=pd.DataFrame(rows); D['ratio_postQC_vs_ena']=D.pairs_postQC/D.ena_read_count
D['hp_rpm_host']=D.hp_k2_clade/D.salmon_mapped*1e6; D['y']=np.log1p(D.hp_rpm_host)
D['hp_rpm_k2only']=D.hp_k2_clade/D.pairs_postQC*1e6; D['y_k2only']=np.log1p(D.hp_rpm_k2only)
D.to_csv(OUT+'per_sample.csv',index=False)
Q=D[D.QC_pass].copy()
def jt(y,g):
    y=np.asarray(y,float); g=np.asarray(g); lv=np.sort(np.unique(g)); J=0.0
    for a in range(len(lv)):
        for b in range(a+1,len(lv)):
            ya=y[g==lv[a]][:,None]; yb=y[g==lv[b]][None,:]; J+=(ya>yb).sum()+0.5*(ya==yb).sum()
    return J
def jt_test(y,g):
    y=np.asarray(y,float); g=np.asarray(g); ng=len(np.unique(g))
    if ng<3: return dict(note='<3 groups; JT not run')
    rng=np.random.default_rng(SEED); Jo=jt(y,g); cnt=0
    for _ in range(NPERM):
        if jt(y,rng.permutation(g))>=Jo: cnt+=1
    n=pd.Series(g).value_counts(); Jmax=sum(n[a]*n[b] for i,a in enumerate(sorted(n.index)) for b in sorted(n.index)[i+1:])
    return dict(J_obs=float(Jo),J_max=float(Jmax),J_null_mean=float(Jmax/2),J_over_Jmax=float(Jo/Jmax),P_one_sided_decreasing=(1+cnt)/(1+NPERM),n_perm=NPERM,seed=SEED,n=int(len(y)),n_groups=int(ng))
def sp1(y,g):
    r=spearmanr(g,y); P2=float(r.pvalue); rho=float(r.statistic)
    return dict(rho=rho,P_two=P2,P_one_sided_rho_lt0=P2/2 if rho<0 else 1-P2/2,n=int(len(y)))
def gsum(df,col,gcol,ocol):
    out=[]
    for (gname,o),s in df.groupby([gcol,ocol]):
        out.append(dict(group=gname,order=int(o),n=int(len(s)),median_rpm=float(s.hp_rpm_host.median()) if col=='y' else float(s.hp_rpm_k2only.median()),
          median_log1p=float(s[col].median()),IQR_log1p=[float(s[col].quantile(.25)),float(s[col].quantile(.75))],
          IQR_rpm=[float((s.hp_rpm_host if col=='y' else s.hp_rpm_k2only).quantile(.25)),float((s.hp_rpm_host if col=='y' else s.hp_rpm_k2only).quantile(.75))],
          n_k2clade_ge10=int((s.hp_k2_clade>=10).sum()),n_k2clade_zero=int((s.hp_k2_clade==0).sum())))
    return sorted(out,key=lambda z:z['order'])
res=dict(label='EXPLORATORY (v0.6e; permanent; never confirmatory)',frozen_sha256_verified=FROZEN,
  qc=dict(n_runs=len(D),n_QC_pass=int(D.QC_pass.sum()),failed=D.loc[~D.QC_pass,['run','group5']].to_dict('records'),libtype=D.libtype.value_counts().to_dict(),
   salmon_pct_median=float(D.salmon_pct.median()),salmon_pct_range=[float(D.salmon_pct.min()),float(D.salmon_pct.max())],ratio_postQC_vs_ena_range=[float(D.ratio_postQC_vs_ena.min()),float(D.ratio_postQC_vs_ena.max())],
   group5_n_after_QC=Q.group5.value_counts().to_dict()))
pri=jt_test(Q.y,Q.group5_order); pri.update(spearman_effect=sp1(Q.y,Q.group5_order),groups=gsum(Q,'y','group5','group5_order'))
pri['verdict']=('本探索性队列中 HP 读段负荷沿组织学分期递减' if pri['P_one_sided_decreasing']<0.05 else '未观察到预声明方向的趋势')
res['E12_primary']=pri
s1=jt_test(Q.y,Q.factor4_order); s1['groups']=gsum(Q,'y','factor4','factor4_order')
s2=dict(group5=sp1(Q.y,Q.group5_order),factor4=sp1(Q.y,Q.factor4_order))
s3=jt_test(Q.y_k2only,Q.group5_order); s3['groups']=gsum(Q,'y_k2only','group5','group5_order'); s3['metric']='log1p(Kraken2 taxid210 clade / fastp post-QC pairs x1e6)'
q4=Q[Q.run!='ERR950183']; s4=jt_test(q4.y,q4.group5_order)
res['E12_sensitivity']={'note':'parallel report only; cannot rewrite primary; no pass/fail',
  'S1_official_factor4_JT':s1,'S2_spearman_ordinal':s2,'S3_kraken2_only_JT':s3,'S4_drop_EA1_ERR950183_JT':s4}
# E13
IM=['CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4','SI','FABP1']
GD=['GKN1','GKN2','TFF2','PGC','PGA3','LIPF','CHIA','MUC6','CLDN18','GAST']
M2=['ATP4A','ATP4B','CBLIF','KCNE2','KCNJ15','KCNJ16','SLC26A7','CCKBR','ESRRG','VSIG1']
M3=['MUC5AC','TFF1','TFF2','MUC6','MUC1','FUT2','FUT3','B3GNT6']
X=np.log2(pd.DataFrame(ex)+1).T; lg={}
def zscore_mean(X,genes,name):
    pres=[g for g in genes if g in X.columns]; miss=[g for g in genes if g not in pres]; lg[name]=dict(present=len(pres),missing=miss)
    if len(pres)/len(genes)<0.7: lg[name]['FAILED_<70%']=True; return None
    Z=(X[pres]-X[pres].mean())/X[pres].std(ddof=1); return Z.mean(axis=1)
S=pd.DataFrame({'M1':zscore_mean(X,IM,'IM')-zscore_mean(X,GD,'GD'),'M2':zscore_mean(X,M2,'M2'),'M3':zscore_mean(X,M3,'M3')})
S=S.join(Q.set_index('run')[['group5','group5_order']]); S.to_csv(OUT+'E13_scores.csv')
e13={'label':'EXPLORATORY secondary','genes':lg,'expected':{'M1':'increase with stage','M2':'decrease with atrophy/IM','M3':'decrease with stage'}}
ps={}
for k in ['M1','M2','M3']:
    r=spearmanr(S.group5_order,S[k]); rng=np.random.default_rng(SEED); bs=[]; v=S[[k,'group5_order']].values
    for _ in range(NBOOT):
        i=rng.integers(0,len(v),len(v)); rr=spearmanr(v[i,1],v[i,0]).statistic
        if np.isfinite(rr): bs.append(rr)
    e13[k]=dict(rho=float(r.statistic),P_two=float(r.pvalue),boot95CI=[float(np.percentile(bs,2.5)),float(np.percentile(bs,97.5))],n_boot_valid=len(bs),n=int(len(S)),
       group_medians={g:float(s[k].median()) for g,s in S.groupby('group5')}); ps[k]=float(r.pvalue)
o=sorted(ps,key=ps.get); m=len(o); prev=1.0; q={}
for i in range(m-1,-1,-1): prev=min(prev,ps[o[i]]*m/(i+1)); q[o[i]]=prev
for k in ps: e13[k]['q_BH']=q[k]; e13[k]['FDR<0.10_descriptive']=bool(q[k]<0.10)
res['E13']=e13
json.dump(res,open(OUT+'E12_E13_results.json','w'),indent=1,ensure_ascii=False,default=float)
print(json.dumps(res,indent=1,ensure_ascii=False,default=float))
print(D[['run','source_name','group5','QC_pass','pairs_postQC','salmon_pct','libtype','bt2_overall_pct','hp_bt2_fragments','hp_k2_clade','hp_rpm_host','hp_rpm_k2only']].to_string())
