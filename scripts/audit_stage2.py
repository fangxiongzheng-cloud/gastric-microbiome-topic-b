# v0.6 §3 measurement audit, stage 2 (sensitivity, NOT confirmatory; no multiplicity correction).
# Primary criterion (frozen v0.6): main effect direction (N>T) does not reverse; P<0.05 not required.
# Implementation choices fixed before looking at audit results (2026-11-10):
#  - every subset count = Kraken2-confirmed (taxid 210 clade) fragments that ALSO pass the filter; denominator = Salmon-mapped host fragments (same as A1)
#  - literal MAPQ>=30: max MAPQ over primary mapped mates in original multi-genome (12-genome) hp.bam
#  - deviation (user-approved t807u): MAPQ>=30 after re-alignment to single genome 26695 (NC_000915.1); labelled 'deviation from protocol (implementation)'
#  - unique/high-confidence (literal): some primary mate with AS>XS or no XS in hp.bam; (deviation): same rule in 26695 re-alignment
#  - low complexity: prinseq-style DUST, 64-bp windows step 32, scaled 0-100, flagged if max window score >7 in either mate
#  - human homology: fragment has any Salmon selective-alignment mapping to GENCODE v46 transcriptome (index salmon_v46)
#  - breadth/evenness on 26695 (deviation alignment, MAPQ>=30), runs with >=100 such fragments
#  - (added 2026-11-10 10:05 after stage-2 crash, before any audit result was seen) 18 runs have no retained k2.out.txt.gz (read-level Kraken2 output); for them read-level subset counts are NaN and those patient pairs drop out of read-level subsets; 'main restricted' reproduces A1 on the same available pairs for comparability
#  - batch: instrument_model from ENA; lane not recoverable (FASTQ deleted, SRA read names) -> reported as unavailable
import pandas as pd,numpy as np,json,subprocess,os,gzip,collections
from scipy.stats import wilcoxon,spearmanr,mannwhitneyu
import statsmodels.formula.api as smf
B='/REPLACE_WITH_PRIVATE_WORKDIR/'; P=B+'results/korea_rna/per_run/'; A=B+'results/audit/'; ST='/REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin/samtools'
D=pd.read_csv(B+'results/core/korea_rna_per_sample.csv')
human=set(open(A+'all/human_hit_readnames.txt').read().split()) if os.path.exists(A+'all/human_hit_readnames.txt') else None
def subtree(rep):
    tx=set(); lvl=None
    for l in open(rep):
        x=l.rstrip('\n').split('\t')
        if len(x)<6: continue
        ind=len(x[5])-len(x[5].lstrip(' '))
        if lvl is not None:
            if ind<=lvl: break
            tx.add(x[4])
        if x[4]=='210': lvl=ind; tx.add('210')
    return tx
def dust(s):
    m=0.
    for st in range(0,max(1,len(s)-63),32):
        w=s[st:st+64]; l=len(w)-2
        if l<3: continue
        c=collections.Counter(w[i:i+3] for i in range(l)); sc=sum(v*(v-1)/2 for v in c.values())/(l-1)
        m=max(m,sc*100/31*(31/((l-1)/2)) if l<62 else sc*100/31)
    return m
def parse(bam):
    R={}
    out=subprocess.run(f"{ST} view -F 0x904 {bam}",shell=True,capture_output=True,text=True).stdout
    for l in out.splitlines():
        x=l.split('\t'); q=x[0]; mq=int(x[4]); tags=dict((t[:2],t[5:]) for t in x[11:])
        u=('XS' not in tags) or (int(tags['AS'])>int(tags['XS']))
        r=R.setdefault(q,[0,False,[]]); r[0]=max(r[0],mq); r[1]=r[1] or u; r[2].append(x[9])
    return R
rows=[]; cov=[]
for run in D.run:
    k2=P+run+'/k2.out.txt.gz'; tx=subtree(P+run+'/k2.report')
    conf=set(); H=parse(P+run+'/hp.bam')
    if not os.path.exists(k2):
        rows.append(dict(run=run,k2out_available=False,bt2_frag=len(H))); continue
    with gzip.open(k2,'rt') as f:
        for l in f:
            x=l.split('\t',3)
            if x[0]=='C' and x[2] in tx: conf.add(x[1])
    G=parse(A+'per_run/'+run+'/hp26695.bam')
    lc={q for q in conf if q in H and any(dust(s)>7 for s in H[q][2])}
    hh={q for q in conf if human is not None and q in human}
    c=dict(run=run,k2out_available=True,k2_conf_recount=len(conf),bt2_frag=len(H),
      lit_mapq30=sum(1 for q in conf if q in H and H[q][0]>=30),lit_unique=sum(1 for q in conf if q in H and H[q][1]),
      dev_mapq30=sum(1 for q in conf if q in G and G[q][0]>=30),dev_unique=sum(1 for q in conf if q in G and G[q][1]),
      dev_mapped=sum(1 for q in conf if q in G),lowcomp=len(lc),humanhit=len(hh) if human is not None else np.nan,
      clean=len(conf-lc-hh),strict=sum(1 for q in conf-lc-hh if q in G and G[q][0]>=30))
    rows.append(c)
    gq=[q for q in G if G[q][0]>=30]
    if len(gq)>=100:
        cv=subprocess.run(f"{ST} view -b -q 30 {A}per_run/{run}/hp26695.bam | {ST} coverage - | tail -1",shell=True,capture_output=True,text=True).stdout.split('\t')
        dep=subprocess.run(f"{ST} view -b -q 30 {A}per_run/{run}/hp26695.bam | {ST} depth -a - | cut -f3",shell=True,capture_output=True,text=True).stdout.split()
        dep=np.array(dep,dtype=float); bins=dep[:len(dep)//10000*10000].reshape(-1,10000).mean(1)
        bases=float(cv[6])*float(cv[2]) if len(cv)>6 else np.nan  # meandepth*length
        exp_b=1-np.exp(-bases/float(cv[2])); cov.append(dict(run=run,n_frag_q30=len(gq),breadth=float(cv[5])/100,expected_breadth_poisson=exp_b,
          breadth_ratio=(float(cv[5])/100)/exp_b,bin10k_cv=float(bins.std()/bins.mean()) if bins.mean()>0 else np.nan,frac_bins_covered=float((bins>0).mean())))
C=pd.DataFrame(rows); D=D.merge(C,on='run'); D.to_csv(A+'audit_per_run.csv',index=False)
CV=pd.DataFrame(cov).merge(D[['run','type','hp_log1p_rpm']],on='run') if cov else pd.DataFrame(); CV.to_csv(A+'audit_breadth_26695.csv',index=False)
def eff(Dx,col,label):
    Dx=Dx.copy(); Dx['y']=np.log1p(Dx[col]/Dx.salmon_mapped*1e6)
    w=Dx.pivot_table(index='pid',columns='type',values='y').dropna(); dd=w.N-w['T']; d2=dd[dd!=0]
    r=dict(subset=label,count_col=col,total_frag=int(Dx[col].sum()),runs_ge10=int((Dx[col]>=10).sum()),pairs=len(w),N_gt_T=int((dd>0).sum()),T_gt_N=int((dd<0).sum()),ties=int((dd==0).sum()))
    if len(d2)>0:
        rk=d2.abs().rank(); rp=rk[d2>0].sum(); tot=len(d2)*(len(d2)+1)/2; r['rank_biserial']=float((rp-(tot-rp))/tot)
        r['P_one_sided_N_gt_T']=float(wilcoxon(w.N,w['T'],alternative='greater').pvalue); r['P_two_sided']=float(wilcoxon(w.N,w['T']).pvalue)
        r['median_fold_rpm1']=float(np.exp(np.median(dd)))
        try:
            M=Dx[Dx.pid.isin(w.index)].copy(); M['logdepth']=np.log10(M.salmon_mapped); M['tum']=(M.type=='T').astype(int)
            m=smf.mixedlm('y~tum+logdepth',M,groups=M.pid).fit(); r['mixed_beta_tumour']=float(m.params['tum']); r['mixed_CI']=[float(v) for v in m.conf_int().loc['tum']]
        except Exception as e: r['mixed_error']=str(e)[:80]
    r['direction_reversed']=bool(r['T_gt_N']>r['N_gt_T']); return r
Dk=D[D.k2out_available==True]
S=[eff(D,'hp_k2_confirmed','main (A1 metric, reproduction)'),eff(D.assign(hp_k2_restr=D.hp_k2_confirmed.where(D.k2out_available==True)),'hp_k2_restr','main restricted to runs with k2.out (comparator for read-level subsets)'),eff(D,'k2_conf_recount','main recount from k2.out'),eff(D,'bt2_frag','bowtie2-only, no Kraken2 filter'),
   eff(D,'lit_mapq30','1 literal v0.6: MAPQ>=30 on multi-genome hp.bam'),eff(D,'lit_unique','2 literal: unique best hit (AS>XS) on multi-genome hp.bam'),
   eff(D,'dev_mapq30','1 DEVIATION FROM PROTOCOL (implementation): MAPQ>=30 after re-alignment to 26695'),eff(D,'dev_unique','2 DEVIATION (implementation): unique best hit on 26695'),
   eff(D,'clean','4 low-complexity (DUST>7) and human-transcriptome-homologous fragments removed'),eff(D,'strict','4+1dev combined: clean AND 26695 MAPQ>=30')]
# 5 Kraken2 concordance
conc=dict(runs_with_k2out=len(Dk),runs_without_k2out=sorted(D.run[D.k2out_available!=True].tolist()),total_bt2=int(Dk.bt2_frag.sum()),total_k2=int(Dk.k2_conf_recount.sum()),overall_concordance=float(Dk.k2_conf_recount.sum()/Dk.bt2_frag.sum()),
  per_run_median_by_type={t:float((g.k2_conf_recount/g.bt2_frag.replace(0,np.nan)).median()) for t,g in Dk.groupby('type')},
  recount_matches_A1_metric_runs=int((Dk.k2_conf_recount==Dk.hp_k2_confirmed).sum()),
  spearman_bt2rpm_vs_k2rpm=float(spearmanr(Dk.bt2_frag/Dk.salmon_mapped,Dk.k2_conf_recount/Dk.salmon_mapped).statistic))
# 6 leave-one-out extremes
w=D.pivot_table(index='pid',columns='type',values='hp_log1p_rpm').dropna(); top=(w.N-w['T']).abs().sort_values(ascending=False).index
loo=[]
for p in top[:10]: e=eff(D[D.pid!=p],'hp_k2_confirmed',f'drop {p}'); loo.append(dict(pid=str(p),rb=e['rank_biserial'],P1=e['P_one_sided_N_gt_T']))
mx=D.groupby('pid').hp_k2_confirmed.max().sort_values(ascending=False).index
S+= [eff(D[~D.pid.isin(top[:5])],'hp_k2_confirmed','6 drop top-5 |delta| patients jointly'),eff(D[~D.pid.isin(top[:10])],'hp_k2_confirmed','6 drop top-10 |delta| patients jointly'),
     eff(D[~D.pid.isin(mx[:10])],'hp_k2_confirmed','6 drop 10 patients with highest HP load')]
# 7 batch
ena=pd.read_csv(B+'data/external/PRJNA1174735_ena.tsv',sep='\t'); inst=ena[ena.run_accession.isin(D.run)].instrument_model.value_counts().to_dict()
out=dict(label='v0.6 §3 measurement audit (SENSITIVITY, not confirmatory; no multiplicity correction)',criterion='direction N>T not reversed; P<0.05 not required',
  subsets=S,kraken2_concordance=conc,loo_single_top10=dict(rows=loo,rb_min=min(x['rb'] for x in loo),rb_max=max(x['rb'] for x in loo),P1_max=max(x['P1'] for x in loo)),
  batch=dict(instrument_model=inst,lane='not recoverable: original FASTQ deleted and SRA-normalised read names carry no flowcell/lane; no adjustment possible'),
  human_homology_check_available=human is not None)
if len(CV):
    out['breadth_26695']=dict(n_runs=len(CV),by_type=CV.groupby('type')[['breadth','breadth_ratio','bin10k_cv','frac_bins_covered']].median().to_dict(),
      spearman_breadthratio_vs_load=float(spearmanr(CV.breadth_ratio,CV.hp_log1p_rpm).statistic),
      MW_breadthratio_N_vs_T_P=float(mannwhitneyu(CV[CV.type=='N'].breadth_ratio,CV[CV.type=='T'].breadth_ratio).pvalue) if CV.type.nunique()==2 else None)
out['any_direction_reversal']=any(s['direction_reversed'] for s in S)
json.dump(out,open(A+'audit_results.json','w'),indent=1,default=float,ensure_ascii=False)
print(pd.DataFrame(S)[['subset','total_frag','runs_ge10','N_gt_T','T_gt_N','ties','rank_biserial','P_one_sided_N_gt_T','median_fold_rpm1','direction_reversed']].to_string())
print(json.dumps({k:out[k] for k in out if k not in ('subsets',)},indent=1,default=float,ensure_ascii=False))
