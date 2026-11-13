# Core step 1+2: assemble matrices and A1 confirmatory (v0.4 A1 star test; v0.6 wording)
import pandas as pd,numpy as np,json,subprocess,os,gzip
from scipy.stats import spearmanr,wilcoxon
B='/REPLACE_WITH_PRIVATE_WORKDIR/'; O=B+'results/korea_rna/'; C=B+'results/core/'; ST='/REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin/samtools'
L=pd.read_csv(B+'data/external/rna_full352.tsv',sep='\t'); rows=[]
def k2clade(f,taxid='210'):
    for l in open(f):
        x=l.rstrip('\n').split('\t')
        if len(x)>4 and x[4]==taxid: return int(x[1])
    return 0
D=pd.read_csv(C+'korea_rna_per_sample.csv')  # assembled by first part of this script (extraction block removed after first run; see git history of this file: extraction used k2 taxid210 clade reads, salmon num_mapped, bt2 fragments, gene TPM=sum transcript TPM)
w=D.pivot_table(index='pid',columns='type',values='hp_log1p_rpm')
res=dict(label='A1 CONFIRMATORY (pre-registered, NOT blind; not independent confirmation)',metric='测序深度归一化的 HP 转录本读段负荷 (Kraken2-confirmed taxid210 reads / Salmon-mapped host fragments x1e6)',
 n_runs=len(D),n_patients=int(D.pid.nunique()),types=D.type.value_counts().to_dict(),libtype=D.libtype.value_counts().to_dict(),
 n_hp_ge10=int((D.hp_k2_confirmed>=10).sum()),n_hp_ge1000=int((D.hp_k2_confirmed>=1000).sum()),detect_rate_ge10=float((D.hp_k2_confirmed>=10).mean()),max_hp=int(D.hp_k2_confirmed.max()))
res['patients_with_only_one_tissue']=int(w.isna().any(axis=1).sum()); w=w.dropna()
res.update(pairs=len(w),N_gt_T=int((w.N>w['T']).sum()),T_gt_N=int((w['T']>w.N).sum()),ties=int((w.N==w['T']).sum()))
dd=w.N-w['T']; pw=wilcoxon(w.N,w['T'],alternative='greater'); pt=wilcoxon(w.N,w['T'])
res['wilcoxon_P_one_sided_N_gt_T']=float(pw.pvalue); res['wilcoxon_P_two_sided']=float(pt.pvalue)
d2=dd[dd!=0]; n=len(d2); rk=d2.abs().rank(); rp=rk[d2>0].sum(); tot=n*(n+1)/2; res['rank_biserial_r_N_vs_T']=float((rp-(tot-rp))/tot)
rng=np.random.default_rng(20261030); bs=[]
v=d2.values
for _ in range(2000):
    x=rng.choice(v,len(v)); rr=pd.Series(np.abs(x)).rank().values; p=rr[x>0].sum(); bs.append((p-(rr.sum()-p))/rr.sum())
res['rank_biserial_95CI_bootstrap']=[float(np.percentile(bs,2.5)),float(np.percentile(bs,97.5))]
W=D.pivot_table(index='pid',columns='type',values='hp_rpm_host').loc[w.index]
res['median_rpm_N']=float(W.N.median()); res['median_rpm_T']=float(W['T'].median()); res['median_fold_N_over_T_rpm_plus1']=float(np.exp(np.median(np.log(W.N+1)-np.log(W['T']+1))))
res['A1_success_P<0.05_and_N>T']=bool(pw.pvalue<0.05 and res['N_gt_T']>res['T_gt_N'])
# secondary: 16S CLR vs log1p RPM
d=D[['s16_clr','hp_log1p_rpm']].dropna(); r=spearmanr(d.s16_clr,d.hp_log1p_rpm); res['spearman_16S_CLR_vs_log1pRPM']=dict(n=len(d),rho=float(r.statistic),P=float(r.pvalue))
# McNemar detection >=10
from statsmodels.stats.contingency_tables import mcnemar
P=D.pivot_table(index='pid',columns='type',values='hp_k2_confirmed').loc[w.index]>=10
tb=[[int((P.N&P['T']).sum()),int((P.N&~P['T']).sum())],[int((~P.N&P['T']).sum()),int((~P.N&~P['T']).sum())]]
res['mcnemar_detect_ge10']=dict(table_NposTpos_NposTneg_NnegTpos_NnegTneg=tb,P=float(mcnemar(tb,exact=True).pvalue))
# depth-adjusted mixed model
import statsmodels.formula.api as smf
M=D[D.pid.isin(w.index)].copy(); M['logdepth']=np.log10(M.salmon_mapped); M['tum']=(M.type=='T').astype(int)
m=smf.mixedlm('hp_log1p_rpm~tum+logdepth',M,groups=M.pid).fit()
res['mixed_depth_adjusted']=dict(beta_tumour=float(m.params['tum']),CI=[float(x) for x in m.conf_int().loc['tum']],P=float(m.pvalues['tum']))
json.dump(res,open(C+'A1_confirmatory.json','w'),indent=1); print(json.dumps(res,indent=1))
