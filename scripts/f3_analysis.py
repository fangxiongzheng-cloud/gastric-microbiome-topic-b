# F3 (Family F3, confirmatory secondary; gated on F2 pass) — PRJNA1061213 differentiation interaction.
# Spec: manuscript/protocol_v0.5_addendum.md §B + plans/analysis_plan_v0.6.md §2 row 4. Decisions logged in logs/freeze_log.md
# (2026-11-11 09:52 CST) before running. Model: OLS Delta ~ poor + stageIII_IV, one-sided on beta_poor (H1 >0), alpha 0.05.
import pandas as pd,numpy as np,json,hashlib,statsmodels.formula.api as smf,statsmodels,sys
from scipy import stats
B='/REPLACE_WITH_PRIVATE_WORKDIR/';O=B+'results/f3/';SEED=20261030
X=B+'data/external/third_cohort/PRJNA1061213/supp/41421_2024_746_MOESM2_ESM.xlsx'
sha=lambda p:hashlib.sha256(open(p,'rb').read()).hexdigest()
F={'v0.5':B+'manuscript/protocol_v0.5_addendum.md','v0.6':B+'plans/analysis_plan_v0.6.md'}
H={k:sha(v) for k,v in F.items()}
assert H['v0.5'].startswith('14dc4ec9') and H['v0.6'].startswith('3243936d'),H
f2=json.load(open(B+'results/f2/F2_results.json')); gate=bool(f2['F2_PRIMARY']['pass_'])
# --- covariates
def tab(s):
    t=pd.read_excel(X,s,header=1); t['pid']=t['No.'].str.extract(r'(\d+)')[0].astype(int); return t.set_index('pid')
t1,t2=tab('Table1'),tab('Table2')
norm=lambda s:s.astype(str).str.strip().str.lower()
gmap={'low':1,'undifferentiated':1,'medium':0,'high':0}
smap={'i':0,'ii':0,'iii':1,'iv':1}
cov=pd.DataFrame({'grade_raw':norm(t1['Grade of differentiation']),'grade_raw_S2':norm(t2['Grade of differentiation']).reindex(t1.index),
                  'tnm_S1':norm(t1['TNM type']),'tnm_S2':norm(t2['TNM type']).reindex(t1.index)})
cov['poor']=cov.grade_raw.map(gmap)  # unknown -> NaN
cov['stage_S1']=cov.tnm_S1.map(smap); cov['stage_S2']=cov.tnm_S2.map(smap)
cov['stage']=np.where(cov.stage_S1==cov.stage_S2,cov.stage_S1,np.nan)   # conflict -> missing (primary)
cov.loc[cov.tnm_S2.isna()|cov.stage_S2.isna(),'stage']=cov.stage_S1   # no S2 counterpart -> S1 (not needed for pairs)
# --- outcome from F2 per-sample tables
def delta(cls):
    d=pd.read_csv(B+f'results/f2/per_sample_{cls}.csv',index_col=0); d=d[d.mindepth>=5000]
    d['lg']=np.log10(d.ra+1e-5)
    w=d.pivot_table(index='pid',columns='type',values=['lg','clr']).dropna()
    return pd.DataFrame({'d_log10':w['lg']['N']-w['lg']['T'],'d_clr':w['clr']['N']-w['clr']['T']})
def fit(df,y,stagecol='stage',cov_type='nonrobust'):
    dd=df[[y,'poor',stagecol]].dropna().rename(columns={stagecol:'st'})
    m=smf.ols(f'{y} ~ poor + st',data=dd).fit(cov_type=cov_type)
    b,se=m.params['poor'],m.bse['poor']; df_=m.df_resid
    if cov_type=='nonrobust': p1=stats.t.sf(b/se,df_); q=stats.t.ppf(0.975,df_)
    else: p1=stats.norm.sf(b/se); q=stats.norm.ppf(0.975)
    return dict(n=int(len(dd)),n_poor=int(dd.poor.sum()),n_wellmod=int((dd.poor==0).sum()),beta_poor=float(b),se=float(se),
                ci95=[float(b-q*se),float(b+q*se)],p_one_sided=float(p1),p_two_sided=float(2*min(p1,1-p1)),
                beta_stageIII_IV=float(m.params['st']),p_stage_two_sided=float(m.pvalues['st']),cov_type=cov_type,r2=float(m.rsquared))
R=dict(spec='F3: OLS Delta(log10 HP N - T) ~ poor(vs well/moderate) + TNM(III-IV vs I-II); one-sided H1 beta_poor>0, alpha=0.05',
       frozen_sha256=H,gate_F2_pass=gate,status='confirmatory secondary (Family F3, 1 test, no correction)' if gate else 'exploratory (F2 not passed)',
       seed=SEED,software=dict(statsmodels=statsmodels.__version__,python=sys.version.split()[0]))
D=delta('silva').join(cov,how='left'); D.to_csv(O+'F3_analysis_table_silva.csv')
npairs=len(D); miss_grade=int(D.poor.isna().sum()); conflict=D.index[(D.stage_S1!=D.stage_S2)&D.stage_S2.notna()].tolist()
R['flow']=dict(pairs_F2=npairs,covariate_rows_matched=int(D.grade_raw.notna().sum()),excluded_grade_unknown=miss_grade,
  excluded_stage_conflict=len([p for p in conflict if not np.isnan(D.loc[p,'poor'])]),stage_conflict_pids=[int(p) for p in conflict],
  grade_S1_vs_S2_mismatch_after_case_norm=int((D.grade_raw!=D.grade_raw_S2).sum()),
  grade_counts=D.grade_raw.value_counts().to_dict(),undifferentiated_pids=[int(p) for p in D.index[D.grade_raw=='undifferentiated']])
P=fit(D,'d_log10'); P['pass']=bool(P['beta_poor']>0 and P['p_one_sided']<0.05); R['PRIMARY']=P
S={}
S['HC3_SE (v0.6 pre-specified)']=fit(D,'d_log10',cov_type='HC3')
S['outcome_dCLR_silva (F2 metric)']=fit(D,'d_clr')
S['outcome_dCLR_silva_HC3']=fit(D,'d_clr',cov_type='HC3')
G=delta('gg2v4').join(cov,how='left')
S['GG2_classifier_log10']=fit(G,'d_log10'); S['GG2_classifier_dCLR']=fit(G,'d_clr')
S['stage_conflict_coded_S1']=fit(D,'d_log10','stage_S1'); S['stage_conflict_coded_S2']=fit(D,'d_log10','stage_S2')
S['exclude_undifferentiated']=fit(D[D.grade_raw!='undifferentiated'],'d_log10')
R['SENSITIVITY']=S
# TNM-stratified (descriptive)
st={}
for lab,v in [('TNM_I_II',0),('TNM_III_IV',1)]:
    dd=D[(D.stage==v)][['d_log10','poor']].dropna(); m=smf.ols('d_log10 ~ poor',data=dd).fit(); q=stats.t.ppf(.975,m.df_resid)
    st[lab]=dict(n=len(dd),n_poor=int(dd.poor.sum()),beta_poor=float(m.params['poor']),ci95=[float(m.params['poor']-q*m.bse['poor']),float(m.params['poor']+q*m.bse['poor'])])
R['TNM_stratified_descriptive']=st
dd=D[['d_log10','poor','stage']].dropna()
R['group_descriptive']={('poor' if g==1 else 'well_moderate'):dict(n=int(len(x)),median_delta_log10=float(x.d_log10.median()),mean_delta_log10=float(x.d_log10.mean()),
    frac_N_gt_T=float((x.d_log10>0).mean())) for g,x in dd.groupby('poor')}
# equivalence check: long format, patient fixed effects, tissue x poor interaction
L=[]
d0=pd.read_csv(B+'results/f2/per_sample_silva.csv',index_col=0); d0=d0[d0.mindepth>=5000]; d0['lg']=np.log10(d0.ra+1e-5)
d0=d0[d0.pid.isin(dd.index)].copy(); d0['normal']=(d0.type=='N').astype(int); d0['poor']=d0.pid.map(dd.poor); d0['st']=d0.pid.map(dd.stage)
ml=smf.ols('lg ~ C(pid) + normal + normal:poor + normal:st',data=d0).fit()
R['long_format_check']=dict(beta_normal_x_poor=float(ml.params['normal:poor']),matches_primary=bool(abs(ml.params['normal:poor']-P['beta_poor'])<1e-8))
json.dump(R,open(O+'F3_results.json','w'),indent=1,default=float); print(json.dumps(R,indent=1,default=float))
