# Step10 E5a (OLS+HC3) + data prep for E5b (Cox in R). v0.6a E5 / v0.6c A. Exploratory. Decisions: freeze_log 2026-11-11 10:20.
import pandas as pd,numpy as np,json,statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
B='/REPLACE_WITH_PRIVATE_WORKDIR/';O=B+'results/step10/E5/'
X=B+'data/external/third_cohort/PRJNA1061213/supp/41421_2024_746_MOESM2_ESM.xlsx'
t=pd.read_excel(X,'Table1',header=1); t['pid']=t['No.'].str.extract(r'(\d+)')[0].astype(int); t=t.set_index('pid')
n=lambda s:s.astype(str).str.strip().str.lower()
CV=pd.DataFrame(index=t.index)
CV['site']=n(t['Tumor location']).map({'gastric antrum':'antrum','gastric body':'body','gastric cardia':'cardia'})
CV['poor']=n(t['Grade of differentiation']).map({'low':1,'undifferentiated':1,'medium':0,'high':0})
CV['stage34']=n(t['TNM type']).map({'i':0,'ii':0,'iii':1,'iv':1})
CV['dMMR']=n(t['MMR']).map({'dmmr':1,'pmmr':0})
CV['age']=pd.to_numeric(t.Age,errors='coerce'); CV['male']=(n(t.Sex)=='male').astype(int)
CV['os_time']=pd.to_numeric(t['Survival period'],errors='coerce'); CV['os_event']=pd.to_numeric(t['Survival state'],errors='coerce')
d=pd.read_csv(B+'results/f2/per_sample_silva.csv',index_col=0); d=d[d.mindepth>=5000]
w=d.pivot_table(index='pid',columns='type',values='clr').dropna()
P=pd.DataFrame({'clr_N':w['N'],'clr_T':w['T']}); P['dHP_TminusN']=P.clr_T-P.clr_N
P.index=P.index.astype(int); M=P.join(CV,how='left')
res={'label':'E5 EXPLORATORY (v0.6a); E5b supplementary only','n_pairs_F2':len(P),'n_pairs_unlinked':int(M.age.isna().sum())}
M.to_csv(O+'E5_analysis_table.csv')
dd=M.dropna(subset=['dHP_TminusN','site','poor','stage34','dMMR','age','male']).copy()
res['E5a_excluded_counts']={c:int(M[c].isna().sum()) for c in ['site','poor','stage34','dMMR']}
m=smf.ols('dHP_TminusN ~ C(site,Treatment("antrum")) + poor + stage34 + dMMR + age + male',data=dd).fit(cov_type='HC3')
ci=m.conf_int()
terms=[k for k in m.params.index if k!='Intercept']
q=multipletests([m.pvalues[k] for k in terms],method='fdr_bh')[1]
res['E5a']=dict(n=len(dd),model='OLS HC3; dHP = CLR_T - CLR_N (SILVA, F2 274 pairs, complete case)',r2=float(m.rsquared),
  coef={k:dict(beta=float(m.params[k]),CI95=[float(ci.loc[k,0]),float(ci.loc[k,1])],P=float(m.pvalues[k]),q_BH=float(qq)) for k,qq in zip(terms,q)},
  n_by_site=dd.site.value_counts().to_dict(),n_poor=int(dd.poor.sum()),n_stage34=int(dd.stage34.sum()),n_dMMR=int(dd.dMMR.sum()))
json.dump(res,open(O+'E5a_results.json','w'),indent=1); print(json.dumps(res,indent=1))
