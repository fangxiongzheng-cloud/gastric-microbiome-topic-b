# v0.6 §2 row 7 "with F2" version of the pre-specified gradient meta (gap G-1, user-approved t889u 2026-11-11).
# Pooling code copied verbatim from scripts/meta4.py (REML tau2 + Hartung-Knapp, Fisher z of rank-biserial r, se=1/sqrt(nonzero-3)).
# Only change: F2 PRJNA1061213 (SILVA primary, 274 pairs) appended as one dataset. Existing meta_gradient_4cohorts* files untouched.
import pandas as pd,numpy as np,json,sys
from scipy.stats import t as tdist
from scipy.optimize import minimize_scalar
M=pd.read_csv('results/korea16s/meta_gradient.csv')
y=json.load(open('results/yu310127/primary_tests.json'))['silva_main']
M=pd.concat([M,pd.DataFrame([dict(study='China PRJNA310127 16S (SILVA)',pairs=y['pairs'],nonzero=y['nonzero'],r=y['rank_biserial_r'],z=np.arctanh(y['rank_biserial_r']),se=1/np.sqrt(y['nonzero']-3),P=y['wilcoxon_clr_P'])])],ignore_index=True)
f=json.load(open('results/f2/F2_results.json'))['silva_main_ge5000']
nz=f['N_gt_T']+f['T_gt_N']
M=pd.concat([M,pd.DataFrame([dict(study='F2 PRJNA1061213 16S (SILVA)',pairs=f['pairs_used'],nonzero=nz,r=f['rank_biserial_r'],z=np.arctanh(f['rank_biserial_r']),se=1/np.sqrt(nz-3),P=f['wilcoxon_two_sided_P'])])],ignore_index=True)
if len(sys.argv)>1: M=M[M.study!='TCGA Spectrum (RNA)']
z,v=M.z.values,M.se.values**2;k=len(z)
def nll(t2):
    w=1/(v+t2);mu=(w*z).sum()/w.sum();return 0.5*(np.log(v+t2).sum()+np.log(w.sum())+(w*(z-mu)**2).sum())
t2=minimize_scalar(nll,bounds=(0,5),method='bounded').x; t2=0 if t2<1e-6 else t2
w=1/(v+t2);mu=(w*z).sum()/w.sum();q=(w*(z-mu)**2).sum()/(k-1);se=np.sqrt(max(q,1e-12)/w.sum());tc=tdist.ppf(.975,k-1)
wf=1/v;Q=(wf*(z-(wf*z).sum()/wf.sum())**2).sum()
res=dict(k=k,total_pairs=int(M.pairs.sum()),method='REML + Hartung-Knapp (Fisher z of rank-biserial r)',r=np.tanh(mu),lo=np.tanh(mu-tc*se),hi=np.tanh(mu+tc*se),P=2*tdist.sf(abs(mu/se),k-1),tau2=t2,Q=Q,I2=max(0,(Q-(k-1))/Q),
  note='v0.6 sec2 row7 with-F2 version; same code as scripts/meta4.py; F2 row P column = two-sided Wilcoxon (descriptive, not used in pooling)')
sfx='_noSpectrum' if len(sys.argv)>1 else '';M.to_csv(f'results/meta_gradient_withF2{sfx}.csv',index=False);json.dump(res,open(f'results/meta_gradient_withF2{sfx}.json','w'),indent=1,default=float);print(M[['study','pairs','nonzero','r']]);print(res)
