import pandas as pd,numpy as np,json
from scipy.stats import t as tdist
from scipy.optimize import minimize_scalar
M=pd.read_csv('results/korea16s/meta_gradient.csv')
y=json.load(open('results/yu310127/primary_tests.json'))['silva_main']
M=pd.concat([M,pd.DataFrame([dict(study='China PRJNA310127 16S (SILVA)',pairs=y['pairs'],nonzero=y['nonzero'],r=y['rank_biserial_r'],z=np.arctanh(y['rank_biserial_r']),se=1/np.sqrt(y['nonzero']-3),P=y['wilcoxon_clr_P'])])],ignore_index=True)
import sys
if len(sys.argv)>1: M=M[M.study!='TCGA Spectrum (RNA)']
z,v=M.z.values,M.se.values**2;k=len(z)
def nll(t2):
    w=1/(v+t2);mu=(w*z).sum()/w.sum();return 0.5*(np.log(v+t2).sum()+np.log(w.sum())+(w*(z-mu)**2).sum())
t2=minimize_scalar(nll,bounds=(0,5),method='bounded').x; t2=0 if t2<1e-6 else t2
w=1/(v+t2);mu=(w*z).sum()/w.sum();q=(w*(z-mu)**2).sum()/(k-1);se=np.sqrt(max(q,1e-12)/w.sum());tc=tdist.ppf(.975,k-1)
wf=1/v;Q=(wf*(z-(wf*z).sum()/wf.sum())**2).sum()
res=dict(k=k,total_pairs=int(M.pairs.sum()),method='REML + Hartung-Knapp (Fisher z of rank-biserial r)',r=np.tanh(mu),lo=np.tanh(mu-tc*se),hi=np.tanh(mu+tc*se),P=2*tdist.sf(abs(mu/se),k-1),tau2=t2,Q=Q,I2=max(0,(Q-(k-1))/Q))
sfx='_noSpectrum' if len(sys.argv)>1 else '';M.to_csv(f'results/meta_gradient_4cohorts{sfx}.csv',index=False);json.dump(res,open(f'results/meta_gradient_4cohorts{sfx}.json','w'),indent=1,default=float);print(M[['study','pairs','r','P']]);print(res)
