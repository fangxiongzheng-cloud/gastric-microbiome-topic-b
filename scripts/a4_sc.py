import scanpy as sc,pandas as pd,numpy as np,json,scipy.sparse as sp
from scipy.stats import kruskal,mannwhitneyu
O='/REPLACE_WITH_PRIVATE_WORKDIR/results/a4_atrophy/'
a=sc.read_h5ad('/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna/work/epithelial.h5ad')
print(a, a.obs.columns.tolist())
X=a.copy(); X.var_names=X.var['sym'].astype(str).values; X.var_names_make_unique()
IM=['CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4']; OX=['ATP4A','ATP4B','CBLIF','PGC','GKN1','GKN2']
dcol=[c for c in ['donor','sample','library'] if c in a.obs][0]
genes=[g for g in IM+OX if g in X.var_names]; print('missing',set(IM+OX)-set(genes))
M=X[:,genes].X; M=M.toarray() if sp.issparse(M) else M
# counts? use per-donor mean of log-normalised expression if data looks normalised
df=pd.DataFrame(M,columns=genes); df['d']=a.obs[dcol].values
P=df.groupby('d',observed=True).mean(); z=(P-P.mean())/P.std()
meta=a.obs[[dcol,'stage','hp']].drop_duplicates(dcol).set_index(dcol)
R=meta.join(pd.DataFrame({'IM':z[[g for g in IM if g in z]].mean(1),'OX':z[[g for g in OX if g in z]].mean(1),'ncells':df.groupby('d',observed=True).size()}))
R.to_csv(O+'sc_donor_scores.csv'); print(R.sort_values(['stage','hp']))
res={'note':'descriptive, n=3 donors per stage x HP; epithelial cells, donor pseudobulk mean of normalised expression'}
for s in ['IM','OX']:
    res[s+'_kruskal_stage']=float(kruskal(*[g[s] for _,g in R.groupby('stage',observed=True)]).pvalue)
    res[s+'_HPpos_vs_neg_MWU_allstages']=float(mannwhitneyu(R[R.hp=='HP+'][s],R[R.hp=='HP-'][s]).pvalue)
    res[s+'_means']=R.groupby(['stage','hp'],observed=True)[s].mean().round(2).to_dict().__repr__()
print(res); json.dump(res,open(O+'sc_exploratory.json','w'),indent=1)
