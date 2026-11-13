import pandas as pd,numpy as np,gzip,json,glob,os
from scipy.stats import spearmanr
B='/REPLACE_WITH_PRIVATE_WORKDIR/'; O=B+'results/a4_atrophy/'
IM=['CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4']; OX=['ATP4A','ATP4B','CBLIF','PGC','GKN1','GKN2']
t2g={}
with gzip.open('/REPLACE_WITH_REFERENCE_ROOT/gencode/gencode.v46.transcripts.fa.gz','rt') as f:
    for l in f:
        if l[0]=='>': x=l[1:].split('|'); t2g[x[0]]=x[5]
L=pd.read_csv(B+'data/external/rna_pilot40.tsv',sep='\t')
ex={}
for _,r in L.iterrows():
    d=B+'results/korea_rna/per_run/'+r.run_accession+'/'
    if not os.path.exists(d+'DONE'): continue
    q=pd.read_csv(d+'salmon/quant.sf.gz',sep='\t',usecols=['Name','TPM'])
    q['g']=q.Name.map(t2g); ex[r.library_name]=q.groupby('g').TPM.sum()
E=np.log2(pd.DataFrame(ex)+1); E.loc[IM+OX].T.to_csv(O+'korea_genes_log2tpm.csv')
s=pd.read_csv(B+'results/korea16s/per_sample_silva.csv',index_col=0); s['lib']=s.pid.astype(str)+'_'+s.type_bs
s=s.set_index('lib')
def score(cols,gs):
    X=E.loc[gs,cols]; z=X.sub(X.mean(1),axis=0).div(X.std(1,ddof=1),axis=0); return z.mean()
N=[c for c in E.columns if c.endswith('_N') and c in s.index]
D=pd.DataFrame({'clr':s.loc[N,'clr'],'reads16s':s.loc[N,'reads'],'IM':score(N,IM),'OX':score(N,OX)})
# ssGSEA-like sensitivity: rank-based (singscore-style mean rank of gene set within sample, all expressed genes)
R=E[N][(E[N]>1).mean(1)>=0.5].rank(); D['IM_rank']=R.loc[[g for g in IM if g in R.index]].mean()/len(R)
D['pid']=[c[:-2] for c in N]
T={c[:-2]:s.loc[c,'clr'] for c in E.columns if c.endswith('_T') and c in s.index}
D['dCLR']=D.pid.map(lambda p: D.loc[p+'_N','clr']-T[p] if p in T else np.nan)
D.to_csv(O+'korea_normal_scores.csv')
def sp(a,b): d=D[[a,b]].dropna(); r=spearmanr(d[a],d[b]); return dict(x=a,y=b,n=len(d),rho=float(r.statistic),P=float(r.pvalue))
res=dict(primary=sp('clr','IM'),secondary=[sp('clr','OX'),sp('dCLR','IM'),sp('clr','IM_rank')],
  genes_missing=[g for g in IM+OX if g not in E.index])
from statsmodels.stats.multitest import multipletests
q=multipletests([x['P'] for x in res['secondary']],method='fdr_bh')[1]
for x,qq in zip(res['secondary'],q): x['q_BH']=float(qq)
res['sens_min1000reads']=sp.__call__ if False else None
d=D[D.reads16s>=1000]; r=spearmanr(d.clr,d.IM); res['sens_16s_ge1000']=dict(n=len(d),rho=float(r.statistic),P=float(r.pvalue))
print(json.dumps(res,indent=1)); json.dump(res,open(O+'korea_tests.json','w'),indent=1)
