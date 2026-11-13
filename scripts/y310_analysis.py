# PRJNA310127 China subset: protocol v0.4 A2 primary test (same rules as Korea P1; depth = min over classifiers)
import pandas as pd, numpy as np, json
from scipy.stats import wilcoxon, binomtest
D='data/external/yu310127/'; O='results/yu310127/'
import os; os.makedirs(O,exist_ok=True)
t=pd.read_csv(D+'exp_table/table.tsv',sep='\t',skiprows=1,index_col=0)
runs=pd.read_csv(D+'runs.tsv',sep='\t',header=None).iloc[:,:2]; runs.columns=['run','alias']
runs['pid']=runs.alias.str[:-1]; runs['type_bs']=runs.alias.str[-1]; runs=runs.set_index('run')
G={}
for n in ['silva','gg2v4']:
    tx=pd.read_csv(D+f'exp_tax_{n}/taxonomy.tsv',sep='\t',index_col=0)['Taxon'].reindex(t.index).fillna('Unassigned')
    bad=tx.str.contains('Mitochondria|Chloroplast|mitochondria|chloroplast')|~tx.str.contains('p__')|tx.str.startswith('Unassigned')|tx.str.contains('Eukaryota')
    lv=tx[~bad].str.split(';').map(lambda x:';'.join([y.strip() for y in x[:6]]))
    g=t[~bad].groupby(lv).sum().T; g.to_csv(O+f'genus_{n}.tsv',sep='\t'); G[n]=g
    print(n,'ASVs',len(t),'removed',int(bad.sum()),'Helicobacter ASVs',int(tx.str.contains('Helicobacter').sum()))
mindepth=pd.concat({k:v.sum(1) for k,v in G.items()},axis=1).min(1)
res={}
for cl,g in G.items():
    tot=g.sum(1); hcols=[c for c in g.columns if c.split(';')[-1].strip().replace('g__','').startswith('Helicobacter')]
    h=g[hcols].sum(1); g2=g.drop(columns=hcols).assign(HELICO=h); lg=np.log(g2+0.5); hc=lg['HELICO']-lg.mean(1)
    df=pd.DataFrame({'reads':tot,'hp':h,'ra':h/tot,'clr':hc}).join(runs[['pid','type_bs']]); df.to_csv(O+f'per_sample_{cl}.csv')
    for label,minr in [('main',1000),('sens5k',5000)]:
        ok=df[mindepth.reindex(df.index)>=minr]
        w=ok.pivot_table(index='pid',columns='type_bs',values=['clr','ra','hp']).dropna(); n=len(w)
        d=w['clr']['N']-w['clr']['T']; s1=wilcoxon(w['clr']['N'],w['clr']['T']); s2=wilcoxon(w['ra']['N'],w['ra']['T'])
        pos=lambda tt:(w['hp'][tt]>=2)&(w['ra'][tt]>=0.001)
        b=int((pos('N')&~pos('T')).sum()); c=int((~pos('N')&pos('T')).sum()); mc=binomtest(b,b+c,0.5).pvalue if b+c>0 else np.nan
        nz=d[d!=0]; rk=nz.abs().rank(); r=(rk[nz>0].sum()-rk[nz<0].sum())/rk.sum()
        res[f'{cl}_{label}']=dict(pairs=n,N_gt_T=int((d>0).sum()),T_gt_N=int((d<0).sum()),ties=int((d==0).sum()),nonzero=int(len(nz)),
          median_ra_N=w['ra']['N'].median(),median_ra_T=w['ra']['T'].median(),median_diff_clr=d.median(),
          wilcoxon_clr_V=s1.statistic,wilcoxon_clr_P=s1.pvalue,wilcoxon_ra_P=s2.pvalue,rank_biserial_r=r,
          prev_N=pos('N').mean(),prev_T=pos('T').mean(),mcnemar_b=b,mcnemar_c=c,mcnemar_P=mc)
json.dump(res,open(O+'primary_tests.json','w'),indent=1,default=float); print(json.dumps(res,indent=1,default=float))
