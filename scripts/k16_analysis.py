# Pre-specified primary test (protocol v0.3 P1)
import pandas as pd, numpy as np, json
from scipy.stats import wilcoxon, binomtest
D='data/external/korea16s/'; O='results/korea16s/'
runs=pd.read_csv(D+'runs16s.tsv',sep='\t',header=None,names=['run','biosample','lib','st','ins','rc','ftp','bytes'])
bs=pd.read_csv('data/external/PRJNA1174735_biosample.tsv',sep='\t')[['biosample','env_medium']]
runs=runs.merge(bs,on='biosample',how='left')
runs['type_bs']=np.where(runs.env_medium.str.startswith('Neoplastic'),'T','N')
runs['type_lib']=runs.lib.str[-1]; runs['pid']=runs.lib.str[:-2]
print('lib vs biosample label agreement', (runs.type_bs==runs.type_lib).mean())
def load(cl):
    g=pd.read_csv(O+f'genus_{cl}.tsv',sep='\t',index_col=0)  # samples x genus counts
    return g
res={}
# common QC set: sample must pass depth in BOTH classifiers (classifier-independent pair set)
_g={cl:pd.read_csv(O+f'genus_{cl}.tsv',sep='\t',index_col=0).sum(1) for cl in ['gg2v4','silva']}
mindepth=pd.concat(_g,axis=1).min(1)
for cl in ['gg2v4','silva']:
    try: g=load(cl)
    except FileNotFoundError: continue
    tot=g.sum(1); hcols=[c for c in g.columns if c.split(';')[-1].strip().replace('g__','').startswith('Helicobacter')]
    h=g[hcols].sum(1)
    g2=g.drop(columns=hcols).assign(HELICO=h)
    lg=np.log(g2+0.5); hc=lg['HELICO']-lg.mean(1)
    df=pd.DataFrame({'reads':tot,'hp':h,'ra':h/tot,'clr':hc}).join(runs.set_index('run')[['pid','type_bs']])
    df.to_csv(O+f'per_sample_{cl}.csv')
    for label,minr in [('main',1000),('sens5k',5000)]:
        ok=df[mindepth.reindex(df.index)>=minr]
        w=ok.pivot_table(index='pid',columns='type_bs',values=['clr','ra','hp'])
        w=w.dropna()
        n=len(w)
        d=w['clr']['N']-w['clr']['T']
        s1=wilcoxon(w['clr']['N'],w['clr']['T']); s2=wilcoxon(w['ra']['N'],w['ra']['T'])
        pos=lambda t:(w['hp'][t]>=2)&(w['ra'][t]>=0.001)
        b=int((pos('N')&~pos('T')).sum()); c=int((~pos('N')&pos('T')).sum())
        mc=binomtest(b,b+c,0.5).pvalue if b+c>0 else np.nan
        nz=d[d!=0]; rk=nz.abs().rank(); r=(rk[nz>0].sum()-rk[nz<0].sum())/rk.sum()
        res[f'{cl}_{label}']=dict(pairs=n,N_gt_T=int((d>0).sum()),T_gt_N=int((d<0).sum()),ties=int((d==0).sum()),
          median_clr_N=w['clr']['N'].median(),median_clr_T=w['clr']['T'].median(),median_diff_clr=d.median(),
          median_ra_N=w['ra']['N'].median(),median_ra_T=w['ra']['T'].median(),mean_ra_N=w['ra']['N'].mean(),mean_ra_T=w['ra']['T'].mean(),
          wilcoxon_clr_V=s1.statistic,wilcoxon_clr_P=s1.pvalue,wilcoxon_ra_P=s2.pvalue,rank_biserial_r=r,
          prev_N=pos('N').mean(),prev_T=pos('T').mean(),mcnemar_b_NposTneg=b,mcnemar_c_NnegTpos=c,mcnemar_P=mc)
    # rarefied sensitivity
    rng=np.random.default_rng(1); ok=df[mindepth.reindex(df.index)>=5000].index
    hr=pd.Series({s:rng.hypergeometric(int(df.hp[s]),int(df.reads[s]-df.hp[s]),5000) for s in ok})
    t=df.loc[ok].assign(hr=hr).pivot_table(index='pid',columns='type_bs',values='hr').dropna()
    res[f'{cl}_rare5k']=dict(pairs=len(t),wilcoxon_P=wilcoxon(t.N,t['T']).pvalue,median_N=t.N.median(),median_T=t['T'].median())
json.dump(res,open(O+'primary_tests.json','w'),indent=1,default=float)
print(json.dumps(res,indent=1,default=float))
