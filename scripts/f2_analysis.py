# F2: PRJNA1061213 one-sided paired Wilcoxon (N>T) on Helicobacter CLR; decisions in results/f2/F2_processing_decisions_20261109.md
import pandas as pd,numpy as np,json,glob
from scipy.stats import wilcoxon
B='/REPLACE_WITH_PRIVATE_WORKDIR/';A=B+'data/external/third_cohort/PRJNA1061213/art/';O=B+'results/f2/'
t=pd.read_csv(A+'exp_table/table.tsv',sep='\t',skiprows=1,index_col=0)
runs=pd.read_csv(B+'data/external/third_cohort/PRJNA1061213_runs.tsv',sep='\t').set_index('run_accession')
st=pd.read_csv(A+'exp_stats/stats.tsv',sep='\t',skiprows=[1],index_col=0); st.to_csv(O+'dada2_stats.csv')
G={};
for n in ['silva','gg2v4']:
    tx=pd.read_csv(A+f'exp_tax_{n}/taxonomy.tsv',sep='\t',index_col=0)['Taxon'].reindex(t.index).fillna('Unassigned')
    bad=tx.str.contains('Mitochondria|Chloroplast|mitochondria|chloroplast')|~tx.str.contains('p__')|tx.str.startswith('Unassigned')|tx.str.contains('Eukaryota')
    lv=tx[~bad].str.split(';').map(lambda x:';'.join([y.strip() for y in x[:6]]))
    G[n]=t[~bad].groupby(lv).sum().T; G[n].to_csv(O+f'genus_{n}.tsv',sep='\t')
mind=pd.concat({n:G[n].sum(1) for n in G},axis=1).reindex(t.columns).min(axis=1).fillna(0)
alias=runs.sample_alias; pid=alias.str.extract(r'(\d+)')[0]; typ=alias.str[:2].map({'CT':'T','NT':'N'})
res={}
def clr_df(g):
    hc=[c for c in g.columns if c.split(';')[-1].strip().replace('g__','').startswith('Helicobacter')]
    h=g[hc].sum(1); g2=g.drop(columns=hc).assign(HELICO=h); lg=np.log(g2+0.5); return pd.DataFrame({'reads':g.sum(1),'hp':h,'ra':h/g.sum(1),'clr':lg['HELICO']-lg.mean(1)}),hc
for n in ['silva','gg2v4']:
    df,hc=clr_df(G[n]); df['pid']=pid.reindex(df.index);df['type']=typ.reindex(df.index);df['mindepth']=mind.reindex(df.index)
    df.to_csv(O+f'per_sample_{n}.csv')
    for lab,thr in [('main_ge5000',5000),('sens_ge1000',1000)]:
        ok=df[df.mindepth>=thr]; w=ok.pivot_table(index='pid',columns='type',values=['clr','ra','hp']).dropna()
        d=w['clr']['N']-w['clr']['T']; nz=d[d!=0]; rk=nz.abs().rank(); r=(rk[nz>0].sum()-rk[nz<0].sum())/rk.sum()
        p1=wilcoxon(w['clr']['N'],w['clr']['T'],alternative='greater'); p2=wilcoxon(w['clr']['N'],w['clr']['T'])
        rng=np.random.default_rng(20261030);v=nz.values;bs=[]
        for _ in range(2000):
            x=rng.choice(v,len(v));rr=pd.Series(np.abs(x)).rank().values;p=rr[x>0].sum();bs.append((p-(rr.sum()-p))/rr.sum())
        npairs_all=int(df.pivot_table(index='pid',columns='type',values='clr').dropna().shape[0])
        res[f'{n}_{lab}']=dict(hp_asv_taxa=hc,pairs_total_in_data=npairs_all,pairs_used=len(w),pairs_dropped_depth=npairs_all-len(w),N_gt_T=int((d>0).sum()),T_gt_N=int((d<0).sum()),ties=int((d==0).sum()),
         wilcoxon_one_sided_P=float(p1.pvalue),wilcoxon_two_sided_P=float(p2.pvalue),rank_biserial_r=float(r),r_CI95_boot=[float(np.percentile(bs,2.5)),float(np.percentile(bs,97.5))],
         median_clr_N=float(w['clr']['N'].median()),median_clr_T=float(w['clr']['T'].median()),median_ra_N=float(w['ra']['N'].median()),median_ra_T=float(w['ra']['T'].median()),
         median_fold_ra_N_over_T_plus1e5=float(np.exp(np.median(np.log(w['ra']['N']+1e-5)-np.log(w['ra']['T']+1e-5)))))
primary=res['silva_main_ge5000']; res['F2_PRIMARY']=dict(spec='SILVA, >=5000 reads (min across classifiers), one-sided paired Wilcoxon N>T',pass_=bool(primary['wilcoxon_one_sided_P']<0.05 and primary['N_gt_T']>primary['T_gt_N']))
res['n_runs_in_table']=int(t.shape[1]); res['median_reads_silva']=float(G['silva'].sum(1).median()); res['samples_lt5000_either']=int((mind<5000).sum())
json.dump(res,open(O+'F2_results.json','w'),indent=1,default=float); print(json.dumps(res,indent=1,default=float))
# genus-level specificity (descriptive): paired effect (rank-biserial r, N>T count) for all genera prevalent >=10% in main set
df,_=clr_df(G['silva']); ok=df.index[mind.reindex(df.index)>=5000]
g=G['silva'].loc[ok]; lg=np.log(g+0.5); clr=lg.sub(lg.mean(1),axis=0); clr['pid']=pid.reindex(clr.index); clr['type']=typ.reindex(clr.index)
rows=[]
for c in g.columns:
    if (g[c]>0).mean()<0.10: continue
    w=clr.pivot_table(index='pid',columns='type',values=c).dropna(); d=w.N-w['T']; nz=d[d!=0]
    if len(nz)<10: continue
    rk=nz.abs().rank(); rows.append(dict(genus=c.split(';')[-1],pairs=len(w),N_gt_T=int((d>0).sum()),T_gt_N=int((d<0).sum()),r=(rk[nz>0].sum()-rk[nz<0].sum())/rk.sum()))
S=pd.DataFrame(rows).sort_values('r',ascending=False).reset_index(drop=True); S['rank_by_r']=S.index+1; S.to_csv(O+'F2_genus_specificity_descriptive.csv',index=False)
print(S.head(8)); print(S[S.genus.str.contains('Helicobacter')])
