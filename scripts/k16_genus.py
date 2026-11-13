import pandas as pd, glob
D='data/external/korea16s/'; O='results/korea16s/'
t=pd.read_csv(D+'exp_table/table.tsv',sep='\t',skiprows=1,index_col=0)
st=pd.concat([pd.read_csv(f+'/stats.tsv',sep='\t',skiprows=[1],index_col=0) for f in glob.glob(D+'exp_stats/*')]); st.to_csv(O+'dada2_stats.csv')
for n in ['silva','gg2v4']:
    tx=pd.read_csv(D+f'exp_tax_{n}/taxonomy.tsv',sep='\t',index_col=0)['Taxon']
    tx=tx.reindex(t.index).fillna('Unassigned')
    bad=tx.str.contains('Mitochondria|Chloroplast|mitochondria|chloroplast')|~tx.str.contains('p__')|tx.str.startswith('Unassigned')|tx.str.contains('Eukaryota')
    tt=t[~bad]
    lv=tx[~bad].str.split(';').map(lambda x:';'.join([y.strip() for y in x[:6]]))
    g=tt.groupby(lv).sum().T
    g.to_csv(O+f'genus_{n}.tsv',sep='\t')
    hp=tx[tx.str.contains('Helicobacter')]
    print(n,'ASVs',len(t),'removed',bad.sum(),'Helicobacter ASVs',len(hp)); print(hp.value_counts().head())
    tx.to_csv(O+f'asv_taxonomy_{n}.tsv',sep='\t')
