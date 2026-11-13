import pandas as pd, numpy as np
D='/REPLACE_WITH_RIGHTS_CLEARED_LEGACY_INPUTS/'
out=[]
def P(*a): s=' '.join(str(x) for x in a); print(s); out.append(s)
m=pd.read_csv(D+'myco/metadata_genus_WIS_overlapping_fungi_bacteria_14494samples.tsv',sep='\t',index_col=0)
c=pd.read_csv(D+'myco/count_data_genus_raw_WIS_overlapping_fungi_bacteria_14494samples.tsv',sep='\t',index_col=0)
s=m[m.investigation.str.contains('STAD')]
P('## NH2022 WIS bacteria genera table: STAD samples',len(s),'cases',s.tcga_case_id.nunique(),'ngenera',c.shape[1])
P(pd.crosstab(s.sample_type,s.experimental_strategy).to_string())
P(s.data_submitting_center_label.value_counts().to_string())
P('plates', s.filename.str.extract(r'TCGA-..-....-...-...-(\d{4})')[0].nunique() if 'filename' in s else 'NA')
pt=s[s.sample_type=='Primary Tumor']
for strat in ['RNA-Seq','WGS']:
  ids=pt[pt.experimental_strategy==strat].index.intersection(c.index); x=c.loc[ids]
  P(f'### PT {strat} n={len(ids)} median total reads(WIS genera)={x.sum(1).median()}')
  for g in ['Fusobacterium','Streptococcus','Prevotella','Veillonella','Peptostreptococcus','Porphyromonas','Lactobacillus','Malassezia']:
    if g in x: P(g,'prev>0=%.2f'%(x[g]>0).mean(),'median=%s'%x[g].median())
f=pd.read_csv(D+'myco/count_data_fungi_decontaminated_raw.tsv',sep='\t',index_col=0)
P('fungi decontam table STAD samples:',len(f.index.intersection(s.index)))
# Poore raw
pm=pd.read_csv(D+'poore2020/Metadata-TCGA-Kraken-17625-Samples.csv',index_col=0)
ps=pm[pm.disease_type.str.contains('Stomach',case=False)]
P('## Poore2020 (retracted) STAD samples',len(ps),'cases',ps.case_uuid.nunique() if 'case_uuid' in ps else '')
P(pd.crosstab(ps.sample_type,ps.experimental_strategy).to_string())
hdr=pd.read_csv(D+'poore2020/Kraken-TCGA-Raw-Data-17625-Samples.csv',nrows=0).columns
hc=[h for h in hdr if 'Helicobacter' in h]; P('Helicobacter cols',hc)
r=pd.read_csv(D+'poore2020/Kraken-TCGA-Raw-Data-17625-Samples.csv',index_col=0,usecols=[hdr[0]]+hc)
r=r.loc[r.index.intersection(ps.index)]
for st in ['Primary Tumor','Solid Tissue Normal']:
  for strat in ['RNA-Seq','WGS']:
    ids=ps[(ps.sample_type==st)&(ps.experimental_strategy==strat)].index.intersection(r.index)
    v=r.loc[ids,hc[-1]]; P(st,strat,'n',len(ids),'Helicobacter prev>0=%.2f median=%s q90=%s'%((v>0).mean(),v.median(),v.quantile(.9)))
open('s0_counts.md','w').write('\n'.join(out))
