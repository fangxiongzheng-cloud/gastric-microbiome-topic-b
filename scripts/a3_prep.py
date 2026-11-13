# Protocol v0.4 A3 — input preparation (Korea 16S, SILVA genus). Seed 20261030.
# A-priori decisions (fixed before any network was computed):
#  * samples: same depth filter as korea16s (min of SILVA/GG2 bacterial reads >=1000), complete T/N pairs only
#  * features: genus-table columns, named as in scripts/k16_da.R
#  * Prevotella: SILVA genera g__Prevotella, g__Prevotella_<n> (e.g. _7, _9) are summed into one feature
#    "Prevotella" (Prevotella sensu lato); Alloprevotella / Paraprevotella / Prevotellaceae-unclassified/uncultured NOT included
#  * prevalence: count>0 in >=20% of N samples AND >=20% of T samples
#  * oral module = Fusobacterium (g__Fusobacterium), Streptococcus (g__Streptococcus),
#    Prevotella (merged as above), Peptostreptococcus (g__Peptostreptococcus); Helicobacter = g__Helicobacter
import pandas as pd, numpy as np, re, json, os
O='results/korea16s/'; R='results/a3_network/'; os.makedirs(R+'input',exist_ok=True)
ps=pd.read_csv(O+'per_sample_silva.csv',index_col=0); g2=pd.read_csv(O+'per_sample_gg2v4.csv',index_col=0)
keep=ps.index[np.minimum(ps.reads,g2.reads.reindex(ps.index))>=1000]
md=ps.loc[keep,['pid','type_bs']]; vc=md.pid.value_counts(); md=md[md.pid.isin(vc[vc==2].index)]
assert (md.groupby('pid').type_bs.nunique()==2).all()
g=pd.read_csv(O+'genus_silva.tsv',sep='\t',index_col=0).loc[md.index]
gen=[c.split(';')[-1].strip() for c in g.columns]
prev_cols=[c for c,x in zip(g.columns,gen) if re.fullmatch(r'g__Prevotella(_\d+)?',x)]
merged={'Prevotella':prev_cols}
gp=g.drop(columns=prev_cols); gp['Prevotella']=g[prev_cols].sum(1)
nm=[]
for c in gp.columns:
    if c=='Prevotella': nm.append(c); continue
    x=re.sub('^g__','',c.split(';')[-1]); nm.append(c if (x=='' or re.search('uncultured|_ge$',x) or not c.split(';')[-1].startswith('g__')) else x)
gp.columns=pd.Index(nm); assert gp.columns.is_unique
pN=(gp[md.type_bs.values=='N']>0).mean(); pT=(gp[md.type_bs.values=='T']>0).mean()
sel=gp.columns[(pN>=0.2)&(pT>=0.2)]
X=gp[sel].astype(int)
ids={t:f'g{i:03d}' for i,t in enumerate(sel)}
pd.DataFrame({'id':[ids[t] for t in sel],'taxon':sel,'prev_N':pN[sel].values,'prev_T':pT[sel].values}).to_csv(R+'input/features.csv',index=False)
md.to_csv(R+'input/samples.csv'); X.to_csv(R+'input/counts_filtered.csv')
module=['Fusobacterium','Streptococcus','Prevotella','Peptostreptococcus']
info=dict(n_pairs=int(len(md)/2),n_samples=int(len(md)),n_features_total=int(gp.shape[1]),n_features_kept=int(len(sel)),
  prevotella_merged_silva_columns=prev_cols,module=module,module_present={m:bool(m in sel) for m in module},
  helicobacter_present=bool('Helicobacter' in sel),
  module_prev={m:dict(N=float(pN.get(m,np.nan)),T=float(pT.get(m,np.nan))) for m in module+['Helicobacter']},
  reads_kept_frac_median=float((X.sum(1)/g.sum(1)).median()))
json.dump(info,open(R+'input/prep_info.json','w'),indent=1); print(json.dumps(info,indent=1))
