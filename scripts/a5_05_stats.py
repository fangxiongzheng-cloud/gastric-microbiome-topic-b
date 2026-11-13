"""A5 step 5: primary test (epithelial->myeloid CellChat strength, HP+ vs HP-, stage-stratified exact permutation) + LR overlap secondary."""
import os, json, itertools, numpy as np, pandas as pd
from scipy.stats import hypergeom
O='/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna'; SEED=20261030
des=pd.read_csv(O+'/sample_design.csv'); rows=[]
for _,r in des.iterrows():
    f=f'{O}/cellchat/{r.donor}_weight.csv'; w=pd.read_csv(f,index_col=0)
    send=[i for i in w.index if i.startswith('Epi_')]
    s=float(w.loc[send,'Myeloid'].sum()) if ('Myeloid' in w.columns and send) else 0.0
    rows.append(dict(donor=r.donor,stage=r.stage,hp=r.hp,strength=s,epi_senders=';'.join(send),has_myeloid='Myeloid' in w.columns))
df=pd.DataFrame(rows); df.to_csv(O+'/primary_donor_strength.csv',index=False)
stages=['GS','IM','GC']
def stat(lbl):  # lbl: dict stage-> boolean array HP+ over df rows of that stage
    return np.mean([df.strength[(df.stage==st).values][lbl[st]].mean()-df.strength[(df.stage==st).values][~lbl[st]].mean() for st in stages])
obs_lbl={st:(df.hp[df.stage==st]=='HP+').values for st in stages}; T=stat(obs_lbl)
combos={st:[np.isin(np.arange(6),c) for c in itertools.combinations(range(6),3)] for st in stages}
null=np.array([stat(dict(zip(stages,c))) for c in itertools.product(*[combos[s] for s in stages])])
p=float((np.abs(null)>=abs(T)-1e-12).mean())
rng=np.random.default_rng(SEED); boots=[]
for _ in range(2000):
    v=[]
    for st in stages:
        a=df.strength[(df.stage==st)&(df.hp=='HP+')].values; b=df.strength[(df.stage==st)&(df.hp=='HP-')].values
        v.append(rng.choice(a,3).mean()-rng.choice(b,3).mean())
    boots.append(np.mean(v))
per=df.groupby(['stage','hp']).strength.agg(['mean','median','min','max']).reset_index()
res=dict(statistic='mean over stages of [mean(HP+) - mean(HP-)] epithelial->Myeloid CellChat aggregateNet weight (p<0.05)',
 T_obs=T,p_two_sided_exact=p,n_perm=len(null),ci95_bootstrap=[float(np.percentile(boots,2.5)),float(np.percentile(boots,97.5))],
 mean_HPpos=float(df.strength[df.hp=='HP+'].mean()),mean_HPneg=float(df.strength[df.hp=='HP-'].mean()),
 per_stage={st:float(df.strength[(df.stage==st)&(df.hp=='HP+')].mean()-df.strength[(df.stage==st)&(df.hp=='HP-')].mean()) for st in stages},
 donors_without_myeloid_or_epi=df.donor[(~df.has_myeloid)|(df.epi_senders=='')].tolist(),
 success=bool(p<0.05))
json.dump(res,open(O+'/primary_test.json','w'),indent=1); per.to_csv(O+'/primary_per_stage.csv',index=False); np.save(O+'/primary_null.npy',null)
print(json.dumps(res,indent=1)); print(df.to_string())
# ---- secondary: LR pairs epi->myeloid, top20 by stratified mean diff in each tool
def key(l,r): return (str(l).upper().replace('_','+')+'^'+str(r).upper().replace('_','+'))
cc=[]; lz=[]
for d in des.donor:
    a=pd.read_csv(f'{O}/cellchat/{d}_lr_all.csv'); a=a[a.source.astype(str).str.startswith('Epi_')&(a.target=='Myeloid')]
    a['k']=[key(l,r) for l,r in zip(a.ligand,a.receptor)]; g=a.groupby('k').prob.sum().rename(d); cc.append(g)
    b=pd.read_csv(f'{O}/liana/{d}.csv'); b=b[b.source.astype(str).str.startswith('Epi_')&(b.target=='Myeloid')]
    b['k']=[key(l,r) for l,r in zip(b.ligand_complex,b.receptor_complex)]; b['sc']=-np.log10(b.magnitude_rank); lz.append(b.groupby('k').sc.max().rename(d))
CC=pd.concat(cc,axis=1).fillna(0); LZ=pd.concat(lz,axis=1).fillna(0)
def sdiff(M):
    out=[]
    for st in stages:
        pos=des.donor[(des.stage==st)&(des.hp=='HP+')]; neg=des.donor[(des.stage==st)&(des.hp=='HP-')]
        out.append(M[pos].mean(1)-M[neg].mean(1))
    return pd.concat(out,axis=1).mean(1)
dc=sdiff(CC); dl=sdiff(LZ); U=sorted(set(dc.index)&set(dl.index))
dcu=dc.loc[U]; dlu=dl.loc[U]; sec={}
for nm,asc in [('HPpos_up',False),('HPneg_up',True)]:
    tc=set(dcu.sort_values(ascending=asc).index[:20]); tl=set(dlu.sort_values(ascending=asc).index[:20]); ov=sorted(tc&tl)
    sec[nm]=dict(universe=len(U),overlap=len(ov),pairs=ov,hypergeom_p=float(hypergeom.sf(len(ov)-1,len(U),20,20)),expected=400/len(U))
pd.DataFrame(dict(cellchat_diff=dc,liana_diff=dl)).to_csv(O+'/lr_epi_myeloid_diffs.csv')
json.dump(sec,open(O+'/secondary_lr_overlap.json','w'),indent=1); print(json.dumps(sec,indent=1))
