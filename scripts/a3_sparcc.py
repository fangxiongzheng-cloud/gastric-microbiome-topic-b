# Protocol v0.4 A3 — primary SparCC networks (FastSpar 1.0.0), bootstrap P, patient-level paired permutation.
# Run: /REPLACE_WITH_ENVIRONMENT_ROOT/py/bin/python scripts/a3_sparcc.py   (<=2 threads; resumable: cached permutations in scratch)
import pandas as pd, numpy as np, subprocess, os, json, itertools, glob, networkx as nx
from networkx.algorithms import community
from statsmodels.stats.multitest import multipletests
SEED=20261030; ITER=100; XIT=10; THR=0.1; NBOOT=1000; NPERM=1000; RHO=0.3; PB=0.01; TH='2'
E='/REPLACE_WITH_ENVIRONMENT_ROOT/envs/a3net/bin/'; S='/REPLACE_WITH_PRIVATE_SCRATCH/scratch_a3/'; R='results/a3_network/'; os.makedirs(S,exist_ok=True)
X=pd.read_csv(R+'input/counts_filtered.csv',index_col=0); md=pd.read_csv(R+'input/samples.csv',index_col=0).loc[X.index]
F=pd.read_csv(R+'input/features.csv'); X.columns=F.id; name=dict(zip(F.id,F.taxon)); rid={v:k for k,v in name.items()}
MOD=[rid[m] for m in ['Fusobacterium','Streptococcus','Prevotella','Peptostreptococcus']]; HP=rid['Helicobacter']
MODPAIRS=list(itertools.combinations(MOD,2)); ids=list(F.id); iu=np.triu_indices(len(ids),1)
def write(sub,fn):
    t=sub.T; t.index.name='#OTU ID'; t.to_csv(fn,sep='\t')
def sparcc(fn,out,seed):
    subprocess.run([E+'fastspar','-c',fn,'-r',out+'.cor','-a',out+'.cov','-i',str(ITER),'-x',str(XIT),'-e',str(THR),'-t',TH,'-s',str(seed),'-y'],check=True,stdout=subprocess.DEVNULL)
    return pd.read_csv(out+'.cor',sep='\t',index_col=0).loc[ids,ids]
def modmean(C): return float(np.mean([C.loc[a,b] for a,b in MODPAIRS]))
def hpor(C): return {name[m]:float(C.loc[HP,m]) for m in MOD}
def graph(C,P=None):
    A=(C.abs()>=RHO)&(~np.eye(len(ids),dtype=bool));
    if P is not None: A&=(P<PB)
    G=nx.Graph(); G.add_nodes_from(ids)
    for i,j in zip(*iu):
        if A.values[i,j]: G.add_edge(ids[i],ids[j],w=float(C.values[i,j]))
    return G
def gmetrics(G):
    m=G.number_of_edges(); pos=sum(1 for *_,d in G.edges(data=True) if d['w']>0)
    Gp=nx.Graph([(u,v) for u,v,d in G.edges(data=True) if d['w']>0]); Gp.add_nodes_from(G.nodes)
    mod=community.modularity(Gp,community.greedy_modularity_communities(Gp)) if Gp.number_of_edges()>0 else np.nan
    lcc=max(nx.connected_components(G),key=len); L=G.subgraph(lcc)
    return dict(n_edges=m,density=nx.density(G),pos_edge_frac=(pos/m if m else np.nan),
      mean_abs_rho=(float(np.mean([abs(d['w']) for *_,d in G.edges(data=True)])) if m else np.nan),
      transitivity=nx.transitivity(G),modularity_pos=mod,lcc_size=len(lcc),
      lcc_avg_path=(nx.average_shortest_path_length(L) if len(lcc)>1 else np.nan))
res={}; obs={}
for grp in ['N','T']:
    fn=S+f'obs_{grp}.tsv'; write(X[md.type_bs=='N' if grp=='N' else md.type_bs=='T'],fn)
    C=sparcc(fn,S+f'obs_{grp}',SEED); obs[grp]=C; C.to_csv(R+f'sparcc_cor_{grp}.csv')
    bd=S+f'boot_{grp}/'; os.makedirs(bd,exist_ok=True)
    if not os.path.exists(bd+'DONE'):
        subprocess.run([E+'fastspar_bootstrap','-c',fn,'-n',str(NBOOT),'-p',bd+'b','-s',str(SEED),'-t',TH],check=True,stdout=subprocess.DEVNULL)
        bs=sorted(glob.glob(bd+'b_*.tsv')); os.makedirs(bd+'cor',exist_ok=True)
        for k,b in enumerate(bs):
            sparcc(b,bd+'cor/tmp',SEED+k+1); os.replace(bd+'cor/tmp.cor',bd+f'cor/cor_{k}.tsv')
        open(bd+'DONE','w').close()
    subprocess.run([E+'fastspar_pvalues','-c',fn,'-r',S+f'obs_{grp}.cor','-p',bd+'cor/cor_','-n',str(NBOOT),'-o',S+f'pval_{grp}.tsv','-t',TH],check=True,stdout=subprocess.DEVNULL)
    P=pd.read_csv(S+f'pval_{grp}.tsv',sep='\t',index_col=0).loc[ids,ids]; P.to_csv(R+f'sparcc_pval_{grp}.csv'); obs[grp+'_P']=P
    rows=[dict(a=name[ids[i]],b=name[ids[j]],rho=C.values[i,j],P_boot=P.values[i,j]) for i,j in zip(*iu)]
    ed=pd.DataFrame(rows); ed['q_BH']=multipletests(ed.P_boot,method='fdr_bh')[1]
    ed['edge']=(ed.rho.abs()>=RHO)&(ed.P_boot<PB); ed.sort_values('rho',key=abs,ascending=False).to_csv(R+f'edges_sparcc_{grp}_all_pairs.csv',index=False)
    ed[ed.edge].sort_values('rho',key=abs,ascending=False).to_csv(R+f'edges_sparcc_{grp}.csv',index=False)
# observed statistics
def stats(CN,CT):
    GN,GT=graph(CN),graph(CT); mN,mT=gmetrics(GN),gmetrics(GT)
    hN,hT=hpor(CN),hpor(CT); dN,dT=dict(GN.degree()),dict(GT.degree())
    return dict(mod=modmean(CT)-modmean(CN), hp={k:hT[k]-hN[k] for k in hN}, hpmean=np.mean(list(hT.values()))-np.mean(list(hN.values())),
      glob={k:mT[k]-mN[k] for k in mN}, deg={k:dT[k]-dN[k] for k in ids}, diff=(CT.values-CN.values)[iu])
o=stats(obs['N'],obs['T'])
# permutations: swap T/N within patient (each patient flipped w.p. 0.5)
pf=S+'perm_cache.npz'; rng=np.random.default_rng(SEED)
pids=md.pid.unique(); flips=rng.integers(0,2,size=(NPERM,len(pids))).astype(bool)
cache=dict(np.load(pf,allow_pickle=True)) if os.path.exists(pf) else {}
done=int(cache.get('n',0)); keys=['mod','hpmean']+[f'hp_{k}' for k in o['hp']]+[f'glob_{k}' for k in o['glob']]+[f'deg_{k}' for k in ids]
arr={k:list(cache.get(k,[])) for k in keys}; diffs=list(cache.get('diffs',[]))
for p in range(done,NPERM):
    fl=dict(zip(pids,flips[p])); lab=np.where([fl[x] for x in md.pid],np.where(md.type_bs=='T','N','T'),md.type_bs)
    Cs={}
    for grp in ['N','T']:
        fn=S+f'perm_{grp}.tsv'; write(X[lab==grp],fn); Cs[grp]=sparcc(fn,S+f'perm_{grp}',SEED+100000+p)
    s=stats(Cs['N'],Cs['T'])
    arr['mod'].append(s['mod']); arr['hpmean'].append(s['hpmean'])
    for k,v in s['hp'].items(): arr[f'hp_{k}'].append(v)
    for k,v in s['glob'].items(): arr[f'glob_{k}'].append(v)
    for k,v in s['deg'].items(): arr[f'deg_{k}'].append(v)
    diffs.append(s['diff'])
    if (p+1)%50==0 or p+1==NPERM:
        np.savez(pf,n=p+1,diffs=np.array(diffs),**{k:np.array(v) for k,v in arr.items()}); print('perm',p+1,flush=True)
def pp(obsv,null):
    null=np.asarray(null,float); null=null[~np.isnan(null)]
    if np.isnan(obsv): return dict(obs=None,P_two=None)
    return dict(obs=float(obsv),null_mean=float(null.mean()),null_sd=float(null.std()),
      P_two=float((1+np.sum(np.abs(null-0)>=abs(obsv)-1e-12))/(1+len(null))),
      P_one_TgtN=float((1+np.sum(null>=obsv-1e-12))/(1+len(null))),n_perm=int(len(null)))
prim=pp(o['mod'],arr['mod'])
prim.update(module_pairs=[(name[a],name[b]) for a,b in MODPAIRS],mean_rho_T=modmean(obs['T']),mean_rho_N=modmean(obs['N']),
  pair_rho_T={f'{name[a]}-{name[b]}':float(obs['T'].loc[a,b]) for a,b in MODPAIRS},
  pair_rho_N={f'{name[a]}-{name[b]}':float(obs['N'].loc[a,b]) for a,b in MODPAIRS},
  success=bool(prim['P_two']<0.05 and o['mod']>0), test='two-sided patient-level paired label permutation (each patient T/N swapped w.p. 0.5), statistic = mean SparCC rho(T) - mean SparCC rho(N) over 6 oral-module pairs')
# secondaries
sec=[]
sec.append(dict(family='Helicobacter-oral',metric='mean rho Helicobacter vs 4 oral genera (T-N)',N=float(np.mean(list(hpor(obs['N']).values()))),T=float(np.mean(list(hpor(obs['T']).values()))),**pp(o['hpmean'],arr['hpmean'])))
for k in o['hp']: sec.append(dict(family='Helicobacter-oral',metric=f'rho Helicobacter-{k} (T-N)',N=hpor(obs['N'])[k],T=hpor(obs['T'])[k],**pp(o['hp'][k],arr[f'hp_{k}'])))
mN,mT=gmetrics(graph(obs['N'])),gmetrics(graph(obs['T']))
for k in o['glob']: sec.append(dict(family='global (|rho|>=0.3 graph)',metric=f'{k} (T-N)',N=mN[k],T=mT[k],**pp(o['glob'][k],arr[f'glob_{k}'])))
dN,dT=dict(graph(obs['N']).degree()),dict(graph(obs['T']).degree())
for k in ids: sec.append(dict(family='degree (|rho|>=0.3 graph)',metric=f'degree {name[k]} (T-N)',N=dN[k],T=dT[k],**pp(o['deg'][k],arr[f'deg_{k}'])))
S2=pd.DataFrame(sec)
for fam in S2.family.unique():
    m=S2.family==fam; v=S2.loc[m,'P_two'].notna(); idx=S2[m][v].index
    S2.loc[idx,'q_BH_family']=multipletests(S2.loc[idx,'P_two'],method='fdr_bh')[1]
S2.to_csv(R+'secondary_tests.csv',index=False)
# differential associations (per-pair permutation)
D=np.array(diffs); od=o['diff']; Pd=(1+(np.abs(D)>=np.abs(od)-1e-12).sum(0))/(1+len(D))
dd=pd.DataFrame(dict(a=[name[ids[i]] for i in iu[0]],b=[name[ids[j]] for j in iu[1]],rho_N=obs['N'].values[iu],rho_T=obs['T'].values[iu],diff_TminusN=od,P_perm=Pd))
dd['q_BH']=multipletests(dd.P_perm,method='fdr_bh')[1]; dd.sort_values('P_perm').to_csv(R+'differential_associations.csv',index=False)
# bootstrap-filtered network metrics (descriptive)
bf={g:gmetrics(graph(obs[g],obs[g+'_P'])) for g in ['N','T']}
degbf={g:{name[k]:v for k,v in graph(obs[g],obs[g+'_P']).degree()} for g in ['N','T']}
pd.DataFrame(degbf).to_csv(R+'degree_bootfiltered.csv')
out=dict(primary=prim,params=dict(tool='FastSpar 1.0.0 (C++ SparCC)',iterations=ITER,exclude_iterations=XIT,exclusion_threshold=THR,n_bootstrap=NBOOT,n_perm=NPERM,edge_rule='|rho|>=0.3 & bootstrap P<0.01',seed=SEED,n_pairs=int(len(pids)),n_features=len(ids)),
  network_metrics_bootfiltered=bf,n_edges_bootfiltered={g:int(bf[g]['n_edges']) for g in bf})
json.dump(out,open(R+'primary_test.json','w'),indent=1,default=float); print(json.dumps(out['primary'],indent=1,default=float))
