# Draft figures for A3 (matplotlib/networkx; no AI images). Okabe-Ito palette.
import pandas as pd, numpy as np, json, networkx as nx, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R='results/a3_network/'; os_=__import__('os'); os_.makedirs(R+'figs',exist_ok=True)
F=pd.read_csv(R+'input/features.csv'); short=lambda t: t if ';' not in t else t.split(';')[-1].split('__')[-1]+' (unc.)'
MOD=['Fusobacterium','Streptococcus','Prevotella','Peptostreptococcus']; OR='#E69F00'; HPc='#D55E00'; POS='#0072B2'; NEG='#CC79A7'
E={g:pd.read_csv(R+f'edges_sparcc_{g}.csv') for g in 'NT'}
G0=nx.Graph(); G0.add_nodes_from(F.taxon)
for g in 'NT':
    for r in E[g].itertuples(): G0.add_edge(r.a,r.b)
conn=[n for n in G0 if G0.degree(n)>0]; iso=[n for n in G0 if G0.degree(n)==0]
pos=nx.kamada_kawai_layout(G0.subgraph(conn)); pos={k:v*0.8 for k,v in pos.items()}
for i,n in enumerate(iso): pos[n]=np.array([1.15*np.cos(2*np.pi*i/len(iso)),1.15*np.sin(2*np.pi*i/len(iso))])
cache=np.load('/REPLACE_WITH_PRIVATE_SCRATCH/scratch_a3/perm_cache.npz'); J=json.load(open(R+'primary_test.json'))['primary']
S=pd.read_csv(R+'secondary_tests.csv')
fig=plt.figure(figsize=(15,10))
for k,g in enumerate('NT'):
    ax=fig.add_subplot(2,3,k+1); G=nx.Graph(); G.add_nodes_from(F.taxon)
    for r in E[g].itertuples(): G.add_edge(r.a,r.b,w=r.rho)
    ec=[POS if d['w']>0 else NEG for *_,d in G.edges(data=True)]; ew=[1+6*(abs(d['w'])-0.3) for *_,d in G.edges(data=True)]
    nc=[OR if n in MOD else HPc if n=='Helicobacter' else '#BBBBBB' for n in G.nodes]; ns=[60+60*G.degree(n) for n in G.nodes]
    nx.draw_networkx_edges(G,pos,ax=ax,edge_color=ec,width=ew,alpha=.8); nx.draw_networkx_nodes(G,pos,ax=ax,node_color=nc,node_size=ns,edgecolors='k',linewidths=.4)
    lab={n:short(n) for n in G.nodes if G.degree(n)>0 or n in MOD or n=='Helicobacter'}
    nx.draw_networkx_labels(G,pos,labels=lab,ax=ax,font_size=5.5)
    ax.set_title(f"{'abcd'[k]}  SparCC network, {'adjacent normal' if g=='N' else 'tumour'} (n=192)\n{G.number_of_edges()} edges (|rho|>=0.3, bootstrap P<0.01)",fontsize=9); ax.axis('off')
ax=fig.add_subplot(2,3,3); ax.hist(cache['mod'],bins=40,color='#999999'); ax.axvline(J['obs'],color=HPc,lw=2)
ax.set_xlabel('mean oral-module SparCC rho, T - N'); ax.set_ylabel('permutations'); ax.set_title(f"c  Primary test: obs {J['obs']:.3f}, two-sided P={J['P_two']:.2f}\n(1000 within-patient T/N swaps)",fontsize=9)
ax=fig.add_subplot(2,3,4); pr=list(J['pair_rho_N']); x=np.arange(len(pr))
ax.bar(x-.2,[J['pair_rho_N'][p] for p in pr],.4,label='Normal',color=POS); ax.bar(x+.2,[J['pair_rho_T'][p] for p in pr],.4,label='Tumour',color=OR)
ax.set_xticks(x); ax.set_xticklabels([p.replace('-','-\n') for p in pr],fontsize=6); ax.axhline(.3,ls=':',c='k',lw=.8); ax.set_ylabel('SparCC rho'); ax.legend(fontsize=7); ax.set_title('d  Oral-module pairwise correlations',fontsize=9)
ax=fig.add_subplot(2,3,5); h=S[S.family=='Helicobacter-oral'].iloc[1:]; x=np.arange(len(h))
ax.bar(x-.2,h.N,.4,color=POS,label='Normal'); ax.bar(x+.2,h['T'],.4,color=OR,label='Tumour')
ax.set_xticks(x); ax.set_xticklabels([m.split('-')[1].split(' ')[0] for m in h.metric],fontsize=7); ax.axhline(-.3,ls=':',c='k',lw=.8); ax.set_ylabel('SparCC rho with Helicobacter')
for i,p in enumerate(h.P_two): ax.text(i,-0.29,f'P={p:.2f}',ha='center',fontsize=6)
ax.legend(fontsize=7); ax.set_title('e  Helicobacter vs oral genera',fontsize=9)
ax=fig.add_subplot(2,3,6); d=S[S.family.str.startswith('degree')].copy(); d['tax']=d.metric.str.replace('degree ','').str.replace(' (T-N)','')
d=d[(d.N>0)|(d['T']>0)].sort_values('T'); y=np.arange(len(d))
ax.scatter(d.N,y,c=POS,label='Normal',s=15); ax.scatter(d['T'],y,c=OR,label='Tumour',s=15)
ax.set_yticks(y); ax.set_yticklabels([short(t) for t in d.tax],fontsize=5.5); ax.set_xlabel('degree'); ax.legend(fontsize=7); ax.set_title('f  Degree centrality (no difference q<0.05)',fontsize=9)
plt.tight_layout(); plt.savefig(R+'figs/FigA3_network_draft.pdf'); plt.savefig(R+'figs/FigA3_network_draft.png',dpi=150); print('ok')
