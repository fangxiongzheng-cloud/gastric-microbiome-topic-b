"""A5 draft figures (matplotlib only)."""
import os, json, numpy as np, pandas as pd, anndata as ad, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O='/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna'; W=O+'/work'; F=O+'/figs'; os.makedirs(F,exist_ok=True)
A=ad.read_h5ad(W+'/annotated.h5ad',backed='r'); ob=A.obs.copy(); U=np.asarray(A.obsm['X_umap'])
cts=sorted(ob.celltype.unique()); cmap=plt.get_cmap('tab20'); col={c:cmap(i%20) for i,c in enumerate(cts)}
fig=plt.figure(figsize=(18,15)); gs=fig.add_gridspec(3,3)
ax=fig.add_subplot(gs[0,0])
for c in cts:
    m=(ob.celltype==c).values; ax.scatter(U[m,0],U[m,1],s=0.3,color=col[c],label=f'{c} ({m.sum()})',rasterized=True)
ax.legend(markerscale=12,fontsize=6,loc='best',frameon=False); ax.set_title('a  UMAP: cell types (Harmony)'); ax.set_xticks([]); ax.set_yticks([])
ax=fig.add_subplot(gs[0,1]); sc_={'GS':'#4daf4a','IM':'#ff7f00','GC':'#e41a1c'}
for st in ['GS','IM','GC']:
    m=(ob.stage==st).values; ax.scatter(U[m,0],U[m,1],s=0.3,color=sc_[st],label=st,rasterized=True)
ax.legend(markerscale=12,frameon=False); ax.set_title('b  UMAP: stage'); ax.set_xticks([]); ax.set_yticks([])
ax=fig.add_subplot(gs[0,2]); ct=pd.crosstab(ob.donor,ob.celltype,normalize='index')
order=[f'{s}_{h}_{i}' for s in ['GS','IM','GC'] for h in ['HPN','HPP'] for i in (1,2,3)]; ct=ct.reindex(order)
bot=np.zeros(len(ct))
for c in cts: ax.bar(range(len(ct)),ct[c].values,bottom=bot,color=col[c],width=0.85); bot+=ct[c].values
ax.set_xticks(range(len(ct))); ax.set_xticklabels(ct.index,rotation=90,fontsize=7); ax.set_ylabel('fraction'); ax.set_title('c  composition per donor')
d=pd.read_csv(O+'/primary_donor_strength.csv'); pr=json.load(open(O+'/primary_test.json')); null=np.load(O+'/primary_null.npy')
ax=fig.add_subplot(gs[1,0])
for i,st in enumerate(['GS','IM','GC']):
    for j,hp in enumerate(['HP-','HP+']):
        v=d.strength[(d.stage==st)&(d.hp==hp)].values; x=i*3+j
        ax.scatter(np.full(len(v),x)+np.linspace(-.12,.12,len(v)),v,color=['#377eb8','#e41a1c'][j],zorder=3); ax.hlines(v.mean(),x-.25,x+.25,color='k')
ax.set_xticks([0,1,3,4,6,7]); ax.set_xticklabels(['GS HP-','GS HP+','IM HP-','IM HP+','GC HP-','GC HP+'],fontsize=8)
ax.set_ylabel('epithelial→myeloid strength (CellChat)'); ax.set_title('d  ★primary: donor-level strength')
ax=fig.add_subplot(gs[1,1]); ax.hist(null,bins=60,color='grey'); ax.axvline(pr['T_obs'],color='r')
ax.set_title(f"e  stratified exact permutation (n={pr['n_perm']})\nT={pr['T_obs']:.3g}, P={pr['p_two_sided_exact']:.3g}"); ax.set_xlabel('mean within-stage HP+ − HP−')
ax=fig.add_subplot(gs[1,2]); lr=pd.read_csv(O+'/lr_epi_myeloid_diffs.csv',index_col=0).fillna(0); sec=json.load(open(O+'/secondary_lr_overlap.json'))
ax.scatter(lr.cellchat_diff,lr.liana_diff,s=6,color='grey')
for p in sec['HPpos_up']['pairs']+sec['HPneg_up']['pairs']:
    if p in lr.index: ax.scatter(lr.cellchat_diff[p],lr.liana_diff[p],s=14,color='r'); ax.annotate(p,(lr.cellchat_diff[p],lr.liana_diff[p]),fontsize=6)
ax.axhline(0,lw=.5,color='k'); ax.axvline(0,lw=.5,color='k'); ax.set_xlabel('CellChat Δprob (HP+−HP−)'); ax.set_ylabel('LIANA Δ(−log10 magnitude rank)')
ax.set_title(f"f  epi→myeloid LR: top-20 overlap {sec['HPpos_up']['overlap']} (HP+↑), {sec['HPneg_up']['overlap']} (HP−↑)",fontsize=9)
tp=O+'/trajectory/pseudotime_cells.csv'
if os.path.exists(tp):
    t=pd.read_csv(tp); ptc=[c for c in t.columns if c.startswith('pt_Lineage')]
    ax=fig.add_subplot(gs[2,0]); ax.scatter(U[:,0],U[:,1],s=0.2,color='lightgrey',rasterized=True)
    pmax=t[ptc].max(1); s_=ax.scatter(t.UMAP1,t.UMAP2,c=pmax,s=0.5,cmap='viridis',rasterized=True); plt.colorbar(s_,ax=ax,label='slingshot pseudotime (max)')
    ax.set_title('g  epithelial trajectory (epithelial-subset UMAP)'); ax.set_xticks([]); ax.set_yticks([])
    ax=fig.add_subplot(gs[2,1]); pg=pd.read_csv(O+'/trajectory/condiments_progressionTest.csv')
    for k,c in enumerate(ptc):
        for hp,cl in [('HP-','#377eb8'),('HP+','#e41a1c')]:
            v=t.loc[(t.hp==hp)&t[c].notna(),c]; h,e=np.histogram(v,bins=50,density=True); ax.plot((e[1:]+e[:-1])/2,h+k*h.max()*0 + k*0.0,color=cl,ls=['-','--',':','-.'][k%4],label=f'{c[3:]} {hp}')
    ax.legend(fontsize=6,frameon=False); ax.set_xlabel('pseudotime'); ax.set_ylabel('density'); ax.set_title('h  HP+ vs HP− along pseudotime (condiments)')
    ax=fig.add_subplot(gs[2,2]); ax.axis('off'); ax.text(0,1,pg.round(4).to_string(index=False),family='monospace',fontsize=7,va='top'); ax.set_title('i  condiments progressionTest')
fig.tight_layout(); fig.savefig(F+'/FigA5_scrna_draft.pdf',dpi=200); fig.savefig(F+'/FigA5_scrna_draft.png',dpi=150)
# supplementary QC
qc=pd.read_csv(O+'/qc_summary.csv'); kn=pd.read_csv(W+'/knee_points.csv')
fig,ax=plt.subplots(1,2,figsize=(13,4.5))
for s,g in kn.groupby('sample'): ax[0].loglog(g['rank'],g.umi,lw=.7)
ax[0].set_xlabel('barcode rank'); ax[0].set_ylabel('UMI'); ax[0].set_title('a  barcode rank plots (18 samples)')
x=np.arange(len(qc)); ax[1].bar(x-.3,qc.called,.3,label='called (OrdMag+EmptyDrops)'); ax[1].bar(x,qc.after_scrublet,.3,label='after QC+Scrublet'); ax[1].bar(x+.3,qc.final,.3,label='final (cap 4000)')
ax[1].set_xticks(x); ax[1].set_xticklabels(qc['sample'],rotation=90,fontsize=7); ax[1].legend(fontsize=7); ax[1].set_title('b  cells per sample')
fig.tight_layout(); fig.savefig(F+'/FigS_A5_qc_draft.pdf'); fig.savefig(F+'/FigS_A5_qc_draft.png',dpi=150)
print('FIGS_DONE')
