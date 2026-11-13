# Step10 E6 (v0.6a, descriptive) + E10 single-cell descriptive part (v0.6d). No HP+/- tests.
import anndata as ad,numpy as np,pandas as pd,scipy.sparse as sp,json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B='/REPLACE_WITH_PRIVATE_WORKDIR/';O=B+'results/step10/E6/'
IM=['CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4','SI','FABP1']; GD=['GKN1','GKN2','TFF2','PGC','PGA3','LIPF','CHIA','MUC6','CLDN18','GAST']
M2=['ATP4A','ATP4B','CBLIF','KCNE2','KCNJ15','KCNJ16','SLC26A7','CCKBR','ESRRG','VSIG1']; M3=['MUC5AC','TFF1','TFF2','MUC6','MUC1','FUT2','FUT3','B3GNT6']
E9=['MUC5AC','TFF1','TFF2','MUC6','MUC1','FUT2','FUT3','B3GNT6','ST3GAL4','A4GNT','CEACAM1','CEACAM5','CEACAM6']
le=pd.read_csv(B+'results/step10/E2/E2_program_leading_edge.tsv',sep='\t',header=None,names=['dir','g'])
a=ad.read_h5ad(B+'results/a5_scrna/work/annotated.h5ad',backed='r')
sym=a.var['sym'].astype(str).values; want=list(dict.fromkeys(IM+GD+M2+M3+E9+le.g.tolist()))
idx=[i for i,s in enumerate(sym) if s in set(want)]; gs=sym[idx]
obs=a.obs[['celltype','donor','stage','hp']].copy()
Xs=[];Cs=[]
for s in range(0,a.n_obs,8000):
    x=a.X[s:s+8000]; c=a.layers['counts'][s:s+8000]
    Xs.append(sp.csr_matrix(x)[:,idx]); Cs.append((sp.csr_matrix(c)[:,idx]>0).astype(np.float32))
X=sp.vstack(Xs).tocsr(); C=sp.vstack(Cs).tocsr()
X=pd.DataFrame(X.toarray(),columns=gs,index=obs.index).T.groupby(level=0).mean().T  # duplicated symbols averaged
C=pd.DataFrame(C.toarray(),columns=gs,index=obs.index).T.groupby(level=0).max().T
grp=obs.groupby(['donor','celltype'],observed=True); keep=[k for k,v in grp.indices.items() if len(v)>=10]
PB=[];PD=[]
for k in keep:
    ii=grp.indices[k]; PB.append(pd.Series(X.iloc[ii].mean().values,index=X.columns,name=k)); PD.append(pd.Series(C.iloc[ii].mean().values,index=C.columns,name=k))
PB=pd.DataFrame(PB); PD=pd.DataFrame(PD); PB.index=pd.MultiIndex.from_tuples(PB.index,names=['donor','celltype']); PD.index=PB.index
mean_ct=PB.groupby(level='celltype').mean(); det_ct=PD.groupby(level='celltype').mean()
mean_ct.T.to_csv(O+'E6_mean_expr_by_celltype_donoravg.csv'); det_ct.T.to_csv(O+'E6_detect_frac_by_celltype_donoravg.csv')
z=(mean_ct-mean_ct.mean())/mean_ct.std(ddof=1).replace(0,np.nan)
up=[g for g in le[le.dir=='UP'].g if g in z.columns]; dn=[g for g in le[le.dir=='DOWN'].g if g in z.columns]
agg=pd.DataFrame({'E2_LE_up(n=%d)'%len(up):z[up].mean(1),'E2_LE_down(n=%d)'%len(dn):z[dn].mean(1)})
agg.to_csv(O+'E6_E2program_aggregate_z_by_celltype.csv')
genes=[g for g in dict.fromkeys(IM+GD+M2+M3) if g in mean_ct.columns]
order=['Epi_Pit','Epi_Neck','Epi_Chief','Epi_Parietal','Epi_IM','Epi_Tumor','Endocrine','T','B','Plasma','Myeloid','Mast','Fibroblast','Endothelial','SMC_pericyte']
order=[o for o in order if o in mean_ct.index]
fig,ax=plt.subplots(figsize=(14,6.5))
for i,g in enumerate(genes):
    zz=z.loc[order,g]; dd=det_ct.loc[order,g]
    ax.scatter([i]*len(order),range(len(order)),s=dd*120+1,c=zz,cmap='viridis',vmin=-2,vmax=2,edgecolors='none')
for j,c in enumerate(agg.columns):
    ax.scatter([len(genes)+1+j]*len(order),range(len(order)),s=60,c=agg.loc[order,c],cmap='viridis',vmin=-2,vmax=2,marker='s')
ax.set_xticks(list(range(len(genes)))+[len(genes)+1,len(genes)+2]); ax.set_xticklabels(genes+list(agg.columns),rotation=90,fontsize=7)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order,fontsize=8); ax.invert_yaxis()
for x0,lab in [(0,'M1 IM'),(len(IM),'M1 GD'),(len(IM)+len(GD),'M2'),(len(IM)+len(GD)+len(M2),'M3')]: ax.axvline(x0-0.5,color='lightgrey',lw=0.5); ax.text(x0,-1.2,lab,fontsize=8)
sm=plt.cm.ScalarMappable(cmap='viridis',norm=plt.Normalize(-2,2)); fig.colorbar(sm,ax=ax,label='z of donor-averaged mean expression (across cell types)')
ax.set_title('E6 (descriptive): GSE249874 donor pseudobulk; dot size = detection fraction',fontsize=9)
fig.tight_layout(); fig.savefig(O+'E6_dotplot.pdf'); fig.savefig(O+'E6_dotplot.png',dpi=200)
# E10 single-cell descriptive: epithelial donor pseudobulk by stage
epi=obs.celltype.astype(str).str.startswith('Epi_'); dg=obs[epi].groupby('donor',observed=True)
E10=pd.DataFrame({d:X.loc[obs.index[epi]].iloc[[obs.index[epi].get_loc(i) for i in v]].mean() for d,v in dg.groups.items()}).T
meta=obs[epi].groupby('donor',observed=True)[['stage','hp']].first()
zz=(E10-E10.mean())/E10.std(ddof=1)
sc=pd.DataFrame({'M1':zz[[g for g in IM if g in zz]].mean(1)-zz[[g for g in GD if g in zz]].mean(1),'M2':zz[[g for g in M2 if g in zz]].mean(1),'M3':zz[[g for g in M3 if g in zz]].mean(1)})
T=pd.concat([meta,sc,E10[[g for g in E9 if g in E10]]],axis=1); T.to_csv(O+'E10_sc_epithelial_donor_pseudobulk.csv')
summ=T.groupby('stage',observed=True)[['M1','M2','M3']+[g for g in E9 if g in T]].median()
summ.to_csv(O+'E10_sc_stage_medians.csv')
json.dump(dict(label='E6 descriptive (v0.6a); E10 single-cell descriptive (v0.6d); no HP+/- tests',n_cells=int(a.n_obs),n_donors=int(obs.donor.nunique()),
  n_donor_celltype_groups=len(keep),genes_missing=[g for g in want if g not in mean_ct.columns][:50],n_LE_up=len(up),n_LE_down=len(dn),
  stage_medians=summ.round(3).to_dict()),open(O+'E6_E10sc_summary.json','w'),indent=1)
print(summ.round(2).to_string()); print(agg.round(2).to_string())
