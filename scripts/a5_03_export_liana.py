"""A5 step 3: export per-donor log-normalised matrices for CellChat (R) and run LIANA rank_aggregate per donor. Resumable."""
import scipy.sparse as _sp
for _c in (_sp.csr_array,_sp.csc_array,_sp.coo_array):
    if not hasattr(_c,'getnnz'): _c.getnnz=lambda self,axis=None: self._getnnz(axis=axis)  # scipy 1.17 removed getnnz; scanpy 1.11.5 still calls it
import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.io as sio, scipy.sparse as sp
import liana as li
SEED=20261030; O='/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna'; W=O+'/work'; os.makedirs(W+'/donor',exist_ok=True); os.makedirs(O+'/liana',exist_ok=True)
# trajectory export (epithelial subset)
E=ad.read_h5ad(W+'/epithelial.h5ad'); os.makedirs(W+'/traj',exist_ok=True)
if not os.path.exists(W+'/traj/obs.csv'):
    ob=E.obs[['sample','stage','hp','donor','epi_type']].copy()
    ob=pd.concat([ob,pd.DataFrame(E.obsm['X_harmony'][:,:10],index=E.obs_names,columns=[f'H{i+1}' for i in range(10)]),
                  pd.DataFrame(E.obsm['X_umap'],index=E.obs_names,columns=['UMAP1','UMAP2'])],axis=1)
    Ec=ad.AnnData(E.layers['counts'],obs=E.obs[['sample']].copy(),var=E.var[['sym']].copy())
    sc.pp.highly_variable_genes(Ec,n_top_genes=2000,flavor='seurat_v3'); h=Ec.var.highly_variable.values  # deviation: no batch_key (loess fails for samples with few epithelial cells)
    sio.mmwrite(W+'/traj/counts_hvg.mtx',sp.csr_matrix(E.layers['counts'][:,h]).T.tocoo()); pd.Series(E.var.sym.values[h]).to_csv(W+'/traj/genes.csv',index=False)
    ob.to_csv(W+'/traj/obs.csv')
del E
A=ad.read_h5ad(W+'/annotated.h5ad'); del A.layers['counts']
A.var_names=A.var.sym.astype(str).values; A.var_names_make_unique()
for d in sorted(A.obs.donor.unique()):
    D=A[A.obs.donor==d].copy(); vc=D.obs.celltype.value_counts(); keep=vc.index[vc>=10]; D=D[D.obs.celltype.isin(keep)].copy()
    base=f'{W}/donor/{d}'
    if not os.path.exists(base+'.mtx'):
        sio.mmwrite(base+'.mtx',sp.csr_matrix(D.X).T.tocoo()); pd.Series(D.var_names).to_csv(base+'_genes.csv',index=False)
        D.obs[['celltype','stage','hp','donor']].to_csv(base+'_meta.csv')
    out=f'{O}/liana/{d}.csv'
    if not os.path.exists(out):
        D.obs['celltype']=D.obs.celltype.astype('category')
        li.mt.rank_aggregate(D,groupby='celltype',resource_name='consensus',expr_prop=0.1,use_raw=False,n_perms=1000,seed=SEED,verbose=False,n_jobs=2)
        D.uns['liana_res'].to_csv(out,index=False); print('liana',d,len(D.uns['liana_res']),flush=True)
print('STEP3_DONE')
