"""A5 step 2: QC (v0.3 S4), Scrublet, cap, Harmony integration, marker-based annotation. See results/a5_scrna/PREDECLARED_a5.md."""
import scipy.sparse as _sp
for _c in (_sp.csr_array,_sp.csc_array,_sp.coo_array):
    if not hasattr(_c,'getnnz'): _c.getnnz=lambda self,axis=None: self._getnnz(axis=axis)  # scipy 1.17 removed getnnz; scanpy 1.11.5 still calls it
import os, json, numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp, harmonypy
sc.settings.n_jobs=2; SEED=20261030; np.random.seed(SEED)
O='/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna'; W=O+'/work'
meta=pd.read_csv('/REPLACE_WITH_PRIVATE_WORKDIR/data/scrna/GSE249874/GSE249874_sample_names_and_accessions.txt.gz',sep='\t')
meta.columns=['library','gsm','title']
def parse(t):
    stage='GC' if 'adenocarcinoma' in t else ('GS' if 'gatritis' in t or 'gastritis' in t else 'IM')
    hp='HP+' if 'positive' in t else 'HP-'; return stage,hp
meta[['stage','hp']]=meta.title.apply(lambda t: pd.Series(parse(t)))
meta['donor']=meta.stage+'_'+meta.hp.str.replace('+','P').str.replace('-','N')+'_'+meta.groupby(['stage','hp']).cumcount().add(1).astype(str)
meta[['library','gsm','stage','hp','donor']].to_csv(O+'/sample_design.csv',index=False)
print(meta[['library','gsm','stage','hp','donor']].to_string(),flush=True)
assert (meta.groupby(['stage','hp']).size()==3).all()
def mad_out(x,n=3,side='both'):
    m=np.median(x); d=np.median(np.abs(x-m))*1.4826
    lo=x<m-n*d; hi=x>m+n*d
    return hi if side=='upper' else (lo|hi)
qcrows=[]; parts=[]
for s in range(1,19):
    f=f'{W}/qc_s{s}.h5ad'
    if not os.path.exists(f):
        A=ad.read_h5ad(f'{W}/raw_s{s}.h5ad'); n0=A.n_obs
        A.var['mt']=A.var.sym.str.startswith('MT-').values if 'sym' in A.var else A.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(A,qc_vars=['mt'],inplace=True,percent_top=None,log1p=True)
        bad=mad_out(A.obs.log1p_n_genes_by_counts.values)|mad_out(A.obs.log1p_total_counts.values)|mad_out(A.obs.pct_counts_mt.values,side='upper')|(A.obs.pct_counts_mt.values>25)
        A=A[~bad].copy(); n1=A.n_obs
        sc.pp.scrublet(A,random_state=SEED); n_dbl=int(A.obs.predicted_doublet.sum()); A=A[~A.obs.predicted_doublet].copy()
        n2=A.n_obs
        if n2>4000: A=A[np.sort(np.random.RandomState(SEED+s).choice(n2,4000,replace=False))].copy()
        r=meta.iloc[s-1]
        for k in ['library','gsm','stage','hp','donor']: A.obs[k]=r[k]
        A.write_h5ad(f,compression='gzip')
        json.dump(dict(sample=f'sample{s}',called=n0,after_mad_mito=n1,doublets=n_dbl,after_scrublet=n2,final=A.n_obs),open(f'{W}/qc_s{s}.json','w'))
    qcrows.append(json.load(open(f'{W}/qc_s{s}.json')))
qc=pd.DataFrame(qcrows); qc.to_csv(O+'/qc_summary.csv',index=False); print(qc.to_string(),flush=True)
if not os.path.exists(W+'/integrated.h5ad'):
    A=ad.concat([ad.read_h5ad(f'{W}/qc_s{s}.h5ad') for s in range(1,19)],join='inner',merge='first')
    A.layers['counts']=A.X.copy()
    sc.pp.highly_variable_genes(A,n_top_genes=2000,flavor='seurat_v3',batch_key='sample',layer='counts')
    sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
    B=A[:,A.var.highly_variable].copy(); sc.pp.scale(B,max_value=10); sc.tl.pca(B,n_comps=30,random_state=SEED)
    ho=harmonypy.run_harmony(B.obsm['X_pca'],A.obs,'sample',random_state=SEED)
    Z=np.asarray(ho.Z_corr); Z=Z.T if Z.shape[0]!=A.n_obs else Z
    A.obsm['X_pca']=B.obsm['X_pca']; A.obsm['X_harmony']=Z; del B
    sc.pp.neighbors(A,n_neighbors=15,use_rep='X_harmony',random_state=SEED)
    sc.tl.leiden(A,resolution=1.0,random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
    sc.tl.umap(A,random_state=SEED)
    A.write_h5ad(W+'/integrated.h5ad',compression='gzip')
A=ad.read_h5ad(W+'/integrated.h5ad')
sym=pd.Series(A.var.sym.values if 'sym' in A.var else A.var_names, index=A.var_names)
MAJ={'Epithelial':['EPCAM','KRT8','KRT18','KRT19'],'T':['CD3D','CD3E','CD2'],'B':['MS4A1','CD79A','CD19'],
 'Plasma':['MZB1','JCHAIN','IGHG1','SDC1'],'Myeloid':['CD68','LYZ','CD14','C1QA','CSF1R','FCER1G'],'Mast':['TPSAB1','CPA3','KIT'],
 'Fibroblast':['DCN','LUM','COL1A1','PDGFRA'],'Endothelial':['PECAM1','VWF','CDH5'],'SMC_pericyte':['ACTA2','MYH11','DES','RGS5'],'Endocrine':['CHGA','CHGB']}
EPI={'Pit':['MUC5AC','TFF1','GKN1','GKN2'],'Neck':['MUC6','TFF2'],'Chief':['PGC','PGA3','LIPF'],'Parietal':['ATP4A','ATP4B','CBLIF'],'IM':['CDX2','MUC2','TFF3','REG4','FABP1','CDH17']}
def g2id(gs): return [i for i in sym.index[sym.isin(gs)]]
def assign(A,sets,key):
    sc_=pd.DataFrame(index=A.obs_names)
    X=A.X.tocsc() if sp.issparse(A.X) else A.X
    # per-gene z across cells, set score = mean z  (computed per cluster means to save memory)
    cl=A.obs[key].astype(str); res={}
    for k,gs in sets.items():
        ids=g2id(gs); idx=[A.var_names.get_loc(i) for i in ids]
        M=np.asarray(X[:,idx].todense()) if sp.issparse(X) else X[:,idx]
        Zs=(M-M.mean(0))/(M.std(0)+1e-9); res[k]=pd.Series(Zs.mean(1),index=A.obs_names)
    S=pd.DataFrame(res); cm=S.groupby(cl.values).mean(); return cm, cm.idxmax(1)
cm,lab=assign(A,MAJ,'leiden'); cm.to_csv(O+'/major_cluster_scores.csv')
A.obs['major']=A.obs.leiden.astype(str).map(lab)
comp=pd.crosstab(A.obs.leiden,A.obs.stage,normalize='index'); comp.to_csv(O+'/cluster_stage_composition.csv')
# epithelial subclustering
E=A[A.obs.major=='Epithelial'].copy()
Ec=ad.AnnData(E.layers['counts'],obs=E.obs[['sample','stage','hp','donor']].copy(),var=E.var[[]].copy())
sc.pp.highly_variable_genes(Ec,n_top_genes=2000,flavor='seurat_v3')  # deviation: no batch_key (loess fails for samples with few epithelial cells); Harmony corrects batch
H=Ec.var.highly_variable.values; del Ec
B=E[:,H].copy(); sc.pp.scale(B,max_value=10); sc.tl.pca(B,n_comps=30,random_state=SEED)
ho=harmonypy.run_harmony(B.obsm['X_pca'],E.obs,'sample',random_state=SEED); Z=np.asarray(ho.Z_corr); Z=Z.T if Z.shape[0]!=E.n_obs else Z
E.obsm['X_harmony']=Z; del B
sc.pp.neighbors(E,n_neighbors=15,use_rep='X_harmony',random_state=SEED)
sc.tl.leiden(E,resolution=1.0,random_state=SEED,flavor='igraph',n_iterations=2,directed=False,key_added='epi_leiden')
sc.tl.umap(E,random_state=SEED)
ecm,elab=assign(E,EPI,'epi_leiden'); gcfrac=E.obs.groupby('epi_leiden').stage.apply(lambda x:(x=='GC').mean())
ecm['frac_GC']=gcfrac; ecm['n']=E.obs.epi_leiden.value_counts()
elab=elab.where(gcfrac.reindex(elab.index)<0.7,'Tumor'); ecm['label']=elab; ecm.to_csv(O+'/epi_cluster_scores.csv')
E.obs['epi_type']=E.obs.epi_leiden.astype(str).map(elab)
A.obs['celltype']=A.obs.major.astype(str); A.obs.loc[E.obs_names,'celltype']='Epi_'+E.obs.epi_type.astype(str)
E.write_h5ad(W+'/epithelial.h5ad',compression='gzip'); A.write_h5ad(W+'/annotated.h5ad',compression='gzip')
ct=pd.crosstab(A.obs.donor,A.obs.celltype); ct.to_csv(O+'/celltype_by_donor.csv'); print(ct.to_string())
print(cm.round(2).to_string()); print(ecm.round(2).to_string())
# consistency check (not used for decisions): IM/CDX2 & parietal by stage
chk=pd.crosstab(A.obs.stage,A.obs.celltype,normalize='index').round(3); chk.to_csv(O+'/celltype_frac_by_stage.csv'); print(chk.to_string())
print('STEP2_DONE')
