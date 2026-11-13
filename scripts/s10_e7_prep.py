# Step10 E7 prep: export sc reference (sparse, symbol-collapsed) + Korea bulk gene counts. See freeze_log 2026-11-11 10:45.
import anndata as ad,numpy as np,pandas as pd,scipy.sparse as sp,gzip,json
B='/REPLACE_WITH_PRIVATE_WORKDIR/';O=B+'results/step10/E7/input/'
a=ad.read_h5ad(B+'results/a5_scrna/work/annotated.h5ad',backed='r')
X=a.layers['counts'][:] if not sp.issparse(a.layers['counts']) else a.layers['counts']
X=sp.csr_matrix(X) if not sp.issparse(X) else X
X=X.tocsc().astype(np.float64)
sym=a.var['sym'].astype(str).values
u,inv=np.unique(sym,return_inverse=True)
G=sp.csr_matrix((np.ones(len(sym)),(np.arange(len(sym)),inv)),shape=(len(sym),len(u)))
Xs=(X@G).tocsc(); Xs.sort_indices()
print('cells x genes',Xs.shape,'nnz',Xs.nnz,'dup symbols collapsed',len(sym)-len(u))
Xs.indptr.astype(np.int64).tofile(O+'sc_indptr.bin'); Xs.indices.astype(np.int32).tofile(O+'sc_indices.bin'); Xs.data.astype(np.float64).tofile(O+'sc_data.bin')
pd.Series(u).to_csv(O+'sc_genes.txt',index=False,header=False)
ob=a.obs[['major','celltype','donor','stage','hp']].copy(); ob['state']=ob.celltype.astype(str).replace({'Epi_Chief':'Epi_Neck_Chief','Epi_Neck':'Epi_Neck_Chief'})
ob.to_csv(O+'sc_obs.csv'); json.dump(dict(n_cells=Xs.shape[0],n_genes=Xs.shape[1],nnz=int(Xs.nnz)),open(O+'sc_dims.json','w'))
print(ob.state.value_counts())
# bulk Korea
t2g={}
with gzip.open('/REPLACE_WITH_REFERENCE_ROOT/gencode/gencode.v46.transcripts.fa.gz','rt') as f:
    for l in f:
        if l[0]=='>': x=l[1:].split('|'); t2g[x[0]]=x[5]
L=pd.read_csv(B+'data/external/rna_full352.tsv',sep='\t'); ex={}
for _,r in L.iterrows():
    q=pd.read_csv(B+f'results/korea_rna/per_run/{r.run_accession}/salmon/quant.sf.gz',sep='\t',usecols=['Name','NumReads'])
    q['g']=q.Name.map(t2g); ex[r.library_name]=q.groupby('g').NumReads.sum()
E=pd.DataFrame(ex).fillna(0).round().astype(np.int64)
E.to_csv(O+'korea_bulk_counts.csv.gz'); print('bulk',E.shape)
