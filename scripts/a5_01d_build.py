"""A5 step 1d: build per-sample raw h5ad of EmptyDrops-called cells (replaces OrdMag-only raw_s*.h5ad, which are moved to work/stale_ordmag/)."""
import os, gzip, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
D='/REPLACE_WITH_PRIVATE_WORKDIR/data/scrna/GSE249874'; W='/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna/work'; T=W+'/ed'; BLK=6794880
feat=pd.read_csv(D+'/GSE249874_raw_feature_features.tsv.gz',sep='\t',header=None,names=['id','sym','type'])
calls={s:pd.read_csv(f'{T}/called_s{s}.csv') for s in range(1,19)}
glob=np.concatenate([calls[s].col.values+(s-1)*BLK for s in range(1,19)]); want=set(glob.tolist())
if not os.path.exists(T+'/called_barcodes.csv'):
    out=[]
    with gzip.open(D+'/GSE249874_raw_feature_barcodes.tsv.gz','rt') as f:
        for i,l in enumerate(f):
            if i in want: out.append((i,l.strip()))
    pd.DataFrame(out,columns=['gcol','barcode']).to_csv(T+'/called_barcodes.csv',index=False)
cb=pd.read_csv(T+'/called_barcodes.csv').set_index('gcol').barcode; rows=[]
for s in range(1,19):
    out=f'{W}/raw_s{s}.h5ad'; c=calls[s].col.values
    if not os.path.exists(out):
        a=np.fromfile(f'{T}/s{s}.bin',np.int32).reshape(-1,3); pos=pd.Series(np.arange(len(c)),index=c)
        k=np.isin(a[:,1],c); a=a[k]
        X=sp.csr_matrix((a[:,2].astype(np.float32),(pos.loc[a[:,1]].values,a[:,0])),shape=(len(c),len(feat)))
        obs=pd.DataFrame(index=[f's{s}_'+b.rsplit('-',1)[0] for b in cb.loc[c+(s-1)*BLK].values]); obs['sample']=f'sample{s}'
        A=ad.AnnData(X=X,obs=obs,var=feat.set_index('id')); A.write_h5ad(out,compression='gzip'); del a,X,A
    rows.append(dict(sample=f'sample{s}',n_called_emptydrops=len(c),median_umi=float(calls[s].umi.median())))
pd.DataFrame(rows).to_csv(W+'/../cellcall_emptydrops_summary.csv',index=False); print(pd.DataFrame(rows).to_string()); print('STEP1D_DONE')
