"""A5 step 1b (deviation: EmptyDrops cell calling). Dump per-sample entries for barcodes with UMI>=10 (for EmptyDrops ambient + testing)."""
import os, numpy as np, pandas as pd
D='/REPLACE_WITH_PRIVATE_WORKDIR/data/scrna/GSE249874'; W='/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna/work'; T=W+'/ed'; os.makedirs(T,exist_ok=True)
MTX=D+'/GSE249874_raw_feature_matrix.mtx.gz'; BLK=6794880
t=np.load(W+'/pass1_totals.npz'); umi=t['umi']; keep=umi>=10
if not os.path.exists(T+'/dump.done'):
    fh={s:open(f'{T}/s{s}.bin','wb') for s in range(1,19)}; done=0
    for c in pd.read_csv(MTX,sep=' ',header=None,skiprows=3,chunksize=20_000_000,dtype=np.int32,engine='c'):
        r=c[0].values-1; j=c[1].values-1; v=c[2].values; k=keep[j]; r,j,v=r[k],j[k],v[k]; ss=j//BLK+1
        for s in np.unique(ss):
            m=ss==s; np.stack([r[m],j[m]-(s-1)*BLK,v[m]],1).astype(np.int32).tofile(fh[s])
        done+=len(c); print('dump',done,flush=True)
    for f in fh.values(): f.close()
    open(T+'/dump.done','w').close()
cc=pd.read_csv(W+'/cellcall_summary.csv'); cc[['sample','ordmag_thr']].to_csv(T+'/ordmag.csv',index=False)
print('STEP1B_DONE')
