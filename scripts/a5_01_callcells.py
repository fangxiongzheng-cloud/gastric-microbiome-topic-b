"""A5 step 1: stream GEO raw (unfiltered) 10x matrix of GSE249874, call cells, write per-sample h5ad.
Pre-declared (before looking at any expression): barcode suffix -k = GEO library 'samplek' (aggr order);
cell call per sample = UMI >= max(500, OrdMag) where OrdMag = 0.1 * 99th pct UMI of top 3000 barcodes
(CellRanger v2 algorithm, expected 3000 cells) AND detected genes >= 200. Memory-safe: chunked streaming.
Resumable: each stage skipped if its output exists."""
import gzip, os, sys, json, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
D='/REPLACE_WITH_PRIVATE_WORKDIR/data/scrna/GSE249874'; O='/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna/work'
os.makedirs(O,exist_ok=True); os.makedirs(O+'/tmp',exist_ok=True)
MTX=D+'/GSE249874_raw_feature_matrix.mtx.gz'; CH=20_000_000
feat=pd.read_csv(D+'/GSE249874_raw_feature_features.tsv.gz',sep='\t',header=None,names=['id','sym','type'])
mito=feat.sym.str.startswith('MT-').values
def header():
    with gzip.open(MTX,'rt') as f:
        n=0
        for l in f:
            n+=1
            if not l.startswith('%'): return n, [int(x) for x in l.split()]
skip,(nr,nc,nnz)=header(); print('mtx',nr,nc,nnz,flush=True); assert nr==len(feat)
def chunks():
    return pd.read_csv(MTX,sep=' ',header=None,skiprows=skip,chunksize=CH,dtype=np.int32,engine='c')
# pass 1
p1=O+'/pass1_totals.npz'
if not os.path.exists(p1):
    # memory: accumulate only over the column range touched by each chunk (mtx is column-sorted in CellRanger output)
    umi=np.zeros(nc,np.int32); ng=np.zeros(nc,np.int32); mt=np.zeros(nc,np.int32); done=0
    for c in chunks():
        r=c[0].values-1; j=c[1].values-1; v=c[2].values; lo=int(j.min()); hi=int(j.max())+1; jj=j-lo
        umi[lo:hi]+=np.bincount(jj,weights=v,minlength=hi-lo).astype(np.int32); ng[lo:hi]+=np.bincount(jj,minlength=hi-lo).astype(np.int32)
        m=mito[r]; mt[lo:hi]+=np.bincount(jj[m],weights=v[m],minlength=hi-lo).astype(np.int32)
        done+=len(c); print('pass1',done,f'{done/nnz:.1%}',flush=True)
    np.savez(p1,umi=umi,ng=ng,mt=mt)
t=np.load(p1); umi,ng,mt=t['umi'],t['ng'],t['mt']
# barcodes file verified (awk) to be 18 contiguous blocks of 6,794,880 (10x v3 whitelist), suffix -1..-18 in order
BLK=6794880; assert nc==18*BLK
samp=(np.arange(nc)//BLK+1).astype(np.int8)
called=np.zeros(nc,bool); rows=[]
for s in range(1,19):
    ix=np.where(samp==s)[0]; u=umi[ix]; top=np.sort(u)[::-1][:3000]
    om=0.1*np.percentile(top,99); thr=max(500,om)
    keep=(u>=thr)&(ng[ix]>=200); called[ix[keep]]=True
    knee=np.sort(u)[::-1]
    rows.append(dict(sample=f'sample{s}',ordmag_thr=round(om,1),thr=round(thr,1),n_called=int(keep.sum()),
        n_umi500=int((u>=500).sum()),median_umi_called=float(np.median(u[keep])) if keep.any() else 0,
        frac_reads_in_called=float(u[keep].sum()/max(u.sum(),1))))
cc=pd.DataFrame(rows); cc.to_csv(O+'/cellcall_summary.csv',index=False); print(cc.to_string(),flush=True)
# knee data for figure (rank vs umi, subsampled)
kn=[]
for s in range(1,19):
    u=np.sort(umi[samp==s])[::-1]; u=u[u>0]; rk=np.unique(np.geomspace(1,len(u),400).astype(int))-1
    kn.append(pd.DataFrame(dict(sample=f'sample{s}',rank=rk+1,umi=u[rk])))
pd.concat(kn).to_csv(O+'/knee_points.csv',index=False)
# pass 2: write entries of called cells to per-sample temp files
newidx=-np.ones(nc,np.int32)
# names of called barcodes: stream file once
cidx=np.where(called)[0]; names={}
if not os.path.exists(O+'/called_barcodes.csv'):
    want=set(cidx.tolist()); out=[]
    with gzip.open(D+'/GSE249874_raw_feature_barcodes.tsv.gz','rt') as f:
        for i,l in enumerate(f):
            if i in want: out.append((i,l.strip()))
    pd.DataFrame(out,columns=['col','barcode']).to_csv(O+'/called_barcodes.csv',index=False)
cb=pd.read_csv(O+'/called_barcodes.csv').set_index('col').barcode
for s in range(1,19):
    ix=np.where(called&(samp==s))[0]; newidx[ix]=np.arange(len(ix))
if not os.path.exists(O+'/tmp/pass2.done'):
    fh={s:open(f'{O}/tmp/s{s}.bin','wb') for s in range(1,19)}; done=0
    for c in chunks():
        r=c[0].values-1; j=c[1].values-1; v=c[2].values; k=called[j]
        r,j,v=r[k],j[k],v[k]; ss=samp[j]
        if len(j)==0:
            done+=len(c); continue
        for s in np.unique(ss):
            m=ss==s; np.stack([r[m],newidx[j[m]].astype(np.int32),v[m]],1).astype(np.int32).tofile(fh[s])
        done+=len(c); print('pass2',done,f'{done/nnz:.1%}',flush=True)
    for f in fh.values(): f.close()
    open(O+'/tmp/pass2.done','w').close()
meta=pd.read_csv(D+'/GSE249874_sample_names_and_accessions.txt.gz',sep='\t')
meta.columns=['library','gsm','title']
for s in range(1,19):
    out=f'{O}/raw_s{s}.h5ad'
    if os.path.exists(out): continue
    ix=np.where(called&(samp==s))[0]; a=np.fromfile(f'{O}/tmp/s{s}.bin',np.int32).reshape(-1,3)
    X=sp.csr_matrix((a[:,2].astype(np.float32),(a[:,1],a[:,0])),shape=(len(ix),nr))
    obs=pd.DataFrame(index=[f's{s}_'+b.rsplit('-',1)[0] for b in cb.loc[ix].values]); obs['sample']=f'sample{s}'
    A=ad.AnnData(X=X,obs=obs,var=feat.set_index('id')); A.var_names_make_unique(); A.write_h5ad(out,compression='gzip')
    print('wrote',out,A.shape,flush=True); del a,X,A
print('STEP1_DONE')
