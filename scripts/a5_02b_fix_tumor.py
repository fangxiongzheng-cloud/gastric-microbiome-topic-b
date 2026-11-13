"""A5 step 2b (deviation, decided after seeing the annotation table but before any CCC/trajectory output):
the cluster-level rule '>=70% GC cells -> Tumor' also labelled gastritis/IM-donor cells as tumour.
Fix: within those clusters only cells from GC donors are 'Tumor'; other cells get the cluster's best normal epithelial lineage
(argmax of the pre-declared Pit/Neck/Chief/Parietal/IM scores). Writes annotated.h5ad/epithelial.h5ad in place (originals kept as *_v1.h5ad)."""
import os, shutil, pandas as pd, anndata as ad
O='/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna'; W=O+'/work'
for f in ['annotated','epithelial']:
    if not os.path.exists(f'{W}/{f}_v1.h5ad'): shutil.copy(f'{W}/{f}.h5ad',f'{W}/{f}_v1.h5ad')
E=ad.read_h5ad(W+'/epithelial_v1.h5ad'); A=ad.read_h5ad(W+'/annotated_v1.h5ad')
ecm=pd.read_csv(O+'/epi_cluster_scores.csv',index_col=0); ecm.index=ecm.index.astype(str)
normal=ecm[['Pit','Neck','Chief','Parietal','IM']].idxmax(1); tumclus=set(ecm.index[ecm.frac_GC>=0.7])
cl=E.obs.epi_leiden.astype(str); isgc=E.obs.stage.values=='GC'
lab=cl.map(normal).values.astype(object); m=cl.isin(tumclus).values & isgc; lab[m]='Tumor'
E.obs['epi_type']=lab; E.obs['epi_type_v1']=cl.map(ecm.label).values
A.obs['celltype']=A.obs.celltype.astype(str); A.obs.loc[E.obs_names,'celltype']='Epi_'+E.obs.epi_type.astype(str)
E.write_h5ad(W+'/epithelial.h5ad',compression='gzip'); A.write_h5ad(W+'/annotated.h5ad',compression='gzip')
ct=pd.crosstab(A.obs.donor,A.obs.celltype); ct.to_csv(O+'/celltype_by_donor.csv'); print(ct.to_string())
pd.crosstab(A.obs.stage,A.obs.celltype,normalize='index').round(3).to_csv(O+'/celltype_frac_by_stage.csv')
print('STEP2B_DONE')
