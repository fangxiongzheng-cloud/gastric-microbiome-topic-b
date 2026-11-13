# v0.6 §4 sensitivity (non-confirmatory; estimation only) + v0.5 H-M sensitivities 1-4 + A.2 attenuation.
# Implementation choices fixed before looking at any adjusted result (2026-11-10 ~11:10 CST):
#  * Cell-composition covariates, approach (a): fixed marker scores (mean z of log2(TPM+1) across all cohort samples), gene lists
#    as locked in addendum v0.6a (epithelial EPCAM,KRT8,KRT18,CDH1; stromal COL1A1,COL1A2,COL3A1,DCN,LUM; immune PTPRC).
#    approach (b): ESTIMATE stromal + immune scores (tidyestimate 1.1.1, RNA-seq mode; ESTIMATE 'purity' formula is
#    Affymetrix-calibrated so it is NOT used), standardised within cohort.
#  * Gradient: mixed model HP ~ tum + logdepth (+covariates) + (1|patient), complete pairs. Korea HP = log1p(RPM) (A1 metric),
#    depth = log10 salmon-mapped; TCGA HP = log10(RPM+1) (H-M metric), depth = log10 homo_reads.
#  * H-M: partial Spearman (rank-residualised) of dHP vs dM_k given d-covariates (T-N); DL meta of Fisher z with se=1/sqrt(n-3-k).
#  * Patient-level bootstrap 1000 (seed 20261030) for §4; A.2 uses 2000 resamples seed 20261101 (v0.5).
import pandas as pd,numpy as np,json,warnings
from scipy.stats import spearmanr,rankdata,norm,pearsonr
import statsmodels.formula.api as smf, statsmodels.api as sm
from multiprocessing import get_context
Pool=lambda n: get_context("fork").Pool(n)
warnings.filterwarnings('ignore')
B='/REPLACE_WITH_PRIVATE_WORKDIR/';CORE=B+'results/core/';O=B+'results/sensitivity/'
IM=['CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4','SI','FABP1']
GD=['GKN1','GKN2','TFF2','PGC','PGA3','LIPF','CHIA','MUC6','CLDN18','GAST']
M2=['ATP4A','ATP4B','CBLIF','KCNE2','KCNJ15','KCNJ16','SLC26A7','CCKBR','ESRRG','VSIG1']
M3=['MUC5AC','TFF1','TFF2','MUC6','MUC1','FUT2','FUT3','B3GNT6']
EPI=['EPCAM','KRT8','KRT18','CDH1'];STR=['COL1A1','COL1A2','COL3A1','DCN','LUM'];IMM=['PTPRC']
def zm(X,g):
    g=[x for x in g if x in X.columns];Z=(X[g]-X[g].mean())/X[g].std(ddof=1);return Z.mean(axis=1)
def scores(X):
    return pd.DataFrame({'M1':zm(X,IM)-zm(X,GD),'M2':zm(X,M2),'M3':zm(X,M3),'epi':zm(X,EPI),'str':zm(X,STR),'imm':zm(X,IMM)})
# ---------- Korea
K=pd.read_csv(CORE+'korea_rna_per_sample.csv').set_index('lib')
need=sorted(set(IM+GD+M2+M3+EPI+STR+IMM+['ATP4A','MUC5AC','CDX2']))
E=pd.read_csv(CORE+'korea_gene_log2tpm.csv.gz',index_col=0);E=E.loc[[g for g in need if g in E.index]].T
SK=scores(E); es=pd.read_csv(O+'korea_estimate.csv').set_index('sample')
SK['est_str']=(es.stromal-es.stromal.mean())/es.stromal.std();SK['est_imm']=(es.immune-es.immune.mean())/es.immune.std()
w=K.pivot_table(index='pid',columns='type',values='hp_log10rpm1').dropna();pk=list(w.index)
DK=pd.DataFrame({'pid':K.pid,'tum':(K.type=='T').astype(int),'hp':K.hp_log1p_rpm,'hp10':K.hp_log10rpm1,'logdepth':np.log10(K.salmon_mapped),'k2':K.hp_k2_confirmed}).join(SK).join(E[['ATP4A','MUC5AC','CDX2']])
DK=DK[DK.pid.isin(pk)].copy()
# ---------- TCGA
s=pd.read_csv(B+'data/external/spectrum2023/stad_helicobacter_per_sample.csv',index_col=0)
s['hp10']=np.log10(s.hpylori_sp_raw/s.homo_reads*1e6+1);s['p']=s.index.str[:12];s['ty']=s.index.str[13:15]
hm=pd.read_csv(CORE+'tcga_stad_hm_genes_tpm.csv');mk=pd.read_csv(O+'tcga_purity_markers_tpm.csv')
e=hm.merge(mk[['barcode','ajcc_pathologic_stage','paper_Lauren.Class']+EPI+STR+IMM],on='barcode');e['k']=e.barcode.str[:16]
e=e.drop_duplicates('k').set_index('k')
gc=[c for c in e.columns if c in need];XT=np.log2(e[gc].astype(float)+1);ST=scores(XT)
et=pd.read_csv(O+'tcga_estimate.csv');et=et[et['sample']!='gene'];et['k']=et['sample'].str.replace('.','-').str[:16];et=et.drop_duplicates('k').set_index('k')
ST['est_str']=(et.stromal-et.stromal.mean())/et.stromal.std();ST['est_imm']=(et.immune-et.immune.mean())/et.immune.std()
Nn=s[s.ty=='11'].groupby('p').head(1);Tt=s[s.ty=='01'].groupby('p').head(1)
pt=sorted(set(Nn.p)&set(Tt.p)&set(e.patient))
ids=[*Nn[Nn.p.isin(pt)].index,*Tt[Tt.p.isin(pt)].index]
DT=pd.DataFrame({'pid':s.loc[ids,'p'],'tum':(s.loc[ids,'ty']=='01').astype(int),'hp':s.loc[ids,'hp10'],'hp10':s.loc[ids,'hp10'],'logdepth':np.log10(s.loc[ids,'homo_reads']),'k2':s.loc[ids,'hpylori_sp_raw']},index=ids).join(ST).join(XT[['ATP4A','MUC5AC','CDX2']])
DT['stage']=DT.pid.map(e.drop_duplicates('patient').set_index('patient').ajcc_pathologic_stage);DT['lauren']=DT.pid.map(e.drop_duplicates('patient').set_index('patient')['paper_Lauren.Class'])
DK['pid']=DK.pid.astype(str);DT['pid']=DT.pid.astype(str)
COH={'korea':DK,'tcga':DT}
COV={'unadj':[],'markers':['epi','str','imm'],'estimate':['est_str','est_imm']}
def fit_tum(D,f):
    for kw in [{},{'method':'powell'}]:
        try:return float(smf.mixedlm(f,D,groups=D.pid).fit(reml=True,**kw).params['tum'])
        except Exception:pass
    return np.nan
def boot_rs(D,seed,nb):
    rng=np.random.default_rng(seed);P=D.pid.unique();g={p:D[D.pid==p] for p in P};out=[]
    for b in range(nb):
        pick=rng.choice(P,len(P),replace=True);out.append(pd.concat([g[p].assign(pid=f'{p}_{i}') for i,p in enumerate(pick)]))
    return out
def _grad(args):
    D,forms=args;return [fit_tum(D,f) for f in forms]
def grad(D,dep=True,nb=1000,seed=20261030):
    base='hp~tum'+('+logdepth' if dep else '');forms={k:base+''.join('+'+c for c in v) for k,v in COV.items()}
    pt={k:fit_tum(D,f) for k,f in forms.items()}
    with Pool(6) as pl:bs=np.array(pl.map(_grad,[(d,list(forms.values())) for d in boot_rs(D,seed,nb)]))
    res={'n_pairs':int(D.pid.nunique()),'formulas':forms}
    for j,k in enumerate(forms):
        res[k]={'beta_tum':pt[k],'CI':list(np.nanpercentile(bs[:,j],[2.5,97.5])),'n_boot_failed':int(np.isnan(bs[:,j]).sum())}
        if k!='unadj':
            d=bs[:,j]-bs[:,0];a=(bs[:,0]-bs[:,j])/bs[:,0]
            res[k]['delta_vs_unadj']=pt[k]-pt['unadj'];res[k]['delta_CI']=list(np.nanpercentile(d,[2.5,97.5]))
            res[k]['prop_attenuation']=(pt['unadj']-pt[k])/pt['unadj'];res[k]['prop_att_CI']=list(np.nanpercentile(a,[2.5,97.5]))
    return res
def deltas(D,cols):
    T=D[D.tum==1].set_index('pid');N=D[D.tum==0].set_index('pid');p=T.index.intersection(N.index)
    return pd.DataFrame({c:T.loc[p,c].values-N.loc[p,c].values for c in cols},index=p)
def pspear(x,y,Z):
    rx,ry=rankdata(x),rankdata(y)
    if Z is None or Z.shape[1]==0:return pearsonr(rx,ry)[0]
    Zr=np.column_stack([np.ones(len(x))]+[rankdata(z) for z in Z.T])
    ex=rx-Zr@np.linalg.lstsq(Zr,rx,rcond=None)[0];ey=ry-Zr@np.linalg.lstsq(Zr,ry,rcond=None)[0];return pearsonr(ex,ey)[0]
def dl(zs,ses):
    zs=np.array(zs);v=np.array(ses)**2;w=1/v;fe=(w*zs).sum()/w.sum();Q=(w*(zs-fe)**2).sum();k=len(zs)
    t2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum()));ws=1/(v+t2);mu=(ws*zs).sum()/ws.sum();se=np.sqrt(1/ws.sum())
    return dict(rho=float(np.tanh(mu)),CI=[float(np.tanh(mu-1.96*se)),float(np.tanh(mu+1.96*se))],P=float(2*norm.sf(abs(mu/se))),tau2=float(t2))
def hm_adj(D,nb=1000,seed=20261030):
    Dd=deltas(D,['hp10','M1','M2','M3','epi','str','imm','est_str','est_imm']);rng=np.random.default_rng(seed);n=len(Dd);out={'n_pairs':n}
    for m in ['M1','M2','M3']:
        out[m]={}
        for k,cv in COV.items():
            Z=Dd[cv].values if cv else None;r=pspear(Dd.hp10.values,Dd[m].values,Z)
            bs=[];bd=[]
            for b in range(nb):
                i=rng.integers(0,n,n);rb=pspear(Dd.hp10.values[i],Dd[m].values[i],Z[i] if Z is not None else None);bs.append(rb)
                if cv:bd.append(rb-pspear(Dd.hp10.values[i],Dd[m].values[i],None))
            o={'rho':float(r),'CI':list(np.nanpercentile(bs,[2.5,97.5])),'z':float(np.arctanh(r)),'se':float(1/np.sqrt(n-3-len(cv)))}
            if cv:o['delta_rho']=float(r-out[m]['unadj']['rho']);o['delta_CI']=list(np.nanpercentile(bd,[2.5,97.5]))
            out[m][k]=o
    return out,Dd
R={'notes':'Non-confirmatory sensitivity (v0.6 §4; v0.5 H-M sens 1-4; A.2). Estimation only; no multiplicity correction; no causal/mediation claims.'}
# §4.1 gradient + H-M adjusted
R['gradient']={c:grad(D) for c,D in COH.items()}
R['hm_adjusted']={};DD={}
for c,D in COH.items():R['hm_adjusted'][c],DD[c]=hm_adj(D)
R['hm_adjusted']['meta']={m:{k:dl([R['hm_adjusted'][c][m][k]['z'] for c in COH],[R['hm_adjusted'][c][m][k]['se'] for c in COH]) for k in COV} for m in ['M1','M2','M3']}
# v0.5 H-M sensitivities
S={}
for c,D in COH.items():
    Dd=DD[c];sc={}
    # 1. multivariable: dHP ~ dscore + mean log depth (+stage, Lauren in TCGA)
    T=D[D.tum==1].set_index('pid');N=D[D.tum==0].set_index('pid');p=Dd.index
    X=Dd.copy();X['mdep']=(T.loc[p,'logdepth'].values+N.loc[p,'logdepth'].values)/2
    if c=='tcga':
        st=T.loc[p,'stage'].astype(str);X['late']=np.where(st.str.contains('III|IV'),1,np.where(st.str.contains('Stage I'),0,np.nan))
        X['lauren']=T.loc[p,'lauren'].values
    for m in ['M1','M2','M3']:
        f=f'hp10~{m}+mdep'+('+late+C(lauren)' if c=='tcga' else '')
        try:
            mm=smf.ols(f,X.dropna(subset=['late','lauren'] if c=='tcga' else None)).fit()
            sc[f'sens1_{m}']=dict(beta=float(mm.params[m]),CI=[float(x) for x in mm.conf_int().loc[m]],P=float(mm.pvalues[m]),n=int(mm.nobs),formula=f)
        except Exception as ex:sc[f'sens1_{m}']={'error':str(ex)}
    if c=='korea':sc['sens1_note']='stage and Lauren not available for Korea PRJNA1174735: unadjustable (depth only).'
    # 2. single genes + leave-one-gene-out
    dg=deltas(D,['hp10','ATP4A','MUC5AC','CDX2'])
    for g in ['ATP4A','MUC5AC','CDX2']:
        r=spearmanr(dg.hp10,dg[g]);sc[f'sens2_gene_{g}']=dict(rho=float(r.statistic),P=float(r.pvalue),n=len(dg))
    X0=E if c=='korea' else XT
    for m,parts in {'M1':(IM,GD),'M2':(M2,None),'M3':(M3,None)}.items():
        rr={}
        for gl in [parts[0]]+([parts[1]] if parts[1] else []):
            for g in gl:
                if g not in X0.columns:continue
                if m=='M1':sco=zm(X0,[x for x in IM if x!=g])-zm(X0,[x for x in GD if x!=g])
                else:sco=zm(X0,[x for x in gl if x!=g])
                tmp=D[['pid','tum','hp10']].assign(sc=sco.reindex(D.index).values);d2=deltas(tmp,['hp10','sc'])
                rr[g]=float(spearmanr(d2.hp10,d2.sc).statistic)
        sc[f'sens2_LOGO_{m}']=dict(min=min(rr.values()),max=max(rr.values()),by_gene_dropped=rr)
    # 3. pairs with HP>=10 reads in normal
    keep=N.loc[p,'k2']>=10;d3=Dd[keep.values]
    for m in ['M1','M2','M3']:
        r=spearmanr(d3.hp10,d3[m]) if len(d3)>3 else None
        sc[f'sens3_{m}']=dict(rho=float(r.statistic),P=float(r.pvalue),n=len(d3)) if r is not None else {'n':len(d3)}
    # 4. unpaired within tissue
    for t,lab in [(0,'N'),(1,'T')]:
        Dt=D[D.tum==t]
        for m in ['M1','M2','M3']:
            r=spearmanr(Dt.hp10,Dt[m]);sc[f'sens4_{lab}_{m}']=dict(rho=float(r.statistic),P=float(r.pvalue),n=len(Dt))
    S[c]=sc
R['hm_v05_sensitivities']=S
# A.2 attenuation (v0.5): Model0 hp10~tum+(1|pid); Model1 + M1+M2+M3; singles. 2000 bootstrap seed 20261101.
A={}
for c,D in COH.items():
    D2=D.assign(hp=D.hp10);forms={'M0':'hp~tum','M1all':'hp~tum+M1+M2+M3','M1only':'hp~tum+M1','M2only':'hp~tum+M2','M3only':'hp~tum+M3'}
    pt={k:fit_tum(D2,f) for k,f in forms.items()}
    with Pool(6) as pl:bs=np.array(pl.map(_grad,[(d,list(forms.values())) for d in boot_rs(D2,20261101,2000)]))
    a={'beta':pt}
    for j,k in enumerate(forms):
        if k=='M0':continue
        att=(bs[:,0]-bs[:,j])/bs[:,0];a[k]=dict(prop_attenuation=(pt['M0']-pt[k])/pt['M0'],CI=list(np.nanpercentile(att,[2.5,97.5])),boot_se=float(np.nanstd(att,ddof=1)))
    A[c]=a
w_=[1/A[c]['M1all']['boot_se']**2 for c in COH];A['pooled_M1all_IV']=float(sum(wi*A[c]['M1all']['prop_attenuation'] for wi,c in zip(w_,COH))/sum(w_))
A['pooled_M1all_IV_CI']=[A['pooled_M1all_IV']-1.96/np.sqrt(sum(w_)),A['pooled_M1all_IV']+1.96/np.sqrt(sum(w_))]
R['A2_attenuation']=A
json.dump(R,open(O+'sensitivity_s4_results.json','w'),indent=1,default=float);print('DONE')
