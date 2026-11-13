# Step10 E7a/E7b/E7c (v0.6b) + E9 (v0.6d) analyses on BayesPrism outputs. EXPLORATORY. Decisions: freeze_log (E7 09:45, E9/E11 entries).
import pandas as pd,numpy as np,json,warnings,os,sys
from scipy.stats import spearmanr,wilcoxon,norm
from scipy import stats
import statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor as vif
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
B='/REPLACE_WITH_PRIVATE_WORKDIR/'; I=os.environ.get('E7DIR',B+'results/step10/E7/'); O7=os.environ.get('E7OUT',I); O9=os.environ.get('E9OUT',B+'results/step10/E9/'); os.makedirs(O9,exist_ok=True)
SEED=20261030; NB=int(os.environ.get('NBOOT',2000))
IM=['CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4','SI','FABP1']; GD=['GKN1','GKN2','TFF2','PGC','PGA3','LIPF','CHIA','MUC6','CLDN18','GAST']
M2=['ATP4A','ATP4B','CBLIF','KCNE2','KCNJ15','KCNJ16','SLC26A7','CCKBR','ESRRG','VSIG1']; M3=['MUC5AC','TFF1','TFF2','MUC6','MUC1','FUT2','FUT3','B3GNT6']
E9G=['MUC5AC','TFF1','TFF2','MUC6','MUC1','FUT2','FUT3','B3GNT6','ST3GAL4','A4GNT','CEACAM1','CEACAM5','CEACAM6']; CIRC={'MUC5AC','TFF1','TFF2','MUC6','A4GNT'}
IMM=['T','B','Plasma','Myeloid','Mast']; STR=['Fibroblast','Endothelial','SMC_pericyte']
def zscore_mean(X,genes,name,log):  # core_hm.py frozen implementation
    pres=[g for g in genes if g in X.columns]; miss=[g for g in genes if g not in pres]; log[name]=dict(present=len(pres),missing=miss)
    if len(pres)/len(genes)<0.7: log[name]['FAILED_<70%']=True; return None
    Z=(X[pres]-X[pres].mean())/X[pres].std(ddof=1); return Z.mean(axis=1)
def mscores(X,log):
    im=zscore_mean(X,IM,'IM',log); gd=zscore_mean(X,GD,'GD',log); m2=zscore_mean(X,M2,'M2',log); m3=zscore_mean(X,M3,'M3',log)
    d={}
    if im is not None and gd is not None: d['M1']=im-gd
    if m2 is not None: d['M2']=m2
    if m3 is not None: d['M3']=m3
    return pd.DataFrame(d, index=X.index)
logit=lambda p:np.log(np.clip(p,1e-4,1-1e-4)/(1-np.clip(p,1e-4,1-1e-4)))
def epi_cpm(f):
    E=pd.read_csv(f,index_col=0); return np.log2(E.div(E.sum(1),axis=0)*1e6+1)
res={'label':'E7 EXPLORATORY (v0.6b); E9 EXPLORATORY (v0.6d)','deconvolution':'BayesPrism 2.2.3 (sparse/chunk-equivalent marker step; see freeze_log)'}
# ---------------- Korea inputs
K=pd.read_csv(B+'results/core/korea_rna_per_sample.csv'); K['pid']=K.pid.astype(str)
w=K.pivot_table(index='pid',columns='type',values='hp_log1p_rpm').dropna(); pids=list(w.index)
K=K[K.pid.isin(pids)].copy().set_index('lib'); K['logdepth']=np.log10(K.salmon_mapped); K['tum']=(K.type=='T').astype(int)
F=pd.read_csv(I+'korea_type_fraction_final.csv',index_col=0); S=pd.read_csv(I+'korea_state_fraction_first.csv',index_col=0)
K['epi']=F.Epithelial.reindex(K.index); K['imm']=F[IMM].sum(1).reindex(K.index); K['str']=F[STR].sum(1).reindex(K.index)
K['gastric_state']=(S['Epi_Pit']+S['Epi_Neck_Chief']).reindex(K.index)
for c in ['epi','imm','str','gastric_state']: K['l_'+c]=logit(K[c])
res['fractions_summary']={t:{'median_N':float(F.loc[[p+'_N' for p in pids],t].median()),'median_T':float(F.loc[[p+'_T' for p in pids],t].median())} for t in F.columns}
st={}
for c in S.columns:
    d=S.loc[[p+'_T' for p in pids],c].values-S.loc[[p+'_N' for p in pids],c].values; st[c]=dict(median_T_minus_N=float(np.median(d)),P_wilcoxon=float(wilcoxon(d).pvalue) if np.any(d!=0) else 1.0)
res['epithelial_state_paired_descriptive']=st
Kd=K.reset_index()
def mixed(df,f):
    for kw in [{},{'method':'powell'}]:
        try:
            m=smf.mixedlm(f,df,groups=df.pid).fit(reml=True,**kw)
            if np.all(np.isfinite(m.params)): return m
        except Exception: pass
    return None
rng=np.random.default_rng(SEED); BI=[rng.choice(pids,len(pids)) for _ in range(NB)]
grp={p:Kd[Kd.pid==p] for p in pids}
def bdf(ix): return pd.concat([grp[p].assign(pid=f'{p}_{j}') for j,p in enumerate(ix)])
forms={'base':'hp_log1p_rpm ~ tum + logdepth','adj_epi':'hp_log1p_rpm ~ tum + logdepth + l_epi','adj_epi_imm_str':'hp_log1p_rpm ~ tum + logdepth + l_epi + l_imm + l_str'}
e7a={}; fits={k:mixed(Kd,f) for k,f in forms.items()}
boots={k:[] for k in forms}
for ix in BI:
    d=bdf(ix)
    for k,f in forms.items(): m=mixed(d,f); boots[k].append(m.params['tum'] if m is not None else np.nan)
for k in forms:
    b=np.array(boots[k]); e7a[k]=dict(beta_tum=float(fits[k].params['tum']),P=float(fits[k].pvalues['tum']),CI95_boot=[float(np.nanpercentile(b,2.5)),float(np.nanpercentile(b,97.5))],n_boot_failed=int(np.isnan(b).sum()))
    if k!='base':
        dlt=np.array(boots[k])-np.array(boots['base']); e7a[k]['change_vs_base']=float(fits[k].params['tum']-fits['base'].params['tum']); e7a[k]['change_CI95_boot']=[float(np.nanpercentile(dlt,2.5)),float(np.nanpercentile(dlt,97.5))]
        e7a[k]['pct_attenuation']=float(100*(1-fits[k].params['tum']/fits['base'].params['tum']))
res['E7a']=e7a
# ---------------- E7b
hm=pd.read_csv(B+'results/core/hm_deltas_korea.csv'); hm['pid']=hm.pid.astype(str); hm=hm.set_index('pid').loc[pids]
D=hm.copy()
for c in ['epi','imm','str']: D['d_'+c]=K.loc[[p+'_T' for p in pids],c].values-K.loc[[p+'_N' for p in pids],c].values
z=lambda s:(s-s.mean())/s.std(ddof=1)
def huber(df,y,xs): return sm.RLM(z(df[y]),sm.add_constant(df[xs].apply(z)),M=sm.robust.norms.HuberT()).fit()
def pspear(df,x,y,cov):
    R=df[[x,y]+cov].rank(); C=sm.add_constant(R[cov]); rx=sm.OLS(R[x],C).fit().resid; ry=sm.OLS(R[y],C).fit().resid
    r=np.corrcoef(rx,ry)[0,1]; n=len(df); k=len(cov); t=r*np.sqrt((n-2-k)/(1-r**2)); return float(r),float(2*stats.t.sf(abs(t),n-2-k))
ep=epi_cpm(I+'korea_epithelial_exp.csv.gz'); lg={}; Me=mscores(ep,lg); res['E7b_model2_gene_presence']=lg
e7b={}; tests=[('E7a','adj_epi',e7a['adj_epi']['P']),('E7a','adj_epi_imm_str',e7a['adj_epi_imm_str']['P'])]
cov=['d_epi','d_imm','d_str']; rng2=np.random.default_rng(SEED); RI=[rng2.integers(0,len(D),len(D)) for _ in range(NB)]
for k in ['M1','M2','M3']:
    x='d'+k; u=huber(D,'dHP',[x]); m=huber(D,'dHP',[x]+cov)
    bb=[huber(D.iloc[ix],'dHP',[x]+cov).params[x] for ix in RI]
    pr,pp=pspear(D,'dHP',x,cov)
    Xv=sm.add_constant(D[[x]+cov].apply(z)).values
    r1=dict(std_beta_unadj=float(u.params[x]),std_beta_adj=float(m.params[x]),P=float(m.pvalues[x]),CI95_boot=[float(np.percentile(bb,2.5)),float(np.percentile(bb,97.5))],
            partial_spearman=pr,partial_spearman_P=pp,VIF={c:float(vif(Xv,i+1)) for i,c in enumerate([x]+cov)})
    if k in Me.columns:
        dm=Me.loc[[p+'_T' for p in pids],k].values-Me.loc[[p+'_N' for p in pids],k].values; r=spearmanr(D.dHP,dm); zz=np.arctanh(r.statistic); se=1/np.sqrt(len(dm)-3)
        r2=dict(rho=float(r.statistic),P=float(r.pvalue),CI95=[float(np.tanh(zz-1.96*se)),float(np.tanh(zz+1.96*se))],n=len(dm))
    else: r2=dict(status='score failed (<70% genes)',P=np.nan)
    e7b[k]={'model1':r1,'model2_epithelial_expression':r2}; tests+=[('E7b_'+k,'model1',r1['P']),('E7b_'+k,'model2',r2['P'])]
ok=[t for t in tests if np.isfinite(t[2])]; q=multipletests([t[2] for t in ok],method='fdr_bh')[1]
res['BH_family_E7']=[dict(item=a,test=b,P=float(p),q=float(qq)) for (a,b,p),qq in zip(ok,q)]
res['E7b']=e7b
# ---------------- E7c TCGA
try:
    td=pd.read_csv(B+'results/core/hm_deltas_tcga.csv').set_index('pid'); et=epi_cpm(I+'tcga_epithelial_exp.csv.gz'); lt={}; Mt=mscores(et,lt); e7c={'gene_presence':lt}
    for k in ['M1','M2','M3']:
        if k not in Mt.columns:
            e7c[k]=dict(status='score failed (<70% genes)',note='direction only'); continue
        dm=Mt.loc[[p+'_T' for p in td.index],k].values-Mt.loc[[p+'_N' for p in td.index],k].values; r=spearmanr(td.dHP,dm); zz=np.arctanh(r.statistic); se=1/np.sqrt(len(dm)-3)
        e7c[k]=dict(rho=float(r.statistic),CI95=[float(np.tanh(zz-1.96*se)),float(np.tanh(zz+1.96*se))],n=len(dm),note='direction only')
    res['E7c_TCGA']=e7c
except Exception as ex: res['E7c_TCGA']=dict(error=str(ex))
json.dump(res,open(O7+'E7_results.json','w'),indent=1,default=float)
D.assign(**{f'd{k}_epi':(Me.loc[[p+'_T' for p in pids],k].values-Me.loc[[p+'_N' for p in pids],k].values) for k in ['M1','M2','M3'] if k in Me.columns}).to_csv(O7+'E7b_korea_deltas.csv')
# ---------------- E9
E=pd.read_csv(B+'results/core/korea_gene_log2tpm.csv.gz',index_col=0).T
r9={}; t9=[]
for g in E9G:
    rr={}
    d=E.loc[[p+'_T' for p in pids],g].values-E.loc[[p+'_N' for p in pids],g].values
    bt=[np.median(d[ix]) for ix in RI]; rr['E9a']=dict(median_diff=float(np.median(d)),CI95_boot=[float(np.percentile(bt,2.5)),float(np.percentile(bt,97.5))],P_wilcoxon=float(wilcoxon(d).pvalue))
    Kd['y']=E.loc[Kd.lib,g].values
    m=mixed(Kd,'y ~ tum + l_epi + logdepth'); dl=(K.loc[[p+'_T' for p in pids],'l_epi'].values-K.loc[[p+'_N' for p in pids],'l_epi'].values)
    rr['E9b']=dict(beta_tum=float(m.params['tum']),CI95=[float(x) for x in m.conf_int().loc['tum']],P=float(m.pvalues['tum']),
        observed_mean_change=float(d.mean()),expected_change_from_epi=float(m.params['l_epi']*dl.mean()))
    ms=mixed(Kd,'y ~ tum + l_epi + l_gastric_state + logdepth')
    rr['E9b_sens_gastric_state']=dict(beta_tum=float(ms.params['tum']),P=float(ms.pvalues['tum']),circular_flag=g in CIRC) if ms is not None else 'fit failed'
    t9.append((g,'E9b',rr['E9b']['P']))
    if g in ep.columns:
        de=ep.loc[[p+'_T' for p in pids],g].values-ep.loc[[p+'_N' for p in pids],g].values
        rr['E9c']=dict(median_diff=float(np.median(de)),P=float(wilcoxon(de).pvalue) if np.any(de!=0) else 1.0)
        r=spearmanr(D.dHP if False else (w.loc[pids,'T']-w.loc[pids,'N']).values,de); rr['E9d']=dict(rho=float(r.statistic),P=float(r.pvalue))
        t9+=[(g,'E9c',rr['E9c']['P']),(g,'E9d',rr['E9d']['P'])]
    else:
        _mk=set(open(B+'results/step10/E7/marker_genes.txt').read().split())
        rr['E9c']=rr['E9d']=('not estimable (marker gene, but not retained in the Korean BayesPrism mixture; consistent with bulk outlier-gene filtering, 16 genes, E7 log/run_nc4.log)' if g in _mk else 'not estimable (not a deconvolution marker gene)')
    r9[g]=rr
q=multipletests([t[2] for t in t9],method='fdr_bh')[1]
for (g,k,_),qq in zip(t9,q): r9[g][k]['q_BH']=float(qq)
for g,rr in r9.items():
    ok=isinstance(rr['E9c'],dict) and rr['E9b']['q_BH']<0.10 and rr['E9c']['q_BH']<0.10 and np.sign(rr['E9b']['beta_tum'])==np.sign(rr['E9c']['median_diff'])
    rr['verdict']=('change within epithelium / beyond epithelial proportion ('+('loss' if rr['E9b']['beta_tum']<0 else 'gain')+')') if ok else 'indistinguishable from epithelial-proportion change'
try:
    Tt=pd.read_csv(B+'results/step10/E3/tcga_selected_log2tpm.csv',index_col=0).T; td=pd.read_csv(B+'results/core/hm_deltas_tcga.csv').set_index('pid')
    et=epi_cpm(I+'tcga_epithelial_exp.csv.gz'); tc={}
    for g in E9G:
        a=float(np.median(Tt.loc[[p+'_T' for p in td.index],g].values-Tt.loc[[p+'_N' for p in td.index],g].values)) if g in Tt.columns else None
        c=float(np.median(et.loc[[p+'_T' for p in td.index],g].values-et.loc[[p+'_N' for p in td.index],g].values)) if g in et.columns else None
        tc[g]=dict(E9a_median_diff=a,E9c_median_diff=c)
    r9['_TCGA_direction_only']=tc
except Exception as ex: r9['_TCGA_direction_only']=dict(error=str(ex))
json.dump({'label':'E9 EXPLORATORY (v0.6d); BH over available E9b/E9c/E9d tests','n_tests_BH':len(t9),'genes':r9},open(O9+'E9_results.json','w'),indent=1,default=float)
# figure: observed vs expected change
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(6,5)); gg=[g for g in E9G]
ob=[r9[g]['E9b']['observed_mean_change'] for g in gg]; ex=[r9[g]['E9b']['expected_change_from_epi'] for g in gg]
ax.scatter(ex,ob,c=['#D55E00' if 'within' in r9[g]['verdict'] else '#999999' for g in gg])
for g,a,b in zip(gg,ex,ob): ax.annotate(g+('*' if g in CIRC else ''),(a,b),fontsize=7)
lim=[min(ex+ob)-0.3,max(ex+ob)+0.3]; ax.plot(lim,lim,ls=':',c='k',lw=0.8); ax.axhline(0,c='lightgrey',lw=0.5); ax.axvline(0,c='lightgrey',lw=0.5)
ax.set_xlabel('expected T−N change from Δ epithelial proportion (log2 TPM)'); ax.set_ylabel('observed mean T−N change (log2 TPM)'); ax.set_title('E9 (exploratory) Korea 174 pairs',fontsize=9)
fig.tight_layout(); fig.savefig(O9+'E9_observed_vs_expected.pdf'); fig.savefig(O9+'E9_observed_vs_expected.png',dpi=200)
print(json.dumps({k:res[k] for k in ['E7a','BH_family_E7']},indent=1,default=float)[:4000]); print('E9 verdicts',{g:r9[g]['verdict'] for g in E9G})
