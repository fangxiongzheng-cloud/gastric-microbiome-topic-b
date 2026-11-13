# Step10: A4 registered primary full-set rerun (PREDECLARED_a4.md) + E10 pseudo-temporal ordering (v0.6d). Decisions: freeze_log 2026-11-11 09:42.
import pandas as pd,numpy as np,json,warnings
from scipy.stats import spearmanr
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
B='/REPLACE_WITH_PRIVATE_WORKDIR/';SEED=20261030;NPERM=10000
IM=['CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4']; OX=['ATP4A','ATP4B','CBLIF','PGC','GKN1','GKN2']
E=pd.read_csv(B+'results/core/korea_gene_log2tpm.csv.gz',index_col=0)
K=pd.read_csv(B+'results/core/korea_rna_per_sample.csv')
Ncols=[c for c in E.columns if c.endswith('_N')]
def score(cols,gs):
    X=E.loc[gs,cols]; z=X.sub(X.mean(1),axis=0).div(X.std(1,ddof=1),axis=0); return z.mean()
# ---------- A4 full
s=pd.read_csv(B+'results/korea16s/per_sample_silva.csv',index_col=0); s['lib']=s.pid.astype(str)+'_'+s.type_bs; s=s.set_index('lib')
N=[c for c in Ncols if c in s.index]
D=pd.DataFrame({'clr':s.loc[N,'clr'],'reads16s':s.loc[N,'reads'],'IM':score(N,IM),'OX':score(N,OX)})
R=E[N][(E[N]>1).mean(1)>=0.5].rank(); D['IM_rank']=R.loc[[g for g in IM if g in R.index]].mean()/len(R)
D['pid']=[c[:-2] for c in N]
T={c[:-2]:s.loc[c,'clr'] for c in E.columns if c.endswith('_T') and c in s.index}
D['dCLR']=D.pid.map(lambda p: D.loc[p+'_N','clr']-T[p] if p in T else np.nan)
def sp(a,b,d=None):
    d=(D if d is None else d)[[a,b]].dropna(); r=spearmanr(d[a],d[b]); return dict(x=a,y=b,n=len(d),rho=float(r.statistic),P=float(r.pvalue))
A=dict(label='A4 registered primary, FULL-SET rerun (PREDECLARED_a4.md; interim n=17 result kept separately)',primary=sp('clr','IM'),prespecified_direction='negative')
A['success']=bool(A['primary']['P']<0.05 and A['primary']['rho']<0)
A['secondary']=[sp('clr','OX'),sp('dCLR','IM'),sp('clr','IM_rank')]
q=multipletests([x['P'] for x in A['secondary']],method='fdr_bh')[1]
for x,qq in zip(A['secondary'],q): x['q_BH']=float(qq)
A['sens_16s_ge1000']=sp('clr','IM',D[D.reads16s>=1000])
D.to_csv(B+'results/step10/A4full/korea_normal_scores_full.csv')
json.dump(A,open(B+'results/step10/A4full/A4_full_results.json','w'),indent=1)
print(json.dumps(A,indent=1))
# ---------- E10 (Korea)
Kp=K.copy(); w=Kp.pivot_table(index='pid',columns='type',values='hp_log1p_rpm').dropna(); pids=set(w.index.astype(str))
Kp['pid']=Kp.pid.astype(str); Kp=Kp[Kp.pid.isin(pids)].copy()
Nall=[p+'_N' for p in sorted(pids)]
imz=score(Nall,IM); oxz=score(Nall,OX)
imt=pd.qcut(imz,3,labels=[0,1,2]).astype(int); oxt=pd.qcut(oxz,3,labels=[2,1,0]).astype(int)  # OX high=0, low=2
Kp['ax1']=np.where(Kp.type=='T',3,Kp.lib.map(imt)); Kp['ax2']=np.where(Kp.type=='T',3,Kp.lib.map(oxt))
Kp['logdepth']=np.log10(Kp.salmon_mapped)
def mixed(df,ax):
    for meth in [None,'powell']:
        try:
            m=smf.mixedlm(f'hp_log1p_rpm~{ax}+logdepth',df,groups=df.pid).fit(reml=True,**({} if meth is None else {'method':meth}))
            if np.isfinite(m.params[ax]): return m
        except Exception: pass
    return None
def boot(df,fun,ax):
    rng=np.random.default_rng(SEED); P=df.pid.unique(); g={p:df[df.pid==p] for p in P}; out=[]
    for i in range(2000):
        ss=[];
        for j,p in enumerate(rng.choice(P,len(P))): x=g[p].copy(); x['pid']=f'{p}_{j}'; ss.append(x)
        m=fun(pd.concat(ss),ax); out.append(m.params[ax] if m is not None else np.nan)
    out=np.array(out); return [float(np.nanpercentile(out,2.5)),float(np.nanpercentile(out,97.5))],int(np.isnan(out).sum())
ols=lambda df,ax: smf.ols(f'hp_log1p_rpm~{ax}+logdepth',df).fit()
def jt(y,g):
    y=np.asarray(y,float); g=np.asarray(g); lv=np.sort(np.unique(g)); J=0.0
    for a in range(len(lv)):
        for b in range(a+1,len(lv)):
            ya=y[g==lv[a]][:,None]; yb=y[g==lv[b]][None,:]; J+=(ya>yb).sum()+0.5*(ya==yb).sum()
    return J
def jt_test(y,g):
    rng=np.random.default_rng(SEED); Jo=jt(y,g); cnt=sum(jt(y,rng.permutation(g))>=Jo for _ in range(NPERM))
    n=pd.Series(g).value_counts(); Jmax=sum(n[a]*n[b] for i,a in enumerate(sorted(n.index)) for b in sorted(n.index)[i+1:])
    return dict(J_obs=float(Jo),J_over_Jmax=float(Jo/Jmax),P_one_sided_decreasing=(1+cnt)/(1+NPERM),n=int(len(y)))
res={'label':'E10 EXPLORATORY (v0.6d): pseudo-temporal ordering from molecular scores; cross-sectional, cannot establish temporal order','n_pairs':len(pids),'n_samples':len(Kp)}
tests=[]
for ax,name in [('ax1','IM_axis'),('ax2','atrophy_OX_axis')]:
    m=mixed(Kp,ax); ci,nf=boot(Kp,mixed,ax)
    r=dict(model='MixedLM REML log1p(RPM) ~ level + log10 depth + (1|pid)',slope=float(m.params[ax]),P=float(m.pvalues[ax]),CI95_boot=ci,boot_failed=nf,
      median_rpm_by_level={int(k):float(v) for k,v in Kp.groupby(ax).hp_rpm_host.median().items()},n_by_level={int(k):int(v) for k,v in Kp[ax].value_counts().items()})
    Nn=Kp[Kp.type=='N']; mo=ols(Nn,ax); cio,_=boot(Nn,lambda d,a: ols(d,a),ax)
    r['normals_only']=dict(model='OLS log1p(RPM) ~ level(0/1/2) + log10 depth (one sample/patient)',slope=float(mo.params[ax]),P=float(mo.pvalues[ax]),CI95_boot=cio,n=len(Nn))
    r['JT_ignore_pairing_descriptive']=jt_test(Kp.hp_log1p_rpm.values,Kp[ax].values)
    res[name]=r; tests+= [(name,'main',r['P']),(name,'normals_only',r['normals_only']['P'])]
q=multipletests([t[2] for t in tests],method='fdr_bh')[1]
for (n_,k,_),qq in zip(tests,q):
    (res[n_] if k=='main' else res[n_]['normals_only'])['q_BH']=float(qq)
Kp[['lib','pid','type','hp_log1p_rpm','logdepth','ax1','ax2']].to_csv(B+'results/step10/E10/E10_korea_levels.csv',index=False)
# TCGA direction only
tn=pd.read_csv(B+'results/a4_atrophy/tcga_normal_scores.csv',index_col=0)
sx=pd.read_csv(B+'data/external/spectrum2023/stad_helicobacter_per_sample.csv',index_col=0)
sx['hp']=np.log10(sx.hpylori_sp_raw/sx.homo_reads*1e6+1); sx['p']=sx.index.str[:12]; sx['ty']=sx.index.str[13:15]
g=pd.read_csv(B+'results/a4_atrophy/tcga_normal_genes.csv',index_col=0)
z=lambda X,gs:(X[gs]-X[gs].mean())/X[gs].std(ddof=1)
tn=tn.dropna(subset=['hp']); gg=g.loc[tn.index]
imT=z(gg,IM).mean(1); oxT=z(gg,OX).mean(1)
nh=sx[sx.ty=='11'].groupby('p').head(1).set_index('p').hp; th=sx[sx.ty=='01'].groupby('p').head(1).set_index('p').hp
TT=pd.DataFrame({'hp_N':nh.reindex(tn.index),'ax1':pd.qcut(imT,3,labels=[0,1,2]).astype(int),'ax2':pd.qcut(oxT,3,labels=[2,1,0]).astype(int),'hp_T':th.reindex(tn.index)})
tc={}
for ax in ['ax1','ax2']:
    long=pd.concat([pd.DataFrame({'pid':TT.index,'lvl':TT[ax],'hp':TT.hp_N}),pd.DataFrame({'pid':TT.index,'lvl':3,'hp':TT.hp_T})]).dropna()
    r=spearmanr(long.lvl,long.hp); rn=spearmanr(TT[ax],TT.hp_N)
    tc[ax]=dict(spearman_all_levels_rho=float(r.statistic),n=len(long),normals_only_rho=float(rn.statistic),n_normals=int(TT.hp_N.notna().sum()),
      n_tumour_with_hp=int(TT.hp_T.notna().sum()),median_hp_by_level={int(k):float(v) for k,v in long.groupby('lvl').hp.median().items()},note='direction only; HP=log10(Spectrum HP/homo reads*1e6+1)')
res['TCGA_direction_only']=tc
json.dump(res,open(B+'results/step10/E10/E10_results.json','w'),indent=1); print(json.dumps(res,indent=1))
