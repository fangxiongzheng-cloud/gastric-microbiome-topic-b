# Step10 E3 (v0.6a): specificity of host-niche associations & competing explanations. EXPLORATORY. Decisions: freeze_log E2/E3 entries.
import pandas as pd,numpy as np,json,warnings
from scipy.stats import spearmanr
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor as vif
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
B='/REPLACE_WITH_PRIVATE_WORKDIR/';O=B+'results/step10/E3/';SEED=20261030
IM=['CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4','SI','FABP1']; GD=['GKN1','GKN2','TFF2','PGC','PGA3','LIPF','CHIA','MUC6','CLDN18','GAST']
M2=['ATP4A','ATP4B','CBLIF','KCNE2','KCNJ15','KCNJ16','SLC26A7','CCKBR','ESRRG','VSIG1']; M3=['MUC5AC','TFF1','TFF2','MUC6','MUC1','FUT2','FUT3','B3GNT6']
IM8=['CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4']
COMP={'prolif':['MKI67','TOP2A','CCNB1','BIRC5','CDK1'],'stroma':['COL1A1','COL1A2','COL3A1','DCN','LUM'],'epi':['EPCAM','KRT8','KRT18','CDH1']}
gm={l.split('\t')[0]:l.rstrip('\n').split('\t')[2:] for l in open('/REPLACE_WITH_REFERENCE_ROOT/msigdb/h.all.v2025.1.Hs.symbols.gmt')}
def zm(X,g):
    g=[x for x in g if x in X.columns]; Z=(X[g]-X[g].mean())/X[g].std(ddof=1); return Z.mean(1)
def mscores(X,drop=()):
    f=lambda gs:[g for g in gs if g not in drop]
    return pd.DataFrame({'M1':zm(X,f(IM))-zm(X,f(GD)),'M2':zm(X,f(M2)),'M3':zm(X,f(M3))})
# ---- Korea
E=pd.read_csv(B+'results/core/korea_gene_log2tpm.csv.gz',index_col=0).T
K=pd.read_csv(B+'results/core/korea_rna_per_sample.csv').set_index('lib')
w=K.pivot_table(index='pid',columns='type',values='hp_log10rpm1').dropna(); pids=[str(p) for p in w.index]
Nl=[p+'_N' for p in pids]; Tl=[p+'_T' for p in pids]
def build(X,hp,Nl,Tl,inf,imm):
    S=mscores(X); 
    for k,g in COMP.items(): S[k]=zm(X,g)
    S['immune_PTPRC']=(X['PTPRC']-X['PTPRC'].mean())/X['PTPRC'].std(ddof=1); S['inflam']=inf.reindex(S.index)
    S['immune_EPIC']=imm.EPIC_immune.reindex(S.index); S['immune_MCP']=imm.MCP_immune.reindex(S.index); S['hp']=hp.reindex(S.index)
    D=pd.DataFrame(S.loc[Tl].values-S.loc[Nl].values,columns=S.columns,index=[x[:-2] for x in Tl]); return S,D
infK=pd.read_csv(O+'inflam_korea.csv',index_col=0).inflam; immK=pd.read_csv(O+'immune_korea.csv',index_col=0)
S,D=build(E,K.hp_log10rpm1,Nl,Tl,infK,immK)
chk=pd.read_csv(B+'results/core/hm_deltas_korea.csv'); chk['pid']=chk.pid.astype(str); chk=chk.set_index('pid').loc[D.index]
assert np.allclose(chk.dHP,D.hp) and np.allclose(chk.dM1,D.M1) and np.allclose(chk.dM3,D.M3),'H-M delta mismatch'
z=lambda s:(s-s.mean())/s.std(ddof=1)
def rlm(df,y,xs):
    X=sm.add_constant(df[xs].apply(z)); return sm.RLM(z(df[y]),X,M=sm.robust.norms.HuberT()).fit()
CV=['prolif','inflam','stroma','immune_PTPRC','epi']
res={'label':'E3 EXPLORATORY (v0.6a)','n_pairs':len(D),'competing_scores':CV,'sens_immune':['immune_EPIC','immune_MCP']}
rng=np.random.default_rng(SEED); BI=[rng.integers(0,len(D),len(D)) for _ in range(2000)]
fam=[]
for k in ['M1','M2','M3']:
    u=rlm(D,'hp',['M'+k[1:]] if False else [k]); b0=u.params[k]
    r={'unadjusted':dict(beta=float(b0),P=float(u.pvalues[k]))}
    for c in CV: m=rlm(D,'hp',[k,c]); r['adj_'+c]=dict(beta=float(m.params[k]),P=float(m.pvalues[k]),pct_change=float(100*(m.params[k]-b0)/abs(b0)))
    m=rlm(D,'hp',[k]+CV); r['adj_all']=dict(beta=float(m.params[k]),P=float(m.pvalues[k]),pct_change=float(100*(m.params[k]-b0)/abs(b0)))
    bb=[];bc=[]
    for ix in BI:
        d=D.iloc[ix]; a=rlm(d,'hp',[k]).params[k]; c=rlm(d,'hp',[k]+CV).params[k]; bb.append(c); bc.append(100*(c-a)/abs(a))
    r['adj_all']['CI95_boot']=[float(np.percentile(bb,2.5)),float(np.percentile(bb,97.5))]; r['adj_all']['pct_change_CI95_boot']=[float(np.percentile(bc,2.5)),float(np.percentile(bc,97.5))]
    Xv=D[[k]+CV].apply(z); r['VIF_full']={c:float(vif(sm.add_constant(Xv).values,i+1)) for i,c in enumerate([k]+CV)}
    for s_ in ['immune_EPIC','immune_MCP']:
        cv2=[c if c!='immune_PTPRC' else s_ for c in CV]; m=rlm(D.dropna(subset=[s_]),'hp',[k]+cv2); r['adj_all_'+s_]=dict(beta=float(m.params[k]),P=float(m.pvalues[k]))
    res['E3a_'+k]=r; fam.append(('E3a_'+k,'adj_all',r['adj_all']['P']))
# E3b
allc=set(sum(COMP.values(),[]))|{'PTPRC'}|set(gm['HALLMARK_INFLAMMATORY_RESPONSE'])
ov={k:sorted(set(g)&allc) for k,g in {'M1':IM+GD,'M2':M2,'M3':M3}.items()}
S2=mscores(E,drop=set(sum(ov.values(),[]))); D2=S2.loc[Tl].values-S2.loc[Nl].values
e3b={'overlap_genes':ov}
for i,k in enumerate(['M1','M2','M3']):
    rr=spearmanr(D.hp,D2[:,i]); e3b[k]={'rho_overlap_removed':float(rr.statistic),'P':float(rr.pvalue),'rho_original':float(spearmanr(D.hp,D[k]).statistic)}
    genes=IM+GD if k=='M1' else (M2 if k=='M2' else M3); lo=[]
    for g in dict.fromkeys(genes):
        if g not in E.columns: continue
        S3=mscores(E,drop={g}); lo.append(spearmanr(D.hp,S3.loc[Tl,k].values-S3.loc[Nl,k].values).statistic)
    e3b[k]['leave_one_gene_out_rho_range']=[float(min(lo)),float(max(lo))]
res['E3b']=e3b
# E3c
hpN=K.loc[Nl,'hp_log10rpm1'].values; e3c={}
for k in ['M1','M2','M3']:
    r1=spearmanr(hpN,D[k]); e3c['i_normalHP_vs_dM_'+k]=dict(rho=float(r1.statistic),P=float(r1.pvalue)); fam.append(('E3c','i_'+k,r1.pvalue))
SN=S.loc[Nl].copy(); SN['IM_A4']=zm(E.loc[Nl],IM8)
for k in ['M2','M3','IM_A4']:
    r2=spearmanr(hpN,SN[k]); e3c['ii_withinNormal_HP_vs_'+k]=dict(rho=float(r2.statistic),P=float(r2.pvalue)); fam.append(('E3c','ii_'+k,r2.pvalue))
res['E3c']=e3c
# E3d
e3d={}; G=E[[g for g in M3 if g in E.columns]]; dG=G.loc[Tl].values-G.loc[Nl].values
for j,g in enumerate(G.columns):
    r=spearmanr(D.hp,dG[:,j]); e3d[g]=dict(rho=float(r.statistic),P=float(r.pvalue)); fam.append(('E3d',g,r.pvalue))
Tall=[l for l in K.index if l.endswith('_T')]; r=spearmanr(K.loc[Tall,'hp_log10rpm1'],E.loc[Tall,'MUC5AC'])
e3d['Oosterlinck_tumour_crosssectional_HP_vs_MUC5AC']=dict(n=len(Tall),rho=float(r.statistic),P=float(r.pvalue)); fam.append(('E3d','Oost',r.pvalue))
res['E3d']=e3d
q=multipletests([f[2] for f in fam],method='fdr_bh')[1]; res['BH_family_E3']=[dict(item=f[0],test=f[1],P=float(f[2]),q=float(qq)) for f,qq in zip(fam,q)]
# TCGA direction
Tt=pd.read_csv(O+'tcga_selected_log2tpm.csv',index_col=0).T; hm=pd.read_csv(B+'results/core/tcga_stad_hm_genes_tpm.csv'); 
td=pd.read_csv(B+'results/core/hm_deltas_tcga.csv').set_index('pid')
infT=pd.read_csv(O+'inflam_tcga.csv',index_col=0).inflam; immT=pd.read_csv(O+'immune_tcga.csv',index_col=0)
Dt=pd.DataFrame(index=td.index); Dt['hp']=td.dHP
for k in ['M1','M2','M3']: Dt[k]=td['d'+k]
for c,g in COMP.items(): s=zm(Tt,g); Dt[c]=[s[p+'_T']-s[p+'_N'] for p in td.index]
s=(Tt.PTPRC-Tt.PTPRC.mean())/Tt.PTPRC.std(ddof=1); Dt['immune_PTPRC']=[s[p+'_T']-s[p+'_N'] for p in td.index]
Dt['inflam']=[infT[p+'_T']-infT[p+'_N'] for p in td.index]
res['TCGA_direction_only']={k:dict(unadj=float(rlm(Dt,'hp',[k]).params[k]),adj_all=float(rlm(Dt,'hp',[k]+CV).params[k])) for k in ['M1','M2','M3']}
res['TCGA_direction_only']['n_pairs']=len(Dt)
D.to_csv(O+'E3_korea_deltas.csv'); json.dump(res,open(O+'E3_results.json','w'),indent=1)
for k in ['M1','M2','M3']: print(k,res['E3a_'+k]['unadjusted'],res['E3a_'+k]['adj_all'],res['E3a_'+k]['VIF_full'])
print(json.dumps({k:res[k] for k in ['E3b','E3c','E3d','TCGA_direction_only']},indent=0)[:3000])
