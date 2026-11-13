import pandas as pd,numpy as np,json
from scipy.stats import wilcoxon,kruskal,chi2_contingency,mannwhitneyu,norm,t as tdist
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
rng=np.random.default_rng(20261030); O='results/korea16s/'
CB={'N':'#0072B2','T':'#D55E00','g':'#999999','y':'#E69F00','gr':'#009E73','p':'#CC79A7'}
plt.rcParams.update({'font.size':8,'pdf.fonttype':42})
# ---------- TCGA dose-response S1/S2 ----------
df=pd.read_csv('results/stage1/hp_calls.csv',index_col=0); T=df[df.tissue=='T'].copy()
cl=pd.read_csv('stad_subtype.csv').set_index('patient'); T=T.join(cl[['Lauren.Class','TNM.Stage']],how='left')
T['sp_log']=np.log1p(T.sp_hp); T['cons']=T.hp_consensus.map({'pos':1,'neg':0}); T['plate']=T.plate.fillna('NA')
st=T['TNM.Stage'].astype(str).str.extract(r'Stage_(IV|III|II|I)')[0]; T['stage']=st.map({'I':1,'II':2,'III':3,'IV':4})
def jt(x,g):
    lv=sorted(np.unique(g)); s=0
    for i in range(len(lv)):
        for j in range(i+1,len(lv)):
            a=x[g==lv[i]];b=x[g==lv[j]]; s+=sum((bb>a).sum()+0.5*(bb==a).sum() for bb in b)
    return s
def jt_perm(x,g,pl,n=1000):
    obs=jt(x,g); c=0
    for _ in range(n):
        g2=g.copy()
        for p in np.unique(pl): ix=np.where(pl==p)[0]; g2[ix]=rng.permutation(g2[ix])
        e=jt(x,g2); c+= e>=obs
    null_mean=None; return obs,c
res=[]
L=T[T['Lauren.Class'].isin(['Intestinal','Diffuse','Mixed'])]
for meas in ['cons','sp_log']:
    d=L[L[meas].notna()]
    if meas=='cons':
        tb=pd.crosstab(d['Lauren.Class'],d.cons); p=chi2_contingency(tb)[1]
        res.append(dict(test='S1 Lauren',measure='consensus pos rate',n=len(d),detail='; '.join(f'{k}: {int(tb.loc[k,1])}/{int(tb.loc[k].sum())}' for k in tb.index),P=p))
        a=d[d['Lauren.Class']=='Intestinal'].cons; b=d[d['Lauren.Class']=='Diffuse'].cons
    else:
        p=kruskal(*[g.sp_log for _,g in d.groupby('Lauren.Class')]).pvalue
        res.append(dict(test='S1 Lauren',measure='Spectrum log1p',n=len(d),detail='; '.join(f'{k}: median {np.expm1(g.sp_log.median()):.0f} reads (n={len(g)})' for k,g in d.groupby('Lauren.Class')),P=p))
    # stage trend within-plate permutation JT
    s=T[T.stage.notna()&T[meas].notna()]
    x=s[meas].values; g=s.stage.values.astype(int); pl=s.plate.values
    obs=jt(x,g); null=[]
    for _ in range(1000):
        g2=g.copy()
        for pp in np.unique(pl): ix=np.where(pl==pp)[0]; g2[ix]=rng.permutation(g2[ix])
        null.append(jt(x,g2))
    null=np.array(null); pv=(np.sum(np.abs(null-null.mean())>=abs(obs-null.mean()))+1)/1001
    det='; '.join(f'St{k}: '+(f'{s[s.stage==k].cons.sum():.0f}/{(s.stage==k).sum()}' if meas=='cons' else f'mean log1p {s[s.stage==k].sp_log.mean():.2f} (n={(s.stage==k).sum()})') for k in [1,2,3,4])
    res.append(dict(test='S2 stage I-IV trend (JT, within-plate perm, two-sided)',measure=meas,n=len(s),detail=det+f'; JT obs-nullmean={obs-null.mean():+.0f}',P=pv))
R=pd.DataFrame(res); R['q_BH']=multipletests(R.P,method='fdr_bh')[1]; R.to_csv(O+'dose_response_TCGA.csv',index=False); print(R.to_string())
# ---------- Korea per-pair & meta ----------
k=pd.read_csv(O+'per_sample_silva.csv',index_col=0); w=k.pivot_table(index='pid',columns='type_bs',values=['clr','ra']).dropna()
def rb(N,Tt):
    d=N-Tt; nz=d[d!=0]; rk=nz.abs().rank(); return (rk[nz>0].sum()-rk[nz<0].sum())/rk.sum(), len(nz)
tc=df.copy(); pairs={}
for col,nm in [('sp_hp','TCGA Spectrum (RNA)'),('tcma_wgs_rpm','TCGA TCMA (WGS)')]:
    pr=pd.concat([df[df.tissue=='T'][col].rename('T'),df[df.tissue=='N'][col].rename('N')],axis=1).dropna(); pairs[nm]=pr
pairs['Korea 16S (SILVA)']=pd.DataFrame({'N':w['clr']['N'],'T':w['clr']['T']})
g2=pd.read_csv(O+'per_sample_gg2v4.csv',index_col=0).pivot_table(index='pid',columns='type_bs',values='clr').dropna()
meta=[]
for nm,pr in pairs.items():
    r,n=rb(pr.N,pr['T']); z=np.arctanh(r); se=1/np.sqrt(max(n-3,1)); meta.append(dict(study=nm,pairs=len(pr),nonzero=n,r=r,z=z,se=se,P=wilcoxon(pr.N,pr['T']).pvalue))
M=pd.DataFrame(meta); wv=1/M.se**2; zf=(wv*M.z).sum()/wv.sum(); Q=(wv*(M.z-zf)**2).sum(); k_=len(M)
tau2=max(0,(Q-(k_-1))/(wv.sum()-(wv**2).sum()/wv.sum())); ws=1/(M.se**2+tau2); zr=(ws*M.z).sum()/ws.sum()
q_hk=(ws*(M.z-zr)**2).sum()/(k_-1); se_hk=np.sqrt(q_hk/ws.sum()); tcrit=tdist.ppf(.975,k_-1)
pooled=dict(r=np.tanh(zr),lo=np.tanh(zr-tcrit*se_hk),hi=np.tanh(zr+tcrit*se_hk),P=2*tdist.sf(abs(zr/se_hk),k_-1),tau2=tau2,Q=Q,I2=max(0,(Q-(k_-1))/Q) if Q>0 else 0)
M.to_csv(O+'meta_gradient.csv',index=False); json.dump(pooled,open(O+'meta_gradient_pooled.json','w'),indent=1); print(M,pooled)
# ---------- Fig 2 ----------
with PdfPages(O+'figs/Fig2_gradient_draft.pdf') as pdf:
    fig,ax=plt.subplots(2,3,figsize=(7.2,5))
    for a,(nm,pr) in zip(ax[0],[('TCGA Spectrum (RNA)',np.log1p(pairs['TCGA Spectrum (RNA)'])),('TCGA TCMA (WGS)',np.log1p(pairs['TCGA TCMA (WGS)'])),('Korea 16S (SILVA)',pairs['Korea 16S (SILVA)'])]):
        for _,r in pr.iterrows(): a.plot([0,1],[r.N,r['T']],color=CB['g'],lw=.4,alpha=.5)
        a.scatter(np.zeros(len(pr)),pr.N,s=5,color=CB['N']); a.scatter(np.ones(len(pr)),pr['T'],s=5,color=CB['T'])
        a.set_xticks([0,1]); a.set_xticklabels(['Adj. normal','Tumour']); a.set_title(f'{nm}\nn={len(pr)} pairs, P={M.set_index("study").loc[nm,"P"]:.1e}',fontsize=7)
        a.set_ylabel('log1p(reads)' if 'TCGA Spec' in nm else ('log1p(RPM)' if 'TCMA' in nm else 'Helicobacter CLR'))
    a=ax[1,0]; d=pairs['Korea 16S (SILVA)']; diff=(d.N-d['T']).sort_values()
    a.bar(range(len(diff)),diff,color=np.where(diff>0,CB['N'],CB['T']),width=1); a.axhline(0,color='k',lw=.5); a.set_title('Korea: per-pair ΔCLR (N−T)',fontsize=7); a.set_xlabel('pairs (sorted)')
    a=ax[1,1]; a.scatter(d.N-d['T'],g2.loc[d.index,'N']-g2.loc[d.index,'T'],s=4,color=CB['gr']); a.set_xlabel('ΔCLR SILVA'); a.set_ylabel('ΔCLR GG2'); a.set_title('Classifier agreement',fontsize=7)
    a=ax[1,2]; y=np.arange(len(M))[::-1]
    lo=np.tanh(M.z-1.96*M.se); hi=np.tanh(M.z+1.96*M.se)
    a.errorbar(M.r,y,xerr=[M.r-lo,hi-M.r],fmt='s',color='k',ms=3); a.errorbar([pooled['r']],[-1],xerr=[[pooled['r']-pooled['lo']],[pooled['hi']-pooled['r']]],fmt='D',color=CB['p'],ms=4)
    a.set_yticks(list(y)+[-1]); a.set_yticklabels(list(M.study)+['RE pooled (HK)']); a.axvline(0,color='k',lw=.5); a.set_xlabel('rank-biserial r (N>T positive)'); a.set_title('Random-effects meta',fontsize=7)
    fig.tight_layout(); pdf.savefig(fig)
# ---------- Fig 3: Korea community ----------
al=pd.read_csv(O+'da/alpha.csv',index_col=0); pc=pd.read_csv(O+'da/pcoa_bray.csv',index_col=0)
C=pd.read_csv(O+'da/da_all_silva.csv'); cons=pd.read_csv(O+'da/da_consensus_silva.csv')
B=pd.read_csv(O+'da/beta_tests.csv'); A=pd.read_csv(O+'da/alpha_tests.csv')
with PdfPages(O+'figs/Fig3_korea_community_draft.pdf') as pdf:
    fig,ax=plt.subplots(2,3,figsize=(7.2,5))
    for a,m,lab in zip(ax[0,:2],['shannon','obs'],['Shannon','Observed genera']):
        ww=al.pivot_table(index='pid',columns='type',values=m)
        for _,r in ww.iterrows(): a.plot([0,1],[r.N,r['T']],color=CB['g'],lw=.3,alpha=.5)
        a.boxplot([ww.N,ww['T']],positions=[0,1],widths=.35,showfliers=False); a.set_xticks([0,1]); a.set_xticklabels(['N','T'])
        a.set_title(f'{lab} (rarefied 5k)\nP={A.set_index("metric").loc[m,"P"]:.1e}',fontsize=7)
    a=ax[0,2]; 
    for t in ['N','T']: s=pc[pc.type==t]; a.scatter(s.PC1,s.PC2,s=4,color=CB[t],label=t,alpha=.7)
    a.legend(fontsize=6); a.set_title(f'Bray PCoA; PERMANOVA (blocked by patient)\nR²={B.R2[0]:.3f}, P={B.P[0]:.3f}',fontsize=7)
    a=ax[1,0]; sig=np.maximum(-np.log10(C.mas_q.clip(1e-20)),0)
    col=[CB['p'] if t in set(cons.taxon) else CB['g'] for t in C.taxon]
    a.scatter(C.mas_coef,sig,s=6,c=col); a.axhline(-np.log10(.05),ls='--',color='k',lw=.5)
    for _,r in cons.iterrows(): a.annotate(r.taxon,(r.mas_coef,-np.log10(max(r.mas_q,1e-20))),fontsize=5)
    a.set_xlabel('MaAsLin2 coef (T vs N)'); a.set_ylabel('-log10 q'); a.set_title('DA; pink = ≥2/3 methods agree',fontsize=7)
    a=ax[1,1]; oral=['Fusobacterium','Streptococcus','Prevotella','Prevotella_7','Peptostreptococcus','Veillonella','Porphyromonas','Lactobacillus']
    o=C[C.taxon.isin(oral)].set_index('taxon').reindex(oral).dropna(how='all')
    mat=o[['anc_lfc','mas_coef','ald_eff']].values; im=a.imshow(mat,cmap='PuOr_r',vmin=-1.6,vmax=1.6,aspect='auto')
    a.set_yticks(range(len(o))); a.set_yticklabels(o.index,fontsize=6); a.set_xticks([0,1,2]); a.set_xticklabels(['ANCOM-BC2','MaAsLin2','ALDEx2'],fontsize=6)
    for i in range(len(o)):
        for j,q in enumerate(['anc_q','mas_q','ald_q']):
            if o.iloc[i][q]<0.05: a.text(j,i,'*',ha='center',va='center')
    a.set_title('Oral taxa effect (T vs N); * q<0.05',fontsize=7); fig.colorbar(im,ax=a,fraction=.05)
    a=ax[1,2]; ra=pd.read_csv(O+'per_sample_silva.csv',index_col=0)
    a.hist([ra[ra.type_bs=='N'].ra,ra[ra.type_bs=='T'].ra],bins=20,color=[CB['N'],CB['T']],label=['N','T']); a.legend(fontsize=6); a.set_xlabel('Helicobacter relative abundance'); a.set_title('Helicobacter dominance',fontsize=7)
    fig.tight_layout(); pdf.savefig(fig)
