# H-M1-3 confirmatory (v0.5 section A, frozen; v0.6 label rule). Primary test only.
import pandas as pd,numpy as np,json
from scipy.stats import spearmanr,norm,t as tdist
from scipy.optimize import minimize_scalar
B='/REPLACE_WITH_PRIVATE_WORKDIR/';C=B+'results/core/'
IM=['CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4','SI','FABP1']
GD=['GKN1','GKN2','TFF2','PGC','PGA3','LIPF','CHIA','MUC6','CLDN18','GAST']
M2=['ATP4A','ATP4B','CBLIF','KCNE2','KCNJ15','KCNJ16','SLC26A7','CCKBR','ESRRG','VSIG1']
M3=['MUC5AC','TFF1','TFF2','MUC6','MUC1','FUT2','FUT3','B3GNT6']
def zscore_mean(X,genes,name,log):
    pres=[g for g in genes if g in X.columns]; miss=[g for g in genes if g not in pres]
    log[name]=dict(present=len(pres),missing=miss)
    if len(pres)/len(genes)<0.7: log[name]['FAILED_<70%']=True; return None
    Z=(X[pres]-X[pres].mean())/X[pres].std(ddof=1); return Z.mean(axis=1)
def cohort(X,hp,pairs,name):
    # X: samples x genes log2(TPM+1), index=sample; pairs: DataFrame pid,N,T sample ids
    lg={}
    im=zscore_mean(X,IM,'IM',lg); gd=zscore_mean(X,GD,'GD',lg)
    S=pd.DataFrame({'M1':im-gd,'M2':zscore_mean(X,M2,'M2',lg),'M3':zscore_mean(X,M3,'M3',lg)})
    out={'genes':lg,'n_pairs':len(pairs)}
    dh=hp[pairs['T']].values-hp[pairs['N']].values
    D=pd.DataFrame({'pid':pairs.pid.values,'dHP':dh})
    for k in ['M1','M2','M3']:
        D['d'+k]=S.loc[pairs['T'],k].values-S.loc[pairs['N'],k].values
        r=spearmanr(D.dHP,D['d'+k]); out[k]=dict(rho=float(r.statistic),P=float(r.pvalue),n=len(D),z=float(np.arctanh(r.statistic)),se=float(1/np.sqrt(len(D)-3)))
    out['dHP_zero_frac']=float((dh==0).mean()); D.to_csv(C+f'hm_deltas_{name}.csv',index=False); return out
def dl(zs,ses):
    zs=np.array(zs);v=np.array(ses)**2;w=1/v;fe=(w*zs).sum()/w.sum();Q=(w*(zs-fe)**2).sum();k=len(zs)
    tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum()))
    ws=1/(v+tau2);mu=(ws*zs).sum()/ws.sum();se=np.sqrt(1/ws.sum());z=mu/se
    return dict(mu_z=float(mu),se=float(se),rho_pooled=float(np.tanh(mu)),CI_rho=[float(np.tanh(mu-1.96*se)),float(np.tanh(mu+1.96*se))],P=float(2*norm.sf(abs(z))),tau2=float(tau2),Q=float(Q),I2=float(max(0,(Q-(k-1))/Q)) if Q>0 else 0.0)
def reml_hksj(zs,ses):
    zs=np.array(zs);v=np.array(ses)**2;k=len(zs)
    def nll(t2):
        w=1/(v+t2);mu=(w*zs).sum()/w.sum();return 0.5*(np.log(v+t2).sum()+np.log(w.sum())+(w*(zs-mu)**2).sum())
    t2=minimize_scalar(nll,bounds=(0,5),method='bounded').x
    w=1/(v+t2);mu=(w*zs).sum()/w.sum()
    q=(w*(zs-mu)**2).sum()/(k-1);se=np.sqrt(max(q,1)*0+q/w.sum()) # HKSJ (unrestricted q; also report truncated variant)
    se_tr=np.sqrt(max(q,1)/w.sum()); P=lambda s:float(2*tdist.sf(abs(mu/s),k-1))
    return dict(tau2_reml=float(t2),mu_z=float(mu),P_HKSJ=P(se),P_HKSJ_max1=P(se_tr))
# ---- Korea
K=pd.read_csv(C+'korea_rna_per_sample.csv'); E=pd.read_csv(C+'korea_gene_log2tpm.csv.gz',index_col=0).T
K=K.set_index('lib'); w=K.pivot_table(index='pid',columns='type',values='hp_log10rpm1').dropna()
kp=pd.DataFrame({'pid':w.index,'N':[p+'_N' for p in w.index],'T':[p+'_T' for p in w.index]})
korea=cohort(E,K.hp_log10rpm1,kp,'korea')
# ---- TCGA
s=pd.read_csv(B+'data/external/spectrum2023/stad_helicobacter_per_sample.csv',index_col=0)
s['hp']=np.log10(s.hpylori_sp_raw/s.homo_reads*1e6+1); s['p']=s.index.str[:12]; s['ty']=s.index.str[13:15]
e=pd.read_csv(C+'tcga_stad_hm_genes_tpm.csv'); e['k']=e.barcode.str[:16]; e=e.set_index('k')
gcols=[c for c in e.columns if c not in ['barcode','patient','sample_type']]
X=np.log2(e[gcols]+1)
# TCGA z-scores computed across all (T+N) samples with expression, as in Korea (all samples in cohort)
Nn=s[s.ty=='11'].groupby('p').head(1).set_index('p'); Tt=s[s.ty=='01'].groupby('p').head(1).set_index('p')
pp=sorted(set(Nn.index)&set(Tt.index)&set(e.patient))
tp=pd.DataFrame({'pid':pp,'N':[Nn.index.get_loc(p) and Nn.loc[p].name if False else None for p in pp]})
nid={p:s.index[(s.p==p)&(s.ty=='11')][0] for p in pp}; tid={p:s.index[(s.p==p)&(s.ty=='01')][0] for p in pp}
tp=pd.DataFrame({'pid':pp,'N':[nid[p] for p in pp],'T':[tid[p] for p in pp]})
tcga=cohort(X,s.hp,tp,'tcga')
res={'korea':korea,'tcga':tcga}
hyp={'M1':'neg','M2':'pos','M3':'two-sided'}
Pm={};meta={}
for k in ['M1','M2','M3']:
    m=dl([korea[k]['z'],tcga[k]['z']],[korea[k]['se'],tcga[k]['se']]); m.update(reml_hksj([korea[k]['z'],tcga[k]['z']],[korea[k]['se'],tcga[k]['se']])); meta[k]=m; Pm[k]=m['P']
order=sorted(Pm,key=Pm.get);hol={};run=0
for i,k in enumerate(order): run=max(run,min(1,(3-i)*Pm[k])); hol[k]=run
for k in meta:
    m=meta[k];m['P_Holm']=hol[k];m['holm_pass']=bool(hol[k]<0.05)
    sgn=lambda x:'neg' if x<0 else 'pos'
    kr,tr=korea[k]['rho'],tcga[k]['rho']; m['cohort_directions']=dict(korea=sgn(kr),tcga=sgn(tr)); m['cohorts_agree_direction']=sgn(kr)==sgn(tr)
    m['pooled_direction']=sgn(m['mu_z'])
    if hyp[k]=='two-sided': m['direction_matches_hypothesis']='n/a (two-sided; interpreted by sign)'; dm=True
    else: dm=(m['pooled_direction']==hyp[k]); m['direction_matches_hypothesis']=bool(dm)
    m['hypothesis']=hyp[k]
    m['verdict']=('Holm通过; 符号为正(ΔHP与Δ黏蛋白生态位评分同向: 生态位丢失伴HP丢失)' if (k=='M3' and m['holm_pass'] and m['pooled_direction']=='pos') else '' ) or ('支持(Holm通过且方向符合)' if (m['holm_pass'] and dm) else '不支持/未通过')
    m['label']=('稳健支持' if (m['holm_pass'] and dm and m['cohorts_agree_direction']) else ('仅合并显著' if (m['holm_pass'] and dm) else '无标签(未通过Holm或方向不符)'))
res['meta']=meta
json.dump(res,open(C+'HM_confirmatory.json','w'),indent=1); print(json.dumps(res,indent=1))
