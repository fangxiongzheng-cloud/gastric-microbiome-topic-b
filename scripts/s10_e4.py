# Step10 E4 (v0.6a): RNA-16S cross-platform calibration, Korean same-BioSample. Exploratory. Decisions: freeze_log 2026-11-11 10:20.
import pandas as pd,numpy as np,json
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
B='/REPLACE_WITH_PRIVATE_WORKDIR/';O=B+'results/step10/E4/';SEED=20261030
D=pd.read_csv(B+'results/core/korea_rna_per_sample.csv')
D=D.dropna(subset=['s16_clr'])
res={'label':'E4 EXPLORATORY (v0.6a)','n_samples_with_both':len(D)}
rng=np.random.default_rng(SEED)
def boot_rho(df,x,y,B_=2000):
    pids=df.pid.unique(); g={p:df[df.pid==p] for p in pids}; out=[]
    for _ in range(B_):
        s=pd.concat([g[p] for p in rng.choice(pids,len(pids))]); out.append(spearmanr(s[x],s[y]).statistic)
    return [float(np.nanpercentile(out,2.5)),float(np.nanpercentile(out,97.5))]
P=[]
for t in ['N','T']:
    d=D[D.type==t]; r=spearmanr(d.hp_log1p_rpm,d.s16_clr)
    res['spearman_'+t]=dict(n=len(d),rho=float(r.statistic),P=float(r.pvalue),CI95_boot=boot_rho(d,'hp_log1p_rpm','s16_clr')); P.append(r.pvalue)
w=D.pivot_table(index='pid',columns='type',values=['hp_log1p_rpm','s16_clr','hp_k2_confirmed']).dropna()
pr=pd.DataFrame({'pid':w.index,'dRNA':w['hp_log1p_rpm']['T']-w['hp_log1p_rpm']['N'],'d16S':w['s16_clr']['T']-w['s16_clr']['N'],
 'maxk2':np.maximum(w['hp_k2_confirmed']['T'],w['hp_k2_confirmed']['N'])}).reset_index(drop=True)
r=spearmanr(pr.dRNA,pr.d16S); P.append(r.pvalue)
bs=[];cc=[];cc2=[]
conc=lambda d:float((np.sign(d.dRNA)==np.sign(d.d16S)).mean())
conc2=lambda d:(lambda e:float((np.sign(e.dRNA)==np.sign(e.d16S)).mean()))(d[(d.dRNA!=0)&(d.d16S!=0)])
for _ in range(2000):
    s=pr.iloc[rng.integers(0,len(pr),len(pr))]; bs.append(spearmanr(s.dRNA,s.d16S).statistic); cc.append(conc(s)); cc2.append(conc2(s))
ci=lambda v:[float(np.nanpercentile(v,2.5)),float(np.nanpercentile(v,97.5))]
nz=pr[(pr.dRNA!=0)&(pr.d16S!=0)]
res['paired_delta']=dict(n_pairs=len(pr),rho=float(r.statistic),P=float(r.pvalue),CI95_boot=ci(bs),
  sign_concordance_all=conc(pr),CI95=ci(cc),sign_concordance_both_nonzero=conc2(pr),n_both_nonzero=len(nz),CI95_nonzero=ci(cc2),
  n_dRNA_zero=int((pr.dRNA==0).sum()),n_d16S_zero=int((pr.d16S==0).sum()),
  direction_table={f'RNA{a}_16S{b}':int(((np.sign(pr.dRNA)==a)&(np.sign(pr.d16S)==b)).sum()) for a in [-1,0,1] for b in [-1,0,1]})
q=multipletests(P,method='fdr_bh')[1]
res['spearman_N']['q_BH']=float(q[0]);res['spearman_T']['q_BH']=float(q[1]);res['paired_delta']['q_BH']=float(q[2])
# strata
D['stratum']=pd.cut(D.hp_k2_confirmed,[-1,9,999,np.inf],labels=['<10','10-999','>=1000'])
D['rna_det']=D.hp_k2_confirmed>=10; D['s16_det']=D.s16_hp_reads>=10; D['det_discord']=D.rna_det!=D.s16_det
D['pr_rna']=D.hp_log1p_rpm.rank(pct=True); D['pr_16s']=D.s16_clr.rank(pct=True); D['rank_discord']=(D.pr_rna-D.pr_16s).abs()>0.5
au=pd.read_csv(B+'results/audit/audit_breadth_26695.csv')[['run','breadth_ratio']]; D=D.merge(au,on='run',how='left')
st={}
for s,d in D.groupby('stratum',observed=True):
    st[str(s)]=dict(n=len(d),n_N=int((d.type=='N').sum()),n_T=int((d.type=='T').sum()),det_discord_rate=float(d.det_discord.mean()),
      n_det_discord=int(d.det_discord.sum()),rank_discord_rate=float(d.rank_discord.mean()),
      s16_detect_rate=float(d.s16_det.mean()),median_s16_clr=float(d.s16_clr.median()))
res['strata_by_rna_k2']=st
pr['stratum']=pd.cut(pr.maxk2,[-1,9,999,np.inf],labels=['<10','10-999','>=1000'])
res['pair_direction_discord_by_max_k2_stratum']={str(s):dict(n=len(d),discord_rate=float((np.sign(d.dRNA)!=np.sign(d.d16S)).mean())) for s,d in pr.groupby('stratum',observed=True)}
desc={}
for k,d in D.groupby('det_discord'):
    desc['discordant' if k else 'concordant']=dict(n=len(d),median_salmon_mapped=float(d.salmon_mapped.median()),median_salmon_pct=float(d.salmon_pct.median()),
       n_audit=int(d.breadth_ratio.notna().sum()),median_breadth_ratio=float(d.breadth_ratio.median()) if d.breadth_ratio.notna().any() else None,
       median_16S_depth=None)
res['detection_discordance_descriptors']=desc
res['note']='Both platforms are relative measures; agreement != accurate bacterial load.'
D.to_csv(O+'E4_per_sample.csv',index=False); pr.to_csv(O+'E4_pairs.csv',index=False)
json.dump(res,open(O+'E4_results.json','w'),indent=1,ensure_ascii=False); print(json.dumps(res,indent=1,ensure_ascii=False))
