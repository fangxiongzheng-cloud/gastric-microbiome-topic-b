# Step10 E8 analysis (v0.6b E8.4). EXPLORATORY / method performance for R1. Decisions: freeze_log E8 entry.
import pandas as pd,numpy as np,json,gzip,glob,os
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests
B='/REPLACE_WITH_PRIVATE_WORKDIR/';O=B+'results/step10/E8/'
LEV=[10,30,100,300,1000,3000,10000,30000]; SEEDS=[20261030,20261031,20261032]
TG=['T1_26695_uniform','T2_26695_rRNA90','T3_F30_offref','T4a_acinonychis','T4b_hepaticus']
def tree(report):
    rows=[l.rstrip('\n').split('\t') for l in open(report)]; clade={'210':set(),'209':set()}; stack=[]
    for r in rows:
        if len(r)<6: continue
        depth=(len(r[5])-len(r[5].lstrip(' ')))//2; tid=r[4]
        while stack and stack[-1][0]>=depth: stack.pop()
        stack.append((depth,tid)); anc={t for _,t in stack}
        for c in clade:
            if c in anc: clade[c].add(tid)
    return clade
def k2clade(report,taxid='210'):
    for l in open(report):
        x=l.rstrip('\n').split('\t')
        if len(x)>4 and x[4]==taxid: return int(x[1])
    return 0
# backgrounds
bg=[]
for d in sorted(glob.glob(O+'bg/*/')):
    if not os.path.exists(d+'DONE'): continue
    r=os.path.basename(d.rstrip('/')); mi=json.load(open(d+'salmon/aux_info/meta_info.json'))
    bg.append(dict(run=r,pairs=int(open(d+'n_pairs.txt').read()),salmon_mapped=mi['num_mapped'],bg_confirmed=k2clade(d+'k2.report'),
                   bg_bt2=int(os.popen(f"/REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin/samtools view -F 0x904 {d}hp.bam | cut -f1 | sort -u | wc -l").read())))
BG=pd.DataFrame(bg); m=pd.read_csv(O+'backgrounds.csv')[['run','type']]; BG=BG.merge(m,on='run'); BG.to_csv(O+'E8_backgrounds_summary.csv',index=False)
# sims
rows=[]
for T in TG:
    for sd in SEEDS:
        d=O+f'sim/{T}_{sd}/'
        if not os.path.exists(d+'DONE'): continue
        cl=tree(d+'k2.report')
        k=pd.read_csv(d+'k2.out.txt.gz',sep='\t',header=None,usecols=[0,1,2],names=['c','name','tid'],dtype=str)
        k['L']=k.name.str.extract(r'^L(\d+)_')[0].astype(int)
        b=pd.read_csv(d+'bt2_primary_name_mapq.tsv',sep='\t',header=None,names=['name','mapq']); b['name']=b.name.str.replace(r'/[12]$','',regex=True)
        bb=b.groupby('name').mapq.max().reset_index(); bb['L']=bb.name.str.extract(r'^L(\d+)_')[0].astype(int)
        for L in LEV:
            kk=k[k.L==L]; x=bb[bb.L==L]
            c210=int(kk.tid.isin(cl['210']).sum()); c209=int((kk.tid.isin(cl['209'])&~kk.tid.isin(cl['210'])).sum())
            uncl=int((kk.tid=='0').sum()); other=len(kk)-c210-c209-uncl
            rows.append(dict(target=T,seed=sd,level=L,bt2_aligned=len(x),mapq30_frac=float((x.mapq>=30).mean()) if len(x) else np.nan,
              k2_input=len(kk),k2_210=c210,k2_209_not210=c209,k2_other=other,k2_unclassified=uncl))
SIM=pd.DataFrame(rows); SIM.to_csv(O+'E8_sim_counts.csv',index=False)
# grid
G=[]
for _,b in BG.iterrows():
    G.append(dict(run=b.run,target='none',seed=0,level=0,sim_conf=0,total=b.bg_confirmed,denom=b.salmon_mapped))
    for _,s in SIM.iterrows(): G.append(dict(run=b.run,target=s.target,seed=s.seed,level=s.level,sim_conf=s.k2_210,total=b.bg_confirmed+s.k2_210,denom=b.salmon_mapped))
G=pd.DataFrame(G); G['spiked_rpm']=G.level/G.denom*1e6; G['obs_rpm_incr']=G.sim_conf/G.denom*1e6; G['total_rpm']=G.total/G.denom*1e6
G.to_csv(O+'E8_grid.csv',index=False)
res={'label':'E8 EXPLORATORY method-performance simulation (v0.6b); simulated reads lack transcription/degradation/library bias - pipeline detection only, not absolute load',
     'backgrounds':BG.to_dict('records'),'denominator_range':[int(BG.salmon_mapped.min()),int(BG.salmon_mapped.max())]}
res['LoB_95pct_confirmed_at_0']=float(np.percentile(BG.bg_confirmed,95)); res['bg_confirmed_counts']=BG.bg_confirmed.tolist()
out={};slopes=[]
for T in TG:
    g=G[G.target==T]; r={}
    det={L:float((g[g.level==L].total>=10).mean()) for L in LEV}; r['detect_rate_ge10_by_level']=det
    lod=next((L for L in LEV if det[L]>=0.95),None); r['LoD_pairs']=lod
    r['LoD_rpm_range']=None if lod is None else [float(g[g.level==lod].spiked_rpm.min()),float(g[g.level==lod].spiked_rpm.max())]
    s=SIM[SIM.target==T]; r['recovery_by_level']={L:dict(median=float((s[s.level==L].k2_210/L).median()),min=float((s[s.level==L].k2_210/L).min()),max=float((s[s.level==L].k2_210/L).max())) for L in LEV}
    r['bt2_aligned_frac_median']=float((s.bt2_aligned/s.level).median()); r['mapq30_frac_median']=float(s.mapq30_frac.median())
    tot=s[['k2_input','k2_210','k2_209_not210','k2_other','k2_unclassified']].sum(); r['kraken_breakdown_of_bt2_aligned']={k:float(tot[k]/tot.k2_input) if tot.k2_input else None for k in ['k2_210','k2_209_not210','k2_other','k2_unclassified']}
    if lod is not None:
        h=g[(g.level>=lod)]; n0=int((h.obs_rpm_incr<=0).sum()); h=h[h.obs_rpm_incr>0]
        X=sm.add_constant(np.log10(h.spiked_rpm)); f=sm.OLS(np.log10(h.obs_rpm_incr),X).fit(); b1=f.params.iloc[1]; se=f.bse.iloc[1]
        tt=(b1-1)/se; p=2*stats.t.sf(abs(tt),f.df_resid); ci=f.conf_int().iloc[1].tolist()
        r['linearity']=dict(n_points=len(h),n_zero_excluded=n0,slope=float(b1),slope_CI95=[float(ci[0]),float(ci[1])],intercept=float(f.params.iloc[0]),R2=float(f.rsquared),P_slope_eq_1=float(p),reference='slope 0.9-1.1 descriptive only')
        slopes.append((T,p))
    if T.startswith('T4'): r['false_positive_frac_by_level']={L:float(s[s.level==L].k2_210.sum()/(3*L)) for L in LEV}; r['fp_rpm_contribution_median_by_level']={L:float(g[g.level==L].obs_rpm_incr.median()) for L in LEV}
    out[T]=r
if slopes:
    q=multipletests([p for _,p in slopes],method='fdr_bh')[1]
    for (T,_),qq in zip(slopes,q): out[T]['linearity']['q_BH']=float(qq)
s1=SIM[SIM.target=='T1_26695_uniform'].set_index(['seed','level']).k2_210; s3=SIM[SIM.target=='T3_F30_offref'].set_index(['seed','level']).k2_210
rat=(s3/s1.replace(0,np.nan)).reset_index(name='ratio'); res['offref_recovery_ratio_T3_over_T1']={'by_level_median':rat.groupby('level').ratio.median().to_dict(),'pooled':float(s3.sum()/s1.sum())}
res['targets']=out
json.dump(res,open(O+'E8_results.json','w'),indent=1,default=float); print(json.dumps(res,indent=1,default=float))
