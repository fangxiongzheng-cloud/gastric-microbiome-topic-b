import pandas as pd,numpy as np
from scipy.stats import mannwhitneyu,fisher_exact,wilcoxon,spearmanr
rng=np.random.default_rng(20261030)
df=pd.read_csv("results/stage1/hp_calls.csv",index_col=0)
T=df[(df.tissue=="T")&df.subtype.notna()].copy()
T["sp_log"]=np.log1p(T.sp_hp); T["sp_rpm_bact"]=np.log1p(T.sp_hp/T.sp_bact*1e6)
T["tcma_log"]=np.log1p(T.tcma_wgs_rpm); T["sp_dec_log"]=np.log1p(T.sp_hp_genus_decontam); T["poore_log"]=np.log1p(T.poore_raw)
T["cons_bin"]=T.hp_consensus.map({"pos":1,"neg":0})
T["plate"]=T.plate.fillna("NA")
measures={"Consensus(pos/neg)":"cons_bin","Spectrum H.pylori raw":"sp_log","Spectrum H.pylori /bact":"sp_rpm_bact","Spectrum Helicobacter decontam":"sp_dec_log","TCMA WGS H.pylori RPM":"tcma_log","Poore raw (retracted;contrast)":"poore_log"}
contr={"EBV vs rest":(["EBV"],None),"MSI vs rest":(["MSI"],None),"EBV+MSI vs GS+CIN":(["EBV","MSI"],["GS","CIN"]),"GS vs rest":(["GS"],None)}
def hedges(a,b):
    na,nb=len(a),len(b); s=np.sqrt(((na-1)*a.var(ddof=1)+(nb-1)*b.var(ddof=1))/(na+nb-2))
    g=(a.mean()-b.mean())/s if s>0 else np.nan; return g*(1-3/(4*(na+nb)-9))
def stat(v,lab): return v[lab].mean()-v[~lab].mean()
out=[]
for mn,col in measures.items():
    for cn,(A,B) in contr.items():
        d=T[T[col].notna()]
        if B: d=d[d.subtype.isin(A+B)]
        lab=d.subtype.isin(A).values; v=d[col].values
        if lab.sum()<3 or (~lab).sum()<3: continue
        obs=stat(pd.Series(v),pd.Series(lab))
        # within-plate permutation
        perm=[];pl=d.plate.values
        for _ in range(1000):
            l2=lab.copy()
            for p in np.unique(pl):
                ix=np.where(pl==p)[0]; l2[ix]=rng.permutation(l2[ix])
            perm.append(v[l2].mean()-v[~l2].mean())
        perm=np.array(perm); pp=(np.sum(np.abs(perm)>=abs(obs))+1)/1001
        if col=="cons_bin":
            tb=[[int(v[lab].sum()),int(lab.sum()-v[lab].sum())],[int(v[~lab].sum()),int((~lab).sum()-v[~lab].sum())]]
            OR,p=fisher_exact(tb); eff=f"pos {tb[0][0]}/{lab.sum()} vs {tb[1][0]}/{(~lab).sum()}; OR={OR:.2f}"; g=np.nan
        else:
            p=mannwhitneyu(v[lab],v[~lab]).pvalue; g=hedges(v[lab],v[~lab]); eff=f"g={g:.2f}"
        out.append(dict(measure=mn,contrast=cn,n_group=int(lab.sum()),n_ref=int((~lab).sum()),effect=eff,hedges_g=g,diff=obs,p_test=p,p_perm_withinplate=pp))
r=pd.DataFrame(out)
from statsmodels.stats.multitest import multipletests
r["q_BH_test"]=multipletests(r.p_test,method="fdr_bh")[1]
r.to_csv("results/stage1/H1_subtype_tests.csv",index=False)
pd.set_option("display.width",250);print(r.drop(columns=["diff"]).round(4).to_string())
# H3 tumour vs normal
N=df[df.tissue=="N"]; Tt=df[df.tissue=="T"]
res=[]
for col,nm in [("sp_hp","Spectrum H.pylori raw"),("tcma_wgs_rpm","TCMA WGS RPM")]:
    pr=pd.concat([Tt[col].rename("T"),N[col].rename("N")],axis=1).dropna()
    w=wilcoxon(pr.N,pr["T"]) if (pr.N!=pr["T"]).any() else None
    res.append(dict(measure=nm,pairs=len(pr),N_gt_T=int((pr.N>pr["T"]).sum()),T_gt_N=int((pr["T"]>pr.N).sum()),ties=int((pr.N==pr["T"]).sum()),det_T=round((pr["T"]>0).mean(),3),det_N=round((pr.N>0).mean(),3),median_T=pr["T"].median(),median_N=pr.N.median(),p_wilcoxon=w.pvalue if w else np.nan))
    # unpaired
    a=Tt[col].dropna();b=N[col].dropna()
    res.append(dict(measure=nm+" (unpaired)",pairs=f"T{len(a)}/N{len(b)}",det_T=round((a>0).mean(),3),det_N=round((b>0).mean(),3),median_T=a.median(),median_N=b.median(),p_wilcoxon=mannwhitneyu(a,b).pvalue))
h3=pd.DataFrame(res);h3.to_csv("results/stage1/H3_tumour_vs_normal.csv",index=False);print(h3.to_string())
