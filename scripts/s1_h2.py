import pandas as pd,numpy as np,statsmodels.formula.api as smf
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
rng=np.random.default_rng(7)
hs=pd.read_csv("results/stage1/host_scores.csv",index_col=0)
df=pd.read_csv("results/stage1/hp_calls.csv",index_col=0);T=df[df.tissue=="T"]
sub=pd.read_csv("data/host/stad_pancan_subtypes.csv");sub["case"]=sub["pan.samplesID"].str[:12];pur=None
d=T.join(hs,how="inner");d=d[d.in380]
d["sp_log"]=np.log1p(d.sp_hp);d["sp_dec_log"]=np.log1p(d.sp_hp_genus_decontam);d["tcma_log"]=np.log1p(d.tcma_wgs_rpm*1e3);d["poore_log"]=np.log1p(d.poore_raw);d["cons"]=d.hp_consensus.map({"pos":1,"neg":0})
d["lbact"]=np.log1p(d.sp_bact);d["plate"]=d.plate.fillna("NA");d["subtype"]=d.subtype.fillna("NA")
feats=[c for c in hs.columns]
clean=lambda s:"".join(ch if ch.isalnum() else "_" for ch in s)
res=[]
for mcol in ["sp_log","sp_dec_log","cons","tcma_log","poore_log"]:
    for f in feats:
        x=d[[mcol,f,"subtype","plate","lbact"]].dropna().rename(columns={f:"y"})
        if len(x)<30: continue
        rho,p=spearmanr(x[mcol],x.y)
        try:
            fit=smf.ols(f"y ~ {mcol} + C(subtype) + lbact",data=x).fit();b=fit.params[mcol];pa=fit.pvalues[mcol]
        except Exception as e: b=pa=np.nan
        # sample-shuffle null for spearman
        nul=[spearmanr(rng.permutation(x[mcol].values),x.y)[0] for _ in range(200)]
        res.append(dict(measure=mcol,feature=f,n=len(x),rho=rho,p_spearman=p,beta_adj=b,p_adj=pa,p_shuffle=(np.sum(np.abs(nul)>=abs(rho))+1)/201))
r=pd.DataFrame(res)
r["q_adj"]=r.groupby("measure").p_adj.transform(lambda p:multipletests(p.fillna(1),method="fdr_bh")[1])
r.to_csv("results/stage1/H2_immunity_tests.csv",index=False)
pd.set_option("display.width",250)
print(r[r.feature.isin(["CTLA4","PDCD1","CD274","LAG3","KEGG_HPYLORI_EPITHELIAL","TCELL_INFLAMED_GEP","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INTERFERON_GAMMA_RESPONSE"])].round(4).to_string())
print("n significant q_adj<0.1 per measure:",r[r.q_adj<0.1].groupby("measure").size().to_dict())
