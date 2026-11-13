import pandas as pd,numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
from scipy.stats import spearmanr
plt.rcParams.update({"font.size":8,"pdf.fonttype":42,"font.family":"DejaVu Sans"})
df=pd.read_csv("results/stage1/hp_calls.csv",index_col=0);T=df[df.tissue=="T"]
# Fig1
f,ax=plt.subplots(1,4,figsize=(14,3.4))
src=pd.Series({"NH2022 RNA-seq\n(host expr.)":380,"Spectrum2023\nRNA Kraken2":int((T.sp_hp.notna()&T.in380).sum()),"TCMA v2 WGS":int((T.tcma_wgs_reads.notna()&T.in380).sum()),"Both":int((T.sp_hp.notna()&T.tcma_wgs_reads.notna()&T.in380).sum()),"PanCanAtlas\nsubtype":int((T.subtype.notna()&T.in380).sum())})
ax[0].barh(src.index[::-1],src.values[::-1],color="#4C72B0");[ax[0].text(v+3,i,str(v),va="center") for i,v in enumerate(src.values[::-1])];ax[0].set_xlabel("Tumour cases (of 380)");ax[0].set_title("a  Data overlap")
j=T[T.tcma_wgs_reads.notna()&T.sp_hp.notna()]
r=spearmanr(j.tcma_wgs_rpm,j.sp_hp);ax[1].scatter(np.log1p(j.tcma_wgs_rpm*1e3),np.log1p(j.sp_hp),s=10,alpha=.7,c="#C44E52");ax[1].set_xlabel("TCMA WGS H. pylori log1p(RPM×1e3)");ax[1].set_ylabel("Spectrum2023 RNA H. pylori log1p(reads)");ax[1].set_title(f"b  TCMA vs Spectrum (n={len(j)}, ρ={r[0]:.2f})")
jp=T[T.tcma_wgs_reads.notna()&T.poore_raw.notna()];r2=spearmanr(jp.tcma_wgs_rpm,jp.poore_raw)
ax[2].scatter(np.log1p(jp.tcma_wgs_rpm*1e3),np.log1p(jp.poore_raw),s=10,alpha=.7,c="grey");ax[2].set_xlabel("TCMA WGS H. pylori log1p(RPM×1e3)");ax[2].set_ylabel("Poore 2020 raw Helicobacter log1p");ax[2].set_title(f"c  TCMA vs Poore raw (n={len(jp)}, ρ={r2[0]:.2f})")
cs=T[T.in380].groupby(["hp_consensus_source","hp_consensus"]).size().unstack(fill_value=0)
cs.plot.bar(stacked=True,ax=ax[3],color={"pos":"#C44E52","neg":"#8C8C8C","discordant":"#DD8452"});ax[3].set_title("d  Consensus H. pylori calls (380)");ax[3].set_xlabel("");ax[3].tick_params(axis='x',rotation=0)
plt.tight_layout();plt.savefig("results/stage1/Fig1_draft.pdf");plt.close()
# Fig2
r=pd.read_csv("results/stage1/H1_subtype_tests.csv");h3=pd.read_csv("results/stage1/H3_tumour_vs_normal.csv")
f=plt.figure(figsize=(14,7));g=f.add_gridspec(2,3)
a=f.add_subplot(g[0,0]);Tc=T[T.hp_consensus.isin(["pos","neg"])&T.subtype.notna()]
pr=Tc.groupby("subtype").hp_consensus.apply(lambda s:(s=="pos").mean()*100).reindex(["CIN","GS","EBV","MSI","HM-SNV"]);n=Tc.subtype.value_counts().reindex(pr.index)
a.bar(pr.index,pr.values,color=["#4C72B0","#55A868","#8172B2","#C44E52","grey"]);[a.text(i,v+1,f"n={n.iloc[i]}",ha="center") for i,v in enumerate(pr.values)];a.set_ylabel("% H. pylori consensus-positive");a.set_title("a  H. pylori positivity by PanCanAtlas subtype")
b=f.add_subplot(g[0,1]);d=T[T.subtype.isin(["CIN","GS","EBV","MSI"])]
b.boxplot([np.log1p(d[d.subtype==s].sp_hp.dropna()) for s in ["CIN","GS","EBV","MSI"]],tick_labels=["CIN","GS","EBV","MSI"],showfliers=False);
for i,s in enumerate(["CIN","GS","EBV","MSI"]): y=np.log1p(d[d.subtype==s].sp_hp.dropna());b.scatter(np.random.normal(i+1,.07,len(y)),y,s=4,alpha=.4)
b.set_ylabel("Spectrum H. pylori log1p(reads)");b.set_title("b  H. pylori load by subtype")
c=f.add_subplot(g[0,2]);rr=r[r.contrast.isin(["EBV vs rest","MSI vs rest","EBV+MSI vs GS+CIN"])&r.hedges_g.notna()].reset_index(drop=True)
y=np.arange(len(rr));c.scatter(rr.hedges_g,y,c=["#C44E52" if p<.05 else "grey" for p in rr.p_perm_withinplate]);c.axvline(0,c="k",lw=.5);c.set_yticks(y);c.set_yticklabels([f"{m[:28]} | {k}" for m,k in zip(rr.measure,rr.contrast)],fontsize=6);c.set_xlabel("Hedges' g (red: within-plate perm P<0.05)");c.set_title("c  Robustness across measures/decontam")
dd=f.add_subplot(g[1,0]);N=df[df.tissue=="N"]
pr2=pd.concat([T.sp_hp.rename("T"),N.sp_hp.rename("N")],axis=1).dropna()
for _,row in pr2.iterrows(): dd.plot([0,1],np.log1p([row["T"],row.N]),c="#C44E52" if row.N>row["T"] else "grey",alpha=.6)
dd.set_xticks([0,1]);dd.set_xticklabels(["Tumour","Adj. normal"]);dd.set_ylabel("Spectrum H. pylori log1p(reads)");dd.set_title(f"d  Paired T vs N (n={len(pr2)}), P={h3.p_wilcoxon[0]:.1e}")
e=f.add_subplot(g[1,1]);pr3=pd.concat([T.tcma_wgs_rpm.rename("T"),N.tcma_wgs_rpm.rename("N")],axis=1).dropna()
for _,row in pr3.iterrows(): e.plot([0,1],np.log1p(1e3*np.array([row["T"],row.N])),c="#C44E52" if row.N>row["T"] else "grey",alpha=.6)
e.set_xticks([0,1]);e.set_xticklabels(["Tumour","Adj. normal"]);e.set_ylabel("TCMA WGS log1p(RPM×1e3)");e.set_title(f"e  Paired T vs N, TCMA (n={len(pr3)}), P={h3.p_wilcoxon[2]:.1e}")
nulls=pd.read_csv("results/stage1/H1_perm_null_MSI.csv") if __import__("os").path.exists("results/stage1/H1_perm_null_MSI.csv") else None
ff=f.add_subplot(g[1,2])
if nulls is not None: ff.hist(nulls.null,bins=40,color="grey");ff.axvline(nulls.obs[0],c="r");ff.set_title("f  Within-plate label permutation null\n(MSI vs rest, Spectrum raw)");ff.set_xlabel("Mean difference log1p")
plt.tight_layout();plt.savefig("results/stage1/Fig2_draft.pdf");plt.close();print("ok")
