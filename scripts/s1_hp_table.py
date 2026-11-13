import pandas as pd,numpy as np
sp=pd.read_csv("data/external/spectrum2023/stad_helicobacter_per_sample.csv",index_col=0)
s1=pd.read_excel("data/external/spectrum2023/spectrum.00904-23-s0002.xlsx",sheet_name="S1",index_col=0)
sp["bact_total"]=s1.drop(columns=["g__Homo"]).sum(1)
s4=pd.read_excel("data/external/spectrum2023/spectrum.00904-23-s0003.xlsx",sheet_name="S4",index_col=0)
sp["bact_total_decontam"]=s4.sum(0).reindex(sp.index)
sp["case"]=sp.index.str[:12];sp["st"]=sp.index.str[13:15]
man=pd.read_csv("stad_star_manifest.csv");man["s16"]=man.cases.str[:16];man["plate"]=man.cases.str[21:25]
pl=man.drop_duplicates("s16").set_index("s16").plate
sp["plate"]=pl.reindex(sp.index).values
sp["tss"]=sp.case.str[5:7]
# TCMA
def tcma(kind):
    T=pd.read_csv(f"data/tcma/v2/{kind}/sample/bacteria.unambiguous.decontam.tissue.sample.reads.txt",sep="\t",index_col=0)
    R=pd.read_csv(f"data/tcma/v2/{kind}/sample/bacteria.unambiguous.decontam.tissue.sample.rpm.txt",sep="\t",index_col=0)
    return T.loc[210] if 210 in T.index else T.loc[209], R.loc[210] if 210 in R.index else R.loc[209], T.sum(0)
w,wr,wt=tcma("WGS"); x,xr,xt=tcma("WXS")
md=pd.read_csv('data/tcma/v2/metadata/metadata.TCMA.sample.txt',sep='\t');stadc=set(md[md.project.astype(str).str.contains('STAD')].case_id.astype(str))|set(md[md.project.astype(str).str.contains('STAD')].bcr_sample_barcode.astype(str).str[:12])
def bycase(s,st): s=s[[c for c in s.index if c[13:15]==st and c[:12] in stadc]]; return pd.Series(s.values,index=[c[:12] for c in s.index]).groupby(level=0).max()
sub=pd.read_csv("data/host/stad_pancan_subtypes.csv");sub["case"]=sub["pan.samplesID"].str[:12]
sub=sub.drop_duplicates("case").set_index("case").Subtype_Selected.str.replace("GI.","")
c380=open("data/host/nh2022_stad_rnaseq_tumour_cases.txt").read().split()
# Poore raw
mp=pd.read_csv("data/host/poore_stad_uuid_barcode_map.csv",index_col=0)
k=pd.read_csv("/REPLACE_WITH_RIGHTS_CLEARED_LEGACY_INPUTS/poore2020/Kraken-TCGA-Raw-Data-17625-Samples.csv",index_col=0,usecols=lambda c:c.startswith("Unnamed") or c.endswith("g__Helicobacter"))
mp["hp"]=k.reindex(mp.index).iloc[:,0]
po=mp[(mp.sample_type=="Primary Tumor")&(mp.experimental_strategy=="RNA-Seq")].groupby("case_barcode").hp.mean()
rows=[]
for st,lab in [("01","T"),("11","N")]:
    t=sp[sp.st==st].groupby("case").agg(sp_hp=("hpylori_sp_raw","max"),sp_hp_genus_decontam=("helico_genus_decontam","max"),sp_bact=("bact_total","max"),plate=("plate","first"),tss=("tss","first"))
    d=t.join([bycase(w,st).rename("tcma_wgs_reads"),bycase(wr,st).rename("tcma_wgs_rpm"),bycase(x,st).rename("tcma_wxs_reads")],how="outer")
    if st=="01": d=d.join(po.rename("poore_raw"),how="left")
    d["tissue"]=lab; rows.append(d)
df=pd.concat(rows); df.index.name="case"
df["subtype"]=sub.reindex(df.index).values
df["in380"]=df.index.isin(c380)
def cons(r):
    a=r.tcma_wgs_reads; b=r.sp_hp
    if pd.notna(a) and pd.notna(b):
        if a>0 and b>=10: return "pos"
        if a==0 and b==0: return "neg"
        return "discordant"
    if pd.notna(b): return "pos" if b>=10 else "neg"
    if pd.notna(a): return "pos" if a>0 else "neg"
    return np.nan
df["hp_consensus"]=df.apply(cons,axis=1)
df["hp_consensus_source"]=np.where(df.tcma_wgs_reads.notna()&df.sp_hp.notna(),"both",np.where(df.sp_hp.notna(),"spectrum_only",np.where(df.tcma_wgs_reads.notna(),"tcma_only","none")))
df.to_csv("results/stage1/hp_calls.csv")
T=df[df.tissue=="T"]
print(T.groupby(["hp_consensus_source","hp_consensus"]).size())
print("T in380 with call:",T[T.in380].hp_consensus.notna().sum(), T[T.in380].hp_consensus.value_counts().to_dict())
print("subtype x call (T, non-discordant):");print(pd.crosstab(T.subtype,T.hp_consensus))
