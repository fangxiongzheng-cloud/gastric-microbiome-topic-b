# A6 primary (TMB) + drivers + signatures. Seed 20261030.
suppressMessages({library(data.table);library(maftools);library(sandwich);library(lmtest)})
set.seed(20261030)
od<-"results/a6_multiomics"
man<-fread("data/a6/manifest_maf.csv")
tum<-man[sample_type=="Primary Tumor"]
hp<-fread("results/stage1/hp_calls.csv")[tissue=="T" & in380==TRUE & hp_consensus %in% c("pos","neg")]
cat("HP set:",table(hp$hp_consensus),"\n")
# read all MAFs
fl<-list.files("data/a6/maf",pattern="maf.gz$",full.names=TRUE)
fid<-sub("__.*","",basename(fl))
cols<-c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Tumor_Sample_Barcode","CONTEXT","t_depth","t_alt_count")
L<-lapply(seq_along(fl),function(i){d<-fread(cmd=paste("zcat",fl[i]),skip="Hugo_Symbol",select=cols,colClasses="character");d[,file_id:=fid[i]];d})
nrow_per<-sapply(L,nrow); cat("files:",length(fl)," zero-variant files:",sum(nrow_per==0),"\n")
M<-rbindlist(L)
# file -> tumour sample barcode (aliquot) from manifest / MAF
fs<-unique(tum[,.(file_id,case,sample,aliquots)])
fs<-fs[file_id %in% fid]
fs[,vial:=substr(sample,14,16)]
# tumour aliquot barcode: from MAF if any variants, else manifest aliquot starting with sample id
tb<-unique(M[,.(file_id,Tumor_Sample_Barcode)])
fs<-merge(fs,tb,by="file_id",all.x=TRUE)
fs[is.na(Tumor_Sample_Barcode),Tumor_Sample_Barcode:=mapply(function(al,s){a<-strsplit(al,";")[[1]];a<-a[startsWith(a,s)];if(length(a))sort(a)[1] else NA_character_},aliquots,sample)]
fs<-fs[substr(sample,14,15)=="01"]
setorder(fs,case,vial,Tumor_Sample_Barcode)
cat("cases with >1 tumour MAF:",sum(duplicated(fs$case)),"\n")
fs<-fs[!duplicated(case)]   # one per case, 01A preferred (sorted)
fs[,wxs_plate:=substr(Tumor_Sample_Barcode,22,25)]
M<-M[file_id %in% fs$file_id]
nonsyn<-c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Site","Translation_Start_Site","Nonstop_Mutation")
cnt<-M[Variant_Classification %in% nonsyn,.N,by=file_id]
fs<-merge(fs,cnt,by="file_id",all.x=TRUE); fs[is.na(N),N:=0L]
fs[,TMB:=N/38]; fs[,logTMB:=ifelse(N>0,log10(TMB),log10((N+0.5)/38))]
cat("zero nonsyn:",sum(fs$N==0),"\n")
cat("NA wxs plate:",sum(is.na(fs$wxs_plate)),"\n");D<-merge(hp[,.(case,hp=hp_consensus,subtype,rna_plate=plate)],fs[,.(case,file_id,Tumor_Sample_Barcode,wxs_plate,N,TMB,logTMB)],by="case")
D[,HP:=factor(hp,levels=c("neg","pos"))]
D[,subtype_u:=factor(ifelse(is.na(subtype)|subtype=="","Unknown",subtype))]
D[,rna_plate_u:=factor(ifelse(is.na(rna_plate)|rna_plate=="","Unknown",rna_plate))]
D[,wxs_plate:=factor(wxs_plate)]
fwrite(D,file.path(od,"a6_tmb_per_case.csv"))
cat("n with MAF:",table(D$HP),"; HP set without MAF:",nrow(hp)-nrow(D),"\n")
fit1<-function(f,dat,lab){m<-lm(f,data=dat);co<-coef(summary(m))["HPpos",];ci<-confint(m)["HPpos",]
  data.table(model=lab,n_pos=sum(dat$HP=="pos"),n_neg=sum(dat$HP=="neg"),beta_log10=co[1],ci_lo=ci[1],ci_hi=ci[2],P=co[4],GMR=10^co[1],GMR_lo=10^ci[1],GMR_hi=10^ci[2],df_resid=m$df.residual)}
res<-list()
res[[1]]<-fit1(logTMB~HP+subtype_u+wxs_plate,D,"PRIMARY: lm log10TMB ~ HP + subtype(Unknown level) + WXS plate")
res[[2]]<-fit1(logTMB~HP+subtype_u+rna_plate_u,D,"sens_a: RNA plate instead of WXS plate")
res[[3]]<-fit1(logTMB~HP+subtype+wxs_plate,D[!is.na(subtype)&subtype!=""],"sens_b: complete-case subtype")
m<-lm(logTMB~HP+subtype_u+wxs_plate,D); ct<-coeftest(m,vcov=vcovHC(m,"HC3"))["HPpos",]; cih<-ct[1]+c(-1,1)*qt(.975,m$df.residual)*ct[2]
res[[4]]<-data.table(model="sens_c: primary model, HC3 SE",n_pos=sum(D$HP=="pos"),n_neg=sum(D$HP=="neg"),beta_log10=ct[1],ci_lo=cih[1],ci_hi=cih[2],P=ct[4],GMR=10^ct[1],GMR_lo=10^cih[1],GMR_hi=10^cih[2],df_resid=m$df.residual)
w<-wilcox.test(TMB~HP,D,conf.int=TRUE)
res[[5]]<-data.table(model="sens_d: unadjusted Wilcoxon (TMB/Mb; estimate=HL shift neg-pos)",n_pos=sum(D$HP=="pos"),n_neg=sum(D$HP=="neg"),beta_log10=unname(w$estimate),ci_lo=w$conf.int[1],ci_hi=w$conf.int[2],P=w$p.value,GMR=NA,GMR_lo=NA,GMR_hi=NA,df_resid=NA)
R<-rbindlist(res); R[,n_in_model:=c(nobs(lm(logTMB~HP+subtype_u+wxs_plate,D)),nobs(lm(logTMB~HP+subtype_u+rna_plate_u,D)),NA,NA,NA)]; print(R); fwrite(R,file.path(od,"a6_primary_tmb.csv"))
desc<-D[,.(n=.N,median_TMB=median(TMB),IQR_lo=quantile(TMB,.25),IQR_hi=quantile(TMB,.75),geo_mean=10^mean(logTMB)),by=HP]; print(desc); fwrite(desc,file.path(od,"a6_tmb_descriptive.csv"))
fwrite(D[,.N,by=.(subtype_u,HP)][order(subtype_u,HP)],file.path(od,"a6_tmb_subtype_by_hp.csv"))
# ---- drivers
drv<-c("TP53","ARID1A","PIK3CA","RHOA","CDH1","KRAS","ERBB2","ERBB3","APC","SMAD4","FBXW7","CTNNB1","KMT2C","KMT2D","PTEN","RNF43","BCOR","ELF3","MUC6","SOX9","CDKN2A","BRAF","ARID2","TGFBR2","ZFHX3","PIK3R1","RASA1","ATM","NRAS")
Mn<-M[Variant_Classification %in% nonsyn & Hugo_Symbol %in% drv & file_id %in% D$file_id]
mut<-unique(Mn[,.(file_id,Hugo_Symbol)])
DR<-rbindlist(lapply(drv,function(g){x<-D$file_id %in% mut[Hugo_Symbol==g]$file_id
  fr<-mean(x); a<-sum(x&D$HP=="pos");b<-sum(D$HP=="pos")-a;c<-sum(x&D$HP=="neg");d<-sum(D$HP=="neg")-c
  if(fr<0.05) return(data.table(gene=g,freq_all=fr,mut_pos=a,n_pos=a+b,mut_neg=c,n_neg=c+d,OR=NA,OR_lo=NA,OR_hi=NA,P=NA,tested=FALSE))
  ft<-fisher.test(matrix(c(a,b,c,d),2,byrow=TRUE))
  data.table(gene=g,freq_all=fr,mut_pos=a,n_pos=a+b,mut_neg=c,n_neg=c+d,OR=unname(ft$estimate),OR_lo=ft$conf.int[1],OR_hi=ft$conf.int[2],P=ft$p.value,tested=TRUE)}))
DR[tested==TRUE,q_BH:=p.adjust(P,"BH")]; setorder(DR,-tested,P); print(DR); fwrite(DR,file.path(od,"a6_driver_fisher.csv"))
# ---- signatures (96-channel from CONTEXT)
S<-M[Variant_Type=="SNP" & file_id %in% D$file_id & nchar(CONTEXT)==11]
S[,ref:=substr(CONTEXT,6,6)]; S<-S[ref==Reference_Allele]
S[,alt:=Tumor_Seq_Allele2]; S[,tri:=substr(CONTEXT,5,7)]
rc<-function(x) chartr("ACGT","TGCA",sapply(lapply(strsplit(x,""),rev),paste,collapse=""))
pur<-S$ref %in% c("A","G"); S[pur,`:=`(tri=rc(tri),alt=chartr("ACGT","TGCA",alt))]; S[,ref:=substr(tri,2,2)]
S<-S[!grepl("N",tri)]
S[,ch:=paste0(substr(tri,1,1),"[",ref,">",alt,"]",substr(tri,3,3))]
suppressMessages(library(MutationalPatterns))
sig<-get_known_signatures(muttype="snv",source="COSMIC_v3.2",genome="GRCh38",incl_poss_artifacts=FALSE)
chs<-rownames(sig); if(is.null(chs)){g<-expand.grid(t=c("A","C","G","T"),f=c("A","C","G","T"),s=c("C>A","C>G","C>T","T>A","T>C","T>G"),stringsAsFactors=FALSE); chs<-paste0(g$f,"[",g$s,"]",g$t); rownames(sig)<-chs}  # MutationalPatterns standard order (substitution > 5' > 3')
cat("sig channels head:",head(chs,3),"\n")
mm<-dcast(S[,.N,by=.(file_id,ch)],ch~file_id,value.var="N",fill=0)
mat<-as.matrix(mm[,-1]); rownames(mat)<-mm$ch
mat<-mat[match(chs,rownames(mat)),,drop=FALSE]; mat[is.na(mat)]<-0; rownames(mat)<-chs
cat("SBS per sample summary:");print(summary(colSums(mat)))
keep<-colSums(mat)>=20; cat("samples with >=20 SBS:",sum(keep),"of",ncol(mat),"\n")
saveRDS(mat,file.path(od,"a6_sbs96_matrix.rds"))
ft<-fit_to_signatures_strict(mat[,keep],as.matrix(sig),max_delta=0.004,method="backwards")
co<-ft$fit_res$contribution; prop<-t(t(co)/colSums(co))
P<-data.table(file_id=colnames(prop),SBS1=prop["SBS1",],SBS5=prop["SBS5",],SBS17=prop["SBS17a",]+prop["SBS17b",],SBS18=prop["SBS18",])
P<-merge(P,D[,.(file_id,case,HP,subtype_u,wxs_plate)],by="file_id")
fwrite(merge(data.table(file_id=colnames(prop),t(prop)),D[,.(file_id,case,HP)],by="file_id"),file.path(od,"a6_signature_proportions_all.csv"))
SG<-rbindlist(lapply(c("SBS1","SBS5","SBS17","SBS18"),function(s){y<-P[[s]];w<-wilcox.test(y~P$HP,conf.int=TRUE)
  m<-lm(y~HP+subtype_u+wxs_plate,P);cc<-coef(summary(m))["HPpos",];ci<-confint(m)["HPpos",]
  data.table(signature=s,n_pos=sum(P$HP=="pos"),n_neg=sum(P$HP=="neg"),mean_pos=mean(y[P$HP=="pos"]),mean_neg=mean(y[P$HP=="neg"]),median_pos=median(y[P$HP=="pos"]),median_neg=median(y[P$HP=="neg"]),
   frac_nonzero_pos=mean(y[P$HP=="pos"]>0),frac_nonzero_neg=mean(y[P$HP=="neg"]>0),P_wilcox=w$p.value,
   adj_beta=cc[1],adj_lo=ci[1],adj_hi=ci[2],adj_P=cc[4])}))
SG[,q_BH:=p.adjust(P_wilcox,"BH")]; SG[,adj_q_BH:=p.adjust(adj_P,"BH")]; print(SG); fwrite(SG,file.path(od,"a6_signatures_tests.csv"))
# QC plot
pdf(file.path(od,"a6_qc_tmb.pdf"),7,4); par(mfrow=c(1,2))
boxplot(logTMB~HP,D,main="log10 TMB by consensus HP",ylab="log10(nonsyn/Mb)");stripchart(logTMB~HP,D,vertical=TRUE,add=TRUE,method="jitter",pch=16,cex=.4,col="grey40")
boxplot(logTMB~HP+subtype_u,D,las=2,cex.axis=.6,main="by subtype x HP");dev.off()
cat("DONE mutation\n")
