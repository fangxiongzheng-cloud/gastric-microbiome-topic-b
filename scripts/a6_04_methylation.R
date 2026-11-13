# A6 methylation: CIMP, MLH1, limma DMP. Seed 20261030.
suppressMessages({library(data.table);library(limma)})
set.seed(20261030); od<-"results/a6_multiomics"
hp<-fread("results/stage1/hp_calls.csv")[tissue=="T" & in380==TRUE & hp_consensus %in% c("pos","neg")]
man<-fread("data/a6/manifest_meth450.csv")[sample_type=="Primary Tumor"]
fl<-list.files("data/a6/meth",pattern="level3betas.txt$",full.names=TRUE); fid<-sub("__.*","",basename(fl))
man<-man[file_id %in% fid & case %in% hp$case]
man[,aliq:=mapply(function(al,s){a<-strsplit(al,";")[[1]];a<-a[startsWith(a,s)];if(length(a))sort(a)[1] else NA_character_},aliquots,sample)]
man[,vial:=substr(sample,14,16)]; setorder(man,case,vial,aliq); cat("dup cases:",sum(duplicated(man$case)),"\n"); man<-man[!duplicated(case)]
man[,meth_plate:=substr(aliq,22,25)]
cat("reading",nrow(man),"beta files\n")
rd<-function(f){x<-fread(f,header=FALSE,col.names=c("probe","b"));setNames(x$b,x$probe)}
b1<-rd(fl[match(man$file_id[1],fid)]); probes<-names(b1)
B<-matrix(NA_real_,length(probes),nrow(man),dimnames=list(probes,man$case))
for(i in seq_len(nrow(man))){x<-rd(fl[match(man$file_id[i],fid)]);B[,i]<-x[probes]}
# purity: PanCanAtlas ABSOLUTE (TCGAbiolinks::Tumor.purity has no STAD rows)
ab<-fread("data/a6/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"); ab[,case:=substr(array,1,12)]; ab<-ab[substr(array,14,15)=="01" & !is.na(purity)][!duplicated(case)]
D<-merge(man[,.(case,file_id,aliq,meth_plate)],hp[,.(case,hp=hp_consensus,subtype,rna_plate=plate)],by="case")
D<-merge(D,ab[,.(case,purity)],by="case",all.x=TRUE)
D[,HP:=factor(hp,levels=c("neg","pos"))]; D[,subtype_u:=factor(ifelse(is.na(subtype)|subtype=="","Unknown",subtype))]; D[,meth_plate:=factor(meth_plate)]
B<-B[,D$case]
cat("n meth:",table(D$HP)," purity missing:",sum(is.na(D$purity)),"\n")
fwrite(D,file.path(od,"a6_meth_samples.csv"))
# annotation
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
an<-as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
an<-an[intersect(rownames(an),rownames(B)),]
genes_of<-function(g){ok<-sapply(seq_len(nrow(an)),function(i){} );NULL}
prom_probes<-function(g){gn<-strsplit(an$UCSC_RefGene_Name,";");gr<-strsplit(an$UCSC_RefGene_Group,";")
  hit<-mapply(function(n,r) any(n==g & r %in% c("TSS1500","TSS200","1stExon")),gn,gr)
  rownames(an)[hit & an$Relation_to_Island=="Island"]}
cimp_genes<-c("CACNA1G","IGF2","NEUROG1","RUNX3","SOCS1")
pp<-lapply(c(cimp_genes,"MLH1"),prom_probes); names(pp)<-c(cimp_genes,"MLH1")
print(sapply(pp,length)); saveRDS(pp,file.path(od,"a6_cimp_mlh1_probes.rds"))
G<-sapply(pp,function(p) colMeans(B[p,,drop=FALSE],na.rm=TRUE))
G<-data.table(case=rownames(G),G)
G[,n_meth:=rowSums(.SD>=0.3,na.rm=TRUE),.SDcols=cimp_genes]; G[,CIMP_high:=n_meth>=3]; G[,CIMP_score:=rowMeans(.SD,na.rm=TRUE),.SDcols=cimp_genes]; G[,MLH1_hyper:=MLH1>=0.3]
G<-merge(G,D[,.(case,HP,subtype_u)],by="case"); fwrite(G,file.path(od,"a6_cimp_mlh1_per_case.csv"))
ftab<-function(x,lab){a<-sum(x&G$HP=="pos");b<-sum(G$HP=="pos")-a;c<-sum(x&G$HP=="neg");d<-sum(G$HP=="neg")-c;f<-fisher.test(matrix(c(a,b,c,d),2,byrow=TRUE))
  data.table(test=lab,method="Fisher",pos=paste0(a,"/",a+b),neg=paste0(c,"/",c+d),estimate=unname(f$estimate),lo=f$conf.int[1],hi=f$conf.int[2],P=f$p.value)}
wt<-function(y,lab){w<-wilcox.test(y~G$HP,conf.int=TRUE);data.table(test=lab,method="Wilcoxon (HL shift neg-pos)",pos=sprintf("median %.3f",median(y[G$HP=="pos"])),neg=sprintf("median %.3f",median(y[G$HP=="neg"])),estimate=unname(w$estimate),lo=w$conf.int[1],hi=w$conf.int[2],P=w$p.value)}
CM<-rbind(ftab(G$CIMP_high,"CIMP-high (>=3/5)"),wt(G$CIMP_score,"CIMP score (mean beta 5 genes)"),wt(G$MLH1,"MLH1 promoter mean beta"),ftab(G$MLH1_hyper,"MLH1 hypermethylated (beta>=0.3)"))
CM[,q_BH:=p.adjust(P,"BH")]; print(CM); fwrite(CM,file.path(od,"a6_cimp_mlh1_tests.csv"))
# subtype-adjusted sensitivity (logistic/linear, subtype_u)
sens<-rbind(
 data.table(test="CIMP-high logistic adj subtype",t(coef(summary(glm(CIMP_high~HP+subtype_u,binomial,G)))["HPpos",])),
 data.table(test="MLH1 beta lm adj subtype",t(coef(summary(lm(MLH1~HP+subtype_u,G)))["HPpos",])),
 data.table(test="CIMP score lm adj subtype",t(coef(summary(lm(CIMP_score~HP+subtype_u,G)))["HPpos",])),use.names=FALSE)
print(sens); fwrite(sens,file.path(od,"a6_cimp_mlh1_sens_adj_subtype.csv"))
# DMP limma
an_all<-as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
sexp<-rownames(an_all)[an_all$chr %in% c("chrX","chrY")]; snpp<-rownames(an_all)[!is.na(an_all$CpG_rs)|!is.na(an_all$SBE_rs)]
keep<-rowMeans(is.na(B))<=0.2 & !(rownames(B) %in% c(sexp,snpp)) & grepl("^cg",rownames(B))
cat("probes kept:",sum(keep),"of",nrow(B),"\n")
Dc<-D[!is.na(purity)]; Bk<-B[keep,Dc$case]
Mv<-log2(pmin(pmax(Bk,0.001),0.999)/(1-pmin(pmax(Bk,0.001),0.999)))
Dc[,meth_plate:=droplevels(meth_plate)]
des<-model.matrix(~HP+subtype_u+purity+meth_plate,Dc); cat("design",dim(des),"rank",qr(des)$rank,"\n")
fit<-eBayes(lmFit(Mv,des)); tt<-topTable(fit,coef="HPpos",number=Inf,sort.by="none")
db<-rowMeans(Bk[,Dc$HP=="pos"],na.rm=TRUE)-rowMeans(Bk[,Dc$HP=="neg"],na.rm=TRUE)
tt$delta_beta<-db[rownames(tt)]; tt$probe<-rownames(tt)
tt$gene<-an_all[rownames(tt),"UCSC_RefGene_Name"]; tt$island<-an_all[rownames(tt),"Relation_to_Island"]
tt$sig<-tt$adj.P.Val<0.05 & abs(tt$delta_beta)>=0.1
setDT(tt); setorder(tt,P.Value)
cat("n DMP (q<0.05):",sum(tt$adj.P.Val<0.05)," with |dB|>=0.1:",sum(tt$sig)," hyper:",sum(tt$sig&tt$delta_beta>0)," hypo:",sum(tt$sig&tt$delta_beta<0),"\n")
cat("min q:",min(tt$adj.P.Val),"\n"); print(head(tt[,.(probe,gene,logFC,delta_beta,P.Value,adj.P.Val)],15))
fwrite(tt,file.path(od,"a6_dmp_limma_all.csv.gz")); fwrite(tt[adj.P.Val<0.05],file.path(od,"a6_dmp_q005.csv"))
fwrite(data.table(n_pos=sum(Dc$HP=="pos"),n_neg=sum(Dc$HP=="neg"),probes=sum(keep),n_q005=sum(tt$adj.P.Val<0.05),n_sig=sum(tt$sig),n_hyper=sum(tt$sig&tt$delta_beta>0),n_hypo=sum(tt$sig&tt$delta_beta<0),min_q=min(tt$adj.P.Val),lambda=median(qchisq(1-tt$P.Value,1))/qchisq(.5,1)),file.path(od,"a6_dmp_summary.csv"))
pdf(file.path(od,"a6_qc_meth.pdf"),8,4);par(mfrow=c(1,3));hist(tt$P.Value,50,main="DMP P-value histogram",xlab="P")
plot(tt$delta_beta,-log10(tt$P.Value),pch=".",xlab="delta beta (pos-neg)",ylab="-log10 P",main="volcano");abline(v=c(-.1,.1),lty=2)
boxplot(MLH1~HP,G,main="MLH1 promoter beta");dev.off()
cat("DONE methylation\n")
