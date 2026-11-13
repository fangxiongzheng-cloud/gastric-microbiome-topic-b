suppressMessages({library(GSVA);library(jsonlite);library(data.table)})
B<-'/REPLACE_WITH_PRIVATE_WORKDIR/'
D<-read.csv(paste0(B,'results/step10/A4full/korea_normal_scores_full.csv'),row.names=1)
E<-as.data.frame(fread(paste0(B,'results/core/korea_gene_log2tpm.csv.gz'))); rownames(E)<-E[[1]]; E<-as.matrix(E[,rownames(D)])
E<-E[apply(E,1,var)>0,]
IM<-c('CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4'); cat('IM genes present',sum(IM%in%rownames(E)),'\n')
sc<-gsva(ssgseaParam(E,list(IM=IM)),verbose=FALSE)
write.csv(data.frame(lib=colnames(sc),IM_ssgsea=as.numeric(sc[1,])),paste0(B,'results/step10/A4full/IM_ssgsea_normals.csv'),row.names=FALSE)
