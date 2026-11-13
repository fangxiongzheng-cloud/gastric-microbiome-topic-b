suppressMessages(library(SummarizedExperiment))
se<-readRDS("data/host/star/stad_star_se.rds"); print(dim(se)); print(assayNames(se))
g<-c("CDX1","CDX2","MUC2","TFF3","VIL1","CDH17","REG4","OLFM4","SI","FABP1","GKN1","GKN2","TFF2","PGC","PGA3","LIPF","CHIA","MUC6","CLDN18","GAST","ATP4A","ATP4B","CBLIF","KCNE2","KCNJ15","KCNJ16","SLC26A7","CCKBR","ESRRG","VSIG1","MUC5AC","TFF1","MUC1","FUT2","FUT3","B3GNT6")
rd<-rowData(se); i<-match(g,rd$gene_name); print(g[is.na(i)])
# PGA3 may be absent/ambiguous; report
tpm<-assay(se,"tpm_unstrand")[i[!is.na(i)],,drop=FALSE]; rownames(tpm)<-g[!is.na(i)]
cd<-as.data.frame(colData(se))[,c("barcode","patient","sample_type")]
out<-cbind(cd,t(tpm)); write.csv(out,"results/core/tcga_stad_hm_genes_tpm.csv",row.names=FALSE)
print(table(cd$sample_type)); print(colnames(out))
tot<-colSums(assay(se,"unstranded")); write.csv(data.frame(barcode=cd$barcode,total_counts=tot),"results/core/tcga_total_counts.csv",row.names=FALSE)
