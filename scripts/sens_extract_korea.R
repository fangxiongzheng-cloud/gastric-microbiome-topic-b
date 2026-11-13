suppressMessages(library(tidyestimate)); B<-"/REPLACE_WITH_PRIVATE_WORKDIR/"
K<-read.csv(gzfile(paste0(B,"results/core/korea_gene_log2tpm.csv.gz")),check.names=FALSE,row.names=NULL)
colnames(K)[1]<-"gene"; K<-K[!duplicated(K$gene),]; rownames(K)<-K$gene; K$gene<-NULL; cat(dim(K),"\n")
f2<-filter_common_genes(K,id="hgnc_symbol",tidy=FALSE,tell_missing=FALSE,find_alias=TRUE)
ek<-estimate_score(f2,is_affymetrix=FALSE); write.csv(ek,paste0(B,"results/sensitivity/korea_estimate.csv"),row.names=FALSE); cat("DONE\n")
