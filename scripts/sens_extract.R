# Sensitivity §4 (v0.6): extract purity marker genes + clinical for TCGA; ESTIMATE (tidyestimate, RNA-seq mode) for TCGA and Korea.
suppressMessages({library(SummarizedExperiment);library(tidyestimate)})
B<-"/REPLACE_WITH_PRIVATE_WORKDIR/"; O<-paste0(B,"results/sensitivity/")
mk<-c("EPCAM","KRT8","KRT18","CDH1","PTPRC","COL1A1","COL1A2","COL3A1","DCN","LUM","ATP4A","MUC5AC","CDX2")
se<-readRDS(paste0(B,"data/host/star/stad_star_se.rds"))
rd<-rowData(se); tpm<-assay(se,"tpm_unstrand"); cd<-as.data.frame(colData(se))
cat(grep("stage|lauren|histolog|primary_diag",colnames(cd),value=TRUE,ignore.case=TRUE),"\n")
i<-match(mk,rd$gene_name); m<-tpm[i,,drop=FALSE]; rownames(m)<-mk
keep<-intersect(c("barcode","patient","sample_type","ajcc_pathologic_stage","primary_diagnosis","paper_Lauren.Class","paper_Lauren Class"),colnames(cd))
write.csv(cbind(cd[,keep,drop=FALSE],t(m)),paste0(O,"tcga_purity_markers_tpm.csv"),row.names=FALSE)
# ESTIMATE TCGA: log2(TPM+1), collapse duplicated symbols by max mean
L<-log2(tpm+1); sym<-rd$gene_name; ok<-!is.na(sym)&!duplicated(sym); L<-L[ok,]; rownames(L)<-sym[ok]
df<-as.data.frame(L,check.names=FALSE)  # note: run 1 passed a gene column too; tidyestimate scored it as a bogus sample "gene" which is dropped downstream
f<-filter_common_genes(df,id="hgnc_symbol",tidy=FALSE,tell_missing=TRUE,find_alias=TRUE)
es<-estimate_score(f,is_affymetrix=FALSE); write.csv(es,paste0(O,"tcga_estimate.csv"),row.names=FALSE)
# Korea
K<-read.csv(gzfile(paste0(B,"results/core/korea_gene_log2tpm.csv.gz")),check.names=FALSE)
colnames(K)[1]<-"gene"; K<-K[!duplicated(K$gene),]
f2<-filter_common_genes(K,id="hgnc_symbol",tidy=FALSE,tell_missing=TRUE,find_alias=TRUE)
ek<-estimate_score(f2,is_affymetrix=FALSE); write.csv(ek,paste0(O,"korea_estimate.csv"),row.names=FALSE)
cat("DONE\n")
