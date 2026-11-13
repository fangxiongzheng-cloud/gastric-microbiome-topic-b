suppressMessages({library(immunedeconv);library(SummarizedExperiment)})
B<-'/REPLACE_WITH_PRIVATE_WORKDIR/';O<-paste0(B,'results/step10/E3/')
E<-read.csv(paste0(B,'results/core/korea_gene_log2tpm.csv.gz'),row.names=1,check.names=FALSE); tpm<-2^as.matrix(E)-1
se<-readRDS(paste0(B,'data/host/star/stad_star_se.rds')); pid<-read.csv(paste0(B,'results/core/hm_deltas_tcga.csv'))$pid
cd<-as.data.frame(colData(se)); s16<-substr(cd$barcode,1,16)
sp<-read.csv(paste0(B,'data/external/spectrum2023/stad_helicobacter_per_sample.csv'),row.names=1)
pick<-function(p,ty){ids<-rownames(sp)[substr(rownames(sp),1,12)==p & substr(rownames(sp),14,15)==ty]; which(s16==substr(ids[1],1,16))[1]}
cols<-c(sapply(pid,pick,ty='11'),sapply(pid,pick,ty='01'))
tt<-assay(se,'tpm_unstrand')[,cols]; colnames(tt)<-c(paste0(pid,'_N'),paste0(pid,'_T')); tt<-rowsum(tt,rowData(se)$gene_name)
imm<-function(X,tag){
  ep<-as.data.frame(deconvolute_epic(X,tumor=TRUE,scale_mrna=TRUE)); 
  epi<-colSums(ep[c('Bcells','CD4_Tcells','CD8_Tcells','Macrophages','NKcells'),,drop=FALSE])
  mc<-deconvolute_mcp_counter(X); keep<-c('T cells','CD8 T cells','Cytotoxic lymphocytes','NK cells','B lineage','Monocytic lineage','Myeloid dendritic cells','Neutrophils')
  z<-t(scale(t(log2(mc[keep,,drop=FALSE]+1)))); 
  out<-data.frame(sample=colnames(X),EPIC_immune=epi[colnames(X)],MCP_immune=colMeans(z)[colnames(X)])
  write.csv(out,paste0(O,'immune_',tag,'.csv'),row.names=FALSE); print(head(out)); print(rownames(ep))}
imm(tpm,'korea'); imm(as.matrix(tt),'tcga')
