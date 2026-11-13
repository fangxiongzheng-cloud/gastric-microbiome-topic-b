suppressMessages({library(singscore);library(edgeR);library(SummarizedExperiment)})
B<-'/REPLACE_WITH_PRIVATE_WORKDIR/';O<-paste0(B,'results/step10/E3/')
gm<-strsplit(readLines('/REPLACE_WITH_REFERENCE_ROOT/msigdb/h.all.v2025.1.Hs.symbols.gmt'),'\t'); hs<-setNames(lapply(gm,function(x) x[-(1:2)]),sapply(gm,'[',1)); inf<-hs[['HALLMARK_INFLAMMATORY_RESPONSE']]
E<-read.csv(paste0(B,'results/core/korea_gene_log2tpm.csv.gz'),row.names=1,check.names=FALSE)
cnt<-read.csv(paste0(B,'results/step10/E7/input/korea_bulk_counts.csv.gz'),row.names=1,check.names=FALSE)
grp<-substr(colnames(cnt),nchar(colnames(cnt)),nchar(colnames(cnt))); keep<-filterByExpr(DGEList(as.matrix(cnt)),group=grp)
X<-as.matrix(E[intersect(rownames(cnt)[keep],rownames(E)),colnames(cnt)]); s<-simpleScore(rankGenes(X),upSet=intersect(inf,rownames(X)))$TotalScore
write.csv(data.frame(sample=colnames(X),inflam=s),paste0(O,'inflam_korea.csv'),row.names=FALSE)
se<-readRDS(paste0(B,'data/host/star/stad_star_se.rds')); pid<-read.csv(paste0(B,'results/core/hm_deltas_tcga.csv'))$pid
cd<-as.data.frame(colData(se)); s16<-substr(cd$barcode,1,16); sp<-read.csv(paste0(B,'data/external/spectrum2023/stad_helicobacter_per_sample.csv'),row.names=1)
pick<-function(p,ty){ids<-rownames(sp)[substr(rownames(sp),1,12)==p & substr(rownames(sp),14,15)==ty]; which(s16==substr(ids[1],1,16))[1]}
cols<-c(sapply(pid,pick,ty='11'),sapply(pid,pick,ty='01')); nm<-c(paste0(pid,'_N'),paste0(pid,'_T'))
ct<-rowsum(assay(se,'unstranded')[,cols],rowData(se)$gene_name); tp<-rowsum(assay(se,'tpm_unstrand')[,cols],rowData(se)$gene_name); colnames(ct)<-colnames(tp)<-nm
k2<-filterByExpr(DGEList(ct),group=rep(c('N','T'),each=length(pid))); Xt<-log2(tp[k2,]+1)
s2<-simpleScore(rankGenes(Xt),upSet=intersect(inf,rownames(Xt)))$TotalScore
write.csv(data.frame(sample=nm,inflam=s2),paste0(O,'inflam_tcga.csv'),row.names=FALSE)
write.csv(log2(tp+1)[intersect(rownames(tp),c('MKI67','TOP2A','CCNB1','BIRC5','CDK1','COL1A1','COL1A2','COL3A1','DCN','LUM','EPCAM','KRT8','KRT18','CDH1','PTPRC','MUC5AC','TFF1','TFF2','MUC6','MUC1','FUT2','FUT3','B3GNT6','ST3GAL4','A4GNT','CEACAM1','CEACAM5','CEACAM6')),],paste0(O,'tcga_selected_log2tpm.csv'))
cat('ok',nrow(X),nrow(Xt),'\n')
