# Step10 E2 (v0.6a): HP-associated host response program in adjacent normals. EXPLORATORY.
suppressMessages({library(limma);library(edgeR);library(fgsea);library(singscore);library(jsonlite);library(SummarizedExperiment)})
B<-'/REPLACE_WITH_PRIVATE_WORKDIR/';O<-paste0(B,'results/step10/E2/'); set.seed(20261030)
K<-read.csv(paste0(B,'results/core/korea_rna_per_sample.csv')); K$pid<-as.character(K$pid)
w<-reshape(K[,c('pid','type','hp_log1p_rpm')],idvar='pid',timevar='type',direction='wide'); pids<-w$pid[complete.cases(w)]
Kn<-K[K$type=='N' & K$pid%in%pids,]; rownames(Kn)<-Kn$lib
cnt<-read.csv(paste0(B,'results/step10/E7/input/korea_bulk_counts.csv.gz'),row.names=1,check.names=FALSE)
E<-read.csv(paste0(B,'results/core/korea_gene_log2tpm.csv.gz'),row.names=1,check.names=FALSE)
zs<-function(X,g,cols){x<-as.matrix(X[g,cols,drop=FALSE]); z<-t(scale(t(x))); colMeans(z)}
Kn$epi<-zs(E,c('EPCAM','KRT8','KRT18','CDH1'),Kn$lib); Kn$imm<-as.numeric(scale(as.numeric(E['PTPRC',Kn$lib])))
Kn$logdepth<-log10(Kn$salmon_mapped)
des<-model.matrix(~hp_log1p_rpm+logdepth+epi+imm,data=Kn)
y<-DGEList(as.matrix(cnt[,Kn$lib])); keep<-filterByExpr(y,des); y<-y[keep,,keep.lib.sizes=FALSE]; y<-calcNormFactors(y)
v<-voom(y,des); f<-eBayes(lmFit(v,des)); tt<-topTable(f,coef='hp_log1p_rpm',n=Inf,sort.by='none'); write.csv(tt,gzfile(paste0(O,'E2a_limma_normals.csv.gz')))
gm<-strsplit(readLines('/REPLACE_WITH_REFERENCE_ROOT/msigdb/h.all.v2025.1.Hs.symbols.gmt'),'\t'); hs<-setNames(lapply(gm,function(x) x[-(1:2)]),sapply(gm,'[',1))
fz<-list(V05_IM=c('CDX1','CDX2','MUC2','TFF3','VIL1','CDH17','REG4','OLFM4','SI','FABP1'),V05_GASTRIC_DIFF=c('GKN1','GKN2','TFF2','PGC','PGA3','LIPF','CHIA','MUC6','CLDN18','GAST'),
 V05_OXYNTIC_M2=c('ATP4A','ATP4B','CBLIF','KCNE2','KCNJ15','KCNJ16','SLC26A7','CCKBR','ESRRG','VSIG1'),V05_MUCIN_M3=c('MUC5AC','TFF1','TFF2','MUC6','MUC1','FUT2','FUT3','B3GNT6'))
sets<-c(hs,fz); st<-setNames(tt$t,rownames(tt))
set.seed(20261030); g<-fgsea(pathways=sets,stats=st,minSize=5,maxSize=500,eps=0,nproc=1)
g<-g[order(g$pval),]; gg<-as.data.frame(g); gg$leadingEdge<-sapply(gg$leadingEdge,paste,collapse=';'); write.csv(gg,paste0(O,'E2a_fgsea.csv'),row.names=FALSE)
res<-list(label='E2 EXPLORATORY (v0.6a)',n_normals=nrow(Kn),n_genes_tested=nrow(tt),n_genes_FDR10=sum(tt$adj.P.Val<0.10),
  n_genes_FDR10_up=sum(tt$adj.P.Val<0.10&tt$logFC>0),n_genes_FDR10_down=sum(tt$adj.P.Val<0.10&tt$logFC<0),
  msigdb_hallmark='h.all.v2025.1.Hs.symbols.gmt',fgsea_version=as.character(packageVersion('fgsea')),n_sets=nrow(g),
  sets_padj_lt_0.10=gg[gg$padj<0.10,c('pathway','NES','pval','padj','size')],frozen_sets=gg[grepl('^V05',gg$pathway),c('pathway','NES','pval','padj','size')])
sig<-g[g$padj<0.10,]
singscore<-function(X,up,down=NULL){ r<-rankGenes(X); if(length(down)) simpleScore(r,upSet=intersect(up,rownames(X)),downSet=intersect(down,rownames(X)))$TotalScore else simpleScore(r,upSet=intersect(up,rownames(X)))$TotalScore }
if(nrow(sig)>0){
  up<-unique(unlist(sig$leadingEdge[sig$NES>0])); dn<-unique(unlist(sig$leadingEdge[sig$NES<0])); both<-intersect(up,dn); up<-setdiff(up,both); dn<-setdiff(dn,both)
  writeLines(c(paste0('UP\t',up),paste0('DOWN\t',dn)),paste0(O,'E2_program_leading_edge.tsv'))
  ge<-rownames(y); X<-as.matrix(E[intersect(ge,rownames(E)),]); sc<-setNames(singscore(X,up,dn),colnames(X))
  Kt<-K[K$pid%in%pids,]; Kt$score<-sc[Kt$lib]
  T_<-Kt[Kt$type=='T',]; c1<-cor.test(T_$score,T_$hp_log1p_rpm,method='spearman',exact=FALSE)
  wN<-Kt[Kt$type=='N',]; wT<-Kt[Kt$type=='T',]; wT<-wT[!duplicated(wT$pid),]; m<-merge(wN,wT,by='pid',suffixes=c('.N','.T'))
  c2<-cor.test(m$score.T-m$score.N,m$hp_log1p_rpm.T-m$hp_log1p_rpm.N,method='spearman',exact=FALSE)
  q<-p.adjust(c(c1$p.value,c2$p.value),'BH')
  res$E2b<-list(n_up=length(up),n_down=length(dn),tumour=list(n=nrow(T_),rho=unname(c1$estimate),P=c1$p.value,q_BH=q[1]),paired_delta=list(n=nrow(m),rho=unname(c2$estimate),P=c2$p.value,q_BH=q[2]),
    normals_score_vs_hp_rho_descriptive=unname(cor(wN$score,wN$hp_log1p_rpm,method='spearman')),note='tumours: the 174 paired T runs (4 unpaired T runs excluded)')
  write.csv(Kt[,c('lib','pid','type','hp_log1p_rpm','score')],paste0(O,'E2_program_scores_korea.csv'),row.names=FALSE)
  # E2c TCGA normals
  se<-readRDS(paste0(B,'data/host/star/stad_star_se.rds')); tn<-read.csv(paste0(B,'results/a4_atrophy/tcga_normal_scores.csv'),row.names=1); tn<-tn[!is.na(tn$hp),]
  cd<-as.data.frame(colData(se)); i<-match(rownames(tn),substr(cd$barcode,1,12)[cd$sample_type=='Solid Tissue Normal']); bc<-cd$barcode[cd$sample_type=='Solid Tissue Normal'][i]
  tp<-assay(se,'tpm_unstrand')[,match(bc,cd$barcode)]; tp<-rowsum(tp,rowData(se)$gene_name); Xt<-log2(tp[intersect(rownames(tp),ge),]+1)
  s2<-singscore(Xt,up,dn); c3<-cor.test(s2,tn$hp,method='spearman',exact=FALSE)
  res$E2c_TCGA_direction_only<-list(n=length(s2),rho=unname(c3$estimate),note='direction only; HP = A4 TCGA measure (tcga_normal_scores.csv hp)')
} else res$E2b<-'not run: no gene set with padj<0.10 (no HP program defined)'
write_json(res,paste0(O,'E2_results.json'),auto_unbox=TRUE,digits=NA,pretty=TRUE); cat(toJSON(res,auto_unbox=TRUE,pretty=TRUE,digits=4))
