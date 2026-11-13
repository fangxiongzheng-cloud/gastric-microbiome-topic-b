# A5 step 6: epithelial trajectory. slingshot (primary), monocle3 (sensitivity), condiments progressionTest. PREDECLARED_a5.md item 8b.
suppressPackageStartupMessages({library(slingshot); library(condiments); library(Matrix); library(SingleCellExperiment)})
set.seed(20261030)
O <- "/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna"; W <- file.path(O,"work","traj"); R <- file.path(O,"trajectory"); dir.create(R,showWarnings=FALSE)
ob <- read.csv(file.path(W,"obs.csv"),row.names=1)
ob <- ob[ob$epi_type %in% c("Pit","Neck","IM","Tumor"),]
ob$lab <- ob$epi_type; ob$lab[ob$epi_type=="Pit" & ob$stage=="GS"] <- "GS_Pit"
Z <- as.matrix(ob[,paste0("H",1:10)])
if (!file.exists(file.path(R,"slingshot.rds"))) {
  sds <- slingshot(Z, clusterLabels=ob$lab, start.clus="GS_Pit"); saveRDS(sds, file.path(R,"slingshot.rds")) } else sds <- readRDS(file.path(R,"slingshot.rds"))
pt <- slingPseudotime(sds); cw <- slingCurveWeights(sds)
lin <- slingLineages(sds); capture.output(print(lin), file=file.path(R,"slingshot_lineages.txt"))
out <- data.frame(cell=rownames(ob), ob[,c("stage","hp","donor","lab","UMAP1","UMAP2")], pt, check.names=FALSE)
colnames(out)[grep("^Lineage",colnames(out))] <- paste0("pt_",colnames(pt))
# condiments: progression of HP+ vs HP- along each lineage
pt0 <- pt; pt0[is.na(pt0)] <- 0
pg <- progressionTest(pt0, cellWeights=cw, conditions=factor(ob$hp), global=TRUE, lineages=TRUE)
pg$q_BH <- p.adjust(pg$p.value,"BH"); write.csv(pg, file.path(R,"condiments_progressionTest.csv"), row.names=FALSE); print(pg)
# monocle3 sensitivity (same cells, given UMAP, root = GS_Pit)
suppressPackageStartupMessages(library(monocle3))
X <- as(readMM(file.path(W,"counts_hvg.mtx")),"CsparseMatrix"); g <- read.csv(file.path(W,"genes.csv"))[[1]]
allob <- read.csv(file.path(W,"obs.csv"),row.names=1); colnames(X) <- rownames(allob); rownames(X) <- make.unique(g)
X <- X[, rownames(ob)]
cds <- new_cell_data_set(X, cell_metadata=ob, gene_metadata=data.frame(gene_short_name=rownames(X),row.names=rownames(X)))
reducedDims(cds)$UMAP <- as.matrix(ob[,c("UMAP1","UMAP2")])
cds <- cluster_cells(cds, reduction_method="UMAP", random_seed=20261030); cds <- learn_graph(cds, use_partition=FALSE)
cds <- order_cells(cds, root_cells=rownames(ob)[ob$lab=="GS_Pit"])
out$pt_monocle3 <- pseudotime(cds)[rownames(ob)]
slmax <- apply(pt,1,function(x) if (all(is.na(x))) NA else max(x,na.rm=TRUE))
ok <- is.finite(out$pt_monocle3) & is.finite(slmax)
ct <- cor.test(slmax[ok], out$pt_monocle3[ok], method="spearman", exact=FALSE)
write.csv(out, file.path(R,"pseudotime_cells.csv"), row.names=FALSE)
jsonlite::write_json(list(n_cells=nrow(ob), n_lineages=length(lin), lineages=lapply(lin,as.character), spearman_slingshot_monocle3=unname(ct$estimate), p=ct$p.value, n_finite=sum(ok),
  versions=list(slingshot=as.character(packageVersion("slingshot")), condiments=as.character(packageVersion("condiments")), monocle3=as.character(packageVersion("monocle3")), R=R.version.string)),
  file.path(R,"trajectory_summary.json"), auto_unbox=TRUE, pretty=TRUE)
cat("STEP6_DONE\n")
