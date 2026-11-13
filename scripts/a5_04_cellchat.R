# A5 step 4: CellChat 2 per donor (see PREDECLARED_a5.md item 7). Resumable: skips donors with existing output.
suppressPackageStartupMessages({library(CellChat); library(Matrix)})
O <- "/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna"; W <- file.path(O,"work","donor"); R <- file.path(O,"cellchat"); dir.create(R,showWarnings=FALSE)
future::plan("sequential"); options(future.globals.maxSize=4e9)
donors <- sub(".mtx$","",list.files(W,pattern="mtx$"))
for (d in donors) {
  if (file.exists(file.path(R,paste0(d,"_weight.csv")))) next
  X <- as(readMM(file.path(W,paste0(d,".mtx"))),"CsparseMatrix"); g <- read.csv(file.path(W,paste0(d,"_genes.csv")))[[1]]
  m <- read.csv(file.path(W,paste0(d,"_meta.csv")),row.names=1); rownames(X) <- g; colnames(X) <- rownames(m)
  m$celltype <- factor(m$celltype)
  cc <- createCellChat(object=X, meta=m, group.by="celltype")
  cc@DB <- subsetDB(CellChatDB.human)   # default: excludes Non-protein Signaling
  cc <- subsetData(cc); cc <- identifyOverExpressedGenes(cc); cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type="triMean", population.size=FALSE, nboot=100, seed.use=20261030)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc); cc <- aggregateNet(cc, thresh=0.05)
  write.csv(cc@net$weight, file.path(R,paste0(d,"_weight_tmp.csv")))
  df <- subsetCommunication(cc, thresh=1); write.csv(df, file.path(R,paste0(d,"_lr_all.csv")), row.names=FALSE)
  saveRDS(cc, file.path(R,paste0(d,".rds")))
  file.rename(file.path(R,paste0(d,"_weight_tmp.csv")), file.path(R,paste0(d,"_weight.csv")))
  cat("cellchat", d, "done\n"); rm(cc,X); gc()
}
cat("STEP4_DONE", as.character(packageVersion("CellChat")), "\n")
