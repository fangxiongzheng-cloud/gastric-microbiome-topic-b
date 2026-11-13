# A5 step 1c (deviation): CellRanger>=3-style cell calling = OrdMag cells (retain) + EmptyDrops rescue (lower=100, niters=10000, FDR<=0.01).
suppressPackageStartupMessages({library(DropletUtils); library(Matrix)})
T <- "/REPLACE_WITH_PRIVATE_WORKDIR/results/a5_scrna/work/ed"; om <- read.csv(file.path(T,"ordmag.csv")); res <- list()
for (s in 1:18) {
  out <- file.path(T,paste0("called_s",s,".csv"))
  if (!file.exists(out)) {
    con <- file(file.path(T,paste0("s",s,".bin")),"rb"); n <- file.info(file.path(T,paste0("s",s,".bin")))$size/4
    a <- matrix(readBin(con,"integer",n=n,size=4),ncol=3,byrow=TRUE); close(con)
    cols <- sort(unique(a[,2])); m <- sparseMatrix(i=a[,1]+1, j=match(a[,2],cols), x=a[,3], dims=c(36601,length(cols))); rm(a); gc()
    set.seed(20261030)
    ed <- emptyDrops(m, lower=100, niters=10000, retain=om$ordmag_thr[s], test.ambient=FALSE)
    tot <- colSums(m); ng <- colSums(m>0)
    cell <- !is.na(ed$FDR) & ed$FDR<=0.01 & ng>=200
    write.csv(data.frame(col=cols[cell], umi=tot[cell], ngenes=ng[cell], fdr=ed$FDR[cell]), out, row.names=FALSE)
    cat("sample",s,"tested",sum(!is.na(ed$FDR)),"called",sum(cell),"\n"); rm(m,ed); gc()
  }
}
cat("STEP1C_DONE\n")
