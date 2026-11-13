# Protocol v0.4 A3 sensitivity: SpiecEasi MB + StARS per group (same samples/features as SparCC). Seed 20261030.
# note: pulsar ncores>1 (mclapply) failed with "object sparseiCov of mode function was not found" in SpiecEasi 1.1.2 + pulsar 0.3.13; run single-core.
# Run with: /REPLACE_WITH_ENVIRONMENT_ROOT/envs/a3net/bin/Rscript scripts/a3_spieceasi.R
suppressMessages({library(SpiecEasi); library(igraph); library(jsonlite); library(Matrix)})
setwd('/REPLACE_WITH_PRIVATE_WORKDIR'); R='results/a3_network/'
X=as.matrix(read.csv(paste0(R,'input/counts_filtered.csv'),row.names=1,check.names=FALSE))
md=read.csv(paste0(R,'input/samples.csv'),row.names=1)[rownames(X),]
mod=c('Fusobacterium','Streptococcus','Prevotella','Peptostreptococcus')
out=list(); edges=list()
for(g in c('N','T')){
  set.seed(20261030)
  se=spiec.easi(X[md$type_bs==g,],method='mb',lambda.min.ratio=1e-2,nlambda=30,sel.criterion='stars',
     pulsar.params=list(rep.num=50,thresh=0.05,ncores=1,seed=20261030))
  B=as.matrix(symBeta(getOptBeta(se),mode='maxabs')); A=as.matrix(getRefit(se)); dimnames(B)=dimnames(A)=list(colnames(X),colnames(X))
  ut=which(upper.tri(A)&A!=0,arr.ind=TRUE)
  ed=data.frame(a=colnames(X)[ut[,1]],b=colnames(X)[ut[,2]],beta=B[ut]); ed=ed[order(-abs(ed$beta)),]
  write.csv(ed,paste0(R,'edges_spieceasi_mb_',g,'.csv'),row.names=FALSE); edges[[g]]=ed
  G=graph_from_adjacency_matrix(A,mode='undirected',diag=FALSE)
  pm=combn(mod,2); modb=apply(pm,2,function(p)B[p[1],p[2]])
  hp=sapply(mod,function(m)B['Helicobacter',m])
  out[[g]]=list(lambda_opt=getOptLambda(se),stars_stability=se$select$stars$summary[getOptInd(se)],n_edges=nrow(ed),
     pos_edge_frac=mean(ed$beta>0),density=edge_density(G),transitivity=transitivity(G,'global'),
     module_pairs=setNames(as.list(modb),apply(pm,2,paste,collapse='-')),module_edges_present=sum(modb!=0),module_mean_beta=mean(modb),
     helicobacter_degree=degree(G)['Helicobacter'][[1]],helicobacter_oral_beta=as.list(hp),
     helicobacter_neg_edges=sum(ed$beta<0 & (ed$a=='Helicobacter'|ed$b=='Helicobacter')),
     top_degree=as.list(sort(degree(G),decreasing=TRUE)[1:8]))
}
out$versions=list(R=R.version.string,SpiecEasi=as.character(packageVersion('SpiecEasi')),pulsar=as.character(packageVersion('pulsar')),huge=as.character(packageVersion('huge')))
write_json(out,paste0(R,'spieceasi_summary.json'),auto_unbox=TRUE,pretty=TRUE,digits=6); cat(toJSON(out,auto_unbox=TRUE,pretty=TRUE))
