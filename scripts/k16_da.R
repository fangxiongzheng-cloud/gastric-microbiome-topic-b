suppressMessages({library(ANCOMBC);library(Maaslin2);library(ALDEx2);library(phyloseq);library(vegan);library(TreeSummarizedExperiment)})
set.seed(1); setwd('/REPLACE_WITH_PRIVATE_WORKDIR'); O='results/korea16s/'
ps=read.csv(paste0(O,'per_sample_silva.csv'),row.names=1); g2=read.csv(paste0(O,'per_sample_gg2v4.csv'),row.names=1)
keep=rownames(ps)[pmin(ps$reads,g2[rownames(ps),'reads'])>=1000]
g=as.matrix(read.delim(paste0(O,'genus_silva.tsv'),row.names=1,check.names=FALSE))[keep,]
md=ps[keep,c('pid','type_bs')]; md$type=factor(md$type_bs,levels=c('N','T'))
# pairs complete
tb=table(md$pid); md=md[md$pid %in% names(tb)[tb==2],]; g=g[rownames(md),]
nm=sub('.*;','',colnames(g)); nm=sub('^g__','',nm); nm=make.unique(ifelse(nm==''|grepl('uncultured|_ge$',nm), colnames(g),nm)); colnames(g)=nm
prev=colMeans(g>0); g=g[,prev>=0.10]; cat('pairs',nrow(md)/2,'genera tested',ncol(g),'\n')
write.csv(data.frame(n_pairs=nrow(md)/2,n_genera=ncol(g)),paste0(O,'da/da_input.csv'))
# ANCOM-BC2
tse=TreeSummarizedExperiment(assays=list(counts=t(g)),colData=S4Vectors::DataFrame(md))
a=ancombc2(data=tse,assay_name='counts',fix_formula='type',rand_formula='(1|pid)',p_adj_method='BH',prv_cut=0,lib_cut=0,struc_zero=FALSE,neg_lb=FALSE,n_cl=6,verbose=FALSE)
ra=a$res; A=data.frame(taxon=ra$taxon,anc_lfc=ra$lfc_typeT,anc_q=ra$q_typeT,anc_pass_ss=ra$passed_ss_typeT)
# MaAsLin2
mo=Maaslin2(input_data=as.data.frame(g),input_metadata=md[,c('pid','type')],output=paste0(O,'da/maaslin2'),fixed_effects='type',random_effects='pid',reference='type,N',normalization='TSS',transform='LOG',min_prevalence=0,plot_heatmap=FALSE,plot_scatter=FALSE,cores=6)
M=mo$results; M=data.frame(taxon=M$feature,mas_coef=M$coef,mas_q=M$qval)
M$taxon=colnames(g)[match(M$taxon,make.names(colnames(g)))]
# ALDEx2 paired: order by pid
o=order(md$pid,md$type); gg=g[o,]; mm=md[o,]
x=aldex.clr(t(gg),as.character(mm$type),mc.samples=128,denom='all',verbose=FALSE)
tt=aldex.ttest(x,paired.test=TRUE); ef=aldex.effect(x,verbose=FALSE)
L=data.frame(taxon=rownames(tt),ald_eff=ef[rownames(tt),'effect'],ald_q=tt$wi.eBH)
R=Reduce(function(a,b)merge(a,b,by='taxon',all=TRUE),list(A,M,L))
# direction: positive = higher in tumour (ALDEx2 effect sign: groups ordered N,T -> effect = T - N? check)
R$anc_sig=R$anc_q<0.05; R$mas_sig=R$mas_q<0.05; R$ald_sig=R$ald_q<0.05
R$n_sig=rowSums(R[,c('anc_sig','mas_sig','ald_sig')],na.rm=TRUE)
write.csv(R,paste0(O,'da/da_all_silva.csv'),row.names=FALSE); save(x,file=paste0(O,'da/aldex_x.rda'))
# diversity (rarefy 5000 for alpha/beta)
gr=g[rowSums(g)>=5000,]; set.seed(1); r=rrarefy(gr,5000); mr=md[rownames(r),]
pp=names(which(table(mr$pid)==2)); r=r[mr$pid%in%pp,]; mr=mr[rownames(r),]
al=data.frame(mr,shannon=diversity(r),obs=rowSums(r>0),simpson=diversity(r,'simpson'))
write.csv(al,paste0(O,'da/alpha.csv'))
w=reshape(al[,c('pid','type','shannon','obs','simpson')],idvar='pid',timevar='type',direction='wide')
out=lapply(c('shannon','obs','simpson'),function(v){t=wilcox.test(w[[paste0(v,'.T')]],w[[paste0(v,'.N')]],paired=TRUE);data.frame(metric=v,pairs=nrow(w),med_N=median(w[[paste0(v,'.N')]]),med_T=median(w[[paste0(v,'.T')]]),P=t$p.value)})
alpha_res=do.call(rbind,out); print(alpha_res)
# beta: Bray & Aitchison; PERMANOVA with permutations restricted within pid
perm=how(nperm=999,blocks=mr$pid)
bc=vegdist(r,'bray'); ai=vegdist(r+1,'robust.aitchison')
b1=adonis2(bc~type,data=mr,permutations=perm); b2=adonis2(ai~type,data=mr,permutations=perm)
# also excluding Helicobacter
rh=r[,!grepl('Helicobacter',colnames(r))]; b3=adonis2(vegdist(rh,'bray')~type,data=mr,permutations=perm)
d1=betadisper(bc,mr$type); pd=permutest(d1,permutations=999)
beta=data.frame(test=c('Bray','robust.Aitchison','Bray_noHelicobacter','Bray_dispersion'),R2=c(b1$R2[1],b2$R2[1],b3$R2[1],NA),P=c(b1$`Pr(>F)`[1],b2$`Pr(>F)`[1],b3$`Pr(>F)`[1],pd$tab$`Pr(>F)`[1]))
print(beta); write.csv(alpha_res,paste0(O,'da/alpha_tests.csv'),row.names=FALSE); write.csv(beta,paste0(O,'da/beta_tests.csv'),row.names=FALSE)
pc=cmdscale(bc,k=2,eig=TRUE); write.csv(data.frame(mr,PC1=pc$points[,1],PC2=pc$points[,2]),paste0(O,'da/pcoa_bray.csv'))
cat('PCoA var',round(pc$eig[1:2]/sum(pc$eig[pc$eig>0]),3),'\n')
