# Step10 E1 (v0.6a): community remodelling after HP removal. EXPLORATORY. Decisions: freeze_log (E1 entry).
suppressMessages({library(vegan);library(phyloseq);library(ANCOMBC);library(ALDEx2);library(Maaslin2);library(metafor);library(jsonlite);library(permute)})
B<-'/REPLACE_WITH_PRIVATE_WORKDIR/'; O<-paste0(B,'results/step10/E1/'); SEED<-20261030
ts<-function(...) cat(format(Sys.time(),'%F %T'),...,'\n')
ORAL<-c('Streptococcus','Veillonella','Prevotella','Fusobacterium','Peptostreptococcus','Parvimonas','Gemella','Porphyromonas','Haemophilus','Neisseria')
load_cohort<-function(coh,db){
  g<-read.delim(paste0(B,'results/',coh,'/genus_',db,'.tsv'),row.names=1,check.names=FALSE)
  m<-read.csv(paste0(B,'results/',coh,'/per_sample_',db,'.csv'),row.names=1)
  if(!'type'%in%names(m)) m$type<-m$type_bs
  m$pid<-as.character(m$pid); m<-m[rownames(m)%in%rownames(g),]
  tot<-rowSums(g[rownames(m),]); m$total<-tot
  pp<-split(rownames(m),m$pid); keep<-unlist(pp[sapply(pp,function(v) length(v)==2 && setequal(m[v,'type'],c('N','T')) && min(tot[v])>=5000)])
  m<-m[keep,]; m<-m[order(m$pid,m$type),]; g<-as.matrix(g[rownames(m),]); list(g=g,m=m)}
prep<-function(g,dropHP=TRUE){
  if(dropHP) g<-g[,!grepl('g__Helicobacter',colnames(g)),drop=FALSE]
  g<-g[,colMeans(g>0)>=0.10,drop=FALSE]; g}
rclr<-function(g){l<-log(g); l[!is.finite(l)]<-NA; mu<-rowMeans(l,na.rm=TRUE); r<-sweep(l,1,mu); r[is.na(r)]<-0; r}
clr05<-function(g){l<-log(g+0.5); sweep(l,1,rowMeans(l))}
e1a<-function(X,m,lab){
  set.seed(SEED); d<-dist(X); ctrl<-how(nperm=9999,blocks=factor(m$pid))
  a<-adonis2(d~type,data=m,permutations=ctrl)
  set.seed(SEED); bd<-betadisper(d,factor(m$type)); pt<-vegan:::permutest.betadisper(bd,permutations=9999)
  list(label=lab,R2=a$R2[1],F=a$F[1],P=a$`Pr(>F)`[1],betadisper_P=pt$tab$`Pr(>F)`[1],mean_dist_to_centroid=as.list(tapply(bd$distances,m$type,mean)))}
res<-list(); ab<-list()
run_cohort<-function(coh,db,full=TRUE){
  ts('cohort',coh,db); L<-load_cohort(coh,db); m<-L$m; g0<-L$g
  gH<-prep(g0,TRUE); gW<-prep(g0,FALSE)
  out<-list(n_pairs=nrow(m)/2,n_genera_noHP=ncol(gH),n_genera_withHP=ncol(gW),nonHP_depth_median=median(rowSums(gH)),nonHP_depth_min=min(rowSums(gH)),
            n_samples_nonHP_depth_lt1000=sum(rowSums(gH)<1000))
  out$E1a<-list(noHP_rCLR=e1a(rclr(gH),m,'HP removed, rCLR (primary)'),noHP_CLR05=e1a(clr05(gH),m,'HP removed, CLR pc0.5'),
               withHP_rCLR=e1a(rclr(gW),m,'with HP, rCLR (comparison)'))
  ts(coh,db,'E1a done')
  # E1c
  dl<-read.csv(paste0(B,'results/',coh,'/per_sample_',db,'.csv'),row.names=1); 
  R<-rclr(gH); nm<-sub('.*;g__','',colnames(R)); idT<-rownames(m)[m$type=='T']; idN<-rownames(m)[m$type=='N']
  stopifnot(all(m[idT,'pid']==m[idN,'pid']))
  dHP<-dl[idT,'clr']-dl[idN,'clr']; e1c<-list()
  for(o in ORAL){ j<-which(nm==o & grepl(';g__',colnames(R)))
    if(length(j)!=1){e1c[[o]]<-list(status=if(length(j)==0)'absent after filter' else 'multiple exact matches'); next}
    dg<-R[idT,j]-R[idN,j]; ct<-suppressWarnings(cor.test(dHP,dg,method='spearman',exact=FALSE)); e1c[[o]]<-list(rho=unname(ct$estimate),P=ct$p.value,n=length(dg))}
  ps<-sapply(e1c,function(x) if(is.null(x$P)) NA else x$P); q<-p.adjust(ps[!is.na(ps)],'BH'); for(k in names(q)) e1c[[k]]$q_BH<-unname(q[k])
  out$E1c<-e1c
  # E1b ANCOMBC2
  set.seed(SEED); ct<-t(gH); meta<-data.frame(row.names=rownames(m),tissue=factor(m$type,levels=c('N','T')),pid=m$pid)
  ps0<-phyloseq(otu_table(ct,taxa_are_rows=TRUE),sample_data(meta))
  t0<-Sys.time(); fit<-ancombc2(data=ps0,fix_formula='tissue',rand_formula='(1|pid)',p_adj_method='BH',prv_cut=0.10,n_cl=2,verbose=FALSE)
  rr<-fit$res; write.csv(rr,paste0(O,'E1b_ancombc2_',coh,'_',db,'.csv'),row.names=FALSE); ts(coh,db,'ANCOMBC2 done',format(Sys.time()-t0))
  out$E1b_ancombc2<-list(n_taxa=nrow(rr),n_q_lt_0.10=sum(rr$q_tissueT<0.10,na.rm=TRUE),n_up_T_q10=sum(rr$q_tissueT<0.10 & rr$lfc_tissueT>0,na.rm=TRUE),n_down_T_q10=sum(rr$q_tissueT<0.10 & rr$lfc_tissueT<0,na.rm=TRUE),
     n_q10_passed_ss=if('passed_ss_tissueT'%in%names(rr)) sum(rr$q_tissueT<0.10 & rr$passed_ss_tissueT,na.rm=TRUE) else NA)
  ab[[paste(coh,db)]]<<-rr
  if(full){
    # ALDEx2 paired
    set.seed(SEED); cond<-as.character(m$type); x<-aldex.clr(ct,cond,mc.samples=128,denom='all',verbose=FALSE)
    tt<-aldex.ttest(x,paired.test=TRUE); ef<-aldex.effect(x,verbose=FALSE); al<-cbind(tt,ef); write.csv(al,paste0(O,'E1b_aldex2_',coh,'_',db,'.csv'))
    out$E1b_aldex2<-list(n_taxa=nrow(al),n_wi_eBH_lt_0.10=sum(al$wi.eBH<0.10,na.rm=TRUE),n_we_eBH_lt_0.10=sum(al$we.eBH<0.10,na.rm=TRUE))
    ts(coh,db,'ALDEx2 done')
    md<-paste0(O,'maaslin2_',coh,'_',db); fitm<-Maaslin2(input_data=as.data.frame(gH),input_metadata=meta,output=md,fixed_effects='tissue',random_effects='pid',
       reference='tissue,N',normalization='TSS',transform='LOG',analysis_method='LM',min_prevalence=0.10,plot_heatmap=FALSE,plot_scatter=FALSE,cores=1)
    rs<-fitm$results; out$E1b_maaslin2<-list(n_taxa=nrow(rs),n_q_lt_0.10=sum(rs$qval<0.10,na.rm=TRUE)); ts(coh,db,'Maaslin2 done')
  }
  out}
for(cc in list(c('f2','silva'),c('korea16s','silva'),c('yu310127','silva'),c('f2','gg2v4'))){
  k<-paste(cc,collapse='_'); res[[k]]<-tryCatch(run_cohort(cc[1],cc[2],full=cc[2]=='silva'),error=function(e){ts('ERROR',k,conditionMessage(e));list(error=conditionMessage(e))})
  write_json(res,paste0(O,'E1_results.json'),auto_unbox=TRUE,digits=NA,pretty=TRUE)}
# meta over 3 SILVA cohorts
A<-lapply(c('f2 silva','korea16s silva','yu310127 silva'),function(k) ab[[k]]); names(A)<-c('f2','korea','yu310127')
if(all(!sapply(A,is.null))){
  common<-Reduce(intersect,lapply(A,function(a) a$taxon)); mt<-list()
  for(tx in common){ y<-sapply(A,function(a) a$lfc_tissueT[a$taxon==tx]); s<-sapply(A,function(a) a$se_tissueT[a$taxon==tx])
    f<-tryCatch(rma(yi=y,sei=s,method='REML'),error=function(e) NULL)
    mt[[tx]]<-data.frame(taxon=tx,lfc_f2=y[1],lfc_korea=y[2],lfc_yu=y[3],pooled=if(is.null(f)) NA else f$b[1],ci_lb=if(is.null(f)) NA else f$ci.lb,ci_ub=if(is.null(f)) NA else f$ci.ub,P=if(is.null(f)) NA else f$pval,I2=if(is.null(f)) NA else f$I2,n_same_sign=max(sum(y>0),sum(y<0)))}
  M<-do.call(rbind,mt); M$q_BH<-p.adjust(M$P,'BH'); M<-M[order(M$P),]; write.csv(M,paste0(O,'E1b_meta_REML_3cohorts.csv'),row.names=FALSE)
  res$E1b_meta<-list(n_common_genera=nrow(M),n_q_lt_0.10=sum(M$q_BH<0.10,na.rm=TRUE),n_all3_same_sign=sum(M$n_same_sign==3),
     top=head(M[,c('taxon','pooled','ci_lb','ci_ub','q_BH','n_same_sign')],15))}
write_json(res,paste0(O,'E1_results.json'),auto_unbox=TRUE,digits=NA,pretty=TRUE); ts('E1 DONE'); file.create(paste0(O,'DONE'))
