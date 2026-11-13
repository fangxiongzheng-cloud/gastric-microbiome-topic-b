suppressMessages({library(SummarizedExperiment);library(GSVA);library(msigdbr);library(immunedeconv);library(DESeq2)})
se<-readRDS("data/host/star/stad_star_se.rds")
se<-se[,se$sample_type=="Primary Tumor"]; se<-se[rowData(se)$gene_type=="protein_coding",]
cnt<-assay(se,"unstranded"); case<-substr(colnames(se),1,12); keep<-!duplicated(case); cnt<-cnt[,keep]; colnames(cnt)<-case[keep]
sym<-rowData(se)$gene_name; 
vst<-vst(round(cnt)); rownames(vst)<-sym; vst<-vst[!duplicated(sym),]
tpm<-assay(se,"tpm_unstrand")[,keep]; colnames(tpm)<-case[keep]; rownames(tpm)<-sym; tpm<-tpm[!duplicated(sym),]
saveRDS(list(vst=vst,tpm=tpm),"results/stage1/host_expr.rds")
h<-msigdbr(species="Homo sapiens",category="H"); k<-msigdbr(species="Homo sapiens",category="C2",subcategory="CP:KEGG_LEGACY")
gs<-c(split(h$gene_symbol,h$gs_name)[c("HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_E2F_TARGETS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")],
 list(KEGG_HPYLORI_EPITHELIAL=k$gene_symbol[k$gs_name=="KEGG_EPITHELIAL_CELL_SIGNALING_IN_HELICOBACTER_PYLORI_INFECTION"],
      TCELL_INFLAMED_GEP=c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6","HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT")))
print(sapply(gs,length))
sc<-gsva(gsvaParam(vst,gs),verbose=FALSE)
mcp<-deconvolute(tpm,"mcp_counter"); m<-as.matrix(mcp[,-1]); rownames(m)<-paste0("MCP_",gsub(" ","_",mcp$cell_type))
ck<-vst[c("CTLA4","PDCD1","CD274","LAG3"),]
out<-rbind(sc,m,ck); write.csv(t(out),"results/stage1/host_scores.csv")
cat("done",ncol(out),"\n")
