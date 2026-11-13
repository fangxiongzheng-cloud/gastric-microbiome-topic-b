# Step10 E5b Cox (v0.6a E5b as approved in v0.6c A). Supplementary only. Exploratory.
suppressMessages({library(survival);library(jsonlite)})
O<-'/REPLACE_WITH_PRIVATE_WORKDIR/results/step10/E5/'
d<-read.csv(paste0(O,'E5_analysis_table.csv'))
d<-d[complete.cases(d[,c('os_time','os_event','clr_N','dHP_TminusN','stage34','age','male')]),]
f<-coxph(Surv(os_time,os_event)~clr_N+dHP_TminusN+stage34+age+male,data=d,ties='efron')
s<-summary(f); z<-cox.zph(f)
co<-s$coefficients; ci<-s$conf.int
q<-p.adjust(co[,5],'BH')
out<-list(label='E5b EXPLORATORY, supplementary only; Cox OS ~ normal HP CLR + dHP(T-N CLR) + TNM(III-IV) + age + sex(male)',
 n=nrow(d),events=sum(d$os_event),
 coef=lapply(setNames(rownames(co),rownames(co)),function(k) list(logHR=co[k,1],HR=co[k,2],HR_CI95=c(ci[k,3],ci[k,4]),P=co[k,5],q_BH=unname(q[k]))),
 concordance=unname(s$concordance[1]),
 schoenfeld=list(global_P=z$table['GLOBAL','p'],per_term=as.list(z$table[-nrow(z$table),'p'])),
 software=paste('survival',packageVersion('survival')))
write_json(out,paste0(O,'E5b_results.json'),auto_unbox=TRUE,digits=NA,pretty=TRUE); cat(toJSON(out,auto_unbox=TRUE,pretty=TRUE,digits=6))
