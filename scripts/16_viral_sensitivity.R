# Assess sensitivity to exclusion of detected HBV and HCV markers.
library(survey)
library(haven)
options(survey.lonely.psu='adjust')
x<-readRDS('outputs/joint/analysis_designs.rds')
full<-x$full
z<-x$target
source('scripts/survey_proportions.R')
res<-flow<-list()
for(j in 1:2){
 period_value<-levels(full$variables$period)[j]
 f<-c('P_HEPC','HEPC_L')[j]
 y<-read_xpt(sprintf('https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/%s.xpt',c(2017,2021)[j],f))
 a<-subset(z,period==period_value & hav_valid %in% 1)
 h<-y$LBXHCR[match(a$variables$SEQN,y$SEQN)]
 hb<-read_xpt(sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/%s.xpt",c(2017,2021)[j],c("P_HEPBD","HEPBD_L")[j]))
 b<-hb$LBDHBG[match(a$variables$SEQN,hb$SEQN)]
 stopifnot(length(b)==nrow(a$variables),length(h)==length(b))
 positive<-b %in% 1 | h %in% 1
 unknown<-!positive & (!b %in% 2 | !h %in% c(2,3))
 flow[[j]]<-data.frame(period=period_value,n=length(b),HBsAg=sum(b %in% 1),RNA=sum(h %in% 1),both=sum(b %in% 1 & h %in% 1),unknown=sum(unknown))
 for(g in c('Baseline','Exclude detected positives','Classifiable negative only')){
  q<-switch(g,Baseline=a,'Exclude detected positives'=a[!positive,],'Classifiable negative only'=a[!positive & !unknown,])
  for(s in c('Negative','Positive')) for(v in c('No doses','Fewer than 2 excluding zero','At least 2','Unknown'))
   res[[length(res)+1]]<-prop(q,with(q$variables,serology==s & dose==v),paste(s,v,sep=' / '),period_value,g)
  q<-subset(q,dose!='Unknown')
  res[[length(res)+1]]<-prop(q,q$variables$hav_joint_indicator,'Secondary indicator (crude)',period_value,g)
 }
}
write.csv(do.call(rbind,flow),'outputs/joint/viral_flow.csv',row.names=FALSE)
write.csv(do.call(rbind,res),'outputs/joint/viral_sensitivity.csv',row.names=FALSE)
print(do.call(rbind,flow))



