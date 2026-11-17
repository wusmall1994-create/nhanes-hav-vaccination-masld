x<-readRDS('outputs/joint/analysis_designs.rds');library(survey);options(survey.lonely.psu='adjust');full<-x$full;target<-x$target;ps<-levels(x$data$period)
source("scripts/survey_proportions.R")
rows<-list()
for(p in ps){period_value<-p;a<-subset(target,period==p & hav_valid %in% 1)
 for(v in c('No doses','Fewer than 2 excluding zero','At least 2','Unknown')){z<-subset(a,dose==v);rows[[length(rows)+1]]<-prop(z,z$variables$hav_susceptible,paste('Negative conditional on',v),p,'Conditional')}
 for(k in c('negative_lt2','known_joint','unknown','negative','positive_lt2')){
  z<-a
  hit<-switch(k,negative_lt2=as.numeric(z$variables$serology=='Negative'&z$variables$dose %in% c('No doses','Fewer than 2 excluding zero')),unknown=as.numeric(z$variables$dose=='Unknown'),negative=z$variables$hav_susceptible,positive_lt2=as.numeric(z$variables$serology=='Positive'&z$variables$dose %in% c('No doses','Fewer than 2 excluding zero')),known_joint=z$variables$hav_joint_indicator)
  rows[[length(rows)+1]]<-prop(z,hit,k,p,'Summary')
 }
}
write.csv(do.call(rbind,rows),'outputs/joint/conditional_summary.csv',row.names=FALSE)
print(do.call(rbind,rows))
