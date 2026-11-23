library(survey)
options(survey.lonely.psu='adjust')
args <- commandArgs(TRUE)
root <- if(length(args)) args[1] else '.'
x <- readRDS(file.path(root,'outputs/joint/analysis_designs.rds'))
full <- update(x$full, core=serology=='Negative' & dose=='No doses')
source(file.path(root,'scripts/survey_proportions.R'))
a <- b <- list()
for(period_value in levels(full$variables$period)) {
 p <- period_value
 z <- subset(full,period==p & vcte_reliable %in% 1 & diabetes_broad %in% 1 & hav_valid %in% 1)
 high <- subset(z,LUXCAPM>=285)
 for(g in c('Observed low alcohol','Lower scenario','Upper scenario','All alcohol statuses')) {
  zz <- switch(g,'Observed low alcohol'=subset(high,alcohol_masld_range %in% 1),
  'Lower scenario'=subset(high,alcohol_masld_range %in% 1 | (alcohol_classifiable %in% 0 & !core)),
  'Upper scenario'=subset(high,alcohol_masld_range %in% 1 | (alcohol_classifiable %in% 0 & core)),high)
  a[[length(a)+1]] <- prop(zz,zz$variables$core,'Negative / No doses',p,g)
 }
 low <- subset(z,LUXCAPM<285 & alcohol_masld_range %in% 1)
 for(s in c('Negative','Positive')) for(v in c('No doses','Fewer than 2 excluding zero','At least 2','Unknown'))
 b[[length(b)+1]] <- prop(low,with(low$variables,serology==s & dose==v),paste(s,v,sep=' / '),p,'CAP <285')
}
a <- do.call(rbind,a); b <- do.call(rbind,b)
stopifnot(identical(a$n[a$domain=='Observed low alcohol'],c(716L,376L)),all(abs(tapply(b$estimate,b$period,sum)-1)<1e-10))
write.csv(a,file.path(root,'outputs/joint/alcohol_eligibility_scenarios.csv'),row.names=FALSE)
write.csv(b,file.path(root,'outputs/joint/cap_low_joint_distribution.csv'),row.names=FALSE)
print(a)


