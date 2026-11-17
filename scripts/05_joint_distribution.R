library(survey)
library(haven)
options(survey.lonely.psu='adjust')
root <- '.'
d <- readRDS(file.path(root,'outputs/harmonized/harmonized_minimal_analysis.rds'))
ps <- c('2017-March 2020','August 2021-August 2023')
d$period <- factor(d$period,levels=ps)
# Verify current-period examination eligibility and phlebotomy weights.
r <- as.data.frame(read_xpt('data/raw/HEPA_L.xpt'))
i <- which(d$period==ps[2]); w <- r$WTPH2YR[match(d$SEQN[i],r$SEQN)]
stopifnot(!any(d$base_eligible[i] %in% 1 & (is.na(w)|w<=0)))
stopifnot(all(abs(d$analysis_weight[i][!is.na(w)&w>0]-w[!is.na(w)&w>0])<1e-6))
d$dose <- ifelse(d$IMQ011 %in% 3,'No doses',ifelse(d$IMQ011 %in% 2,'Fewer than 2 excluding zero',ifelse(d$IMQ011 %in% 1,'At least 2','Unknown')))
d$serology <- ifelse(d$hav_susceptible %in% 1,'Negative',ifelse(d$hav_susceptible %in% 0,'Positive',NA))
d$birthplace <- factor(ifelse(d$DMDBORN4 %in% 1,'US-born',ifelse(d$DMDBORN4 %in% 2,'Foreign-born',NA)),levels=c('US-born','Foreign-born'))
d$race <- factor(ifelse(d$RIDRETH3 %in% c(1,2),'Hispanic',ifelse(d$RIDRETH3 %in% 3,'Non-Hispanic White',ifelse(d$RIDRETH3 %in% 4,'Non-Hispanic Black',ifelse(d$RIDRETH3 %in% 6,'Non-Hispanic Asian','Other/multiracial')))))
d$diagnosed <- ifelse(d$DIQ010 %in% c(1,2,3),as.numeric(d$DIQ010==1),NA)
d$age10 <- d$RIDAGEYR/10
late <- 'data/raw'
for (m in c('HUQ_L','HIQ_L')) {
 x <- as.data.frame(read_xpt(file.path(late,paste0(m,'.xpt'))))
 for(v in intersect(c('HUQ030','HIQ011'),names(x))) {d[[v]]<-NA_real_;d[[v]][i]<-x[[v]][match(d$SEQN[i],x$SEQN)]}
}
d$usual_care <- ifelse(d$HUQ030 %in% c(1,2,3),as.numeric(d$HUQ030 %in% c(1,3)),NA)
d$insured <- ifelse(d$HIQ011 %in% c(1,2),as.numeric(d$HIQ011==1),NA)
d$stratum <- interaction(d$period,d$SDMVSTRA,drop=TRUE)
d$psu <- interaction(d$period,d$SDMVSTRA,d$SDMVPSU,drop=TRUE)
u <- subset(d,adult20 %in% 1 & !is.na(analysis_weight) & analysis_weight>0)
full <- svydesign(ids=~psu,strata=~stratum,weights=~analysis_weight,nest=TRUE,data=u)
target <- subset(full,base_eligible %in% 1 & target_masld_diabetes %in% 1)
stopifnot(identical(as.integer(table(target$variables$period)),c(755L,383L)))
out <- 'outputs/joint';dir.create(out,recursive=TRUE,showWarnings=FALSE)
source("scripts/survey_proportions.R")
rows<-birth<-selection<-service<-list()
for(p in ps) {
 period_value<-p
 a<-subset(target,period==p & hav_valid %in% 1)
 for(s in c('Negative','Positive')) for(v in c('No doses','Fewer than 2 excluding zero','At least 2','Unknown')) {
  rows[[length(rows)+1]]<-prop(a,a$variables$serology==s&a$variables$dose==v,paste(s,v,sep=' / '),p)
 }
 # Additional marginal and conditional measures make denominators explicit.
 for(b in c('US-born','Foreign-born')) {
  z<-subset(a,birthplace==b)
  specs<-list('Anti-HAV negative'=z$variables$hav_susceptible,'Unknown history'=as.numeric(z$variables$dose=='Unknown'))
  for(k in names(specs)) birth[[length(birth)+1]]<-prop(z,specs[[k]],k,p,b)
  z<-subset(z,dose!='Unknown')
  birth[[length(birth)+1]]<-prop(z,as.numeric(z$variables$dose=='At least 2'),'At least 2 reported doses among known histories',p,b)
  birth[[length(birth)+1]]<-prop(z,z$variables$hav_joint_indicator,'Negative and fewer than 2 among known histories',p,b)
  z<-subset(a,birthplace==b)
  for(s in c('Negative','Positive')) for(v in c('No doses','Fewer than 2 excluding zero','At least 2','Unknown')) birth[[length(birth)+1]]<-prop(z,as.numeric(z$variables$serology==s & z$variables$dose==v),paste(s,v,sep=' / '),p,b)
 }
 # Expanded CAP/diabetes domain describes a different phenotype; not a bias correction.
 base<-subset(full,period==p & vcte_reliable %in% 1 & diabetes_broad %in% 1 & LUXCAPM>=285 & hav_valid %in% 1)
 for(g in c('All alcohol statuses','Unclassifiable','Below threshold')) {
  z<-base
  if(g=='Unclassifiable') z<-subset(base,alcohol_classifiable %in% 0)
  if(g=='Below threshold') z<-subset(base,alcohol_masld_range %in% 1)
  for(s in c('Negative','Positive')) for(v in c('No doses','Fewer than 2 excluding zero','At least 2','Unknown')) selection[[length(selection)+1]]<-prop(z,as.numeric(z$variables$serology==s & z$variables$dose==v),paste(s,v,sep=' / '),p,g)
 }
 # Physician-diagnosed diabetes available in both periods; service variables only later.
 for(v in c('diagnosed',if(p==ps[2]) c('usual_care','insured'))) for(level in 1:0) {
  z<-a[a$variables[[v]] %in% level,]
  if(nrow(z$variables)==0) next
  service[[length(service)+1]]<-prop(z,as.numeric(z$variables$dose=='Unknown'),'Unknown history',p,paste(v,level,sep='='))
  zz<-subset(z,dose!='Unknown')
  service[[length(service)+1]]<-prop(zz,zz$variables$hav_joint_indicator,'Negative and fewer than 2 among known histories',p,paste(v,level,sep='='))
 }
}
write.csv(do.call(rbind,rows),file.path(out,'joint_distribution.csv'),row.names=FALSE)
write.csv(do.call(rbind,birth),file.path(out,'birthplace_components.csv'),row.names=FALSE)
write.csv(do.call(rbind,selection),file.path(out,'alcohol_joint_distribution.csv'),row.names=FALSE)
write.csv(do.call(rbind,service),file.path(out,'service_descriptions.csv'),row.names=FALSE)
sparse<-with(subset(target$variables,hav_valid %in% 1),as.data.frame(table(period,birthplace,race,dose=='Unknown')))
write.csv(sparse,file.path(out,'sparse_birth_race.csv'),row.names=FALSE)
# Period-specific simplified models avoid asserting one pooled national prevalence.
models<-list()
for(p in ps) {
 a<-subset(target,period==p & !is.na(hav_joint_indicator) & !is.na(birthplace) & !is.na(race))
 for(m in c('Birthplace only','Age and birthplace','Age birthplace and race')) {
  form<-switch(m,'Birthplace only'=hav_joint_indicator~birthplace,'Age and birthplace'=hav_joint_indicator~age10+birthplace,'Age birthplace and race'=hav_joint_indicator~age10+birthplace+race)
  fit<-svyglm(form,design=a,family=quasibinomial())
  stopifnot(fit$converged,all(is.finite(coef(fit))))
  b<-coef(fit)['birthplaceForeign-born'];se<-sqrt(vcov(fit)['birthplaceForeign-born','birthplaceForeign-born']);cr<-qt(.975,fit$df.residual)
  models[[length(models)+1]]<-data.frame(period=p,model=m,n=nrow(a$variables),events=sum(a$variables$hav_joint_indicator),OR=exp(b),ci_low=exp(b-cr*se),ci_high=exp(b+cr*se),df=fit$df.residual,converged=fit$converged)
 }
}
write.csv(do.call(rbind,models),file.path(out,'simplified_models.csv'),row.names=FALSE)
saveRDS(list(data=d,full=full,target=target),file.path(out,'analysis_designs.rds'))
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
print(do.call(rbind,rows));print(do.call(rbind,models));print(do.call(rbind,service))
