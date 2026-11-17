library(survey)
library(splines)
options(survey.lonely.psu="adjust")
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)>=1) args[[1]] else "."
outdir <- file.path(root,"outputs","masld_hav_hbv_exploratory")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d <- readRDS(file.path(root,"outputs/harmonized/harmonized_minimal_analysis.rds"))
d <- d[d$adult20 %in% 1 & !is.na(d$analysis_weight) & d$analysis_weight>0,]
d$period <- factor(d$period)
d$stack_stratum <- interaction(d$period,d$SDMVSTRA,drop=TRUE)
d$stack_psu <- interaction(d$period,d$SDMVSTRA,d$SDMVPSU,drop=TRUE)
d$age10 <- (d$RIDAGEYR-60)/10
d$bmi5 <- (d$BMXBMI-35)/5
d$hba1c <- d$LBXGH-7.5
d$birthplace <- factor(ifelse(d$DMDBORN4 %in% 1,"US-born",ifelse(d$DMDBORN4 %in% 2,"Foreign-born",NA)),levels=c("US-born","Foreign-born"))
d$sex <- factor(ifelse(d$RIAGENDR==1,"Male","Female"),levels=c("Male","Female"))
race <- c("1"="Hispanic","2"="Hispanic","3"="Non-Hispanic White","4"="Non-Hispanic Black","6"="Non-Hispanic Asian","7"="Other/multiracial")
d$race <- factor(unname(race[as.character(d$RIDRETH3)]),levels=c("Non-Hispanic White","Hispanic","Non-Hispanic Black","Non-Hispanic Asian","Other/multiracial"))
d$fibrosis <- factor(ifelse(d$LUXSMED>=8,"LSM >=8 kPa","LSM <8 kPa"),levels=c("LSM <8 kPa","LSM >=8 kPa"))
des <- svydesign(ids=~stack_psu,strata=~stack_stratum,weights=~analysis_weight,nest=TRUE,data=d)
target <- subset(des,base_eligible==1 & target_masld_diabetes==1)
covars <- c("age10","birthplace","sex","race","bmi5","hba1c","fibrosis","period")
target$variables$complete_model <- complete.cases(target$variables[,c("hav_joint_indicator",covars)])
cc <- subset(target,complete_model)
formula <- hav_joint_indicator ~ age10 + birthplace + sex + race + bmi5 + hba1c + fibrosis + period
fit <- svyglm(formula,design=cc,family=quasibinomial())
stopifnot(fit$converged,all(is.finite(coef(fit))),min(eigen(vcov(fit),symmetric=TRUE)$values)>0)
b<-coef(fit);se<-sqrt(diag(vcov(fit)));crit<-qt(.975,fit$df.residual)
or<-data.frame(term=names(b),estimate=exp(b),ci_low=exp(b-crit*se),ci_high=exp(b+crit*se))
write.csv(or,file.path(outdir,"adjusted_odds_ratios.csv"),row.names=FALSE)
print(or)
diag <- data.frame(n=nrow(cc$variables),events=sum(cc$variables$hav_joint_indicator==1),nonevents=sum(cc$variables$hav_joint_indicator==0),parameters=length(b),design_df=degf(des),residual_df=fit$df.residual,converged=fit$converged,condition_number=kappa(vcov(fit)))
write.csv(diag,file.path(outdir,"model_diagnostics.csv"),row.names=FALSE)
audit<-list()
for(p in levels(d$period)){
 z<-subset(target$variables,period==p)
 audit[[length(audit)+1]]<-data.frame(period=p,target=nrow(z),outcome_available=sum(!is.na(z$hav_joint_indicator)),complete_model=sum(z$complete_model),events=sum(z$hav_joint_indicator[z$complete_model]==1),missing_covariates=sum(!is.na(z$hav_joint_indicator)&!z$complete_model))
}
write.csv(do.call(rbind,audit),file.path(outdir,"complete_case_flow.csv"),row.names=FALSE)
print(do.call(rbind,audit))
# Counterfactual covariate margins over the same complete-case empirical population.
standard<-cc$variables
w<-standard$analysis_weight/sum(standard$analysis_weight)
margin<-function(variable,level){
 nd<-standard
 nd[[variable]]<-factor(level,levels=levels(standard[[variable]]))
 X<-model.matrix(delete.response(terms(fit)),nd,contrasts.arg=fit$contrasts,xlev=fit$xlevels)
 X<-X[,names(b),drop=FALSE]
 prob<-plogis(drop(X%*%b))
 list(estimate=sum(w*prob),gradient=colSums(X*as.numeric(w*prob*(1-prob))))
}
margins<-list()
for(v in c("birthplace","sex","race","fibrosis","period")){
 ref<-margin(v,levels(standard[[v]])[1])
 for(level in levels(standard[[v]])){
  m<-margin(v,level);g<-m$gradient
  sem<-sqrt(drop(t(g)%*%vcov(fit)%*%g))
  gd<-g-ref$gradient
  sed<-sqrt(max(0,drop(t(gd)%*%vcov(fit)%*%gd)))
  margins[[length(margins)+1]]<-data.frame(variable=v,level=level,estimate=m$estimate,ci_low=max(0,m$estimate-crit*sem),ci_high=min(1,m$estimate+crit*sem),difference=m$estimate-ref$estimate,diff_low=m$estimate-ref$estimate-crit*sed,diff_high=m$estimate-ref$estimate+crit*sed)
 }
}
write.csv(do.call(rbind,margins),file.path(outdir,"adjusted_prevalence_margins.csv"),row.names=FALSE)
# Check age functional form without selecting covariates by statistical significance.
spline_fit<-svyglm(update(formula,.~.-age10+ns(RIDAGEYR,df=3)),design=cc,family=quasibinomial())
stopifnot(spline_fit$converged)
sb<-coef(spline_fit);ss<-sqrt(diag(vcov(spline_fit)));sc<-qt(.975,spline_fit$df.residual)
write.csv(data.frame(term=names(sb),estimate=exp(sb),ci_low=exp(sb-sc*ss),ci_high=exp(sb+sc*ss)),file.path(outdir,"age_spline_sensitivity.csv"),row.names=FALSE)
# Outcome count audit and descriptive subgroup estimates, with period-specific df.
subgroups<-list()
for(p in levels(d$period)){
 pd<-subset(target,period==p)
 pdf<-degf(subset(des,period==p))
 for(v in c("birthplace","race")){
  for(level in levels(d[[v]])){
   sd<-pd[pd$variables[[v]] %in% level,]
   y<-sd$variables$hav_joint_indicator
   n<-sum(!is.na(y));events<-sum(y==1,na.rm=TRUE)
   if(n>0 && events>0 && events<n){
    pr<-svyciprop(~hav_joint_indicator,sd,method="logit",df=pdf,na.rm=TRUE)
    lim<-confint(pr)
    subgroups[[length(subgroups)+1]]<-data.frame(period=p,variable=v,level=level,n=n,events=events,nonevents=n-events,estimate=coef(pr)[1],ci_low=lim[1],ci_high=lim[2],flag=if(n<30 || events<10 || n-events<10)"Sparse cell" else "")
   } else if(n>0) {
    subgroups[[length(subgroups)+1]]<-data.frame(period=p,variable=v,level=level,n=n,events=events,nonevents=n-events,estimate=events/n,ci_low=NA_real_,ci_high=NA_real_,flag="Zero outcome cell; CI not estimable")
   }
  }
 }
}
write.csv(do.call(rbind,subgroups),file.path(outdir,"subgroup_prevalence.csv"),row.names=FALSE)
print(diag)
print(do.call(rbind,margins))
capture.output(sessionInfo(),file=file.path(outdir,"model_session_info.txt"))

# Primary comparator display uses the same joint age-birthplace standard as Table S8.
d$capgroup <- factor(ifelse(d$LUXCAPM>=285,"CAP >=285","CAP <285"),levels=c("CAP <285","CAP >=285"))
des <- update(des,capgroup=d$capgroup)
diabetes_design <- subset(des,base_eligible==1 & diabetes_broad==1 & !is.na(birthplace))
standard_cap <- diabetes_design$variables
standard_cap$standard_weight <- ave(standard_cap$analysis_weight,standard_cap$period,FUN=function(w)w/sum(w)/2)
cap_design <- subset(diabetes_design,alcohol_masld_range==1)
cap_rows <- list()
for(outcome in c("hav_joint_indicator","hav_susceptible","hbv_susceptible")){
 f <- svyglm(reformulate("period*capgroup + ns(RIDAGEYR,df=3) + birthplace",outcome),design=cap_design,family=quasibinomial())
 bc <- coef(f);vc <- vcov(f);tc <- qt(.975,f$df.residual)
 cap_margin<-function(p,g){
  nd<-standard_cap
  nd$period<-factor(p,levels=levels(d$period))
  nd$capgroup<-factor(g,levels=levels(d$capgroup))
  X<-model.matrix(delete.response(terms(f)),nd,contrasts.arg=f$contrasts,xlev=f$xlevels)[,names(bc),drop=FALSE]
  pr<-plogis(drop(X%*%bc));w<-nd$standard_weight
  list(estimate=sum(w*pr),gradient=colSums(X*as.numeric(w*pr*(1-pr))))
 }
 for(p in levels(d$period)){
  lo<-cap_margin(p,"CAP <285");hi<-cap_margin(p,"CAP >=285")
  for(g in levels(d$capgroup)){
   mm<-if(g=="CAP <285")lo else hi
   se<-sqrt(drop(t(mm$gradient)%*%vc%*%mm$gradient))
   z<-subset(cap_design$variables,period==p & capgroup==g)
   cap_rows[[length(cap_rows)+1]]<-data.frame(period=p,outcome=outcome,group=g,n=sum(!is.na(z[[outcome]])),estimate=mm$estimate,ci_low=max(0,mm$estimate-tc*se),ci_high=min(1,mm$estimate+tc*se))
  }
  gd<-hi$gradient-lo$gradient;sd<-sqrt(drop(t(gd)%*%vc%*%gd))
  cap_rows[[length(cap_rows)+1]]<-data.frame(period=p,outcome=outcome,group="Difference",n=NA,estimate=hi$estimate-lo$estimate,ci_low=hi$estimate-lo$estimate-tc*sd,ci_high=hi$estimate-lo$estimate+tc*sd)
 }
}
write.csv(do.call(rbind,cap_rows),file.path(outdir,"main_comparator_table.csv"),row.names=FALSE)
print(with(cc$variables,table(race,birthplace)))
print(kappa(cov2cor(vcov(fit))))
