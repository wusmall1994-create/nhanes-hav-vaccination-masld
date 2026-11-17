library(survey)
options(survey.lonely.psu="adjust")
x <- readRDS("outputs/joint/analysis_designs.rds")
full <- x$full
target <- x$target
source("scripts/survey_proportions.R")
out <- list()
add <- function(des,hit,label,scenario) {
 out[[length(out)+1]] <<- prop(des,hit,label,period_value,scenario)
}
for (period_value in levels(x$data$period)) {
 p <- period_value
 a <- subset(target,period==p)
 add(a, !a$variables$hav_valid %in% 1,"Missing HAV serology","Primary target")
 base <- subset(full,period==p & vcte_reliable %in% 1 & alcohol_masld_range %in% 1 & diabetes_broad %in% 1 & hav_valid %in% 1)
 for (scenario in c("CAP 285","CAP 274","CAP 248","Diagnosed diabetes")) {
  threshold <- switch(scenario,"CAP 274"=274,"CAP 248"=248,285)
  z <- subset(base,LUXCAPM>=threshold)
  if(scenario=="Diagnosed diabetes") z<-subset(z,DIQ010 %in% 1)
  v<-z$variables
  for(s in c("Negative","Positive")) for(k in c("No doses","Fewer than 2 excluding zero","At least 2","Unknown")) add(z,v$serology==s & v$dose==k,paste(s,k,sep=" / "),scenario)
  add(z,v$dose=="Unknown","Unknown history (either result)",scenario)
 }
}
# Direct age standardization, fixed equal-period mixture of valid-serology target age distributions.
a <- subset(target,hav_valid %in% 1)
a <- update(a, age_band=cut(RIDAGEYR,c(19,39,59,Inf)))
ps <- levels(x$data$period)
ref <- Reduce("+",lapply(ps,function(p) coef(svymean(~age_band,subset(a,period==p)))))/2
for(period_value in ps){
 p<-period_value
 z<-subset(a,period==p)
 for(s in c("Negative","Positive","Either")) for(k in c("No doses","Fewer than 2 excluding zero","At least 2","Unknown")){
  if(s=="Either" && k!="Unknown") next
  hit<-as.numeric((s=="Either"|z$variables$serology==s)&z$variables$dose==k)
  zz<-update(z,hit=hit)
  lev<-levels(z$variables$age_band)
  for(j in seq_along(lev)){
   zz$variables[[paste0("num",j)]]<-hit*(z$variables$age_band==lev[j])
   zz$variables[[paste0("den",j)]]<-as.numeric(z$variables$age_band==lev[j])
  }
  totals<-svytotal(~num1+den1+num2+den2+num3+den3,zz)
  q<-svycontrast(totals,bquote(.(ref[1])*num1/den1 + .(ref[2])*num2/den2 + .(ref[3])*num3/den3))
  est<-as.numeric(coef(q));se<-as.numeric(SE(q));df<-degf(subset(full,period==p))
  ci<-if(est>0&&est<1) plogis(qlogis(est)+c(-1,1)*qt(.975,df)*se/(est*(1-est))) else c(NA,NA)
  label<-if(s=="Either") "Unknown history (either result)" else paste(s,k,sep=" / ")
  out[[length(out)+1]]<-data.frame(period=p,domain="Age standardized",category=label,n=length(hit),events=sum(hit),estimate=est,se=se,ci_low=ci[1],ci_high=ci[2],df=df,sparse=length(hit)<30|sum(hit)<10|sum(1-hit)<10)
 }
}
result<-do.call(rbind,out)
stopifnot(all(abs(aggregate(estimate~period+domain,subset(result,category %in% unique(read.csv("outputs/joint/joint_distribution.csv")$category)),sum)$estimate-1)<1e-9))
write.csv(result,"outputs/joint/joint_robustness.csv",row.names=FALSE)
write.csv(data.frame(age=c("20-39","40-59","60+"),weight=unname(ref)),"outputs/joint/joint_age_standard.csv",row.names=FALSE)
print(subset(result,category %in% c("Missing HAV serology","Negative / No doses","Negative / At least 2","Unknown history (either result)"))[,1:8],row.names=FALSE)
