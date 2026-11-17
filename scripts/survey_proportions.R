# Pointwise survey-logit CIs, period-specific adult-universe df; no CI for boundary cells.
prop <- function(des,hit,label,period,domain='All valid HAV') {
 des <- update(des,hit=as.numeric(hit)); y<-des$variables$hit;n<-sum(!is.na(y));events<-sum(y==1,na.rm=TRUE)
 df<-degf(subset(full,period==period_value))
 a<-svymean(~hit,des,na.rm=TRUE); est<-as.numeric(coef(a));se<-as.numeric(SE(a));ci<-c(NA_real_,NA_real_)
 if(events>0 & events<n) ci<-as.numeric(confint(svyciprop(~hit,des,method='logit',df=df,na.rm=TRUE)))
 data.frame(period=period,domain=domain,category=label,n=n,events=events,estimate=est,se=se,ci_low=ci[1],ci_high=ci[2],df=df,sparse=n<30|events<10|(n-events)<10)
}
