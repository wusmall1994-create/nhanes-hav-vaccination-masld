library(survey)
options(survey.lonely.psu="adjust")
x <- readRDS("outputs/joint/analysis_designs.rds")
full <- x$full
source("scripts/survey_proportions.R")
out <- lapply(levels(x$data$period),function(p){
 period_value <<- p
 a <- subset(x$target,period==p & hav_valid %in% 1 & dose=="No doses")
 prop(a,a$variables$serology=="Positive","Positive conditional on no doses",p)
})
z<-do.call(rbind,out)
stopifnot(all(z$n==c(458,252)),all(z$events==c(225,94)))
write.csv(z,"outputs/joint/zero_dose_positive.csv",row.names=FALSE)
print(z)
