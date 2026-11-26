# Estimate classification measures against observed total anti-HAV.
library(survey)
options(survey.lonely.psu="adjust")
x <- readRDS("outputs/joint/analysis_designs.rds")
full <- x$full
source("scripts/survey_proportions.R")
rows <- list()
for (period_value in levels(full$variables$period)) {
  a <- subset(x$target, period == period_value & hav_valid %in% 1 & dose != "Unknown")
  stopifnot(nrow(a$variables) == c(616,324)[match(period_value,levels(full$variables$period))])
  for (threshold in c("At least two doses","Any dose")) {
    yes <- if (threshold == "Any dose") a$variables$dose != "No doses" else a$variables$dose == "At least 2"
    a <- update(a, reported=yes, positive=serology=="Positive")
    for (metric in c("Sensitivity","Specificity","PPV","NPV")) {
      z <- switch(metric,Sensitivity=subset(a,positive),Specificity=subset(a,!positive),PPV=subset(a,reported),NPV=subset(a,!reported))
      hit <- with(z$variables,switch(metric,Sensitivity=reported,Specificity=!reported,PPV=positive,NPV=!positive))
      rows[[length(rows)+1]] <- prop(z,hit,metric,period_value,threshold)
    }
  }
}
results <- do.call(rbind,rows)
write.csv(results,"outputs/joint/dose_history_performance.csv",row.names=FALSE)
print(results[,c("period","domain","category","events","n","estimate","ci_low","ci_high")],row.names=FALSE)

