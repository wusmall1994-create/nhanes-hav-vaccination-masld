library(survey)
options(survey.lonely.psu = "adjust")
x <- readRDS("outputs/joint/analysis_designs.rds")
full <- x$full
source("scripts/survey_proportions.R")
target <- subset(x$target, !is.na(hav_joint_indicator))
target <- update(target, age_band = cut(RIDAGEYR, c(19, 39, 59, Inf),
                     labels = c("20-39", "40-59", "60+")))
rows <- list()
for (period_value in levels(x$data$period)) {
  for (age in levels(target$variables$age_band)) {
    z <- subset(target, period == period_value & age_band == age)
    rows[[length(rows) + 1]] <- prop(z, z$variables$hav_joint_indicator,
       "Negative anti-HAV and fewer than two reported doses", period_value, age)
  }
}
write.csv(do.call(rbind, rows), "outputs/joint/known_history_age_groups.csv",
          row.names = FALSE)
