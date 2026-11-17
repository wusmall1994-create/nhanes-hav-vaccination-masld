library(survey)
library(splines)
options(survey.lonely.psu = "adjust")
x <- readRDS("outputs/joint/analysis_designs.rds")
design <- x$full
v <- design$variables
v$capgroup <- factor(ifelse(v$LUXCAPM >= 285, "CAP >=285", "CAP <285"),
                     levels = c("CAP <285", "CAP >=285"))
v$negative_zero <- as.numeric(v$serology == "Negative" & v$dose == "No doses")
v$negative_two <- as.numeric(v$serology == "Negative" & v$dose == "At least 2")
v$unknown_history <- as.numeric(v$dose == "Unknown")
design$variables <- v
domain <- subset(design, base_eligible %in% 1 & diabetes_broad %in% 1 &
                   alcohol_masld_range %in% 1 & hav_valid %in% 1 &
                   !is.na(birthplace))
standard <- domain$variables
standard$standard_weight <- ave(standard$analysis_weight, standard$period,
                                FUN = function(w) w / sum(w) / 2)
stopifnot(abs(sum(standard$standard_weight) - 1) < 1e-12)
results <- list()
diagnostics <- list()
for (outcome in c("negative_zero", "negative_two", "unknown_history")) {
  fit <- svyglm(reformulate("period * capgroup + ns(RIDAGEYR, df = 3) + birthplace",
                           outcome),
                design = domain, family = quasibinomial())
  coefficients <- coef(fit)
  covariance <- vcov(fit)
  critical <- qt(.975, fit$df.residual)
  margin <- function(period, group) {
    newdata <- standard
    newdata$period <- factor(period, levels = levels(v$period))
    newdata$capgroup <- factor(group, levels = levels(v$capgroup))
    matrix <- model.matrix(delete.response(terms(fit)), newdata,
                           contrasts.arg = fit$contrasts, xlev = fit$xlevels)
    matrix <- matrix[, names(coefficients), drop = FALSE]
    predictions <- plogis(drop(matrix %*% coefficients))
    w <- newdata$standard_weight
    list(estimate = sum(w * predictions),
         gradient = colSums(matrix * as.numeric(w * predictions * (1 - predictions))))
  }
  diagnostics[[outcome]] <- data.frame(outcome = outcome, n = nrow(domain$variables),
      parameters = length(coefficients), residual_df = fit$df.residual,
      converged = fit$converged, standard_n = nrow(standard))
  for (period in levels(v$period)) {
    low <- margin(period, "CAP <285")
    high <- margin(period, "CAP >=285")
    for (group in levels(v$capgroup)) {
      estimate <- if (group == "CAP <285") low else high
      se <- sqrt(drop(t(estimate$gradient) %*% covariance %*% estimate$gradient))
      z <- domain$variables[domain$variables$period == period &
                              domain$variables$capgroup == group, ]
      events <- sum(z[[outcome]])
      n <- nrow(z)
      results[[length(results) + 1]] <- data.frame(period = period, outcome = outcome,
          group = group, n = n, events = events, estimate = estimate$estimate,
          ci_low = max(0, estimate$estimate - critical * se),
          ci_high = min(1, estimate$estimate + critical * se),
          sparse = n < 30 | events < 10 | n - events < 10)
    }
    gradient <- high$gradient - low$gradient
    se <- sqrt(drop(t(gradient) %*% covariance %*% gradient))
    difference <- high$estimate - low$estimate
    results[[length(results) + 1]] <- data.frame(period = period, outcome = outcome,
        group = "Difference", n = NA, events = NA, estimate = difference,
        ci_low = difference - critical * se, ci_high = difference + critical * se,
        sparse = NA)
  }
}
write.csv(do.call(rbind, results), "outputs/joint/cap_joint_comparison.csv", row.names = FALSE)
write.csv(do.call(rbind, diagnostics), "outputs/joint/cap_joint_diagnostics.csv", row.names = FALSE)
write.csv(aggregate(standard_weight ~ period + birthplace, standard, sum),
          "outputs/joint/cap_joint_standard.csv", row.names = FALSE)
print(do.call(rbind, results), row.names = FALSE)
