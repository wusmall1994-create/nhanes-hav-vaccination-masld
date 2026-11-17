library(survey)
options(survey.lonely.psu = "adjust")
x <- readRDS("outputs/joint/analysis_designs.rds")
source("scripts/survey_proportions.R")
full <- x$full
target <- subset(x$target, hav_valid %in% 1)
target <- update(target, age_band = cut(RIDAGEYR, c(19, 39, 59, Inf),
                                         labels = c("20-39", "40-59", "60+")))
# Census P25-1130 counts from the NCI 2000 US standard population.
source_counts <- data.frame(
  age = c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
          "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+"),
  count = c(18257225, 17722067, 19511370, 22179956, 22479229, 19805793,
            17224359, 13307234, 10654272, 9409940, 8725574, 7414559,
            4900234, 4259173),
  band = c(rep("20-39", 4), rep("40-59", 4), rep("60+", 6)))
standard <- aggregate(count ~ band, source_counts, sum)
standard <- standard[match(levels(target$variables$age_band), standard$band), ]
standard$weight <- standard$count / sum(standard$count)
stopifnot(abs(sum(standard$weight) - 1) < 1e-12)
ref <- standard$weight
rows <- list()
cells <- list()
for (period_value in levels(x$data$period)) {
  z <- subset(target, period == period_value)
  for (s in c("Negative", "Positive", "Either")) {
    for (dose in c("No doses", "Fewer than 2 excluding zero", "At least 2", "Unknown")) {
      if (s == "Either" && dose != "Unknown") next
      hit <- as.numeric((s == "Either" | z$variables$serology == s) &
                          z$variables$dose == dose)
      label <- if (s == "Either") "Unknown history (either result)" else paste(s, dose, sep = " / ")
      zz <- update(z, hit = hit)
      for (j in seq_along(ref)) {
        in_band <- z$variables$age_band == standard$band[j]
        zz$variables[[paste0("num", j)]] <- hit * in_band
        zz$variables[[paste0("den", j)]] <- as.numeric(in_band)
        cell <- subset(zz, age_band == standard$band[j])
        cells[[length(cells) + 1]] <- prop(cell, cell$variables$hit, label,
                                           period_value, standard$band[j])
      }
      totals <- svytotal(~num1 + den1 + num2 + den2 + num3 + den3, zz)
      contrast <- svycontrast(totals, bquote(.(ref[1]) * num1/den1 +
                            .(ref[2]) * num2/den2 + .(ref[3]) * num3/den3))
      estimate <- as.numeric(coef(contrast))
      se <- as.numeric(SE(contrast))
      df <- degf(subset(full, period == period_value))
      ci <- if (estimate > 0 && estimate < 1) {
        plogis(qlogis(estimate) + c(-1, 1) * qt(.975, df) * se /
                 (estimate * (1 - estimate)))
      } else c(NA, NA)
      rows[[length(rows) + 1]] <- data.frame(period = period_value,
        category = label, n = length(hit), events = sum(hit), estimate = estimate,
        se = se, ci_low = ci[1], ci_high = ci[2], df = df,
        sparse = length(hit) < 30 | sum(hit) < 10 | sum(1 - hit) < 10)
    }
  }
}
result <- do.call(rbind, rows)
sums <- aggregate(estimate ~ period,
                  subset(result, category != "Unknown history (either result)"), sum)
stopifnot(all(abs(sums$estimate - 1) < 1e-10))
write.csv(result, "outputs/joint/us2000_joint_distribution.csv", row.names = FALSE)
write.csv(do.call(rbind, cells), "outputs/joint/us2000_age_cells.csv", row.names = FALSE)
write.csv(standard, "outputs/joint/us2000_age_standard.csv", row.names = FALSE)
print(standard, row.names = FALSE)
print(subset(result, category %in% c("Negative / No doses", "Negative / At least 2",
                                    "Unknown history (either result)")), row.names = FALSE)
