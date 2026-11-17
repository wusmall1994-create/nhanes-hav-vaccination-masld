# Final publication outputs and deterministic numerical checks.
read_result <- function(file) read.csv(file.path("outputs", file), check.names = FALSE)
destination <- file.path("outputs", "final_tables")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
joint <- read_result("joint/joint_distribution.csv")
conditional <- read_result("joint/conditional_summary.csv")
cap <- read_result("joint/cap_joint_comparison.csv")
external <- read_result("joint/us2000_joint_distribution.csv")
age_cells <- read_result("joint/us2000_age_cells.csv")
age_indicator <- read_result("joint/known_history_age_groups.csv")
periods <- unique(joint$period)
for (period in periods) {
  z <- joint[joint$period == period, ]
  stopifnot(nrow(z) == 8L, abs(sum(z$estimate) - 1) < 1e-10,
            sum(z$events) == unique(z$n),
            all(z$ci_low <= z$estimate), all(z$estimate <= z$ci_high))
  z <- external[external$period == period &
                  external$category != "Unknown history (either result)", ]
  stopifnot(nrow(z) == 8L, abs(sum(z$estimate) - 1) < 1e-10)
  for (dose in c("No doses", "Fewer than 2 excluding zero", "At least 2", "Unknown")) {
    negative <- joint[joint$period == period &
                        joint$category == paste("Negative", dose, sep = " / "), ]
    positive <- joint[joint$period == period &
                        joint$category == paste("Positive", dose, sep = " / "), ]
    dose_estimate <- conditional[conditional$period == period &
        conditional$category == paste("Negative conditional on", dose), ]
    stopifnot(nrow(dose_estimate) == 1L,
              dose_estimate$n == negative$events + positive$events,
              dose_estimate$events == negative$events,
              abs(dose_estimate$estimate -
                 negative$estimate / (negative$estimate + positive$estimate)) < 1e-10)
  }
  for (outcome in unique(cap$outcome)) {
    z <- cap[cap$period == period & cap$outcome == outcome, ]
    estimates <- setNames(z$estimate, z$group)
    stopifnot(abs(estimates["Difference"] -
       (estimates["CAP >=285"] - estimates["CAP <285"])) < 1e-12)
  }
}
stopifnot(nrow(age_indicator) == 6L,
          identical(as.integer(age_indicator$n), c(30L, 217L, 369L, 25L, 90L, 209L)))
empirical <- read_result("joint/joint_age_standard.csv")
standard2000 <- read_result("joint/us2000_age_standard.csv")
stopifnot(abs(sum(empirical$weight) - 1) < 1e-10,
          abs(sum(standard2000$weight) - 1) < 1e-10,
          all(age_cells$sparse ==
              (age_cells$n < 30 | age_cells$events < 10 |
                 age_cells$n - age_cells$events < 10)))
publication_proportions <- function(data) {
  data$estimate_percent <- 100 * data$estimate
  data$ci_low_percent <- 100 * data$ci_low
  data$ci_high_percent <- 100 * data$ci_high
  data
}
write.csv(read_result("secondary/target_characteristics_by_period.csv"),
          file.path(destination, "Table_1.csv"), row.names = FALSE)
write.csv(publication_proportions(joint),
          file.path(destination, "Table_2.csv"), row.names = FALSE)
dose_conditional <- conditional[
  grepl("^Negative conditional on ", conditional$category), ]
stopifnot(nrow(dose_conditional) == 8L)
write.csv(publication_proportions(dose_conditional),
          file.path(destination, "Table_3.csv"), row.names = FALSE)
write.csv(publication_proportions(cap),
          file.path(destination, "Table_4.csv"), row.names = FALSE)
index <- data.frame(
  table = paste0("S", 1:20),
  source = c(
    "masld_hav_hbv/variable_dictionary.csv; harmonized/participant_flow.csv",
    "secondary/age_standardized_prevalence.csv; secondary/age_adjusted_period_differences.csv",
    "secondary/age_standardized_prevalence.csv; secondary/age_adjusted_period_differences.csv",
    "secondary/phenotype_sensitivity_temporal.csv",
    "secondary/model_diagnostics.csv",
    "secondary/age_cell_counts.csv",
    "secondary/selection_profiles.csv",
    "secondary/age_standardized_prevalence.csv; secondary/age_adjusted_group_differences.csv; secondary/additional_checks.csv; masld_hav_hbv_exploratory/main_comparator_table.csv",
    "secondary/target_alcohol_selection.csv; secondary/target_alcohol_skip.csv; secondary/alcohol_differences.csv",
    "masld_hav_hbv_exploratory/adjusted_odds_ratios.csv; masld_hav_hbv_exploratory/age_spline_sensitivity.csv; masld_hav_hbv_exploratory/complete_case_flow.csv",
    "masld_hav_hbv_exploratory/adjusted_prevalence_margins.csv",
    "masld_hav_hbv_exploratory/subgroup_prevalence.csv",
    "secondary/age_standardized_prevalence.csv; secondary/age_adjusted_period_differences.csv",
    "joint/known_history_age_groups.csv",
    "secondary/age_standardized_prevalence.csv; secondary/age_adjusted_group_differences.csv",
    "joint/birthplace_components.csv",
    "joint/simplified_models.csv; joint/sparse_birth_race.csv",
    "joint/service_descriptions.csv",
    "joint/alcohol_joint_distribution.csv",
    "joint/joint_robustness.csv; joint/joint_age_standard.csv; joint/us2000_joint_distribution.csv; joint/us2000_age_cells.csv; joint/us2000_age_standard.csv"))
for (source in unlist(strsplit(index$source, "; ", fixed = TRUE))) {
  stopifnot(file.exists(file.path("outputs", source)))
}
write.csv(index, file.path(destination, "Supplementary_table_sources.csv"),
          row.names = FALSE)
writeLines("All numerical checks passed.", file.path(destination, "validation.txt"))
