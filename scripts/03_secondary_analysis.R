#!/usr/bin/env Rscript

# Age-standardized secondary analysis of HAV/HBV serologic negativity among
# US adults with diabetes, with and without probable MASLD.
#
# Primary standardization model:
#   outcome ~ survey period * MASLD status + natural spline(age, df = 3)
#
# Predictive margins are standardized to a common empirical age distribution:
# an equal mixture of the two survey-period age distributions among all
# base-eligible adults with diabetes. Linear-age and categorical-age models are
# retained as sensitivity analyses. The categorical model is not primary
# because the youngest current-period cells are sparse.

suppressPackageStartupMessages({
  library(survey)
  library(splines)
})

options(survey.lonely.psu = "adjust")

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else "."
input_file <- file.path(root, "outputs", "harmonized", "harmonized_minimal_analysis.rds")
outdir <- file.path(root, "outputs", "secondary")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Missing harmonized dataset: ", input_file)

d <- readRDS(input_file)
period_levels <- c("2017-March 2020", "August 2021-August 2023")
group_levels <- c("Diabetes only", "MASLD + diabetes")
d$period <- factor(as.character(d$period), levels = period_levels)
d$group <- ifelse(
  d$target_masld_diabetes == 1, "MASLD + diabetes",
  ifelse(d$diabetes_only == 1, "Diabetes only", NA_character_)
)
d$group <- factor(d$group, levels = group_levels)
d$age_group <- cut(
  d$RIDAGEYR, breaks = c(19, 39, 59, Inf),
  labels = c("20-39", "40-59", "60+"), right = TRUE
)
d$female <- ifelse(d$RIAGENDR %in% c(1, 2), as.integer(d$RIAGENDR == 2), NA_real_)
d$race_hispanic <- ifelse(d$RIDRETH3 %in% c(1, 2, 3, 4, 6, 7), as.integer(d$RIDRETH3 %in% c(1, 2)), NA_real_)
d$race_white <- ifelse(d$RIDRETH3 %in% c(1, 2, 3, 4, 6, 7), as.integer(d$RIDRETH3 == 3), NA_real_)
d$race_black <- ifelse(d$RIDRETH3 %in% c(1, 2, 3, 4, 6, 7), as.integer(d$RIDRETH3 == 4), NA_real_)
d$race_asian <- ifelse(d$RIDRETH3 %in% c(1, 2, 3, 4, 6, 7), as.integer(d$RIDRETH3 == 6), NA_real_)
d$race_other <- ifelse(d$RIDRETH3 %in% c(1, 2, 3, 4, 6, 7), as.integer(d$RIDRETH3 == 7), NA_real_)
d$foreign_born <- ifelse(d$DMDBORN4 %in% c(1, 2), as.integer(d$DMDBORN4 == 2), NA_real_)
d$obesity <- ifelse(!is.na(d$BMXBMI), as.integer(d$BMXBMI >= 30), NA_real_)
d$physician_diabetes <- ifelse(d$DIQ010 %in% c(1, 2, 3), as.integer(d$DIQ010 == 1), NA_real_)
d$fibrosis_8 <- ifelse(!is.na(d$LUXSMED), as.integer(d$LUXSMED >= 8), NA_real_)
d$age_20_39 <- ifelse(!is.na(d$age_group), as.integer(d$age_group == "20-39"), NA_real_)
d$age_40_59 <- ifelse(!is.na(d$age_group), as.integer(d$age_group == "40-59"), NA_real_)
d$age_60_plus <- ifelse(!is.na(d$age_group), as.integer(d$age_group == "60+"), NA_real_)

# Nest masked design identifiers within survey period. Cycle-specific weights
# are retained: WTMECPRP-derived weights in the early period and WTPH2YR in the
# current period.
d$stack_stratum <- interaction(d$period, d$SDMVSTRA, drop = TRUE)
d$stack_psu <- interaction(d$period, d$SDMVSTRA, d$SDMVPSU, drop = TRUE)
universe <- d$adult20 == 1 & !is.na(d$analysis_weight) & d$analysis_weight > 0 &
  !is.na(d$SDMVSTRA) & !is.na(d$SDMVPSU)
du <- d[universe, ]
stacked_design <- svydesign(
  ids = ~stack_psu, strata = ~stack_stratum, weights = ~analysis_weight,
  nest = TRUE, data = du
)
analysis_design <- subset(stacked_design, base_eligible == 1 & !is.na(group))

# Descriptive characteristics of the combined MASLD-diabetes phenotype in each
# survey period. These are period-specific weighted descriptions, not tests of
# temporal change.
characteristic_specs <- list(
  list(characteristic = "Age", level = "Mean", variable = "RIDAGEYR", type = "continuous", unit = "years"),
  list(characteristic = "Age group", level = "20-39", variable = "age_20_39", type = "binary", unit = ""),
  list(characteristic = "Age group", level = "40-59", variable = "age_40_59", type = "binary", unit = ""),
  list(characteristic = "Age group", level = "60+", variable = "age_60_plus", type = "binary", unit = ""),
  list(characteristic = "Sex", level = "Female", variable = "female", type = "binary", unit = ""),
  list(characteristic = "Race and ethnicity", level = "Non-Hispanic White", variable = "race_white", type = "binary", unit = ""),
  list(characteristic = "Race and ethnicity", level = "Hispanic", variable = "race_hispanic", type = "binary", unit = ""),
  list(characteristic = "Race and ethnicity", level = "Non-Hispanic Black", variable = "race_black", type = "binary", unit = ""),
  list(characteristic = "Race and ethnicity", level = "Non-Hispanic Asian", variable = "race_asian", type = "binary", unit = ""),
  list(characteristic = "Race and ethnicity", level = "Other/multiracial", variable = "race_other", type = "binary", unit = ""),
  list(characteristic = "Nativity", level = "Foreign-born", variable = "foreign_born", type = "binary", unit = ""),
  list(characteristic = "Body mass index", level = "Mean", variable = "BMXBMI", type = "continuous", unit = "kg/m2"),
  list(characteristic = "Body mass index", level = ">=30 kg/m2", variable = "obesity", type = "binary", unit = ""),
  list(characteristic = "Hemoglobin A1c", level = "Mean", variable = "LBXGH", type = "continuous", unit = "%"),
  list(characteristic = "Physician-diagnosed diabetes", level = "Yes", variable = "physician_diabetes", type = "binary", unit = ""),
  list(characteristic = "Controlled attenuation parameter", level = "Mean", variable = "LUXCAPM", type = "continuous", unit = "dB/m"),
  list(characteristic = "Liver stiffness measurement", level = "Mean", variable = "LUXSMED", type = "continuous", unit = "kPa"),
  list(characteristic = "Liver stiffness measurement", level = ">=8 kPa", variable = "fibrosis_8", type = "binary", unit = "")
)

characteristic_rows <- list()
kc0 <- 1L
for (p in period_levels) {
  td <- subset(analysis_design, period == p & group == "MASLD + diabetes")
  for (spec in characteristic_specs) {
    y <- td$variables[[spec$variable]]
    n <- sum(!is.na(y))
    if (spec$type == "continuous") {
      est <- svymean(as.formula(paste0("~", spec$variable)), td, na.rm = TRUE)
      value <- unname(coef(est)[1])
      se <- unname(SE(est)[1])
      ci <- as.numeric(confint(est)[1, ])
      events <- NA_integer_
    } else {
      est <- svyciprop(as.formula(paste0("~", spec$variable)), td, method = "logit", na.rm = TRUE)
      value <- unname(coef(est)[1])
      se <- unname(SE(est)[1])
      ci <- as.numeric(confint(est)[1, ])
      events <- sum(y == 1, na.rm = TRUE)
    }
    characteristic_rows[[kc0]] <- data.frame(
      period = p, characteristic = spec$characteristic, level = spec$level,
      variable = spec$variable, type = spec$type, unit = spec$unit,
      n = n, events = events, estimate = value, se = se,
      ci_low = ci[1], ci_high = ci[2], stringsAsFactors = FALSE
    )
    kc0 <- kc0 + 1L
  }
}
target_characteristics <- do.call(rbind, characteristic_rows)
write.csv(target_characteristics, file.path(outdir, "target_characteristics_by_period.csv"), row.names = FALSE, na = "")

design_info <- data.frame(
  metric = c(
    "Stacked adult universe n", "Base-eligible diabetes n", "Stacked strata",
    "Stacked PSU pairs", "Stacked design df"
  ),
  value = c(
    nrow(stacked_design$variables), nrow(analysis_design$variables),
    length(unique(stacked_design$variables$stack_stratum)),
    length(unique(stacked_design$variables$stack_psu)), degf(stacked_design)
  ),
  stringsAsFactors = FALSE
)
write.csv(design_info, file.path(outdir, "stacked_design_summary.csv"), row.names = FALSE)

endpoint_labels <- c(
  hav_susceptible = "Total anti-HAV negativity",
  hbv_susceptible = "HBV triple-panel negativity",
  hav_joint_indicator = "HAV joint indicator",
  dual_susceptible = "Both serologic negativity definitions",
  either_susceptible = "Either serologic negativity definition"
)

# Cell-size audit before choosing a standardization estimator.
cell_rows <- list()
k <- 1L
for (ep in names(endpoint_labels)) {
  for (p in period_levels) {
    for (g in group_levels) {
      for (a in levels(d$age_group)) {
        keep <- d$base_eligible == 1 & d$period == p & d$group == g & d$age_group == a
        y <- d[[ep]][keep]
        n <- sum(!is.na(y))
        events <- sum(y == 1, na.rm = TRUE)
        nonevents <- sum(y == 0, na.rm = TRUE)
        flags <- character(0)
        if (n < 30) flags <- c(flags, "n<30")
        if (events < 10) flags <- c(flags, "events<10")
        if (nonevents < 10) flags <- c(flags, "nonevents<10")
        cell_rows[[k]] <- data.frame(
          period = p, group = g, age_group = a, endpoint = endpoint_labels[[ep]],
          endpoint_variable = ep, n = n, events = events, nonevents = nonevents,
          reliability = if (length(flags) == 0L) "Meets count screen" else paste("Use caution:", paste(flags, collapse = "; ")),
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
}
age_cells <- do.call(rbind, cell_rows)
write.csv(age_cells, file.path(outdir, "age_cell_counts.csv"), row.names = FALSE)

# Survey-weighted age composition by period and phenotype.
age_composition_rows <- list()
k <- 1L
for (p in period_levels) {
  for (g in group_levels) {
    sd <- subset(analysis_design, period == p & group == g)
    mean_age <- svymean(~RIDAGEYR, sd, na.rm = TRUE)
    for (a in levels(d$age_group)) {
      indicator_name <- paste0("age_ind_", gsub("[^0-9A-Za-z]", "_", a))
      sd$variables[[indicator_name]] <- as.integer(sd$variables$age_group == a)
      prop <- svymean(as.formula(paste0("~", indicator_name)), sd, na.rm = TRUE)
      age_composition_rows[[k]] <- data.frame(
        period = p, group = g, age_group = a,
        n = sum(sd$variables$age_group == a, na.rm = TRUE),
        weighted_proportion = unname(coef(prop)[1]),
        weighted_proportion_se = unname(SE(prop)[1]),
        weighted_mean_age = unname(coef(mean_age)[1]),
        weighted_mean_age_se = unname(SE(mean_age)[1]),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
}
age_composition <- do.call(rbind, age_composition_rows)
write.csv(age_composition, file.path(outdir, "age_composition.csv"), row.names = FALSE)

# Build a common empirical age standard. Each survey period contributes 50%
# after normalizing its survey weights among all base-eligible adults with
# diabetes, regardless of MASLD status.
standard <- analysis_design$variables[, c("RIDAGEYR", "period", "analysis_weight")]
standard$standard_weight <- NA_real_
for (p in period_levels) {
  idx <- standard$period == p
  standard$standard_weight[idx] <- 0.5 * standard$analysis_weight[idx] / sum(standard$analysis_weight[idx])
}
stopifnot(abs(sum(standard$standard_weight) - 1) < 1e-10)
standard$age_group <- cut(
  standard$RIDAGEYR, breaks = c(19, 39, 59, Inf),
  labels = c("20-39", "40-59", "60+"), right = TRUE
)
standard_summary <- data.frame(
  component = c("Weighted mean age", levels(standard$age_group)),
  value = c(
    weighted.mean(standard$RIDAGEYR, standard$standard_weight),
    vapply(levels(standard$age_group), function(a) sum(standard$standard_weight[standard$age_group == a]), numeric(1))
  ),
  stringsAsFactors = FALSE
)
write.csv(standard_summary, file.path(outdir, "common_age_standard.csv"), row.names = FALSE)

model_specs <- list(
  primary_spline = ~ period * group + ns(RIDAGEYR, df = 3),
  sensitivity_linear = ~ period * group + RIDAGEYR,
  sensitivity_categorical = ~ period * group + age_group
)

make_prediction_matrix <- function(fit, newdata) {
  tt <- delete.response(terms(fit))
  X <- model.matrix(tt, data = newdata, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  beta_names <- names(coef(fit))
  missing_cols <- setdiff(beta_names, colnames(X))
  if (length(missing_cols) > 0L) stop("Prediction matrix missing: ", paste(missing_cols, collapse = ", "))
  X[, beta_names, drop = FALSE]
}

standardized_margin <- function(fit, period_value, group_value) {
  nd <- standard
  nd$period <- factor(period_value, levels = period_levels)
  nd$group <- factor(group_value, levels = group_levels)
  nd$age_group <- factor(nd$age_group, levels = c("20-39", "40-59", "60+"))
  X <- make_prediction_matrix(fit, nd)
  beta <- coef(fit)
  eta <- drop(X %*% beta)
  prob <- plogis(eta)
  w <- nd$standard_weight
  estimate <- sum(w * prob)
  gradient <- colSums(X * as.numeric(w * prob * (1 - prob)))
  names(gradient) <- colnames(X)
  list(estimate = estimate, gradient = gradient)
}

contrast_from_gradients <- function(fit, estimate, gradient, label, contrast_type, endpoint, spec) {
  v <- vcov(fit)
  gradient <- gradient[names(coef(fit))]
  se <- sqrt(drop(t(gradient) %*% v %*% gradient))
  df <- fit$df.residual
  crit <- qt(0.975, df = df)
  data.frame(
    model_spec = spec, endpoint = endpoint_labels[[endpoint]], endpoint_variable = endpoint,
    contrast_type = contrast_type, contrast = label, estimate = estimate, se = se,
    ci_low = estimate - crit * se, ci_high = estimate + crit * se,
    p_value = if (contrast_type == "standardized prevalence") NA_real_ else 2 * pt(-abs(estimate / se), df = df),
    model_residual_df = df, stringsAsFactors = FALSE
  )
}

margin_rows <- list()
period_contrast_rows <- list()
group_contrast_rows <- list()
interaction_rows <- list()
model_rows <- list()
km <- kp <- kg <- ki <- kc <- 1L

for (spec_name in names(model_specs)) {
  for (ep in names(endpoint_labels)) {
    formula <- as.formula(paste(ep, paste(deparse(model_specs[[spec_name]]), collapse = " ")))
    fit <- svyglm(formula, design = analysis_design, family = quasibinomial())
    if (!fit$converged || any(!is.finite(coef(fit)))) stop("Model failed: ", spec_name, " / ", ep)

    model_rows[[kc]] <- data.frame(
      model_spec = spec_name, endpoint = endpoint_labels[[ep]], endpoint_variable = ep,
      observations = nrow(model.frame(fit)), coefficients = length(coef(fit)),
      residual_df = fit$df.residual, converged = fit$converged,
      stringsAsFactors = FALSE
    )
    kc <- kc + 1L

    margins <- list()
    for (p in period_levels) {
      for (g in group_levels) {
        key <- paste(p, g, sep = "||")
        margins[[key]] <- standardized_margin(fit, p, g)
        rr <- contrast_from_gradients(
          fit, margins[[key]]$estimate, margins[[key]]$gradient,
          paste(p, g, sep = " | "), "standardized prevalence", ep, spec_name
        )
        rr$period <- p
        rr$group <- g
        rr$ci_low <- pmax(0, rr$ci_low)
        rr$ci_high <- pmin(1, rr$ci_high)
        margin_rows[[km]] <- rr
        km <- km + 1L
      }
    }

    for (g in group_levels) {
      early <- margins[[paste(period_levels[1], g, sep = "||")]]
      current <- margins[[paste(period_levels[2], g, sep = "||")]]
      period_contrast_rows[[kp]] <- contrast_from_gradients(
        fit, current$estimate - early$estimate, current$gradient - early$gradient,
        paste("Current minus early period in", g), "period difference", ep, spec_name
      )
      period_contrast_rows[[kp]]$group <- g
      kp <- kp + 1L
    }

    for (p in period_levels) {
      diab <- margins[[paste(p, "Diabetes only", sep = "||")]]
      both <- margins[[paste(p, "MASLD + diabetes", sep = "||")]]
      group_contrast_rows[[kg]] <- contrast_from_gradients(
        fit, both$estimate - diab$estimate, both$gradient - diab$gradient,
        paste("MASLD + diabetes minus diabetes only in", p), "group difference", ep, spec_name
      )
      group_contrast_rows[[kg]]$period <- p
      kg <- kg + 1L
    }

    e_d <- margins[[paste(period_levels[1], "Diabetes only", sep = "||")]]
    c_d <- margins[[paste(period_levels[2], "Diabetes only", sep = "||")]]
    e_m <- margins[[paste(period_levels[1], "MASLD + diabetes", sep = "||")]]
    c_m <- margins[[paste(period_levels[2], "MASLD + diabetes", sep = "||")]]
    did_est <- (c_m$estimate - e_m$estimate) - (c_d$estimate - e_d$estimate)
    did_grad <- (c_m$gradient - e_m$gradient) - (c_d$gradient - e_d$gradient)
    interaction_rows[[ki]] <- contrast_from_gradients(
      fit, did_est, did_grad,
      "Difference in period changes: MASLD + diabetes minus diabetes only",
      "period-by-MASLD interaction on absolute prevalence scale", ep, spec_name
    )
    ki <- ki + 1L
  }
}

margins_out <- do.call(rbind, margin_rows)
period_out <- do.call(rbind, period_contrast_rows)
group_out <- do.call(rbind, group_contrast_rows)
interaction_out <- do.call(rbind, interaction_rows)
models_out <- do.call(rbind, model_rows)

write.csv(margins_out, file.path(outdir, "age_standardized_prevalence.csv"), row.names = FALSE, na = "")
write.csv(period_out, file.path(outdir, "age_adjusted_period_differences.csv"), row.names = FALSE, na = "")
write.csv(group_out, file.path(outdir, "age_adjusted_group_differences.csv"), row.names = FALSE, na = "")
write.csv(interaction_out, file.path(outdir, "period_by_masld_interactions.csv"), row.names = FALSE, na = "")
write.csv(models_out, file.path(outdir, "model_diagnostics.csv"), row.names = FALSE, na = "")

# Phenotype-definition sensitivity analyses. Each model is limited to adults
# meeting the variant diabetes definition and the variant probable-MASLD CAP
# threshold. All predictive margins retain the same primary common age
# standard so that changes are comparable across definitions.
variant_specs <- list(
  "CAP >=285; broad diabetes" = list(cutoff = 285, diabetes = "diabetes_broad"),
  "CAP >=274; broad diabetes" = list(cutoff = 274, diabetes = "diabetes_broad"),
  "CAP >=248; broad diabetes" = list(cutoff = 248, diabetes = "diabetes_broad"),
  "CAP >=285; physician-diagnosed diabetes" = list(cutoff = 285, diabetes = "physician_diabetes")
)

standardized_period_margin <- function(fit, period_value) {
  nd <- standard
  nd$period <- factor(period_value, levels = period_levels)
  X <- make_prediction_matrix(fit, nd)
  beta <- coef(fit)
  prob <- plogis(drop(X %*% beta))
  w <- nd$standard_weight
  gradient <- colSums(X * as.numeric(w * prob * (1 - prob)))
  names(gradient) <- colnames(X)
  list(estimate = sum(w * prob), gradient = gradient)
}

variant_rows <- list()
kv <- 1L
for (variant_name in names(variant_specs)) {
  vs <- variant_specs[[variant_name]]
  diabetes_flag <- stacked_design$variables[[vs$diabetes]] == 1
  target_flag <- diabetes_flag & stacked_design$variables$LUXCAPM >= vs$cutoff &
    stacked_design$variables$alcohol_masld_range == 1 & stacked_design$variables$base_eligible == 1
  vd <- stacked_design[target_flag, ]
  for (ep in names(endpoint_labels)) {
    fit <- svyglm(
      as.formula(paste0(ep, " ~ period + ns(RIDAGEYR, df = 3)")),
      design = vd, family = quasibinomial()
    )
    early <- standardized_period_margin(fit, period_levels[1])
    current <- standardized_period_margin(fit, period_levels[2])
    v <- vcov(fit)
    df <- fit$df.residual
    crit <- qt(0.975, df = df)
    for (p in period_levels) {
      mm <- if (p == period_levels[1]) early else current
      se <- sqrt(drop(t(mm$gradient) %*% v %*% mm$gradient))
      variant_rows[[kv]] <- data.frame(
        variant = variant_name, endpoint = endpoint_labels[[ep]], endpoint_variable = ep,
        result_type = "standardized prevalence", period = p,
        n = sum(vd$variables$period == p & !is.na(vd$variables[[ep]])),
        estimate = mm$estimate, se = se,
        ci_low = max(0, mm$estimate - crit * se),
        ci_high = min(1, mm$estimate + crit * se),
        p_value = NA_real_, residual_df = df, stringsAsFactors = FALSE
      )
      kv <- kv + 1L
    }
    diff <- current$estimate - early$estimate
    grad <- current$gradient - early$gradient
    se <- sqrt(drop(t(grad) %*% v %*% grad))
    variant_rows[[kv]] <- data.frame(
      variant = variant_name, endpoint = endpoint_labels[[ep]], endpoint_variable = ep,
      result_type = "period difference", period = "Current minus early",
      n = NA_integer_, estimate = diff, se = se,
      ci_low = diff - crit * se, ci_high = diff + crit * se,
      p_value = 2 * pt(-abs(diff / se), df = df), residual_df = df,
      stringsAsFactors = FALSE
    )
    kv <- kv + 1L
  }
}
variant_out <- do.call(rbind, variant_rows)
write.csv(variant_out, file.path(outdir, "phenotype_sensitivity_temporal.csv"), row.names = FALSE, na = "")

fmt_pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f%%", 100 * x))
fmt_ci <- function(low, high) ifelse(is.na(low) | is.na(high), "NA", sprintf("%.1f%% to %.1f%%", 100 * low, 100 * high))
primary_margins <- margins_out[margins_out$model_spec == "primary_spline", ]
primary_period <- period_out[period_out$model_spec == "primary_spline", ]
primary_interaction <- interaction_out[interaction_out$model_spec == "primary_spline", ]
primary_vars <- c("hav_susceptible", "hbv_susceptible", "hav_joint_indicator")

report <- c(
  "# Formal temporal analysis: age-standardized HAV/HBV serologic negativity",
  "",
  "## Estimand and model",
  "",
  "The primary estimates are survey-weighted predictive margins from a logistic model containing survey period, MASLD status among adults with diabetes, their interaction, and a natural spline for age with 3 degrees of freedom. All margins use the same empirical age standard, formed as an equal mixture of the two survey-period age distributions among base-eligible adults with diabetes.",
  "",
  sprintf("The stacked design contains %d design degrees of freedom. The common standard has a mean age of %.1f years.", degf(stacked_design), standard_summary$value[1]),
  "",
  "## Primary age-standardized estimates",
  ""
)

for (ep in primary_vars) {
  for (g in group_levels) {
    z <- primary_margins[primary_margins$endpoint_variable == ep & primary_margins$group == g, ]
    z <- z[match(period_levels, z$period), ]
    dd <- primary_period[primary_period$endpoint_variable == ep & primary_period$group == g, ]
    report <- c(report, sprintf(
      "- %s, %s: %s (%s) in 2017-March 2020 and %s (%s) in August 2021-August 2023; adjusted change %s (%s), P=%.3f.",
      endpoint_labels[[ep]], g,
      fmt_pct(z$estimate[1]), fmt_ci(z$ci_low[1], z$ci_high[1]),
      fmt_pct(z$estimate[2]), fmt_ci(z$ci_low[2], z$ci_high[2]),
      fmt_pct(dd$estimate), fmt_ci(dd$ci_low, dd$ci_high), dd$p_value
    ))
  }
}

report <- c(report, "", "## Period-by-MASLD interaction", "")
for (ep in primary_vars) {
  z <- primary_interaction[primary_interaction$endpoint_variable == ep, ]
  report <- c(report, sprintf(
    "- %s: difference in adjusted period changes %s (%s), P=%.3f.",
    endpoint_labels[[ep]], fmt_pct(z$estimate), fmt_ci(z$ci_low, z$ci_high), z$p_value
  ))
}

flagged_cells <- age_cells[age_cells$reliability != "Meets count screen", ]
report <- c(
  report,
  "",
  "## Robustness and reporting decision",
  "",
  sprintf("- %d of %d age-specific cells were flagged by the count screen; this is why direct age-group standardization is not the primary estimator.", nrow(flagged_cells), nrow(age_cells)),
  "- Linear-age and categorical-age specifications are provided as sensitivity analyses in the CSV outputs.",
  "- Period-by-MASLD interactions are exploratory and should not be interpreted as causal policy effects.",
  ""
)

writeLines(report, file.path(outdir, "formal_temporal_report.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(outdir, "session_info.txt"))
cat(paste(report, collapse = "\n"), "\n")

library(haven)
fd <- file.path(root,"outputs","harmonized")
flow <- read.csv(file.path(fd,"sample_flow.csv"))
flow$invalid_weight <- flow$alcohol_missing <- flow$module_missing <- flow$response_missing <- flow$above_threshold_comparator <- 0
profiles <- cross <- list()
for (i in seq_along(period_levels)) {
  p <- period_levels[i]
  raw <- read_xpt(file.path(root,c("data/raw_2017_2020/P_ALQ.xpt","data/raw/ALQ_L.xpt")[i]))
  z <- subset(d, period==p & adult20==1 & vcte_reliable==1)
  valid <- !is.na(z$analysis_weight) & z$analysis_weight>0
  missing <- z$alcohol_classifiable==0
  included <- valid & !missing
  flow$invalid_weight[i] <- sum(!valid)
  flow$alcohol_missing[i] <- sum(valid & missing)
  flow$module_missing[i] <- sum(valid & missing & !(z$SEQN %in% raw$SEQN))
  flow$response_missing[i] <- sum(valid & missing & z$SEQN %in% raw$SEQN)
  flow$above_threshold_comparator[i] <- sum(included & z$diabetes_only %in% 1 & z$alcohol_masld_range %in% 0)
  stopifnot(sum(included)==flow$base_eligible[i])
  for (kept in c(TRUE,FALSE)){
    x <- z[included==kept,]
    profiles[[length(profiles)+1]] <- data.frame(
      period=p,included=kept,n=nrow(x),age=mean(x$RIDAGEYR),
      women_n=sum(x$female==1,na.rm=TRUE),women=mean(x$female,na.rm=TRUE),
      foreign_n=sum(x$foreign_born==1,na.rm=TRUE),foreign=mean(x$foreign_born,na.rm=TRUE),
      diabetes_n=sum(x$diabetes_broad==1,na.rm=TRUE),diabetes=mean(x$diabetes_broad,na.rm=TRUE),
      cap_n=sum(x$LUXCAPM>=285),cap=mean(x$LUXCAPM>=285))
  }
  h <- subset(analysis_design, period==p & group=="MASLD + diabetes" & hav_valid==1)
  for (code in c(3,2,1,9,7,NA)) {
    v <- h$variables
    dose <- if (is.na(code)) is.na(v$IMQ011) else v$IMQ011 %in% code
    for (neg in 1:0) {
      hit <- as.numeric(dose & v$hav_susceptible==neg)
      est <- svymean(~hit, update(h,hit=hit))
      cross[[length(cross)+1]] <- data.frame(period=p,code=code,negative=neg,n=sum(hit),N=length(hit),estimate=coef(est)[1],se=SE(est)[1])
    }
  }
  h <- update(h, unknown=as.integer(!IMQ011 %in% c(1,2,3)))
  est <- svymean(~unknown,h)
  cat(p,"weighted unknown",coef(est)[1],"included weight*:",sum(z$analysis_weight[included]),"\n")
}
write.csv(flow,file.path(fd,"participant_flow.csv"),row.names=FALSE)
write.csv(do.call(rbind,profiles),file.path(outdir,"selection_profiles.csv"),row.names=FALSE)
write.csv(do.call(rbind,cross),file.path(outdir,"hav_dose_cross.csv"),row.names=FALSE)
model_add <- function(a, std, mode, outcome, form) {
  standard <<- std
  fit <- svyglm(form,design=a,family=quasibinomial())
  stopifnot(fit$converged, all(is.finite(coef(fit))))
  rows <- list()
  target_margins <- lapply(period_levels, function(p) standardized_margin(fit,p,group_levels[2]))
  names(target_margins) <- period_levels
  for (p in period_levels) {
    hi <- target_margins[[p]]
    lo <- standardized_margin(fit,p,group_levels[1])
    for (type in c("target","group difference")) {
      diff <- type=="group difference"
      row <- contrast_from_gradients(fit,hi$estimate-if(diff) lo$estimate else 0,
        hi$gradient-if(diff) lo$gradient else 0,p,type,outcome,mode)
      rows[[length(rows)+1]] <- row
    }
  }
  if (mode == "age and nativity") {
    early <- target_margins[[1]]
    late <- target_margins[[2]]
    rows[[length(rows)+1]] <- contrast_from_gradients(fit,late$estimate-early$estimate,
      late$gradient-early$gradient,"Later minus earlier","period difference",outcome,mode)
  }
  do.call(rbind,rows)
}
std_age <- standard
a <- update(analysis_design,
  lower=ifelse(hav_valid==1,as.integer(hav_joint_indicator %in% 1),NA_real_),
  upper=ifelse(hav_valid==1,as.integer(hav_joint_indicator %in% 1 | (hav_susceptible==1 & !IMQ011 %in% c(1,2,3))),NA_real_))
out <- list(
  model_add(a,std_age,"lower","hav_joint_indicator",lower ~ period*group + ns(RIDAGEYR,df=3)),
  model_add(a,std_age,"upper","hav_joint_indicator",upper ~ period*group + ns(RIDAGEYR,df=3)))
nat <- subset(analysis_design,!is.na(foreign_born))
std_nat <- nat$variables[,c("RIDAGEYR","period","analysis_weight","foreign_born")]
std_nat$age_group <- cut(std_nat$RIDAGEYR, breaks=c(19,39,59,Inf), labels=c("20-39","40-59","60+"))
std_nat$standard_weight <- ave(std_nat$analysis_weight,std_nat$period,FUN=function(w) w/sum(w)/2)
cap <- subset(nat,alcohol_masld_range==1)
for (ep in c("hav_joint_indicator","hav_susceptible","hbv_susceptible")){
  out[[length(out)+1]] <- model_add(nat,std_nat,"age and nativity",ep,reformulate("period*group + ns(RIDAGEYR,df=3) + foreign_born",ep))
  out[[length(out)+1]] <- model_add(cap,std_nat,"age and nativity; low alcohol",ep,reformulate("period*group + ns(RIDAGEYR,df=3) + foreign_born",ep))
  out[[length(out)+1]] <- model_add(subset(analysis_design,alcohol_masld_range==1),std_age,"age; low alcohol",ep,reformulate("period*group + ns(RIDAGEYR,df=3)",ep))
}
write.csv(do.call(rbind,out),file.path(outdir,"additional_checks.csv"),row.names=FALSE)


# Selection audit before applying the alcohol eligibility criterion.
d <- readRDS(file.path(root,"outputs/harmonized/harmonized_minimal_analysis.rds"))
d <- d[d$adult20 %in% 1 & !is.na(d$analysis_weight) & d$analysis_weight>0,]
d$status <- ifelse(d$alcohol_classifiable==0,"Unclassifiable",ifelse(d$alcohol_masld_range==1,"Below threshold","At or above threshold"))
d$subgroup <- as.integer(d$vcte_reliable==1 & d$diabetes_broad==1 & d$LUXCAPM>=285)
d$anydose <- ifelse(d$IMQ011 %in% 1:3,as.integer(d$IMQ011 %in% 1:2),NA)
d$twodose <- ifelse(d$IMQ011 %in% 1:3,as.integer(d$IMQ011==1),NA)
d$history_unknown <- as.integer(!d$IMQ011 %in% 1:3)
d$skip_reason <- ifelse(d$ALQ111 %in% 2,"Lifetime abstainer",
 ifelse(d$ALQ121 %in% 0,"No drinks in past year",
 ifelse(is.na(d$ALQ111),"ALQ111 missing",
 ifelse(d$ALQ111 %in% c(7,9),"ALQ111 refused/unknown",
 ifelse(is.na(d$ALQ121),"ALQ121 missing",
 ifelse(d$ALQ121 %in% c(77,99),"ALQ121 refused/unknown",ifelse(d$alcohol_classifiable==1,"Drinking quantity available","Quantity missing/refused/unknown")))))))
res <- skips <- differences <- list()
for(p in unique(d$period)){
 u <- d[d$period==p,]
 des <- svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~analysis_weight,nest=TRUE,data=u)
 sub <- subset(des,subgroup==1)
 for(g in c("Below threshold","At or above threshold","Unclassifiable","All classifiable")){
   hit <- as.integer(if (g=="All classifiable") sub$variables$status!="Unclassifiable" else sub$variables$status==g)
   pr <- svymean(~hit,update(sub,hit=hit))
   sd <- sub[hit==1,]
   for(v in c("hav_susceptible","anydose","twodose","history_unknown","hav_joint_indicator")){
     y <- sd$variables[[v]]
     est <- svymean(reformulate(v),sd,na.rm=TRUE)
     crit <- qt(.975,degf(des))
     res[[length(res)+1]] <- data.frame(period=p,status=g,n=nrow(sd$variables),subgroup_n=nrow(sub$variables),
       weighted_group=coef(pr)[1],variable=v,valid=sum(!is.na(y)),events=sum(y==1,na.rm=TRUE),
       estimate=coef(est)[1],ci_low=max(0,coef(est)[1]-crit*SE(est)[1]),ci_high=min(1,coef(est)[1]+crit*SE(est)[1]))
   }
   for(r in unique(sd$variables$skip_reason)){
     skips[[length(skips)+1]] <- data.frame(period=p,status=g,reason=r,n=sum(sd$variables$skip_reason==r))
   }
 }
 # Estimate within-design contrasts, retaining covariance between domains.
 for (comparator in c("Below threshold","All classifiable")) {
   comparison <- if(comparator=="Below threshold") subset(sub,status %in% c("Below threshold","Unclassifiable")) else sub
   comparison <- update(comparison,missing_alcohol=as.integer(status=="Unclassifiable"))
   for(v in c("hav_susceptible","hav_joint_indicator")){
     means <- svyby(reformulate(v),~missing_alcohol,comparison,svymean,na.rm=TRUE,covmat=TRUE)
     contrast <- svycontrast(means,c(-1,1))
     estimate <- as.numeric(coef(contrast))
     se <- as.numeric(SE(contrast))
     critical <- qt(.975,degf(des))
     differences[[length(differences)+1]] <- data.frame(
       period=p,comparator=comparator,variable=v,
       estimate=estimate,se=se,ci_low=estimate-critical*se,
       ci_high=estimate+critical*se,design_df=degf(des))
   }
 }

}
write.csv(do.call(rbind,res),file.path(outdir,"target_alcohol_selection.csv"),row.names=FALSE)
write.csv(do.call(rbind,skips),file.path(outdir,"target_alcohol_skip.csv"),row.names=FALSE)
write.csv(do.call(rbind,differences),file.path(outdir,"alcohol_differences.csv"),row.names=FALSE)
