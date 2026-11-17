#!/usr/bin/env Rscript

# Data harmonization for the two-period analysis.
#
# This script harmonizes the variables required from NHANES 2017-March 2020
# and August 2021-August 2023, applies common phenotype and outcome definitions,
# and writes the intermediate analysis object used by the formal models.

suppressPackageStartupMessages({
  library(haven)
  library(survey)
})

options(survey.lonely.psu = "adjust")

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else "."
old_raw <- file.path(root, "data", "raw_2017_2020")
current_file <- file.path(root, "data", "processed", "masld_hav_hbv_analysis.csv")
outdir <- file.path(root, "outputs", "harmonized")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required_old <- c(
  "P_DEMO.xpt", "P_LUX.xpt", "P_DIQ.xpt", "P_GHB.xpt", "P_ALQ.xpt",
  "P_BMX.xpt", "P_HEPA.xpt", "P_HEPB_S.xpt", "P_HEPBD.xpt", "P_IMQ.xpt"
)
missing_old <- required_old[!file.exists(file.path(old_raw, required_old))]
if (length(missing_old) > 0L) {
  stop("Missing early-period files: ", paste(missing_old, collapse = ", "))
}
if (!file.exists(current_file)) stop("Missing current processed file: ", current_file)

read_cols <- function(file, columns) {
  dat <- as.data.frame(read_xpt(file))
  missing <- setdiff(columns, names(dat))
  if (length(missing) > 0L) stop(basename(file), " missing columns: ", paste(missing, collapse = ", "))
  dat[, columns, drop = FALSE]
}

merge_module <- function(base, file, columns) {
  other <- read_cols(file, c("SEQN", columns))
  if (anyDuplicated(other$SEQN)) stop(basename(file), " is not participant-level")
  merge(base, other, by = "SEQN", all.x = TRUE, sort = FALSE)
}

derive_old_period <- function() {
  demo_columns <- c(
    "SEQN", "RIDSTATR", "RIAGENDR", "RIDAGEYR", "RIDRETH3", "DMDBORN4",
    "WTMECPRP", "SDMVSTRA", "SDMVPSU"
  )
  d <- read_cols(file.path(old_raw, "P_DEMO.xpt"), demo_columns)
  modules <- list(
    P_LUX = c("LUAXSTAT", "LUANMVGP", "LUXSMED", "LUXSIQRM", "LUXCAPM"),
    P_DIQ = c("DIQ010"),
    P_GHB = c("LBXGH"),
    P_ALQ = c("ALQ111", "ALQ121", "ALQ130"),
    P_BMX = c("BMXBMI"),
    P_HEPA = c("LBXHA"),
    P_HEPB_S = c("LBXHBS"),
    P_HEPBD = c("LBXHBC", "LBDHBG"),
    P_IMQ = c("IMQ011")
  )
  for (nm in names(modules)) {
    d <- merge_module(d, file.path(old_raw, paste0(nm, ".xpt")), modules[[nm]])
  }

  d$period <- "2017-March 2020"
  d$analysis_weight <- d$WTMECPRP
  d$adult20 <- as.integer(d$RIDAGEYR >= 20)
  d$vcte_complete <- as.integer(d$LUAXSTAT == 1)
  d$vcte_reliable <- as.integer(
    d$LUAXSTAT == 1 & d$LUANMVGP >= 10 & d$LUXSIQRM <= 30 &
      !is.na(d$LUXCAPM) & !is.na(d$LUXSMED)
  )
  d$diabetes_broad <- as.integer(d$DIQ010 == 1 | (!is.na(d$LBXGH) & d$LBXGH >= 6.5))

  days_per_year <- c(
    `0` = 0, `1` = 365, `2` = 300, `3` = 182, `4` = 104, `5` = 52,
    `6` = 30, `7` = 12, `8` = 9, `9` = 4.5, `10` = 1.5
  )
  d$ALQ130[d$ALQ130 %in% c(777, 999)] <- NA_real_
  d$drink_days_year <- unname(days_per_year[as.character(d$ALQ121)])
  d$alcohol_g_day <- d$drink_days_year * d$ALQ130 * 14 / 365
  nondrinker <- d$ALQ111 == 2 | d$ALQ121 == 0
  d$alcohol_g_day[nondrinker %in% TRUE] <- 0
  d$alcohol_classifiable <- as.integer(!is.na(d$alcohol_g_day))
  female <- d$RIAGENDR == 2
  d$alcohol_masld_range <- ifelse(
    is.na(d$alcohol_g_day), NA_integer_,
    as.integer(ifelse(female, d$alcohol_g_day < 20, d$alcohol_g_day < 30))
  )

  # In adults with diabetes, diabetes itself satisfies the MASLD
  # cardiometabolic criterion. Other CMRFs are not needed for this
  # analysis of the combined phenotype and diabetes-only comparator.
  d$masld_in_diabetes <- as.integer(
    d$diabetes_broad == 1 & d$LUXCAPM >= 285 & d$alcohol_masld_range == 1
  )
  d$target_masld_diabetes <- as.integer(d$masld_in_diabetes == 1)
  d$diabetes_only <- as.integer(d$diabetes_broad == 1 & d$masld_in_diabetes == 0)

  d$hav_valid <- as.integer(d$LBXHA %in% c(1, 2))
  d$hav_susceptible <- ifelse(d$LBXHA %in% c(1, 2), as.integer(d$LBXHA == 2), NA_real_)
  d$havax_incomplete_none <- ifelse(
    d$IMQ011 %in% c(1, 2, 3), as.integer(d$IMQ011 %in% c(2, 3)), NA_real_
  )
  d$hav_joint_indicator <- ifelse(
    !is.na(d$hav_susceptible) & !is.na(d$havax_incomplete_none),
    as.integer(d$hav_susceptible == 1 & d$havax_incomplete_none == 1), NA_real_
  )

  hbv_valid <- d$LBDHBG %in% c(1, 2) & d$LBXHBC %in% c(1, 2) & d$LBXHBS %in% c(1, 2)
  d$hbv_valid <- as.integer(hbv_valid)
  d$hbv_susceptible <- ifelse(
    hbv_valid,
    as.integer(d$LBDHBG == 2 & d$LBXHBC == 2 & d$LBXHBS == 2),
    NA_real_
  )
  d$dual_susceptible <- ifelse(
    !is.na(d$hav_susceptible) & !is.na(d$hbv_susceptible),
    as.integer(d$hav_susceptible == 1 & d$hbv_susceptible == 1), NA_real_
  )
  d$either_susceptible <- ifelse(
    !is.na(d$hav_susceptible) & !is.na(d$hbv_susceptible),
    as.integer(d$hav_susceptible == 1 | d$hbv_susceptible == 1), NA_real_
  )
  d$base_eligible <- as.integer(
    d$adult20 == 1 & d$vcte_reliable == 1 & d$alcohol_classifiable == 1 &
      !is.na(d$analysis_weight) & d$analysis_weight > 0
  )
  d
}

derive_current_period <- function() {
  d <- read.csv(current_file, check.names = FALSE)
  # Recover integer questionnaire codes directly from the original XPT.
  # Read alcohol values directly from the source transport file.
  alcohol <- as.data.frame(read_xpt(file.path(root, "data", "raw", "ALQ_L.xpt")))
  idx <- match(d$SEQN, alcohol$SEQN)
  for (variable in c("ALQ111", "ALQ121", "ALQ130")) d[[variable]] <- alcohol[[variable]][idx]
  d$period <- "August 2021-August 2023"
  # Read sampling weights from the original laboratory module to preserve zeros.
  weights <- as.data.frame(read_xpt(file.path(root, "data", "raw", "HEPA_L.xpt")))
  d$WTPH2YR <- weights$WTPH2YR[match(d$SEQN, weights$SEQN)]
  d$analysis_weight <- d$WTPH2YR
  # Apply the same alcohol cleaning to the late-period processed input.
  # The source CSV retains the NHANES questionnaire variables prior to recoding.
  days_per_year <- c(
    `0` = 0, `1` = 365, `2` = 300, `3` = 182, `4` = 104, `5` = 52,
    `6` = 30, `7` = 12, `8` = 9, `9` = 4.5, `10` = 1.5
  )
  d$ALQ130[d$ALQ130 %in% c(777, 999)] <- NA_real_
  d$alcohol_g_day <- unname(days_per_year[as.character(d$ALQ121)]) * d$ALQ130 * 14 / 365
  d$alcohol_g_day[d$ALQ111 %in% 2 | d$ALQ121 %in% 0] <- 0
  d$alcohol_classifiable <- as.integer(!is.na(d$alcohol_g_day))
  d$alcohol_masld_range <- as.integer(d$alcohol_g_day < ifelse(d$RIAGENDR == 2, 20, 30))
  d$base_eligible <- as.integer(
    d$adult20 == 1 & d$vcte_reliable == 1 & d$alcohol_classifiable == 1 &
      !is.na(d$analysis_weight) & d$analysis_weight > 0
  )
  d$masld <- as.integer(
    !is.na(d$LUXCAPM) & d$LUXCAPM >= 285 & d$cmrf_any == 1 & d$alcohol_masld_range %in% 1
  )
  d$target_masld_diabetes <- as.integer(d$masld == 1 & d$diabetes_broad == 1)
  d$diabetes_only <- as.integer(d$diabetes_broad == 1 & d$masld == 0)
  needed <- c(
    "SEQN", "period", "RIDAGEYR", "SDMVSTRA", "SDMVPSU", "analysis_weight",
    "adult20", "vcte_complete", "vcte_reliable", "alcohol_classifiable",
    "base_eligible", "diabetes_broad", "target_masld_diabetes", "diabetes_only",
    "hav_valid", "hbv_valid", "hav_susceptible", "hbv_susceptible",
    "hav_joint_indicator", "dual_susceptible", "either_susceptible"
  )
  missing <- setdiff(needed, names(d))
  if (length(missing) > 0L) stop("Current processed file missing columns: ", paste(missing, collapse = ", "))
  d
}

old <- derive_old_period()
current <- derive_current_period()

common <- c(
  "SEQN", "period", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDBORN4",
  "DIQ010", "LBXGH", "BMXBMI", "LUXCAPM", "LUXSMED", "IMQ011", "ALQ111", "ALQ121", "ALQ130",
  "SDMVSTRA", "SDMVPSU", "analysis_weight",
  "adult20", "vcte_complete", "vcte_reliable", "alcohol_classifiable", "alcohol_masld_range",
  "base_eligible", "diabetes_broad", "target_masld_diabetes", "diabetes_only",
  "hav_valid", "hbv_valid", "hav_susceptible", "hbv_susceptible",
  "hav_joint_indicator", "dual_susceptible", "either_susceptible"
)
old <- old[, common]
current <- current[, common]
all_periods <- rbind(old, current)

period_order <- c("2017-March 2020", "August 2021-August 2023")
all_periods$period <- factor(all_periods$period, levels = period_order)

flow_for_period <- function(d) {
  target <- d$base_eligible == 1 & d$target_masld_diabetes == 1
  diab_only <- d$base_eligible == 1 & d$diabetes_only == 1
  data.frame(
    period = as.character(d$period[1]),
    total = nrow(d),
    age_20_plus = sum(d$adult20 == 1, na.rm = TRUE),
    complete_vcte = sum(d$adult20 == 1 & d$vcte_complete == 1, na.rm = TRUE),
    reliable_vcte = sum(d$adult20 == 1 & d$vcte_reliable == 1, na.rm = TRUE),
    base_eligible = sum(d$base_eligible == 1, na.rm = TRUE),
    masld_diabetes = sum(target, na.rm = TRUE),
    masld_diabetes_hav = sum(target & !is.na(d$hav_susceptible), na.rm = TRUE),
    masld_diabetes_hbv = sum(target & !is.na(d$hbv_susceptible), na.rm = TRUE),
    masld_diabetes_joint_indicator = sum(target & !is.na(d$hav_joint_indicator), na.rm = TRUE),
    diabetes_only = sum(diab_only, na.rm = TRUE),
    diabetes_only_hav = sum(diab_only & !is.na(d$hav_susceptible), na.rm = TRUE),
    diabetes_only_hbv = sum(diab_only & !is.na(d$hbv_susceptible), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

flow <- do.call(rbind, lapply(period_order, function(p) flow_for_period(all_periods[all_periods$period == p, ])))
write.csv(flow, file.path(outdir, "sample_flow.csv"), row.names = FALSE, na = "")

make_design <- function(d) {
  universe <- d$adult20 == 1 & !is.na(d$analysis_weight) & d$analysis_weight > 0 &
    !is.na(d$SDMVSTRA) & !is.na(d$SDMVPSU)
  du <- d[universe, ]
  svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~analysis_weight,
    nest = TRUE, data = du
  )
}

designs <- setNames(lapply(period_order, function(p) {
  make_design(all_periods[all_periods$period == p, ])
}), period_order)

design_summary <- data.frame(
  period = period_order,
  adult_weighted_universe_n = vapply(designs, function(x) nrow(x$variables), numeric(1)),
  strata = vapply(designs, function(x) length(unique(x$variables$SDMVSTRA)), numeric(1)),
  psu_pairs = vapply(designs, function(x) nrow(unique(x$variables[c("SDMVSTRA", "SDMVPSU")])), numeric(1)),
  design_df = vapply(designs, degf, numeric(1)),
  stringsAsFactors = FALSE
)
write.csv(design_summary, file.path(outdir, "design_summary.csv"), row.names = FALSE)

endpoint_labels <- c(
  hav_susceptible = "HAV serologic susceptibility",
  hbv_susceptible = "HBV triple-panel susceptibility",
  hav_joint_indicator = "HAV joint indicator",
  dual_susceptible = "Susceptible to both HAV and HBV",
  either_susceptible = "Susceptible to HAV or HBV"
)
group_definitions <- c(
  "MASLD + diabetes" = "target_masld_diabetes",
  "Diabetes only" = "diabetes_only"
)

estimate_one <- function(design, period, group_label, group_var, endpoint, endpoint_label) {
  keep <- design$variables$base_eligible == 1 & design$variables[[group_var]] == 1
  subdesign <- design[keep, ]
  y <- subdesign$variables[[endpoint]]
  denominator <- sum(!is.na(y))
  events <- sum(y == 1, na.rm = TRUE)
  nonevents <- sum(y == 0, na.rm = TRUE)
  if (denominator == 0L || events == 0L || nonevents == 0L) {
    return(data.frame(
      period = period, group = group_label, endpoint = endpoint_label,
      endpoint_variable = endpoint, n = denominator, events = events, nonevents = nonevents,
      prevalence = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      rse = NA_real_, reliability = "Not estimable", stringsAsFactors = FALSE
    ))
  }
  estimate <- svyciprop(as.formula(paste0("~", endpoint)), subdesign, method = "logit", na.rm = TRUE)
  p <- unname(coef(estimate)[1])
  se <- unname(SE(estimate)[1])
  ci <- as.numeric(confint(estimate)[1, ])
  rse <- se / p
  reasons <- character(0)
  if (denominator < 30) reasons <- c(reasons, "n<30")
  if (events < 10) reasons <- c(reasons, "events<10")
  if (nonevents < 10) reasons <- c(reasons, "nonevents<10")
  if (is.finite(rse) && rse > 0.30) reasons <- c(reasons, "RSE>30%")
  reliability <- if (length(reasons) == 0L) "Meets screening criteria" else paste("Use caution:", paste(reasons, collapse = "; "))
  data.frame(
    period = period, group = group_label, endpoint = endpoint_label,
    endpoint_variable = endpoint, n = denominator, events = events, nonevents = nonevents,
    prevalence = p, se = se, ci_low = ci[1], ci_high = ci[2], rse = rse,
    reliability = reliability, stringsAsFactors = FALSE
  )
}

estimate_rows <- list()
k <- 1L
for (p in period_order) {
  for (g in names(group_definitions)) {
    for (ep in names(endpoint_labels)) {
      estimate_rows[[k]] <- estimate_one(
        designs[[p]], p, g, group_definitions[[g]], ep, endpoint_labels[[ep]]
      )
      k <- k + 1L
    }
  }
}
estimates <- do.call(rbind, estimate_rows)
write.csv(estimates, file.path(outdir, "core_endpoint_estimates.csv"), row.names = FALSE, na = "")

# Crude independent-period percentage-point differences. These are a
# crude time-period comparison.
difference_rows <- list()
k <- 1L
for (g in names(group_definitions)) {
  for (ep in names(endpoint_labels)) {
    z <- estimates[estimates$group == g & estimates$endpoint_variable == ep, ]
    z <- z[match(period_order, z$period), ]
    min_df <- min(design_summary$design_df)
    if (nrow(z) == 2L && all(is.finite(z$prevalence)) && all(is.finite(z$se))) {
      diff <- z$prevalence[2] - z$prevalence[1]
      se_diff <- sqrt(sum(z$se^2))
      crit <- qt(0.975, df = min_df)
      low <- diff - crit * se_diff
      high <- diff + crit * se_diff
      p_value <- 2 * pt(-abs(diff / se_diff), df = min_df)
    } else {
      diff <- se_diff <- low <- high <- p_value <- NA_real_
    }
    difference_rows[[k]] <- data.frame(
      group = g, endpoint = endpoint_labels[[ep]], endpoint_variable = ep,
      early_prevalence = z$prevalence[1], current_prevalence = z$prevalence[2],
      absolute_change = diff, se_change = se_diff, ci_low = low, ci_high = high,
      p_value = p_value, comparison_df = min_df,
      reliability = if (all(z$reliability == "Meets screening criteria")) {
        "Both period estimates meet screening criteria"
      } else {
        paste(z$period, z$reliability, collapse = " | ")
      },
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
differences <- do.call(rbind, difference_rows)
write.csv(differences, file.path(outdir, "crude_period_differences.csv"), row.names = FALSE, na = "")

capture.output(sessionInfo(), file = file.path(outdir, "session_info.txt"))
saveRDS(all_periods, file.path(outdir, "harmonized_minimal_analysis.rds"))
