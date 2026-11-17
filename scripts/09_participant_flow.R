#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
theme_set(theme_classic(base_family="Arial",base_size=10))
args <- commandArgs(trailingOnly = TRUE)
root <- "."
feas_dir <- file.path(root, "outputs", "harmonized")
formal_dir <- file.path(root, "outputs", "secondary")
fig_dir <- "outputs/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

flow <- read.csv(file.path(feas_dir, "participant_flow.csv"), check.names = FALSE)

draw_flow <- function() {
  periods <- c("2017-March 2020", "August 2021-August 2023")
  x_positions <- c(0.26, 0.74)
  y_positions <- c(0.83, 0.685, 0.54, 0.395, 0.25, 0.105)
  labels <- lapply(seq_len(2), function(i) c(
    sprintf("NHANES participants\nn = %s", format(flow$total[i], big.mark = ",")),
    sprintf("Age 20 years or older\nn = %s", format(flow$age_20_plus[i], big.mark = ",")),
    sprintf("Reliable VCTE\nn = %s", format(flow$reliable_vcte[i], big.mark = ",")),
    sprintf("Base-eligible adults: n = %s\nInvalid weights excluded: %s\nUnclassifiable alcohol intake excluded: %s", format(flow$base_eligible[i], big.mark = ","), flow$invalid_weight[i], flow$alcohol_missing[i]),
    sprintf("Probable MASLD + diabetes\nn = %s", format(flow$masld_diabetes[i], big.mark = ",")),
    sprintf(
      "Primary joint distribution: %s\nComplete HBV panel: %s\nJoint indicator with known history: %s",
      format(flow$masld_diabetes_hav[i], big.mark = ","),
      format(flow$masld_diabetes_hbv[i], big.mark = ","),
      format(flow$masld_diabetes_joint_indicator[i], big.mark = ",")
    )
  ))

  grid.newpage()
  for (i in seq_along(periods)) {
    grid.text(periods[i], x = x_positions[i], y = 0.95,
              gp = gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold", col = "#16324F"))
    for (j in seq_along(y_positions)) {
      fill <- if (j == 5) "#DDEBF7" else if (j == 6) "#FFF2CC" else "#F6F8FA"
      grid.roundrect(x = x_positions[i], y = y_positions[j], width = 0.39, height = 0.105,
                     r = unit(0.02, "snpc"), gp = gpar(fill = fill, col = "#5B6770", lwd = 1.1))
      grid.text(labels[[i]][j], x = x_positions[i], y = y_positions[j],
                gp = gpar(fontfamily = "Arial", fontsize = if (j == 6) 8.5 else 9.3, col = "#202124"))
      if (j < length(y_positions)) {
        grid.lines(x = unit(c(x_positions[i], x_positions[i]), "npc"),
                   y = unit(c(y_positions[j] - 0.055, y_positions[j + 1] + 0.055), "npc"),
                   arrow = arrow(type = "closed", length = unit(0.11, "inches")),
                   gp = gpar(col = "#6B7280", lwd = 1.0))
      }
    }
  }
}

png(file.path(fig_dir, "Figure_1.png"), width = 3000, height = 3300, res = 400, bg = "white")
draw_flow()
dev.off()

svglite::svglite(file.path(fig_dir, 'Figure_1.svg'), width=7.5, height=8.25)
draw_flow()
dev.off()

tiff(file.path(fig_dir, "Figure_1.tiff"), width = 7.5, height = 8.25, units = "in", res = 600,
     compression = "lzw", bg = "white")
draw_flow()
dev.off()
