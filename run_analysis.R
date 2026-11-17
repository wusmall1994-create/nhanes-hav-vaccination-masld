#!/usr/bin/env Rscript
# Run from the repository root after downloading and preparing NHANES data.
scripts <- sort(list.files("scripts", pattern = "^[0-9]{2}_.*[.]R$", full.names = TRUE))
stopifnot(length(scripts) == 9L)
rscript <- file.path(R.home("bin"), "Rscript")
for (script in scripts) {
  message("Running ", script)
  status <- system2(rscript, shQuote(script))
  if (status != 0L) stop("Analysis failed: ", script)
}
