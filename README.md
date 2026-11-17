# Joint distribution of hepatitis A serology and vaccination in NHANES

Analysis code for the repeated cross-sectional study of US adults with diabetes and probable MASLD. The primary analysis reports the eight-cell distribution of total anti-HAV and self-reported hepatitis A vaccine doses separately for NHANES 2017閳ユ彈arch 2020 and August 2021閳ユ弬ugust 2023. Unknown dose histories are retained.

## Reproduce the analysis

Run from the repository root with Python 3 (dependencies in `requirements.txt`) and R with packages `haven`, `survey`, `ggplot2`, and `svglite`. The figure scripts use Arial.

```sh
python scripts/00_download_nhanes.py
python scripts/01_data_preparation.py
Rscript run_analysis.R
```

The downloader retrieves public CDC/NCHS transport files and records URLs and SHA-256 hashes. Inputs, generated datasets and outputs are local. `outputs/joint/` contains the primary distribution, dose-conditional estimates, age standardization, phenotype sensitivity analyses and exploratory descriptions. `outputs/secondary/` and `outputs/masld_hav_hbv_exploratory/` contain the secondary estimates and model outputs. Figures are written to `outputs/figures/`.

The shared-low-alcohol CAP comparison in Table 4 uses all adults with valid HAV serology and known birthplace, retaining unknown vaccination histories. Three survey-weighted logistic models include period, CAP group, their interaction, a natural age spline (3 df), and birthplace. Predictive margins use the equally weighted period-specific age/birthplace distributions in this same eligible comparison sample. Absolute differences and 95% CIs retain model covariance and use 33 residual design degrees of freedom. The analysis writes cap_joint_comparison.csv, model diagnostics and standard-population summaries under outputs/joint/.

The external age sensitivity uses the 2000 US standard population (Census P25-1130) from https://seer.cancer.gov/stdpopulations/stdpop.19ages.html, aggregated to ages 20–39, 40–59 and 60+ and renormalized among adults aged 20+. Source counts and the resulting fixed weights are explicit in scripts/12_us2000_age_standard.R. It writes the full joint distribution, age-specific denominators, event counts and intervals, and standard weights. The internal empirical standard supports comparisons within this study. The external standard examines dependence on a younger reference age structure. These estimates describe the study phenotype and retain all valid HAV results, including unknown histories.
