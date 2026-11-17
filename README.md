# Joint distribution of hepatitis A serology and vaccination in NHANES

Analysis code for the repeated cross-sectional study of US adults with diabetes and probable MASLD. The primary analysis reports the eight-cell distribution of total anti-HAV and self-reported hepatitis A vaccine doses separately for NHANES 2017–March 2020 and August 2021–August 2023. Unknown dose histories are retained.

## Reproduce the analysis

Run from the repository root with Python 3 (dependencies in `requirements.txt`) and R with packages `haven`, `survey`, `ggplot2`, and `svglite`. The figure scripts use Arial.

```sh
python scripts/00_download_nhanes.py
python scripts/01_data_preparation.py
Rscript run_analysis.R
```

The downloader retrieves public CDC/NCHS transport files and records URLs and SHA-256 hashes. Inputs, generated datasets and outputs are local. `outputs/joint/` contains the primary distribution, dose-conditional estimates, age standardization, phenotype sensitivity analyses and exploratory descriptions. `outputs/secondary/` and `outputs/masld_hav_hbv_exploratory/` contain the secondary estimates and model outputs. Figures are written to `outputs/figures/`.
