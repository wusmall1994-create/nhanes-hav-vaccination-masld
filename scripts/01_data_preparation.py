"""Prepare an auditable NHANES 2021-2023 analysis dataset.

This script never modifies raw XPT files. It merges participant-level modules,
derives prespecified MASLD/diabetes/viral-serology variables, and writes a flat
CSV for the R survey-analysis script.
"""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
OUT = ROOT / "outputs" / "masld_hav_hbv"


def read_xpt(name: str, cols: list[str] | None = None) -> pd.DataFrame:
    df = pd.read_sas(RAW / f"{name}.xpt", format="xport")
    if cols is not None:
        missing = sorted(set(cols) - set(df.columns))
        if missing:
            raise KeyError(f"{name}: missing columns {missing}")
        df = df[cols]
    return df


def merge_one(base: pd.DataFrame, name: str, cols: list[str]) -> pd.DataFrame:
    other = read_xpt(name, ["SEQN", *cols])
    if other["SEQN"].duplicated().any():
        raise ValueError(f"{name} is not participant-level")
    return base.merge(other, on="SEQN", how="left", validate="one_to_one")


def yes_no(series: pd.Series, yes=1, no=2) -> pd.Series:
    out = pd.Series(np.nan, index=series.index, dtype="float64")
    out.loc[series == yes] = 1
    out.loc[series == no] = 0
    return out


def main() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    demo_cols = [
        "SEQN", "RIDSTATR", "RIAGENDR", "RIDAGEYR", "RIDRETH3",
        "DMDBORN4", "DMDEDUC2", "DMDMARTZ", "WTMEC2YR", "SDMVSTRA",
        "SDMVPSU", "INDFMPIR",
    ]
    d = read_xpt("DEMO_L", demo_cols)

    modules = {
        "LUX_L": ["LUAXSTAT", "LUAPNME", "LUANMVGP", "LUXSMED", "LUXSIQRM", "LUXCAPM"],
        "HEPA_L": ["WTPH2YR", "LBXHA"],
        "HEPB_S_L": ["LBXHBS"],
        "HEPBD_L": ["LBXHBC", "LBDHBG"],
        "HEPC_L": ["LBXHCR", "LBDHCI"],
        "IMQ_L": ["IMQ011"],
        "DIQ_L": ["DIQ010", "DID040", "DIQ050"],
        "GHB_L": ["LBXGH"],
        "ALQ_L": ["ALQ111", "ALQ121", "ALQ130"],
        "BMX_L": ["BMDSTATS", "BMXBMI", "BMXWAIST"],
        "HIQ_L": ["HIQ011"],
        "HUQ_L": ["HUQ010", "HUQ030", "HUQ042"],
        "SMQ_L": ["SMQ020", "SMQ040"],
        "BIOPRO_L": ["LBXSCR", "LBXSTR"],
        "ALB_CR_L": ["URDACT"],
        "BPXO_L": ["BPXOSY1", "BPXODI1", "BPXOSY2", "BPXODI2", "BPXOSY3", "BPXODI3"],
        "BPQ_L": ["BPQ020", "BPQ150"],
        "HDL_L": ["LBDHDD"],
    }
    for name, cols in modules.items():
        d = merge_one(d, name, cols)

    # Analysis eligibility and VCTE quality.
    d["adult20"] = (d.RIDAGEYR >= 20).astype(int)
    d["vcte_complete"] = (d.LUAXSTAT == 1).astype(int)
    d["vcte_reliable"] = (
        (d.LUAXSTAT == 1)
        & (d.LUANMVGP >= 10)
        & (d.LUXSIQRM <= 30)
        & d.LUXCAPM.notna()
        & d.LUXSMED.notna()
    ).astype(int)

    # Diabetes: diagnosed (not gestational-only per DIQ010 wording) or HbA1c >=6.5%.
    d["diabetes_diagnosed"] = (d.DIQ010 == 1).astype(int)
    d["diabetes_broad"] = ((d.DIQ010 == 1) | (d.LBXGH >= 6.5)).astype(int)
    d["diabetes_duration"] = np.where(
        (d.DIQ010 == 1) & d.DID040.between(1, 80),
        np.maximum(d.RIDAGEYR - d.DID040, 0),
        np.nan,
    )
    d["insulin_use"] = yes_no(d.DIQ050)
    d["hba1c_ge7"] = np.where(d.LBXGH.notna(), (d.LBXGH >= 7).astype(int), np.nan)

    # Average alcohol in g/day. Frequency categories are converted to midpoint days/year.
    days_per_year = {0: 0, 1: 365, 2: 300, 3: 182, 4: 104, 5: 52,
                     6: 30, 7: 12, 8: 9, 9: 4.5, 10: 1.5}
    d["drink_days_year"] = d.ALQ121.map(days_per_year)
    d["alcohol_g_day"] = d["drink_days_year"] * d.ALQ130 * 14 / 365
    nondrinker = (d.ALQ111 == 2) | (d.ALQ121 == 0)
    d.loc[nondrinker, "alcohol_g_day"] = 0
    d["alcohol_classifiable"] = d.alcohol_g_day.notna().astype(int)
    female = d.RIAGENDR == 2
    d["alcohol_masld_range"] = np.where(
        d.alcohol_g_day.notna(),
        np.where(female, d.alcohol_g_day < 20, d.alcohol_g_day < 30).astype(int),
        np.nan,
    )
    d["metald_range"] = np.where(
        d.alcohol_g_day.notna(),
        np.where(female,
                 d.alcohol_g_day.between(20, 50, inclusive="left"),
                 d.alcohol_g_day.between(30, 60, inclusive="left")).astype(int),
        np.nan,
    )

    # Cardiometabolic risk factors from the 2023 MASLD consensus.
    asian = d.RIDRETH3 == 6
    overweight = np.where(asian, d.BMXBMI >= 23, d.BMXBMI >= 25)
    waist_high = np.where(female, d.BMXWAIST >= 80,
                          np.where(asian, d.BMXWAIST >= 90, d.BMXWAIST >= 94))
    d["cmrf_adiposity"] = (overweight | waist_high).astype(int)
    d["cmrf_dysglycemia"] = ((d.LBXGH >= 5.7) | (d.DIQ010.isin([1, 3]))).astype(int)
    bp_cols = ["BPXOSY1", "BPXOSY2", "BPXOSY3"]
    dbp_cols = ["BPXODI1", "BPXODI2", "BPXODI3"]
    d["mean_sbp"] = d[bp_cols].mean(axis=1, skipna=True)
    d["mean_dbp"] = d[dbp_cols].mean(axis=1, skipna=True)
    d["cmrf_bp"] = ((d.mean_sbp >= 130) | (d.mean_dbp >= 85) | (d.BPQ020 == 1)).astype(int)
    d["cmrf_tg"] = (d.LBXSTR >= 150).astype(int)
    d["cmrf_hdl"] = np.where(female, d.LBDHDD <= 50, d.LBDHDD <= 40).astype(int)
    d["cmrf_any"] = d[["cmrf_adiposity", "cmrf_dysglycemia", "cmrf_bp", "cmrf_tg", "cmrf_hdl"]].max(axis=1)

    # Active viral hepatitis. HCV RNA values 2 and 3 both mean no detectable RNA.
    d["active_hbv"] = np.where(d.LBDHBG.isin([1, 2]), (d.LBDHBG == 1).astype(int), np.nan)
    d["active_hcv"] = np.where(d.LBXHCR.isin([1, 2, 3]), (d.LBXHCR == 1).astype(int), np.nan)
    d["no_active_b_or_c"] = np.where(
        d.active_hbv.notna() & d.active_hcv.notna(),
        ((d.active_hbv == 0) & (d.active_hcv == 0)).astype(int),
        np.nan,
    )

    # Main and sensitivity SLD/MASLD definitions. MASLD is an affirmative
    # diagnosis and may coexist with viral hepatitis; active HBV/HCV status is
    # retained for description but is not an exclusion criterion.
    for cutoff in (248, 274, 285):
        d[f"sld_{cutoff}"] = (d.LUXCAPM >= cutoff).astype(int)
        d[f"masld_{cutoff}"] = (
            (d[f"sld_{cutoff}"] == 1)
            & (d.cmrf_any == 1)
            & (d.alcohol_masld_range == 1)
        ).astype(int)
    d["masld"] = d.masld_285
    d["target_masld_diabetes"] = ((d.masld == 1) & (d.diabetes_broad == 1)).astype(int)
    d["target_dx_diabetes"] = ((d.masld == 1) & (d.diabetes_diagnosed == 1)).astype(int)
    d["target_sld_diabetes_all_alcohol"] = (
        (d.sld_285 == 1) & (d.cmrf_any == 1) & (d.diabetes_broad == 1)
    ).astype(int)
    d["phenotype4"] = np.select(
        [
            (d.masld == 0) & (d.diabetes_broad == 0),
            (d.masld == 1) & (d.diabetes_broad == 0),
            (d.masld == 0) & (d.diabetes_broad == 1),
            (d.masld == 1) & (d.diabetes_broad == 1),
        ],
        ["Neither", "MASLD only", "Diabetes only", "MASLD + diabetes"],
        default="Unclassified",
    )

    # HAV outcomes and self-reported vaccination.
    d["hav_valid"] = d.LBXHA.isin([1, 2]).astype(int)
    d["hav_susceptible"] = np.where(d.LBXHA.isin([1, 2]), (d.LBXHA == 2).astype(int), np.nan)
    d["havax_complete"] = np.where(d.IMQ011.isin([1, 2, 3]), (d.IMQ011 == 1).astype(int), np.nan)
    d["havax_incomplete_none"] = np.where(d.IMQ011.isin([1, 2, 3]), d.IMQ011.isin([2, 3]).astype(int), np.nan)
    d["hav_joint_indicator"] = np.where(
        d.hav_susceptible.notna() & d.havax_incomplete_none.notna(),
        ((d.hav_susceptible == 1) & (d.havax_incomplete_none == 1)).astype(int),
        np.nan,
    )

    # HBV serologic profiles.
    hbv_valid = d.LBDHBG.isin([1, 2]) & d.LBXHBC.isin([1, 2]) & d.LBXHBS.isin([1, 2])
    d["hbv_valid"] = hbv_valid.astype(int)
    d["hbv_profile"] = np.select(
        [
            hbv_valid & (d.LBDHBG == 2) & (d.LBXHBC == 2) & (d.LBXHBS == 2),
            hbv_valid & (d.LBDHBG == 2) & (d.LBXHBC == 2) & (d.LBXHBS == 1),
            hbv_valid & (d.LBDHBG == 2) & (d.LBXHBC == 1) & (d.LBXHBS == 1),
            hbv_valid & (d.LBDHBG == 2) & (d.LBXHBC == 1) & (d.LBXHBS == 2),
            hbv_valid & (d.LBDHBG == 1),
        ],
        ["Susceptible", "Vaccine-pattern protection", "Resolved infection", "Isolated anti-HBc", "Current infection"],
        default="Missing/indeterminate",
    )
    d["hbv_susceptible"] = np.where(hbv_valid, (d.hbv_profile == "Susceptible").astype(int), np.nan)
    d["dual_susceptible"] = np.where(
        d.hav_susceptible.notna() & d.hbv_susceptible.notna(),
        ((d.hav_susceptible == 1) & (d.hbv_susceptible == 1)).astype(int),
        np.nan,
    )
    d["either_susceptible"] = np.where(
        d.hav_susceptible.notna() & d.hbv_susceptible.notna(),
        ((d.hav_susceptible == 1) | (d.hbv_susceptible == 1)).astype(int),
        np.nan,
    )

    # Liver, kidney and access-to-care phenotypes.
    d["significant_fibrosis"] = np.where(d.LUXSMED.notna(), (d.LUXSMED >= 8).astype(int), np.nan)
    d["advanced_fibrosis"] = np.where(d.LUXSMED.notna(), (d.LUXSMED >= 12).astype(int), np.nan)
    sex_female = (d.RIAGENDR == 2).astype(int)
    kappa = np.where(sex_female == 1, 0.7, 0.9)
    alpha = np.where(sex_female == 1, -0.241, -0.302)
    scr_ratio = d.LBXSCR / kappa
    d["egfr_2021"] = 142 * np.minimum(scr_ratio, 1) ** alpha * np.maximum(scr_ratio, 1) ** -1.2 * 0.9938 ** d.RIDAGEYR * np.where(sex_female == 1, 1.012, 1)
    d["ckd"] = np.where(
        d.egfr_2021.notna() | d.URDACT.notna(),
        ((d.egfr_2021 < 60) | (d.URDACT >= 30)).astype(int),
        np.nan,
    )
    d["insured"] = yes_no(d.HIQ011)
    d["usual_care"] = np.where(d.HUQ030.isin([1, 2, 3]), d.HUQ030.isin([1, 3]).astype(int), np.nan)
    d["general_health_fairpoor"] = np.where(d.HUQ010.isin([1, 2, 3, 4, 5]), d.HUQ010.isin([4, 5]).astype(int), np.nan)
    d["current_smoker"] = np.where(
        d.SMQ020.isin([1, 2]),
        ((d.SMQ020 == 1) & d.SMQ040.isin([1, 2])).astype(int),
        np.nan,
    )

    # Labels used by survey tables/models.
    d["age_group"] = pd.cut(d.RIDAGEYR, [19, 39, 59, np.inf], labels=["20-39", "40-59", "60+"], right=True).astype("object")
    d["sex"] = d.RIAGENDR.map({1: "Male", 2: "Female"})
    d["race5"] = d.RIDRETH3.map({1: "Hispanic", 2: "Hispanic", 3: "Non-Hispanic White", 4: "Non-Hispanic Black", 6: "Non-Hispanic Asian", 7: "Other/multiracial"})
    d["race4_model"] = d.RIDRETH3.map({1: "Hispanic", 2: "Hispanic", 3: "Non-Hispanic White", 4: "Non-Hispanic Black", 6: "Other", 7: "Other"})
    d["foreign_born"] = np.where(d.DMDBORN4.isin([1, 2]), (d.DMDBORN4 == 2).astype(int), np.nan)
    d["education"] = d.DMDEDUC2.map({1: "<High school", 2: "<High school", 3: "High school/GED", 4: "Some college", 5: "College graduate+"})
    d["pir_group"] = pd.cut(d.INDFMPIR, [-np.inf, 1.3, 3.5, np.inf], labels=["<1.3", "1.3-3.5", ">3.5"], right=False).astype("object")
    d["bmi_group"] = pd.cut(d.BMXBMI, [-np.inf, 25, 30, np.inf], labels=["<25", "25-29.9", "30+"], right=False).astype("object")

    # Define the eligible universe before domain analyses.
    d["base_eligible"] = (
        (d.adult20 == 1)
        & (d.vcte_reliable == 1)
        & (d.WTPH2YR > 0)
        & (d.alcohol_classifiable == 1)
    ).astype(int)

    # Keep raw codes plus derived analysis columns for auditability.
    output_csv = PROCESSED / "masld_hav_hbv_analysis.csv"
    d.to_csv(output_csv, index=False)

    flow = {
        "NHANES total": int(len(d)),
        "Age >=20": int((d.adult20 == 1).sum()),
        "Age >=20 with complete VCTE": int(((d.adult20 == 1) & (d.vcte_complete == 1)).sum()),
        "Age >=20 with reliable VCTE": int(((d.adult20 == 1) & (d.vcte_reliable == 1)).sum()),
        "Base eligible": int((d.base_eligible == 1).sum()),
        "MASLD + diabetes": int(((d.base_eligible == 1) & (d.target_masld_diabetes == 1)).sum()),
        "Target with valid HAV": int(((d.base_eligible == 1) & (d.target_masld_diabetes == 1) & (d.hav_valid == 1)).sum()),
        "Target with valid HBV triple panel": int(((d.base_eligible == 1) & (d.target_masld_diabetes == 1) & (d.hbv_valid == 1)).sum()),
        "Target with both HAV and HBV": int(((d.base_eligible == 1) & (d.target_masld_diabetes == 1) & (d.hav_valid == 1) & (d.hbv_valid == 1)).sum()),
    }
    with open(OUT / "flow_counts.json", "w", encoding="utf-8") as f:
        json.dump(flow, f, ensure_ascii=False, indent=2)

    qc = {
        "rows": int(len(d)),
        "unique_seqn": int(d.SEQN.nunique()),
        "duplicate_seqn": int(d.SEQN.duplicated().sum()),
        "strata": int(d.SDMVSTRA.nunique()),
        "stratum_psu_pairs": int(d[["SDMVSTRA", "SDMVPSU"]].drop_duplicates().shape[0]),
        "weight_positive": int((d.WTPH2YR > 0).sum()),
        "raw_files": sorted(p.name for p in RAW.glob("*.xpt")),
        "definitions": {
            "primary_steatosis": "CAP >=285 dB/m",
            "reliable_vcte": "complete, >=10 valid measures, IQR/median <=30%, nonmissing CAP and LSM",
            "diabetes": "DIQ010=1 or HbA1c >=6.5%",
            "MASLD_alcohol": "<20 g/day women, <30 g/day men",
            "HAV_susceptible": "total anti-HAV negative",
            "HBV_susceptible": "HBsAg-, anti-HBc-, anti-HBs-",
        },
    }
    with open(OUT / "data_qc.json", "w", encoding="utf-8") as f:
        json.dump(qc, f, ensure_ascii=False, indent=2)

    dictionary = [
        ("WTPH2YR", "HEPA_L/HEPBD_L", "Phlebotomy 2-year weight", "Survey weight"),
        ("SDMVSTRA", "DEMO_L", "Masked variance stratum", "Survey design"),
        ("SDMVPSU", "DEMO_L", "Masked variance PSU", "Survey design"),
        ("masld", "Derived", "CAP >=285, >=1 CMRF, alcohol <20 g/day women/<30 men; viral hepatitis may coexist", "Primary exposure"),
        ("diabetes_broad", "DIQ_L/GHB_L", "Doctor-diagnosed diabetes or HbA1c >=6.5%", "Primary exposure"),
        ("hav_susceptible", "HEPA_L", "Total anti-HAV negative (LBXHA=2)", "Primary outcome"),
        ("hbv_susceptible", "HEPBD_L/HEPB_S_L", "HBsAg-, total anti-HBc-, anti-HBs-", "Primary outcome"),
        ("dual_susceptible", "Derived", "Both HAV and HBV serologically susceptible", "Secondary outcome"),
        ("either_susceptible", "Derived", "HAV or HBV serologically susceptible", "Secondary outcome"),
        ("hav_joint_indicator", "HEPA_L/IMQ_L", "Anti-HAV negative plus incomplete/no self-reported HepA vaccination", "Secondary outcome"),
        ("significant_fibrosis", "LUX_L", "LSM >=8 kPa", "Clinical correlate"),
        ("ckd", "BIOPRO_L/ALB_CR_L", "2021 CKD-EPI eGFR <60 or UACR >=30 mg/g", "Clinical correlate"),
        ("insured", "HIQ_L", "Covered by health insurance", "Access correlate"),
        ("usual_care", "HUQ_L", "Has one or more usual places for health care", "Access correlate"),
    ]
    pd.DataFrame(dictionary, columns=["variable", "source_file", "definition", "role"]).to_csv(
        OUT / "variable_dictionary.csv", index=False
    )

    manifest_rows = []
    for path in sorted(RAW.glob("*.xpt")):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        stem = path.stem
        manifest_rows.append({
            "file": path.name,
            "bytes": path.stat().st_size,
            "sha256": digest,
            "source_url": f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/{stem}.xpt",
        })
    pd.DataFrame(manifest_rows).to_csv(OUT / "source_manifest.csv", index=False)

    print(json.dumps(flow, ensure_ascii=False, indent=2))
    print(f"Wrote {output_csv}")


if __name__ == "__main__":
    main()
