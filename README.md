# Depression-Independent Thalamocortical Correlates of Post-Stroke Fatigue

Analysis code for:

> **Depression-Independent Thalamocortical Correlates of Post-Stroke Fatigue**  

---

## Overview

This repository contains R analysis code for a longitudinal resting-state fMRI study examining neural correlates of fatigue in stroke survivors and age-matched controls. The primary contribution is a depression dissociation approach: fatigue-connectivity associations are tested both in depression-adjusted models and after residualizing SF-36 Vitality scores with respect to depressive symptom burden (GDSS).

## Repository Structure

```
├── 01_FC_ALFF_fALFF_compiler.mat  # Compile imaging features
├── 02_merge_and_demographics.R    # Merge REDCap + imaging data; demographics table
├── 03_mixed_models_FC.R           # Mixed-effects models: FC ~ vitality + GDSS
├── 04_mixed_models_ALFF_fALFF.R   # Mixed-effects models: ALFF/fALFF ~ vitality + GDSS
├── 05_parcel_group_tvals_network.R        # Parcel-wise stroke vs control group differences (Network)
├── 06_parcel_group_tvals_ROI.R        # Parcel-wise stroke vs control group differences (ROI)
├── 07_figures.R                   # Figures: scatter plots, network heatmap, FC node strength surface maps
├── 08_lesion_site_sensitivity.R   # Lesion site moderation analyses
├── 09_lesion_size_moderation.R    # Lesion volume moderation analyses
├── 10_supplementary_tables.R      # Supplementary Excel tables (S1–S5)
└── README.md
```

---

## Requirements

**R version:** 4.3 or later

**Required packages:**

| Package | Version | Purpose |
|---------|---------|---------|
| lme4 | ≥ 1.1 | Linear mixed-effects models |
| lmerTest | ≥ 3.1 | p-values for mixed models |
| broom.mixed | ≥ 0.2 | Tidy model outputs |
| dplyr | ≥ 1.1 | Data manipulation |
| readr | ≥ 2.1 | CSV I/O |
| readxl | ≥ 1.4 | REDCap Excel import |
| stringr | ≥ 1.5 | String processing |
| ggplot2 | ≥ 3.4 | Plotting |
| cowplot | ≥ 1.1 | Figure assembly |
| patchwork | ≥ 1.1 | Panel layout |
| ggseg | ≥ 1.6 | Brain surface visualization |
| sf | ≥ 1.0 | Spatial features for ggseg |
| openxlsx | ≥ 4.2 | Supplementary tables |
| scales | ≥ 1.2 | Color scaling |

Install all at once:
```r
install.packages(c(
  "lme4", "lmerTest", "broom.mixed", "dplyr", "readr", "readxl",
  "stringr", "ggplot2", "cowplot", "patchwork", "sf", "openxlsx", "scales"
))

# ggseg requires separate installation
install.packages("ggseg", repos = c("https://ggseg.r-universe.dev",
                                     "https://cloud.r-project.org"))
```

---

## Data

Participant data are not publicly available due to confidentiality requirements.

**Input files required:**
| File | Description |
|------|-------------|
| `data/FCS_Demographics_and_Behavior.xlsx` | REDCap export with demographics, SF-36 Vitality, and GDSS scores |
| `data/regionFC_and_ChaCo_alff_falff_features_long.csv` | Long-format imaging features (FC, ALFF, fALFF, ChaCo) |

**Expected data format for imaging CSV:**

| Column | Description |
|--------|-------------|
| study_id | Participant identifier |
| event_name | REDCap event (e.g. acute_arm_1, 3month_arm_1) |
| feature_type | Feature family (e.g. NET_BETWEEN_FC, ALFF_ROI) |
| feature_name | Specific feature (e.g. LIM-FP, ctx-lh-insula) |
| value | Feature value (Fisher-z for FC, raw for ALFF/fALFF) |

---

## Usage

-Run scripts in order

**Before running:** edit the `CONFIG` section at the top of each script to set your local data and output paths.

---

## MRI Preprocessing

Preprocessing was performed using:
- **fMRIPrep** v23.2.1 (Esteban et al., 2019; doi:10.1038/s41592-018-0235-4)
- **FreeSurfer** v7.3.2 (Dale et al., 1999)
- **ICA-AROMA** via fMRIPost-AROMA v0.0.12 (Pruim et al., 2015; doi:10.1016/j.neuroimage.2015.02.064)
- Spatial smoothing: 6 mm FWHM (SUSAN; Smith & Brady, 1997)
- Standard space: MNI152NLin6Asym (Evans et al., 2012)

Functional network assignment followed the 7-network parcellation of Yeo et al. (2011; doi:10.1152/jn.00338.2011).

---

## License

MIT License. See `LICENSE` for details.

---

## Contact

Joan Stilling 
Weill Cornell Medicine
