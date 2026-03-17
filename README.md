# SANT vs Splenic Lymphoma Analysis

This repository contains an R-based exploratory analysis workflow for differentiating **SANT** and **splenic lymphoma** using clinicopathologic and imaging features.

## Project Overview

- Input dataset: `standard_data.csv`
- Main analysis script: `analysis.R`
- Core outputs:
  - Baseline comparison table (`Table1_Output_Corrected.csv`)
  - Main figures (`figures/main/`)
  - Supplementary figures (`figures/supplementary/`)
  - ROC/OR/bootstrap summaries (`figures/*.csv`)

## How to Run

From the repository root:

```bash
Rscript analysis.R
```

The script will:

1. Install/load required R packages (`table1`, `ggplot2`, `ggpubr`, `pROC`, `dplyr`, etc.).
2. Clean and recode variables (including missing/unknown handling).
3. Generate Table 1 and diagnostic/comparison figures.
4. Export main and supplementary figure PDFs under `figures/`.

## Figure Structure (Current)

- `figures/main/`
  - `Figure2_Significant_Comparisons.pdf`
  - `Figure3_Key_Proportions.pdf`
  - `Figure4_Diagnostic_Framework.pdf`
- `figures/supplementary/`
  - `Supplementary_FigureS1_ROC_Comparison.pdf`
  - `Supplementary_FigureS2_OR_Bootstrap.pdf`
  - `Supplementary_FigureS3_Correlations.pdf`
  - `Supplementary_FigureS4_Clinical_Utility_Exploratory.pdf`

## Notes

- This is a small-sample retrospective-style analysis and should be interpreted as exploratory.
- Model-based utility panels are intentionally placed in supplementary figures to reduce overclaiming risk.
- If Arial is unavailable on your system, plotting automatically falls back to `sans`.

