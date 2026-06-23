[README.md](https://github.com/user-attachments/files/29236355/README.md)
# R Scripts for Reproducing the Statistical Analysis

**Manuscript:** Asrat Z, Fekadu M, Woldu Z, Demissew S, Opgenoorth L, Miehe G, Zech W. *Distribution and Regeneration Potential of the Soil Seed Bank of Ericaceous Vegetation in the Bale Mountains.* Journal of Vegetation Science (2026).

**Zenodo DOI:** https://doi.org/10.5281/zenodo.10251925

---

## Overview

This repository contains all R scripts needed to reproduce the statistical analyses, ordinations, and figures presented in the manuscript. The scripts are organised into four numbered files that should be run in sequence.

## Directory Structure

```
R_scripts/
├── README.md                          ← This file
├── 01_diversity_ordination.R          ← Diversity indices, NMDS, CCA, PERMANOVA, ANOSIM
├── 02_glmm_models.R                   ← GLMMs for Erica seed density & germination (glmmTMB)
├── 03_anova_posthoc.R                 ← ANOVA, Tukey HSD, normality tests, Figure 4/5/6 regressions
├── 04_figure_generation.R             ← ggplot2 figures (Figs. 2, 4, 5, 6)
└── data/                              ← CSV exports from the original .xls dataset
    ├── Plot_data.csv
    ├── TOTAL_SPECIES_LIST.csv
    ├── SOIL_seedbank-species_list.csv
    ├── Aboveground-SPP.csv
    ├── Soil-Seedbank_3_column_data.csv
    ├── Aboveground(stand)_3_column_dat.csv
    ├── SOIL_seed-germ_dat.csv
    ├── Environmental_data.csv
    └── Germination_count_per_10_days_.csv
```

## Software Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| **R** | 4.3.0+ | Statistical computing |
| RStudio (recommended) | 2023.06+ | IDE |

## R Package Dependencies

Install all required packages before running the scripts:

```r
install.packages(c("vegan", "glmmTMB", "lme4", "MuMIn", "DHARMa",
                   "dplyr", "tidyr", "readr", "ggplot2", "patchwork", "pwr"))
```

| Package | Version | Used in Script | Purpose |
|---------|---------|-----------------|---------|
| `vegan` | 2.6-4 | 01, 04 | Diversity indices, NMDS, CCA, PERMANOVA, ANOSIM |
| `glmmTMB` | 1.1.8+ | 02 | Generalized Linear Mixed Models (negative binomial) |
| `lme4` | 1.1-34+ | 02 | Mixed-effects models support |
| `MuMIn` | 1.47.5+ | 02 | AICc model selection |
| `DHARMa` | 0.4.6+ | 02 | Residual diagnostics for GLMMs |
| `dplyr` | 1.1.2+ | All | Data manipulation |
| `tidyr` | 1.3.0+ | 01, 04 | Reshaping data (pivot_wider) |
| `readr` | 2.1.4+ | All | Reading CSV files |
| `ggplot2` | 3.4.2+ | 04 | Figure generation |
| `patchwork` | 1.1.2+ | 04 | Combining multi-panel figures |
| `pwr` | 1.3-2+ | 02 | Post-hoc power analysis |

## How to Reproduce the Analysis

### Step 1: Set up the working directory

```r
# Set the working directory to the R_scripts folder
setwd("/path/to/R_scripts")

# Verify data files are present
list.files("data/")
```

### Step 2: Run the scripts in order

```r
# Script 1: Diversity and ordination analysis
source("01_diversity_ordination.R")
# Output: alpha_diversity.csv, alpha_diversity_summary_by_PES.csv

# Script 2: GLMM models
source("02_glmm_models.R")
# Output: glmm_erica_cover_coefficients.csv, model_comparison CSVs

# Script 3: ANOVA and post-hoc tests
source("03_anova_posthoc.R")
# Output: figure4_regression_results.csv

# Script 4: Figure generation
source("04_figure_generation.R")
# Output: PNG figures in figures/ directory (600 DPI)
```

### Step 3: Check outputs

```r
# Generated files
list.files("output/")
list.files("figures/")
```

## Script Details

### 01_diversity_ordination.R

**What it does:**
- Computes alpha diversity (species richness S, Shannon H', Pielou's evenness J') for standing vegetation and germinated seed bank
- Calculates Sørensen similarity between standing vegetation and seed bank
- Performs Non-metric Multidimensional Scaling (NMDS) with Bray-Curtis dissimilarity
- Performs Canonical Correspondence Analysis (CCA) with environmental variables
- Runs PERMANOVA and ANOSIM tests for compositional differences among Principal Experimental Sites (PES)

**Key outputs:**
- Alpha diversity table (per-plot and summary by PES)
- NMDS stress value and ANOSIM R/p statistics
- CCA eigenvalues and constrained inertia
- PERMANOVA F-statistic and R²

**Reproduces:** Results §3.1 (vegetation and seed bank composition), §3.3 (environmental and disturbance effects), Table 1

### 02_glmm_models.R

**What it does:**
- Fits negative binomial GLMMs for Erica cover (proxy for seed density) and germination potential
- Fixed effects: elevation, disturbance intensity, stoniness, moss cover
- Random effect: PES (principal experimental site)
- Performs model selection via AICc (dredge-style comparison of candidate models)
- Compares linear vs. quadratic polynomial regressions for elevation
- Conducts post-hoc power analysis using the `pwr` package

**Key outputs:**
- GLMM coefficient tables (estimate, SE, z, p)
- AICc model comparison tables
- Post-hoc power estimates

**Reproduces:** Results §3.2 (seed bank density and germination), Table 3, post-hoc power analysis

### 03_anova_posthoc.R

**What it does:**
- Performs Shapiro-Wilk normality tests on all response variables
- Runs one-way ANOVA comparing means across PES for all response variables
- Performs Tukey's HSD post-hoc pairwise comparisons
- Runs Welch's t-tests comparing standing vegetation vs seed bank diversity
- Fits linear regressions for Figure 4 (Erica cover vs moss, Erica cover vs disturbance, Erica seedlings vs disturbance, Erica germination vs disturbance)
- Fits cubic regressions for Figure 5 (elevational patterns)

**Key outputs:**
- ANOVA F and p-values
- Tukey HSD adjusted p-values for all pairwise comparisons
- R² and p-values for all four Figure 4 panels
- Figure 4 regression results table

**Reproduces:** Table 1, Table 2, Table 4, Figure 4 statistics, Figure 5 statistics, Figure 6 statistics

### 04_figure_generation.R

**What it does:**
- Generates Figure 2 (temporal emergence patterns) using ggplot2 with LOESS smoothing
- Generates Figure 4 (4-panel environmental drivers) with regression fits and R²/p annotations
- Generates Figure 5 (4-panel elevational patterns) with polynomial fits
- Generates Figure 6 (2-panel seed density and germination vs gradients)
- All figures saved at 600 DPI in PNG format

**Key outputs:**
- `Figure2_temporal_emergence.png`
- `Figure4_environmental_drivers.png` (or 4 individual panels)
- `Figure5_elevational_patterns.png` (or 4 individual panels)
- `Figure6_seed_density_germination.png` (or 2 individual panels)

**Note:** Figure 1 (study area map) and Figure 3 (CCA/NMDS biplots) are generated separately — Figure 1 requires GIS software, and Figure 3 ordination biplots are produced interactively in R using `plot()` on the ordination objects from Script 01.

## Data Description

| File | Description | Rows × Cols |
|------|-------------|-------------|
| `Plot_data.csv` | Plot metadata (plot number, elevation, site, GPS coordinates, locality) | 50 × 6 |
| `TOTAL_SPECIES_LIST.csv` | Complete species list with family, habit, life style, endemicity | 113 × 6 |
| `Aboveground-SPP.csv` | Species list observed in standing vegetation | 88 × 2 |
| `SOIL_seedbank-species_list.csv` | Species list observed in soil seed bank germination | 73 × 2 |
| `Soil-Seedbank_3_column_data.csv` | Plot × species abundance for seed bank | 310 × 4 |
| `Aboveground(stand)_3_column_dat.csv` | Plot × species abundance for standing vegetation | 919 × 4 |
| `SOIL_seed-germ_dat.csv` | Per-plot germination data (moss, lichen, Erica germination, total germination) | 50 × 5 |
| `Environmental_data.csv` | Per-plot environmental variables (elevation, slope, soil, rock, deadwood, litter, Erica tree/seedling, shrub cover, moss, lichen, disturbance) | 50 × 13 |
| `Germination_count_per_10_days_.csv` | Per-plot germination counts at 10-day intervals (14 intervals) | 50 × 17 |

## Random Seed

All stochastic analyses use `set.seed(42)` for full reproducibility.

## Citation

If you use these scripts or data, please cite:

> Asrat Z, Fekadu M, Woldu Z, Demissew S, Opgenoorth L, Miehe G, Zech W. Distribution and Regeneration Potential of the Soil Seed Bank of Ericaceous Vegetation in the Bale Mountains. *Journal of Vegetation Science* (2026). doi: [pending]

## License

- **Data:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **Code:** MIT License

## Contact

**Corresponding author:** Zerihun Asrat (zerihunasrat@wku.edu.et)

**Data repository:** https://doi.org/10.5281/zenodo.20694321
