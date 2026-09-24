# Replication Materials

## Scale-dependent predictors in spatial criminology: random forest evidence from street segments to census tracts

---

## Overview

This directory contains all data, code, and manuscript files required to replicate the analyses reported in the paper. All models are estimated in R using the `ranger` package with a fixed random seed (`set.seed(42)`).

The replication procedure has two paths:

- **Full replication** — re-runs all models from the raw data files. Recommended for complete verification. Scripts 02 and 03 were run on an HPC cluster and have approximately 2–4 hours total runtime; replicators running locally will need to update the file paths at the top of each script (see note below).
- **Results verification only** — loads the pre-computed results files and verifies that the numbers in the manuscript tables match. Suitable for reviewers who wish to confirm reported results without re-running the models.

---

## Directory Structure

```
replication/
├── README.md                                      # this file
├── data/
│   ├── segments_final_RF_ready_aligned.csv        # segment-level data (n = 14,188; row order matches gpkg geometry)
│   ├── segments_aggregated_DA_2021.csv            # DA-level aggregated data (n = 1,009)
│   ├── segments_aggregated_CT_2021.csv            # CT-level aggregated data (n = 127)
│   └── W_segs_touches.rds                         # pre-built network-topology spatial weights (segment level)
├── code/
│   ├── 01_RF_Model_Data_Analysis.R                # local data inspection and predictor verification
│   ├── 02_RF_analysis_SLAGS.R                     # segment-level RF models (cluster script)
│   ├── 03_RF_analysis_CTs_DAs_SLAGS.R             # DA- and CT-level RF models (cluster script)
│   ├── 04_ALE_Plots.R                             # ALE plots (Figures 3 and 4)
│   └── 05_NB_ZINB_comparison.R                    # negative binomial comparison (Appendix S1)
├── results/
│   ├── RF_count_SLAGS/
│   │   ├── RF_summary.csv                         # segment-level R² results (Table 2)
│   │   ├── varimp_placeOnly_[CRIME].csv           # segment-level variable importance, place-only models
│   │   ├── varimp_withLags_[CRIME].csv            # segment-level variable importance, place+lags models
│   │   └── oob_convergence_TFV_2025.png           # convergence diagnostic
│   ├── RF_count_DA_SLAGS/
│   │   ├── RF_summary_DA.csv                      # DA-level R² results (Table 2)
│   │   ├── varimp_DA_placeOnly_[CRIME].csv
│   │   └── varimp_DA_withLags_[CRIME].csv
│   ├── RF_count_CT_SLAGS/
│   │   ├── RF_summary_CT.csv                      # CT-level R² results (Table 2)
│   │   ├── varimp_CT_placeOnly_[CRIME].csv
│   │   └── varimp_CT_withLags_[CRIME].csv
│   ├── NB_comparison/
│   │   ├── NB_ZINB_full_coefficients.csv          # all NB/ZINB coefficients (Appendix S1)
│   │   └── NB_ZINB_summary.csv                    # summary: preferred model, W×y β and p (Table S1)
│   ├── Figure3_ALE_crimemix2025_segment_TOB.png   # Figure 3
│   └── Figure4_ALE_MedIncome_CT_BNER.png          # Figure 4
└── manuscript/
    └── [manuscript file]
```

---

## Software Requirements

```r
install.packages(c("ranger", "randomForest", "dplyr", "ggplot2", "iml", "scales",
                   "MASS", "pscl", "spdep", "sf"))
```

| Package | Version tested | Purpose |
|---|---|---|
| `ranger` | 0.16+ | Random forest estimation (all models) |
| `randomForest` | 4.7+ | OOB convergence diagnostic only |
| `dplyr` | 2.0+ | Data manipulation |
| `ggplot2` | 3.4+ | Variable importance plots and Figures 1–2 |
| `iml` | 0.11+ | ALE plots (Figures 3 and 4) |
| `scales` | 1.3+ | Y-axis formatting in Figures 1–2 |
| `MASS` | 7.3+ | Negative binomial regression (`glm.nb`) |
| `pscl` | 1.5+ | Zero-inflated NB (`zeroinfl`) and Vuong test |
| `spdep` | 1.3+ | Spatial weights and spatial lag of outcome |
| `sf` | 1.0+ | Boundary polygon geometry |

R version 4.3 or later is recommended. The segment-level models benefit from at least 8 GB of RAM.

---

## Key Model Settings

All RF models use the following settings, fixed in the code and must not be changed for results to replicate exactly:

| Setting | Value |
|---|---|
| Random seed | `set.seed(42)` |
| Number of trees | 1000 |
| Computational threads | 1 (`num.threads = 1`) |
| Variable importance | Permutation-based |
| Missing value handling | Median imputation |
| Outcome year | 2025 |
| Crime types modelled | TOB, BNEC, TOV, BNER, Theft, Mischief, TFV |

> **Important:** `ranger` produces non-identical results across repeated fits when run with multiple threads, even with the same seed. 

---

## Note on File Paths

Scripts 02 and 03 were run on the HPC cluster and use Linux paths (e.g., `/home/yourname/data/`, `/home/yourname/results/`). Replicators running locally must update the `data_path` and `results_path` variables at the top of each script to point to their local copies of the data files before running. Script 01 uses Windows paths and is provided for reference only — it does not fit any models.

---

## Full Replication Procedure

Run the scripts in the following order.

### Step 1 — Inspect data and verify predictor groups (optional)

```r
source("code/01_RF_Model_Data_Analysis.R")
```

Loads `segments_final_RF_ready_aligned.csv`, constructs predictor groups (network, diversity, business, census, spatial context variables), applies median imputation, and prints a variable group summary to the console. This script does not fit any models — it is provided for data inspection and to document the predictor construction logic used in Script 02.

**Update required:** change `data_dir` at the top of the script to your local data directory.

**Expected console output:** variable group counts 600+ predictors across place-only specification.

---

### Step 2 — Segment-level random forest models

```r
source("code/02_RF_analysis_SLAGS.R")
```

Fits 14 random forest regression models (7 crime types × 2 specifications) at the street segment level (n = 14,188). Begins with an OOB convergence diagnostic for TFV to confirm the tree count is sufficient, then fits two models per crime type:

- **Place-characteristics-only:** network + diversity + business + census + spatial context predictors 
- **Place + lags:** above plus 22 years of lagged crime counts, 2003–2024 

**Update required:** change `data_path` and `results_path` at the top of the script to your local paths.

**Output files written to** `results/RF_count_SLAGS/`:
- `RF_summary.csv` — OOB R² for all 14 models (populates Table 2, segment columns)
- `varimp_placeOnly_[CRIME].csv` — permutation importance for all predictors, place-only model
- `varimp_withLags_[CRIME].csv` — permutation importance for all predictors, place+lags model
- `oob_convergence_TFV_2025.png` — convergence diagnostic plot
- `varimp_placeOnly_[CRIME].png`, `varimp_withLags_[CRIME].png` — importance bar charts

**Expected runtime:** approximately 2–3 hours single-threaded.

---

### Step 3 — DA- and CT-level random forest models

```r
source("code/03_RF_analysis_CTs_DAs_SLAGS.R")
```

Fits 28 random forest regression models (7 crime types × 2 specifications × 2 scales: DA and CT) using the same pipeline as Script 02 applied to the aggregated DA and CT files. The random seed is reset to `set.seed(42)` at the start of each scale.

**Update required:** change `data_path` and `results_path` within the `LEVELS` list at the top of the script to your local paths.

**Output files written to** `results/RF_count_DA_SLAGS/` and `results/RF_count_CT_SLAGS/`:
- `RF_summary_DA.csv`, `RF_summary_CT.csv` — OOB R² (populates Table 2, DA and CT columns)
- `varimp_DA_placeOnly_[CRIME].csv`, `varimp_DA_withLags_[CRIME].csv`
- `varimp_CT_placeOnly_[CRIME].csv`, `varimp_CT_withLags_[CRIME].csv`

**Expected runtime:** approximately 20–40 minutes.

---

### Step 4 — ALE plots (Figures 3 and 4)

```r
source("code/04_ALE_Plots.R")
```

Refits two models independently of Scripts 02–03:

- **Figure 3:** `diversity_25` at the segment level predicting TOB (rank 1 in the place-only model)
- **Figure 4:** leading census-category predictor at the CT level predicting BNER (`CT_med_income_ln_2006`)

The script prints the top overall predictor and top census-category predictor to the console before plotting. Confirm these match the manuscript before treating the plots as final — this is a fresh model fit, not a saved object.

**Update required:** update input data paths and output PNG paths at the top of the script.

**Output files:**
- `results/Figure3_ALE_crimemix2025_segment_TOB.png` — Figure 3
- `results/Figure4_ALE_MedIncome_CT_BNER.png` — Figure 4

**Expected runtime:** approximately 10–20 minutes.

---

### Step 5 — Negative binomial comparison (Appendix S1)

```r
source("code/05_NB_ZINB_comparison.R")
```

Estimates negative binomial (NB) and zero-inflated NB (ZINB) regression models using the top-five predictors identified by each RF place-only model, at all three spatial scales and for all seven crime types. A spatial lag of the outcome (W×y) is included at every scale: network-topology weights (pre-built from segment adjacency via `st_touches()` on the gpkg geometries) at the street segment scale, and first-order Queen's-contiguity weights at the DA and CT scales. The Vuong test selects between NB and ZINB per model (NB preferred throughout). Predictors are standardised to unit variance prior to estimation.

**Update required:** update `data_path`, `varimp_path`, `weights_path`, `BOUNDARY_GPKG`, and `OUT_DIR` within the `LEVELS` list at the top of the script.

**Output files written to** `results/NB_comparison/`:
- `NB_ZINB_full_coefficients.csv` — all coefficients for all models
- `NB_ZINB_summary.csv` — per-model summary: preferred model, top-predictor significance counts, W×y β and p (populates Table S1)

**Expected runtime:** approximately 5–10 minutes.

---

## Results Verification Only (Without Re-Running Models)

To verify that the pre-computed results match the manuscript tables without re-running the models:

```r
# Load pre-computed results
seg  <- read.csv("results/RF_count_SLAGS/RF_summary.csv")
da   <- read.csv("results/RF_count_DA_SLAGS/RF_summary_DA.csv")
ct   <- read.csv("results/RF_count_CT_SLAGS/RF_summary_CT.csv")
nb   <- read.csv("results/NB_comparison/NB_ZINB_summary.csv")

# Table 2 — OOB R² at all three scales
print(seg[, c("crime_type", "rsq_place_only", "rsq_with_lags")])
print(da[,  c("crime_type", "rsq_place_only", "rsq_with_lags")])
print(ct[,  c("crime_type", "rsq_place_only", "rsq_with_lags")])

# Table S1 — NB comparison summary
print(nb[, c("scale", "crime", "n_top_sig_p05", "slag_estimate", "slag_p_value", "slag_sig")])

# Tables 3–4 — variable importance (example: TOB place-only, segment level)
vi <- read.csv("results/RF_count_SLAGS/varimp_placeOnly_TOB.csv")
print(head(vi[order(-vi$importance), ], 10))
```

All pre-computed results files were generated using the scripts and settings documented above.

---

## Notes on Reproducibility

- **Tree-count convergence** was verified for all models via OOB error diagnostics at each spatial scale prior to publication. The convergence diagnostic in Script 02 (OOB MSE by number of trees for TFV) can be run independently to confirm this.
- **CT-level category stability check** (Table 5) was conducted by refitting the CT place-only model 500 times with independent random seeds and recording which predictor category obtained the higher importance share in each refit. Results are provided as pre-computed summaries.
- **Spatial weights:** The segment-level network-topology weights (`W_segs_touches.rds`) were constructed using `st_touches()` on the segment geometries from the Vancouver road network GeoPackage, identifying segments that share an intersection node. The segment data file (`segments_final_RF_ready_aligned.csv`) has been reordered to match the row order of the GeoPackage geometry, which is required for the pre-built weights to align correctly with the data. DA- and CT-level Queen's-contiguity weights are constructed on the fly in Scripts 03 and 05 from the boundary GeoPackage.
- **Data sources:** crime data were obtained from the Vancouver Police Department Open Data Portal; census sociodemographic variables from Statistics Canada (2001–2021 Censuses). The data files included here are pre-processed derivatives. Raw data sources and processing steps are described in full in the Data and Methods section of the manuscript.

---

## Contact

For questions about the replication materials, please contact the corresponding author.
