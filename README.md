# Urban Emission Policies and Transport CO2 — a Double Machine Learning Analysis

Empirical project for the course *Machine Learning in Econometrics*. The
analysis estimates the causal effect of two urban transport policies —
**congestion pricing (CP)** and **low-emission zones (LEZ)** — on city-level
road-transport CO2 emissions, using **Double/Debiased Machine Learning (DML)**
for the partially linear model with city fixed effects (`DoubleML` package,
`mlr3` learners).

The full task description is in [`instructions.pdf`](instructions.pdf); the
research questions addressed are:

1. **Q1 — Average effects:** How do CP and LEZ each affect transport CO2, and
   what is the joint effect of running both?
2. **Q2 — Heterogeneity:** Do the effects vary with city characteristics
   (population, GDP, transit quality, ...)?
3. **Q3 — Interaction:** Are the two policies super- or sub-additive when
   combined?

**Data:** an anonymized panel of European cities (city-year observations) in
`Data/urban_emissions_panel.RData` / `.csv`; all variables are documented in
[`Data/variable_descriptions.csv`](Data/variable_descriptions.csv).

**Methodology details** (estimator, control matrix, learners, clustering,
repeated cross-fitting) are documented in
[`README_codebase.md`](README_codebase.md).

## Project Structure

```
Empirical_project_ML_Econometrics/
├── Code/
│   ├── 00_setup.R                        # Packages, output folder, DML wrapper (dml_plr),
│   │                                     #   control-matrix builder (build_W), plot helpers
│   ├── 01_data_preparation.R             # Outcome, treatments, control sets; writes prepared_data.rds
│   ├── 02_main_and_joint_effects.R       # Q1 & Q3: average effects + policy-regime parametrization
│   ├── 03a_heterogeneity_primary.R       # Q2: pre-specified treatment-by-characteristic interactions
│   ├── 03b_heterogeneity_all_controls.R  # Q2 robustness: ALL controls interacted (manual, appendix)
│   ├── 04_sensitivity_analysis.R         # Learner choice, control sets, announcement placebo
│   ├── 05a_run_all_sequential.R          # Runner: full pipeline (00 -> 01 -> 02 -> 03a -> 04)
│   └── 05b_run_all_parallel.R            # Same pipeline, 02/03a/04 as parallel background processes
│
├── Data/
│   ├── urban_emissions_panel.RData       # Analysis dataset (also as .csv)
│   ├── urban_emissions_panel.csv
│   └── variable_descriptions.csv         # Codebook for all variables
│
├── output/                               # Created by the runners
│   ├── prepared_data.rds
│   ├── tables/
│   │   ├── tab_main_effects.csv            # from 02
│   │   ├── tab_joint_effects_regimes.csv   # from 02
│   │   ├── tab_heterogeneity_primary.csv   # from 03a
│   │   ├── tab_sensitivity.csv             # from 04
│   │   └── tab_sensitivity_placebo.csv     # from 04
│   └── figures/
│       ├── fig_main_effects.png
│       ├── fig_joint_effects_regimes.png
│       ├── fig_heterogeneity_primary_cp.png
│       ├── fig_heterogeneity_primary_lez.png
│       ├── fig_sensitivity_learner.png
│       └── fig_sensitivity_controls.png
│
├── output_heterogeneity_all_controls/    # Precomputed results of 03b (slowest run)
│   ├── tables/tab_heterogeneity_full.csv
│   └── figures/fig_heterogeneity_full_{cp,lez}.png
│
├── instructions.pdf                      # Course assignment / task description
├── README_codebase.md                    # Script-by-script and methodology documentation
├── README.md
├── .gitignore
└── Empirical_project_ML_Econometrics.Rproj
```

Output filenames mirror the script that produces them, with a `tab_` / `fig_`
prefix distinguishing the table from the figure of the same result (e.g.
`tab_main_effects.csv` and `fig_main_effects.png` both come from
`02_main_and_joint_effects.R`).

## Requirements / Getting Started

- R (≥ 4.2 recommended)
- All required packages are installed and loaded automatically by the setup
  script, which every analysis script sources first:

```r
# --- Installs and loads all required packages ---
source("Code/00_setup.R")
```

Packages used: `DoubleML`, `mlr3`, `mlr3learners`, `data.table`, `glmnet`,
`ranger`, `hdm`, `ggplot2`, `dplyr`, `tidyr` (plus `future` for parallel
cross-fitting). Figures use the Times New Roman font, which is assumed to be
available on the system.

### Running the analysis

Open `Empirical_project_ML_Econometrics.Rproj` (so the working directory is
the project root) and run either:

```r
source("Code/05a_run_all_sequential.R")   # runs everywhere, one script after another
```

or, on a machine with many free cores (~15 for full speed):

```r
source("Code/05b_run_all_parallel.R")     # runs 02/03a/04 side by side, identical results
```

All results are written to `output/tables/` (csv) and `output/figures/` (png).
Both runners produce identical numbers (fixed seeds, `n_rep = 5` repeated
cross-fitting). The full pipeline is compute-intensive — the runners print
elapsed time as they go. To write to a different folder, set
`out_dir <- "my_folder"` before sourcing.

The robustness appendix `Code/03b_heterogeneity_all_controls.R` is **not**
part of the pipeline (it is by far the slowest run). Its precomputed results
ship in `output_heterogeneity_all_controls/`; to reproduce them:

```r
out_dir <- "output"                  # any folder containing prepared_data.rds
source("Code/03b_heterogeneity_all_controls.R")
```
