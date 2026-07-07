# Codebase: Urban emission policies (CP & LEZ) - DML analysis

Run everything with `source("05_run_all.R")` from the project root.
Results are written to `output/tables/` and `output/figures/`.

| Script | Report section | Content |
|---|---|---|
| `00_setup.R` | Proposed estimator | Packages, `build_W`, cross-fitted `dml_plm`, `lincom`, and plotting helpers |
| `01_data_preparation.R` | Data preparation | Outcome, treatments, controls, centered heterogeneity variables, and descriptives |
| `02_main_effects_dml.R` | Results (Q1, Q3) | Average CP and LEZ effects, their interaction, and the total joint effect |
| `03_heterogeneity_dml.R` | Results (Q2) | CATE heterogeneity through treatment-by-characteristic interactions |
| `04_sensitivity_analysis.R` | Sensitivity | Nuisance learners (plugin lasso / CV lasso / ridge / random forest), control sets, and cross-fitting fold counts |
| `05_run_all.R` | Runner | Runs the complete analysis |

## Key modeling choices

- The outcome is log transport CO2, so coefficients are log-point effects.
- Controls exclude likely mediators in the baseline specification.
- `build_W()` adds year and country effects, pairwise interactions, and selected squared terms.
- Nuisance functions use the plugin-lasso (`hdm::rlasso`, theory-based lambda) as baseline; CV-based learners appear as sensitivity comparisons.
- Cross-fitting folds are assigned by city, keeping all years of one city together.
- Inference comes from `summary()`/`confint()` on the final residual-on-residual OLS, as in the lecture; a cluster-robust alternative (`vcov_cluster`) is sketched in `00_setup.R` but not active.
- Heterogeneity follows `Code.R`: policy-by-centered-characteristic interactions approximate the CATE.
