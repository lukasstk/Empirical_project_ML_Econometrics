# Codebase: Urban emission policies (CP & LEZ) - DML analysis

Run everything with `source("05_run_all.R")` from the project root.
Results are written to `output/tables/` and `output/figures/`.

| Script | Report section | Content |
|---|---|---|
| `00_setup.R` | Proposed estimator | Packages, `build_W`, the `dml_plr` wrapper around the `DoubleML` package, and plotting helpers |
| `01_data_preparation.R` | Data preparation | Outcome, treatments, controls, centered heterogeneity variables, and descriptives |
| `02_main_effects_dml.R` | Results (Q1, Q3) | Average CP and LEZ effects and their interaction; a second run with mutually exclusive regime dummies gives the total joint effect directly |
| `03_heterogeneity_dml.R` | Results (Q2) | CATE heterogeneity through treatment-by-characteristic interactions |
| `04_sensitivity_analysis.R` | Sensitivity | Nuisance learners (CV lasso 1se/min / ridge / random forest), control sets, fold counts, and the announced-but-not-active placebo test |
| `05_run_all.R` | Runner | Runs the complete analysis |

## Key modeling choices

- The outcome is log transport CO2, so coefficients are log-point effects.
- Controls exclude likely mediators in the baseline specification.
- `build_W()` adds year and country effects, pairwise interactions, and selected squared terms.
- Estimation runs entirely on the `DoubleML` package (partialling-out score, DML2); nuisance functions use the CV-lasso with the lambda.1se rule as baseline, other learners appear as sensitivity comparisons.
- Cross-fitting folds are assigned by city (via `cluster_cols`), keeping all years of one city together; standard errors are cluster-robust at the city level.
- The main estimates use `n_rep = 5` repeated cross-fitting splits, so fold-assignment noise is included in the reported confidence intervals.
- Heterogeneity follows `Code.R`: policy-by-centered-characteristic interactions approximate the CATE.
