# Codebase: Urban emission policies (CP & LEZ) - DML analysis

Run everything with `source("05_run_all.R")` from the project root.
Results are written to a date-stamped folder, `output_YYYY-MM-DD/tables/` and
`output_YYYY-MM-DD/figures/`, so reruns never overwrite earlier results (the
original run lives in `output/`). To pick a folder name yourself, set
`out_dir <- "my_folder"` in the console before sourcing (see `00_setup.R`).

| Script | Report section | Content |
|---|---|---|
| `00_setup.R` | Proposed estimator | Packages, `build_W`, the `dml_plr` wrapper around the `DoubleML` package, and plotting helpers |
| `01_data_preparation.R` | Data preparation | Outcome, treatments, controls, centered heterogeneity variables, and descriptives |
| `02_main_effects_dml.R` | Results (Q1, Q3) | Average CP and LEZ effects and their interaction; a second run with mutually exclusive regime dummies gives the total joint effect directly |
| `03_heterogeneity_dml.R` | Results (Q2) | CATE heterogeneity through treatment-by-characteristic interactions |
| `04_sensitivity_analysis.R` | Sensitivity | Nuisance learners (rlasso / ridge / random forest), control sets, excluding the COVID years, and the announced-but-not-active placebo test |
| `05_run_all.R` | Runner | Runs the complete analysis |

## Key modeling choices

- The outcome is log transport CO2, so coefficients are log-point effects.
- Controls exclude likely mediators in the baseline specification.
- `build_W()` adds year and country effects, pairwise interactions, and selected squared terms.
- Estimation runs on the `DoubleML` package (partialling-out score, DML2). With several treatment columns, DoubleML's default (`use_other_treat_as_covariate = TRUE`) estimates each coefficient in turn, folding the other treatments into that run's covariates - exactly the lecture's "One-By-One Double LASSO" procedure for a vector of target coefficients. Nuisance functions use the plugin lasso (`hdm::rlasso`) as baseline; ridge and a random forest appear as learner-choice sensitivity comparisons in 04, chosen because each tests a specific assumption (sparsity; linear-after-selection form) rather than just a different lambda-selection rule.
- Cross-fitting folds are assigned by city (via `cluster_cols`), keeping all years of one city together; standard errors are cluster-robust at the city level.
- The main estimates use `n_rep = 5` repeated cross-fitting splits, so fold-assignment noise is included in the reported confidence intervals.
- Heterogeneity follows the lecture's wage-gap example: policy-by-centered-characteristic interactions approximate the CATE.
