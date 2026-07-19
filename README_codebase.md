# Codebase: Double Machine Learning for Urban Climate Policy Evaluation (CP & LEZ)

Run everything with `source("Code/05a_run_all_sequential.R")` from the
project root, or with `source("Code/05b_run_all_parallel.R")` on a machine
with enough cores (~15) to run scripts 02/03a/04 side by side — both produce
identical results. All scripts live in `Code/` and are written to be sourced
from the project root; the data lives in `Data/`. Results are
written to `output/tables/` (csv) and `output/figures/` (png); reruns
overwrite this folder. To write elsewhere instead, set
`out_dir <- "my_folder"` in the console before sourcing (see `00_setup.R`).

Output filenames mirror the script that writes them, with a `tab_`/`fig_`
prefix for tables vs. figures over the same result (e.g. `tab_main_effects.csv`
and `fig_main_effects.png` both come from `02_main_and_joint_effects.R`).

| Script | Report section | Content |
|---|---|---|
| `00_setup.R` | Proposed estimator | Packages, `build_W`, the `dml_plr` wrapper around the `DoubleML` package, and plotting helpers |
| `01_data_preparation.R` | Data preparation | Outcome, treatments, and control sets |
| `02_main_and_joint_effects.R` | Results (Q1, Q3) | Two parametrizations of the same model: (a) CP, LEZ, and their interaction (Q3 directly); (b) mutually exclusive regimes only CP / only LEZ / both vs. neither, so the total joint effect gets a proper CI. Prints a cross-parametrization consistency check |
| `03a_heterogeneity_primary.R` | Results (Q2) | CATE heterogeneity through treatment-by-characteristic interactions with a small pre-specified set of characteristics; `03b_heterogeneity_all_controls.R` (manual, appendix) repeats this with **all** baseline controls |
| `04_sensitivity_analysis.R` | Sensitivity | Nuisance learners (rlasso / ridge / random forest, all on the same control matrix), control sets, and the announced-but-not-active placebo test |
| `05a_run_all_sequential.R` | Runner | Runs the complete analysis sequentially |
| `05b_run_all_parallel.R` | Runner | Same pipeline, 02/03a/04 as parallel background processes |

## Key modeling choices

- The outcome is log transport CO2, so coefficients are log-point effects.
- Controls exclude likely mediators in the baseline specification.
- `build_W()` adds year and country effects, pairwise interactions, and selected squared terms.
- Estimation runs on the `DoubleML` package (partialling-out score, DML2). With several treatment columns, DoubleML's default (`use_other_treat_as_covariate = TRUE`) estimates each coefficient in turn, folding the other treatments into that run's covariates - exactly the lecture's "One-By-One Double LASSO" procedure for a vector of target coefficients. Nuisance functions use the plugin lasso (`hdm::rlasso`) as baseline; ridge and a random forest appear as learner-choice sensitivity comparisons in 04, chosen because each tests a specific assumption (sparsity; linear-after-selection form) rather than just a different lambda-selection rule.
- Cross-fitting folds are assigned by city (via `cluster_cols`), keeping all years of one city together; standard errors are cluster-robust at the city level.
- The main estimates use `n_rep = 5` repeated cross-fitting splits, so fold-assignment noise is included in the reported confidence intervals.
- Heterogeneity follows the lecture's wage-gap example: policy-by-centered-characteristic interactions approximate the CATE. The primary spec (03a) uses a small set of characteristics chosen ex ante; the manual companion `03b_heterogeneity_all_controls.R` interacts all baseline controls (including the noise variables, as a selection check) as a robustness appendix.
