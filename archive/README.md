# Archive: estimation-method evolution (why the code looks the way it does)

This folder preserves results from the two retired estimation setups, so the
current baseline can be compared against them and the design decisions can be
reconstructed later (e.g. for the report's methods section or a defense).

The three generations, in order:

| generation | folder | code | learner | inference |
|---|---|---|---|---|
| 1. hand-rolled DML | `2026-07-07_hand_rolled_pipeline/` | `git show 72ba19b:00_setup.R` (fn `dml_plm`) | cv.glmnet, shared Y-fit, joint OLS | classical lm SEs (NOT clustered) |
| 2. DoubleML + CV-lasso | `2026-07-08_doubleml_cv_lasso/` | commit ff0f478 | cv.glmnet lambda.1se, per-column | cluster-robust (package) |
| 3. DoubleML + plugin lasso (CURRENT) | `2026-07-08_doubleml_rlasso_baseline/` | working tree 2026-07-08 | hdm::rlasso plugin lambda via custom mlr3 learner | cluster-robust (package) |

## Headline comparison: total effect of both policies vs no policy

| generation | estimate | std.error | runtime (script 02) |
|---|---|---|---|
| 1. hand-rolled (joint OLS + lincom)  | -0.785 | **0.041** | ~10 min |
| 2. DoubleML CV-lasso (n_rep=1)       | -0.637 | 0.129 | ~12 min (parallel) / 80 min (n_rep=5, serial) |
| 3. DoubleML rlasso (n_rep=1)         | -0.776 | 0.150 | **~3 min** |

## Why generation 1 was retired

- **Its default SEs are invalid for this panel.** The joint OLS used classical
  `summary(lm)` SEs, which treat 3,100 city-years as independent; with ~183
  cities and strong within-city correlation they are overconfident by a factor
  of ~4 (SE 0.041 vs ~0.15 clustered; p = 1.5e-82 is not a credible p-value).
  A hand-rolled `vcov_cluster` existed but was optional/commented out.
  Note the POINT estimate (-0.785) agrees almost perfectly with generation 3
  (-0.776) - the old code was not wrong about the effect, only about the
  uncertainty.
- Per-coefficient orthogonality: DoubleML moves the other treatment columns
  into the nuisances (flexible, ML-based) instead of handling treatment
  correlation linearly in the final OLS.
- Package = reference implementation of Chernozhukov et al. (2018); no
  hand-rolled estimation code to defend line by line.

**When to revisit generation 1:** if a CI for a linear combination of
coefficients is ever needed (the joint OLS provides coefficient covariances,
DoubleML does not - the regime parametrization in 02 exists as the
workaround), or as a didactic reference for what DoubleML does internally.

## Why generation 2 was retired as baseline (kept as sensitivity row)

- Runtime: the inner 10-fold CV of cv.glmnet makes every nuisance fit ~11
  lasso path fits; with per-column refits and n_rep=5 this pushed script 02
  to ~80 minutes.
- The plugin lambda (Belloni/Chernozhukov/Hansen, as in the lecture's
  rlassoEffects) needs one fit per nuisance, ~10x faster, and ties the
  project back to the course material.
- CV-lasso (lambda.1se and lambda.min) remains as learner-sensitivity rows in
  04_sensitivity_analysis.R - the plugin-vs-CV comparison doubles as a check
  on the custom rlasso learner wrapper and on the sparsity assumption behind
  the plugin rule.

## Interesting observations from the CV-vs-plugin comparison (02, n_rep=1, seed 42)

Same seed, same city-level folds, so differences are purely the lambda rule:

| term | CV-lasso | rlasso | note |
|---|---|---|---|
| cp_active | -0.169 | -0.180 | stable; SE shrinks 0.138 -> 0.094 (stronger regularization leaves more identifying variation in D) |
| lez_active | +0.460 | +0.289 | drops ~1/3: the LEZ coefficient is sensitive to how flexibly controls are partialled out - consistent with the suspicion that it partly reflects selection/confounding (see placebo check in 04), not a causal "+X%" |
| cp_x_lez | -1.111 | -0.944 | stable in sign/size; synergy story unchanged |
| both (run 2) | -0.637 | -0.776 | cross-parametrization consistency check theta_cp+theta_lez+theta_int vs both: gap 0.18 (CV) vs 0.06 (rlasso) - internal consistency improved |

All signs, orderings, and the significance pattern are identical across the
three generations: CP alone small/uncertain, LEZ alone positive (suspicious,
see above), strong negative synergy, total effect of both around -50%.

## File inventory

- `2026-07-07_hand_rolled_pipeline/`: original output of generation 1
  (both-policies lincom table; heterogeneity and sensitivity tables from the
  2026-07-07 20:24-21:22 run). The code: `git show 72ba19b` (or
  `git show ff0f478^:00_setup.R`).
- `2026-07-08_doubleml_cv_lasso/`: script-02 tables of generation 2
  (n_rep=1 dev run; values transcribed from the console at 4 decimals).
- `2026-07-08_doubleml_rlasso_baseline/`: script-02 tables of the first
  generation-3 run (n_rep=1), for comparison against later n_rep=5 output.
