# =============================================================================
# 00_setup.R
# Packages, output folders, and helper functions used by all other scripts.
#
# Project: Effects of congestion pricing (CP) and low-emission zones (LEZ)
#          on urban road-transport CO2 emissions.
# Method:  Double/Debiased Machine Learning (DML) for the partially linear
#          model, estimated with the DoubleML package (Bach/Chernozhukov/
#          Kurz/Spindler) - the reference implementation of Chernozhukov
#          et al. (2018) - with mlr3 nuisance learners. Cross-fitting folds
#          are assigned at the CITY level and standard errors are
#          cluster-robust at the city level (cluster_cols below).
# =============================================================================

packages <- c("DoubleML", "mlr3", "mlr3learners", "data.table",
              "glmnet", "ranger", "ggplot2", "dplyr", "tidyr")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))
lgr::get_logger("mlr3")$set_threshold("warn")   # silence per-fold fitting logs

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

# Save a results table as csv and print it to the console
save_table <- function(tab, name) {
  write.csv(tab, file.path("output/tables", paste0(name, ".csv")),
            row.names = FALSE)
  cat("\n====", name, "====\n")
  print(tab, digits = 4)
}

# Convert a log-point estimate into an exact percentage effect:
# 100 * (exp(beta) - 1), as in the wage-gap interpretation in Code.R
pct <- function(x) 100 * (exp(x) - 1)

# -----------------------------------------------------------------------------
# build_W: construct the high-dimensional control matrix
#
# Mirrors the construction in Code.R:
#   - model.matrix() with (controls)^2 -> all main effects + pairwise
#     interactions (flexible controls without estimating the final model yet)
#   - remove constant columns (no information, numerical problems for lasso)
#   - demean every column
# Additions for this panel application:
#   - year dummies (common shocks: fuel price trends, EU regulation, COVID...)
#   - country dummies (time-constant country differences)
#   - squared terms for selected continuous controls via I(x^2)
# -----------------------------------------------------------------------------
build_W <- function(df, ctrl_vars, sq_vars = NULL, interactions = TRUE) {
  rhs <- paste(ctrl_vars, collapse = " + ")
  if (interactions) rhs <- paste0("(", rhs, ")^2")
  if (!is.null(sq_vars)) {
    rhs <- paste(rhs, "+", paste(paste0("I(", sq_vars, "^2)"), collapse = " + "))
  }
  f <- as.formula(paste("~ -1 + factor(year) + country_id +", rhs))
  W <- model.matrix(f, data = df)
  W <- W[, apply(W, 2, var) > 0, drop = FALSE]        # drop constant columns
  # syntactic, unique column names (":" -> "."), consistent with D in dml_plr
  colnames(W) <- make.names(colnames(W), unique = TRUE)
  scale(W, center = TRUE, scale = FALSE)               # demean each column
}

# -----------------------------------------------------------------------------
# make_learner: mlr3 nuisance learners for E[Y|W] and E[D|W]
#
#   "lasso"     - cv.glmnet, alpha = 1, lambda by 10-fold CV with the
#                 one-standard-error rule (lambda.1se): the most strongly
#                 regularized model within one SE of the CV-error minimum,
#                 as recommended in the lecture. BASELINE learner.
#   "lasso_min" - as above but lambda at the CV minimum (sensitivity check)
#   "ridge"     - cv.glmnet, alpha = 0, lambda.1se
#   "rf"        - random forest (ranger); pass a W without interactions,
#                 since the forest captures nonlinearities itself
# -----------------------------------------------------------------------------
make_learner <- function(learner = c("lasso", "lasso_min", "ridge", "rf")) {
  learner <- match.arg(learner)
  switch(learner,
         lasso     = lrn("regr.cv_glmnet", alpha = 1, s = "lambda.1se"),
         lasso_min = lrn("regr.cv_glmnet", alpha = 1, s = "lambda.min"),
         ridge     = lrn("regr.cv_glmnet", alpha = 0, s = "lambda.1se"),
         rf        = lrn("regr.ranger", num.trees = 500, min.node.size = 5))
}

# -----------------------------------------------------------------------------
# dml_plr: cross-fitted DML for the partially linear model (via DoubleML)
#
#   Y = D %*% theta + g(W) + e,   E[e | D, W] = 0
#
# Thin wrapper around DoubleMLPLR that only wires the data together and
# collects the results in a data frame; all estimation, cross-fitting, and
# inference is done by the package. What the package does internally
# (Neyman-orthogonal "partialling out" score, DML2 procedure):
#   1. Split observations into n_folds folds. Because the data enter as
#      DoubleMLClusterData with cluster_cols = city, folds are assigned at
#      the CITY level (all years of a city stay in one fold - no leakage
#      through serial correlation) and SEs are CLUSTER-ROBUST at the city
#      level.
#   2. For each fold: learn E[Y|W] and E[D|W] on the other folds, predict
#      on the held-out fold, form out-of-fold residuals.
#   3. Solve the orthogonal moment condition on the pooled residuals (DML2).
#   4. With n_rep > 1, repeat everything on n_rep independent fold splits;
#      the reported estimate is the MEDIAN across splits and the reported
#      variance includes the dispersion across splits,
#      median(se_r^2 + (theta_r - theta_median)^2), so split noise is part
#      of the confidence interval (no separate seed analysis needed).
#
# D can contain several treatment columns (cp_active, lez_active, their
# interaction, treatment x characteristic interactions...). DoubleML then
# estimates each coefficient IN TURN, moving the other treatment columns
# into the controls for that run ("use other treat as covariate").
# -----------------------------------------------------------------------------
dml_plr <- function(Y, D, W, cluster, learner = "lasso",
                    n_folds = 5, n_rep = 1, seed = 42, level = 0.95) {
  D <- as.matrix(D)
  # syntactic, unique names: consistent with build_W, required by mlr3
  colnames(D) <- make.names(colnames(D), unique = TRUE)
  stopifnot(nrow(D) == length(Y), nrow(W) == length(Y),
            length(cluster) == length(Y))

  df <- data.table::as.data.table(cbind(
    data.frame(Y_out = Y, city_cl = as.character(cluster)), D, W))
  dml_data <- DoubleMLClusterData$new(
    df,
    y_col        = "Y_out",
    d_cols       = colnames(D),
    x_cols       = colnames(W),
    cluster_cols = "city_cl")

  ml <- make_learner(learner)
  set.seed(seed)                       # reproducible fold assignment
  model <- DoubleMLPLR$new(dml_data,
                           ml_l = ml$clone(),   # nuisance E[Y|W]
                           ml_m = ml$clone(),   # nuisance E[D|W]
                           n_folds = n_folds,
                           n_rep   = n_rep,
                           score   = "partialling out")
  model$fit()

  ci <- model$confint(level = level)
  list(results = data.frame(term      = names(model$coef),
                            estimate  = model$coef,
                            std.error = model$se,
                            p.value   = model$pval,
                            conf.low  = ci[, 1],
                            conf.high = ci[, 2],
                            row.names = NULL),
       model = model)
}

# -----------------------------------------------------------------------------
# plot_effects: coefficient plot with confidence intervals
# -----------------------------------------------------------------------------
plot_effects <- function(res, title, file = NULL) {
  p <- ggplot(res, aes(x = estimate, y = reorder(term, estimate))) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    labs(x = "Effect on log transport CO2 (95% CI)", y = NULL,
         title = title) +
    theme_minimal(base_size = 12)
  if (!is.null(file)) {
    ggsave(file.path("output/figures", file), p,
           width = 8, height = max(3, 0.4 * nrow(res) + 1.5))
  }
  p
}
