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
#          are assigned at the city level, with standard errors clustered
#          at the city level too (cluster_cols below). With several
#          treatment columns, DoubleML's default (use_other_treat_as_covariate
#          = TRUE, left unchanged here) estimates each coefficient in turn,
#          folding the other treatments into that run's covariates - this is
#          exactly "One-By-One Double LASSO" from the lecture (see dml_plr
#          below for the precise correspondence).
# =============================================================================

packages <- c("DoubleML", "mlr3", "mlr3learners", "data.table",
              "glmnet", "ranger", "hdm", "ggplot2", "dplyr", "tidyr")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))
lgr::get_logger("mlr3")$set_threshold("warn")   # silence per-fold fitting logs

# Parallel cross-fitting: DoubleML trains the per-fold nuisance models via
# mlr3::resample(), which distributes folds across the workers declared here
# (future is an mlr3 dependency). n_folds = 5 folds run at once, so more than
# 5 workers buys nothing except for the 10-fold sensitivity spec.
future::plan("multisession", workers = 5)

# All results are written below out_dir. Date-stamped so a rerun never
# overwrites the results of an earlier run (the original results live in
# "output/"). The exists() guard keeps one consistent folder when several
# scripts source() this file in the same session, and also lets you pick a
# name manually before running, e.g.:  out_dir <- "output_nrep5"
if (!exists("out_dir")) out_dir <- paste0("output_", format(Sys.Date()))
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
cat("Results are written to:", out_dir, "\n")

# Save a results table as csv and print it to the console
save_table <- function(tab, name) {
  write.csv(tab, file.path(out_dir, "tables", paste0(name, ".csv")),
            row.names = FALSE)
  cat("\n====", name, "====\n")
  print(tab, digits = 4)
}

# Convert a log-point estimate into an exact percentage effect:
# 100 * (exp(beta) - 1), as in the lecture's wage-gap interpretation.
pct <- function(x) 100 * (exp(x) - 1)

# -----------------------------------------------------------------------------
# build_W: construct the high-dimensional control matrix
#
# Same construction as the lecture's control matrix:
#   - model.matrix() with (controls)^2 -> main effects + pairwise
#     interactions, flexible controls without touching the final model
#   - drop constant columns (no information, breaks lasso numerically)
#   - demean every column
# Added for this panel:
#   - year dummies (common shocks: fuel prices, EU regulation, COVID...)
#   - country dummies (time-constant country differences)
#   - squared terms for a few continuous controls via I(x^2)
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
# LearnerRegrRlasso: custom mlr3 learner wrapping hdm::rlasso
#
# Lasso with the theory-based plugin penalty of Belloni/Chernozhukov/Hansen,
# the same lambda rule as rlassoEffects() in the lecture. lambda comes from
# a closed-form formula, so each nuisance fit is a single lasso fit instead
# of the ~11 path fits of cv.glmnet's inner 10-fold CV (~10x faster).
# Cross-fitting is kept anyway even though the plugin theory would allow
# full-sample fitting, so inference doesn't rest on the exact sparsity
# conditions full-sample validity would require.
# No mlr3 package wraps hdm, hence this minimal learner (train/predict only).
# -----------------------------------------------------------------------------
LearnerRegrRlasso <- R6::R6Class("LearnerRegrRlasso",
  inherit = mlr3::LearnerRegr,
  public = list(
    initialize = function() {
      super$initialize(
        id            = "regr.rlasso",
        feature_types = c("numeric", "integer"),
        predict_types = "response",
        param_set     = paradox::ps(
          # post = FALSE: use the lasso coefficients directly for prediction,
          # no post-lasso OLS refit (matches the lecture's baseline rlasso call)
          post = paradox::p_lgl(default = FALSE, tags = "train")),
        packages      = "hdm",
        label         = "Lasso with plugin lambda (hdm::rlasso)"
      )
    }
  ),
  private = list(
    .train = function(task) {
      pv <- self$param_set$get_values(tags = "train")
      hdm::rlasso(x = as.matrix(task$data(cols = task$feature_names)),
                  y = task$data(cols = task$target_names)[[1]],
                  post = isTRUE(pv$post))
    },
    .predict = function(task) {
      X <- as.matrix(task$data(cols = task$feature_names))
      list(response = as.numeric(predict(self$model, newdata = X)))
    }
  )
)

# -----------------------------------------------------------------------------
# make_learner: mlr3 nuisance learners for E[Y|W] and E[D|W]
#
#   "rlasso" - hdm::rlasso with the plugin lambda (see LearnerRegrRlasso
#              above). BASELINE, used in 02/03: single fit per nuisance,
#              no inner CV; same penalty rule as the lecture's
#              rlassoEffects().
#   "ridge"  - mlr3learners' built-in regr.cv_glmnet (alpha = 0, lambda by
#              10-fold CV, one-standard-error rule) - no custom wrapper
#              needed, unlike rlasso, since glmnet already ships an mlr3
#              learner. Not sparse (see 04): checks whether the result
#              depends on rlasso's sparsity assumption specifically.
#   "rf"     - random forest (ranger); pass a W without interactions,
#              since the forest captures nonlinearities itself. Checks
#              whether the result depends on the linear-after-selection
#              functional form.
# Only rlasso is used in 02/03; ridge and rf appear as learner-choice
# sensitivity rows in 04 (CV-tuned lasso variants were deliberately left
# out there - see 04's comment on why).
# -----------------------------------------------------------------------------
make_learner <- function(learner = c("rlasso", "ridge", "rf")) {
  learner <- match.arg(learner)
  switch(learner,
         rlasso = LearnerRegrRlasso$new(),
         ridge  = lrn("regr.cv_glmnet", alpha = 0, s = "lambda.1se"),
         # num.threads = 1: folds already run in parallel (future workers,
         # see plan() above); ranger's own threading on top would oversubscribe
         rf     = lrn("regr.ranger", num.trees = 500, min.node.size = 5,
                      num.threads = 1))
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
#      the city level (all years of a city stay in one fold, so there's no
#      leakage through serial correlation) and SEs are cluster-robust at
#      the city level.
#   2. For each fold: learn E[Y|W] and E[D|W] on the other folds, predict
#      on the held-out fold, form out-of-fold residuals.
#   3. Solve the orthogonal moment condition on the pooled residuals (DML2).
#   4. With n_rep > 1, repeat everything on n_rep independent fold splits;
#      the reported estimate is the median across splits, and the reported
#      variance includes the dispersion across splits,
#      median(se_r^2 + (theta_r - theta_median)^2), so split noise is
#      already part of the confidence interval (no separate seed analysis
#      needed).
#
# D can contain several treatment columns (cp_active, lez_active, their
# interaction, treatment x characteristic interactions...). With
# use_other_treat_as_covariate = TRUE (DoubleMLData's own default, left
# unchanged), DoubleML estimates each coefficient IN TURN, moving the other
# treatment columns into that run's covariates alongside W. This is exactly
# "One-By-One Double LASSO" from the lecture (Lecture 10, "Inference on
# Multiple Coefficients"):
#   for each l = 1,...,q: run Double LASSO of Y on D_l, treating
#   D_(-l)'beta_{1,-l} + W'beta_2 - the OTHER treatments together with the
#   genuine controls, as ONE nuisance function - then collect the q
#   estimates with marginal confidence intervals. hdm::rlassoEffects()
#   (the lecture's function for this, confirmed by its own source: it
#   loops over each target column and folds every other column, including
#   the other targets, into one combined control matrix) implements the
#   identical procedure.
# -----------------------------------------------------------------------------
dml_plr <- function(Y, D, W, cluster, learner = "rlasso",
                    n_folds = 5, n_rep = 5, seed = 42, level = 0.95) {
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
# Shared plot styling: large-font theme (report/print size) plus a color rule
# for coefficient plots - grey if the 95% CI includes 0 (not significant),
# red if significant and negative (lower emissions), steelblue if significant
# and positive (higher emissions). The facet_wrap variant of this theme
# (theme.adjusted.facet) is only needed by the sensitivity figures, so
# it's defined locally in 04_sensitivity_analysis.R instead of here.
# -----------------------------------------------------------------------------
theme.adjusted <- theme(
  text = element_text(family = "Times New Roman"),
  axis.text.x = element_text(angle = 0, hjust = 0.5, margin = margin(t = 5), size = 22),
  axis.title.x = element_text(margin = margin(t = 20), size = 32),
  axis.text.y = element_text(hjust = 1, margin = margin(r = 10), size = 22, angle = 0),
  axis.title.y = element_text(margin = margin(r = 20), size = 32),
  title = element_text(color = "black"),
  plot.title = element_text(size = 28, color = "black", face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 17, color = "black", face = "italic"),
  panel.grid.major = element_line(color = "darkgray", linewidth = 0.2),
  panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
  plot.background = element_rect(fill = "white", color = NA)
)

classify_colour <- function(estimate, conf.low, conf.high) {
  ifelse(conf.low <= 0 & conf.high >= 0, "grey50",
         ifelse(estimate < 0, "red", "steelblue"))
}

# -----------------------------------------------------------------------------
# pretty_term: human-readable axis/facet labels for the raw D column names
# (only_cp, cp_x_lez, cp_x_log_population, ...). CP = congestion pricing,
# LEZ = low-emission zone - spelled out once in each plot's title, kept as
# the abbreviation everywhere else (axis text, facet strips).
# -----------------------------------------------------------------------------
pretty_term <- function(term) {
  base <- c(
    cp_active  = "CP",
    lez_active = "LEZ",
    cp_x_lez   = "CP:LEZ",
    only_cp    = "CP only",
    only_lez   = "LEZ only",
    both       = "CP & LEZ (both)"
  )
  het <- c(
    log_population        = "population (log)",
    log_gdp_pc            = "GDP per capita (log)",
    log_pop_density       = "population density (log)",
    public_transit_score  = "public transit score",
    political_green       = "green party vote share",
    fuel_price             = "fuel price",
    tourism_intensity      = "tourism intensity",
    industry_logistics     = "logistics industry share",
    coastal                = "coastal city"
  )
  out <- unname(base[term])
  cp_hit  <- is.na(out) & grepl("^cp_x_",  term)
  lez_hit <- is.na(out) & grepl("^lez_x_", term)
  out[cp_hit]  <- paste0("CP:",  het[sub("^cp_x_",  "", term[cp_hit])])
  out[lez_hit] <- paste0("LEZ:", het[sub("^lez_x_", "", term[lez_hit])])
  ifelse(is.na(out), term, out)
}

# -----------------------------------------------------------------------------
# plot_effects: coefficient plot with confidence intervals
# -----------------------------------------------------------------------------
plot_effects <- function(res, title, file = NULL, subtitle = NULL,
                         x_breaks = scales::pretty_breaks(n = 10)) {
  res$sig_colour  <- classify_colour(res$estimate, res$conf.low, res$conf.high)
  res$term_label  <- pretty_term(res$term)
  p <- ggplot(res, aes(x = estimate, y = reorder(term_label, estimate))) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7, colour = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = sig_colour),
                   height = 0.25, linewidth = 1) +
    geom_point(aes(colour = sig_colour), size = 3.5) +
    scale_colour_identity() +
    # check.overlap = FALSE: always show every category's label, even if
    # adjacent labels visually touch - otherwise ggplot silently drops
    # every other one when the plot is short on vertical space, which is
    # why some plots showed all rows and others only every second one.
    scale_y_discrete(guide = guide_axis(check.overlap = FALSE)) +
    # x_breaks default targets ~10 "nice" breaks; pass
    # scales::breaks_width(1) or breaks_width(0.5) per call for plots
    # whose narrow range makes that target-based default too dense/uneven.
    scale_x_continuous(breaks = x_breaks) +
    labs(x = "Effect on log transport CO2", y = NULL,
         title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme.adjusted
  if (!is.null(file)) {
    # Fixed size for every figure in the report, so they all sit at the
    # same \includegraphics width/aspect ratio without per-figure tuning.
    ggsave(file.path(out_dir, "figures", file), p,
           width = 10, height = 4, dpi = 300)
  }
  p
}
