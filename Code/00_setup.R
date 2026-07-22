# =============================================================================
# 00_setup.R
# Packages, output folders, and helper functions used by all other scripts.
#
# Project: Effects of congestion pricing (CP) and low-emission zones (LEZ)
#          on urban road-transport CO2 emissions.
# Method:  Double/Debiased Machine Learning (DML) for the partially linear
#          model with city fixed effects (DoubleML package, mlr3 learners).
# =============================================================================

packages <- c("DoubleML", "mlr3", "mlr3learners", "data.table",
              "glmnet", "ranger", "hdm", "ggplot2", "dplyr", "tidyr")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))
lgr::get_logger("mlr3")$set_threshold("warn")   # silence per-fold fitting logs

# Parallel cross-fitting: nuisance fits are distributed across 5 workers
# (one per fold).
future::plan("multisession", workers = 5)

# Output folder; set out_dir before sourcing to pick a different name.
if (!exists("out_dir")) out_dir <- "output"
ensure_out_dirs <- function(dir) {
  dir.create(file.path(dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
}
ensure_out_dirs(out_dir)
cat("Results are written to:", out_dir, "\n")

# Save a results table as csv and print it to the console
save_table <- function(tab, name) {
  write.csv(tab, file.path(out_dir, "tables", paste0(name, ".csv")),
            row.names = FALSE)
  cat("\n====", name, "====\n")
  print(tab, digits = 4)
}

# Log-point estimate -> exact percentage effect: 100 * (exp(beta) - 1)
pct <- function(x) 100 * (exp(x) - 1)

# ---- demean_by_city_fwl: within-transformation (city fixed effects) -------------
# Subtract each city's own mean over its years from every column.
# rowsum(M, id) and table(id) both sort groups alphabetically by default,
# so they'd already line up - but n_i is re-indexed by rownames(sums)
# explicitly, so correctness never depends on that coincidence.
# sums / n_i then recycles row-wise (n_i[i] divides row i of sums).
demean_by_city_fwl <- function(M, id) {
  M    <- as.matrix(M)
  id   <- as.character(id)
  sums <- rowsum(M, id)                        # one row per city (sorted)
  n_i  <- as.vector(table(id)[rownames(sums)]) # years per city, same order
  means <- sums / n_i
  M - means[id, , drop = FALSE]
}

# ---- build_W: construct the high-dimensional control matrix -----------------
# Main effects + all pairwise interactions + squared terms for sq_vars +
# year and country dummies; constant columns dropped; every column
# within-city demeaned (city FE), then zero-variance columns dropped again.
build_W <- function(df, ctrl_vars, sq_vars = NULL, interactions = TRUE) {
  rhs <- paste(ctrl_vars, collapse = " + ")
  if (interactions) rhs <- paste0("(", rhs, ")^2")
  if (!is.null(sq_vars)) {
    rhs <- paste(rhs, "+", paste(paste0("I(", sq_vars, "^2)"), collapse = " + "))
  }
  f <- as.formula(paste("~ -1 + factor(year) + country_id +", rhs))
  W <- model.matrix(f, data = df)
  W <- W[, apply(W, 2, var) > 0, drop = FALSE]        # drop constant columns
  # syntactic, unique column names (":" -> ".") so mlr3 accepts them
  colnames(W) <- make.names(colnames(W), unique = TRUE)
  W <- demean_by_city_fwl(W, df$city_id)                  # city fixed effects
  W[, apply(W, 2, var) > 1e-12, drop = FALSE]
}

# ---- LearnerRegrRlasso: mlr3 learner wrapping hdm::rlasso --------------------
# Lasso with the theory-based plugin penalty of Belloni/Chernozhukov/Hansen.
# DoubleMLPLR requires ml_l/ml_m to be mlr3::Learner objects - hdm::rlasso is
# a standalone function, not one. This class just adapts it: .train()/.predict()
# delegate to hdm::rlasso()/predict.rlasso() under the hood, so DoubleML sees
# an ordinary learner while the actual estimation is unchanged.
LearnerRegrRlasso <- R6::R6Class("LearnerRegrRlasso",
  inherit = mlr3::LearnerRegr,
  public = list(
    initialize = function() {
      super$initialize(
        id            = "regr.rlasso",
        feature_types = c("numeric", "integer"),
        predict_types = "response",
        param_set     = paradox::ps(
          # post = FALSE: no post-lasso OLS refit
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

# ---- make_learner: mlr3 nuisance learners for E[Y|W] and E[D|W] --------------
#   "rlasso" - plugin-lambda lasso (baseline, used in 02/03)
#   "ridge"  - cv_glmnet with alpha = 0, lambda.1se (sensitivity, 04)
#   "rf"     - random forest via ranger (sensitivity, 04)
make_learner <- function(learner = c("rlasso", "ridge", "rf")) {
  learner <- match.arg(learner)
  switch(learner,
         rlasso = LearnerRegrRlasso$new(),
         ridge  = lrn("regr.cv_glmnet", alpha = 0, s = "lambda.1se"),
         # num.threads = 1: folds already run in parallel (future workers)
         rf     = lrn("regr.ranger", num.trees = 500, min.node.size = 5,
                      num.threads = 1))
}

# ---- dml_plr: cross-fitted DML for the partially linear model ----------------
#   Y = D %*% theta + g(W) + e,   E[e | D, W] = 0
# Y and D get the same within-city demeaning as W in build_W. Folds are
# assigned at the city level and SEs are cluster-robust at the city level
# (DoubleMLClusterData). With several treatment columns, each coefficient is
# estimated in turn with the other treatments folded into that run's
# covariates (DoubleML default). n_rep > 1 repeats the cross-fitting on
# independent fold splits and folds the across-split dispersion into the SEs.
dml_plr <- function(Y, D, W, cluster, learner = "rlasso",
                    n_folds = 5, n_rep = 5, seed = 42, level = 0.95) {
  D <- as.matrix(D)
  # syntactic, unique names: consistent with build_W, required by mlr3
  colnames(D) <- make.names(colnames(D), unique = TRUE)
  stopifnot(nrow(D) == length(Y), nrow(W) == length(Y),
            length(cluster) == length(Y))

  # City fixed effects: within-demean Y and D per city
  Y <- as.vector(demean_by_city_fwl(matrix(Y, ncol = 1), cluster))
  D <- demean_by_city_fwl(D, cluster)

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

# ---- Shared plot styling ------------------------------------------------------
theme.adjusted <- theme(
  text = element_text(family = "Times New Roman"),
  axis.text.x = element_text(angle = 0, hjust = 0.5, margin = margin(t = 5), size = 15),
  axis.title.x = element_text(margin = margin(t = 20), size = 16),
  axis.text.y = element_text(hjust = 1, margin = margin(r = 10), size = 15, angle = 0),
  axis.title.y = element_text(margin = margin(r = 20), size = 20),
  title = element_text(color = "black"),
  panel.grid.major = element_line(color = "darkgray", linewidth = 0.2),
  panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
  plot.background = element_rect(fill = "white", color = NA)
)

# CI colour rule, mapped to the meaning of the effect: grey = CI includes 0,
# steelblue = significant emission reduction (negative coefficient),
# red = significant emission increase (positive coefficient)
classify_colour <- function(estimate, conf.low, conf.high) {
  ifelse(conf.low <= 0 & conf.high >= 0, "grey50",
         ifelse(estimate < 0, "steelblue", "red"))
}

# Significance stars, R model convention: *** p<0.001, ** p<0.01, * p<0.05
sig_stars <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

# ---- pretty_term: human-readable labels for the raw D column names ----------
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
    log_tourism_intensity  = "tourism intensity (log)",
    log_elevation          = "elevation (log)",
    industry_logistics     = "logistics industry share",
    coastal                = "coastal city"
  )
  out <- unname(base[term])
  cp_hit  <- is.na(out) & grepl("^cp_x_",  term)
  lez_hit <- is.na(out) & grepl("^lez_x_", term)
  # unmapped characteristics fall back to their raw variable name
  cp_z  <- sub("^cp_x_",  "", term[cp_hit])
  lez_z <- sub("^lez_x_", "", term[lez_hit])
  out[cp_hit]  <- paste0("CP:",  ifelse(is.na(het[cp_z]),  cp_z,  het[cp_z]))
  out[lez_hit] <- paste0("LEZ:", ifelse(is.na(het[lez_z]), lez_z, het[lez_z]))
  ifelse(is.na(out), term, out)
}

# ---- plot_effects: coefficient plot with confidence intervals ---------------
# No plot titles: the report supplies figure captions.
plot_effects <- function(res, file = NULL,
                         x_breaks = scales::pretty_breaks(n = 10),
                         height = 4) {
  res$sig_colour  <- classify_colour(res$estimate, res$conf.low, res$conf.high)
  res$stars       <- sig_stars(res$p.value)
  res$term_label  <- pretty_term(res$term)
  # row order follows the input table (reversed: first row on top)
  res$term_label  <- factor(res$term_label, levels = rev(unique(res$term_label)))
  p <- ggplot(res, aes(x = estimate, y = term_label)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7, colour = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = sig_colour),
                   height = 0.25, linewidth = 1) +
    geom_point(aes(colour = sig_colour), size = 3.5) +
    # significance stars just above each significant point
    geom_text(aes(label = stars, colour = sig_colour), nudge_y = 0.25,
              size = 5, family = "Times New Roman", fontface = "bold") +
    scale_colour_identity() +
    # always draw every y label, even when they visually touch
    scale_y_discrete(guide = guide_axis(check.overlap = FALSE)) +
    scale_x_continuous(breaks = x_breaks) +
    labs(x = "Effect on log transport CO2", y = NULL) +
    theme_minimal(base_size = 12) +
    theme.adjusted
  if (!is.null(file)) {
    ggsave(file.path(out_dir, "figures", file), p,
           width = 10, height = height, dpi = 300)
  }
  p
}
