# =============================================================================
# 00_setup.R
# Packages, output folders, and helper functions used by all other scripts.
#
# Project: Effects of congestion pricing (CP) and low-emission zones (LEZ)
#          on urban road-transport CO2 emissions.
# Method:  Double/Debiased Machine Learning (DML) with cross-fitting,
#          following the partialling-out logic from the wage-gap example
#          (Code.R), extended to panel data with two treatments.
# =============================================================================

packages <- c("glmnet", "ranger", "hdm",
              "ggplot2", "dplyr", "tidyr")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

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
  # syntactic, unique column names (":" -> "."), consistent with D in dml_plm
  colnames(W) <- make.names(colnames(W), unique = TRUE)
  scale(W, center = TRUE, scale = FALSE)               # demean each column
}

# -----------------------------------------------------------------------------
# vcov_cluster: cluster-robust covariance matrix for an lm fit
# (Liang-Zeger / "clustered standard errors"), computed by hand.
# NOT USED BY DEFAULT - optional extension beyond the course material,
# kept for reference; see the commented line in dml_plm.
#
#   V = c * (X'X)^-1 [ sum_g s_g s_g' ] (X'X)^-1,   s_g = sum_{i in g} x_i e_i
#
# Intuition: the middle ("meat") sums the score contributions x_i * e_i
# WITHIN each cluster first, so errors of the same city may be arbitrarily
# correlated; only different cities are assumed independent. The outer
# factors ("bread") are the usual OLS (X'X)^-1 - hence "sandwich" estimator.
# c = G/(G-1) * (n-1)/(n-k) is the standard small-sample correction
# (G = number of clusters). With independent errors this reduces to a
# heteroskedasticity-robust estimator; the classical lm vcov would instead
# use sigma^2 (X'X)^-1.
# -----------------------------------------------------------------------------
vcov_cluster <- function(fit, cluster) {
  X <- model.matrix(fit)
  e <- residuals(fit)
  n <- nrow(X)
  k <- ncol(X)
  G <- length(unique(cluster))
  S     <- rowsum(X * e, group = as.character(cluster))  # G x k score sums
  meat  <- crossprod(S)                                  # sum_g s_g s_g'
  bread <- solve(crossprod(X))                           # (X'X)^-1
  adj   <- G / (G - 1) * (n - 1) / (n - k)
  adj * bread %*% meat %*% bread
}

# -----------------------------------------------------------------------------
# dml_plm: cross-fitted DML estimator for the partially linear model
#
#   Y = D %*% theta + g(W) + e,   E[e | D, W] = 0
#
# D can contain several treatment columns at once (here: cp_active,
# lez_active, their interaction, and possibly treatment x characteristic
# interactions for the heterogeneity analysis).
#
# Steps (Neyman-orthogonal partialling out with K-fold cross-fitting):
#   1. Split CITIES (not city-years!) into K folds, so that all years of a
#      city stay in the same fold. This prevents leakage through serial
#      correlation within a city.
#   2. For each fold k: learn E[Y|W] and E[D_j|W] on the other folds,
#      predict on fold k, and form out-of-fold residuals.
#   3. Final OLS of residualized Y on residualized D.
#   4. Inference via summary() and confint() on the final OLS, as in the
#      lecture examples (a cluster-robust alternative is sketched in a
#      comment below; see vcov_cluster).
#
# learner:
#   "lasso"  - cv.glmnet, alpha = 1, lambda chosen by cross-validation
#              (10-fold CV inside every training step -> 11 path fits per call)
#   "rlasso" - hdm::rlasso with the theory-based plugin lambda, as used by
#              rlassoEffects() in the lecture (Code.R). No inner CV: lambda
#              comes from a closed-form rule, so each call is a single fit
#              (~11x fewer lasso fits than "lasso"). The outer cross-fitting
#              is unchanged.
#   "ridge"  - cv.glmnet, alpha = 0
#   "rf"     - random forest (ranger); pass a W without interactions,
#              since the forest captures nonlinearities itself
# lambda_rule: "lambda.1se" (default) - one-standard-error rule: the most
#              strongly regularized model within one SE of the CV-error
#              minimum, as recommended in the lecture. "lambda.min" uses the
#              CV minimum directly. Only relevant for "lasso"/"ridge".
# -----------------------------------------------------------------------------
dml_plm <- function(Y, D, W, cluster, learner = c("lasso", "rlasso", "ridge", "rf"),
                    K = 5, seed = 1, lambda_rule = "lambda.1se",
                    level = 0.95) {
  learner <- match.arg(learner)
  D <- as.matrix(D)
  # D names must stay syntactic: data.frame()/lm() below mangle non-syntactic
  # names, which would break the coefficient lookup via colnames(D)
  colnames(D) <- make.names(colnames(D), unique = TRUE)
  n <- length(Y)
  stopifnot(nrow(D) == n, nrow(W) == n, length(cluster) == n)

  # IMPORTANT: FOLDS ARE ASSIGNED AT THE CITY LEVEL, SO ALL YEARS OF A CITY
  # STAY IN THE SAME FOLD. THIS PREVENTS WITHIN-CITY LEAKAGE ACROSS FOLDS.
  set.seed(seed)
  ids <- unique(as.character(cluster))
  id_fold <- sample(rep_len(1:K, length(ids)))
  names(id_fold) <- ids
  fold <- id_fold[as.character(cluster)]

  # Standard observation-level alternative (each city-year is assigned separately):
  # fold <- sample(rep_len(1:K, n))

  fit_predict <- function(W_tr, y_tr, W_te) {
    if (learner %in% c("lasso", "ridge")) {
      fit <- cv.glmnet(W_tr, y_tr, alpha = ifelse(learner == "lasso", 1, 0))
      as.numeric(predict(fit, newx = W_te, s = lambda_rule))
    } else if (learner == "rlasso") {
      # Plugin penalty instead of CV: lambda from the theoretical rule of
      # Belloni/Chernozhukov/Hansen (same choice as rlassoEffects in Code.R).
      # post = FALSE: use the lasso coefficients directly for prediction,
      # without the post-lasso OLS refit (as in the first rlasso call in Code.R).
      fit <- rlasso(x = W_tr, y = y_tr, post = FALSE)
      as.numeric(predict(fit, newdata = W_te))
    } else {
      fit <- ranger(y = y_tr, x = as.data.frame(W_tr),
                    num.trees = 500, min.node.size = 5, seed = seed)
      predict(fit, data = as.data.frame(W_te))$predictions
    }
  }

  # out-of-fold residuals for Y and every treatment column
  Y_res <- rep(NA_real_, n)
  D_res <- matrix(NA_real_, n, ncol(D), dimnames = list(NULL, colnames(D)))
  for (k in 1:K) {
    cat(sprintf("  fold %d/%d: Y", k, K)); flush.console()   # <- before the Y fit
    tr <- fold != k
    te <- fold == k
    Y_res[te] <- Y[te] -
      fit_predict(W[tr, , drop = FALSE], Y[tr], W[te, , drop = FALSE])
    for (j in seq_len(ncol(D))) {
      cat(" |", colnames(D)[j]); flush.console()             # <- before each D column fit
      D_res[te, j] <- D[te, j] -
        fit_predict(W[tr, , drop = FALSE], D[tr, j], W[te, , drop = FALSE])
    }
    cat(" âœ“\n")                                              # <- fold finished
  }

  # final partialling-out regression (residualized Y on residualized D)
  dat <- data.frame(Y_res = Y_res, D_res)
  fit <- lm(Y_res ~ ., data = dat)
  # NOTE (beyond course material): the classical standard errors below
  # assume independent errors. With repeated observations of the same city
  # one can instead use the cluster-robust covariance matrix (vcov_cluster
  # above). To do so, replace the summary()/confint() block with:
  # V    <- vcov_cluster(fit, cluster)
  # keep <- colnames(D)
  # b    <- coef(fit)[keep]
  # se   <- sqrt(diag(V))[keep]
  # z    <- qnorm(1 - (1 - level) / 2)
  # results <- data.frame(term      = keep,
  #                       estimate  = b,
  #                       std.error = se,
  #                       p.value   = 2 * pnorm(-abs(b / se)),
  #                       conf.low  = b - z * se,
  #                       conf.high = b + z * se,
  #                       row.names = NULL)

  # Inference directly from summary() and confint(), as in Code.R:
  ct <- summary(fit)$coefficients      # Estimate | Std. Error | t | Pr(>|t|)
  ci <- confint(fit, level = level)
  keep <- colnames(D)
  results <- data.frame(term      = keep,
                        estimate  = ct[keep, 1],
                        std.error = ct[keep, 2],
                        p.value   = ct[keep, 4],
                        conf.low  = ci[keep, 1],
                        conf.high = ci[keep, 2],
                        row.names = NULL)
  list(results = results, fit = fit,
       Y_res = Y_res, D_res = D_res, fold = fold)
}

# -----------------------------------------------------------------------------
# lincom: linear combination of DML coefficients with delta-method SE
# Example: total effect of running both policies =
#          theta_cp + theta_lez + theta_interaction
# L is a named vector of weights, e.g. c(cp_active = 1, lez_active = 1, ...)
# -----------------------------------------------------------------------------
lincom <- function(dml_out, L, level = 0.95) {
  cf <- coef(dml_out$fit)
  # The variance of a SUM of coefficients needs their covariances too:
  #   Var(L'b) = L' V L, with V the full coefficient covariance matrix.
  # summary() only reports the diagonal of V (the SEs), so this is the one
  # place where vcov() is required.
  V  <- vcov(dml_out$fit)
  Lf <- setNames(rep(0, length(cf)), names(cf))
  stopifnot(all(names(L) %in% names(Lf)))
  Lf[names(L)] <- L
  est <- sum(Lf * cf)
  se  <- sqrt(drop(t(Lf) %*% V %*% Lf))
  z   <- qnorm(1 - (1 - level) / 2)
  data.frame(estimate = est, std.error = se,
             conf.low = est - z * se, conf.high = est + z * se,
             p.value  = 2 * pnorm(-abs(est / se)))
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

