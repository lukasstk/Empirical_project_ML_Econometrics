# =============================================================================
# 04_sensitivity_analysis.R
# Robustness of the main results to the modeling choices:
#
#   (1) ML learner for the nuisance functions:
#       plugin lasso (baseline, as in 02/03) / CV lasso / ridge / random forest
#   (2) Control set: baseline / minimal (no interactions) / + mediators
#   (3) Number of cross-fitting folds: K = 2 / 5 (baseline) / 10
#
# All checks reuse the dml_plm estimator from the lecture material; no
# additional methods are introduced.
# =============================================================================

source("00_setup.R")
prep <- readRDS("output/prepared_data.rds")
data    <- prep$data

Y <- data$log_transport_co2
D <- cbind(cp_active  = data$cp_active,
           lez_active = data$lez_active,
           cp_x_lez   = data$cp_x_lez)

W_base <- build_W(data, prep$ctrl_baseline, prep$sq_vars)
# Main-effects-only matrix (no interactions):
#   used for the random forest (captures nonlinearities itself) and for
#   the "minimal controls" specification
W_main <- build_W(data, prep$ctrl_baseline, interactions = FALSE)
# Extended matrix including potential mediators (bad controls)
W_ext  <- build_W(data, c(prep$ctrl_baseline, prep$ctrl_mediators), prep$sq_vars)

sens <- list()
run_spec <- function(label, Yv, Dv, Wv, ...) {
  out <- dml_plm(Yv, Dv, Wv, cluster = data$city_id, seed = 42, ...)
  cbind(out$results, spec = label)
}

# ---- (1) Learner choice ------------------------------------------------------
# The baseline learner is the plugin-lasso, as in 02 and 03. The CV-lasso,
# ridge and random forest rows check that the results do not hinge on that
# choice of nuisance learner.
# The CV-based rows use the one-standard-error rule (lambda.1se): among all
# lambdas whose CV error is within one standard error of the minimum, take
# the largest (most strongly regularized), as recommended in the lecture to
# avoid overly complex models despite CV.
sens$plugin <- run_spec("plugin lasso (baseline)",  Y, D, W_base, learner = "rlasso")
sens$lasso  <- run_spec("lasso, CV lambda.1se",     Y, D, W_base, learner = "lasso")
sens$ridge  <- run_spec("ridge, CV lambda.1se",     Y, D, W_base, learner = "ridge")
sens$rf     <- run_spec("random forest",            Y, D, W_main, learner = "rf")

# ---- (2) Control sets --------------------------------------------------------
# Same learner as the baseline, so each row changes exactly one thing.
sens$min   <- run_spec("minimal controls (main effects only)",
                       Y, D, W_main, learner = "rlasso")
sens$ext   <- run_spec("+ mediators (pm25, fleet shares) [bad controls]",
                       Y, D, W_ext,  learner = "rlasso")

# ---- (3) Number of cross-fitting folds ---------------------------------------
# K trades off training data per nuisance fit (K = 2: only half the cities
# to learn from) against computation. The DML theory allows any fixed K, so
# the estimates should be stable; this justifies the baseline choice K = 5.
sens$k2  <- run_spec("cross-fitting: K = 2",  Y, D, W_base, learner = "rlasso", K = 2)
sens$k10 <- run_spec("cross-fitting: K = 10", Y, D, W_base, learner = "rlasso", K = 10)

sens_tab <- do.call(rbind, sens)
rownames(sens_tab) <- NULL
save_table(sens_tab, "tab_sensitivity_specs")

p_sens <- ggplot(sens_tab, aes(x = estimate, y = spec)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  facet_wrap(~ term, scales = "free_x") +
  labs(x = "Effect on log transport CO2 (95% CI)", y = NULL,
       title = "Sensitivity of the main effects to modeling choices") +
  theme_minimal(base_size = 11)
ggsave("output/figures/fig_sensitivity.png", p_sens, width = 11, height = 6)

cat("\nSensitivity analysis done. Check that:\n",
    " - estimates are similar across learners, control sets, and fold counts\n",
    " - the +mediators row shrinks the LEZ effect (expected: bad controls)\n")
