# =============================================================================
# 04_sensitivity_analysis.R
# Robustness of the main results to the modeling choices:
#
#   (1) ML learner for the nuisance functions:
#       CV-lasso lambda.1se (baseline, as in 02/03) / CV-lasso lambda.min /
#       ridge / random forest
#   (2) Control set: baseline / minimal (no interactions) / + mediators
#   (3) Number of cross-fitting folds: n_folds = 2 / 5 (baseline) / 10
#   (4) Placebo treatment: policy announced (major topic of public
#       discussion) but not yet implemented - tests identification itself
#
# All checks use the same DoubleML estimator as 02/03 (dml_plr, 00_setup.R).
# Split noise from the random fold assignment needs no separate check here:
# the main estimates in 02 use n_rep = 5 repeated cross-fitting, which
# includes the across-split dispersion in the reported CIs. The rows below
# use n_rep = 1 to keep the runtime manageable; each row changes exactly
# one thing relative to the baseline.
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
  out <- dml_plr(Yv, Dv, Wv, cluster = data$city_id, seed = 42, ...)
  cbind(out$results, spec = label)
}

# ---- (1) Learner choice ------------------------------------------------------
# The baseline learner is the CV-lasso with the one-standard-error rule
# (lambda.1se): among all lambdas whose CV error is within one standard
# error of the minimum, take the largest (most strongly regularized), as
# recommended in the lecture to avoid overly complex models despite CV.
# The other rows check that the results do not hinge on that choice:
# lambda.min (less regularization), ridge (no selection), random forest
# (fully nonparametric).
sens$lasso  <- run_spec("CV lasso, lambda.1se (baseline)", Y, D, W_base,
                        learner = "lasso")
sens$lmin   <- run_spec("CV lasso, lambda.min",            Y, D, W_base,
                        learner = "lasso_min")
sens$ridge  <- run_spec("ridge, lambda.1se",               Y, D, W_base,
                        learner = "ridge")
sens$rf     <- run_spec("random forest",                   Y, D, W_main,
                        learner = "rf")

# ---- (2) Control sets --------------------------------------------------------
# Same learner as the baseline, so each row changes exactly one thing.
sens$min   <- run_spec("minimal controls (main effects only)",
                       Y, D, W_main, learner = "lasso")
sens$ext   <- run_spec("+ mediators (pm25, fleet shares) [bad controls]",
                       Y, D, W_ext,  learner = "lasso")

# ---- (3) Number of cross-fitting folds ---------------------------------------
# n_folds trades off training data per nuisance fit (n_folds = 2: only half
# the cities to learn from) against computation. The DML theory allows any
# fixed number of folds, so the estimates should be stable; this justifies
# the baseline choice n_folds = 5.
sens$k2  <- run_spec("cross-fitting: 2 folds",  Y, D, W_base,
                     learner = "lasso", n_folds = 2)
sens$k10 <- run_spec("cross-fitting: 10 folds", Y, D, W_base,
                     learner = "lasso", n_folds = 10)

# ---- (4) Placebo treatment: announced but not yet active ---------------------
# Blocks (1)-(3) vary the ESTIMATOR (learner, controls, folds). This check
# targets the IDENTIFYING ASSUMPTION (conditional ignorability), which no
# estimator choice can fix. Between the year a policy first became a major
# topic of local public discussion and its implementation year (~1.5 years
# on average in this data), the policy cannot mechanically affect emissions:
# no charge is collected, no vehicle is banned. A placebo indicator for
# exactly these city-years should therefore have a coefficient close to zero.
#
# A clearly negative placebo coefficient would mean emissions fall BEFORE the
# policy exists, which has two competing readings:
#   - anticipation: households and firms adapt early (e.g. replacing a diesel
#     car before the LEZ starts). The decline is still policy-caused, but the
#     pre-implementation comparison years are then partly treated, so the
#     main estimates UNDERSTATE the total policy effect.
#   - selection / awareness: the public debate reflects (and further fuels)
#     general environmental awareness in the city - people behave more
#     mindfully and reduce driving and emissions regardless of the future
#     policy, and greener local governments typically push other measures
#     (transit, cycling, parking) at the same time. This decline is NOT
#     caused by the policy instruments, so the main estimates OVERSTATE
#     their effect.
# The data cannot fully separate the two (every city that discussed a policy
# eventually implemented it), but a near-zero placebo rules BOTH out, and a
# negative one bounds the main estimates from one side or the other.
#
# The placebo indicators enter JOINTLY with the actual treatment indicators,
# so their coefficients are net of the active-policy effects.
plc <- function(announce, impl) {
  as.integer(announce > 0 & data$year >= announce &
               (impl == 0 | data$year < impl))
}
D_plac <- cbind(D,
                cp_pre  = plc(data$cp_announce_year,  data$cp_impl_year),
                lez_pre = plc(data$lez_announce_year, data$lez_impl_year))
placebo <- dml_plr(Y, D_plac, W_base, cluster = data$city_id,
                   learner = "lasso", seed = 42)
save_table(placebo$results, "tab_placebo_announcement")

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
    " - the +mediators row shrinks the LEZ effect (expected: bad controls)\n",
    " - cp_pre / lez_pre in tab_placebo_announcement are close to zero\n",
    "   (identification check; see block (4) comment for what a negative\n",
    "   coefficient would mean: anticipation vs. selection/awareness)\n")
