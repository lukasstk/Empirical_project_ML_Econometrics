# =============================================================================
# 02_main_effects_dml.R
# Questions 1 & 3: average policy effects and super-/sub-additivity.
#
# Two DML runs on the same partially linear model, differing only in how the
# treatment regimes are parametrized:
#
#   Run 1 (regime parametrization):  D = (only CP, only LEZ, both)
#     mutually exclusive regime dummies against "no policy".
#     - only_cp / only_lez: effect of each policy ALONE              -> question 1
#     - both: total effect of running both policies vs. neither, as a single
#       coefficient with its own directly computed CI - no linear combination
#       of estimates (and hence no coefficient covariances) needed
#                                                                -> questions 1 & 3
#
#   Run 2 (interaction parametrization):  D = (CP, LEZ, CP*LEZ)
#     - theta_int: deviation of the joint effect from the sum of the
#       individual effects (< 0: super-additive synergy)        -> question 3
#
# Estimated with cross-fitted DML (DoubleML package, dml_plr in 00_setup.R):
# each treatment coefficient is estimated in turn, with the other treatment
# columns folded into that run's nuisance/control set alongside W - this is
# DoubleML's default (use_other_treat_as_covariate = TRUE) and exactly the
# lecture's "One-By-One Double LASSO" procedure for a vector of target
# coefficients (see dml_plr in 00_setup.R for the precise correspondence).
# Plugin-lasso nuisances (hdm::rlasso, lambda from the Belloni/Chernozhukov/
# Hansen formula - one fit per nuisance instead of an inner 10-fold CV;
# CV-tuned lasso, ridge, and random forest are checked as alternative
# learners in 04), city-level folds, cluster-robust SEs. Repeated
# cross-fitting (n_rep below: 1 for development, 5 for the final run) puts
# split noise into the CIs.
# =============================================================================

source("00_setup.R")
prep <- readRDS(file.path(out_dir, "prepared_data.rds"))
data    <- prep$data

# Repeated cross-fitting splits: 5 puts split noise into the CIs (see
# dml_plr, 00_setup.R); cheap with the rlasso learner. Set to 1 only if a
# quick development run is needed.
n_rep <- 5

# ---- (1) Build outcome, treatments, controls --------------------------------
Y <- data$log_transport_co2
D <- cbind(cp_active  = data$cp_active,
           lez_active = data$lez_active,
           cp_x_lez   = data$cp_x_lez)

# High-dimensional control matrix:
# main effects + all pairwise interactions + squared terms
# + year dummies + country dummies (see build_W in 00_setup.R)
W <- build_W(data, prep$ctrl_baseline, prep$sq_vars)
cat("Control matrix W:", nrow(W), "x", ncol(W), "\n")

# ---- (2) Run 1: regime parametrization (single- and joint-policy effects) ---
D_regime <- cbind(only_cp  = data$cp_active  * (1 - data$lez_active),
                  only_lez = data$lez_active * (1 - data$cp_active),
                  both     = data$cp_active  * data$lez_active)
regime <- dml_plr(Y, D_regime, W, cluster = data$city_id,
                  learner = "rlasso", n_folds = 5, n_rep = n_rep, seed = 42)

res_regime <- regime$results
res_regime$pct_effect    <- pct(res_regime$estimate)
res_regime$pct_conf.low  <- pct(res_regime$conf.low)
res_regime$pct_conf.high <- pct(res_regime$conf.high)
save_table(res_regime, "tab_policy_regimes")

# ---- (3) Run 2: interaction parametrization (super-/sub-additivity) ---------
main <- dml_plr(Y, D, W, cluster = data$city_id,
                learner = "rlasso", n_folds = 5, n_rep = n_rep, seed = 42)

res_main <- main$results
res_main$pct_effect    <- pct(res_main$estimate)
res_main$pct_conf.low  <- pct(res_main$conf.low)
res_main$pct_conf.high <- pct(res_main$conf.high)
save_table(res_main, "tab_main_effects_dml")

# Consistency check across parametrizations (should hold approximately):
#   only_cp  ~ theta_cp,  only_lez ~ theta_lez,
#   both     ~ theta_cp + theta_lez + theta_int
# See tab_policy_regimes.csv ("both" row) vs. tab_main_effects_dml.csv
# (cp_active + lez_active + cp_x_lez) for the report.

# ---- (4) Coefficient plots ---------------------------------------------------
# x_breaks: fixed 0.5-step ticks, not the default ~10-target pretty_breaks -
# these two plots' narrower estimate range made that default too dense/uneven.
plot_effects(res_regime,
             "Policy regimes vs. no policy: congestion pricing (CP) and low-emission zone (LEZ)",
             "fig_policy_regimes.png",
             x_breaks = scales::breaks_width(0.5))
plot_effects(res_main,
             "Average policy effects: congestion pricing (CP) and low-emission zone (LEZ)",
             "fig_main_effects.png",
             x_breaks = scales::breaks_width(0.5))
