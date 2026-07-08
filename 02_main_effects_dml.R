# =============================================================================
# 02_main_effects_dml.R
# Questions 1 & 3: average policy effects and super-/sub-additivity.
#
# Two DML runs on the same partially linear model, differing only in how the
# treatment regimes are parametrized:
#
#   Run 1 (interaction parametrization):  D = (CP, LEZ, CP*LEZ)
#     - theta_cp / theta_lez: effect of each policy ALONE (the interaction
#       term is 0 when only one policy is active)               -> question 1
#     - theta_int: deviation of the joint effect from the sum of the
#       individual effects (< 0: super-additive synergy)        -> question 3
#
#   Run 2 (regime parametrization):  D = (only CP, only LEZ, both)
#     mutually exclusive regime dummies against "no policy". The coefficient
#     on `both` IS the total effect of running both policies vs. neither,
#     as a single coefficient with its own package-computed CI - no linear
#     combination of estimates (and hence no coefficient covariances) needed.
#
# Estimated with cross-fitted DML (DoubleML package), CV-lasso nuisances,
# city-level folds, cluster-robust SEs; n_rep = 5 repeated cross-fitting
# splits, so split noise is included in the confidence intervals
# (see dml_plr in 00_setup.R).
# =============================================================================

source("00_setup.R")
prep <- readRDS("output/prepared_data.rds")
data    <- prep$data

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

# ---- (2) Run 1: interaction parametrization ---------------------------------
main <- dml_plr(Y, D, W, cluster = data$city_id,
                learner = "lasso", n_folds = 5, n_rep = 5, seed = 42)

res_main <- main$results
res_main$pct_effect    <- pct(res_main$estimate)
res_main$pct_conf.low  <- pct(res_main$conf.low)
res_main$pct_conf.high <- pct(res_main$conf.high)
save_table(res_main, "tab_main_effects_dml")

# ---- (3) Run 2: regime parametrization (total effect of both policies) ------
D_regime <- cbind(only_cp  = data$cp_active  * (1 - data$lez_active),
                  only_lez = data$lez_active * (1 - data$cp_active),
                  both     = data$cp_active  * data$lez_active)
regime <- dml_plr(Y, D_regime, W, cluster = data$city_id,
                  learner = "lasso", n_folds = 5, n_rep = 5, seed = 42)

res_regime <- regime$results
res_regime$pct_effect    <- pct(res_regime$estimate)
res_regime$pct_conf.low  <- pct(res_regime$conf.low)
res_regime$pct_conf.high <- pct(res_regime$conf.high)
save_table(res_regime, "tab_policy_regimes")

# Consistency check across parametrizations (should hold approximately):
#   only_cp  ~ theta_cp,  only_lez ~ theta_lez,
#   both     ~ theta_cp + theta_lez + theta_int
both_row <- res_regime[res_regime$term == "both", ]

# ---- (4) Coefficient plots ---------------------------------------------------
plot_effects(res_main,
             "Average policy effects (cross-fitted DML, 95% CI)",
             "fig_main_effects.png")
plot_effects(res_regime,
             "Policy regimes vs. no policy (cross-fitted DML, 95% CI)",
             "fig_policy_regimes.png")

# ---- (5) Console interpretation helper --------------------------------------
cat("\n--- Interpretation (log points -> percent: 100*(exp(b)-1)) ---\n")
for (i in seq_len(nrow(res_main))) {
  cat(sprintf("%-12s % .4f  ->  %+.1f%% emissions [%.1f%%, %.1f%%]\n",
              res_main$term[i], res_main$estimate[i], res_main$pct_effect[i],
              res_main$pct_conf.low[i], res_main$pct_conf.high[i]))
}
cat(sprintf("%-12s % .4f  ->  %+.1f%% emissions [%.1f%%, %.1f%%]\n",
            "both(total)", both_row$estimate, both_row$pct_effect,
            both_row$pct_conf.low, both_row$pct_conf.high))
cat("theta_int < 0 => joint implementation MORE effective than additive.\n")
