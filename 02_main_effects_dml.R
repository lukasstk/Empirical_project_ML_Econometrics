# =============================================================================
# 02_main_effects_dml.R
# Question 1: Were the policies effective in reducing emissions?
# Question 3: Is the joint effect larger than the sum of individual effects?
#
# Partially linear model with two treatments and their interaction:
#   log(CO2_it) = theta_cp * CP_it + theta_lez * LEZ_it
#                 + theta_int * (CP_it * LEZ_it) + g(W_it) + e_it
#
#   theta_cp  : effect of congestion pricing alone
#   theta_lez : effect of a low-emission zone alone
#   theta_int : deviation of the joint effect from additivity (question 3);
#               theta_int < 0  -> super-additive (synergy),
#               theta_int > 0  -> sub-additive
#
# Estimated with cross-fitted DML, following the partialling-out procedure
# from the lecture and Code.R, extended to two treatments and their
# interaction. Nuisance functions use the plugin-lasso (hdm::rlasso,
# theory-based lambda - the same penalty rule as rlassoEffects() in the
# lecture); the outer K-fold city-level cross-fitting is kept. The
# cross-validated lasso appears as a comparison in 04_sensitivity_analysis.R.
# =============================================================================

source("00_setup.R")
prep <- readRDS("output/prepared_data.rds")
data    <- prep$data

# ---- (1) Build outcome, treatments, controls --------------------------------
Y <- data$log_transport_co2
D <- cbind(cp_active = data$cp_active,
           lez_active = data$lez_active,
           cp_x_lez   = data$cp_x_lez)

# High-dimensional control matrix:
# main effects + all pairwise interactions + squared terms
# + year dummies + country dummies (see build_W in 00_setup.R)
W <- build_W(data, prep$ctrl_baseline, prep$sq_vars)
cat("Control matrix W:", nrow(W), "x", ncol(W), "\n")

# ---- (2) Main estimate: cross-fitted DML with plugin-lasso ------------------
main <- dml_plm(Y, D, W, cluster = data$city_id,
                learner = "rlasso", K = 5, seed = 42)

res_main <- main$results
res_main$pct_effect   <- pct(res_main$estimate)
res_main$pct_conf.low <- pct(res_main$conf.low)
res_main$pct_conf.high<- pct(res_main$conf.high)
save_table(res_main, "tab_main_effects_dml")

# Total effect of implementing BOTH policies (linear combination):
# theta_cp + theta_lez + theta_int, with delta-method CI
both <- lincom(main, c(cp_active = 1, lez_active = 1, cp_x_lez = 1))
both <- cbind(term = "both policies (total)", both,
              pct_effect = pct(both$estimate))
save_table(both, "tab_both_policies_total")

# ---- (3) Coefficient plot ----------------------------------------------------
plot_effects(res_main,
             "Average policy effects (cross-fitted DML, 95% CI)",
             "fig_main_effects.png")

# ---- (4) Console interpretation helper --------------------------------------
cat("\n--- Interpretation (log points -> percent: 100*(exp(b)-1)) ---\n")
for (i in seq_len(nrow(res_main))) {
  cat(sprintf("%-12s % .4f  ->  %+.1f%% emissions [%.1f%%, %.1f%%]\n",
              res_main$term[i], res_main$estimate[i], res_main$pct_effect[i],
              res_main$pct_conf.low[i], res_main$pct_conf.high[i]))
}
cat(sprintf("%-12s % .4f  ->  %+.1f%% emissions\n",
            "both(total)", both$estimate, both$pct_effect))
cat("theta_int < 0 => joint implementation MORE effective than additive.\n")
