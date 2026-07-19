# =============================================================================
# 02_main_and_joint_effects.R
# Questions 1 & 3: average policy effects and super-/sub-additivity.
# Two parametrizations of the same partially linear model with city fixed
# effects (same Y, same W, coefficients are linear recombinations of each
# other):
#   (a) interaction parametrization D = (CP, LEZ, CP*LEZ):
#       theta_cp / theta_lez answer Q1 for each policy alone and the
#       interaction theta_int answers Q3 (super-/sub-additivity) directly.
#   (b) regime parametrization D = (only CP, only LEZ, both), base category
#       "neither policy": theta_both answers Q1's "both policies vs. no
#       policy" with a proper SE/CI, which (a) can only deliver as a
#       CI-less sum of coefficients (cross-coefficient covariances are not
#       reported by DoubleML).
# =============================================================================

source("Code/00_setup.R")
prep <- readRDS(file.path(out_dir, "prepared_data.rds"))
data    <- prep$data

# Repeated cross-fitting splits (5 = final run, 1 = quick development run)
n_rep <- 5

# ---- (1) Build outcome, treatments, controls --------------------------------
Y <- data$log_transport_co2
D <- cbind(cp_active  = data$cp_active,
           lez_active = data$lez_active,
           cp_x_lez   = data$cp_x_lez)

W <- build_W(data, prep$ctrl_baseline, prep$sq_vars)
cat("Control matrix W:", nrow(W), "x", ncol(W), "\n")

# ---- (2) DML run (a): interaction parametrization ---------------------------
main <- dml_plr(Y, D, W, cluster = data$city_id,
                learner = "rlasso", n_folds = 5, n_rep = n_rep, seed = 42)

res_main <- main$results
res_main$pct_effect    <- pct(res_main$estimate)
res_main$pct_conf.low  <- pct(res_main$conf.low)
res_main$pct_conf.high <- pct(res_main$conf.high)
save_table(res_main, "tab_main_effects")

# ---- (3) DML run (b): regime parametrization (question 1, joint effect) -----
# Mutually exclusive policy regimes; "both" is the total effect of running
# both policies vs. neither, with its own SE/CI.
D_reg <- cbind(only_cp  = data$cp_active  * (1 - data$lez_active),
               only_lez = data$lez_active * (1 - data$cp_active),
               both     = data$cp_active  * data$lez_active)

regimes <- dml_plr(Y, D_reg, W, cluster = data$city_id,
                   learner = "rlasso", n_folds = 5, n_rep = n_rep, seed = 42)

res_reg <- regimes$results
res_reg$pct_effect    <- pct(res_reg$estimate)
res_reg$pct_conf.low  <- pct(res_reg$conf.low)
res_reg$pct_conf.high <- pct(res_reg$conf.high)
save_table(res_reg, "tab_joint_effects_regimes")

# Consistency across parametrizations (equal up to cross-fitting noise):
#   theta_cp + theta_lez + theta_int  (a)  =  theta_both  (b)
#   theta_both - theta_only_cp - theta_only_lez  (b)  =  theta_int  (a)
both_hat <- sum(res_main$estimate)
cat(sprintf("\nBoth vs. neither -- (a) reconstructed sum: %.4f | (b) direct with CI: %.4f [%.4f, %.4f]\n",
            both_hat,
            res_reg$estimate[res_reg$term == "both"],
            res_reg$conf.low[res_reg$term == "both"],
            res_reg$conf.high[res_reg$term == "both"]))

# ---- (4) Coefficient plots (one per parametrization) ------------------------
plot_effects(res_main, "fig_main_effects.png",
             x_breaks = scales::breaks_width(0.25))
plot_effects(res_reg,  "fig_joint_effects_regimes.png",
             x_breaks = scales::breaks_width(0.25))
