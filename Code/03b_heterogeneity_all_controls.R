# =============================================================================
# 03b_heterogeneity_all_controls.R
# Question 2 (robustness specification): treatment interacted with ALL
# baseline controls, including the placebo variables. Specification check
# for 03a_heterogeneity_primary.R - read jointly, not row by row.
#
# Not part of the 05 pipeline; run manually on an out_dir that already
# contains prepared_data.rds, e.g.:
#   Rscript -e 'out_dir <- "output"; source("Code/03b_heterogeneity_all_controls.R")'
# =============================================================================

source("Code/00_setup.R")
prep <- readRDS(file.path(out_dir, "prepared_data.rds"))
data    <- prep$data

Y <- data$log_transport_co2
W <- build_W(data, prep$ctrl_baseline, prep$sq_vars)

# Interaction set: all baseline controls except log_area_km2 (exact identity
# with log_population - log_pop_density); factor expanded to dummies,
# constant columns dropped, every column centered.
het_ctrl <- setdiff(prep$ctrl_baseline, "log_area_km2")
Z <- model.matrix(~ ., data = data[het_ctrl])[, -1, drop = FALSE]
Z <- Z[, apply(Z, 2, var) > 0, drop = FALSE]
Zc <- scale(Z, center = TRUE, scale = FALSE)
colnames(Zc) <- make.names(colnames(Zc), unique = TRUE)
cat("Heterogeneity characteristics Z:", ncol(Zc),
    "columns (log_area_km2 excluded)\n")

# ---- CATE interaction approach with cross-fitted DML ------------------------
D_cp_int  <- data$cp_active  * Zc
colnames(D_cp_int)  <- paste0("cp_x_",  colnames(Zc))
D_lez_int <- data$lez_active * Zc
colnames(D_lez_int) <- paste0("lez_x_", colnames(Zc))

D_het <- cbind(cp_active  = data$cp_active,
               lez_active = data$lez_active,
               cp_x_lez   = data$cp_x_lez,
               D_cp_int, D_lez_int)

# Slowest run in the project; drop to n_rep = 1 for a quick development check.
het_full <- dml_plr(Y, D_het, W, cluster = data$city_id,
                    learner = "rlasso", n_folds = 5, n_rep = 5, seed = 42)
save_table(het_full$results, "tab_heterogeneity_full")

# ---- Coefficient plots (one per policy) --------------------------------------
res_cp  <- het_full$results[grepl("^cp_x_",  het_full$results$term) &
                            het_full$results$term != "cp_x_lez", ]
res_lez <- het_full$results[grepl("^lez_x_", het_full$results$term), ]
fig_height <- max(4, 0.3 * nrow(res_cp) + 1.5)

plot_effects(res_cp,  "fig_heterogeneity_full_cp.png",  height = fig_height)
plot_effects(res_lez, "fig_heterogeneity_full_lez.png", height = fig_height)

# ---- Significance-count diagnostic --------------------------------------------
# Naive count of nominal 5% hits vs. the number expected by chance under a
# global null (not a formal joint test).
sig_count <- function(res, label) {
  n   <- nrow(res)
  hit <- sum(res$p.value < 0.05, na.rm = TRUE)
  cat(sprintf("%s: %d/%d terms significant at 5%% (%.1f expected by chance under a global null)\n",
              label, hit, n, 0.05 * n))
}
sig_count(res_cp,  "CP interactions")
sig_count(res_lez, "LEZ interactions")

cat("\nFull heterogeneity spec done (", ncol(Zc), "characteristics).\n")
