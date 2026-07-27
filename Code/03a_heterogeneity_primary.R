# =============================================================================
# 03a_heterogeneity_primary.R
# Question 2 (primary specification): treatment interacted with a small,
# pre-specified set of city characteristics.
# Companion: 03b_heterogeneity_all_controls.R (all baseline controls
# interacted; robustness appendix, run manually - not part of the 05
# pipeline).
# =============================================================================

source("Code/00_setup.R")
prep <- readRDS(file.path(out_dir, "prepared_data.rds"))
data    <- prep$data

Y <- data$log_transport_co2
W <- build_W(data, prep$ctrl_baseline, prep$sq_vars)

# Pre-specified heterogeneity characteristics, chosen ex ante.
# log_area_km2 excluded (exact identity with log_population - log_pop_density).
het_vars_primary <- c("log_population", "log_gdp_pc", "log_pop_density",
                      "public_transit_score", "political_green", "fuel_price",
                      "log_tourism_intensity", "logistics_activity", "coastal")

Zc <- scale(as.matrix(data[het_vars_primary]), center = TRUE, scale = FALSE)
colnames(Zc) <- make.names(colnames(Zc), unique = TRUE)

# ---- Linear CATE approximation (interaction) with cross-fitted DML ----------
D_cp_int  <- data$cp_active  * Zc
colnames(D_cp_int)  <- paste0("cp_x_",  colnames(Zc))
D_lez_int <- data$lez_active * Zc
colnames(D_lez_int) <- paste0("lez_x_", colnames(Zc))

D_het <- cbind(cp_active  = data$cp_active,
               lez_active = data$lez_active,
               cp_x_lez   = data$cp_x_lez,
               D_cp_int, D_lez_int)

het_primary <- dml_plr(Y, D_het, W, cluster = data$city_id,
                       learner = "rlasso", n_folds = 5, n_rep = 5, seed = 42)
save_table(het_primary$results, "tab_heterogeneity_primary")

# ---- Coefficient plots (one per policy) --------------------------------------
res_cp  <- het_primary$results[grepl("^cp_x_",  het_primary$results$term) &
                               het_primary$results$term != "cp_x_lez", ]
res_lez <- het_primary$results[grepl("^lez_x_", het_primary$results$term), ]

plot_effects(res_cp,  "fig_heterogeneity_primary_cp.png")
plot_effects(res_lez, "fig_heterogeneity_primary_lez.png")

cat("\nPrimary heterogeneity spec done (", length(het_vars_primary),
    "pre-specified characteristics).\n")
