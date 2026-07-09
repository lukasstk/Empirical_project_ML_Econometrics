# =============================================================================
# 03_heterogeneity_dml.R
# Question 2: Which type of cities benefit the most/least from the policies?
#
# Interaction approach, as in the lecture's wage-gap heterogeneity example:
# add treatment x characteristic interactions as additional target
# parameters. Characteristics are centered, so the main treatment
# coefficients stay the effect for an average city, and each interaction
# coefficient says how the effect shifts per unit above average.
#
#   log(CO2) = theta_cp*CP + theta_lez*LEZ + theta_int*CP*LEZ
#              + CP*(Z - mean(Z))'gamma_cp + LEZ*(Z - mean(Z))'gamma_lez
#              + g(W) + e
#
# Estimated with cross-fitted DML (DoubleML package, dml_plr in 00_setup.R),
# same procedure as 02 ("One-By-One Double LASSO"): each of the 3 + 2x9
# treatment columns is estimated in turn, with the remaining treatment
# columns folded into that run's nuisance/control set alongside W.
# Plugin-lasso nuisances, city-level folds, cluster-robust SEs. n_rep = 5,
# same as 02; 21 separate treatment columns is why this run is slower.
#
# A "type" of city is only defined through its observable characteristics
# (size, wealth, density, transit quality, ...), so "which type benefits
# most" is the same question as "how does the effect vary with each
# characteristic" - exactly what these interactions measure. Same
# construction as the lecture's `female` x education/experience
# interactions, just with CP/LEZ as treatments and city traits as
# characteristics.
#
# Reading the table: because characteristics are centered, the first three
# rows (cp_active, lez_active, cp_x_lez) are the effects for an average
# city, not the same numbers as in 02 - report the average effects from
# 02, not from here.
# =============================================================================

source("00_setup.R")
prep <- readRDS(file.path(out_dir, "prepared_data.rds"))
data    <- prep$data

Y <- data$log_transport_co2
W <- build_W(data, prep$ctrl_baseline, prep$sq_vars)

het_vars_c <- paste0(prep$het_vars, "_c")   # centered characteristics
Zc <- as.matrix(data[het_vars_c])

# ---- CATE interaction approach with cross-fitted DML ------------------------
# Treatment matrix: 3 base treatments + CP x Z and LEZ x Z interactions.
D_cp_int  <- data$cp_active  * Zc
colnames(D_cp_int)  <- paste0("cp_x_",  prep$het_vars)
D_lez_int <- data$lez_active * Zc
colnames(D_lez_int) <- paste0("lez_x_", prep$het_vars)

D_het <- cbind(cp_active  = data$cp_active,
               lez_active = data$lez_active,
               cp_x_lez   = data$cp_x_lez,
               D_cp_int, D_lez_int)

het <- dml_plr(Y, D_het, W, cluster = data$city_id,
               learner = "rlasso", n_folds = 5, n_rep = 5, seed = 42)
save_table(het$results, "tab_heterogeneity_dml")

# Coefficient plots, separately for CP and LEZ interactions
res_cp  <- het$results[grepl("^cp_x_",  het$results$term) &
                       het$results$term != "cp_x_lez", ]
res_lez <- het$results[grepl("^lez_x_", het$results$term), ]

plot_effects(res_cp,
             "Heterogeneity of the congestion pricing (CP) effect",
             "fig_het_cp.png")
plot_effects(res_lez,
             "Heterogeneity of the low-emission zone (LEZ) effect",
             "fig_het_lez.png")

