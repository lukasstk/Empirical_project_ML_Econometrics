# =============================================================================
# 04_sensitivity_analysis.R
# Robustness of the main results to the modeling choices:
#   (1) ML learner for the nuisance functions: rlasso (baseline) / ridge /
#       random forest
#   (2) Control set: baseline / minimal (no interactions) / + mediators
#   (3) Placebo treatment: policy announced but not yet implemented
# All checks use the same estimator as 02/03 (dml_plr, 00_setup.R); each
# row changes exactly one thing relative to the baseline.
# =============================================================================

source("Code/00_setup.R")
prep <- readRDS(file.path(out_dir, "prepared_data.rds"))
data    <- prep$data

Y <- data$log_transport_co2
D <- cbind(cp_active  = data$cp_active,
           lez_active = data$lez_active,
           cp_x_lez   = data$cp_x_lez)

W_base <- build_W(data, prep$ctrl_baseline, prep$sq_vars)
# Main-effects-only matrix (no interactions): minimal-controls spec
W_main <- build_W(data, prep$ctrl_baseline, interactions = FALSE)
# Baseline control set extended by the mediators
W_ext  <- build_W(data, c(prep$ctrl_baseline, prep$ctrl_mediators), prep$sq_vars)

sens <- list()
run_spec <- function(label, Yv, Dv, Wv, ...) {
  out <- dml_plr(Yv, Dv, Wv, cluster = data$city_id, seed = 42, ...)
  cbind(out$results, spec = label)
}

# ---- (1) Learner choice ------------------------------------------------------
sens$rlasso    <- run_spec("rlasso (baseline)", Y, D, W_base,
                           learner = "rlasso")
sens$ridge     <- run_spec("CV ridge (1se)",     Y, D, W_base,
                           learner = "ridge")
# Same W as the other learners, so this row changes only the learner.
# Slowest row of the table: ranger on the full 600+ column matrix.
sens$rf        <- run_spec("random forest",      Y, D, W_base,
                           learner = "rf")

# ---- (2) Control sets --------------------------------------------------------
sens$min   <- run_spec("minimal controls", Y, D, W_main, learner = "rlasso")
sens$ext   <- run_spec("with added mediators",
                       Y, D, W_ext,  learner = "rlasso")

# ---- (3) Placebo treatment: announced but not yet active ---------------------
# Indicator for city-years between announcement and implementation; enters
# jointly with the actual treatment indicators. The table does double duty:
# cp_pre / lez_pre are the placebo test, and because DoubleML folds the
# other d_cols into the covariates, the cp_active / lez_active rows are the
# main effects adjusted for the announcement period (anticipation).
plc <- function(announce, impl) {
  as.integer(announce > 0 & data$year >= announce &
               (impl == 0 | data$year < impl))
}
D_plac <- cbind(D,
                cp_pre  = plc(data$cp_announce_year,  data$cp_impl_year),
                lez_pre = plc(data$lez_announce_year, data$lez_impl_year))
placebo <- dml_plr(Y, D_plac, W_base, cluster = data$city_id,
                   learner = "rlasso", seed = 42)
save_table(placebo$results, "tab_sensitivity_placebo")

# ---- Combined sensitivity table -------------------------------------------
sens_tab <- do.call(rbind, sens)
rownames(sens_tab) <- NULL
save_table(sens_tab, "tab_sensitivity")

# ---- Sensitivity figures (one per check group) ----------------------------
# Pin spec order to the definition order above (reversed: baseline on top)
sens_tab$spec <- factor(sens_tab$spec, levels = rev(unique(sens_tab$spec)))

sens_tab$sig_colour <- classify_colour(sens_tab$estimate, sens_tab$conf.low, sens_tab$conf.high)
# Full names for the facet strips; explicit levels fix panel order to
# CP, LEZ, CP:LEZ top-to-bottom
facet_titles <- c(
  cp_active  = "Congestion pricing (CP)",
  lez_active = "Low-emission zone (LEZ)",
  cp_x_lez   = "CP:LEZ"
)
sens_tab$facet_label <- factor(facet_titles[sens_tab$term], levels = facet_titles)

# Split sens_tab into one plot per check group
spec_group_levels <- c("Learner sensitivity", "Control-set sensitivity")
spec_groups <- c(
  "rlasso (baseline)"     = "Learner sensitivity",
  "CV ridge (1se)"        = "Learner sensitivity",
  "random forest"         = "Learner sensitivity",
  "minimal controls"      = "Control-set sensitivity",
  "with added mediators"  = "Control-set sensitivity"
)
sens_tab$spec_group <- factor(spec_groups[as.character(sens_tab$spec)],
                              levels = spec_group_levels)

# facet_wrap variant of theme.adjusted (00_setup.R), only used here
theme.adjusted.facet <- theme(
  text = element_text(family = "Times New Roman"),
  axis.text.x = element_text(angle = 0, hjust = 0.5, margin = margin(t = 5), size = 12),
  axis.title.x = element_text(hjust = 0.39, margin = margin(t = 20), size = 12),
  axis.text.y = element_text(hjust = 1, margin = margin(r = 10), size = 12, angle = 0),
  axis.title.y = element_blank(),
  title = element_text(color = "black"),
  plot.background = element_rect(fill = "white", color = NA),
  # facet strip replaced by an in-panel top-left label (see plot_sens)
  strip.text = element_blank(),
  strip.background = element_blank()
)

# One figure per check group: 3 term-facets stacked (CP, LEZ, CP:LEZ), term
# name as a bold in-panel corner label instead of the default facet strip.
# No plot titles: the report supplies figure captions.
plot_sens <- function(tab, file) {
  tab$spec  <- droplevels(tab$spec)
  tab$stars <- sig_stars(tab$p.value)
  # must stay a factor with the same levels as tab$facet_label, otherwise
  # ggplot sorts the panels alphabetically
  facet_name_labels <- data.frame(
    facet_label = factor(levels(tab$facet_label),
                         levels = levels(tab$facet_label)))
  p <- ggplot(tab, aes(x = estimate, y = spec)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7, colour = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = sig_colour),
                   height = 0.25, linewidth = 0.6) +
    geom_point(aes(colour = sig_colour), size = 1.5) +
    # significance stars just above each significant point
    geom_text(aes(label = stars, colour = sig_colour), nudge_y = 0.3,
              size = 3.5, family = "Times New Roman", fontface = "bold") +
    scale_colour_identity() +
    scale_y_discrete(guide = guide_axis(check.overlap = FALSE)) +
    # in-panel label: same x start in every panel so left edges line up
    geom_text(data = facet_name_labels, aes(x = 0.6, y = Inf, label = facet_label),
             inherit.aes = FALSE, hjust = 0, vjust = 1.4,
             fontface = "bold", family = "Times New Roman", size = 3) +
    facet_wrap(~ facet_label, ncol = 1) +
    # same fixed scale for every facet/plot, so panels are comparable
    scale_x_continuous(limits = c(-1.5, 1), breaks = seq(-1.5, 1, by = 0.5)) +
    labs(x = "Effect on log transport CO2", y = NULL) +
    theme_minimal(base_size = 11) +
    theme.adjusted.facet

  # fixed height so all three figures share the same size/aspect ratio
  ggsave(file.path(out_dir, "figures", file), p,
         width = 10, height = 8, dpi = 300)
  p
}

plot_sens(sens_tab[sens_tab$spec_group == "Learner sensitivity", ],
         "fig_sensitivity_learner.png")
plot_sens(sens_tab[sens_tab$spec_group == "Control-set sensitivity", ],
         "fig_sensitivity_controls.png")
