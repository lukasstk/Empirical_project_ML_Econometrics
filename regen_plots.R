# =============================================================================
# regen_plots.R
# Regenerates all coefficient plots from already-computed result tables,
# without re-running any DML fitting. Use this to check plot styling/margin
# changes in 00_setup.R / 04_sensitivity_analysis.R quickly.
#
# Usage: set out_dir to the run you want to (re)plot, then source this file.
#   out_dir <- "output_2026-07-08"
#   source("regen_plots.R")
# Figures are overwritten in <out_dir>/figures/.
# =============================================================================

if (!exists("out_dir")) out_dir <- "output_2026-07-08"
source("00_setup.R")

# Local copies of theme.adjusted and plot_effects() (both normally defined
# in 00_setup.R) - redefined here so you can tweak values and re-source
# just this file, without touching 00_setup.R. These override the versions
# loaded above for the rest of this script.
theme.adjusted <- theme(
  text = element_text(family = "Times New Roman"),
  axis.text.x = element_text(angle = 0, hjust = 0.5, margin = margin(t = 5), size = 15),
  axis.title.x = element_text(margin = margin(t = 20), size = 20),
  axis.text.y = element_text(hjust = 1, margin = margin(r = 10), size = 15, angle = 0),
  axis.title.y = element_text(margin = margin(r = 20), size = 20),
  title = element_text(color = "black"),
  plot.title = element_text(size = 18, color = "black", face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 17, color = "black", face = "italic"),
  panel.grid.major = element_line(color = "darkgray", linewidth = 0.2),
  panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
  plot.background = element_rect(fill = "white", color = NA)
)

# Local copy of pretty_term() (normally in 00_setup.R) - same reason as
# theme.adjusted/plot_effects() above: edit labels here, re-source this
# file only.
pretty_term <- function(term) {
  base <- c(
    cp_active  = "CP",
    lez_active = "LEZ",
    cp_x_lez   = "CP:LEZ",
    only_cp    = "CP only",
    only_lez   = "LEZ only",
    both       = "CP & LEZ (both)"
  )
  het <- c(
    log_population        = "population (log)",
    log_gdp_pc            = "GDP per capita (log)",
    log_pop_density       = "population density (log)",
    public_transit_score  = "public transit score",
    political_green       = "green party vote share",
    fuel_price             = "fuel price",
    tourism_intensity      = "tourism intensity",
    industry_logistics     = "logistics industry share",
    coastal                = "coastal city"
  )
  out <- unname(base[term])
  cp_hit  <- is.na(out) & grepl("^cp_x_",  term)
  lez_hit <- is.na(out) & grepl("^lez_x_", term)
  out[cp_hit]  <- paste0("CP:",  het[sub("^cp_x_",  "", term[cp_hit])])
  out[lez_hit] <- paste0("LEZ:", het[sub("^lez_x_", "", term[lez_hit])])
  ifelse(is.na(out), term, out)
}

plot_effects <- function(res, title, file = NULL, subtitle = NULL,
                         x_breaks = scales::pretty_breaks(n = 10)) {
  res$sig_colour <- classify_colour(res$estimate, res$conf.low, res$conf.high)
  res$term_label <- pretty_term(res$term)
  p <- ggplot(res, aes(x = estimate, y = reorder(term_label, estimate))) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7, colour = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = sig_colour),
                   height = 0.25, linewidth = 1) +
    geom_point(aes(colour = sig_colour), size = 3.5) +
    scale_colour_identity() +
    scale_y_discrete(guide = guide_axis(check.overlap = FALSE)) +
    scale_x_continuous(breaks = x_breaks) +
    labs(x = "Effect on log transport CO2", y = NULL,
         title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme.adjusted
  if (!is.null(file)) {
    ggsave(file.path(out_dir, "figures", file), p,
           width = 10, height = 4, dpi = 300)
  }
  p
}

tab_dir <- file.path(out_dir, "tables")

res_main   <- read.csv(file.path(tab_dir, "tab_main_effects_dml.csv"))
res_regime <- read.csv(file.path(tab_dir, "tab_policy_regimes.csv"))
het        <- read.csv(file.path(tab_dir, "tab_heterogeneity_dml.csv"))
sens_tab   <- read.csv(file.path(tab_dir, "tab_sensitivity_specs.csv"))

plot_effects(res_regime,
             "Policy regimes vs. no policy: congestion pricing (CP) and low-emission zone (LEZ)",
             "fig_policy_regimes.png",
             x_breaks = scales::breaks_width(0.25))
plot_effects(res_main,
             "Average policy effects: congestion pricing (CP) and low-emission zone (LEZ)",
             "fig_main_effects.png",
             x_breaks = scales::breaks_width(0.25))

res_cp  <- het[grepl("^cp_x_",  het$term) & het$term != "cp_x_lez", ]
res_lez <- het[grepl("^lez_x_", het$term), ]
plot_effects(res_cp,
             "Heterogeneity of the congestion pricing (CP) effect",
             "fig_het_cp.png")
plot_effects(res_lez,
             "Heterogeneity of the low-emission zone (LEZ) effect",
             "fig_het_lez.png")

# theme.adjusted.facet lives in 04_sensitivity_analysis.R now, not
# 00_setup.R - duplicated here since this script doesn't source 04.
theme.adjusted.facet <- theme(
  text = element_text(family = "Times New Roman"),
  axis.text.x = element_text(angle = 0, hjust = 0.5, margin = margin(t = 5), size = 12),
  axis.title.x = element_text(hjust = 0.39, margin = margin(t = 20), size = 15),
  axis.text.y = element_text(hjust = 1, margin = margin(r = 10), size = 12, angle = 0),
  axis.title.y = element_blank(),
  title = element_text(color = "black"),
  plot.title = element_text(size = 15, color = "black", face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 10, color = "black", face = "italic"),
  plot.background = element_rect(fill = "white", color = NA),
  strip.text = element_blank(),
  strip.background = element_blank()
)

# Cached tab_sensitivity_specs.csv predates both the label-shortening and
# the CV-lasso-variant removal in 04_sensitivity_analysis.R - match on the
# original long labels (not the current short text) or this patch silently
# matches nothing, and drop the rows that no longer exist in the current
# code (CV lasso 1se/min, cross-fitting fold checks) so this preview
# matches what the updated pipeline will actually produce.
sens_tab <- sens_tab[!sens_tab$spec %in% c(
  "CV lasso, lambda.1se", "CV lasso, lambda.min",
  "cross-fitting: 2 folds", "cross-fitting: 10 folds"), ]
sens_tab$spec[sens_tab$spec == "plugin lasso, hdm rlasso (baseline)"] <- "rlasso (baseline)"
sens_tab$spec[sens_tab$spec == "ridge, lambda.1se"] <- "CV ridge (1se)"
sens_tab$spec[sens_tab$spec == "minimal controls (main effects only)"] <- "minimal controls"
sens_tab$spec[sens_tab$spec == "+ mediators (pm25, fleet shares) [bad controls]"] <- "with added mediators"

# The COVID-exclusion row (block (3) in 04_sensitivity_analysis.R) doesn't
# exist in this cached CSV at all. Duplicate the baseline row under that
# spec name purely to preview the 3-group grid layout - these numbers are
# NOT real COVID-exclusion results, only a placeholder for layout testing.
covid_placeholder <- sens_tab[sens_tab$spec == "rlasso (baseline)", ]
covid_placeholder$spec <- "excl. COVID (2020-21)"
sens_tab <- rbind(sens_tab, covid_placeholder)

# spec has no explicit plot order otherwise, so ggplot falls back to
# alphabetical - pin it to the order specs first appear in the table
# instead, reversed so the y-axis (which draws its first level at the
# bottom) puts baseline at the top.
sens_tab$spec <- factor(sens_tab$spec, levels = rev(unique(sens_tab$spec)))

sens_tab$sig_colour <- classify_colour(sens_tab$estimate, sens_tab$conf.low, sens_tab$conf.high)
# Full names for the facet strips (not pretty_term()'s short "CP"/"LEZ" -
# each facet spells out its own term, so no shared subtitle is needed).
facet_titles <- c(
  cp_active  = "Congestion pricing (CP)",
  lez_active = "Low-emission zone (LEZ)",
  cp_x_lez   = "CP:LEZ"
)
# Explicit factor levels so panel order is always CP, LEZ, CP:LEZ
# top-to-bottom, consistently across all three sensitivity figures below.
sens_tab$facet_label <- factor(facet_titles[sens_tab$term], levels = facet_titles)

# Three different sensitivity checks: rlasso/ridge/random forest vary the
# nuisance LEARNER; minimal controls/with added mediators vary the
# CONTROL SET (W); excl. COVID varies the SAMPLE. Used below to split
# sens_tab into three separate plots, one per group.
spec_group_levels <- c("Learner sensitivity", "Control-set sensitivity",
                       "Sample robustness")
spec_groups <- c(
  "rlasso (baseline)"     = "Learner sensitivity",
  "CV ridge (1se)"        = "Learner sensitivity",
  "random forest"         = "Learner sensitivity",
  "minimal controls"      = "Control-set sensitivity",
  "with added mediators"  = "Control-set sensitivity",
  "excl. COVID (2020-21)" = "Sample robustness"
)
sens_tab$spec_group <- factor(spec_groups[as.character(sens_tab$spec)],
                              levels = spec_group_levels)

# Three separate figures - one per spec_group, so no group label/strip is
# needed. The 3 term-facets are always CP first, LEZ middle, the CP:LEZ
# interaction last (facet_titles' definition order above), stacked
# (ncol = 1) rather than side by side; the term name is a bold corner
# label per panel instead of the default facet strip.
plot_sens <- function(tab, label, file) {
  tab$spec <- droplevels(tab$spec)
  facet_name_labels <- data.frame(facet_label = levels(tab$facet_label))
  p <- ggplot(tab, aes(x = estimate, y = spec)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7, colour = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = sig_colour),
                   height = 0.25, linewidth = 0.6) +
    geom_point(aes(colour = sig_colour), size = 1.5) +
    scale_colour_identity() +
    scale_y_discrete(guide = guide_axis(check.overlap = FALSE)) +
    # x = 0.5, hjust = 0 (left-justified starting AT x = 0.5): every label
    # starts at the same x, so left edges line up in a column across the 3
    # stacked panels. This theme's text size = 3 (vs. 6 in the real
    # 04_sensitivity_analysis.R plot), so labels are narrower here and can
    # sit further right without overflowing past x = 1 - I can't render
    # this to confirm exactly, so nudge x left if it clips, or right if
    # there's visible slack.
    geom_text(data = facet_name_labels, aes(x = 0.6, y = Inf, label = facet_label),
             inherit.aes = FALSE, hjust = 0, vjust = 1.4,
             fontface = "bold", family = "Times New Roman", size = 3) +
    facet_wrap(~ facet_label, ncol = 1) +
    # Same fixed scale for every facet/plot (not free_x any more): the
    # x = 0 reference line then sits at the same, centered position
    # everywhere, and panels are directly comparable to each other.
    scale_x_continuous(limits = c(-1.5, 1), breaks = seq(-1.5, 1, by = 0.5)) +
    labs(x = "Effect on log transport CO2", y = NULL,
         title = paste0("Sensitivity analysis of coefficients: ", label)) +
    theme_minimal(base_size = 11) +
    theme.adjusted.facet

  # Fixed height for every sensitivity figure (not row-count-dependent any
  # more), so all three sit at the same size/aspect ratio.
  ggsave(file.path(out_dir, "figures", file), p,
         width = 10, height = 8, dpi = 300)
  p
}

plot_sens(sens_tab[sens_tab$spec_group == "Learner sensitivity", ],
         "Different learners", "fig_sensitivity_learner.png")
plot_sens(sens_tab[sens_tab$spec_group == "Control-set sensitivity", ],
         "Different Control sets (W)", "fig_sensitivity_controls.png")
plot_sens(sens_tab[sens_tab$spec_group == "Sample robustness", ],
         "Non-Covid Sample", "fig_sensitivity_sample.png")

cat("\nDone. Figures written to", file.path(out_dir, "figures"), "\n")
