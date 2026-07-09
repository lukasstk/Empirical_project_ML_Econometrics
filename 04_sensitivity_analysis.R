# =============================================================================
# 04_sensitivity_analysis.R
# Robustness of the main results to the modeling choices:
#
#   (1) ML learner for the nuisance functions:
#       plugin lasso, hdm::rlasso (baseline, as in 02/03) / ridge / random
#       forest - each asks whether the result depends on the sparsity
#       assumption or the linear-after-selection form; CV-tuned lasso
#       (lambda.1se/min) was left out since that's just a lambda-selection
#       comparison, not specific to this DML setup
#   (2) Control set: baseline / minimal (no interactions) / + mediators
#   (3) Sample: excluding the COVID years (2020-2021), an emissions shock
#       unrelated to CP/LEZ policy that year dummies may not fully absorb
#   (4) Placebo treatment: policy announced (major topic of public
#       discussion) but not yet implemented - tests identification itself
#
# All checks use the same DoubleML estimator as 02/03 (dml_plr, 00_setup.R).
# All rows use n_rep = 5 repeated cross-fitting splits (the dml_plr default),
# so the across-split dispersion is included in every reported CI; each row
# changes exactly one thing relative to the baseline.
# =============================================================================

source("00_setup.R")
prep <- readRDS(file.path(out_dir, "prepared_data.rds"))
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
# The baseline learner is the plugin-penalty lasso (hdm::rlasso): lambda
# from the theory-based formula of Belloni/Chernozhukov/Hansen, as used by
# rlassoEffects() in the lecture - one fit per nuisance, no inner CV.
#
# Deliberately NOT included: CV-tuned lasso (lambda.1se / lambda.min). That
# comparison is about which cross-validation rule picks lambda - a generic
# concern for any lasso application, not something specific to this DML
# setup or this data. Ridge and random forest are kept instead, since they
# ask a substantive question each: does the result depend on the sparsity
# assumption specifically (ridge), or on the linear-after-selection
# functional form (random forest) - not just "which tuning rule."
#
# Spec labels are kept short on purpose (full detail is in the comments
# above/below) - they double as facet y-axis text in the sensitivity
# figures below, and long labels there force cramped panels and
# truncated/overlapping text.
sens$rlasso    <- run_spec("rlasso (baseline)", Y, D, W_base,
                           learner = "rlasso")
# Ridge shrinks every coefficient but never selects a sparse support. DML's
# general theory (Chernozhukov et al. 2018) only needs the product of both
# nuisances' convergence rates to vanish faster than 1/sqrt(n) - it is not
# itself tied to sparsity - but the rate argument covered in the lecture
# (Belloni/Chernozhukov/Hansen) establishes this specifically via the
# sparsity assumption. Ridge would need its own separate rate condition
# (e.g. bounded effective dimension) that isn't established here, so it's
# included as a check on the point estimate only, not a learner with the
# same covered inference guarantee as the baseline.
sens$ridge     <- run_spec("CV ridge (1se)",     Y, D, W_base,
                           learner = "ridge")
# Same caveat for the forest: no explicit sparse support either, and its
# own convergence-rate theory (e.g. Wager/Athey) is a separate literature
# not covered in the course. Kept for the same reason - a check on the
# point estimate under a learner whose rate isn't established by the
# sparsity argument.
sens$rf        <- run_spec("random forest",      Y, D, W_main,
                           learner = "rf")

# ---- (2) Control sets --------------------------------------------------------
# Same learner as the baseline, so each row changes exactly one thing.
sens$min   <- run_spec("minimal controls", Y, D, W_main, learner = "rlasso")
# Note: W_ext extends the BASELINE control set (with interactions) plus
# mediators, not the minimal one above - "with added mediators" deliberately
# doesn't say "minimal", since it builds on the baseline spec, not on
# minimal controls.
sens$ext   <- run_spec("with added mediators",
                       Y, D, W_ext,  learner = "rlasso")

# ---- (3) Sample: exclude COVID years (2020-2021) ------------------------------
# Lockdowns and reduced traffic caused a large emissions drop unrelated to
# CP/LEZ policy in 2020-2021. Year dummies in W absorb the common shock,
# but if COVID hit treated and untreated cities differently, the estimate
# could still be biased. Y, D, W are rebuilt from the filtered data (not
# just subsetted) so W is re-centered on the reduced sample, same as the
# construction at the top of this script.
data_no_covid <- data[!(data$year %in% c(2020, 2021)), ]
Y_nc <- data_no_covid$log_transport_co2
D_nc <- cbind(cp_active  = data_no_covid$cp_active,
             lez_active = data_no_covid$lez_active,
             cp_x_lez   = data_no_covid$cp_x_lez)
W_nc <- build_W(data_no_covid, prep$ctrl_baseline, prep$sq_vars)
no_covid <- dml_plr(Y_nc, D_nc, W_nc, cluster = data_no_covid$city_id,
                    learner = "rlasso", seed = 42)
sens$no_covid <- cbind(no_covid$results, spec = "excl. COVID (2020-21)")

# ---- (4) Placebo treatment: announced but not yet active ---------------------
# Blocks (1)-(3) vary the estimator or the sample. This check targets the
# identifying assumption (conditional ignorability) instead, which no
# estimator or sample choice can fix. Between the year a policy first became
# a major topic of local public discussion and its implementation year
# (~1.5 years on average here), it cannot mechanically affect emissions -
# no charge is collected, no vehicle is banned - so a placebo indicator for
# exactly these city-years should have a coefficient close to zero.
#
# A clearly negative placebo coefficient means emissions fall before the
# policy exists, with two competing readings:
#   - anticipation: households and firms adapt early (e.g. replacing a
#     diesel car before the LEZ starts). Still policy-caused, but the
#     pre-implementation years are then partly treated, so the main
#     estimates understate the total policy effect.
#   - selection / awareness: the public debate reflects (and fuels) general
#     environmental awareness - people drive and emit less regardless of
#     the future policy, and greener local governments tend to push other
#     measures (transit, cycling, parking) at the same time. Not caused by
#     the policy instruments, so the main estimates overstate their effect.
# The data can't fully separate the two (every city that discussed a policy
# eventually implemented it), but a near-zero placebo rules both out, and a
# negative one bounds the main estimates from one side or the other.
#
# The placebo indicators enter jointly with the actual treatment indicators,
# so their coefficients are net of the active-policy effects.
plc <- function(announce, impl) {
  as.integer(announce > 0 & data$year >= announce &
               (impl == 0 | data$year < impl))
}
D_plac <- cbind(D,
                cp_pre  = plc(data$cp_announce_year,  data$cp_impl_year),
                lez_pre = plc(data$lez_announce_year, data$lez_impl_year))
placebo <- dml_plr(Y, D_plac, W_base, cluster = data$city_id,
                   learner = "rlasso", seed = 42)
save_table(placebo$results, "tab_placebo_announcement")

sens_tab <- do.call(rbind, sens)
rownames(sens_tab) <- NULL
save_table(sens_tab, "tab_sensitivity_specs")

# spec has no explicit plot order otherwise, so ggplot falls back to
# alphabetical - pin it to the definition order above instead (baseline
# first, "with added mediators" last), reversed so the y-axis (which draws
# its first level at the bottom) puts baseline at the top.
sens_tab$spec <- factor(sens_tab$spec, levels = rev(unique(sens_tab$spec)))

sens_tab$sig_colour <- classify_colour(sens_tab$estimate, sens_tab$conf.low, sens_tab$conf.high)
# Full names for the facet strips (not pretty_term()'s short "CP"/"LEZ" -
# each facet spells out its own term, so no shared subtitle is needed).
facet_titles <- c(
  cp_active  = "Congestion pricing (CP)",
  lez_active = "Low-emission zone (LEZ)",
  cp_x_lez   = "CP:LEZ interaction"
)
# Explicit factor levels (not alphabetical default) so panel order is
# always CP, LEZ, CP:LEZ top-to-bottom, consistently across all three
# sensitivity figures below.
sens_tab$facet_label <- factor(facet_titles[sens_tab$term], levels = facet_titles)

# Three different sensitivity checks (see comments above run_spec() blocks
# (1), (2), (3)): rlasso/ridge/random forest vary the nuisance LEARNER;
# minimal controls/with added mediators vary the CONTROL SET (W); excl.
# COVID varies the SAMPLE. Used below to split sens_tab into three
# separate plots, one per group.
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

# facet_wrap variant of theme.adjusted (00_setup.R) - only used here, so
# it's local to this script rather than in the shared setup file.
theme.adjusted.facet <- theme(
  text = element_text(family = "Times New Roman"),
  axis.text.x = element_text(angle = 0, hjust = 0.5, margin = margin(t = 5), size = 22.5),
  axis.title.x = element_text(hjust = 0.39, margin = margin(t = 20), size = 25),
  axis.text.y = element_text(hjust = 1, margin = margin(r = 10), size = 20, angle = 0),
  axis.title.y = element_blank(),
  title = element_text(color = "black"),
  plot.title = element_text(size = 25, color = "black", face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 17, color = "black", face = "italic"),
  plot.background = element_rect(fill = "white", color = NA),
  # Default facet strip is replaced by an in-panel top-left label (see
  # plot_sens() below) - a full-width strip above each of the 3 stacked
  # panels looked odd, so it's hidden here.
  strip.text = element_blank(),
  strip.background = element_blank()
)

# Three separate figures - learner-choice, control-set, and sample (excl.
# COVID) checks test genuinely different things (see the run_spec() block
# comments above), so each gets its own full title and its own file. No
# group label/strip needed any more since each figure now covers exactly
# one spec_group. The 3 term-facets are always CP first, LEZ middle, the
# CP:LEZ interaction last (facet_titles' definition order above), stacked
# (ncol = 1) rather than side by side; the term name is a bold corner
# label per panel instead of the default facet strip.
plot_sens <- function(tab, label, file) {
  tab$spec <- droplevels(tab$spec)
  facet_name_labels <- data.frame(facet_label = levels(tab$facet_label))
  p <- ggplot(tab, aes(x = estimate, y = spec)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7, colour = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = sig_colour),
                   height = 0.25, linewidth = 1) +
    geom_point(aes(colour = sig_colour), size = 3.5) +
    scale_colour_identity() +
    scale_y_discrete(guide = guide_axis(check.overlap = FALSE)) +
    # x = 0.25, hjust = 0 (left-justified starting AT x = 0.25): every
    # label starts at the same x, so left edges line up in a column across
    # the 3 stacked panels. Pushed as far right as a rough estimate allows
    # without the longest label ("Low-emission zone (LEZ)", ~23 chars at
    # this text size/panel width) overflowing past the x = 1 panel edge -
    # I can't render this to confirm exactly, so nudge x left if it clips,
    # or right if there's visible slack.
    geom_text(data = facet_name_labels, aes(x = 0.25, y = Inf, label = facet_label),
             inherit.aes = FALSE, hjust = 0, vjust = 1.4,
             fontface = "bold", family = "Times New Roman", size = 6) +
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
         "Control set (W)", "fig_sensitivity_controls.png")
plot_sens(sens_tab[sens_tab$spec_group == "Sample robustness", ],
         "Sample", "fig_sensitivity_sample.png")
