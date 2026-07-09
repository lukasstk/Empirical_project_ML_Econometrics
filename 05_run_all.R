# =============================================================================
# 05_run_all.R
# Master script: runs the full analysis pipeline in order.
# Outputs land in a date-stamped folder, <out_dir>/tables (csv) and
# <out_dir>/figures (png), e.g. output_2026-07-08/ - so a rerun never
# overwrites earlier results (the original run lives in output/). Set
# out_dir before running to choose the folder name yourself (00_setup.R).
#
# The whole pipeline runs on the DoubleML package (see 00_setup.R): "One-By-
# One Double LASSO", each treatment coefficient estimated in turn with the
# other treatments folded into that run's nuisance set (DoubleML's default
# use_other_treat_as_covariate = TRUE). The baseline nuisance learner is the
# plugin-penalty lasso (hdm::rlasso, no inner CV -> one fit per nuisance);
# ridge and a random forest appear as learner-choice sensitivity rows in 04.
# Cross-fitting folds run on 5 parallel workers (future::plan in
# 00_setup.R). Expect a few minutes for 02/03 and most of the total time in
# 04 (random forest row).
# =============================================================================

t0 <- Sys.time()
run <- function(script) {
  cat("\n============================================================\n")
  cat(">>", script, "\n")
  cat("============================================================\n")
  source(script, echo = FALSE)
  cat(">> done after", round(difftime(Sys.time(), t0, units = "mins"), 1),
      "minutes total\n")
}

run("01_data_preparation.R")
run("02_main_effects_dml.R")
run("03_heterogeneity_dml.R")
run("04_sensitivity_analysis.R")

cat("\nAll scripts finished. See", file.path(out_dir, "tables"), "and",
    file.path(out_dir, "figures"), "\n")
