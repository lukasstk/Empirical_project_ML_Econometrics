# =============================================================================
# 05_run_all.R
# Master script: runs the full analysis pipeline in order.
# Outputs land in output/tables (csv) and output/figures (png).
#
# The baseline nuisance learner is the plugin-lasso (hdm::rlasso, no inner
# cross-validation), so most of the remaining runtime comes from the
# CV-lasso / ridge / random-forest comparison rows in the sensitivity script.
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

cat("\nAll scripts finished. See output/tables and output/figures.\n")
