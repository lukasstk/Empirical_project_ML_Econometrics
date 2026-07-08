# =============================================================================
# 05_run_all.R
# Master script: runs the full analysis pipeline in order.
# Outputs land in output/tables (csv) and output/figures (png).
#
# The whole pipeline runs on the DoubleML package (see 00_setup.R). The
# baseline nuisance learner is the CV-lasso, so every nuisance fit runs an
# inner 10-fold cross-validation; expect a total runtime in the tens of
# minutes (02 uses n_rep = 5 repetitions, 03 has 21 treatment columns).
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
