# =============================================================================
# 05a_run_all_sequential.R
# Master script: runs the full analysis pipeline sequentially, in order.
# Outputs land in output/tables (csv) and output/figures (png); set
# out_dir before running to choose a different folder name (00_setup.R).
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

run("Code/00_setup.R")
run("Code/01_data_preparation.R")
run("Code/02_main_and_joint_effects.R")
run("Code/03a_heterogeneity_primary.R")  # Code/03b_heterogeneity_all_controls.R
                                         # (robustness appendix) is run manually
run("Code/04_sensitivity_analysis.R")

cat("\nAll scripts finished. See", file.path(out_dir, "tables"), "and",
    file.path(out_dir, "figures"), "\n")
