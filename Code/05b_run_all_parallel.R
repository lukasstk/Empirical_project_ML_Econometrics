# =============================================================================
# 05b_run_all_parallel.R
# Script-level parallel version of 05a_run_all_sequential.R: runs 00+01
# sequentially,
# then launches 02/03/04 (which are independent, each reading only
# prepared_data.rds) as background Rscript processes on the same out_dir.
# Console output of the children goes to .run_logs/ in the project root
# (gitignored; kept out of out_dir so the output folder stays clean),
# completion is signalled via sentinel files. Full speed needs ~15 free
# cores (3 children x 5 future workers each).
# =============================================================================

t0 <- Sys.time()

# ---- (1) Setup + data preparation, sequentially -------------------------------
source("Code/00_setup.R")      # fixes out_dir for this run
source("Code/01_data_preparation.R")

n_cores <- parallel::detectCores()
cat("\nDetected", n_cores, "cores;",
    "3 children x 5 future workers each wants ~15.\n")
if (n_cores < 15) cat("Fewer than 15 cores: children will time-slice",
                      "(still correct, just less than the full 3x speedup).\n")

# ---- (2) Launch 02/03/04 as background processes ------------------------------
# Code/03b_heterogeneity_all_controls.R (robustness appendix) is not in this
# list; run it manually on the same out_dir when needed.
scripts <- c("Code/02_main_and_joint_effects.R",
             "Code/03a_heterogeneity_primary.R",
             "Code/04_sensitivity_analysis.R")
# Logs live outside out_dir so the output folder holds only results.
# Deliberately a RELATIVE path: absolute Windows temp paths contain
# backslashes, which turn into (invalid) escape sequences when embedded in
# the child R command below - the children then can never write their
# sentinels and the wait loop hangs at 0/3 forever.
log_dir <- ".run_logs"
dir.create(log_dir, showWarnings = FALSE)
cat("Child logs go to:", log_dir, "\n")

# basename(): scripts live in Code/, but sentinel/log files are flat names
sentinel <- function(s) file.path(log_dir, paste0(".done_", basename(s)))
for (s in scripts) unlink(sentinel(s))   # stale sentinels from earlier runs

for (s in scripts) {
  # Each child gets out_dir injected before 00_setup.R runs, sources its
  # script, and writes OK/FAIL to its sentinel.
  expr <- sprintf(
    "out_dir <- '%s'; ok <- tryCatch({ source('%s', echo = FALSE); TRUE },
       error = function(e) { message(conditionMessage(e)); FALSE });
     writeLines(if (ok) 'OK' else 'FAIL', '%s')",
    out_dir, s, sentinel(s))
  system2("Rscript",
          args   = c("-e", shQuote(expr)),
          stdout = file.path(log_dir, paste0(basename(s), ".log")),
          stderr = file.path(log_dir, paste0(basename(s), ".log")),
          wait   = FALSE)
  cat(">> launched", s, "in the background\n")
}

# ---- (3) Wait for all three, then report --------------------------------------
repeat {
  done <- vapply(scripts, function(s) file.exists(sentinel(s)), logical(1))
  cat(sprintf("\r%d/%d scripts finished (%.1f min elapsed) ",
              sum(done), length(scripts),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (all(done)) break
  Sys.sleep(10)
}
cat("\n\n")

for (s in scripts) {
  status <- readLines(sentinel(s), warn = FALSE)[1]
  cat(sprintf("%-28s %s\n", s, status))
  if (!identical(status, "OK"))
    cat("   -> see", file.path(log_dir, paste0(basename(s), ".log")), "\n")
}

cat("\nAll scripts finished after",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "minutes. See",
    file.path(out_dir, "tables"), "and", file.path(out_dir, "figures"), "\n")
