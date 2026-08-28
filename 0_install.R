# ==============================================================================
# 0_install.R  —  package verification (container environment)
# ==============================================================================
# Inside the Singularity container all packages are pre-installed in the
# Dockerfile. This script only verifies they are present and warns if any are
# missing. No conda environment or custom libPaths is required.
# ==============================================================================

cat("=== 0_install.R: package verification (container mode) ===\n")
cat(sprintf("R: %s\n", R.version.string))
cat(sprintf("libPaths[1]: %s\n", .libPaths()[1]))

critical <- c(
  "conflicted",
  "dplyr", "tidyr", "lubridate", "readxl", "data.table", "writexl",
  "tibble", "stringr", "purrr", "readr", "forcats",
  "terra", "sf",
  "ctmm",
  "ncdf4",
  "mixtools", "circular", "suncalc",
  "brms", "rstan", "tidybayes", "loo", "bayesplot", "posterior",
  "ggplot2", "patchwork", "scales", "ggridges",
  "parallel", "future", "future.apply"
)

optional <- c("sp", "raster", "gridExtra", "corrplot")

ok        <- character(0)
fail_crit <- character(0)
fail_opt  <- character(0)

for (p in critical) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("[OK] %s\n", p))
    ok <- c(ok, p)
  } else {
    cat(sprintf("[MISSING] %s\n", p))
    fail_crit <- c(fail_crit, p)
  }
}

for (p in optional) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("[OK] %s (optional)\n", p))
    ok <- c(ok, p)
  } else {
    cat(sprintf("[WARN] %s missing (optional — pipeline still runs)\n", p))
    fail_opt <- c(fail_opt, p)
  }
}

cat("\n=== SUMMARY ===\n")
cat(sprintf("OK (critical) : %d / %d\n", sum(critical %in% ok), length(critical)))
cat(sprintf("OK (optional) : %d / %d\n", sum(optional %in% ok), length(optional)))

if (length(fail_opt) > 0) {
  cat(sprintf("[WARN] Missing optional: %s\n", paste(fail_opt, collapse = ", ")))
}

if (length(fail_crit) > 0) {
  cat(sprintf("\n[ERROR] Missing required packages: %s\n", paste(fail_crit, collapse = ", ")))
  cat("Rebuild the Dockerfile and update the container image.\n")
  quit(status = 1)
}

cat("\nAll required packages present. 0_install.R SUCCESSFUL.\n")
