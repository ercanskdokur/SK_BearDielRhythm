# ==============================================================================
# 09c_LOBO_extended.R  —  Grouped bear-CV: H1/H2/H3/H4 (for model selection)
# ==============================================================================
# 09b_LOBO.R did grouped-CV for H0 only. But MODEL SELECTION rests on obs-LOO
# (stacking H2=0.909) and obs-LOO is OPTIMISTIC (H0 inflated by +3584 elpd).
# H2's advantage over H1 in obs-LOO is only ~760 -> below the optimism gap.
#
# This script runs grouped 10-fold bear-CV for H1, H2, H3, H4 as well; the model
# ranking is re-evaluated on GROUPED-CV elpd. The ranking may DIFFER from obs-LOO
#
# TIME: each model ~100-112h (same budget as 09b). A SLURM array (0-3) is
# recommended for parallel runs, or one at a time via the MODEL env var.
#
# ENV: MODEL = "H1" | "H2" | "H3" | "H4"  (default: all in sequence — long!)
#
# OUTPUTS: Table_V3_07_LOBO_extended.csv (appended each run)
#          output/models/kfold_<MODEL>.rds
#          -> Table S20 (intensity) / S22 (duration)  [canonical combiner is 09e]
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("09c_LOBO_extended")
log_step("=== 09c_LOBO_extended: START ===")

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE)
options(brms.backend = "rstan")
options(future.globals.maxSize = 8 * 1024^3)
stan_cache <- file.path(tmp_dir, "stan_cache")
dir.create(stan_cache, recursive = TRUE, showWarnings = FALSE)

K_FOLDS <- as.integer(Sys.getenv("K_FOLDS", unset = "10"))
which_models <- Sys.getenv("MODEL", unset = "")
todo <- if (nchar(which_models) > 0) which_models else c("H1","H2","H3","H4")
log_step("LOBO models to run: %s", paste(todo, collapse = ", "))

out_csv <- file.path(tbl_path, paste0(tag_out("Table_V3_07_LOBO"), ".csv"))

run_lobo <- function(tag) {
  fp <- file.path(mod_path, sprintf("%s_%s.rds", tag, RESP_KEY()))
  if (!file.exists(fp)) { log_step("MISSING: %s", fp); return(NULL) }
  h <- readRDS(fp); fit <- h$fit
  n_bears <- dplyr::n_distinct(h$data$BearID_f)
  log_step("%s: grouped %d-fold CV, %d bears, n=%d", tag, K_FOLDS, n_bears, h$n)

  t0 <- Sys.time()
  kf <- tryCatch(
    brms::kfold(fit, folds = "grouped", group = "BearID_f", K = K_FOLDS,
                chains = 2, cores = 2, iter = 1500, warmup = 1000,
                save_fits = FALSE, recompile = FALSE),
    error = function(e) { log_step("%s kfold FAIL: %s", tag, e$message); NULL })
  log_step("%s kfold time: %.1f min", tag,
           as.numeric(difftime(Sys.time(), t0, units = "mins")))
  if (is.null(kf)) return(NULL)

  loo_std <- loo::loo(fit, cores = 1)
  row <- data.frame(
    Model = tag,
    elpd_obsLOO   = round(loo_std$estimates["elpd_loo","Estimate"], 1),
    se_obsLOO     = round(loo_std$estimates["elpd_loo","SE"], 1),
    elpd_groupCV  = round(kf$estimates["elpd_kfold","Estimate"], 1),
    se_groupCV    = round(kf$estimates["elpd_kfold","SE"], 1),
    optimism      = round(loo_std$estimates["elpd_loo","Estimate"] -
                         kf$estimates["elpd_kfold","Estimate"], 1)
  )
  save_rds_safe(kf, file.path(mod_path, sprintf("kfold_%s.rds", tag)))
  # A SLURM array (0-3) runs concurrently -> appending to a shared file creates a
  # race/corruption. Instead EACH MODEL writes its own part file (atomic, race-free);
  # the merge below reads the part files.
  write.csv(row, file.path(tbl_path, sprintf("Table_V3_07_LOBO_%s_%s.csv", tag, RESP_KEY())),
            row.names = FALSE)
  cat("\n=== ", tag, " LOBO ===\n"); print(row, row.names = FALSE)
  row
}

for (tag in todo) run_lobo(tag)

# ---- Summary: grouped-CV ranking vs obs-LOO ranking ------------------------
# Merge the part files (each model wrote its own). The last array task to finish
# writes the full table; earlier tasks that see an incomplete set still write with
# a FULL OVERWRITE (not append) -> no corruption.
parts <- list.files(tbl_path, pattern = "^Table_V3_07_LOBO_H[0-9]+\\.csv$",
                    full.names = TRUE)
if (length(parts) > 0) {
  res <- do.call(rbind, lapply(parts, utils::read.csv))
  res <- res[!duplicated(res$Model, fromLast = TRUE), ]
  res$rank_obsLOO  <- rank(-res$elpd_obsLOO)
  res$rank_groupCV <- rank(-res$elpd_groupCV)
  res <- res[order(res$rank_groupCV), ]
  write.csv(res, out_csv, row.names = FALSE)
  cat("\n=== RANKING: obs-LOO vs grouped-CV ===\n")
  print(res[, c("Model","elpd_obsLOO","rank_obsLOO","elpd_groupCV","rank_groupCV","optimism")],
        row.names = FALSE)
  log_step("Combined table written (%d models). If rankings differ: model selection should follow grouped-CV.",
           nrow(res))
}

log_step("=== 09c_LOBO_extended DONE ===")
sink(type = "message"); sink()
