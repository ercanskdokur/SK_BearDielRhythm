# ==============================================================================
# 09e_lobo_combine.R  —  Combine per-fold elpd -> grouped-CV ranking
# ==============================================================================
# Reads the models/lobo/<MODEL>_<RESP>_fold*.rds files written by 09d_lobo_fold.R.
# If ALL K folds are present for a model: concatenate the pointwise elpd (each
# observation is held out EXACTLY once -> the total spans all N observations),
# elpd_kfold=sum, se=sqrt(N)*sd(pointwise). Compare with obs-LOO -> OPTIMISM.
#
# OUTPUT: tables/Table_V3_07_LOBO_<RESP>.csv
#         -> Table S20 (intensity) / S22 (duration)
# ENV: RESP=intensity|duration | K_FOLDS=10
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log(sprintf("09e_lobo_combine_%s", RESP_KEY()))
RESP <- RESP_KEY()
K    <- as.integer(Sys.getenv("K_FOLDS", unset = "10"))
lobo_dir <- file.path(mod_path, "lobo")
log_step("=== 09e combine | RESP=%s K=%d ===", RESP, K)

hyps <- c("H0","H1","H2","H3","H4","H5")
rows <- list()
for (M in hyps) {
  ffs <- file.path(lobo_dir, sprintf("%s_%s_fold%d.rds", M, RESP, seq_len(K)))
  have <- file.exists(ffs)
  if (!all(have)) { log_step("%s: %d/%d folds present -> SKIP (incomplete)", M, sum(have), K); next }
  fl <- lapply(ffs, readRDS)
  pw <- unlist(lapply(fl, `[[`, "elpd_pointwise"))
  N  <- length(pw)
  ndiv <- sum(vapply(fl, function(x) x$divergent, numeric(1)))
  mrh  <- max(vapply(fl, function(x) x$max_rhat, numeric(1)), na.rm = TRUE)
  elpd_kfold <- sum(pw)
  se_kfold   <- sqrt(N) * stats::sd(pw)

  # obs-LOO (same fit)
  h <- readRDS(file.path(mod_path, sprintf("%s_%s.rds", M, RESP)))
  lo <- tryCatch(loo::loo(h$fit, cores = 1), error = function(e) NULL)
  elpd_obs <- if (is.null(lo)) NA_real_ else lo$estimates["elpd_loo","Estimate"]
  se_obs   <- if (is.null(lo)) NA_real_ else lo$estimates["elpd_loo","SE"]

  rows[[M]] <- data.frame(
    Model = M, n_obs = N, n_folds = K,
    elpd_obsLOO  = round(elpd_obs, 1), se_obsLOO = round(se_obs, 1),
    elpd_groupCV = round(elpd_kfold, 1), se_groupCV = round(se_kfold, 1),
    optimism = round(elpd_obs - elpd_kfold, 1),
    fold_divergent = ndiv, fold_max_rhat = round(mrh, 4)
  )
  log_step("%s: obs-LOO=%.1f  grouped-CV=%.1f  optimism=%.1f (div=%d, Rhat=%.3f)",
           M, elpd_obs, elpd_kfold, elpd_obs - elpd_kfold, ndiv, mrh)
}

if (!length(rows)) { log_step("No model has all K folds — combine skipped"); quit(save="no") }
res <- do.call(rbind, rows); rownames(res) <- NULL
res$rank_obsLOO  <- rank(-res$elpd_obsLOO)
res$rank_groupCV <- rank(-res$elpd_groupCV)
# difference from the grouped-CV best
best <- which.max(res$elpd_groupCV)
res$dCV_vs_best <- round(res$elpd_groupCV - res$elpd_groupCV[best], 1)
res <- res[order(res$rank_groupCV), ]

out_csv <- file.path(tbl_path, paste0(tag_out("Table_V3_07_LOBO"), ".csv"))
# -> Table S20 (intensity) / S22 (duration): grouped leave-bears-out CV
write.csv(res, out_csv, row.names = FALSE)
cat("\n=== RANKING: obs-LOO vs grouped bear-CV (", RESP, ") ===\n")
print(res[, c("Model","elpd_obsLOO","rank_obsLOO","elpd_groupCV","se_groupCV",
              "rank_groupCV","dCV_vs_best","optimism")], row.names = FALSE)
log_step("=== 09e combine DONE -> %s ===", basename(out_csv))
sink(type = "message"); try(sink(), silent = TRUE)
