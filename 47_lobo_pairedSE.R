# ==============================================================================
# 47_lobo_pairedSE.R  —  Grouped bear-CV: PAIRED elpd difference + SE  [B2]
# ==============================================================================
# Call:  Rscript 47_lobo_pairedSE.R
#
# WHY (audit B2):
#   Table_V3_07 gave dCV_vs_best WITHOUT an SE; the per-model SE is ~355 but because
#   the folds are SHARED, the PAIRED SE of the difference is far smaller and that is
#   what matters. The saved fold objects carry per-observation elpd (elpd_pointwise)
#   -> all models use the same foldmap/observations -> the paired difference can be
#   computed:
#     elpd_diff = sum(elpd_best_i - elpd_M_i),  se_diff = sd(diff)*sqrt(N).
#
# OUTPUT: tables/Table_V3_07_LOBO_pairedSE_<RESP>.csv  (per response)
#           -> Table S21 (intensity) / S23 (duration)
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("47_lobo_pairedSE")
log_step("=== 47_lobo_pairedSE: START ===")

lobo_dir <- file.path(mod_path, "lobo")
models <- c("H0","H1","H2","H3","H4","H5")

for (RESP in c("intensity","duration")) {
  log_step("--- RESP=%s ---", RESP)
  concat <- list(); bref <- NULL; ok <- TRUE
  for (M in models) {
    fps <- file.path(lobo_dir, sprintf("%s_%s_fold%d.rds", M, RESP, 1:10))
    if (!all(file.exists(fps))) { log_step("!! %s: missing fold, skipping", M); ok <- FALSE; break }
    parts <- lapply(fps, function(fp) { x <- readRDS(fp)
      list(elpd = x$elpd_pointwise, bear = x$bearid) })
    elpd <- unlist(lapply(parts, `[[`, "elpd"))
    bear <- unlist(lapply(parts, `[[`, "bear"))
    concat[[M]] <- elpd
    if (is.null(bref)) bref <- bear
    else if (!identical(bear, bref)) log_step("!! %s: bearid alignment DIFFERS (caution)", M)
  }
  if (!ok) next
  lens <- sapply(concat, length)
  if (length(unique(lens)) != 1) { log_step("!! lengths differ: %s", paste(lens, collapse=",")); next }
  N <- lens[1]
  P <- do.call(cbind, concat)                      # observation x model
  elpd_tot <- colSums(P)
  best <- names(which.max(elpd_tot))
  log_step("Best grouped-CV model: %s (elpd=%.1f)", best, elpd_tot[best])
  res <- lapply(models, function(M) {
    diff <- P[, best] - P[, M]
    data.frame(Model = M, elpd_groupCV = round(elpd_tot[M], 1),
               elpd_diff_vs_best = round(sum(diff), 1),
               se_diff = round(stats::sd(diff) * sqrt(N), 1),
               n_obs = N, best = best)
  })
  out <- do.call(rbind, res)
  out <- out[order(-out$elpd_groupCV), ]
  out$z_diff <- ifelse(out$se_diff > 0, round(out$elpd_diff_vs_best / out$se_diff, 2), NA)
  # -> Table S21 (intensity) / S23 (duration): LOBO paired-SE model comparison
  utils::write.csv(out, file.path(tbl_path, sprintf("Table_V3_07_LOBO_pairedSE_%s.csv", RESP)),
                   row.names = FALSE)
  cat(sprintf("\n=== LOBO paired-SE (%s) ===\n", RESP)); print(out, row.names = FALSE)
}
log_step("=== 47_lobo_pairedSE DONE ===")
sink(type = "message"); sink()
