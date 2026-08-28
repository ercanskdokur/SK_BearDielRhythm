# ==============================================================================
# 06b_loo_extended.R  —  Extended LOO: H0/H1/H2/H3/H4/H5
# ==============================================================================
# Compares the hypothesis models (H3 hour x human, H4 hour x temperature, H5
# plasticity) by LOO on the same 50K rows. Because all models inherit H0's
# subsample, elpd is COMPARABLE.
#
# CRITICAL NOTE: observation-level LOO is OPTIMISTIC because of hourly
# autocorrelation, so this table alone is NOT sufficient for model selection; it
# must be read together with the grouped leave-bears-out CV (09c_LOBO_extended.R).
#
# Produces (SI, one file per RESP -> intensity / duration):
#   Table_V3_04_LOO_Compare   -> Table S14 (intensity) / S15 (duration)
#   Table_V3_05_ParetoK       -> Table S16 (intensity) / S17 (duration)
#   Table_V3_06_StackingWeights -> Table S18 (intensity) / S19 (duration)
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("06b_loo_extended")
log_step("=== 06b_loo_extended: START ===")

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))

RESP <- RESP_KEY()
# The expected family depends on RESP: intensity=zi-Beta, duration=beta-binom.
# (H0-H5 within the same RESP share the family -> LOO is valid.)
EXPECTED_FAM <- if (identical(RESP, "duration")) "beta_binomial" else "zero_inflated_beta"
hyps <- c("H0","H1","H2","H3","H4","H5")
model_files <- stats::setNames(sprintf("%s_%s.rds", hyps, RESP), hyps)
log_step("RESP=%s -> models: %s", RESP, paste(model_files, collapse=", "))

fits <- list()
for (nm in names(model_files)) {
  fp <- file.path(mod_path, model_files[nm])
  if (!file.exists(fp)) { log_step("MISSING: %s", model_files[nm]); next }
  h <- readRDS(fp)
  if (!identical(h$family, EXPECTED_FAM)) { log_step("SKIP %s (family=%s, expected=%s)", nm, h$family, EXPECTED_FAM); next }
  if (!identical(h$n, NULL) && exists("N_REF")) {
    if (h$n != N_REF) log_step("WARNING: %s n=%d != reference %d (LOO comparison biased!)", nm, h$n, N_REF)
  }
  fits[[nm]] <- h$fit
  if (!exists("N_REF")) N_REF <- h$n
  log_step("Loaded %s (n=%d)", nm, h$n)
}
if (length(fits) < 2) stop("Not enough models")

# ---- n-consistency check (REQUIRED for LOO validity) -----------------------
ns <- sapply(fits, function(f) stats::nobs(f))
log_step("Model nobs: %s", paste(names(ns), ns, sep = "=", collapse = " | "))
if (length(unique(ns)) > 1) {
  log_step("!!! WARNING: models have different n -> elpd NOT COMPARABLE. Same subsample required.")
}

# ---- LOO + Pareto-k --------------------------------------------------------
loo_list <- list(); pk_rows <- list()
for (nm in names(fits)) {
  log_step("LOO: %s ...", nm)
  l <- tryCatch(loo::loo(fits[[nm]], cores = 1),
                error = function(e) { log_step("  ERR: %s", e$message); NULL })
  if (is.null(l)) next
  loo_list[[nm]] <- l
  pk <- l$diagnostics$pareto_k
  pk_rows[[nm]] <- data.frame(
    Model = nm, n_obs = length(pk),
    pk_le_07 = sum(pk <= 0.7), pk_07_1 = sum(pk > 0.7 & pk <= 1), pk_gt_1 = sum(pk > 1),
    elpd_loo = round(l$estimates["elpd_loo","Estimate"], 1),
    se_elpd  = round(l$estimates["elpd_loo","SE"], 1),
    p_loo    = round(l$estimates["p_loo","Estimate"], 1)
  )
}
pk_df <- do.call(rbind, pk_rows); rownames(pk_df) <- NULL
# -> Table S16 (intensity) / S17 (duration): Pareto-k reliability
write.csv(pk_df, file.path(tbl_path, paste0(tag_out("Table_V3_05_ParetoK"), ".csv")), row.names = FALSE)
cat("\n=== Pareto-k (extended) ===\n"); print(pk_df, row.names = FALSE)

# ---- loo_compare -----------------------------------------------------------
cmp <- loo::loo_compare(loo_list)
cmp_df <- as.data.frame(cmp); cmp_df$Model <- rownames(cmp_df)
cmp_df <- cmp_df[, c("Model","elpd_diff","se_diff","elpd_loo","se_elpd_loo","p_loo","se_p_loo")]
num <- sapply(cmp_df, is.numeric); cmp_df[num] <- lapply(cmp_df[num], round, 2)
rownames(cmp_df) <- NULL
# -> Table S14 (intensity) / S15 (duration): LOO model comparison (H0-H5)
write.csv(cmp_df, file.path(tbl_path, paste0(tag_out("Table_V3_04_LOO_Compare"), ".csv")), row.names = FALSE)
cat("\n=== loo_compare (extended) ===\n"); print(cmp_df, row.names = FALSE)

# ---- Stacking weights ------------------------------------------------------
w <- tryCatch(loo::loo_model_weights(loo_list, method = "stacking"),
              error = function(e) { log_step("stacking ERR: %s", e$message); NULL })
if (!is.null(w)) {
  w_df <- data.frame(Model = names(w), Stack_w = round(as.numeric(w), 3))
  w_df <- w_df[order(-w_df$Stack_w), ]
  # -> Table S18 (intensity) / S19 (duration): stacking weights
  write.csv(w_df, file.path(tbl_path, paste0(tag_out("Table_V3_06_StackingWeights"), ".csv")), row.names = FALSE)
  cat("\n=== Stacking weights (extended) ===\n"); print(w_df, row.names = FALSE)
}

save_rds_safe(loo_list, file.path(mod_path, "loo_objects_extended.rds"))
log_step("=== 06b_loo_extended DONE ===")
sink(type = "message"); sink()
