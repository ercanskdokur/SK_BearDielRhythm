# ==============================================================================
# 06c_dumpspec_compare.R  —  RAW vs RESID dump-parametrization comparison
# ==============================================================================
# Main deliverable of the "24 fits": dump distance (d2GarbageDump) is collinear
# with LandPC1 (VIF ~12.3). To address this, every hypothesis was fitted twice:
#   raw   : H{h}_{resp}.rds            (d2GarbageDump original)
#   resid : H{h}_{resp}_dumpresid.rds  (residual orthogonal to LandPC1+LandPC2)
#
# raw and resid use the SAME response + SAME rows (resid only re-expresses the
# dump column) => elpd is DIRECTLY comparable; loo_compare is valid.
#
# For each hypothesis this script:
#   (1) LOO(raw), LOO(resid) -> elpd_diff (resid - raw), se_diff
#   (2) compares the dump coefficient (fixef d2GarbageDump_km_sc) between the two
#       models -> was collinearity inflating the coefficient?
#
# Produces (SI, one file per RESP -> intensity / duration):
#   Table_V3_07_DumpSpec_LOO  -> Table S24 (intensity) / S25 (duration)
#   Table_V3_08_DumpSpec_Coef -> Table S26 (intensity) / S27 (duration)
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("06c_dumpspec_compare")
log_step("=== 06c_dumpspec_compare: START ===")

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))

RESP <- RESP_KEY()
hyps <- c("H0","H1","H2","H3","H4","H5")
dump_par <- "d2GarbageDump_km_sc"   # fixef row name
log_step("RESP=%s", RESP)

loo_rows  <- list()
coef_rows <- list()

for (h in hyps) {
  fp_raw <- file.path(mod_path, sprintf("%s_%s.rds", h, RESP))
  fp_res <- file.path(mod_path, sprintf("%s_%s_dumpresid.rds", h, RESP))
  if (!file.exists(fp_raw) || !file.exists(fp_res)) {
    log_step("SKIP %s: raw=%s resid=%s", h, file.exists(fp_raw), file.exists(fp_res))
    next
  }
  mr <- readRDS(fp_raw); ms <- readRDS(fp_res)
  fr <- mr$fit; fs <- ms$fit

  # --- n-consistency (REQUIRED for elpd comparison) -------------------------
  nr <- stats::nobs(fr); nss <- stats::nobs(fs)
  if (nr != nss) { log_step("!!! %s: n_raw=%d != n_resid=%d -> elpd NOT COMPARABLE, skip", h, nr, nss); next }

  # --- LOO ------------------------------------------------------------------
  lr <- tryCatch(loo::loo(fr, cores = 1), error = function(e){log_step("  LOO raw ERR %s: %s", h, e$message); NULL})
  ls <- tryCatch(loo::loo(fs, cores = 1), error = function(e){log_step("  LOO resid ERR %s: %s", h, e$message); NULL})
  if (is.null(lr) || is.null(ls)) next

  cmp <- loo::loo_compare(list(raw = lr, resid = ls))
  # cmp: best model is the first row (its elpd_diff = 0); the other's diff is signed
  best <- rownames(cmp)[1]
  other <- rownames(cmp)[2]
  elpd_diff <- cmp[other, "elpd_diff"]   # other - best (negative)
  se_diff   <- cmp[other, "se_diff"]
  # standardize: resid - raw
  if (best == "raw") { d_rs <- elpd_diff; s_rs <- se_diff } else { d_rs <- -elpd_diff; s_rs <- se_diff }

  loo_rows[[h]] <- data.frame(
    Model = h,
    elpd_raw   = round(lr$estimates["elpd_loo","Estimate"], 1),
    elpd_resid = round(ls$estimates["elpd_loo","Estimate"], 1),
    elpd_diff_resid_minus_raw = round(d_rs, 2),
    se_diff = round(s_rs, 2),
    favored = ifelse(abs(d_rs) < 2*s_rs, "equivalent", best),
    p_loo_raw   = round(lr$estimates["p_loo","Estimate"], 1),
    p_loo_resid = round(ls$estimates["p_loo","Estimate"], 1),
    n = nr
  )

  # --- dump coefficient raw vs resid ----------------------------------------
  frx <- brms::fixef(fr); fsx <- brms::fixef(fs)
  if (dump_par %in% rownames(frx) && dump_par %in% rownames(fsx)) {
    coef_rows[[h]] <- data.frame(
      Model = h,
      raw_est = round(frx[dump_par, "Estimate"], 4),
      raw_lo  = round(frx[dump_par, "Q2.5"],   4),
      raw_hi  = round(frx[dump_par, "Q97.5"],  4),
      resid_est = round(fsx[dump_par, "Estimate"], 4),
      resid_lo  = round(fsx[dump_par, "Q2.5"],   4),
      resid_hi  = round(fsx[dump_par, "Q97.5"],  4),
      shrink_pct = round(100*(1 - abs(fsx[dump_par,"Estimate"])/abs(frx[dump_par,"Estimate"])), 1)
    )
  } else {
    log_step("%s: dump coefficient not in fixef (raw=%s resid=%s)", h,
             dump_par %in% rownames(frx), dump_par %in% rownames(fsx))
  }
  log_step("%s: elpd raw=%.1f resid=%.1f diff(r-r)=%.2f+-%.2f", h,
           lr$estimates["elpd_loo","Estimate"], ls$estimates["elpd_loo","Estimate"], d_rs, s_rs)
}

if (length(loo_rows)) {
  loo_df <- do.call(rbind, loo_rows); rownames(loo_df) <- NULL
  # -> Table S24 (intensity) / S25 (duration): dump specification (raw vs resid) LOO
  write.csv(loo_df, file.path(tbl_path, paste0(tag_out("Table_V3_07_DumpSpec_LOO"), ".csv")), row.names = FALSE)
  cat("\n=== DumpSpec LOO (raw vs resid) ===\n"); print(loo_df, row.names = FALSE)
} else log_step("WARNING: no hypothesis had a raw+resid pair")

if (length(coef_rows)) {
  coef_df <- do.call(rbind, coef_rows); rownames(coef_df) <- NULL
  # -> Table S26 (intensity) / S27 (duration): dump coefficient, raw vs resid
  write.csv(coef_df, file.path(tbl_path, paste0(tag_out("Table_V3_08_DumpSpec_Coef"), ".csv")), row.names = FALSE)
  cat("\n=== dump coefficient raw vs resid ===\n"); print(coef_df, row.names = FALSE)
}

log_step("=== 06c_dumpspec_compare DONE ===")
sink(type = "message"); sink()
