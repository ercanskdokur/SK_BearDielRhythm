# ==============================================================================
# 48_ar_compare.R  —  Autocorrelation robustness: H3 base vs +bearday-RI vs +AR(1)  [D2]
# ==============================================================================
# Call:  Rscript 48_ar_compare.R   (AFTER 12b + 12c finish)
#
# WHY (audit D2):
#   This table shows
#   whether H3's FIXED effects (dump/road/temperature tertiles + covariates) change
#   when autocorrelation is accounted for -> a direct answer to the autocorrelation
#   objection.
#
# OUTPUT: tables/Table_V3_16_AR_FixedCompare_<RESP>.csv
#   NOTE: Table_V3_16 is NOT a numbered SI table (referee-defense output, referenced
#   in text). Do not confuse it with Table S16 (Pareto-k, produced by 06b).
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("48_ar_compare")
log_step("=== 48_ar_compare: START ===")

getfit <- function(fp) {
  if (!file.exists(fp)) { log_step("!! missing: %s", basename(fp)); return(NULL) }
  x <- readRDS(fp); if (is.list(x) && !is.null(x$fit)) x$fit else x
}
fe_df <- function(fit, tag) {
  if (is.null(fit)) return(NULL)
  fe <- summary(fit)$fixed
  data.frame(Parameter = rownames(fe),
             est = round(fe[, "Estimate"], 4),
             lo = round(fe[, "l-95% CI"], 4), hi = round(fe[, "u-95% CI"], 4),
             stringsAsFactors = FALSE) |> setNames(c("Parameter",
             paste0(c("est_","lo_","hi_"), tag)))
}

for (RESP in c("intensity", "duration")) {
  log_step("--- RESP=%s ---", RESP)
  base <- getfit(file.path(mod_path, sprintf("H3_%s.rds", RESP)))
  brid <- getfit(file.path(mod_path, sprintf("H3_beardayRI_%s.rds", RESP)))
  ar1  <- getfit(file.path(mod_path, sprintf("H3_thinAR1_%s.rds", RESP)))
  tabs <- Filter(Negate(is.null), list(fe_df(base, "base"), fe_df(brid, "beardayRI"),
                                       fe_df(ar1, "AR1")))
  if (length(tabs) == 0) { log_step("!! %s: no model present, skipping", RESP); next }
  out <- Reduce(function(a, b) merge(a, b, by = "Parameter", all = TRUE), tabs)
  # add the ar(1) parameter too (if present)
  if (!is.null(ar1)) {
    cr <- tryCatch(summary(ar1)$cor_pars, error = function(e) NULL)
    if (!is.null(cr)) { cat("\n=== AR(1) cor_pars (", RESP, ") ===\n"); print(cr) }
  }
  # NOT a numbered SI table (see header)
  utils::write.csv(out, file.path(tbl_path, sprintf("Table_V3_16_AR_FixedCompare_%s.csv", RESP)),
                   row.names = FALSE)
  cat(sprintf("\n=== AR FIXED COMPARE (%s) ===\n", RESP)); print(out, row.names = FALSE)
}
log_step("=== 48_ar_compare DONE ===")
sink(type = "message"); sink()
