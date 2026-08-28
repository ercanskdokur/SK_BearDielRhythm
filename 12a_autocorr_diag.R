# ==============================================================================
# 12a_autocorr_diag.R  —  Temporal-autocorrelation DIAGNOSTIC (light, short)
#   (1) H3 residual lag-1 ACF (within bear-day series) -> how big is the problem?
#   (2) verify whether brms ar() works with zero_inflated_beta via a SMALL test fit
#       -> so we don't waste the heavy thinned-AR(1) job.
# OUTPUT: output/tables/Table_V3_09_AutocorrDiag.csv + log
#         -> Table S31 (SI)
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("12a_autocorr_diag")
log_step("=== 12a_autocorr_diag: START ===")

options(brms.backend = "rstan")

h3 <- load_model("H3")
fit <- h3$fit
d   <- h3$data
log_step("H3 loaded: n=%d", nrow(d))

# ---- (1) RESIDUAL LAG-1 ACF (bear-day series) ------------------------------
log_step("Computing residuals (ndraws=300)...")
res <- residuals(fit, summary = TRUE, ndraws = 300)[, "Estimate"]
d$.res <- res
d$.sid <- paste(d$BearID_f, d$Year_f, d$DOY, sep = "_")

split_ix <- split(seq_len(nrow(d)), d$.sid)
acf1 <- vapply(split_ix, function(ix) {
  sub <- d[ix, , drop = FALSE]
  sub <- sub[order(sub$Hour_block), ]
  r <- sub$.res
  if (length(r) < 4 || stats::sd(r, na.rm = TRUE) == 0) return(NA_real_)
  suppressWarnings(stats::cor(r[-length(r)], r[-1], use = "complete.obs"))
}, numeric(1))

len <- vapply(split_ix, length, integer(1))
keep <- !is.na(acf1) & len >= 4
mean_acf1 <- stats::weighted.mean(acf1[keep], len[keep])
med_acf1  <- stats::median(acf1[keep])
log_step("Bear-day series: %d (>=4 hours: %d). Lag-1 ACF: mean=%.4f, median=%.4f",
         length(split_ix), sum(keep), mean_acf1, med_acf1)

# ---- (2) brms ar() SUPPORT TEST (small subset, 1 chain, 200 iter) ----------
set.seed(42)
sids <- unique(d$.sid)
sub_sids <- sample(sids, min(150, length(sids)))
sub <- d[d$.sid %in% sub_sids, , drop = FALSE]
sub$.sid <- factor(sub$.sid)
log_step("ar() test: subset n=%d, %d series", nrow(sub), nlevels(sub$.sid))

ar_status <- "UNKNOWN"; ar_msg <- ""
test <- tryCatch({
  #   this test fit
  #   silently errored and wrote brms_ar_zibeta="FAIL" (a FALSE negative). The
  #   report concluded "ar() cannot be fit with zi_beta" from that; but 12c fits
  #   ar() SUCCESSFULLY on the real response. The correct response is Act_intensity
  #   (the H3 intensity model).
  bform_ar <- brms::bf(
    Act_intensity ~ DumpProx_f + RoadProx_f +
      ar(time = Hour_block, gr = .sid, p = 1),
    zi ~ 1)
  f <- brms::brm(bform_ar, data = sub, family = brms::zero_inflated_beta(),
                 chains = 1, iter = 200, warmup = 100, refresh = 0, silent = 2,
                 seed = 1)
  "OK"
}, error = function(e) { ar_msg <<- conditionMessage(e); "FAIL" })
ar_status <- test
log_step("brms ar() zi_beta support: %s %s", ar_status, ar_msg)

out <- data.frame(
  metric = c("mean_resid_lag1_acf", "median_resid_lag1_acf",
             "n_series", "n_series_used", "brms_ar_zibeta"),
  value  = c(round(mean_acf1, 4), round(med_acf1, 4),
             length(split_ix), sum(keep), ar_status))
# -> Table S31 (SI): residual temporal-autocorrelation diagnostics
write.csv(out, file.path(tbl_path, "Table_V3_09_AutocorrDiag.csv"), row.names = FALSE)
cat("\n=== AUTOCORR DIAG ===\n"); print(out, row.names = FALSE)
if (nzchar(ar_msg)) cat("\nar() error msg:\n", ar_msg, "\n")

log_step("=== 12a_autocorr_diag DONE ===")
cat("DONE 12a\n")
sink(type = "message"); sink()
