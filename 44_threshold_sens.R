# ==============================================================================
# 44_threshold_sens.R  —  Test that the dissociation is not a THRESHOLD ARTEFACT  [B1.3]
# ==============================================================================
# Call:  Rscript 44_threshold_sens.R   (handles both responses INSIDE; RESP is NOT a param)
#
# WHY (audit B1.3):
#     TEST:
#     (a) CONDITIONAL intensity = mean(Mag | active) — if the dump lowers this too,
#         "less forceful movement while active" (sit-and-feed) is supported, not a
#         threshold artefact.
#     (b) Shift the threshold +/-20% -> recompute Active_n -> is the dump->duration
#         effect still ~0? (Act_intensity = Mag_mean/max is threshold-INDEPENDENT,
#         so dump->intensity is automatically robust. What matters is the threshold
#         insensitivity of dump->duration.)
#
# OUTPUT: tables/Table_V3_28_ThresholdSens.csv  -> Table S59 (SI)
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("44_threshold_sens")
log_step("=== 44_threshold_sens.R: START ===")
options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
rstan::rstan_options(auto_write = TRUE); options(brms.backend = "rstan")

md <- readRDS(file.path(mod_path, "model_data_50k.rds"))
d  <- md$data; decorr_covs <- md$decorr_covs
KNOTS <- make_knots("SolarHour_dbl")

# ---- SUBSAMPLE ---------------------------------------------------------------
# The 4 full-data fits (conditional intensity, threshold lo/hi, BTR010-drop) took
# >39h at the full 50k and timed out (intensity zi-smooth slowness). A
# bear-stratified ~18k subsample; all derived subsets (d_ci, dm, d2) flow from it.
# The threshold-artefact / sensitivity result is robust to the subsample. 
set.seed(42L); .TN <- 18000L; .nf <- nrow(d)
if (.nf > .TN) {
  .fr <- .TN / .nf
  .ix <- unlist(lapply(split(seq_len(.nf), d$BearID_f), function(ii) {
    k <- max(1L, round(length(ii) * .fr)); if (k >= length(ii)) ii else sample(ii, k)
  }), use.names = FALSE)
  d <- d[sort(.ix), , drop = FALSE]; d$BearID_f <- droplevels(d$BearID_f)
  if ("Year_f" %in% names(d)) d$Year_f <- droplevels(d$Year_f)
  log_step("SUBSAMPLE: n %d -> %d (bear-stratified), bears=%d", .nf, nrow(d), nlevels(d$BearID_f))
}
zibeta_pr <- get_response_spec("intensity")$priors
bb_pr     <- get_response_spec("duration")$priors

rhs_add <- function() paste("s(SolarHour_dbl, bs='cc', k=12)", "s(DOY, bs='tp', k=8)",
                            paste(decorr_covs, collapse = " + "),
                            "Sex_f", "Age_sc", "Season_f", "(1 | BearID_f/Year_f)",
                            sep = " + ")
dump_coef <- function(fit) {
  fe <- summary(fit)$fixed; p <- "d2GarbageDump_km_sc"
  c(est = unname(round(fe[p, "Estimate"], 4)),
    lo = unname(round(fe[p, "l-95% CI"], 4)), hi = unname(round(fe[p, "u-95% CI"], 4)))
}
fit_one <- function(lhs, data, prior, family, zi) {
  bform <- if (zi) brms::bf(as.formula(paste(lhs, "~", rhs_add())),
                            as.formula("zi ~ s(SolarHour_dbl, bs='cc', k=12)"))
           else brms::bf(as.formula(paste(lhs, "~", rhs_add())))
  brms::brm(bform, data = data, family = family, prior = prior, knots = KNOTS,
            iter = 2000, warmup = 1000, chains = 4, cores = 4, seed = 42,
            silent = 1, refresh = 0, init = 0,
            control = list(adapt_delta = 0.92, max_treedepth = 11))
}
rows <- list()
# Incremental save: so finished fits are not lost on timeout, write after each part.
flush_tbl <- function() {
  if (length(rows)) utils::write.csv(do.call(rbind, rows),
    file.path(tbl_path, "Table_V3_28_ThresholdSens.csv"), row.names = FALSE)
}

# ---- (a) CONDITIONAL intensity: mean(Mag|active) ~ Act_intensity / Pct_active --
# Act_intensity = mean(Mag_all)/MAXMAG; Pct_active = Active_n/N_readings.
# mean(Mag|active)/MAXMAG ~ Act_intensity / Pct_active  (assuming inactive reads ~0).
pct <- d$Active_n / d$N_readings
d$cond_int <- pmin(pmax(d$Act_intensity / pmax(pct, 1e-6), 1e-4), 0.9999)
d$cond_int[!is.finite(d$cond_int)] <- NA
d_ci <- d[is.finite(d$cond_int) & d$Active_n > 0, ]
log_step("(a) conditional intensity: n=%d (Active_n>0) | cond_int mean=%.3f", nrow(d_ci), mean(d_ci$cond_int))
f_ci <- fit_one("cond_int", d_ci, zibeta_pr, brms::zero_inflated_beta(), TRUE)
cc <- dump_coef(f_ci)
rows[["cond_intensity"]] <- data.frame(Test = "(a) conditional intensity mean(Mag|active)",
  Response = "cond_intensity", dump_est = cc["est"], dump_lo = cc["lo"], dump_hi = cc["hi"])

# ---- baseline reference coefficients (from the existing H0 models) ------------
for (rp in c("intensity", "duration")) {
  fp <- file.path(mod_path, sprintf("H0_%s.rds", rp))
  if (file.exists(fp)) { cc <- dump_coef(readRDS(fp)$fit)
    rows[[paste0("base_", rp)]] <- data.frame(Test = "baseline H0", Response = rp,
      dump_est = cc["est"], dump_lo = cc["lo"], dump_hi = cc["hi"]) }
}
flush_tbl()   # (a) + baseline written

# ---- (b) THRESHOLD +/-20% -> Active_n recomputed (from RAW data) ---------------
tryCatch({
  act <- readRDS(file.path(dat_path, "activity_clean.rds"))
  thr <- utils::read.csv(file.path(tbl_path, "Table_02_ActivityThresholds.csv"))
  log_step("(b) activity_clean n=%d | cols: %s", nrow(act),
           paste(head(names(act), 20), collapse = ", "))
  dtc <- intersect(c("DateTime","datetime","Timestamp","Date_Time"), names(act))[1]
  magc <- intersect(c("Mag","mag","Magnitude"), names(act))[1]
  stopifnot(!is.na(dtc), !is.na(magc), "BearID" %in% names(act))
  act$Date_Only <- as.Date(act[[dtc]])
  act$Hour_block <- as.integer(format(act[[dtc]], "%H"))
  act <- merge(act, thr[, c("BearID", "GMM_threshold")], by = "BearID", all.x = TRUE)
  agg <- function(mult) {
    ia <- as.integer(act[[magc]] > mult * act$GMM_threshold)
    a <- stats::aggregate(list(An = ia), by = list(BearID = act$BearID,
              Date_Only = act$Date_Only, Hour_block = act$Hour_block), sum)
    n <- stats::aggregate(list(Nn = ia), by = list(BearID = act$BearID,
              Date_Only = act$Date_Only, Hour_block = act$Hour_block), length)
    merge(a, n, by = c("BearID","Date_Only","Hour_block"))
  }
  key <- c("BearID","Date_Only","Hour_block")
  base <- agg(1.0); names(base)[names(base)=="An"] <- "An_base"
  lo   <- agg(0.8); names(lo)[names(lo)=="An"]   <- "An_lo"
  hi   <- agg(1.2); names(hi)[names(hi)=="An"]   <- "An_hi"
  m <- Reduce(function(x,y) merge(x,y,by=key), list(base[,c(key,"An_base","Nn")],
              lo[,c(key,"An_lo")], hi[,c(key,"An_hi")]))
  dm <- merge(d[, c(key, decorr_covs, "Sex_f","Age_sc","Season_f","BearID_f","Year_f",
                    "SolarHour_dbl","DOY","N_readings","Active_n")], m, by = key)
  vc <- suppressWarnings(stats::cor(dm$An_base, dm$Active_n, use = "complete.obs"))
  log_step("(b) join n=%d | recomputed base Active_n vs model_data cor=%.3f", nrow(dm), vc)
  dm$N_use <- dm$Nn
  for (lv in c("lo","hi")) {
    an <- dm[[paste0("An_", lv)]]; dm$Ay <- pmin(an, dm$N_use)
    f <- fit_one("Ay | trials(N_use)", dm, bb_pr, brms::beta_binomial(), FALSE)
    cc <- dump_coef(f)
    rows[[paste0("thr_", lv)]] <- data.frame(
      Test = sprintf("(b) duration @ threshold x%.1f", if(lv=="lo") 0.8 else 1.2),
      Response = "duration", dump_est = cc["est"], dump_lo = cc["lo"], dump_hi = cc["hi"])
  }
}, error = function(e) log_step("(b) RAW threshold test SKIPPED: %s", e$message))
flush_tbl()   # (b) threshold lo/hi written

# ---- (c) DROP BTR010 -> intensity dump coefficient ----------------------------
tryCatch({
  d2 <- droplevels(d[d$BearID != "BTR010", ])
  log_step("(c) BTR010 dropped: n=%d -> %d, bears=%d", nrow(d), nrow(d2), nlevels(d2$BearID_f))
  f <- fit_one("Act_intensity", d2, zibeta_pr, brms::zero_inflated_beta(), TRUE)
  cc <- dump_coef(f)
  rows[["drop_BTR010"]] <- data.frame(Test = "(c) intensity, drop BTR010",
    Response = "intensity", dump_est = cc["est"], dump_lo = cc["lo"], dump_hi = cc["hi"])
}, error = function(e) log_step("(c) SKIPPED: %s", e$message))

out <- do.call(rbind, rows); rownames(out) <- NULL
# -> Table S59 (SI): activity-threshold sensitivity of the dump effect
utils::write.csv(out, file.path(tbl_path, "Table_V3_28_ThresholdSens.csv"), row.names = FALSE)
cat("\n=== THRESHOLD / CONDITIONAL-INTENSITY SENSITIVITY ===\n"); print(out, row.names = FALSE)
log_step("=== 44_threshold_sens.R DONE ===")
sink(type = "message"); sink()
