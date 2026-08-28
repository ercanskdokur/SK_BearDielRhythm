# ==============================================================================
# 05_prep_modeldata.R  —  SHARED model data: 50K subsample + tertiles
# ==============================================================================
# WHY A SEPARATE SCRIPT:
#   There are TWO responses (intensity and duration) and both MUST use the SAME
#   rows, because:
#     - LOO/elpd can only be compared on an identical observation set;
#     - the "dump does not change duration but lowers intensity" contrast is only
#       meaningful on the same rows — different subsamples could create an
#       artificial difference.
#   So the subsample is produced ONCE, in ONE place, and written to disk.
#
# Output: models/model_data_50k.rds
#          $data       : 50K rows, two responses + tertiles + SolarHour
#          $breaks     : tertile cut points (from all data — a fixed reference)
#          $decorr_covs: model covariates (the tertile-cut ones REMOVED)
#          $time_var, $knots
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("05_prep_modeldata")
log_step("=== 05_prep_modeldata.R: START ===")

N_TARGET <- as.integer(Sys.getenv("N_TARGET", "50000"))
TIME_VAR <- get_time_var()
log_step("Diel axis: %s | target n: %d", TIME_VAR, N_TARGET)

# ---- DATA --------------------------------------------------------------------
meta <- readRDS(file.path(mod_path, "decorr_meta.rds"))
decorr_covs <- meta$decorr_covs

d_full <- readRDS(file.path(dat_path, "hourly_decorr.rds"))
d_full <- dplyr::filter(d_full, Season %in% c("Mating", "Hyperphagia"))

# ---- TEMPERATURE: daily -> hourly --------------------------------------------
# Use temp_hourly_sc (ERA5-Land + lapse rate) instead of temp_mean_sc (daily,
# zero within-day variance) for the hour-resolved diel models.
if (!"temp_hourly_sc" %in% names(d_full)) {
  stop("temp_hourly_sc missing — run 03b_era5_hourly.R first.\n",
       "  (Continuing with daily temp_mean leaves H4 UNTESTABLE.)")
}
decorr_covs <- setdiff(decorr_covs, "temp_mean_sc")
decorr_covs <- unique(c(decorr_covs, "temp_hourly_sc"))
log_step("Covariate set (temp_mean_sc -> temp_hourly_sc): %s",
         paste(decorr_covs, collapse = ", "))

need <- c(decorr_covs, "Act_intensity", "Active_n", "N_readings",
          "Sex_f", "Age_sc", "BearID_f", "Year_f", "Season_f", TIME_VAR, "DOY")
need <- unique(need)
if (!"DOY" %in% names(d_full)) {
  d_full$DOY <- as.integer(lubridate::yday(d_full$Date_Only))
}
miss <- setdiff(need, names(d_full))
if (length(miss) > 0) stop("Missing column(s): ", paste(miss, collapse = ", "))

n0 <- nrow(d_full)
d_full <- tidyr::drop_na(d_full, dplyr::any_of(need))
log_step("NA filter: %s -> %s rows, %d bears",
         format(n0, big.mark = ","), format(nrow(d_full), big.mark = ","),
         dplyr::n_distinct(d_full$BearID))

# ---- TERTILES: from ALL data, ONCE -------------------------------------------
# Near = LOW distance (close to humans). cut() assigns labels in increasing x
# order, so the lowest-distance tertile becomes "Near" — correct.
make_tertile <- function(x, labels, brk = NULL) {
  qs <- if (is.null(brk)) stats::quantile(x, c(0, 1/3, 2/3, 1), na.rm = TRUE) else brk
  qs[1] <- -Inf; qs[length(qs)] <- Inf
  cut(x, breaks = unique(qs), labels = labels, include.lowest = TRUE)
}
pick <- function(raw, sc) if (raw %in% names(d_full)) raw else sc
dump_src <- pick("d2GarbageDump_km", "d2GarbageDump_km_sc")
road_src <- pick("d2Roads_km",       "d2Roads_km_sc")
temp_src <- if ("temp_hourly" %in% names(d_full)) "temp_hourly" else "temp_hourly_sc"

brk <- list(
  dump = stats::quantile(d_full[[dump_src]], c(0, 1/3, 2/3, 1), na.rm = TRUE),
  road = stats::quantile(d_full[[road_src]], c(0, 1/3, 2/3, 1), na.rm = TRUE),
  temp = stats::quantile(d_full[[temp_src]], c(0, 1/3, 2/3, 1), na.rm = TRUE)
)
log_step("Tertile cut points (from ALL data — a FIXED reference for every model/script):")
log_step("  dump (%s): %s", dump_src, paste(round(brk$dump, 3), collapse = " | "))
log_step("  road (%s): %s", road_src, paste(round(brk$road, 3), collapse = " | "))
log_step("  temp (%s): %s  <- HOURLY (ERA5), NOT daily", temp_src,
         paste(round(brk$temp, 2), collapse = " | "))

d_full$DumpProx_f <- make_tertile(d_full[[dump_src]], c("Near","Mid","Far"), brk$dump)
d_full$RoadProx_f <- make_tertile(d_full[[road_src]], c("Near","Mid","Far"), brk$road)
d_full$TempTert_f <- make_tertile(d_full[[temp_src]], c("Cool","Mild","Warm"), brk$temp)
d_full$DumpProx_o <- ordered(d_full$DumpProx_f, levels = c("Near","Mid","Far"))
d_full$RoadProx_o <- ordered(d_full$RoadProx_f, levels = c("Near","Mid","Far"))
d_full$TempTert_o <- ordered(d_full$TempTert_f, levels = c("Cool","Mild","Warm"))
d_full <- d_full[stats::complete.cases(
  d_full[, c("DumpProx_f","RoadProx_f","TempTert_f")]), ]

for (v in c("DumpProx_f","RoadProx_f","TempTert_f")) {
  log_step("%s: %s", v, paste(names(table(d_full[[v]])),
           as.integer(table(d_full[[v]])), sep = "=", collapse = " | "))
}

# ---- 50K SUBSAMPLE -----------------------------------------------------------
# NOTE: group_by(BearID_f) %>% slice_sample(prop=p) takes the SAME proportion
# from each bear, so it PRESERVES the bear composition (it does not balance it).
# This is deliberate: keep the sample population-representative.
set.seed(42)
if (nrow(d_full) > N_TARGET) {
  p <- N_TARGET / nrow(d_full)
  d_mod <- d_full %>%
    dplyr::group_by(BearID_f) %>%
    dplyr::slice_sample(prop = p) %>%
    dplyr::ungroup()
  log_step("Subsample: %s -> %s rows (prop=%.4f per bear)",
           format(nrow(d_full), big.mark = ","),
           format(nrow(d_mod), big.mark = ","), p)
} else {
  d_mod <- d_full
  log_step("Subsample unnecessary: n=%d <= %d", nrow(d_mod), N_TARGET)
}

# ---- RESPONSE SANITY CHECKS --------------------------------------------------
log_step("--- INTENSITY response (Act_intensity = mean(Mag)/%.1f) ---", MAXMAG)
log_step("  summary: %s", paste(round(summary(d_mod$Act_intensity), 4), collapse = " "))
log_step("  zero fraction: %.4f | at 0.9999: %d",
         mean(d_mod$Act_intensity == 0), sum(d_mod$Act_intensity >= 0.9999))
if (mean(d_mod$Act_intensity == 0) < 0.01) {
  log_step("  !! zero fraction very low — zero_inflated_beta may be unnecessary")
}
log_step("--- DURATION response (Active_n | trials(N_readings)) ---")
log_step("  Active_n summary: %s", paste(summary(d_mod$Active_n), collapse = " "))
log_step("  N_readings: %s", paste(names(table(d_mod$N_readings)),
         as.integer(table(d_mod$N_readings)), sep = "=", collapse = " | "))
pa <- d_mod$Active_n / d_mod$N_readings
log_step("  proportion: %.2f%% zeros, %.2f%% ones  <- U-shaped; zi-Beta unusable, beta-binom required",
         100 * mean(pa == 0), 100 * mean(pa == 1))

# ---- DIEL AXIS CHECK ---------------------------------------------------------
log_step("--- Diel axis: %s ---", TIME_VAR)
log_step("  range: %.3f - %.3f", min(d_mod[[TIME_VAR]]), max(d_mod[[TIME_VAR]]))
log_step("  cc knots will be set explicitly to c(0,24) (old code did not ->")
log_step("    with data 0..23 mgcv compressed the 24-hour cycle into 23 hours)")

# ---- SAVE --------------------------------------------------------------------
save_rds_safe(list(
  data        = d_mod,
  breaks      = brk,
  sources     = list(dump = dump_src, road = road_src, temp = temp_src),
  decorr_covs = decorr_covs,
  time_var    = TIME_VAR,
  knots       = make_knots(TIME_VAR),
  n           = nrow(d_mod),
  n_bears     = dplyr::n_distinct(d_mod$BearID),
  seed        = 42
), file.path(mod_path, "model_data_50k.rds"))

log_step("SAVED: model_data_50k.rds — n=%d, bears=%d",
         nrow(d_mod), dplyr::n_distinct(d_mod$BearID))
log_step("All 05* scripts now read this file — the subsample is NOT regenerated.")
log_step("=== 05_prep_modeldata.R DONE ===")
sink(type = "message"); sink()
