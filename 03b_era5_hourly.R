# ==============================================================================
# 03b_era5_hourly.R  —  HOURLY temperature (ERA5-Land) + lapse-rate + Parton-Logan
# ==============================================================================
# FIX: ERA5-Land hourly 2 m temperature (~9 km, hourly reanalysis).
#   The critical advantage of ERA5-Land: it carries REAL weather variation
#   (fronts, cloud, wind). That is information INDEPENDENT of the solar cycle
#   -> not collinear with s(SolarHour) -> the thermal effect is truly identifiable.
#
# SENSITIVITY: Parton-Logan diurnal reconstruction (from E-OBS tn/tx).
#   CAUTION: the P-L within-day shape is a DETERMINISTIC function of the hour ->
#   heavily collinear with s(SolarHour). Identification comes ONLY from the same
#   hour having different tn/tx on different days. Hence it is SENSITIVITY, not
#   primary.
#
# LAPSE-RATE: the study area spans 1300-2900 m. Neither E-OBS (0.1 deg ~11 km) nor
#   ERA5-Land (~9 km) resolves this elevation range. A -6.5 degC/km correction is
#   applied to both sources (the Elevation raster is already extracted).
#
# PREREQUISITE — CDS API (once, by the user):
#   1) https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=download
#      -> ACCEPT the "Terms of use" (otherwise 403)
#   2) https://cds.climate.copernicus.eu/profile -> Personal Access Token
#   3) ~/.cdsapirc:
#        url: https://cds.climate.copernicus.eu/api
#        key: <TOKEN>
#      chmod 600 ~/.cdsapirc ; pip install --user "cdsapi>=0.7.2"
#
# OUTPUTS:
#   data/rasters/ERA5_Land/era5_land_t2m_<YEAR>.nc   (downloaded raw)
#   output/data/era5_hourly_bear.rds                 (bear x hour temperature)
#   output/data/hourly_full.rds                      (UPDATED: temp_hourly added)
#   output/tables/Table_V2_18_TempSourceCompare.csv  (ERA5 vs P-L vs E-OBS daily)
#     NOTE: Table_V2_18 and FigV2_07 below are method-justification outputs; they
#     are NOT numbered SI tables/figures.
#   output/figures/FigV2_07_TempDiel.{png,pdf}
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("03b_era5_hourly")
log_step("=== 03b_era5_hourly.R: START ===")

DO_DOWNLOAD <- as.logical(Sys.getenv("DO_DOWNLOAD", "TRUE"))
era5_dir <- file.path(gis_path, "ERA5_Land")
dir.create(era5_dir, recursive = TRUE, showWarnings = FALSE)

# ---- PYTHON PATH -------------------------------------------------------------
# cdsapi needs Python 3.7+ (cads_api_client uses `from __future__ import
# annotations`) and a working ssl module (so pip can install). A bare
# `import cdsapi` test can HIDE the version problem because the failing import is
# lazy -> the Client() call is what actually fails. Point PYBIN at a Python >=3.7
# with ssl and cdsapi installed. Verified end-to-end (a real ERA5-Land request
# returned "successful" and the .nc downloaded).
PYBIN <- Sys.getenv("PYBIN", unset = "python3")
if (!file.exists(PYBIN)) stop("Python not found: ", PYBIN)
py_ver <- system2(PYBIN, "--version", stdout = TRUE, stderr = TRUE)
log_step("Python: %s (%s)", PYBIN, paste(py_ver, collapse = " "))

# cdsapi + .cdsapirc pre-check: fail BEFORE an hours-long job
if (DO_DOWNLOAD) {
  chk <- system2(PYBIN, c("-c", shQuote(paste(
    "import cdsapi",
    "c = cdsapi.Client()",
    "print('CDS_OK')", sep = "; "))), stdout = TRUE, stderr = TRUE)
  if (!any(grepl("CDS_OK", chk))) {
    stop("CDS API not ready:\n  ", paste(chk, collapse = "\n  "),
         "\n\nFIX:\n",
         "  1) https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=download\n",
         "     -> ACCEPT 'Terms of use' (otherwise the download returns 403)\n",
         "  2) https://cds.climate.copernicus.eu/profile -> Personal Access Token\n",
         "  3) create ~/.cdsapirc:\n",
         "       url: https://cds.climate.copernicus.eu/api\n",
         "       key: <TOKEN>\n",
         "     chmod 600 ~/.cdsapirc\n",
         "  4) python3 -m pip install --user 'cdsapi>=0.7.2'\n")
  }
  log_step("CDS API ready (cdsapi + ~/.cdsapirc verified)")
}

hourly_full <- readRDS(file.path(dat_path, "hourly_full.rds"))
log_step("hourly_full: %s rows", format(nrow(hourly_full), big.mark = ","))

yrs <- sort(unique(hourly_full$Year))
log_step("Years: %s", paste(yrs, collapse = ", "))

# ==============================================================================
# 1. ERA5-LAND DOWNLOAD (via python cdsapi)
# ==============================================================================
# Small bbox + Apr-Nov + hourly -> ~50-150 MB per year. Depending on the CDS queue
# it can take 5-40 min per year. Downloaded year by year so that a single failed
# year does not waste the rest (skip-existing).
bbox_n <- bbox_lat[2] + 0.2; bbox_s <- bbox_lat[1] - 0.2
bbox_w <- bbox_lon[1] - 0.2; bbox_e <- bbox_lon[2] + 0.2
log_step("ERA5 bbox: N=%.1f W=%.1f S=%.1f E=%.1f", bbox_n, bbox_w, bbox_s, bbox_e)

era5_fp <- function(y) file.path(era5_dir, sprintf("era5_land_t2m_%d.nc", y))

if (DO_DOWNLOAD) {
  py_tmpl <- '
import cdsapi, sys
c = cdsapi.Client()
c.retrieve(
    "reanalysis-era5-land",
    {
        "variable": ["2m_temperature"],
        "year": ["%d"],
        "month": ["04","05","06","07","08","09","10","11"],
        "day": [%s],
        "time": [%s],
        "area": [%.2f, %.2f, %.2f, %.2f],
        "data_format": "netcdf",
        "download_format": "unarchived",
    },
    "%s",
)
print("OK %d")
'
  # NOTE: these already contain quotes -> use [%s] in the template, NOT ["%s"].
  # (["%s"] would produce [""01",...,"31""] -> SyntaxError: leading zeros in
  #  decimal integer literals, which silently failed every year.)
  days  <- paste(sprintf("\"%02d\"", 1:31), collapse = ",")
  hours <- paste(sprintf("\"%02d:00\"", 0:23), collapse = ",")

  for (y in yrs) {
    fp <- era5_fp(y)
    if (file.exists(fp) && file.info(fp)$size > 1e5) {
      log_step("  %d: exists (%.1f MB) — skipped", y, file.info(fp)$size / 1024^2)
      next
    }
    py <- sprintf(py_tmpl, y, days, hours, bbox_n, bbox_w, bbox_s, bbox_e, fp, y)
    py_file <- file.path(tmp_dir, sprintf("era5_%d.py", y))
    writeLines(py, py_file)
    log_step("  %d: downloading from CDS... (may be long due to queue)", y)
    t0 <- Sys.time()
    rc <- system2(PYBIN, py_file, stdout = TRUE, stderr = TRUE)
    if (!file.exists(fp)) {
      log_step("  !! %d DOWNLOAD FAILED. cdsapi output:\n%s", y,
               paste(utils::tail(rc, 15), collapse = "\n"))
      log_step("  !! Check: (a) does ~/.cdsapirc exist? (b) ERA5-Land 'Terms of use' accepted?")
    } else {
      log_step("  %d: OK (%.1f MB, %.1f min)", y, file.info(fp)$size / 1024^2,
               as.numeric(difftime(Sys.time(), t0, units = "mins")))
    }
  }
} else {
  log_step("DO_DOWNLOAD=FALSE — existing .nc files will be used")
}

have <- yrs[file.exists(vapply(yrs, era5_fp, character(1)))]
if (length(have) == 0) {
  stop("No ERA5 file present. Check the CDS setup (~/.cdsapirc + licence).")
}
log_step("Available ERA5 years: %s", paste(have, collapse = ", "))

# ---- READING THE ERA5 TIME AXIS ----------------------------------------------
# PITFALL: for netCDF produced by the new CDS, terra::time() returns ALL NA:
#   > terra::time(r)
#   [1] NA NA NA NA
# but the time information lives in the LAYER NAMES:
#   > names(r)
#   "t2m_valid_time=1560384000" "t2m_valid_time=1560405600" ...
# These are Unix epoch seconds (1560384000 = 2019-06-13 00:00:00 UTC; consecutive
# differences of 21600 s = 6 h -> verified).
# So terra::time() is tried first, and if it is NA the layer name is parsed.
era5_times <- function(r) {
  tt <- terra::time(r)
  if (!all(is.na(tt))) {
    return(lubridate::with_tz(as.POSIXct(tt), TZ_DATA))
  }
  nm <- names(r)
  ep <- suppressWarnings(as.numeric(sub(".*valid_time=", "", nm)))
  if (all(is.na(ep))) {
    stop("ERA5 time axis could not be read from terra::time() nor layer names.\n",
         "  Layer names: ", paste(utils::head(nm, 3), collapse = " | "))
  }
  as.POSIXct(ep, origin = "1970-01-01", tz = "UTC")
}

# ==============================================================================
# 2. BEAR x HOUR TEMPERATURE EXTRACTION
# ==============================================================================
# hourly_full stores mean Longitude/Latitude for each (BearID, Date_Only,
# Hour_block) (kept in 03) -> bilinear sampling from the ERA5 grid.
era5_bear <- with_ckpt("03b_era5_extract", {
  key <- hourly_full %>%
    dplyr::filter(!is.na(Longitude), !is.na(Latitude)) %>%
    dplyr::select(BearID, Date_Only, Hour_block, Year, Longitude, Latitude) %>%
    dplyr::distinct()
  log_step("Bear-hours to extract: %s", format(nrow(key), big.mark = ","))

  out <- vector("list", length(have))
  for (i in seq_along(have)) {
    y <- have[i]
    ky <- key[key$Year == y, , drop = FALSE]
    if (nrow(ky) == 0) next
    r <- terra::rast(era5_fp(y))          # time-dimensioned SpatRaster (in K)
    # ERA5 timestamps are UTC; the pipeline is also UTC (see 01) -> direct match.
    # No time-zone conversion — this is the point that stops the shift re-entering.
    tt <- era5_times(r)
    log_step("  %d: %d layers, time %s .. %s (UTC)", y, terra::nlyr(r),
             format(min(tt), "%Y-%m-%d %H:%M"), format(max(tt), "%Y-%m-%d %H:%M"))
    want <- as.POSIXct(paste(ky$Date_Only, sprintf("%02d:00:00", ky$Hour_block)),
                       tz = TZ_DATA)
    idx <- match(want, tt)
    ok <- !is.na(idx)
    log_step("  %d: %s bear-hours, ERA5 layer match %.1f%%",
             y, format(nrow(ky), big.mark = ","), 100 * mean(ok))
    if (!any(ok)) next
    pts <- terra::vect(as.matrix(ky[ok, c("Longitude", "Latitude")]),
                       type = "points", crs = "EPSG:4326")
    # Batch extraction per unique layer (one-by-one would be far too slow)
    vals <- rep(NA_real_, nrow(ky))
    ui <- unique(idx[ok])
    for (k in ui) {
      sel <- which(ok & idx == k)
      pk <- terra::vect(as.matrix(ky[sel, c("Longitude", "Latitude")]),
                        type = "points", crs = "EPSG:4326")
      vv <- terra::extract(r[[k]], pk, method = "bilinear")
      vals[sel] <- vv[[2]]
    }
    ky$t2m_K <- vals
    out[[i]] <- ky
  }
  dplyr::bind_rows(out)
})

era5_bear$temp_era5_raw <- era5_bear$t2m_K - 273.15    # K -> degC
log_step("ERA5 raw temperature: %s", paste(round(summary(era5_bear$temp_era5_raw), 2), collapse = " "))

# ---- LAPSE-RATE CORRECTION ---------------------------------------------------
# The gap between the ERA5-Land grid cell's MEAN elevation and the bear's ACTUAL
# elevation can reach 10 degC over a 1300-2900 m area.
# With no ERA5-Land geopotential, we approximate the grid elevation by aggregating
# the Elevation raster to the ERA5 resolution.
elev_r <- terra::rast(file.path(gis_path, raster_files[["Elevation"]]))
r1 <- terra::rast(era5_fp(have[1]))[[1]]
elev_ll   <- terra::project(elev_r, terra::crs(r1))
elev_grid <- terra::resample(elev_ll, r1, method = "average")
pts_all <- terra::vect(as.matrix(era5_bear[, c("Longitude", "Latitude")]),
                       type = "points", crs = "EPSG:4326")
era5_bear$Elev_grid <- terra::extract(elev_grid, pts_all)[[2]]
era5_bear$Elev_bear <- terra::extract(elev_ll,   pts_all)[[2]]
era5_bear$temp_hourly <- lapse_correct(era5_bear$temp_era5_raw,
                                       era5_bear$Elev_bear, era5_bear$Elev_grid)
log_step("Elevation gap (bear - grid): median %.0f m, range %.0f..%.0f m",
         stats::median(era5_bear$Elev_bear - era5_bear$Elev_grid, na.rm = TRUE),
         min(era5_bear$Elev_bear - era5_bear$Elev_grid, na.rm = TRUE),
         max(era5_bear$Elev_bear - era5_bear$Elev_grid, na.rm = TRUE))
log_step("After lapse correction: %s",
         paste(round(summary(era5_bear$temp_hourly), 2), collapse = " "))
save_rds_safe(era5_bear, file.path(dat_path, "era5_hourly_bear.rds"))

# ==============================================================================
# 3. PARTON-LOGAN SENSITIVITY (from E-OBS tn/tx)
# ==============================================================================
# 04_decorrelate.R dropped temp_max because it correlated with temp_mean (r=0.85).
# Here BOTH tn/tx are NEEDED — not dropped — to build the within-day curve.
pl_ok <- all(c("temp_min", "temp_max", "SR", "SS") %in% names(hourly_full))
if (pl_ok) {
  hourly_full$temp_pl <- parton_logan(
    hour_dec = hourly_full$Hour_block + 0.5,
    tn = hourly_full$temp_min, tx = hourly_full$temp_max,
    SR = hourly_full$SR,       SS = hourly_full$SS
  )
  log_step("Parton-Logan: %s", paste(round(summary(hourly_full$temp_pl), 2), collapse = " "))
} else {
  log_step("!! Parton-Logan skipped — missing: %s",
           paste(setdiff(c("temp_min","temp_max","SR","SS"), names(hourly_full)), collapse = ", "))
  hourly_full$temp_pl <- NA_real_
}

# ==============================================================================
# 4. JOIN onto hourly_full + z-score
# ==============================================================================
hourly_full <- hourly_full %>%
  dplyr::left_join(
    era5_bear %>% dplyr::select(BearID, Date_Only, Hour_block, temp_hourly, temp_era5_raw),
    by = c("BearID", "Date_Only", "Hour_block")
  )
na_pct <- 100 * mean(is.na(hourly_full$temp_hourly))
log_step("temp_hourly NA: %.2f%%", na_pct)
if (na_pct > 10) log_step("!! WARNING: high temp_hourly NA rate — check ERA5 coverage")

hourly_full$temp_hourly_sc <- as.numeric(scale(hourly_full$temp_hourly))
hourly_full$temp_pl_sc     <- as.numeric(scale(hourly_full$temp_pl))
save_rds_safe(hourly_full, file.path(dat_path, "hourly_full.rds"))

# ==============================================================================
# 5. SOURCE COMPARISON — this table IS the justification of the method
# ==============================================================================
# It must show: the daily E-OBS within-day variance is ZERO; ERA5 and P-L are not.
# And ERA5's within-day residual variance (holding the hour fixed) must be LARGER
# than P-L's — because P-L is deterministic while ERA5 carries real weather.
resid_sd <- function(x, hour) {
  ok <- !is.na(x)
  if (sum(ok) < 100) return(NA_real_)
  stats::sd(stats::residuals(stats::lm(x[ok] ~ factor(hour[ok]))))
}
cmp <- data.frame(
  Source = c("E-OBS daily (old)", "Parton-Logan (tn/tx)", "ERA5-Land hourly (new)"),
  Mean = c(mean(hourly_full$temp_mean, na.rm = TRUE),
           mean(hourly_full$temp_pl, na.rm = TRUE),
           mean(hourly_full$temp_hourly, na.rm = TRUE)),
  SD_total = c(stats::sd(hourly_full$temp_mean, na.rm = TRUE),
               stats::sd(hourly_full$temp_pl, na.rm = TRUE),
               stats::sd(hourly_full$temp_hourly, na.rm = TRUE)),
  SD_within_day = c(
    mean(tapply(hourly_full$temp_mean,   paste(hourly_full$BearID, hourly_full$Date_Only), stats::sd), na.rm = TRUE),
    mean(tapply(hourly_full$temp_pl,     paste(hourly_full$BearID, hourly_full$Date_Only), stats::sd), na.rm = TRUE),
    mean(tapply(hourly_full$temp_hourly, paste(hourly_full$BearID, hourly_full$Date_Only), stats::sd), na.rm = TRUE)
  ),
  SD_resid_given_hour = c(
    resid_sd(hourly_full$temp_mean,   hourly_full$Hour_block),
    resid_sd(hourly_full$temp_pl,     hourly_full$Hour_block),
    resid_sd(hourly_full$temp_hourly, hourly_full$Hour_block)
  )
)
cmp[, -1] <- round(cmp[, -1], 3)
cmp$Note <- c("within-day variance = 0 -> H4 NOT TESTABLE",
              "within-day variance present but a deterministic function of hour -> partial identification",
              "within-day variance + REAL weather residual -> H4 a genuine rival")
# NOTE: method-justification table, NOT a numbered SI table.
write.csv(cmp, file.path(tbl_path, "Table_V2_18_TempSourceCompare.csv"), row.names = FALSE)
cat("\n=== TEMPERATURE SOURCE COMPARISON ===\n"); print(cmp, row.names = FALSE)

if (all(!is.na(hourly_full$temp_hourly)) || sum(!is.na(hourly_full$temp_hourly)) > 1000) {
  cr <- stats::cor(hourly_full$temp_hourly, hourly_full$temp_pl, use = "complete.obs")
  log_step("cor(ERA5, Parton-Logan) = %.3f", cr)
  cr2 <- stats::cor(hourly_full$temp_hourly, hourly_full$temp_mean, use = "complete.obs")
  log_step("cor(ERA5, E-OBS daily) = %.3f", cr2)
}

# ==============================================================================
# 6. FIGURE: within-day temperature profile, three sources
# ==============================================================================
pdat <- hourly_full %>%
  dplyr::select(Hour_block, temp_mean, temp_pl, temp_hourly) %>%
  tidyr::pivot_longer(-Hour_block, names_to = "src", values_to = "T") %>%
  dplyr::filter(!is.na(T)) %>%
  dplyr::mutate(src = dplyr::recode(src,
    temp_mean = "E-OBS daily (old)", temp_pl = "Parton-Logan",
    temp_hourly = "ERA5-Land hourly")) %>%
  dplyr::group_by(src, Hour_block) %>%
  dplyr::summarise(m = mean(T), lo = stats::quantile(T, .1),
                   hi = stats::quantile(T, .9), .groups = "drop")

p <- ggplot2::ggplot(pdat, ggplot2::aes(Hour_block, m, colour = src, fill = src)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = .15, colour = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_colour_manual(values = c("E-OBS daily (old)" = "#999999",
                                          "Parton-Logan" = "#E69F00",
                                          "ERA5-Land hourly" = "#0072B2"), name = NULL) +
  ggplot2::scale_fill_manual(values = c("E-OBS daily (old)" = "#999999",
                                        "Parton-Logan" = "#E69F00",
                                        "ERA5-Land hourly" = "#0072B2"), name = NULL) +
  ggplot2::labs(x = "Hour of day (UTC)", y = "2 m temperature (°C)",
                title = "Within-day temperature: why the thermal control needed fixing",
                subtitle = "E-OBS daily is flat by construction — it cannot represent midday-heat avoidance") +
  theme_bear()
save_fig(p, "FigV2_07_TempDiel", w = 9, h = 5)

log_step("=== 03b_era5_hourly.R DONE ===")
sink(type = "message"); sink()
