# ==============================================================================
# 00_setup.R  —  shared configuration, paths, constants and helper functions
# ==============================================================================
# HPC cluster (SLURM); set PROJECT_ROOT to your project directory.
# R: 4.2.0
# This script is source()d by every other script.
# ==============================================================================

# ---- 1. lib paths & scratch ---------------------------------------------------
# In the container (Singularity) environment packages live in the system
# library, so no libPaths override is needed.

tmp_dir <- Sys.getenv("TMPDIR", unset = "/path/to/scratch/tmp")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir)

# ---- 2. PATHS -----------------------------------------------------------------
# Project root holds the analysis subtree BrownBearDielAct/. The raw input data
# (data/) is shared and never modified; all intermediate and final outputs are
# written under the project subtree so a clean run is fully reproducible.
base_path  <- Sys.getenv("PROJECT_ROOT", unset = "/path/to/project")
proj_path  <- file.path(base_path, "BrownBearDielAct")   # analysis root
data_path  <- file.path(base_path, "data")               # raw input (shared)
code_path  <- file.path(proj_path, "scripts")

act_path   <- file.path(data_path, "Activity")
pos_path   <- file.path(data_path, "Position")
meta_path  <- file.path(data_path, "Bear_metadata.xlsx")
gis_path   <- file.path(data_path, "rasters")
shp_path   <- file.path(data_path, "shapefiles")
eobs_path  <- file.path(data_path, "rasters", "E_OBS_Climate")

tbl_path   <- file.path(proj_path, "tables")
fig_path   <- file.path(proj_path, "figures")
mod_path   <- file.path(proj_path, "models")
dat_path   <- file.path(proj_path, "data")        # intermediate outputs (not raw input)
log_path   <- file.path(proj_path, "logs")
ckpt_path  <- file.path(proj_path, "checkpoints")

for (d in c(tbl_path, fig_path, mod_path, dat_path, log_path, ckpt_path)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- 3. PACKAGES (only load — install via 0_install.R) ------------------------
# NOTE: the 'conflicted' package is deliberately NOT loaded. Inside brms::brm()
# it raised an error for every ambiguous function (area, autocor, cluster,
# col_factor, ddirichlet, ...). R's natural "last attached wins" rule already
# gives the priority we want: brms (16) > mixtools (13), brms > ctmm (11),
# terra (8) > tidyverse (2), future (29) > ctmm (11). If 'conflicted' is ever
# loaded, the 'safe_prefer' stub below is a no-op (try-silent), so it does no harm.
pkgs <- c(
  "tidyverse", "lubridate", "readxl", "data.table", "writexl", "tibble",
  "terra", "sf", "sp",
  "ctmm",
  "ncdf4",
  "mixtools", "circular", "suncalc",
  "brms", "tidybayes", "loo", "bayesplot", "posterior",
  "ggplot2", "gridExtra", "grid", "patchwork", "scales", "ggridges",
  "corrplot",
  "parallel", "future", "future.apply"
)

missing_now <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_now) > 0) {
  stop(sprintf("Missing package(s): %s\n  Run 0_install.R first.",
               paste(missing_now, collapse = ", ")))
}

suppressPackageStartupMessages({
  for (p in pkgs) library(p, character.only = TRUE)
})

# ---- 4. NAMESPACE CONFLICTS — explicit pinning --------------------------------
# 'conflicted' is disabled: it raised errors inside brm()'s internal calls.
# R's natural attachment order (last attached wins) already provides the
# priority we want. The 'safe_prefer' calls below are NO-OPs — kept as legacy
# but with no effect. Defensive detach guard:
options(dplyr.summarise.inform = FALSE)
options(conflicts.policy = list(warn = FALSE))

# Defensive: if 'conflicted' was somehow loaded (another script, .Rprofile,
# etc.) detach it so it does not intercept brm()'s internal calls.
if ("package:conflicted" %in% search()) {
  try(detach("package:conflicted", unload = TRUE, force = TRUE), silent = TRUE)
}

# No-op stub — the legacy safe_prefer calls have no effect.
safe_prefer <- function(name, pkg) invisible(NULL)

suppressMessages(suppressWarnings({
  # ---- dplyr: core data manipulation ----
  for (fn in c("filter","select","mutate","summarise","summarize",
               "rename","arrange","lag","lead","first","last","count",
               "between","desc","collapse","slice","intersect","union",
               "setdiff","combine","n","across","cur_data")) {
    safe_prefer(fn, "dplyr")
  }
  # ---- lubridate: date/time ----
  for (fn in c("hour","minute","second","month","year","day","week",
               "wday","yday","date","days","weeks","months","years",
               "hours","minutes","seconds")) {
    safe_prefer(fn, "lubridate")
  }
  # ---- terra: spatial (conflicts with sf, raster, sp, mgcv) ----
  for (fn in c("extract","crs","project","values","origin","time","area")) {
    safe_prefer(fn, "terra")
  }
  # ---- brms: ar, autocor (conflict with stats::ar, ctmm::autocor) ----
  for (fn in c("ar", "autocor")) {
    safe_prefer(fn, "brms")
  }
  # ---- future/parallel: cluster (conflicts with ctmm::cluster) ----
  safe_prefer("cluster", "future")
  # ---- ctmm proactive: known conflicts (distance, projection, raster, extent, meta, zoom) ----
  # ctmm also exports these names; if they clash with another package, let terra/sp win.
  for (fn in c("distance", "projection", "extent", "meta", "zoom", "units")) {
    safe_prefer(fn, "terra")
    safe_prefer(fn, "sp")
  }
  # ---- ggplot2: annotate (conflicts with NLP or base) ----
  safe_prefer("annotate", "ggplot2")
  safe_prefer("Position", "ggplot2")
  # ---- purrr ----
  for (fn in c("discard","keep","transpose","flatten","compact","pluck",
               "set_names","modify")) {
    safe_prefer(fn, "purrr")
  }
  # ---- tidyr ----
  for (fn in c("expand","pack","unpack","extract","fill","complete")) {
    safe_prefer(fn, "tidyr")
  }
  # ---- terra::extract MUST win (override tidyr::extract) ----
  safe_prefer("extract", "terra")
  # ---- stringr ----
  for (fn in c("str_detect","str_replace","str_extract")) {
    safe_prefer(fn, "stringr")
  }
  # ---- forcats ----
  safe_prefer("fct_recode", "forcats")
  # ---- readr: all col_* (conflict with vroom) ----
  for (fn in c("col_character","col_integer","col_double","col_logical",
               "col_factor","col_date","col_datetime","col_time",
               "col_number","col_skip","col_guess","col_big_integer",
               "cols","cols_only","read_csv","read_tsv","read_delim")) {
    safe_prefer(fn, "readr")
  }
  # ---- base operators (conflict with posterior/terra) ----
  for (fn in c("%in%", "Recall", "body<-", "kronecker")) {
    safe_prefer(fn, "base")
  }
  # ---- utils: shadows ctmm/terra head/tail ----
  for (fn in c("head", "tail")) {
    safe_prefer(fn, "utils")
  }
  # ---- stats: the posterior package shadows sd/var/median/quantile/mad/cor ----
  for (fn in c("sd", "var", "median", "mad", "quantile", "cor",
               "weighted.mean", "aggregate", "lag")) {
    safe_prefer(fn, "stats")
  }
}))

# ---- 5. R 4.2.0 S4 dispatch workaround ----------------------------------------
# ctmm calls setGeneric("nrow"); under R 4.2.0 this raises a stopifnot finalize
# error. Define global shadows — ctmm uses its own namespace internally.
nrow <- function(x) base::nrow(x)
ncol <- function(x) base::ncol(x)

# ---- 6. RNG seed --------------------------------------------------------------
set.seed(42)

# ---- 6. terra / raster tmpdir hijack ------------------------------------------
terra::terraOptions(tempdir = tmp_dir, memfrac = 0.6, progress = 0)
if (requireNamespace("raster", quietly = TRUE)) {
  raster::rasterOptions(tmpdir = tmp_dir)
}

# ---- 7. STUDY AREA ------------------------------------------------------------
study_lat  <- 40.33
study_lon  <- 42.59
study_crs  <- "EPSG:32637"
bbox_lat   <- c(40.0, 40.8)
bbox_lon   <- c(42.0, 43.2)

# ---- 7b. TIME ZONE ------------------------------------------------------------
# In the raw Vectronic files the UTC and LMT columns are IDENTICAL (the collar
# LMT offset was never set) -> LMT == UTC. The whole pipeline runs in UTC; sun
# times are also computed in UTC. Local time is used for presentation only.
TZ_DATA  <- "UTC"
TZ_LOCAL <- "Europe/Istanbul"   # for figures/labels only

# ---- 7c. ACCELEROMETER SCALE --------------------------------------------------
# Mag = sqrt(X^2+Y^2+Z^2), each axis 0..255 (a >=255 saturation QC exists in 01)
# -> theoretical maximum 255*sqrt(3) = 441.67. The response is normalized by
# this: Act_intensity = mean(Mag)/MAXMAG in [0,1]. This replaces the old
# arbitrary "Activity/5000" constant (observed max mean(Mag) = 421.4 -> 0.954,
# no ceiling problem) and removes any dependence on N_readings.
MAXMAG <- 255 * sqrt(3)

# ---- 7d. SOLAR TIME -----------------------------------------------------------
# The diel axis is sun-relative time, NOT wall-clock time. Rationale:
#   (a) At the study site sunset shifts ~3 h across the active season -> a
#       wall-clock axis smears the crepuscular peaks.
#   (b) A sun-anchored axis is immune to time-zone errors.
#   (c) The 0 and 24 boundary knots of the cc basis represent the SAME instant
#       biologically.
#
# mode = "sunset" : hours elapsed since sunset [0,24) (same as script 31)
# mode = "double" : double-anchored — stretches day and
#                   night SEPARATELY so both sunrise and sunset align. More
#                   correct for a two-peaked (crepuscular) species; because day
#                   length changes ~10->14 h from April to November, dawn still
#                   drifts under the "sunset" mode.
solar_hour <- function(hour_dec, SR, SS, mode = c("double", "sunset")) {
  mode <- match.arg(mode)
  if (mode == "sunset") return((hour_dec - SS) %% 24)
  # --- double-anchored ---
  # Target: SR -> 6, SS -> 18 (i.e. the standard 12 h day / 12 h night template)
  daylen   <- (SS - SR) %% 24
  nightlen <- 24 - daylen
  h <- hour_dec
  is_day <- ((h - SR) %% 24) < daylen
  out <- numeric(length(h))
  # day: [SR, SS] -> [6, 18]
  out[is_day] <- 6 + 12 * (((h[is_day] - SR[is_day]) %% 24) / daylen[is_day])
  # night: [SS, SR+24] -> [18, 30] -> mod 24
  out[!is_day] <- (18 + 12 * (((h[!is_day] - SS[!is_day]) %% 24) / nightlen[!is_day])) %% 24
  out
}

# ---- 7e. WITHIN-DAY TEMPERATURE: diurnal-curve reconstruction ----------------
# A sensitivity alternative to the ERA5-Land hourly data. E-OBS provides only
# daily tn/tx; this function reconstructs the within-day curve from them.
# WARNING: the within-day shape is a deterministic function of the hour -> it is
# collinear with s(Hour). Identification comes ONLY from the same hour having
# different tn/tx on different days. ERA5-Land should be the primary thermal
# control.
parton_logan <- function(hour_dec, tn, tx, SR, SS, a = 1.86, b = 2.20, c_ = -0.17) {
  daylen   <- (SS - SR) %% 24
  nightlen <- 24 - daylen
  Tsset    <- tn + (tx - tn) * sin(pi * (SS - SR) / (daylen + 2 * a))  # T at sunset
  h <- hour_dec
  is_day <- (h >= SR) & (h <= SS)
  out <- numeric(length(h))
  out[is_day] <- tn[is_day] + (tx[is_day] - tn[is_day]) *
    sin(pi * (h[is_day] - SR[is_day]) / (daylen[is_day] + 2 * a))
  hn <- ifelse(h < SR, h + 24 - SS, h - SS)   # hours elapsed since sunset
  out[!is_day] <- tn[!is_day] + (Tsset[!is_day] - tn[!is_day]) *
    exp(-b * hn[!is_day] / nightlen[!is_day])
  out
}

# ---- 7f. ELEVATION LAPSE-RATE CORRECTION --------------------------------------
# The E-OBS 0.1 degree (~11 km) and ERA5-Land ~9 km grids cannot resolve the
# 1300-2900 m study area. Correct to each individual's true elevation with the
# standard environmental lapse rate. T_bear = T_grid - lapse * (Elev_bear - Elev_grid)
LAPSE_RATE <- 0.0065   # degC / m
lapse_correct <- function(t_grid, elev_bear, elev_grid, lapse = LAPSE_RATE) {
  t_grid - lapse * (elev_bear - elev_grid)
}

# ---- 7g. RESPONSE-VARIABLE SPECIFICATION --------------------------------------
# Two complementary hourly responses are kept as primary because they capture
# DIFFERENT biological dimensions and dissociate at the dump:
#
#   intensity : Act_intensity = mean(Mag)/MAXMAG   [0,1], ~6.6% zeros
#               -> zero_inflated_beta.  "How INTENSELY is it moving?"
#   duration  : Active_n | trials(N_readings)      0..12 integer
#               -> beta_binomial.       "For what FRACTION of the hour is it active?"
#
# Why zi-Beta CANNOT be used for duration: Pct_active is 45.16% zeros + 25.33%
# ones (U-shaped). The Beta distribution can inflate zeros but not ones. The
# beta-binomial handles both boundaries naturally (this is why H0bb comes out
# clean on every statistic in Table S10/S11).
#
# Usage: each 05* script reads the `RESP` env variable (intensity|duration).
get_response_spec <- function(resp = Sys.getenv("RESP", unset = "intensity")) {
  if (!resp %in% c("intensity", "duration")) {
    stop("RESP must be 'intensity' or 'duration', got: ", resp)
  }
  if (resp == "intensity") {
    list(
      key = "intensity",
      lhs = "Act_intensity",
      family = brms::zero_inflated_beta(),
      family_name = "zero_inflated_beta",
      has_zi = TRUE,
      label = "Movement intensity (mean Mag / sensor max)",
      priors = c(
        brms::prior(normal(0, 1),         class = "b"),
        brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
        brms::prior(student_t(3, 0, 2.5), class = "sd"),
        brms::prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "zi"),
        brms::prior(gamma(0.01, 0.01),    class = "phi")
      )
    )
  } else {
    list(
      key = "duration",
      lhs = "Active_n | trials(N_readings)",
      family = brms::beta_binomial(),
      family_name = "beta_binomial",
      has_zi = FALSE,
      label = "Active duration (proportion of readings active)",
      priors = c(
        brms::prior(normal(0, 1),         class = "b"),
        brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
        brms::prior(student_t(3, 0, 2.5), class = "sd"),
        brms::prior(gamma(0.01, 0.01),    class = "phi")
      )
    )
  }
}

# ---- 7g2. MODEL PATH + OUTPUT NAMING ------------------------------------------
# Every hypothesis is fitted TWICE (intensity + duration). Model files are named
# <HYP>_<RESP>.rds; the derived scripts (06b/07b/09*/10*/12*/14/16) read the RESP
# env, open the correct model and write outputs with the SAME suffix. Without
# this the two responses would SILENTLY overwrite each other's tables/figures.
RESP_KEY <- function() Sys.getenv("RESP", unset = "intensity")

model_fp <- function(hyp, resp = RESP_KEY()) {
  file.path(mod_path, sprintf("%s_%s.rds", hyp, resp))
}
# Appends the response suffix to an output name:
#   tag_out("Table_V3_04_LOO") -> "Table_V3_04_LOO_intensity"
tag_out <- function(base, resp = RESP_KEY()) sprintf("%s_%s", base, resp)

# Safely load a model; otherwise a meaningful error.
load_model <- function(hyp, resp = RESP_KEY()) {
  fp <- model_fp(hyp, resp)
  if (!file.exists(fp)) {
    stop(sprintf("Model not found: %s\n  Run first: HYP=%s RESP=%s Rscript 05_fit_models.R",
                 fp, hyp, resp))
  }
  readRDS(fp)
}

# ---- 7h. DIEL AXIS ------------------------------------------------------------
# Primary axis is SolarHour_dbl (double-anchored solar time, produced in 02).
# Can be changed via the TIME_AXIS env: solar_double | solar_sunset | clock
get_time_var <- function(axis = Sys.getenv("TIME_AXIS", unset = "solar_double")) {
  switch(axis,
    solar_double = "SolarHour_dbl",
    solar_sunset = "SolarHour_ss",
    clock        = "Hour_block",
    stop("TIME_AXIS must be solar_double | solar_sunset | clock, got: ", axis))
}


# ---- 8. SEASONS / AGE ---------------------------------------------------------
season_defs <- list(
  Mating      = list(label = "Mating (Apr-Jun)",      short = "Mating",
                     start_md = "04-01", end_md = "06-30", months = 4:6),
  Hyperphagia = list(label = "Hyperphagia (Jul-Nov)", short = "Hyperphagia",
                     start_md = "07-01", end_md = "11-30", months = 7:11),
  Hibernation = list(label = "Hibernation (Dec-Mar)", short = "Hibernation",
                     start_md = "12-01", end_md = "03-31", months = c(12, 1, 2, 3))
)

assign_season <- function(dates) {
  m <- lubridate::month(dates)
  dplyr::case_when(
    m %in% 4:6           ~ "Mating",
    m %in% 7:11          ~ "Hyperphagia",
    m %in% c(12, 1, 2, 3)~ "Hibernation",
    TRUE                 ~ NA_character_
  )
}

assign_age_category <- function(age) {
  dplyr::case_when(
    age >= 5 ~ "Adult",
    age >= 2 ~ "Subadult",
    TRUE     ~ NA_character_
  )
}

# ---- 9. RASTER REGISTRY -------------------------------------------------------
# Environmental raster set (12 layers):
raster_files <- c(
  Elevation        = "Elevation.tif",
  Slope            = "Slope.tif",
  Aspect           = "Aspect.tif",
  Roughness        = "Roughness.tif",
  d2Builtareas     = "d2Builtareas.tif",
  d2Crops          = "d2Crops.tif",
  d2Forest         = "d2Forest.tif",
  d2GarbageDump    = "d2GarbageDump.tif",
  d2ProtectedAreas = "d2ProtectedAreas.tif",
  d2Rangeland      = "d2Rangeland.tif",
  d2Roads          = "d2Roads.tif",
  d2Water          = "d2Water.tif"
)

# Distance rasters (converted to km)
dist_vars <- c("d2Builtareas", "d2Crops", "d2Forest", "d2GarbageDump",
               "d2ProtectedAreas", "d2Rangeland", "d2Roads", "d2Water")

# Topographic rasters (no transform)
topo_vars <- c("Elevation", "Slope", "Aspect", "Roughness")

# Climate variables (from E-OBS NetCDF):
climate_vars <- c("precipitation", "temp_mean", "temp_min", "temp_max")

# _km version used in the models
dist_vars_km <- paste0(dist_vars, "_km")

# All environmental covariates (in the models):
all_covs <- c(topo_vars, dist_vars_km, climate_vars)

# ---- 10. COLOR PALETTE --------------------------------------------------------
pal <- list(
  female      = "#882255",
  male        = "#117733",
  mating      = "#44AA99",
  hyperphagia = "#CC6677",
  hibernation = "#332288",
  subadult    = "#88CCEE",
  adult       = "#AA4499",
  multi       = c("#44AA99", "#CC6677", "#332288", "#DDCC77",
                  "#882255", "#117733", "#88CCEE", "#AA4499"),
  gradient_low  = "#F7F7F7",
  gradient_mid  = "#FDAE61",
  gradient_high = "#D73027",
  day   = "#DDCC77",
  night = "#332288"
)

theme_bear <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = "white", color = NA),
      panel.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major  = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor  = ggplot2::element_blank(),
      axis.line         = ggplot2::element_line(color = "grey30"),
      axis.ticks        = ggplot2::element_line(color = "grey30"),
      axis.text         = ggplot2::element_text(color = "grey20"),
      axis.title        = ggplot2::element_text(color = "grey10", face = "bold"),
      plot.title        = ggplot2::element_text(color = "grey10", face = "bold",
                                                size = base_size + 2, hjust = 0),
      plot.subtitle     = ggplot2::element_text(color = "grey40", size = base_size - 1),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key        = ggplot2::element_rect(fill = "white", color = NA),
      strip.background  = ggplot2::element_rect(fill = "grey92"),
      strip.text        = ggplot2::element_text(face = "bold"),
      plot.margin       = ggplot2::margin(10, 10, 10, 10)
    )
}

# ---- 11. HELPERS --------------------------------------------------------------

# Logger: writes to both stdout and a log file (sink split=TRUE).
# On error the sink is cleaned up; output is always flushed.
init_log <- function(script_name) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  lf <- file.path(log_path, sprintf("%s_%s.log", script_name, ts))
  con <- file(lf, open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  cat(sprintf("=== %s | log: %s | %s ===\n",
              script_name, lf, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("R: %s | host: %s | pid: %d\n",
              R.version.string, Sys.info()[["nodename"]], Sys.getpid()))
  # On error, print a stack trace and clean up the sink.
  options(error = function() {
    cat("\n\n!!! ERROR !!!\n")
    cat(geterrmessage(), "\n")
    cat("\n--- Traceback ---\n")
    try(print(sys.calls()), silent = TRUE)
    sink(type = "message"); try(sink(), silent = TRUE)
    quit(status = 1, save = "no")
  })
  invisible(lf)
}

# Safe RDS write: write to a temp file first, then rename — no partial writes.
save_rds_safe <- function(obj, fp) {
  tmp <- paste0(fp, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, fp)
  cat(sprintf("  [RDS] %s (%.1f MB)\n",
              basename(fp), file.info(fp)$size / 1024^2))
  invisible(fp)
}

log_step <- function(msg, ...) {
  cat(sprintf("[%s] %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              sprintf(msg, ...)))
  flush.console()
}

# Save a figure — both PNG and PDF.
save_fig <- function(p, name_no_ext, fig_dir = fig_path,
                     w = 8, h = 5, dpi = 300) {
  png_path <- file.path(fig_dir, paste0(name_no_ext, ".png"))
  pdf_path <- file.path(fig_dir, paste0(name_no_ext, ".pdf"))
  ggplot2::ggsave(png_path, plot = p, width = w, height = h,
                  dpi = dpi, bg = "white")
  ggplot2::ggsave(pdf_path, plot = p, width = w, height = h,
                  bg = "white", device = "pdf")
  cat(sprintf("  Saved: %s.{png,pdf}\n", name_no_ext))
}

# Checkpoint helper: load if the file exists, otherwise evaluate expr and save.
with_ckpt <- function(name, expr, force = FALSE) {
  fp <- file.path(ckpt_path, paste0(name, ".rds"))
  force_env <- as.logical(Sys.getenv("FORCE_RECOMPUTE", "FALSE"))
  if (file.exists(fp) && !force && !isTRUE(force_env)) {
    log_step("ckpt HIT  : %s", name)
    return(readRDS(fp))
  }
  log_step("ckpt MISS : %s — computing...", name)
  res <- force(expr)
  saveRDS(res, fp)
  log_step("ckpt SAVE : %s", basename(fp))
  res
}

# Vectronic activity file reader
read_activity_file <- function(fp) {
  cn <- c("No", "CollarID", "UTC_Date", "UTC_Time", "Date", "Time",
          "Origin", "SCTS_Date", "SCTS_Time", "Act_Mode", "DT",
          "ActivityX", "ActivityY", "ActivityZ", "Temp", "AnimalID", "GroupID")
  tryCatch({
    df <- readr::read_table(fp, skip = 3, col_names = cn,
                            na = c("NA", "N/A", "", " "), show_col_types = FALSE)
    df <- df %>%
      dplyr::mutate(dplyr::across(c(ActivityX, ActivityY, ActivityZ), as.numeric))
    df
  }, error = function(e) {
    warning(sprintf("Activity read failed: %s — %s", basename(fp), e$message))
    NULL
  })
}

# Vectronic GPS file reader — field parsing via regex
read_position_file <- function(fp) {
  tryCatch({
    lines <- readLines(fp, warn = FALSE, encoding = "latin1")
    lines <- iconv(lines, from = "latin1", to = "UTF-8", sub = "")
    lines <- gsub("\xEF\xBB\xBF", "", lines)
    lines <- gsub("\r", "", lines)

    data_lines <- lines[-(1:3)]
    data_lines <- data_lines[nchar(trimws(data_lines)) > 0]
    if (length(data_lines) == 0) return(NULL)

    results <- vector("list", length(data_lines))
    for (i in seq_along(data_lines)) {
      fields <- strsplit(trimws(data_lines[i]), "\\s+")[[1]]
      nf <- length(fields)
      if (nf < 14) next

      lmt_date <- fields[5]
      lmt_time <- fields[6]
      origin   <- fields[7]
      lat <- NA_real_; lon <- NA_real_

      for (j in 10:min(nf, 20)) {
        val_str <- fields[j]
        if (grepl("^[0-9]+,[0-9]{5,}$", val_str)) {
          val_num <- as.numeric(gsub(",", ".", val_str))
          if (is.na(lat) && !is.na(val_num)) {
            lat <- val_num
          } else if (!is.na(lat) && is.na(lon) && !is.na(val_num)) {
            lon <- val_num
            break
          }
        }
      }

      dop <- NA_real_
      if (!is.na(lon)) {
        for (j in 10:min(nf, 20)) {
          if (grepl("^[0-9]+,[0-9]{5,}$", fields[j])) {
            val_num <- as.numeric(gsub(",", ".", fields[j]))
            if (abs(val_num - lon) < 0.0001) {
              if (j + 2 <= nf) {
                dop <- suppressWarnings(as.numeric(gsub(",", ".", fields[j + 2])))
              }
              break
            }
          }
        }
      }

      results[[i]] <- data.frame(
        Date = lmt_date, Time = lmt_time, Origin = origin,
        Latitude = lat, Longitude = lon, DOP = dop,
        stringsAsFactors = FALSE
      )
    }

    df <- dplyr::bind_rows(results)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df
  }, error = function(e) {
    warning(sprintf("Position read failed: %s — %s", basename(fp), e$message))
    NULL
  })
}

# Activity threshold via a 2-component GMM
get_gmm_threshold <- function(mag_vector) {
  clean <- mag_vector[mag_vector > 0 & mag_vector < 255]
  if (length(clean) < 50) return(stats::median(clean, na.rm = TRUE))
  tryCatch({
    m   <- mixtools::normalmixEM(log(clean), k = 2, maxrestarts = 50, verb = FALSE)
    idx <- order(m$mu)
    mu1 <- m$mu[idx[1]]; s1 <- m$sigma[idx[1]]; l1 <- m$lambda[idx[1]]
    mu2 <- m$mu[idx[2]]; s2 <- m$sigma[idx[2]]; l2 <- m$lambda[idx[2]]
    f   <- function(x) l1 * stats::dnorm(x, mu1, s1) - l2 * stats::dnorm(x, mu2, s2)
    root <- tryCatch(stats::uniroot(f, c(mu1, mu2))$root,
                     error = function(e) mean(c(mu1, mu2)))
    exp(root)
  }, error = function(e) stats::median(clean, na.rm = TRUE))
}

# E-OBS bilinear interpolation
get_eobs_value <- function(lon, lat, date, nc, var_name, nc_lon, nc_lat, nc_dates) {
  time_idx <- which(nc_dates == as.Date(date))
  if (length(time_idx) == 0) return(NA_real_)

  lon_idx_low  <- suppressWarnings(max(which(nc_lon <= lon)))
  lon_idx_high <- suppressWarnings(min(which(nc_lon >= lon)))
  lat_idx_low  <- suppressWarnings(max(which(nc_lat <= lat)))
  lat_idx_high <- suppressWarnings(min(which(nc_lat >= lat)))

  if (!is.finite(lon_idx_low)  || !is.finite(lon_idx_high) ||
      !is.finite(lat_idx_low)  || !is.finite(lat_idx_high)) return(NA_real_)

  if (lon_idx_low == lon_idx_high && lat_idx_low == lat_idx_high) {
    val <- ncdf4::ncvar_get(nc, var_name,
                            start = c(lon_idx_low, lat_idx_low, time_idx),
                            count = c(1, 1, 1))
    return(val)
  }

  corners <- expand.grid(lon_idx = c(lon_idx_low, lon_idx_high),
                         lat_idx = c(lat_idx_low, lat_idx_high))
  values <- numeric(4)
  for (i in 1:4) {
    values[i] <- ncdf4::ncvar_get(nc, var_name,
                                  start = c(corners$lon_idx[i],
                                            corners$lat_idx[i],
                                            time_idx),
                                  count = c(1, 1, 1))
  }

  if (all(is.na(values))) return(NA_real_)

  lon_w <- if (lon_idx_low == lon_idx_high) 0 else
    (lon - nc_lon[lon_idx_low]) / (nc_lon[lon_idx_high] - nc_lon[lon_idx_low])
  lat_w <- if (lat_idx_low == lat_idx_high) 0 else
    (lat - nc_lat[lat_idx_low]) / (nc_lat[lat_idx_high] - nc_lat[lat_idx_low])

  if (all(!is.na(values))) {
    values[1] * (1 - lon_w) * (1 - lat_w) +
      values[2] * lon_w       * (1 - lat_w) +
      values[4] * lon_w       * lat_w       +
      values[3] * (1 - lon_w) * lat_w
  } else {
    mean(values, na.rm = TRUE)
  }
}

cat("=== 00_setup.R ready ===\n")
