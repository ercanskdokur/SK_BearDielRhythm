# ==============================================================================
# 03_covariate_extraction.R
# ==============================================================================
# 1) ctmm GPS interpolation (30 min)
# 2) GIS raster values (12 layers: 4 topographic + 8 distance)
# 3) E-OBS climate data (4 daily NetCDFs; bilinear interpolation)
# 4) Join with the hourly / daily / bout response tables
#
# Note on tables: this script writes an intermediate correlation matrix over the
# E-OBS climate set. The FINAL SI Table S6 (Kendall correlation) is produced by
# corr_era5.R over the modelling covariate set (precipitation + ERA5 temp_hourly).
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("03_covariate_extraction")
log_step("=== 03_covariate_extraction.R: START ===")

`%||%` <- function(a, b) if (!is.null(a)) a else b

pos_all          <- readRDS(file.path(dat_path, "position_clean.rds"))
hourly_intensity <- readRDS(file.path(dat_path, "hourly_intensity.rds"))
daily_active     <- readRDS(file.path(dat_path, "daily_active.rds"))
active_bouts     <- readRDS(file.path(dat_path, "active_bouts.rds"))
meta             <- readRDS(file.path(dat_path, "metadata.rds"))

# ==============================================================================
# 1. LOAD RASTERS (12-layer set)
# ==============================================================================
log_step("Loading rasters (%d layers)...", length(raster_files))

raster_list <- list()
for (nm in names(raster_files)) {
  fp <- file.path(gis_path, raster_files[nm])
  if (file.exists(fp)) {
    r <- terra::rast(fp)
    names(r) <- nm
    raster_list[[nm]] <- r
  } else {
    log_step("  WARNING: raster not found -> %s", fp)
  }
}

ref_raster <- raster_list[[1]]
rcrs <- terra::crs(ref_raster, describe = TRUE)
log_step("  Raster CRS: %s (EPSG:%s)", rcrs$name, rcrs$code)
log_step("  Loaded: %d / %d", length(raster_list), length(raster_files))

for (nm in names(raster_list)) {
  e <- as.vector(terra::ext(raster_list[[nm]]))
  log_step("    %s: ext=[%.0f,%.0f,%.0f,%.0f] res=%.1f",
           nm, e[1], e[2], e[3], e[4], terra::res(raster_list[[nm]])[1])
}

# ==============================================================================
# CKPT: GPS INTERPOLATION (slowest step)
# ==============================================================================
pos_interp <- with_ckpt("03a_pos_interp", {
  log_step("ctmm GPS interpolation (30 min)...")
  pos_sf <- sf::st_as_sf(pos_all, coords = c("Longitude", "Latitude"),
                         crs = 4326) %>%
    sf::st_transform(as.numeric(gsub("EPSG:", "", study_crs)))
  pos_all$X_utm <- sf::st_coordinates(pos_sf)[, 1]
  pos_all$Y_utm <- sf::st_coordinates(pos_sf)[, 2]

  bears <- unique(pos_all$BearID)
  log_step("  n bears: %d", length(bears))
  interp_list <- list()

  for (b in seq_along(bears)) {
    bear <- bears[b]
    df_bear <- pos_all %>%
      filter(BearID == bear) %>%
      arrange(DateTime) %>%
      filter(!duplicated(DateTime))

    if (nrow(df_bear) < 30) {
      log_step("    %s: too few fixes (%d) -> raw data", bear, nrow(df_bear))
      interp_list[[bear]] <- df_bear %>%
        select(BearID, DateTime, Date_Only, Year, Month,
               Longitude, Latitude, X_utm, Y_utm,
               Sex, Current_Age, Age_Category, Season, BearYear)
      next
    }

    tryCatch({
      telem_df <- data.frame(
        individual.local.identifier = bear,
        timestamp     = format(df_bear$DateTime, "%Y-%m-%d %H:%M:%S"),
        location.long = df_bear$Longitude,
        location.lat  = df_bear$Latitude
      )
      # TZ: switched to UTC in 01 (collar LMT == UTC). ctmm must also use UTC,
      # otherwise the interpolation time axis shifts by 3 hours.
      telem <- ctmm::as.telemetry(telem_df, timeformat = "%Y-%m-%d %H:%M:%S",
                                  timezone = "UTC")
      guess <- ctmm::ctmm.guess(telem, interactive = FALSE)
      fit <- tryCatch(
        ctmm::ctmm.select(telem, guess, verbose = FALSE, cores = 1),
        error = function(e) ctmm::ctmm.fit(telem, guess)
      )
      if (is.list(fit) && !inherits(fit, "ctmm")) fit <- fit[[1]]

      t_start <- min(df_bear$DateTime)
      t_end   <- max(df_bear$DateTime)
      t_seq   <- seq(t_start, t_end, by = "30 min")
      pred <- predict(telem, fit, t = as.numeric(t_seq))

      # ctmm telemetry is in projected (x,y meters) coords; back-project to lon/lat
      telem_proj <- tryCatch(attr(telem, "projection"),
                             error = function(e) NULL)
      if (is.null(telem_proj) || !nzchar(telem_proj)) {
        telem_proj <- tryCatch(ctmm::projection(telem),
                               error = function(e) NULL)
      }

      pred_x <- pred$x %||% pred[["x"]]
      pred_y <- pred$y %||% pred[["y"]]
      if (is.null(pred_x) || is.null(pred_y)) {
        stop("x/y not found in ctmm predict output")
      }

      tmp_xy <- data.frame(
        x = as.numeric(pred_x),
        y = as.numeric(pred_y)
      )
      pred_sf_proj <- sf::st_as_sf(tmp_xy, coords = c("x", "y"),
                                   crs = telem_proj)
      pred_sf_ll <- sf::st_transform(pred_sf_proj, 4326)
      pred_ll <- sf::st_coordinates(pred_sf_ll)

      pred_df <- data.frame(
        BearID    = bear,
        DateTime  = as.POSIXct(as.numeric(pred$t),
                               origin = "1970-01-01", tz = "UTC"),   # aligned with 01
        Longitude = pred_ll[, 1],
        Latitude  = pred_ll[, 2],
        stringsAsFactors = FALSE
      )

      # Drop large gaps (>12h)
      df_bear <- df_bear %>% arrange(DateTime)
      rt <- df_bear$DateTime
      gs <- c(); ge <- c()
      for (j in 2:length(rt)) {
        gh <- as.numeric(difftime(rt[j], rt[j - 1], units = "hours"))
        if (gh > 12) { gs <- c(gs, rt[j - 1]); ge <- c(ge, rt[j]) }
      }
      if (length(gs) > 0) {
        in_gap <- rep(FALSE, nrow(pred_df))
        for (g in seq_along(gs)) {
          in_gap <- in_gap | (pred_df$DateTime > gs[g] & pred_df$DateTime < ge[g])
        }
        pred_df <- pred_df[!in_gap, ]
      }

      pred_sf <- sf::st_as_sf(pred_df, coords = c("Longitude", "Latitude"),
                              crs = 4326) %>%
        sf::st_transform(as.numeric(gsub("EPSG:", "", study_crs)))
      pred_df$X_utm <- sf::st_coordinates(pred_sf)[, 1]
      pred_df$Y_utm <- sf::st_coordinates(pred_sf)[, 2]
      pred_df$Date_Only    <- as.Date(pred_df$DateTime)
      pred_df$Year         <- lubridate::year(pred_df$DateTime)
      pred_df$Month        <- lubridate::month(pred_df$DateTime)
      pred_df$Sex          <- df_bear$Sex[1]
      pred_df$Current_Age  <- df_bear$Current_Age[1]
      pred_df$Age_Category <- df_bear$Age_Category[1]
      pred_df$Season       <- assign_season(pred_df$DateTime)
      pred_df$BearYear     <- paste(bear, pred_df$Year, sep = "_")

      interp_list[[bear]] <- pred_df
      if (b %% 5 == 0) log_step("    %d/%d (%s, %d fixes)",
                                 b, length(bears), bear, nrow(pred_df))
    }, error = function(e) {
      log_step("    %s: ctmm error (%s) -> raw data", bear, e$message)
      interp_list[[bear]] <<- df_bear %>%
        select(BearID, DateTime, Date_Only, Year, Month,
               Longitude, Latitude, X_utm, Y_utm,
               Sex, Current_Age, Age_Category, Season, BearYear)
    })
  }

  pp <- bind_rows(interp_list)
  pp %>%
    filter(Month >= 4, Month <= 11,
           Latitude  >= bbox_lat[1], Latitude  <= bbox_lat[2],
           Longitude >= bbox_lon[1], Longitude <= bbox_lon[2])
})
log_step("  pos_interp: %s", format(nrow(pos_interp), big.mark = ","))

# ==============================================================================
# CKPT: GIS RASTER VALUES
# ==============================================================================
pos_interp <- with_ckpt("03b_pos_with_gis", {
  log_step("Extracting GIS raster values...")
  pp <- pos_interp
  pos_vect <- terra::vect(pp, geom = c("X_utm", "Y_utm"), crs = study_crs)

  for (nm in names(raster_list)) {
    r <- raster_list[[nm]]
    r_crs <- terra::crs(r)
    if (terra::crs(pos_vect) != r_crs) {
      pv <- terra::project(pos_vect, r_crs)
    } else {
      pv <- pos_vect
    }
    vals <- terra::extract(r, pv)
    pp[[nm]] <- vals[[nm]]
    na_n <- sum(is.na(pp[[nm]]))
    log_step("    %s: %d NA (%.1f%%)", nm, na_n, na_n / nrow(pp) * 100)
  }

  # Convert distance rasters to km
  for (v in dist_vars) {
    if (v %in% names(pp)) pp[[paste0(v, "_km")]] <- pp[[v]] / 1000
  }
  pp
})

# ==============================================================================
# CKPT: E-OBS CLIMATE DATA
# ==============================================================================
daily_pos <- with_ckpt("03c_eobs_daily", {
  log_step("E-OBS NetCDF climate data (path: %s)", eobs_path)

  eobs_vars <- data.frame(
    file = c("rr_ens_mean_0.1deg_reg_v31.0e.nc",
             "tg_ens_mean_0.1deg_reg_v31.0e.nc",
             "tn_ens_mean_0.1deg_reg_v31.0e.nc",
             "tx_ens_mean_0.1deg_reg_v31.0e.nc"),
    var_name = c("rr", "tg", "tn", "tx"),
    col_name = c("precipitation", "temp_mean", "temp_min", "temp_max"),
    stringsAsFactors = FALSE
  )

  # Daily position (BearID, Date_Only, mean lon/lat)
  dp <- pos_interp %>%
    group_by(BearID, Date_Only) %>%
    summarise(
      mean_lon = mean(Longitude, na.rm = TRUE),
      mean_lat = mean(Latitude,  na.rm = TRUE),
      .groups = "drop"
    )
  log_step("  bear-days: %d", nrow(dp))

  for (v in seq_len(nrow(eobs_vars))) {
    info <- eobs_vars[v, ]
    log_step("  [%d/%d] %s (%s)", v, nrow(eobs_vars), info$col_name, info$var_name)
    nc_file <- file.path(eobs_path, info$file)
    if (!file.exists(nc_file)) {
      log_step("    FILE MISSING: %s", nc_file)
      dp[[info$col_name]] <- NA_real_
      next
    }
    nc <- ncdf4::nc_open(nc_file)
    nc_lon <- ncdf4::ncvar_get(nc, "longitude")
    nc_lat <- ncdf4::ncvar_get(nc, "latitude")
    nc_t   <- ncdf4::ncvar_get(nc, "time")
    nc_d   <- as.Date("1950-01-01") + nc_t

    dp[[info$col_name]] <- NA_real_
    for (i in seq_len(nrow(dp))) {
      dp[[info$col_name]][i] <- get_eobs_value(
        lon = dp$mean_lon[i], lat = dp$mean_lat[i],
        date = dp$Date_Only[i],
        nc = nc, var_name = info$var_name,
        nc_lon = nc_lon, nc_lat = nc_lat, nc_dates = nc_d
      )
      if (i %% 2000 == 0) log_step("    %d / %d", i, nrow(dp))
    }
    ncdf4::nc_close(nc)
    nan <- sum(is.na(dp[[info$col_name]]))
    log_step("    done, NA: %d (%.1f%%)", nan, nan / nrow(dp) * 100)
  }
  dp
})

# ==============================================================================
# 5. RESPONSE-COVARIATE JOIN
# ==============================================================================
log_step("Response x covariate join...")
raster_names <- names(raster_files)

# Hourly: mean over BearID x Date_Only x Hour_block
# Longitude/Latitude are also kept (needed for the spatial grid in later scripts)
hourly_covs <- pos_interp %>%
  mutate(Hour_block = lubridate::hour(DateTime)) %>%
  group_by(BearID, Date_Only, Hour_block) %>%
  summarise(
    Longitude = mean(Longitude, na.rm = TRUE),
    Latitude  = mean(Latitude,  na.rm = TRUE),
    across(all_of(c(raster_names, dist_vars_km)), ~ mean(.x, na.rm = TRUE)),
    N_fixes_hour = n(),
    .groups = "drop"
  ) %>%
  left_join(daily_pos %>% select(BearID, Date_Only,
                                 precipitation, temp_mean, temp_min, temp_max),
            by = c("BearID", "Date_Only"))

hourly_full <- hourly_intensity %>%
  left_join(hourly_covs, by = c("BearID", "Date_Only", "Hour_block"))
log_step("  hourly_full: %s rows x %d cols",
         format(nrow(hourly_full), big.mark = ","), ncol(hourly_full))

# Daily: mean over BearID x Date_Only (+ Longitude/Latitude for spatial)
daily_covs <- pos_interp %>%
  group_by(BearID, Date_Only) %>%
  summarise(
    Longitude = mean(Longitude, na.rm = TRUE),
    Latitude  = mean(Latitude,  na.rm = TRUE),
    across(all_of(c(raster_names, dist_vars_km)), ~ mean(.x, na.rm = TRUE)),
    N_fixes_day = n(),
    .groups = "drop"
  ) %>%
  left_join(daily_pos %>% select(BearID, Date_Only,
                                 precipitation, temp_mean, temp_min, temp_max),
            by = c("BearID", "Date_Only"))

daily_full <- daily_active %>%
  left_join(daily_covs, by = c("BearID", "Date_Only"))
log_step("  daily_full: %s rows x %d cols",
         format(nrow(daily_full), big.mark = ","), ncol(daily_full))

# Bout: mean of the fixes falling within the bout window
bout_covs <- active_bouts %>%
  select(BearID, Date_Only, Bout_start, Bout_end, Bout_length_min)

bout_gis_list <- list()
for (bear in unique(bout_covs$BearID)) {
  bear_pos   <- pos_interp %>% filter(BearID == bear)
  bear_bouts <- bout_covs  %>% filter(BearID == bear)
  if (nrow(bear_pos) == 0 || nrow(bear_bouts) == 0) next

  for (i in seq_len(nrow(bear_bouts))) {
    bo <- bear_bouts[i, ]
    bo_mid <- bo$Bout_start + (bo$Bout_end - bo$Bout_start) / 2
    bf <- bear_pos %>%
      filter(DateTime >= bo$Bout_start - 900,
             DateTime <= bo$Bout_end   + 900)
    if (nrow(bf) == 0) {
      td <- abs(as.numeric(difftime(bear_pos$DateTime, bo_mid, units = "mins")))
      cl <- which.min(td)
      if (td[cl] <= 60) bf <- bear_pos[cl, ] else next
    }
    gis_avg <- bf %>%
      summarise(across(all_of(c(raster_names, dist_vars_km)),
                       ~ mean(.x, na.rm = TRUE)))
    gis_avg$BearID     <- bear
    gis_avg$Bout_start <- bo$Bout_start
    bout_gis_list[[length(bout_gis_list) + 1]] <- gis_avg
  }
  if (which(unique(bout_covs$BearID) == bear) %% 10 == 0) {
    log_step("    bout cov %d / %d",
             which(unique(bout_covs$BearID) == bear),
             length(unique(bout_covs$BearID)))
  }
}
bout_gis <- bind_rows(bout_gis_list)

bout_full <- active_bouts %>%
  left_join(bout_gis, by = c("BearID", "Bout_start")) %>%
  left_join(daily_pos %>% select(BearID, Date_Only,
                                 precipitation, temp_mean, temp_min, temp_max),
            by = c("BearID", "Date_Only"))
log_step("  bout_full: %s rows x %d cols",
         format(nrow(bout_full), big.mark = ","), ncol(bout_full))

# ==============================================================================
# 6. CORRELATION ANALYSIS (Kendall, |r|>0.7 flag)
# ==============================================================================
# NOTE: this is an INTERMEDIATE correlation over the E-OBS climate set. The final
# SI Table S6 (and Fig S3b) use the modelling covariate set (precipitation +
# ERA5 temp_hourly) and are produced by corr_era5.R.
log_step("Correlation (Kendall) analysis...")
covariate_cols <- c(topo_vars, dist_vars_km, climate_vars)

cor_data <- daily_full %>%
  select(any_of(covariate_cols)) %>%
  tidyr::drop_na()

if (nrow(cor_data) > 100) {
  cor_matrix <- cor(cor_data, method = "kendall", use = "pairwise.complete.obs")
  # Intermediate (E-OBS climate); superseded by corr_era5.R -> Table S6.
  write.csv(as.data.frame(cor_matrix),
            file.path(tbl_path, "Table_04_CorrelationMatrix.csv"))

  png(file.path(fig_path, "Fig_00_CorrelationMatrix.png"),
      width = 10, height = 10, units = "in", res = 300)
  corrplot::corrplot(cor_matrix, method = "color", type = "lower",
                     tl.col = "black", tl.srt = 45, tl.cex = 0.8,
                     addCoef.col = "black", number.cex = 0.6,
                     col = colorRampPalette(c("#332288", "white", "#CC6677"))(100),
                     title = "Kendall's Rank Correlation - Environmental Covariates",
                     mar = c(0, 0, 2, 0))
  dev.off()

  pdf(file.path(fig_path, "Fig_00_CorrelationMatrix.pdf"),
      width = 10, height = 10)
  corrplot::corrplot(cor_matrix, method = "color", type = "lower",
                     tl.col = "black", tl.srt = 45, tl.cex = 0.8,
                     addCoef.col = "black", number.cex = 0.6,
                     col = colorRampPalette(c("#332288", "white", "#CC6677"))(100),
                     title = "Kendall's Rank Correlation - Environmental Covariates",
                     mar = c(0, 0, 2, 0))
  dev.off()

  hc <- which(abs(cor_matrix) > 0.7 & abs(cor_matrix) < 1, arr.ind = TRUE)
  if (nrow(hc) > 0) {
    log_step("  HIGH CORRELATION (|r|>0.7):")
    for (i in seq_len(nrow(hc))) {
      log_step("    %s ~ %s : r=%.3f",
               rownames(cor_matrix)[hc[i, 1]],
               colnames(cor_matrix)[hc[i, 2]],
               cor_matrix[hc[i, 1], hc[i, 2]])
    }
  } else {
    log_step("  No high correlation.")
  }
}

# ==============================================================================
# 7. SAVE RDS
# ==============================================================================
log_step("Saving RDS...")
save_rds_safe(hourly_full, file.path(dat_path, "hourly_full.rds"))
save_rds_safe(daily_full,  file.path(dat_path, "daily_full.rds"))
save_rds_safe(bout_full,   file.path(dat_path, "bout_full.rds"))
save_rds_safe(pos_interp,  file.path(dat_path, "position_interpolated.rds"))
save_rds_safe(daily_pos,   file.path(dat_path, "daily_positions_climate.rds"))

write.csv(head(hourly_full, 1000), file.path(dat_path, "hourly_full_preview.csv"), row.names = FALSE)
write.csv(head(daily_full,  1000), file.path(dat_path, "daily_full_preview.csv"),  row.names = FALSE)
write.csv(head(bout_full,   1000), file.path(dat_path, "bout_full_preview.csv"),   row.names = FALSE)

cat("\n========================================\n")
cat("  COVARIATE EXTRACTION SUMMARY\n")
cat("========================================\n")
cat(sprintf("pos_interp:  %s\n", format(nrow(pos_interp), big.mark = ",")))
cat(sprintf("hourly_full: %s rows x %d cols\n",
            format(nrow(hourly_full), big.mark = ","), ncol(hourly_full)))
cat(sprintf("daily_full:  %s rows x %d cols\n",
            format(nrow(daily_full), big.mark = ","), ncol(daily_full)))
cat(sprintf("bout_full:   %s rows x %d cols\n",
            format(nrow(bout_full), big.mark = ","), ncol(bout_full)))
cat(sprintf("\nCovariates: %s\n", paste(covariate_cols, collapse = ", ")))

log_step("=== 03_covariate_extraction.R: DONE ===")
sink(type = "message"); sink()
