# ==============================================================================
# 02_response_variables.R
# ==============================================================================
# Response variables:
#   1. Hourly intensity (per-hour movement magnitude)
#   2. Daily active minutes
#   3. Active bout length (consecutive active period)
#
# Produces (SI): Table S4 (Table_02_ActivityThresholds — per-bear GMM thresholds),
#                Table S5 (Table_03_DescriptiveStats, sheet Hourly_bySexSeason —
#                          descriptive activity by sex and season).
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("02_response_variables")
log_step("=== 02_response_variables.R: START ===")

act_all  <- readRDS(file.path(dat_path, "activity_clean.rds"))
meta     <- readRDS(file.path(dat_path, "metadata.rds"))
daily_qc <- readRDS(file.path(dat_path, "daily_qc.rds"))

# ==============================================================================
# CKPT: GMM activity thresholds (per bear)
# ==============================================================================
threshold_table <- with_ckpt("02a_thresholds", {
  log_step("GMM thresholds per bear...")
  act_all %>%
    group_by(BearID) %>%
    summarise(
      N_records     = n(),
      Mag_mean      = mean(Mag, na.rm = TRUE),
      Mag_median    = stats::median(Mag, na.rm = TRUE),
      GMM_threshold = get_gmm_threshold(Mag),
      .groups = "drop"
    )
})
log_step("Threshold range: %.1f - %.1f (median: %.1f)",
         min(threshold_table$GMM_threshold),
         max(threshold_table$GMM_threshold),
         stats::median(threshold_table$GMM_threshold))

act_all <- act_all %>%
  left_join(threshold_table %>% select(BearID, GMM_threshold), by = "BearID") %>%
  mutate(Is_Active = as.integer(Mag > GMM_threshold))

# -> Table S4 (SI): per-bear GMM activity thresholds
write.csv(threshold_table,
          file.path(tbl_path, "Table_02_ActivityThresholds.csv"),
          row.names = FALSE)

# ==============================================================================
# Sun times (per date)
# ==============================================================================
# ---- TIME ZONE ---------------------------------------------------------------
# Sun times are computed in TZ_DATA (UTC), the same zone as the activity axis
# (collar LMT == UTC — see 01_data_preparation.R). Using the local zone here
# while the activity DateTime is effectively UTC would produce a ~3-hour shift,
# mis-assigning the four t_period classes and propagating that error to the
# nocturnality index (11), the conflict windows (15) and every diel curve.
log_step("Computing sun times... (tz=%s — aligned with 01)", TZ_DATA)
sun_times <- suncalc::getSunlightTimes(
  date = unique(as.Date(act_all$Date_Only)),
  lat  = study_lat, lon = study_lon, tz = TZ_DATA,
  keep = c("sunrise", "sunset", "dawn", "dusk")
) %>%
  mutate(
    SR   = lubridate::hour(sunrise) + lubridate::minute(sunrise) / 60,
    SS   = lubridate::hour(sunset)  + lubridate::minute(sunset)  / 60,
    Dawn = lubridate::hour(dawn)    + lubridate::minute(dawn)    / 60,
    Dusk = lubridate::hour(dusk)    + lubridate::minute(dusk)    / 60,
    Date_Only = date,
    # Night length — needed as an offset in the nocturnality model: the night
    # fraction depends mechanically on night length, which shifts from ~10 to
    # ~14 h across April-November at 40N; Season_f (2 levels) cannot control it.
    night_len_h = 24 - ((SS - SR) %% 24),
    day_len_h   = (SS - SR) %% 24
  ) %>%
  select(Date_Only, SR, SS, Dawn, Dusk, night_len_h, day_len_h)

log_step("  Day-length range: %.2f - %.2f h (swing %.2f)",
         min(sun_times$day_len_h), max(sun_times$day_len_h),
         diff(range(sun_times$day_len_h)))
log_step("  Sunset range: %.2f - %.2f (swing %.2f h)",
         min(sun_times$SS), max(sun_times$SS), diff(range(sun_times$SS)))

# ==============================================================================
# 1) HOURLY INTENSITY
# ==============================================================================
log_step("Computing hourly intensity...")
hourly_intensity <- with_ckpt("02b_hourly", {
  act_all %>%
    mutate(Hour_block = Hour) %>%
    group_by(BearID, BearYear, Date_Only, Year, Month, Hour_block,
             Sex, Current_Age, Age_Category, Season) %>%
    summarise(
      # --- TWO RESPONSES ------------------------------------------------------
      # Two complementary responses are kept because they measure DIFFERENT
      # biological dimensions and dissociate at the dump:
      #   H0   (zi-Beta, intensity): dump distance +, pd=1.00 -> moves LESS nearby
      #   H0bb (beta-binom, duration): dump distance -0.020 [-0.096, 0.056] -> flat
      # => Bears near the dump are active for the SAME duration but move LESS
      #    intensely: a sit-and-feed signature.
      Activity   = sum(Mag, na.rm = TRUE),        # legacy response — for back-comparison
      Mag_mean   = mean(Mag, na.rm = TRUE),       # INTENSITY -> zi-Beta response
      Active_n   = sum(Is_Active, na.rm = TRUE),  # DURATION  -> beta-binomial numerator
      N_readings = n(),                           #           -> beta-binomial trials()
      Pct_active = mean(Is_Active, na.rm = TRUE), # descriptive (=Active_n/N_readings)
      .groups = "drop"
    ) %>%
    filter(N_readings >= 6) %>%
    mutate(
      # INTENSITY response: [0,1], no arbitrary constant. MAXMAG = 255*sqrt(3) =
      # 441.67, the accelerometer's theoretical maximum (each axis 0..255). The
      # old "Activity/5000" was arbitrary and grew mechanically with N_readings.
      # Observed max mean(Mag) = 421.4 -> 0.954: no ceiling/censoring problem.
      # Observed zero fraction ~6.6% -> zero_inflated_beta appropriate.
      Act_intensity = pmin(Mag_mean / MAXMAG, 0.9999),
      Active_n      = as.integer(Active_n),
      N_readings    = as.integer(N_readings)
    ) %>%
    left_join(sun_times, by = "Date_Only") %>%
    mutate(
      Hour_mid = Hour_block + 0.5,
      # t_period: both Hour_mid and SR/SS/Dawn/Dusk are in the same zone (UTC).
      t_period = case_when(
        Hour_mid >= Dawn & Hour_mid < SR   ~ "morning_crepuscular",
        Hour_mid >= SR   & Hour_mid < SS   ~ "diurnal",
        Hour_mid >= SS   & Hour_mid < Dusk ~ "evening_crepuscular",
        TRUE ~ "nocturnal"
      ),
      # --- SOLAR AXIS (primary diel axis) -------------------------------------
      # The wall-clock axis had two problems: (a) sunset drifts ~3 h across the
      # season, smearing crepuscular peaks; (b) s(Hour_block, bs='cc') without
      # explicit knots let mgcv place boundary knots at 0 and 23, compressing
      # the 24-hour cycle into 23 hours. On the solar axis 0 and 24 are the same
      # instant biologically, so knots = c(0,24) is natural.
      SolarHour_dbl = solar_hour(Hour_mid, SR, SS, mode = "double"),  # primary
      SolarHour_ss  = solar_hour(Hour_mid, SR, SS, mode = "sunset")   # comparator (cf. 31)
    )
})
log_step("  Hourly: %s records, %d bears",
         format(nrow(hourly_intensity), big.mark = ","),
         n_distinct(hourly_intensity$BearID))

# ==============================================================================
# 2) DAILY ACTIVE MINUTES
# ==============================================================================
log_step("Computing daily active minutes...")
daily_active <- with_ckpt("02c_daily", {
  act_all %>%
    group_by(BearID, BearYear, Date_Only, Year, Month,
             Sex, Current_Age, Age_Category, Season) %>%
    summarise(
      Active_bins      = sum(Is_Active, na.rm = TRUE),
      Total_bins       = n(),
      Active_minutes   = Active_bins * 5,
      Inactive_minutes = (Total_bins - Active_bins) * 5,
      Pct_active       = Active_bins / Total_bins,
      Mean_Mag         = mean(Mag, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(Total_bins >= 144)
})
log_step("  Daily: %s records, %d bears",
         format(nrow(daily_active), big.mark = ","),
         n_distinct(daily_active$BearID))

# Extra table split by time-period
daily_active_period <- act_all %>%
  left_join(sun_times, by = "Date_Only") %>%
  mutate(
    Hour_dec = Hour + lubridate::minute(DateTime) / 60,
    t_period = case_when(
      Hour_dec >= Dawn & Hour_dec < SR   ~ "morning_crepuscular",
      Hour_dec >= SR   & Hour_dec < SS   ~ "diurnal",
      Hour_dec >= SS   & Hour_dec < Dusk ~ "evening_crepuscular",
      TRUE ~ "nocturnal"
    )
  ) %>%
  group_by(BearID, BearYear, Date_Only, Year, Month,
           Sex, Current_Age, Age_Category, Season, t_period) %>%
  summarise(
    Active_minutes = sum(Is_Active, na.rm = TRUE) * 5,
    Total_minutes  = n() * 5,
    .groups = "drop"
  )

# ==============================================================================
# 3) ACTIVE BOUT LENGTH
# ==============================================================================
log_step("Computing active bouts...")
active_bouts <- with_ckpt("02d_bouts", {
  compute_bouts <- function(df_bd) {
    if (nrow(df_bd) == 0) return(NULL)
    df_bd <- df_bd %>% arrange(DateTime)
    av <- df_bd$Is_Active
    bid <- cumsum(c(1, diff(av) != 0))
    df_bd$bout_id <- bid
    ab <- df_bd %>%
      filter(Is_Active == 1) %>%
      group_by(bout_id) %>%
      summarise(
        Bout_start      = min(DateTime),
        Bout_end        = max(DateTime),
        Bout_length_min = n() * 5,
        Mean_Mag        = mean(Mag, na.rm = TRUE),
        .groups = "drop"
      )
    if (nrow(ab) == 0) return(NULL)
    ab %>%
      mutate(
        BearID       = df_bd$BearID[1],
        BearYear     = df_bd$BearYear[1],
        Date_Only    = df_bd$Date_Only[1],
        Year         = df_bd$Year[1],
        Month        = df_bd$Month[1],
        Sex          = df_bd$Sex[1],
        Current_Age  = df_bd$Current_Age[1],
        Age_Category = df_bd$Age_Category[1],
        Season       = df_bd$Season[1]
      )
  }

  bear_days <- act_all %>% distinct(BearID, Date_Only) %>% arrange(BearID, Date_Only)
  log_step("  bouts for %d bear-days...", nrow(bear_days))

  bout_list <- vector("list", nrow(bear_days))
  for (i in seq_len(nrow(bear_days))) {
    bd <- bear_days[i, ]
    df_sub <- act_all %>% filter(BearID == bd$BearID, Date_Only == bd$Date_Only)
    bo <- compute_bouts(df_sub)
    if (!is.null(bo)) bout_list[[i]] <- bo
    if (i %% 5000 == 0) log_step("    bout %d / %d", i, nrow(bear_days))
  }
  bind_rows(bout_list)
})

active_bouts <- active_bouts %>%
  mutate(Bout_hour = lubridate::hour(Bout_start) + lubridate::minute(Bout_start) / 60) %>%
  left_join(sun_times, by = "Date_Only") %>%
  mutate(
    t_period = case_when(
      Bout_hour >= Dawn & Bout_hour < SR   ~ "morning_crepuscular",
      Bout_hour >= SR   & Bout_hour < SS   ~ "diurnal",
      Bout_hour >= SS   & Bout_hour < Dusk ~ "evening_crepuscular",
      TRUE ~ "nocturnal"
    )
  )

log_step("  Bouts: %s records, %d bears, mean %.1f min",
         format(nrow(active_bouts), big.mark = ","),
         n_distinct(active_bouts$BearID),
         mean(active_bouts$Bout_length_min))

# ==============================================================================
# 5. DESCRIPTIVE STATS
# ==============================================================================
log_step("Descriptive statistics...")

desc_hourly <- hourly_intensity %>%
  group_by(Sex, Season) %>%
  summarise(N = n(), N_bears = n_distinct(BearID),
            Mean_activity   = mean(Activity, na.rm = TRUE),
            SD_activity     = sd(Activity, na.rm = TRUE),
            Median_activity = stats::median(Activity, na.rm = TRUE),
            .groups = "drop")

desc_daily <- daily_active %>%
  group_by(Sex, Season) %>%
  summarise(N = n(), N_bears = n_distinct(BearID),
            Mean_active_min   = mean(Active_minutes, na.rm = TRUE),
            SD_active_min     = sd(Active_minutes, na.rm = TRUE),
            Median_active_min = stats::median(Active_minutes, na.rm = TRUE),
            .groups = "drop")

desc_bouts <- active_bouts %>%
  group_by(Sex, Season) %>%
  summarise(N = n(), N_bears = n_distinct(BearID),
            Mean_bout_min   = mean(Bout_length_min, na.rm = TRUE),
            SD_bout_min     = sd(Bout_length_min, na.rm = TRUE),
            Median_bout_min = stats::median(Bout_length_min, na.rm = TRUE),
            .groups = "drop")

desc_hourly_age <- hourly_intensity %>%
  group_by(Age_Category, Sex, Season) %>%
  summarise(N = n(), N_bears = n_distinct(BearID),
            Mean_activity = mean(Activity, na.rm = TRUE),
            SD_activity   = sd(Activity, na.rm = TRUE),
            .groups = "drop")

# -> Table S5 (SI): sheet Hourly_bySexSeason = descriptive activity by sex and season
writexl::write_xlsx(list(
  Hourly_bySexSeason     = desc_hourly,
  Daily_bySexSeason      = desc_daily,
  Bouts_bySexSeason      = desc_bouts,
  Hourly_byAgeSexSeason  = desc_hourly_age
), file.path(tbl_path, "Table_03_DescriptiveStats.xlsx"))

cat("\n========================================\n")
cat("  RESPONSE VARIABLES SUMMARY\n")
cat("========================================\n")
cat("Hourly:\n");   print(desc_hourly)
cat("\nDaily:\n");  print(desc_daily)
cat("\nBouts:\n"); print(desc_bouts)

# ==============================================================================
# 6. SAVE RDS
# ==============================================================================
log_step("Saving RDS...")
save_rds_safe(hourly_intensity,    file.path(dat_path, "hourly_intensity.rds"))
save_rds_safe(daily_active,        file.path(dat_path, "daily_active.rds"))
save_rds_safe(daily_active_period, file.path(dat_path, "daily_active_period.rds"))
save_rds_safe(active_bouts,        file.path(dat_path, "active_bouts.rds"))
save_rds_safe(threshold_table,     file.path(dat_path, "gmm_thresholds.rds"))
save_rds_safe(sun_times,           file.path(dat_path, "sun_times.rds"))

log_step("=== 02_response_variables.R: DONE ===")
sink(type = "message"); sink()
