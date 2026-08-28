# ==============================================================================
# 01_data_preparation.R
# ==============================================================================
# Reads the Activity + Position files, joins them with the metadata, and applies
# a stepwise filter (age > 2, 72 h post-capture, Apr-Nov, bounding box, QC).
# Output: activity_clean.rds, position_clean.rds, bearyear_summary.rds + tables.
#
# Produces (SI): Table S1  (Table_01a_FilterSteps    — data-filtering steps),
#                Table S2  (Table_01b_FinalBearYear   — final bear-year inventory),
#                Table S3  (Table_01c_GPSIntervals    — per-bear GPS fix intervals).
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("01_data_preparation")
log_step("=== 01_data_preparation.R: START ===")

# ==============================================================================
# CKPT 1: METADATA + FILE LISTS
# ==============================================================================
ckpt1 <- with_ckpt("01a_metadata_and_files", {
  log_step("Reading metadata...")
  metadata_raw <- readxl::read_excel(meta_path)

  metadata <- metadata_raw %>%
    mutate(
      Code         = as.character(Code),
      Capture_Date = as.Date(Date, format = "%d/%m/%Y"),
      Capture_Year = lubridate::year(Capture_Date),
      Sex          = as.character(Sex),
      Capture_Age  = suppressWarnings(as.numeric(`Age (Year)`))
    ) %>%
    filter(!is.na(Capture_Age), !is.na(Capture_Date)) %>%
    select(Code, Name, Sex, Capture_Age, Capture_Year, Capture_Date,
           Mass = `Mass (kg)`, CollarID = `Collar ID`) %>%
    arrange(Code, Capture_Date)

  log_step("  Metadata: %d records, %d unique bears",
           nrow(metadata), n_distinct(metadata$Code))

  meta_first <- metadata %>%
    group_by(Code) %>%
    slice_min(Capture_Date, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(Code, Sex, Capture_Age, Capture_Year, Capture_Date)

  act_files <- list.files(act_path, pattern = "\\.txt$", ignore.case = TRUE)
  pos_files <- list.files(pos_path, pattern = "\\.txt$", ignore.case = TRUE)
  act_bears <- gsub("\\.txt$", "", act_files, ignore.case = TRUE)
  pos_bears <- gsub("\\.txt$", "", pos_files, ignore.case = TRUE)

  # BTR0093 -> BTR093 normalization
  if ("BTR0093" %in% pos_bears && !"BTR093" %in% pos_bears) {
    src <- file.path(pos_path, "BTR0093.txt")
    dst <- file.path(pos_path, "BTR093.txt")
    if (file.exists(src) && !file.exists(dst)) {
      file.copy(src, dst)
      log_step("  Position: copied BTR0093.txt -> BTR093.txt")
    }
    pos_bears <- gsub("BTR0093", "BTR093", pos_bears)
  }
  overlap_bears <- intersect(act_bears, pos_bears)
  log_step("  Activity: %d | Position: %d | Overlap: %d",
           length(act_bears), length(pos_bears), length(overlap_bears))

  list(metadata = metadata, meta_first = meta_first,
       act_bears = act_bears, pos_bears = pos_bears,
       overlap_bears = overlap_bears)
})

metadata      <- ckpt1$metadata
meta_first    <- ckpt1$meta_first
overlap_bears <- ckpt1$overlap_bears

# ==============================================================================
# CKPT 2: RAW ACTIVITY DATA
# ==============================================================================
act_all <- with_ckpt("01b_activity_raw", {
  log_step("Reading activity (%d bears)...", length(overlap_bears))
  act_list <- list()
  for (bear in overlap_bears) {
    fp <- file.path(act_path, paste0(bear, ".txt"))
    df <- read_activity_file(fp)
    if (is.null(df)) {
      log_step("  WARNING: could not read activity -> %s", bear); next
    }
    df$BearID <- bear
    act_list[[bear]] <- df
    if (which(overlap_bears == bear) %% 10 == 0) {
      log_step("  read %d / %d", which(overlap_bears == bear), length(overlap_bears))
    }
  }
  bind_rows(act_list)
})
log_step("Total raw activity: %s", format(nrow(act_all), big.mark = ","))

# ==============================================================================
# CKPT 3: RAW POSITION DATA
# ==============================================================================
pos_all <- with_ckpt("01c_position_raw", {
  log_step("Reading position (%d bears)...", length(overlap_bears))
  pos_list <- list()
  for (bear in overlap_bears) {
    fp <- file.path(pos_path, paste0(bear, ".txt"))
    if (!file.exists(fp)) {
      alt_names <- c(paste0("BTR0", gsub("BTR", "", bear)), bear)
      for (alt in alt_names) {
        fp_alt <- file.path(pos_path, paste0(alt, ".txt"))
        if (file.exists(fp_alt)) { fp <- fp_alt; break }
      }
    }
    if (!file.exists(fp)) { log_step("  WARNING: position not found -> %s", bear); next }
    df <- read_position_file(fp)
    if (is.null(df)) { log_step("  WARNING: could not read position -> %s", bear); next }
    df$BearID <- bear
    pos_list[[bear]] <- df
    if (which(overlap_bears == bear) %% 10 == 0) {
      log_step("  read %d / %d", which(overlap_bears == bear), length(overlap_bears))
    }
  }
  bind_rows(pos_list)
})
log_step("Total raw position: %s", format(nrow(pos_all), big.mark = ","))

# ==============================================================================
# 4. DateTime parse + metadata join + bear-year
# ==============================================================================
# ---- TIME ZONE: critical detail --------------------------------------------
# The raw Vectronic files store UTC and LMT columns that are IDENTICAL: the
# collars' LMT offset was never programmed, so LMT == UTC. read_activity_file
# reads columns 5-6 (LMT), so DateTime is effectively UTC. Labelling it as a
# local zone would create a ~3-hour systematic shift relative to the sun times
# (sunrise/sunset/dawn/dusk) computed in 02, mis-timing t_period, the
# nocturnality index, the conflict windows and every diel curve. We therefore
# keep everything in UTC; the sun times in 02 are also computed in UTC and the
# diel axis is a sun-relative SolarHour, so the time zone becomes irrelevant.
# Local time is used for presentation only (lubridate::with_tz).
# ------------------------------------------------------------------------------
TZ_DATA <- "UTC"   # collar LMT == UTC (verified)
log_step("DateTime parse + metadata join... (tz=%s — LMT==UTC verified)", TZ_DATA)

act_all <- act_all %>%
  mutate(
    DateTime  = lubridate::dmy_hms(paste(Date, Time), tz = TZ_DATA),
    Date_Only = as.Date(DateTime),
    Year      = lubridate::year(DateTime),
    Month     = lubridate::month(DateTime),
    Hour      = lubridate::hour(DateTime)
  ) %>%
  filter(!is.na(DateTime)) %>%
  left_join(meta_first, by = c("BearID" = "Code")) %>%
  mutate(
    Current_Age  = Capture_Age + (Year - Capture_Year),
    Age_Category = assign_age_category(Current_Age),
    Season       = assign_season(DateTime),
    BearYear     = paste(BearID, Year, sep = "_")
  )

pos_all <- pos_all %>%
  mutate(
    DateTime  = lubridate::dmy_hms(paste(Date, Time), tz = TZ_DATA),   # UTC (see above)
    Date_Only = as.Date(DateTime),
    Year      = lubridate::year(DateTime),
    Month     = lubridate::month(DateTime)
  ) %>%
  filter(!is.na(DateTime)) %>%
  left_join(meta_first, by = c("BearID" = "Code")) %>%
  mutate(
    Current_Age  = Capture_Age + (Year - Capture_Year),
    Age_Category = assign_age_category(Current_Age),
    Season       = assign_season(DateTime),
    BearYear     = paste(BearID, Year, sep = "_")
  )

# ==============================================================================
# 5. FILTERING (step by step)
# ==============================================================================
filter_tracker <- data.frame()

record_filter <- function(step_name, df_act, df_pos) {
  ac <- df_act %>% group_by(BearID) %>% summarise(N_act = n(), .groups = "drop")
  pc <- df_pos %>% group_by(BearID) %>% summarise(N_pos = n(), .groups = "drop")
  combined <- full_join(ac, pc, by = "BearID")
  combined$Step       <- step_name
  combined$Total_act  <- nrow(df_act)
  combined$Total_pos  <- nrow(df_pos)
  combined$N_bears_act <- n_distinct(df_act$BearID)
  combined$N_bears_pos <- n_distinct(df_pos$BearID)
  combined
}

filter_tracker <- bind_rows(filter_tracker, record_filter("0_Raw", act_all, pos_all))

log_step("Filter 1: has metadata...");
act_all <- act_all %>% filter(!is.na(Sex), !is.na(Capture_Age))
pos_all <- pos_all %>% filter(!is.na(Sex), !is.na(Capture_Age))
filter_tracker <- bind_rows(filter_tracker, record_filter("1_WithMetadata", act_all, pos_all))

log_step("Filter 2: age > 2...")
act_all <- act_all %>% filter(Current_Age > 2)
pos_all <- pos_all %>% filter(Current_Age > 2)
filter_tracker <- bind_rows(filter_tracker, record_filter("2_Age_gt2", act_all, pos_all))

log_step("Filter 3: drop first 72 h post-capture...")
capture_dates_list <- metadata %>%
  select(Code, Capture_Date) %>%
  filter(!is.na(Capture_Date)) %>%
  group_by(Code) %>%
  summarise(cap_dates = list(sort(unique(Capture_Date))), .groups = "drop") %>%
  tibble::deframe()

filter_post_capture <- function(df) {
  keep <- rep(TRUE, nrow(df))
  for (bear in unique(df$BearID)) {
    caps <- capture_dates_list[[bear]]
    if (is.null(caps)) next
    idx <- which(df$BearID == bear)
    bt <- df$DateTime[idx]
    for (cap_date in as.list(caps)) {
      cap_time <- as.POSIXct(paste(cap_date, "00:00:00"), tz = TZ_DATA)  # aligned with UTC
      hd <- as.numeric(difftime(bt, cap_time, units = "hours"))
      keep[idx] <- keep[idx] & !(hd >= 0 & hd < 72)
    }
  }
  df[keep, ]
}

act_all <- filter_post_capture(act_all)
pos_all <- filter_post_capture(pos_all)
filter_tracker <- bind_rows(filter_tracker, record_filter("3_Post72h", act_all, pos_all))

log_step("Filter 4: active season (Apr-Nov)...")
act_all <- act_all %>% filter(Month >= 4 & Month <= 11)
pos_all <- pos_all %>% filter(Month >= 4 & Month <= 11)
filter_tracker <- bind_rows(filter_tracker, record_filter("4_AprNov", act_all, pos_all))

log_step("Filter 5: bounding box (lat:%.1f-%.1f, lon:%.1f-%.1f)...",
         bbox_lat[1], bbox_lat[2], bbox_lon[1], bbox_lon[2])
n_before <- nrow(pos_all)
pos_all <- pos_all %>%
  filter(!is.na(Latitude), !is.na(Longitude),
         Latitude  >= bbox_lat[1], Latitude  <= bbox_lat[2],
         Longitude >= bbox_lon[1], Longitude <= bbox_lon[2])
log_step("  bbox: dropped %d records (%.2f%%)",
         n_before - nrow(pos_all), (n_before - nrow(pos_all)) / n_before * 100)
filter_tracker <- bind_rows(filter_tracker, record_filter("5_BoundingBox", act_all, pos_all))

log_step("Filter 6: activity QC (origin, saturation, daily completeness)...")
act_all <- act_all %>% filter(Origin == "Collar" | is.na(Origin))
act_all <- act_all %>%
  mutate(
    ActivityX = tidyr::replace_na(ActivityX, 0),
    ActivityY = tidyr::replace_na(ActivityY, 0),
    ActivityZ = tidyr::replace_na(ActivityZ, 0),
    Mag = sqrt(ActivityX^2 + ActivityY^2 + ActivityZ^2)
  )
daily_qc <- act_all %>%
  group_by(BearID, Date_Only) %>%
  summarise(
    n_rec = n(),
    pct_saturated = sum(ActivityX >= 255 | ActivityY >= 255 | ActivityZ >= 255,
                        na.rm = TRUE) / n(),
    .groups = "drop"
  ) %>%
  filter(n_rec >= 50, pct_saturated < 0.10)
act_all <- act_all %>% semi_join(daily_qc, by = c("BearID", "Date_Only"))
filter_tracker <- bind_rows(filter_tracker, record_filter("6_ActivityQC", act_all, pos_all))

log_step("Filter 7: GPS validity (DOP)...")
pos_all <- pos_all %>%
  filter(!is.na(Latitude), !is.na(Longitude),
         Origin == "Collar" | is.na(Origin))
if ("DOP" %in% names(pos_all)) {
  pos_all <- pos_all %>% filter(is.na(DOP) | DOP <= 20)
}
filter_tracker <- bind_rows(filter_tracker, record_filter("7_ValidGPS", act_all, pos_all))

# Common bears
common_bears <- intersect(unique(act_all$BearID), unique(pos_all$BearID))
act_all <- act_all %>% filter(BearID %in% common_bears)
pos_all <- pos_all %>% filter(BearID %in% common_bears)
filter_tracker <- bind_rows(filter_tracker, record_filter("8_Final", act_all, pos_all))

# ==============================================================================
# 6. Bear-year summary + GPS interval
# ==============================================================================
log_step("Bear-year summary...")
bearyear_summary <- act_all %>%
  group_by(BearID, Year, Sex, Current_Age, Age_Category, BearYear) %>%
  summarise(
    N_activity_records = n(),
    N_active_days      = n_distinct(Date_Only),
    Date_start         = min(Date_Only),
    Date_end           = max(Date_Only),
    Seasons            = paste(unique(Season), collapse = ", "),
    .groups = "drop"
  ) %>%
  left_join(
    pos_all %>% group_by(BearID, Year) %>%
      summarise(N_positions = n(),
                N_pos_days = n_distinct(Date_Only), .groups = "drop"),
    by = c("BearID", "Year")
  ) %>%
  arrange(BearID, Year)

log_step("Final: %d bears, %d bear-years",
         n_distinct(bearyear_summary$BearID), nrow(bearyear_summary))

log_step("GPS fix-interval analysis...")
gps_interval_summary <- pos_all %>%
  arrange(BearID, DateTime) %>%
  group_by(BearID) %>%
  mutate(dt_min = as.numeric(difftime(DateTime, lag(DateTime), units = "mins"))) %>%
  ungroup() %>%
  filter(!is.na(dt_min), dt_min > 0, dt_min < 1440) %>%
  group_by(BearID) %>%
  summarise(
    median_interval_min = stats::median(dt_min, na.rm = TRUE),
    mean_interval_min   = mean(dt_min, na.rm = TRUE),
    q25_interval        = stats::quantile(dt_min, 0.25, na.rm = TRUE),
    q75_interval        = stats::quantile(dt_min, 0.75, na.rm = TRUE),
    N_fixes = n(),
    .groups = "drop"
  )

# ==============================================================================
# 7. WRITE TABLES
# ==============================================================================
log_step("Writing tables...")

filter_steps_summary <- filter_tracker %>%
  group_by(Step) %>%
  summarise(
    N_bears_activity      = first(N_bears_act),
    N_bears_position      = first(N_bears_pos),
    Total_activity_records = first(Total_act),
    Total_position_records = first(Total_pos),
    .groups = "drop"
  ) %>%
  arrange(Step)

bear_filter_wide <- filter_tracker %>%
  select(BearID, Step, N_act) %>%
  tidyr::pivot_wider(names_from = Step, values_from = N_act, values_fill = 0) %>%
  left_join(meta_first %>% select(Code, Sex, Capture_Age),
            by = c("BearID" = "Code")) %>%
  arrange(BearID)

final_table <- bearyear_summary %>%
  mutate(Duration_days = as.numeric(Date_end - Date_start)) %>%
  select(BearID, Year, Sex, Age = Current_Age, Age_Category,
         N_activity_records, N_active_days, N_positions,
         Date_start, Date_end, Duration_days, Seasons)

writexl::write_xlsx(list(
  FilterSteps      = filter_steps_summary,
  BearFilterDetail = bear_filter_wide,
  FinalBearYear    = final_table,
  GPSIntervals     = gps_interval_summary
), file.path(tbl_path, "Table_01_DataFiltering.xlsx"))

# -> Table S1 (SI): data-filtering steps (bears/records retained per QC stage)
write.csv(filter_steps_summary, file.path(tbl_path, "Table_01a_FilterSteps.csv"), row.names = FALSE)
# -> Table S2 (SI): final analytic bear-year inventory
write.csv(final_table,          file.path(tbl_path, "Table_01b_FinalBearYear.csv"), row.names = FALSE)
# -> Table S3 (SI): per-bear GPS fix-interval summary
write.csv(gps_interval_summary, file.path(tbl_path, "Table_01c_GPSIntervals.csv"), row.names = FALSE)

cat("\n========================================\n")
cat("  DATA PREPARATION SUMMARY\n")
cat("========================================\n")
print(filter_steps_summary)
cat(sprintf("\nFinal: %d bears, %d bear-years, F=%d M=%d\n",
            n_distinct(act_all$BearID), nrow(bearyear_summary),
            sum(final_table$Sex == "F"), sum(final_table$Sex == "M")))

# ==============================================================================
# 8. SAVE RDS
# ==============================================================================
log_step("Saving RDS...")
save_rds_safe(act_all,              file.path(dat_path, "activity_clean.rds"))
save_rds_safe(pos_all,              file.path(dat_path, "position_clean.rds"))
save_rds_safe(meta_first,           file.path(dat_path, "metadata.rds"))
save_rds_safe(bearyear_summary,     file.path(dat_path, "bearyear_summary.rds"))
save_rds_safe(gps_interval_summary, file.path(dat_path, "gps_intervals.rds"))
save_rds_safe(daily_qc,             file.path(dat_path, "daily_qc.rds"))

log_step("=== 01_data_preparation.R: DONE ===")
sink(type = "message"); sink()
