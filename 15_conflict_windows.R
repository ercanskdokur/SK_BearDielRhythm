# ==============================================================================
# 15_conflict_windows.R  —  Human-bear conflict-risk windows (v2)
# ==============================================================================
# GOAL (applied): conflict risk is high in the twilight windows where both humans
#   and bears are active. METRIC = the share of activity that falls in the
#   conflict window (timing-normalised; cleaned of the overall activity level),
#   morning/evening SEPARATE, STRATIFIED by human-pressure (dump/road) proximity.
#
# WINDOWS (solar time): morning = SR-1..SR+4 ; evening = SS-4..SS+1
#
# OUTPUTS: Table_V3_12_ConflictWindows.csv (share metric: disturbance x prox x window)
#            -> Table S34 (SI)
#          Table_V3_13_ConflictPerBear.csv (per-individual; for phenotypes)
#            -> Table S35 (SI)
#          FigV3_06_ConflictRisk.{png,pdf}
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("15_conflict_windows")
log_step("=== 15_conflict_windows (v2): START ===")

d <- readRDS(file.path(dat_path, "hourly_decorr.rds"))
d <- dplyr::filter(d, Season %in% c("Mating", "Hyperphagia"))
need <- c("SR","SS","Hour_mid","Pct_active","BearID_f","d2GarbageDump_km","d2Roads_km")
miss <- setdiff(need, names(d)); if (length(miss)>0) stop("Missing: ", paste(miss, collapse=", "))

# The tertile cuts come from the FIXED reference in model_data (computed once from
# all data in 05_prep) -> "Near" here means the SAME distance as in the
# models/figures. The old 15 cut its own quantiles and shifted the "Near"
# definition.
brk <- tryCatch(readRDS(file.path(mod_path, "model_data_50k.rds"))$breaks,
                error = function(e) NULL)
if (is.null(brk)) log_step("!! model_data breaks missing — falling back to local quantiles")

# ---- Window flags -----------------------------------------------------------
d$in_morning  <- d$Hour_mid >= (d$SR - 1) & d$Hour_mid < (d$SR + 4)
d$in_evening  <- d$Hour_mid >= (d$SS - 4) & d$Hour_mid < (d$SS + 1)
d$in_conflict <- d$in_morning | d$in_evening
log_step("Window observation share: conflict=%.3f (morning=%.3f, evening=%.3f)",
         mean(d$in_conflict), mean(d$in_morning), mean(d$in_evening))

make_tertile <- function(x, labels, fixed = NULL) {
  qs <- if (!is.null(fixed)) fixed else stats::quantile(x, c(0,1/3,2/3,1), na.rm = TRUE)
  qs[1] <- -Inf; qs[length(qs)] <- Inf
  cut(x, breaks = unique(qs), labels = labels, include.lowest = TRUE)
}
d$DumpProx <- make_tertile(d$d2GarbageDump_km, c("Near","Mid","Far"), brk$dump)
d$RoadProx <- make_tertile(d$d2Roads_km,       c("Near","Mid","Far"), brk$road)

# ---- (1) TIMING-NORMALISED SHARE: activity share falling in the window ------
# Per disturbance x proximity group: sum(activity_window)/sum(activity_total)
frac_by <- function(prox_col, dist_label) {
  df <- d[!is.na(d[[prox_col]]), ]
  g <- split(seq_len(nrow(df)), df[[prox_col]])
  do.call(rbind, lapply(names(g), function(lv) {
    ix <- g[[lv]]; tot <- sum(df$Pct_active[ix], na.rm = TRUE)
    data.frame(Disturbance = dist_label, Proximity = lv,
               Morning = sum(df$Pct_active[ix][df$in_morning[ix]], na.rm=TRUE)/tot,
               Evening = sum(df$Pct_active[ix][df$in_evening[ix]], na.rm=TRUE)/tot,
               Conflict_total = sum(df$Pct_active[ix][df$in_conflict[ix]], na.rm=TRUE)/tot)
  }))
}
tab <- rbind(frac_by("DumpProx","Dump"), frac_by("RoadProx","Road"))
tab[,3:5] <- round(tab[,3:5], 4)
# The UNIFORM-activity NULL is added to the table. Windows span 10/24 h (nominal
# 0.417); the empirical null is the share of observation-hours falling in the
# window. This lets one read "Near below chance (avoidance), Far at chance, Mid
# above"; it also shows the raw shares are NON-monotonic.
base_m <- mean(d$in_morning); base_e <- mean(d$in_evening); base_c <- mean(d$in_conflict)
tab$Null_morning  <- round(base_m, 4)
tab$Null_evening  <- round(base_e, 4)
tab$Null_conflict <- round(base_c, 4)
log_step("NULL (uniform activity): morning=%.4f evening=%.4f conflict=%.4f | nominal 10/24=%.4f",
         base_m, base_e, base_c, 10/24)
# -> Table S34 (SI): conflict-window activity share by proximity
write.csv(tab, file.path(tbl_path, "Table_V3_12_ConflictWindows.csv"), row.names = FALSE)
cat("\n=== CONFLICT-WINDOW ACTIVITY FRACTION (timing-normalized; +NULL) ===\n"); print(tab, row.names = FALSE)

# ---- (2) PER-BEAR metrics (for phenotypes; unchanged) ----------------------
perbear <- d %>% dplyr::group_by(BearID_f) %>%
  dplyr::summarise(
    conflict_activity = mean(Pct_active[in_conflict], na.rm = TRUE),
    overall_activity  = mean(Pct_active, na.rm = TRUE),
    conflict_overlap  = sum(Pct_active[in_conflict], na.rm=TRUE)/sum(Pct_active, na.rm=TRUE),
    evening_overlap   = sum(Pct_active[in_evening],  na.rm=TRUE)/sum(Pct_active, na.rm=TRUE),
    n_obs = dplyr::n(), .groups = "drop")
# -> Table S35 (SI): per-bear conflict overlap
write.csv(perbear, file.path(tbl_path, "Table_V3_13_ConflictPerBear.csv"), row.names = FALSE)
log_step("Per-bear conflict metrics: %d bears", nrow(perbear))

# ---- (3) FIGURE: morning/evening share x proximity, dump/road facet --------
long <- rbind(
  data.frame(Disturbance=tab$Disturbance, Proximity=tab$Proximity,
             Window="Morning", Fraction=tab$Morning),
  data.frame(Disturbance=tab$Disturbance, Proximity=tab$Proximity,
             Window="Evening", Fraction=tab$Evening))
long$Proximity <- factor(long$Proximity, levels = c("Near","Mid","Far"))
long$Window <- factor(long$Window, levels = c("Morning","Evening"))
# Baseline: window-hour share (expected if activity were uniform)
base_m <- mean(d$in_morning); base_e <- mean(d$in_evening)
hl <- data.frame(Window=c("Morning","Evening"), y=c(base_m, base_e))

p <- ggplot2::ggplot(long, ggplot2::aes(x = Proximity, y = Fraction,
                                        color = Window, group = Window)) +
  ggplot2::geom_hline(data = hl, ggplot2::aes(yintercept = y, color = Window),
                      linetype = "dashed", alpha = 0.6) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::facet_wrap(~ Disturbance) +
  ggplot2::scale_color_manual(values = c(Morning = pal$mating, Evening = pal$hyperphagia)) +
  ggplot2::labs(x = "Proximity to disturbance",
                y = "Share of daily activity in window",
                color = "Conflict window",
                title = "Human-bear conflict-window activity share",
                subtitle = "Timing-normalized; dashed = share expected if activity were uniform") +
  theme_bear() +
  ggplot2::theme(plot.title    = ggplot2::element_text(size = 12),
                 plot.subtitle = ggplot2::element_text(size = 9.5, colour = "grey40"))
save_fig(p, "FigV3_06_ConflictRisk", w = 8.5, h = 4.8)

log_step("=== 15_conflict_windows (v2) DONE ===")
cat("DONE 15\n")
sink(type = "message"); sink()
