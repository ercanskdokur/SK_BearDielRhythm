# ==============================================================================
# 11_nocturnality_index.R  —  Nocturnality index + human pressure
# ==============================================================================
# FRAMING: human presence is known to increase nocturnality in mammals. A
#   direct, interpretable response: the fraction of per-bear-day activity that
#   falls at NIGHT.
#
# This is the SCALAR, easy-to-communicate complement of H3 (hour x human-pressure
# smooth), and it rests on solar time (t_period: derived from SR/SS/Dawn/Dusk) —
# resolving the clock-hour confounder up front.
#
# RESPONSE (beta-binomial; NO arbitrary scale):
#   Night_active | trials(Total_active)
#   Night_active  = count of ACTIVE readings in night+twilight hours
#   Total_active  = count of active readings across all hours (per bear-day)
#
# MODEL:
#   Night_active | trials(Total_active) ~
#     d2GarbageDump_km_sc + d2Roads_km_sc + d2ProtectedAreas_km_sc +
#     LandPC1_sc + LandPC2_sc + temp_mean_sc + Season_f + Sex_f + Age_sc +
#     (1 | BearID_f/Year_f)
#
# OUTPUTS:
#   output/data/nocturnality_daily.rds
#   output/models/Noct_betabinom.rds
#   output/tables/Table_V3_08_Nocturnality_Coef.csv
#     -> Table S28 (strict night)  /  Table S30 (wide night, sensitivity)
#   output/figures/FigV3_02_Nocturnality_vs_Dump.{png,pdf}
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("11_nocturnality_index")
log_step("=== 11_nocturnality_index: START ===")

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE)
options(brms.backend = "rstan")

# ---- DATA: hourly_decorr (t_period + Pct_active + N_readings + covariates) --
meta <- readRDS(file.path(mod_path, "decorr_meta.rds"))
decorr_covs <- meta$decorr_covs
d <- readRDS(file.path(dat_path, "hourly_decorr.rds"))
d <- dplyr::filter(d, Season %in% c("Mating", "Hyperphagia"))
need <- c("t_period", "Active_n", "N_readings", "BearID_f", "Year_f",
          "Season_f", "Date_Only", "night_len_h")
miss <- setdiff(need, names(d))
if (length(miss) > 0) stop("Missing column(s): ", paste(miss, collapse = ", "))

# Hourly active-reading count
d$active_reads <- d$Active_n   # produced directly in 02 (no round(Pct*N) needed)

# ---- NIGHT DEFINITION --------------------------------------------------------
NIGHT_DEF <- Sys.getenv("NIGHT_DEF", unset = "strict")
if (!NIGHT_DEF %in% c("strict", "wide")) stop("NIGHT_DEF: strict|wide")
d$is_night <- if (NIGHT_DEF == "strict") {
  d$t_period == "nocturnal"
} else {
  d$t_period %in% c("nocturnal", "evening_crepuscular", "morning_crepuscular")
}
log_step("NIGHT DEFINITION: %s (%s)", NIGHT_DEF,
         ifelse(NIGHT_DEF == "strict", "nocturnal only (strict night definition)",
                "night + twilight — OLD definition, sensitivity"))
log_step("t_period distribution: %s",
         paste(names(table(d$t_period)), as.integer(table(d$t_period)),
               sep = "=", collapse = " | "))
log_step("is_night fraction: %.3f", mean(d$is_night))


# ---- AGGREGATE TO BEAR-DAY --------------------------------------------------
cov_sc <- intersect(c("d2GarbageDump_km_sc","d2Roads_km_sc","d2ProtectedAreas_km_sc",
                      "LandPC1_sc","LandPC2_sc","temp_mean_sc","precipitation_sc",
                      "Elevation_sc","Slope_sc","Age_sc"), names(d))
noct <- d %>%
  dplyr::group_by(BearID_f, Year_f, Season_f, Date_Only) %>%
  dplyr::summarise(
    Total_active = sum(active_reads, na.rm = TRUE),
    Night_active = sum(active_reads[is_night], na.rm = TRUE),
    Sex_f = dplyr::first(Sex_f),
    # night length is constant within a bear-day -> first() suffices. NEEDED for the offset.
    night_len_h = dplyr::first(night_len_h),
    dplyr::across(dplyr::all_of(cov_sc), ~ mean(.x, na.rm = TRUE)),
    n_hours = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(Total_active >= 5, n_hours >= 12)   # enough activity + day coverage
noct$DOY <- as.integer(lubridate::yday(noct$Date_Only))
noct$Total_active <- as.integer(noct$Total_active)
noct$Night_active <- as.integer(pmin(noct$Night_active, noct$Total_active))
log_step("Nocturnality data: %d bear-days, %d bears, mean night-fraction=%.3f",
         nrow(noct), dplyr::n_distinct(noct$BearID_f),
         mean(noct$Night_active / noct$Total_active))
save_rds_safe(noct, file.path(dat_path, "nocturnality_daily.rds"))

# ---- MODEL (beta-binomial) -------------------------------------------------
# ---- TWO FIXES --------------------------------------------------------------
# NIGHT LENGTH was UNCONTROLLED. The night FRACTION mechanically depends on
#     night LENGTH; at 40N it varies from ~10 to ~14 h across Apr->Nov. Because
#     dump/road proximity correlates with season (bears use the dump in certain
#     periods), leaving it uncontrolled would let the dump coefficient ABSORB the
#     photoperiod. FIX: log(night length) offset. Expected night-fraction is
#     proportional to night length -> on the logit scale the offset is the right
#     control.
# temp_mean_sc -> temp_hourly_sc is NOT applied (this is a bear-day aggregation;
#     the daily mean is already the right scale here) — daily temperature is
#     APPROPRIATE for this model.
fixed_terms <- paste(intersect(
  c("d2GarbageDump_km_sc","d2Roads_km_sc","d2ProtectedAreas_km_sc",
    "LandPC1_sc","LandPC2_sc","temp_mean_sc","Sex_f","Age_sc","Season_f"),
  names(noct)), collapse = " + ")

extra <- c()
if ("night_len_h" %in% names(noct)) {
  extra <- c(extra, "offset(log(night_len_h))")
  log_step("Night-length offset ADDED (range %.2f - %.2f h)",
           min(noct$night_len_h), max(noct$night_len_h))
} else {
  log_step("!! night_len_h missing — photoperiod left uncontrolled (re-run 02)")
}
if ("DOY" %in% names(noct)) {
  extra <- c(extra, "s(DOY, bs = 'tp', k = 5)")
  log_step("s(DOY, bs='tp', k=5) ADDED (was computed but unused before)")
}

fml <- brms::bf(as.formula(paste(
  "Night_active | trials(Total_active) ~", fixed_terms,
  if (length(extra)) paste("+", paste(extra, collapse = " + ")) else "",
  "+ (1 | BearID_f/Year_f)")))
log_step("Noct FORMULA: %s", paste(deparse(fml$formula), collapse = " "))

priors_bb <- c(
  brms::prior(normal(0, 1),         class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(gamma(0.01, 0.01),    class = "phi")
)

t0 <- Sys.time()
.noct_fp <- file.path(mod_path, sprintf("Noct_betabinom%s.rds",
                                        if (NIGHT_DEF == "strict") "" else "_wide"))
if (file.exists(.noct_fp) && !nzchar(Sys.getenv("FORCE_REFIT"))) {
  log_step("Reusing saved nocturnality fit (no refit; force with FORCE_REFIT): %s",
           basename(.noct_fp))
  .o <- readRDS(.noct_fp); fit_noct <- .o$fit
  if (!is.null(.o$data)) noct <- .o$data
} else {
  fit_noct <- brms::brm(
    formula = fml, data = noct, family = brms::beta_binomial(),
    prior = priors_bb,
    iter = 3000, warmup = 1000, chains = 4, cores = 4,
    seed = 42, silent = 1, refresh = 250, init = 0,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    save_pars = brms::save_pars(all = TRUE)
  )
  log_step("Noct fit time: %.1f min", as.numeric(difftime(Sys.time(), t0, "mins")))
}

s <- summary(fit_noct); fe <- s$fixed
cat("\n=== NOCTURNALITY FIXED EFFECTS ===\n"); print(fe)
coef_df <- data.frame(
  Parameter = rownames(fe),
  Estimate = round(fe[,"Estimate"],4), Q2.5 = round(fe[,"l-95% CI"],4),
  Q97.5 = round(fe[,"u-95% CI"],4), Rhat = round(fe[,"Rhat"],4),
  ESS_bulk = round(fe[,"Bulk_ESS"]),
  Direction_sig = (sign(fe[,"l-95% CI"]) == sign(fe[,"u-95% CI"]))
)
# strict (primary) and wide (sensitivity) must NOT overwrite each other
# -> Table S28 (strict) / Table S30 (wide): nocturnality model, linear dump
sfx <- if (NIGHT_DEF == "strict") "" else "_wide"
write.csv(coef_df, file.path(tbl_path, sprintf("Table_V3_08_Nocturnality_Coef%s.csv", sfx)),
          row.names = FALSE)
save_rds_safe(list(fit = fit_noct, family = "beta_binomial",
                   night_def = NIGHT_DEF, data = noct),
              file.path(mod_path, sprintf("Noct_betabinom%s.rds", sfx)))

# ---- FIGURE: nocturnality ~ dump distance ----------------------------------
ce <- tryCatch(brms::conditional_effects(fit_noct, effects = "d2GarbageDump_km_sc",
                                         method = "posterior_epred"),
               error = function(e) NULL)
if (!is.null(ce)) {
  dpl <- ce[["d2GarbageDump_km_sc"]]
  # epred returns a proportion (trials=1 reference) — interpret on a 0-1 scale
  p <- ggplot2::ggplot(dpl, ggplot2::aes(x = d2GarbageDump_km_sc, y = estimate__)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower__, ymax = upper__),
                         alpha = 0.2, fill = pal$hyperphagia) +
    ggplot2::geom_line(color = pal$hyperphagia, linewidth = 1) +
    ggplot2::labs(x = "Distance to garbage dump (z-score)",
                  y = "Night activity proportion (nocturnality)",
                  title = "Nocturnality ~ distance to garbage dump",
                  subtitle = "Positive slope = Nocturnality rises with distance") +
    theme_bear() +
    ggplot2::theme(plot.title    = ggplot2::element_text(size = 12),
                   plot.subtitle = ggplot2::element_text(size = 9.5, colour = "grey40"))
  save_fig(p, sprintf("FigV3_02_Nocturnality_vs_Dump%s", sfx), w = 8, h = 5)
}

log_step("=== 11_nocturnality_index DONE ===")
sink(type = "message"); sink()
