# ==============================================================================
# 13_nocturnality_moon.R  —  MOON-ILLUMINATION control for the nocturnality model
# ==============================================================================
# WHY: a standard confounder for nocturnality is moonlight; full-moon nights are
#   bright and animals may adjust night activity. If the human-pressure (dump/road)
#   coefficients SURVIVE controlling for moonlight, the finding is not a moonlight
#   artefact.
#
# APPROACH: to the bear-day data produced by 11_nocturnality_index.R
#   (nocturnality_daily.rds) add a moon-illumination fraction
#   (suncalc::getMoonIllumination, 0=new moon, 1=full moon); refit the
#   beta-binomial; compare dump/road coefficients with vs without moon.
#
# OUTPUTS: output/tables/Table_V3_10_Nocturnality_Moon_Coef.csv   -> Table S32 (SI)
#          output/tables/Table_V3_11_Moon_DumpRoad_Compare.csv    -> Table S33 (SI)
#          output/models/Noct_moon_betabinom.rds
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("13_nocturnality_moon")
log_step("=== 13_nocturnality_moon: START ===")

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE)
options(brms.backend = "rstan")

noct <- readRDS(file.path(dat_path, "nocturnality_daily.rds"))
log_step("Nocturnality data: %d bear-days", nrow(noct))

# ---- MOON ILLUMINATION FRACTION (per day; location-independent) ------------
noct$Date_Only <- as.Date(noct$Date_Only)
mi <- suncalc::getMoonIllumination(date = noct$Date_Only)
noct$moon_fraction <- mi$fraction
noct$moon_sc <- as.numeric(scale(noct$moon_fraction))
log_step("Moon illumination: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]",
         mean(noct$moon_fraction), sd(noct$moon_fraction),
         min(noct$moon_fraction), max(noct$moon_fraction))

# ---- MODEL: original + moon_sc ---------------------------------------------
fixed_terms <- paste(intersect(
  c("d2GarbageDump_km_sc","d2Roads_km_sc","d2ProtectedAreas_km_sc",
    "LandPC1_sc","LandPC2_sc","temp_mean_sc","Sex_f","Age_sc","Season_f",
    "moon_sc"),
  names(noct)), collapse = " + ")

extra <- c()
if ("night_len_h" %in% names(noct)) {
  extra <- c(extra, "offset(log(night_len_h))")
  log_step("Night-length offset ADDED (aligned with 11) — range %.2f - %.2f h",
           min(noct$night_len_h), max(noct$night_len_h))
} else {
  log_step("!! night_len_h missing — cannot align with 11 (regenerate nocturnality_daily)")
}
if ("DOY" %in% names(noct)) {
  extra <- c(extra, "s(DOY, bs = 'tp', k = 5)")
  log_step("s(DOY, bs='tp', k=5) ADDED (aligned with 11)")
}

fml <- brms::bf(as.formula(paste(
  "Night_active | trials(Total_active) ~", fixed_terms,
  if (length(extra)) paste("+", paste(extra, collapse = " + ")) else "",
  "+ (1 | BearID_f/Year_f)")))
log_step("Noct+Moon FORMULA (aligned with 11 +moon): %s",
         paste(deparse(fml$formula), collapse = " "))

priors_bb <- c(
  brms::prior(normal(0, 1),         class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(gamma(0.01, 0.01),    class = "phi")
)

t0 <- Sys.time()
fit <- brms::brm(fml, data = noct, family = brms::beta_binomial(),
                 prior = priors_bb, iter = 3000, warmup = 1000,
                 chains = 4, cores = 4, seed = 42, silent = 1, refresh = 250,
                 init = 0, control = list(adapt_delta = 0.95, max_treedepth = 12),
                 save_pars = brms::save_pars(all = TRUE))
log_step("Noct+Moon fit time: %.1f min", as.numeric(difftime(Sys.time(), t0, "mins")))

s <- summary(fit); fe <- s$fixed
coef_df <- data.frame(
  Parameter = rownames(fe),
  Estimate = round(fe[,"Estimate"],4), Q2.5 = round(fe[,"l-95% CI"],4),
  Q97.5 = round(fe[,"u-95% CI"],4), Rhat = round(fe[,"Rhat"],4),
  ESS_bulk = round(fe[,"Bulk_ESS"]),
  Direction_sig = (sign(fe[,"l-95% CI"]) == sign(fe[,"u-95% CI"])))
# -> Table S32 (SI): nocturnality + lunar-illumination model
write.csv(coef_df, file.path(tbl_path, "Table_V3_10_Nocturnality_Moon_Coef.csv"),
          row.names = FALSE)
cat("\n=== NOCT+MOON FIXED EFFECTS ===\n"); print(coef_df, row.names = FALSE)
save_rds_safe(list(fit = fit, family = "beta_binomial", data = noct),
              file.path(mod_path, "Noct_moon_betabinom.rds"))

# ---- COMPARISON: with-moon vs without-moon (original model) ----------------
orig_fp <- file.path(tbl_path, "Table_V3_08_Nocturnality_Coef.csv")
if (file.exists(orig_fp)) {
  o <- utils::read.csv(orig_fp)
  pick <- c("d2GarbageDump_km_sc","d2Roads_km_sc","d2ProtectedAreas_km_sc")
  cmp <- merge(
    o[o$Parameter %in% pick, c("Parameter","Estimate","Q2.5","Q97.5")],
    coef_df[coef_df$Parameter %in% pick, c("Parameter","Estimate","Q2.5","Q97.5")],
    by = "Parameter", suffixes = c("_noMoon","_withMoon"))
  # -> Table S33 (SI): moon control, dump/road coefficient comparison
  write.csv(cmp, file.path(tbl_path, "Table_V3_11_Moon_DumpRoad_Compare.csv"),
            row.names = FALSE)
  cat("\n=== DUMP/ROAD: without-moon vs with-moon ===\n"); print(cmp, row.names = FALSE)
}

log_step("=== 13_nocturnality_moon DONE ===")
cat("DONE 13\n")
sink(type = "message"); sink()
