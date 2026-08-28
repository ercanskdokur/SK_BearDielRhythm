# ==============================================================================
# 05_fit_models.R  —  H0..H5 x {intensity, duration}, single parametric script
# ==============================================================================
# Call:  HYP=H3 RESP=intensity Rscript 05_fit_models.R
# SLURM: --array=0-11  (6 hypotheses x 2 responses; see HYP_GRID below)
# Output: models/<HYP>_<RESP>.rds  (models only; no SI table)
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))

# ---- PARAMETERS --------------------------------------------------------------
HYP_GRID <- expand.grid(hyp = c("H0","H1","H2","H3","H4","H5"),
                        resp = c("intensity","duration"),
                        stringsAsFactors = FALSE)

aid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
if (nzchar(aid)) {
  i <- as.integer(aid) + 1L
  if (i < 1 || i > nrow(HYP_GRID)) stop("array id out of range: ", aid)
  HYP  <- HYP_GRID$hyp[i]
  Sys.setenv(RESP = HYP_GRID$resp[i])
} else {
  HYP <- Sys.getenv("HYP", unset = "H0")
}
RS <- get_response_spec()
if (!HYP %in% HYP_GRID$hyp) stop("invalid HYP: ", HYP)

init_log(sprintf("05_fit_%s_%s_%s", HYP, RS$key,
                 Sys.getenv("DUMPSPEC", unset = "raw")))
log_step("=== 05_fit_models.R | HYP=%s | RESP=%s (%s) | DUMPSPEC=%s ===",
         HYP, RS$key, RS$family_name, Sys.getenv("DUMPSPEC", unset = "raw"))

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE)
options(brms.backend = "rstan")
dir.create(file.path(tmp_dir, "stan_cache"), recursive = TRUE, showWarnings = FALSE)

# ---- SHARED DATA (05_prep_modeldata.R) ---------------------------------------
# CRITICAL: both responses use the SAME 50K rows. A different subsample could
# artificially create the "dump does not change duration but lowers intensity"
# contrast.
md_fp <- file.path(mod_path, "model_data_50k.rds")
if (!file.exists(md_fp)) stop("model_data_50k.rds missing — run 05_prep_modeldata.R first")
md <- readRDS(md_fp)
d_mod       <- md$data
decorr_covs <- md$decorr_covs

# ---- DUMPSPEC: collinearity sensitivity of the dump covariate ----------------
# FINDING (04_decorrelate): d2GarbageDump_km_sc VIF=12.31 (PROBLEMATIC),
#   R2_on_others = 0.9150 -> the focal variable is 91.5% explained by the others.
#   Dominant co-var: LandPC1 (loadings 0.43-0.45, all equal positive = a
#   "remoteness from everything" axis). Ecological reason: the dump is ~6 km west
#   of the town.
#
# DECISION: fit BOTH specifications and compare by LOO.
#   DUMPSPEC=raw   -> d2GarbageDump_km_sc raw ("per 1 km from the dump")
#   DUMPSPEC=resid -> residualized against LandPC1+LandPC2, then rescaled
#                     ("dump proximity BEYOND what general remoteness explains")
#
# CRITICAL: the formula, response and ROWS stay identical; only the CONTENT of
#   the dump COLUMN changes -> elpd is DIRECTLY comparable (same response, same n).
DUMPSPEC <- Sys.getenv("DUMPSPEC", unset = "raw")
if (!DUMPSPEC %in% c("raw", "resid")) stop("invalid DUMPSPEC: ", DUMPSPEC)

if (identical(DUMPSPEC, "resid")) {
  stopifnot(all(c("d2GarbageDump_km_sc", "LandPC1_sc", "LandPC2_sc") %in% names(d_mod)))
  aux <- stats::lm(d2GarbageDump_km_sc ~ LandPC1_sc + LandPC2_sc, data = d_mod)
  r2_removed <- summary(aux)$r.squared
  d_mod$d2GarbageDump_km_sc <- as.numeric(scale(stats::residuals(aux)))
  log_step("DUMPSPEC=resid | dump ~ LandPC1+LandPC2 R2=%.4f -> residualize + rescale",
           r2_removed)
  log_step("  NOTE: the coefficient is NOT on the raw km scale. Interpretation: dump")
  log_step("  proximity beyond general remoteness. Over-adjustment risk is reported.")
} else {
  log_step("DUMPSPEC=raw | dump covariate RAW (VIF~12.3, collinearity reported explicitly)")
}

# The diel axis is read from the TIME_AXIS env, NOT from md. The data carries all
# three axis columns (SolarHour_dbl / SolarHour_ss / Hour_block, produced in 02),
# so no re-prep is needed for the axis sensitivity:
#   TIME_AXIS=solar_double (default, primary)
#   TIME_AXIS=clock        (sensitivity: wall-clock — the old design's axis)
#   TIME_AXIS=solar_sunset (sensitivity: sunset-anchored only)
TIME_VAR <- get_time_var()
if (!TIME_VAR %in% names(d_mod)) {
  stop("Axis column missing: ", TIME_VAR, " — rerun 02_response_variables.R")
}
KNOTS <- make_knots(TIME_VAR)
AXIS_TAG <- Sys.getenv("TIME_AXIS", unset = "solar_double")
log_step("Data: n=%d, bears=%d | diel axis=%s (%s) | knots=c(%s)",
         nrow(d_mod), md$n_bears, TIME_VAR, AXIS_TAG, paste(KNOTS[[1]], collapse = ","))

# ---- SHARED FORMULA PARTS ----------------------------------------------------
re_term  <- "(1 | BearID_f/Year_f)"
doy_term <- "s(DOY, bs = 'tp', k = 8)"                       # BUG-2 fix (cc -> tp)
smooth   <- function(by = NULL) {
  if (is.null(by)) sprintf("s(%s, bs = 'cc', k = 12)", TIME_VAR)
  else sprintf("s(%s, by = %s, bs = 'cc', k = 12)", TIME_VAR, by)
}
# BUG-3 fix: the CONTINUOUS form of the tertile-cut variables is removed
covs_minus <- function(...) paste(setdiff(decorr_covs, c(...)), collapse = " + ")
base_covs  <- paste(decorr_covs, collapse = " + ")

# ---- HYPOTHESIS FORMULAS -----------------------------------------------------
rhs <- switch(HYP,

  # H0 — main effects (baseline)
  H0 = paste(smooth(), doy_term, base_covs, "Sex_f", "Age_sc", "Season_f", re_term,
             sep = " + "),

  # H1 — subsidy and activity budget: a Season x dump-distance interaction.
  # Same formula, run for BOTH responses; the finding is the CONTRAST of the
  # interaction across the two responses.
  H1 = paste(smooth(), doy_term,
             covs_minus("d2GarbageDump_km_sc"),
             "Sex_f", "Age_sc",
             "Season_f * d2GarbageDump_km_sc",
             re_term, sep = " + "),

  # H2 — diel shape ~ sex x season
  H2 = paste(smooth(), smooth("SexSeason_o"), doy_term,
             base_covs, "Sex_f", "Age_sc", "Season_f", "SexSeason_f",
             re_term, sep = " + "),

  # H3 — MAIN HYPOTHESIS: anthropogenic diel modulation, with an HOURLY thermal control
  H3 = paste(smooth(), smooth("DumpProx_o"), smooth("RoadProx_o"), smooth("TempTert_o"),
             "DumpProx_f", "RoadProx_f", "TempTert_f",
             doy_term,
             covs_minus("d2GarbageDump_km_sc", "d2Roads_km_sc", "temp_hourly_sc"),
             "Sex_f", "Age_sc", "Season_f",                   # BUG-4 fix
             re_term, sep = " + "),

  # H4 — thermal counter-hypothesis (now a REAL rival: hourly ERA5 temperature)
  H4 = paste(smooth(), smooth("TempTert_o"), "TempTert_f", doy_term,
             covs_minus("temp_hourly_sc"),
             "Sex_f", "Age_sc", "Season_f", re_term, sep = " + "),

  # H5 — individual plasticity (random slopes)
  # DECISION CRITERION: H5 vs H0 by LOO. A variance component whose CrI excludes
  # zero is NOT evidence (SDs are bounded >=0 and their CrIs almost never include 0).
  H5 = paste(smooth(), doy_term, base_covs, "Sex_f", "Age_sc", "Season_f",
             "(1 + d2GarbageDump_km_sc + d2Roads_km_sc | BearID_f)",
             "(1 | BearID_f:Year_f)", sep = " + ")
)

# Sex x season interaction factor for H2
if (HYP == "H2") {
  d_mod$SexSeason_f <- droplevels(interaction(d_mod$Sex_f, d_mod$Season_f, sep = "."))
  d_mod$SexSeason_o <- ordered(d_mod$SexSeason_f, levels = levels(d_mod$SexSeason_f))
  log_step("SexSeason: %s", paste(names(table(d_mod$SexSeason_f)),
           as.integer(table(d_mod$SexSeason_f)), sep = "=", collapse = " | "))
}

fml_str <- paste(RS$lhs, "~", rhs)
log_step("FORMULA: %s", fml_str)

# zi part only for zero_inflated_beta
bform <- if (RS$has_zi) {
  brms::bf(as.formula(fml_str), as.formula(sprintf("zi ~ %s", smooth())))
} else {
  brms::bf(as.formula(fml_str))
}

# ---- FIT ---------------------------------------------------------------------
# H3 is heaviest (4 smooths) -> higher adapt_delta/treedepth
heavy <- HYP %in% c("H3", "H5")
t0 <- Sys.time()
log_step("brms::brm start — %s/%s, n=%d, %s", HYP, RS$key, nrow(d_mod),
         ifelse(heavy, "heavy profile", "standard profile"))

fit <- brms::brm(
  formula = bform,
  data    = d_mod,
  family  = RS$family,
  prior   = RS$priors,
  knots   = KNOTS,                      # BUG-1 fix: cc boundary knots 0 and 24
  iter    = 3000, warmup = 1000,
  chains  = 4, cores = 4,
  seed    = 42, silent = 1, refresh = 250, init = 0,
  control = list(adapt_delta   = if (heavy) 0.98 else 0.95,
                 max_treedepth = if (heavy) 15 else 13),
  save_pars = brms::save_pars(all = TRUE)
)
elapsed <- difftime(Sys.time(), t0, units = "mins")
log_step("fit time: %.1f minutes", as.numeric(elapsed))

# ---- DIAGNOSTICS -------------------------------------------------------------
s <- summary(fit); fe <- s$fixed
cat("\n=== FIXED EFFECTS ===\n"); print(fe)
if (!is.null(s$splines)) { cat("\n=== SPLINES ===\n"); print(s$splines) }
if (!is.null(s$random))  { cat("\n=== RANDOM ===\n");  print(s$random) }

np       <- brms::nuts_params(fit)
ndiv     <- as.integer(sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE))
rhat_max <- max(fe[, "Rhat"], na.rm = TRUE)
ess_min  <- as.integer(min(fe[, "Bulk_ESS"], na.rm = TRUE))

# ---- SAVE --------------------------------------------------------------------
# The primary axis (solar_double) is written without a suffix; sensitivity axes
# go to a separate file so they do NOT overwrite the main models.
ax_sfx <- if (AXIS_TAG == "solar_double") "" else paste0("_", AXIS_TAG)
# DUMPSPEC suffix: the primary (raw) is unsuffixed so downstream scripts run
# unchanged; the sensitivity version (resid) goes to a SEPARATE file.
ds_sfx <- if (DUMPSPEC == "raw") "" else paste0("_dump", DUMPSPEC)
out_fp <- file.path(mod_path, sprintf("%s_%s%s%s.rds", HYP, RS$key, ax_sfx, ds_sfx))
save_rds_safe(list(fit = fit, hyp = HYP, response = RS$key,
                   family = RS$family_name, formula = fml_str,
                   time_var = TIME_VAR, time_axis = AXIS_TAG, knots = KNOTS,
                   dumpspec = DUMPSPEC,
                   decorr_covs = decorr_covs,
                   n = nrow(d_mod), n_bears = md$n_bears,
                   elapsed_min = as.numeric(elapsed),
                   rhat_max = rhat_max, ess_min = ess_min, divergent = ndiv,
                   data = d_mod),
              out_fp)

log_step("%s/%s: Rhat_max=%.4f | ESS_bulk_min=%d | divergent=%d",
         HYP, RS$key, rhat_max, ess_min, ndiv)
if (rhat_max > 1.01) log_step("!! WARNING: Rhat_max>1.01 — convergence suspect")
if (ndiv > 0)        log_step("!! WARNING: %d divergent transition(s)", ndiv)

# ---- FOCAL COEFFICIENTS ------------------------------------------------------
# NOTE: the variables are DISTANCES (d2GarbageDump = distance from the dump).
# "Proximity" is the OPPOSITE. We log this explicitly to avoid a sign misreading.
focal <- grep("d2GarbageDump|d2Roads|Season_f|DumpProx|RoadProx|TempTert",
              rownames(fe), value = TRUE)
if (length(focal) > 0) {
  cat("\n=== FOCAL COEFFICIENTS (NOTE: variables are DISTANCE, not proximity) ===\n")
  print(fe[focal, , drop = FALSE])
  cat("Reading rule: a POSITIVE d2X coefficient => response RISES with distance from X\n")
  cat("                                          => bears CLOSE to X show a LOWER response\n")
}

log_step("=== 05_fit_models.R (%s/%s) DONE ===", HYP, RS$key)
sink(type = "message"); sink()
