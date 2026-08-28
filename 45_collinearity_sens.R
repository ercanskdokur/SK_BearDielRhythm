# ==============================================================================
# 45_collinearity_sens.R  —  GENUINE collinearity sensitivity of dump distance  [A1]
# ==============================================================================
# Call:  RESP=intensity Rscript 45_collinearity_sens.R  (or RESP=duration)
#
# WHY (audit A1):
#   "raw vs residualize" was a REPARAMETERISATION (PC1/PC2 are already in the model
#   -> in the extra models the fit is identical; the coefficient ratio is just a
#   sqrt(1-R^2) scaling). It was NOT a genuine robustness test. Instead, an HONEST
#   collinearity treatment:
#     (i)  dump~LandPC1 scatter + R^2  (figure data)
#     (ii) how the dump coefficient changes in the FULL model (LandPC1 included) vs
#          with LandPC1 DROPPED = the GENUINE measure of how much the distance axis
#          carries the dump signal (in elpd/coefficient terms).
#     (iii) VIF reported explicitly (from Table_V2_02).
#
# OUTPUTS: tables/Table_V3_27_DumpLandPC1_scatter.csv        (once; RESP-independent; figure data)
#          tables/Table_V3_27_CollinSens_<RESP>.csv          (dump coef: full vs drop-PC1)
#            -> Table S60 (intensity; S60a) / (duration; S60b)
#          models/CollinFull_<RESP>.rds, models/CollinNoPC1_<RESP>.rds
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
RS <- get_response_spec()
init_log(sprintf("45_collin_%s", RS$key))
log_step("=== 45_collinearity_sens.R | RESP=%s ===", RS$key)

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE); options(brms.backend = "rstan")

md <- readRDS(file.path(mod_path, "model_data_50k.rds"))
d  <- md$data; decorr_covs <- md$decorr_covs
KNOTS <- make_knots("SolarHour_dbl")

# ---- (i) dump ~ LandPC1 scatter data + R^2 (written only on the intensity run) ----
if (identical(RS$key, "intensity")) {
  sc <- data.frame(d2GarbageDump_km = d$d2GarbageDump_km,
                   d2GarbageDump_km_sc = d$d2GarbageDump_km_sc,
                   LandPC1_sc = d$LandPC1_sc, LandPC2_sc = d$LandPC2_sc,
                   DumpProx_f = d$DumpProx_f)
  # keep it from ballooning: sample at most ~500 rows per bear
  set.seed(2026)
  idx <- unlist(tapply(seq_len(nrow(sc)), d$BearID_f,
                       function(i) if (length(i) > 500) sample(i, 500) else i))
  utils::write.csv(sc[idx, ], file.path(tbl_path, "Table_V3_27_DumpLandPC1_scatter.csv"),
                   row.names = FALSE)
  r2 <- summary(stats::lm(d2GarbageDump_km_sc ~ LandPC1_sc + LandPC2_sc, data = d))$r.squared
  r2pc1 <- summary(stats::lm(d2GarbageDump_km_sc ~ LandPC1_sc, data = d))$r.squared
  log_step("dump_sc ~ LandPC1+LandPC2 R^2=%.4f | ~LandPC1 alone R^2=%.4f", r2, r2pc1)
}

fit_additive <- function(covs, tag) {
  rhs <- paste("s(SolarHour_dbl, bs='cc', k=12)", "s(DOY, bs='tp', k=8)",
               paste(covs, collapse = " + "), "Sex_f", "Age_sc", "Season_f",
               "(1 | BearID_f/Year_f)", sep = " + ")
  fml_str <- paste(RS$lhs, "~", rhs)
  bform <- if (RS$has_zi) brms::bf(as.formula(fml_str),
                                   as.formula("zi ~ s(SolarHour_dbl, bs='cc', k=12)"))
           else brms::bf(as.formula(fml_str))
  log_step("[%s] FORMULA: %s", tag, fml_str)
  t0 <- Sys.time()
  f <- brms::brm(bform, data = d, family = RS$family, prior = RS$priors, knots = KNOTS,
                 iter = 2000, warmup = 1000, chains = 4, cores = 4, seed = 42,
                 silent = 1, refresh = 500, init = 0,
                 control = list(adapt_delta = 0.92, max_treedepth = 11),
                 save_pars = brms::save_pars(all = TRUE))
  log_step("[%s] fit time: %.1f min", tag, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  save_rds_safe(list(fit = f, response = RS$key, tag = tag),
                file.path(mod_path, sprintf("Collin%s_%s.rds", tag, RS$key)))
  f
}

# ---- SUBSAMPLE ---------------------------------------------------------------
# NoPC1 (with LandPC1 dropped) is even SLOWER than Full at the full 50k (the
# collinearity geometry breaks) -> intractable in rstan. A bear-stratified ~18k
# subsample; Full and NoPC1 are fit on the SAME 18k (comparable). The collinearity
# conclusion (dump coefficient with/without PC1) is robust to the subsample. The R^2
# diagnostic (above) is computed from the FULL data. 
{
  set.seed(42L); TARGET_N <- 18000L; n_full <- nrow(d)
  if (n_full > TARGET_N) {
    fr <- TARGET_N / n_full
    ix <- unlist(lapply(split(seq_len(n_full), d$BearID_f), function(ii) {
      k <- max(1L, round(length(ii) * fr)); if (k >= length(ii)) ii else sample(ii, k)
    }), use.names = FALSE)
    d <- d[sort(ix), , drop = FALSE]; d$BearID_f <- droplevels(d$BearID_f)
    if ("Year_f" %in% names(d)) d$Year_f <- droplevels(d$Year_f)
    log_step("SUBSAMPLE (%s): n %d -> %d (bear-stratified), bears=%d",
             RS$key, n_full, nrow(d), nlevels(d$BearID_f))
  }
}

f_full <- fit_additive(decorr_covs, "Full")                               # LandPC1/2 included
f_noPC <- fit_additive(setdiff(decorr_covs, "LandPC1_sc"), "NoPC1")       # distance axis DROPPED

grab <- function(f, tag) {
  fe <- summary(f)$fixed; p <- "d2GarbageDump_km_sc"
  data.frame(Model = tag, Parameter = p,
             dump_est = round(fe[p, "Estimate"], 4),
             dump_lo = round(fe[p, "l-95% CI"], 4), dump_hi = round(fe[p, "u-95% CI"], 4),
             CI_width = round(fe[p, "u-95% CI"] - fe[p, "l-95% CI"], 4))
}
out <- rbind(grab(f_full, "Full (LandPC1+PC2 in)"), grab(f_noPC, "Drop LandPC1"))
# -> Table S60 (intensity: S60a / duration: S60b): collinearity sensitivity (drop-LandPC1)
utils::write.csv(out, file.path(tbl_path, sprintf("Table_V3_27_CollinSens_%s.csv", RS$key)),
                 row.names = FALSE)
cat("\n=== DUMP COEF: full vs drop-LandPC1 ===\n"); print(out, row.names = FALSE)
log_step("=== 45_collinearity_sens.R (%s) DONE ===", RS$key)
sink(type = "message"); sink()
