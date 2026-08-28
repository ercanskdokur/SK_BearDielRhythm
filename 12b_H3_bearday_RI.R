# ==============================================================================
# 12b_H3_bearday_RI.R  —  H3 + bear-day random intercept (autocorrelation robust) [D2]
# ==============================================================================
# Call:  RESP=intensity Rscript 12b_H3_bearday_RI.R  (or RESP=duration)
# Within-day hours are not independent. Adding (1|BearDay) to H3 absorbs
# day-level correlation; if the fixed effects/diel smooths stay the same as the
# original H3, the results are robust.
#
# OUTPUT: models/H3_beardayRI_<RESP>.rds  (feeds 48_ar_compare)
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
RS <- get_response_spec()
init_log(sprintf("12b_H3_bearday_RI_%s", RS$key))
log_step("=== 12b_H3_bearday_RI | RESP=%s ===", RS$key)

options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
rstan::rstan_options(auto_write = TRUE); options(brms.backend = "rstan")

h3 <- load_model("H3"); d_mod <- h3$data
decorr_covs <- if (!is.null(h3$decorr_covs)) h3$decorr_covs else
               readRDS(file.path(mod_path, "model_data_50k.rds"))$decorr_covs
d_mod$BearDay_f <- factor(paste(d_mod$BearID_f, d_mod$Year_f, d_mod$DOY, sep = "_"))
TIME_VAR <- "SolarHour_dbl"; KNOTS <- make_knots(TIME_VAR)
log_step("n=%d, bear-day groups=%d", nrow(d_mod), nlevels(d_mod$BearDay_f))

sm <- function(by = NULL) {
  if (is.null(by)) sprintf("s(%s, bs='cc', k=12)", TIME_VAR)
  else sprintf("s(%s, by=%s, bs='cc', k=12)", TIME_VAR, by)
}
covs_minus <- setdiff(decorr_covs, c("d2GarbageDump_km_sc","d2Roads_km_sc","temp_hourly_sc"))
h3_rhs <- paste(sm(), sm("DumpProx_o"), sm("RoadProx_o"), sm("TempTert_o"),
                "DumpProx_f", "RoadProx_f", "TempTert_f",
                "s(DOY, bs='tp', k=8)", paste(covs_minus, collapse = " + "),
                "Sex_f", "Age_sc", "Season_f", sep = " + ")
fml_str <- paste(RS$lhs, "~", h3_rhs, "+ (1 | BearID_f/Year_f) + (1 | BearDay_f)")
log_step("FORMULA: %s", fml_str)
bform <- if (RS$has_zi) { brms::bf(as.formula(fml_str), as.formula(sm2 <- "zi ~ s(SolarHour_dbl, bs='cc', k=12)")) } else { brms::bf(as.formula(fml_str)) }

t0 <- Sys.time()
fit <- brms::brm(bform, data = d_mod, family = RS$family, prior = RS$priors, knots = KNOTS,
                 iter = 2000, warmup = 1000, chains = 4, cores = 4, seed = 42,
                 silent = 1, refresh = 250, init = 0,
                 control = list(adapt_delta = 0.92, max_treedepth = 11),
                 save_pars = brms::save_pars(all = TRUE))
log_step("12b fit time: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

cat("\n=== 12b FIXED EFFECTS ===\n"); print(round(brms::fixef(fit), 4))
save_rds_safe(list(fit = fit, response = RS$key, family = RS$family_name,
                   formula = fml_str, data = d_mod),
              file.path(mod_path, sprintf("H3_beardayRI_%s.rds", RS$key)))
log_step("=== 12b_H3_bearday_RI (%s) DONE ===", RS$key)
sink(type = "message"); sink()
