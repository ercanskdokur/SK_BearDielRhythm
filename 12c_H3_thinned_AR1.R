# ==============================================================================
# 12c_H3_thinned_AR1.R  —  H3 + literal AR(1), on thinned data  [D2]
# ==============================================================================
# Call:  RESP=intensity Rscript 12c_H3_thinned_AR1.R  (or RESP=duration)
# brms ar() adds a latent AR term; because the full 50K is heavy, we thin to
# even hours within a day and add AR(1). If the fixed effects are consistent with
# the original H3, autocorrelation does not change the results. The AR time index
# is Hour_block (within-day ordering).
#
# OUTPUT: models/H3_thinAR1_<RESP>.rds  (feeds 48_ar_compare)
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
RS <- get_response_spec()
init_log(sprintf("12c_H3_thinned_AR1_%s", RS$key))
log_step("=== 12c_H3_thinned_AR1 | RESP=%s ===", RS$key)

options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
rstan::rstan_options(auto_write = TRUE); options(brms.backend = "rstan")

h3 <- load_model("H3"); d_mod <- h3$data
decorr_covs <- if (!is.null(h3$decorr_covs)) h3$decorr_covs else
               readRDS(file.path(mod_path, "model_data_50k.rds"))$decorr_covs
d_mod$BearDay_f <- factor(paste(d_mod$BearID_f, d_mod$Year_f, d_mod$DOY, sep = "_"))
d_thin <- d_mod[d_mod$Hour_block %% 2 == 0, ]
d_thin$BearDay_f <- droplevels(d_thin$BearDay_f)
log_step("Thinning: n %d -> %d, bear-day=%d", nrow(d_mod), nrow(d_thin),
         nlevels(d_thin$BearDay_f))
TIME_VAR <- "SolarHour_dbl"; KNOTS <- make_knots(TIME_VAR)

sm <- function(by = NULL) {
  if (is.null(by)) sprintf("s(%s, bs='cc', k=12)", TIME_VAR)
  else sprintf("s(%s, by=%s, bs='cc', k=12)", TIME_VAR, by)
}
covs_minus <- setdiff(decorr_covs, c("d2GarbageDump_km_sc","d2Roads_km_sc","temp_hourly_sc"))
h3_rhs <- paste(sm(), sm("DumpProx_o"), sm("RoadProx_o"), sm("TempTert_o"),
                "DumpProx_f", "RoadProx_f", "TempTert_f",
                "s(DOY, bs='tp', k=8)", paste(covs_minus, collapse = " + "),
                "Sex_f", "Age_sc", "Season_f", sep = " + ")
fml_str <- paste(RS$lhs, "~", h3_rhs,
                 "+ (1 | BearID_f/Year_f) + ar(time = Hour_block, gr = BearDay_f, p = 1)")
log_step("FORMULA: %s", fml_str)
bform <- if (RS$has_zi) { brms::bf(as.formula(fml_str), as.formula("zi ~ s(SolarHour_dbl, bs='cc', k=12)")) } else { brms::bf(as.formula(fml_str)) }

t0 <- Sys.time()
fit <- brms::brm(bform, data = d_thin, family = RS$family, prior = RS$priors, knots = KNOTS,
                 iter = 2000, warmup = 1000, chains = 4, cores = 4, seed = 42,
                 silent = 1, refresh = 250, init = 0,
                 control = list(adapt_delta = 0.97, max_treedepth = 12),
                 save_pars = brms::save_pars(all = TRUE))
log_step("12c fit time: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

cat("\n=== 12c FIXED EFFECTS ===\n"); print(round(brms::fixef(fit), 4))
ar_est <- tryCatch(brms::posterior_summary(fit, variable = "ar[1]"), error = function(e) NULL)
if (!is.null(ar_est)) { cat("\n=== AR(1) coefficient ===\n"); print(round(ar_est, 4)) }
save_rds_safe(list(fit = fit, response = RS$key, family = RS$family_name,
                   formula = fml_str, data = d_thin),
              file.path(mod_path, sprintf("H3_thinAR1_%s.rds", RS$key)))
log_step("=== 12c_H3_thinned_AR1 (%s) DONE ===", RS$key)
sink(type = "message"); sink()
