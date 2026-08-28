# ==============================================================================
# 42_tensor_H3.R  —  CONTINUOUS (varying-coefficient) version of H3  [B4 + B5]
# ==============================================================================
# Call:  RESP=intensity Rscript 42_tensor_H3.R   (or RESP=duration)
#
# WHY (audit B4/B5):
#   The tertile H3 (`s(time, by=DumpProx_o)` etc.) produced three problems:
#   (B4) the Mid class was not consistently INTERMEDIATE -> the tertile
#        discretisation misses the gradient (conflict windows, night-fraction,
#        diel curves all non-monotonic).
#   (B5) temperature entered only as a 3-level tertile, and the CONTINUOUS
#        temp_hourly_sc was dropped from H3 -> coarse thermal control + hour x
#        tertile collinearity.
#   FIX: no cut-points. The diel curve varies with the CONTINUOUS value of
#     dump/road/temperature (varying-coefficient cyclic smooth), and temperature
#     STAYS a CONTINUOUS main effect. This gives a genuine dose-response surface
#     and removes the three arbitrary cut points.
#
# FORMULA (the continuous counterpart of tertile H3):
#   y ~ s(time,cc,12)
#     + s(time, by=d2GarbageDump_km_sc, cc,12)   # diel shape ~ dump distance (continuous)
#     + s(time, by=d2Roads_km_sc, cc,12)         # diel shape ~ road distance (continuous)
#     + s(time, by=temp_hourly_sc, cc,12)        # diel shape ~ hourly temperature (continuous)
#     + s(DOY,tp,8)
#     + <decorr_covs: dump/road/temp_hourly STAY as main effects> + Sex+Age+Season
#     + (1|BearID/Year)
#   Note: in mgcv s(time,by=x) does NOT include x's MAIN effect; x is already linear
#   in decorr_covs -> B5's "keep the continuous temperature main effect" is satisfied.
#
# OUTPUTS: models/H3cont_<RESP>.rds
#          tables/Table_V3_25_TensorH3_smooths_<RESP>.csv  -> Table S52 (intensity) / S53 (duration)
#          tables/Table_V3_25_TensorH3_LOO_<RESP>.csv       (elpd: cont vs tertile H3; NOT a numbered SI table)
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
RS <- get_response_spec()
init_log(sprintf("42_tensor_H3_%s", RS$key))
log_step("=== 42_tensor_H3.R | RESP=%s (%s) ===", RS$key, RS$family_name)

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE); options(brms.backend = "rstan")

md <- readRDS(file.path(mod_path, "model_data_50k.rds"))
d  <- md$data
decorr_covs <- md$decorr_covs
TIME_VAR <- "SolarHour_dbl"; KNOTS <- make_knots(TIME_VAR)
log_step("Data n=%d, bears=%d | decorr_covs: %s", nrow(d), md$n_bears,
         paste(decorr_covs, collapse = ", "))

# ---- SUBSAMPLE ---------------------------------------------------------------
# At the full 50k this tensor varying-coefficient model is intractable in rstan
# (~140 s/iter, 4 continuous by-smooths + zi submodel). Only for this robustness
# analysis a bear-stratified ~18k subsample is taken; it is enough to estimate the
# dose-response surface. The primary H0-H5 models are fit at the FULL 50k and are
# not touched.
set.seed(42L)
TARGET_N <- 18000L
n_full <- nrow(d)
SUBSAMPLED <- FALSE
if (n_full > TARGET_N) {
  frac <- TARGET_N / n_full
  idx <- unlist(lapply(split(seq_len(n_full), d$BearID_f), function(ii) {
    k <- max(1L, round(length(ii) * frac))
    if (k >= length(ii)) ii else sample(ii, k)
  }), use.names = FALSE)
  d <- d[sort(idx), , drop = FALSE]
  d$BearID_f <- droplevels(d$BearID_f)
  if ("Year_f" %in% names(d)) d$Year_f <- droplevels(d$Year_f)
  SUBSAMPLED <- TRUE
  log_step("SUBSAMPLE: n %d -> %d (bear-stratified, frac=%.3f), bears=%d",
           n_full, nrow(d), frac, nlevels(d$BearID_f))
}

sm  <- function(by = NULL) {
  if (is.null(by)) sprintf("s(%s, bs='cc', k=12)", TIME_VAR)
  else sprintf("s(%s, by=%s, bs='cc', k=12)", TIME_VAR, by)
}
base_covs <- paste(decorr_covs, collapse = " + ")   # dump/road/temp_hourly main effects included

rhs <- paste(sm(),
             sm("d2GarbageDump_km_sc"), sm("d2Roads_km_sc"), sm("temp_hourly_sc"),
             "s(DOY, bs='tp', k=8)",
             base_covs, "Sex_f", "Age_sc", "Season_f",
             "(1 | BearID_f/Year_f)", sep = " + ")
fml_str <- paste(RS$lhs, "~", rhs)
log_step("FORMULA: %s", fml_str)

# zi submodel: a zi diel smooth (zi~s(time)) in the tensor added another cyclic
# smooth and slowed the model ~3x. The tensor's focus is how the MEAN diel curve
# varies with dump/road/temperature (B4/B5); zi is a nuisance here.
# -> zi INTERCEPT-only. The primary H0-H5 and tertile H3 keep zi~s(time).
bform <- if (RS$has_zi) { brms::bf(as.formula(fml_str), as.formula("zi ~ 1")) } else { brms::bf(as.formula(fml_str)) }

t0 <- Sys.time()
fit <- brms::brm(bform, data = d, family = RS$family, prior = RS$priors, knots = KNOTS,
                 iter = 2000, warmup = 1000, chains = 4, cores = 4, seed = 42,
                 silent = 1, refresh = 250, init = 0,
                 control = list(adapt_delta = 0.92, max_treedepth = 11),
                 save_pars = brms::save_pars(all = TRUE))
log_step("fit time: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

s <- summary(fit)
np <- brms::nuts_params(fit)
ndiv <- as.integer(sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE))
rhat_max <- max(s$fixed[, "Rhat"], na.rm = TRUE)
log_step("Rhat_max=%.4f | divergent=%d", rhat_max, ndiv)

# ---- smooth (varying-coefficient) summary ----------------------------------
if (!is.null(s$splines)) {
  sp <- as.data.frame(s$splines)
  sp$Term <- rownames(sp)
  # -> Table S52 (intensity) / S53 (duration): tensor-interaction smooth SDs
  utils::write.csv(sp[, c("Term", setdiff(names(sp), "Term"))],
                   file.path(tbl_path, sprintf("Table_V3_25_TensorH3_smooths_%s.csv", RS$key)),
                   row.names = FALSE)
  cat("\n=== SMOOTHS ===\n"); print(sp)
}

save_rds_safe(list(fit = fit, hyp = "H3cont", response = RS$key,
                   family = RS$family_name, formula = fml_str, time_var = TIME_VAR,
                   knots = KNOTS, decorr_covs = decorr_covs, n = nrow(d),
                   n_bears = md$n_bears, rhat_max = rhat_max, divergent = ndiv,
                   data = d),
              file.path(mod_path, sprintf("H3cont_%s.rds", RS$key)))

# ---- LOO: continuous H3 vs tertile H3 --------------------------------------
tryCatch({
  loo_c <- brms::loo(fit)
  h3t_fp <- file.path(mod_path, sprintf("H3_%s.rds", RS$key))
  if (!SUBSAMPLED && file.exists(h3t_fp)) {
    h3t <- readRDS(h3t_fp)$fit
    loo_t <- brms::loo(h3t)
    cmpl <- loo::loo_compare(list(H3cont = loo_c, H3tertile = loo_t))
    cmpdf <- as.data.frame(cmpl); cmpdf$model <- rownames(cmpdf)
    # LOO cross-check (NOT a numbered SI table)
    utils::write.csv(cmpdf, file.path(tbl_path,
                     sprintf("Table_V3_25_TensorH3_LOO_%s.csv", RS$key)), row.names = FALSE)
    cat("\n=== LOO cont vs tertile ===\n"); print(cmpl)
  } else if (SUBSAMPLED) {
    # elpd is not comparable with the full-50k tertile H3 (different n) -> SKIPPED.
    log_step("LOO vs tertile SKIPPED: tensor was subsampled (n=%d), tertile H3 is full 50k.", nrow(d))
  }
}, error = function(e) log_step("LOO ERR: %s", e$message))

log_step("=== 42_tensor_H3.R (%s) DONE ===", RS$key)
sink(type = "message"); sink()
