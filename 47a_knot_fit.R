# ==============================================================================
# 47a_knot_fit.R  —  Knot-sensitivity: SINGLE fit (SLURM array task)
# ==============================================================================
# 8 tasks = {H3,H4} x k{8,10,12,15}. Each task rebuilds the SAME fixed 18k
# bear-stratified subsample (seed=42) -> the fits share the same rows. Intensity
# (zi-beta).
#   sbatch --array=0-7 ...
# OUTPUT: models/knotsens_{HYP}_k{K}_intensity.rds  (fit + loo)  [feeds 47b]
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
suppressWarnings(suppressMessages({library(dplyr); library(loo)}))

GRID <- expand.grid(HYP = c("H3","H4"), K = c(8L,10L,12L,15L), stringsAsFactors = FALSE)
aid  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
stopifnot(aid >= 0, aid < nrow(GRID))
HYP <- GRID$HYP[aid + 1L]; K <- GRID$K[aid + 1L]

init_log(sprintf("47a_knot_fit_%s_k%d", HYP, K))
log_step("=== 47a_knot_fit | HYP=%s | k=%d | intensity ===", HYP, K)

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE)
options(brms.backend = "rstan")

RS <- get_response_spec("intensity"); stopifnot(RS$has_zi)
md <- readRDS(file.path(mod_path, "model_data_50k.rds"))
d_full <- md$data; decorr_covs <- md$decorr_covs
TIME_VAR <- get_time_var(); KNOTS <- make_knots(TIME_VAR)

set.seed(42L); TARGET_N <- 18000L
d18 <- d_full %>% dplyr::group_by(BearID_f) %>%
  dplyr::slice_sample(prop = TARGET_N / nrow(d_full)) %>% dplyr::ungroup()
log_step("Subsample: %s -> %s rows, %d bears",
         format(nrow(d_full), big.mark=","), format(nrow(d18), big.mark=","),
         dplyr::n_distinct(d18$BearID_f))

re_term  <- "(1 | BearID_f/Year_f)"
doy_term <- "s(DOY, bs = 'tp', k = 8)"
smooth_k <- function(K, by = NULL) {
  if (is.null(by)) sprintf("s(%s, bs = 'cc', k = %d)", TIME_VAR, K)
  else sprintf("s(%s, by = %s, bs = 'cc', k = %d)", TIME_VAR, by, K)
}
covs_minus <- function(...) paste(setdiff(decorr_covs, c(...)), collapse = " + ")
build_rhs <- function(HYP, K) {
  if (HYP == "H3")
    paste(smooth_k(K), smooth_k(K,"DumpProx_o"), smooth_k(K,"RoadProx_o"), smooth_k(K,"TempTert_o"),
          "DumpProx_f","RoadProx_f","TempTert_f", doy_term,
          covs_minus("d2GarbageDump_km_sc","d2Roads_km_sc","temp_hourly_sc"),
          "Sex_f","Age_sc","Season_f", re_term, sep = " + ")
  else
    paste(smooth_k(K), smooth_k(K,"TempTert_o"), "TempTert_f", doy_term,
          covs_minus("temp_hourly_sc"), "Sex_f","Age_sc","Season_f", re_term, sep = " + ")
}

fml <- paste(RS$lhs, "~", build_rhs(HYP, K))
bform <- brms::bf(as.formula(fml), as.formula(sprintf("zi ~ %s", smooth_k(K))))
log_step("FORMULA: %s", fml)

t0 <- Sys.time()
fit <- brms::brm(bform, data = d18, family = RS$family, prior = RS$priors, knots = KNOTS,
                 iter = 2000, warmup = 1000, chains = 4, cores = slurm_cpus,
                 seed = 42, silent = 1, refresh = 250, init = 0,
                 control = list(adapt_delta = 0.92, max_treedepth = 11),
                 save_pars = brms::save_pars(all = TRUE))
log_step("fit time: %.1f min", as.numeric(difftime(Sys.time(), t0, units="mins")))
fit <- brms::add_criterion(fit, "loo")
saveRDS(list(fit = fit, HYP = HYP, K = K), file.path(mod_path, sprintf("knotsens_%s_k%d_intensity.rds", HYP, K)))
log_step("=== 47a DONE: knotsens_%s_k%d_intensity.rds ===", HYP, K)
sink(type = "message"); sink()
