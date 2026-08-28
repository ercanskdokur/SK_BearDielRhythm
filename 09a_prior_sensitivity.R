# ==============================================================================
# 09a_prior_sensitivity.R  —  Prior-sensitivity analysis (priorsense + multi-tier)
# ==============================================================================
# METHODS:
#   (1) priorsense::powerscale_sensitivity — analytic power-scaling
#   (2) Manual 3-tier: weak / moderate / tight priors fit + posterior comparison
#
# Runs on H0_hourly (zero_inflated_beta + decorrelated covariates). Because the
# decorrelated set is already well-behaved, sensitivity should typically be small.
#
# OUTPUTS:
#   Table_V2_14_PriorSenseAnalytical.csv   (analytic priorsense; NOT a numbered SI table)
#   Table_V2_15_PriorSensTier.csv          -> Table S57 (SI)
#   Table_V2_16_PriorSensConv.csv          -> Table S58 (SI)
#   FigV2_06_PriorSensitivity.{png,pdf}
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("09a_prior_sensitivity")
log_step("=== 09a_prior_sensitivity: START ===")

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE)
options(brms.backend = "rstan")
stan_cache <- file.path(tmp_dir, "stan_cache")
dir.create(stan_cache, recursive = TRUE, showWarnings = FALSE)

# ---- (1) priorsense power-scaling analytic test ----------------------------
ZIB <- "zero_inflated_beta"
h0 <- load_model("H0")
if (!identical(h0$family, ZIB)) {
  stop(sprintf("H0 family wrong: %s (expected zi_beta)", h0$family))
}
fit_default <- h0$fit
d_mod <- h0$data
SCALE_CONST <- h0$scale_const
log_step("H0 loaded: n=%d, family=%s, SCALE=%.0f",
         nrow(d_mod), h0$family, SCALE_CONST)
# The old `Activity_scaled` response no longer exists; H0 was fit with
# Act_intensity. That block produced a numeric(0) assignment -> REMOVED.
# priorsense operates on the fitted h0$fit; the manual tier uses Act_intensity.
log_step("priorsense::powerscale_sensitivity on H0 (Act_intensity)...")

# priorsense >=1.0.0: powerscale_sensitivity(fit)
ps_result <- tryCatch({
  priorsense::powerscale_sensitivity(fit_default)
}, error = function(e) {
  log_step("priorsense FAIL: %s", e$message)
  NULL
})

if (!is.null(ps_result)) {
  cat("\n=== priorsense::powerscale_sensitivity ===\n")
  tryCatch({
    print(ps_result)
    # Defensive: priorsense output structure is version-dependent; write whatever
    # it returns.
    ps_df <- tryCatch(as.data.frame(ps_result$sensitivity),
                       error = function(e) as.data.frame(ps_result))
    if (!"variable" %in% names(ps_df) && !is.null(rownames(ps_df))) {
      ps_df$variable <- rownames(ps_df); rownames(ps_df) <- NULL
    }
    # NOT a numbered SI table (analytic diagnostic)
    write.csv(ps_df,
              file.path(tbl_path, "Table_V2_14_PriorSenseAnalytical.csv"),
              row.names = FALSE)
    log_step("priorsense table written (cols: %s)",
             paste(names(ps_df), collapse = ","))
  }, error = function(e) {
    log_step("priorsense table FAIL: %s — continuing", e$message)
  })
} else {
  log_step("priorsense skipped — using the manual multi-tier sensitivity")
}

# ---- (2) MANUAL 3-TIER PRIOR FIT -------------------------------------------
log_step("Manual 3-tier prior fit...")
decorr_covs <- readRDS(file.path(mod_path, "model_data_50k.rds"))$decorr_covs
fixed_terms <- paste(c(decorr_covs, "Sex_f", "Age_sc"), collapse = " + ")
# Response Act_intensity, diel axis SolarHour_dbl, DOY tp (not cc — data filtered
# to Apr-Nov).
fml_str <- paste(
  "Act_intensity ~ s(SolarHour_dbl, bs='cc', k=12) + s(DOY, bs='tp', k=8) +",
  fixed_terms, "+ (1 | BearID_f/Year_f)"
)
bform <- brms::bf(as.formula(fml_str),
                  zi ~ s(SolarHour_dbl, bs = 'cc', k = 12))
KNOTS <- make_knots("SolarHour_dbl")

# 3 tiers — zi_beta priors (nuisance phi/zi_Intercept kept FIXED across tiers,
# only b/Intercept/sd width varied — a genuine "fixed effect" sensitivity)
priors_list <- list(
  "Tier1_weak" = c(
    brms::prior(normal(0, 5),         class = "b"),
    brms::prior(student_t(3, 0, 5),   class = "Intercept"),
    brms::prior(student_t(3, 0, 5),   class = "sd"),
    brms::prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "zi"),
    brms::prior(gamma(0.01, 0.01),    class = "phi")
  ),
  "Tier2_moderate" = c(   # same as H0
    brms::prior(normal(0, 1),         class = "b"),
    brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
    brms::prior(student_t(3, 0, 2.5), class = "sd"),
    brms::prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "zi"),
    brms::prior(gamma(0.01, 0.01),    class = "phi")
  ),
  "Tier3_tight" = c(
    brms::prior(normal(0, 0.3),       class = "b"),
    brms::prior(student_t(3, 0, 1),   class = "Intercept"),
    brms::prior(student_t(3, 0, 1),   class = "sd"),
    brms::prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "zi"),
    brms::prior(gamma(0.01, 0.01),    class = "phi")
  )
)

# H0 fit = Tier2 (already computed)
fits_tier <- list(Tier2_moderate = fit_default)

for (nm in c("Tier1_weak", "Tier3_tight")) {
  rds_fp <- file.path(mod_path, sprintf("fit_prior_%s.rds", nm))
  if (file.exists(rds_fp) && file.info(rds_fp)$size > 1e6) {
    f <- tryCatch(readRDS(rds_fp),
                  error = function(e) { log_step("  %s READ FAIL: %s", nm, e$message); NULL })
    # Family check: skip old hurdle_gamma caches
    fam_ok <- !is.null(f) &&
              tryCatch(identical(f$family$family, ZIB), error = function(e) FALSE)
    if (fam_ok) {
      log_step("CHECKPOINT (%s, family=%s): %s already present (%.1f MB), refit skipped.",
               ZIB, ZIB, nm, file.info(rds_fp)$size / 1024^2)
      fits_tier[[nm]] <- f; next
    }
    log_step("CHECKPOINT but family mismatch — %s will be REFIT", nm)
  }
  log_step("Fit %s...", nm)
  t0 <- Sys.time()
  f <- tryCatch(
    brms::brm(
      formula = bform,
      data    = d_mod,
      family  = brms::zero_inflated_beta(),
      prior   = priors_list[[nm]],
      knots   = KNOTS,
      iter    = 2000, warmup = 1000,
      chains  = 4, cores = 4,
      seed    = 42, silent = 1, refresh = 250,
      init    = 0,
      control = list(adapt_delta = 0.92, max_treedepth = 11),
      save_pars = brms::save_pars(all = TRUE)
    ),
    error = function(e) { log_step("  %s FAIL: %s", nm, e$message); NULL }
  )
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_step("  %s time: %.1f min", nm, dt)
  if (!is.null(f)) {
    save_rds_safe(f, rds_fp)
    fits_tier[[nm]] <- f
  }
}

# ---- POSTERIOR COEFFICIENT TABLE ACROSS TIERS ------------------------------
cmp_rows <- list()
for (nm in names(fits_tier)) {
  s <- summary(fits_tier[[nm]])$fixed
  cmp_rows[[nm]] <- data.frame(
    Tier      = nm,
    Term      = rownames(s),
    Estimate  = round(s[, "Estimate"], 4),
    SE        = round(s[, "Est.Error"], 4),
    Q2.5      = round(s[, "l-95% CI"], 4),
    Q97.5     = round(s[, "u-95% CI"], 4),
    Rhat      = round(s[, "Rhat"], 4),
    ESS_bulk  = round(s[, "Bulk_ESS"]),
    stringsAsFactors = FALSE
  )
}
cmp <- do.call(rbind, cmp_rows); rownames(cmp) <- NULL
# -> Table S57 (SI): prior-sensitivity, tier-wise coefficients
write.csv(cmp, file.path(tbl_path, "Table_V2_15_PriorSensTier.csv"), row.names = FALSE)

# ---- COMPARISON FIGURE -----------------------------------------------------
cmp_plot <- cmp %>% dplyr::filter(Term != "Intercept", !grepl("^s", Term))
p_sens <- ggplot2::ggplot(cmp_plot,
                          ggplot2::aes(x = Estimate, y = Term, color = Tier)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_pointrange(ggplot2::aes(xmin = Q2.5, xmax = Q97.5),
                           position = ggplot2::position_dodge(width = 0.6),
                           size = 0.4) +
  ggplot2::scale_color_manual(values = c(Tier1_weak = "#88CCEE",
                                         Tier2_moderate = "#332288",
                                         Tier3_tight = "#CC6677")) +
  ggplot2::labs(x = "Posterior mean (95% CrI)", y = NULL,
                title = "Prior sensitivity: 3-tier comparison (decorrelated covariates)",
                subtitle = sprintf("H0_hourly, zero_inflated_beta, n=%d", nrow(d_mod))) +
  theme_bear() +
  ggplot2::theme(legend.position = "top")
save_fig(p_sens, "FigV2_06_PriorSensitivity", w = 10, h = 7)

# Convergence summary
conv_rows <- list()
for (nm in names(fits_tier)) {
  s <- summary(fits_tier[[nm]])$fixed
  conv_rows[[nm]] <- data.frame(
    Tier         = nm,
    Rhat_max     = round(max(s[, "Rhat"], na.rm = TRUE), 4),
    ESS_bulk_min = round(min(s[, "Bulk_ESS"], na.rm = TRUE)),
    ESS_tail_min = round(min(s[, "Tail_ESS"], na.rm = TRUE)),
    converged    = all(s[, "Rhat"] < 1.05, na.rm = TRUE)
  )
}
conv_df <- do.call(rbind, conv_rows); rownames(conv_df) <- NULL
# -> Table S58 (SI): prior-sensitivity, tier-wise convergence
write.csv(conv_df,
          file.path(tbl_path, "Table_V2_16_PriorSensConv.csv"),
          row.names = FALSE)
cat("\n=== Convergence per tier ===\n"); print(conv_df, row.names = FALSE)

log_step("=== 09a_prior_sensitivity DONE ===")
sink(type = "message"); sink()
