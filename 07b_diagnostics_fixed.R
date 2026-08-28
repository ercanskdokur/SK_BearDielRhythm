# ==============================================================================
# 07b_diagnostics_fixed.R  —  CORRECTED PPC + latent-scale VPC
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("07b_diagnostics_fixed")
log_step("=== 07b_diagnostics_fixed: START ===")

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))

# ---- Load models -----------------------------------------------------------
model_files <- c(
  stats::setNames(sprintf("%s_%s.rds", c("H0","H1","H2","H3","H4","H5"), RESP_KEY()),
                  c("H0","H1","H2","H3","H4","H5"))
)
fits <- list(); fams <- character()
for (nm in names(model_files)) {
  fp <- file.path(mod_path, model_files[nm])
  if (!file.exists(fp)) { log_step("MISSING (skip): %s", model_files[nm]); next }
  h <- readRDS(fp)
  fits[[nm]] <- h$fit; fams[nm] <- h$family
  log_step("Loaded %s (family=%s, n=%d)", nm, h$family, h$n)
}
if (length(fits) == 0) stop("No model loaded")

# ============================================================================
# 1) CORRECTED PPC STATISTICS
# ============================================================================
ppc_rows <- list()
for (nm in names(fits)) {
  fit <- fits[[nm]]
  log_step("PPC(fixed): %s", nm)
  y    <- as.numeric(brms::get_y(fit))
  yrep <- tryCatch(brms::posterior_predict(fit, ndraws = 500),
                   error = function(e) NULL)
  if (is.null(yrep)) next

  # beta_binomial -> convert to proportions (variable trials); others are [0,1].
  # NOTE: fams[nm] returns NAME-ATTRIBUTED, so identical(c(H0="beta_binomial"),
  #   "beta_binomial") is ALWAYS FALSE -> the beta_binomial branch never ran and
  #   duration PPC was computed on raw counts (0..12), with prop_at_upper actually
  #   being P(Active_n>=1). Use fams[[nm]] (drop the name).
  is_bb <- identical(fams[[nm]], "beta_binomial")
  if (is_bb) {
    tr <- fit$data$N_readings
    y_use    <- y / tr
    yrep_use <- sweep(yrep, 2, tr, "/")
    upper_obs <- mean(y == tr)
    upper_rep <- rowMeans(sweep(yrep, 2, tr, "==") == 1)  # eq trials
  } else {
    y_use <- y; yrep_use <- yrep
    upper_obs <- mean(y >= 0.999)
    upper_rep <- rowMeans(yrep >= 0.999)
  }

  stat_fns <- list(
    mean      = function(v) mean(v, na.rm = TRUE),
    sd        = function(v) stats::sd(v, na.rm = TRUE),
    prop_zero = function(v) mean(v == 0, na.rm = TRUE),   # <-- FIXED (true zero)
    q05       = function(v) unname(stats::quantile(v, 0.05, na.rm = TRUE)),
    q95       = function(v) unname(stats::quantile(v, 0.95, na.rm = TRUE))
  )

  for (st in names(stat_fns)) {
    obs <- stat_fns[[st]](y_use)
    rep <- apply(yrep_use, 1, stat_fns[[st]])
    pm  <- mean(rep, na.rm = TRUE); rsd <- stats::sd(rep, na.rm = TRUE)
    bayes_p <- mean(rep >= obs, na.rm = TRUE)
    rel_err <- if (abs(obs) > 1e-9) (obs - pm) / abs(obs) else (obs - pm)
    z_disc  <- if (rsd > 1e-12) (obs - pm) / rsd else NA_real_
    # MISFIT only if BOTH statistical (extreme p) AND material (>10% relative error)
    flag <- if ((bayes_p < 0.025 | bayes_p > 0.975) & abs(rel_err) > 0.10)
              "MISFIT" else if (bayes_p < 0.025 | bayes_p > 0.975)
              "small_n-detectable" else "OK"
    ppc_rows[[paste(nm, st)]] <- data.frame(
      Model = nm, Statistic = st,
      Observed = round(obs, 4), PostMean = round(pm, 4),
      Rel_err_pct = round(100 * rel_err, 1),
      Z_discrep = round(z_disc, 1),
      Bayes_p = round(bayes_p, 3), Flag = flag
    )
  }
  # upper-bound calibration (clipping/pile-up check)
  ppc_rows[[paste(nm, "prop_upper")]] <- data.frame(
    Model = nm, Statistic = "prop_at_upper",
    Observed = round(upper_obs, 4), PostMean = round(mean(upper_rep), 4),
    Rel_err_pct = round(100 * ((upper_obs - mean(upper_rep)) /
                                 max(abs(upper_obs), 1e-9)), 1),
    Z_discrep = NA_real_, Bayes_p = round(mean(upper_rep >= upper_obs), 3),
    Flag = if (abs(upper_obs - mean(upper_rep)) > 0.02) "check_clipping" else "OK"
  )
}
ppc_df <- do.call(rbind, ppc_rows); rownames(ppc_df) <- NULL
# -> Table S10 (intensity) / S11 (duration): posterior predictive checks
write.csv(ppc_df, file.path(tbl_path, paste0(tag_out("Table_V3_02_PPC_Stats"), ".csv")),
          row.names = FALSE)
cat("\n=== PPC (fixed) ===\n"); print(ppc_df, row.names = FALSE)

# ============================================================================
# 2) LATENT-SCALE VPC
# ============================================================================
# ---- NEW: family-specific + simulation-based CROSS-CHECK --------------------
obs_level_var <- function(fit, family_name) {
  dr <- posterior::as_draws_matrix(fit)
  mu <- mean(brms::posterior_epred(fit, re_formula = NA, ndraws = 200))
  if (family_name == "zero_inflated_beta") {
    # zi-Beta: observation-level logit-scale variance of the Beta component.
    # NOTE: mu here is the zi-mixed epred; rescale by the zi probability to get the
    # Beta component's mu.
    zi <- if ("zi" %in% colnames(dr)) mean(dr[, "zi"]) else 0
    mu_beta <- min(max(mu / max(1 - zi, 1e-6), 1e-4), 1 - 1e-4)
    phi <- mean(dr[, "phi"])
    list(var = trigamma(mu_beta * phi) + trigamma((1 - mu_beta) * phi),
         label = sprintf("Beta-logit trigamma (mu=%.3f, phi=%.1f)", mu_beta, phi))
  } else if (family_name == "beta_binomial") {
    tr_mean <- mean(fit$data$N_readings, na.rm = TRUE)
    mu_p <- min(max(mu / tr_mean, 1e-4), 1 - 1e-4)
    phi <- mean(dr[, "phi"])
    list(var = pi^2 / 3 + trigamma(mu_p * phi) + trigamma((1 - mu_p) * phi),
         label = sprintf("BetaBinom-logit (p=%.3f, phi=%.2f)", mu_p, phi))
  } else {
    list(var = NA_real_, label = paste("UNKNOWN family:", family_name))
  }
}

vpc_rows <- list()
for (nm in names(fits)) {
  fit <- fits[[nm]]
  fam <- tryCatch(fit$family$family, error = function(e) NA_character_)
  vc  <- tryCatch(brms::VarCorr(fit), error = function(e) NULL)
  if (is.null(vc)) next
  re_var <- numeric()
  for (g in names(vc)) {
    if (!is.null(vc[[g]]$sd) && nrow(vc[[g]]$sd) > 0) {
      sd_int <- vc[[g]]$sd[grep("Intercept", rownames(vc[[g]]$sd))[1], "Estimate"]
      if (!is.na(sd_int)) re_var[g] <- sd_int^2
    }
  }
  if (length(re_var) == 0) next

  olv <- obs_level_var(fit, fam)
  log_step("%s (%s): observation-level variance = %.4f  [%s]", nm, fam, olv$var, olv$label)
  log_step("   (old code used pi^2/3 = %.4f here — WRONG)", pi^2/3)

  re_total     <- sum(re_var)
  denom_latent <- re_total + olv$var
  for (g in names(re_var)) {
    vpc_rows[[paste(nm, g)]] <- data.frame(
      Model = nm, Family = fam, Group = g,
      Var = round(re_var[g], 4),
      Var_dist = round(olv$var, 4),
      VPC_latent_pct = round(100 * re_var[g] / denom_latent, 2),
      VPC_pi2_3_OLD = round(100 * re_var[g] / (re_total + pi^2/3), 2)  # for transparency
    )
  }
  vpc_rows[[paste(nm, "Distribution")]] <- data.frame(
    Model = nm, Family = fam, Group = olv$label,
    Var = round(olv$var, 4), Var_dist = round(olv$var, 4),
    VPC_latent_pct = round(100 * olv$var / denom_latent, 2),
    VPC_pi2_3_OLD = NA_real_
  )
}
vpc_df <- do.call(rbind, vpc_rows); rownames(vpc_df) <- NULL
# -> Table S12 (intensity) / S13 (duration): latent-scale variance partition (VPC)
write.csv(vpc_df, file.path(tbl_path, paste0(tag_out("Table_V3_03_VPC_latent"), ".csv")), row.names = FALSE)
cat("\n=== VPC latent-scale (family-specific) ===\n")
print(vpc_df, row.names = FALSE)
cat("\nNOTE: the 'VPC_pi2_3_OLD' column is for TRANSPARENCY only — it shows what the\n")
cat("     old (wrong) constant produced. It does NOT go into the report.\n")

# ---- SIMULATION-BASED CROSS-CHECK ------------------------------------------
# performance::variance_decomposition() works via posterior predictive simulation
# and does NOT depend on the family formula -> the best way to verify the analytic
# result. If performance is unavailable in the container this block is silently
# skipped (the analytic result is still written).
if (!requireNamespace("performance", quietly = TRUE)) {
  log_step("!! package 'performance' unavailable — simulation cross-check SKIPPED.")
  log_step("   To install (user-lib, outside the container): install.packages('performance')")
  log_step("   The analytic VPC (Table_V3_03) was still written.")
} else {
  vd_rows <- list()
  for (nm in names(fits)) {
    vd <- tryCatch(performance::variance_decomposition(fits[[nm]], ci = 0.95),
                   error = function(e) {
                     log_step("  var_decomp %s error: %s", nm, conditionMessage(e)); NULL })
    if (is.null(vd)) next
    vd_rows[[nm]] <- data.frame(Model = nm,
                                ICC_sim = round(vd$ICC_decomposed, 4),
                                CI_low  = round(vd$ICC_CI[1], 4),
                                CI_high = round(vd$ICC_CI[2], 4))
    log_step("%s: simulation ICC = %.4f [%.4f, %.4f]", nm,
             vd$ICC_decomposed, vd$ICC_CI[1], vd$ICC_CI[2])
  }
  if (length(vd_rows) > 0) {
    vd_df <- do.call(rbind, vd_rows); rownames(vd_df) <- NULL
    # cross-check output; NOT a numbered SI table
    write.csv(vd_df, file.path(tbl_path, paste0(tag_out("Table_V3_03b_VPC_simulated"), ".csv")), row.names = FALSE)
    cat("\n=== VPC — simulation-based cross-check ===\n"); print(vd_df, row.names = FALSE)
    cat("\nIf analytic (V3_03) and simulation (V3_03b) DISAGREE, trust the simulation.\n")
  }
}

# ============================================================================
# 3) ZERO-CALIBRATION FIGURE (observed vs model zero fraction)
# ============================================================================
zr <- ppc_df[ppc_df$Statistic == "prop_zero", ]
if (nrow(zr) > 0) {
  zr$Model <- factor(zr$Model, levels = c("H0","H1","H2","H3","H4","H5"))
  ok6 <- c(H0="#999999", H1="#E69F00", H2="#56B4E9",
           H3="#009E73", H4="#D55E00", H5="#CC79A7")
  p <- ggplot2::ggplot(zr, ggplot2::aes(x = Observed, y = PostMean,
                                        colour = Model, label = Model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey60") +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_colour_manual(values = ok6, name = "Hypothesis") +
    ggplot2::labs(x = "Observed zero proportion", y = "Model-predicted zero proportion",
                  title = "Zero-inflation calibration (correct prop_zero)",
                  subtitle = "Near line = good fit; far = zi component under/over-estimated") +
    theme_bear()
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(size = 3.4, show.legend = FALSE,
              max.overlaps = Inf, box.padding = 0.7, point.padding = 0.5,
              min.segment.length = 0, seed = 42, segment.colour = "grey55")
  } else {
    p <- p + ggplot2::geom_text(vjust = -0.9, hjust = -0.2, size = 3.4,
                                show.legend = FALSE)
  }
  save_fig(p, tag_out("FigV3_01_PPC_zero_calibration"), w = 7.5, h = 6)
}

log_step("=== 07b_diagnostics_fixed DONE ===")
sink(type = "message"); sink()
