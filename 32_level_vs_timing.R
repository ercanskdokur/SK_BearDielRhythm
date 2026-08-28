# ==============================================================================
# 32_level_vs_timing.R  —  Is the dump effect LEVEL or TIMING?
# ==============================================================================
#
# QUESTION: does the "bears near the dump" finding come from a REDISTRIBUTION of
#   activity (timing / temporal niche) or simply from being LESS active (level /
#   energy saving)?
#
# METHOD — G-computation (marginal):
#   For each tertile, subsample real data rows, FORCE DumpProx to the target level,
#   take posterior_epred per hour, average over rows. The covariate distribution
#   reflects the real data (no median/reference "imaginary bear"). re_formula = NA
#   -> population-level curve.
#
# INFERENCE — per posterior draw:
#   AUC_g      = area under the curve       [LEVEL]
#   norm_g(t)  = curve / AUC                [SHAPE]
#   night_frac = night share of norm curve  [SHAPE summary]
#
# OUTPUTS:
#   tables/Table_V3_23_level_vs_timing_<RESP>.csv     -> Table S48 (intensity) / S49 (duration)
#   tables/Table_V3_24_normalized_curves_<RESP>.csv   -> Table S50 (intensity) / S51 (duration)
#   figures/FigV3_10_level_vs_timing_<RESP>.{png,pdf}
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))

RS <- get_response_spec()
init_log(sprintf("32_level_vs_timing_%s", RS$key))
log_step("=== 32_level_vs_timing | RESP=%s ===", RS$key)

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))

N_REP   <- as.integer(Sys.getenv("N_REP",   "1500"))  # G-comp row subsample
N_DRAWS <- as.integer(Sys.getenv("N_DRAWS", "400"))   # posterior draws

# ---- MODEL -------------------------------------------------------------------
obj <- load_model("H3")
fit <- obj$fit
d   <- obj$data
TIME_VAR <- obj$time_var
log_step("H3 | n=%d | diel axis=%s | family=%s", nrow(d), TIME_VAR, obj$family)
if (!"DumpProx_o" %in% names(d)) stop("DumpProx_o missing — this model is not in the H3 template")

# ---- Night hours: derive FROM THE DATA --------------------------------------
# No fixed assumption. t_period is now CORRECT (time zone moved to UTC in 01/02).
if ("t_period" %in% names(d)) {
  nt <- tapply(d$t_period == "nocturnal", round(d[[TIME_VAR]]), mean, na.rm = TRUE)
  night_hours <- as.numeric(names(nt))[which(nt > 0.5)]
  log_step("Night hours (on the %s axis, from t_period, >50%% nocturnal): %s",
           TIME_VAR, paste(night_hours, collapse = ","))
} else {
  night_hours <- c(21:23, 0:4)
  log_step("!! no t_period — default night hours: %s", paste(night_hours, collapse = ","))
}

# ---- G-computation -----------------------------------------------------------
set.seed(42)
base <- d[sample.int(nrow(d), min(N_REP, nrow(d))), , drop = FALSE]
hours <- seq(0, 24, length.out = 49)[-49]
dt_h  <- 24 / length(hours)
lvls  <- c("Near", "Mid", "Far")
log_step("G-computation: %d rows x %d hours x 3 tertiles, %d draws",
         nrow(base), length(hours), N_DRAWS)

cube <- array(NA_real_, c(N_DRAWS, length(hours), length(lvls)),
              dimnames = list(NULL, NULL, lvls))
t0 <- Sys.time()
for (hi in seq_along(hours)) {
  nd_l <- lapply(lvls, function(lv) {
    x <- base
    x[[TIME_VAR]] <- hours[hi]
    x$DumpProx_f  <- lv
    x$DumpProx_o  <- factor(lv, levels = lvls, ordered = TRUE)
    x$.lvl <- lv
    x
  })
  nd <- do.call(rbind, nd_l)
  ep <- brms::posterior_epred(fit, newdata = nd, re_formula = NA, ndraws = N_DRAWS)
  # beta-binomial epred returns COUNTS -> convert to a proportion so the two
  # responses are comparable
  if (!RS$has_zi) ep <- sweep(ep, 2, nd$N_readings, "/")
  for (li in seq_along(lvls)) {
    cube[, hi, li] <- rowMeans(ep[, nd$.lvl == lvls[li], drop = FALSE])
  }
  if (hi %% 8 == 0) log_step("  hour %d/%d (%.1f min)", hi, length(hours),
                             as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
log_step("posterior_epred total: %.1f min",
         as.numeric(difftime(Sys.time(), t0, units = "mins")))

# ---- LEVEL: AUC --------------------------------------------------------------
auc <- apply(cube, c(1, 3), sum) * dt_h          # draws x level

# ---- SHAPE: normalised curve + night share -----------------------------------
norm_cube <- cube
for (li in seq_along(lvls)) norm_cube[, , li] <- cube[, , li] / auc[, li]
night_idx  <- which(round(hours) %in% night_hours)
night_frac <- apply(norm_cube[, night_idx, , drop = FALSE], c(1, 3), sum) * dt_h

qsum <- function(x) c(mean = mean(x), q025 = unname(stats::quantile(x, .025)),
                      q50 = unname(stats::quantile(x, .5)),
                      q975 = unname(stats::quantile(x, .975)))
res <- do.call(rbind, lapply(seq_along(lvls), function(li) {
  data.frame(metric = c("AUC_total_activity", "night_fraction_of_curve"),
             dump = lvls[li],
             rbind(qsum(auc[, li]), qsum(night_frac[, li])), row.names = NULL)
}))

# ---- CONTRASTS: Near - Far ---------------------------------------------------
d_auc   <- auc[, "Near"] - auc[, "Far"]
d_night <- night_frac[, "Near"] - night_frac[, "Far"]
# whole shape difference: max absolute difference between the normalised curves
d_shape <- apply(abs(norm_cube[, , "Near"] - norm_cube[, , "Far"]), 1, max)
pd <- function(x) max(mean(x > 0), mean(x < 0))
excl0 <- function(x) !(stats::quantile(x, .025) < 0 & stats::quantile(x, .975) > 0)

contr <- data.frame(
  metric = c("AUC_total_activity", "night_fraction_of_curve"),
  dump = "Near - Far",
  rbind(qsum(d_auc), qsum(d_night)),
  pd = c(pd(d_auc), pd(d_night)),
  excludes_zero = c(excl0(d_auc), excl0(d_night)),
  row.names = NULL
)
auc_pct <- 100 * d_auc / auc[, "Far"]
log_step("AUC Near-Far: %.5f [%.5f, %.5f]  (%%%.1f [%%%.1f, %%%.1f]), pd=%.3f",
         mean(d_auc), stats::quantile(d_auc,.025), stats::quantile(d_auc,.975),
         mean(auc_pct), stats::quantile(auc_pct,.025), stats::quantile(auc_pct,.975),
         pd(d_auc))
log_step("night_frac Near-Far: %.5f [%.5f, %.5f], pd=%.3f",
         mean(d_night), stats::quantile(d_night,.025), stats::quantile(d_night,.975),
         pd(d_night))
log_step("Normalised curve max diff: %.5f [%.5f, %.5f]",
         mean(d_shape), stats::quantile(d_shape,.025), stats::quantile(d_shape,.975))

out <- rbind(res, contr[, names(res)])
out$pd <- c(rep(NA_real_, nrow(res)), contr$pd)
out$response <- RS$key; out$n_draws <- N_DRAWS; out$n_rep <- N_REP
out$time_axis <- TIME_VAR

# ---- VERDICT -----------------------------------------------------------------
lvl_sig   <- contr$excludes_zero[1]
shape_sig <- contr$excludes_zero[2]
verdict <- if (lvl_sig && !shape_sig) {
  "LEVEL: the dump changes total activity, not its distribution"
} else if (!lvl_sig && shape_sig) {
  "TIMING: the dump changes the WITHIN-DAY distribution of activity (temporal-niche thesis)"
} else if (lvl_sig && shape_sig) {
  "BOTH: total activity AND distribution change — report SEPARATELY"
} else {
  "BOTH UNCERTAIN: neither the level nor the shape contrast excludes 0"
}
log_step("VERDICT (%s): %s", RS$key, verdict)
log_step("NOTE: cross-check this verdict against the measurement-side")
log_step("     intensity-vs-duration contrast. They should say the same thing.")
out$verdict <- verdict

# -> Table S48 (intensity) / S49 (duration): level-vs-timing decomposition
write.csv(out, file.path(tbl_path, sprintf("Table_V3_23_level_vs_timing_%s.csv", RS$key)),
          row.names = FALSE)

curve_df <- do.call(rbind, lapply(seq_along(lvls), function(li) {
  raw <- cube[, , li]; nrm <- norm_cube[, , li]
  data.frame(dump = lvls[li], hour = hours,
             raw_mean = colMeans(raw),
             raw_lo = apply(raw, 2, stats::quantile, .025),
             raw_hi = apply(raw, 2, stats::quantile, .975),
             norm_mean = colMeans(nrm),
             norm_lo = apply(nrm, 2, stats::quantile, .025),
             norm_hi = apply(nrm, 2, stats::quantile, .975))
}))
# -> Table S50 (intensity) / S51 (duration): normalised diel curves
write.csv(curve_df,
          file.path(tbl_path, sprintf("Table_V3_24_normalized_curves_%s.csv", RS$key)),
          row.names = FALSE)

# ---- FIGURE: left = raw (level+shape), right = normalised (shape only) -------
ok_pal <- c(Near = "#D55E00", Mid = "#009E73", Far = "#0072B2")
curve_df$dump <- factor(curve_df$dump, levels = lvls)

p_raw <- ggplot2::ggplot(curve_df, ggplot2::aes(hour, raw_mean, colour = dump, fill = dump)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = raw_lo, ymax = raw_hi), alpha = .18, colour = NA) +
  ggplot2::geom_line(linewidth = .9) +
  ggplot2::scale_colour_manual(values = ok_pal, name = "Dump distance") +
  ggplot2::scale_fill_manual(values = ok_pal, name = "Dump distance") +
  ggplot2::labs(title = "Raw diel curves", subtitle = "Level + shape combined",
                x = "Solar time",
                y = if (RS$key == "duration") "Active duration" else "Movement intensity") +
  ggplot2::theme_minimal(base_size = 11)

p_norm <- ggplot2::ggplot(curve_df, ggplot2::aes(hour, norm_mean, colour = dump, fill = dump)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = norm_lo, ymax = norm_hi), alpha = .18, colour = NA) +
  ggplot2::geom_line(linewidth = .9) +
  ggplot2::scale_colour_manual(values = ok_pal, name = "Dump distance") +
  ggplot2::scale_fill_manual(values = ok_pal, name = "Dump distance") +
  ggplot2::labs(title = "Normalised diel curves",
                subtitle = "Shape only (each curve integrates to 1)",
                x = "Solar time", y = "Proportion of daily activity") +
  ggplot2::theme_minimal(base_size = 11)

p <- patchwork::wrap_plots(p_raw, p_norm, nrow = 1, guides = "collect") +
  patchwork::plot_annotation(
    title = sprintf("Is the dump effect level or timing? (%s)", RS$key)
  ) & ggplot2::theme(legend.position = "bottom")
save_fig(p, sprintf("FigV3_10_level_vs_timing_%s", RS$key), w = 11, h = 4.6)

log_step("=== 32_level_vs_timing DONE ===")
sink(type = "message"); sink()
