# ==============================================================================
# 31_axis_comparison.R  —  Diel-axis sensitivity: solar vs wall-clock
# ==============================================================================
#
# QUESTION (the original, still-valid question from today/31):
#   Is H3's evening peak a real biological moment, or a smear produced by averaging
#   the season-drifting sunset on the wall-clock axis?
#   At this latitude (~40N, per Section 2.1) sunset drifts ~3 HOURS across the season.
#   Wall-clock 18:00 is before sunset in June (broad daylight) and after it in
#   November (dark) — the wall-clock axis treats these as the SAME point.
#
# EXPECTATION: on the solar axis the peak should SHARPEN (smaller half-maximum
#   width), the amplitude should GROW and the dump Near-Far difference should
#   INCREASE. If not, wall-clock smearing matters less than assumed.
#
# PREREQUISITE (two fits):
#   TIME_AXIS=solar_double HYP=H3 RESP=intensity Rscript 05_fit_models.R  -> H3_intensity.rds
#   TIME_AXIS=clock        HYP=H3 RESP=intensity Rscript 05_fit_models.R  -> H3_intensity_clock.rds
#
# OUTPUTS:
#   tables/Table_V3_22_AxisComparison_<RESP>.csv
#     -> Table S46 (intensity) / S47 (duration)
#   figures/FigV3_11_AxisComparison_<RESP>.{png,pdf}
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))

RS <- get_response_spec()
init_log(sprintf("31_axis_comparison_%s", RS$key))
log_step("=== 31_axis_comparison | RESP=%s ===", RS$key)

N_DRAWS <- as.integer(Sys.getenv("N_DRAWS", "600"))

# ---- LOAD THE TWO MODELS -----------------------------------------------------
fps <- c(solar_double = file.path(mod_path, sprintf("H3_%s.rds", RS$key)),
         clock        = file.path(mod_path, sprintf("H3_%s_clock.rds", RS$key)))
miss <- fps[!file.exists(fps)]
if (length(miss) > 0) {
  stop("Missing fit(s):\n  ", paste(miss, collapse = "\n  "),
       "\n\nRun these first:\n",
       sprintf("  TIME_AXIS=solar_double HYP=H3 RESP=%s Rscript 05_fit_models.R\n", RS$key),
       sprintf("  TIME_AXIS=clock        HYP=H3 RESP=%s Rscript 05_fit_models.R\n", RS$key))
}

# ---- DIEL CURVE + METRICS PER AXIS -------------------------------------------
# G-computation: subsample real data rows, FORCE DumpProx, marginalise over rows.
# Does NOT create an "imaginary bear with median covariates".
curve_metrics <- function(fp, axis_name) {
  obj <- readRDS(fp)
  fit <- obj$fit; d <- obj$data; TV <- obj$time_var
  log_step("--- %s | %s | n=%d ---", axis_name, TV, nrow(d))

  set.seed(42)
  base <- d[sample.int(nrow(d), min(1200, nrow(d))), , drop = FALSE]
  hours <- seq(0, 24, length.out = 49)[-49]
  dt_h  <- 24 / length(hours)
  lvls  <- c("Near", "Mid", "Far")

  cube <- array(NA_real_, c(N_DRAWS, length(hours), length(lvls)),
                dimnames = list(NULL, NULL, lvls))
  for (hi in seq_along(hours)) {
    nd_l <- lapply(lvls, function(lv) {
      x <- base
      x[[TV]] <- hours[hi]
      x$DumpProx_f <- lv
      x$DumpProx_o <- factor(lv, levels = lvls, ordered = TRUE)
      x$.lvl <- lv
      x
    })
    nd <- do.call(rbind, nd_l)
    ep <- brms::posterior_epred(fit, newdata = nd, re_formula = NA, ndraws = N_DRAWS)
    if (!RS$has_zi) ep <- sweep(ep, 2, nd$N_readings, "/")   # to proportion
    for (li in seq_along(lvls)) {
      cube[, hi, li] <- rowMeans(ep[, nd$.lvl == lvls[li], drop = FALSE])
    }
  }

  out <- do.call(rbind, lapply(seq_along(lvls), function(li) {
    m     <- colMeans(cube[, , li])
    amp   <- max(m) - min(m)
    half  <- min(m) + 0.5 * amp
    data.frame(
      axis = axis_name, dump = lvls[li],
      peak_time      = round(hours[which.max(m)], 2),
      amplitude      = round(amp, 5),
      # half-maximum width: SMALL = SHARP peak. This is the main metric.
      halfmax_width_h = round(sum(m >= half) * dt_h, 2),
      auc            = round(sum(m) * dt_h, 5)
    )
  }))
  # max difference between the Near-Far curves = MAGNITUDE of the anthropogenic signal
  out$dump_maxdiff_near_far <- round(max(abs(colMeans(cube[, , "Near"]) -
                                             colMeans(cube[, , "Far"]))), 5)
  attr(out, "curves") <- data.frame(
    axis = axis_name,
    hour = rep(hours, times = 3),
    dump = rep(lvls, each = length(hours)),
    est  = c(colMeans(cube[, , 1]), colMeans(cube[, , 2]), colMeans(cube[, , 3]))
  )
  out
}

res <- list(); curves <- list()
for (ax in names(fps)) {
  m <- curve_metrics(fps[[ax]], ax)
  res[[ax]] <- m; curves[[ax]] <- attr(m, "curves")
  for (i in seq_len(nrow(m))) {
    log_step("  %s/%s: peak=%.2f, amplitude=%.5f, half-max width=%.2f h",
             ax, m$dump[i], m$peak_time[i], m$amplitude[i], m$halfmax_width_h[i])
  }
  log_step("  %s: dump Near-Far max diff = %.5f", ax, m$dump_maxdiff_near_far[1])
}
cmp <- do.call(rbind, res); rownames(cmp) <- NULL
cmp$response <- RS$key

# ---- VERDICT -----------------------------------------------------------------
w_solar <- mean(cmp$halfmax_width_h[cmp$axis == "solar_double"])
w_clock <- mean(cmp$halfmax_width_h[cmp$axis == "clock"])
d_solar <- cmp$dump_maxdiff_near_far[cmp$axis == "solar_double"][1]
d_clock <- cmp$dump_maxdiff_near_far[cmp$axis == "clock"][1]

verdict <- if (w_solar < w_clock * 0.9 && d_solar > d_clock * 1.1) {
  "SOLAR AXIS REQUIRED: wall-clock both smears the peak and weakens the dump signal"
} else if (w_solar < w_clock * 0.9) {
  "SOLAR AXIS sharpens the peak but the dump-signal magnitude is similar"
} else if (d_solar > d_clock * 1.1) {
  "SOLAR AXIS grows the dump signal but the peak sharpness is similar"
} else {
  "AXIS DIFFERENCE SMALL: wall-clock smearing is unimportant in this data (reportable null)"
}
log_step("Mean half-max width: solar=%.2f h, clock=%.2f h (%.1f%% diff)",
         w_solar, w_clock, 100 * (w_clock - w_solar) / w_clock)
log_step("Dump Near-Far max diff: solar=%.5f, clock=%.5f (%.1f%% diff)",
         d_solar, d_clock, 100 * (d_solar - d_clock) / d_clock)
log_step("VERDICT: %s", verdict)
cmp$verdict <- verdict

# -> Table S46 (intensity) / S47 (duration): solar-vs-clock axis comparison
write.csv(cmp, file.path(tbl_path, sprintf("Table_V3_22_AxisComparison_%s.csv", RS$key)),
          row.names = FALSE)
cat("\n=== AXIS COMPARISON ===\n"); print(cmp, row.names = FALSE)

# ---- FIGURE ------------------------------------------------------------------
cdf <- do.call(rbind, curves)
cdf$dump <- factor(cdf$dump, levels = c("Near","Mid","Far"))
cdf$axis <- factor(cdf$axis, levels = c("clock","solar_double"),
                   labels = c("Wall-clock hour (old)", "Double-anchored solar time (new)"))
ok <- c(Near = "#D55E00", Mid = "#009E73", Far = "#0072B2")

p <- ggplot2::ggplot(cdf, ggplot2::aes(hour, est, colour = dump)) +
  ggplot2::geom_line(linewidth = .9) +
  ggplot2::facet_wrap(~ axis, scales = "free_x") +
  ggplot2::scale_colour_manual(values = ok, name = "Dump distance") +
  ggplot2::labs(x = "Diel axis", y = RS$label,
                title = "Does the wall-clock axis blur the crepuscular peak?",
                subtitle = paste0("Sunset drifts ~3 h across Apr-Nov; the clock axis ",
                                  "averages over that drift | ", verdict)) +
  theme_bear()
save_fig(p, sprintf("FigV3_11_AxisComparison_%s", RS$key), w = 11, h = 4.5)

log_step("=== 31_axis_comparison DONE ===")
sink(type = "message"); sink()
