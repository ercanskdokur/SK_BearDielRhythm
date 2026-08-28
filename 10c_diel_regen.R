# ==============================================================================
# 10c_diel_regen.R  —  Model-based diel curves: SOLAR-AXIS + MARGINAL regeneration
# ==============================================================================
# FIX — the G-computation template from script 32 (level_vs_timing):
#   - diel axis = obj$time_var (= SolarHour_dbl); only the TARGET grouping factor
#     and the hour are forced, all OTHER covariates are kept at their REAL observed
#     values (no median "imaginary bear") -> a true population marginal.
#   - re_formula = NA -> population level.
#   - beta_binomial epred returns COUNTS -> divide by N_readings to get a
#     PROPORTION (RS$has_zi==FALSE).
#   - y = true unit (RS$label); x = solar time (SR=6, SS=18), no clock hours.
#
# OUTPUTS (one per RESP, tag_out appends RESP):
#   FigV3_03_H3_DielByDumpProx_<RESP>.{png,pdf}
#   FigV3_04_H3_DielByRoadProx_<RESP>.{png,pdf}
#   FigV3_05_H4_DielByTemp_<RESP>.{png,pdf}
#   FigV2_M2_H2_DielBySexSeason_<RESP>.{png,pdf}
#   tables/Table_V3_30_DielCurves_<fig>_<RESP>.csv   (numeric curves; NOT a numbered SI table)
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))

RS <- get_response_spec()
init_log(sprintf("10c_diel_regen_%s", RS$key))
log_step("=== 10c_diel_regen | RESP=%s | has_zi=%s ===", RS$key, RS$has_zi)

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))

N_REP   <- as.integer(Sys.getenv("N_REP",   "1500"))   # G-comp row subsample
N_DRAWS <- as.integer(Sys.getenv("N_DRAWS", "400"))    # posterior draws
OK <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")

# ---- Derive night hours FROM THE DATA (same as script 32) --------------------
derive_night <- function(d, tv) {
  if ("t_period" %in% names(d)) {
    nt <- tapply(d$t_period == "nocturnal", round(d[[tv]]), mean, na.rm = TRUE)
    as.numeric(names(nt))[which(nt > 0.5)]
  } else c(21:23, 0:4)
}

# ---- Marginal diel curve (G-comp) for one model + one grouping factor --------
# group_f: parametric/label factor name (e.g. DumpProx_f); group_o: by-smooth ordered
# name (e.g. DumpProx_o). group_levels: levels to plot. Only these + the hour are forced.
gcomp_diel <- function(hyp, group_f, group_o, group_levels) {
  obj <- load_model(hyp); fit <- obj$fit; d <- obj$data; tv <- obj$time_var
  log_step("[%s] n=%d axis=%s family=%s group=%s (%s)", hyp, nrow(d), tv, obj$family,
           group_f, paste(group_levels, collapse="/"))
  have_f <- group_f %in% names(d)
  have_o <- group_o %in% names(d)
  if (!have_f && !have_o) stop(sprintf("neither %s nor %s is in the %s data", group_f, group_o, hyp))

  set.seed(42L)
  base  <- d[sample.int(nrow(d), min(N_REP, nrow(d))), , drop = FALSE]
  hours <- seq(0, 24, length.out = 49)[-49]
  cube  <- array(NA_real_, c(N_DRAWS, length(hours), length(group_levels)),
                 dimnames = list(NULL, NULL, group_levels))
  t0 <- Sys.time()
  for (hi in seq_along(hours)) {
    nd_l <- lapply(group_levels, function(lv) {
      x <- base
      x[[tv]] <- hours[hi]
      if (have_f) x[[group_f]] <- factor(lv,  levels = group_levels)
      if (have_o) x[[group_o]] <- ordered(lv, levels = group_levels)
      x$.lvl <- lv
      x
    })
    nd <- do.call(rbind, nd_l)
    ep <- brms::posterior_epred(fit, newdata = nd, re_formula = NA, ndraws = N_DRAWS)
    if (!RS$has_zi) ep <- sweep(ep, 2, nd$N_readings, "/")   # beta_binom: count -> proportion
    for (li in seq_along(group_levels))
      cube[, hi, li] <- rowMeans(ep[, nd$.lvl == group_levels[li], drop = FALSE])
    if (hi %% 12 == 0)
      log_step("  hour %d/%d (%.1f min)", hi, length(hours),
               as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  log_step("  epred total: %.1f min", as.numeric(difftime(Sys.time(), t0, units="mins")))

  curve <- do.call(rbind, lapply(seq_along(group_levels), function(li) {
    m <- cube[, , li]
    data.frame(group = group_levels[li], hour = hours,
               est = colMeans(m),
               lo  = apply(m, 2, stats::quantile, .025),
               hi  = apply(m, 2, stats::quantile, .975))
  }))
  curve$group <- factor(curve$group, levels = group_levels)
  list(curve = curve, night = derive_night(d, tv), tv = tv)
}

# ---- Plot: solar axis + night shading + real-unit y --------------------------
plot_diel <- function(gc, pal, fig_name, subtitle, legend_lab) {
  cur <- gc$curve; nh <- gc$night
  # split the night hours into contiguous blocks (typical on the solar axis: [0,6) and (18,24])
  night_rects <- data.frame(xmin = numeric(0), xmax = numeric(0))
  if (length(nh)) {
    sh <- sort(nh); brk <- c(0, which(diff(sh) > 1), length(sh))
    for (b in seq_len(length(brk) - 1)) {
      seg <- sh[(brk[b] + 1):brk[b + 1]]
      night_rects <- rbind(night_rects,
                           data.frame(xmin = min(seg) - .5, xmax = max(seg) + .5))
    }
  }
  p <- ggplot2::ggplot(cur, ggplot2::aes(hour, est, colour = group, fill = group))
  if (nrow(night_rects))
    p <- p + ggplot2::geom_rect(data = night_rects, inherit.aes = FALSE,
             ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
             alpha = 0.08, fill = "#2C3E50")
  p <- p +
    ggplot2::geom_vline(xintercept = c(6, 18), linetype = "dashed",
                        colour = "grey55", linewidth = .35) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = .16, colour = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_manual(values = pal, name = legend_lab) +
    ggplot2::scale_fill_manual(values = pal, name = legend_lab) +
    ggplot2::scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
      labels = c("0", "sunrise (6)", "12", "sunset (18)", "24"),
      limits = c(0, 24), expand = c(0, 0)) +
    ggplot2::labs(x = "Solar time  (double-anchored: sunrise=6, sunset=18)",
                  y = RS$label, title = "Marginal diel activity curve",
                  subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
  save_fig(p, fig_name, w = 8, h = 5)
  write.csv(cur, file.path(tbl_path, sprintf("Table_V3_30_DielCurves_%s.csv",
            sub("^Fig[^_]*_[0-9A-Za-z]+_", "", fig_name))), row.names = FALSE)
  log_step("Drawn + table: %s", fig_name)
}

tri  <- c(Near = "#D55E00", Mid = "#009E73", Far = "#0072B2")
ther <- c(Cool = "#0072B2", Mild = "#009E73", Warm = "#D55E00")

# ---- FigV3_03: H3 x dump proximity -------------------------------------------
g03 <- gcomp_diel("H3", "DumpProx_f", "DumpProx_o", c("Near","Mid","Far"))
plot_diel(g03, tri, tag_out("FigV3_03_H3_DielByDumpProx"),
          "H3, marginal over all other covariates — by garbage-dump proximity",
          "Dump proximity")

# ---- FigV3_04: H3 x road proximity -------------------------------------------
g04 <- gcomp_diel("H3", "RoadProx_f", "RoadProx_o", c("Near","Mid","Far"))
plot_diel(g04, tri, tag_out("FigV3_04_H3_DielByRoadProx"),
          "H3, marginal over all other covariates — by road proximity",
          "Road proximity")

# ---- FigV3_05: H4 x thermal regime -------------------------------------------
g05 <- gcomp_diel("H4", "TempTert_f", "TempTert_o", c("Cool","Mild","Warm"))
plot_diel(g05, ther, tag_out("FigV3_05_H4_DielByTemp"),
          "H4, marginal over all other covariates — by hourly thermal regime",
          "Thermal regime")

# ---- FigV2_M2: H2 x sex x season ----------------------------------------
obj_h2 <- load_model("H2")
ss_var_f <- if ("SexSeason_f" %in% names(obj_h2$data)) "SexSeason_f" else NA
ss_var_o <- if ("SexSeason_o" %in% names(obj_h2$data)) "SexSeason_o" else NA
if (!is.na(ss_var_o) || !is.na(ss_var_f)) {
  ss_levels <- levels(obj_h2$data[[ if (!is.na(ss_var_o)) ss_var_o else ss_var_f ]])
  ss_pal <- stats::setNames(OK[seq_along(ss_levels)], ss_levels)
  gM2 <- gcomp_diel("H2",
                    if (!is.na(ss_var_f)) ss_var_f else "SexSeason_f",
                    if (!is.na(ss_var_o)) ss_var_o else "SexSeason_o",
                    ss_levels)
  plot_diel(gM2, ss_pal, tag_out("FigV2_M2_H2_DielBySexSeason"),
            "H2, marginal over all other covariates — by sex x season",
            "Sex . Season")
} else {
  log_step("!! No SexSeason factor in the H2 data — FigV2_M2 SKIPPED")
}

log_step("=== 10c_diel_regen DONE ===")
sink(type = "message"); sink()
