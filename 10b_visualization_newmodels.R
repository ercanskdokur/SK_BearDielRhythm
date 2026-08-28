# ==============================================================================
# 10b_visualization_newmodels.R  —  H3/H4 diel curves (human pressure & thermal)
# ==============================================================================
# Plots the diel activity curve per tertile, on a common axis, from H3 (hour x
# human pressure) and H4 (hour x temperature). GOAL: show the nocturnality shift
# VISUALLY — "do individuals near the dump push peak activity into the
# night/twilight hours?"
#
# OUTPUTS:
#   FigV3_03_H3_DielByDumpProx.{png,pdf}
#   FigV3_04_H3_DielByRoadProx.{png,pdf}
#   FigV3_05_H4_DielByTemp.{png,pdf}
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("10b_visualization_newmodels")
log_step("=== 10b_visualization_newmodels: START ===")

# Helper: build a diel-curve grid for a model + by-group factor, then plot it
plot_diel_by <- function(model_fp, group_var, group_levels, group_pal,
                         fig_name, subtitle, legend_lab = group_var) {
  if (!file.exists(model_fp)) { log_step("MISSING: %s", model_fp); return(invisible()) }
  h <- readRDS(model_fp); fit <- h$fit; d <- h$data
  tv <- if (!is.null(h$time_var)) h$time_var else "SolarHour_dbl"  # model time variable (cyclic spline; the grid must carry it)
  is_dur <- "N_readings" %in% names(d)  # beta-binom (duration) -> trials(N_readings) needed

  # Fix other _sc covariates at their median; factors at their reference level
  sc_vars <- grep("_sc$", names(d), value = TRUE)
  base_row <- as.data.frame(lapply(d[, sc_vars, drop = FALSE],
                                   function(x) stats::median(x, na.rm = TRUE)))
  base_row$DOY <- 200L
  base_row$Sex_f <- factor("F", levels = levels(d$Sex_f))
  if ("Season_f" %in% names(d))
    base_row$Season_f <- factor("Hyperphagia", levels = levels(d$Season_f))

  hours <- 0:23
  grid <- do.call(rbind, lapply(group_levels, function(gl) {
    g <- base_row[rep(1, length(hours)), , drop = FALSE]
    g$Hour_block <- hours
    if (!identical(tv, "Hour_block")) g[[tv]] <- hours   # model expects SolarHour_dbl (BUG FIX: was missing -> epred failed)
    if (is_dur) g$N_readings <- stats::median(d$N_readings, na.rm = TRUE)
    g[[group_var]] <- factor(gl, levels = group_levels)
    # ordered counterpart (for the by-smooth)
    og <- paste0(sub("_f$", "", group_var), "_o")
    if (og %in% names(d)) g[[og]] <- ordered(gl, levels = group_levels)
    # fix the other by-factors at the reference (H3 has dump/road/temp together)
    for (other in c("DumpProx_f","RoadProx_f","TempTert_f")) {
      if (other %in% names(d) && other != group_var) {
        lv <- levels(d[[other]])
        g[[other]] <- factor(lv[1], levels = lv)
        oo <- paste0(sub("_f$","",other), "_o")
        if (oo %in% names(d)) g[[oo]] <- ordered(lv[1], levels = lv)
      }
    }
    g
  }))

  ep <- tryCatch(brms::posterior_epred(fit, newdata = grid, re_formula = NA,
                                       ndraws = 500),
                 error = function(e) { log_step("epred ERR: %s", e$message); NULL })
  if (is.null(ep)) return(invisible())
  if (is_dur) ep <- ep / stats::median(d$N_readings, na.rm = TRUE)  # beta-binom count -> proportion
  grid$est <- apply(ep, 2, mean)
  grid$lo  <- apply(ep, 2, stats::quantile, 0.025)
  grid$hi  <- apply(ep, 2, stats::quantile, 0.975)

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = Hour_block, y = est,
                                          color = .data[[group_var]],
                                          fill = .data[[group_var]])) +
    ggplot2::annotate("rect", xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf,
                      alpha = 0.06, fill = pal$night) +
    ggplot2::annotate("rect", xmin = 20, xmax = 24, ymin = -Inf, ymax = Inf,
                      alpha = 0.06, fill = pal$night) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = group_pal) +
    ggplot2::scale_fill_manual(values = group_pal) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 4)) +
    ggplot2::labs(x = "Solar time",
                  y = "Predicted activity (model scale)",
                  color = legend_lab, fill = legend_lab,
                  title = "Diel activity curve", subtitle = subtitle) +
    theme_bear()
  save_fig(p, fig_name, w = 8, h = 5)
  log_step("Drawn: %s", fig_name)
}

tri_pal <- c(pal$hyperphagia, pal$mating, pal$hibernation)

# H3 — dump proximity
plot_diel_by(model_fp("H3"), "DumpProx_f",
             c("Near","Mid","Far"), tri_pal, tag_out("FigV3_03_H3_DielByDumpProx"),
             "By garbage-dump proximity (temperature-controlled)",
             legend_lab = "Dump proximity")
# H3 — road proximity
plot_diel_by(model_fp("H3"), "RoadProx_f",
             c("Near","Mid","Far"), tri_pal, tag_out("FigV3_04_H3_DielByRoadProx"),
             "By road proximity (temperature-controlled)",
             legend_lab = "Road proximity")
# H4 — thermal regime
plot_diel_by(model_fp("H4"), "TempTert_f",
             c("Cool","Mild","Warm"), c(pal$mating, pal$day, pal$gradient_high),
             tag_out("FigV3_05_H4_DielByTemp"), "By thermal regime (heat avoidance)",
             legend_lab = "Thermal regime")

log_step("=== 10b_visualization_newmodels DONE ===")
sink(type = "message"); sink()
