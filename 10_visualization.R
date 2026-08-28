# ==============================================================================
# 10_visualization.R  —  Manuscript-ready figures (Fig 1-5)
# ==============================================================================
# Five main figures for the manuscript, from the new framework:
#   Fig 1: Diel rhythm by Sex x Season (descriptive — Fig03 style)
#   Fig 2: H2 estimated smooths — Hour x Sex x Season (model-based)
#   Fig 3: H1 interaction plot — Season x Dump distance (marginal effects)
#   Fig 4: Coefficient forest plot — H0 with P(direction)
#   Fig 5: Variance partition + Conservation scenarios
#
# OUTPUTS: output/figures/FigV2_M1_DielSexSeason.{png,pdf}
#                        FigV2_M2_H2_HourSmooths.{png,pdf}
#                        FigV2_M3_H1_DumpInteraction.{png,pdf}
#                        FigV2_M4_CoefficientForest.{png,pdf}
#                        FigV2_M5_VPC_Scenarios.{png,pdf}
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("10_visualization")
log_step("=== 10_visualization (manuscript figs): START ===")

# ---- DATA ----
hourly_full <- readRDS(file.path(dat_path, "hourly_decorr.rds"))
h0_h <- load_model("H0")
fit_h0_hourly <- h0_h$fit

fp_h1_h <- model_fp("H1")
fp_h2_h <- model_fp("H2")

# Sun reference
sun_times_ref <- suncalc::getSunlightTimes(
  date = seq(as.Date("2020-04-01"), as.Date("2020-11-30"), by = "day"),
  lat = study_lat, lon = study_lon, tz = "Europe/Istanbul",
  keep = c("sunrise", "sunset")
)
mean_sr <- mean(lubridate::hour(sun_times_ref$sunrise) +
                lubridate::minute(sun_times_ref$sunrise) / 60)
mean_ss <- mean(lubridate::hour(sun_times_ref$sunset) +
                lubridate::minute(sun_times_ref$sunset) / 60)

# ==============================================================================
# FIG 1 — Diel rhythm by Sex x Season (descriptive)
# ==============================================================================
log_step("Fig M1: Diel by Sex x Season (descriptive)")
fig1_data <- hourly_full %>%
  dplyr::filter(Season %in% c("Mating","Hyperphagia")) %>%
  dplyr::group_by(Sex, Season, Hour_block) %>%
  dplyr::summarise(Mean_act = mean(Activity, na.rm = TRUE),
                   SE = stats::sd(Activity, na.rm = TRUE) / sqrt(dplyr::n()),
                   .groups = "drop")

fig1 <- ggplot2::ggplot(fig1_data,
                        ggplot2::aes(x = Hour_block, y = Mean_act,
                                     color = Sex, fill = Sex)) +
  ggplot2::annotate("rect", xmin = 0, xmax = mean_sr,
                    ymin = -Inf, ymax = Inf,
                    fill = pal$night, alpha = 0.08) +
  ggplot2::annotate("rect", xmin = mean_ss, xmax = 24,
                    ymin = -Inf, ymax = Inf,
                    fill = pal$night, alpha = 0.08) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Mean_act - SE,
                                    ymax = Mean_act + SE),
                       alpha = 0.2, color = NA) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_vline(xintercept = c(mean_sr, mean_ss),
                      linetype = "dashed", color = "grey50") +
  ggplot2::facet_wrap(~ Season) +
  ggplot2::scale_color_manual(values = c(F = pal$female, M = pal$male),
                              labels = c(F = "Female", M = "Male")) +
  ggplot2::scale_fill_manual(values  = c(F = pal$female, M = pal$male),
                             labels = c(F = "Female", M = "Male")) +
  ggplot2::scale_x_continuous(breaks = seq(0, 23, 6),
                              labels = sprintf("%02d:00", seq(0, 23, 6))) +
  ggplot2::labs(x = "Hour of day", y = "Mean hourly activity intensity",
                title = "Diel activity rhythm by sex and season",
                subtitle = "Sarıkamış brown bears, NE Türkiye") +
  theme_bear()
save_fig(fig1, tag_out("FigV2_M1_DielSexSeason"), w = 12, h = 5)

# ==============================================================================
# FIG 2 — H2 estimated smooths (model-based)
# ==============================================================================
if (file.exists(fp_h2_h)) {
  log_step("Fig M2: H2 conditional smooths")
  h2 <- readRDS(fp_h2_h)
  fit_h2 <- h2$fit

  ce_h2 <- tryCatch(
    brms::conditional_smooths(fit_h2),  # get all smooths
    error = function(e) { log_step("H2 ce FAIL: %s", e$message); NULL }
  )

  if (!is.null(ce_h2) && length(ce_h2) > 0) {
    # Find the Hour_block smooth (the one with Sex_Season_f variable)
    sm <- NULL
    for (i in seq_along(ce_h2)) {
      cand <- ce_h2[[i]]
      if (!is.null(cand) && "Sex_Season_f" %in% names(cand) && "Hour_block" %in% names(cand)) {
        sm <- cand; break
      }
    }
    if (!is.null(sm)) {
      sm <- sm %>%
        dplyr::mutate(
          Sex = ifelse(grepl("^F_", as.character(Sex_Season_f)), "Female", "Male"),
          Season = ifelse(grepl("Mating", as.character(Sex_Season_f)),
                          "Mating", "Hyperphagia"),
          Season = factor(Season, levels = c("Mating","Hyperphagia"))
        )

      fig2 <- ggplot2::ggplot(sm,
                              ggplot2::aes(x = Hour_block, y = estimate__,
                                           color = Sex, fill = Sex)) +
        ggplot2::annotate("rect", xmin = 0, xmax = mean_sr,
                          ymin = -Inf, ymax = Inf,
                          fill = pal$night, alpha = 0.08) +
        ggplot2::annotate("rect", xmin = mean_ss, xmax = 24,
                          ymin = -Inf, ymax = Inf,
                          fill = pal$night, alpha = 0.08) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower__, ymax = upper__),
                             alpha = 0.25, color = NA) +
        ggplot2::geom_line(linewidth = 1.1) +
        ggplot2::facet_wrap(~ Season) +
        ggplot2::scale_color_manual(values = c(Female = pal$female, Male = pal$male)) +
        ggplot2::scale_fill_manual(values  = c(Female = pal$female, Male = pal$male)) +
        ggplot2::scale_x_continuous(breaks = seq(0, 23, 6),
                                    labels = sprintf("%02d:00", seq(0, 23, 6))) +
        ggplot2::geom_vline(xintercept = c(mean_sr, mean_ss),
                            linetype = "dashed", color = "grey50") +
        ggplot2::labs(x = "Hour of day", y = "Estimated smooth (95% CrI)",
                      title = "H2: Diel rhythm smooth by sex x season",
                      subtitle = "brms zero_inflated_beta | s(Hour, by=Sex_Season, cc, k=12)") +
        theme_bear()
      save_fig(fig2, tag_out("FigV2_M2_H2_HourSmooths"), w = 12, h = 5)
    }
  }
}

# ==============================================================================
# FIG 3 — H1 Season x Dump interaction
# ==============================================================================
if (file.exists(fp_h1_h)) {
  log_step("Fig M3: H1 Season x Dump")
  h1 <- readRDS(fp_h1_h)
  fit_h1 <- h1$fit

  ce_h1 <- tryCatch(
    brms::conditional_effects(fit_h1,
                              effects = "d2GarbageDump_km_sc:Season_f",
                              plot = FALSE),
    error = function(e) { log_step("H1 ce FAIL: %s", e$message); NULL }
  )

  if (!is.null(ce_h1) && length(ce_h1) > 0) {
    d <- ce_h1[[1]]
    fig3 <- ggplot2::ggplot(d,
                            ggplot2::aes(x = d2GarbageDump_km_sc,
                                         y = estimate__,
                                         color = Season_f, fill = Season_f)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower__, ymax = upper__),
                           alpha = 0.25, color = NA) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::scale_color_manual(values = c(Mating = pal$mating,
                                              Hyperphagia = pal$hyperphagia)) +
      ggplot2::scale_fill_manual(values  = c(Mating = pal$mating,
                                              Hyperphagia = pal$hyperphagia)) +
      ggplot2::labs(x = "Distance to garbage dump (z-scored)",
                    y = "Predicted hourly activity",
                    title = "H1: Season x dump distance interaction",
                    subtitle = "Predicted from H1 model",
                    color = "Season", fill = "Season") +
      theme_bear()
    save_fig(fig3, tag_out("FigV2_M3_H1_DumpInteraction"), w = 10, h = 6)
  }
}

# ==============================================================================
# FIG 4 — Coefficient forest (H0_hourly)
# ==============================================================================
log_step("Fig M4: H0 coefficient forest")
fe <- summary(fit_h0_hourly)$fixed
draws <- brms::as_draws_matrix(fit_h0_hourly, variable = "^b_", regex = TRUE)
p_gt_0 <- colMeans(draws > 0)

fe_df <- data.frame(
  Parameter = rownames(fe),
  Estimate  = fe[, "Estimate"],
  Q2.5      = fe[, "l-95% CI"],
  Q97.5     = fe[, "u-95% CI"],
  stringsAsFactors = FALSE
)
fe_df$P_gt_0 <- p_gt_0[paste0("b_", fe_df$Parameter)]
fe_df$Strong <- (fe_df$P_gt_0 > 0.95 | fe_df$P_gt_0 < 0.05)
fe_df$Direction <- ifelse(fe_df$Estimate > 0, "Positive", "Negative")
fe_df <- fe_df %>%
  dplyr::filter(Parameter != "Intercept") %>%
  dplyr::mutate(Parameter = factor(Parameter,
                                   levels = Parameter[order(Estimate)]))

fig4 <- ggplot2::ggplot(fe_df,
                        ggplot2::aes(x = Estimate, y = Parameter,
                                     xmin = Q2.5, xmax = Q97.5,
                                     color = Strong)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_pointrange(size = 0.5, linewidth = 0.8) +
  ggplot2::geom_text(ggplot2::aes(x = Q97.5, label = sprintf("%.2f", P_gt_0)),
                     hjust = -0.3, size = 3, color = "grey30") +
  ggplot2::scale_color_manual(values = c(`TRUE` = pal$female, `FALSE` = "grey60"),
                              labels = c(`TRUE` = "P(direction) > 0.95",
                                         `FALSE` = "P(direction) <= 0.95")) +
  ggplot2::labs(x = "Posterior coefficient (z-scored predictor)",
                y = NULL,
                title = "H0 reference model: posterior coefficients",
                subtitle = sprintf("Hourly_Intensity | zero_inflated_beta | n=%d | numbers = P(coef > 0)",
                                   h0_h$n),
                color = NULL) +
  ggplot2::xlim(NA, max(fe_df$Q97.5) * 1.2) +
  theme_bear() +
  ggplot2::theme(legend.position = "bottom")
save_fig(fig4, tag_out("FigV2_M4_CoefficientForest"), w = 10, h = 7)

# ==============================================================================
# FIG 5 — Variance partition + Dump scenarios
# ==============================================================================
log_step("Fig M5: VPC + Dump scenarios")


vpc_fp  <- file.path(tbl_path, paste0(tag_out("Table_V3_03_VPC_latent"), ".csv"))
scen_fp <- file.path(tbl_path, paste0(tag_out("Table_V2_13_DumpScenarios"), ".csv"))

if (file.exists(vpc_fp)) {
  vpc <- read.csv(vpc_fp)
  # 07b format: Model | Family | Group | Var | Var_dist | VPC_latent_pct | VPC_pi2_3_OLD
  # The VPC_pi2_3_OLD column is NOT plotted — it is in the table for transparency only.
  vpc <- vpc[, intersect(c("Model","Group","VPC_latent_pct"), names(vpc))]

  p5a <- ggplot2::ggplot(vpc,
                         ggplot2::aes(x = Model, y = VPC_latent_pct, fill = Group)) +
    ggplot2::geom_col(position = "stack", color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = setNames(pal$multi[1:length(unique(vpc$Group))],
                                                 unique(vpc$Group))) +
    ggplot2::labs(x = NULL, y = "Variance partition (%)",
                  title = "(A) Latent-scale variance partition",
                  subtitle = paste0("Family-specific observation-level variance",
                                    " | ", RESP_KEY())) +
    theme_bear() +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (file.exists(scen_fp)) {
    scen <- read.csv(scen_fp)
    p5b <- ggplot2::ggplot(scen,
                           ggplot2::aes(x = Dump_quantile, y = Predicted_Activity_mean,
                                        ymin = Q2.5, ymax = Q97.5,
                                        color = Season, group = Season)) +
      ggplot2::geom_pointrange(position = ggplot2::position_dodge(0.3),
                               size = 0.6) +
      ggplot2::geom_line(position = ggplot2::position_dodge(0.3),
                         linewidth = 0.9) +
      ggplot2::scale_color_manual(values = c(Mating = pal$mating,
                                              Hyperphagia = pal$hyperphagia)) +
      ggplot2::labs(x = "Dump distance quantile (z-scored)",
                    y = "Predicted hourly activity",
                    title = "(B) Predicted activity under dump-distance scenarios",
                    subtitle = "Season x dump (H1_hourly)") +
      theme_bear()

    fig5 <- patchwork::wrap_plots(p5a, p5b, ncol = 2, widths = c(1, 1))
    save_fig(fig5, tag_out("FigV2_M5_VPC_Scenarios"), w = 14, h = 6)
  } else {
    save_fig(p5a, tag_out("FigV2_M5_VPC_Scenarios"), w = 8, h = 6)
  }
}

log_step("=== 10_visualization DONE ===")
sink(type = "message"); sink()
