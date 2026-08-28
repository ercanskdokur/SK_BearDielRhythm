# ==============================================================================
# 14_individual_specialization.R  —  Individual plasticity / specialization
# ==============================================================================
# QUESTION: do individuals differ consistently in their RESPONSE to human pressure
#   (dump/road) (a behavioural syndrome/specialization), or are they all similar?
#   H5's random-slope BLUPs are a direct measure of this. It is also linked to the
#   consistency of the individual's actual dump-proximity use (from position data).
#
# OUTPUTS: Table_V3_14_IndivSlopes.csv (per-individual slopes = phenotypes)
#            -> Table S36 (intensity) / S37 (duration)
#          Table_V3_15_SlopeSD.csv     (group-level sd = among-individual variation)
#            -> Table S38 (intensity) / S39 (duration)
#          FigV3_07_ReactionNorm.{png,pdf}
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("14_individual_specialization")
log_step("=== 14_individual_specialization: START ===")
options(brms.backend = "rstan")

h5 <- load_model("H5")
fit <- h5$fit; d <- h5$data
fe  <- brms::fixef(fit)

# ---- (1) Per-individual random slopes (BLUP) = phenotypes -------------------
re <- brms::ranef(fit)$BearID_f   # [bear, stat, parameter]
pars <- dimnames(re)[[3]]
log_step("ranef parameters: %s", paste(pars, collapse = ", "))
get_par <- function(p) if (p %in% pars) re[, "Estimate", p] else rep(NA_real_, dim(re)[1])

slopes <- data.frame(
  BearID_f = dimnames(re)[[1]],
  re_Intercept = get_par("Intercept"),
  re_dump_slope = get_par("d2GarbageDump_km_sc"),
  re_road_slope = get_par("d2Roads_km_sc"))
# Total (fixed + individual) slope = the individual's actual response
slopes$total_dump_slope <- fe["d2GarbageDump_km_sc","Estimate"] + slopes$re_dump_slope
slopes$total_road_slope <- fe["d2Roads_km_sc","Estimate"] + slopes$re_road_slope

# ---- Position-based dump-proximity use + consistency -----------------------
if ("d2GarbageDump_km" %in% names(d)) {
  expo <- d %>% dplyr::group_by(BearID_f) %>%
    dplyr::summarise(dump_exposure_mean = mean(d2GarbageDump_km, na.rm = TRUE),
                     dump_exposure_cv   = stats::sd(d2GarbageDump_km, na.rm=TRUE) /
                                          mean(d2GarbageDump_km, na.rm = TRUE),
                     .groups = "drop")
  slopes <- merge(slopes, expo, by = "BearID_f", all.x = TRUE)
}
# -> Table S36 (intensity) / S37 (duration): individual reaction-norm slopes
write.csv(slopes, file.path(tbl_path, paste0(tag_out("Table_V3_14_IndivSlopes"), ".csv")), row.names = FALSE)
log_step("Individual slopes: %d bears", nrow(slopes))

# ---- (2) Group-level SD (is there among-individual variation?) --------------
vc <- brms::VarCorr(fit)$BearID_f$sd   # [parameter, stat]
sd_df <- data.frame(
  Parameter = rownames(vc),
  SD_Estimate = round(vc[,"Estimate"],4),
  SD_Q2.5 = round(vc[,"Q2.5"],4), SD_Q97.5 = round(vc[,"Q97.5"],4),
  Credibly_nonzero = vc[,"Q2.5"] > 0)
# -> Table S38 (intensity) / S39 (duration): random-slope standard deviations
write.csv(sd_df, file.path(tbl_path, paste0(tag_out("Table_V3_15_SlopeSD"), ".csv")), row.names = FALSE)
cat("\n=== AMONG-INDIVIDUAL SLOPE SD (plasticity) ===\n"); print(sd_df, row.names = FALSE)

# ---- (3) FIGURE: dump-response reaction norm (one line per bear) ------------
dump_rng <- range(d$d2GarbageDump_km_sc, na.rm = TRUE)
xx <- seq(dump_rng[1], dump_rng[2], length.out = 30)
norm_df <- do.call(rbind, lapply(seq_len(nrow(slopes)), function(i) {
  data.frame(BearID_f = slopes$BearID_f[i], x = xx,
             y = slopes$re_Intercept[i] + fe["Intercept","Estimate"] +
                 slopes$total_dump_slope[i] * xx)
}))
pop_df <- data.frame(x = xx, y = fe["Intercept","Estimate"] +
                       fe["d2GarbageDump_km_sc","Estimate"] * xx)
p <- ggplot2::ggplot(norm_df, ggplot2::aes(x = x, y = y, group = BearID_f)) +
  ggplot2::geom_line(alpha = 0.35, color = pal$mating, linewidth = 0.5) +
  ggplot2::geom_line(data = pop_df, ggplot2::aes(x = x, y = y),
                     inherit.aes = FALSE, color = pal$hyperphagia, linewidth = 1.3) +
  ggplot2::labs(x = "Distance to garbage dump (z-score)",
                y = "Latent activity (logit scale)",
                title = "Individual reaction norms to dump proximity",
                subtitle = "Thin lines = individual bears; thick = population mean (H5)") +
  theme_bear() +
  ggplot2::theme(plot.title    = ggplot2::element_text(size = 12),
                 plot.subtitle = ggplot2::element_text(size = 9.5, colour = "grey40"))
save_fig(p, tag_out("FigV3_07_ReactionNorm"), w = 7, h = 5)

log_step("=== 14_individual_specialization DONE ===")
cat("DONE 14\n")
sink(type = "message"); sink()
