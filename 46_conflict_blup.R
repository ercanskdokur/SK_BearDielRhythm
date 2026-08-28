# ==============================================================================
# 46_conflict_blup.R  —  Per-bear conflict overlap: PARTIAL POOLING (BLUP)  [B8]
# ==============================================================================
# Call:  Rscript 46_conflict_blup.R
#
#
# OUTPUTS: tables/Table_V3_29_ConflictOverlap_BLUP.csv  -> Table S56 (SI)
#          figures/FigV3_11_ConflictOverlap_BLUP.{png,pdf}
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("46_conflict_blup")
log_step("=== 46_conflict_blup: START ===")
options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
rstan::rstan_options(auto_write = TRUE); options(brms.backend = "rstan")

d <- readRDS(file.path(dat_path, "hourly_decorr.rds"))
d <- dplyr::filter(d, Season %in% c("Mating", "Hyperphagia"))
log_step("hourly_decorr n=%d | cols: %s", nrow(d), paste(head(names(d),25), collapse=", "))

# Window flags (identical to 15)
d$in_conflict <- (d$Hour_mid >= (d$SR - 1) & d$Hour_mid < (d$SR + 4)) |
                 (d$Hour_mid >= (d$SS - 4) & d$Hour_mid < (d$SS + 1))

# Active-reading COUNTS (beta-binomial needs integers)
if (all(c("Active_n","N_readings") %in% names(d))) {
  d$an <- d$Active_n; d$nn <- d$N_readings
} else {
  nn <- if ("N_readings" %in% names(d)) d$N_readings else 12L
  d$nn <- nn; d$an <- round(d$Pct_active * nn)
}
d$an <- pmin(d$an, d$nn)

perbear <- d %>% dplyr::group_by(BearID_f) %>%
  dplyr::summarise(active_in = sum(an[in_conflict], na.rm = TRUE),
                   active_tot = sum(an, na.rm = TRUE),
                   n_obs = dplyr::n(), .groups = "drop")
perbear <- perbear[perbear$active_tot > 0, ]
perbear$active_in <- pmin(perbear$active_in, perbear$active_tot)
log_step("BLUP data: %d bears | active_tot range %d..%d",
         nrow(perbear), min(perbear$active_tot), max(perbear$active_tot))

fit <- brms::brm(
  active_in | trials(active_tot) ~ 1 + (1 | BearID_f),
  data = perbear, family = brms::beta_binomial(),
  prior = c(brms::prior(student_t(3,0,2.5), class="Intercept"),
            brms::prior(student_t(3,0,2.5), class="sd"),
            brms::prior(gamma(0.01,0.01), class="phi")),
  iter = 3000, warmup = 1000, chains = 4, cores = 4, seed = 42,
  silent = 1, refresh = 0, init = 0, control = list(adapt_delta = 0.95))

# Per-bear partially-pooled overlap (trials=1 -> proportion)
nd <- data.frame(BearID_f = perbear$BearID_f, active_tot = 1)
ep <- brms::posterior_epred(fit, newdata = nd, re_formula = NULL)
raw_overlap <- perbear$active_in / perbear$active_tot
out <- data.frame(
  BearID_f = perbear$BearID_f, n_obs = perbear$n_obs,
  active_tot = perbear$active_tot,
  raw_overlap = round(raw_overlap, 4),
  blup_overlap = round(apply(ep, 2, mean), 4),
  blup_lo = round(apply(ep, 2, stats::quantile, 0.025), 4),
  blup_hi = round(apply(ep, 2, stats::quantile, 0.975), 4))
out <- out[order(out$blup_overlap), ]
pop_mean <- round(brms::posterior_epred(fit,
              newdata = data.frame(active_tot = 1), re_formula = NA) |> mean(), 4)
log_step("Population mean overlap=%.4f | raw range %.3f..%.3f -> BLUP range %.3f..%.3f",
         pop_mean, min(raw_overlap), max(raw_overlap), min(out$blup_overlap), max(out$blup_overlap))
# -> Table S56 (SI): per-bear conflict overlap, partially pooled (BLUP)
utils::write.csv(out, file.path(tbl_path, "Table_V3_29_ConflictOverlap_BLUP.csv"), row.names = FALSE)
cat("\n=== CONFLICT OVERLAP (raw vs partially-pooled BLUP) ===\n"); print(out, row.names = FALSE)

# Caterpillar
out$BearID_f <- factor(out$BearID_f, levels = out$BearID_f)
p <- ggplot2::ggplot(out, ggplot2::aes(x = blup_overlap, y = BearID_f)) +
  ggplot2::geom_vline(xintercept = pop_mean, linetype = "dashed", color = "grey50") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = blup_lo, xmax = blup_hi), height = 0, color = "grey60") +
  ggplot2::geom_point(ggplot2::aes(size = n_obs), color = pal$mating) +
  ggplot2::geom_point(ggplot2::aes(x = raw_overlap), shape = 4, color = "grey40", alpha = 0.7) +
  ggplot2::scale_size_continuous(range = c(1, 4), name = "bear-hours") +
  ggplot2::labs(x = "Conflict-window overlap (partially pooled)", y = NULL,
                title = "Per-bear conflict-window overlap, shrunk",
                subtitle = "Point = BLUP (size ~ n); x = raw mean; dashed = population mean") +
  theme_bear()
save_fig(p, "FigV3_11_ConflictOverlap_BLUP", w = 7.5, h = 8)

log_step("=== 46_conflict_blup DONE ===")
sink(type = "message"); sink()
