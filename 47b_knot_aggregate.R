# ==============================================================================
# 47b_knot_aggregate.R  —  Knot-sensitivity: aggregate the 8 fits, build a table
# ==============================================================================
# Runs AFTER the 47a array finishes. For each k:
#   (1) H3-vs-H4 LOO elpd difference (same rows -> paired)
#   (2) H3 Near-Far night-fraction contrast (G-computation; script 32 method)
# OUTPUT: tables/Table_V3_25_KnotSensitivity_intensity.csv  -> Table S61 (SI)
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
suppressWarnings(suppressMessages({library(dplyr); library(loo)}))

init_log("47b_knot_aggregate")
log_step("=== 47b_knot_aggregate | intensity ===")
slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))

# Rebuild d18 with the SAME seed (for G-comp) + night hours
md <- readRDS(file.path(mod_path, "model_data_50k.rds"))
d_full <- md$data; TIME_VAR <- get_time_var()
set.seed(42L); TARGET_N <- 18000L
d18 <- d_full %>% dplyr::group_by(BearID_f) %>%
  dplyr::slice_sample(prop = TARGET_N / nrow(d_full)) %>% dplyr::ungroup()
nt <- tapply(d18$t_period == "nocturnal", round(d18[[TIME_VAR]]), mean, na.rm = TRUE)
NIGHT_HOURS <- as.numeric(names(nt))[which(nt > 0.5)]
log_step("d18=%d rows | night hours: %s", nrow(d18), paste(NIGHT_HOURS, collapse=","))

gcomp_nightfrac <- function(fit, n_rep = 1000L, n_draws = 300L) {
  set.seed(42L)
  base  <- d18[sample.int(nrow(d18), min(n_rep, nrow(d18))), , drop = FALSE]
  hours <- seq(0, 24, length.out = 49)[-49]; dt <- 24 / length(hours)
  lvls  <- c("Near","Far")
  cube  <- array(NA_real_, c(n_draws, length(hours), 2), dimnames = list(NULL,NULL,lvls))
  for (hi in seq_along(hours)) {
    nd <- do.call(rbind, lapply(lvls, function(lv) {
      x <- base; x[[TIME_VAR]] <- hours[hi]
      x$DumpProx_f <- lv
      x$DumpProx_o <- factor(lv, levels = c("Near","Mid","Far"), ordered = TRUE)
      x$.lvl <- lv; x
    }))
    ep <- brms::posterior_epred(fit, newdata = nd, re_formula = NA, ndraws = n_draws)
    for (li in 1:2) cube[, hi, li] <- rowMeans(ep[, nd$.lvl == lvls[li], drop = FALSE])
  }
  auc <- apply(cube, c(1,3), sum) * dt
  ni  <- which(round(hours) %in% NIGHT_HOURS)
  nf  <- sapply(1:2, function(li) (apply(cube[, ni, li, drop=FALSE], 1, sum) * dt) / auc[, li])
  dN  <- nf[,1] - nf[,2]
  c(nf_near = mean(nf[,1]), nf_far = mean(nf[,2]), contrast = mean(dN),
    lo = unname(stats::quantile(dN,.025)), hi = unname(stats::quantile(dN,.975)),
    pd = max(mean(dN>0), mean(dN<0)))
}

Kset <- c(8L,10L,12L,15L)
rows <- list()
for (K in Kset) {
  f3 <- readRDS(file.path(mod_path, sprintf("knotsens_H3_k%d_intensity.rds", K)))$fit
  f4 <- readRDS(file.path(mod_path, sprintf("knotsens_H4_k%d_intensity.rds", K)))$fit
  lc <- loo::loo_compare(f3$criteria$loo, f4$criteria$loo)
  e3 <- f3$criteria$loo$estimates["elpd_loo","Estimate"]
  e4 <- f4$criteria$loo$estimates["elpd_loo","Estimate"]
  nf <- gcomp_nightfrac(f3)
  rows[[as.character(K)]] <- data.frame(
    k = K, elpd_H3 = e3, elpd_H4 = e4,
    elpd_diff_H3_minus_H4 = e3 - e4, se_diff = lc[2,"se_diff"],
    nightfrac_near = nf["nf_near"], nightfrac_far = nf["nf_far"],
    nightfrac_contrast_NearMinusFar = nf["contrast"],
    contrast_lo = nf["lo"], contrast_hi = nf["hi"], contrast_pd = nf["pd"],
    row.names = NULL)
  log_step("k=%d | elpd H3=%.1f H4=%.1f diff=%.1f (se=%.1f) | night Near-Far=%.4f [%.4f,%.4f] pd=%.3f",
           K, e3, e4, e3 - e4, lc[2,"se_diff"], nf["contrast"], nf["lo"], nf["hi"], nf["pd"])
}
out <- do.call(rbind, rows)
out$best_model <- ifelse(out$elpd_H3 >= out$elpd_H4, "H3", "H4")  # definitive: higher elpd
out$response <- "intensity"; out$n_rows <- nrow(d18); out$time_axis <- TIME_VAR
# -> Table S61 (SI): diel-spline knot-count sensitivity
write.csv(out, file.path(tbl_path, "Table_V3_25_KnotSensitivity_intensity.csv"), row.names = FALSE)
cat("\n=== KNOT SENSITIVITY SUMMARY ===\n"); print(out, row.names = FALSE)
log_step("=== 47b DONE -> Table_V3_25_KnotSensitivity_intensity.csv ===")
sink(type = "message"); sink()
