# ==============================================================================
# 11b_nocturnality_tertile.R  —  DIAGNOSTIC: nocturnality index by dump TERTILE
# ==============================================================================
#   This script refits the index with a TERTILE and produces:
#     (a) the true (possibly non-monotonic) Near/Mid/Far pattern,
#     (b) the MARGINAL night-fraction per tertile (G-comp) — a scale directly
#         comparable with the H3 night-fraction,
#   so we can tell "linear artefact" from "genuine reversed finding".
#
# SAME as script 11: strict night (nocturnal only), bear-day aggregation,
#   offset(log(night_len)), s(DOY), same covariates, beta_binomial, temp_mean
#   (the correct scale at the bear-day level — same as the 11:117 note).
# ONLY DIFFERENCE: d2GarbageDump_km_sc (linear) -> DumpProx_f (Near/Mid/Far tertile).
#
# OUTPUTS:
#   tables/Table_V3_08b_NoctTertile_Coef.csv       -> Table S29 (SI) (Near ref; Mid/Far contrasts)
#   tables/Table_V3_08b_NoctTertile_Marginal.csv   (per-tertile marginal night-fraction; supports S29)
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("11b_nocturnality_tertile")
log_step("=== 11b_nocturnality_tertile: START ===")

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE)
options(brms.backend = "rstan")

# ---- DATA + NIGHT DEFINITION (identical to script 11) ------------------------
d <- readRDS(file.path(dat_path, "hourly_decorr.rds"))
d <- dplyr::filter(d, Season %in% c("Mating", "Hyperphagia"))
d$active_reads <- d$Active_n
d$is_night <- d$t_period == "nocturnal"     # STRICT night — script 11 primary
log_step("is_night fraction (strict): %.3f", mean(d$is_night))

cov_sc <- intersect(c("d2GarbageDump_km_sc","d2Roads_km_sc","d2ProtectedAreas_km_sc",
                      "LandPC1_sc","LandPC2_sc","temp_mean_sc","Age_sc"), names(d))
noct <- d %>%
  dplyr::group_by(BearID_f, Year_f, Season_f, Date_Only) %>%
  dplyr::summarise(
    Total_active = sum(active_reads, na.rm = TRUE),
    Night_active = sum(active_reads[is_night], na.rm = TRUE),
    Sex_f = dplyr::first(Sex_f),
    night_len_h = dplyr::first(night_len_h),
    dplyr::across(dplyr::all_of(cov_sc), ~ mean(.x, na.rm = TRUE)),
    n_hours = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(Total_active >= 5, n_hours >= 12)
noct$DOY <- as.integer(lubridate::yday(noct$Date_Only))
noct$Total_active <- as.integer(noct$Total_active)
noct$Night_active <- as.integer(pmin(noct$Night_active, noct$Total_active))

# ---- DUMP TERTILE: tertile of the bear-day mean distance ---------------------
# The H3 tertiles are at the observation level; here we cut the bear-day MEAN
# distance by its own tertile (Near = smallest distance = close to the dump). Same
# label scheme.
qd <- stats::quantile(noct$d2GarbageDump_km_sc, c(0, 1/3, 2/3, 1), na.rm = TRUE)
qd[1] <- -Inf; qd[length(qd)] <- Inf
noct$DumpProx_f <- cut(noct$d2GarbageDump_km_sc, breaks = unique(qd),
                       labels = c("Near","Mid","Far"), include.lowest = TRUE)
log_step("Dump tertile distribution (bear-day): %s",
         paste(names(table(noct$DumpProx_f)), as.integer(table(noct$DumpProx_f)),
               sep = "=", collapse = " | "))
log_step("Tertile distance ranges (z): Near<=%.3f | Mid<=%.3f | Far",
         qd[2], qd[3])
log_step("Nocturnality data: %d bear-days, %d bears, mean night-fraction=%.3f",
         nrow(noct), dplyr::n_distinct(noct$BearID_f),
         mean(noct$Night_active / noct$Total_active))

# ---- MODEL: dump TERTILE, rest identical to script 11 ------------------------
fixed_terms <- paste(c("DumpProx_f",
  intersect(c("d2Roads_km_sc","d2ProtectedAreas_km_sc","LandPC1_sc","LandPC2_sc",
              "temp_mean_sc","Sex_f","Age_sc","Season_f"), names(noct))),
  collapse = " + ")
fml <- brms::bf(as.formula(paste(
  "Night_active | trials(Total_active) ~", fixed_terms,
  "+ offset(log(night_len_h)) + s(DOY, bs = 'tp', k = 5) + (1 | BearID_f/Year_f)")))
log_step("FORMULA: %s", paste(deparse(fml$formula), collapse = " "))

priors_bb <- c(
  brms::prior(normal(0, 1),         class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(gamma(0.01, 0.01),    class = "phi"))

t0 <- Sys.time()
fit <- brms::brm(formula = fml, data = noct, family = brms::beta_binomial(),
  prior = priors_bb, iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 42, silent = 1, refresh = 250, init = 0,
  control = list(adapt_delta = 0.95, max_treedepth = 11),
  save_pars = brms::save_pars(all = TRUE))
log_step("fit time: %.1f min", as.numeric(difftime(Sys.time(), t0, "mins")))

s <- summary(fit); fe <- s$fixed
cat("\n=== NOCT TERTILE FIXED EFFECTS ===\n"); print(fe)
coef_df <- data.frame(
  Parameter = rownames(fe),
  Estimate = round(fe[,"Estimate"],4), Q2.5 = round(fe[,"l-95% CI"],4),
  Q97.5 = round(fe[,"u-95% CI"],4), Rhat = round(fe[,"Rhat"],4),
  ESS_bulk = round(fe[,"Bulk_ESS"]),
  Direction_sig = (sign(fe[,"l-95% CI"]) == sign(fe[,"u-95% CI"])))
# -> Table S29 (SI): nocturnality model, dump tertile
write.csv(coef_df, file.path(tbl_path, "Table_V3_08b_NoctTertile_Coef.csv"),
          row.names = FALSE)

# ---- MARGINAL night-fraction per tertile (G-comp) ----------------------------
# posterior_epred beta_binomial returns COUNTS (E[Night_active]); divide by
# Total_active to get a PROPORTION -> same scale as the script 32 night_fraction.
set.seed(42L)
base <- noct[sample.int(nrow(noct), min(2000L, nrow(noct))), , drop = FALSE]
lvls <- c("Near","Mid","Far")
marg <- lapply(lvls, function(lv) {
  nd <- base; nd$DumpProx_f <- factor(lv, levels = lvls)
  ep <- brms::posterior_epred(fit, newdata = nd, re_formula = NA, ndraws = 400)
  p  <- rowMeans(sweep(ep, 2, nd$Total_active, "/"))   # per-draw marginal proportion
  c(mean = mean(p), q025 = unname(stats::quantile(p, .025)),
    q975 = unname(stats::quantile(p, .975)))
})
marg_df <- data.frame(DumpProx = lvls, do.call(rbind, marg), row.names = NULL)
# Near - Far contrast
ndN <- base; ndN$DumpProx_f <- factor("Near", levels = lvls)
ndF <- base; ndF$DumpProx_f <- factor("Far",  levels = lvls)
pN <- rowMeans(sweep(brms::posterior_epred(fit, newdata=ndN, re_formula=NA, ndraws=400), 2, ndN$Total_active, "/"))
pF <- rowMeans(sweep(brms::posterior_epred(fit, newdata=ndF, re_formula=NA, ndraws=400), 2, ndF$Total_active, "/"))
dNF <- pN - pF
marg_df <- rbind(marg_df,
  data.frame(DumpProx = "Near - Far", mean = mean(dNF),
             q025 = unname(stats::quantile(dNF,.025)),
             q975 = unname(stats::quantile(dNF,.975))))
marg_df$pd <- c(rep(NA_real_, 3), max(mean(dNF>0), mean(dNF<0)))
write.csv(marg_df, file.path(tbl_path, "Table_V3_08b_NoctTertile_Marginal.csv"),
          row.names = FALSE)
cat("\n=== MARGINAL NIGHT-FRACTION (tertile) ===\n"); print(marg_df, row.names = FALSE)

log_step("VERDICT: Near marginal=%.3f, Far marginal=%.3f, Near-Far=%.4f [%.4f,%.4f]",
         marg_df$mean[1], marg_df$mean[3], mean(dNF),
         stats::quantile(dNF,.025), stats::quantile(dNF,.975))
log_step("  Near>Far -> ALIGNED with H3/29 (near=more nocturnal); Near<Far ->")
log_step("  aligned with the linear index (far=more nocturnal) -> genuine reversed finding.")
log_step("=== 11b_nocturnality_tertile DONE ===")
sink(type = "message"); sink()
