# ==============================================================================
# 29_dump_partitioning.R  —  INTRA-SPECIFIC temporal partitioning among dump users
# ==============================================================================
#
# QUESTION: do bears using the dump separate their timing to AVOID one another?
#   In brown bears, dominant males may exclude females/subadults both spatially
#   and TEMPORALLY. If so, the dump is not just a "resource" but a stage for
#   intra-specific competition -> H3's population-level curve hides a heterogeneous
#   average.
#
# OUTPUTS:
#   tables/Table_V3_19_DumpPartition_smooth_<RESP>.csv   (hourly prediction per class)
#     -> Table S40 (intensity) / S41 (duration)
#   tables/Table_V3_20_DumpPartition_summary_<RESP>.csv  (peak/amplitude/night-fraction)
#     -> Table S42 (intensity) / S43 (duration)
#   figures/FigV3_08_DumpPartitioning_<RESP>.{png,pdf}
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))

RS <- get_response_spec()
init_log(sprintf("29_dump_partitioning_%s", RS$key))
log_step("=== 29_dump_partitioning | RESP=%s ===", RS$key)

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE)
options(brms.backend = "rstan")

NEAR_KM <- as.numeric(Sys.getenv("NEAR_KM", "2"))   # the old version also used 2 km

# ---- DATA: shared model data -------------------------------------------------
md <- readRDS(file.path(mod_path, "model_data_50k.rds"))
d  <- md$data
TIME_VAR <- md$time_var
KNOTS    <- md$knots

dump_col <- if ("d2GarbageDump_km" %in% names(d)) "d2GarbageDump_km" else NULL
if (is.null(dump_col)) stop("d2GarbageDump_km missing — the raw km column is required (for the threshold)")

n0 <- nrow(d)
d <- d[d[[dump_col]] < NEAR_KM, , drop = FALSE]
log_step("Dump-near filter (<%.1f km): %s -> %s hours, %d bears",
         NEAR_KM, format(n0, big.mark = ","), format(nrow(d), big.mark = ","),
         dplyr::n_distinct(d$BearID))
if (nrow(d) < 3000) log_step("!! WARNING: low n (%d) — power limited", nrow(d))

# REPORT the demographic distribution (this was the old version's weak point)
demo <- d %>%
  dplyr::distinct(BearID, Sex_f, Age_Category) %>%
  dplyr::count(Sex_f, Age_Category)
cat("\n=== DUMP-NEAR BEAR COMPOSITION ===\n"); print(demo)
log_step("Age range: %.1f - %.1f (used continuously, NOT classified)",
         min(d$Current_Age, na.rm = TRUE), max(d$Current_Age, na.rm = TRUE))

d$Sex_o <- ordered(d$Sex_f, levels = levels(d$Sex_f))

# ---- FORMULA -----------------------------------------------------------------
# NOTE: brms does NOT support ti() (only t2). The earlier `ti(...)` gave a terms_sm
# error. `t2(..., full = TRUE)` is the interaction component of a main-effect-
# inclusive tensor, orthogonalised against the existing s(SolarHour)+s(Age_sc) ->
# the "does age change the diel SHAPE" question is preserved and runs under brms.
fml_str <- paste(
  RS$lhs, "~",
  sprintf("s(%s, bs = 'cc', k = 12)", TIME_VAR),                       # common diel curve
  sprintf("+ s(%s, by = Sex_o, bs = 'cc', k = 12)", TIME_VAR),         # sex DEVIATION
  "+ s(Age_sc, bs = 'tp', k = 4)",                                     # age main effect
  sprintf("+ t2(%s, Age_sc, bs = c('cc','tp'), k = c(12, 4), full = TRUE)", TIME_VAR), # age x diel
  "+ Sex_f + Season_f + s(DOY, bs = 'tp', k = 5)",
  "+ (1 | BearID_f/Year_f)"
)
log_step("FORMULA: %s", fml_str)

bform <- if (RS$has_zi) {
  brms::bf(as.formula(fml_str),
           as.formula(sprintf("zi ~ s(%s, bs = 'cc', k = 12)", TIME_VAR)))
} else brms::bf(as.formula(fml_str))

saved_fp <- file.path(mod_path, sprintf("DumpPartition_%s.rds", RS$key))
if (file.exists(saved_fp) && !nzchar(Sys.getenv("FORCE_REFIT"))) {
  log_step("Reusing saved fit (no refit; force with FORCE_REFIT): %s",
           basename(saved_fp))
  .obj <- readRDS(saved_fp); fit <- .obj$fit
  if (!is.null(.obj$data))     d        <- .obj$data
  if (!is.null(.obj$time_var)) TIME_VAR <- .obj$time_var
  if (!is.null(.obj$near_km))  NEAR_KM  <- .obj$near_km
} else {
  t0 <- Sys.time()
  fit <- brms::brm(formula = bform, data = d, family = RS$family, prior = RS$priors,
                   knots = KNOTS, iter = 3000, warmup = 1000, chains = 4, cores = 4,
                   seed = 42, silent = 1, refresh = 250, init = 0,
                   control = list(adapt_delta = 0.97, max_treedepth = 14),
                   save_pars = brms::save_pars(all = TRUE))
  log_step("fit time: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  save_rds_safe(list(fit = fit, response = RS$key, near_km = NEAR_KM,
                     formula = fml_str, time_var = TIME_VAR, data = d), saved_fp)
}

s <- summary(fit); cat("\n=== FIXED ===\n"); print(s$fixed)
if (!is.null(s$splines)) { cat("\n=== SPLINES ===\n"); print(s$splines) }

# ==============================================================================
ok <- c("#D55E00", "#0072B2")   # Okabe-Ito (F, M) — colour-blind friendly
sexlv <- levels(d$Sex_f)

# ---- COMPOSITION (REQUIRED for the caption) ----------------------------------
nb <- dplyr::n_distinct(d$BearID)
by <- dplyr::n_distinct(paste(d$BearID, d$Year))
sex_n <- d %>% dplyr::distinct(BearID, Sex_f) %>% dplyr::count(Sex_f)
nF <- sum(sex_n$n[sex_n$Sex_f == sexlv[1]]); nM <- sum(sex_n$n[sex_n$Sex_f == sexlv[2]])
age_rng <- range(d$Current_Age, na.rm = TRUE)
log_step("COMPOSITION (<%.0f km): %d individuals (%d F, %d M), %d bear-years, %d bear-hours, age %.0f-%.0f y",
         NEAR_KM, nb, nF, nM, by, nrow(d), age_rng[1], age_rng[2])

# year<->Age_sc linear mapping (derive from data; do not trust the scale attribute)
age_map   <- stats::lm(Age_sc ~ Current_Age, data = d)
age_sc_of <- function(yr) as.numeric(stats::predict(age_map, newdata = data.frame(Current_Age = yr)))
hours <- seq(0, 24, length.out = 97)[-97]

# night hours FROM THE DATA (t_period) — no fixed assumption
night_hours <- if ("t_period" %in% names(d)) {
  tb <- tapply(d$t_period == "nocturnal", round(d[[TIME_VAR]]), mean, na.rm = TRUE)
  as.numeric(names(tb))[which(tb > 0.5)]
} else c(21:23, 0:4)
log_step("Night hours (%s, from t_period): %s", TIME_VAR, paste(night_hours, collapse = ","))

pred_grid <- function(age_vec) {
  g <- expand.grid(t = hours, Sex_f = sexlv, Age_year = age_vec, stringsAsFactors = FALSE)
  g[[TIME_VAR]] <- g$t
  g$Sex_o   <- ordered(g$Sex_f, levels = sexlv)
  g$Age_sc  <- age_sc_of(g$Age_year)
  g$Season_f <- levels(d$Season_f)[which.max(table(d$Season_f))]
  g$DOY <- stats::median(d$DOY, na.rm = TRUE)
  if (!RS$has_zi) g$N_readings <- stats::median(d$N_readings, na.rm = TRUE)
  g
}
epred_prop <- function(g) {
  e <- brms::posterior_epred(fit, newdata = g, re_formula = NA, ndraws = 800)
  if (!RS$has_zi) e <- e / stats::median(d$N_readings, na.rm = TRUE)   # to proportion
  e
}

# ---- (a) TOP: age-gradient curves (full-year sequence) -----------------------
age_years <- seq(ceiling(age_rng[1]), floor(age_rng[2]), by = 1)
gc_ <- pred_grid(age_years); ec_ <- epred_prop(gc_)
curve_df <- data.frame(hour = gc_$t, Sex = gc_$Sex_f, Age = gc_$Age_year,
                       est = colMeans(ec_),
                       lo = apply(ec_, 2, stats::quantile, 0.025),
                       hi = apply(ec_, 2, stats::quantile, 0.975))
# -> Table S40 (intensity) / S41 (duration): dump-partitioning diel curves by age x sex
write.csv(curve_df,
          file.path(tbl_path, sprintf("Table_V3_19_DumpPartition_smooth_%s.csv", RS$key)),
          row.names = FALSE)

# ---- (b,c) BOTTOM: amplitude & night-fraction, CONTINUOUS function of age (CrI) --
age_fine <- seq(age_rng[1], age_rng[2], length.out = 40)
gp_ <- pred_grid(age_fine); ep_ <- epred_prop(gp_)
prof <- do.call(rbind, lapply(
  split(seq_len(nrow(gp_)), list(gp_$Sex_f, gp_$Age_year)), function(ix) {
    sub <- ep_[, ix, drop = FALSE]; h <- gp_$t[ix]         # draws x hour
    amp <- apply(sub, 1, function(r) max(r) - min(r))
    nf  <- apply(sub, 1, function(r) sum(r[round(h) %in% night_hours]) / sum(r))
    pk  <- h[which.max(colMeans(sub))]
    data.frame(Sex = gp_$Sex_f[ix][1], Age = gp_$Age_year[ix][1],
               peak_hour_solar = round(pk, 2),
               amp = mean(amp), amp_lo = stats::quantile(amp, .025), amp_hi = stats::quantile(amp, .975),
               nf  = mean(nf),  nf_lo  = stats::quantile(nf, .025),  nf_hi  = stats::quantile(nf, .975))
  }))
rownames(prof) <- NULL
prof$response <- RS$key; prof$near_km <- NEAR_KM; prof$time_axis <- TIME_VAR
# -> Table S42 (intensity) / S43 (duration): dump-partitioning age profile
write.csv(prof,
          file.path(tbl_path, sprintf("Table_V3_20_DumpPartition_summary_%s.csv", RS$key)),
          row.names = FALSE)

# text anchor: amplitude & night-fraction at the youngest / median / oldest age
anch <- prof[prof$Age %in% c(min(prof$Age), stats::median(prof$Age), max(prof$Age)), ]
cat("\n=== AGE PROFILE (three age anchors) ===\n")
print(anch[, c("Sex","Age","peak_hour_solar","amp","nf")], row.names = FALSE)

# sex contrast at the median age: are males and females active at the SAME hour?
if (length(sexlv) == 2) {
  amed <- age_years[which.min(abs(age_years - stats::median(d$Current_Age, na.rm = TRUE)))]
  gm <- pred_grid(amed); em <- epred_prop(gm)
  iF <- which(gm$Sex_f == sexlv[1]); iM <- which(gm$Sex_f == sexlv[2])
  dif <- em[, iM, drop = FALSE] - em[, iF, drop = FALSE]
  pd  <- apply(dif, 2, function(x) max(mean(x > 0), mean(x < 0)))
  hrs_sig <- gm$t[iF][pd > 0.95]
  log_step("Sex difference (age=%d) pd>0.95 solar hours: %s", amed,
           if (length(hrs_sig)) paste(round(hrs_sig, 1), collapse = ", ") else "NONE")
  log_step("  -> %s", if (length(hrs_sig) >= 4) "EVIDENCE of temporal partitioning"
           else "NO clear temporal partitioning")
}

# ---- FIGURE (2 rows: curves / [amplitude | night-fraction]) ------------------
y_short <- trimws(sub("\\s*\\(.*$", "", RS$label))   # "Active duration" / "Movement intensity"
p_top <- ggplot2::ggplot(curve_df, ggplot2::aes(hour, est, group = Age, colour = Age)) +
  ggplot2::geom_line(linewidth = .55) +
  ggplot2::facet_wrap(~ Sex, nrow = 1) +
  ggplot2::scale_colour_viridis_c(option = "viridis", name = "Age (y)") +
  ggplot2::labs(x = sprintf("Solar time (%s; 0 = anchor)", TIME_VAR), y = y_short) +
  theme_bear() +
  ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 12))
p_amp <- ggplot2::ggplot(prof, ggplot2::aes(Age, amp, colour = Sex, fill = Sex)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = amp_lo, ymax = amp_hi), alpha = .15, colour = NA) +
  ggplot2::geom_line(linewidth = .9) +
  ggplot2::scale_colour_manual(values = ok) + ggplot2::scale_fill_manual(values = ok) +
  ggplot2::labs(x = "Age (years)", y = "Diel amplitude") + theme_bear()
p_nf <- ggplot2::ggplot(prof, ggplot2::aes(Age, nf, colour = Sex, fill = Sex)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = nf_lo, ymax = nf_hi), alpha = .15, colour = NA) +
  ggplot2::geom_line(linewidth = .9) +
  ggplot2::scale_colour_manual(values = ok) + ggplot2::scale_fill_manual(values = ok) +
  ggplot2::labs(x = "Age (years)", y = "Night fraction") + theme_bear()
p <- p_top / (p_amp | p_nf) +
  patchwork::plot_layout(heights = c(1.05, 1)) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = sprintf("Intra-specific temporal partitioning at the dump (< %.0f km)", NEAR_KM),
    subtitle = sprintf("Age modelled continuously | %d bears (%d F, %d M), %d bear-years, ages %.0f-%.0f y | %s",
                       nb, nF, nM, by, age_rng[1], age_rng[2], RS$key),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 11, colour = "grey40")))
save_fig(p, sprintf("FigV3_08_DumpPartitioning_%s", RS$key), w = 11, h = 8.5)

log_step("=== 29_dump_partitioning DONE ===")
sink(type = "message"); sink()
