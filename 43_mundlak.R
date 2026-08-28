# ==============================================================================
# 43_mundlak.R  —  Within- vs between-individual decomposition (Mundlak device)  [B6]
# ==============================================================================
# Call:  RESP=intensity Rscript 43_mundlak.R   (or RESP=duration)
#
# WHY (audit B6):
#   DumpProx/d2GarbageDump is time-varying (positions move). The Near/Far contrast
#   blends TWO different things:
#     - WITHIN  : a given bear APPROACHING the dump (plasticity)
#     - BETWEEN : some bears LIVING near the dump (habitat selection / distribution)
#   The report language ("bears near the dump shift to night") implies WITHIN, but
#   the model does not separate them. The standard partial fix = Mundlak: enter each
#   variable's INDIVIDUAL MEAN (between) and its WITHIN-DEVIATION as separate terms.
#   The difference in coefficients = direct evidence on plasticity vs distribution.
#
# OUTPUTS: models/Mundlak_<RESP>.rds
#          tables/Table_V3_26_Mundlak_<RESP>.csv
#            -> Table S54 (intensity) / S55 (duration)  (between vs within, dump & road)
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
RS <- get_response_spec()
init_log(sprintf("43_mundlak_%s", RS$key))
log_step("=== 43_mundlak.R | RESP=%s (%s) ===", RS$key, RS$family_name)

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE); options(brms.backend = "rstan")

md <- readRDS(file.path(mod_path, "model_data_50k.rds"))
d  <- md$data
decorr_covs <- md$decorr_covs
TIME_VAR <- "SolarHour_dbl"; KNOTS <- make_knots(TIME_VAR)

# ---- Mundlak decomposition: per-individual mean + within-deviation for dump & road
grp <- d$BearID_f
bm  <- function(x) ave(x, grp, FUN = function(z) mean(z, na.rm = TRUE))
d$dump_bm <- bm(d$d2GarbageDump_km_sc); d$dump_wc <- d$d2GarbageDump_km_sc - d$dump_bm
d$road_bm <- bm(d$d2Roads_km_sc);       d$road_wc <- d$d2Roads_km_sc       - d$road_bm
log_step("Mundlak: dump between SD=%.3f within SD=%.3f | road between SD=%.3f within SD=%.3f",
         sd(d$dump_bm), sd(d$dump_wc), sd(d$road_bm), sd(d$road_wc))

# drop the continuous _sc dump/road, replace with bm+wc; keep the rest as is
other_covs <- setdiff(decorr_covs, c("d2GarbageDump_km_sc", "d2Roads_km_sc"))
rhs <- paste("s(SolarHour_dbl, bs='cc', k=12)", "s(DOY, bs='tp', k=8)",
             "dump_bm", "dump_wc", "road_bm", "road_wc",
             paste(other_covs, collapse = " + "),
             "Sex_f", "Age_sc", "Season_f", "(1 | BearID_f/Year_f)", sep = " + ")
fml_str <- paste(RS$lhs, "~", rhs)
log_step("FORMULA: %s", fml_str)

bform <- if (RS$has_zi) { brms::bf(as.formula(fml_str), as.formula("zi ~ s(SolarHour_dbl, bs='cc', k=12)")) } else { brms::bf(as.formula(fml_str)) }

t0 <- Sys.time()
fit <- brms::brm(bform, data = d, family = RS$family, prior = RS$priors, knots = KNOTS,
                 iter = 2000, warmup = 1000, chains = 4, cores = 4, seed = 42,
                 silent = 1, refresh = 250, init = 0,
                 control = list(adapt_delta = 0.92, max_treedepth = 11),
                 save_pars = brms::save_pars(all = TRUE))
log_step("fit time: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

s <- summary(fit); fe <- s$fixed
pick <- c("dump_bm", "dump_wc", "road_bm", "road_wc")
out <- data.frame(
  Parameter = pick,
  Component = c("between (distribution)", "within (plasticity)",
                "between (distribution)", "within (plasticity)"),
  Estimate = round(fe[pick, "Estimate"], 4),
  Q2.5 = round(fe[pick, "l-95% CI"], 4), Q97.5 = round(fe[pick, "u-95% CI"], 4),
  Rhat = round(fe[pick, "Rhat"], 4), ESS_bulk = round(fe[pick, "Bulk_ESS"]),
  Direction_sig = (sign(fe[pick, "l-95% CI"]) == sign(fe[pick, "u-95% CI"])))
# -> Table S54 (intensity) / S55 (duration): within-vs-between (Mundlak) decomposition
utils::write.csv(out, file.path(tbl_path, sprintf("Table_V3_26_Mundlak_%s.csv", RS$key)),
                 row.names = FALSE)
cat("\n=== MUNDLAK (between vs within) ===\n"); print(out, row.names = FALSE)

save_rds_safe(list(fit = fit, response = RS$key, family = RS$family_name,
                   formula = fml_str, data = d),
              file.path(mod_path, sprintf("Mundlak_%s.rds", RS$key)))
log_step("=== 43_mundlak.R (%s) DONE ===", RS$key)
sink(type = "message"); sink()
