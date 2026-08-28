# ==============================================================================
# 09d_lobo_fold.R  —  SINGLE FOLD grouped bear-CV (SLURM array; restart-safe)
# ==============================================================================
#
# LOGIC (same definition as brms::kfold):
#   - Bears are split deterministically into K groups (foldmap; all models use the
#     SAME fold -> paired/fair comparison). Held-out = this fold's bears.
#   - The model is REFIT on the remaining bears (update; iter=1500).
#   - On the held-out observations, the log pointwise predictive density is via
#     log_lik with allow_new_levels=TRUE, sample_new_levels="uncertainty" ->
#     predicting an UNSEEN bear (random effects drawn from the estimated sd). This
#     is the "new bear" CV.
#   - Output: held-out per-observation elpd (pointwise) -> combine sums it.
#
# ENV: MODEL=H0..H5 | RESP=intensity|duration | FOLD=1..K | K_FOLDS=10
# OUTPUT: models/lobo/<MODEL>_<RESP>_fold<FOLD>.rds  (SKIPPED if present)
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))

RESP  <- RESP_KEY()
MODEL <- Sys.getenv("MODEL", unset = "")
FOLD  <- as.integer(Sys.getenv("FOLD", unset = "0"))
K     <- as.integer(Sys.getenv("K_FOLDS", unset = "10"))
if (!nchar(MODEL) || FOLD < 1) stop("MODEL and FOLD (>=1) env are required")

init_log(sprintf("09d_lobo_%s_%s_f%d", MODEL, RESP, FOLD))
log_step("=== 09d fold | MODEL=%s RESP=%s FOLD=%d/%d ===", MODEL, RESP, FOLD, K)

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))
options(mc.cores = max(1L, slurm_cpus))
rstan::rstan_options(auto_write = TRUE)
options(brms.backend = "rstan")

lobo_dir <- file.path(mod_path, "lobo")
dir.create(lobo_dir, recursive = TRUE, showWarnings = FALSE)
out_fp <- file.path(lobo_dir, sprintf("%s_%s_fold%d.rds", MODEL, RESP, FOLD))
if (file.exists(out_fp)) { log_step("ALREADY EXISTS -> skip: %s", basename(out_fp)); quit(save = "no", status = 0) }

# ---- Model + data ------------------------------------------------------------
fp <- file.path(mod_path, sprintf("%s_%s.rds", MODEL, RESP))
if (!file.exists(fp)) stop("Model missing: ", fp)
h  <- readRDS(fp); fit <- h$fit; d <- h$data; tv <- h$time_var
log_step("Model loaded: n=%d, axis=%s, family=%s", nrow(d), tv, h$family)

# ---- Foldmap (deterministic, shared per RESP) --------------------------------
# The first task to run writes it; the others read it. tryCatch re-reads against a
# possible race. All models use the same RESP foldmap.
fm_fp <- file.path(lobo_dir, sprintf("foldmap_%s_K%d.rds", RESP, K))
make_fm <- function() {
  bears <- sort(unique(as.character(d$BearID_f)))
  set.seed(2026)
  fm <- data.frame(BearID_f = bears,
                   fold = sample(rep_len(seq_len(K), length(bears))),
                   stringsAsFactors = FALSE)
  fm
}
if (!file.exists(fm_fp)) {
  fm <- make_fm()
  tryCatch(saveRDS(fm, fm_fp), error = function(e) NULL)
} else {
  fm <- tryCatch(readRDS(fm_fp), error = function(e) make_fm())
}
heldout <- fm$BearID_f[fm$fold == FOLD]
if (!length(heldout)) stop("No held-out bear for fold ", FOLD, " (K too large?)")
test_idx <- as.character(d$BearID_f) %in% heldout
train <- d[!test_idx, , drop = FALSE]
test  <- d[ test_idx, , drop = FALSE]
log_step("Held-out bears (fold %d): %s | train n=%d, test n=%d",
         FOLD, paste(heldout, collapse = ","), nrow(train), nrow(test))

# ---- REFIT the fold model (update: formula/prior/family preserved) -----------
t0 <- Sys.time()
knots_arg <- stats::setNames(list(c(0, 24)), tv)
fit_tr <- update(fit, newdata = train,
                 iter = 1500, warmup = 1000, chains = 2, cores = max(1L, slurm_cpus),
                 knots = knots_arg, seed = 4200 + FOLD,
                 control = list(adapt_delta = 0.97, max_treedepth = 12),
                 refresh = 250, silent = 1, init = 0,
                 save_pars = brms::save_pars(all = TRUE))
log_step("fold fit time: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

# convergence summary (diagnostic)
np <- brms::nuts_params(fit_tr)
n_div <- sum(np$Value[np$Parameter == "divergent__"])
rh <- tryCatch(max(brms::rhat(fit_tr), na.rm = TRUE), error = function(e) NA_real_)
log_step("fold diag: divergent=%d, max_Rhat=%.4f", n_div, rh)

# ---- Held-out log pointwise predictive density -------------------------------
# UNSEEN bear: random effects drawn from the estimated sd (uncertainty).
ll <- brms::log_lik(fit_tr, newdata = test, allow_new_levels = TRUE,
                    sample_new_levels = "uncertainty")
logmeanexp <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
elpd_i <- apply(ll, 2, logmeanexp)          # lpd per held-out observation
log_step("fold elpd = %.1f (n_test=%d, mean/obs=%.3f)",
         sum(elpd_i), length(elpd_i), mean(elpd_i))

saveRDS(list(model = MODEL, resp = RESP, fold = FOLD, K = K,
             heldout_bears = heldout, n_test = length(elpd_i),
             elpd_pointwise = elpd_i, elpd_fold = sum(elpd_i),
             bearid = as.character(test$BearID_f),
             divergent = n_div, max_rhat = rh), out_fp)
log_step("=== 09d fold %d SAVED: %s ===", FOLD, basename(out_fp))
sink(type = "message"); try(sink(), silent = TRUE)
