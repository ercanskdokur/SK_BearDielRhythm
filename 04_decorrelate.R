# ==============================================================================
# 04_decorrelate.R  —  PCA-based decorrelation of collinear covariates
# ==============================================================================
# GOAL: resolve catastrophic multicollinearity (VIF > 10, cond.no = 4003):
#   - 5 land-cover distance variables (Built/Crops/Forest/Rangeland/Water, r=0.81-0.99)
#     -> PCA -> 2 orthogonal PCs (LandPC1, LandPC2)
#   - temp_mean & temp_max (r=0.85) -> keep only temp_mean
#   - Slope & Roughness (r=0.83)    -> keep only Slope
#
# Produces (SI): Table S7 (Table_V2_01_PCA_Loadings   — PCA loadings),
#                Table S8 (Table_V2_02_VIF_After       — VIF, canonical set),
#                Table S9 (Table_V2_02b_FocalCollinearity — focal collinearity),
#                Fig S3a  (FigV2_01_PCA_Biplot          — PCA biplot).
# Also writes the modelling-ready decorrelated data frames and the PCA object
# (pca_landcover.rds) for transforming new data.
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
init_log("04_decorrelate")
log_step("=== 04_decorrelate.R: START ===")

# ---- 1. LOAD DATA ------------------------------------------------------------
hourly_full <- readRDS(file.path(dat_path, "hourly_full.rds"))
daily_full  <- readRDS(file.path(dat_path, "daily_full.rds"))
bout_full   <- readRDS(file.path(dat_path, "bout_full.rds"))

log_step("hourly: n=%d, daily: n=%d, bout: n=%d",
         nrow(hourly_full), nrow(daily_full), nrow(bout_full))

# ---- 2. COVARIATE SELECTION --------------------------------------------------
# Decorrelated set:
#   - LandPC1, LandPC2          (5 -> 2; PCA on Built/Crops/Forest/Rangeland/Water)
#   - d2Roads_km                (VIF ~1.98, OK)
#   - d2GarbageDump_km          (VIF ~27, but conceptually critical; keep)
#   - d2ProtectedAreas_km       (VIF ~17.6; re-checked after PCA removes its co-vars)
#   - temp_mean                 (drop temp_max; keep precipitation)
#   - precipitation             (VIF 1.16, OK)
#   - Elevation                 (VIF ~2)
#   - Slope                     (drop Roughness)
#   - Sex_f, Age_sc             (demographic — needed for H2)
#
# Land-cover PCA scope: 5 distance variables
land_cov_vars <- c("d2Builtareas_km", "d2Crops_km", "d2Forest_km",
                   "d2Rangeland_km", "d2Water_km")

# ---- 3. TRAIN PCA (on hourly_full — the largest data) -----------------------
# Drop rows with NA, then z-score, then prcomp
log_step("Training PCA: %d rows, 5 variables", nrow(hourly_full))
pca_data <- hourly_full %>%
  dplyr::select(dplyr::all_of(land_cov_vars)) %>%
  tidyr::drop_na()
log_step("  Complete cases: %d (%.1f%%)",
         nrow(pca_data), 100 * nrow(pca_data) / nrow(hourly_full))

pca_landcover <- stats::prcomp(pca_data, center = TRUE, scale. = TRUE)

# Variance explained
ve <- (pca_landcover$sdev^2) / sum(pca_landcover$sdev^2)
log_step("  Variance explained: PC1=%.1f%% PC2=%.1f%% PC3=%.1f%% PC4=%.1f%% PC5=%.1f%%",
         100*ve[1], 100*ve[2], 100*ve[3], 100*ve[4], 100*ve[5])

# Loading table
loadings <- as.data.frame(pca_landcover$rotation)
loadings$Variable <- rownames(loadings)
loadings <- loadings[, c("Variable", paste0("PC", 1:5))]
loadings_long <- loadings %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("PC"), ~ round(., 3)))
# -> Table S7 (SI): PCA loadings of the land-cover distance variables
write.csv(loadings_long,
          file.path(tbl_path, "Table_V2_01_PCA_Loadings.csv"),
          row.names = FALSE)
log_step("  PCA loadings written")
cat("\nPCA Loadings:\n")
print(loadings_long)

# Determine how many PCs to keep — Kaiser criterion (eigenvalue > 1) or 80% variance
cum_ve <- cumsum(ve)
n_keep_kaiser <- sum(pca_landcover$sdev^2 > 1)
n_keep_80pct  <- which(cum_ve >= 0.80)[1]
n_keep <- max(2, min(n_keep_kaiser, n_keep_80pct))
log_step("  Kaiser criterion: %d PC | 80%% variance: %d PC | chosen: %d",
         n_keep_kaiser, n_keep_80pct, n_keep)
# Force a minimum of 2 PCs (we need at least 2 to capture both the
# Built-Crops-Range axis AND a Forest-Water axis)
if (n_keep < 2) n_keep <- 2

# Save the PCA object so new data can be transformed later
save_rds_safe(list(pca = pca_landcover, n_keep = n_keep, vars = land_cov_vars),
              file.path(mod_path, "pca_landcover.rds"))

# ---- 4. ADD PCA SCORES TO ALL TABLES ----------------------------------------
apply_pca <- function(df, pca_obj, vars, n_keep) {
  # Which rows can enter the PCA
  ok <- stats::complete.cases(df[, vars, drop = FALSE])
  scores <- matrix(NA_real_, nrow = nrow(df), ncol = n_keep,
                   dimnames = list(NULL, paste0("LandPC", 1:n_keep)))
  if (any(ok)) {
    # predict() expects raw data; prcomp transforms internally
    new_scores <- stats::predict(pca_obj, newdata = df[ok, vars, drop = FALSE])
    scores[ok, ] <- new_scores[, 1:n_keep, drop = FALSE]
  }
  for (k in seq_len(n_keep)) {
    df[[paste0("LandPC", k, "_sc")]] <- as.numeric(scale(scores[, k]))
  }
  df
}

hourly_full <- apply_pca(hourly_full, pca_landcover, land_cov_vars, n_keep)
daily_full  <- apply_pca(daily_full,  pca_landcover, land_cov_vars, n_keep)
bout_full   <- apply_pca(bout_full,   pca_landcover, land_cov_vars, n_keep)

# ---- 5. ENSURE ADDITIONAL Z-SCORES ------------------------------------------
# The _sc versions mostly exist already; here we just add any that are missing.
extra_vars <- c("d2Roads_km", "d2GarbageDump_km", "d2ProtectedAreas_km",
                "temp_mean", "precipitation", "Elevation", "Slope")

ensure_z <- function(df) {
  for (v in extra_vars) {
    z <- paste0(v, "_sc")
    if (!z %in% names(df) && v %in% names(df) && is.numeric(df[[v]])) {
      df[[z]] <- as.numeric(scale(df[[v]]))
    }
  }
  # Age_sc
  if ("Current_Age" %in% names(df) && !"Age_sc" %in% names(df)) {
    df$Age_sc <- as.numeric(scale(df$Current_Age))
  }
  if ("Sex" %in% names(df) && !"Sex_f" %in% names(df)) {
    df$Sex_f <- factor(df$Sex, levels = c("F", "M"))
  }
  if (!"BearID_f" %in% names(df) && "BearID" %in% names(df)) {
    df$BearID_f <- factor(df$BearID)
  }
  if (!"Year_f" %in% names(df) && "Year" %in% names(df)) {
    df$Year_f <- factor(df$Year)
  }
  if (!"Season_f" %in% names(df) && "Season" %in% names(df)) {
    df$Season_f <- factor(df$Season, levels = c("Mating", "Hyperphagia"))
  }
  df
}

hourly_full <- ensure_z(hourly_full)
daily_full  <- ensure_z(daily_full)
bout_full   <- ensure_z(bout_full)

# ---- 6. RECOMPUTE VIF (decorrelated set) ------------------------------------
decorr_covs <- c(paste0("LandPC", seq_len(n_keep), "_sc"),
                 "d2Roads_km_sc", "d2GarbageDump_km_sc", "d2ProtectedAreas_km_sc",
                 "temp_mean_sc", "precipitation_sc",
                 "Elevation_sc", "Slope_sc")

vif_check <- function(df, covs, label) {
  X <- as.matrix(df[, intersect(covs, names(df)), drop = FALSE])

  # --- ORDER MATTERS ---------------------------------------------------------
  # Clean rows FIRST (complete.cases), THEN check for constant columns. Doing it
  # the other way (dropping columns that contain any NA first) drops every column
  # in hourly_full, producing an EMPTY VIF table.
  X <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(X) < 10) {
    log_step("  VIF[%s]: n=%d after complete.cases — skipping", label, nrow(X))
    return(NULL)
  }
  keep <- apply(X, 2, function(cc) stats::sd(cc) > 0)
  if (any(!keep)) {
    log_step("  VIF[%s]: constant column(s) dropped: %s",
             label, paste(colnames(X)[!keep], collapse = ", "))
  }
  X <- X[, keep, drop = FALSE]
  if (ncol(X) < 2) return(NULL)
  log_step("  VIF[%s]: n=%d rows x %d covariates", label, nrow(X), ncol(X))
  # ---------------------------------------------------------------------------

  R <- stats::cor(X)
  vifs <- diag(solve(R))
  kappa_val <- kappa(R)
  data.frame(
    Dataset = label,
    Covariate = names(vifs),
    VIF = round(vifs, 2),
    Flag = ifelse(vifs > 10, "PROBLEMATIC",
                  ifelse(vifs > 5, "CAUTION", "OK")),
    kappa = round(kappa_val, 1),
    N_rows = nrow(X)
  )
}

vif_h <- vif_check(hourly_full, decorr_covs, "hourly")
vif_d <- vif_check(daily_full,  decorr_covs, "daily")
vif_b <- vif_check(bout_full,   decorr_covs, "bout")
vif_all <- dplyr::bind_rows(vif_h, vif_d, vif_b)

if (is.null(vif_all) || nrow(vif_all) == 0) {
  stop("VIF table came out EMPTY — the ordering bug may have returned; stopping.")
}
# -> Table S8 (SI): variance-inflation factors (canonical set)
write.csv(vif_all, file.path(tbl_path, "Table_V2_02_VIF_After.csv"), row.names = FALSE)
cat("\n=== VIF (decorrelated set) ===\n")
print(vif_all, row.names = FALSE)

# ---- 6b. FOCAL-VARIABLE COLLINEARITY CHECK ----------------------------------
# The focal variable d2GarbageDump_km has VIF ~27 but is conceptually critical,
# so it is kept and its collinearity is reported explicitly. The PCA loadings
# (PC1: 0.452/0.446/0.453/0.453/0.432 — all ~equal positive) show LandPC1 is a
# pure "remoteness from everything" axis; because dumps sit near settlements,
# d2GarbageDump is expected to load onto LandPC1.
focal_report <- function(df, label) {
  need <- c("d2GarbageDump_km_sc", "LandPC1_sc", "LandPC2_sc",
            "d2Roads_km_sc", "d2ProtectedAreas_km_sc")
  need <- intersect(need, names(df))
  Z <- df[, need, drop = FALSE]
  Z <- Z[stats::complete.cases(Z), , drop = FALSE]
  if (nrow(Z) < 10 || !"d2GarbageDump_km_sc" %in% names(Z)) return(NULL)
  r2 <- summary(stats::lm(d2GarbageDump_km_sc ~ ., data = as.data.frame(Z)))$r.squared
  data.frame(Dataset = label,
             Focal = "d2GarbageDump_km_sc",
             R2_on_others = round(r2, 4),
             VIF_implied = round(1 / (1 - r2), 2),
             Verdict = ifelse(1/(1-r2) > 5,
                              "CAUTION: consider residualizing against LandPC1", "OK"))
}
focal_all <- dplyr::bind_rows(focal_report(hourly_full, "hourly"),
                              focal_report(daily_full,  "daily"))
if (!is.null(focal_all) && nrow(focal_all) > 0) {
  # -> Table S9 (SI): focal collinearity of the dump-distance variable
  write.csv(focal_all, file.path(tbl_path, "Table_V2_02b_FocalCollinearity.csv"),
            row.names = FALSE)
  cat("\n=== FOCAL VARIABLE (d2GarbageDump) COLLINEARITY ===\n")
  print(focal_all, row.names = FALSE)
}

# ---- 7. PCA BIPLOT FIGURE ---------------------------------------------------
pca_scores_df <- data.frame(
  PC1 = pca_landcover$x[, 1],
  PC2 = pca_landcover$x[, 2]
)
loadings_df <- as.data.frame(pca_landcover$rotation[, 1:2])
loadings_df$Variable <- gsub("_km$", "", rownames(loadings_df))
# Arrow scaling
arrow_scale <- max(abs(c(pca_scores_df$PC1, pca_scores_df$PC2))) /
               max(abs(c(loadings_df$PC1, loadings_df$PC2))) * 0.7

p_biplot <- ggplot2::ggplot(pca_scores_df, ggplot2::aes(x = PC1, y = PC2)) +
  ggplot2::geom_point(alpha = 0.05, size = 0.4, color = "grey50") +
  ggplot2::geom_segment(data = loadings_df,
                        ggplot2::aes(x = 0, y = 0,
                                     xend = PC1 * arrow_scale,
                                     yend = PC2 * arrow_scale),
                        arrow = grid::arrow(length = ggplot2::unit(0.3, "cm")),
                        color = pal$hyperphagia, linewidth = 1.0) +
  ggplot2::geom_text(data = loadings_df,
                     ggplot2::aes(x = PC1 * arrow_scale * 1.1,
                                  y = PC2 * arrow_scale * 1.1,
                                  label = Variable),
                     color = pal$hyperphagia, fontface = "bold",
                     size = 4) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::labs(x = sprintf("PC1 (%.1f%% var)", 100 * ve[1]),
                y = sprintf("PC2 (%.1f%% var)", 100 * ve[2]),
                title = "PCA biplot – 5 land-cover distance variables",
                subtitle = "Used to construct LandPC1, LandPC2") +
  theme_bear()
# -> Fig S3a (SI): PCA biplot of the land-cover distance variables
save_fig(p_biplot, "FigV2_01_PCA_Biplot", w = 9, h = 7)

# ---- 8. SAVE DECORRELATED DATA FRAMES ---------------------------------------
save_rds_safe(hourly_full, file.path(dat_path, "hourly_decorr.rds"))
save_rds_safe(daily_full,  file.path(dat_path, "daily_decorr.rds"))
save_rds_safe(bout_full,   file.path(dat_path, "bout_decorr.rds"))

# Also save the decorrelated covariate list (for downstream scripts)
save_rds_safe(list(n_keep = n_keep,
                   decorr_covs = decorr_covs,
                   land_cov_vars = land_cov_vars),
              file.path(mod_path, "decorr_meta.rds"))

cat("\n========================================\n")
cat("  DECORRELATE SUMMARY\n")
cat("========================================\n")
cat(sprintf("PCs retained        : %d\n", n_keep))
cat(sprintf("Decorrelated covars : %d\n", length(decorr_covs)))
cat(sprintf("Max VIF (hourly)    : %.2f\n",
            ifelse(is.null(vif_h), NA, max(vif_h$VIF, na.rm = TRUE))))
cat(sprintf("Kappa (hourly)      : %.1f\n",
            ifelse(is.null(vif_h), NA, vif_h$kappa[1])))
cat("Before: Max VIF = 776, Kappa = 4003 (ILL-CONDITIONED)\n")

log_step("=== 04_decorrelate.R: DONE ===")
sink(type = "message"); sink()
