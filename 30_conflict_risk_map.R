# ==============================================================================
# 30_conflict_risk_map.R  —  Spatial human-bear conflict-risk map
# ==============================================================================
# LOGIC — risk = EXPOSURE x OVERLAP:
#   Conflict risk in a cell is the product of two things:
#     (a) the bear's tendency to be in that cell            -> "exposure"
#     (b) when active there, coinciding with human hours     -> "temporal overlap"
#   H3 gives (b): from the cell's dump/road distance it predicts the diel curve and
#   we compute the SHARE of the curve falling in the conflict window (timing-
#   normalised; cleaned of the overall activity level -> no "level vs timing" mix).
#   For (a) we use a smoothed GPS use-intensity surface (focal mean, not a true KDE).
#
# IMPORTANT LIMIT: this is not a PROBABILITY map, but a "given H3, if a bear is
#   active here, what share of its activity falls in the human window" map. Real
#   conflict also depends on human density; we do NOT multiply by a separate human
#   layer because we have no human-use data. 
#
# WINDOWS (solar time, SAME definition as 15_conflict_windows.R):
#   morning = SR-1 .. SR+4 ; evening = SS-4 .. SS+1
#
# OUTPUTS:
#   tables/Table_V3_21_ConflictHotspots_<RESP>.csv
#     -> Table S44 (intensity) / S45 (duration)
#   figures/FigV3_09_ConflictRiskMap_<RESP>.{png,pdf}
#   models/conflict_risk_rasters_<RESP>.tif  (3 layers: morning, evening, combined)
# ==============================================================================

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))


suppressMessages({ library(rstan); library(brms) })

RS <- get_response_spec()
init_log(sprintf("30_conflict_risk_map_%s", RS$key))
log_step("=== 30_conflict_risk_map | RESP=%s ===", RS$key)

slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = max(1L, slurm_cpus))
GRID_RES <- as.numeric(Sys.getenv("GRID_RES", "500"))   # metres
N_DRAWS  <- as.integer(Sys.getenv("N_DRAWS", "300"))

# ---- MODEL -------------------------------------------------------------------
h3 <- load_model("H3")
fit <- h3$fit
d   <- h3$data
TIME_VAR <- h3$time_var
log_step("H3 loaded: n=%d, axis=%s", nrow(d), TIME_VAR)

# ---- 1. PREDICTION GRID ------------------------------------------------------
# GRID_RES-metre cells within the study-area polygon; a covariate per cell.
sa <- sf::st_read(file.path(shp_path, "StudyAreaSk.shp"), quiet = TRUE)
sa <- sf::st_transform(sa, study_crs)
tmpl <- terra::rast(terra::vect(sa), resolution = GRID_RES, crs = study_crs)
grid <- terra::as.points(terra::init(tmpl, 1))
grid <- terra::mask(terra::init(tmpl, 1), terra::vect(sa))
gp   <- terra::as.points(grid)
log_step("Grid: %d cells (%.0f m)", nrow(gp), GRID_RES)

gdf <- as.data.frame(terra::geom(gp)[, c("x", "y")])
names(gdf) <- c("X_utm", "Y_utm")

# Extract covariates from the rasters
for (nm in names(raster_files)) {
  fp <- file.path(gis_path, raster_files[nm])
  if (!file.exists(fp)) next
  r <- terra::rast(fp)
  # NOTE: as.matrix(gdf) was taken from a GROWING gdf inside the loop (from round 2
  #   on, >2 columns -> "[vect] not an appropriate matrix"). Reuse the fixed grid
  #   points gp (study_crs); they are row-aligned with gdf.
  pv <- terra::project(gp, terra::crs(r))
  gdf[[nm]] <- terra::extract(r, pv)[[2]]
}
for (v in dist_vars) if (v %in% names(gdf)) gdf[[paste0(v, "_km")]] <- gdf[[v]] / 1000
gdf <- gdf[stats::complete.cases(gdf[, intersect(dist_vars_km, names(gdf))]), ]
log_step("Cells with covariates: %d", nrow(gdf))

# ---- 2. MAP CELLS TO H3'S TERTILES -------------------------------------------
# CRITICAL: the tertile boundaries were cut from ALL data in 05_prep_modeldata.R
# and stored in model_data_50k.rds. If we don't use the same boundaries, the map
# and the model use different "Near" definitions.
md  <- readRDS(file.path(mod_path, "model_data_50k.rds"))
brk <- md$breaks
cut_ref <- function(x, b, labs) {
  b[1] <- -Inf; b[length(b)] <- Inf
  cut(x, breaks = unique(b), labels = labs, include.lowest = TRUE)
}
gdf$DumpProx_f <- cut_ref(gdf$d2GarbageDump_km, brk$dump, c("Near","Mid","Far"))
gdf$RoadProx_f <- cut_ref(gdf$d2Roads_km,       brk$road, c("Near","Mid","Far"))
gdf$DumpProx_o <- ordered(gdf$DumpProx_f, levels = c("Near","Mid","Far"))
gdf$RoadProx_o <- ordered(gdf$RoadProx_f, levels = c("Near","Mid","Far"))
log_step("Cell tertile distribution — dump: %s",
         paste(names(table(gdf$DumpProx_f)), as.integer(table(gdf$DumpProx_f)),
               sep = "=", collapse = " | "))

# LandPC scores: same transform as the PCA trained in 04
pcaobj <- readRDS(file.path(mod_path, "pca_landcover.rds"))
ok <- stats::complete.cases(gdf[, pcaobj$vars, drop = FALSE])
sc <- matrix(NA_real_, nrow(gdf), pcaobj$n_keep)
sc[ok, ] <- stats::predict(pcaobj$pca, newdata = gdf[ok, pcaobj$vars, drop = FALSE])[, 1:pcaobj$n_keep]
for (k in seq_len(pcaobj$n_keep)) gdf[[paste0("LandPC", k, "_sc")]] <- as.numeric(scale(sc[, k]))


for (v in md$decorr_covs) {
  if (!v %in% names(gdf)) gdf[[v]] <- stats::median(d[[v]], na.rm = TRUE)
}
gdf$TempTert_f <- factor("Mild", levels = c("Cool","Mild","Warm"))
gdf$TempTert_o <- ordered("Mild", levels = c("Cool","Mild","Warm"))
gdf$Sex_f    <- levels(d$Sex_f)[1]
gdf$Age_sc   <- stats::median(d$Age_sc, na.rm = TRUE)
gdf$Season_f <- levels(d$Season_f)[which.max(table(d$Season_f))]
gdf$DOY      <- stats::median(d$DOY, na.rm = TRUE)
if (!RS$has_zi) gdf$N_readings <- stats::median(d$N_readings, na.rm = TRUE)

# ---- 3. CONFLICT WINDOWS (on the solar axis) ---------------------------------
# SAME definition as 15_conflict_windows.R: morning = SR-1..SR+4, evening = SS-4..SS+1.
# Because the solar axis is double-anchored (SR->6, SS->18) the windows are FIXED:
#   morning = 5..10 ; evening = 14..19
if (TIME_VAR == "SolarHour_dbl") {
  win_morning <- c(5, 10); win_evening <- c(14, 19)
  log_step("Windows (double-anchored solar): morning=%s, evening=%s",
           paste(win_morning, collapse="-"), paste(win_evening, collapse="-"))
} else {
  stop("30 was written only for the SolarHour_dbl axis (TIME_AXIS=solar_double). Got: ",
       TIME_VAR)
}

hours <- seq(0, 24, length.out = 49)[-49]
dt_h  <- 24 / length(hours)
in_win <- function(h, w) h >= w[1] & h < w[2]

# ---- 4. DIEL CURVE PER CELL -> WINDOW SHARE ----------------------------------
# Chunk the cells for memory.
CH <- 2000
chunks <- split(seq_len(nrow(gdf)), ceiling(seq_len(nrow(gdf)) / CH))
res <- vector("list", length(chunks))
t0 <- Sys.time()
for (ci in seq_along(chunks)) {
  ix <- chunks[[ci]]
  cur <- matrix(NA_real_, length(ix), length(hours))
  for (hi in seq_along(hours)) {
    nd <- gdf[ix, , drop = FALSE]
    nd[[TIME_VAR]] <- hours[hi]
    ep <- brms::posterior_epred(fit, newdata = nd, re_formula = NA, ndraws = N_DRAWS)
    cur[, hi] <- colMeans(ep)
  }
  auc <- rowSums(cur) * dt_h
  res[[ci]] <- data.frame(
    idx = ix,
    frac_morning = rowSums(cur[, in_win(hours, win_morning), drop = FALSE]) * dt_h / auc,
    frac_evening = rowSums(cur[, in_win(hours, win_evening), drop = FALSE]) * dt_h / auc,
    auc = auc
  )
  log_step("  chunk %d/%d (%.1f min)", ci, length(chunks),
           as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
rr <- do.call(rbind, res)
gdf$frac_morning <- NA_real_; gdf$frac_evening <- NA_real_; gdf$auc <- NA_real_
gdf$frac_morning[rr$idx] <- rr$frac_morning
gdf$frac_evening[rr$idx] <- rr$frac_evening
gdf$auc[rr$idx]          <- rr$auc
gdf$frac_conflict <- gdf$frac_morning + gdf$frac_evening

# ---- 5. EXPOSURE LAYER (smoothed GPS use-intensity) --------------------------
# Temporal overlap alone is not risk: a high overlap in a cell the bear NEVER goes
# to is meaningless. A simple smoothed use-intensity surface from GPS positions
# (focal-mean smoothing of the fix count; not a true kernel density estimate).
pos <- readRDS(file.path(dat_path, "position_interpolated.rds"))
pv  <- terra::vect(as.matrix(pos[, c("X_utm","Y_utm")]), type = "points", crs = study_crs)
use <- terra::rasterize(pv, tmpl, fun = "length", background = 0)  # "count" -> masked by dplyr::count; length is the base fn
use <- terra::focal(use, w = 5, fun = "mean", na.policy = "omit")   # smooth
use_v <- terra::extract(use, terra::vect(as.matrix(gdf[, c("X_utm","Y_utm")]),
                                         type = "points", crs = study_crs))[[2]]
gdf$use_rel <- use_v / max(use_v, na.rm = TRUE)
gdf$risk <- gdf$frac_conflict * gdf$use_rel     # overlap x exposure
log_step("Risk: %s", paste(round(summary(gdf$risk), 4), collapse = " "))

# ---- 6. RASTER + HOTSPOT TABLE -----------------------------------------------
mk <- function(col) {
  r <- terra::rast(tmpl)
  cells <- terra::cellFromXY(r, as.matrix(gdf[, c("X_utm","Y_utm")]))
  r[cells] <- gdf[[col]]
  names(r) <- col
  r
}
stk <- c(mk("frac_morning"), mk("frac_evening"), mk("risk"))
terra::writeRaster(stk, file.path(mod_path, sprintf("conflict_risk_rasters_%s.tif", RS$key)),
                   overwrite = TRUE)

# Sort over UNIQUE CELLS (X,Y): a repeated coordinate was inflating the "top-10".
gdf_u <- gdf[!duplicated(gdf[, c("X_utm","Y_utm")]), ]
hot <- gdf_u[order(-gdf_u$risk), ][seq_len(min(50, nrow(gdf_u))), ]
hot_out <- data.frame(
  rank = seq_len(nrow(hot)),
  X_utm = round(hot$X_utm), Y_utm = round(hot$Y_utm),
  d2GarbageDump_km = round(hot$d2GarbageDump_km, 3),
  d2Roads_km = round(hot$d2Roads_km, 3),
  frac_morning = round(hot$frac_morning, 4),
  frac_evening = round(hot$frac_evening, 4),
  use_rel = round(hot$use_rel, 4),
  risk = round(hot$risk, 5),
  response = RS$key
)
# -> Table S44 (intensity) / S45 (duration): spatial conflict hotspots (top cells)
write.csv(hot_out,
          file.path(tbl_path, sprintf("Table_V3_21_ConflictHotspots_%s.csv", RS$key)),
          row.names = FALSE)
cat("\n=== TOP 10 RISKIEST CELLS ===\n"); print(utils::head(hot_out, 10), row.names = FALSE)

# ---- 7. FIGURE ---------------------------------------------------------------
pdat <- gdf[!is.na(gdf$risk), ]
dump_shp <- tryCatch(sf::st_transform(sf::st_read(
  file.path(shp_path, "GarbageDumpSites.shp"), quiet = TRUE), study_crs), error = function(e) NULL)
road_shp <- tryCatch(sf::st_transform(sf::st_read(
  file.path(shp_path, "MainRoads.shp"), quiet = TRUE), study_crs), error = function(e) NULL)

p <- ggplot2::ggplot() +
  ggplot2::geom_raster(data = pdat, ggplot2::aes(X_utm, Y_utm, fill = risk)) +
  ggplot2::scale_fill_viridis_c(option = "inferno", name = "Conflict\nrisk") +
  { if (!is.null(road_shp)) ggplot2::geom_sf(data = road_shp, colour = "grey85",
                                             linewidth = .4, alpha = .8) } +
  { if (!is.null(dump_shp)) ggplot2::geom_sf(data = dump_shp, fill = NA,
                                             colour = "cyan", linewidth = .7) } +
  ggplot2::coord_sf(datum = sf::st_crs(study_crs)) +
  ggplot2::labs(
    title = "Human–bear conflict risk: temporal overlap x space use",
    caption = "Cyan = garbage dump; grey = main roads. NOT a probability; human density not included.",
    x = NULL, y = NULL) +
  theme_bear()
save_fig(p, sprintf("FigV3_09_ConflictRiskMap_%s", RS$key), w = 9, h = 7)

log_step("=== 30_conflict_risk_map DONE ===")
sink(type = "message"); sink()
