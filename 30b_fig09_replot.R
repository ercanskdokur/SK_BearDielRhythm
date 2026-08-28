# ==============================================================================
# 30b_fig09_replot.R  —  FigV3_09 (conflict-risk map) replot ONLY
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
suppressMessages({ library(terra); library(sf) })
RS <- get_response_spec()
init_log(sprintf("30b_fig09_replot_%s", RS$key))
log_step("=== 30b_fig09_replot | RESP=%s ===", RS$key)

rp <- file.path(mod_path, sprintf("conflict_risk_rasters_%s.tif", RS$key))
if (!file.exists(rp)) stop("Missing: ", rp)
r <- terra::rast(rp)
if (!"risk" %in% names(r)) stop("no risk layer; layers: ", paste(names(r), collapse=", "))
pdat <- as.data.frame(r[["risk"]], xy = TRUE, na.rm = TRUE)
names(pdat) <- c("X_utm", "Y_utm", "risk")
log_step("Raster cells: %d, risk range [%.4f, %.4f]",
         nrow(pdat), min(pdat$risk), max(pdat$risk))

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
  theme_bear() +
  ggplot2::theme(plot.title    = ggplot2::element_text(size = 12),
                 plot.subtitle = ggplot2::element_text(size = 9.5, colour = "grey40"))
save_fig(p, sprintf("FigV3_09_ConflictRiskMap_%s", RS$key), w = 9, h = 7)
log_step("=== 30b_fig09_replot DONE: FigV3_09_ConflictRiskMap_%s ===", RS$key)
sink(type = "message"); sink()
