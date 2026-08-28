# ==============================================================================
# corr_era5.R  —  Environmental covariate correlation matrix (ERA5 hourly temp)
# ==============================================================================
# Kendall correlation among all environmental covariates, using the hourly ERA5
# temperature (temp_hourly) rather than the daily E-OBS mean.
#
# Produces:
#   Table_04b_CorrelationMatrix_ERA5.csv  -> Table S6 (SI)
#   Fig_00b_CorrelationMatrix_ERA5.{png,pdf} -> Fig S3b (SI)
# ==============================================================================
suppressWarnings(suppressMessages({library(corrplot); library(dplyr); library(tidyr)}))
B <- file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct")
d <- readRDS(file.path(B, "data", "hourly_full.rds"))
cat("hourly_full rows=", nrow(d), " cols=", ncol(d), "\n")
cat("name matches:", paste(grep("temp|precip|d2|Elev|Slope|Aspect|Rough", names(d), value = TRUE, ignore.case = TRUE), collapse = ", "), "\n")
topo  <- c("Elevation", "Slope", "Aspect", "Roughness")
distk <- paste0(c("d2Builtareas", "d2Crops", "d2Forest", "d2GarbageDump", "d2ProtectedAreas", "d2Rangeland", "d2Roads", "d2Water"), "_km")
want  <- c(topo, distk, "precipitation", "temp_hourly")
have  <- intersect(want, names(d)); miss <- setdiff(want, names(d))
cat("SELECTED:", paste(have, collapse = ", "), "\n")
if (length(miss)) cat("MISSING :", paste(miss, collapse = ", "), "\n")
cd <- d %>% dplyr::select(dplyr::all_of(have)) %>% tidyr::drop_na()
cat("complete rows=", nrow(cd), "\n")
set.seed(42); n <- nrow(cd); ss <- if (n > 15000) cd[sample(n, 15000), ] else cd
cm <- cor(ss, method = "kendall", use = "pairwise.complete.obs")
# -> Table S6 (SI): environmental covariate correlation matrix (Kendall tau)
write.csv(as.data.frame(cm), file.path(B, "tables", "Table_04b_CorrelationMatrix_ERA5.csv"))
ttl <- "Kendall correlation - Environmental covariates"
# -> Fig S3b (SI)
mk <- function(dev) { fp <- file.path(B, "figures", paste0("Fig_00b_CorrelationMatrix_ERA5.", dev))
  if (dev == "png") png(fp, width = 10, height = 10, units = "in", res = 300) else pdf(fp, width = 10, height = 10)
  corrplot::corrplot(cm, method = "color", type = "lower", tl.col = "black", tl.srt = 45, tl.cex = 0.8,
    addCoef.col = "black", number.cex = 0.6,
    col = colorRampPalette(c("#332288", "white", "#CC6677"))(100), title = ttl, mar = c(0, 0, 2, 0)); dev.off() }
mk("png"); mk("pdf")
cat("\n--- temp_hourly Kendall tau vs others ---\n"); print(round(sort(cm[, "temp_hourly"]), 3))
cat("DONE n_subsample=", nrow(ss), "\n")
