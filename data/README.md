# `data/` — input datasets

The pipeline reads all raw inputs from a `data/` folder whose location is set by
`PROJECT_ROOT` in `00_setup.R` (`data_path <- file.path(PROJECT_ROOT, "data")`).
This file documents the **full expected layout**, which parts ship with the
repository, and where to obtain the parts that do not (large or third-party
licensed datasets are not redistributed here).

```
data/
├── Activity/                  # per-bear accelerometer .txt   (private; see example_data/)
├── Position/                  # per-bear GPS .txt             (private; see example_data/)
├── Bear_metadata.xlsx         # bear metadata                 (private; see example_data/)
├── rasters/                   # environmental predictors
│   ├── *.tif                  # 12 predictor rasters          [SHIPPED as EnvironmentalPreds/]
│   ├── E_OBS_Climate/*.nc     # E-OBS daily climate           [DOWNLOAD — see below]
│   └── ERA5_Land/*.nc         # ERA5-Land hourly 2 m temp     [DOWNLOAD / auto — see below]
├── shapefiles/                # StudyAreaSk, MainRoads, SK_AMNP, GarbageDumpSites  [SHIPPED]
├── esri_lulc/*.tif            # ESRI 10 m land cover tiles     [DOWNLOAD — see below]
├── wdpa/{ARM,GEO}/*.shp       # WDPA protected areas           [DOWNLOAD — see below]
├── ne/*.shp                   # Natural Earth country borders  [DOWNLOAD — see below]
└── dump_symbol.png            # small legend icon              [SHIPPED]
```

## Shipped with the repository

| Repo folder | Copy to | Contents |
|---|---|---|
| `EnvironmentalPreds/*.tif` | `data/rasters/` | 8 distance rasters (`d2Builtareas, d2Crops, d2Forest, d2GarbageDump, d2ProtectedAreas, d2Rangeland, d2Roads, d2Water`) + 4 topographic (`Elevation, Slope, Aspect, Roughness`) |
| `shapefiles/` | `data/shapefiles/` | `StudyAreaSk`, `MainRoads`, `SK_AMNP`, `GarbageDumpSites` |
| `example_data/` | see `example_data/README.md` | synthetic `BTR_example` activity/position/metadata |

The 12 `.tif` predictors are the **direct inputs** to the modelling pipeline
(`03_covariate_extraction.R`) and to `plot_predictors.py`. They were derived from
the source datasets below: the `d2*` rasters are Euclidean distances to the
corresponding source vectors / land-cover classes, and the topographic rasters
come from a DEM.

## Download yourself (not redistributed here)

These are large and/or carry third-party licences. Obtain them from the original
providers and place them at the paths shown above. Versions listed are those
used in this study (encoded in the original filenames).

| Dataset | Used by | Version / files | Source |
|---|---|---|---|
| **ERA5-Land** hourly 2 m temperature (primary thermal covariate) | `03b_era5_hourly.R` | `ERA5_Land/era5_land_t2m_<YEAR>.nc`, variable `2m_temperature`, ~9 km hourly | Copernicus Climate Data Store, dataset `reanalysis-era5-land`: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=download — **`03b` downloads this automatically via the Python `cdsapi`**; requires a `~/.cdsapirc` key and acceptance of the ERA5-Land licence. |
| **E-OBS** daily gridded climate (old/sensitivity precipitation + daily temperature) | `03_covariate_extraction.R` | `E_OBS_Climate/{rr,tg,tn,tx}_ens_mean_0.1deg_reg_v31.0e.nc` (precip, mean, min, max), 0.1° | ECA&D E-OBS ensemble gridded dataset: https://www.ecad.eu (also on the Copernicus CDS). |
| **ESRI 10 m Annual Land Cover** (figure background + upstream land-cover distance rasters) | `plot_studyarea.py` | `esri_lulc/{37S,37T,38S,38T}_20220101-20230101.tif` (2022 epoch, UTM 37/38 S/T) | Esri / Impact Observatory 10 m Annual LULC, via Esri Living Atlas: https://livingatlas.arcgis.com |
| **WDPA** protected areas | `plot_studyarea.py` (and upstream `d2ProtectedAreas` together with `SK_AMNP`) | `wdpa/ARM/...`, `wdpa/GEO/...` `WDPA_WDOECM_Aug2026_Public_<ISO>_shp-polygons.shp` (Aug 2026 release; Armenia + Georgia) | World Database on Protected Areas (UNEP-WCMC & IUCN): https://www.protectedplanet.net |
| **Natural Earth** country borders | `plot_studyarea.py` | `ne/ne_10m_admin_0_boundary_lines_land.shp`, 1:10 m | Natural Earth (public domain): https://www.naturalearthdata.com |

## Which script reads what

- `03_covariate_extraction.R` — `rasters/*.tif` (all 12) + `E_OBS_Climate/*.nc`.
- `03b_era5_hourly.R` — `ERA5_Land/*.nc` (self-downloads if absent).
- `corr_era5.R` — covariate correlations incl. the ERA5 hourly temperature (Table S6).
- `plot_predictors.py` — `rasters/*.tif` (predictor panel figure).
- `plot_studyarea.py` — `shapefiles/*`, `esri_lulc/*.tif`, `wdpa/{ARM,GEO}/*.shp`, `ne/*.shp`, `dump_symbol.png`.
- `30_conflict_risk_map.R`, `30b_fig09_replot.R` — `shapefiles/StudyAreaSk`, `GarbageDumpSites`, `MainRoads`.

> Scripts from `03` onward therefore require these external layers. Stages `01`–`02`
> run from `example_data/` alone (no rasters, no cluster). Please respect each
> provider's licence and attribution terms when redistributing any derived product.
