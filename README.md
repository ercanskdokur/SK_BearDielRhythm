# SK_BearDielRhythm

Analysis code for a study of **brown bear (*Ursus arctos*) diel activity in a
human-dominated landscape of northeastern Türkiye** (Sarıkamış region). The
pipeline fits Bayesian hierarchical generalized additive mixed models (GAMMs,
`brms`/`rstan`) to accelerometer- and GPS-derived activity, and asks whether
bears respond to human infrastructure (a garbage dump, roads) by **re-timing and
re-shaping** their activity rather than by reducing it.

> This repository accompanies the manuscript:
> **[Authors]. [Title]. [Journal] ([Year]). [DOI / preprint link].**
> _(fill in on submission/acceptance)_

## What the pipeline does

- Ingests raw Vectronic activity + GPS files, applies QC filtering, and derives
  two complementary responses: **movement intensity** (zero-inflated Beta) and
  **active duration** (beta-binomial).
- Extracts environmental covariates (distances to dump/roads/land-cover, terrain)
  and **ERA5-Land hourly temperature**; decorrelates them (PCA + VIF).
- Fits a hypothesis set **H0–H5** and compares them with LOO / grouped
  leave-bears-out cross-validation (LOBO).
- Derived analyses: nocturnality index, lunar control, conflict windows and a
  spatial conflict-risk surface, individual reaction norms, intra-specific
  dump partitioning, and a suite of robustness/sensitivity checks.

Each script header documents its purpose and, where applicable, **which
Supplementary Table (S*xx*) it produces**.

## Repository structure

```
.
├── 00_setup.R, 0_install.R        # shared config/paths/helpers; package install
├── 01_… 48_…                      # numbered analysis pipeline (R)
├── plot_predictors.py             # predictor-raster figure
├── plot_studyarea.py              # study-area land-cover map
├── sub_*.sh                       # SLURM submit scripts (HPC)
├── Dockerfile                     # the analysis container recipe (R 4.2.0)
├── EnvironmentalPreds/            # 12 GeoTIFF predictor rasters
├── shapefiles/                    # study area, roads, protected area, dump sites
├── example_data/                  # synthetic demo inputs (BTR_example) + README
└── data/README.md                 # where to obtain the external datasets
```

## Requirements

The pipeline is designed to run **only inside the analysis container** (R 4.2.0
with `brms`, `rstan`, `terra`, `sf`, `ctmm`, …). Build it from the `Dockerfile`:

```bash
docker build -t skbear-r:4.2.0 .
# On an HPC cluster, convert to Singularity/Apptainer:
apptainer build skbear_r.sif docker://<your-registry>/skbear-r:4.2.0
```

## Quick start (runnable demo, no cluster / no GIS)

Stages `01`–`02` run offline on the bundled **synthetic** example bear:

```bash
# Expected layout under PROJECT_ROOT:
#   $PROJECT_ROOT/BrownBearDielAct/scripts/   <- the scripts in this repo
#   $PROJECT_ROOT/data/                       <- inputs
# Place the example inputs (see example_data/README.md):
#   example_data/Activity/BTR_example.txt   -> data/Activity/BTR_example.txt
#   example_data/Position/BTR_example.txt   -> data/Position/BTR_example.txt
#   example_data/Bear_metadata.xlsx         -> data/Bear_metadata.xlsx

PROJECT_ROOT=/path/to/project Rscript 01_data_preparation.R
PROJECT_ROOT=/path/to/project Rscript 02_response_variables.R
```

This reproduces the ingestion → QC → response-variable flow. **Stages `03`
onward require the external covariate rasters and ERA5-Land temperature** (see
[`data/README.md`](data/README.md)) and were run on a SLURM cluster inside the
container (`sub_*.sh`).

## Data availability

- **`example_data/`** — small, fully synthetic inputs for one fictional bear
  (`BTR_example`); **not real animal data**, provided to demonstrate the input
  formats and run stages 01–02.
- **`EnvironmentalPreds/` + `shapefiles/`** — the environmental predictors and
  study-area vectors used by the models and figures.
- **External datasets** (ERA5-Land, E-OBS, ESRI land cover, WDPA, Natural Earth)
  are **not redistributed**; obtain them from the original providers as
  documented in [`data/README.md`](data/README.md).
- Real bear GPS/accelerometer data are sensitive (a legally protected species)
  and are available from the authors on reasonable request.

## Citation

If you use this code, please cite the manuscript above. _(Update once available.)_

## License

Released under the **MIT License** — see [`LICENSE`](LICENSE). Note that the
third-party datasets referenced in `data/README.md` retain their own licences
and terms.

## Contact

Maintainer: **Ercan Sıkdokur** ([@ercanskdokur](https://github.com/ercanskdokur)).
For questions, use GitHub or the corresponding-author details in the manuscript.
