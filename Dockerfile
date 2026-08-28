# ==============================================================================
# Dockerfile — SkBear R 4.2.0 Pipeline
# ==============================================================================
# Base: rocker/geospatial:4.2.0
#   - Ubuntu 20.04 LTS, R 4.2.0
#   - includes GDAL 3.4+, GEOS, PROJ, libudunits2
#   - terra, sf, sp, raster already installed
#
# This is the canonical build recipe for the analysis container. The pipeline is
# designed to run ONLY inside this container. Reproduce the image from this file:
#
# Build (Docker):
#   docker build -t skbear-r:4.2.0 .
# Convert to a Singularity/Apptainer image on the HPC cluster:
#   apptainer build skbear_r.sif docker://<your-registry>/skbear-r:4.2.0
# (or push the local image to your own registry namespace first).
# ==============================================================================

FROM rocker/geospatial:4.2.0

LABEL description="SkBear bear activity analysis pipeline — R 4.2.0"

# ------------------------------------------------------------------------------
# 1. Extra system dependencies (rocker/geospatial already ships the spatial libs)
# ------------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    # ctmm dependencies
    libfftw3-dev \
    libgmp-dev \
    libmpfr-dev \
    # NetCDF (ncdf4)
    libnetcdf-dev \
    libhdf5-dev \
    # readxl / writexl
    libxslt1-dev \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------------------
# 2. Stan / rstan C++14 toolchain
# ------------------------------------------------------------------------------
RUN mkdir -p /etc/R && \
    printf 'CXX14 = g++ -std=c++14 -fPIC\nCXX14FLAGS = -O3 -fPIC\nCXX14PICFLAGS = -fPIC\n' \
    >> /etc/R/Makevars.site

# ------------------------------------------------------------------------------
# 3. R packages (rocker/geospatial already provides terra, sf, sp, raster)
# ------------------------------------------------------------------------------

# 3a. Core data wrangling
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/focal/latest')); \
  install.packages(c( \
    'conflicted', \
    'tidyverse', \
    'lubridate', 'readxl', 'data.table', 'writexl', \
    'tibble', 'stringr', 'purrr', 'readr', 'forcats' \
  ), dependencies = TRUE, Ncpus = 4L); \
  if (!'tidyverse' %in% rownames(installed.packages())) stop('tidyverse not installed')"

# 3b. NetCDF + circadian
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/focal/latest')); \
  install.packages(c('ncdf4', 'mixtools', 'circular', 'suncalc'), \
    dependencies = TRUE, Ncpus = 4L); \
  if (!'mixtools' %in% rownames(installed.packages())) stop('mixtools not installed')"

# 3c. rstan (longest step — its own layer)
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/focal/latest')); \
  install.packages('rstan', dependencies = TRUE, Ncpus = 4L); \
  if (!'rstan' %in% rownames(installed.packages())) stop('rstan not installed')"

# 3d. brms + Bayesian tooling
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/focal/latest')); \
  install.packages(c('brms', 'tidybayes', 'loo', 'bayesplot', 'posterior'), \
    dependencies = TRUE, Ncpus = 4L); \
  if (!'brms' %in% rownames(installed.packages())) stop('brms not installed')"

# 3e. GPS movement model — dependencies=NA: does not reinstall existing pkgs (terra/sf)
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/focal/latest')); \
  install.packages('ctmm', dependencies = NA, Ncpus = 4L); \
  if (!'ctmm' %in% rownames(installed.packages())) stop('ctmm not installed')"

# 3f. Visualization + parallel
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/focal/latest')); \
  install.packages(c( \
    'ggplot2', 'patchwork', 'scales', 'ggridges', \
    'gridExtra', 'corrplot', \
    'future', 'future.apply' \
  ), dependencies = TRUE, Ncpus = 4L)"

# ------------------------------------------------------------------------------
# 4. Verify installation — build FAILS if anything is missing
# ------------------------------------------------------------------------------
RUN Rscript -e "\
  pkgs <- c( \
    'conflicted','tidyverse','lubridate','readxl','data.table','writexl', \
    'terra','sf','sp','raster','ctmm','ncdf4', \
    'mixtools','circular','suncalc', \
    'rstan','brms','tidybayes','loo','bayesplot','posterior', \
    'ggplot2','patchwork','scales','ggridges','gridExtra','corrplot', \
    'future','future.apply' \
  ); \
  inst <- rownames(installed.packages(lib.loc = .libPaths())); \
  missing <- pkgs[!pkgs %in% inst]; \
  if (length(missing) > 0) { \
    cat('MISSING PACKAGES:', paste(missing, collapse=', '), '\n'); \
    quit(status = 1) \
  } else { \
    cat('=== All packages OK ===\n') \
  }"

# ------------------------------------------------------------------------------
# 5. Working directory
# ------------------------------------------------------------------------------
WORKDIR /workspace

CMD ["R", "--no-save"]
