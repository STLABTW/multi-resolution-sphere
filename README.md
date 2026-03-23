# Multi-Resolution Model on the Sphere

Spherical multi-resolution regression applied to global sea-surface temperature (SST) data.

## Overview

`fullmodel-max.R` fits the model to the SST annual-maximum dataset and produces:
- predicted SST map (Mollweide projection)
- prediction standard-error map
- raw-data map

The method uses multi-resolution thin-plate splines (MRTS) on the sphere combined with a locally-supported Matérn covariance (Wendland tapering).

## Files

| File | Description |
|------|-------------|
| `fullmodel-max.R` | Main analysis script for the annual-maximum dataset |
| `fn_pcc_test_pre.R` | C++ kernel functions (compiled via Rcpp) and MRTS helpers |
| `fn_0610.R` | Variogram objective functions for parameter estimation |
| `effectivefn.R` | Matérn effective-range calculator |
| `integral_table2.rds` | Pre-computed lookup table for the kernel integrals |
| `data_sst_max_20240419.csv` | SST annual-maximum dataset (tracked via Git LFS) |

## Data

`data_sst_max_20240419.csv` contains ~6.4 million ocean-grid observations with columns:
`latitude`, `longitude`, `temperature` (°C).

The file is stored in this repository via **Git LFS**. Clone with LFS support:
```bash
git lfs install
git clone <repo-url>
```

## Dependencies (R packages)

```r
install.packages(c(
  "icosa", "fields", "sf", "rnaturalearth", "rnaturalearthdata",
  "pracma", "raster", "maps", "SparseM", "ggplot2",
  "RSpectra", "Rcpp", "RcppArmadillo", "RcppEigen",
  "RcppNumerical", "GpGp", "MASS", "Matrix"
))
```

## Usage

Set the working directory to the repo root and run:
```r
source("fullmodel-max.R")
```

Output PNG files are saved to `realdata/`.
