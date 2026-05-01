# mrtsSphere

[![CRAN status](https://www.r-pkg.org/badges/version/mrtsSphere)](https://cran.r-project.org/package=mrtsSphere)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/mrtsSphere)](https://cran.r-project.org/package=mrtsSphere)
[![License: GPL v2+](https://img.shields.io/badge/License-GPL%20v2%2B-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)

Multi-resolution thin-plate spline (MRTS) basis functions on the sphere
for large-scale spatial regression and prediction. R implementation of
the method in:

> Huang, H.-Y., Huang, H.-C., and Ing, C.-K. (2025).
> **Multi-Resolution Spatial Methods on the Sphere: Efficient Prediction
> for Global Data.** *Environmetrics*.
> DOI: [10.1002/env.70092](https://doi.org/10.1002/env.70092)

The basis is constructed from the eigen-decomposition of a centered
spherical kernel and is evaluated on the prediction grid via a
parallel C++ routine (Rcpp + optional OpenMP).

## Installation

From CRAN (recommended):

```r
install.packages("mrtsSphere")
```

Development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("STLABTW/multi-resolution-sphere")
```

The package compiles C++ code on installation; you need a working
toolchain (Xcode CLT on macOS, Rtools on Windows, `r-base-dev` on
Linux). OpenMP is optional — without it the package still works,
single-threaded.

## Quick start

```r
library(mrtsSphere)

# Build a 20 x 10 grid of (lat, lon) locations on the sphere.
n_lon  <- 20
n_lat  <- 10
lon_seq <- seq(-180, 176, length.out = n_lon)
lat_seq <- seq( -90,  87, length.out = n_lat)
grid <- as.matrix(expand.grid(lat = lat_seq, lon = lon_seq))

# Pick 100 knots at random.
set.seed(1)
knots <- grid[sample(nrow(grid), 100), ]

# 10 multi-resolution basis functions evaluated on the full grid.
res <- mrts_sphere(knots, k = 10, X = grid)
dim(res$mrts)   # 200 x 10
```

For a full multi-resolution mixed model (MRMM) fit — basis +
GLS-on-tapered-Matérn-residuals + kriging predictions with prediction
SE — use `mrmm_fit()`:

```r
# z is the vector of observations at the rows of `obs_loc` (lat, lon).
fit <- mrmm_fit(z = z, loc = obs_loc, knots = obs_loc, k = 50,
                taper_range  = 0.2,   # radians on unit sphere
                matern_range = 0.07,
                sill         = 1.6,
                nugget       = 0.02,
                smoothness   = 0.5)

pred <- predict(fit, newdata = grid, se.fit = TRUE)
pred$fit      # mean prediction at each grid row
pred$se.fit   # prediction standard error
```

The covariance parameters above are user-supplied. A future release
will add `mrmm_estimate_params()`; for now you can estimate them using
`geoR`, `gstat`, or the variogram-based scripts under
`inst/paper/`.

A longer worked example that simulates a spherical Gaussian random
field with `fields` and recovers it through both OLS and MRMM ships
with the package at:

```r
system.file("articles", "mrtsSphere.Rmd", package = "mrtsSphere")
```

(Render with `rmarkdown::render()` if you have pandoc installed.)

## Reproducing the paper

The original analysis scripts from the paper are bundled under
`inst/paper/`. After installing the package:

```r
file.path(system.file("paper", package = "mrtsSphere"), "fullmodel-max.R")
```

The SST input dataset (`data_sst_max_20240419.csv`) is tracked via
Git LFS in [this repository](https://github.com/STLABTW/multi-resolution-sphere).

## Citation

```r
citation("mrtsSphere")
```

## License

GPL (>= 2). See <https://www.gnu.org/licenses/gpl-2.0> for the full
license text.
