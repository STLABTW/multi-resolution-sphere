# Paper-companion scripts

These scripts reproduce the analyses in:

> Multi-resolution approximations of Gaussian processes for large
> spatial datasets on the sphere. *Environmetrics*, 2025.
> [doi:10.1002/env.70092](https://doi.org/10.1002/env.70092)

They are bundled with the `mrtsSphere` package for reference but are
**not** part of the package itself — they are not tested by
`R CMD check` and are not automatically loaded.

## Files

| File | Description |
|------|-------------|
| `fullmodel-max.R` | Main analysis script for the SST annual-maximum dataset |
| `fn_pcc_test_pre.R` | Original C++ kernel functions and helpers (now in `src/mrts.cpp`) |
| `fn_0610.R` | Variogram objective functions for parameter estimation |
| `effectivefn.R` | Matérn effective-range calculator |

## Data

`data_sst_max_20240419.csv` is tracked via Git LFS in the repository
root. To run the paper scripts, fetch it with:

```bash
git lfs install
git lfs pull
```

## Running

```r
# from the package source tree
setwd(system.file("paper", package = "mrtsSphere"))
source("fullmodel-max.R")
```
