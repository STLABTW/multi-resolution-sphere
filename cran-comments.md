# cran-comments.md

## Submission (0.2.0, feature release)

This is a feature release of the CRAN-resident package `mrtsSphere`.

New in 0.2.0:

- `mrmm_fit()` fits the multi-resolution mixed model (MRMM) by GLS
  with a tapered Matérn residual covariance, given user-supplied
  covariance parameters.
- `predict.mrmm()` returns kriging-style predictions and prediction
  standard errors.
- `tapered_matern_sphere()` is exported for users who want direct
  access to the covariance.

A future release will add a parameter-estimation routine
(`mrmm_estimate_params()`); the present release deliberately keeps the
covariance parameters as user inputs to limit scope.

The package implements the multi-resolution thin-plate spline (MRTS)
basis on the sphere from Huang, Huang, and Ing (2025), *Environmetrics*,
<doi:10.1002/env.70092>.

## Test environments

- macOS (Apple silicon), R 4.5.0 — `R CMD check --as-cran`
  - 0 errors | 0 warnings | 0 notes
- win-builder R-release (R 4.6.0)
  - 0 errors | 0 warnings | 1 note (see below)
- win-builder R-devel (R-devel, 2026-04-26 r89963)
  - 0 errors | 0 warnings | 1 note (see below)

## Notes from win-builder

Both Windows builds report:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Hao-Yun Huang <hhuscout@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  Environmetrics (17:52)
  Ing (16:23)
```

- "New submission" — this is the package's first CRAN upload.
- "Environmetrics" — the name of the journal in which the underlying
  paper was published (Wiley journal *Environmetrics*).
- "Ing" — surname of co-author Ching-Kang Ing.

Both flagged words are intentional and correctly spelled.

## Notes for reviewers

- The package compiles 'C++' code via 'Rcpp' / 'RcppEigen' /
  'RcppNumerical'. 'OpenMP' is optional and guarded with `#ifdef
  _OPENMP`; the standard `SHLIB_OPENMP_*FLAGS` are used in
  `src/Makevars` and `src/Makevars.win` so the package builds on systems
  without OpenMP support (e.g. the default macOS toolchain).
- A precomputed kernel-integral lookup table is shipped as internal data
  in `R/sysdata.rda` (~ 0.5 MB). It is regenerable from
  `data-raw/prepare_sysdata.R`.
- The `inst/paper/` directory contains the original analysis scripts
  used to produce the figures in the paper. They are *not* part of the
  package and are not run by `R CMD check`.

## Reverse dependencies

There are no reverse dependencies (first submission).
