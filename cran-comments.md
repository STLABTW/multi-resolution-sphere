# cran-comments.md

## Submission

This is the first CRAN submission of `mrtsSphere`.

The package implements the multi-resolution thin-plate spline (MRTS)
basis on the sphere from Huang, Huang, and Ing (2025), *Environmetrics*,
<doi:10.1002/env.70092>.

## Test environments

- macOS (Apple silicon), R 4.5.0 — `R CMD check --as-cran`
  - 0 errors | 0 warnings | 0 notes (run after first GitHub push)
- (To run before submission: `devtools::check_win_devel()` and
  `rhub::rhub_check()` for Linux/Windows/macOS-cran flavours.)

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
