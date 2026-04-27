#' mrtsSphere: Multi-Resolution Thin-Plate Splines on the Sphere
#'
#' Constructs multi-resolution thin-plate spline basis functions on the
#' sphere for use in spatial regression and large-scale spatial prediction
#' problems.
#'
#' The main user-facing function is [mrts_sphere()].
#'
#' @references
#' Multi-resolution approximations of Gaussian processes for large
#' spatial datasets on the sphere. *Environmetrics*, 2025.
#' \doi{10.1002/env.70092}
#'
#' @keywords internal
#' @useDynLib mrtsSphere, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RSpectra eigs_sym
"_PACKAGE"
