#' Pairwise great-circle distances on the unit sphere
#'
#' @param loc1,loc2 Two-column matrices of `(latitude, longitude)` in
#'   degrees. If `loc2` is `NULL`, distances are computed within `loc1`.
#'
#' @return An `nrow(loc1) x nrow(loc2)` matrix of central angles in
#'   radians (equivalently, great-circle distances on the unit sphere).
#'
#' @keywords internal
gcdist_sphere <- function(loc1, loc2 = NULL) {
    if (is.null(loc2)) loc2 <- loc1
    rad <- pi / 180
    lat1 <- loc1[, 1] * rad; lon1 <- loc1[, 2] * rad
    lat2 <- loc2[, 1] * rad; lon2 <- loc2[, 2] * rad
    cos_d <- outer(sin(lat1), sin(lat2)) +
        outer(cos(lat1), cos(lat2)) * cos(outer(lon1, lon2, "-"))
    cos_d[cos_d >  1] <-  1
    cos_d[cos_d < -1] <- -1
    acos(cos_d)
}

# Matérn correlation as a function of central angle d (radians).
# rho(d) = (2^(1-nu) / gamma(nu)) * (d/range)^nu * BesselK(d/range, nu)
# Reduces to exp(-d/range) when nu == 0.5.
matern_corr_internal <- function(d, range, nu) {
    if (isTRUE(all.equal(nu, 0.5))) {
        return(exp(-d / range))
    }
    x <- d / range
    out <- x
    zero <- (x == 0)
    out[zero] <- 1
    if (any(!zero)) {
        out[!zero] <- (2^(1 - nu) / gamma(nu)) *
            x[!zero]^nu * besselK(x[!zero], nu)
    }
    out
}

#' Tapered Matérn covariance on the sphere
#'
#' Builds a covariance matrix of the form
#'
#' \deqn{\mathrm{sill} \cdot \rho_{\mathrm{Mat}}(d;\, c, \nu) \cdot
#'       W(d;\, a) \;+\; \mathrm{nugget} \cdot I}{
#'   sill * Matern(d; c, nu) * Wendland(d; a) + nugget * I}
#'
#' where \eqn{d} is the great-circle distance on the unit sphere
#' (radians), \eqn{\rho_\mathrm{Mat}} is the Matérn correlation with
#' range \eqn{c} and smoothness \eqn{\nu}, and \eqn{W} is the
#' Wendland \eqn{C^2} taper of compact support \eqn{a} (degree
#' \eqn{k = 2}, dimension \eqn{= 2}). The nugget term is added only when
#' the matrix is symmetric, i.e. when `loc2 = NULL`.
#'
#' @param loc1 Two-column numeric matrix of `(latitude, longitude)`
#'   locations in degrees.
#' @param loc2 Optional second set of locations. If `NULL`, the
#'   covariance is computed within `loc1` (symmetric output) and the
#'   nugget is added on the diagonal.
#' @param taper_range Compact-support range \eqn{a} of the Wendland
#'   taper, in radians on the unit sphere.
#' @param matern_range Range parameter \eqn{c} of the Matérn correlation.
#' @param smoothness Smoothness \eqn{\nu} of the Matérn correlation.
#'   Default `0.5` (exponential).
#' @param sill Marginal variance of the structured component. Must be
#'   non-negative.
#' @param nugget Nugget variance, added on the diagonal when
#'   `loc2 = NULL`. Must be non-negative.
#'
#' @return An `nrow(loc1) x nrow(loc2)` numeric covariance matrix.
#'
#' @examples
#' set.seed(1)
#' loc <- cbind(runif(20, -60, 60), runif(20, -180, 180))
#' Sigma <- tapered_matern_sphere(loc,
#'     taper_range = 0.4,
#'     matern_range = 0.1,
#'     smoothness = 0.5,
#'     sill = 1.0,
#'     nugget = 0.05)
#' isSymmetric(Sigma)
#'
#' @export
tapered_matern_sphere <- function(loc1, loc2 = NULL,
                                  taper_range, matern_range,
                                  smoothness = 0.5,
                                  sill = 1, nugget = 0) {
    loc1 <- as.matrix(loc1)
    if (ncol(loc1) != 2L) {
        stop("`loc1` must be a 2-column matrix of (lat, lon).")
    }
    symmetric <- is.null(loc2)
    if (symmetric) {
        loc2_eff <- loc1
    } else {
        loc2_eff <- as.matrix(loc2)
        if (ncol(loc2_eff) != 2L) {
            stop("`loc2` must be a 2-column matrix of (lat, lon).")
        }
    }
    if (!(taper_range > 0)) stop("`taper_range` must be > 0.")
    if (!(matern_range > 0)) stop("`matern_range` must be > 0.")
    if (!(smoothness > 0))   stop("`smoothness` must be > 0.")
    if (sill < 0)   stop("`sill` must be non-negative.")
    if (nugget < 0) stop("`nugget` must be non-negative.")

    d <- gcdist_sphere(loc1, loc2_eff)

    # Wendland C^2 taper, dimension 2:  (1 - r)^6 * (35/3 r^2 + 6 r + 1)
    r <- d / taper_range
    taper <- (1 - r)^6 * ((35 / 3) * r^2 + 6 * r + 1)
    taper[r > 1] <- 0

    matern <- matern_corr_internal(d, matern_range, smoothness)

    Sigma <- sill * matern * taper

    if (symmetric && nugget > 0) {
        diag(Sigma) <- diag(Sigma) + nugget
    }
    Sigma
}
