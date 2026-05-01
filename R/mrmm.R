#' Fit a multi-resolution mixed model (MRMM) on the sphere
#'
#' Combines the [mrts_sphere()] basis with a tapered Matérn residual
#' covariance and obtains the basis coefficients by generalised least
#' squares (GLS).
#'
#' Given training observations \eqn{z} at locations \eqn{s_1,\ldots,s_n}
#' on the sphere, the model is
#' \deqn{z(s) = B(s)\,\beta + \varepsilon(s),\qquad
#'       \mathrm{Cov}(\varepsilon(s_i),\varepsilon(s_j)) =
#'       \Sigma_{ij},}{
#'   z(s) = B(s) beta + epsilon(s),  Cov(epsilon_i, epsilon_j) = Sigma_ij,}
#' where \eqn{B(s)} is the MRTS basis at location \eqn{s} and \eqn{\Sigma}
#' is the tapered-Matérn covariance built by
#' [tapered_matern_sphere()] using the user-supplied parameters. The
#' GLS estimate
#' \deqn{\hat{\beta} = (B^\top \Sigma^{-1} B)^{-1} B^\top \Sigma^{-1} z}{
#'   beta_hat = (B' Sigma^{-1} B)^{-1} B' Sigma^{-1} z}
#' is computed via a Cholesky factorisation of \eqn{\Sigma}. Predictions
#' (mean and standard error) at new locations are obtained with
#' [predict.mrmm()].
#'
#' Covariance parameters must be supplied. A separate parameter-estimation
#' routine (`mrmm_estimate_params()`) is planned for a future release.
#'
#' @param z Numeric vector of observations at the training locations.
#' @param loc Two-column numeric matrix of `(latitude, longitude)` in
#'   degrees giving the training locations. Must have `length(z)` rows.
#' @param knots Two-column numeric matrix of `(latitude, longitude)` in
#'   degrees giving the knot locations of the MRTS basis.
#' @param k Integer. Number of MRTS basis functions to use
#'   (`2 <= k <= nrow(knots)`).
#' @param taper_range Wendland compact-support range \eqn{a} (radians on
#'   the unit sphere).
#' @param matern_range Matérn range parameter \eqn{c}.
#' @param sill Sill of the structured residual component
#'   (\eqn{\sigma_y^2}).
#' @param nugget Nugget variance (\eqn{\sigma^2_\varepsilon}); diagonal
#'   added to the covariance.
#' @param smoothness Matérn smoothness \eqn{\nu}. Default `0.5`.
#'
#' @return An object of class `"mrmm"`: a list with components
#' \describe{
#'   \item{`coefficients`}{GLS estimate \eqn{\hat\beta}, length `k`.}
#'   \item{`fitted`, `residuals`}{Fitted values and residuals at the
#'     training locations.}
#'   \item{`coef_cov`}{The matrix \eqn{(B^\top\Sigma^{-1}B)^{-1}}, used by
#'     [predict.mrmm()] for standard errors.}
#'   \item{`sigma_chol`}{Upper-triangular Cholesky factor of \eqn{\Sigma}
#'     at the training locations.}
#'   \item{`basis_train`, `train_loc`, `train_z`}{Stored training
#'     basis, locations, and observations.}
#'   \item{`knots`, `k`, `params`, `call`}{Inputs preserved on the fit.}
#' }
#'
#' @seealso [predict.mrmm()], [mrts_sphere()],
#'   [tapered_matern_sphere()].
#'
#' @references
#' Huang, H.-Y., Huang, H.-C., and Ing, C.-K. (2025). Multi-Resolution
#' Spatial Methods on the Sphere: Efficient Prediction for Global Data.
#' *Environmetrics*. \doi{10.1002/env.70092}
#'
#' @examples
#' set.seed(1)
#' n <- 80
#' loc <- cbind(runif(n, -60, 60), runif(n, -180, 180))
#' z   <- 2 + sin(loc[, 1] * pi / 180) + rnorm(n, sd = 0.3)
#'
#' knots <- loc[sample.int(n, 30), ]
#' fit <- mrmm_fit(z = z, loc = loc, knots = knots, k = 8,
#'                 taper_range  = 1.0,
#'                 matern_range = 0.2,
#'                 sill         = 0.5,
#'                 nugget       = 0.1,
#'                 smoothness   = 0.5)
#' fit
#' length(fit$coefficients)   # k
#'
#' @export
mrmm_fit <- function(z, loc, knots, k,
                     taper_range, matern_range,
                     sill, nugget, smoothness = 0.5) {
    z <- as.numeric(z)
    loc <- as.matrix(loc)
    knots <- as.matrix(knots)

    if (ncol(loc) != 2L)   stop("`loc` must be a 2-column matrix of (lat, lon).")
    if (ncol(knots) != 2L) stop("`knots` must be a 2-column matrix of (lat, lon).")
    if (length(z) != nrow(loc)) {
        stop("length(z) must equal nrow(loc).")
    }
    k <- as.integer(k)
    if (k < 2L || k > nrow(knots)) {
        stop("`k` must be an integer in [2, nrow(knots)].")
    }

    # MRTS basis at training locations.
    B <- mrts_sphere(knot = knots, k = k, X = loc)$mrts

    # Tapered Matérn residual covariance at training locations.
    Sigma <- tapered_matern_sphere(loc,
        taper_range  = taper_range,
        matern_range = matern_range,
        smoothness   = smoothness,
        sill         = sill,
        nugget       = nugget)

    R <- chol(Sigma)
    Sinv_B <- backsolve(R, forwardsolve(t(R), B))
    Sinv_z <- backsolve(R, forwardsolve(t(R), z))

    BtSinvB  <- crossprod(B, Sinv_B)
    BtSinvz  <- crossprod(B, Sinv_z)
    coef_cov <- solve(BtSinvB)
    beta_hat <- drop(coef_cov %*% BtSinvz)

    fitted <- as.numeric(B %*% beta_hat)
    resid  <- z - fitted

    structure(
        list(
            coefficients = beta_hat,
            fitted       = fitted,
            residuals    = resid,
            coef_cov     = coef_cov,
            sigma_chol   = R,
            basis_train  = B,
            train_loc    = loc,
            train_z      = z,
            knots        = knots,
            k            = k,
            params       = list(
                taper_range  = taper_range,
                matern_range = matern_range,
                sill         = sill,
                nugget       = nugget,
                smoothness   = smoothness
            ),
            call         = match.call()
        ),
        class = "mrmm"
    )
}

#' @export
print.mrmm <- function(x, ...) {
    cat("Multi-resolution mixed model (mrtsSphere)\n")
    cat(sprintf("  Training points : %d\n", length(x$train_z)))
    cat(sprintf("  Knots           : %d\n", nrow(x$knots)))
    cat(sprintf("  Basis rank (k)  : %d\n", x$k))
    cat("Tapered-Matern parameters:\n")
    p <- x$params
    cat(sprintf("  taper_range  = %g (rad)\n", p$taper_range))
    cat(sprintf("  matern_range = %g\n",        p$matern_range))
    cat(sprintf("  sill         = %g\n",        p$sill))
    cat(sprintf("  nugget       = %g\n",        p$nugget))
    cat(sprintf("  smoothness   = %g\n",        p$smoothness))
    invisible(x)
}
