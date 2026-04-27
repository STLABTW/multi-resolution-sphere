#' Multi-resolution thin-plate spline basis on the sphere
#'
#' Builds a set of `k` multi-resolution thin-plate spline (MRTS) basis
#' functions on the sphere from a set of knot locations, and evaluates
#' them at the prediction locations `X`.
#'
#' The first basis function is constant (`sqrt(1/n)`); the remaining
#' `k - 1` basis functions are obtained from the eigen-decomposition of
#' the centered knot kernel matrix, following the construction described
#' in the reference.
#'
#' @param knot Numeric matrix with two columns giving knot locations as
#'   `(latitude, longitude)` in degrees.
#' @param k Integer. Number of basis functions to construct (the rank of
#'   the basis). Must satisfy `2 <= k <= nrow(knot)`.
#' @param X Numeric matrix with two columns giving prediction locations
#'   as `(latitude, longitude)` in degrees, where the basis is evaluated.
#'
#' @return A list with one element:
#' \describe{
#'   \item{`mrts`}{An `nrow(X) x k` numeric matrix whose columns are the
#'     basis functions evaluated at the rows of `X`.}
#' }
#'
#' @references
#' Multi-resolution approximations of Gaussian processes for large
#' spatial datasets on the sphere. *Environmetrics*, 2025.
#' \doi{10.1002/env.70092}
#'
#' @examples
#' ## Build a small global grid in (lat, lon) degrees.
#' n_lon <- 12
#' n_lat <- 8
#' lon_seq <- seq(-180, 150, length.out = n_lon)
#' lat_seq <- seq( -80,  80, length.out = n_lat)
#' grid <- as.matrix(expand.grid(lat = lat_seq, lon = lon_seq))
#'
#' ## Pick 30 knots and evaluate the MRTS basis at every grid point.
#' set.seed(1)
#' knots <- grid[sample(nrow(grid), 30), ]
#' res <- mrts_sphere(knots, k = 5, X = grid)
#' dim(res$mrts)   # nrow(grid) x k
#'
#' ## Recovering a simulated spherical exponential field with the basis.
#' if (requireNamespace("fields", quietly = TRUE)) {
#'   # Great-circle distance (km) -> exponential covariance.
#'   d_grid    <- fields::rdist.earth(grid[, 2:1], miles = FALSE)
#'   cov_field <- exp(-d_grid / 2000)
#'   y_true    <- as.numeric(t(chol(cov_field + diag(1e-8, nrow(grid)))) %*%
#'                           rnorm(nrow(grid)))
#'
#'   # Noisy observations at the knot locations.
#'   obs_idx <- match(data.frame(t(knots)), data.frame(t(grid)))
#'   z_obs   <- y_true[obs_idx] + rnorm(nrow(knots), sd = 0.3)
#'
#'   # Project into the MRTS basis (least squares) and predict on the grid.
#'   B_obs    <- res$mrts[obs_idx, , drop = FALSE]
#'   beta_hat <- qr.solve(B_obs, z_obs)
#'   y_hat    <- res$mrts %*% beta_hat
#'   sqrt(mean((y_hat - y_true)^2))   # RMSE
#' }
#'
#' @export
mrts_sphere <- function(knot, k, X) {
    if (!is.matrix(knot) || ncol(knot) != 2L) {
        stop("`knot` must be a 2-column numeric matrix of (lat, lon).")
    }
    if (!is.matrix(X) || ncol(X) != 2L) {
        stop("`X` must be a 2-column numeric matrix of (lat, lon).")
    }
    n <- nrow(knot)
    N <- nrow(X)
    k <- as.integer(k)
    if (k < 2L || k > n) {
        stop("`k` must be an integer in [2, nrow(knot)].")
    }

    onev <- rep(1 / n, n)
    K <- cpp_K(knot[, 1], knot[, 2], n)
    Q <- diag(1, n, n) - (1 / n)
    eiK <- RSpectra::eigs_sym(Q %*% K %*% Q, k)

    eiKvecmval <- matrix(0, n, k)
    for (i in seq_len(k)) {
        eiKvecmval[, i] <- eiK$vectors[, i] / eiK$values[i]
    }

    dm_train <- cpp_Kmatrix6(
        KK = k,
        X = knot,
        ggrids = X,
        Konev = as.numeric(K %*% onev),
        eiKvecmval = eiKvecmval,
        n = n,
        N = N,
        aa_grid = .integral_table$aa_grid,
        tab_vals = .integral_table$tab_vals
    )

    list(mrts = dm_train)
}
