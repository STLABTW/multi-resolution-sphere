#' Predict from a multi-resolution mixed model fit
#'
#' Computes kriging-style predictions and (optionally) prediction
#' standard errors at new spherical locations from an [mrmm_fit()]
#' object.
#'
#' Given the GLS coefficient \eqn{\hat\beta} stored in `object`, the
#' prediction at a new location \eqn{s^*} is
#' \deqn{\hat z(s^*) = B(s^*)\hat\beta + \Sigma(s^*, S_{\mathrm{train}})\,
#'                     \Sigma_{\mathrm{train}}^{-1}\,(z - B\hat\beta).}{
#'   z_hat(s*) = B(s*) beta_hat + Sigma(s*, S_train) Sigma_train^{-1} (z - B beta_hat).}
#' The prediction variance is the closed-form expression that arises in
#' kriging with linear-trend basis (sometimes called "universal
#' kriging"); see the reference.
#'
#' @param object An object of class `"mrmm"` from [mrmm_fit()].
#' @param newdata Two-column numeric matrix of `(latitude, longitude)`
#'   in degrees giving the prediction locations.
#' @param se.fit Logical; if `TRUE` (default), also return the
#'   prediction standard error at each new location.
#' @param ... Currently unused.
#'
#' @return A list with
#' \describe{
#'   \item{`fit`}{Numeric vector of predicted values, one per row of
#'     `newdata`.}
#'   \item{`se.fit`}{Numeric vector of prediction standard errors (only
#'     if `se.fit = TRUE`).}
#' }
#'
#' @seealso [mrmm_fit()].
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
#' knots <- loc[sample.int(n, 30), ]
#' fit <- mrmm_fit(z = z, loc = loc, knots = knots, k = 8,
#'                 taper_range = 1.0, matern_range = 0.2,
#'                 sill = 0.5, nugget = 0.1)
#'
#' new_loc <- cbind(seq(-50, 50, length.out = 5),
#'                  seq(-100, 100, length.out = 5))
#' p <- predict(fit, newdata = new_loc, se.fit = TRUE)
#' p$fit
#' p$se.fit
#'
#' @export
predict.mrmm <- function(object, newdata, se.fit = TRUE, ...) {
    newdata <- as.matrix(newdata)
    if (ncol(newdata) != 2L) {
        stop("`newdata` must be a 2-column matrix of (lat, lon).")
    }
    p <- object$params

    B_new <- mrts_sphere(knot = object$knots, k = object$k, X = newdata)$mrts

    C_cross <- tapered_matern_sphere(newdata, object$train_loc,
        taper_range  = p$taper_range,
        matern_range = p$matern_range,
        smoothness   = p$smoothness,
        sill         = p$sill,
        nugget       = 0)

    R <- object$sigma_chol
    Sinv_resid <- backsolve(R, forwardsolve(t(R), object$residuals))
    mean_pred <- as.numeric(B_new %*% object$coefficients + C_cross %*% Sinv_resid)

    if (!isTRUE(se.fit)) {
        return(list(fit = mean_pred))
    }

    Sinv_B  <- backsolve(R, forwardsolve(t(R), object$basis_train))
    Sinv_Ct <- backsolve(R, forwardsolve(t(R), t(C_cross)))

    f2 <- B_new - C_cross %*% Sinv_B

    term1 <- rowSums((f2 %*% object$coef_cov) * f2)
    term2 <- rowSums(C_cross * t(Sinv_Ct))

    var_pred <- term1 - term2 + p$sill
    var_pred[var_pred < 0] <- 0

    list(
        fit    = mean_pred,
        se.fit = sqrt(var_pred)
    )
}
