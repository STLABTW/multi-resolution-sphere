test_that("mrmm_fit returns the expected S3 object and shapes", {
    set.seed(1)
    n <- 60
    loc <- cbind(runif(n, -60, 60), runif(n, -180, 180))
    z   <- 1 + 0.5 * sin(loc[, 1] * pi / 180) + rnorm(n, sd = 0.2)
    knots <- loc[sample.int(n, 25), ]

    fit <- mrmm_fit(z = z, loc = loc, knots = knots, k = 6,
                    taper_range = 1.0, matern_range = 0.2,
                    sill = 0.4, nugget = 0.1)

    expect_s3_class(fit, "mrmm")
    expect_length(fit$coefficients, 6L)
    expect_length(fit$fitted, n)
    expect_length(fit$residuals, n)
    expect_equal(dim(fit$coef_cov), c(6L, 6L))
})

test_that("predict.mrmm returns mean and non-negative SEs", {
    set.seed(2)
    n <- 50
    loc <- cbind(runif(n, -50, 50), runif(n, -160, 160))
    z   <- rnorm(n)
    knots <- loc[sample.int(n, 20), ]
    fit <- mrmm_fit(z, loc, knots, k = 5,
                    taper_range = 1.2, matern_range = 0.3,
                    sill = 0.5, nugget = 0.2)

    new_loc <- cbind(seq(-40, 40, length.out = 8),
                     seq(-150, 150, length.out = 8))

    p <- predict(fit, newdata = new_loc, se.fit = TRUE)
    expect_named(p, c("fit", "se.fit"))
    expect_length(p$fit, nrow(new_loc))
    expect_length(p$se.fit, nrow(new_loc))
    expect_true(all(p$se.fit >= 0))

    p_mean_only <- predict(fit, newdata = new_loc, se.fit = FALSE)
    expect_named(p_mean_only, "fit")
})

test_that("with sill = 0 (Sigma = nugget * I) GLS reduces to OLS", {
    set.seed(3)
    n <- 40
    loc <- cbind(runif(n, -50, 50), runif(n, -150, 150))
    z   <- rnorm(n)
    knots <- loc[sample.int(n, 18), ]
    k <- 5

    fit_gls <- mrmm_fit(z, loc, knots, k = k,
                        taper_range = 1.0, matern_range = 0.2,
                        sill = 0, nugget = 0.5)

    B   <- mrts_sphere(knots, k = k, X = loc)$mrts
    ols <- as.numeric(qr.solve(crossprod(B), crossprod(B, z)))
    expect_equal(as.numeric(fit_gls$coefficients), ols, tolerance = 1e-8)
})

test_that("tapered_matern_sphere is symmetric, PSD-ish, and respects support", {
    set.seed(4)
    loc <- cbind(runif(15, -45, 45), runif(15, -120, 120))
    Sigma <- tapered_matern_sphere(loc,
        taper_range = 0.5, matern_range = 0.1,
        smoothness = 0.5, sill = 1.0, nugget = 0.05)
    expect_true(isSymmetric(Sigma))
    expect_true(min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) > 0)

    # With a tiny taper, off-diagonals of a fresh draw vanish.
    Sigma_tiny <- tapered_matern_sphere(loc,
        taper_range = 1e-6, matern_range = 0.1,
        sill = 1.0, nugget = 0.0)
    expect_equal(Sigma_tiny[upper.tri(Sigma_tiny)],
                 rep(0, sum(upper.tri(Sigma_tiny))))
})
