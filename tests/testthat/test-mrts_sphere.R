test_that("mrts_sphere returns expected shape and first basis is constant", {
    set.seed(42)
    lon_seq <- seq(-180, 150, length.out = 12)
    lat_seq <- seq(-80, 80, length.out = 8)
    grid <- as.matrix(expand.grid(lat = lat_seq, lon = lon_seq))

    n_obs <- 25
    knots <- grid[sample(nrow(grid), n_obs), ]

    k <- 5
    res <- mrts_sphere(knots, k = k, X = grid)

    expect_type(res, "list")
    expect_named(res, "mrts")
    expect_equal(dim(res$mrts), c(nrow(grid), k))

    # First basis function is the constant sqrt(1/n_knot).
    expect_equal(res$mrts[, 1], rep(sqrt(1 / n_obs), nrow(grid)))
})

test_that("mrts_sphere validates inputs", {
    set.seed(1)
    knots <- matrix(runif(20, -80, 80), ncol = 2)
    grid  <- matrix(runif(40, -80, 80), ncol = 2)

    expect_error(mrts_sphere(as.numeric(knots), 3, grid), "knot")
    expect_error(mrts_sphere(knots, 3, as.numeric(grid)), "X")
    expect_error(mrts_sphere(knots, 1L, grid),  "k.*\\[2,")
    expect_error(mrts_sphere(knots, nrow(knots) + 1L, grid), "k.*\\[2,")
})
