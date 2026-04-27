# You need this for the gamma function and the Bessel function
# No external packages needed

# Matérn correlation function
matern_cor <- function(h, range, nu) {
  h <- as.numeric(h)
  h[h == 0] <- 1e-10 # avoid division by zero
  part1 <- (2^(1 - nu)) / gamma(nu)
  part2 <- (h / range)^nu
  part3 <- besselK(h / range, nu)
  return(part1 * part2 * part3)
}

# Effective range calculator: distance where corr drops to target
matern_effective_range <- function(range, nu, target_cor = 0.05) {
  # f(d) = |corr(d) - target|
  f <- function(d) abs(matern_cor(d, range, nu) - target_cor)
  out <- optimize(f, interval = c(1e-8, 50))  # Search over a reasonable interval
  return(out$minimum)
}

matern_eff_range <- function(range, nu, corlevel = 0.05) {
  # Solve for h: matern_cor(h) = corlevel
  # For exponential: h = -range * log(corlevel)
  if (nu == 0.5) {
    h_eff <- -range * log(corlevel)
  } else {
    matern_cor <- function(h) {
      x <- h / range
      factor <- 1 / (2^(nu-1) * gamma(nu))
      factor * (x^nu) * besselK(x, nu)
    }
    f <- function(h) matern_cor(h) - corlevel
    h_eff <- tryCatch(uniroot(f, lower = 1e-8, upper = 10 * range)$root, error = function(e) NA)
  }
  return(h_eff)
}

# Example usage:
range <- 1
nu <- 0.5   # exponential
matern_effective_range(range, nu)
# Try for nu=1.5 or other range values!
