# matern_cov_sphere <- function(u, c_hat, nu) {
#   x <- u / c_hat
#   # For numerical stability, cap extremely small x
#   x[x < 1e-8] <- 1e-8
#   factor <- (2^(nu-1)) / gamma(nu)
#   cov <- factor * (x^nu) * besselK(x, nu)
#   cov[u == 0] <- 1 # 保證變異數為1
#   return(cov)
# }
matern_cov_sphere <- function(u, c_hat, nu) {
  x <- u / c_hat
  # For numerical stability, cap extremely small x
  x[x < 1e-8] <- 1e-8
  factor <- (1/2^(nu-1)) / gamma(nu)
  cov <- factor * (x^nu) * besselK(x, nu)
  cov[u == 0] <- 1 # 保證變異數為1
  return(cov)
}

library(Matrix)

mle_negloglik_matern_X <- function(param, y, X, distm) {
  sigma2 <- param[1]
  phi    <- param[2]
  nugget <- param[3]
  nu     <- param[4]
  
  n <- length(y)
  
  # Covariance matrix
  K <- sigma2 * matern_cov_sphere(distm, phi, nu)
  diag(K) <- diag(K) + nugget
  
  chol_K <- try(chol(K), silent=TRUE)
  if (inherits(chol_K, "try-error")) return(1e10)
  
  # Generalized least squares estimate for beta:
  y_ <- as.numeric(y)
  X_ <- as.matrix(X)
  Kinv_y <- backsolve(chol_K, forwardsolve(t(chol_K), y_))
  Kinv_X <- backsolve(chol_K, forwardsolve(t(chol_K), X_))
  beta_hat <- solve(t(X_) %*% Kinv_X, t(X_) %*% Kinv_y)
  
  resid <- y_ - X_ %*% beta_hat
  Kinv_resid <- backsolve(chol_K, forwardsolve(t(chol_K), resid))
  logdetK <- 2 * sum(log(diag(chol_K)))
  
  ll <- 0.5 * (t(resid) %*% Kinv_resid + logdetK + n * log(2 * pi))
  return(as.numeric(ll))
}



fun = function(param) {
  vy_hat = param[1]
  c_hat = param[2]
  sigma2_hat = param[3]
  nu = param[4]   # 新增 smoothness

  cov_vu <- matern_cov_sphere(v$u, c_hat, nu)
  # 將 v$u 和 xdt 都套用 Matern covariance
  vr = (vy_hat + sigma2_hat) - vy_hat*cov_vu*
    (1-v$u/a)^6*(35/3*(v$u/a)^2+6*v$u/a+1)*as.numeric(v$u<=a)
  r1 = sum(v$n * ((v$v - vr) / vr)^2)

  vr2 = (vy_hat + sigma2_hat) -
    vy_hat*matern_cov_sphere(xdt, c_hat, nu)*
    (1-xdt/a)^6*(35/3*(xdt/a)^2+6*xdt/a+1)*as.numeric(xdt<=a)
  v2 = apply(res[,1,]*res[,2,], MARGIN = 1, FUN = sum)
  n2 = apply(res[, 2, ], MARGIN = 1, FUN = sum)
  fr = v2 / n2
  fr[1] = 0
  r2 = sum(n2 * ((fr - vr2) / vr2)^2)

  return(r1 + r2)
}

fun_exp_no_a = function(param)
{
  vy_hat = param[1] ; c_hat = param[2] ; a = a ; sigma2_hat = param[3]
  vr = (vy_hat+sigma2_hat)-vy_hat*exp(-v$u/c_hat)
  r1 = sum(v$n*((v$v-vr)/vr)^2)
  vr2 = (vy_hat+sigma2_hat)-vy_hat*exp(-xdt/c_hat)
  v2 = apply(res[,1,]*res[,2,], MARGIN = 1, FUN = sum)
  n2 = apply(res[,2,], MARGIN = 1, FUN = sum)
  fr = v2/n2 ; fr[1] = 0
  r2 = sum(n2*((fr-vr2)/vr2)^2)
  return(r1+r2)
}

fun_no_a = function(param) {
  vy_hat = param[1]
  c_hat = param[2]
  sigma2_hat = param[3]
  nu = param[4]   # 新增 smoothness
  
  cov_vu <- matern_cov_sphere(v$u, c_hat, nu)
  # 將 v$u 和 xdt 都套用 Matern covariance
  vr = (vy_hat + sigma2_hat) - vy_hat*cov_vu
  r1 = sum(v$n * ((v$v - vr) / vr)^2)
  
  vr2 = (vy_hat + sigma2_hat) -
    vy_hat*matern_cov_sphere(xdt, c_hat, nu)
  v2 = apply(res[,1,]*res[,2,], MARGIN = 1, FUN = sum)
  n2 = apply(res[,2,], MARGIN = 1, FUN = sum)
  fr = v2 / n2
  fr[1] = 0
  r2 = sum(n2 * ((fr - vr2) / vr2)^2)
  
  return(r1 + r2)
}

fun_a = function(param) {
  vy_hat = param[1]
  c_hat = param[2]
  sigma2_hat = param[3]
  nu = param[4]   # 新增 smoothness
  a = param[5]
  
  cov_vu <- matern_cov_sphere(v$u, c_hat, nu)
  # 將 v$u 和 xdt 都套用 Matern covariance
  vr = (vy_hat + sigma2_hat) - vy_hat*cov_vu*
    (1-v$u/a)^6*(35/3*(v$u/a)^2+6*v$u/a+1)*as.numeric(v$u<=a)
  r1 = sum(v$n * ((v$v - vr) / vr)^2)
  
  vr2 = (vy_hat + sigma2_hat) -
    vy_hat*matern_cov_sphere(xdt, c_hat, nu)*
    (1-xdt/a)^6*(35/3*(xdt/a)^2+6*xdt/a+1)*as.numeric(xdt<=a)
  v2 = apply(res[,1,]*res[,2,], MARGIN = 1, FUN = sum)
  n2 = apply(res[, 2, ], MARGIN = 1, FUN = sum)
  fr = v2 / n2
  fr[1] = 0
  r2 = sum(n2 * ((fr - vr2) / vr2)^2)
  
  return(r1 + r2)
}


# fun = function(param) {
#   vy_hat = param[1]
#   c_hat = param[2]
#   sigma2_hat = param[3]
#   nu = param[4]
#   
#   cov_vu <- matern_cov_sphere(v$u, c_hat, nu)
#   kernel <- (1-v$u/a)^6*(35/3*(v$u/a)^2+6*v$u/a+1)*as.numeric(v$u<=a)
#   vr = (vy_hat + sigma2_hat) - vy_hat * cov_vu * kernel
#   if(any(!is.finite(vr))) { print("vr not finite"); print(param); return(Inf) }
#   if(any(vr == 0)) { print("vr==0"); print(param); return(Inf) }
#   
#   r1 = sum(v$n * ((v$v - vr) / vr)^2)
#   
#   cov_xdt <- matern_cov_sphere(xdt, c_hat, nu)
#   kernel2 <- (1-xdt/a)^6*(35/3*(xdt/a)^2+6*xdt/a+1)*as.numeric(xdt<=a)
#   vr2 = (vy_hat + sigma2_hat) - vy_hat * cov_xdt * kernel2
#   if(any(!is.finite(vr2))) { print("vr2 not finite"); print(param); return(Inf) }
#   if(any(vr2 == 0)) { print("vr2==0"); print(param); return(Inf) }
#   
#   v2 = apply(res[, 1, ], MARGIN = 1, FUN = sum)
#   n2 = apply(res[, 2, ], MARGIN = 1, FUN = sum)
#   fr = v2 / n2
#   fr[1] = 0
#   r2 = sum(n2 * ((fr - vr2) / vr2)^2)
#   
#   val = r1 + r2
#   if(!is.finite(val)) { print("Objective not finite"); print(param); return(Inf) }
#   return(val)
# }