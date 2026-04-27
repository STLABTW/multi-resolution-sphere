# ==============================================================================
# fullmodel-max.R
# Multi-resolution model on the sphere — SST annual maximum dataset
#
# Depends on:  fn_pcc_test_pre.R   (C++ kernels, mrts_sphere, compiled via Rcpp)
#              fn_0610.R            (variogram objective functions)
#              effectivefn.R        (matern_eff_range)
#              data_sst_max_20240419.csv
#              integral_table2.rds
# ==============================================================================

# --- Model parameters (from variogram-based estimation) -----------------------
a          <- 0.1973   # Wendland compactness (angular, radians on unit sphere)
KK <- bb   <- 50       # number of MRTS basis functions used
c_hat      <- 0.0652   # Matérn range parameter
vy_hat     <- 1.6372   # partial sill
sigma2_hat <- 0.0152   # nugget
nu_hat     <- 0.5      # Matérn smoothness

nd         <- 10000    # training sample size

# --- Libraries ----------------------------------------------------------------
library(icosa)
library(fields)            # rdist.earth, Wendland, stationary.cov
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(pracma)
library(raster)
library(maps)
library(SparseM)
library(ggplot2)

# --- Source helpers (also compiles C++ via Rcpp) ------------------------------
source("fn_pcc_test_pre.R")
source("fn_0610.R")
source("effectivefn.R")

# --- Load data; remove land points --------------------------------------------
filename  <- "data_sst_max_20240419.csv"
y2        <- read.csv(filename)   # columns: latitude, longitude, temperature

pts_sf    <- st_as_sf(y2, coords = c("longitude", "latitude"), crs = 4326)
land_poly <- ne_countries(scale = "medium", returnclass = "sf")
is_land   <- lengths(st_intersects(pts_sf, land_poly)) > 0
y2        <- y2[!is_land, , drop = FALSE]

# --- Build icosahedral control-point grid ------------------------------------
gLow  <- hexagrid(tessellation = c(3, 3, 4))
y     <- y2$temperature
grids <- cbind(y2$latitude, y2$longitude)
cells <- locate(gLow, grids[, 2:1])
X     <- cart2sph(gLow@faceCenters[unique(cells), ])[, 2:1] / pi * 180
rownames(X) <- NULL
n     <- nrow(X)

# --- Kernel matrix K and spectral decomposition (bigK eigenpairs) -------------
bigK       <- 2000
onev       <- rep(1/n, n)
K          <- cpp_K(X[, 1], X[, 2], n)
Q          <- diag(1, n, n) - (1/n)
eiK        <- eigs_sym(Q %*% K %*% Q, bigK)
rm(Q)

eiKvecmval <- matrix(0, n, bigK)
for (i in 1:bigK) eiKvecmval[, i] <- eiK$vector[, i] / eiK$values[i]

# --- Pre-computed integral table ----------------------------------------------
tab <- readRDS("integral_table2.rds")

# --- Prediction grid (icosahedral face centres) -------------------------------
ggrids      <- gg <- cart2sph(gLow@faceCenters)[, 2:1] / pi * 180
rownames(ggrids) <- rownames(gg) <- NULL
NN          <- nrow(gg)

# --- Training sample ----------------------------------------------------------
set.seed(728)
pickn2  <- sample(seq_len(nrow(y2)), nd)
y_t     <- y2$temperature[pickn2]
grids_t <- cbind(y2$latitude[pickn2], y2$longitude[pickn2])
N_t     <- length(y_t)

# --- Design matrix for training locations ------------------------------------
dm_train <- cpp_Kmatrix6(bb, X, grids_t, K %*% onev, eiKvecmval, n, N_t,
                          tab$aa_grid, tab$tab_vals)

# --- Training covariance and inverse -----------------------------------------
cov_exp_hat <- vy_hat *
  Wendland(rdist.earth(grids_t[, 2:1], R = 1), theta = a, k = 2, dimension = 2) *
  stationary.cov(grids_t[, 2:1], theta = c_hat, Dist.args = list(R = 1),
                 Distance = "rdist.earth", Covariance = "Matern",
                 smoothness = nu_hat) +
  diag(sigma2_hat, N_t, N_t)

incov_exp_hat <- SparseM::solve(cov_exp_hat)

# --- GLS coefficient estimate ------------------------------------------------
infincov <- SparseM::solve(t(dm_train[, 1:bb]) %*% incov_exp_hat %*% dm_train[, 1:bb])
beta_hat <- infincov %*% t(dm_train[, 1:bb]) %*% incov_exp_hat %*% y_t

# --- Prediction at grid locations --------------------------------------------
Kmatrix   <- cpp_Kmatrix6(bigK, X, gg, K %*% onev, eiKvecmval, n, NN,
                           tab$aa_grid, tab$tab_vals)

cross_cov <- vy_hat *
  Wendland(rdist.earth(gg[, 2:1], grids_t[, 2:1], R = 1), theta = a, k = 2, dimension = 2) *
  stationary.cov(gg[, 2:1], grids_t[, 2:1], theta = c_hat, Dist.args = list(R = 1),
                 Distance = "rdist.earth", Covariance = "Matern",
                 smoothness = nu_hat)

rema     <- cross_cov %*% incov_exp_hat
khat_exp <- Kmatrix[, 1:KK] %*% beta_hat[1:KK] +
            rema %*% (y_t - dm_train[, 1:KK] %*% beta_hat[1:KK])

# --- Prediction standard error -----------------------------------------------
c2    <- SparseM::solve(t(dm_train[, 1:KK]) %*% incov_exp_hat %*% dm_train[, 1:KK])
f2    <- Kmatrix[, 1:KK] - rema %*% dm_train[, 1:KK]
resid <- numeric(NN)
for (i in 1:NN) {
  resid[i] <- t(f2[i, ]) %*% c2 %*% f2[i, ] -
              t(cross_cov[i, ]) %*% incov_exp_hat %*% cross_cov[i, ] + vy_hat
}
resid <- sqrt(resid)

# --- Clamp values and mask land points ----------------------------------------
ran       <- range(c(y, y_t))
psu_k_exp <- clamp(khat_exp[, 1], ran[1], ran[2])
residc    <- clamp(resid, 0, 2)
for (i in 1:NN) {
  if (!is.na(map.where("world", ggrids[i, 2], ggrids[i, 1])))
    psu_k_exp[i] <- residc[i] <- NA
}

# --- Plotting -----------------------------------------------------------------
dir_name  <- "realdata"
dir.create(dir_name, showWarnings = FALSE)
data_name <- tools::file_path_sans_ext(basename(filename))

world      <- ne_countries(scale = "medium", returnclass = "sf")
world_moll <- st_transform(world, crs = "+proj=moll")

make_sf_moll <- function(lon, lat, val) {
  df   <- data.frame(x = lon, y = lat, region = val)
  df_sf <- st_as_sf(df, coords = c("x", "y"), crs = 4326)
  st_transform(df_sf, crs = "+proj=moll")
}

base_map <- list(
  geom_sf(data = world_moll, fill = "gray95", color = "gray60", linewidth = 0.2),
  geom_sf(data = world,      fill = NA,       color = "gray50", linewidth = 0.1),
  coord_sf(crs = "+proj=moll", expand = FALSE),
  labs(color = ""),
  theme_void(),
  theme(legend.text        = element_text(size = 10),
        legend.title       = element_text(size = 10),
        legend.key.height  = unit(1.5, "cm"),
        legend.key.width   = unit(0.4, "cm"))
)

# Prediction plot
df_moll  <- make_sf_moll(ggrids[, 2], ggrids[, 1], psu_k_exp)
plotfile <- file.path(dir_name, sprintf("pred_real_exp_all_%s.png", data_name))
png(plotfile, width = 1000, height = 500, res = 150)
print(ggplot() +
  base_map +
  geom_sf(data = df_moll, aes(color = region), size = 1.2, shape = 16) +
  scale_color_gradientn(
    colours  = c("#00AAAA", "#FFFFBB", "#FF3333"),
    limits   = ran,
    name     = expression(""*degree*"C"),
    na.value = "transparent"))
dev.off()
cat("Saved:", plotfile, "\n")

# Prediction standard error plot
df_moll  <- make_sf_moll(ggrids[, 2], ggrids[, 1], residc)
plotfile <- file.path(dir_name, sprintf("pred_real_resid_tp_all_%s.png", data_name))
png(plotfile, width = 1000, height = 500, res = 150)
print(ggplot() +
  base_map +
  geom_sf(data = df_moll, aes(color = region), size = 1.1, shape = 16) +
  scale_color_gradientn(
    colours  = c("#00AAAA", "#FFFFBB", "#FF3333"),
    limits   = c(0, 2.1),
    na.value = "transparent"))
dev.off()
cat("Saved:", plotfile, "\n")

# Raw data plot
df_moll  <- make_sf_moll(y2$longitude, y2$latitude, y2$temperature)
plotfile <- file.path(dir_name, sprintf("realdata_%s.png", data_name))
png(plotfile, width = 1000, height = 500, res = 150)
print(ggplot() +
  geom_sf(data = world_moll, fill = "gray90", color = "black", linewidth = 0.2) +
  geom_sf(data = world,      fill = NA,       color = "black", linewidth = 0.1) +
  coord_sf(crs = "+proj=moll", expand = FALSE) +
  geom_sf(data = df_moll, aes(color = region), size = 0.5, shape = 16) +
  scale_color_gradientn(
    colours  = c("#00AAAA", "#FFFFBB", "#FF3333"),
    limits   = ran,
    name     = expression(""*degree*"C"),
    na.value = "transparent") +
  labs(color = "") +
  theme_void() +
  theme(legend.text        = element_text(size = 10),
        legend.title       = element_text(size = 10),
        legend.key.height  = unit(1.5, "cm"),
        legend.key.width   = unit(0.4, "cm")))
dev.off()
cat("Saved:", plotfile, "\n")
