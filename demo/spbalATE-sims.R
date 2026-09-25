### spbalATE-sims.R
### Dupont-style simulations with an unmeasured spatial confounder, in three
### settings: (i) raster (gridded) data, (ii) point data in a box, and
### (iii) point data in a non-rectangular polygon. For each, spBalance()
### balances a spatial smooth and spbalATE() estimates the ATE (true value 3)
### with iid, HAC and spatial block bootstrap standard errors.
###
### Run with: demo("spbalATE-sims", package = "spbalance")
### Needs the fields, spbootstrap, terra, and sf packages.

library(spbalance)
library(fields)
library(spbootstrap)
library(terra)
library(sf)

set.seed(2026)

# settings shared by all three cases
n.boot   <- 100                                  # bootstrap replicates
ps.form  <- exposure ~ s(x, y, k = 100)          # balance a spatial smooth
lambdas  <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 5, 10)
opt      <- list(tol = 1e-7, max.iter = 1000, alpha = 0.5, beta = 0.5)
par.fun  <- if (.Platform$OS.type == "unix") {  # evaluate replicates in parallel
  function(X, FUN) parallel::mclapply(X, FUN, mc.cores = max(1, parallel::detectCores() - 1))
} else lapply


#------------------------------------------------------------------------#
#--- Dupont-style data-generating process                             ---#
#------------------------------------------------------------------------#

# Two independent Matern Gaussian processes on [0, 10] x [0, 10], simulated
# on a fine grid by circulant embedding (using 'fields'), then evaluated
# anywhere by bilinear interpolation:
#   z  : unmeasured confounder (affects exposure and outcome)
#   z' : affects the outcome only
gp.grid <- list(x = seq(0, 10, length.out = 128), y = seq(0, 10, length.out = 128))
gp.sim  <- circulantEmbeddingSetup(
  gp.grid, M = c(512, 512),
  cov.args = list(Covariance = "Matern", aRange = 1, smoothness = 1)
)
simGP <- function() list(x = gp.grid$x, y = gp.grid$y, z = circulantEmbedding(gp.sim))

simDupont <- function(locs, beta = 3) {
  # exposure ~ Bernoulli(expit(2 * standardised z)); outcome = beta * exposure - z - z' + e
  z  <- interp.surface(simGP(), locs)
  zp <- interp.surface(simGP(), locs)
  exposure <- rbinom(nrow(locs), 1, plogis(2 * (z - mean(z)) / sd(z)))
  data.frame(x = locs[, 1], y = locs[, 2], exposure = exposure,
             outcome = beta * exposure - z - zp + rnorm(nrow(locs)))
}

fitATE <- function(df, boot) {
  fit <- spBalance(ps.form, data = df, lambda = lambdas, tuning = "coefvar", opt.params = opt)
  spbalATE(fit, df, outcome = "outcome", estimator = c("HT", "Hajek"),
           se = c("iid", "hac", "boot"), boot = c(boot, list(apply.fun = par.fun)))
}


#------------------------------------------------------------------------#
#--- (i) Raster data: a 30 x 30 grid                                  ---#
#------------------------------------------------------------------------#

r <- rast(nrows = 30, ncols = 30, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = "local")
sim.r <- simDupont(crds(r))
r$exposure <- sim.r$exposure
r$outcome  <- sim.r$outcome

# spBalance() works with data frames, so convert (one row per cell, with x, y)
grid.df <- as.data.frame(r, xy = TRUE)

# fitted spatial CBPS model using 'spBalance'
fit.grid <- spBalance(ps.form, data = grid.df, lambda = lambdas, tuning = "coefvar", opt.params = opt)

# The fitted object includes:
#   (i) fit$pi.hat: estimated propensity scores
#   (ii) fit$lambda: the lambda value selected by 'coefvar'
# among other outputs.

# The spBalance output can be used to estimate an ATE via the spbalATE() function.
# There are several inputs to be aware of:
#   outcome: you must specify the name of the outcome variable in your data frame
#   estimator: some combination of "HT", "Hajek", or "aIPTW". The latter requires
#         specification of a `mu.formula`
#   se: one or more of "iid", "hac", and "boot". The bootstrap estimator
#         is recommended, however note that it takes the most time to compute.
#   boot: If "boot" is included as an se argument, it requires additional
#         information. This includes
#             - n.boot: number of replicates
#             - block.l: block side length; if NULL (default) it is chosen
#                   using a procedure suggested by Nordman and Lahiri (2007)
#             - type: type of data, one of "grid" (data on a complete regular
#                   grid), "point" (irregular points contained in a rectangle
#                   of size b.x \by b.y), or "polygon" (points contained in
#                   a polygon boundary)
#             - apply.fun: type of 'apply' function used. default is 'lapply',
#                   but allows for parallelization using 'mclapply'.


# generate ATE estimates and standard errors
#   note: demo takes ~2 minutes using 7 cores on my local machine
ate.grid <- spbalATE(fit = fit.grid, data = grid.df, outcome = "outcome",
  estimator = c("HT", "Hajek"), se = c("iid", "hac", "boot"),
  boot = list(n.boot = n.boot, apply.fun = par.fun)
)

# resulting point estimates and
ate.grid



#------------------------------------------------------------------------#
#--- (ii) Point data in a box                                         ---#
#------------------------------------------------------------------------#

# Repeat the above procedure with irregular point locations, uniformly
# sampled within a [0,10]x[0,10] box.

box.df <- simDupont(cbind(runif(900, 0, 10), runif(900, 0, 10)))


# spatial CBPS
fit.box <- spBalance(ps.form, data = box.df, lambda = lambdas, tuning = "coefvar", opt.params = opt)

# estimate ATE
#   the arguments to boot change, since data are now irregularly spaced on a
#     rectangular domain. if block length is specified, it should be in
#     coordinate units.
#   demo takes ~2 minutes using 7 cores on my local machine
ate.box <- spbalATE(fit = fit.box, data = box.df, outcome = "outcome",
  estimator = c("HT", "Hajek"), se = c("iid", "hac", "boot"),
  boot = list(n.boot = n.boot, type = "point", box.x = c(0,10), box.y = c(0,10),
              apply.fun = par.fun)
)

ate.box


#------------------------------------------------------------------------#
#--- (iii) Point data in a non-rectangular polygon                    ---#
#------------------------------------------------------------------------#

study.area <- st_sfc(st_polygon(list(rbind(
  c(0.5, 1), c(6, 0.3), c(9.7, 3), c(8.5, 9.5), c(4, 8), c(1, 9.8), c(0.3, 5), c(0.5, 1)
))))
poly.pts <- st_coordinates(st_sample(study.area, 900))
poly.df  <- simDupont(poly.pts)

# points within a polygon
plot(study.area)
points(poly.pts, pch = 20, cex = 0.25)

# estimate spatial CBPS
fit.poly <- spBalance(ps.form, data = poly.df, lambda = lambdas, tuning = "coefvar", opt.params = opt)

# estimate ATE
#   once again, boot should be of type "polygon" and polygon boundary should be
#     supplied as an 'sf' object.
#   demo takes ~2 minutes using 7 cores on my local machine
ate.poly <- spbalATE(fit = fit.poly, data = poly.df, outcome = "outcome",
  estimator = c("HT", "Hajek"), se = c("iid", "hac", "boot"),
  boot = list(n.boot = n.boot, type = "polygon", boundary = study.area,
              apply.fun = par.fun)
)

ate.poly

#------------------------------------------------------------------------#
#--- Summary                                                          ---#
#------------------------------------------------------------------------#

# the three data settings, colored by exposure
op <- par(mfrow = c(1, 3), mar = c(3, 3, 2, 1))
plot(r$exposure, main = "(i) raster", legend = FALSE, col = c("grey85", "firebrick"))
plot(box.df$x, box.df$y, pch = 16, cex = 0.6, main = "(ii) box", asp = 1,
     col = c("grey60", "firebrick")[box.df$exposure + 1], xlab = "", ylab = "")
plot(study.area, main = "(iii) polygon", border = "grey30")
points(poly.df$x, poly.df$y, pch = 16, cex = 0.6,
       col = c("grey60", "firebrick")[poly.df$exposure + 1])
par(op)

# summary of results across sims
summ <- do.call(rbind, lapply(
  list(raster = ate.grid, box = ate.box, polygon = ate.poly),
  function(a) a$estimates
))
summ$setting <- rep(c("raster", "box", "polygon"), each = 2)
cat("True ATE = 3. One simulated dataset per setting, so these illustrate the\n",
    "interface rather than the estimators' properties.\n\n", sep = "")
print(summ[, c("setting", "estimator", "estimate", "se.iid", "se.hac", "se.boot")],
      digits = 3, row.names = FALSE)
cat(sprintf("\nRun times: raster %.0f s, box %.0f s, polygon %.0f s\n",
            t1[["elapsed"]], t2[["elapsed"]], t3[["elapsed"]]))



